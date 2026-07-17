// Full-image smoke gate for the trace kernel (workstream deliverable 3): a
// 64x64 Kerr render with no disk, background shaded by escape direction, must
// dispatch on Lavapipe and produce a finite, non-constant field whose
// horizon-shadow fraction sits in a broad sane band. This is a behavioural
// smoke, not a pixel-parity gate; the physics it exercises (Kerr-Schild metric,
// Christoffels, Cartesian RK4) is pinned by kernel_parity_test.cpp.

#include "sirius/backend/device.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <span>
#include <string>
#include <vector>

namespace {

using sirius::backend::BufferHandle;
using sirius::backend::BufferUsage;
using sirius::backend::ComputeDevice;
using sirius::backend::CreateVulkanDevice;
using sirius::backend::EnumerateVulkanDevices;

std::vector<std::uint32_t> LoadSpirv(const std::string& path) {
    std::ifstream file(path, std::ios::binary | std::ios::ate);
    if (!file) {
        return {};
    }
    const auto size = static_cast<std::size_t>(file.tellg());
    std::vector<std::uint32_t> words(size / sizeof(std::uint32_t));
    file.seekg(0);
    file.read(reinterpret_cast<char*>(words.data()), static_cast<std::streamsize>(size));
    return words;
}

TEST(KernelTrace, KerrRenderIsFiniteNonConstantWithBoundedShadow) {
#ifndef SIRIUS_KERNEL_DIR
    GTEST_SKIP() << "kernels not compiled (slangc absent at configure time)";
#else
    const auto devices = EnumerateVulkanDevices();
    ASSERT_TRUE(devices.has_value()) << devices.error().Description();
    if (devices->empty()) {
        GTEST_SKIP() << "no Vulkan device present";
    }
    auto device = CreateVulkanDevice(0);
    ASSERT_TRUE(device.has_value()) << device.error().Description();

    const auto spirv = LoadSpirv(std::string(SIRIUS_KERNEL_DIR) + "/trace.spv");
    ASSERT_FALSE(spirv.empty()) << "trace.spv missing or empty";
    const auto kernel = (*device)->LoadKernel(spirv);
    ASSERT_TRUE(kernel.has_value()) << kernel.error().Description();

    constexpr std::uint32_t kWidth = 64;
    constexpr std::uint32_t kHeight = 64;

    // Kerr M=1, a=0.9 viewed down the spin axis from 40M. fov and distance are
    // chosen so the capture cross-section covers a modest image fraction.
    std::vector<float> params(32, 0.0f);
    params[0] = kWidth;
    params[1] = kHeight;
    params[2] = 0.0f;    // metricId (Kerr-Schild family)
    params[3] = 1.0f;    // M
    params[4] = 0.9f;    // a
    params[5] = 0.0f;    // Q
    params[6] = 0.0f;    // Lambda
    params[7] = 0.0f; params[8] = 0.0f; params[9] = -40.0f;   // camera position
    params[10] = 0.0f; params[11] = 0.0f; params[12] = 1.0f;  // forward
    params[13] = 1.0f; params[14] = 0.0f; params[15] = 0.0f;  // right
    params[16] = 0.0f; params[17] = 1.0f; params[18] = 0.0f;  // up
    params[19] = 0.6f;   // fov (radians)
    params[20] = 1.0f;   // aspect
    params[21] = 2000.0f;  // maxSteps
    params[22] = 0.06f;    // stepScale
    params[23] = 0.02f;    // minStep
    params[24] = 2.0f;     // maxStep
    params[25] = 100.0f;   // escapeRadius
    params[26] = 1.12f;    // captureFactor
    params[27] = 0.0f;     // disk disabled

    std::vector<float> radiance(kWidth * kHeight * 4, 0.0f);

    const auto rbuf =
        (*device)->CreateBuffer(radiance.size() * sizeof(float), BufferUsage::kStorage);
    const auto pbuf =
        (*device)->CreateBuffer(params.size() * sizeof(float), BufferUsage::kStorage);
    ASSERT_TRUE(rbuf && pbuf);
    ASSERT_TRUE((*device)->WriteBuffer(*rbuf, std::as_bytes(std::span<const float>(radiance))));
    ASSERT_TRUE((*device)->WriteBuffer(*pbuf, std::as_bytes(std::span<const float>(params))));

    const BufferHandle binding[] = {*rbuf, *pbuf};
    const auto dispatched =
        (*device)->Dispatch(*kernel, binding, (kWidth + 7) / 8, (kHeight + 7) / 8, 1);
    ASSERT_TRUE(dispatched.has_value()) << dispatched.error().Description();

    ASSERT_TRUE(
        (*device)->ReadBuffer(*rbuf, std::as_writable_bytes(std::span<float>(radiance))));

    // Finiteness across every channel.
    for (float value : radiance) {
        ASSERT_TRUE(std::isfinite(value)) << "non-finite radiance sample";
    }

    // Non-constancy and shadow fraction from per-pixel luminance.
    float min_lum = std::numeric_limits<float>::max();
    float max_lum = 0.0f;
    std::size_t shadow_pixels = 0;
    for (std::uint32_t p = 0; p < kWidth * kHeight; ++p) {
        const float r = radiance[p * 4 + 0];
        const float g = radiance[p * 4 + 1];
        const float b = radiance[p * 4 + 2];
        const float lum = 0.2126f * r + 0.7152f * g + 0.0722f * b;
        min_lum = std::min(min_lum, lum);
        max_lum = std::max(max_lum, lum);
        if (r < 1e-4f && g < 1e-4f && b < 1e-4f) {
            ++shadow_pixels;
        }
    }
    EXPECT_GT(max_lum - min_lum, 1e-3f) << "radiance field is constant";

    const double fraction =
        static_cast<double>(shadow_pixels) / static_cast<double>(kWidth * kHeight);
    std::cout << "[ trace    ] shadow fraction=" << fraction << " luminance range=["
              << min_lum << ", " << max_lum << "]\n";
    EXPECT_GE(fraction, 0.005) << "shadow fraction too small: " << fraction;
    EXPECT_LE(fraction, 0.30) << "shadow fraction too large: " << fraction;
#endif
}

}  // namespace
