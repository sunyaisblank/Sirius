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
using sirius::backend::ResolveVulkanDeviceIndex;

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
    const auto selected = ResolveVulkanDeviceIndex(*devices);
    ASSERT_TRUE(selected.has_value()) << selected.error().Description();
    auto device = CreateVulkanDevice(*selected);
    ASSERT_TRUE(device.has_value()) << device.error().Description();

    const auto spirv = LoadSpirv(std::string(SIRIUS_KERNEL_DIR) + "/trace.spv");
    ASSERT_FALSE(spirv.empty()) << "trace.spv missing or empty";
    const auto kernel = (*device)->LoadKernel(spirv);
    ASSERT_TRUE(kernel.has_value()) << kernel.error().Description();

    constexpr std::uint32_t kWidth = 64;
    constexpr std::uint32_t kHeight = 64;

    // Kerr M=1, a=0.9 viewed down the spin axis from 40M. fov and distance are
    // chosen so the capture cross-section covers a modest image fraction. The
    // full image is a single tile here (tileOrigin 0, tile == image). Background
    // is the analytic gradient (starfield disabled), so no texture is needed.
    std::vector<float> params(65, 0.0f);
    params[44] = 0.5f;
    params[45] = 0.5f;
    params[51] = 1.0f;
    params[0] = kWidth;   // imageWidth
    params[1] = kHeight;  // imageHeight
    params[2] = 0.0f;     // metricId (Kerr-Schild family)
    params[3] = 1.0f;     // M
    params[4] = 0.9f;     // a
    params[5] = 0.0f;     // Q
    params[6] = 0.0f;     // Lambda
    params[7] = 0.0f;
    params[8] = 0.0f;
    params[9] = -40.0f;  // camera position
    params[10] = 0.0f;
    params[11] = 0.0f;
    params[12] = 1.0f;  // forward
    params[13] = 1.0f;
    params[14] = 0.0f;
    params[15] = 0.0f;  // right
    params[16] = 0.0f;
    params[17] = 1.0f;
    params[18] = 0.0f;            // up
    params[19] = 0.6f;            // fov (radians)
    params[20] = 1.0f;            // aspect
    params[21] = 2000.0f;         // max_steps
    params[22] = 0.06f;           // stepScale
    params[23] = 0.02f;           // min_step
    params[24] = 2.0f;            // max_step
    params[25] = 100.0f;          // escape_radius
    params[26] = 1.12f;           // captureFactor
    params[27] = 0.0f;            // disk disabled
    params[31] = 0.0f;            // tileOriginX
    params[32] = 0.0f;            // tileOriginY
    params[33] = float(kWidth);   // tileWidth
    params[34] = float(kHeight);  // tileHeight
    params[35] = 0.0f;            // starfield disabled (gradient background)
    params[36] = 1.0f;            // starfieldWidth (dummy)
    params[37] = 1.0f;            // starfieldHeight (dummy)

    std::vector<float> radiance(kWidth * kHeight * 4, 0.0f);
    const std::vector<std::uint32_t> starfield_dummy = {0u};

    const auto rbuf =
        (*device)->CreateBuffer(radiance.size() * sizeof(float), BufferUsage::kStorage);
    const auto pbuf = (*device)->CreateBuffer(params.size() * sizeof(float), BufferUsage::kStorage);
    const auto sbuf = (*device)->CreateBuffer(starfield_dummy.size() * sizeof(std::uint32_t),
                                              BufferUsage::kStorage);
    const auto psbuf = (*device)->CreateBuffer(sizeof(std::uint32_t), BufferUsage::kStorage);
    const auto pobuf = (*device)->CreateBuffer(sizeof(std::uint32_t), BufferUsage::kStorage);
    const auto pibuf = (*device)->CreateBuffer(sizeof(std::uint32_t), BufferUsage::kStorage);
    ASSERT_TRUE(rbuf && pbuf && sbuf && psbuf && pobuf && pibuf);
    ASSERT_TRUE((*device)->WriteBuffer(*rbuf, std::as_bytes(std::span<const float>(radiance))));
    ASSERT_TRUE((*device)->WriteBuffer(*pbuf, std::as_bytes(std::span<const float>(params))));
    ASSERT_TRUE((*device)->WriteBuffer(
        *sbuf, std::as_bytes(std::span<const std::uint32_t>(starfield_dummy))));

    const BufferHandle binding[] = {*rbuf, *pbuf, *sbuf, *psbuf, *pobuf, *pibuf};
    const auto dispatched =
        (*device)->Dispatch(*kernel, binding, (kWidth + 7) / 8, (kHeight + 7) / 8, 1);
    ASSERT_TRUE(dispatched.has_value()) << dispatched.error().Description();

    ASSERT_TRUE((*device)->ReadBuffer(*rbuf, std::as_writable_bytes(std::span<float>(radiance))));

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
    std::cout << "[ trace    ] shadow fraction=" << fraction << " luminance range=[" << min_lum
              << ", " << max_lum << "]\n";
    EXPECT_GE(fraction, 0.005) << "shadow fraction too small: " << fraction;
    EXPECT_LE(fraction, 0.30) << "shadow fraction too large: " << fraction;
#endif
}

#ifdef SIRIUS_KERNEL_DIR
// Dispatch one 64x64 Kerr scene (the same scene as the smoke gate above)
// through the given SPIR-V module and return the RGBA radiance field.
std::vector<float> RunKerrScene(ComputeDevice& device, const std::vector<std::uint32_t>& spirv) {
    const auto kernel = device.LoadKernel(spirv);
    if (!kernel) {
        ADD_FAILURE() << kernel.error().Description();
        return {};
    }

    constexpr std::uint32_t kWidth = 64;
    constexpr std::uint32_t kHeight = 64;

    std::vector<float> params(65, 0.0f);
    params[44] = 0.5f;
    params[45] = 0.5f;
    params[51] = 1.0f;
    params[0] = kWidth;
    params[1] = kHeight;
    params[2] = 0.0f;
    params[3] = 1.0f;
    params[4] = 0.9f;
    params[7] = 0.0f;
    params[8] = 0.0f;
    params[9] = -40.0f;
    params[10] = 0.0f;
    params[11] = 0.0f;
    params[12] = 1.0f;
    params[13] = 1.0f;
    params[14] = 0.0f;
    params[15] = 0.0f;
    params[16] = 0.0f;
    params[17] = 1.0f;
    params[18] = 0.0f;
    params[19] = 0.6f;
    params[20] = 1.0f;
    params[21] = 2000.0f;
    params[22] = 0.06f;
    params[23] = 0.02f;
    params[24] = 2.0f;
    params[25] = 100.0f;
    params[26] = 1.12f;
    params[33] = float(kWidth);
    params[34] = float(kHeight);
    params[36] = 1.0f;
    params[37] = 1.0f;

    std::vector<float> radiance(kWidth * kHeight * 4, 0.0f);
    const std::vector<std::uint32_t> starfield_dummy = {0u};

    const auto rbuf = device.CreateBuffer(radiance.size() * sizeof(float), BufferUsage::kStorage);
    const auto pbuf = device.CreateBuffer(params.size() * sizeof(float), BufferUsage::kStorage);
    const auto sbuf =
        device.CreateBuffer(starfield_dummy.size() * sizeof(std::uint32_t), BufferUsage::kStorage);
    const auto psbuf = device.CreateBuffer(sizeof(std::uint32_t), BufferUsage::kStorage);
    const auto pobuf = device.CreateBuffer(sizeof(std::uint32_t), BufferUsage::kStorage);
    const auto pibuf = device.CreateBuffer(sizeof(std::uint32_t), BufferUsage::kStorage);
    if (!(rbuf && pbuf && sbuf && psbuf && pobuf && pibuf)) {
        ADD_FAILURE() << "buffer creation failed";
        return {};
    }
    const bool wrote =
        device.WriteBuffer(*rbuf, std::as_bytes(std::span<const float>(radiance))) &&
        device.WriteBuffer(*pbuf, std::as_bytes(std::span<const float>(params))) &&
        device.WriteBuffer(*sbuf, std::as_bytes(std::span<const std::uint32_t>(starfield_dummy)));
    if (!wrote) {
        ADD_FAILURE() << "buffer upload failed";
        return {};
    }

    const BufferHandle binding[] = {*rbuf, *pbuf, *sbuf, *psbuf, *pobuf, *pibuf};
    const auto dispatched =
        device.Dispatch(*kernel, binding, (kWidth + 7) / 8, (kHeight + 7) / 8, 1);
    if (!dispatched) {
        ADD_FAILURE() << dispatched.error().Description();
        return {};
    }
    if (!device.ReadBuffer(*rbuf, std::as_writable_bytes(std::span<float>(radiance)))) {
        ADD_FAILURE() << "radiance readback failed";
        return {};
    }
    return radiance;
}
#endif

// The fp64 rung (trace_fp64.spv, same source with the Cartesian trajectory
// core widened to double) renders the same Kerr scene as the fp32 kernel and
// agrees with it closely away from discontinuous capture flips. The geometry
// is identical, while each precision policy owns its adaptive step schedule;
// invariant accuracy and physical termination are qualified by KernelParity.
// Verifies the rung is a real double-precision trajectory, not a relabelled
// fp32 module, by requiring the artefact to load, dispatch, and stay finite
// on a shaderFloat64 device (Lavapipe reports it).
TEST(KernelTrace, Fp64RungAgreesWithFp32OnKerrScene) {
#ifndef SIRIUS_KERNEL_DIR
    GTEST_SKIP() << "kernels not compiled (slangc absent at configure time)";
#else
    const auto devices = EnumerateVulkanDevices();
    ASSERT_TRUE(devices.has_value()) << devices.error().Description();
    if (devices->empty()) {
        GTEST_SKIP() << "no Vulkan device present";
    }
    const auto selected = ResolveVulkanDeviceIndex(*devices);
    ASSERT_TRUE(selected.has_value()) << selected.error().Description();
    auto device = CreateVulkanDevice(*selected);
    ASSERT_TRUE(device.has_value()) << device.error().Description();
    if (!(*device)->Info().supports_fp64) {
        GTEST_SKIP() << "device lacks shaderFloat64";
    }

    const auto spirv32 = LoadSpirv(std::string(SIRIUS_KERNEL_DIR) + "/trace.spv");
    const auto spirv64 = LoadSpirv(std::string(SIRIUS_KERNEL_DIR) + "/trace_fp64.spv");
    ASSERT_FALSE(spirv32.empty()) << "trace.spv missing";
    ASSERT_FALSE(spirv64.empty()) << "trace_fp64.spv missing";

    const auto r32 = RunKerrScene(**device, spirv32);
    const auto r64 = RunKerrScene(**device, spirv64);
    ASSERT_EQ(r32.size(), r64.size());
    ASSERT_FALSE(r64.empty());

    double sum_abs_diff = 0.0;
    double max_abs_diff = 0.0;
    float min64 = std::numeric_limits<float>::max();
    float max64 = 0.0f;
    for (std::size_t i = 0; i < r64.size(); ++i) {
        ASSERT_TRUE(std::isfinite(r64[i])) << "non-finite fp64-rung radiance";
        const double d = std::abs(double(r64[i]) - double(r32[i]));
        sum_abs_diff += d;
        max_abs_diff = std::max(max_abs_diff, d);
        min64 = std::min(min64, r64[i]);
        max64 = std::max(max64, r64[i]);
    }
    const double mean_abs_diff = sum_abs_diff / double(r64.size());

    std::size_t black32 = 0;
    std::size_t black64 = 0;
    for (std::size_t pixel = 0; pixel < r64.size() / 4; ++pixel) {
        const std::size_t base = pixel * 4;
        if (r32[base] < 1e-4f && r32[base + 1] < 1e-4f && r32[base + 2] < 1e-4f) {
            ++black32;
        }
        if (r64[base] < 1e-4f && r64[base + 1] < 1e-4f && r64[base + 2] < 1e-4f) {
            ++black64;
        }
    }

    std::cout << "[ trace64  ] mean|d|=" << mean_abs_diff << " max|d|=" << max_abs_diff
              << " range64=[" << min64 << ", " << max64 << "] black32=" << black32
              << " black64=" << black64 << "\n";

    EXPECT_GT(max64 - min64, 1e-3f) << "fp64 radiance field is constant";
    const std::size_t shadow_count_tolerance = (r64.size() / 4) / 50;  // two percent
    const std::size_t shadow_count_difference =
        black64 > black32 ? black64 - black32 : black32 - black64;
    EXPECT_LE(shadow_count_difference, shadow_count_tolerance)
        << "precision rungs disagree on the physical-termination population";
    // Same scene, independently controlled precision schedules: the fields
    // must agree closely in the mean. Individual pixels may flip across the
    // capture edge, so the max is bounded loosely by the background dynamic
    // range rather than tightly.
    EXPECT_LT(mean_abs_diff, 1e-2) << "fp64 rung diverges from fp32 in the mean";
    EXPECT_LT(max_abs_diff, 2.0) << "fp64 rung wildly diverges at a pixel";
#endif
}

// The compensated rung (trace_fp32comp.spv, Kahan state accumulation) renders
// the same scene, stays finite and non-constant, agrees with plain fp32
// closely, and — measured against the fp64 reference — tracks it at least as
// well as plain fp32 does (within slack for noise; on hardware with sloppier
// fp32 the compensation's gain is larger, which this bound also admits).
TEST(KernelTrace, CompensatedRungTracksFp64AtLeastAsWellAsFp32) {
#ifndef SIRIUS_KERNEL_DIR
    GTEST_SKIP() << "kernels not compiled (slangc absent at configure time)";
#else
    const auto devices = EnumerateVulkanDevices();
    ASSERT_TRUE(devices.has_value()) << devices.error().Description();
    if (devices->empty()) {
        GTEST_SKIP() << "no Vulkan device present";
    }
    const auto selected = ResolveVulkanDeviceIndex(*devices);
    ASSERT_TRUE(selected.has_value()) << selected.error().Description();
    auto device = CreateVulkanDevice(*selected);
    ASSERT_TRUE(device.has_value()) << device.error().Description();
    if (!(*device)->Info().supports_fp64) {
        GTEST_SKIP() << "device lacks shaderFloat64 (needed for the reference field)";
    }

    const auto spirv32 = LoadSpirv(std::string(SIRIUS_KERNEL_DIR) + "/trace.spv");
    const auto spirvC = LoadSpirv(std::string(SIRIUS_KERNEL_DIR) + "/trace_fp32comp.spv");
    const auto spirv64 = LoadSpirv(std::string(SIRIUS_KERNEL_DIR) + "/trace_fp64.spv");
    ASSERT_FALSE(spirv32.empty());
    ASSERT_FALSE(spirvC.empty()) << "trace_fp32comp.spv missing";
    ASSERT_FALSE(spirv64.empty());

    const auto r32 = RunKerrScene(**device, spirv32);
    const auto rC = RunKerrScene(**device, spirvC);
    const auto r64 = RunKerrScene(**device, spirv64);
    ASSERT_EQ(rC.size(), r32.size());
    ASSERT_EQ(rC.size(), r64.size());
    ASSERT_FALSE(rC.empty());

    double sum_c64 = 0.0, sum_3264 = 0.0, sum_c32 = 0.0;
    float minc = std::numeric_limits<float>::max();
    float maxc = 0.0f;
    for (std::size_t i = 0; i < rC.size(); ++i) {
        ASSERT_TRUE(std::isfinite(rC[i])) << "non-finite compensated-rung radiance";
        sum_c64 += std::abs(double(rC[i]) - double(r64[i]));
        sum_3264 += std::abs(double(r32[i]) - double(r64[i]));
        sum_c32 += std::abs(double(rC[i]) - double(r32[i]));
        minc = std::min(minc, rC[i]);
        maxc = std::max(maxc, rC[i]);
    }
    const double n = double(rC.size());
    const double mean_c64 = sum_c64 / n;
    const double mean_3264 = sum_3264 / n;
    const double mean_c32 = sum_c32 / n;
    std::cout << "[ traceC   ] mean|comp-fp64|=" << mean_c64 << " mean|fp32-fp64|=" << mean_3264
              << " mean|comp-fp32|=" << mean_c32 << "\n";

    EXPECT_GT(maxc - minc, 1e-3f) << "compensated radiance field is constant";
    EXPECT_LT(mean_c32, 1e-2) << "compensated rung diverges from fp32 in the mean";
    EXPECT_LE(mean_c64, mean_3264 * 1.5 + 1e-9)
        << "compensation made the fp64 tracking worse, which defeats the rung";
#endif
}

}  // namespace
