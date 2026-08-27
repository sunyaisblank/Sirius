// Ray-bundle wiring in the trace kernel (specification P2, Vulkan side).
//
// trace.slang propagates two geodesic-deviation vectors through its render loop
// when the beam flag (params[43]) is set, using gr_deviation's full-Riemann
// Jacobi acceleration (DNGR eq A.19), and carries the beam's transverse
// expansion in the radiance alpha channel. This suite runs on Lavapipe and pins:
//   - beams off leaves alpha == 1 everywhere (the default render is unmoved);
//   - beams on yields a finite area-expansion field that grows toward the
//     shadow, the lensing signature the CPU bundle also shows.
//
// The kernel and CPU use independent live integrators, so this is behavioural
// parity rather than a bitwise claim. The exact radial/circular congruence
// tolerance belongs to the double-precision covariant oracle; the two live
// paths are additionally joined at the ellipse-filtered point-catalogue output.

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
using sirius::backend::KernelHandle;
using sirius::backend::ResolveVulkanDeviceIndex;

std::vector<std::uint32_t> LoadSpirv(const std::string& path) {
    std::ifstream file(path, std::ios::binary | std::ios::ate);
    if (!file) return {};
    const auto size = static_cast<std::size_t>(file.tellg());
    std::vector<std::uint32_t> words(size / sizeof(std::uint32_t));
    file.seekg(0);
    file.read(reinterpret_cast<char*>(words.data()), static_cast<std::streamsize>(size));
    return words;
}

// Kerr a=0.9 down the spin axis from 40M, single tile, gradient background.
std::vector<float> BaseParams(std::uint32_t w, std::uint32_t h) {
    std::vector<float> params(72, 0.0f);
    params[46] = 0.5f;
    params[47] = 0.5f;
    params[53] = 1.0f;
    params[0] = float(w);
    params[1] = float(h);
    params[2] = 0.0f;  // Kerr-Schild family.
    params[3] = 1.0f;  // M
    params[4] = 0.9f;  // a
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
    params[18] = 0.0f;      // up
    params[19] = 0.6f;      // fov
    params[20] = 1.0f;      // aspect
    params[21] = 2000.0f;   // max_steps
    params[22] = 0.06f;     // stepScale
    params[23] = 0.02f;     // min_step
    params[24] = 2.0f;      // max_step
    params[25] = 100.0f;    // escape_radius
    params[26] = 1.12f;     // captureFactor
    params[27] = 0.0f;      // disk disabled
    params[34] = float(w);  // tileWidth
    params[35] = float(h);  // tileHeight
    params[37] = 1.0f;      // dummy starfield dims
    params[38] = 1.0f;
    return params;
}

std::vector<float> Dispatch(ComputeDevice& device, KernelHandle kernel,
                            const std::vector<float>& params, std::uint32_t w, std::uint32_t h) {
    std::vector<float> radiance(w * h * 4, 0.0f);
    const std::vector<std::uint32_t> star_dummy = {0u};
    const auto rbuf = device.CreateBuffer(radiance.size() * sizeof(float), BufferUsage::kStorage);
    const auto pbuf = device.CreateBuffer(params.size() * sizeof(float), BufferUsage::kStorage);
    const auto sbuf =
        device.CreateBuffer(star_dummy.size() * sizeof(std::uint32_t), BufferUsage::kStorage);
    const auto psbuf =
        device.CreateBuffer(star_dummy.size() * sizeof(std::uint32_t), BufferUsage::kStorage);
    const auto pobuf =
        device.CreateBuffer(star_dummy.size() * sizeof(std::uint32_t), BufferUsage::kStorage);
    const auto pibuf =
        device.CreateBuffer(star_dummy.size() * sizeof(std::uint32_t), BufferUsage::kStorage);
    EXPECT_TRUE(rbuf && pbuf && sbuf && psbuf && pobuf && pibuf);
    EXPECT_TRUE(device.WriteBuffer(*rbuf, std::as_bytes(std::span<const float>(radiance))));
    EXPECT_TRUE(device.WriteBuffer(*pbuf, std::as_bytes(std::span<const float>(params))));
    EXPECT_TRUE(
        device.WriteBuffer(*sbuf, std::as_bytes(std::span<const std::uint32_t>(star_dummy))));
    const BufferHandle binding[] = {*rbuf, *pbuf, *sbuf, *psbuf, *pobuf, *pibuf};
    EXPECT_TRUE(device.Dispatch(kernel, binding, (w + 7) / 8, (h + 7) / 8, 1).has_value());
    EXPECT_TRUE(device.ReadBuffer(*rbuf, std::as_writable_bytes(std::span<float>(radiance))));
    return radiance;
}

TEST(KernelBeam, BeamFlagWiresDeviationWithoutMovingDefault) {
#ifndef SIRIUS_KERNEL_DIR
    GTEST_SKIP() << "kernels not compiled (slangc absent at configure time)";
#else
    const auto devices = EnumerateVulkanDevices();
    ASSERT_TRUE(devices.has_value()) << devices.error().Description();
    if (devices->empty()) GTEST_SKIP() << "no Vulkan device present";
    const auto selected = ResolveVulkanDeviceIndex(*devices);
    ASSERT_TRUE(selected.has_value()) << selected.error().Description();
    auto device = CreateVulkanDevice(*selected);
    ASSERT_TRUE(device.has_value()) << device.error().Description();
    const auto spirv = LoadSpirv(std::string(SIRIUS_KERNEL_DIR) + "/trace.spv");
    ASSERT_FALSE(spirv.empty()) << "trace.spv missing";
    const auto kernel = (*device)->LoadKernel(spirv);
    ASSERT_TRUE(kernel.has_value()) << kernel.error().Description();

    constexpr std::uint32_t kW = 64, kH = 64;

    // Beams off: the default render, alpha must be exactly 1 everywhere.
    std::vector<float> off = Dispatch(**device, *kernel, BaseParams(kW, kH), kW, kH);
    std::vector<float> off_rgb = off;
    for (std::uint32_t p = 0; p < kW * kH; ++p) {
        ASSERT_FLOAT_EQ(off[p * 4 + 3], 1.0f) << "beams-off alpha moved at pixel " << p;
    }

    // Beams on.
    std::vector<float> p_on = BaseParams(kW, kH);
    p_on[43] = 1.0f;
    std::vector<float> on = Dispatch(**device, *kernel, p_on, kW, kH);

    // RGB is untouched by the beam wiring; only alpha changes.
    for (std::uint32_t p = 0; p < kW * kH; ++p) {
        for (int c = 0; c < 3; ++c) {
            ASSERT_FLOAT_EQ(on[p * 4 + c], off_rgb[p * 4 + c])
                << "beam wiring perturbed RGB at pixel " << p << " channel " << c;
        }
    }

    // Beam expansion (alpha) is finite; measure it against proximity to the
    // shadow. Escaped pixels near the shadow are strongly lensed (large
    // expansion); edge pixels sit near the seed value.
    float center_max = 0.0f;  // Near the shadow (image centre).
    float edge_expansion = 0.0f;
    int edge_count = 0;
    for (std::uint32_t y = 0; y < kH; ++y) {
        for (std::uint32_t x = 0; x < kW; ++x) {
            std::uint32_t p = y * kW + x;
            float alpha = on[p * 4 + 3];
            ASSERT_TRUE(std::isfinite(alpha)) << "non-finite beam expansion";
            bool shadow = on[p * 4 + 0] < 1e-4f && on[p * 4 + 1] < 1e-4f && on[p * 4 + 2] < 1e-4f;
            float dx = float(x) - kW / 2.0f, dy = float(y) - kH / 2.0f;
            float rad = std::sqrt(dx * dx + dy * dy);
            if (!shadow && rad < 18.0f) center_max = std::max(center_max, alpha);
            if (rad > 28.0f) {
                edge_expansion += alpha;
                edge_count++;
            }
        }
    }
    float edge_mean = (edge_count > 0) ? edge_expansion / edge_count : 0.0f;
    std::cout << "[kernel beam] near-shadow max expansion=" << center_max
              << " edge mean expansion=" << edge_mean << "\n";

    // The bundle expands toward the shadow far more than at the frame edge - the
    // lensing signature, consistent with the CPU bundle's defocusing.
    EXPECT_GT(center_max, 2.0f * std::max(edge_mean, 1e-6f));
    EXPECT_GT(center_max, 1.0f);
#endif
}

}  // namespace
