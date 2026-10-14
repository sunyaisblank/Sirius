// Postconditions verified here: device enumeration reports rather than
// throws; a Slang-compiled kernel dispatched through the Vulkan adapter
// reproduces the CPU computation of f(r) = 1 - 2M/r within fp32 rounding.
// On machines with no Vulkan ICD both tests skip, cleanly and loudly.

#include "sirius/backend/device.h"

#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <vector>

namespace {

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

TEST(VulkanBackend, EnumerationReportsInsteadOfThrowing) {
    const auto devices = EnumerateVulkanDevices();
    ASSERT_TRUE(devices.has_value()) << devices.error().Description();
    if (devices->empty()) {
        GTEST_SKIP() << "no Vulkan device present";
    }
    for (const auto& info : *devices) {
        EXPECT_FALSE(info.name.empty());
    }
}

TEST(VulkanBackend, SlangKernelMatchesCpuReference) {
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

    const auto spirv = LoadSpirv(std::string(SIRIUS_KERNEL_DIR) + "/smoke.spv");
    ASSERT_FALSE(spirv.empty()) << "smoke.spv missing or empty";

    const auto kernel = (*device)->LoadKernel(spirv);
    ASSERT_TRUE(kernel.has_value()) << kernel.error().Description();

    constexpr std::uint32_t kCount = 4096;
    constexpr float kMass = 0.5f;
    std::vector<float> radii(kCount);
    for (std::uint32_t i = 0; i < kCount; ++i) {
        radii[i] = 1.0f + 0.01f * static_cast<float>(i);
    }
    const std::vector<float> params = {kMass, static_cast<float>(kCount)};

    const auto radii_buffer =
        (*device)->CreateBuffer(radii.size() * sizeof(float), BufferUsage::kStorage);
    const auto factors_buffer =
        (*device)->CreateBuffer(radii.size() * sizeof(float), BufferUsage::kStorage);
    const auto params_buffer =
        (*device)->CreateBuffer(params.size() * sizeof(float), BufferUsage::kStorage);
    ASSERT_TRUE(radii_buffer && factors_buffer && params_buffer);

    ASSERT_TRUE((*device)->WriteBuffer(
        *radii_buffer, std::as_bytes(std::span<const float>(radii))));
    ASSERT_TRUE((*device)->WriteBuffer(
        *params_buffer, std::as_bytes(std::span<const float>(params))));

    const sirius::backend::BufferHandle bindings[] = {*radii_buffer, *factors_buffer,
                                                      *params_buffer};
    const auto dispatched = (*device)->Dispatch(*kernel, bindings, (kCount + 63) / 64, 1, 1);
    ASSERT_TRUE(dispatched.has_value()) << dispatched.error().Description();

    std::vector<float> factors(kCount);
    ASSERT_TRUE(
        (*device)->ReadBuffer(*factors_buffer, std::as_writable_bytes(std::span(factors))));

    float max_difference = 0.0f;
    for (std::uint32_t i = 0; i < kCount; ++i) {
        const float reference = 1.0f - 2.0f * kMass / radii[i];
        max_difference = std::max(max_difference, std::abs(factors[i] - reference));
    }
    EXPECT_LE(max_difference, 1e-6f) << "kernel diverges from CPU reference";
#endif
}

}  // namespace
