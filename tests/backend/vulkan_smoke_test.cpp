// Postconditions verified here: device enumeration reports rather than
// throws; a Slang-compiled kernel dispatched through the Vulkan adapter
// reproduces the CPU computation of f(r) = 1 - 2M/r within fp32 rounding.
// On machines with no Vulkan ICD both tests skip, cleanly and loudly.

#include "sirius/backend/device.h"
#include "sirius/backend/vulkan/vulkan_device.h"

#include <gtest/gtest.h>

#include "support/scoped_environment.h"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <optional>
#include <span>
#include <thread>
#include <vector>

namespace {

using sirius::backend::BufferUsage;
using sirius::backend::ComputeDevice;
using sirius::backend::CreateVulkanDevice;
using sirius::backend::EnumerateVulkanDevices;
using sirius::backend::ResolveVulkanDeviceIndex;
using sirius::test::ScopedEnvironmentVariable;

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

TEST(VulkanBackend, KernelPrecisionDeclinesUnsupportedFloat64AndMalformedInstructions) {
    using sirius::backend::ValidateVulkanKernelPrecision;
    // Independent SPIR-V wire fixtures: Shader capability followed by Float64.
    const std::vector<std::uint32_t> fp32 = {0x07230203u, 0x00010500u, 0u, 1u, 0u, 0x00020011u, 1u};
    auto fp64 = fp32;
    fp64.insert(fp64.end(), {0x00020011u, 10u});
    EXPECT_TRUE(ValidateVulkanKernelPrecision(fp32, false));
    EXPECT_TRUE(ValidateVulkanKernelPrecision(fp32, true));
    EXPECT_TRUE(ValidateVulkanKernelPrecision(fp64, true));
    const auto refused = ValidateVulkanKernelPrecision(fp64, false);
    ASSERT_FALSE(refused.has_value());
    EXPECT_EQ(refused.error().domain(), sirius::base::ErrorDomain::kKernel);
    EXPECT_NE(refused.error().detail().find("shaderFloat64"), std::string::npos);

    // The real adapter must decline before touching an uninitialised Vulkan
    // handle. This makes removal of the production call a failing control.
    sirius::backend::VulkanDevice unopened;
    const auto unbound = unopened.LoadKernel(fp64);
    ASSERT_FALSE(unbound.has_value());
    EXPECT_EQ(unbound.error().detail(), refused.error().detail());
    for (std::size_t length = 0; length <= 5; ++length) {
        EXPECT_FALSE(unopened.LoadKernel(std::span(fp32).first(length)).has_value());
    }
    auto bad_magic = fp32;
    bad_magic[0] = 0;
    auto bad_schema = fp32;
    bad_schema[4] = 1;
    auto zero_extent = fp32;
    zero_extent[5] = 17u;
    auto overflowing_extent = fp32;
    overflowing_extent[5] = 0xffff0011u;
    auto missing_capability = fp32;
    missing_capability.resize(6);
    missing_capability[5] = 0x00010011u;
    auto extra_capability_operand = fp64;
    extra_capability_operand[7] = 0x00030011u;
    extra_capability_operand.push_back(0u);
    for (const auto& malformed : {bad_magic, bad_schema, zero_extent, overflowing_extent,
                                  missing_capability, extra_capability_operand}) {
        const auto result = unopened.LoadKernel(malformed);
        ASSERT_FALSE(result.has_value());
        EXPECT_EQ(result.error().operation(), "validate shader module");
    }
}

sirius::base::Expected<float> DispatchSmokeKernel(std::size_t device_index,
                                                  std::span<const std::uint32_t> spirv) {
    auto device = CreateVulkanDevice(device_index);
    if (!device) return std::unexpected(device.error());

    auto kernel = (*device)->LoadKernel(spirv);
    if (!kernel) return std::unexpected(kernel.error());

    constexpr std::uint32_t kCount = 4096;
    constexpr float kMass = 0.5f;
    std::vector<float> radii(kCount);
    for (std::uint32_t i = 0; i < kCount; ++i) {
        radii[i] = 1.0f + 0.01f * static_cast<float>(i);
    }
    const std::vector<float> params = {kMass, static_cast<float>(kCount)};

    auto radii_buffer =
        (*device)->CreateBuffer(radii.size() * sizeof(float), BufferUsage::kStorage);
    if (!radii_buffer) return std::unexpected(radii_buffer.error());
    auto factors_buffer =
        (*device)->CreateBuffer(radii.size() * sizeof(float), BufferUsage::kStorage);
    if (!factors_buffer) return std::unexpected(factors_buffer.error());
    auto params_buffer =
        (*device)->CreateBuffer(params.size() * sizeof(float), BufferUsage::kStorage);
    if (!params_buffer) return std::unexpected(params_buffer.error());

    if (auto written =
            (*device)->WriteBuffer(*radii_buffer, std::as_bytes(std::span<const float>(radii)));
        !written) {
        return std::unexpected(written.error());
    }
    if (auto written =
            (*device)->WriteBuffer(*params_buffer, std::as_bytes(std::span<const float>(params)));
        !written) {
        return std::unexpected(written.error());
    }

    const sirius::backend::BufferHandle bindings[] = {*radii_buffer, *factors_buffer,
                                                      *params_buffer};
    if (auto dispatched = (*device)->Dispatch(*kernel, bindings, (kCount + 63) / 64, 1, 1);
        !dispatched) {
        return std::unexpected(dispatched.error());
    }

    std::vector<float> factors(kCount);
    if (auto read =
            (*device)->ReadBuffer(*factors_buffer, std::as_writable_bytes(std::span(factors)));
        !read) {
        return std::unexpected(read.error());
    }

    float max_difference = 0.0f;
    for (std::uint32_t i = 0; i < kCount; ++i) {
        const float reference = 1.0f - 2.0f * kMass / radii[i];
        max_difference = std::max(max_difference, std::abs(factors[i] - reference));
    }
    return max_difference;
}

TEST(VulkanBackend, EnumerationReportsInsteadOfThrowing) {
    const auto devices = EnumerateVulkanDevices();
    ASSERT_TRUE(devices.has_value()) << devices.error().Description();
    if (devices->empty()) {
        GTEST_SKIP() << "no Vulkan device present";
    }
    for (const auto& info : *devices) {
        EXPECT_FALSE(info.name.empty());
        EXPECT_FALSE(info.driver_name.empty())
            << "external attestation cannot distinguish Dozen/MoltenVK/native drivers";
        EXPECT_GT(info.api_version, 0u);
        EXPECT_GT(info.vendor_id, 0u);
        EXPECT_GT(info.render_memory_bytes, 0u)
            << "the adapter cannot govern buffers without an allocatable render heap";
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

    const auto selected = ResolveVulkanDeviceIndex(*devices);
    ASSERT_TRUE(selected.has_value()) << selected.error().Description();
    const auto spirv = LoadSpirv(std::string(SIRIUS_KERNEL_DIR) + "/smoke.spv");
    ASSERT_FALSE(spirv.empty()) << "smoke.spv missing or empty";

    const auto max_difference = DispatchSmokeKernel(*selected, spirv);
    ASSERT_TRUE(max_difference.has_value()) << max_difference.error().Description();
    EXPECT_LE(*max_difference, 1e-6f) << "kernel diverges from CPU reference";
#endif
}

TEST(VulkanBackend, WorkerThreadDispatchTearsDownSafely) {
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
    const auto spirv = LoadSpirv(std::string(SIRIUS_KERNEL_DIR) + "/smoke.spv");
    ASSERT_FALSE(spirv.empty()) << "smoke.spv missing or empty";

    std::optional<sirius::base::Expected<float> > result;
    std::thread worker([&] { result.emplace(DispatchSmokeKernel(*selected, spirv)); });
    worker.join();

    ASSERT_TRUE(result.has_value());
    ASSERT_TRUE(result->has_value()) << result->error().Description();
    EXPECT_LE(**result, 1e-6f) << "worker-thread kernel diverges from CPU reference";
#endif
}

TEST(VulkanBackend, DeviceSelectionIsStrictAndRangeChecked) {
    const std::vector<sirius::backend::DeviceInfo> devices = {
        {.name = "device zero"},
        {.name = "device one"},
    };

    {
        ScopedEnvironmentVariable selector("SIRIUS_VULKAN_DEVICE", nullptr);
        const auto selected = ResolveVulkanDeviceIndex(devices);
        ASSERT_TRUE(selected.has_value()) << selected.error().Description();
        EXPECT_EQ(*selected, 0u);
    }
    {
        ScopedEnvironmentVariable selector("SIRIUS_VULKAN_DEVICE", "1");
        const auto selected = ResolveVulkanDeviceIndex(devices);
        ASSERT_TRUE(selected.has_value()) << selected.error().Description();
        EXPECT_EQ(*selected, 1u);
    }
    for (const char* invalid : {"-1", "1tail", " 1", "2"}) {
        ScopedEnvironmentVariable selector("SIRIUS_VULKAN_DEVICE", invalid);
        const auto selected = ResolveVulkanDeviceIndex(devices);
        EXPECT_FALSE(selected.has_value()) << invalid;
    }
}

}  // namespace
