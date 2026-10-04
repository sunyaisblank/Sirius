// CPU controls cover portability negotiation and precision refusal. Device
// checks cover enumeration, shader parity and worker-thread teardown. Missing
// devices skip the device checks; strict qualification rejects those skips.

#include "sirius/backend/device.h"
#include "sirius/backend/vulkan/vulkan_device.h"
#include "sirius/backend/vulkan/vulkan_portability.h"

#include <gtest/gtest.h>

#include "support/scoped_environment.h"

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
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

// CPU-only loader boundary. The production creation functions call these
// entry points; extension names are copied while their create-info is alive.
struct PortabilityProbe {
    inline static thread_local PortabilityProbe* current = nullptr;
    std::vector<const char*> advertised;
    std::vector<std::string> enabled;
    VkResult count_result = VK_SUCCESS;
    VkResult list_result = VK_SUCCESS;
    VkResult create_result = VK_SUCCESS;
    int creates = 0;
    VkInstanceCreateFlags flags = 0;
    const VkApplicationInfo* application = nullptr;
    const VkPhysicalDeviceFeatures* features = nullptr;
    const VkDeviceQueueCreateInfo* queues = nullptr;
    std::uint32_t queue_count = 0;
    VkPhysicalDevice physical = VK_NULL_HANDLE;

    PortabilityProbe() { current = this; }
    ~PortabilityProbe() { current = nullptr; }

    static VkResult Enumerate(std::uint32_t* count, VkExtensionProperties* properties) {
        if (properties == nullptr) {
            *count = static_cast<std::uint32_t>(current->advertised.size());
            return current->count_result;
        }
        const auto written =
            std::min(*count, static_cast<std::uint32_t>(current->advertised.size()));
        for (std::uint32_t i = 0; i < written; ++i) {
            properties[i] = {};
            std::strncpy(properties[i].extensionName, current->advertised[i],
                         VK_MAX_EXTENSION_NAME_SIZE - 1);
        }
        *count = written;
        return current->list_result;
    }
    static VKAPI_ATTR VkResult VKAPI_CALL InstanceExtensions(const char* layer,
                                                             std::uint32_t* count,
                                                             VkExtensionProperties* properties) {
        EXPECT_EQ(layer, nullptr);
        return Enumerate(count, properties);
    }
    static VKAPI_ATTR VkResult VKAPI_CALL DeviceExtensions(VkPhysicalDevice physical,
                                                           const char* layer, std::uint32_t* count,
                                                           VkExtensionProperties* properties) {
        EXPECT_EQ(layer, nullptr);
        current->physical = physical;
        return Enumerate(count, properties);
    }
    static VKAPI_ATTR VkResult VKAPI_CALL Instance(const VkInstanceCreateInfo* info,
                                                   const VkAllocationCallbacks* allocator,
                                                   VkInstance* result) {
        EXPECT_EQ(allocator, nullptr);
        ++current->creates;
        current->flags = info->flags;
        current->application = info->pApplicationInfo;
        for (std::uint32_t i = 0; i < info->enabledExtensionCount; ++i) {
            current->enabled.emplace_back(info->ppEnabledExtensionNames[i]);
        }
        *result = VK_NULL_HANDLE;  // No actual loader/device is touched by these tests.
        return current->create_result;
    }
    static VKAPI_ATTR VkResult VKAPI_CALL Device(VkPhysicalDevice physical,
                                                 const VkDeviceCreateInfo* info,
                                                 const VkAllocationCallbacks* allocator,
                                                 VkDevice* result) {
        EXPECT_EQ(allocator, nullptr);
        EXPECT_EQ(physical, current->physical);
        ++current->creates;
        current->features = info->pEnabledFeatures;
        current->queues = info->pQueueCreateInfos;
        current->queue_count = info->queueCreateInfoCount;
        for (std::uint32_t i = 0; i < info->enabledExtensionCount; ++i) {
            current->enabled.emplace_back(info->ppEnabledExtensionNames[i]);
        }
        *result = VK_NULL_HANDLE;
        return current->create_result;
    }
};

TEST(VulkanBackend, PortabilityInstanceOptInRequiresAdvertisedExtension) {
    for (const bool advertised : {false, true}) {
        PortabilityProbe probe;
        probe.advertised = {"VK_EXT_debug_utils"};
        if (advertised) probe.advertised.push_back("VK_KHR_portability_enumeration");
        const VkApplicationInfo application{.sType = VK_STRUCTURE_TYPE_APPLICATION_INFO};
        const char* existing[] = {"VK_EXT_debug_utils"};
        const VkInstanceCreateInfo info{
            .sType = VK_STRUCTURE_TYPE_INSTANCE_CREATE_INFO,
            .pApplicationInfo = &application,
            .enabledExtensionCount = 1,
            .ppEnabledExtensionNames = existing,
        };
        const auto result = sirius::backend::detail::CreateInstanceWithPortability(
            info, PortabilityProbe::InstanceExtensions, PortabilityProbe::Instance);
        ASSERT_TRUE(result.has_value()) << result.error().Description();
        EXPECT_EQ(probe.creates, 1);
        EXPECT_EQ(probe.application, &application);
        EXPECT_EQ(probe.flags,
                  advertised
                      ? VkInstanceCreateFlags{VK_INSTANCE_CREATE_ENUMERATE_PORTABILITY_BIT_KHR}
                      : 0u);
        const std::vector<std::string> expected =
            advertised
                ? std::vector<std::string>{"VK_EXT_debug_utils", "VK_KHR_portability_enumeration"}
                : std::vector<std::string>{"VK_EXT_debug_utils"};
        EXPECT_EQ(probe.enabled, expected);
    }
}

TEST(VulkanBackend, PortabilityDeviceEnablesSubsetWithoutChangingPrecisionOrQueues) {
    for (const bool advertised : {false, true}) {
        for (const VkBool32 fp64 : {VK_FALSE, VK_TRUE}) {
            PortabilityProbe probe;
            if (advertised) probe.advertised = {"VK_KHR_portability_subset"};
            const VkPhysicalDeviceFeatures features{.shaderFloat64 = fp64};
            const VkDeviceQueueCreateInfo queue{
                .sType = VK_STRUCTURE_TYPE_DEVICE_QUEUE_CREATE_INFO,
                .queueFamilyIndex = 3,
                .queueCount = 1,
            };
            const VkDeviceCreateInfo info{
                .sType = VK_STRUCTURE_TYPE_DEVICE_CREATE_INFO,
                .queueCreateInfoCount = 1,
                .pQueueCreateInfos = &queue,
                .pEnabledFeatures = &features,
            };
            const auto result = sirius::backend::detail::CreateDeviceWithPortability(
                VK_NULL_HANDLE, info, PortabilityProbe::DeviceExtensions, PortabilityProbe::Device);
            ASSERT_TRUE(result.has_value()) << result.error().Description();
            EXPECT_EQ(probe.creates, 1);
            EXPECT_EQ(probe.features, &features);
            EXPECT_EQ(probe.features->shaderFloat64, fp64);
            EXPECT_EQ(probe.queues, &queue);
            EXPECT_EQ(probe.queue_count, 1u);
            const std::vector<std::string> expected =
                advertised ? std::vector<std::string>{"VK_KHR_portability_subset"}
                           : std::vector<std::string>{};
            EXPECT_EQ(probe.enabled, expected);
        }
    }
}

TEST(VulkanBackend, PortabilityEnumerationErrorsDeclineBeforeCreation) {
    for (const bool count_failure : {false, true}) {
        for (const VkResult failure : {VK_ERROR_INITIALIZATION_FAILED, VK_INCOMPLETE}) {
            PortabilityProbe probe;
            probe.advertised = {"VK_KHR_portability_enumeration", "VK_KHR_portability_subset"};
            (count_failure ? probe.count_result : probe.list_result) = failure;
            const auto instance = sirius::backend::detail::CreateInstanceWithPortability(
                {.sType = VK_STRUCTURE_TYPE_INSTANCE_CREATE_INFO},
                PortabilityProbe::InstanceExtensions, PortabilityProbe::Instance);
            ASSERT_FALSE(instance.has_value());
            EXPECT_EQ(instance.error().operation(), "enumerate Vulkan instance extensions");
            const auto device = sirius::backend::detail::CreateDeviceWithPortability(
                VK_NULL_HANDLE, {.sType = VK_STRUCTURE_TYPE_DEVICE_CREATE_INFO},
                PortabilityProbe::DeviceExtensions, PortabilityProbe::Device);
            ASSERT_FALSE(device.has_value());
            EXPECT_EQ(device.error().operation(), "enumerate Vulkan device extensions");
            EXPECT_EQ(probe.creates, 0);
        }
    }
}

TEST(VulkanBackend, PortabilityCreationFailuresRemainExplicit) {
    PortabilityProbe probe;
    probe.create_result = VK_ERROR_EXTENSION_NOT_PRESENT;
    const auto instance = sirius::backend::detail::CreateInstanceWithPortability(
        {.sType = VK_STRUCTURE_TYPE_INSTANCE_CREATE_INFO}, PortabilityProbe::InstanceExtensions,
        PortabilityProbe::Instance);
    ASSERT_FALSE(instance.has_value());
    EXPECT_EQ(instance.error().operation(), "create Vulkan instance");
    const auto device = sirius::backend::detail::CreateDeviceWithPortability(
        VK_NULL_HANDLE, {.sType = VK_STRUCTURE_TYPE_DEVICE_CREATE_INFO},
        PortabilityProbe::DeviceExtensions, PortabilityProbe::Device);
    ASSERT_FALSE(device.has_value());
    EXPECT_EQ(device.error().operation(), "create Vulkan logical device");
    EXPECT_EQ(probe.creates, 2);
}

TEST(VulkanBackend, PipelineCacheImportRequiresExactDeviceIdentityAndHeader) {
    using sirius::backend::detail::VulkanPipelineCacheDataCompatible;
    // Fixed version-one wire header, not a host-structure reinterpretation.
    const std::array<std::byte, 32> header{
        std::byte{32}, std::byte{0},  std::byte{0},    std::byte{0},    std::byte{1},
        std::byte{0},  std::byte{0},  std::byte{0},    std::byte{0x02}, std::byte{0x10},
        std::byte{0},  std::byte{0},  std::byte{0xbf}, std::byte{0x15}, std::byte{0},
        std::byte{0},  std::byte{0},  std::byte{1},    std::byte{2},    std::byte{3},
        std::byte{4},  std::byte{5},  std::byte{6},    std::byte{7},    std::byte{8},
        std::byte{9},  std::byte{10}, std::byte{11},   std::byte{12},   std::byte{13},
        std::byte{14}, std::byte{15},
    };
    const VkPhysicalDeviceProperties properties{
        .vendorID = 0x1002,
        .deviceID = 0x15bf,
        .pipelineCacheUUID = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15},
    };
    ASSERT_TRUE(VulkanPipelineCacheDataCompatible(header, properties));
    for (std::size_t length = 0; length < header.size(); ++length)
        EXPECT_FALSE(
            VulkanPipelineCacheDataCompatible(std::span(header).first(length), properties));
    for (const std::size_t changed : {0, 4, 8, 12, 16, 31}) {
        auto incompatible = header;
        incompatible[changed] ^= std::byte{1};
        EXPECT_FALSE(VulkanPipelineCacheDataCompatible(incompatible, properties)) << changed;
    }
    auto big_endian = header;
    big_endian[0] = std::byte{0};
    big_endian[3] = std::byte{32};
    EXPECT_FALSE(VulkanPipelineCacheDataCompatible(big_endian, properties));
    auto extended_header = header;
    extended_header[0] = std::byte{33};
    EXPECT_FALSE(VulkanPipelineCacheDataCompatible(extended_header, properties));
}

TEST(VulkanBackend, PipelineCacheImportRejectsOversizeWithoutTruncation) {
    using sirius::backend::detail::kVulkanPipelineCacheBlobLimit;
    using sirius::backend::detail::VulkanPipelineCacheDataCompatible;
    const VkPhysicalDeviceProperties properties{};
    std::vector<std::byte> data(kVulkanPipelineCacheBlobLimit + 1);
    data[0] = std::byte{32};
    data[4] = std::byte{1};
    EXPECT_TRUE(VulkanPipelineCacheDataCompatible(
        std::span(data).first(kVulkanPipelineCacheBlobLimit), properties));
    EXPECT_FALSE(VulkanPipelineCacheDataCompatible(data, properties));
    // Compatibility controls touch only the guard, never submit fabricated
    // opaque driver data. Production imports exclusively Vulkan-exported bytes.
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
    // Both dispatches execute the same idempotent shader. A second call must
    // reuse this device's pipeline while still measuring actual queue work.
    sirius::backend::DispatchTiming timing;
    float max_difference = 0.0f;
    for (int invocation = 0; invocation < 2; ++invocation) {
        // A cached call that omits execution must not inherit the first result.
        const std::vector<float> unwritten(kCount, -42.0f);
        if (auto written = (*device)->WriteBuffer(*factors_buffer,
                                                  std::as_bytes(std::span<const float>(unwritten)));
            !written) {
            return std::unexpected(written.error());
        }
        if (auto dispatched =
                (*device)->Dispatch(*kernel, bindings, (kCount + 63) / 64, 1, 1, &timing);
            !dispatched) {
            return std::unexpected(dispatched.error());
        }
        if (!std::isfinite(timing.submit_wait_ms) || timing.submit_wait_ms <= 0.0 ||
            timing.submit_wait_ms > 1000.0) {
            return sirius::base::Fail(sirius::base::ErrorDomain::kDevice, "measure smoke dispatch",
                                      "actual submission timing is invalid or exceeds the stop");
        }
        const std::string prefix = invocation == 0 ? "dispatch_created_" : "dispatch_cached_";
        ::testing::Test::RecordProperty(prefix + "pipeline_setup_ms",
                                        std::to_string(timing.pipeline_setup_ms));
        ::testing::Test::RecordProperty(prefix + "command_setup_ms",
                                        std::to_string(timing.command_setup_ms));
        ::testing::Test::RecordProperty(prefix + "submit_wait_ms",
                                        std::to_string(timing.submit_wait_ms));
        ::testing::Test::RecordProperty(prefix + "cleanup_ms", std::to_string(timing.cleanup_ms));
        ::testing::Test::RecordProperty(prefix + "total_ms", std::to_string(timing.total_ms));
        EXPECT_EQ(timing.pipeline_created, invocation == 0);
        EXPECT_TRUE(std::isfinite(timing.pipeline_setup_ms));
        EXPECT_TRUE(std::isfinite(timing.command_setup_ms));
        EXPECT_TRUE(std::isfinite(timing.submit_wait_ms));
        EXPECT_TRUE(std::isfinite(timing.cleanup_ms));
        EXPECT_TRUE(std::isfinite(timing.total_ms));
        EXPECT_GE(timing.pipeline_setup_ms, 0.0);
        EXPECT_GE(timing.command_setup_ms, 0.0);
        EXPECT_GT(timing.submit_wait_ms, 0.0);
        EXPECT_GE(timing.cleanup_ms, 0.0);
        EXPECT_NEAR(timing.total_ms,
                    timing.pipeline_setup_ms + timing.command_setup_ms + timing.submit_wait_ms +
                        timing.cleanup_ms,
                    1e-9 * std::max(1.0, timing.total_ms));
        std::vector<float> factors(kCount);
        if (auto read =
                (*device)->ReadBuffer(*factors_buffer, std::as_writable_bytes(std::span(factors)));
            !read) {
            return std::unexpected(read.error());
        }
        for (std::uint32_t i = 0; i < kCount; ++i) {
            if (!std::isfinite(factors[i])) {
                return sirius::base::Fail(sirius::base::ErrorDomain::kDevice, "read smoke dispatch",
                                          "nonfinite shader output");
            }
            const float reference = 1.0f - 2.0f * kMass / radii[i];
            max_difference = std::max(max_difference, std::abs(factors[i] - reference));
        }
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

TEST(VulkanBackend, HostMemoryPreferencePreservesCompatibleHeapAndCoherence) {
    using sirius::backend::detail::VulkanHostMemoryType;
    constexpr auto coherent =
        VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT;
    constexpr auto cached = coherent | VK_MEMORY_PROPERTY_HOST_CACHED_BIT;
    VkPhysicalDeviceMemoryProperties properties{};
    properties.memoryHeapCount = 5;
    properties.memoryHeaps[0].size = 4096;
    properties.memoryHeaps[1].size = 4096;
    properties.memoryHeaps[2].size = 2048;
    properties.memoryHeaps[3].size = 16384;
    properties.memoryHeaps[4].size = 8192;
    properties.memoryTypeCount = 8;
    properties.memoryTypes[0] = {coherent, 0};
    properties.memoryTypes[1] = {cached, 1};
    properties.memoryTypes[2] = {cached, 0};
    properties.memoryTypes[3] = {cached, 2};
    properties.memoryTypes[4] = {
        VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_CACHED_BIT, 3};
    properties.memoryTypes[5] = {coherent, 4};
    properties.memoryTypes[6] = {cached, 4};
    properties.memoryTypes[7] = {cached | VK_MEMORY_PROPERTY_DEVICE_COHERENT_BIT_AMD, 4};

    EXPECT_EQ(VulkanHostMemoryType(properties, 0xffu), 6u);  // Largest coherent heap.
    EXPECT_EQ(VulkanHostMemoryType(properties, 0x0fu), 2u);  // Cache on the original heap.
    EXPECT_EQ(VulkanHostMemoryType(properties, 0x03u), 0u);  // Equal-size different heap.
    EXPECT_EQ(VulkanHostMemoryType(properties, 0x09u), 0u);  // Smaller cached heap.
    EXPECT_EQ(VulkanHostMemoryType(properties, 0x11u), 0u);  // Larger noncoherent heap.
    EXPECT_EQ(VulkanHostMemoryType(properties, 0xa0u), 5u);  // Do not add device-coherent flags.
    EXPECT_EQ(VulkanHostMemoryType(properties, 0x01u), 0u);  // Cached type incompatible.
    EXPECT_EQ(VulkanHostMemoryType(properties, 0x04u), 2u);  // Already cached.
    EXPECT_FALSE(VulkanHostMemoryType(properties, 0x10u));   // No coherent type.
    EXPECT_FALSE(VulkanHostMemoryType(properties, 0u));      // No compatible type.
}

TEST(VulkanBackend, BufferAllocationLimitCountsActualResidencyAndPreservesExistingBuffers) {
    const auto devices = EnumerateVulkanDevices();
    ASSERT_TRUE(devices.has_value()) << devices.error().Description();
    if (devices->empty()) GTEST_SKIP() << "no Vulkan device present";
    const auto selected = ResolveVulkanDeviceIndex(*devices);
    ASSERT_TRUE(selected.has_value()) << selected.error().Description();
    auto opened = CreateVulkanDevice(*selected);
    ASSERT_TRUE(opened.has_value()) << opened.error().Description();
    auto& device = **opened;
    EXPECT_EQ(device.BufferAllocationBytes(), 0u);
    ASSERT_TRUE(device.SetBufferAllocationLimit(0).has_value());
    EXPECT_FALSE(device.RequiredBufferAllocationBytes(0, BufferUsage::kStorage));
    const auto required = device.RequiredBufferAllocationBytes(1024, BufferUsage::kStorage);
    ASSERT_TRUE(required) << required.error().Description();
    EXPECT_GE(*required, 1024U);
    EXPECT_EQ(device.BufferAllocationBytes(), 0U);
    EXPECT_FALSE(device.CreateBuffer(1, BufferUsage::kStorage).has_value());
    EXPECT_EQ(device.BufferAllocationBytes(), 0u);
    ASSERT_TRUE(device.SetBufferAllocationLimit(1024 * 1024).has_value());
    const auto buffer = device.CreateBuffer(1024, BufferUsage::kStorage);
    ASSERT_TRUE(buffer.has_value()) << buffer.error().Description();
    auto resident = device.BufferAllocationBytes();
    EXPECT_GE(resident, 1024u);
    EXPECT_EQ(resident, *required);
    const std::array<std::uint32_t, 4> sent{0x12345678u, 0xffffffffu, 0u, 0xabcdef01u};
    ASSERT_TRUE(device.WriteBuffer(*buffer, std::as_bytes(std::span(sent))).has_value());
    ASSERT_TRUE(device.SetBufferAllocationLimit(resident + 1).has_value());
    EXPECT_FALSE(device.CreateBuffer(2, BufferUsage::kStorage).has_value());
    EXPECT_EQ(device.BufferAllocationBytes(), resident);
    EXPECT_FALSE(device.SetBufferAllocationLimit(resident - 1).has_value());
    EXPECT_FALSE(device.CreateBuffer(2, BufferUsage::kStorage).has_value());
    EXPECT_EQ(device.BufferAllocationBytes(), resident);
    std::array<std::uint32_t, 4> received{};
    ASSERT_TRUE(
        device.ReadBuffer(*buffer, std::as_writable_bytes(std::span(received))).has_value());
    EXPECT_EQ(received, sent);
    // Observe this driver's allocation requirement for a one-byte buffer, then
    // make the same request with one fewer byte of actual residency available.
    ASSERT_TRUE(device.SetBufferAllocationLimit(1024 * 1024).has_value());
    const auto tiny_required = device.RequiredBufferAllocationBytes(1, BufferUsage::kStorage);
    ASSERT_TRUE(tiny_required) << tiny_required.error().Description();
    EXPECT_EQ(device.BufferAllocationBytes(), resident);
    const auto tiny = device.CreateBuffer(1, BufferUsage::kStorage);
    ASSERT_TRUE(tiny.has_value());
    const auto tiny_allocation = device.BufferAllocationBytes() - resident;
    EXPECT_GE(tiny_allocation, 1u);
    EXPECT_EQ(tiny_allocation, *tiny_required);
    RecordProperty("one_byte_buffer_actual_allocation", std::to_string(tiny_allocation));
    resident = device.BufferAllocationBytes();
    ASSERT_TRUE(device.SetBufferAllocationLimit(resident + tiny_allocation - 1).has_value());
    const auto padded = device.CreateBuffer(1, BufferUsage::kStorage);
    EXPECT_FALSE(padded.has_value());
    if (!padded && tiny_allocation > 1) {
        EXPECT_NE(padded.error().detail().find("actual Vulkan allocation"), std::string::npos);
    }
    EXPECT_EQ(device.BufferAllocationBytes(), resident);
    ASSERT_TRUE(
        device.ReadBuffer(*buffer, std::as_writable_bytes(std::span(received))).has_value());
    EXPECT_EQ(received, sent);
    // The exact full residency cap rejects another allocation without damage.
    ASSERT_TRUE(device.SetBufferAllocationLimit(resident).has_value());
    EXPECT_FALSE(device.CreateBuffer(1, BufferUsage::kStorage).has_value());
    EXPECT_EQ(device.BufferAllocationBytes(), resident);
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
        const auto absent = ResolveVulkanDeviceIndex({});
        ASSERT_FALSE(absent);
        EXPECT_EQ(absent.error().domain(), sirius::base::ErrorDomain::kDevice);
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
        const auto absent = ResolveVulkanDeviceIndex({});
        ASSERT_FALSE(absent);
        EXPECT_EQ(absent.error().domain(), sirius::base::ErrorDomain::kConfiguration);
    }
    ScopedEnvironmentVariable selector("SIRIUS_VULKAN_DEVICE", "0");
    const auto absent = ResolveVulkanDeviceIndex({});
    ASSERT_FALSE(absent);
    EXPECT_EQ(absent.error().domain(), sirius::base::ErrorDomain::kConfiguration);
}

}  // namespace
