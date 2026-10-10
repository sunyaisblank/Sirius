#include "support/test_resource.h"

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
#include <format>
#include <fstream>
#include <limits>
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
    bool supports_fma32 = false;
    unsigned feature_queries = 0;
    const void* create_chain = nullptr;
    std::optional<sirius::backend::detail::ShaderFmaFeatures> fma;

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
        current->create_chain = info->pNext;
        if (info->pNext && static_cast<const VkBaseInStructure*>(info->pNext)->sType ==
                               sirius::backend::detail::kShaderFmaFeaturesType) {
            current->fma =
                *static_cast<const sirius::backend::detail::ShaderFmaFeatures*>(info->pNext);
        }
        for (std::uint32_t i = 0; i < info->enabledExtensionCount; ++i) {
            current->enabled.emplace_back(info->ppEnabledExtensionNames[i]);
        }
        *result = VK_NULL_HANDLE;
        return current->create_result;
    }
    static VKAPI_ATTR void VKAPI_CALL Features(VkPhysicalDevice physical,
                                               VkPhysicalDeviceFeatures2* features) {
        EXPECT_EQ(physical, current->physical);
        EXPECT_EQ(features->sType, VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_FEATURES_2);
        ASSERT_NE(features->pNext, nullptr);
        auto* fma = static_cast<sirius::backend::detail::ShaderFmaFeatures*>(features->pNext);
        EXPECT_EQ(fma->sType, sirius::backend::detail::kShaderFmaFeaturesType);
        EXPECT_EQ(fma->pNext, nullptr);
        ++current->feature_queries;
        fma->shaderFmaFloat16 = VK_TRUE;
        fma->shaderFmaFloat32 = current->supports_fma32 ? VK_TRUE : VK_FALSE;
        fma->shaderFmaFloat64 = VK_TRUE;
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

TEST(VulkanBackend, Fma32DeviceAdmissionRequiresAdvertisedFeatureAndEveryFloatControl) {
    for (unsigned mask = 0; mask < 32; ++mask) {
        for (const VkBool32 fp64 : {VK_FALSE, VK_TRUE}) {
            SCOPED_TRACE(mask);
            PortabilityProbe probe;
            probe.advertised = {"VK_KHR_portability_subset"};
            if (mask & 1u) probe.advertised.push_back("VK_KHR_shader_fma");
            probe.supports_fma32 = (mask & 2u) != 0;
            sirius::backend::DeviceInfo device{
                .preserves_fp32_denormals = (mask & 4u) != 0,
                .rounds_fp32_to_nearest = (mask & 8u) != 0,
                .preserves_fp32_signed_zero_inf_nan = (mask & 16u) != 0,
                .fma_fp32_enabled = true,
            };
            const VkPhysicalDeviceFeatures features{.shaderFloat64 = fp64};
            const VkDeviceQueueCreateInfo queue{.sType = VK_STRUCTURE_TYPE_DEVICE_QUEUE_CREATE_INFO,
                                                .queueCount = 1};
            const VkPhysicalDeviceTimelineSemaphoreFeatures caller_chain{
                .sType = VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_TIMELINE_SEMAPHORE_FEATURES};
            const VkDeviceCreateInfo info{.sType = VK_STRUCTURE_TYPE_DEVICE_CREATE_INFO,
                                          .pNext = &caller_chain,
                                          .queueCreateInfoCount = 1,
                                          .pQueueCreateInfos = &queue,
                                          .pEnabledFeatures = &features};
            const auto result = sirius::backend::detail::CreateDeviceWithPortability(
                VK_NULL_HANDLE, info, PortabilityProbe::DeviceExtensions, PortabilityProbe::Device,
                &device, PortabilityProbe::Features);
            ASSERT_TRUE(result) << result.error().Description();
            const bool admitted = mask == 31u;
            EXPECT_EQ(device.fma_fp32_enabled, admitted);
            EXPECT_EQ(probe.feature_queries, (mask & 29u) == 29u ? 1u : 0u);
            EXPECT_EQ(probe.fma.has_value(), admitted);
            EXPECT_EQ(probe.features, &features);
            EXPECT_EQ(probe.features->shaderFloat64, fp64);
            EXPECT_EQ(probe.queues, &queue);
            if (admitted) {
                EXPECT_EQ(probe.fma->pNext, &caller_chain);
                EXPECT_EQ(probe.fma->shaderFmaFloat16, VK_FALSE);
                EXPECT_EQ(probe.fma->shaderFmaFloat32, VK_TRUE);
                EXPECT_EQ(probe.fma->shaderFmaFloat64, VK_FALSE);
                EXPECT_EQ(probe.enabled, (std::vector<std::string>{"VK_KHR_portability_subset",
                                                                   "VK_KHR_shader_fma"}));
            } else {
                EXPECT_EQ(probe.create_chain, &caller_chain);
                EXPECT_EQ(probe.enabled, (std::vector<std::string>{"VK_KHR_portability_subset"}));
            }
        }
    }
    for (const bool enumeration_failure : {false, true}) {
        PortabilityProbe probe;
        probe.advertised = {"VK_KHR_shader_fma"};
        probe.supports_fma32 = true;
        if (enumeration_failure)
            probe.list_result = VK_INCOMPLETE;
        else
            probe.create_result = VK_ERROR_INITIALIZATION_FAILED;
        sirius::backend::DeviceInfo device{.preserves_fp32_denormals = true,
                                           .rounds_fp32_to_nearest = true,
                                           .preserves_fp32_signed_zero_inf_nan = true,
                                           .fma_fp32_enabled = true};
        const auto result = sirius::backend::detail::CreateDeviceWithPortability(
            VK_NULL_HANDLE, {.sType = VK_STRUCTURE_TYPE_DEVICE_CREATE_INFO},
            PortabilityProbe::DeviceExtensions, PortabilityProbe::Device, &device,
            PortabilityProbe::Features);
        EXPECT_FALSE(result);
        EXPECT_FALSE(device.fma_fp32_enabled);
        EXPECT_EQ(probe.creates, enumeration_failure ? 0 : 1);
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

TEST(VulkanBackend, KernelPrecisionDeclinesFmaBeforeDriverWorkUnlessBinary32WasEnabled) {
    using sirius::backend::ValidateVulkanKernelPrecision;
    const std::vector<std::uint32_t> fma = {
        0x07230203u, 0x00010500u, 0u,  8u,          0u, 0x00020011u, 1u, 0x00020011u, 6030u,
        0x00030016u, 1u,          32u, 0x0006114bu, 1u, 2u,          3u, 4u,          5u};
    EXPECT_TRUE(ValidateVulkanKernelPrecision(fma, false, true));
    EXPECT_TRUE(ValidateVulkanKernelPrecision(fma, true, true));
    const auto refused = ValidateVulkanKernelPrecision(fma, true);
    ASSERT_FALSE(refused);
    EXPECT_NE(refused.error().detail().find("shaderFmaFloat32"), std::string::npos);
    sirius::backend::VulkanDevice unopened;
    const auto unbound = unopened.LoadKernel(fma);
    ASSERT_FALSE(unbound);
    EXPECT_EQ(unbound.error().detail(), refused.error().detail());
    auto undeclared = fma;
    undeclared.erase(undeclared.begin() + 7, undeclared.begin() + 9);
    EXPECT_FALSE(ValidateVulkanKernelPrecision(undeclared, true));
    auto bad_extent = fma;
    bad_extent[12] = 0x0005114bu;
    EXPECT_FALSE(ValidateVulkanKernelPrecision(bad_extent, true, true));
    for (const auto width : {16u, 64u}) {
        auto unsupported = fma;
        unsupported[11] = width;
        EXPECT_FALSE(ValidateVulkanKernelPrecision(unsupported, true, true));
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
#ifndef SIRIUS_TEST_HAS_KERNELS
    GTEST_SKIP() << "kernels not compiled (slangc absent at configure time)";
#else
    const auto devices = EnumerateVulkanDevices();
    ASSERT_TRUE(devices.has_value()) << devices.error().Description();
    if (devices->empty()) {
        GTEST_SKIP() << "no Vulkan device present";
    }

    const auto selected = ResolveVulkanDeviceIndex(*devices);
    ASSERT_TRUE(selected.has_value()) << selected.error().Description();
    const auto spirv = LoadSpirv(sirius::test::ResourcePath("kernels/smoke.spv"));
    ASSERT_FALSE(spirv.empty()) << "smoke.spv missing or empty";

    const auto max_difference = DispatchSmokeKernel(*selected, spirv);
    ASSERT_TRUE(max_difference.has_value()) << max_difference.error().Description();
    EXPECT_LE(*max_difference, 1e-6f) << "kernel diverges from CPU reference";
#endif
}

TEST(VulkanBackend, WorkerThreadDispatchTearsDownSafely) {
#ifndef SIRIUS_TEST_HAS_KERNELS
    GTEST_SKIP() << "kernels not compiled (slangc absent at configure time)";
#else
    const auto devices = EnumerateVulkanDevices();
    ASSERT_TRUE(devices.has_value()) << devices.error().Description();
    if (devices->empty()) {
        GTEST_SKIP() << "no Vulkan device present";
    }

    const auto selected = ResolveVulkanDeviceIndex(*devices);
    ASSERT_TRUE(selected.has_value()) << selected.error().Description();
    const auto spirv = LoadSpirv(sirius::test::ResourcePath("kernels/smoke.spv"));
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

TEST(VulkanBackend, IdenticalKernelWordsReusePipelineAcrossBuffers) {
#ifndef SIRIUS_TEST_HAS_KERNELS
    GTEST_SKIP() << "kernels not compiled (slangc absent at configure time)";
#else
    const auto inventory = EnumerateVulkanDevices();
    ASSERT_TRUE(inventory) << inventory.error().Description();
    if (inventory->empty()) {
        GTEST_SKIP() << "no Vulkan device present";
    }
    const auto selected = ResolveVulkanDeviceIndex(*inventory);
    ASSERT_TRUE(selected) << selected.error().Description();
    auto opened = CreateVulkanDevice(*selected);
    ASSERT_TRUE(opened) << opened.error().Description();
    auto& device = **opened;
    auto words = LoadSpirv(sirius::test::ResourcePath("kernels/smoke.spv"));
    ASSERT_GE(words.size(), 5U);
    const auto unchanged = words;
    const auto first = device.LoadKernel(words);
    ASSERT_TRUE(first) << first.error().Description();
    const auto again = device.LoadKernel(unchanged);
    ASSERT_TRUE(again) << again.error().Description();
    EXPECT_EQ(again->value, first->value);

    // Change the smoke factor from 1-2M/r to 1-4M/r. A different complete
    // module must not inherit the original handle or executable pipeline.
    std::uint32_t float_type = 0;
    std::size_t changed_constants = 0;
    for (std::size_t offset = 5; offset < words.size();) {
        const auto extent = words[offset] >> 16;
        const auto opcode = words[offset] & 0xffffU;
        ASSERT_GT(extent, 0U);
        ASSERT_LE(extent, words.size() - offset);
        if (opcode == 22U && extent == 3U && words[offset + 2] == 32U)
            float_type = words[offset + 1];
        if (opcode == 43U && extent == 4U && words[offset + 1] == float_type &&
            words[offset + 3] == 0x40000000U) {
            words[offset + 3] = 0x40800000U;
            ++changed_constants;
        }
        offset += extent;
    }
    ASSERT_EQ(changed_constants, 1U);
    const auto variant = device.LoadKernel(words);
    ASSERT_TRUE(variant) << variant.error().Description();
    EXPECT_NE(variant->value, first->value);
    const auto after_mutation = device.LoadKernel(unchanged);
    ASSERT_TRUE(after_mutation) << after_mutation.error().Description();
    EXPECT_EQ(after_mutation->value, first->value);

    constexpr std::array<float, 4> radii{1, 2, 4, 8};
    constexpr std::array<float, 3> masses{.5f, .75f, .5f};
    constexpr std::array<std::array<float, 4>, 3> expected{
        {{0, .5f, .75f, .875f}, {-.5f, .25f, .625f, .8125f}, {-1, 0, .5f, .75f}}};
    std::array<sirius::backend::BufferHandle, 3> prior_outputs{};
    std::array<std::array<float, 4>, 3> saved_outputs{};
    for (std::size_t invocation = 0; invocation < expected.size(); ++invocation) {
        const std::array<float, 2> params{masses[invocation], 4};
        constexpr std::array<float, 4> unwritten{-42, -42, -42, -42};
        const auto in = device.CreateBuffer(sizeof(radii), BufferUsage::kStorage);
        const auto out = device.CreateBuffer(sizeof(unwritten), BufferUsage::kStorage);
        const auto parameters = device.CreateBuffer(sizeof(params), BufferUsage::kStorage);
        ASSERT_TRUE(in);
        ASSERT_TRUE(out);
        ASSERT_TRUE(parameters);
        ASSERT_TRUE(device.WriteBuffer(*in, std::as_bytes(std::span(radii))));
        ASSERT_TRUE(device.WriteBuffer(*out, std::as_bytes(std::span(unwritten))));
        ASSERT_TRUE(device.WriteBuffer(*parameters, std::as_bytes(std::span(params))));
        const std::array bindings{*in, *out, *parameters};
        sirius::backend::DispatchTiming timing;
        ASSERT_TRUE(device.Dispatch(invocation == 2   ? *variant
                                    : invocation == 0 ? *first
                                                      : *again,
                                    bindings, 1, 1, 1, &timing));
        EXPECT_EQ(timing.pipeline_created, invocation != 1);
        prior_outputs[invocation] = *out;
        // Rebinding a cached pipeline must write this output, leaving every
        // previous output unchanged. Keep the original smoke tolerance:
        // Vulkan division need not round each exactly representable quotient.
        for (std::size_t prior = 0; prior <= invocation; ++prior) {
            std::array<float, 4> actual{};
            ASSERT_TRUE(
                device.ReadBuffer(prior_outputs[prior], std::as_writable_bytes(std::span(actual))));
            if (prior == invocation) {
                for (std::size_t value = 0; value < actual.size(); ++value)
                    EXPECT_NEAR(actual[value], expected[prior][value], 1e-6f);
                saved_outputs[prior] = actual;
            } else {
                EXPECT_EQ(std::memcmp(actual.data(), saved_outputs[prior].data(), sizeof(actual)),
                          0);
            }
        }
    }
#endif
}

TEST(VulkanBackend, IndependentPairCompletesDistinctKernelsAndRejectsSharedBuffers) {
#ifndef SIRIUS_TEST_HAS_KERNELS
    GTEST_SKIP() << "kernels not compiled (slangc absent at configure time)";
#else
    const auto inventory = EnumerateVulkanDevices();
    ASSERT_TRUE(inventory) << inventory.error().Description();
    if (inventory->empty()) GTEST_SKIP() << "no Vulkan device present";
    const auto selected = ResolveVulkanDeviceIndex(*inventory);
    ASSERT_TRUE(selected) << selected.error().Description();
    auto opened = CreateVulkanDevice(*selected);
    ASSERT_TRUE(opened) << opened.error().Description();
    auto& device = **opened;
    ASSERT_TRUE(device.SupportsIndependentPair());
    auto* vulkan = dynamic_cast<sirius::backend::VulkanDevice*>(&device);
    ASSERT_NE(vulkan, nullptr);
    const auto timestamp_properties = vulkan->TimestampProperties();
    const bool timestamps_supported =
        timestamp_properties.valid_bits >= 36 && timestamp_properties.valid_bits <= 64 &&
        std::isfinite(timestamp_properties.period_ns) && timestamp_properties.period_ns > 0;
    EXPECT_FALSE(vulkan->LastDispatchTimestamp());  // Default dispatch has no queries.
    RecordProperty("timestamp_valid_bits", std::to_string(timestamp_properties.valid_bits));
    RecordProperty("timestamp_period_ns", std::format("{:.17g}", timestamp_properties.period_ns));
    auto words = LoadSpirv(sirius::test::ResourcePath("kernels/smoke.spv"));
    ASSERT_GE(words.size(), 5U);
    const auto original = device.LoadKernel(words);
    ASSERT_TRUE(original);
    std::uint32_t float_type = 0;
    std::size_t changed = 0;
    for (std::size_t offset = 5; offset < words.size();) {
        const auto extent = words[offset] >> 16, opcode = words[offset] & 0xffffU;
        ASSERT_GT(extent, 0U);
        ASSERT_LE(extent, words.size() - offset);
        if (opcode == 22U && extent == 3U && words[offset + 2] == 32U)
            float_type = words[offset + 1];
        if (opcode == 43U && extent == 4U && words[offset + 1] == float_type &&
            words[offset + 3] == 0x40000000U) {
            words[offset + 3] = 0x40800000U;
            ++changed;
        }
        offset += extent;
    }
    ASSERT_EQ(changed, 1U);
    const auto variant = device.LoadKernel(words);
    ASSERT_TRUE(variant);
    std::array<std::array<sirius::backend::BufferHandle, 3>, 2> bindings{};
    std::array<float, 194> radii{}, untouched{};
    for (std::size_t i = 0; i < radii.size(); ++i) radii[i] = float(4U << (i % 3));
    untouched.fill(-42);
    constexpr std::array<unsigned, 2> counts{67, 131};
    constexpr std::array<float, 2> masses{.5f, .75f};
    for (std::size_t i = 0; i < bindings.size(); ++i) {
        const std::array<float, 2> parameters{masses[i], float(counts[i])};
        const auto in = device.CreateBuffer(sizeof(radii), BufferUsage::kStorage);
        const auto out = device.CreateBuffer(sizeof(untouched), BufferUsage::kStorage);
        const auto params = device.CreateBuffer(sizeof(parameters), BufferUsage::kStorage);
        ASSERT_TRUE(in);
        ASSERT_TRUE(out);
        ASSERT_TRUE(params);
        bindings[i] = {*in, *out, *params};
        ASSERT_TRUE(device.WriteBuffer(*in, std::as_bytes(std::span(radii))));
        ASSERT_TRUE(device.WriteBuffer(*params, std::as_bytes(std::span(parameters))));
    }
    const auto resident = device.BufferAllocationBytes();
    std::array<sirius::backend::ComputeDispatch, 2> commands{
        sirius::backend::ComputeDispatch{*original, bindings[0], 2, 1, 1},
        sirius::backend::ComputeDispatch{*variant, bindings[1], 3, 1, 1}};
    sirius::backend::IndependentPairTiming timing;
    const auto check_timestamp = [&](const std::string& prefix,
                                     const sirius::backend::DispatchTiming& host, bool enabled) {
        const auto observation = vulkan->LastDispatchTimestamp();
        if (!enabled) {
            EXPECT_FALSE(observation);
            return;
        }
        ASSERT_TRUE(observation);
        EXPECT_EQ(observation->query_result, VK_SUCCESS);
        EXPECT_NE(observation->availability[0], 0U);
        EXPECT_NE(observation->availability[1], 0U);
        ASSERT_TRUE(observation->device_span_ms);
        EXPECT_TRUE(std::isfinite(*observation->device_span_ms));
        EXPECT_GE(*observation->device_span_ms, 0);
        EXPECT_NEAR(observation->host_submit_ms + observation->host_wait_ms, host.submit_wait_ms,
                    1e-9 * std::max(1., host.submit_wait_ms));
        RecordProperty(prefix + "_device_span_ms",
                       std::format("{:.17g}", *observation->device_span_ms));
        RecordProperty(prefix + "_host_submit_ms",
                       std::format("{:.17g}", observation->host_submit_ms));
        RecordProperty(prefix + "_host_wait_ms", std::format("{:.17g}", observation->host_wait_ms));
        RecordProperty(prefix + "_begin_tick", std::to_string(observation->ticks[0]));
        RecordProperty(prefix + "_end_tick", std::to_string(observation->ticks[1]));
    };
    for (unsigned invocation = 0; invocation < 4; ++invocation) {
        const bool requested = invocation == 1 || invocation == 3;
        const auto configured = vulkan->SetDispatchTimestampsEnabled(requested);
        if (!requested || timestamps_supported) {
            ASSERT_TRUE(configured) << configured.error().Description();
        } else {
            EXPECT_FALSE(configured);  // Unsupported diagnostics do not disable ordinary work.
        }
        for (const auto& bound : bindings)
            ASSERT_TRUE(device.WriteBuffer(bound[1], std::as_bytes(std::span(untouched))));
        if (invocation == 1) {
            sirius::backend::DispatchTiming single;
            const auto& first_command = commands[0];
            ASSERT_TRUE(device.Dispatch(first_command.kernel, first_command.buffers,
                                        first_command.groups_x, 1, 1, &single));
            check_timestamp("single", single, timestamps_supported);
        }
        ASSERT_TRUE(device.DispatchIndependentPair(commands, &timing));
        check_timestamp("pair_" + std::to_string(invocation), timing.combined,
                        requested && timestamps_supported);
        EXPECT_EQ(timing.pipeline_creations, invocation == 0 ? 2U : 0U);
        EXPECT_EQ(timing.combined.pipeline_created, invocation == 0);
        EXPECT_TRUE(std::isfinite(timing.combined.submit_wait_ms));
        EXPECT_GT(timing.combined.submit_wait_ms, 0);
        EXPECT_NEAR(timing.combined.total_ms,
                    timing.combined.pipeline_setup_ms + timing.combined.command_setup_ms +
                        timing.combined.submit_wait_ms + timing.combined.cleanup_ms,
                    1e-9 * std::max(1., timing.combined.total_ms));
        EXPECT_EQ(device.BufferAllocationBytes(), resident);
        for (std::size_t row = 0; row < bindings.size(); ++row) {
            std::array<float, 194> actual{};
            ASSERT_TRUE(
                device.ReadBuffer(bindings[row][1], std::as_writable_bytes(std::span(actual))));
            for (std::size_t i = 0; i < actual.size(); ++i) {
                if (i < counts[row])
                    EXPECT_NEAR(actual[i], 1 - (row == 0 ? 2 : 4) * masses[row] / radii[i], 1e-6f);
                else
                    EXPECT_EQ(actual[i], untouched[i]);
            }
        }
        std::swap(commands[0], commands[1]);
    }
    for (const auto& bound : bindings)
        ASSERT_TRUE(device.WriteBuffer(bound[1], std::as_bytes(std::span(untouched))));
    auto overlapping = bindings[1];
    overlapping[1] = bindings[0][1];
    commands[1].buffers = overlapping;
    const auto refused = device.DispatchIndependentPair(commands, &timing);
    ASSERT_FALSE(refused);
    EXPECT_FALSE(vulkan->LastDispatchTimestamp());  // Refusal cannot retain a prior sample.
    EXPECT_EQ(refused.error().operation(), "dispatch independent compute pair");
    EXPECT_EQ(timing.combined.total_ms, 0);
    EXPECT_EQ(timing.pipeline_creations, 0U);
    ASSERT_TRUE(vulkan->SetDispatchTimestampsEnabled(false));
    for (const auto& bound : bindings) {
        std::array<float, 194> actual{};
        ASSERT_TRUE(device.ReadBuffer(bound[1], std::as_writable_bytes(std::span(actual))));
        EXPECT_EQ(actual, untouched);
    }
    // A refused pair must leave the completed primary usable for a different
    // single stream and then a pair, without stale bindings or query commands.
    commands[1].buffers = bindings[1];
    sirius::backend::DispatchTiming recovered_single;
    ASSERT_TRUE(device.Dispatch(commands[0].kernel, commands[0].buffers, commands[0].groups_x, 1, 1,
                                &recovered_single));
    EXPECT_FALSE(vulkan->LastDispatchTimestamp());
    for (std::size_t row = 0; row < bindings.size(); ++row) {
        std::array<float, 194> actual{};
        ASSERT_TRUE(device.ReadBuffer(bindings[row][1], std::as_writable_bytes(std::span(actual))));
        for (std::size_t i = 0; i < actual.size(); ++i)
            if (row == 0 && i < counts[row])
                EXPECT_NEAR(actual[i], 1 - 2 * masses[row] / radii[i], 1e-6f);
            else
                EXPECT_EQ(actual[i], untouched[i]);
    }
    ASSERT_TRUE(device.DispatchIndependentPair(commands, &timing));
    EXPECT_FALSE(vulkan->LastDispatchTimestamp());
    EXPECT_EQ(timing.pipeline_creations, 0U);
    EXPECT_EQ(device.BufferAllocationBytes(), resident);
    for (std::size_t row = 0; row < bindings.size(); ++row) {
        std::array<float, 194> actual{};
        ASSERT_TRUE(device.ReadBuffer(bindings[row][1], std::as_writable_bytes(std::span(actual))));
        for (std::size_t i = 0; i < actual.size(); ++i)
            if (i < counts[row])
                EXPECT_NEAR(actual[i], 1 - (row == 0 ? 2 : 4) * masses[row] / radii[i], 1e-6f);
            else
                EXPECT_EQ(actual[i], untouched[i]);
    }
#endif
}

TEST(VulkanBackend, TimestampSpansBoundWrapAndInvalidCounters) {
    using sirius::backend::detail::VulkanTimestampSpanMs;
    const auto simple = VulkanTimestampSpanMs(100, 110, 36, 2, .1);
    ASSERT_TRUE(simple);
    EXPECT_NEAR(*simple, .000020, 1e-20);
    // Independent fixed counter witnesses: 30 ticks through a 36-bit wrap,
    // and four through a 64-bit wrap. Upper undefined bits are deliberately set.
    const auto wrapped =
        VulkanTimestampSpanMs(0xffffffffffffffecULL, 0xabcde0000000000aULL, 36, 2, .1);
    ASSERT_TRUE(wrapped);
    EXPECT_NEAR(*wrapped, .000060, 1e-20);
    const auto full = VulkanTimestampSpanMs(0xfffffffffffffffdULL, 1, 64, 2, .1);
    ASSERT_TRUE(full);
    EXPECT_NEAR(*full, .000008, 1e-20);
    EXPECT_FALSE(VulkanTimestampSpanMs(100, 110, 36, 2, 137438.953472));
    EXPECT_FALSE(VulkanTimestampSpanMs(100, 110, 36, 2, 200000));
    for (const auto bits : {0U, 35U, 65U})
        EXPECT_FALSE(VulkanTimestampSpanMs(100, 110, bits, 2, .1));
    for (const auto period : {0., -1., std::numeric_limits<double>::infinity(),
                              std::numeric_limits<double>::quiet_NaN()})
        EXPECT_FALSE(VulkanTimestampSpanMs(100, 110, 36, period, .1));
    for (const auto completion :
         {-1., std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()})
        EXPECT_FALSE(VulkanTimestampSpanMs(100, 110, 36, 2, completion));
}

TEST(VulkanBackend, SubmissionErrorsDoNotImplyIdleCompletion) {
    using sirius::backend::detail::VulkanSubmissionNeedsCompletion;
    // Fixed Vulkan lifetime decisions, not injected device or physics outcomes.
    for (const auto waited : {VK_SUCCESS, VK_ERROR_DEVICE_LOST, VK_ERROR_OUT_OF_HOST_MEMORY,
                              VK_ERROR_OUT_OF_DEVICE_MEMORY}) {
        SCOPED_TRACE(waited);
        EXPECT_TRUE(VulkanSubmissionNeedsCompletion(VK_ERROR_DEVICE_LOST, waited));
        EXPECT_FALSE(VulkanSubmissionNeedsCompletion(VK_ERROR_OUT_OF_HOST_MEMORY, waited));
        EXPECT_FALSE(VulkanSubmissionNeedsCompletion(VK_ERROR_OUT_OF_DEVICE_MEMORY, waited));
    }
    EXPECT_FALSE(VulkanSubmissionNeedsCompletion(VK_SUCCESS, VK_SUCCESS));
    EXPECT_FALSE(VulkanSubmissionNeedsCompletion(VK_SUCCESS, VK_ERROR_DEVICE_LOST));
    EXPECT_TRUE(VulkanSubmissionNeedsCompletion(VK_SUCCESS, VK_ERROR_OUT_OF_HOST_MEMORY));
    EXPECT_TRUE(VulkanSubmissionNeedsCompletion(VK_SUCCESS, VK_ERROR_OUT_OF_DEVICE_MEMORY));
    EXPECT_TRUE(VulkanSubmissionNeedsCompletion(VK_ERROR_UNKNOWN, VK_SUCCESS));
    EXPECT_TRUE(VulkanSubmissionNeedsCompletion(VK_SUCCESS, VK_ERROR_UNKNOWN));
    EXPECT_TRUE(VulkanSubmissionNeedsCompletion(VK_SUCCESS, VK_TIMEOUT));
    EXPECT_TRUE(VulkanSubmissionNeedsCompletion(VK_SUCCESS, VK_NOT_READY));
}

}  // namespace
