#include "sirius/backend/vulkan/vulkan_device.h"

#include "sirius/backend/vulkan/vulkan_portability.h"
#include "sirius/base/contracts.h"

#include <algorithm>
#include <charconv>
#include <chrono>
#include <cstdlib>
#include <cstring>
#include <format>
#include <mutex>
#include <new>
#include <string_view>
#include <utility>

#if defined(__linux__)
#include <dlfcn.h>
#endif

namespace sirius::backend {

namespace {

using base::ErrorDomain;
using base::Expected;
using base::Fail;

// Headers on the measured toolchain are Vulkan 1.3 while the loader is 1.4;
// requesting 1.3 is compatible with both (specification section 1.7 evidence).
constexpr std::uint32_t kApiVersion = VK_MAKE_API_VERSION(0, 1, 3, 0);

VkBufferCreateInfo BufferCreateInfo(std::uint64_t size_bytes, BufferUsage usage) {
    const VkBufferUsageFlags usage_flags = usage == BufferUsage::kStorage
                                               ? VK_BUFFER_USAGE_STORAGE_BUFFER_BIT
                                               : VK_BUFFER_USAGE_UNIFORM_BUFFER_BIT;
    return {
        .sType = VK_STRUCTURE_TYPE_BUFFER_CREATE_INFO,
        .size = size_bytes,
        .usage = usage_flags,
        .sharingMode = VK_SHARING_MODE_EXCLUSIVE,
    };
}

// One bounded serialized blob survives fresh render devices. Vulkan objects and
// buffers retain their existing per-device lifetime. Concurrent device exports
// replace this slot; they cannot accumulate entries for different adapters.
struct ProcessPipelineCache {
    std::mutex mutex;
    std::vector<std::byte> data;
};

ProcessPipelineCache& PipelineCacheStore() {
    static ProcessPipelineCache cache;
    return cache;
}

[[nodiscard]] std::string VkResultText(VkResult result) {
    return std::format("VkResult {}", static_cast<int>(result));
}

[[nodiscard]] DeviceKind ClassifyDevice(const VkPhysicalDeviceProperties& properties) {
    switch (properties.deviceType) {
        case VK_PHYSICAL_DEVICE_TYPE_INTEGRATED_GPU:
            return DeviceKind::kIntegratedGpu;
        case VK_PHYSICAL_DEVICE_TYPE_DISCRETE_GPU:
            return DeviceKind::kDiscreteGpu;
        case VK_PHYSICAL_DEVICE_TYPE_CPU:
            return DeviceKind::kSoftware;
        default:
            return DeviceKind::kOther;
    }
}

[[nodiscard]] DeviceInfo DescribeDevice(VkPhysicalDevice physical) {
    VkPhysicalDeviceFloatControlsProperties float_controls{
        .sType = VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_FLOAT_CONTROLS_PROPERTIES,
    };
    VkPhysicalDeviceDriverProperties driver{
        .sType = VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_DRIVER_PROPERTIES,
        .pNext = &float_controls,
    };
    VkPhysicalDeviceProperties2 properties2{
        .sType = VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_PROPERTIES_2,
        .pNext = &driver,
    };
    vkGetPhysicalDeviceProperties2(physical, &properties2);
    const VkPhysicalDeviceProperties& properties = properties2.properties;
    VkPhysicalDeviceFeatures features{};
    vkGetPhysicalDeviceFeatures(physical, &features);
    VkPhysicalDeviceMemoryProperties memory{};
    vkGetPhysicalDeviceMemoryProperties(physical, &memory);

    std::uint64_t device_local = 0;
    for (std::uint32_t i = 0; i < memory.memoryHeapCount; ++i) {
        if ((memory.memoryHeaps[i].flags & VK_MEMORY_HEAP_DEVICE_LOCAL_BIT) != 0) {
            device_local = std::max(device_local, memory.memoryHeaps[i].size);
        }
    }

    std::uint64_t render_memory = 0;
    if (const auto type = detail::VulkanHostMemoryType(memory, ~std::uint32_t{0})) {
        render_memory = memory.memoryHeaps[memory.memoryTypes[*type].heapIndex].size;
    }

    return DeviceInfo{
        .name = properties.deviceName,
        .driver_name = driver.driverName,
        .driver_info = driver.driverInfo,
        .kind = ClassifyDevice(properties),
        .vendor_id = properties.vendorID,
        .device_id = properties.deviceID,
        .api_version = properties.apiVersion,
        .driver_id = static_cast<std::uint32_t>(driver.driverID),
        .device_local_bytes = device_local,
        .render_memory_bytes = render_memory,
        .supports_fp64 = features.shaderFloat64 == VK_TRUE,
        .preserves_fp32_denormals = float_controls.shaderDenormPreserveFloat32 == VK_TRUE,
        .rounds_fp32_to_nearest = float_controls.shaderRoundingModeRTEFloat32 == VK_TRUE,
        .rounds_fp64_to_nearest = float_controls.shaderRoundingModeRTEFloat64 == VK_TRUE,
        .preserves_fp32_signed_zero_inf_nan =
            float_controls.shaderSignedZeroInfNanPreserveFloat32 == VK_TRUE,
    };
}

[[nodiscard]] Expected<void> RetainDozenThreadRuntime(const DeviceInfo& info) {
#if defined(__linux__)
    if (info.driver_id != static_cast<std::uint32_t>(VK_DRIVER_ID_MESA_DOZEN)) {
        return {};
    }

    // WSL's D3D12 runtime registers a pthread TLS destructor from
    // libd3d12core.so. Dozen normally dlcloses the libd3d12 wrapper and its
    // core with the Vulkan instance, which can leave a render worker calling
    // unmapped code as the thread exits. Process-lifetime references to both
    // mappings keep the registered destructor valid. The raw handles are
    // intentionally never dlclosed: the required lifetime is the process,
    // not a VulkanDevice instance.
    struct RuntimeLease {
        void* core = nullptr;
        void* wrapper = nullptr;
        std::string error;
    };
    static const RuntimeLease lease = [] {
        dlerror();
        void* core = dlopen("libd3d12core.so", RTLD_NOW | RTLD_LOCAL);
        const char* core_error = core == nullptr ? dlerror() : nullptr;
        if (core == nullptr) {
            return RuntimeLease{.error = core_error == nullptr ? "could not load libd3d12core.so"
                                                               : core_error};
        }
        dlerror();
        void* wrapper = dlopen("libd3d12.so", RTLD_NOW | RTLD_LOCAL);
        const char* wrapper_error = wrapper == nullptr ? dlerror() : nullptr;
        return RuntimeLease{
            .core = core,
            .wrapper = wrapper,
            .error = wrapper_error == nullptr ? "could not load libd3d12.so" : wrapper_error,
        };
    }();
    if (lease.core == nullptr || lease.wrapper == nullptr) {
        return Fail(ErrorDomain::kDevice, "retain Dozen D3D12 thread runtime", lease.error);
    }
#else
    (void)info;
#endif
    return {};
}

[[nodiscard]] Expected<VkInstance> CreateInstance() {
    const VkApplicationInfo app_info{
        .sType = VK_STRUCTURE_TYPE_APPLICATION_INFO,
        .pApplicationName = "sirius",
        .applicationVersion = 1,
        .pEngineName = "sirius",
        .engineVersion = 1,
        .apiVersion = kApiVersion,
    };
    const VkInstanceCreateInfo create_info{
        .sType = VK_STRUCTURE_TYPE_INSTANCE_CREATE_INFO,
        .pApplicationInfo = &app_info,
    };
    return detail::CreateInstanceWithPortability(create_info);
}

[[nodiscard]] Expected<std::vector<VkPhysicalDevice>> ListPhysicalDevices(VkInstance instance) {
    std::uint32_t count = 0;
    if (const VkResult r = vkEnumeratePhysicalDevices(instance, &count, nullptr); r != VK_SUCCESS) {
        return Fail(ErrorDomain::kDevice, "enumerate physical devices", VkResultText(r));
    }
    std::vector<VkPhysicalDevice> devices(count);
    if (count > 0) {
        if (const VkResult r = vkEnumeratePhysicalDevices(instance, &count, devices.data());
            r != VK_SUCCESS) {
            return Fail(ErrorDomain::kDevice, "enumerate physical devices", VkResultText(r));
        }
    }
    return devices;
}

}  // namespace

namespace detail {

Expected<VkInstance> CreateInstanceWithPortability(
    VkInstanceCreateInfo create_info, PFN_vkEnumerateInstanceExtensionProperties enumerate,
    PFN_vkCreateInstance create) {
    std::uint32_t count = 0;
    if (const VkResult r = enumerate(nullptr, &count, nullptr); r != VK_SUCCESS) {
        return Fail(ErrorDomain::kDevice, "enumerate Vulkan instance extensions", VkResultText(r));
    }
    std::vector<VkExtensionProperties> advertised(count);
    if (count > 0) {
        if (const VkResult r = enumerate(nullptr, &count, advertised.data()); r != VK_SUCCESS) {
            // A changing list (VK_INCOMPLETE) is not evidence of absence. Decline
            // instead of silently losing the portability opt-in.
            return Fail(ErrorDomain::kDevice, "enumerate Vulkan instance extensions",
                        VkResultText(r));
        }
        advertised.resize(count);
    }
    std::vector<const char*> enabled;
    for (std::uint32_t i = 0; i < create_info.enabledExtensionCount; ++i) {
        enabled.push_back(create_info.ppEnabledExtensionNames[i]);
    }
    if (std::any_of(advertised.begin(), advertised.end(), [](const auto& extension) {
            return std::strcmp(extension.extensionName,
                               VK_KHR_PORTABILITY_ENUMERATION_EXTENSION_NAME) == 0;
        })) {
        enabled.push_back(VK_KHR_PORTABILITY_ENUMERATION_EXTENSION_NAME);
        create_info.flags |= VK_INSTANCE_CREATE_ENUMERATE_PORTABILITY_BIT_KHR;
    }
    create_info.enabledExtensionCount = static_cast<std::uint32_t>(enabled.size());
    create_info.ppEnabledExtensionNames = enabled.empty() ? nullptr : enabled.data();
    VkInstance instance = VK_NULL_HANDLE;
    if (const VkResult r = create(&create_info, nullptr, &instance); r != VK_SUCCESS) {
        return Fail(ErrorDomain::kDevice, "create Vulkan instance", VkResultText(r));
    }
    return instance;
}

Expected<VkDevice> CreateDeviceWithPortability(VkPhysicalDevice physical,
                                               VkDeviceCreateInfo create_info,
                                               PFN_vkEnumerateDeviceExtensionProperties enumerate,
                                               PFN_vkCreateDevice create, DeviceInfo* device_info,
                                               PFN_vkGetPhysicalDeviceFeatures2 get_features) {
    if (device_info) device_info->fma_fp32_enabled = false;
    std::uint32_t count = 0;
    if (const VkResult r = enumerate(physical, nullptr, &count, nullptr); r != VK_SUCCESS) {
        return Fail(ErrorDomain::kDevice, "enumerate Vulkan device extensions", VkResultText(r));
    }
    std::vector<VkExtensionProperties> advertised(count);
    if (count > 0) {
        if (const VkResult r = enumerate(physical, nullptr, &count, advertised.data());
            r != VK_SUCCESS) {
            return Fail(ErrorDomain::kDevice, "enumerate Vulkan device extensions",
                        VkResultText(r));
        }
        advertised.resize(count);
    }
    std::vector<const char*> enabled;
    for (std::uint32_t i = 0; i < create_info.enabledExtensionCount; ++i) {
        enabled.push_back(create_info.ppEnabledExtensionNames[i]);
    }
    // The name is stable even in SDKs where the provisional subset structures
    // are hidden behind VK_ENABLE_BETA_EXTENSIONS. No subset features are used
    // by this storage-buffer compute adapter; Vulkan 1.3 meets its 1.1 dependency.
    constexpr const char* kPortabilitySubset = "VK_KHR_portability_subset";
    if (std::any_of(advertised.begin(), advertised.end(), [](const auto& extension) {
            return std::strcmp(extension.extensionName, kPortabilitySubset) == 0;
        })) {
        enabled.push_back(kPortabilitySubset);
    }
    ShaderFmaFeatures fma{.sType = kShaderFmaFeaturesType};
    constexpr const char* kShaderFma = "VK_KHR_shader_fma";
    const bool eligible =
        device_info && device_info->preserves_fp32_denormals &&
        device_info->rounds_fp32_to_nearest && device_info->preserves_fp32_signed_zero_inf_nan &&
        std::any_of(advertised.begin(), advertised.end(), [](const auto& extension) {
            return std::strcmp(extension.extensionName, kShaderFma) == 0;
        });
    bool enable_fma = false;
    if (eligible) {
        VkPhysicalDeviceFeatures2 features{.sType = VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_FEATURES_2,
                                           .pNext = &fma};
        get_features(physical, &features);
        enable_fma = fma.shaderFmaFloat32 == VK_TRUE;
        if (enable_fma) {
            // Preserve the caller's chain and legacy Float64 admission. Only
            // binary32 FMA is needed by the optional retained Transport module.
            fma.pNext = const_cast<void*>(create_info.pNext);
            fma.shaderFmaFloat16 = VK_FALSE;
            fma.shaderFmaFloat64 = VK_FALSE;
            create_info.pNext = &fma;
            enabled.push_back(kShaderFma);
        }
    }
    create_info.enabledExtensionCount = static_cast<std::uint32_t>(enabled.size());
    create_info.ppEnabledExtensionNames = enabled.empty() ? nullptr : enabled.data();
    VkDevice device = VK_NULL_HANDLE;
    if (const VkResult r = create(physical, &create_info, nullptr, &device); r != VK_SUCCESS) {
        return Fail(ErrorDomain::kDevice, "create Vulkan logical device", VkResultText(r));
    }
    if (device_info) device_info->fma_fp32_enabled = enable_fma;
    return device;
}

}  // namespace detail

std::optional<std::uint32_t> detail::VulkanHostMemoryType(
    const VkPhysicalDeviceMemoryProperties& properties, std::uint32_t memory_type_bits) {
    constexpr VkMemoryPropertyFlags kRequired =
        VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT;
    std::optional<std::uint32_t> selected;
    VkDeviceSize selected_heap_size = 0;
    for (std::uint32_t i = 0; i < properties.memoryTypeCount; ++i) {
        const auto& type = properties.memoryTypes[i];
        const auto heap_size = properties.memoryHeaps[type.heapIndex].size;
        if ((memory_type_bits & (1u << i)) != 0 && (type.propertyFlags & kRequired) == kRequired &&
            heap_size > selected_heap_size) {
            selected = i;
            selected_heap_size = heap_size;
        }
    }
    if (!selected) return std::nullopt;

    const auto& original = properties.memoryTypes[*selected];
    const auto preferred_flags = original.propertyFlags | VK_MEMORY_PROPERTY_HOST_CACHED_BIT;
    // Uncached coherent memory is commonly write-combined; reading the retained
    // output then costs much more than cached host memory. Keep the original
    // heap/capacity, and add only host caching, never feature-dependent flags.
    for (std::uint32_t i = 0; i < properties.memoryTypeCount; ++i) {
        const auto& type = properties.memoryTypes[i];
        if ((memory_type_bits & (1u << i)) != 0 && type.heapIndex == original.heapIndex &&
            type.propertyFlags == preferred_flags) {
            return i;
        }
    }
    return selected;
}

bool detail::VulkanPipelineCacheDataCompatible(std::span<const std::byte> data,
                                               const VkPhysicalDeviceProperties& properties) {
    // Version-one fields are tightly packed little-endian bytes, independent of
    // the host's byte order or C structure packing. Only Vulkan-exported bytes
    // reach the driver; this guard discards incompatible optimization data.
    // https://docs.vulkan.org/refpages/latest/refpages/source/VkPipelineCacheHeaderVersionOne.html
    constexpr std::size_t kHeaderBytes = 32;
    if (data.size() < kHeaderBytes || data.size() > kVulkanPipelineCacheBlobLimit) return false;
    const auto word = [&](std::size_t offset) {
        std::uint32_t value = 0;
        for (std::size_t byte = 0; byte < 4; ++byte)
            value |= std::uint32_t(std::to_integer<unsigned char>(data[offset + byte]))
                     << (8 * byte);
        return value;
    };
    return word(0) == kHeaderBytes &&
           word(4) == static_cast<std::uint32_t>(VK_PIPELINE_CACHE_HEADER_VERSION_ONE) &&
           word(8) == properties.vendorID && word(12) == properties.deviceID &&
           std::memcmp(data.data() + 16, properties.pipelineCacheUUID, VK_UUID_SIZE) == 0;
}

Expected<std::vector<DeviceInfo>> EnumerateVulkanDevices() {
    auto instance = CreateInstance();
    if (!instance) {
        // No loader or no ICD is a decline: the caller falls back to the
        // CPU tracer with a clear message, never a fabricated device.
        return std::vector<DeviceInfo>{};
    }
    auto devices = ListPhysicalDevices(*instance);
    if (!devices) {
        vkDestroyInstance(*instance, nullptr);
        return std::unexpected(devices.error());
    }
    std::vector<DeviceInfo> infos;
    infos.reserve(devices->size());
    for (VkPhysicalDevice physical : *devices) {
        infos.push_back(DescribeDevice(physical));
    }
    vkDestroyInstance(*instance, nullptr);
    return infos;
}

Expected<std::size_t> ResolveVulkanDeviceIndex(std::span<const DeviceInfo> devices) {
    const char* raw = std::getenv("SIRIUS_VULKAN_DEVICE");
    if (raw == nullptr || *raw == '\0') {
        if (devices.empty()) {
            return Fail(ErrorDomain::kDevice, "select Vulkan device", "no devices were enumerated");
        }
        return std::size_t{0};
    }

    const std::string_view value(raw);
    std::size_t index = 0;
    const auto [end, error] = std::from_chars(value.data(), value.data() + value.size(), index);
    if (error != std::errc{} || end != value.data() + value.size()) {
        return Fail(ErrorDomain::kConfiguration, "select Vulkan device",
                    std::format("SIRIUS_VULKAN_DEVICE='{}' is not a zero-based integer", value));
    }
    if (index >= devices.size()) {
        return Fail(ErrorDomain::kConfiguration, "select Vulkan device",
                    std::format("SIRIUS_VULKAN_DEVICE={} is out of range for {} device(s)", index,
                                devices.size()));
    }
    return index;
}

Expected<std::unique_ptr<ComputeDevice>> CreateVulkanDevice(std::size_t index) {
    auto instance = CreateInstance();
    if (!instance) {
        return std::unexpected(instance.error());
    }
    auto device = std::make_unique<VulkanDevice>();
    device->instance_ = *instance;

    auto physicals = ListPhysicalDevices(device->instance_);
    if (!physicals) {
        return std::unexpected(physicals.error());
    }
    if (index >= physicals->size()) {
        return Fail(ErrorDomain::kDevice, "open Vulkan device",
                    std::format("index {} out of range ({} devices)", index, physicals->size()));
    }
    device->physical_ = (*physicals)[index];
    device->info_ = DescribeDevice(device->physical_);
    if (auto retained = RetainDozenThreadRuntime(device->info_); !retained) {
        return std::unexpected(retained.error());
    }

    // Queue family with compute.
    std::uint32_t family_count = 0;
    vkGetPhysicalDeviceQueueFamilyProperties(device->physical_, &family_count, nullptr);
    std::vector<VkQueueFamilyProperties> families(family_count);
    vkGetPhysicalDeviceQueueFamilyProperties(device->physical_, &family_count, families.data());
    std::uint32_t family = family_count;
    for (std::uint32_t i = 0; i < family_count; ++i) {
        if ((families[i].queueFlags & VK_QUEUE_COMPUTE_BIT) != 0) {
            family = i;
            break;
        }
    }
    if (family == family_count) {
        return Fail(ErrorDomain::kDevice, "open Vulkan device",
                    std::format("'{}' has no compute queue", device->info_.name));
    }
    device->queue_family_ = family;

    const float priority = 1.0f;
    const VkDeviceQueueCreateInfo queue_info{
        .sType = VK_STRUCTURE_TYPE_DEVICE_QUEUE_CREATE_INFO,
        .queueFamilyIndex = family,
        .queueCount = 1,
        .pQueuePriorities = &priority,
    };
    // fp64 is enabled whenever the hardware offers it so the precision
    // ladder can select double kernels at run time.
    VkPhysicalDeviceFeatures enabled{};
    enabled.shaderFloat64 = device->info_.supports_fp64 ? VK_TRUE : VK_FALSE;
    const VkDeviceCreateInfo device_info{
        .sType = VK_STRUCTURE_TYPE_DEVICE_CREATE_INFO,
        .queueCreateInfoCount = 1,
        .pQueueCreateInfos = &queue_info,
        .pEnabledFeatures = &enabled,
    };
    auto logical = detail::CreateDeviceWithPortability(device->physical_, device_info,
                                                       vkEnumerateDeviceExtensionProperties,
                                                       vkCreateDevice, &device->info_);
    if (!logical) {
        return std::unexpected(logical.error());
    }
    device->device_ = *logical;
    vkGetDeviceQueue(device->device_, family, 0, &device->queue_);

    device->InitialisePipelineCache();

    const VkCommandPoolCreateInfo pool_info{
        .sType = VK_STRUCTURE_TYPE_COMMAND_POOL_CREATE_INFO,
        .flags = VK_COMMAND_POOL_CREATE_RESET_COMMAND_BUFFER_BIT,
        .queueFamilyIndex = family,
    };
    if (const VkResult r =
            vkCreateCommandPool(device->device_, &pool_info, nullptr, &device->command_pool_);
        r != VK_SUCCESS) {
        return Fail(ErrorDomain::kDevice, "create command pool", VkResultText(r));
    }

    // Sized for tile-grained dispatch: a handful of sets in flight, each
    // with a handful of bindings; revisited with the memory governor.
    const VkDescriptorPoolSize pool_sizes[] = {
        {.type = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, .descriptorCount = 256},
        {.type = VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER, .descriptorCount = 64},
    };
    const VkDescriptorPoolCreateInfo descriptor_pool_info{
        .sType = VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO,
        .flags = VK_DESCRIPTOR_POOL_CREATE_FREE_DESCRIPTOR_SET_BIT,
        .maxSets = 64,
        .poolSizeCount = 2,
        .pPoolSizes = pool_sizes,
    };
    if (const VkResult r = vkCreateDescriptorPool(device->device_, &descriptor_pool_info, nullptr,
                                                  &device->descriptor_pool_);
        r != VK_SUCCESS) {
        return Fail(ErrorDomain::kDevice, "create descriptor pool", VkResultText(r));
    }

    return Expected<std::unique_ptr<ComputeDevice>>{std::move(device)};
}

VulkanDevice::~VulkanDevice() {
    if (device_ != VK_NULL_HANDLE) {
        vkDeviceWaitIdle(device_);
        // Export while the owning device is alive and no operations can modify
        // its cache. A failed optimization snapshot must not escape a destructor.
        if (pipeline_cache_modified_) {
            try {
                (void)SnapshotPipelineCache();
            } catch (const std::bad_alloc&) {
            }
        }
        for (auto& [key, pipeline] : pipelines_) {
            vkDestroyPipeline(device_, pipeline.pipeline, nullptr);
            vkDestroyPipelineLayout(device_, pipeline.layout, nullptr);
            vkDestroyDescriptorSetLayout(device_, pipeline.set_layout, nullptr);
        }
        for (const Kernel& kernel : kernels_) {
            vkDestroyShaderModule(device_, kernel.module, nullptr);
        }
        for (Buffer& buffer : buffers_) {
            vkDestroyBuffer(device_, buffer.buffer, nullptr);
            vkFreeMemory(device_, buffer.memory, nullptr);
        }
        if (pipeline_cache_ != VK_NULL_HANDLE)
            vkDestroyPipelineCache(device_, pipeline_cache_, nullptr);
        vkDestroyDescriptorPool(device_, descriptor_pool_, nullptr);
        vkDestroyCommandPool(device_, command_pool_, nullptr);
        vkDestroyDevice(device_, nullptr);
    }
    if (instance_ != VK_NULL_HANDLE) {
        vkDestroyInstance(instance_, nullptr);
    }
}

void VulkanDevice::InitialisePipelineCache() {
    vkGetPhysicalDeviceProperties(physical_, &pipeline_cache_properties_);
    std::vector<std::byte> seed;
    try {
        auto& stored = PipelineCacheStore();
        std::lock_guard lock(stored.mutex);
        if (detail::VulkanPipelineCacheDataCompatible(stored.data, pipeline_cache_properties_))
            seed = stored.data;
        else if (!stored.data.empty())
            pipeline_cache_stats_.import_discarded = true;
    } catch (const std::bad_alloc&) {
        // Import is optional; an empty cache remains usable if the host cannot
        // afford the bounded serialized copy.
        seed.clear();
        pipeline_cache_stats_.import_discarded = true;
    }
    // Default flags preserve Vulkan's internally synchronized pipeline creation.
    // https://docs.vulkan.org/refpages/latest/refpages/source/vkCreatePipelineCache.html
    VkPipelineCacheCreateInfo cache_info{
        .sType = VK_STRUCTURE_TYPE_PIPELINE_CACHE_CREATE_INFO,
        .initialDataSize = seed.size(),
        .pInitialData = seed.empty() ? nullptr : seed.data(),
    };
    auto result = vkCreatePipelineCache(device_, &cache_info, nullptr, &pipeline_cache_);
    if (result != VK_SUCCESS && !seed.empty()) {
        // Cache import cannot add a product admission requirement. Retry empty,
        // then retain the original uncached pipeline path if allocation fails.
        pipeline_cache_ = VK_NULL_HANDLE;
        seed.clear();
        pipeline_cache_stats_.import_discarded = true;
        cache_info.initialDataSize = 0;
        cache_info.pInitialData = nullptr;
        result = vkCreatePipelineCache(device_, &cache_info, nullptr, &pipeline_cache_);
    }
    pipeline_cache_stats_.creation_result = result;
    pipeline_cache_stats_.enabled = result == VK_SUCCESS;
    if (result == VK_SUCCESS)
        pipeline_cache_stats_.imported_bytes = seed.size();
    else
        pipeline_cache_ = VK_NULL_HANDLE;
}

Expected<VulkanPipelineCacheStats> VulkanDevice::SnapshotPipelineCache() {
    if (pipeline_cache_ == VK_NULL_HANDLE)
        return Fail(ErrorDomain::kDevice, "snapshot pipeline cache", "cache is not initialized");
    std::size_t available = 0;
    if (const auto result = vkGetPipelineCacheData(device_, pipeline_cache_, &available, nullptr);
        result != VK_SUCCESS)
        return Fail(ErrorDomain::kDevice, "size pipeline cache", VkResultText(result));
    pipeline_cache_stats_.available_bytes = available;
    if (!pipeline_cache_modified_) return pipeline_cache_stats_;
    pipeline_cache_stats_.exported_bytes = 0;
    pipeline_cache_stats_.export_discarded = false;
    const auto discard = [&] {
        pipeline_cache_stats_.export_discarded = true;
        pipeline_cache_modified_ = false;
        return pipeline_cache_stats_;
    };
    // Never truncate a blob or allocate an unbounded serialization buffer. The
    // driver's internal cache is separate from this bound and explicit buffers.
    if (available < 32 || available > detail::kVulkanPipelineCacheBlobLimit) return discard();
    std::vector<std::byte> data;
    try {
        data.resize(available);
    } catch (const std::bad_alloc&) {
        return discard();
    }
    std::size_t written = available;
    const auto result = vkGetPipelineCacheData(device_, pipeline_cache_, &written, data.data());
    if (result == VK_INCOMPLETE || written > data.size()) return discard();
    if (result != VK_SUCCESS)
        return Fail(ErrorDomain::kDevice, "export pipeline cache", VkResultText(result));
    data.resize(written);
    if (!detail::VulkanPipelineCacheDataCompatible(data, pipeline_cache_properties_))
        return discard();
    auto& stored = PipelineCacheStore();
    std::lock_guard lock(stored.mutex);
    stored.data = std::move(data);
    pipeline_cache_stats_.exported_bytes = written;
    pipeline_cache_modified_ = false;
    return pipeline_cache_stats_;
}

Expected<void> ValidateVulkanKernelPrecision(std::span<const std::uint32_t> spirv,
                                             bool supports_fp64, bool fma_fp32_enabled) {
    // SPIR-V binary encoding: five-word header, high 16 bits instruction word
    // count, low 16 bits opcode; OpCapability=17 and Float64=10. Khronos authority:
    // https://github.com/KhronosGroup/SPIRV-Headers/blob/main/include/spirv/unified1/spirv.hpp11
    constexpr std::uint32_t kMagic = 0x07230203u;
    constexpr std::uint32_t kOpCapability = 17u;
    constexpr std::uint32_t kFloat64 = 10u;
    constexpr std::uint32_t kFmaKHR = 6030u;
    constexpr std::uint32_t kOpFmaKHR = 4427u;
    if (spirv.size() <= 5 || spirv[0] != kMagic || spirv[4] != 0) {
        return Fail(ErrorDomain::kKernel, "validate shader module", "malformed SPIR-V header");
    }
    bool needs_fp64 = false;
    bool needs_fma = false;
    std::map<std::uint32_t, std::uint32_t> float_widths;
    std::vector<std::uint32_t> fma_types;
    for (std::size_t offset = 5; offset < spirv.size();) {
        const std::uint32_t word_count = spirv[offset] >> 16;
        const std::uint32_t opcode = spirv[offset] & 0xffffu;
        if (word_count == 0 || word_count > spirv.size() - offset) {
            return Fail(ErrorDomain::kKernel, "validate shader module",
                        "malformed SPIR-V instruction extent");
        }
        if (opcode == kOpCapability) {
            if (word_count != 2) {
                return Fail(ErrorDomain::kKernel, "validate shader module",
                            "malformed SPIR-V OpCapability");
            }
            needs_fp64 = needs_fp64 || spirv[offset + 1] == kFloat64;
            needs_fma = needs_fma || spirv[offset + 1] == kFmaKHR;
        }
        needs_fma = needs_fma || opcode == kOpFmaKHR;
        if (opcode == 22u && word_count == 3)
            float_widths.emplace(spirv[offset + 1], spirv[offset + 2]);
        if (opcode == kOpFmaKHR) {
            if (word_count != 6)
                return Fail(ErrorDomain::kKernel, "validate shader module",
                            "malformed SPIR-V OpFmaKHR");
            fma_types.push_back(spirv[offset + 1]);
        }
        offset += word_count;
    }
    if (needs_fp64 && !supports_fp64) {
        return Fail(ErrorDomain::kKernel, "load shader precision",
                    "SPIR-V Float64 requested but the device lacks shaderFloat64");
    }
    if (needs_fma && !fma_fp32_enabled)
        return Fail(ErrorDomain::kKernel, "load shader precision",
                    "SPIR-V FMAKHR requested but the logical device has not enabled "
                    "shaderFmaFloat32 with required controls");
    for (const auto type : fma_types) {
        const auto width = float_widths.find(type);
        if (width == float_widths.end() || width->second != 32u)
            return Fail(ErrorDomain::kKernel, "load shader precision",
                        "only binary32 OpFmaKHR is enabled on this logical device");
    }
    return {};
}

Expected<KernelHandle> VulkanDevice::LoadKernel(std::span<const std::uint32_t> spirv) {
    if (auto precision =
            ValidateVulkanKernelPrecision(spirv, info_.supports_fp64, info_.fma_fp32_enabled);
        !precision) {
        return std::unexpected(precision.error());
    }
    // Identical modules share their device-lived pipelines even when a new
    // retained compute instance owns different input and output buffers.
    for (std::size_t index = 0; index < kernels_.size(); ++index) {
        if (std::ranges::equal(kernels_[index].words, spirv))
            return KernelHandle{static_cast<std::uint32_t>(index)};
    }
    Kernel kernel;
    try {
        kernel.words.assign(spirv.begin(), spirv.end());
        // Complete host allocations before acquiring the Vulkan object.
        kernels_.reserve(kernels_.size() + 1);
    } catch (const std::bad_alloc&) {
        return Fail(ErrorDomain::kKernel, "load shader module",
                    "host shader storage allocation failed");
    }
    const VkShaderModuleCreateInfo create_info{
        .sType = VK_STRUCTURE_TYPE_SHADER_MODULE_CREATE_INFO,
        .codeSize = kernel.words.size() * sizeof(std::uint32_t),
        .pCode = kernel.words.data(),
    };
    if (const VkResult r = vkCreateShaderModule(device_, &create_info, nullptr, &kernel.module);
        r != VK_SUCCESS) {
        return Fail(ErrorDomain::kKernel, "create shader module", VkResultText(r));
    }
    kernels_.push_back(std::move(kernel));
    return KernelHandle{static_cast<std::uint32_t>(kernels_.size() - 1)};
}

Expected<void> VulkanDevice::SetBufferAllocationLimit(std::uint64_t bytes) {
    if (bytes < buffer_allocation_bytes_) {
        return Fail(ErrorDomain::kDevice, "set buffer allocation limit",
                    "limit is below the actual resident buffer allocations");
    }
    buffer_allocation_limit_ = bytes;
    return {};
}

Expected<BufferHandle> VulkanDevice::CreateBuffer(std::uint64_t size_bytes, BufferUsage usage) {
    SIRIUS_PRE(size_bytes > 0);
    if (size_bytes > buffer_allocation_limit_ - buffer_allocation_bytes_) {
        return Fail(ErrorDomain::kDevice, "allocate buffer memory",
                    "requested buffer exceeds the remaining explicit allocation budget");
    }
    const auto buffer_info = BufferCreateInfo(size_bytes, usage);
    Buffer buffer{.size_bytes = size_bytes, .usage = usage};
    if (const VkResult r = vkCreateBuffer(device_, &buffer_info, nullptr, &buffer.buffer);
        r != VK_SUCCESS) {
        return Fail(ErrorDomain::kDevice, "create buffer", VkResultText(r));
    }

    VkMemoryRequirements requirements{};
    vkGetBufferMemoryRequirements(device_, buffer.buffer, &requirements);
    if (requirements.size > buffer_allocation_limit_ - buffer_allocation_bytes_) {
        vkDestroyBuffer(device_, buffer.buffer, nullptr);
        return Fail(ErrorDomain::kDevice, "allocate buffer memory",
                    "actual Vulkan allocation exceeds the remaining explicit allocation budget");
    }
    VkPhysicalDeviceMemoryProperties memory_properties{};
    vkGetPhysicalDeviceMemoryProperties(physical_, &memory_properties);

    const auto type_index =
        detail::VulkanHostMemoryType(memory_properties, requirements.memoryTypeBits);
    if (!type_index) {
        vkDestroyBuffer(device_, buffer.buffer, nullptr);
        return Fail(ErrorDomain::kDevice, "allocate buffer memory",
                    "no host-visible coherent memory type");
    }

    const VkMemoryAllocateInfo allocate_info{
        .sType = VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO,
        .allocationSize = requirements.size,
        .memoryTypeIndex = *type_index,
    };
    if (const VkResult r = vkAllocateMemory(device_, &allocate_info, nullptr, &buffer.memory);
        r != VK_SUCCESS) {
        vkDestroyBuffer(device_, buffer.buffer, nullptr);
        return Fail(ErrorDomain::kDevice, "allocate buffer memory", VkResultText(r));
    }
    if (const VkResult r = vkBindBufferMemory(device_, buffer.buffer, buffer.memory, 0);
        r != VK_SUCCESS) {
        vkDestroyBuffer(device_, buffer.buffer, nullptr);
        vkFreeMemory(device_, buffer.memory, nullptr);
        return Fail(ErrorDomain::kDevice, "bind buffer memory", VkResultText(r));
    }

    buffers_.push_back(buffer);
    buffer_allocation_bytes_ += requirements.size;
    return BufferHandle{static_cast<std::uint32_t>(buffers_.size() - 1)};
}

Expected<std::uint64_t> VulkanDevice::RequiredBufferAllocationBytes(std::uint64_t size_bytes,
                                                                    BufferUsage usage) {
    if (size_bytes == 0)
        return Fail(ErrorDomain::kDevice, "query buffer allocation",
                    "buffer size must be positive");
    const auto buffer_info = BufferCreateInfo(size_bytes, usage);
    VkBuffer buffer = VK_NULL_HANDLE;
    if (const VkResult r = vkCreateBuffer(device_, &buffer_info, nullptr, &buffer); r != VK_SUCCESS)
        return Fail(ErrorDomain::kDevice, "query buffer allocation", VkResultText(r));
    VkMemoryRequirements requirements{};
    vkGetBufferMemoryRequirements(device_, buffer, &requirements);
    vkDestroyBuffer(device_, buffer, nullptr);
    if (requirements.size < size_bytes)
        return Fail(ErrorDomain::kDevice, "query buffer allocation",
                    "driver allocation requirement is smaller than the requested buffer");
    return requirements.size;
}

Expected<void> VulkanDevice::WriteBuffer(BufferHandle handle, std::span<const std::byte> data) {
    SIRIUS_PRE(handle.value < buffers_.size());
    Buffer& buffer = buffers_[handle.value];
    SIRIUS_PRE(data.size_bytes() <= buffer.size_bytes);
    void* mapped = nullptr;
    if (const VkResult r = vkMapMemory(device_, buffer.memory, 0, data.size_bytes(), 0, &mapped);
        r != VK_SUCCESS) {
        return Fail(ErrorDomain::kDevice, "map buffer for write", VkResultText(r));
    }
    std::memcpy(mapped, data.data(), data.size_bytes());
    vkUnmapMemory(device_, buffer.memory);
    return {};
}

Expected<void> VulkanDevice::ReadBuffer(BufferHandle handle, std::span<std::byte> out) {
    SIRIUS_PRE(handle.value < buffers_.size());
    Buffer& buffer = buffers_[handle.value];
    SIRIUS_PRE(out.size_bytes() <= buffer.size_bytes);
    void* mapped = nullptr;
    if (const VkResult r = vkMapMemory(device_, buffer.memory, 0, out.size_bytes(), 0, &mapped);
        r != VK_SUCCESS) {
        return Fail(ErrorDomain::kDevice, "map buffer for read", VkResultText(r));
    }
    std::memcpy(out.data(), mapped, out.size_bytes());
    vkUnmapMemory(device_, buffer.memory);
    return {};
}

Expected<VulkanDevice::Pipeline*> VulkanDevice::GetOrCreatePipeline(
    KernelHandle kernel, std::span<const BufferHandle> buffers, bool* created) {
    if (created != nullptr) *created = false;
    PipelineKey key{.kernel = kernel.value, .bindings = {}};
    key.bindings.reserve(buffers.size());
    for (const BufferHandle handle : buffers) {
        key.bindings.push_back(buffers_[handle.value].usage);
    }
    if (const auto found = pipelines_.find(key); found != pipelines_.end()) {
        return &found->second;
    }

    std::vector<VkDescriptorSetLayoutBinding> bindings;
    bindings.reserve(buffers.size());
    for (std::uint32_t i = 0; i < buffers.size(); ++i) {
        bindings.push_back(VkDescriptorSetLayoutBinding{
            .binding = i,
            .descriptorType = key.bindings[i] == BufferUsage::kStorage
                                  ? VK_DESCRIPTOR_TYPE_STORAGE_BUFFER
                                  : VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER,
            .descriptorCount = 1,
            .stageFlags = VK_SHADER_STAGE_COMPUTE_BIT,
        });
    }

    Pipeline pipeline{};
    const VkDescriptorSetLayoutCreateInfo layout_info{
        .sType = VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_CREATE_INFO,
        .bindingCount = static_cast<std::uint32_t>(bindings.size()),
        .pBindings = bindings.data(),
    };
    if (const VkResult r =
            vkCreateDescriptorSetLayout(device_, &layout_info, nullptr, &pipeline.set_layout);
        r != VK_SUCCESS) {
        return Fail(ErrorDomain::kKernel, "create descriptor set layout", VkResultText(r));
    }
    const VkPipelineLayoutCreateInfo pipeline_layout_info{
        .sType = VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO,
        .setLayoutCount = 1,
        .pSetLayouts = &pipeline.set_layout,
    };
    if (const VkResult r =
            vkCreatePipelineLayout(device_, &pipeline_layout_info, nullptr, &pipeline.layout);
        r != VK_SUCCESS) {
        vkDestroyDescriptorSetLayout(device_, pipeline.set_layout, nullptr);
        return Fail(ErrorDomain::kKernel, "create pipeline layout", VkResultText(r));
    }
    // slangc names every SPIR-V compute entry point "main".
    const VkComputePipelineCreateInfo pipeline_info{
        .sType = VK_STRUCTURE_TYPE_COMPUTE_PIPELINE_CREATE_INFO,
        .stage =
            VkPipelineShaderStageCreateInfo{
                .sType = VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO,
                .stage = VK_SHADER_STAGE_COMPUTE_BIT,
                .module = kernels_[kernel.value].module,
                .pName = "main",
            },
        .layout = pipeline.layout,
    };
    if (const VkResult r = vkCreateComputePipelines(device_, pipeline_cache_, 1, &pipeline_info,
                                                    nullptr, &pipeline.pipeline);
        r != VK_SUCCESS) {
        vkDestroyPipelineLayout(device_, pipeline.layout, nullptr);
        vkDestroyDescriptorSetLayout(device_, pipeline.set_layout, nullptr);
        return Fail(ErrorDomain::kKernel, "create compute pipeline", VkResultText(r));
    }

    pipeline_cache_modified_ = pipeline_cache_ != VK_NULL_HANDLE;

    auto [inserted, _] = pipelines_.emplace(std::move(key), pipeline);
    if (created != nullptr) *created = true;
    return &inserted->second;
}

Expected<void> VulkanDevice::Dispatch(KernelHandle kernel, std::span<const BufferHandle> buffers,
                                      std::uint32_t groups_x, std::uint32_t groups_y,
                                      std::uint32_t groups_z, DispatchTiming* timing) {
    const ComputeDispatch command{kernel, buffers, groups_x, groups_y, groups_z};
    return DispatchCommands(std::span(&command, 1), timing, nullptr);
}

Expected<void> VulkanDevice::DispatchIndependentPair(const std::array<ComputeDispatch, 2>& commands,
                                                     IndependentPairTiming* timing) {
    if (timing) *timing = {};
    for (const auto first : commands[0].buffers)
        for (const auto second : commands[1].buffers)
            if (first.value == second.value)
                return Fail(ErrorDomain::kDevice, "dispatch independent compute pair",
                            "commands share a bound buffer");
    return DispatchCommands(commands, timing ? &timing->combined : nullptr,
                            timing ? &timing->pipeline_creations : nullptr);
}

Expected<void> VulkanDevice::DispatchCommands(std::span<const ComputeDispatch> commands,
                                              DispatchTiming* timing,
                                              std::uint32_t* pipeline_creations) {
    SIRIUS_PRE(!commands.empty() && commands.size() <= 2);
    for (const auto& item : commands) {
        SIRIUS_PRE(item.kernel.value < kernels_.size());
        SIRIUS_PRE(item.groups_x > 0 && item.groups_y > 0 && item.groups_z > 0);
        for (const BufferHandle handle : item.buffers) SIRIUS_PRE(handle.value < buffers_.size());
    }

    if (timing != nullptr) *timing = {};
    if (pipeline_creations) *pipeline_creations = 0;
    const auto dispatch_start = std::chrono::steady_clock::now();
    std::array<Pipeline*, 2> pipelines{};
    std::array<VkDescriptorSetLayout, 2> layouts{};
    for (std::size_t i = 0; i < commands.size(); ++i) {
        bool created = false;
        auto pipeline = GetOrCreatePipeline(commands[i].kernel, commands[i].buffers, &created);
        if (timing) timing->pipeline_created = timing->pipeline_created || created;
        if (pipeline_creations) *pipeline_creations += created ? 1 : 0;
        if (!pipeline) return std::unexpected(pipeline.error());
        pipelines[i] = *pipeline;
        layouts[i] = (*pipeline)->set_layout;
    }
    const auto pipeline_end = std::chrono::steady_clock::now();

    const VkDescriptorSetAllocateInfo set_info{
        .sType = VK_STRUCTURE_TYPE_DESCRIPTOR_SET_ALLOCATE_INFO,
        .descriptorPool = descriptor_pool_,
        .descriptorSetCount = static_cast<std::uint32_t>(commands.size()),
        .pSetLayouts = layouts.data(),
    };
    std::array<VkDescriptorSet, 2> sets{};
    if (const VkResult r = vkAllocateDescriptorSets(device_, &set_info, sets.data());
        r != VK_SUCCESS) {
        return Fail(ErrorDomain::kDevice, "allocate descriptor set", VkResultText(r));
    }
    const auto free_sets = [&] {
        vkFreeDescriptorSets(device_, descriptor_pool_, set_info.descriptorSetCount, sets.data());
    };

    for (std::size_t command_index = 0; command_index < commands.size(); ++command_index) {
        const auto buffers = commands[command_index].buffers;
        std::vector<VkDescriptorBufferInfo> buffer_infos(buffers.size());
        std::vector<VkWriteDescriptorSet> writes(buffers.size());
        for (std::uint32_t i = 0; i < buffers.size(); ++i) {
            const Buffer& buffer = buffers_[buffers[i].value];
            buffer_infos[i] = VkDescriptorBufferInfo{
                .buffer = buffer.buffer,
                .offset = 0,
                .range = buffer.size_bytes,
            };
            writes[i] = VkWriteDescriptorSet{
                .sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET,
                .dstSet = sets[command_index],
                .dstBinding = i,
                .descriptorCount = 1,
                .descriptorType = buffer.usage == BufferUsage::kStorage
                                      ? VK_DESCRIPTOR_TYPE_STORAGE_BUFFER
                                      : VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER,
                .pBufferInfo = &buffer_infos[i],
            };
        }
        vkUpdateDescriptorSets(device_, static_cast<std::uint32_t>(writes.size()), writes.data(), 0,
                               nullptr);
    }

    const VkCommandBufferAllocateInfo command_info{
        .sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_ALLOCATE_INFO,
        .commandPool = command_pool_,
        .level = VK_COMMAND_BUFFER_LEVEL_PRIMARY,
        .commandBufferCount = 1,
    };
    VkCommandBuffer command = VK_NULL_HANDLE;
    if (const VkResult r = vkAllocateCommandBuffers(device_, &command_info, &command);
        r != VK_SUCCESS) {
        free_sets();
        return Fail(ErrorDomain::kDevice, "allocate command buffer", VkResultText(r));
    }

    const VkCommandBufferBeginInfo begin_info{
        .sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO,
        .flags = VK_COMMAND_BUFFER_USAGE_ONE_TIME_SUBMIT_BIT,
    };
    if (const auto result = vkBeginCommandBuffer(command, &begin_info); result != VK_SUCCESS) {
        vkFreeCommandBuffers(device_, command_pool_, 1, &command);
        free_sets();
        return Fail(ErrorDomain::kDevice, "begin compute command buffer", VkResultText(result));
    }
    for (std::size_t i = 0; i < commands.size(); ++i) {
        const auto& item = commands[i];
        const VkMemoryBarrier before{
            .sType = VK_STRUCTURE_TYPE_MEMORY_BARRIER,
            .srcAccessMask = VK_ACCESS_SHADER_WRITE_BIT | VK_ACCESS_HOST_WRITE_BIT,
            .dstAccessMask = VK_ACCESS_SHADER_READ_BIT | VK_ACCESS_SHADER_WRITE_BIT,
        };
        vkCmdPipelineBarrier(
            command, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT | VK_PIPELINE_STAGE_HOST_BIT,
            VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 1, &before, 0, nullptr, 0, nullptr);
        vkCmdBindPipeline(command, VK_PIPELINE_BIND_POINT_COMPUTE, pipelines[i]->pipeline);
        vkCmdBindDescriptorSets(command, VK_PIPELINE_BIND_POINT_COMPUTE, pipelines[i]->layout, 0, 1,
                                &sets[i], 0, nullptr);
        vkCmdDispatch(command, item.groups_x, item.groups_y, item.groups_z);
        const VkMemoryBarrier after{
            .sType = VK_STRUCTURE_TYPE_MEMORY_BARRIER,
            .srcAccessMask = VK_ACCESS_SHADER_WRITE_BIT,
            .dstAccessMask = VK_ACCESS_SHADER_READ_BIT | VK_ACCESS_HOST_READ_BIT,
        };
        vkCmdPipelineBarrier(command, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,
                             VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT | VK_PIPELINE_STAGE_HOST_BIT, 0,
                             1, &after, 0, nullptr, 0, nullptr);
    }
    if (const auto result = vkEndCommandBuffer(command); result != VK_SUCCESS) {
        vkFreeCommandBuffers(device_, command_pool_, 1, &command);
        free_sets();
        return Fail(ErrorDomain::kDevice, "end compute command buffer", VkResultText(result));
    }

    const VkSubmitInfo submit_info{
        .sType = VK_STRUCTURE_TYPE_SUBMIT_INFO,
        .commandBufferCount = 1,
        .pCommandBuffers = &command,
    };
    const auto submit_start = std::chrono::steady_clock::now();
    VkResult submit_result = vkQueueSubmit(queue_, 1, &submit_info, VK_NULL_HANDLE);
    if (submit_result == VK_SUCCESS) {
        submit_result = vkQueueWaitIdle(queue_);
    }

    const auto submit_end = std::chrono::steady_clock::now();
    vkFreeCommandBuffers(device_, command_pool_, 1, &command);
    free_sets();
    const auto dispatch_end = std::chrono::steady_clock::now();
    if (timing != nullptr) {
        const auto milliseconds = [](auto duration) {
            return std::chrono::duration<double, std::milli>(duration).count();
        };
        timing->pipeline_setup_ms = milliseconds(pipeline_end - dispatch_start);
        timing->command_setup_ms = milliseconds(submit_start - pipeline_end);
        timing->submit_wait_ms = milliseconds(submit_end - submit_start);
        timing->cleanup_ms = milliseconds(dispatch_end - submit_end);
        timing->total_ms = milliseconds(dispatch_end - dispatch_start);
    }

    if (submit_result != VK_SUCCESS) {
        return Fail(ErrorDomain::kDevice, "submit compute dispatch", VkResultText(submit_result));
    }
    return {};
}

}  // namespace sirius::backend
