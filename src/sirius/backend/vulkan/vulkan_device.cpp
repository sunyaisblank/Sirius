#include "sirius/backend/vulkan/vulkan_device.h"

#include "sirius/base/contracts.h"

#include <charconv>
#include <cstdlib>
#include <cstring>
#include <format>
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
    VkPhysicalDeviceDriverProperties driver{
        .sType = VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_DRIVER_PROPERTIES,
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

    constexpr VkMemoryPropertyFlags kRenderMemory =
        VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT;
    std::uint64_t render_memory = 0;
    for (std::uint32_t i = 0; i < memory.memoryTypeCount; ++i) {
        if ((memory.memoryTypes[i].propertyFlags & kRenderMemory) != kRenderMemory) {
            continue;
        }
        const std::uint32_t heap = memory.memoryTypes[i].heapIndex;
        render_memory = std::max(render_memory, memory.memoryHeaps[heap].size);
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
    VkInstance instance = VK_NULL_HANDLE;
    if (const VkResult r = vkCreateInstance(&create_info, nullptr, &instance); r != VK_SUCCESS) {
        return Fail(ErrorDomain::kDevice, "create Vulkan instance", VkResultText(r));
    }
    return instance;
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
    if (devices.empty()) {
        return Fail(ErrorDomain::kDevice, "select Vulkan device", "no devices were enumerated");
    }
    const char* raw = std::getenv("SIRIUS_VULKAN_DEVICE");
    if (raw == nullptr || *raw == '\0') {
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
    if (const VkResult r =
            vkCreateDevice(device->physical_, &device_info, nullptr, &device->device_);
        r != VK_SUCCESS) {
        return Fail(ErrorDomain::kDevice, "create Vulkan logical device", VkResultText(r));
    }
    vkGetDeviceQueue(device->device_, family, 0, &device->queue_);

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
        for (auto& [key, pipeline] : pipelines_) {
            vkDestroyPipeline(device_, pipeline.pipeline, nullptr);
            vkDestroyPipelineLayout(device_, pipeline.layout, nullptr);
            vkDestroyDescriptorSetLayout(device_, pipeline.set_layout, nullptr);
        }
        for (VkShaderModule module : kernels_) {
            vkDestroyShaderModule(device_, module, nullptr);
        }
        for (Buffer& buffer : buffers_) {
            vkDestroyBuffer(device_, buffer.buffer, nullptr);
            vkFreeMemory(device_, buffer.memory, nullptr);
        }
        vkDestroyDescriptorPool(device_, descriptor_pool_, nullptr);
        vkDestroyCommandPool(device_, command_pool_, nullptr);
        vkDestroyDevice(device_, nullptr);
    }
    if (instance_ != VK_NULL_HANDLE) {
        vkDestroyInstance(instance_, nullptr);
    }
}

Expected<KernelHandle> VulkanDevice::LoadKernel(std::span<const std::uint32_t> spirv) {
    SIRIUS_PRE(!spirv.empty());
    const VkShaderModuleCreateInfo create_info{
        .sType = VK_STRUCTURE_TYPE_SHADER_MODULE_CREATE_INFO,
        .codeSize = spirv.size_bytes(),
        .pCode = spirv.data(),
    };
    VkShaderModule module = VK_NULL_HANDLE;
    if (const VkResult r = vkCreateShaderModule(device_, &create_info, nullptr, &module);
        r != VK_SUCCESS) {
        return Fail(ErrorDomain::kKernel, "create shader module", VkResultText(r));
    }
    kernels_.push_back(module);
    return KernelHandle{static_cast<std::uint32_t>(kernels_.size() - 1)};
}

Expected<BufferHandle> VulkanDevice::CreateBuffer(std::uint64_t size_bytes, BufferUsage usage) {
    SIRIUS_PRE(size_bytes > 0);
    const VkBufferUsageFlags usage_flags = usage == BufferUsage::kStorage
                                               ? VK_BUFFER_USAGE_STORAGE_BUFFER_BIT
                                               : VK_BUFFER_USAGE_UNIFORM_BUFFER_BIT;
    const VkBufferCreateInfo buffer_info{
        .sType = VK_STRUCTURE_TYPE_BUFFER_CREATE_INFO,
        .size = size_bytes,
        .usage = usage_flags,
        .sharingMode = VK_SHARING_MODE_EXCLUSIVE,
    };
    Buffer buffer{.size_bytes = size_bytes, .usage = usage};
    if (const VkResult r = vkCreateBuffer(device_, &buffer_info, nullptr, &buffer.buffer);
        r != VK_SUCCESS) {
        return Fail(ErrorDomain::kDevice, "create buffer", VkResultText(r));
    }

    VkMemoryRequirements requirements{};
    vkGetBufferMemoryRequirements(device_, buffer.buffer, &requirements);
    VkPhysicalDeviceMemoryProperties memory_properties{};
    vkGetPhysicalDeviceMemoryProperties(physical_, &memory_properties);

    // Host-visible coherent memory in this increment; device-local staging
    // is governor work, where the budget model decides placement.
    constexpr VkMemoryPropertyFlags kWanted =
        VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT;
    std::uint32_t type_index = memory_properties.memoryTypeCount;
    std::uint64_t selected_heap_size = 0;
    for (std::uint32_t i = 0; i < memory_properties.memoryTypeCount; ++i) {
        const bool allowed = (requirements.memoryTypeBits & (1u << i)) != 0;
        const bool suitable = (memory_properties.memoryTypes[i].propertyFlags & kWanted) == kWanted;
        const std::uint64_t heap_size =
            memory_properties.memoryHeaps[memory_properties.memoryTypes[i].heapIndex].size;
        if (allowed && suitable && heap_size > selected_heap_size) {
            type_index = i;
            selected_heap_size = heap_size;
        }
    }
    if (type_index == memory_properties.memoryTypeCount) {
        vkDestroyBuffer(device_, buffer.buffer, nullptr);
        return Fail(ErrorDomain::kDevice, "allocate buffer memory",
                    "no host-visible coherent memory type");
    }

    const VkMemoryAllocateInfo allocate_info{
        .sType = VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO,
        .allocationSize = requirements.size,
        .memoryTypeIndex = type_index,
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
    return BufferHandle{static_cast<std::uint32_t>(buffers_.size() - 1)};
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
    KernelHandle kernel, std::span<const BufferHandle> buffers) {
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
                .module = kernels_[kernel.value],
                .pName = "main",
            },
        .layout = pipeline.layout,
    };
    if (const VkResult r = vkCreateComputePipelines(device_, VK_NULL_HANDLE, 1, &pipeline_info,
                                                    nullptr, &pipeline.pipeline);
        r != VK_SUCCESS) {
        vkDestroyPipelineLayout(device_, pipeline.layout, nullptr);
        vkDestroyDescriptorSetLayout(device_, pipeline.set_layout, nullptr);
        return Fail(ErrorDomain::kKernel, "create compute pipeline", VkResultText(r));
    }

    auto [inserted, _] = pipelines_.emplace(std::move(key), pipeline);
    return &inserted->second;
}

Expected<void> VulkanDevice::Dispatch(KernelHandle kernel, std::span<const BufferHandle> buffers,
                                      std::uint32_t groups_x, std::uint32_t groups_y,
                                      std::uint32_t groups_z) {
    SIRIUS_PRE(kernel.value < kernels_.size());
    SIRIUS_PRE(groups_x > 0 && groups_y > 0 && groups_z > 0);
    for (const BufferHandle handle : buffers) {
        SIRIUS_PRE(handle.value < buffers_.size());
    }

    auto pipeline = GetOrCreatePipeline(kernel, buffers);
    if (!pipeline) {
        return std::unexpected(pipeline.error());
    }

    const VkDescriptorSetAllocateInfo set_info{
        .sType = VK_STRUCTURE_TYPE_DESCRIPTOR_SET_ALLOCATE_INFO,
        .descriptorPool = descriptor_pool_,
        .descriptorSetCount = 1,
        .pSetLayouts = &(*pipeline)->set_layout,
    };
    VkDescriptorSet set = VK_NULL_HANDLE;
    if (const VkResult r = vkAllocateDescriptorSets(device_, &set_info, &set); r != VK_SUCCESS) {
        return Fail(ErrorDomain::kDevice, "allocate descriptor set", VkResultText(r));
    }

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
            .dstSet = set,
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

    const VkCommandBufferAllocateInfo command_info{
        .sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_ALLOCATE_INFO,
        .commandPool = command_pool_,
        .level = VK_COMMAND_BUFFER_LEVEL_PRIMARY,
        .commandBufferCount = 1,
    };
    VkCommandBuffer command = VK_NULL_HANDLE;
    if (const VkResult r = vkAllocateCommandBuffers(device_, &command_info, &command);
        r != VK_SUCCESS) {
        vkFreeDescriptorSets(device_, descriptor_pool_, 1, &set);
        return Fail(ErrorDomain::kDevice, "allocate command buffer", VkResultText(r));
    }

    const VkCommandBufferBeginInfo begin_info{
        .sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO,
        .flags = VK_COMMAND_BUFFER_USAGE_ONE_TIME_SUBMIT_BIT,
    };
    vkBeginCommandBuffer(command, &begin_info);
    vkCmdBindPipeline(command, VK_PIPELINE_BIND_POINT_COMPUTE, (*pipeline)->pipeline);
    vkCmdBindDescriptorSets(command, VK_PIPELINE_BIND_POINT_COMPUTE, (*pipeline)->layout, 0, 1,
                            &set, 0, nullptr);
    vkCmdDispatch(command, groups_x, groups_y, groups_z);
    vkEndCommandBuffer(command);

    const VkSubmitInfo submit_info{
        .sType = VK_STRUCTURE_TYPE_SUBMIT_INFO,
        .commandBufferCount = 1,
        .pCommandBuffers = &command,
    };
    VkResult submit_result = vkQueueSubmit(queue_, 1, &submit_info, VK_NULL_HANDLE);
    if (submit_result == VK_SUCCESS) {
        submit_result = vkQueueWaitIdle(queue_);
    }

    vkFreeCommandBuffers(device_, command_pool_, 1, &command);
    vkFreeDescriptorSets(device_, descriptor_pool_, 1, &set);

    if (submit_result != VK_SUCCESS) {
        return Fail(ErrorDomain::kDevice, "submit compute dispatch", VkResultText(submit_result));
    }
    return {};
}

}  // namespace sirius::backend
