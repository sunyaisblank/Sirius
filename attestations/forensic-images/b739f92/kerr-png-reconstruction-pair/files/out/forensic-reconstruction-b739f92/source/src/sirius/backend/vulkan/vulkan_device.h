#pragma once

// Vulkan compute adapter behind the ComputeDevice seam. Deliberately
// synchronous and host-visible-only in this first increment: correctness
// and Lavapipe testability first; device-local staging and submission
// overlap arrive with the memory governor (specification programme 4),
// where profiling can justify them.

#include <cstdint>
#include <map>
#include <vector>

#include <vulkan/vulkan.h>

#include "sirius/backend/device.h"

namespace sirius::backend {

class VulkanDevice final : public ComputeDevice {
  public:
    // Use CreateVulkanDevice(); this is public only for std::make_unique.
    VulkanDevice() = default;
    ~VulkanDevice() override;

    VulkanDevice(const VulkanDevice&) = delete;
    VulkanDevice& operator=(const VulkanDevice&) = delete;
    VulkanDevice(VulkanDevice&&) = delete;
    VulkanDevice& operator=(VulkanDevice&&) = delete;

    [[nodiscard]] const DeviceInfo& Info() const noexcept override { return info_; }

    [[nodiscard]] base::Expected<KernelHandle> LoadKernel(
        std::span<const std::uint32_t> spirv) override;

    [[nodiscard]] base::Expected<BufferHandle> CreateBuffer(std::uint64_t size_bytes,
                                                            BufferUsage usage) override;

    [[nodiscard]] base::Expected<void> WriteBuffer(BufferHandle buffer,
                                                   std::span<const std::byte> data) override;

    [[nodiscard]] base::Expected<void> ReadBuffer(BufferHandle buffer,
                                                  std::span<std::byte> out) override;

    [[nodiscard]] base::Expected<void> Dispatch(KernelHandle kernel,
                                                std::span<const BufferHandle> buffers,
                                                std::uint32_t groups_x, std::uint32_t groups_y,
                                                std::uint32_t groups_z) override;

  private:
    friend base::Expected<std::unique_ptr<ComputeDevice>> CreateVulkanDevice(std::size_t index);

    struct Buffer {
        VkBuffer buffer = VK_NULL_HANDLE;
        VkDeviceMemory memory = VK_NULL_HANDLE;
        std::uint64_t size_bytes = 0;
        BufferUsage usage = BufferUsage::kStorage;
    };

    // Pipelines are cached per kernel and binding signature; the signature
    // is the ordered list of buffer usages, which fixes the descriptor
    // layout.
    struct PipelineKey {
        std::uint32_t kernel = 0;
        std::vector<BufferUsage> bindings;
        auto operator<=>(const PipelineKey&) const = default;
    };

    struct Pipeline {
        VkDescriptorSetLayout set_layout = VK_NULL_HANDLE;
        VkPipelineLayout layout = VK_NULL_HANDLE;
        VkPipeline pipeline = VK_NULL_HANDLE;
    };

    [[nodiscard]] base::Expected<Pipeline*> GetOrCreatePipeline(
        KernelHandle kernel, std::span<const BufferHandle> buffers);

    VkInstance instance_ = VK_NULL_HANDLE;
    VkPhysicalDevice physical_ = VK_NULL_HANDLE;
    VkDevice device_ = VK_NULL_HANDLE;
    VkQueue queue_ = VK_NULL_HANDLE;
    std::uint32_t queue_family_ = 0;
    VkCommandPool command_pool_ = VK_NULL_HANDLE;
    VkDescriptorPool descriptor_pool_ = VK_NULL_HANDLE;
    DeviceInfo info_;

    std::vector<VkShaderModule> kernels_;
    std::vector<Buffer> buffers_;
    std::map<PipelineKey, Pipeline> pipelines_;
};

}  // namespace sirius::backend
