#pragma once

// Vulkan compute adapter behind the ComputeDevice seam. Deliberately
// synchronous and host-visible-only in this first increment: correctness
// and Lavapipe testability first; device-local staging and submission
// overlap arrive with the memory governor (specification programme 4),
// where profiling can justify them.

#include "sirius/backend/device.h"

#include <vulkan/vulkan.h>

#include <cstdint>
#include <limits>
#include <map>
#include <optional>
#include <vector>

namespace sirius::backend {

// Check instruction framing and device-dependent Float64/FMA requirements before
// asking the driver to create a module. Full SPIR-V validity remains spirv-val's gate.
[[nodiscard]] base::Expected<void> ValidateVulkanKernelPrecision(
    std::span<const std::uint32_t> spirv, bool supports_fp64, bool fma_fp32_enabled = false);

namespace detail {
// Keep the largest compatible coherent host heap and its existing tie order;
// prefer the same memory properties with host caching on that selected heap.
[[nodiscard]] std::optional<std::uint32_t> VulkanHostMemoryType(
    const VkPhysicalDeviceMemoryProperties& properties, std::uint32_t memory_type_bits);
// Bound the one process-lived serialized blob, not driver-internal cache memory.
inline constexpr std::size_t kVulkanPipelineCacheBlobLimit = 32 * 1024 * 1024;
[[nodiscard]] bool VulkanPipelineCacheDataCompatible(std::span<const std::byte> data,
                                                     const VkPhysicalDeviceProperties& properties);
// The enclosing actual host completion interval must exclude a second wrap.
// Upper undefined counter bits are ignored; availability is checked by the caller.
[[nodiscard]] std::optional<double> VulkanTimestampSpanMs(std::uint64_t begin, std::uint64_t end,
                                                          std::uint32_t valid_bits,
                                                          double period_ns,
                                                          double host_completion_ms);
// A submit reporting device loss is not an idle-completion guarantee. Waited
// applies only when submitted was VK_SUCCESS. Allocation errors leave submitted
// state unchanged; an unknown error cannot prove that resources are unused.
[[nodiscard]] bool VulkanSubmissionNeedsCompletion(VkResult submitted, VkResult waited);
}  // namespace detail

struct VulkanPipelineCacheStats {
    bool enabled = false;
    VkResult creation_result = VK_SUCCESS;
    std::size_t imported_bytes = 0;
    std::size_t available_bytes = 0;
    std::size_t exported_bytes = 0;
    bool export_discarded = false;
    bool import_discarded = false;
};

struct VulkanTimestampProperties {
    std::uint32_t valid_bits = 0;  // Selected compute queue; zero means unsupported.
    double period_ns = 0;
};

// One undivided device span, including original barriers and scheduling between
// the two markers. This is not isolated shader time or calibrated host time.
// Query failures/unavailable/ambiguous counters have no device_span_ms value.
struct VulkanDispatchTimestamp {
    VulkanTimestampProperties properties;
    VkResult query_result = VK_NOT_READY;
    std::array<std::uint64_t, 2> ticks{}, availability{};
    std::optional<double> device_span_ms;
    double host_submit_ms = 0, host_wait_ms = 0;
};

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

    [[nodiscard]] VulkanPipelineCacheStats PipelineCacheStatistics() const noexcept {
        return pipeline_cache_stats_;
    }
    // Caller serializes this with other device operations, as with Dispatch.
    // Publish a complete bounded snapshot for later fresh devices; no disk I/O.
    [[nodiscard]] base::Expected<VulkanPipelineCacheStats> SnapshotPipelineCache();

    [[nodiscard]] VulkanTimestampProperties TimestampProperties() const noexcept {
        return timestamp_properties_;
    }
    // Optional diagnostics, disabled by default. Serialize with all device work.
    // Reuses exactly two queries; no buffer allocation or extra submission/wait.
    // Setup refusal leaves ordinary dispatch usable. Query failure never changes
    // a successfully completed numerical dispatch or its host governor timing.
    [[nodiscard]] base::Expected<void> SetDispatchTimestampsEnabled(bool enabled);
    [[nodiscard]] const std::optional<VulkanDispatchTimestamp>& LastDispatchTimestamp()
        const noexcept {
        return last_dispatch_timestamp_;
    }

    [[nodiscard]] base::Expected<KernelHandle> LoadKernel(
        std::span<const std::uint32_t> spirv) override;

    [[nodiscard]] base::Expected<BufferHandle> CreateBuffer(std::uint64_t size_bytes,
                                                            BufferUsage usage) override;

    [[nodiscard]] base::Expected<std::uint64_t> RequiredBufferAllocationBytes(
        std::uint64_t size_bytes, BufferUsage usage) override;

    [[nodiscard]] base::Expected<void> SetBufferAllocationLimit(std::uint64_t bytes) override;
    [[nodiscard]] std::uint64_t BufferAllocationBytes() const noexcept override {
        return buffer_allocation_bytes_;
    }

    [[nodiscard]] base::Expected<void> WriteBuffer(BufferHandle buffer,
                                                   std::span<const std::byte> data) override;

    [[nodiscard]] base::Expected<void> ReadBuffer(BufferHandle buffer,
                                                  std::span<std::byte> out) override;

    [[nodiscard]] base::Expected<void> Dispatch(KernelHandle kernel,
                                                std::span<const BufferHandle> buffers,
                                                std::uint32_t groups_x, std::uint32_t groups_y,
                                                std::uint32_t groups_z,
                                                DispatchTiming* timing = nullptr) override;

    [[nodiscard]] bool SupportsIndependentPair() const noexcept override { return true; }
    [[nodiscard]] base::Expected<void> DispatchIndependentPair(
        const std::array<ComputeDispatch, 2>& commands,
        IndependentPairTiming* timing = nullptr) override;

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
        KernelHandle kernel, std::span<const BufferHandle> buffers, bool* created);
    [[nodiscard]] base::Expected<void> DispatchCommands(std::span<const ComputeDispatch> commands,
                                                        DispatchTiming* timing,
                                                        std::uint32_t* pipeline_creations);
    void InitialisePipelineCache();

    // At most one accepted submission can lack a confirmed completion. Further
    // work/query reuse refuses; destruction waits before releasing its resources.
    struct PendingDispatch {
        VkCommandBuffer command = VK_NULL_HANDLE;
        std::array<VkDescriptorSet, 2> sets{};
        std::uint32_t set_count = 0;
    };
    std::optional<PendingDispatch> pending_dispatch_;

    VkInstance instance_ = VK_NULL_HANDLE;
    VkPhysicalDevice physical_ = VK_NULL_HANDLE;
    VkDevice device_ = VK_NULL_HANDLE;
    VkQueue queue_ = VK_NULL_HANDLE;
    std::uint32_t queue_family_ = 0;
    VkCommandPool command_pool_ = VK_NULL_HANDLE;
    VkDescriptorPool descriptor_pool_ = VK_NULL_HANDLE;
    VkPipelineCache pipeline_cache_ = VK_NULL_HANDLE;
    VkQueryPool timestamp_pool_ = VK_NULL_HANDLE;
    VulkanTimestampProperties timestamp_properties_;
    std::optional<VulkanDispatchTimestamp> last_dispatch_timestamp_;
    VkPhysicalDeviceProperties pipeline_cache_properties_{};
    VulkanPipelineCacheStats pipeline_cache_stats_;
    bool pipeline_cache_modified_ = false;
    DeviceInfo info_;

    struct Kernel {
        std::vector<std::uint32_t> words;
        VkShaderModule module = VK_NULL_HANDLE;
    };
    std::vector<Kernel> kernels_;
    std::vector<Buffer> buffers_;
    std::uint64_t buffer_allocation_bytes_ = 0;
    std::uint64_t buffer_allocation_limit_ = std::numeric_limits<std::uint64_t>::max();
    std::map<PipelineKey, Pipeline> pipelines_;
};

}  // namespace sirius::backend
