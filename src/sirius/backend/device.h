#pragma once

// The device abstraction of docs/ARCHITECTURE.md section 4: enumerate,
// report capabilities, load kernels, move buffers, dispatch. Everything
// above this seam is backend-agnostic; every vendor API is an adapter
// behind it. The first adapter is Vulkan compute (vulkan/), which alone
// reaches AMD, Intel, and NVIDIA silicon plus Lavapipe software fallback.

#include "sirius/base/error.h"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <span>
#include <string>
#include <vector>

namespace sirius::backend {

enum class DeviceKind {
    kSoftware,       // CPU implementation behind a GPU API (Lavapipe, SwiftShader)
    kIntegratedGpu,  // shares system memory; the 2 GB budget class
    kDiscreteGpu,
    kOther,
};

[[nodiscard]] constexpr const char* ToString(DeviceKind kind) noexcept {
    switch (kind) {
        case DeviceKind::kSoftware:
            return "software";
        case DeviceKind::kIntegratedGpu:
            return "integrated";
        case DeviceKind::kDiscreteGpu:
            return "discrete";
        case DeviceKind::kOther:
            return "other";
    }
    return "unknown";
}

struct DeviceInfo {
    std::string name;
    DeviceKind kind = DeviceKind::kOther;
    // Largest device-local heap; the memory governor derives tile budgets
    // from this, never from full-frame hope.
    std::uint64_t device_local_bytes = 0;
    // Whether fp64 kernels are available (precision ladder rung one).
    bool supports_fp64 = false;
};

// Opaque per-device handles; values are indices into the owning device's
// tables and are meaningless across devices.
struct KernelHandle {
    std::uint32_t value = 0;
};

struct BufferHandle {
    std::uint32_t value = 0;
};

enum class BufferUsage {
    kStorage,  // read-write structured data
    kUniform,  // small read-only parameter blocks
};

// One compute device. Synchronous by design at this seam: a Dispatch
// returns when results are readable. Tile-level parallelism lives above
// (the scheduler overlaps tiles, not intra-tile commands), which keeps
// every adapter simple enough to verify against the CPU reference.
class ComputeDevice {
  public:
    virtual ~ComputeDevice() = default;

    [[nodiscard]] virtual const DeviceInfo& Info() const noexcept = 0;

    // SPIR-V module with a single compute entry point named "main"
    // (slangc renames every entry point to "main" on SPIR-V emission).
    [[nodiscard]] virtual base::Expected<KernelHandle> LoadKernel(
        std::span<const std::uint32_t> spirv) = 0;

    [[nodiscard]] virtual base::Expected<BufferHandle> CreateBuffer(std::uint64_t size_bytes,
                                                                    BufferUsage usage) = 0;

    [[nodiscard]] virtual base::Expected<void> WriteBuffer(BufferHandle buffer,
                                                           std::span<const std::byte> data) = 0;

    [[nodiscard]] virtual base::Expected<void> ReadBuffer(BufferHandle buffer,
                                                          std::span<std::byte> out) = 0;

    // Binds `buffers` to descriptor set 0, bindings 0..N-1 in order, and
    // dispatches the given workgroup counts.
    [[nodiscard]] virtual base::Expected<void> Dispatch(KernelHandle kernel,
                                                        std::span<const BufferHandle> buffers,
                                                        std::uint32_t groups_x,
                                                        std::uint32_t groups_y,
                                                        std::uint32_t groups_z) = 0;
};

// Enumerates Vulkan-visible devices (empty vector when no loader or ICD is
// present; that is a decline, not an error).
[[nodiscard]] base::Expected<std::vector<DeviceInfo>> EnumerateVulkanDevices();

// Opens device `index` as enumerated by EnumerateVulkanDevices.
[[nodiscard]] base::Expected<std::unique_ptr<ComputeDevice>> CreateVulkanDevice(std::size_t index);

}  // namespace sirius::backend
