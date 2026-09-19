#pragma once

// Vulkan render path: dispatches retained stages or the legacy trace shader
// through backend::ComputeDevice within independent work/residency bounds and writes linear
// radiance into the session's display buffer. Tonemapping and grading stay
// host-side (the display pipeline owns them); this path returns linear radiance
// per pixel exactly as the CPU tracer does, so the two backends feed the same
// output writers.
//
// The renderer opens the selected device and applies the memory governor and
// precision ladder. Retained transport shares the host source owner; legacy
// shaders upload their source resources. It declines (a base::Error) for any
// metric or scene semantics outside the Vulkan render path and when a requested
// precision rung is unsupported, never substituting a different render.

#include "sirius/base/error.h"
#include "sirius/render/dispatch_governor.h"
#include "sirius/render/memory_governor.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <string>

namespace sirius::render {

struct SessionConfig;
class DisplayBuffer;

// The precision-ladder rung the render ran on (recorded in metadata/logs).
enum class PrecisionRung {
    Fp32,      // Retained binary32 Kerr-family transport; scalar legacy metrics.
    Fp32Comp,  // Same retained Kerr-family path; compensated legacy metrics.
    Fp64,      // Retained Kerr-family transport with exact binary64 products;
               // binary64 legacy metrics. Requires shaderFloat64.
};

// Independent residency and submission caps. Expensive workloads admit at most
// 256 active fp32 beam/catalogue trajectories in 64x4 bands; fp64 and heavy
// compensated fp32 retain 64x1 pending wider physical evidence. Ordinary fp32
// and compensated fp32 use 512x16 independently of their residency tile.
// These are hard work bounds,
// not duration guarantees; feedback cannot preempt an individual trajectory.
struct VulkanDispatchLimits {
    int tile_edge_cap = kMaxTileEdge;
    int max_band_width = kOrdinaryMaxBandWidth;
    int max_band_rows = kOrdinaryMaxBandRows;
    std::int64_t max_pixels = kOrdinaryMaxPixels;
    double default_target_ms = kDefaultDispatchTargetMs;
};

[[nodiscard]] VulkanDispatchLimits ResolveVulkanDispatchLimits(PrecisionRung precision,
                                                               bool ray_bundles,
                                                               bool point_starfield);

// What the Vulkan render produced, for logging and the parity/governor tests.
struct VulkanRenderStats {
    std::string device_name;
    std::size_t device_index = 0;
    std::string metric_name;
    TilePlan tile_plan;
    std::uint64_t explicit_buffer_allocation_bytes = 0;
    std::uint64_t continuation_capacity = 0;
    // Actual Init/Advance/Finalize submissions; separate from zero-ray driver initialization.
    std::array<std::int64_t, 3> continuation_dispatches{};
    std::array<double, 3> maximum_continuation_ms{};
    bool retained_intervals = false;
    std::uint64_t camera_batches = 0;
    std::uint64_t accepted_intervals = 0;
    // Film camera, joint RK, projection, dense sampling, initialization, smooth ray camera.
    std::array<std::int64_t, 6> retained_stage_dispatches{};
    PrecisionRung precision = PrecisionRung::Fp32;
    bool starfield_uploaded = false;
    bool point_catalogue_uploaded = false;
    int tiles_rendered = 0;
    int work_tile_edge = 0;            // Host publication work, independent of device residency.
    std::int64_t band_dispatches = 0;  // governed ray submissions, excluding initialization
    double dispatch_seconds = 0.0;
    double maximum_dispatch_ms = 0.0;
    std::int64_t maximum_dispatch_rays = 0;
    std::int64_t dispatch_target_overshoots = 0;
    std::int64_t dispatch_subdivisions = 0;
    int dispatch_fallbacks = 0;
    int initialization_dispatches = 0;  // software only, zero active rays
    double initialization_seconds = 0.0;
    double initialization_submit_wait_ms = 0.0;
    double seconds = 0.0;  // complete render wall time, including initialization
};

// Checks the scene features represented by the selected Vulkan route.
// Auto-selection and the dispatch boundary both use this contract.
[[nodiscard]] base::Expected<void> ValidateVulkanRenderConfig(const SessionConfig& config);

// Renders `config`'s scene on the selected Vulkan device into `display` (which the
// caller has already sized to config.width x config.height). `on_tile` reports
// progress as (tiles_done, tiles_total). Preconditions: the display buffer is
// initialised to the config resolution. Postcondition on success: every pixel of
// `display` holds finite linear radiance; on failure nothing is partially
// committed to the caller beyond the error return.
[[nodiscard]] base::Expected<VulkanRenderStats> RenderVulkanToDisplay(
    const SessionConfig& config, DisplayBuffer& display,
    const std::function<void(int tiles_done, int tiles_total)>& on_tile = {},
    const std::function<bool()>& should_cancel = {});

}  // namespace sirius::render
