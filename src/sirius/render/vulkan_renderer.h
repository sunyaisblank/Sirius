#pragma once

// Vulkan render path: dispatches the Slang `trace` compute kernel through a
// backend::ComputeDevice, one governed tile at a time, and writes linear
// radiance into the session's display buffer. Tonemapping and grading stay
// host-side (the display pipeline owns them); this path returns linear radiance
// per pixel exactly as the CPU tracer does, so the two backends feed the same
// output writers.
//
// The renderer is self-contained: it opens the device, applies the memory
// governor and the precision ladder, uploads the starfield when it fits the
// budget, and dispatches the kernel. It declines loudly (a base::Error) for any
// metric or scene semantics outside the Vulkan render path and when a requested
// precision rung is unsupported, never substituting a different render.

#include "sirius/base/error.h"
#include "sirius/render/memory_governor.h"

#include <cstdint>
#include <functional>
#include <string>

namespace sirius::render {

struct SessionConfig;
class DisplayBuffer;

// The precision-ladder rung the render ran on (recorded in metadata/logs).
enum class PrecisionRung {
    Fp32,      // plain single precision, the default (trace.spv)
    Fp32Comp,  // fp32 with Kahan-compensated state accumulation
               // (trace_fp32comp.spv); SIRIUS_PRECISION=fp32-comp, any device
    Fp64,      // double-precision trajectory core (trace_fp64.spv); selected by
               // SIRIUS_PRECISION=fp64 on devices reporting shaderFloat64
};

// Maximum memory-tile edge and row-band height for one Vulkan workload. fp64
// trajectories and the wider ray-bundle/point-catalogue kernel use the strict
// physical-watchdog footprint; ordinary fp32 paths retain the memory governor's
// normal tile cap and adaptive row growth.
struct VulkanDispatchLimits {
    int tile_edge_cap = kMaxTileEdge;
    int max_band_rows = kMaxTileEdge;
};

[[nodiscard]] VulkanDispatchLimits ResolveVulkanDispatchLimits(PrecisionRung precision,
                                                               bool ray_bundles,
                                                               bool point_starfield);

// What the Vulkan render produced, for logging and the parity/governor tests.
struct VulkanRenderStats {
    std::string device_name;
    std::string metric_name;
    TilePlan tile_plan;
    PrecisionRung precision = PrecisionRung::Fp32;
    bool starfield_uploaded = false;
    bool point_catalogue_uploaded = false;
    int tiles_rendered = 0;
    int band_dispatches = 0;  // compute submissions; >= tiles_rendered under banding
    double seconds = 0.0;
};

// Checks the scene features the current one-sample Vulkan kernel represents.
// Auto-selection and the dispatch boundary both use this contract.
[[nodiscard]] base::Expected<void> ValidateVulkanRenderConfig(const SessionConfig& config);

// Renders `config`'s scene on the first Vulkan device into `display` (which the
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
