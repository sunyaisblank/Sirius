#pragma once

// VRAM budget governor (docs/ARCHITECTURE.md section 6, specification programme
// 4). The 780M-class target has a 2 GB budget that a naive full-frame HDR
// pipeline exhausts, so device residency is bounded by construction: the
// governor derives the largest square tile whose device-resident working set
// fits a fixed fraction of the reported budget, and the Vulkan render path
// dispatches at that tile size. Full-frame buffers stay host-side; the only
// per-tile device buffer the trace kernel allocates is the RGBA32F radiance
// readback (ray state lives in registers), so the per-pixel working set is
// 16 bytes and any persistent device buffers (params, an uploaded starfield)
// are passed as a fixed overhead. A budget too small to hold the fixed overhead
// plus a minimal tile declines loudly rather than over-committing the device.

#include "sirius/base/error.h"

#include <cstdint>

namespace sirius::render {

// Fraction of the reported budget the governor is allowed to place resident.
// Half leaves headroom for descriptor pools, command buffers, staging the
// driver performs implicitly, and the fragmentation a single heap suffers under
// repeated allocation; it is deliberately conservative for the 2 GB target.
inline constexpr double kResidencyFraction = 0.5;

// Per-tile device working set, bytes per pixel: one RGBA32F radiance readback
// texel. The trace kernel keeps geodesic state (position, momentum, deviation
// vectors, accumulators) in registers, so it is not device-resident and does
// not enter this figure; when the ray-bundle path (programme 5) makes that
// state a buffer, the constant grows and every consumer follows.
inline constexpr std::uint64_t kTileWorkingSetBytesPerPixel = 16;

// Smallest tile edge worth dispatching. Below the 8x8 workgroup a tile wastes
// the launch; a budget that cannot fit an 8x8 tile plus the fixed overhead is a
// decline, not a 1x1 crawl.
inline constexpr int kMinTileEdge = 8;

// Sanity cap on the tile edge independent of budget: beyond this a single tile
// is large enough that per-tile scheduling buys nothing, and the readback
// staging grows without bound.
inline constexpr int kMaxTileEdge = 4096;

// The tile size the scheduler receives on the Vulkan path, with the budget
// arithmetic that produced it (recorded for the render metadata and tests).
struct TilePlan {
    int tile_edge = 0;                       // Square tile edge in pixels.
    std::uint64_t budget_bytes = 0;          // Reported (or overridden) budget.
    std::uint64_t usable_bytes = 0;          // kResidencyFraction * budget.
    std::uint64_t fixed_overhead_bytes = 0;  // Persistent device buffers.
    std::uint64_t tile_working_set_bytes = 0;  // tile_edge^2 * bytes-per-pixel.
};

// Derives the largest tile that fits kResidencyFraction of `budget_bytes` after
// reserving `fixed_overhead_bytes` of persistent device residency, capped at the
// image extent and kMaxTileEdge. Preconditions: image dimensions positive.
// Declines (kDevice error) when the budget cannot seat the fixed overhead plus a
// kMinTileEdge tile.
[[nodiscard]] base::Expected<TilePlan> DeriveTilePlan(std::uint64_t budget_bytes, int image_width,
                                                      int image_height,
                                                      std::uint64_t fixed_overhead_bytes);

// Resolves the budget to plan against: the SIRIUS_MEMORY_BUDGET_MB override when
// set (for constrained-budget testing), else the device-local heap size. A zero
// device budget with no override is itself a decline signal the caller reports.
[[nodiscard]] std::uint64_t ResolveBudgetBytes(std::uint64_t device_local_bytes);

}  // namespace sirius::render
