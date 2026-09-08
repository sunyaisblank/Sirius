#pragma once

// The memory governor owns residency. This controller chooses the area of each
// synchronous compute submission within independent, non-bypassable shape caps.
// Timing feedback can shrink subsequent work; it cannot preempt an individual
// trajectory or guarantee duration on an unmeasured device/scene.

#include "sirius/base/error.h"
#include "sirius/render/pixel_sampling.h"

#include <algorithm>
#include <array>
#include <cstdint>
#include <functional>

namespace sirius::render {

inline constexpr double kDefaultDispatchTargetMs = 250.0;
// The physical Dozen heavy-fp32 128x128 render retained exact output with this
// target: 288 submissions, maximum 777.5 ms, no minimum-work fallback. A lower
// target repeatedly shrank useful bands despite the fixed submission cost.
// Other precision/workload profiles retain the original default above.
inline constexpr double kHeavyFp32DispatchTargetMs = 750.0;

inline constexpr int kInitialBandRows = 1;
// Ordinary fp32/compensated Kerr disk renders also need submission bounds:
// residency-derived growth reached 512x128 and 1948.7 ms on Radeon/Dozen.
// Full 512x512, four-sample probes at 512x16 stayed below 379 ms on both
// rungs; fp32 preserved the original PNG bytes. This is measured work sizing,
// not a duration guarantee for other scenes or devices. Residency stays separate.
inline constexpr int kOrdinaryMaxBandWidth = 512;
inline constexpr int kOrdinaryMaxBandRows = 16;
inline constexpr std::int64_t kOrdinaryMaxPixels = 8192;
inline constexpr int kWatchdogSafeMaxTileEdge = 64;
inline constexpr int kWatchdogSafeMaxBandWidth = 64;
// On the physical Radeon/Dozen fp32 beam/catalogue probe, 64x4 completed and
// preserved output bytes, while a 64x8 sample took 1888 ms. Keep the accepted
// cap at 256 active pixels. This observation is not an all-scene time guarantee.
inline constexpr int kWatchdogSafeMaxBandRows = 4;
inline constexpr std::int64_t kWatchdogSafeMaxPixels = 256;
// Wider fp64 and compensated-fp32 submissions have not earned that evidence.
inline constexpr int kConservativeMaxBandRows = 1;
inline constexpr std::int64_t kConservativeMaxPixels = 64;
// Stop submitting an irreducible band after this observed submit/wait duration.
// This leaves margin below the target route's approximately two-second watchdog;
// the first unexpectedly slow submission still cannot be preempted by the host.
inline constexpr double kDispatchStopMs = 1000.0;
inline constexpr double kBandGrowthCap = 2.0;

// Absolute image rectangle. Subdivision changes only private dispatch packing;
// each child must restart its samples before its pixels can be published.
struct DispatchRegion {
    int x;
    int y;
    int width;
    int height;

    [[nodiscard]] std::int64_t Pixels() const { return static_cast<std::int64_t>(width) * height; }
};

// Preconditions: positive dimensions and more than one pixel. Split rows first,
// then columns, preserving the exact disjoint union even for odd-sized tails.
[[nodiscard]] std::array<DispatchRegion, 2> SplitDispatchRegion(const DispatchRegion& region);

class BandController {
  public:
    BandController(int tile_edge, double target_ms)
        : BandController(tile_edge, target_ms, tile_edge) {}
    BandController(int tile_edge, double target_ms, int max_band_rows)
        : BandController(tile_edge, target_ms, max_band_rows,
                         static_cast<std::int64_t>(std::max(1, tile_edge)) *
                             std::clamp(max_band_rows, 1, std::max(1, tile_edge))) {}
    BandController(int band_width, double target_ms, int max_band_rows, std::int64_t max_pixels)
        : max_pixels_(std::max<std::int64_t>(1, max_pixels)),
          max_rows_(std::max(1, max_band_rows)),
          target_ms_(target_ms),
          pixels_(std::min<std::int64_t>(std::max(1, band_width), max_pixels_)),
          safety_pixel_cap_(max_pixels_) {}

    // Preconditions: positive dimensions and band_width <= the hard area cap.
    // Both area and row bounds apply even when adaptation is disabled. Safety
    // fallback keeps one row and requires subdivision within the sticky safety area cap.
    [[nodiscard]] int NextRows(int remaining_rows, int band_width) const;

    // False means invalid timing or one pixel exceeded kDispatchStopMs: the caller
    // must decline further submissions. Zero duration activates one-pixel subdivision;
    // a reducible overshoot halves the actual area into a sticky, non-growing safety
    // cap. Ordinary under-target observations
    // double actual area; over-target ones halve it. This reaches useful whole
    // bands despite a fixed submission-time floor.
    // Tail feedback uses only the area actually dispatched.
    bool Record(std::int64_t dispatched_pixels, double measured_ms);

    [[nodiscard]] bool Enabled() const { return target_ms_ > 0.0; }
    [[nodiscard]] double TargetMs() const { return target_ms_; }
    [[nodiscard]] bool SafetyFallback() const { return safety_fallback_; }
    [[nodiscard]] std::int64_t SafetyPixelCap() const { return safety_pixel_cap_; }

  private:
    std::int64_t max_pixels_;
    int max_rows_;
    double target_ms_;
    std::int64_t pixels_;
    std::int64_t safety_pixel_cap_;
    bool safety_fallback_ = false;
};

// Execute one private logical band within the caller's hard shape/area caps.
// Subdivision only shrinks this region. Submission observes the original camera
// sample order; completion is called only after every sample of a final region.
// A reducible safety fallback discards that region's partial accumulation and
// replays each disjoint child from sample zero. No device timing is synthesized.
[[nodiscard]] base::Expected<void> ExecuteDispatchRegions(
    const DispatchRegion& region, int samples_per_pixel, BandController& bands,
    const std::function<base::Expected<double>(const DispatchRegion&, const CameraSample&, int)>&
        submit,
    const std::function<base::Expected<void>(const DispatchRegion&)>& completed,
    const std::function<bool()>& should_cancel = {});

// Resolves the per-dispatch wall-time target: the supplied profile default unless
// SIRIUS_DISPATCH_TARGET_MS overrides it. Zero disables adaptation but retains hard caps; a
// negative, non-finite, overflowing, or unparseable value is a loud error, never a silent default
// (the same contract SIRIUS_PRECISION carries).
[[nodiscard]] base::Expected<double> ResolveDispatchTargetMs(
    double default_target_ms = kDefaultDispatchTargetMs);

}  // namespace sirius::render
