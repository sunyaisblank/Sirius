#pragma once

// Dispatch-duration governor (companion to the memory governor). The memory
// governor bounds what a tile makes device-resident; this governor bounds how
// long any single compute dispatch keeps the device busy. The two are separate
// authorities because the constraints are physically different: residency is a
// byte budget the device reports, duration is an operating-system watchdog
// (Windows TDR, ~2 s) that removes the device when one dispatch overruns it.
// Software Vulkan has no watchdog, which is why only physical-GPU validation
// exposed the constraint (E3 runbook, 2026-07-28: a 2864x2864 governed tile on
// the Radeon 780M died with D3D12 device removal).
//
// The mechanism is row banding: a tile is submitted as a sequence of shorter
// tiles ("bands") of the same width. The trace kernel derives every pixel from
// absolute image coordinates, so banding is value-transparent by construction;
// it changes submission granularity and nothing else. Band height adapts to a
// wall-time target per dispatch, because a static height cannot serve both the
// fp32 and fp64 rungs (an order of magnitude apart on consumer silicon) nor
// devices of unknown throughput.
//
// The controller's state is a band AREA in pixels, not a row count, and every
// measurement is attributed to the pixels the dispatch actually covered. Both
// choices are load-bearing (adversarial review of 2026-07-28, two confirmed
// findings): a row-count controller misreads a tile's truncated tail band as
// a fast full-height band and ratchets up by the growth cap at every tile
// boundary, and it carries a height learned on a narrow last-column tile into
// the next full-width tile, where the same rows cost proportionally more —
// the whole-tile dispatch shape the governor exists to prevent. Area
// accounting closes both routes. The residual risk is a genuine cost cliff
// between adjacent regions of the image: one dispatch after the cliff may
// overshoot before feedback arrives, bounded by kBandGrowthCap times the
// target per step; a lower SIRIUS_DISPATCH_TARGET_MS buys margin at the cost
// of submission overhead.

#include "sirius/base/error.h"

#include <algorithm>
#include <cstdint>

namespace sirius::render {

// Wall-time target per dispatch, milliseconds. A quarter second keeps an
// eight-fold estimation error inside a 2 s watchdog while wasting little
// submission overhead on the rungs that could sustain longer dispatches.
inline constexpr double kDefaultDispatchTargetMs = 250.0;

// Rows in the first band, at governed-tile width, before any measurement
// exists. One row is the smallest dispatch this banding scheme can represent
// and is therefore the only safe cold start for an unknown workload. The
// physical fp64 path demonstrated that an eight-row estimate which was safe
// before a kernel-cost increase can cross the operating-system watchdog before
// the controller receives its first measurement. Capped growth recovers useful
// throughput within a few bands.
inline constexpr int kInitialBandRows = 1;

// The row-band governor cannot submit less than one full tile row. Keep that
// irreducible fp64 dispatch at or below the 64-pixel width already exercised by
// the physical fp64 runtime gate; wider tiles can still exceed the WSL2/D3D12
// watchdog even at one row after the widened trajectory kernel's cost grows.
inline constexpr int kFp64MaxTileEdge = 64;

// Per-step growth bound. Growth is damped so one spuriously fast measurement
// (a band of sky, a warm cache) cannot balloon the next dispatch past the
// watchdog; shrinking is deliberately unbounded because overshooting the
// target risks device removal while undershooting only costs overhead.
inline constexpr double kBandGrowthCap = 2.0;

// Adaptive band-area controller. One instance covers the bands within a single
// tile; the renderer starts a fresh controller at every tile so a band learned
// on a cheap sky region cannot cross a spatial cost cliff at the next tile.
// NextRows converts the learned area to rows at the width of the band being
// planned, keeping partial-width tiles safe. Pure arithmetic; no device handles
// cross this boundary, so it is unit-testable without Vulkan.
//
// Invariants: NextRows always returns a value in [1, remaining_rows]; with the
// governor disabled (target_ms <= 0) it returns remaining_rows unchanged, which
// restores the historical one-dispatch-per-tile behaviour.
class BandController {
  public:
    // tile_edge bounds the learned area at one governed tile (the radiance
    // buffer's capacity); target_ms <= 0 disables banding (the escape hatch
    // SIRIUS_DISPATCH_TARGET_MS=0 documents).
    BandController(int tile_edge, double target_ms)
        : max_pixels_(static_cast<std::int64_t>(std::max(1, tile_edge)) * std::max(1, tile_edge)),
          target_ms_(target_ms),
          pixels_(std::min<std::int64_t>(
              static_cast<std::int64_t>(kInitialBandRows) * std::max(1, tile_edge), max_pixels_)) {}

    // Rows the next dispatch should cover for a band `band_width` pixels wide,
    // given the rows still unrendered in the current tile. Preconditions:
    // remaining_rows > 0, band_width > 0. Postcondition: 1 <= result <=
    // remaining_rows, unless banding is disabled, in which case result ==
    // remaining_rows.
    [[nodiscard]] int NextRows(int remaining_rows, int band_width) const;

    // Feeds one dispatch's outcome back: the pixels it ACTUALLY covered (which
    // a tail band makes smaller than the controller asked for) and its wall
    // time. Growth is capped at kBandGrowthCap per step; shrink is
    // proportional and uncapped (floor one pixel). Disabled controllers
    // ignore measurements.
    void Record(std::int64_t dispatched_pixels, double measured_ms);

    [[nodiscard]] bool Enabled() const { return target_ms_ > 0.0; }
    [[nodiscard]] double TargetMs() const { return target_ms_; }

  private:
    std::int64_t max_pixels_;
    double target_ms_;
    std::int64_t pixels_;  // learned band area
};

// Resolves the per-dispatch wall-time target: kDefaultDispatchTargetMs unless
// SIRIUS_DISPATCH_TARGET_MS overrides it. Zero disables banding; a negative,
// non-finite, overflowing, or unparseable value is a loud error, never a
// silent default (the same contract SIRIUS_PRECISION carries).
[[nodiscard]] base::Expected<double> ResolveDispatchTargetMs();

}  // namespace sirius::render
