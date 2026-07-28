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

#include "sirius/base/error.h"

#include <algorithm>

namespace sirius::render {

// Wall-time target per dispatch, milliseconds. A quarter second keeps an
// eight-fold estimation error inside a 2 s watchdog while wasting little
// submission overhead on the rungs that could sustain longer dispatches.
inline constexpr double kDefaultDispatchTargetMs = 250.0;

// Rows in the first band, before any measurement exists. Small enough that
// even the fp64 rung on integrated silicon finishes well inside the target's
// order of magnitude; the controller grows from here within a few bands.
inline constexpr int kInitialBandRows = 8;

// Per-step growth bound. Growth is damped so one spuriously fast measurement
// (a band of sky, a warm cache) cannot balloon the next dispatch past the
// watchdog; shrinking is deliberately unbounded because overshooting the
// target risks device removal while undershooting only costs overhead.
inline constexpr double kBandGrowthCap = 2.0;

// Adaptive band-height controller. One instance persists across a render so
// the learned throughput carries from tile to tile. Pure arithmetic; no device
// handles cross this boundary, so it is unit-testable without Vulkan.
//
// Invariants: NextRows always returns a value in [1, remaining_rows]; with the
// governor disabled (target_ms <= 0) it returns remaining_rows unchanged, which
// restores the historical one-dispatch-per-tile behaviour.
class BandController {
  public:
    // max_rows caps a band at the governed tile edge; target_ms <= 0 disables
    // banding (the escape hatch SIRIUS_DISPATCH_TARGET_MS=0 documents).
    BandController(int max_rows, double target_ms)
        : max_rows_(std::max(1, max_rows)),
          target_ms_(target_ms),
          rows_(std::clamp(kInitialBandRows, 1, std::max(1, max_rows))) {}

    // Rows the next dispatch should cover, given the rows still unrendered in
    // the current tile. Postcondition: 1 <= result <= remaining_rows, unless
    // banding is disabled, in which case result == remaining_rows.
    [[nodiscard]] int NextRows(int remaining_rows) const;

    // Feeds one dispatch's measured wall time back into the controller. Growth
    // is capped at kBandGrowthCap per step; shrink is proportional and
    // uncapped (floor one row). Disabled controllers ignore measurements.
    void Record(double measured_ms);

    [[nodiscard]] bool Enabled() const { return target_ms_ > 0.0; }
    [[nodiscard]] double TargetMs() const { return target_ms_; }

  private:
    int max_rows_;
    double target_ms_;
    int rows_;
};

// Resolves the per-dispatch wall-time target: kDefaultDispatchTargetMs unless
// SIRIUS_DISPATCH_TARGET_MS overrides it. Zero disables banding; a negative or
// unparseable value is a loud error, never a silent default (the same contract
// SIRIUS_PRECISION carries).
[[nodiscard]] base::Expected<double> ResolveDispatchTargetMs();

}  // namespace sirius::render
