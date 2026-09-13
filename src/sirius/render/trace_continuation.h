#pragma once

// Shared host layout for gr_trace_continuation.slang. These records belong to
// the active dispatch region, independently of radiance tile residency.
#include "sirius/base/error.h"

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <string>

namespace sirius::render {

enum class TraceAction : std::uint32_t { Initialise = 0, Advance = 1, Finalize = 2 };
enum class TracePhase : std::uint32_t { Running = 1, Terminal = 2, Finalized = 3 };
enum class TraceTermination : std::uint32_t {
    Running = 0,
    Capture = 1,
    Disk = 2,
    ObserverEscape = 3,
    WorkLimit = 4,
    MinimumStep = 5,
    UnrepresentedEvent = 6,
    InvalidInput = 7,
    IdentityMismatch = 8,
    Arithmetic = 9,
    OppositeEscape = 10,
};

template <typename Real>
struct alignas(4 * sizeof(Real)) TraceContinuationRecord {
    // x, k, four X columns, four covariant V columns, compensated x/k, and
    // (adaptive step, accepted affine length, last interval, beam seed).
    std::array<std::array<Real, 4>, 13> physical;
    static constexpr std::size_t kPositionColumns = 2;
    static constexpr std::size_t kCovariantColumns = 6;
    static constexpr std::size_t kIntegration = 12;
    std::array<float, 4> source;
    std::array<float, 4> volume;
    // phase, termination, attempts, accepted attempts.
    std::array<std::uint32_t, 4> control;
    // absolute pixel x/y, original sample index, epoch.
    std::array<std::uint32_t, 4> identity;
};

using FloatTraceContinuation = TraceContinuationRecord<float>;
using DoubleTraceContinuation = TraceContinuationRecord<double>;
static_assert(sizeof(FloatTraceContinuation) == 272);
static_assert(sizeof(DoubleTraceContinuation) == 480);
static_assert(offsetof(FloatTraceContinuation, control) == 240);
static_assert(offsetof(DoubleTraceContinuation, control) == 448);
static_assert(offsetof(FloatTraceContinuation, identity) == 256);
static_assert(offsetof(DoubleTraceContinuation, identity) == 464);

inline constexpr std::uint32_t kTraceParameterCount = 70;
inline constexpr std::size_t kTraceActionParameter = 68;
inline constexpr std::size_t kTraceEpochParameter = 69;

[[nodiscard]] constexpr bool IsPhysicalTraceTermination(std::uint32_t value) {
    return value == static_cast<std::uint32_t>(TraceTermination::Capture) ||
           value == static_cast<std::uint32_t>(TraceTermination::Disk) ||
           value == static_cast<std::uint32_t>(TraceTermination::ObserverEscape) ||
           value == static_cast<std::uint32_t>(TraceTermination::OppositeEscape);
}

// Validate the actual readback before the host advances a sample or publishes
// radiance. The previous record is required for Advance and Finalize; Init is
// the only operation that may discard an old generation.
template <typename Real>
[[nodiscard]] base::Expected<bool> ValidateTraceContinuation(
    const TraceContinuationRecord<Real>& current, const TraceContinuationRecord<Real>* previous,
    TraceAction action, const std::array<std::uint32_t, 4>& identity,
    std::uint32_t maximum_attempts) {
    const auto invalid = [](const std::string& why) -> base::Expected<bool> {
        return base::Fail(base::ErrorDomain::kKernel, "validate trace continuation", why);
    };
    if (current.identity != identity) return invalid("pixel, sample or epoch mismatch");
    for (const auto& vector : current.physical) {
        for (const auto value : vector) {
            if (!std::isfinite(value)) return invalid("non-finite physical state");
        }
    }
    for (const auto value : current.source) {
        if (!std::isfinite(value)) return invalid("non-finite source state");
    }
    for (const auto value : current.volume) {
        if (!std::isfinite(value)) return invalid("non-finite volume state");
    }
    const auto [phase, reason, attempts, accepts] = current.control;
    const auto running = static_cast<std::uint32_t>(TracePhase::Running);
    const auto terminal = static_cast<std::uint32_t>(TracePhase::Terminal);
    const auto finalized = static_cast<std::uint32_t>(TracePhase::Finalized);
    if (attempts > maximum_attempts || accepts > attempts) {
        return invalid("invalid attempt counters");
    }
    if (phase != running && phase != terminal && phase != finalized) {
        return invalid("invalid phase");
    }
    if ((phase == running && reason != 0) ||
        (phase != running && !IsPhysicalTraceTermination(reason))) {
        return invalid("trace termination " + std::to_string(reason) + " cannot publish radiance");
    }
    if (action == TraceAction::Initialise) {
        if (phase == finalized || attempts != 0 || accepts != 0) {
            return invalid("initialization advanced or finalized a ray");
        }
    } else {
        if (previous == nullptr || previous->identity != identity) {
            return invalid("missing matching previous state");
        }
        const auto [old_phase, old_reason, old_attempts, old_accepts] = previous->control;
        const bool physical_unchanged =
            std::memcmp(&current, previous, offsetof(TraceContinuationRecord<Real>, control)) == 0;
        if (action == TraceAction::Finalize) {
            if (old_phase != terminal || phase != finalized || reason != old_reason ||
                attempts != old_attempts || accepts != old_accepts || !physical_unchanged) {
                return invalid("finalization changed a trajectory or did not finalize");
            }
        } else if (action == TraceAction::Advance) {
            if (old_phase == terminal) {
                if (current.control != previous->control || !physical_unchanged) {
                    return invalid("terminal state advanced");
                }
            } else if (old_phase != running || phase == finalized || attempts < old_attempts ||
                       attempts - old_attempts > 1 || accepts < old_accepts ||
                       accepts - old_accepts > attempts - old_attempts ||
                       (phase == running && attempts == old_attempts)) {
                return invalid("advance violated the one-attempt transition");
            }
            if (accepts == old_accepts &&
                (std::memcmp(current.physical.data(), previous->physical.data(),
                             TraceContinuationRecord<Real>::kIntegration *
                                 sizeof(current.physical[0])) != 0 ||
                 std::memcmp(&current.physical[TraceContinuationRecord<Real>::kIntegration][1],
                             &previous->physical[TraceContinuationRecord<Real>::kIntegration][1],
                             3 * sizeof(Real)) != 0 ||
                 std::memcmp(current.source.data(), previous->source.data(),
                             sizeof(current.source)) != 0 ||
                 std::memcmp(current.volume.data(), previous->volume.data(),
                             sizeof(current.volume)) != 0)) {
                return invalid("unaccepted attempt changed committed state or source effects");
            }
        } else {
            return invalid("invalid host action");
        }
    }
    return phase != running;
}

}  // namespace sirius::render
