#pragma once

#include "sirius/core/point_source_response.h"
#include "sirius/core/starfield.h"

#include <array>
#include <cstddef>
#include <expected>
#include <functional>
#include <span>
#include <vector>

namespace sirius::render {

inline constexpr int kPointDetectorBlockEdge = 32;
inline constexpr std::size_t kPointDetectorBatchCapacity =
    kPointDetectorBlockEdge * kPointDetectorBlockEdge;

// A scalar sampler uses one original sample's z, q = q0 + L z, |z| <= 4.
// A batch sampler uses a common film q at the same pupil. Subdivision never
// changes the original L or parent Gaussian.
using DetectorCoordinate = std::array<double, 2>;

struct PointDetectorProbe {
    bool visible = false;
    std::array<double, 3> direction{};
    core::AngularMatrix2 source_derivative{};  // Source radians per sampler coordinate.
    double camera_over_source_frequency = 1;
    double transmission = 1;
    std::size_t inner_attempts = 0;
    std::size_t tail_attempts = 0;
};

enum class PointDetectorFailure {
    InvalidInput,
    ProjectionUnavailable,
    TraceFailed,
    Cancelled,
    WorkLimit,
    Unresolved,
    Arithmetic
};

struct PointDetectorPolicy {
    double absolute_rgb_error = 1e-7;
    double relative_rgb_error = 1e-3;
    double geometry_error = 0.0025;
    double maximum_linearization_residual = 0.25;
    double root_error = 2e-6;
    unsigned minimum_depth = 1;
    unsigned maximum_depth = 16;
    std::size_t maximum_probes = 16384;
    std::size_t maximum_cells = 4096;
    std::size_t maximum_candidate_visits = 131072;
    unsigned maximum_newton_steps = 16;
};

struct PointDetectorStatistics {
    std::size_t probes = 0;
    std::size_t probe_requests = 0;
    std::size_t probe_cache_comparisons = 0;
    std::size_t cells = 0;
    std::size_t candidate_visits = 0;
    std::size_t newton_steps = 0;
    std::size_t roots = 0;
    std::size_t reserved_cache_bytes = 0;
    std::size_t inner_attempts = 0;
    std::size_t tail_attempts = 0;
    double maximum_root_coordinate_error = 0;
};

struct PointDetectorResult {
    std::array<double, 3> rgb{};
    std::array<double, 3> estimated_error{};
    PointDetectorStatistics statistics;
};

struct PointDetectorError {
    PointDetectorFailure reason;
    PointDetectorStatistics statistics;
};

using PointDetectorSampler = std::function<std::expected<PointDetectorProbe, PointDetectorFailure>(
    const DetectorCoordinate&)>;

// Several original camera samples at the SAME pupil, expressed in one smooth
// film chart. Each keeps q = centre + chart_from_standard*z, |z| <= 4.
// Their kernels are neither merged nor replaced by the discovery envelope.
struct PointDetectorFootprint {
    DetectorCoordinate centre{};
    core::AngularMatrix2 chart_from_standard{{{1, 0}, {0, 1}}};
};

struct PointDetectorBatchResult {
    struct Sample {
        std::array<double, 3> rgb{};
        std::array<double, 3> estimated_error{};
    };
    std::vector<Sample> samples;
    PointDetectorStatistics statistics;
};

// Finite adaptive estimator, not an enclosure of every possible source map.
// The callback must trace the actual continuous camera family, including at
// candidate images. Failure has no RGB payload and cannot publish a partial sum.
[[nodiscard]] std::expected<PointDetectorResult, PointDetectorError> EvaluatePointDetector(
    const core::StarfieldSpatialIndex& catalogue, double brightness_scale,
    const PointDetectorSampler& sample, const std::function<bool()>& cancelled,
    const PointDetectorPolicy& policy = {});

// Shared finite image discovery for up to kPointDetectorBatchCapacity samples. The
// callback traces actual chart coordinates and returns derivatives per chart
// unit. Root error is tightened for the narrowest original Gaussian; entirely
// invisible regions retain its spatial sampling depth. Every returned sample
// must independently satisfy the original radiance allowance. Failure carries
// no partial batch; except for cancellation, callers may retry the smaller
// original packets when their common discovery envelope cannot be resolved.
[[nodiscard]] std::expected<PointDetectorBatchResult, PointDetectorError>
EvaluatePointDetectorBatch(const core::StarfieldSpatialIndex& catalogue, double brightness_scale,
                           std::span<const PointDetectorFootprint> footprints,
                           const PointDetectorSampler& sample,
                           const std::function<bool()>& cancelled,
                           const PointDetectorPolicy& policy = {});

// Try common discovery, bisect declined regions along their widest film axis,
// and finish unresolved leaves with the original scalar detector. The shared
// attempts consume at most max(1024, 128*sample_count) extra probes in total;
// exhaustion of that scheduling budget never relaxes an original leaf's policy.
// Cancellation or any failed original footprint withholds the entire result.
[[nodiscard]] std::expected<PointDetectorBatchResult, PointDetectorError>
EvaluatePointDetectorGroup(const core::StarfieldSpatialIndex& catalogue, double brightness_scale,
                           std::span<const PointDetectorFootprint> footprints,
                           const PointDetectorSampler& sample,
                           const std::function<bool()>& cancelled,
                           const PointDetectorPolicy& policy = {});

}  // namespace sirius::render
