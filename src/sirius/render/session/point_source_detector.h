#pragma once

#include "sirius/core/point_source_response.h"
#include "sirius/core/starfield.h"

#include <array>
#include <cstddef>
#include <expected>
#include <functional>

namespace sirius::render {

// All coordinates belong to one original camera sample and fixed pupil.
// q = q0 + L z, |z| <= 4. Subdivision never changes L or the parent Gaussian.
using DetectorCoordinate = std::array<double, 2>;

struct PointDetectorProbe {
    bool visible = false;
    std::array<double, 3> direction{};
    core::AngularMatrix2 source_derivative{};  // Source radians per original z.
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

// Finite adaptive estimator, not an enclosure of every possible source map.
// The callback must trace the actual continuous camera family, including at
// candidate images. Failure has no RGB payload and cannot publish a partial sum.
[[nodiscard]] std::expected<PointDetectorResult, PointDetectorError> EvaluatePointDetector(
    const core::StarfieldSpatialIndex& catalogue, double brightness_scale,
    const PointDetectorSampler& sample, const std::function<bool()>& cancelled,
    const PointDetectorPolicy& policy = {});

}  // namespace sirius::render
