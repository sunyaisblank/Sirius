#pragma once

// Numerical event localisation for a finite outer trace boundary. This is a
// core accepted-segment operation so a backend never depends upward on render
// orchestration merely to impose a causal-patch boundary condition.

#include "sirius/base/contracts.h"
#include "sirius/core/tensor.h"

#include <cmath>
#include <optional>
#include <utility>

namespace sirius::core {

struct AcceptedTraceSegmentSample {
    Vec4 position;
    Vec4 tangent;
};

[[nodiscard]] inline AcceptedTraceSegmentSample SampleAcceptedTraceSegment(
    const Vec4& start_position, const Vec4& start_tangent, const Vec4& end_position,
    const Vec4& end_tangent, double affine_length, double fraction) {
    SIRIUS_PRE(std::isfinite(affine_length) && affine_length > 0.0);
    SIRIUS_PRE(std::isfinite(fraction) && fraction >= 0.0 && fraction <= 1.0);
    const double s2 = fraction * fraction;
    const double s3 = s2 * fraction;
    const double h00 = 2.0 * s3 - 3.0 * s2 + 1.0;
    const double h10 = s3 - 2.0 * s2 + fraction;
    const double h01 = -2.0 * s3 + 3.0 * s2;
    const double h11 = s3 - s2;
    const double dh00 = 6.0 * s2 - 6.0 * fraction;
    const double dh10 = 3.0 * s2 - 4.0 * fraction + 1.0;
    const double dh01 = -6.0 * s2 + 6.0 * fraction;
    const double dh11 = 3.0 * s2 - 2.0 * fraction;

    return {
        start_position * h00 + start_tangent * (affine_length * h10) + end_position * h01 +
            end_tangent * (affine_length * h11),
        start_position * (dh00 / affine_length) + start_tangent * dh10 +
            end_position * (dh01 / affine_length) + end_tangent * dh11,
    };
}

// Locate the first sampled crossing of a finite outer causal boundary on an
// accepted RK segment. The integrator may evaluate a trial endpoint outside
// the patch, but the caller never needs to impose its boundary condition there.
[[nodiscard]] inline std::optional<AcceptedTraceSegmentSample> FindCausalBoundaryEvent(
    const Vec4& start_position, const Vec4& start_tangent, const Vec4& end_position,
    const Vec4& end_tangent, double affine_length, double boundary_radius) {
    if (!(std::isfinite(affine_length) && affine_length > 0.0 && std::isfinite(boundary_radius) &&
          boundary_radius > 0.0)) {
        return std::nullopt;
    }
    const auto radius_at = [&](double fraction) {
        const auto sample = SampleAcceptedTraceSegment(start_position, start_tangent, end_position,
                                                       end_tangent, affine_length, fraction);
        return std::pair{sample,
                         std::hypot(sample.position(1), sample.position(2), sample.position(3))};
    };
    const auto [start_sample, start_radius] = radius_at(0.0);
    if (!std::isfinite(start_radius)) return std::nullopt;
    if (start_radius >= boundary_radius) return std::nullopt;

    constexpr int kBoundaryBrackets = 64;
    double inside_fraction = 0.0;
    for (int interval = 1; interval <= kBoundaryBrackets; ++interval) {
        const double outside_fraction = static_cast<double>(interval) / kBoundaryBrackets;
        const double outside_radius = radius_at(outside_fraction).second;
        if (!std::isfinite(outside_radius)) return std::nullopt;
        if (outside_radius < boundary_radius) {
            inside_fraction = outside_fraction;
            continue;
        }

        double lower = inside_fraction;
        double upper = outside_fraction;
        for (int iteration = 0; iteration < 64; ++iteration) {
            const double midpoint = lower + 0.5 * (upper - lower);
            if (radius_at(midpoint).second < boundary_radius) {
                lower = midpoint;
            } else {
                upper = midpoint;
            }
        }
        auto boundary_sample = radius_at(upper).first;
        const double sampled_radius = std::hypot(
            boundary_sample.position(1), boundary_sample.position(2), boundary_sample.position(3));
        if (!(std::isfinite(sampled_radius) && sampled_radius > 0.0)) return std::nullopt;
        const double inward_scale = std::nextafter(boundary_radius, 0.0) / sampled_radius;
        for (int component = 1; component < 4; ++component) {
            boundary_sample.position(component) *= inward_scale;
        }
        return boundary_sample;
    }
    return std::nullopt;
}

}  // namespace sirius::core
