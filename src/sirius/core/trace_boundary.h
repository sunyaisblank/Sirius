#pragma once

// Numerical event localisation for a finite outer trace boundary. This is a
// core accepted-segment operation so a backend never depends upward on render
// orchestration merely to impose a causal-patch boundary condition.

#include "sirius/base/contracts.h"
#include "sirius/core/tensor.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <optional>
#include <utility>

namespace sirius::core {

struct AcceptedTraceSegmentSample {
    Vec4 position;
    Vec4 tangent;
    double fraction;
};

namespace detail {

struct UnitIntervalPolynomialRoots {
    std::array<double, 6> values{};
    int count = 0;
};

[[nodiscard]] inline double EvaluatePolynomial(const std::array<double, 7>& coefficients,
                                               int degree, double argument) {
    double value = coefficients[degree];
    for (int power = degree - 1; power >= 0; --power) {
        value = value * argument + coefficients[power];
    }
    return value;
}

inline void AppendUnitRoot(UnitIntervalPolynomialRoots& roots, double root,
                           double fraction_tolerance) {
    root = std::clamp(root, 0.0, 1.0);
    if (roots.count > 0 &&
        std::abs(root - roots.values[static_cast<std::size_t>(roots.count - 1)]) <=
            fraction_tolerance) {
        return;
    }
    if (roots.count < static_cast<int>(roots.values.size())) {
        roots.values[static_cast<std::size_t>(roots.count++)] = root;
    }
}

// Isolate every real root of a degree-at-most-six polynomial on [0,1].  The
// derivative roots partition the interval into monotone pieces; sign changes
// give transverse roots and a zero at a derivative root gives a tangent
// contact.  The recursion terminates at a linear polynomial and therefore does
// not assume that a finite sample grid happens to hit a narrow crossing.
[[nodiscard]] inline UnitIntervalPolynomialRoots FindPolynomialRootsOnUnitInterval(
    const std::array<double, 7>& coefficients, int degree) {
    UnitIntervalPolynomialRoots roots;
    double scale = 0.0;
    for (int power = 0; power <= degree; ++power) {
        scale = std::max(scale, std::abs(coefficients[power]));
    }
    const double value_tolerance = 256.0 * std::numeric_limits<double>::epsilon() *
                                   std::max(scale, std::numeric_limits<double>::min());
    const double fraction_tolerance = 512.0 * std::numeric_limits<double>::epsilon();

    while (degree > 0 && std::abs(coefficients[degree]) <= value_tolerance) --degree;
    if (degree == 0) return roots;
    if (degree == 1) {
        const double root = -coefficients[0] / coefficients[1];
        if (std::isfinite(root) && root >= -fraction_tolerance &&
            root <= 1.0 + fraction_tolerance) {
            AppendUnitRoot(roots, root, fraction_tolerance);
        }
        return roots;
    }

    std::array<double, 7> derivative{};
    for (int power = 1; power <= degree; ++power) {
        derivative[power - 1] = static_cast<double>(power) * coefficients[power];
    }
    const auto critical = FindPolynomialRootsOnUnitInterval(derivative, degree - 1);

    std::array<double, 8> boundaries{};
    int boundary_count = 1;
    boundaries[0] = 0.0;
    for (int index = 0; index < critical.count; ++index) {
        const double point = critical.values[static_cast<std::size_t>(index)];
        if (point > fraction_tolerance && point < 1.0 - fraction_tolerance) {
            boundaries[static_cast<std::size_t>(boundary_count++)] = point;
        }
    }
    boundaries[static_cast<std::size_t>(boundary_count++)] = 1.0;

    for (int interval = 0; interval + 1 < boundary_count; ++interval) {
        const double lower = boundaries[static_cast<std::size_t>(interval)];
        const double upper = boundaries[static_cast<std::size_t>(interval + 1)];
        const double lower_value = EvaluatePolynomial(coefficients, degree, lower);
        const double upper_value = EvaluatePolynomial(coefficients, degree, upper);

        if (std::abs(lower_value) <= value_tolerance) {
            AppendUnitRoot(roots, lower, fraction_tolerance);
        }
        if ((lower_value < -value_tolerance && upper_value > value_tolerance) ||
            (lower_value > value_tolerance && upper_value < -value_tolerance)) {
            double root_lower = lower;
            double root_upper = upper;
            double root_lower_value = lower_value;
            for (int iteration = 0; iteration < 64; ++iteration) {
                const double midpoint = root_lower + 0.5 * (root_upper - root_lower);
                const double midpoint_value = EvaluatePolynomial(coefficients, degree, midpoint);
                if ((root_lower_value < 0.0) == (midpoint_value < 0.0)) {
                    root_lower = midpoint;
                    root_lower_value = midpoint_value;
                } else {
                    root_upper = midpoint;
                }
            }
            AppendUnitRoot(roots, root_lower + 0.5 * (root_upper - root_lower), fraction_tolerance);
        }
    }
    const double final_value = EvaluatePolynomial(coefficients, degree, 1.0);
    if (std::abs(final_value) <= value_tolerance) {
        AppendUnitRoot(roots, 1.0, fraction_tolerance);
    }
    return roots;
}

[[nodiscard]] inline std::array<double, 7> SphericalSegmentPolynomial(
    const Vec4& start_position, const Vec4& start_tangent, const Vec4& end_position,
    const Vec4& end_tangent, double affine_length, double boundary_radius) {
    // Each spatial Hermite component is a cubic a*s^3+b*s^2+c*s+d.
    double spatial[3][4]{};
    for (int axis = 0; axis < 3; ++axis) {
        const double start = start_position(axis + 1);
        const double finish = end_position(axis + 1);
        const double start_delta = affine_length * start_tangent(axis + 1);
        const double finish_delta = affine_length * end_tangent(axis + 1);
        spatial[axis][0] = start;
        spatial[axis][1] = start_delta;
        spatial[axis][2] = -3.0 * start + 3.0 * finish - 2.0 * start_delta - finish_delta;
        spatial[axis][3] = 2.0 * start - 2.0 * finish + start_delta + finish_delta;
    }

    std::array<double, 7> coefficients{};
    coefficients[0] = -boundary_radius * boundary_radius;
    for (const auto& component : spatial) {
        for (int left = 0; left <= 3; ++left) {
            for (int right = 0; right <= 3; ++right) {
                coefficients[static_cast<std::size_t>(left + right)] +=
                    component[left] * component[right];
            }
        }
    }
    return coefficients;
}

}  // namespace detail

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
        fraction,
    };
}

// Locate the observer-nearest contact with a spherical capture surface on one
// accepted cubic-Hermite segment.  Squared radius minus R^2 is a sextic.  Its
// complete derivative-root partition above detects entries, exits, and tangent
// contacts even when both segment endpoints lie outside the sphere.
[[nodiscard]] inline std::optional<AcceptedTraceSegmentSample> FindSphericalCaptureEvent(
    const Vec4& start_position, const Vec4& start_tangent, const Vec4& end_position,
    const Vec4& end_tangent, double affine_length, double capture_radius) {
    if (!(std::isfinite(affine_length) && affine_length > 0.0 && std::isfinite(capture_radius) &&
          capture_radius > 0.0)) {
        return std::nullopt;
    }
    for (int component = 0; component < 4; ++component) {
        if (!std::isfinite(start_position(component)) || !std::isfinite(start_tangent(component)) ||
            !std::isfinite(end_position(component)) || !std::isfinite(end_tangent(component))) {
            return std::nullopt;
        }
    }

    const auto coefficients = detail::SphericalSegmentPolynomial(
        start_position, start_tangent, end_position, end_tangent, affine_length, capture_radius);
    const auto roots = detail::FindPolynomialRootsOnUnitInterval(coefficients, 6);
    if (roots.count == 0) return std::nullopt;

    const double fraction = roots.values[0];
    auto event = SampleAcceptedTraceSegment(start_position, start_tangent, end_position,
                                            end_tangent, affine_length, fraction);
    const double sampled_radius =
        std::hypot(event.position(1), event.position(2), event.position(3));
    if (!(std::isfinite(sampled_radius) && sampled_radius > 0.0)) return std::nullopt;
    const double surface_scale = capture_radius / sampled_radius;
    for (int component = 1; component < 4; ++component) {
        event.position(component) *= surface_scale;
    }
    return event;
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
