#pragma once

#include "sirius/core/trace_boundary.h"

namespace sirius::core {

// Hairer's DP continuous extension, expressed in retained increments:
// x(s) = origin + s * (increment + (1-s) * (a + s * (b + (1-s) * c))).
// This rounded position curve locates events; its tangent does not replace a
// separately admitted physical tangent supplied by the retained sampler.
struct DopriPositionSegment {
    Vec4 origin;
    Vec4 increment;
    Vec4 a;
    Vec4 b;
    Vec4 c;
    double interval;
    double parameter_limit = 1.0;

    [[nodiscard]] Vec4 Displacement(double fraction) const {
        SIRIUS_PRE(std::isfinite(fraction) && fraction >= 0.0 && fraction <= parameter_limit);
        return (increment + (a + (b + c * (1.0 - fraction)) * fraction) * (1.0 - fraction)) *
               fraction;
    }

    [[nodiscard]] AcceptedTraceSegmentSample Sample(double fraction) const {
        SIRIUS_PRE(std::isfinite(interval) && interval > 0.0);
        const auto displacement = Displacement(fraction);
        const Vec4 derivative = increment + a * (1.0 - 2.0 * fraction) +
                                b * (fraction * (2.0 - 3.0 * fraction)) +
                                c * (2.0 * fraction * (1.0 - fraction) * (1.0 - 2.0 * fraction));
        return {origin + displacement, derivative / interval, fraction};
    }

    // Fractions remain relative to the original accepted step. No coefficient
    // or interval is changed when an earlier event limits the search.
    [[nodiscard]] DopriPositionSegment Restricted(double limit) const {
        SIRIUS_PRE(std::isfinite(limit) && limit >= 0.0 && limit <= parameter_limit);
        auto restricted = *this;
        restricted.parameter_limit = limit;
        return restricted;
    }

    [[nodiscard]] bool IsFinite() const {
        if (!(std::isfinite(interval) && interval > 0.0 && std::isfinite(parameter_limit) &&
              parameter_limit >= 0.0 && parameter_limit <= 1.0)) {
            return false;
        }
        for (int component = 0; component < 4; ++component) {
            for (const auto* value : {&origin, &increment, &a, &b, &c}) {
                if (!std::isfinite((*value)(component))) return false;
            }
        }
        return true;
    }
};

namespace detail {

[[nodiscard]] inline std::array<Vec4, 5> DopriPositionCoefficients(
    const DopriPositionSegment& segment) {
    return {segment.origin, segment.increment + segment.a, segment.b + segment.c - segment.a,
            -segment.b - segment.c * 2.0, segment.c};
}

template <std::size_t Size>
[[nodiscard]] inline PolynomialRoots<Size - 1> FindDopriRoots(
    const std::array<double, Size>& coefficients, double limit) {
    PolynomialRoots<Size - 1> roots;
    if (!(std::isfinite(limit) && limit > 0.0 && limit <= 1.0)) return roots;
    // Isolate on a unit interval, then restore original step fractions.
    // Scaling the coefficients also avoids overflow in the derivative chain.
    double scale = 0.0;
    for (double coefficient : coefficients) {
        if (!std::isfinite(coefficient)) return roots;
        scale = std::max(scale, std::abs(coefficient));
    }
    if (scale == 0.0) return roots;  // A coplanar curve has no isolated roots.
    std::array<double, Size> restricted{};
    double power = 1.0;
    for (std::size_t index = 0; index < Size; ++index) {
        restricted[index] = (coefficients[index] / scale) * power;
        power *= limit;
    }
    roots = FindPolynomialRootsOnUnitInterval(restricted, static_cast<int>(Size - 1));
    for (int index = 0; index < roots.count; ++index) {
        roots.values[static_cast<std::size_t>(index)] *= limit;
    }
    return roots;
}

[[nodiscard]] inline std::optional<AcceptedTraceSegmentSample> FindDopriQuadricEvent(
    const DopriPositionSegment& segment, const std::array<double, 3>& axes,
    SphericalBoundarySense sense) {
    if (!segment.IsFinite()) return std::nullopt;
    switch (sense) {
        case SphericalBoundarySense::AnyContact:
        case SphericalBoundarySense::IncreasingRadius:
        case SphericalBoundarySense::DecreasingRadius:
            break;
        default:
            return std::nullopt;
    }
    for (double axis : axes) {
        if (!(std::isfinite(axis) && axis > 0.0)) return std::nullopt;
    }
    const auto position = DopriPositionCoefficients(segment);
    std::array<std::array<double, 5>, 3> spatial{};
    int dominant_axis = 0;
    for (int axis = 0; axis < 3; ++axis) {
        for (int power = 0; power <= 4; ++power) {
            spatial[axis][power] = position[power](axis + 1) / axes[axis];
        }
        if (std::abs(spatial[axis][0]) > std::abs(spatial[dominant_axis][0])) dominant_axis = axis;
    }
    std::array<double, 9> coefficients{};
    // Factor the largest origin square minus one to retain a small radial
    // offset at a large origin, rather than subtracting two large squares.
    const double origin = std::abs(spatial[dominant_axis][0]);
    coefficients[0] = (origin - 1.0) * (origin + 1.0);
    for (int axis = 0; axis < 3; ++axis) {
        for (int left = 0; left <= 4; ++left) {
            for (int right = 0; right <= 4; ++right) {
                if (left == 0 && right == 0 && axis == dominant_axis) continue;
                coefficients[static_cast<std::size_t>(left + right)] +=
                    spatial[axis][left] * spatial[axis][right];
            }
        }
    }
    const auto roots = FindDopriRoots(coefficients, segment.parameter_limit);
    std::array<double, 9> derivative{};
    double derivative_scale = 0.0;
    for (int power = 1; power <= 8; ++power) {
        derivative[power - 1] = static_cast<double>(power) * coefficients[power];
        derivative_scale += std::abs(derivative[power - 1]);
    }
    if (!std::isfinite(derivative_scale)) return std::nullopt;
    const double direction_tolerance =
        512.0 * std::numeric_limits<double>::epsilon() * derivative_scale;
    for (int index = 0; index < roots.count; ++index) {
        const double fraction = roots.values[static_cast<std::size_t>(index)];
        const double direction = EvaluatePolynomial(derivative, 7, fraction);
        if (sense == SphericalBoundarySense::IncreasingRadius &&
            !(direction > direction_tolerance)) {
            continue;
        }
        if (sense == SphericalBoundarySense::DecreasingRadius &&
            !(direction < -direction_tolerance)) {
            continue;
        }
        auto sample = segment.Sample(fraction);
        bool finite = true;
        for (int component = 0; component < 4; ++component) {
            finite = finite && std::isfinite(sample.position(component)) &&
                     std::isfinite(sample.tangent(component));
        }
        if (finite) return sample;
    }
    return std::nullopt;
}

}  // namespace detail

// The squared spherical/ellipsoidal radius of a quartic curve has degree eight.
// Complete derivative-root partitioning includes same-side entries/exits and
// tangent contacts; directional modes exclude numerically stationary contacts.
[[nodiscard]] inline std::optional<AcceptedTraceSegmentSample> FindSphericalBoundaryEvent(
    const DopriPositionSegment& segment, double radius, SphericalBoundarySense sense) {
    return detail::FindDopriQuadricEvent(segment, {radius, radius, radius}, sense);
}

[[nodiscard]] inline std::optional<AcceptedTraceSegmentSample> FindKerrEllipsoidBoundaryEvent(
    const DopriPositionSegment& segment, double radius, double spin, SphericalBoundarySense sense) {
    if (!std::isfinite(spin)) return std::nullopt;
    const double equatorial_axis = std::hypot(radius, spin);
    return detail::FindDopriQuadricEvent(segment, {equatorial_axis, equatorial_axis, radius},
                                         sense);
}

using DopriDiskPlaneRoots = detail::PolynomialRoots<4>;

// Equatorial disk-plane contacts, including isolated tangencies. A curve
// contained in the plane has no isolated event. Results use original fractions.
[[nodiscard]] inline DopriDiskPlaneRoots FindDiskPlaneRoots(const DopriPositionSegment& segment) {
    if (!segment.IsFinite()) return {};
    const auto position = detail::DopriPositionCoefficients(segment);
    std::array<double, 5> coefficients{};
    for (int power = 0; power <= 4; ++power) coefficients[power] = position[power](3);
    return detail::FindDopriRoots(coefficients, segment.parameter_limit);
}

}  // namespace sirius::core
