#pragma once

#include "sirius/core/tensor.h"

#include "radial_quadrature.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <expected>
#include <numeric>

namespace sirius::core::relativity {

enum class SphericalKerrSchildChart { Ingoing, Outgoing };

// Derivatives of Cartesian coordinates and coordinate tangent, with respect
// to one launch parameter at the input event.
struct SphericalTailVariation {
    Vec4 position;
    Vec4 tangent;
};

struct SphericalTailResult {
    Vec4 position;
    Vec4 tangent;
    double affine_length;
    std::array<SphericalTailVariation, 4> event_variations;
    std::array<double, 4> affine_variations;
    unsigned panels;
    // Analytic Simpson bounds for the integrals and their impact-squared
    // derivatives at the represented scalar invariants. The complete event
    // and coordinate transforms have separate error obligations.
    std::array<RadialIntegralBound, 6> radial_integral_bounds;
};

namespace finite_tail_detail {
using Jet = Dual<double>;
using Vector = std::array<Jet, 3>;

inline bool Finite(const Jet& x) { return std::isfinite(x.real) && std::isfinite(x.dual); }

inline Jet Dot(const Vector& a, const Vector& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }

struct Momentum {
    Vector normal;
    Vector angular;
    Jet radius;
    Jet frequency;
    Jet impact_squared;
    Jet coordinate_time;
};

inline std::expected<Momentum, FiniteTailFailure> Extract(const Vec4& position, const Vec4& tangent,
                                                          const SphericalTailVariation& column,
                                                          double mass, double orientation,
                                                          double sources, double boundary,
                                                          double tolerance) {
    const Vector x{Jet(position(1), column.position(1)), Jet(position(2), column.position(2)),
                   Jet(position(3), column.position(3))};
    const Vector k{Jet(tangent(1), column.tangent(1)), Jet(tangent(2), column.tangent(2)),
                   Jet(tangent(3), column.tangent(3))};
    const Jet radius = sqrt(Dot(x, x));
    // For r>3M the centrifugal potential (1-2M/r)/r² decreases strictly.
    // An outward solution admitted here cannot turn anywhere before the
    // requested finite sphere and cannot re-enter the source support.
    if (!Finite(radius) || !(radius.real > std::max(3.0 * mass, sources)) ||
        !(radius.real < boundary))
        return std::unexpected(FiniteTailFailure::Domain);
    const Jet h = Jet(2.0 * mass) / radius;
    const Jet f = Jet(1.0) - h;
    Vector n;
    for (unsigned i = 0; i < 3; ++i) n[i] = x[i] / radius;
    const Jet radial = Dot(n, k);
    const Jet frequency = h * radial * orientation - f * Jet(tangent(0), column.tangent(0));
    if (!Finite(frequency) || !Finite(radial) || !(frequency.real > 0.0) || !(radial.real > 0.0))
        return std::unexpected(FiniteTailFailure::Domain);
    Vector transverse;
    for (unsigned i = 0; i < 3; ++i) transverse[i] = k[i] - radial * n[i];
    const Jet temporal = frequency * frequency;
    const Jet longitudinal = radial * radial;
    const Jet angular = f * Dot(transverse, transverse);
    const Jet constraint = temporal - longitudinal - angular;
    const double scale = temporal.real + longitudinal.real + angular.real;
    const double variation_scale =
        std::abs(temporal.dual) + std::abs(longitudinal.dual) + std::abs(angular.dual);
    // Small relative null defects can still move the outward radial root
    // substantially near a turn. Bound the defect against radial speed².
    if (!Finite(constraint) || !std::isfinite(scale) || !(scale > 0.0) ||
        !std::isfinite(variation_scale) ||
        std::abs(constraint.real) > tolerance * longitudinal.real ||
        std::abs(constraint.dual) > tolerance * std::max(scale, variation_scale))
        return std::unexpected(FiniteTailFailure::Constraint);
    Vector impact;
    for (unsigned i = 0; i < 3; ++i) impact[i] = transverse[i] * radius / frequency;
    const Jet impact_squared = Dot(impact, impact);
    if (!Finite(impact_squared)) return std::unexpected(FiniteTailFailure::Arithmetic);
    return Momentum{n,         impact,         radius,
                    frequency, impact_squared, Jet(position(0), column.position(0))};
}

struct Rotation {
    Jet cosine;
    Jet scaled_sine;
};

inline Rotation AngularRotation(Jet impact_squared, Jet angular_integral) {
    const Jet square = impact_squared * angular_integral * angular_integral;
    // Evaluate cos(sqrt(z)) and sinc(sqrt(z)) as functions of z near zero.
    // This preserves the first variation of an exactly radial launch.
    if (square.real < 1e-6) {
        const Jet cosine =
            Jet(1.0) + square * (Jet(-0.5) + square * (Jet(1.0 / 24.0) - square / 720.0));
        const Jet sinc =
            Jet(1.0) + square * (Jet(-1.0 / 6.0) + square * (Jet(1.0 / 120.0) - square / 5040.0));
        return {cosine, angular_integral * sinc};
    }
    const Jet angle = sqrt(square);
    return {cos(angle), angular_integral * sin(angle) / angle};
}
}  // namespace finite_tail_detail

// Private continuation candidate for exactly spherical, uncharged vacuum
// Kerr-Schild geometry. The caller supplies the complete radius of all source
// support. The output variations are at the common finite radius, not at a
// common affine time. Affine variation is returned separately for that change
// of parameter. The six radial quadratures carry their own enclosures.
[[nodiscard]] inline std::expected<SphericalTailResult, FiniteTailFailure> ContinueSphericalTail(
    const Vec4& position, const Vec4& tangent, const std::array<SphericalTailVariation, 4>& columns,
    double mass, SphericalKerrSchildChart chart, double source_support_radius, double boundary,
    const SphericalTailControl& control = {}) {
    using namespace finite_tail_detail;
    if (!RoundingEnvironment() || !std::isfinite(mass) || mass < 0.0 ||
        !std::isfinite(source_support_radius) || source_support_radius < 0.0 ||
        !std::isfinite(boundary) || !std::isfinite(control.absolute_tolerance) ||
        !(control.absolute_tolerance > 0.0) || !std::isfinite(control.relative_tolerance) ||
        control.relative_tolerance < 0.0 || !std::isfinite(control.constraint_tolerance) ||
        control.constraint_tolerance < 0.0)
        return std::unexpected(FiniteTailFailure::Input);
    for (int i = 0; i < 4; ++i) {
        if (!std::isfinite(position(i)) || !std::isfinite(tangent(i)))
            return std::unexpected(FiniteTailFailure::Input);
        for (const auto& column : columns)
            if (!std::isfinite(column.position(i)) || !std::isfinite(column.tangent(i)))
                return std::unexpected(FiniteTailFailure::Input);
    }
    const double orientation = chart == SphericalKerrSchildChart::Ingoing ? 1.0 : -1.0;
    std::array<Momentum, 4> momenta;
    for (unsigned column = 0; column < columns.size(); ++column) {
        const auto extracted =
            Extract(position, tangent, columns[column], mass, orientation, source_support_radius,
                    boundary, control.constraint_tolerance);
        if (!extracted) return std::unexpected(extracted.error());
        momenta[column] = *extracted;
    }
    const double start = momenta[0].radius.real;
    const double b2 = momenta[0].impact_squared.real;
    const auto a = RadialAtLaunch(start, mass, orientation, b2);
    Quadrature quadrature{mass, orientation, b2, start, boundary, control};
    const auto integrated = quadrature.Run();
    if (!integrated) return std::unexpected(integrated.error());

    SphericalTailResult result{};
    result.panels = quadrature.panels;
    std::array<double, 6> estimates;
    for (unsigned i = 0; i < estimates.size(); ++i) {
        const auto enclosure = (*integrated)[i];
        estimates[i] = std::midpoint(enclosure.lo, enclosure.hi);
        const double budget =
            control.absolute_tolerance + control.relative_tolerance * std::abs(estimates[i]);
        if (MidpointError(enclosure) > budget)
            return std::unexpected(FiniteTailFailure::Arithmetic);
        result.radial_integral_bounds[i] = {enclosure.lo, enclosure.hi};
    }
    const double h = 2.0 * mass / boundary;
    const double f = 1.0 - h;
    for (unsigned index = 0; index < momenta.size(); ++index) {
        const auto& input = momenta[index];
        const auto integral = [&](unsigned i) {
            // Leibniz rule for a moving launch radius and fixed target radius.
            return Jet(estimates[i],
                       estimates[i + 3] * input.impact_squared.dual - a[i] * input.radius.dual);
        };
        const Jet affine = integral(0) / input.frequency;
        const Jet time = input.coordinate_time + integral(2);
        const Jet radial =
            input.frequency * sqrt(Jet(1.0) - input.impact_squared * (f / boundary / boundary));
        const Jet temporal = (radial * (orientation * h) - input.frequency) / f;
        const auto rotation = AngularRotation(input.impact_squared, integral(1));
        std::array<Jet, 4> x, k;
        x[0] = time;
        k[0] = temporal;
        for (unsigned i = 0; i < 3; ++i) {
            const Jet n =
                input.normal[i] * rotation.cosine + input.angular[i] * rotation.scaled_sine;
            const Jet angular = input.angular[i] * rotation.cosine -
                                input.normal[i] * input.impact_squared * rotation.scaled_sine;
            x[i + 1] = n * boundary;
            k[i + 1] = n * radial + angular * input.frequency / boundary;
        }
        if (!Finite(affine) || !(affine.real > 0.0) || !Finite(radial) || !(radial.real > 0.0))
            return std::unexpected(FiniteTailFailure::Arithmetic);
        if (index == 0) result.affine_length = affine.real;
        result.affine_variations[index] = affine.dual;
        for (unsigned i = 0; i < 4; ++i) {
            if (!Finite(x[i]) || !Finite(k[i]))
                return std::unexpected(FiniteTailFailure::Arithmetic);
            if (index == 0) {
                result.position(i) = x[i].real;
                result.tangent(i) = k[i].real;
            }
            result.event_variations[index].position(i) = x[i].dual;
            result.event_variations[index].tangent(i) = k[i].dual;
        }
    }
    return result;
}
}  // namespace sirius::core::relativity
