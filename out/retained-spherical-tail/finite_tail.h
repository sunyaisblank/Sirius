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

struct SphericalTailEventBounds {
    std::array<SphericalTailBound, 4> position;
    std::array<SphericalTailBound, 4> tangent;
    SphericalTailBound affine_length;
};

struct SphericalTailResult {
    Vec4 position;
    Vec4 tangent;
    double affine_length;
    std::array<SphericalTailVariation, 4> event_variations;
    std::array<double, 4> affine_variations;
    unsigned panels;
    std::array<SphericalTailBound, 6> radial_integral_bounds;
    SphericalTailEventBounds event_bounds;
    std::array<SphericalTailEventBounds, 4> variation_bounds;
};

namespace finite_tail_detail {
struct Jet {
    Range real;
    Range dual;

    Jet() = default;
    explicit Jet(double value) : real(value) {}
    explicit Jet(Range value, Range derivative = {}) : real(value), dual(derivative) {}
};

inline Jet operator+(Jet a, Jet b) { return Jet(a.real + b.real, a.dual + b.dual); }
inline Jet operator-(Jet a, Jet b) { return Jet(a.real - b.real, a.dual - b.dual); }
inline Jet operator*(Jet a, Jet b) {
    return Jet(a.real * b.real, a.dual * b.real + a.real * b.dual);
}
inline Jet operator/(Jet a, Jet b) {
    const Range quotient = a.real / b.real;
    return Jet(quotient, (a.dual - quotient * b.dual) / b.real);
}
inline Jet sqrt(Jet a) {
    const Range root = sqrt(a.real);
    return Jet(root, a.dual / (Range(2.0) * root));
}
inline bool Finite(Jet a) { return a.real.Finite() && a.dual.Finite(); }

using Vector = std::array<Jet, 3>;
using Four = std::array<Range, 4>;
using Columns = std::array<std::array<Jet, 4>, 4>;

inline Jet Dot(const Vector& a, const Vector& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }

inline Jet NormSquared(const Vector& a) {
    auto value = Square(a[0].real) + Square(a[1].real) + Square(a[2].real);
    value.lo = std::max(0.0, value.lo);
    const auto derivative =
        Range(2.0) * (a[0].real * a[0].dual + a[1].real * a[1].dual + a[2].real * a[2].dual);
    return Jet(value, derivative);
}

struct Momentum {
    Vector normal;
    Vector angular;
    Jet radius;
    Jet frequency;
    Jet impact_squared;
    Jet coordinate_time;
};

inline std::expected<Momentum, FiniteTailFailure> Extract(const std::array<Jet, 4>& position,
                                                          const std::array<Jet, 4>& tangent,
                                                          double mass, double orientation,
                                                          double sources, double boundary,
                                                          double tolerance) {
    const Vector x{position[1], position[2], position[3]};
    const Vector k{tangent[1], tangent[2], tangent[3]};
    const Jet radius = sqrt(NormSquared(x));
    // For r>3M the centrifugal potential (1-2M/r)/r² decreases strictly.
    // An outward solution admitted here cannot turn anywhere before the
    // requested finite sphere and cannot re-enter the source support.
    if (!Finite(radius) || !(radius.real.lo > std::max((Range(3.0) * Range(mass)).hi, sources)) ||
        !(radius.real.hi < boundary))
        return std::unexpected(FiniteTailFailure::Domain);
    const Jet h = Jet(2.0) * Jet(mass) / radius;
    const Jet f = Jet(1.0) - h;
    Vector n;
    for (unsigned i = 0; i < 3; ++i) n[i] = x[i] / radius;
    const Jet radial = Dot(n, k);
    const Jet frequency = h * radial * Jet(orientation) - f * tangent[0];
    if (!Finite(frequency) || !Finite(radial) || !(frequency.real.lo > 0.0) ||
        !(radial.real.lo > 0.0))
        return std::unexpected(FiniteTailFailure::Domain);
    Vector transverse;
    for (unsigned i = 0; i < 3; ++i) transverse[i] = k[i] - radial * n[i];
    const Jet temporal = frequency * frequency;
    const Jet longitudinal = radial * radial;
    const Jet angular = f * NormSquared(transverse);
    const Jet constraint = temporal - longitudinal - angular;
    const auto scale = temporal.real + longitudinal.real + angular.real;
    const auto variation_scale = Abs(temporal.dual) + Abs(longitudinal.dual) + Abs(angular.dual);
    // Small relative null defects can still move the outward radial root
    // substantially near a turn. Bound the defect against radial speed².
    if (!Finite(constraint) || !scale.Finite() || !variation_scale.Finite() || !(scale.lo > 0.0) ||
        Magnitude(constraint.real) > (Range(tolerance) * longitudinal.real).lo ||
        Magnitude(constraint.dual) >
            (Range(tolerance) * Range(std::max(scale.lo, variation_scale.lo))).lo)
        return std::unexpected(FiniteTailFailure::Constraint);
    Vector impact;
    for (unsigned i = 0; i < 3; ++i) impact[i] = transverse[i] * radius / frequency;
    const Jet impact_squared = NormSquared(impact);
    if (!Finite(impact_squared)) return std::unexpected(FiniteTailFailure::Arithmetic);
    return Momentum{n, impact, radius, frequency, impact_squared, position[0]};
}

struct Rotation {
    Jet cosine;
    Jet scaled_sine;
};

inline Range RotationSeries(Range z, unsigned factorial_offset, bool derivative = false) {
    // On [0, 4] the tail starting at term 12 is alternating with decreasing
    // magnitudes, including for the differentiated sinc series. Horner's
    // leading interval includes both zero and the first omitted coefficient.
    constexpr unsigned terms = 12;
    std::array<Range, terms + 1> coefficients;
    Range denominator(1.0);
    unsigned factorial = 1;
    for (unsigned i = 0; i < coefficients.size(); ++i) {
        while (factorial < 2 * i + factorial_offset) {
            ++factorial;
            denominator = denominator * Range(factorial);
        }
        const double sign = (i + static_cast<unsigned>(derivative)) % 2 == 0 ? 1.0 : -1.0;
        coefficients[i] = Range(sign * (derivative ? i + 1.0 : 1.0)) / denominator;
    }
    Range result(std::min(0.0, coefficients.back().lo), std::max(0.0, coefficients.back().hi));
    for (unsigned i = terms; i > 0; --i) result = coefficients[i - 1] + z * result;
    return result;
}

inline std::expected<Rotation, FiniteTailFailure> AngularRotation(Jet impact_squared,
                                                                  Jet angular_integral) {
    Jet square = impact_squared * angular_integral * angular_integral;
    square.real.lo = std::max(0.0, square.real.lo);
    if (!Finite(square) || square.real.hi > 4.0) return std::unexpected(FiniteTailFailure::Domain);
    const Range cosine = RotationSeries(square.real, 0);
    const Range sinc = RotationSeries(square.real, 1);
    const Range sinc_derivative = RotationSeries(square.real, 3, true);
    return Rotation{Jet(cosine, Range(-0.5) * sinc * square.dual),
                    angular_integral * Jet(sinc, sinc_derivative * square.dual)};
}

inline SphericalTailBound Bound(Range value) { return {value.lo, value.hi}; }

inline bool Within(Range value, double absolute, double relative) {
    if (!value.Finite()) return false;
    const double estimate = std::midpoint(value.lo, value.hi);
    const auto budget = Range(absolute) + Range(relative) * Range(std::abs(estimate));
    return MidpointError(value) <= budget.lo;
}

inline std::expected<SphericalTailResult, FiniteTailFailure> Continue(
    const Columns& position, const Columns& tangent, double mass, SphericalKerrSchildChart chart,
    double source_support_radius, double boundary, const SphericalTailControl& control) {
    if (Cancelled(control)) return std::unexpected(FiniteTailFailure::Cancelled);
    if (!RoundingEnvironment() || !std::isfinite(mass) || mass < 0.0 ||
        !std::isfinite(source_support_radius) || source_support_radius < 0.0 ||
        !std::isfinite(boundary) || !std::isfinite(control.absolute_tolerance) ||
        !(control.absolute_tolerance > 0.0) || !std::isfinite(control.relative_tolerance) ||
        control.relative_tolerance < 0.0 || !std::isfinite(control.constraint_tolerance) ||
        control.constraint_tolerance < 0.0 || !std::isfinite(control.event_absolute_tolerance) ||
        !(control.event_absolute_tolerance > 0.0) ||
        !std::isfinite(control.event_relative_tolerance) || control.event_relative_tolerance < 0.0)
        return std::unexpected(FiniteTailFailure::Input);
    for (unsigned column = 0; column < position.size(); ++column)
        for (unsigned i = 0; i < position[column].size(); ++i)
            if (!Finite(position[column][i]) || !Finite(tangent[column][i]))
                return std::unexpected(FiniteTailFailure::Input);
    const double orientation = chart == SphericalKerrSchildChart::Ingoing ? 1.0 : -1.0;
    std::array<Momentum, 4> momenta;
    for (unsigned column = 0; column < momenta.size(); ++column) {
        const auto extracted =
            Extract(position[column], tangent[column], mass, orientation, source_support_radius,
                    boundary, control.constraint_tolerance);
        if (!extracted) return std::unexpected(extracted.error());
        momenta[column] = *extracted;
    }
    const Range start = momenta[0].radius.real;
    const Range b2 = momenta[0].impact_squared.real;
    const auto a = RadialAtLaunch(start, mass, orientation, b2);
    Quadrature quadrature{mass, orientation, b2, start, boundary, control};
    const auto integrated = quadrature.Run();
    if (!integrated) return std::unexpected(integrated.error());

    SphericalTailResult result{};
    result.panels = quadrature.panels;
    for (unsigned i = 0; i < integrated->size(); ++i) {
        const auto enclosure = (*integrated)[i];
        if (!Within(enclosure, control.absolute_tolerance, control.relative_tolerance))
            return std::unexpected(FiniteTailFailure::Arithmetic);
        result.radial_integral_bounds[i] = Bound(enclosure);
    }
    const Jet target(boundary);
    const Jet h = Jet(2.0) * Jet(mass) / target;
    const Jet f = Jet(1.0) - h;
    const auto within = [&](Jet value) {
        return Within(value.real, control.event_absolute_tolerance,
                      control.event_relative_tolerance) &&
               Within(value.dual, control.event_absolute_tolerance,
                      control.event_relative_tolerance);
    };
    const auto midpoint = [](Range value) { return std::midpoint(value.lo, value.hi); };
    for (unsigned index = 0; index < momenta.size(); ++index) {
        const auto& input = momenta[index];
        const auto integral = [&](unsigned i) {
            // Leibniz rule for a moving launch radius and fixed target radius.
            return Jet((*integrated)[i],
                       (*integrated)[i + 3] * input.impact_squared.dual - a[i] * input.radius.dual);
        };
        const Jet affine = integral(0) / input.frequency;
        const Jet time = input.coordinate_time + integral(2);
        const Jet radial =
            input.frequency * sqrt(Jet(1.0) - input.impact_squared * (f / target / target));
        const Jet temporal = (radial * Jet(orientation) * h - input.frequency) / f;
        const auto rotation = AngularRotation(input.impact_squared, integral(1));
        if (!rotation) return std::unexpected(rotation.error());
        std::array<Jet, 4> x, k;
        x[0] = time;
        k[0] = temporal;
        for (unsigned i = 0; i < 3; ++i) {
            const Jet n =
                input.normal[i] * rotation->cosine + input.angular[i] * rotation->scaled_sine;
            const Jet angular = input.angular[i] * rotation->cosine -
                                input.normal[i] * input.impact_squared * rotation->scaled_sine;
            x[i + 1] = n * target;
            k[i + 1] = n * radial + angular * input.frequency / target;
        }
        if (!(affine.real.lo > 0.0) || !Finite(radial) || !(radial.real.lo > 0.0))
            return std::unexpected(FiniteTailFailure::Arithmetic);
        if (!within(affine)) return std::unexpected(FiniteTailFailure::EventAccuracy);
        if (index == 0) {
            result.affine_length = midpoint(affine.real);
            result.event_bounds.affine_length = Bound(affine.real);
        }
        result.affine_variations[index] = midpoint(affine.dual);
        result.variation_bounds[index].affine_length = Bound(affine.dual);
        for (unsigned i = 0; i < 4; ++i) {
            if (!within(x[i]) || !within(k[i]))
                return std::unexpected(FiniteTailFailure::EventAccuracy);
            if (index == 0) {
                result.position(i) = midpoint(x[i].real);
                result.tangent(i) = midpoint(k[i].real);
                result.event_bounds.position[i] = Bound(x[i].real);
                result.event_bounds.tangent[i] = Bound(k[i].real);
            }
            result.event_variations[index].position(i) = midpoint(x[i].dual);
            result.event_variations[index].tangent(i) = midpoint(k[i].dual);
            result.variation_bounds[index].position[i] = Bound(x[i].dual);
            result.variation_bounds[index].tangent[i] = Bound(k[i].dual);
        }
    }
    if (Cancelled(control)) return std::unexpected(FiniteTailFailure::Cancelled);
    return result;
}
}  // namespace finite_tail_detail

// Encloses the finite-radius event and launch derivatives of the outward null
// ray defined by the input direction, Killing frequency, and angular momentum.
// The null residual guard limits the difference from the input radial speed.
// The caller supplies the complete radius of all source support.
[[nodiscard]] inline std::expected<SphericalTailResult, FiniteTailFailure> ContinueSphericalTail(
    const Vec4& position, const Vec4& tangent, const std::array<SphericalTailVariation, 4>& columns,
    double mass, SphericalKerrSchildChart chart, double source_support_radius, double boundary,
    const SphericalTailControl& control = {}) {
    using namespace finite_tail_detail;
    Columns x, k;
    for (unsigned column = 0; column < columns.size(); ++column)
        for (unsigned i = 0; i < x[column].size(); ++i) {
            x[column][i] = Jet(Range(position(i)), Range(columns[column].position(i)));
            k[column][i] = Jet(Range(tangent(i)), Range(columns[column].tangent(i)));
        }
    return Continue(x, k, mass, chart, source_support_radius, boundary, control);
}
}  // namespace sirius::core::relativity
