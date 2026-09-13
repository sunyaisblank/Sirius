#pragma once

#include "sirius/core/source_sky_map.h"
#include "sirius/core/trace_boundary.h"

#include <algorithm>
#include <expected>

namespace sirius::core::relativity {

// Vacuum Kerr continuation from an outward exterior event to radial infinity.
// These tolerances concern the tail only; the inner trace and the detector
// integration have separate errors. Callers must budget all three.
struct KerrInfinityConfig {
    double absolute_tolerance = 1e-10;
    double relative_tolerance = 1e-9;
    double input_constraint_tolerance = 1e-8;
    unsigned maximum_attempts = 4096;
};

enum class KerrInfinityFailure { InvalidInput, RadialUnresolved, WorkLimit, Arithmetic };

struct KerrInfinityResult {
    SourceSkyAngularMap map;
    double frequency;
    std::array<double, 2> frequency_derivative;
    double maximum_local_error_ratio;
    unsigned accepted_steps;
    unsigned attempted_steps;
};

namespace kerr_infinity_detail {
using Jet = Dual<double>;
using Vector = std::array<Jet, 3>;
using State = std::array<std::array<Jet, 6>, 2>;

inline Jet Dot(const Vector& a, const Vector& b) {
    Jet result;
    for (unsigned i = 0; i < 3; ++i) result += a[i] * b[i];
    return result;
}
inline Vector SpinCross(const Vector& a) { return {-a[1], a[0], Jet{}}; }
inline bool Finite(const Jet& a) { return std::isfinite(a.real) && std::isfinite(a.dual); }

struct Constants {
    Jet inverse_radius;
    Jet energy;
    Jet angular_momentum;
    Jet carter_combination;
};

// n is the oblate angular unit vector in the public ingoing KS orientation.
// v is p_theta e_theta + (L/sin(theta)) e_phi, represented without polar angles.
// The covector projection below is regular at both axes.
inline bool InitialState(double mass, double spin, const Metric4d& metric,
                         const Tensor<Dual<double>, 4, 4, 4>& metric_derivatives,
                         const Vec4& position, const Vec4& tangent, const Vec4& displacement,
                         const Vec4& coordinate_variation, double seed, Constants& constants,
                         std::array<Jet, 6>& state, Jet& frequency, double constraint_tolerance) {
    const double scale = std::max(
        {std::abs(position(1)), std::abs(position(2)), std::abs(position(3)), std::abs(spin)});
    if (!(scale > 0.0) || !std::isfinite(scale)) return false;
    Vector x, k;
    std::array<Jet, 4> momentum{};
    std::array<Jet, 4> four_tangent;
    for (unsigned mu = 0; mu < 4; ++mu)
        four_tangent[mu] = Jet(tangent(mu), coordinate_variation(mu) / seed);
    for (unsigned mu = 0; mu < 4; ++mu) {
        for (unsigned nu = 0; nu < 4; ++nu) {
            Jet value(metric(mu, nu).real);
            for (unsigned axis = 0; axis < 4; ++axis)
                value.dual += metric_derivatives(axis, mu, nu).real * displacement(axis) / seed;
            momentum[mu] += value * four_tangent[nu];
        }
    }
    frequency = momentum[0];  // -E for the positive-frequency past tangent.
    if (!Finite(frequency) || !(frequency.real > 0.0)) return false;
    Jet null_residual;
    double null_scale = 0.0;
    double null_derivative_scale = 0.0;
    for (unsigned mu = 0; mu < 4; ++mu) {
        const auto term = momentum[mu] * four_tangent[mu];
        null_residual += term;
        null_scale += std::abs(term.real);
        null_derivative_scale += std::abs(momentum[mu].dual * four_tangent[mu].real) +
                                 std::abs(momentum[mu].real * four_tangent[mu].dual);
    }
    if (!Finite(null_residual) || !(null_scale > 0.0) ||
        std::abs(null_residual.real) > constraint_tolerance * null_scale ||
        std::abs(null_residual.dual) >
            constraint_tolerance * std::max(null_scale, null_derivative_scale))
        return false;
    // Remove affine-frequency scale, including its derivative, before forming
    // squared constants. Signed E=-1 is retained in the frame-dragging terms.
    for (auto& value : momentum) value /= frequency;
    for (auto& value : four_tangent) value /= frequency;
    for (unsigned i = 0; i < 3; ++i) {
        x[i] = Jet(position(i + 1), displacement(i + 1) / seed);
        k[i] = four_tangent[i + 1];
    }
    Vector scaled;
    for (unsigned i = 0; i < 3; ++i) scaled[i] = x[i] / scale;
    const double a_scaled = spin / scale;
    const Jet difference = Dot(scaled, scaled) - Jet(a_scaled * a_scaled);
    const Jet discriminant =
        sqrt(difference * difference + scaled[2] * scaled[2] * (4.0 * a_scaled * a_scaled));
    const Jet radius = sqrt((difference + discriminant) * 0.5) * scale;
    if (!Finite(radius) ||
        !(radius.real > mass + std::sqrt(std::max(0.0, mass * mass - spin * spin))))
        return false;
    const Jet r2 = radius * radius;
    const Jet denominator = r2 + Jet(spin * spin);
    const Vector n{(radius * x[0] + x[1] * spin) / denominator,
                   (radius * x[1] - x[0] * spin) / denominator, x[2] / radius};
    const Vector p{momentum[1], momentum[2], momentum[3]};
    const Vector w{radius * p[0] + p[1] * spin, radius * p[1] - p[0] * spin, radius * p[2]};
    const Jet projection = Dot(n, w);
    Vector v;
    for (unsigned i = 0; i < 3; ++i) v[i] = w[i] - n[i] * projection;
    // Outgoing oblate radial branch; checking a Euclidean radius is insufficient.
    const Jet radial_derivative =
        (radius * (x[0] * k[0] + x[1] * k[1]) + denominator * x[2] * k[2] / radius) /
        (discriminant * (scale * scale));
    if (!Finite(radial_derivative) || !(radial_derivative.real > 0.0)) return false;
    constants.inverse_radius = Jet(1.0) / radius;
    constants.energy = -momentum[0];
    constants.angular_momentum = n[0] * v[1] - n[1] * v[0];
    constants.carter_combination =
        Dot(v, v) - constants.energy * constants.angular_momentum * (2.0 * spin) +
        constants.energy * constants.energy * (Jet(1.0) - n[2] * n[2]) * (spin * spin);
    for (unsigned i = 0; i < 3; ++i) {
        state[i] = n[i];
        state[i + 3] = v[i];
        if (!Finite(n[i]) || !Finite(v[i])) return false;
    }
    return Finite(constants.carter_combination);
}

// Check all extrema of H on the complete proposed radial interval, not just
// the launch point. A grazing radial root cannot be traversed by u as time.
inline bool OutwardInterval(double mass, double spin, const Constants& c) {
    const double u = c.inverse_radius.real;
    const double e = c.energy.real;
    const double q = c.carter_combination.real;
    const double z = spin * spin * e - spin * c.angular_momentum.real;
    const std::array<double, 7> polynomial{e * e,
                                           0.0,
                                           (2.0 * e * z - q) * u * u,
                                           2.0 * mass * q * u * u * u,
                                           (z * z - spin * spin * q) * u * u * u * u,
                                           0.0,
                                           0.0};
    std::array<double, 7> derivative{};
    double scale = 0.0;
    for (unsigned i = 0; i <= 4; ++i) {
        if (!std::isfinite(polynomial[i])) return false;
        scale += std::abs(polynomial[i]);
        if (i) derivative[i - 1] = static_cast<double>(i) * polynomial[i];
    }
    const auto extrema = sirius::core::detail::FindPolynomialRootsOnUnitInterval(derivative, 3);
    const auto positive = [&](double fraction) {
        return sirius::core::detail::EvaluatePolynomial(polynomial, 4, fraction) >
               128.0 * std::numeric_limits<double>::epsilon() * scale;
    };
    if (!positive(0.0) || !positive(1.0)) return false;
    for (int i = 0; i < extrema.count; ++i)
        if (!positive(extrema.values[static_cast<unsigned>(i)])) return false;
    return true;
}

inline bool Derivative(double mass, double spin, const std::array<Constants, 2>& constants,
                       double fraction, const State& state, State& result) {
    for (unsigned column = 0; column < 2; ++column) {
        const auto& c = constants[column];
        // Carry the varying lower limit: u(s)=u0(1-s), including du0/dangle.
        const Jet u = c.inverse_radius * (1.0 - fraction);
        const Jet u2 = u * u;
        const Jet d = Jet(1.0) - u * (2.0 * mass) + u2 * (spin * spin);
        const Jet p = c.energy + (c.energy * (spin * spin) - c.angular_momentum * spin) * u2;
        const Jet h = p * p - d * u2 * c.carter_combination;
        if (!(h.real > 0.0) || !(d.real > 0.0)) return false;
        const Jet root = sqrt(h);
        const Jet drag =
            (c.energy * u * (2.0 * mass * spin) - c.angular_momentum * u2 * (spin * spin)) / d;
        const Jet twist = drag / root + Jet(spin) / d;
        const Vector n{state[column][0], state[column][1], state[column][2]};
        const Vector v{state[column][3], state[column][4], state[column][5]};
        const Vector cross_n = SpinCross(n), cross_v = SpinCross(v);
        const Jet v2 = Dot(v, v);
        for (unsigned i = 0; i < 3; ++i) {
            const Jet force = c.energy * c.energy * n[2] * (Jet(i == 2 ? 1.0 : 0.0) - n[2] * n[i]) *
                              (spin * spin);
            result[column][i] = c.inverse_radius * (v[i] / root + twist * cross_n[i]);
            result[column][i + 3] =
                c.inverse_radius * ((force - v2 * n[i]) / root + twist * cross_v[i]);
            if (!Finite(result[column][i]) || !Finite(result[column][i + 3])) return false;
        }
    }
    return true;
}
}  // namespace kerr_infinity_detail

// Inputs use the public ingoing Cartesian chart. Displacements are at the
// actual handoff event; coordinate variations are K=V-Gamma(k,X). Both contain
// the original launch seed. Spin is dimensional a (angular momentum / mass),
// not the dimensionless a/M. The metric/dg must be those of vacuum Kerr at x.
[[nodiscard]] inline std::expected<KerrInfinityResult, KerrInfinityFailure> TraceKerrInfinity(
    double mass, double spin, const Metric4d& metric,
    const Tensor<Dual<double>, 4, 4, 4>& metric_derivatives, const Vec4& position,
    const Vec4& tangent, const std::array<Vec4, 2>& endpoint_displacements,
    const std::array<Vec4, 2>& coordinate_variations, double angular_seed,
    const KerrInfinityConfig& config = {}) {
    using namespace kerr_infinity_detail;
    if (!std::isfinite(mass) || mass < 0.0 || !std::isfinite(spin) ||
        (mass > 0.0 && std::abs(spin) > mass) || !IsFinite(position) || !IsFinite(tangent) ||
        !std::isfinite(angular_seed) || !(angular_seed > 0.0) ||
        !std::isfinite(config.absolute_tolerance) || !(config.absolute_tolerance > 0.0) ||
        !std::isfinite(config.relative_tolerance) || config.relative_tolerance < 0.0 ||
        !std::isfinite(config.input_constraint_tolerance) ||
        config.input_constraint_tolerance < 0.0)
        return std::unexpected(KerrInfinityFailure::InvalidInput);
    State state{};
    std::array<Constants, 2> constants;
    std::array<Jet, 2> frequencies;
    for (unsigned column = 0; column < 2; ++column) {
        if (!IsFinite(endpoint_displacements[column]) || !IsFinite(coordinate_variations[column]) ||
            !InitialState(mass, spin, metric, metric_derivatives, position, tangent,
                          endpoint_displacements[column], coordinate_variations[column],
                          angular_seed, constants[column], state[column], frequencies[column],
                          config.input_constraint_tolerance))
            return std::unexpected(KerrInfinityFailure::InvalidInput);
    }
    if (!OutwardInterval(mass, spin, constants[0]))
        return std::unexpected(KerrInfinityFailure::RadialUnresolved);
    // Dormand-Prince 5(4), with a single accepted step shared by both columns.
    constexpr std::array<double, 7> nodes{0, 1.0 / 5, 3.0 / 10, 4.0 / 5, 8.0 / 9, 1, 1};
    constexpr double tableau[7][6] = {
        {},
        {1.0 / 5},
        {3.0 / 40, 9.0 / 40},
        {44.0 / 45, -56.0 / 15, 32.0 / 9},
        {19372.0 / 6561, -25360.0 / 2187, 64448.0 / 6561, -212.0 / 729},
        {9017.0 / 3168, -355.0 / 33, 46732.0 / 5247, 49.0 / 176, -5103.0 / 18656},
        {35.0 / 384, 0, 500.0 / 1113, 125.0 / 192, -2187.0 / 6784, 11.0 / 84}};
    constexpr std::array<double, 7> fourth{
        5179.0 / 57600, 0, 7571.0 / 16695, 393.0 / 640, -92097.0 / 339200, 187.0 / 2100, 1.0 / 40};
    double fraction = 0.0, step = 0.125;
    KerrInfinityResult result{};
    while (fraction < 1.0) {
        if (result.attempted_steps++ >= config.maximum_attempts)
            return std::unexpected(KerrInfinityFailure::WorkLimit);
        step = std::min(step, 1.0 - fraction);
        if (!(fraction + step > fraction)) return std::unexpected(KerrInfinityFailure::Arithmetic);
        std::array<State, 7> stages{};
        State fifth{};
        for (unsigned stage = 0; stage < 7; ++stage) {
            State trial = state;
            for (unsigned previous = 0; previous < stage; ++previous)
                for (unsigned column = 0; column < 2; ++column)
                    for (unsigned i = 0; i < 6; ++i)
                        trial[column][i] +=
                            stages[previous][column][i] * (step * tableau[stage][previous]);
            if (!Derivative(mass, spin, constants, fraction + nodes[stage] * step, trial,
                            stages[stage]))
                return std::unexpected(KerrInfinityFailure::Arithmetic);
            if (stage == 6) fifth = trial;
        }
        double error = 0.0;
        for (unsigned column = 0; column < 2; ++column) {
            for (unsigned i = 0; i < 6; ++i) {
                Jet lower = state[column][i];
                for (unsigned stage = 0; stage < 7; ++stage)
                    lower += stages[stage][column][i] * (step * fourth[stage]);
                if (!Finite(lower) || !Finite(fifth[column][i]))
                    return std::unexpected(KerrInfinityFailure::Arithmetic);
                const auto normalized = [&](double old, double high, double low) {
                    const double difference = std::abs(high - low);
                    const double scale =
                        config.absolute_tolerance +
                        config.relative_tolerance * std::max(std::abs(old), std::abs(high));
                    if (!std::isfinite(difference) || !std::isfinite(scale) || !(scale > 0.0))
                        return std::numeric_limits<double>::infinity();
                    const double ratio = difference / scale;
                    return std::isfinite(ratio) ? ratio : std::numeric_limits<double>::infinity();
                };
                error = std::max(
                    {error, normalized(state[column][i].real, fifth[column][i].real, lower.real),
                     normalized(state[column][i].dual, fifth[column][i].dual, lower.dual)});
            }
        }
        if (!std::isfinite(error)) return std::unexpected(KerrInfinityFailure::Arithmetic);
        if (error <= 1.0) {
            state = fifth;
            fraction += step;
            ++result.accepted_steps;
            result.maximum_local_error_ratio = std::max(result.maximum_local_error_ratio, error);
        }
        step *= error == 0.0 ? 5.0 : std::clamp(0.9 * std::pow(error, -0.2), 0.1, 5.0);
    }
    std::optional<CelestialTangentBasis<double>> basis;
    for (unsigned column = 0; column < 2; ++column) {
        Vector n{state[column][0], state[column][1], state[column][2]};
        const Jet norm = sqrt(Dot(n, n));
        if (!Finite(norm) || !(norm.real > 0.0))
            return std::unexpected(KerrInfinityFailure::Arithmetic);
        for (auto& value : n) value /= norm;
        if (column == 0) {
            for (unsigned i = 0; i < 3; ++i) result.map.direction[i] = n[i].real;
            basis = MakeCelestialTangentBasis(result.map.direction);
            if (!basis) return std::unexpected(KerrInfinityFailure::Arithmetic);
        }
        for (unsigned i = 0; i < 3; ++i) {
            result.map.jacobian[0][column] += basis->first[i] * n[i].dual;
            result.map.jacobian[1][column] += basis->second[i] * n[i].dual;
        }
        result.frequency_derivative[column] = frequencies[column].dual;
    }
    result.frequency = frequencies[0].real;
    result.map.determinant = result.map.jacobian[0][0] * result.map.jacobian[1][1] -
                             result.map.jacobian[0][1] * result.map.jacobian[1][0];
    if (!std::isfinite(result.map.determinant))
        return std::unexpected(KerrInfinityFailure::Arithmetic);
    return result;
}

}  // namespace sirius::core::relativity
