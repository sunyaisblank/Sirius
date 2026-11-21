#pragma once

#include "sirius/core/kerr_infinity.h"

namespace sirius::core::relativity::finite_tail_detail {

using FilmVector = std::array<Range, 2>;
using FilmMatrix = std::array<FilmVector, 2>;

// Interval first variation over a rectangle of continuous film coordinates.
// Values and partial derivatives start from the same represented pupil.
struct FilmJet {
    Range value;
    FilmVector derivative{};

    FilmJet() = default;
    explicit FilmJet(double constant) : value(constant) {}
    explicit FilmJet(Range domain) : value(domain) {}
    FilmJet(Range domain, FilmVector gradient) : value(domain), derivative(gradient) {}
};

inline Range Enclosure(const FilmJet& a) { return a.value; }

inline FilmJet operator+(const FilmJet& a, const FilmJet& b) {
    return {a.value + b.value,
            {a.derivative[0] + b.derivative[0], a.derivative[1] + b.derivative[1]}};
}
inline FilmJet operator-(const FilmJet& a, const FilmJet& b) {
    return {a.value - b.value,
            {a.derivative[0] - b.derivative[0], a.derivative[1] - b.derivative[1]}};
}
inline FilmJet operator*(const FilmJet& a, const FilmJet& b) {
    return {a.value * b.value,
            {a.derivative[0] * b.value + a.value * b.derivative[0],
             a.derivative[1] * b.value + a.value * b.derivative[1]}};
}
inline FilmJet operator/(const FilmJet& a, const FilmJet& b) {
    const Range result = a.value / b.value;
    return {result,
            {(a.derivative[0] - result * b.derivative[0]) / b.value,
             (a.derivative[1] - result * b.derivative[1]) / b.value}};
}
inline FilmJet sqrt(const FilmJet& a) {
    const Range result = sqrt(a.value), scale = Range(2.0) * result;
    return {result, {a.derivative[0] / scale, a.derivative[1] / scale}};
}
inline FilmJet Square(const FilmJet& a) {
    const Range scale = Range(2.0) * a.value;
    return {Square(a.value), {scale * a.derivative[0], scale * a.derivative[1]}};
}

struct StationaryChart {
    FilmVector coordinate;  // u_min, H(u_min)/E².
    FilmMatrix jacobian;
    Range lambda, eta;  // L/E, Q/E².
    Range inverse_start, initial_mu, energy, mu_rate, nx, ny;
};

inline std::optional<StationaryChart> MakeStationaryChart(
    const kerr_infinity_detail::DetectorRadialBounds<FilmJet>& state, double spin) {
    using kerr_infinity_detail::EvaluateRadial;
    const auto normalized = [&] {
        auto c = state.coefficients;
        for (auto& coefficient : c) coefficient = coefficient / state.coefficients[0];
        c[0] = FilmJet(1.0);
        return c;
    }();
    const auto& c = normalized;
    const Range a = Range(4.0) * c[4].value, b = Range(3.0) * c[3].value,
                d = Range(2.0) * c[2].value;
    if (!(d.hi < 0.0 && b.lo > 0.0 && state.inward)) return {};
    const Range u = Range(-2.0) * d / (b + sqrt(Square(b) - Range(4.0) * a * d));
    if (!(u.Finite() && u.lo > state.start.hi && u.hi < state.horizon.lo)) return {};
    const std::array slope{FilmJet(0.0), FilmJet(2.0) * c[2], FilmJet(3.0) * c[3],
                           FilmJet(4.0) * c[4]};
    const Range curvature =
        Range(2.0) * c[2].value + u * (Range(6.0) * c[3].value + u * Range(12.0) * c[4].value);
    if (!(curvature.lo > 0.0)) return {};
    // H_u is zero at the stationary radius. The partial of H
    // gives the derivative of the minimum without a radial-root term.
    const FilmJet h = EvaluateRadial(c, FilmJet(u));
    const FilmJet slope_u = EvaluateRadial(slope, FilmJet(u));
    const FilmVector stationary_gradient{(Range(0.0) - slope_u.derivative[0]) / curvature,
                                         (Range(0.0) - slope_u.derivative[1]) / curvature};
    const Range positive_energy = Range(0.0) - state.energy.value;
    const Range lambda = (Range(0.0) - state.angular.value) / positive_energy;
    const Range eta =
        (state.carter.value - Square(state.angular.value - Range(spin) * state.energy.value)) /
        Square(positive_energy);
    if (!h.value.Finite() || !lambda.Finite() || !eta.Finite() ||
        !std::ranges::all_of(h.derivative, &Range::Finite) ||
        !std::ranges::all_of(stationary_gradient, &Range::Finite))
        return {};
    return StationaryChart{{u, h.value},
                           {stationary_gradient, h.derivative},
                           lambda,
                           eta,
                           state.start,
                           state.mu.value,
                           state.energy.value,
                           state.mu_rate.value,
                           state.nx.value,
                           state.ny.value};
}

template <typename Scalar>
struct InvertedConstants {
    Scalar lambda, eta, normalized_carter;
};

template <typename Scalar>
inline std::optional<InvertedConstants<Scalar>> InvertStationaryScalars(
    double mass, double spin, const std::array<Scalar, 2>& target) {
    const auto& [u, minimum] = target;
    const Scalar a(spin), a2 = Square(a), twice_mass(2.0 * mass);
    if (!(Enclosure(a).lo > 0.0 && Enclosure(u).lo > 0.0) || !Enclosure(minimum).Finite())
        return {};
    const Scalar u2 = Square(u), delta = Scalar(1.0) + u * (a2 * u - twice_mass);
    const Scalar k = Scalar(2.0) + u * (Scalar(2.0) * a2 * u - twice_mass) / delta;
    const Scalar coefficient = Scalar(4.0) - k;
    // Write P=1+(a²-aλ)u². H=P²-Du²C and H_u=0 give
    // (4-k)P²-4P+kH=0. The branch continuous from P=4/(4-k)
    // gives the exterior spherical photon constants at H=0.
    const Scalar p = (Scalar(4.0) + sqrt(Scalar(16.0) - Scalar(4.0) * coefficient * k * minimum)) /
                     (Scalar(2.0) * coefficient);
    const Scalar normalized_carter = (Square(p) - minimum) / (delta * u2);
    const Scalar lambda = (a2 - (p - Scalar(1.0)) / u2) / a;
    const Scalar eta = normalized_carter - Square(lambda - a);
    if (!(Enclosure(p).lo > 0.0 && Enclosure(eta).lo > 0.0)) return {};
    return InvertedConstants<Scalar>{lambda, eta, normalized_carter};
}

inline std::optional<InvertedConstants<Range>> InvertStationaryChart(double mass, double spin,
                                                                     const FilmVector& target) {
    return InvertStationaryScalars(mass, spin, target);
}

struct ScaledTurningPair {
    Range s, turn_offset, inner_offset, gap_over_s;
    double contraction;
};

inline std::optional<ScaledTurningPair> EncloseScaledTurningPair(double mass, double spin, Range u,
                                                                 Range s) {
    using kerr_infinity_detail::RadialCoefficients;
    if (!(s.lo >= 0.0) || !s.Finite()) return {};
    const Range s2 = Square(s);
    const auto constants = InvertStationaryChart(mass, spin, {u, Range(0.0) - s2});
    if (!constants) return {};
    const auto c =
        RadialCoefficients(mass, spin, Range(1.0), constants->lambda, constants->normalized_carter);
    const Range k = c[2] + u * (Range(3.0) * c[3] + u * Range(6.0) * c[4]);
    const Range l = c[3] + u * Range(4.0) * c[4];
    const Range base = Range(1.0) / sqrt(k);
    if (!base.Finite()) return {};

    // Write u_turn = u_min - s*t_minus, u_inner = u_min + s*t_plus,
    // and H(u_min) = -s². The exact translated polynomial has no
    // linear term. The two rescaled root equations remain regular at s=0.
    double maximum_contraction = 0.0;
    const Range broad = Widen(0.5 * base.lo, 2.0 * base.hi);
    const auto root = [&](double sign) -> std::optional<Range> {
        const Range linear = Range(sign) * l * s, quadratic = c[4] * s2;
        const auto factor = [&](Range t) { return k + t * (linear + quadratic * t); };
        const Range enclosure = factor(broad);
        const Range speed = sqrt(enclosure);
        const Range gradient =
            (linear + Range(2.0) * quadratic * broad) / (Range(2.0) * enclosure * speed);
        const double contraction = Magnitude(gradient);
        if (!gradient.Finite() || !(contraction < 1.0)) return {};
        maximum_contraction = std::max(maximum_contraction, contraction);
        Range result = Range(1.0) / speed;
        if (!(result.Finite() && broad.lo < result.lo && result.hi < broad.hi)) return {};
        for (unsigned i = 0; i < 5; ++i) {
            const Range next = Range(1.0) / sqrt(factor(result));
            if (!next.Finite()) return {};
            result = Range(std::max(result.lo, next.lo), std::min(result.hi, next.hi));
        }
        return result.Finite() ? std::optional(result) : std::nullopt;
    };
    const auto turn = root(-1.0), inner = root(1.0);
    if (!turn || !inner) return {};
    return ScaledTurningPair{s, *turn, *inner, *turn + *inner, maximum_contraction};
}

struct ChartProof {
    FilmMatrix preconditioner;
    FilmMatrix residual;
    double contraction;
};

inline std::optional<ChartProof> ProveInjectiveChart(const FilmMatrix& jacobian) {
    const auto mid = [](Range x) { return std::midpoint(x.lo, x.hi); };
    const auto& [first, second] = jacobian;
    const double a = mid(first[0]), b = mid(first[1]), c = mid(second[0]), d = mid(second[1]);
    const double determinant = a * d - b * c;
    if (!std::isfinite(determinant) || determinant == 0.0) return {};
    const FilmMatrix inverse{FilmVector{Range(d / determinant), Range(-b / determinant)},
                             FilmVector{Range(-c / determinant), Range(a / determinant)}};
    ChartProof proof{inverse, {}, 0.0};
    for (unsigned i = 0; i < 2; ++i) {
        for (unsigned j = 0; j < 2; ++j) {
            proof.residual[i][j] = Range(i == j ? 1.0 : 0.0) - (inverse[i][0] * jacobian[0][j] +
                                                                inverse[i][1] * jacobian[1][j]);
            if (!proof.residual[i][j].Finite()) return {};
        }
        const Range row = Abs(proof.residual[i][0]) + Abs(proof.residual[i][1]);
        proof.contraction = std::max(proof.contraction, row.hi);
    }
    if (!(proof.contraction < 1.0)) return {};
    return proof;
}

}  // namespace sirius::core::relativity::finite_tail_detail
