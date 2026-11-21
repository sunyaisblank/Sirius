#pragma once

#include "sirius/core/kerr_detector/stationary_chart.h"

namespace sirius::core::relativity::finite_tail_detail {

inline Range LogReduced(Range mantissa) {
    const auto z = (mantissa - Range(1.0)) / (mantissa + Range(1.0));
    const auto square = finite_tail_detail::Square(z);
    Range term = z, sum{};
    for (unsigned n = 1; n < 40; n += 2) {
        sum = sum + term / Range(n);
        term = term * square;
    }
    const double remainder =
        (Range(2.0) * finite_tail_detail::Abs(term) / (Range(41.0) * (Range(1.0) - square))).hi;
    return Range(2.0) * sum + Range(-remainder, remainder);
}

inline Range LogPoint(double value) {
    if (!(value > 0.0) || !std::isfinite(value)) return finite_tail_detail::InvalidRange();
    if (value == 1.0) return Range(0.0);
    int exponent;
    const double mantissa = std::frexp(value, &exponent) * 2;
    const static Range log_two = LogReduced(Range(2.0));
    return LogReduced(Range(mantissa)) + Range(exponent - 1) * log_two;
}

inline FilmJet Log(const FilmJet& a) {
    if (!a.value.Finite() || !(a.value.lo > 0.0)) return FilmJet(InvalidRange());
    return {Range(LogPoint(a.value.lo).lo, LogPoint(a.value.hi).hi),
            {a.derivative[0] / a.value, a.derivative[1] / a.value}};
}

struct RemainderDerivatives {
    Range value;
    Range radial;
    Range log_scale;
};

struct RegularizedRadialPhase {
    Range leading_factor;  // tau = -2 leading_factor log(s) + time_remainder.
    Range limiting_drag;   // phi_drag = limiting_drag * tau + spin_remainder.
    FilmVector leading_gradient, drag_gradient;
    RemainderDerivatives time_remainder, spin_remainder;
    Range time_log_derivative;  // s d(tau)/ds.
    double spin_derivative_constant, spin_derivative_log_coefficient;
};

inline std::optional<FilmJet> DifferentiateScaledTurningOffset(const FilmJet& k,
                                                               const FilmJet& linear,
                                                               const FilmJet& quadratic,
                                                               Range offset) {
    const FilmJet t(offset);
    const FilmJet polynomial = k + t * (linear + quadratic * t);
    const FilmJet fixed_point = FilmJet(1.0) / sqrt(polynomial);
    const Range map_slope = (Range(0.0) - linear.value - Range(2.0) * quadratic.value * offset) /
                            (Range(2.0) * polynomial.value * sqrt(polynomial.value));
    const Range denominator = Range(1.0) - map_slope;
    FilmJet result(offset);
    for (unsigned i = 0; i < 2; ++i) result.derivative[i] = fixed_point.derivative[i] / denominator;
    if (!std::ranges::all_of(result.derivative, &Range::Finite)) return {};
    return result;
}

inline Range ScatteringWeight(Range span, Range gap, Range x) {
    const auto endpoint = [](double width, double separation, double fraction) {
        if (fraction == 0.0) return Range(0.0);
        const Range w(width), t(fraction);
        return Range(2.0) * w * t / sqrt(Range(1.0) + Range(separation) / (w * Square(t)));
    };
    // The weight increases with leg span and x, and decreases with the gap.
    return Range(endpoint(span.lo, gap.hi, x.lo).lo, endpoint(span.hi, gap.lo, x.hi).hi);
}

inline Range GapFraction(Range gap, Range radial_distance) {
    const auto value = [](double delta, double q2, double at_zero) {
        if (q2 == 0.0) return Range(delta == 0.0 ? at_zero : 1.0);
        if (delta == 0.0) return Range(0.0);
        return Range(delta) / (Range(delta) + Range(q2));
    };
    return Range(std::max(0.0, value(gap.lo, radial_distance.hi, 0.0).lo),
                 std::min(1.0, value(gap.hi, radial_distance.lo, 1.0).hi));
}

inline void EncloseWeightedPhasePanel(RemainderDerivatives& result, const FilmJet& integrand,
                                      const std::array<Range, 2>& weight_rates, Range s,
                                      Range integrated_weight) {
    result.value = result.value + integrated_weight * integrand.value;
    result.radial = result.radial + integrated_weight * (integrand.derivative[0] +
                                                         weight_rates[0] * integrand.value);
    result.log_scale = result.log_scale + integrated_weight * (s * integrand.derivative[1] +
                                                               weight_rates[1] * integrand.value);
}

inline std::optional<RegularizedRadialPhase> EncloseRegularizedRadialPhase(double mass, double spin,
                                                                           Range centre,
                                                                           Range scale, Range start,
                                                                           unsigned panels) {
    using kerr_infinity_detail::EvaluateRadial;
    using kerr_infinity_detail::RadialCoefficients;
    if (panels == 0 || panels > 1024 || !RoundingEnvironment()) return {};
    const auto turning = EncloseScaledTurningPair(mass, spin, centre, scale);
    if (!turning || !start.Finite() || !(scale.hi > 0.0 && scale.hi < 0.25)) return {};
    const FilmJet u(centre, {Range(1.0), Range(0.0)});
    const FilmJet s(scale, {Range(0.0), Range(1.0)});
    const FilmJet s2 = Square(s);
    const auto constants = InvertStationaryScalars(mass, spin, std::array{u, FilmJet(0.0) - s2});
    if (!constants) return {};
    const auto c = RadialCoefficients(mass, spin, FilmJet(1.0), constants->lambda,
                                      constants->normalized_carter);
    const FilmJet k = c[2] + u * (FilmJet(3.0) * c[3] + u * FilmJet(6.0) * c[4]);
    const FilmJet l = c[3] + u * FilmJet(4.0) * c[4];
    const FilmJet curve = c[4] * s2;
    const auto minus =
        DifferentiateScaledTurningOffset(k, FilmJet(0.0) - l * s, curve, turning->turn_offset);
    const auto plus = DifferentiateScaledTurningOffset(k, l * s, curve, turning->inner_offset);
    if (!minus || !plus) return {};
    const FilmJet gap_coefficient = *minus + *plus;
    const FilmJet turn = u - s * *minus, inner = u + s * *plus;
    const FilmJet gap_s = s * gap_coefficient;
    const Range gap(std::max(0.0, gap_s.value.lo), gap_s.value.hi);
    const std::array quartic_factor{FilmJet(1.0) / (turn * inner), c[3] + (turn + inner) * c[4],
                                    c[4]};
    const FilmJet gt = sqrt(EvaluateRadial(quartic_factor, turn)), leading = FilmJet(1.0) / gt;
    const FilmJet a(spin), a2 = Square(a), two_mass(2.0 * mass);
    const auto d = [&](const FilmJet& v) { return FilmJet(1.0) + v * (a2 * v - two_mass); };
    const FilmJet d1 = a * two_mass, d2 = FilmJet(0.0) - a2 * constants->lambda;
    const FilmJet drag = turn * (d1 + turn * d2) / d(turn);
    if (!leading.value.Finite() || !drag.value.Finite()) return {};

    RemainderDerivatives time{}, spin_angle{};
    FilmJet exact_time{};
    Range spin_constant{}, spin_log{};
    for (const Range end : {start, Range(0.0)}) {
        // q^2=span*x^2, u=turn-q^2, x in [0,1].
        const FilmJet span = turn - FilmJet(end);
        const FilmJet length = sqrt(span);
        if (!length.value.Finite()) return {};
        exact_time = exact_time +
                     FilmJet(2.0) * leading *
                         (Log(length + sqrt(span + gap_s)) - FilmJet(0.5) * Log(gap_coefficient));

        const std::array gap_rate{
            gap_coefficient.derivative[0] / gap_coefficient.value,
            Range(1.0) + scale * gap_coefficient.derivative[1] / gap_coefficient.value};
        const Range plain_span_s = span.derivative[1] / span.value;
        const std::array span_rate{span.derivative[0] / span.value, scale * plain_span_s};
        double largest_drag_change = 0.0;
        for (unsigned i = 0; i < panels; ++i) {
            const Range left = Range(i) / Range(panels), right = Range(i + 1) / Range(panels);
            const Range x(left.lo, right.hi), x2 = Square(x);
            const FilmJet position = turn * FilmJet(Range(1.0) - x2) + FilmJet(end * x2);
            const FilmJet g = sqrt(EvaluateRadial(quartic_factor, position));
            if (!g.value.Finite()) return {};
            const FilmJet sum_radius = position + turn;
            // Divide both residual integrands by q² algebraically. This
            // preserves the zero of the numerator at the outer turn.
            const FilmJet time_change =
                (quartic_factor[1] + quartic_factor[2] * sum_radius) / (g * gt * (g + gt));
            const FilmJet drag_change =
                (drag * (a2 * sum_radius - two_mass) - d1 - d2 * sum_radius) / (g * d(position));
            const Range fraction = GapFraction(gap, span.value * x2);
            std::array<Range, 2> weight_rate;
            for (unsigned axis = 0; axis < 2; ++axis) {
                weight_rate[axis] =
                    span_rate[axis] + Range(0.5) * fraction * (span_rate[axis] - gap_rate[axis]);
            }
            const Range weight = ScatteringWeight(span.value, gap, x) / Range(panels);
            EncloseWeightedPhasePanel(time, time_change, weight_rate, scale, weight);
            EncloseWeightedPhasePanel(spin_angle, drag_change, weight_rate, scale, weight);
            spin_constant =
                spin_constant + weight * (Abs(drag_change.derivative[1]) +
                                          Range(1.5) * Abs(plain_span_s) * Abs(drag_change.value));
            largest_drag_change = std::max(largest_drag_change, Magnitude(drag_change.value));
        }
        // Integral w*r dx <= Delta*log(1+span/Delta). Factor out
        // Delta=s*g to bound the scale derivative by s*(C+D*|log s|).
        const Range log_factor =
            Range(0.5) * Abs(gap_rate[1]) * Range(largest_drag_change) * gap_coefficient.value;
        const Range log_span = Log(span + gap_s).value - Log(gap_coefficient).value;
        if (!log_span.Finite()) return {};
        spin_constant = spin_constant + log_factor * Range(std::max(0.0, log_span.hi));
        spin_log = spin_log + log_factor;
    }
    time.value = time.value + exact_time.value;
    time.radial = time.radial + exact_time.derivative[0];
    time.log_scale = time.log_scale + scale * exact_time.derivative[1];

    // s*log(s) extends to zero and decreases on [0,s_max<1/e].
    const Range weighted_log((Range(scale.hi) * LogPoint(scale.hi)).lo, 0.0);
    const Range transverse =
        Range(-2.0) * (leading.value + weighted_log * leading.derivative[1]) + time.log_scale;
    for (const Range check :
         {time.value, time.radial, time.log_scale, spin_angle.value, spin_angle.radial,
          spin_angle.log_scale, transverse, spin_constant, spin_log}) {
        if (!check.Finite()) return {};
    }
    return RegularizedRadialPhase{leading.value, drag.value, leading.derivative, drag.derivative,
                                  time,          spin_angle, transverse,         spin_constant.hi,
                                  spin_log.hi};
}

struct LatitudeRootProof {
    Range half_period;
    double period_h_bound, endpoint_h_bound, winding_s2_bound, slope_upper;
};

inline std::optional<LatitudeRootProof> EncloseLatitudeRootMonotonicity(
    double mass, double spin, Range centre, Range scale, Range launch_mu,
    double catalogue_clearance, const RegularizedRadialPhase& radial) {
    // Here the second jet axis is h=-s^2, removing the explicit s^2
    // from each polar derivative. The launch pupil is fixed.
    const Range s2 = Square(scale);
    const FilmJet h(Range(0.0) - s2, {Range(0.0), Range(1.0)});
    const auto invariants = InvertStationaryScalars(mass, spin, std::array{FilmJet(centre), h});
    if (!invariants) return {};
    const FilmJet a2 = Square(FilmJet(spin)), eta = invariants->eta;
    const FilmJet b = a2 - eta - Square(invariants->lambda);
    const FilmJet omega2 = FilmJet(0.5) * (sqrt(Square(b) + FilmJet(4.0) * a2 * eta) - b);
    const FilmJet amplitude2 = eta / omega2;
    const Range amplitude = sqrt(amplitude2.value);
    const FilmJet frequency2 = omega2 + a2 * amplitude2 * FilmJet(Range(0.0, 1.0));
    const FilmJet time_density = FilmJet(1.0) / sqrt(frequency2);
    constexpr Range pi{0x1.921fb54442d18p+1, 0x1.921fb54442d19p+1};
    const Range period = pi * time_density.value;
    const Range period_h = pi * Abs(time_density.derivative[1]);
    const Range zero_turn_gap = amplitude - Abs(launch_mu);
    const double clearance = std::min(zero_turn_gap.lo, catalogue_clearance);
    if (!(clearance > 0.0 && amplitude.lo > 0.0 && period.lo > 0.0)) return {};

    // At an admissible source (or the launch), |mu| lies at least
    // clearance below A. Differentiate asin(|mu|/A) with respect to h.
    // Its denominator is 2 A^2 sqrt(A^2-mu^2).
    const Range angle_h = amplitude * Abs(amplitude2.derivative[1]) /
                          (Range(2.0) * amplitude2.value * sqrt(amplitude * Range(clearance)));
    const Range partial_h = period_h + angle_h * time_density.value;

    // Any latitude root after m complete polar half-librations satisfies
    // m T <= tau. Multiplying by s^2 and bounding s^2 log(1/s) at
    // s_max encloses m s^2 for every possible winding at its roots.
    const Range max_s2(Square(Range(scale.hi)));
    const Range upper_at_scale =
        Range(-2.0) * radial.leading_factor * LogPoint(scale.hi) + radial.time_remainder.value;
    const Range winding_s2 = max_s2 * upper_at_scale / period;
    const Range correction = Range(2.0) * (winding_s2 * period_h + Range(2.0) * max_s2 * partial_h);
    const Range crossing_slope = radial.time_log_derivative + correction;
    if (!(winding_s2.Finite() && crossing_slope.Finite() && crossing_slope.hi < 0.0)) return {};
    return LatitudeRootProof{period, period_h.hi, partial_h.hi, winding_s2.hi, crossing_slope.hi};
}

}  // namespace sirius::core::relativity::finite_tail_detail
