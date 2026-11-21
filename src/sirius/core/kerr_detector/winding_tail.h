#pragma once

#include "sirius/core/kerr_detector/radial_remainder.h"

namespace sirius::core::relativity::finite_tail_detail {

inline Range SinePoint(double angle) {
    using finite_tail_detail::Square;
    const Range x(angle), x2 = Square(x);
    Range term = x, sum = x;
    for (unsigned n = 1; n <= 12; ++n) {
        term = term * x2 / Range((2 * n) * (2 * n + 1));
        sum = n % 2 != 0 ? sum - term : sum + term;
    }
    // For 0 <= angle <= pi, the Taylor terms decrease after the
    // first term. Degree 27 gives the outward omitted-term bound.
    const double remainder = (term * x2 / Range(26 * 27)).hi;
    return Range((sum - Range(remainder)).lo, sum.hi);
}

inline constexpr Range kPiBounds{0x1.921fb54442d18p+1, 0x1.921fb54442d19p+1};

inline Range SineLobe(Range angle) {
    const double lower = std::max(0.0, angle.lo);
    const double upper = std::min(kPiBounds.hi, angle.hi);
    if (!(lower <= upper)) return finite_tail_detail::InvalidRange();
    const Range first = SinePoint(lower), last = SinePoint(upper);
    const bool includes_peak = lower <= kPiBounds.hi * 0.5 && upper >= kPiBounds.lo * 0.5;
    return Range(std::min(first.lo, last.lo),
                 includes_peak ? 1.0 : std::min(1.0, std::max(first.hi, last.hi)));
}

struct WindingTail {
    unsigned first_winding;
    Range half_azimuth_radial;
    double azimuth_slope_per_winding_lower, azimuth_slope_per_winding_upper;
    double source_polar_speed_lower, wrap_count_slope;
    double inverse_determinant_coefficient, tail_prefactor;
    int tail_exponent;
};

inline Range IntegerPower(Range factor, unsigned exponent) {
    Range result(1.0);
    while (exponent != 0) {
        if (exponent % 2 != 0) result = result * factor;
        exponent /= 2;
        if (exponent != 0) factor = Square(factor);
    }
    return result;
}

inline std::optional<WindingTail> EncloseWindingTail(
    double mass, double spin, Range centre, Range scale, Range launch_mu,
    double catalogue_clearance, Range film_determinant, const RegularizedRadialPhase& radial,
    const LatitudeRootProof& latitude, unsigned first_requested) {
    const FilmJet u(centre, {Range(1.0), Range(0.0)}), s(scale, {Range(0.0), Range(1.0)});
    const auto constants =
        InvertStationaryScalars(mass, spin, std::array{u, FilmJet(0.0) - Square(s)});
    if (!constants || !(catalogue_clearance > 0.0 && film_determinant.lo > 0.0)) return {};
    const FilmJet a2 = Square(FilmJet(spin)), &eta = constants->eta;
    const FilmJet p = a2 - eta - Square(constants->lambda);
    const FilmJet omega2 = (sqrt(Square(p) + FilmJet(4.0) * a2 * eta) - p) * FilmJet(0.5);
    const FilmJet amplitude2 = eta / omega2;
    const FilmJet drag(radial.limiting_drag, radial.drag_gradient);
    const auto integrands = [&](Range sine2) {
        const FilmJet z2(sine2);
        const FilmJet t = FilmJet(1.0) / sqrt(omega2 + a2 * amplitude2 * z2);
        const FilmJet phi = (constants->lambda / (FilmJet(1.0) - amplitude2 * z2) + drag) * t;
        return std::array{t, phi};
    };
    constexpr unsigned panels = 256;
    const Range width = kPiBounds / Range(panels);
    FilmJet period{}, half_azimuth{};
    for (unsigned i = 0; i < panels; ++i) {
        const auto [t, phi] = integrands(Square(SineLobe(Range(i, i + 1) * width)));
        period = period + t * FilmJet(width);
        half_azimuth = half_azimuth + phi * FilmJet(width);
    }
    if (!(half_azimuth.derivative[0].lo > 0.0 && period.value.lo > 0.0)) return {};
    const auto [time_density, angle_density] = integrands({0.0, 1.0});
    const Range amplitude = sqrt(amplitude2.value);
    const Range initial_distance = amplitude2.value - Square(Abs(launch_mu));
    const Range final_distance = amplitude * Range(catalogue_clearance);
    if (!(initial_distance.lo > 0.0 && final_distance.lo > 0.0)) return {};
    const Range turn_gradient =
        (Abs(launch_mu) / sqrt(initial_distance) + amplitude / sqrt(final_distance)) /
        (Range(2.0) * amplitude2.value);
    const auto endpoint_bound = [&](const FilmJet& integrand, unsigned axis) {
        return Range(2.0) * kPiBounds * Abs(integrand.derivative[axis]) +
               Abs(integrand.value) * turn_gradient * Abs(amplitude2.derivative[axis]);
    };
    const Range period_radial = Abs(period.derivative[0]);
    const Range time_end_radial = endpoint_bound(time_density, 0);
    const Range az_end_radial = endpoint_bound(angle_density, 0);
    const Range az_end_scale = endpoint_bound(angle_density, 1);
    const Range az_half_scale = Abs(half_azimuth.derivative[1]);
    const Range spin_radial = Abs(radial.spin_remainder.radial);
    const Range log_slope(-latitude.slope_upper);
    if (!(log_slope.lo > 0.0 && radial.time_remainder.value.lo > 0.0)) return {};

    // At a latitude root, m*T <= tau <= (m+2)*T. The second
    // inequality bounds log(1/s), and the first gives s <= gain*base^m.
    const Range phase_coefficient = Range(2.0) * radial.leading_factor;
    const Range multiplicative(8.0), base(19.0 / 64.0);
    const Range log_ratio = Range(0.0) - LogPoint(base.lo);
    if (!((radial.time_remainder.value / phase_coefficient).hi < LogPoint(multiplicative.lo).lo &&
          (period.value / phase_coefficient).lo > log_ratio.hi))
        return {};
    const Range log_depth_linear = period.value / phase_coefficient;
    const Range log_depth_constant = Range(2.0) * log_depth_linear;

    // Bound |d log(s)/du| along the unique latitude root by B*m+C.
    const Range normalized_phase_gradient = Range(2.0) * Abs(radial.leading_gradient[0]);
    const Range bw = (normalized_phase_gradient * log_depth_linear + period_radial) / log_slope;
    const Range bc = (normalized_phase_gradient * log_depth_constant +
                      Abs(radial.time_remainder.radial) + time_end_radial) /
                     log_slope;
    const Range aw =
        az_half_scale + Range(radial.spin_derivative_log_coefficient) * log_depth_linear;
    const Range ac = az_end_scale + Range(radial.spin_derivative_constant) +
                     Range(radial.spin_derivative_log_coefficient) * log_depth_constant;

    const Range edge_time =
        Range(-2.0) * radial.leading_factor * LogPoint(scale.hi) + radial.time_remainder.value;
    unsigned first = std::max(8U, first_requested);
    Range error{};
    for (; first <= 1024; ++first) {
        const Range m(first);
        // For any m >= first, (B*m+C)(D*m+E) base^m decreases.
        // Each factor grows by at most (m+1)/m.
        const Range ratio = base * Square((m + Range(1.0)) / m);
        if (!(ratio.hi < 1.0 && (m * period.value).lo > edge_time.hi)) continue;
        const Range shrinking =
            multiplicative * IntegerPower(base, first) * (bw * m + bc) * (aw * m + ac);
        error = az_end_radial + spin_radial + shrinking;
        if ((m * (half_azimuth.derivative[0] - Range(1.0))).lo > error.hi) break;
    }
    if (first > 1024) return {};
    const Range slope_lower = half_azimuth.derivative[0] - error / Range(first);
    const Range slope_upper = half_azimuth.derivative[0] + error / Range(first);

    // At a latitude root, d(mu)/dlog(s) is the polar velocity times
    // the latitude timing residual derivative. The source sphere area
    // element is dmu*dphi. Chain the original standard-film Jacobian
    // through h=-s², then bound the Gaussian by one and its area
    // normalization 2*pi*(1-exp(-8)) from below by one.
    const Range polar_speed = sqrt(final_distance * omega2.value);
    const Range film_gain = Range(2.0) / (Range(polar_speed.lo) * log_slope *
                                          Range(slope_lower.lo) * Range(film_determinant.lo));
    const Range wraps =
        Range(slope_upper.hi) * (Range(centre.hi) - Range(centre.lo)) / (Range(2.0) * kPiBounds);
    // For each of the two polar legs the number of azimuth roots is
    // <= 1 + wraps*m. Each response <= film_gain * s² / m.
    const Range base2 = Square(base);
    // base⁴<1/128, so base^(2*m) <= 2^(-floor(7*m/2)).
    // Keep the binary exponent separate to avoid intermediate underflow.
    if (!(Square(base2).hi < 1.0 / 128.0)) return {};
    const Range prefactor = Range(2.0) * Square(multiplicative) * Range(film_gain.hi) *
                            (Range(1.0) / Range(first) + wraps) / (Range(1.0) - base2);
    if (!(prefactor.Finite() && prefactor.hi > 0.0 && polar_speed.lo > 0.0)) return {};
    return WindingTail{first,
                       half_azimuth.derivative[0],
                       slope_lower.lo,
                       slope_upper.hi,
                       polar_speed.lo,
                       wraps.hi,
                       film_gain.hi,
                       prefactor.hi,
                       -((7 * int(first)) / 2)};
}

}  // namespace sirius::core::relativity::finite_tail_detail
