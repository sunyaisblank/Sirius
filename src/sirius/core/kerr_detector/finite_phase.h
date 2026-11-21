#pragma once

#include "sirius/core/kerr_detector/winding_tail.h"

namespace sirius::core::relativity::finite_tail_detail {

using PhasePair = std::array<Range, 2>;

template <typename Function>
PhasePair PhaseGaussPanel(Range left, Range right, const Function& integrand) {
    const Range width = right - left, middle = (left + right) * Range(.5), half = width * Range(.5);
    const Range split = Range(2.0) * sqrt(Range(6.0) / Range(5.0));
    const Range far = half * sqrt((Range(3.0) + split) / Range(7.0));
    const Range near = half * sqrt((Range(3.0) - split) / Range(7.0));
    const Range w0 = (Range(18.0) - sqrt(Range(30.0))) / Range(72.0);
    const Range w1 = (Range(18.0) + sqrt(Range(30.0))) / Range(72.0);

    const auto a = integrand(middle - far), b = integrand(middle + far);
    const auto c = integrand(middle - near), d = integrand(middle + near);
    using Series = TaylorSeries<8>;
    const auto all = integrand(Series::Variable(left.lo, right.hi));
    const Range factor = Square(Square(Square(width))) * width / Range(44100.0);
    PhasePair result;
    for (unsigned i = 0; i < 2; ++i) {
        const double error = (factor * Abs(all[i].coefficient[8])).hi;
        result[i] = width * ((a[i] + b[i]) * w0 + (c[i] + d[i]) * w1) + Range(-error, error);
    }
    return result;
}

struct FiniteRadialPhase {
    Range leading, drag_at_turn, time_residual, spin_residual;
};

inline Range LogBounds(Range value) {
    if (!value.Finite() || !(value.lo > 0.0)) return InvalidRange();
    return Range(LogPoint(value.lo).lo, LogPoint(value.hi).hi);
}

inline std::optional<FiniteRadialPhase> EvaluateFiniteRadialPhase(double mass, double spin,
                                                                  Range centre, Range scale,
                                                                  Range start) {
    using kerr_infinity_detail::EvaluateRadial;
    using kerr_infinity_detail::RadialCoefficients;
    if (!RoundingEnvironment()) return {};
    const auto roots = EncloseScaledTurningPair(mass, spin, centre, scale);
    const auto inv = InvertStationaryChart(mass, spin, {centre, Range(0.0) - Square(scale)});
    if (!roots || !inv) return {};
    const Range outer = centre - scale * roots->turn_offset;
    const Range inner = centre + scale * roots->inner_offset;
    const Range gd = scale * roots->gap_over_s;
    const Range gap(std::max(0.0, gd.lo), gd.hi);
    const auto c = RadialCoefficients(mass, spin, Range(1.0), inv->lambda, inv->normalized_carter);
    const std::array cofactor{Range(1.0) / (inner * outer), c[3] + (inner + outer) * c[4], c[4]};
    const Range gt = sqrt(EvaluateRadial(cofactor, outer));
    const Range a(spin), a2 = Square(a), two_m(2.0 * mass), d1 = a * two_m,
                         d2 = Range(0.0) - a2 * inv->lambda;
    const Range dt = outer * (d1 + outer * d2) / (Range(1.0) + outer * (a2 * outer - two_m));
    const Range leading = Range(1.0) / gt;
    if (!dt.Finite() || !leading.Finite()) return {};
    const auto divided = [&](auto q2) {
        using N = decltype(q2);
        const N u = N(outer) - q2, sum = u + N(outer);
        const N g = sqrt((N(cofactor[2]) * u + N(cofactor[1])) * u + N(cofactor[0]));
        const N d = N(1.0) + u * (N(a2) * u - N(two_m));
        return std::array{
            (N(cofactor[1]) + N(cofactor[2]) * sum) / (g * N(gt) * (g + N(gt))),
            (N(dt) * (N(a2) * sum - N(two_m)) - N(d1) - N(d2) * sum) / (g * d),
        };
    };
    Range time{}, azimuth{};
    for (Range end : {start, Range(0.0)}) {
        const Range span = outer - end, highest = sqrt(span);
        const Range log_edge =
            LogBounds(highest + sqrt(span + gap)) - Range(.5) * LogBounds(roots->gap_over_s);
        time = time + Range(2.0) * leading * log_edge;
        if (!highest.Finite() || !time.Finite()) return {};
        // The singular part was integrated exactly. Since q²/sqrt(q²+gap)
        // <= q, the remaining contribution below cut is enclosed by cut².
        constexpr double cut = 0x1p-7;
        const auto near = divided(Range(0.0, cut * cut));
        time = time + Range(0.0, cut * cut) * near[0];
        azimuth = azimuth + Range(0.0, cut * cut) * near[1];
        if (highest.lo < cut) return {};
        double left = cut;
        while (left < highest.lo) {
            const double base = std::scalbn(1.0, std::ilogb(left));
            const double stop = std::min(2.0 * base, highest.lo);
            const auto integrand = [&](auto x) {
                using N = decltype(x);
                const N x2 = x * x;
                const N ds = N(2.0 * base * base) * x2 / sqrt(x2 + N(gap / Range(base * base)));
                const auto f = divided(N(base * base) * x2);
                return std::array{ds * f[0], ds * f[1]};
            };
            constexpr unsigned subdivisions = 8;
            const Range first(left / base), last(stop / base);
            const Range stride = (last - first) / Range(subdivisions);
            for (unsigned i = 0; i < subdivisions; ++i) {
                const auto result = PhaseGaussPanel(first + Range(i) * stride,
                                                    first + Range(i + 1) * stride, integrand);
                time = time + result[0];
                azimuth = azimuth + result[1];
            }
            left = stop;
        }
        const Range far(highest.lo, highest.hi);
        const auto end_density = divided(Square(far));
        const Range missed = Range(0.0, (Range(highest.hi) - Range(highest.lo)).hi > 0.0
                                            ? (Range(highest.hi) - Range(highest.lo)).hi
                                            : 0.0);
        const Range factor = Range(2.0) * far * missed;
        time = time + factor * end_density[0];
        azimuth = azimuth + factor * end_density[1];
    }
    if (!time.Finite() || !azimuth.Finite()) return {};
    return FiniteRadialPhase{leading, dt, time, azimuth};
}

}  // namespace sirius::core::relativity::finite_tail_detail
