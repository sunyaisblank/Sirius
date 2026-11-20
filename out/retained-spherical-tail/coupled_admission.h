#pragma once

#include "coupled_handoff.h"

namespace sirius::core::relativity {
namespace finite_tail_detail {

// Use the lower scale obtained by treating the represented midpoint as the
// larger magnitude. This bounds the original symmetric error for every value
// in the enclosure, including an endpoint on the opposite side of zero.
inline Range CoupledComponentError(double value, SphericalTailBound bound, Range lower_scale) {
    const Range enclosure = Enclosure(bound);
    if (!enclosure.Finite() || !(enclosure.lo <= value && value <= enclosure.hi))
        return InvalidRange();
    return Abs(Range(value) - enclosure) / lower_scale;
}
}  // namespace finite_tail_detail

// Conservative upper bound for Geodesic::CoupledStateError against the exact
// terminal event, with the live central RMS and per-column scales. The launch
// projection and the subsequent source maps have separate error obligations.
[[nodiscard]] inline double SphericalCoupledStateErrorBound(const SphericalCoupledTail& tail,
                                                            const IntegratorConfig& config,
                                                            const Rk45CoupledState& control) {
    using namespace finite_tail_detail;
    const auto positive = [](double value) { return std::isfinite(value) && value > 0.0; };
    if (!RoundingEnvironment() || !positive(config.abs_tolerance) ||
        !positive(config.rel_tolerance) || !positive(control.length_scale) ||
        !positive(control.frequency_scale) || !positive(control.tolerance))
        return kInfinity;
    const Range length(control.length_scale), frequency(control.frequency_scale);
    const Range tolerance(control.tolerance), absolute(config.abs_tolerance),
        relative(config.rel_tolerance);
    Range sum;
    for (unsigned i = 0; i < tail.event_bounds.position.size(); ++i) {
        const Range x_scale = absolute * length + relative * Abs(Range(tail.position(i)));
        const Range k_scale = absolute * frequency + relative * Abs(Range(tail.tangent(i)));
        const auto x =
            CoupledComponentError(tail.position(i), tail.event_bounds.position[i], x_scale);
        const auto k =
            CoupledComponentError(tail.tangent(i), tail.event_bounds.tangent[i], k_scale);
        sum = sum + Square(x) + Square(k);
    }
    if (!sum.Finite()) return kInfinity;
    const auto mean = sum / Range(8.0);
    if (!mean.Finite()) return kInfinity;
    double ratio = Widen(0.0, std::sqrt(mean.hi)).hi;
    for (unsigned i = 0; i < tail.variations.size(); ++i) {
        const auto& launch = control.variations[i];
        if (!positive(control.column_scale[i]) || !IsFinite(launch.displacement) ||
            !IsFinite(launch.derivative))
            return kInfinity;
        const Range column(control.column_scale[i]);
        const auto& bounds = tail.variation_bounds[i];
        for (unsigned j = 0; j < bounds.displacement.size(); ++j) {
            const double displacement = tail.variations[i].displacement(j);
            const double derivative = tail.variations[i].derivative(j);
            const Range x_scale = tolerance * (length * column + Abs(Range(displacement)));
            const Range v_scale = tolerance * (frequency * column + Abs(Range(derivative)));
            const auto x = CoupledComponentError(displacement, bounds.displacement[j], x_scale);
            const auto v = CoupledComponentError(derivative, bounds.derivative[j], v_scale);
            ratio = std::max({ratio, x.hi, v.hi});
        }
    }
    return ratio;
}
}  // namespace sirius::core::relativity
