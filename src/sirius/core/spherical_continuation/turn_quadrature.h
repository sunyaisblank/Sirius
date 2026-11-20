#pragma once

#include "radial_quadrature.h"

#include <optional>

namespace sirius::core::relativity::finite_tail_detail {

inline std::optional<Range> VacuumTurnRange(Range start, Range impact_squared, double mass,
                                            double source_support) {
    if (mass == 0.0) {
        const auto turn = sqrt(impact_squared);
        if (turn.Finite() && turn.lo > source_support && turn.hi < start.lo) return turn;
        return {};
    }
    const Range two_mass = Range(2.0) * Range(mass);
    const double photon_sphere = (Range(3.0) * Range(mass)).hi;
    double lower = std::max(source_support, std::nextafter(photon_sphere, kInfinity));
    double upper = start.lo;
    if (!(lower > photon_sphere && upper > lower)) return {};
    const auto potential = [&](double radius) {
        const Range r(radius);
        return Range(1.0) - (Range(1.0) - two_mass / r) * impact_squared / Square(r);
    };
    // The radial potential is strictly increasing outside 3M for b² > 0.
    if (!(potential(lower).hi < 0.0 && potential(upper).lo > 0.0)) return {};
    for (unsigned i = 0; i < std::numeric_limits<double>::digits; ++i) {
        const double middle = std::midpoint(lower, upper);
        if (!(lower < middle && middle < upper)) break;
        const Range value = potential(middle);
        if (value.hi < 0.0)
            lower = middle;
        else if (value.lo > 0.0)
            upper = middle;
        else {
            // Continue tightening both sides when the interval uncertainty at
            // the midpoint prevents a sign decision.
            const double left = std::midpoint(lower, middle);
            const double right = std::midpoint(middle, upper);
            const bool raise_lower = left > lower && potential(left).hi < 0.0;
            const bool lower_upper = right < upper && potential(right).lo > 0.0;
            if (!raise_lower && !lower_upper) break;
            if (raise_lower) lower = left;
            if (lower_upper) upper = right;
        }
    }
    return Range(lower, upper);
}

struct TurnQuadrature {
    using Bounds = Quadrature::Bounds;
    using Result = Quadrature::Result;

    double mass;
    double orientation;
    Range impact_squared;
    Range start_radius;
    double end_radius;
    Range turn;
    const SphericalTailControl& control;
    unsigned panels = 0;

    Result Run() {
        // Differentiate the simple root r_t³ − b² r_t + 2Mb² = 0.
        turn_derivative_ =
            (turn - Range(2.0) * Range(mass)) / (Range(3.0) * Square(turn) - impact_squared);
        if (!turn_derivative_.Finite()) return std::unexpected(FiniteTailFailure::Domain);
        return Integrate(0.0, 1.0, control.absolute_tolerance);
    }

  private:
    Range turn_derivative_;

    template <typename Number>
    std::array<Number, 6> Evaluate(Number parameter) const {
        using D = RadialJet<Number>;
        const D b2{Number(impact_squared), Number(1.0)};
        const D periapsis{Number(turn), Number(turn_derivative_)};
        const D u(parameter), u_squared(parameter * parameter);
        const D three_turn = D(3.0) * periapsis;
        const D root_factor = three_turn * periapsis - b2;
        const D two_mass = D(2.0) * D(mass);
        const auto leg = [&](Range end, double sign) {
            // In r = r_t + y², dr = 2y dy and the same y multiplies
            // the radial speed. Cancelling it leaves a regular kernel at
            // the turn, with y = u sqrt(end - r_t).
            const D span = D{Number(end)} - periapsis;
            const D y_squared = span * u_squared;
            const D radius = periapsis + y_squared;
            const D radius_squared = radius * radius;
            const D factor = root_factor + y_squared * (three_turn + y_squared);
            const D measure = D(2.0) * sqrt(span) / sqrt(factor / (radius_squared * radius));
            const D h = two_mass / radius;
            const D lapse = D(1.0) - h;
            const D time = (D(2.0 * sign) * span * u * h - measure) / lapse;
            return std::array{measure, measure / radius_squared, time};
        };
        const auto inward = leg(start_radius, -orientation);
        const auto outward = leg(Range(end_radius), orientation);
        const auto affine = inward[0] + outward[0];
        const auto angle = inward[1] + outward[1];
        const auto time = inward[2] + outward[2];
        return {affine.real, angle.real, time.real, affine.dual, angle.dual, time.dual};
    }

    Result Integrate(double start, double end, double abs_budget) {
        if (Cancelled(control)) return std::unexpected(FiniteTailFailure::Cancelled);
        if (panels >= control.maximum_panels) return std::unexpected(FiniteTailFailure::WorkLimit);
        ++panels;

        const double middle = std::midpoint(start, end);
        if (!(start < middle && middle < end))
            return std::unexpected(FiniteTailFailure::Arithmetic);
        const Range width = Range(end) - Range(start);
        const Range half = width * Range(0.5);
        const Range center = Range(start) + half;
        const Range radical = Range(2.0) * sqrt(Range(6.0) / Range(5.0));
        const Range far = half * sqrt((Range(3.0) + radical) / Range(7.0));
        const Range near = half * sqrt((Range(3.0) - radical) / Range(7.0));
        const auto far_left = Evaluate(center - far), far_right = Evaluate(center + far);
        const auto near_left = Evaluate(center - near), near_right = Evaluate(center + near);
        const Range weight_shift = sqrt(Range(30.0));
        const Range far_weight = (Range(18.0) - weight_shift) / Range(72.0);
        const Range near_weight = (Range(18.0) + weight_shift) / Range(72.0);

        using OrderEight = TaylorSeries<8>;
        const auto derivatives = Evaluate(OrderEight::Variable(start, end));
        const Range width2 = width * width;
        const Range width4 = width2 * width2;
        // Four-node Gauss-Legendre. With f^(8)/8! the remainder
        // coefficient is (4!)^4 / (9 (8!)²) = 1/44100.
        const Range remainder_factor = width4 * width4 * width / Range(44100.0);
        Bounds result;
        bool accepted = true;
        for (unsigned i = 0; i < result.size(); ++i) {
            const auto term = derivatives[i].coefficient[8];
            const double error = (remainder_factor * Range(Magnitude(term))).hi;
            const auto estimate = width * ((far_left[i] + far_right[i]) * far_weight +
                                           (near_left[i] + near_right[i]) * near_weight);
            result[i] = estimate + Range(-error, error);
            const double value = std::midpoint(result[i].lo, result[i].hi);
            accepted = accepted && result[i].Finite() &&
                       MidpointError(result[i]) <=
                           std::max(abs_budget, control.relative_tolerance * std::abs(value));
        }
        if (accepted) return result;
        const auto left = Integrate(start, middle, abs_budget * 0.5);
        if (!left) return left;
        const auto right = Integrate(middle, end, abs_budget * 0.5);
        if (!right) return right;
        for (unsigned i = 0; i < result.size(); ++i) {
            result[i] = (*left)[i] + (*right)[i];
            if (!result[i].Finite()) return std::unexpected(FiniteTailFailure::Arithmetic);
        }
        return result;
    }

  public:
    TurnQuadrature(double mass, double orientation, Range impact_squared, Range start_radius,
                   double end_radius, Range turn, const SphericalTailControl& control)
        : mass(mass),
          orientation(orientation),
          impact_squared(impact_squared),
          start_radius(start_radius),
          end_radius(end_radius),
          turn(turn),
          control(control) {}
};
}  // namespace sirius::core::relativity::finite_tail_detail
