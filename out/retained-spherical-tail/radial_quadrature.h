#pragma once

#include <algorithm>
#include <array>
#include <cfenv>
#include <cmath>
#include <expected>
#include <limits>
#include <numeric>

namespace sirius::core::relativity {

enum class FiniteTailFailure { Input, Domain, Constraint, WorkLimit, Arithmetic };

struct SphericalTailControl {
    double absolute_tolerance = 1e-13;
    double relative_tolerance = 1e-11;
    double constraint_tolerance = 1e-8;
    unsigned maximum_panels = 4096;
};

struct RadialIntegralBound {
    double lower = 0;
    double upper = 0;
};

namespace finite_tail_detail {

inline bool RoundingEnvironment() {
    if (!std::numeric_limits<double>::is_iec559 || std::fegetround() != FE_TONEAREST) return false;
    const volatile double subnormal = std::numeric_limits<double>::denorm_min();
    const volatile double normal = std::numeric_limits<double>::min();
    return subnormal != 0.0 && normal * 0.5 != 0.0;
}

// Binary64 intervals enclose arithmetic on represented inputs in the default
// IEEE environment. Each floating operation is widened by one adjacent word.
struct Range {
    double lo = 0;
    double hi = 0;

    constexpr Range() = default;
    constexpr explicit Range(double value) : lo(value), hi(value) {}
    constexpr Range(double lower, double upper) : lo(lower), hi(upper) {}

    bool Finite() const { return std::isfinite(lo) && std::isfinite(hi) && lo <= hi; }
};

inline constexpr double kInfinity = std::numeric_limits<double>::infinity();
inline Range Widen(double low, double high) {
    return {std::nextafter(low, -kInfinity), std::nextafter(high, kInfinity)};
}
inline Range InvalidRange() { return {-kInfinity, kInfinity}; }
inline Range operator+(Range a, Range b) {
    if (!a.Finite() || !b.Finite()) return InvalidRange();
    return Widen(a.lo + b.lo, a.hi + b.hi);
}
inline Range operator-(Range a, Range b) {
    if (!a.Finite() || !b.Finite()) return InvalidRange();
    return Widen(a.lo - b.hi, a.hi - b.lo);
}
inline Range operator*(Range a, Range b) {
    if (!a.Finite() || !b.Finite()) return InvalidRange();
    const auto values = {a.lo * b.lo, a.lo * b.hi, a.hi * b.lo, a.hi * b.hi};
    const auto [lo, hi] = std::minmax(values);
    return Widen(lo, hi);
}
inline Range operator/(Range a, Range b) {
    if (!b.Finite() || b.lo <= 0.0) return InvalidRange();
    return a * Widen(1.0 / b.hi, 1.0 / b.lo);
}
inline Range sqrt(Range a) {
    if (!a.Finite() || a.lo <= 0.0) return InvalidRange();
    return Widen(std::sqrt(a.lo), std::sqrt(a.hi));
}

inline double MidpointError(Range enclosure) {
    const double estimate = std::midpoint(enclosure.lo, enclosure.hi);
    return std::max((Range(estimate) - Range(enclosure.lo)).hi,
                    (Range(enclosure.hi) - Range(estimate)).hi);
}

// Coefficients enclose f^(i)(r)/i! for every r in the panel. Only fourth
// derivatives are needed by the single-panel Simpson remainder.
struct Taylor {
    std::array<Range, 5> coefficient{};

    Taylor() = default;
    explicit Taylor(double value) { coefficient[0] = Range(value); }
    static Taylor Variable(double start, double end) {
        Taylor result;
        result.coefficient[0] = {start, end};
        result.coefficient[1] = Range(1.0);
        return result;
    }
};

inline Taylor operator+(const Taylor& a, const Taylor& b) {
    Taylor result;
    for (unsigned n = 0; n < a.coefficient.size(); ++n)
        result.coefficient[n] = a.coefficient[n] + b.coefficient[n];
    return result;
}

inline Taylor operator-(const Taylor& a, const Taylor& b) {
    Taylor result;
    for (unsigned n = 0; n < a.coefficient.size(); ++n)
        result.coefficient[n] = a.coefficient[n] - b.coefficient[n];
    return result;
}
inline Taylor operator*(const Taylor& a, const Taylor& b) {
    Taylor result;
    for (unsigned n = 0; n < a.coefficient.size(); ++n)
        for (unsigned j = 0; j <= n; ++j)
            result.coefficient[n] = result.coefficient[n] + a.coefficient[j] * b.coefficient[n - j];
    return result;
}
inline Taylor operator/(const Taylor& a, const Taylor& b) {
    Taylor result;
    for (unsigned n = 0; n < a.coefficient.size(); ++n) {
        Range rest = a.coefficient[n];
        for (unsigned j = 1; j <= n; ++j)
            rest = rest - b.coefficient[j] * result.coefficient[n - j];
        result.coefficient[n] = rest / b.coefficient[0];
    }
    return result;
}
inline Taylor sqrt(const Taylor& a) {
    Taylor result;
    result.coefficient[0] = sqrt(a.coefficient[0]);
    const Range denominator = Range(2.0) * result.coefficient[0];
    for (unsigned n = 1; n < a.coefficient.size(); ++n) {
        Range rest = a.coefficient[n];
        for (unsigned j = 1; j < n; ++j)
            rest = rest - result.coefficient[j] * result.coefficient[n - j];
        result.coefficient[n] = rest / denominator;
    }
    return result;
}

inline std::array<double, 3> RadialAtLaunch(double r, double mass, double orientation,
                                            double impact_squared) {
    const double h = 2.0 * mass / r;
    const double f = 1.0 - h;
    const double inverse_speed = 1.0 / std::sqrt(1.0 - f * (impact_squared / r / r));
    return {inverse_speed, inverse_speed / r / r, (orientation * h - inverse_speed) / f};
}

struct Quadrature {
    using Bounds = std::array<Range, 6>;
    using Result = std::expected<Bounds, FiniteTailFailure>;

    double mass;
    double orientation;
    double impact_squared;
    double start_radius;
    double end_radius;
    const SphericalTailControl& control;
    unsigned panels = 0;

    template <typename Number>
    std::array<Number, 6> Evaluate(Number parameter) const {
        using std::sqrt;
        // Subtract the exactly integrable flat part of affine length and
        // coordinate time before changing variables to u=1/r. Rationalize
        // (1/sqrt(D)-1)/u² to avoid subtractive cancellation near u=0.
        const Number inverse_end = Number(1.0) / Number(end_radius);
        const Number span = Number(1.0) / Number(start_radius) - inverse_end;
        const Number u = inverse_end + parameter * span;
        const Number u_squared = u * u;
        const Number two_mass = Number(2.0) * Number(mass);
        const Number f = Number(1.0) - two_mass * u;
        const Number speed = sqrt(Number(1.0) - f * (Number(impact_squared) * u_squared));
        const Number inverse_speed = Number(1.0) / speed;
        const Number residual = Number(impact_squared) / (speed * (Number(1.0) + speed));
        const Number half_inverse_cubed =
            Number(0.5) * inverse_speed * inverse_speed * inverse_speed;
        Number time = Number(0.0) - residual;
        if (orientation < 0.0 && mass > 0.0) time = time - Number(2.0) * two_mass / (u * f);
        // Frequency-scaled affine length, angle per unit impact, and
        // Kerr-Schild time, followed by their impact-squared partials.
        return {span * f * residual,
                span * inverse_speed,
                span * time,
                span * f * half_inverse_cubed,
                span * f * u_squared * half_inverse_cubed,
                Number(0.0) - span * half_inverse_cubed};
    }

    Result Run() {
        auto result =
            Integrate(0, 1, Evaluate(Range(0.0)), Evaluate(Range(1.0)), control.absolute_tolerance);
        if (!result) return result;
        const Range radial_span = Range(end_radius) - Range(start_radius);
        (*result)[0] = (*result)[0] + radial_span;
        (*result)[2] = (*result)[2] - radial_span;
        for (const auto enclosure : *result)
            if (!enclosure.Finite()) return std::unexpected(FiniteTailFailure::Arithmetic);
        return result;
    }

  private:
    Result Integrate(double start, double end, const Bounds& start_values, const Bounds& end_values,
                     double abs_budget) {
        if (panels >= control.maximum_panels) return std::unexpected(FiniteTailFailure::WorkLimit);
        ++panels;
        const Range width = Range(end) - Range(start);
        const double middle = std::midpoint(start, end);
        if (!(start < middle && middle < end))
            return std::unexpected(FiniteTailFailure::Arithmetic);
        const Range exact_middle = Range(start) + width * Range(0.5);
        const auto mid_values = Evaluate(exact_middle);
        const auto derivatives = Evaluate(Taylor::Variable(start, end));
        const Range width2 = width * width;
        // Single-panel Simpson error <= width^5 max|f^(4)| / 2880.
        // Taylor coefficient 4 is scaled by 4!, hence 2880/24 = 120.
        const Range width_factor = width2 * width2 * width / Range(120.0);
        Bounds result;
        bool accepted = true;
        for (unsigned i = 0; i < result.size(); ++i) {
            const Range fourth = derivatives[i].coefficient[4];
            if (!fourth.Finite()) {
                accepted = false;
                break;
            }
            const double magnitude = std::max(std::abs(fourth.lo), std::abs(fourth.hi));
            const double remainder = (width_factor * Range(magnitude)).hi;
            const Range estimate = (start_values[i] + Range(4.0) * mid_values[i] + end_values[i]) *
                                   (width / Range(6.0));
            result[i] = estimate + Range(-remainder, remainder);
            if (!result[i].Finite()) {
                accepted = false;
                break;
            }
            const double midpoint = std::midpoint(result[i].lo, result[i].hi);
            accepted = accepted &&
                       MidpointError(result[i]) <=
                           std::max(abs_budget, control.relative_tolerance * std::abs(midpoint));
        }
        if (accepted) return result;
        // The rounded midpoint defines adjacent child domains. Its point value
        // stands apart from the interval for the parent's exact midpoint.
        const auto child_values = Evaluate(Range(middle));
        const auto left = Integrate(start, middle, start_values, child_values, abs_budget * 0.5);
        if (!left) return std::unexpected(left.error());
        const auto right = Integrate(middle, end, child_values, end_values, abs_budget * 0.5);
        if (!right) return std::unexpected(right.error());
        for (unsigned i = 0; i < result.size(); ++i) {
            result[i] = (*left)[i] + (*right)[i];
            if (!result[i].Finite()) return std::unexpected(FiniteTailFailure::Arithmetic);
        }
        return result;
    }
};
}  // namespace finite_tail_detail
}  // namespace sirius::core::relativity
