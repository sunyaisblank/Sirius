#pragma once

#include "sirius/core/metrics/kerr_schild_family.h"
#include "sirius/core/source_sky_map.h"
#include "sirius/core/spherical_continuation/radial_quadrature.h"
#include "sirius/core/trace_boundary.h"
#include "sirius/core/twofold.h"

#include <algorithm>
#include <cfloat>
#include <expected>

namespace sirius::core::relativity {

// Vacuum Kerr continuation from an exterior event to radial infinity.
// These tolerances concern the tail only; the inner trace and the detector
// integration have separate errors. Callers must budget all three.
struct KerrInfinityConfig {
    double absolute_tolerance = 1e-10;
    double relative_tolerance = 1e-9;
    double input_constraint_tolerance = 1e-8;
    unsigned maximum_attempts = 4096;
    bool allow_radial_turn = false;
    // Keep any skipped inward segment outside an explicitly governed radius.
    double minimum_turn_radius = 0.0;
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
    bool inward = false;
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
    // Use the oblate radial branch; the Euclidean radial sign can disagree.
    const Jet radial_derivative =
        (radius * (x[0] * k[0] + x[1] * k[1]) + denominator * x[2] * k[2] / radius) /
        (discriminant * (scale * scale));
    if (!Finite(radial_derivative) || radial_derivative.real == 0.0) return false;
    constants.inward = radial_derivative.real < 0.0;
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

template <typename Scalar>
std::array<Scalar, 5> RadialCoefficients(double mass, double spin, Scalar energy,
                                         Scalar angular_momentum, Scalar carter) {
    const Scalar a(spin);
    const auto a2 = a * a;
    const auto z = a2 * energy - a * angular_momentum;
    return {energy * energy, Scalar{}, Scalar(2) * energy * z - carter,
            Scalar(2) * Scalar(mass) * carter, z * z - a2 * carter};
}

template <typename Scalar, std::size_t N>
Scalar EvaluateRadial(const std::array<Scalar, N>& coefficients, Scalar u) {
    auto value = coefficients.back();
    for (std::size_t i = N - 1; i > 0; --i) value = value * u + coefficients[i - 1];
    return value;
}

inline bool ProvePositiveRadialBernstein(
    std::array<finite_tail_detail::Range, 5> coefficients, unsigned& budget) {
    using finite_tail_detail::Range;
    if (budget == 0) return false;
    --budget;
    if (!std::ranges::all_of(coefficients, &Range::Finite)) return false;
    if (std::ranges::all_of(coefficients, [](Range c) { return c.lo > 0.0; })) return true;
    if (coefficients.front().hi <= 0.0 || coefficients.back().hi <= 0.0) return false;

    // Each half lies in the convex hull of its Bernstein coefficients.
    // De Casteljau subdivision at the exact parameter 1/2 tightens that hull
    // without losing the cancellation already represented by the polynomial.
    std::array<Range, 5> left, right;
    left.front() = coefficients.front();
    right.back() = coefficients.back();
    for (std::size_t remaining = coefficients.size() - 1; remaining != 0; --remaining) {
        for (std::size_t i = 0; i < remaining; ++i)
            coefficients[i] = (coefficients[i] + coefficients[i + 1]) * Range(0.5);
        left[coefficients.size() - remaining] = coefficients.front();
        right[remaining - 1] = coefficients[remaining - 1];
    }
    return ProvePositiveRadialBernstein(left, budget) &&
           ProvePositiveRadialBernstein(right, budget);
}

inline bool ProvePositiveRadialPotential(
    const std::array<finite_tail_detail::Range, 5>& coefficients,
    finite_tail_detail::Range interval, unsigned& budget) {
    using finite_tail_detail::Range;
    if (budget == 0 || !interval.Finite()) return false;
    if (EvaluateRadial(coefficients, interval).lo > 0.0) {
        --budget;
        return true;
    }

    // Compose H(lower + width*t) by Horner in the power basis, then
    // convert the resulting quartic on t in [0,1] to Bernstein form.
    const Range lower(interval.lo), width = Range(interval.hi) - lower;
    std::array<Range, 5> power{};
    power[0] = coefficients.back();
    for (std::size_t degree = 1; degree < power.size(); ++degree) {
        for (std::size_t i = degree; i != 0; --i)
            power[i] = lower * power[i] + width * power[i - 1];
        power[0] = lower * power[0] + coefficients[power.size() - 1 - degree];
    }
    const auto& [p0, p1, p2, p3, p4] = power;
    return ProvePositiveRadialBernstein(
        {p0, p0 + p1 * Range(0.25), p0 + p1 * Range(0.5) + p2 / Range(6.0),
         p0 + p1 * Range(0.75) + p2 * Range(0.5) + p3 * Range(0.25),
         p0 + p1 + p2 + p3 + p4},
        budget);
}

template <typename Scalar>
std::array<Scalar, 4> DeflateRadialTurn(const std::array<Scalar, 5>& coefficients, Scalar root) {
    // H(u) = (root - u) G(u). Horner deflation avoids dividing a near-zero
    // numerator by root-u when evaluating G at the turn.
    std::array<Scalar, 4> result;
    result[3] = Scalar{} - coefficients[4];
    for (std::size_t i = 3; i > 0; --i) result[i - 1] = result[i] * root - coefficients[i];
    return result;
}

// A short, nonoverlapping expansion for the sign of the represented radial
// potential. The words are ordered from least to most significant. Adding
// one word uses Shewchuk's grow-expansion with the existing error-free sums.
class RadialExpansion {
  public:
    explicit RadialExpansion(double value = 0.0) { Add(value); }

    int Sign() const {
        if (!valid_ || size_ == 0) return 0;
        return words_[size_ - 1] > 0.0 ? 1 : -1;
    }

    friend RadialExpansion operator+(RadialExpansion left, const RadialExpansion& right) {
        left.valid_ = left.valid_ && right.valid_;
        for (std::size_t i = 0; i < right.size_; ++i) left.Add(right.words_[i]);
        return left;
    }

    friend RadialExpansion operator-(RadialExpansion left, const RadialExpansion& right) {
        left.valid_ = left.valid_ && right.valid_;
        for (std::size_t i = 0; i < right.size_; ++i) left.Add(-right.words_[i]);
        return left;
    }

    friend RadialExpansion operator*(const RadialExpansion& left, const RadialExpansion& right) {
        RadialExpansion result;
        result.valid_ = left.valid_ && right.valid_;
        for (std::size_t i = 0; result.valid_ && i < left.size_; ++i) {
            for (std::size_t j = 0; result.valid_ && j < right.size_; ++j) {
                const auto product = Twofold::Product(left.words_[i], right.words_[j]);
                // Keep both the rounded product and its FMA residual away
                // from exponent limits. The residual then fits in binary64.
                const double magnitude = std::abs(product.hi);
                if (!(magnitude >= 0x1p-900 && magnitude <= 0x1p900)) {
                    result.valid_ = false;
                    break;
                }
                result.Add(product.lo);
                result.Add(product.hi);
            }
        }
        return result;
    }

  private:
    void Add(double value) {
        if (!valid_ || value == 0.0) return;
        if (!std::isfinite(value)) {
            valid_ = false;
            return;
        }
        std::size_t written = 0;
        for (std::size_t i = 0; i < size_; ++i) {
            const auto sum = Twofold::Sum(value, words_[i]);
            if (!std::isfinite(sum.hi) || !std::isfinite(sum.lo)) {
                valid_ = false;
                return;
            }
            if (sum.lo != 0.0) words_[written++] = sum.lo;
            value = sum.hi;
        }
        if (value != 0.0) {
            if (written == words_.size()) {
                valid_ = false;
                return;
            }
            words_[written++] = value;
        }
        size_ = written;
    }

    std::array<double, 32> words_{};
    std::size_t size_ = 0;
    bool valid_ = true;
};

inline int ExactRadialSign(double mass, double spin, const Constants& c, double inverse_radius) {
    // The enclosing solver has already checked rounding and gradual
    // underflow. Error-free transforms also require binary64 evaluation.
    if (FLT_EVAL_METHOD != 0) return 0;
    const RadialExpansion m(mass), a(spin), e(c.energy.real), l(c.angular_momentum.real);
    const RadialExpansion q(c.carter_combination.real), u(inverse_radius);
    const RadialExpansion one(1.0), two(2.0);
    const auto au = a * u;
    const auto p = e + (a * e - l) * (au * u);
    const auto d = one - (two * m) * u + au * au;
    return (p * p - d * (u * u) * q).Sign();
}

struct RadialTurn {
    Jet root, start, travel;
    std::array<Jet, 4> factor;
};

inline std::optional<RadialTurn> AdmitRadialTurn(double mass, double spin, const Constants& c,
                                                 double minimum_radius) {
    using finite_tail_detail::Range;
    if (!(mass > 0.0) || !finite_tail_detail::RoundingEnvironment()) return {};

    const auto interval_coefficients =
        RadialCoefficients(mass, spin, Range(c.energy.real), Range(c.angular_momentum.real),
                           Range(c.carter_combination.real));
    const auto coefficients =
        RadialCoefficients(mass, spin, c.energy, c.angular_momentum, c.carter_combination);
    double lower = c.inverse_radius.real;
    if (!(lower > 0.0) || !(EvaluateRadial(interval_coefficients, Range(lower)).lo > 0.0))
        return {};

    const Range M(mass), a(spin), absolute_spin(std::abs(spin));
    const Range outer_horizon = absolute_spin.IsPoint(mass)
                                    ? M
                                    : M + sqrt((M - absolute_spin) * (M + absolute_spin));
    const Range horizon_inverse = Range(1.0) / outer_horizon;
    if (!horizon_inverse.Finite()) return {};
    const double limit = minimum_radius > 0.0
                             ? std::min(horizon_inverse.lo,
                                        (Range(1.0) / Range(minimum_radius)).lo)
                             : horizon_inverse.lo;
    if (!(lower < limit)) return {};

    const auto sign = [&](double inverse_radius) {
        const Range bound = EvaluateRadial(interval_coefficients, Range(inverse_radius));
        if (bound.lo > 0.0) return 1;
        if (bound.hi < 0.0) return -1;
        return ExactRadialSign(mass, spin, c, inverse_radius);
    };
    const auto negative = [&](double candidate) {
        return lower < candidate && sign(candidate) < 0;
    };
    // Search only as far as the physical horizon or the configured
    // cutoff. A sign-definite turn inside the ergosphere is still regular.
    double upper = std::min((Range(0.25) / M).lo, limit);
    if (!negative(upper)) {
        upper = limit;
        if (!negative(upper)) {
            // The forbidden interval may be wholly inside the search range.
            // Roots of H'(u)/u suggest trial points. Interval signs and the
            // deflated factor below certify the bracket independently.
            const std::array derivative{2.0 * coefficients[2].real,
                                        3.0 * coefficients[3].real * limit,
                                        4.0 * coefficients[4].real * limit * limit};
            if (!std::ranges::all_of(derivative, [](double x) { return std::isfinite(x); }))
                return {};
            const auto extrema =
                sirius::core::detail::FindPolynomialRootsOnUnitInterval(derivative, 2);
            upper = lower;
            for (int i = 0; i < extrema.count; ++i) {
                const double candidate = extrema.values[static_cast<unsigned>(i)] * limit;
                if (negative(candidate)) {
                    upper = candidate;
                    break;
                }
            }
            if (!(lower < upper)) return {};
        }
    }

    unsigned sign_evaluations = 60;
    while (sign_evaluations != 0) {
        --sign_evaluations;
        const double mid = std::midpoint(lower, upper);
        if (!(lower < mid && mid < upper)) break;
        const int value = sign(mid);
        if (value > 0) {
            lower = mid;
        } else if (value < 0) {
            upper = mid;
        } else {
            // The sign at the centre can become ambiguous before the
            // bracket is narrow enough. Signs on either side can still
            // shrink it without inferring a sign for the ambiguous point.
            const auto tighten = [&](double proven, bool positive) {
                double unknown = mid;
                for (unsigned i = 0; i < 4 && sign_evaluations != 0; ++i) {
                    const double candidate = std::midpoint(proven, unknown);
                    if (candidate == proven || candidate == unknown) break;
                    --sign_evaluations;
                    const int candidate_sign = sign(candidate);
                    if (positive ? candidate_sign > 0 : candidate_sign < 0)
                        proven = candidate;
                    else
                        unknown = candidate;
                }
                return proven;
            };
            const double left = tighten(lower, true);
            const double right = tighten(upper, false);
            if (left == lower && right == upper) break;
            lower = left;
            upper = right;
        }
    }
    if (!(lower > c.inverse_radius.real) ||
        !(upper - lower <= 128 * std::numeric_limits<double>::epsilon() * upper))
        return {};

    // Both factors of the Kerr radial denominator are positive through
    // the bracket when the larger horizon radius is still below the turn.
    if (!((Range(1.0) - outer_horizon * Range(upper)).lo > 0.0)) return {};

    // The deflated factor must stay positive from infinity through the
    // entire bracket. This rules out an earlier or a multiple root.
    const auto factor = DeflateRadialTurn(interval_coefficients, Range(lower, upper));
    unsigned panels = 16;
    if (!ProvePositiveRadialPotential({factor[0], factor[1], factor[2], factor[3], Range(0.0)},
                                     Range(0.0, upper), panels))
        return {};

    const double root = std::midpoint(lower, upper);
    const Jet fixed_root(root);
    const auto slope = EvaluateRadial(DeflateRadialTurn(coefficients, fixed_root), fixed_root);
    // H'(root) = -G(root). Implicit differentiation carries the moving
    // turning point into both columns instead of holding it fixed.
    const Jet varied_root(root, EvaluateRadial(coefficients, fixed_root).dual / slope.real);
    const Jet start = -sqrt(varied_root - c.inverse_radius);
    const Jet travel = sqrt(varied_root) - start;
    if (!Finite(varied_root) || !Finite(start) || !Finite(travel)) return {};
    return RadialTurn{varied_root, start, travel, DeflateRadialTurn(coefficients, varied_root)};
}

inline bool TurnDerivative(double mass, double spin, const std::array<Constants, 2>& constants,
                           const std::array<RadialTurn, 2>& turns, double fraction,
                           const State& state, State& result) {
    for (unsigned column = 0; column < 2; ++column) {
        const auto& c = constants[column];
        const auto& turn = turns[column];
        // q increases through zero. du/dq = -2q, H = q^2 G, so the
        // Mino-time rate dτ/dq = 2/sqrt(G) on both radial branches.
        const Jet q = turn.start + turn.travel * fraction;
        const Jet u = fraction == 1.0 ? Jet{} : turn.root - q * q;
        const Jet u2 = u * u;
        const Jet d = Jet(1.0) - u * (2.0 * mass) + u2 * (spin * spin);
        const Jet factor = EvaluateRadial(turn.factor, u);
        if (!(factor.real > 0.0) || !(d.real > 0.0)) return false;
        const Jet advance = turn.travel * 2.0 / sqrt(factor);
        const Jet drag =
            (c.energy * u * (2.0 * mass * spin) - c.angular_momentum * u2 * (spin * spin)) / d;
        const Jet twist = advance * drag + turn.travel * q * (2.0 * spin) / d;
        const Vector n{state[column][0], state[column][1], state[column][2]};
        const Vector v{state[column][3], state[column][4], state[column][5]};
        const Vector cross_n = SpinCross(n), cross_v = SpinCross(v);
        const Jet v2 = Dot(v, v);
        for (unsigned i = 0; i < 3; ++i) {
            const Jet force = c.energy * c.energy * n[2] * (Jet(i == 2 ? 1.0 : 0.0) - n[2] * n[i]) *
                              (spin * spin);
            result[column][i] = advance * v[i] + twist * cross_n[i];
            result[column][i + 3] = advance * (force - v2 * n[i]) + twist * cross_v[i];
            if (!Finite(result[column][i]) || !Finite(result[column][i + 3])) return false;
        }
    }
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

// A dark vacuum detector sample needs only capture classification. Enclose
// the past null launch and prove the radial potential is strictly positive all
// the way to the outer horizon. No terminal state is computed here.
namespace kerr_infinity_detail {
struct DetectorRadialBounds {
    std::array<finite_tail_detail::Range, 5> coefficients;
    finite_tail_detail::Range start;
    finite_tail_detail::Range horizon;
    bool inward;
};

[[nodiscard]] inline std::optional<DetectorRadialBounds> EncloseDetectorRadialState(
    const KerrSchildParams& parameters, const std::array<finite_tail_detail::Range, 4>& position,
    const std::array<finite_tail_detail::Range, 4>& tangent) {
    using finite_tail_detail::Range;
    const auto [mass, spin, charge, cosmological_constant] = parameters;
    if (charge != 0.0 || cosmological_constant != 0.0 || !(mass > std::abs(spin)) ||
        !std::isfinite(mass) || !finite_tail_detail::RoundingEnvironment() ||
        !std::ranges::all_of(position, &Range::Finite) || !std::ranges::all_of(tangent, &Range::Finite))
        return std::nullopt;

    using finite_tail_detail::Square;
    const Range M(mass), a(spin);
    const Range x = position[1], y = position[2], z = position[3];
    const Range reduced = Square(x) + Square(y) + Square(z) - Square(a);
    const Range discriminant = sqrt(Square(reduced) + Range(4.0) * Square(a) * Square(z));
    const Range r = sqrt((reduced + discriminant) * Range(0.5));
    if (!(r.lo > (Range(2.0) * M).hi)) return std::nullopt;

    const Range kx = tangent[1], ky = tangent[2], kz = tangent[3];
    // The denominator of the oblate radius gradient is positive. Its
    // numerator fixes the radial branch before looking for any roots.
    const Range r2_plus_a2 = r * r + a * a;
    const Range radial = r * (x * kx + y * ky) + r2_plus_a2 * z * kz / r;
    if (!(radial.hi < 0.0 || radial.lo > 0.0)) return std::nullopt;

    const auto ell = metric_detail::KerrSchildNullCovector(x, y, z, r, spin);
    const auto H = metric_detail::KerrSchildFactor(r, z, parameters);
    const Range spatial = ell[1] * kx + ell[2] * ky + ell[3] * kz;
    const Range B = H * spatial;
    const Range C = Square(kx) + Square(ky) + Square(kz) + H * Square(spatial);
    const Range minus_gtt = Range(1.0) - H;
    const Range null_time = (B - sqrt(Square(B) + minus_gtt * C)) / minus_gtt;
    const Range time = tangent[0];
    if (!null_time.Finite() || !(null_time.hi < 0.0) || !(time.hi < 0.0)) return std::nullopt;
    // Point inputs must be close to null. The regional caller constructs
    // correlated null tangents from a certified observer frame and unit rays.
    const bool point = std::ranges::all_of(position, [](Range v) { return v.lo == v.hi; }) &&
                       std::ranges::all_of(tangent, [](Range v) { return v.lo == v.hi; });
    if (point && finite_tail_detail::Magnitude(null_time - time) >
                     -time.hi * KerrInfinityConfig{}.input_constraint_tolerance)
        return std::nullopt;
    const Range kt(std::min(time.lo, null_time.lo), std::max(time.hi, null_time.hi));
    const Range correction = H * (kt + spatial);
    const Range energy = kt - correction;
    if (!(energy.hi < 0.0)) return std::nullopt;
    const Range px = kx + correction * ell[1];
    const Range py = ky + correction * ell[2];
    const Range pz = kz + correction * ell[3];
    const Range angular = x * py - y * px;

    // Express the Carter combination as the squared norm of the angular
    // covector less the axial spin term. This form is regular on the axis.
    const std::array w{r * px + a * py, r * py - a * px, r * pz};
    const Range projected = ell[1] * w[0] + ell[2] * w[1] + ell[3] * w[2];
    const Range vx = w[0] - ell[1] * projected + a * energy * ell[2];
    const Range vy = w[1] - ell[2] * projected - a * energy * ell[1];
    const Range vz = w[2] - ell[3] * projected;
    const Range carter = Square(vx) + Square(vy) + Square(vz);
    const auto coefficients = RadialCoefficients(mass, spin, energy, angular, carter);
    const auto start = Range(1.0) / r;
    const auto horizon = Range(1.0) / (M + sqrt(M * M - a * a));
    if (!start.Finite() || !horizon.Finite() || !(start.lo > 0.0 && start.hi < horizon.lo))
        return std::nullopt;

    return DetectorRadialBounds{coefficients, start, horizon, radial.hi < 0.0};
}
}  // namespace kerr_infinity_detail

[[nodiscard]] inline bool CertifyKerrDetectorCapture(
    const KerrSchildParams& parameters, const std::array<finite_tail_detail::Range, 4>& position,
    const std::array<finite_tail_detail::Range, 4>& tangent, unsigned maximum_panels) {
    const auto bounds = kerr_infinity_detail::EncloseDetectorRadialState(parameters, position, tangent);
    if (!bounds || !bounds->inward) return false;
    unsigned budget = std::min(maximum_panels, 64U);
    return kerr_infinity_detail::ProvePositiveRadialPotential(
        bounds->coefficients, {bounds->start.lo, bounds->horizon.hi}, budget);
}

// The shared launch can publish an escape only if the whole launch enclosure
// is on the escape branch. The subsequent continuation certifies the
// represented turn and transports both angular columns.
[[nodiscard]] inline bool CertifyKerrDetectorEscape(
    const KerrSchildParams& parameters, const std::array<finite_tail_detail::Range, 4>& position,
    const std::array<finite_tail_detail::Range, 4>& tangent, double minimum_radius,
    unsigned maximum_panels) {
    using finite_tail_detail::Range;
    using namespace kerr_infinity_detail;
    const auto bounds = EncloseDetectorRadialState(parameters, position, tangent);
    if (!bounds) return false;
    if (!bounds->inward) {
        unsigned budget = std::min(maximum_panels, 64U);
        return ProvePositiveRadialPotential(bounds->coefficients, {0.0, bounds->start.hi}, budget);
    }

    const auto& c = bounds->coefficients;
    const double limit = minimum_radius > 0.0
                             ? std::min(bounds->horizon.lo, (Range(1.0) / Range(minimum_radius)).lo)
                             : bounds->horizon.lo;
    if (!(limit > bounds->start.hi) || !(EvaluateRadial(c, bounds->start).lo > 0.0)) return false;
    const auto centre = [](Range r) { return std::midpoint(r.lo, r.hi); };
    const std::array derivative{2.0 * centre(c[2]), 3.0 * centre(c[3]) * limit,
                                4.0 * centre(c[4]) * limit * limit};
    if (!std::ranges::all_of(derivative, [](double v) { return std::isfinite(v); })) return false;
    const auto extrema = sirius::core::detail::FindPolynomialRootsOnUnitInterval(derivative, 2);
    for (int i = 0; i < extrema.count; ++i) {
        const double probe = extrema.values[static_cast<std::size_t>(i)] * limit;
        if (probe > bounds->start.hi && EvaluateRadial(c, Range(probe)).hi < 0.0) return true;
    }
    return EvaluateRadial(c, Range(limit)).hi < 0.0;
}

[[nodiscard]] inline bool CertifyKerrDetectorCapture(const KerrSchildParams& parameters,
                                                     const Vec4& position, const Vec4& tangent,
                                                     unsigned maximum_panels) {
    using finite_tail_detail::Range;
    std::array<Range, 4> position_bounds, tangent_bounds;
    for (int i = 0; i < 4; ++i) {
        position_bounds[i] = Range(position(i));
        tangent_bounds[i] = Range(tangent(i));
    }
    return CertifyKerrDetectorCapture(parameters, position_bounds, tangent_bounds, maximum_panels);
}

// Inputs use the public ingoing Cartesian chart. Displacements are at the
// actual handoff event; coordinate variations are K=V-Gamma(k,X). Both contain
// the original launch seed. Spin is dimensional a (angular momentum / mass),
// not the dimensionless a/M. The metric/dg must be those of vacuum Kerr at x.
// The optional turn radius is written only with a complete map.
[[nodiscard]] inline std::expected<KerrInfinityResult, KerrInfinityFailure> TraceKerrInfinity(
    double mass, double spin, const Metric4d& metric,
    const Tensor<Dual<double>, 4, 4, 4>& metric_derivatives, const Vec4& position,
    const Vec4& tangent, const std::array<Vec4, 2>& endpoint_displacements,
    const std::array<Vec4, 2>& coordinate_variations, double angular_seed,
    const KerrInfinityConfig& config = {}, double* turn_radius = nullptr) {
    using namespace kerr_infinity_detail;
    if (!std::isfinite(mass) || mass < 0.0 || !std::isfinite(spin) ||
        (mass > 0.0 && std::abs(spin) > mass) || !IsFinite(position) || !IsFinite(tangent) ||
        !std::isfinite(angular_seed) || !(angular_seed > 0.0) ||
        !std::isfinite(config.absolute_tolerance) || !(config.absolute_tolerance > 0.0) ||
        !std::isfinite(config.relative_tolerance) || config.relative_tolerance < 0.0 ||
        !std::isfinite(config.input_constraint_tolerance) ||
        config.input_constraint_tolerance < 0.0 || !std::isfinite(config.minimum_turn_radius) ||
        config.minimum_turn_radius < 0.0)
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
    std::array<RadialTurn, 2> turns;
    const bool through_turn = constants[0].inward;
    if (through_turn) {
        if (!config.allow_radial_turn) return std::unexpected(KerrInfinityFailure::InvalidInput);
        for (unsigned column = 0; column < constants.size(); ++column) {
            const auto turn =
                AdmitRadialTurn(mass, spin, constants[column], config.minimum_turn_radius);
            if (!turn) return std::unexpected(KerrInfinityFailure::RadialUnresolved);
            turns[column] = *turn;
        }
    } else if (!OutwardInterval(mass, spin, constants[0])) {
        return std::unexpected(KerrInfinityFailure::RadialUnresolved);
    }
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
            const double stage_fraction = fraction + nodes[stage] * step;
            const bool derived = through_turn ? TurnDerivative(mass, spin, constants, turns,
                                                               stage_fraction, trial, stages[stage])
                                              : Derivative(mass, spin, constants, stage_fraction,
                                                           trial, stages[stage]);
            if (!derived) return std::unexpected(KerrInfinityFailure::Arithmetic);
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
    if (turn_radius) *turn_radius = through_turn ? 1.0 / turns[0].root.real : 0.0;
    return result;
}

}  // namespace sirius::core::relativity
