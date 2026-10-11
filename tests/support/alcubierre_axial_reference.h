#pragma once

// Independent finite axial null reference for the declared constant-velocity
// Alcubierre tanh profile. No live metric, derivative, connection, launch,
// integrator, boundary interpolant or sky helper supplies an expected value.
// Clark, Hiscock & Larson, gr-qc/9907019, equations (1)--(3), (22)--(25):
// https://arxiv.org/html/gr-qc/9907019
//
// Put q=x-v*t, w=v*(f(|q|)-1), and s=+/-1 for dx/dt=v*f+s.
// The metric becomes -dt^2+(dq-w*dt)^2. Its conserved negative Killing contraction gives
// C=-g(k,partial_t+v*partial_x)=(1+s*w)*k^t and k^q=s*C (equivalently k^x-v*k^t is constant).
// A zero-beta Eulerian past ray with local sky direction -s has k^t_0=-1,
// hence C=-D(q0), D=1+s*w; it does NOT set f(q0)=1. Then
// dt/dq=s/D, lambda=(q-q0)/(s*C), k^x=(v*f+s)*C/D.
// The finite Eulerian sky is exactly (0,-s,0,0), not an infinity limit.
//
// This finite reference covers escaping axial branches whose D is globally
// positive and bounded away from zero, with any coordinate turning points
// strictly inside the selected sphere. Trapped/horizon/near-D-zero branches
// are excluded from THIS finite witness, not from the product specification.
// Coarse/fine agreement is a required empirical reference refinement check;
// it is not a rigorous quadrature error bound.

#include <algorithm>
#include <array>
#include <cmath>
#include <numbers>
#include <stdexcept>

namespace sirius::test::alcubierre_axial_reference {

using Scalar = long double;
using Four = std::array<Scalar, 4>;

struct Input {
    Scalar velocity, radius, sigma;
    Scalar initial_time, initial_x, escape_radius;
    int branch;  // s; the camera and terminal Eulerian sky direction is -s.
};

enum class Fate { Escape };

struct State {
    Fate fate = Fate::Escape;
    Four position{}, tangent{}, sky{}, initial_position{}, initial_tangent{};
    Scalar affine = 0;
    Scalar initial_q = 0, terminal_q = 0, killing_constant = 0;
    Scalar denominator_lower_bound = 0;
    unsigned coordinate_turning_points = 0;
};

inline Scalar LogCosh(Scalar value) {
    const Scalar magnitude = std::abs(value);
    return magnitude + std::log1p(std::exp(-2 * magnitude)) - std::numbers::ln2_v<Scalar>;
}

inline Scalar Shape(const Input& input, Scalar q) {
    // tanh(b+a)-tanh(b-a) divided by 2*tanh(a) equals
    // cosh(a)^2/[cosh(b+a)*cosh(b-a)]. Its logarithmic form avoids
    // cancellation of far-field tanhs and overflow of far-field coshs.
    const Scalar a = input.sigma * input.radius, b = input.sigma * std::abs(q);
    return std::exp(2 * LogCosh(a) - LogCosh(b + a) - LogCosh(b - a));
}

inline Scalar Denominator(const Input& input, Scalar q) {
    return 1 + input.branch * input.velocity * (Shape(input, q) - 1);
}

inline void Validate(const Input& input, unsigned panels_per_segment) {
    if (!std::isfinite(input.velocity) || std::abs(input.velocity) > 10 ||
        !std::isfinite(input.radius) || !(input.radius > 0) || input.radius > 1000 ||
        !std::isfinite(input.sigma) || !(input.sigma > 0) || input.sigma > 1000 ||
        !(input.sigma * input.radius >= 0.1L && input.sigma * input.radius <= 100) ||
        !std::isfinite(input.initial_time) || !std::isfinite(input.initial_x) ||
        !std::isfinite(input.escape_radius) || !(input.escape_radius > std::abs(input.initial_x)) ||
        (input.branch != -1 && input.branch != 1) || panels_per_segment == 0 ||
        panels_per_segment % 2 != 0 || panels_per_segment > 8192) {
        throw std::invalid_argument("finite axial reference input domain");
    }
    // f ranges between zero and one, so this is a global D bound. The selected
    // six cases have lower bound >=0.4. A zero/near-zero D requires a different
    // causal witness rather than repairing a denominator or skipping a region.
    if (!(std::min(Scalar(1), 1 - input.branch * input.velocity) >= 0.25L))
        throw std::invalid_argument("trapped or near-D-zero axial branch excluded from witness");
}

template <typename Function>
Scalar Simpson(Function function, Scalar lower, Scalar upper, unsigned panels) {
    if (lower == upper) return 0;
    const Scalar step = (upper - lower) / panels;
    Scalar sum = 0;
    for (unsigned i = 0; i <= panels; ++i) {
        const Scalar value = function(lower + i * step);
        if (!std::isfinite(value)) throw std::runtime_error("nonfinite axial quadrature field");
        const unsigned weight = (i == 0 || i == panels) ? 1 : (i % 2 == 0 ? 2 : 4);
        sum += weight * value;
    }
    const Scalar result = step * sum / 3;
    if (!std::isfinite(result)) throw std::runtime_error("nonfinite axial quadrature sum");
    return result;
}

inline Scalar TimeAt(const Input& input, Scalar q, unsigned panels_per_segment) {
    const Scalar q0 = input.initial_x - input.velocity * input.initial_time;
    Scalar lower = std::min(q0, q), upper = std::max(q0, q);
    const Scalar exterior_denominator = 1 - input.branch * input.velocity;
    // Split the entire integration interval at the walls and their tails.
    // There is no profile/tail truncation. Integrating only the departure from
    // the exactly flat comoving slope keeps long exterior intervals well conditioned.
    const Scalar tail = input.radius + 16 / input.sigma;
    const std::array<Scalar, 5> breaks{-tail, -input.radius, 0, input.radius, tail};
    const auto correction = [&](Scalar event_q) {
        const Scalar denominator = Denominator(input, event_q);
        if (!std::isfinite(denominator) || !(denominator >= 0.25L))
            throw std::runtime_error("axial reference denominator left declared branch");
        return 1 / denominator - 1 / exterior_denominator;
    };
    Scalar integral = 0;
    for (const Scalar split : breaks) {
        if (split > lower && split < upper) {
            integral += Simpson(correction, lower, split, panels_per_segment);
            lower = split;
        }
    }
    integral += Simpson(correction, lower, upper, panels_per_segment);
    if (q < q0) integral = -integral;
    return input.initial_time + input.branch * ((q - q0) / exterior_denominator + integral);
}

inline void RequireFinite(const State& state) {
    for (const auto* vector : {&state.position, &state.tangent, &state.sky, &state.initial_position,
                               &state.initial_tangent}) {
        for (const Scalar value : *vector) {
            if (!std::isfinite(value)) throw std::runtime_error("nonfinite axial reference vector");
        }
    }
    for (const Scalar value : {state.affine, state.initial_q, state.terminal_q,
                               state.killing_constant, state.denominator_lower_bound}) {
        if (!std::isfinite(value)) throw std::runtime_error("nonfinite axial reference scalar");
    }
}

inline State Orbit(const Input& input, unsigned panels_per_segment) {
    Validate(input, panels_per_segment);
    const Scalar q0 = input.initial_x - input.velocity * input.initial_time;
    const Scalar d0 = Denominator(input, q0), c = -d0;
    if (!std::isfinite(q0) || !std::isfinite(d0) || !(d0 >= 0.25L))
        throw std::invalid_argument("nonfinite or excluded axial launch branch");
    const Scalar travel_sign = -input.branch;
    const auto x_at_q = [&](Scalar q) {
        return q + input.velocity * TimeAt(input, q, panels_per_segment);
    };
    unsigned turning_points = 0;
    // k^x=0 iff f=-s/v. For the safe superluminal branch there are two
    // coordinate turns. Both must remain interior, otherwise this reference
    // does not pretend its eventual far root is the first outward event.
    if (input.branch * input.velocity < -1) {
        Scalar lower = 0, upper = input.radius + 16 / input.sigma;
        const Scalar turning_shape = -Scalar(input.branch) / input.velocity;
        if (!(Shape(input, upper) < turning_shape))
            throw std::runtime_error("axial coordinate-turn bracket unavailable");
        for (unsigned i = 0; i < 96; ++i) {
            const Scalar middle = (lower + upper) / 2;
            if (Shape(input, middle) > turning_shape)
                lower = middle;
            else
                upper = middle;
        }
        const Scalar turning_q = (lower + upper) / 2;
        for (const Scalar q : {-turning_q, turning_q}) {
            if (travel_sign * (q - q0) > 0) {
                const Scalar x = x_at_q(q);
                if (!std::isfinite(x) || !(std::abs(x) < input.escape_radius))
                    throw std::invalid_argument("earlier coordinate turn outside reference sphere");
                ++turning_points;
            }
        }
    }

    const Scalar target = travel_sign * input.escape_radius;
    const auto boundary_residual = [&](Scalar travel) {
        const Scalar q = q0 + travel_sign * travel;
        const Scalar result = travel_sign * (x_at_q(q) - target);
        if (!std::isfinite(result)) throw std::runtime_error("nonfinite axial event residual");
        return result;
    };
    Scalar lower = 0, upper = 2 * (input.escape_radius + std::abs(input.initial_x));
    bool bracketed = boundary_residual(upper) > 0;
    for (unsigned i = 0; i < 16 && !bracketed; ++i) {
        upper *= 2;
        bracketed = boundary_residual(upper) > 0;
    }
    if (!bracketed) throw std::runtime_error("finite axial escape bracket unavailable");
    for (unsigned i = 0; i < 96; ++i) {
        const Scalar middle = (lower + upper) / 2;
        if (boundary_residual(middle) > 0)
            upper = middle;
        else
            lower = middle;
    }

    const Scalar travel = (lower + upper) / 2, q = q0 + travel_sign * travel;
    const Scalar time = TimeAt(input, q, panels_per_segment);
    const Scalar f = Shape(input, q), kt = c / Denominator(input, q);
    State state;
    state.position = {time, q + input.velocity * time, 0, 0};
    state.tangent = {kt, (input.velocity * f + input.branch) * kt, 0, 0};
    state.sky = {0, travel_sign, 0, 0};
    state.initial_position = {input.initial_time, input.initial_x, 0, 0};
    state.initial_tangent = {-1, -input.velocity * Shape(input, q0) - input.branch, 0, 0};
    state.affine = travel / d0;
    state.initial_q = q0;
    state.terminal_q = q;
    state.killing_constant = c;
    state.denominator_lower_bound = std::min(Scalar(1), 1 - input.branch * input.velocity);
    state.coordinate_turning_points = turning_points;
    RequireFinite(state);
    const Scalar scale = std::max(input.radius, 1 / input.sigma);
    if (!(state.affine > 0) || !(state.tangent[0] < 0) ||
        !(std::abs(state.position[1] - target) / scale < 1.0e-12L) ||
        !(travel_sign * state.tangent[1] > 0)) {
        throw std::runtime_error("axial reference did not establish outward terminal event");
    }
    if (input.velocity == 0) {
        const Scalar exact_affine = input.escape_radius - travel_sign * input.initial_x;
        const Four exact_position{input.initial_time - exact_affine, target, 0, 0};
        const Four exact_tangent{-1, travel_sign, 0, 0};
        if (!(std::abs(state.affine - exact_affine) / scale < 1.0e-12L))
            throw std::runtime_error("flat axial affine self-check failed");
        for (unsigned axis = 0; axis < 4; ++axis) {
            if (!(std::abs(state.position[axis] - exact_position[axis]) / scale < 1.0e-12L) ||
                state.tangent[axis] != exact_tangent[axis] ||
                state.sky[axis] != (axis == 1 ? travel_sign : 0))
                throw std::runtime_error("flat axial state self-check failed");
        }
    }
    return state;
}

struct RefinedState {
    State coarse, fine;
    Scalar maximum_normalized_gap = 0;
    unsigned coarse_panels_per_segment = 0, fine_panels_per_segment = 0;
};

inline RefinedState Refine(const Input& input, unsigned coarse_panels_per_segment = 512,
                           unsigned fine_panels_per_segment = 1024) {
    if (!(fine_panels_per_segment > coarse_panels_per_segment))
        throw std::invalid_argument("axial reference refinement must increase panel count");
    RefinedState result{Orbit(input, coarse_panels_per_segment),
                        Orbit(input, fine_panels_per_segment), 0, coarse_panels_per_segment,
                        fine_panels_per_segment};
    // Check every field before reducing: std::max can hide a coarse NaN.
    RequireFinite(result.coarse);
    RequireFinite(result.fine);
    const Scalar scale = std::max(input.radius, 1 / input.sigma);
    const auto include = [&](Scalar gap) {
        if (!std::isfinite(gap)) throw std::runtime_error("nonfinite axial reference refinement");
        result.maximum_normalized_gap = std::max(result.maximum_normalized_gap, gap);
    };
    include(std::abs(result.coarse.affine - result.fine.affine) / scale);
    include(std::abs(result.coarse.terminal_q - result.fine.terminal_q) / scale);
    include(std::abs(result.coarse.killing_constant - result.fine.killing_constant));
    include(std::abs(result.coarse.denominator_lower_bound - result.fine.denominator_lower_bound));
    for (unsigned axis = 0; axis < 4; ++axis) {
        include(std::abs(result.coarse.position[axis] - result.fine.position[axis]) / scale);
        include(
            std::abs(result.coarse.initial_position[axis] - result.fine.initial_position[axis]) /
            scale);
        include(std::abs(result.coarse.tangent[axis] - result.fine.tangent[axis]));
        include(std::abs(result.coarse.initial_tangent[axis] - result.fine.initial_tangent[axis]));
        include(std::abs(result.coarse.sky[axis] - result.fine.sky[axis]));
    }
    if (result.coarse.coordinate_turning_points != result.fine.coordinate_turning_points ||
        !(result.maximum_normalized_gap < 1.0e-6L)) {
        throw std::runtime_error("axial reference refinement exceeds declared uncertainty");
    }
    return result;
}

}  // namespace sirius::test::alcubierre_axial_reference
