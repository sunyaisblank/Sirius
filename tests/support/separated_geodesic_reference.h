#pragma once

// An independent, long-double Carter reference for complete vacuum Kerr rays.
// LaunchCameraRay supplies a measured initial event/tangent only. This file
// uses no product metric, connection, chart map, stepper, event interpolant,
// Jacobi state or infinity continuation to form its expected answers.
//
// With d lambda/d chi=Sigma, r''=R'(r)/2 and theta''=Theta'(theta)/2
// continue smoothly through radial and polar turns without a square-root
// sign switch. The first outward Euclidean sphere or the past horizon is
// localized by independently integrating partial RK4 steps. Captures use a
// regular outgoing chart; exterior events return to the public ingoing chart.
// Carter potentials/conventions: Gralla & Lupsasca, PRD 101, 044032 (2020),
// https://arxiv.org/abs/1910.12881, equations 3--9. Exact extremality is not
// inferred from the paper's elliptic formulas; the differential equations and
// regular-chart limits below are evaluated directly.

#include "sirius/core/camera_launch.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <optional>
#include <stdexcept>

namespace sirius::test::separated_reference {
using Scalar = long double;
using Three = std::array<Scalar, 3>;
using Four = std::array<Scalar, 4>;
using Matrix = std::array<std::array<Scalar, 2>, 2>;

inline Scalar Dot(const Three& a, const Three& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
inline Three Unit(Three a) {
    const Scalar norm = std::sqrt(Dot(a, a));
    if (!(norm > 0)) throw std::runtime_error("reference zero direction");
    for (auto& value : a) value /= norm;
    return a;
}
inline std::array<Three, 2> Basis(const Three& n) {
    unsigned least = 0;
    for (unsigned i = 1; i < 3; ++i)
        if (std::abs(n[i]) < std::abs(n[least])) least = i;
    Three first{};
    for (unsigned i = 0; i < 3; ++i) first[i] = (i == least ? 1 : 0) - n[i] * n[least];
    first = Unit(first);
    return {first, Three{n[1] * first[2] - n[2] * first[1], n[2] * first[0] - n[0] * first[2],
                         n[0] * first[1] - n[1] * first[0]}};
}

struct Geometry {
    Scalar radius, theta, azimuth, h;
    Four ell;
    std::array<std::array<Scalar, 4>, 4> metric{};
    std::array<Four, 3> spatial{};

    Scalar Inner(const Four& a, const Four& b) const {
        Scalar value = 0;
        for (unsigned mu = 0; mu < 4; ++mu)
            for (unsigned nu = 0; nu < 4; ++nu) value += metric[mu][nu] * a[mu] * b[nu];
        return value;
    }
    Three Sky(const Four& k) const {
        const Scalar frequency = -k[0] / std::sqrt(1 + h);
        if (!(frequency > 0)) throw std::runtime_error("reference cone orientation");
        Three n{};
        for (unsigned axis = 0; axis < 3; ++axis) n[axis] = Inner(k, spatial[axis]) / frequency;
        return Unit(n);
    }
    std::array<Four, 2> Screen(const Four& k) const {
        const auto basis = Basis(Sky(k));
        std::array<Four, 2> screen{};
        for (unsigned col = 0; col < 2; ++col)
            for (unsigned axis = 0; axis < 3; ++axis)
                for (unsigned mu = 0; mu < 4; ++mu)
                    screen[col][mu] += basis[col][axis] * spatial[axis][mu];
        return screen;
    }
};

inline Geometry At(const Four& x, Scalar mass, Scalar spin, bool outgoing = false) {
    const Scalar rho2 = x[1] * x[1] + x[2] * x[2] + x[3] * x[3];
    const Scalar reduced = rho2 - spin * spin;
    const Scalar r2 = (reduced + std::sqrt(reduced * reduced + 4 * spin * spin * x[3] * x[3])) / 2;
    Geometry g;
    g.radius = std::sqrt(r2);
    if (!(g.radius > 0)) throw std::runtime_error("reference singular radial chart");
    g.theta = std::acos(std::clamp(x[3] / g.radius, Scalar(-1), Scalar(1)));
    g.azimuth = std::atan2(x[2], x[1]);
    const Scalar handed = outgoing ? -spin : spin;
    g.ell = {outgoing ? -1.L : 1.L, (g.radius * x[1] + handed * x[2]) / (r2 + spin * spin),
             (g.radius * x[2] - handed * x[1]) / (r2 + spin * spin), x[3] / g.radius};
    const Scalar sigma = r2 + spin * spin * std::cos(g.theta) * std::cos(g.theta);
    g.h = 2 * mass * g.radius / sigma;
    for (unsigned mu = 0; mu < 4; ++mu)
        for (unsigned nu = 0; nu < 4; ++nu)
            g.metric[mu][nu] =
                (mu == nu ? (mu == 0 ? -1.L : 1.L) : 0.L) + g.h * g.ell[mu] * g.ell[nu];
    // Fixed Cartesian coordinate axes, Gram-Schmidt in the spatial metric.
    // Their temporal components vanish and u_mu=(-1/sqrt(1+h),0,0,0).
    for (unsigned axis = 0; axis < 3; ++axis) {
        g.spatial[axis][axis + 1] = 1;
        for (unsigned prior = 0; prior < axis; ++prior) {
            const Scalar projection = g.Inner(g.spatial[axis], g.spatial[prior]);
            for (unsigned mu = 0; mu < 4; ++mu)
                g.spatial[axis][mu] -= projection * g.spatial[prior][mu];
        }
        const Scalar length = std::sqrt(g.Inner(g.spatial[axis], g.spatial[axis]));
        for (auto& value : g.spatial[axis]) value /= length;
    }
    return g;
}

enum class Fate { Escape, Capture, Disk };
struct Result {
    Fate fate;
    Four x, k;
    Three finite_direction, infinity_direction;
    Scalar affine, energy, angular_momentum, carter;
    unsigned radial_turns = 0, polar_turns = 0, steps = 0;
    bool outgoing = false;
};

namespace detail {
// r, dr/dchi, theta, dtheta/dchi, outgoing Cartesian azimuth, t_out, lambda.
using State = std::array<Scalar, 7>;
struct Constants {
    Scalar mass, spin, energy, angular, carter;
};
inline Scalar Combination(const Constants& c) {
    return c.carter + (c.angular - c.spin * c.energy) * (c.angular - c.spin * c.energy);
}
inline Scalar Radius(const State& y, Scalar a) {
    return std::sqrt(y[0] * y[0] + a * a * std::sin(y[2]) * std::sin(y[2]));
}
// Time/Cartesian-azimuth shifts from ingoing to outgoing, with the same
// explicitly declared exterior gauge (zero shifts at r=2r+) as the public API.
inline std::array<Scalar, 2> Shift(Scalar r, const Constants& c) {
    const Scalar plus = c.mass + std::sqrt(c.mass * c.mass - c.spin * c.spin);
    const Scalar minus = c.mass - std::sqrt(c.mass * c.mass - c.spin * c.spin);
    const Scalar anchor = 2 * plus;
    Scalar time, inverse_delta;
    if (plus == minus) {
        inverse_delta = -1 / (r - plus) + 1 / (anchor - plus);
        time = 2 * c.mass * std::log((r - plus) / (anchor - plus)) +
               (c.mass * c.mass + c.spin * c.spin) * inverse_delta;
    } else {
        const Scalar lp = std::log((r - plus) / (anchor - plus));
        const Scalar lm = std::log((r - minus) / (anchor - minus));
        inverse_delta = (lp - lm) / (plus - minus);
        time = 2 * c.mass * (plus * lp - minus * lm) / (plus - minus);
    }
    const Scalar angle = c.spin == 0 ? 0
                                     : -2 * c.spin * inverse_delta +
                                           2 * (std::atan2(r, c.spin) - std::atan2(anchor, c.spin));
    return {-2 * time, angle};
}
inline State Derivative(const State& y, const Constants& c) {
    const Scalar r = y[0], v = y[1], theta = y[2], a = c.spin;
    const Scalar sine = std::sin(theta), cosine = std::cos(theta);
    if (!(std::abs(sine) > 1e-5L)) throw std::runtime_error("reference polar chart conditioning");
    const Scalar delta = r * r - 2 * c.mass * r + a * a;
    const Scalar P = c.energy * (r * r + a * a) - a * c.angular;
    const Scalar combination = Combination(c);
    // (P-v)/Delta = K/(P+v) on the inward past branch. It removes the
    // cancelling horizon pole instead of clamping Delta or either coordinate.
    const Scalar regular = v < 0 && P < 0 ? combination / (P + v) : (P - v) / delta;
    const Scalar azimuth =
        c.angular / (sine * sine) - a * c.energy + a * regular + a * v / (r * r + a * a);
    const Scalar time = P + 2 * c.mass * r * regular + a * (c.angular - a * c.energy * sine * sine);
    return {v,
            2 * c.energy * r * P - (r - c.mass) * combination,
            y[3],
            -a * a * c.energy * c.energy * sine * cosine +
                c.angular * c.angular * cosine / (sine * sine * sine),
            azimuth,
            time,
            r * r + a * a * cosine * cosine};
}
template <std::size_t N, typename Function>
inline std::array<Scalar, N> Rk4(const std::array<Scalar, N>& y, Scalar h, Function f) {
    const auto add = [](auto a, const auto& b, Scalar factor) {
        for (unsigned i = 0; i < N; ++i) a[i] += factor * b[i];
        return a;
    };
    const auto k1 = f(y), k2 = f(add(y, k1, h / 2));
    const auto k3 = f(add(y, k2, h / 2)), k4 = f(add(y, k3, h));
    auto result = y;
    for (unsigned i = 0; i < N; ++i) result[i] += h * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) / 6;
    return result;
}
inline Three Infinity(const State& y, const Constants& c, unsigned steps) {
    // Independent scalar Carter tail in u=1/r, not the product's sphere-vector
    // equations or adaptive method. The complete inner turn has already run.
    std::array<Scalar, 3> state{y[2], y[3], y[4] - std::atan2(c.spin, y[0])};
    const Scalar initial = 1 / y[0], h = -initial / steps;
    for (unsigned step = 0; step < steps; ++step) {
        const Scalar u = initial + h * step;
        const auto f = [&](Scalar at, const std::array<Scalar, 3>& q) {
            const Scalar sine = std::sin(q[0]), cosine = std::cos(q[0]);
            const Scalar D = 1 - 2 * c.mass * at + c.spin * c.spin * at * at;
            const Scalar P = c.energy + (c.spin * c.spin * c.energy - c.spin * c.angular) * at * at;
            const Scalar H = P * P - D * at * at * Combination(c);
            if (!(H > 0) || !(D > 0) || !(std::abs(sine) > 1e-5L))
                throw std::runtime_error("reference unresolved outward tail");
            const Scalar root = std::sqrt(H);
            return std::array<Scalar, 3>{
                -q[1] / root,
                (c.spin * c.spin * c.energy * c.energy * sine * cosine -
                 c.angular * c.angular * cosine / (sine * sine * sine)) /
                    root,
                -(c.angular / (sine * sine) - c.spin * c.energy + c.spin * P / D) / root -
                    c.spin / D};
        };
        const auto add = [](auto q, const auto& d, Scalar factor) {
            for (unsigned i = 0; i < 3; ++i) q[i] += factor * d[i];
            return q;
        };
        const auto k1 = f(u, state), k2 = f(u + h / 2, add(state, k1, h / 2));
        const auto k3 = f(u + h / 2, add(state, k2, h / 2)), k4 = f(u + h, add(state, k3, h));
        for (unsigned i = 0; i < 3; ++i)
            state[i] += h * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) / 6;
    }
    return {std::sin(state[0]) * std::cos(state[2]), std::sin(state[0]) * std::sin(state[2]),
            std::cos(state[0])};
}
}  // namespace detail

inline Result Trace(const core::CameraLaunch& launch, double mass, double spin, double outer_radius,
                    Scalar refinement = .002L, double disk_inner = 0, double disk_outer = 0) {
    const Scalar scale = mass > 0 ? mass : 1;
    Four x{}, k{};
    for (unsigned i = 0; i < 4; ++i) {
        x[i] = launch.position(static_cast<int>(i)) / scale;
        k[i] = launch.tangent(static_cast<int>(i));
    }
    const Scalar M = mass / scale, a = spin / scale, outer = outer_radius / scale;
    const auto initial = At(x, M, a);
    if (mass == 0 && spin == 0) {
        const Three position{x[1], x[2], x[3]}, velocity{k[1], k[2], k[3]};
        const Scalar vv = Dot(velocity, velocity), along = Dot(position, velocity);
        const Scalar affine =
            (-along + std::sqrt(along * along + vv * (outer * outer - Dot(position, position)))) /
            vv;
        for (unsigned i = 0; i < 4; ++i) x[i] += affine * k[i];
        const auto direction = Unit(velocity);
        for (auto& value : x) value *= scale;
        return {Fate::Escape, x, k, direction, direction, affine * scale, k[0], 0, 0};
    }
    Four p{};
    for (unsigned mu = 0; mu < 4; ++mu)
        for (unsigned nu = 0; nu < 4; ++nu) p[mu] += initial.metric[mu][nu] * k[nu];
    const Scalar r = initial.radius, sine = std::sin(initial.theta),
                 cosine = std::cos(initial.theta);
    const Scalar radial = (r * (x[1] * k[1] + x[2] * k[2]) + (r * r + a * a) * x[3] * k[3] / r) /
                          (r * r + a * a * cosine * cosine);
    const Scalar polar = (cosine * radial - k[3]) / (r * sine);
    const Scalar sigma = r * r + a * a * cosine * cosine;
    const Scalar energy = -p[0], angular = x[1] * p[2] - x[2] * p[1];
    const Scalar carter =
        sigma * sigma * polar * polar +
        cosine * cosine * (angular * angular / (sine * sine) - a * a * energy * energy);
    const detail::Constants c{M, a, energy, angular, carter};
    const auto shift = detail::Shift(r, c);
    detail::State y{r,
                    sigma * radial,
                    initial.theta,
                    sigma * polar,
                    initial.azimuth + shift[1],
                    x[0] + shift[0],
                    0};
    const Scalar horizon = M + std::sqrt(M * M - a * a);
    unsigned radial_turns = 0, polar_turns = 0, attempts = 0;
    Fate fate = Fate::Escape;
    for (; attempts < 200000; ++attempts) {
        const auto derivative = detail::Derivative(y, c);
        const Scalar rate =
            std::max({1.L, std::abs(y[1]) / y[0], std::abs(y[3]), std::abs(derivative[4])});
        const Scalar h = refinement / rate;
        const auto advance = [&](Scalar fraction) {
            return detail::Rk4(y, h * fraction,
                               [&](const auto& state) { return detail::Derivative(state, c); });
        };
        auto next = advance(1);
        bool captured = next[0] <= horizon;
        bool escaped = detail::Radius(next, a) >= outer && next[1] > 0;
        Scalar terminal_fraction = 1;
        if (captured || escaped) {
            Scalar lo = 0, hi = 1;
            for (unsigned iteration = 0; iteration < 52; ++iteration) {
                const Scalar mid = (lo + hi) / 2;
                const auto candidate = advance(mid);
                if (captured ? candidate[0] <= horizon : detail::Radius(candidate, a) >= outer)
                    hi = mid;
                else
                    lo = mid;
            }
            terminal_fraction = (lo + hi) / 2;
            next = advance(terminal_fraction);
            fate = captured ? Fate::Capture : Fate::Escape;
        }
        bool disk = false;
        const Scalar equator = std::acos(-1.L) / 2;
        if (disk_inner > 0 && disk_outer > disk_inner &&
            (y[2] - equator) * (next[2] - equator) <= 0 && y[2] != next[2]) {
            Scalar lo = 0, hi = terminal_fraction;
            const bool from_north = y[2] < equator;
            for (unsigned iteration = 0; iteration < 52; ++iteration) {
                const Scalar mid = (lo + hi) / 2;
                const auto candidate = advance(mid);
                if (from_north ? candidate[2] >= equator : candidate[2] <= equator)
                    hi = mid;
                else
                    lo = mid;
            }
            const auto event = advance((lo + hi) / 2);
            if (event[0] >= disk_inner / scale && event[0] <= disk_outer / scale) {
                next = event;
                captured = escaped = false;
                disk = true;
                fate = Fate::Disk;
            }
        }
        if (y[1] * next[1] < 0) ++radial_turns;
        if (y[3] * next[3] < 0) ++polar_turns;
        y = next;
        if (captured || escaped || disk) break;
    }
    if (attempts == 200000) throw std::runtime_error("reference work bound exhausted");
    const auto derivative = detail::Derivative(y, c);
    const Scalar rr = y[0], st = std::sin(y[2]), ct = std::cos(y[2]);
    const Scalar rho = std::sqrt(rr * rr + a * a);
    Scalar azimuth = y[4], time = y[5], azimuth_rate = derivative[4], time_rate = derivative[5];
    Three infinity{};
    if (fate != Fate::Capture) {
        const auto terminal_shift = detail::Shift(rr, c);
        azimuth -= terminal_shift[1];
        time -= terminal_shift[0];
        const Scalar delta = rr * rr - 2 * M * rr + a * a;
        time_rate += 4 * M * rr / delta * y[1];
        azimuth_rate += (2 * a / delta - 2 * a / (rr * rr + a * a)) * y[1];
        auto incoming = y;
        incoming[4] = azimuth;
        if (fate == Fate::Escape)
            infinity = detail::Infinity(incoming, c,
                                        refinement < .001L   ? 8192
                                        : refinement < .002L ? 4096
                                                             : 2048);
    }
    x = {time, rho * st * std::cos(azimuth), rho * st * std::sin(azimuth), rr * ct};
    const Scalar terminal_sigma = rr * rr + a * a * ct * ct;
    const Scalar transverse_rate = rr / rho * y[1] * st + rho * ct * y[3];
    k = {time_rate / terminal_sigma,
         (transverse_rate * std::cos(azimuth) - rho * st * std::sin(azimuth) * azimuth_rate) /
             terminal_sigma,
         (transverse_rate * std::sin(azimuth) + rho * st * std::cos(azimuth) * azimuth_rate) /
             terminal_sigma,
         (ct * y[1] - rr * st * y[3]) / terminal_sigma};
    const auto terminal = At(x, M, a, fate == Fate::Capture);
    const auto finite = terminal.Sky(k);
    for (auto& value : x) value *= scale;
    return {fate,
            x,
            k,
            finite,
            infinity,
            y[6] * scale,
            energy,
            angular * scale,
            carter * scale * scale,
            radial_turns,
            polar_turns,
            attempts + 1,
            fate == Fate::Capture};
}
}  // namespace sirius::test::separated_reference
