#pragma once

// Independent finite spherical null orbits in one declared XY plane. Lengths
// and affine parameter are in a caller-declared natural unit; Lambda is in its
// inverse square. No product metric, launch, connection, chart, integrator,
// event interpolant or sky helper supplies an expected answer.
// Lebedev & Lake, arXiv:1609.05183, equations 1--14:
// https://arxiv.org/pdf/1609.05183. With f=1-2M/r-Lambda*r^2/3,
// C=p_t>0 for a past ray, L=x*p_y-y*p_x, and W^2=C^2-f*L^2/r^2.
// Affine/angle integrands are 1/W and L/(r^2 W), at finite endpoints.

#include <array>
#include <cmath>
#include <numbers>
#include <stdexcept>

namespace sirius::test::kottler_reference {

using Scalar = long double;
using Four = std::array<Scalar, 4>;
using Three = std::array<Scalar, 3>;

struct Horizons {
    Scalar capture, cosmological;
};

inline Horizons Roots(Scalar mass, Scalar lambda) {
    if (!std::isfinite(mass) || !std::isfinite(lambda) || !(mass >= 0) || !(lambda > 0) ||
        !(9 * lambda * mass * mass < 1))
        throw std::invalid_argument("reference sub-Nariai spherical domain");
    if (mass == 0) return {0, std::sqrt(3 / lambda)};
    // Independent trigonometric cubic solution, rather than the product's
    // scale-safe bisection authority. This finite corpus is away from Nariai.
    const Scalar phase = std::acos(-3 * mass * std::sqrt(lambda)) / 3;
    const Scalar factor = 2 / std::sqrt(lambda);
    return {factor * std::cos(phase + 4 * std::numbers::pi_v<Scalar> / 3),
            factor * std::cos(phase)};
}

inline Scalar Lapse(Scalar radius, Scalar mass, Scalar lambda) {
    return 1 - 2 * mass / radius - lambda * radius * radius / 3;
}

template <typename Function>
Three Integrate(Function function, Scalar start, Scalar end, unsigned panels) {
    if (panels == 0 || panels % 2 != 0) throw std::invalid_argument("reference panel count");
    const Scalar step = (end - start) / panels;
    const auto first = function(start), last = function(end);
    Three sum{};
    for (unsigned field = 0; field < 3; ++field) sum[field] = first[field] + last[field];
    for (unsigned i = 1; i < panels; ++i) {
        const auto value = function(start + i * step);
        for (unsigned field = 0; field < 3; ++field)
            sum[field] += (i % 2 == 0 ? 2 : 4) * value[field];
    }
    for (auto& value : sum) value *= step / 3;
    return sum;
}

// Declared outgoing-time gauge t_out=t_in-2*integral_anchor^r (1-f)/f dr,
// anchor=(r_capture+r_cosmological)/2. This independently quadratures the
// smooth exterior interval, rather than invoking the production chart map or
// its partial-fraction logarithms. It is only needed for captures.
inline Scalar OutgoingShift(Scalar radius, Scalar mass, Scalar lambda, unsigned panels) {
    const auto roots = Roots(mass, lambda);
    if (!(mass > 0) || !(radius > roots.capture) || !(radius < roots.cosmological))
        throw std::invalid_argument("reference outgoing gauge domain");
    const Scalar anchor = (roots.capture + roots.cosmological) / 2;
    const auto integral = Integrate(
        [=](Scalar r) {
            const Scalar f = Lapse(r, mass, lambda);
            return Three{(1 - f) / f, 0, 0};
        },
        anchor, radius, panels);
    return -2 * integral[0];
}

// Public finite sky uses an ingoing Eulerian observer with fixed Cartesian
// x/y/z coordinate seeds, Gram-Schmidt in I+H*n*n^T. A rotated spherical
// radial/azimuth frame would have a different published angle.
inline Four Sky(const Four& position, const Four& tangent, Scalar mass, Scalar lambda) {
    const Scalar r = std::hypot(position[1], position[2]);
    const Scalar h = 1 - Lapse(r, mass, lambda);
    const Three radial{position[1] / r, position[2] / r, 0};
    const auto dot = [](const Three& a, const Three& b) {
        return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    };
    const auto inner = [&](const Three& a, const Three& b) {
        return dot(a, b) + h * dot(a, radial) * dot(b, radial);
    };
    std::array<Three, 3> basis{};
    for (unsigned axis = 0; axis < 3; ++axis) {
        basis[axis][axis] = 1;
        for (unsigned prior = 0; prior < axis; ++prior) {
            const Scalar projection = inner(basis[axis], basis[prior]);
            for (unsigned i = 0; i < 3; ++i) basis[axis][i] -= projection * basis[prior][i];
        }
        const Scalar length = std::sqrt(inner(basis[axis], basis[axis]));
        for (auto& value : basis[axis]) value /= length;
    }
    const Three k{tangent[1], tangent[2], tangent[3]};
    const Scalar frequency = -tangent[0] / std::sqrt(1 + h);
    if (!(frequency > 0)) throw std::runtime_error("reference past Eulerian sky cone");
    Four sky{};
    for (unsigned axis = 0; axis < 3; ++axis)
        sky[axis + 1] =
            (dot(k, basis[axis]) + h * (tangent[0] + dot(k, radial)) * dot(basis[axis], radial)) /
            frequency;
    const Scalar norm = std::sqrt(sky[1] * sky[1] + sky[2] * sky[2] + sky[3] * sky[3]);
    for (unsigned axis = 1; axis < 4; ++axis) sky[axis] /= norm;
    return sky;
}

enum class Fate { Capture, Escape };
struct State {
    Fate fate;
    Four position, tangent, sky;
    Scalar affine, angle, killing_magnitude, angular_momentum, turning_radius;
};

// Camera components are local radial and azimuth directions at (r,pi/2,0),
// hence azimuth is +Cartesian y, L=r*n_phi; theta is zero. The components may
// be measured binary32 pinhole output, but launch normalization and the
// Eulerian metric/tetrad are independently derived here. Camera local
// frequency 1 does NOT set the Killing magnitude C to 1.
inline State Orbit(Scalar mass, Scalar lambda, Scalar launch_radius, Scalar escape_radius,
                   Scalar radial_direction, Scalar azimuth_direction, unsigned panels) {
    const auto roots = Roots(mass, lambda);
    if (!std::isfinite(launch_radius) || !std::isfinite(escape_radius) ||
        !(launch_radius > roots.capture) || !(launch_radius < roots.cosmological) ||
        !(escape_radius > launch_radius) || !(escape_radius <= roots.cosmological * (1 + 1.0e-14L)))
        throw std::invalid_argument("reference finite causal endpoint domain");
    const Scalar norm = std::hypot(radial_direction, azimuth_direction);
    if (!std::isfinite(norm) || !(norm > 0))
        throw std::invalid_argument("reference camera direction");
    const Scalar n = radial_direction / norm, transverse = azimuth_direction / norm;
    const Scalar h0 = 1 - Lapse(launch_radius, mass, lambda);
    const Scalar c = (1 + h0 * n) / std::sqrt(1 + h0);
    const Scalar angular = launch_radius * transverse;
    const Scalar l2 = angular * angular, a2 = c * c + lambda * l2 / 3;
    const Scalar radial_start = (h0 + n) / std::sqrt(1 + h0);
    if (!(c > 0) || radial_start == 0)
        throw std::invalid_argument("reference unresolved launch branch");
    const bool inward = radial_start < 0;
    const bool reflection = inward && (mass == 0 ? l2 > 0 : a2 - l2 / (27 * mass * mass) < 0);
    const bool capture = inward && !reflection;
    if (capture && mass == 0)
        throw std::invalid_argument("reference radial de Sitter origin not in this finite corpus");
    Scalar turn = 0;
    if (reflection) {
        if (mass == 0) {
            turn = std::abs(angular) / std::sqrt(a2);
        } else {
            Scalar lower = 3 * mass, upper = launch_radius;
            for (unsigned i = 0; i < 128; ++i) {
                const Scalar middle = (lower + upper) / 2;
                const Scalar potential =
                    a2 - l2 / (middle * middle) + 2 * mass * l2 / (middle * middle * middle);
                if (potential > 0)
                    upper = middle;
                else
                    lower = middle;
            }
            turn = (lower + upper) / 2;
        }
    }
    const Scalar terminal = capture ? roots.capture : escape_radius;
    // The regular past branch is outgoing-chart inward capture or ingoing-chart
    // outward escape. Rationalizing removes the cancelling horizon pole:
    // k^t=-W-L^2/[r^2(C+W)]. No lapse clamp is needed at either physical horizon.
    const auto time_rate = [=](Scalar r, Scalar w, bool inward_ingoing) {
        if (!inward_ingoing) return -w - l2 / (r * r * (c + w));
        const Scalar f = Lapse(r, mass, lambda);
        return (-c - (1 - f) * w) / f;
    };
    const auto direct = [=](Scalar r) {
        const Scalar w = std::sqrt(a2 - l2 / (r * r) + 2 * mass * l2 / (r * r * r));
        return Three{1 / w, angular / (r * r * w), time_rate(r, w, false) / w};
    };
    Three integral{};
    if (reflection) {
        // r=turn+s^2. Factor the polynomial at its independently found simple
        // root; dr/W=2/sqrt((a2*(r^2+r*turn+turn^2)-L^2)/r^3), including s=0.
        const auto leg = [&](Scalar end, bool inward_ingoing) {
            return Integrate(
                [=](Scalar s) {
                    const Scalar r = turn + s * s;
                    const Scalar quotient =
                        (a2 * (r * r + r * turn + turn * turn) - l2) / (r * r * r);
                    const Scalar root = std::sqrt(quotient), measure = 2 / root;
                    const Scalar w = s * root;
                    return Three{measure, angular * measure / (r * r),
                                 time_rate(r, w, inward_ingoing) * measure};
                },
                0, std::sqrt(end - turn), panels);
        };
        const auto inner = leg(launch_radius, true), outer = leg(escape_radius, false);
        for (unsigned field = 0; field < 3; ++field) integral[field] = inner[field] + outer[field];
    } else {
        integral = capture ? Integrate(direct, terminal, launch_radius, panels)
                           : Integrate(direct, launch_radius, terminal, panels);
    }
    const Scalar w = std::sqrt(a2 - l2 / (terminal * terminal) +
                               2 * mass * l2 / (terminal * terminal * terminal));
    const Scalar radial = capture ? -w : w;
    const Scalar cosine = std::cos(integral[1]), sine = std::sin(integral[1]);
    State result;
    result.fate = capture ? Fate::Capture : Fate::Escape;
    result.position = {
        integral[2] + (capture ? OutgoingShift(launch_radius, mass, lambda, panels) : 0),
        terminal * cosine, terminal * sine, 0};
    result.tangent = {time_rate(terminal, w, false), radial * cosine - angular * sine / terminal,
                      radial * sine + angular * cosine / terminal, 0};
    result.sky = capture ? Four{} : Sky(result.position, result.tangent, mass, lambda);
    result.affine = integral[0];
    result.angle = integral[1];
    result.killing_magnitude = c;
    result.angular_momentum = angular;
    result.turning_radius = turn;
    return result;
}

}  // namespace sirius::test::kottler_reference
