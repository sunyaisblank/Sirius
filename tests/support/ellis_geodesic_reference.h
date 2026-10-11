#pragma once

// Independent finite Ellis null orbits in units of the throat radius. No
// product metric, connection, launch, stepper, event or sky helper is used.
// ds²=-dt²+dell²+(ell²+1)dOmega², unit energy magnitude, J=(ell²+1)phi'.
// The represented past-directed affine tangent has t'=-1.
// Nakajima & Asada, PRD85 107501 (2012), equations1--9:
// https://arxiv.org/pdf/1204.3710. Finite endpoints below are derived from
// those conserved quantities; the published infinity bending angle is not
// substituted for the tracer's finite directional boundary.

#include <array>
#include <cmath>
#include <stdexcept>

namespace sirius::test::ellis_reference {

using Scalar = long double;
using Four = std::array<Scalar, 4>;

struct State {
    Four position;
    Four tangent;
    Four sky;
    Scalar affine;
    Scalar angle;
};

template <typename Function>
Scalar Integrate(Function function, Scalar start, Scalar end, unsigned panels) {
    if (panels == 0 || panels % 2 != 0) throw std::invalid_argument("reference panel count");
    const Scalar step = (end - start) / panels;
    Scalar sum = function(start) + function(end);
    for (unsigned i = 1; i < panels; ++i) sum += (i % 2 == 0 ? 2 : 4) * function(start + i * step);
    return sum * step / 3;
}

// rho is isotropic radius/b0; impact is signed J/b0. Both rays start inward
// at phi=0. Reflection uses ell=sqrt(J²-1)/cos(u), removing its turning
// singularity. Transmission uses ell=tan(u), regular across the throat.
inline State Orbit(Scalar launch_rho, Scalar terminal_rho, Scalar impact, bool opposite,
                   unsigned panels) {
    const Scalar j = std::abs(impact);
    const Scalar launch_l = launch_rho - 1 / (4 * launch_rho);
    const Scalar terminal_l = terminal_rho - 1 / (4 * terminal_rho);
    if (!std::isfinite(launch_rho) || !std::isfinite(terminal_rho) || !std::isfinite(j) ||
        !(launch_rho > 0.5L) || !(terminal_rho > 0) || j == 1 ||
        !(j < launch_rho + 1 / (4 * launch_rho)))
        throw std::invalid_argument("reference finite Ellis orbit domain");
    Scalar angle = 0, affine = 0;
    const bool reflection = j > 1;
    if (reflection) {
        const Scalar turn = std::sqrt(j * j - 1);
        if (!(launch_l > turn) || !(terminal_l > turn) || opposite)
            throw std::invalid_argument("reference reflecting endpoint");
        const Scalar u0 = std::acos(turn / launch_l), u1 = std::acos(turn / terminal_l);
        const auto angular = [j](Scalar u) {
            const Scalar s = std::sin(u) / j;
            return 1 / std::sqrt(1 - s * s);
        };
        const auto length = [j](Scalar u) {
            const Scalar s = std::sin(u) / j, c = std::cos(u);
            return j * std::sqrt(1 - s * s) / (c * c);
        };
        angle = Integrate(angular, 0, u0, panels) + Integrate(angular, 0, u1, panels);
        affine = Integrate(length, 0, u0, panels) + Integrate(length, 0, u1, panels);
    } else {
        if (terminal_l > 0 || opposite != (terminal_l < 0))
            throw std::invalid_argument("reference transmitting endpoint");
        const Scalar u0 = std::atan(terminal_l), u1 = std::atan(launch_l);
        const auto angular = [j](Scalar u) {
            const Scalar c = std::cos(u);
            return j / std::sqrt(1 - j * j * c * c);
        };
        const auto length = [j](Scalar u) {
            const Scalar c = std::cos(u);
            return 1 / (c * c * std::sqrt(1 - j * j * c * c));
        };
        angle = Integrate(angular, u0, u1, panels);
        affine = Integrate(length, u0, u1, panels);
    }
    angle = std::copysign(angle, impact);
    const Scalar radius = terminal_rho + 1 / (4 * terminal_rho);
    const Scalar radial = (reflection ? 1 : -1) * std::sqrt(1 - j * j / (radius * radius));
    const Scalar angular = impact / radius;
    const Scalar conformal = radius / terminal_rho;
    const Scalar c = std::cos(angle), s = std::sin(angle);
    // On the opposite end inversion reflects the radial component, while
    // keeping the tangential component: the independent I-2nn^T Jacobian.
    const Scalar sky_radial = opposite ? -radial : radial;
    return {{-affine, terminal_rho * c, terminal_rho * s, 0},
            {-1, (radial * c - angular * s) / conformal, (radial * s + angular * c) / conformal, 0},
            {0, sky_radial * c - angular * s, sky_radial * s + angular * c, 0},
            affine,
            angle};
}

}  // namespace sirius::test::ellis_reference
