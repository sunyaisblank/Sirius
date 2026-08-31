#pragma once

// Exact equatorial circular-orbit primitives for the asymptotically flat,
// uncharged Schwarzschild/Kerr family.  CPU metric, disk, frequency-transfer,
// and ray-tracer consumers share these functions so their ISCO and emitter
// angular velocity cannot become separate physical authorities.  The
// double-precision Boyer-Lindquist oracle and the Slang device implementation
// remain independent comparators.

#include <cmath>
#include <optional>

namespace sirius::core::relativity {

// Bardeen-Press-Teukolsky equatorial ISCO.  mass and spin are dimensional
// geometric quantities with |spin| <= mass.  Positive spin selects the
// co-rotating branch and negative spin the counter-rotating branch for the
// disk's fixed positive-phi orbital orientation.
[[nodiscard]] inline std::optional<double> TryKerrIscoRadius(double mass, double spin) {
    if (!std::isfinite(mass) || !(mass > 0.0) || !std::isfinite(spin) || std::abs(spin) > mass) {
        return std::nullopt;
    }

    const double dimensionless_spin = spin / mass;
    const double spin_magnitude = std::abs(dimensionless_spin);
    if (spin_magnitude == 0.0) return 6.0 * mass;

    const double z1 = 1.0 + std::cbrt(1.0 - spin_magnitude * spin_magnitude) *
                                (std::cbrt(1.0 + spin_magnitude) + std::cbrt(1.0 - spin_magnitude));
    const double z2 = std::sqrt(3.0 * spin_magnitude * spin_magnitude + z1 * z1);
    const double radicand = (3.0 - z1) * (3.0 + z1 + 2.0 * z2);
    if (!std::isfinite(z1) || !std::isfinite(z2) || !std::isfinite(radicand) || radicand < 0.0) {
        return std::nullopt;
    }

    const double branch = dimensionless_spin >= 0.0 ? -1.0 : 1.0;
    const double radius = mass * (3.0 + z2 + branch * std::sqrt(radicand));
    if (!std::isfinite(radius) || !(radius > 0.0)) return std::nullopt;
    return radius;
}

// Coordinate angular velocity of the same fixed positive-phi equatorial
// circular-geodesic branch,
//
//   Omega = sqrt(M) / (r^(3/2) + a sqrt(M)).
//
// The caller still owns the stable/timelike orbit domain (normally r >= ISCO).
// This primitive rejects an undefined or reversed denominator rather than
// manufacturing a finite orbital rate.
[[nodiscard]] inline std::optional<double> TryKerrCircularOrbitAngularVelocity(double mass,
                                                                               double spin,
                                                                               double radius) {
    if (!std::isfinite(mass) || !(mass > 0.0) || !std::isfinite(spin) || std::abs(spin) > mass ||
        !std::isfinite(radius) || !(radius > 0.0)) {
        return std::nullopt;
    }

    const double sqrt_mass = std::sqrt(mass);
    const double denominator = std::pow(radius, 1.5) + spin * sqrt_mass;
    if (!std::isfinite(denominator) || !(denominator > 0.0)) return std::nullopt;

    const double angular_velocity = sqrt_mass / denominator;
    if (!std::isfinite(angular_velocity) || !(angular_velocity > 0.0)) return std::nullopt;
    return angular_velocity;
}

}  // namespace sirius::core::relativity
