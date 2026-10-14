#pragma once

// Common interface for accretion disk models: temperature and density profiles,
// geometry bounds for ray intersection, kinematics, and relativistic factors.
// Ported from PHAD000A.h.

#include <cmath>
#include <numbers>
#include <optional>
#include <string>

namespace sirius::core {

// Accretion disk model interface.
class IDiskModel {
  public:
    virtual ~IDiskModel() = default;

    // Model name (e.g. "Novikov-Thorne").
    virtual const char* ModelName() const = 0;

    // Inner radius (typically the ISCO), in units of M.
    virtual double InnerRadius() const = 0;

    // Outer radius, in units of M.
    virtual double OuterRadius() const = 0;

    // Half-thickness H(r) [M]; 0 for an infinitely thin disk.
    virtual double HalfThickness(double r) const = 0;

    // Whether a point (r, theta) is inside the disk volume.
    virtual bool IsInDisk(double r, double theta) const = 0;

    // Temperature [K] at (r, z); z is the height above the midplane [M].
    virtual double Temperature(double r, double z = 0) const = 0;

    // Radiative flux F [W/m^2] at the surface.
    virtual double Flux(double r) const = 0;

    // Angular velocity Omega of the disk material at r.
    virtual double AngularVelocity(double r) const = 0;

    // Radial (accretion) velocity; negative for infall.
    virtual double RadialVelocity(double r) const {
        (void)r;
        return 0.0;  // Default: no radial motion.
    }

    // Set the black hole mass (solar masses or geometric) and spin a* = a/M.
    virtual void SetBlackHoleParameters(double mass, double spin) = 0;

    virtual double BlackHoleMass() const = 0;
    virtual double BlackHoleSpin() const = 0;
};

// Height z to polar angle theta = pi/2 - atan(z/r) for a thin disk. The
// cylindrical radius must be positive: the axis has no owned azimuthal disk
// chart and may not be silently projected onto the equatorial plane.
[[nodiscard]] inline std::optional<double> ZToTheta(double r, double z) {
    if (!std::isfinite(r) || !std::isfinite(z) || !(r > 0.0)) return std::nullopt;
    const double theta = std::numbers::pi / 2.0 - std::atan(z / r);
    if (!std::isfinite(theta) || !(theta > 0.0) || !(theta < std::numbers::pi)) {
        return std::nullopt;
    }
    return theta;
}

// Polar angle theta to height z = r cot(theta) for a thin disk. Both polar
// axes are coordinate singularities and decline instead of impersonating the
// finite equatorial height z=0.
[[nodiscard]] inline std::optional<double> ThetaToZ(double r, double theta) {
    if (!std::isfinite(r) || !(r > 0.0) || !std::isfinite(theta) || !(theta > 0.0) ||
        !(theta < std::numbers::pi)) {
        return std::nullopt;
    }
    const double height = r * std::cos(theta) / std::sin(theta);
    if (!std::isfinite(height)) return std::nullopt;
    return height;
}

}  // namespace sirius::core
