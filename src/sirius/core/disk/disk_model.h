#pragma once

// Common interface for accretion disk models: temperature and density profiles,
// geometry bounds for ray intersection, kinematics, and relativistic factors.
// Ported from PHAD000A.h.

#include <cmath>
#include <string>

namespace sirius::core {

// Result of a ray/disk intersection test.
struct DiskIntersection {
    bool hit = false;            // True when the ray intersects the disk.
    double r = 0.0;              // Radial coordinate at intersection [M].
    double theta = 0.0;          // Polar angle at intersection [rad].
    double phi = 0.0;            // Azimuthal angle at intersection [rad].
    double z = 0.0;              // Height above the midplane [M].
    double temperature = 0.0;    // Local temperature [K].
    double density = 0.0;        // Local density [kg/m^3 or dimensionless].
    double optical_depth = 0.0;  // Optical depth to the surface.
};

// Emission properties at a disk point.
struct DiskEmission {
    double radiance = 0.0;     // Specific intensity I_nu [W/(m^2 sr Hz)].
    double temperature = 0.0;  // Effective temperature [K].
    double redshift = 0.0;     // Gravitational + Doppler redshift factor g.
    double beaming = 1.0;      // Relativistic beaming factor delta^4.
    double flux = 0.0;         // Radiative flux F [W/m^2].
};

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

    // Density at (r, z); units depend on the model.
    virtual double Density(double r, double z = 0) const = 0;

    // Radiative flux F [W/m^2] at the surface.
    virtual double Flux(double r) const = 0;

    // Angular velocity Omega of the disk material at r.
    virtual double AngularVelocity(double r) const = 0;

    // Radial (accretion) velocity; negative for infall.
    virtual double RadialVelocity(double r) const {
        (void)r;
        return 0.0;  // Default: no radial motion.
    }

    // Doppler + gravitational redshift factor g = nu_obs/nu_emit (g < 1 redshift,
    // g > 1 blueshift). Subclasses override with the full relativistic form.
    virtual double RedshiftFactor(double r, double phi, double observer_theta) const {
        (void)r;
        (void)phi;
        (void)observer_theta;
        return 1.0;
    }

    // Relativistic beaming factor delta = [gamma (1 - beta . n)]^{-1}.
    virtual double BeamingFactor(double r, double phi, double observer_theta) const {
        (void)r;
        (void)phi;
        (void)observer_theta;
        return 1.0;
    }

    // Set the black hole mass (solar masses or geometric) and spin a* = a/M.
    virtual void SetBlackHoleParameters(double mass, double spin) = 0;

    virtual double BlackHoleMass() const = 0;
    virtual double BlackHoleSpin() const = 0;
};

// Height z to polar angle theta = pi/2 - atan(z/r) for a thin disk.
inline double ZToTheta(double r, double z) {
    if (r <= 0) return M_PI / 2.0;
    return M_PI / 2.0 - std::atan(z / r);
}

// Polar angle theta to height z = r cot(theta) for a thin disk.
inline double ThetaToZ(double r, double theta) {
    double sin_theta = std::sin(theta);
    if (std::abs(sin_theta) < 1e-10) return 0;
    return r * std::cos(theta) / sin_theta;
}

}  // namespace sirius::core
