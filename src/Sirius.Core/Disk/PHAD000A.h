// =============================================================================
// PHAD000A.h - Accretion Disk Model Interface
// Component ID: PHAD000A (Physics/Accretion Disk/Interface)
// =============================================================================
//
// PURPOSE:
// Defines the common interface for all accretion disk models in Sirius.
// This allows the rendering pipeline to work with different disk models
// (Novikov-Thorne, Polish Donut, MAD disk) through a unified interface.
//
// IMPLEMENTATIONS:
// - PHAD001A: Novikov-Thorne thin disk (relativistic, radiatively efficient)
// - PHAD002A: Polish Donut thick disk (geometrically thick torus)
// - PHAD003A: Generic parametric disk (configurable temperature/density profiles)
//
// MATHEMATICAL FOUNDATION:
// Accretion disk models provide:
// - Temperature profile T(r, z) for blackbody emission
// - Density profile ρ(r, z) for optical depth calculation
// - Geometry bounds for ray intersection testing
//
// TESTS: TSPH005A.cpp (disk model tests)
// =============================================================================

#ifndef PHAD000A_H
#define PHAD000A_H

#include <string>
#include <cmath>

namespace Sirius {

// =============================================================================
// Disk Intersection Result
// =============================================================================
struct DiskIntersection {
    bool hit = false;           // True if ray intersects disk
    double r = 0.0;             // Radial coordinate at intersection [M]
    double theta = 0.0;         // Polar angle at intersection [rad]
    double phi = 0.0;           // Azimuthal angle at intersection [rad]
    double z = 0.0;             // Height above midplane [M]
    double temperature = 0.0;   // Local temperature [K]
    double density = 0.0;       // Local density [kg/m³ or dimensionless]
    double optical_depth = 0.0; // Optical depth to surface
};

// =============================================================================
// Disk Emission Properties
// =============================================================================
struct DiskEmission {
    double radiance = 0.0;      // Specific intensity I_ν [W/(m²·sr·Hz)]
    double temperature = 0.0;   // Effective temperature [K]
    double redshift = 0.0;      // Gravitational + Doppler redshift factor g
    double beaming = 1.0;       // Relativistic beaming factor δ⁴
    double flux = 0.0;          // Radiative flux F [W/m²]
};

// =============================================================================
// IDiskModel - Accretion Disk Model Interface
// =============================================================================
class IDiskModel {
public:
    virtual ~IDiskModel() = default;

    // =========================================================================
    // Identification
    // =========================================================================

    /// @brief Get model name (e.g., "Novikov-Thorne", "Polish Donut")
    virtual const char* modelName() const = 0;

    // =========================================================================
    // Geometry
    // =========================================================================

    /// @brief Inner radius of disk (typically ISCO)
    /// @return Inner radius in units of M
    virtual double innerRadius() const = 0;

    /// @brief Outer radius of disk
    /// @return Outer radius in units of M
    virtual double outerRadius() const = 0;

    /// @brief Half-thickness at given radius
    /// @param r Radial coordinate [M]
    /// @return Half-thickness H(r) [M], 0 for infinitely thin disk
    virtual double halfThickness(double r) const = 0;

    /// @brief Check if a point is within the disk
    /// @param r Radial coordinate [M]
    /// @param theta Polar angle from spin axis [rad]
    /// @return true if point is inside disk volume
    virtual bool isInDisk(double r, double theta) const = 0;

    // =========================================================================
    // Thermodynamic Properties
    // =========================================================================

    /// @brief Temperature at given location
    /// @param r Radial coordinate [M]
    /// @param z Height above midplane [M] (optional, default 0)
    /// @return Temperature in Kelvin
    virtual double temperature(double r, double z = 0) const = 0;

    /// @brief Density at given location
    /// @param r Radial coordinate [M]
    /// @param z Height above midplane [M] (optional, default 0)
    /// @return Density (units depend on model - typically dimensionless or kg/m³)
    virtual double density(double r, double z = 0) const = 0;

    /// @brief Radiative flux at disk surface
    /// @param r Radial coordinate [M]
    /// @return Flux F [W/m²] in physical units
    virtual double flux(double r) const = 0;

    // =========================================================================
    // Kinematics
    // =========================================================================

    /// @brief Angular velocity of disk material (Keplerian or sub-Keplerian)
    /// @param r Radial coordinate [M]
    /// @return Angular velocity Ω [rad/s or 1/M]
    virtual double angularVelocity(double r) const = 0;

    /// @brief Radial velocity (accretion flow)
    /// @param r Radial coordinate [M]
    /// @return Radial velocity v_r (negative for infall)
    virtual double radialVelocity(double r) const {
        (void)r;
        return 0.0;  // Default: no radial motion
    }

    // =========================================================================
    // Relativistic Effects
    // =========================================================================

    /// @brief Compute Doppler + gravitational redshift factor g = ν_obs/ν_emit
    /// @param r Radial coordinate [M]
    /// @param phi Azimuthal angle [rad]
    /// @param observer_theta Observer inclination from spin axis [rad]
    /// @return Redshift factor g (g < 1 means redshift, g > 1 means blueshift)
    virtual double redshiftFactor(double r, double phi, double observer_theta) const {
        // Default implementation: simple Doppler approximation
        // Subclasses should override with full relativistic calculation
        (void)r; (void)phi; (void)observer_theta;
        return 1.0;
    }

    /// @brief Relativistic beaming factor δ = [γ(1 - β·n̂)]^(-1)
    /// @param r Radial coordinate [M]
    /// @param phi Azimuthal angle [rad]
    /// @param observer_theta Observer inclination [rad]
    /// @return Beaming factor (>1 for approaching, <1 for receding)
    virtual double beamingFactor(double r, double phi, double observer_theta) const {
        (void)r; (void)phi; (void)observer_theta;
        return 1.0;
    }

    // =========================================================================
    // Configuration
    // =========================================================================

    /// @brief Set black hole parameters
    /// @param mass Mass in solar masses or geometric units
    /// @param spin Dimensionless spin a* = a/M
    virtual void setBlackHoleParameters(double mass, double spin) = 0;

    /// @brief Get black hole mass
    virtual double blackHoleMass() const = 0;

    /// @brief Get black hole spin
    virtual double blackHoleSpin() const = 0;
};

// =============================================================================
// Utility Functions
// =============================================================================

/// @brief Convert height z to polar angle theta (for thin disk)
/// @param r Radial coordinate [M]
/// @param z Height above midplane [M]
/// @return Polar angle θ = π/2 - atan(z/r)
inline double zToTheta(double r, double z) {
    if (r <= 0) return M_PI / 2.0;
    return M_PI / 2.0 - std::atan(z / r);
}

/// @brief Convert polar angle theta to height z (for thin disk)
/// @param r Radial coordinate [M]
/// @param theta Polar angle [rad]
/// @return Height z = r * cot(θ) relative to equatorial plane
inline double thetaToZ(double r, double theta) {
    double sin_theta = std::sin(theta);
    if (std::abs(sin_theta) < 1e-10) return 0;
    return r * std::cos(theta) / sin_theta;
}

} // namespace Sirius

#endif // PHAD000A_H
