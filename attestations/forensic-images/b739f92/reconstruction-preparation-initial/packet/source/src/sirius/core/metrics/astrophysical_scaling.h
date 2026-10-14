#pragma once

// Astrophysically-scaled supermassive black hole parameters: conversion between
// geometric units (G=c=1) and SI, plus horizon, ISCO, ergosphere, and Eddington
// scalings. References: EHT Collaboration (2019) ApJL 875, L1-L6; Thorne (1974)
// for the a* <= 0.998 spin limit. Ported from PHSM001A.h.

#include <algorithm>
#include <cmath>

namespace sirius::core {

// Local SI constant set for the scaling relations. Kept separate from the
// constants authority because the legacy code carried its own copy; the values
// are behaviour-preserving and use domain notation.
namespace physical_constants {
constexpr double G = 6.67430e-11;            // Gravitational constant [m^3/(kg s^2)]
constexpr double c = 2.99792458e8;           // Speed of light [m/s]
constexpr double M_sun = 1.98892e30;         // Solar mass [kg]
constexpr double sigma_T = 6.6524587e-29;    // Thomson cross-section [m^2]
constexpr double m_p = 1.67262192e-27;       // Proton mass [kg]
constexpr double pc = 3.08567758e16;         // Parsec [m]
constexpr double Mpc = 3.08567758e22;        // Megaparsec [m]
constexpr double year = 3.15576e7;           // Year [s]
constexpr double k_B = 1.380649e-23;         // Boltzmann constant [J/K]
constexpr double h = 6.62607015e-34;         // Planck constant [J s]
constexpr double sigma_SB = 5.670374419e-8;  // Stefan-Boltzmann [W/(m^2 K^4)]
                                             // (CODATA 2018, matches
                                             // constants::physical::kStefanBoltzmann)
}  // namespace physical_constants

// Supermassive black hole parameters with derived geometric and SI scalings.
struct SmbhParams {
    // Primary parameters (user-specified)
    double mass_solar = 1.0e6;      // Mass in solar masses [M_sun]
    float spin_parameter = 0.0f;    // Dimensionless spin a* = a/M in [-0.998, 0.998]
    float inclination_deg = 45.0f;  // Observer inclination angle [degrees]
    float distance_mpc = 10.0f;     // Distance in megaparsecs

    // Derived quantities (computed via ComputeDerived())
    double m_kg = 0.0;         // Mass in kilograms
    double m_geometric = 0.0;  // GM/c^2 in meters (gravitational radius)
    double r_g = 0.0;          // Gravitational radius GM/c^2 [m]
    double r_s = 0.0;          // Schwarzschild radius 2GM/c^2 [m]
    double r_isco = 0.0;       // ISCO radius in meters
    double r_isco_m = 0.0;     // ISCO radius in units of M
    double r_horizon = 0.0;    // Event horizon radius in meters
    double r_horizon_m = 0.0;  // Event horizon radius in units of M
    double distance_m = 0.0;   // Distance in meters
    double l_edd = 0.0;        // Eddington luminosity [W]
    double t_g = 0.0;          // Gravitational time unit GM/c^3 [s]

    // Invariants (enforced by clamping, not assertion):
    //   mass_solar in [10^5, 10^11] for SMBH class
    //   |spin_parameter| <= 0.998 (Thorne limit)
    //   distance_mpc > 0
    //   inclination_deg in [0, 90]

    // Compute all derived quantities from the primary parameters.
    void ComputeDerived() {
        using namespace physical_constants;

        // Clamp spin to Thorne limit
        if (spin_parameter > 0.998f) spin_parameter = 0.998f;
        if (spin_parameter < -0.998f) spin_parameter = -0.998f;

        // Mass conversions
        m_kg = mass_solar * M_sun;
        m_geometric = G * m_kg / (c * c);  // GM/c^2 in meters
        r_g = m_geometric;
        r_s = 2.0 * m_geometric;

        // Distance
        distance_m = distance_mpc * Mpc;

        // ISCO radius (depends on spin)
        r_isco_m = ComputeIsco(spin_parameter);
        r_isco = r_isco_m * m_geometric;

        // Event horizon radius
        // r_+ = M + sqrt(M^2 - a^2) where a = a* * M
        float a = std::abs(spin_parameter);
        r_horizon_m = 1.0 + std::sqrt(1.0 - a * a);
        r_horizon = r_horizon_m * m_geometric;

        // Eddington luminosity: L_Edd = 4 pi G M m_p c / sigma_T
        l_edd = 4.0 * M_PI * G * m_kg * m_p * c / sigma_T;

        // Gravitational time unit: GM/c^3
        t_g = m_geometric / c;
    }

    // ISCO radius in units of M:
    //   r_ISCO = M {3 + Z2 -/+ sqrt[(3 - Z1)(3 + Z1 + 2 Z2)]}
    // where -/+ selects prograde/retrograde orbits.
    static double ComputeIsco(float a_star) {
        double a = std::abs(static_cast<double>(a_star));
        if (a > 1.0) a = 1.0;

        double Z1 = 1.0 + std::cbrt(1.0 - a * a) * (std::cbrt(1.0 + a) + std::cbrt(1.0 - a));
        double Z2 = std::sqrt(3.0 * a * a + Z1 * Z1);

        // Guard against floating-point errors making (3 - Z1) slightly negative
        double factor1 = std::max(0.0, 3.0 - Z1);
        double factor2 = 3.0 + Z1 + 2.0 * Z2;
        double sqrt_term = std::sqrt(factor1 * factor2);

        double r_isco;
        if (a_star >= 0.0f) {
            // Prograde orbit
            r_isco = 3.0 + Z2 - sqrt_term;
        } else {
            // Retrograde orbit
            r_isco = 3.0 + Z2 + sqrt_term;
        }

        return r_isco;  // In units of M
    }

    // Inner horizon radius in units of M: r_- = M - sqrt(M^2 - a^2).
    static double ComputeInnerHorizon(float a_star) {
        float a = std::abs(a_star);
        if (a > 1.0f) a = 1.0f;
        return 1.0 - std::sqrt(1.0 - a * a);
    }

    // Ergosphere radius at polar angle theta: r_ergo = M + sqrt(M^2 - a^2 cos^2 theta).
    double ComputeErgosphereRadius(double theta_rad) const {
        double a = std::abs(static_cast<double>(spin_parameter));
        double cos_theta = std::cos(theta_rad);
        return m_geometric * (1.0 + std::sqrt(1.0 - a * a * cos_theta * cos_theta));
    }

    // Convert a radius from M units to meters.
    double RadiusToMeters(double r_M) const { return r_M * m_geometric; }

    // Convert a radius from meters to M units.
    double RadiusFromMeters(double r_m) const { return r_m / m_geometric; }

    // Angular size of r_g as seen from the observer [radians].
    double AngularSizeOfRg() const { return m_geometric / distance_m; }

    // Angular size of the shadow (~sqrt(27) M for Schwarzschild; Kerr is
    // similar but depends on spin and inclination).
    double ShadowAngularSize() const { return std::sqrt(27.0) * AngularSizeOfRg(); }

    // Preset configurations for known SMBHs.

    // M87* configuration (EHT 2019).
    static SmbhParams M87Star() {
        SmbhParams p;
        p.mass_solar = 6.5e9;
        p.spin_parameter = 0.90f;
        p.inclination_deg = 17.0f;
        p.distance_mpc = 16.8f;
        p.ComputeDerived();
        return p;
    }

    // Sagittarius A* configuration.
    static SmbhParams SgrAStar() {
        SmbhParams p;
        p.mass_solar = 4.0e6;
        p.spin_parameter = 0.50f;
        p.inclination_deg = 45.0f;
        p.distance_mpc = 0.0082f;  // 8.2 kpc
        p.ComputeDerived();
        return p;
    }

    // Gargantua-like configuration (Interstellar).
    static SmbhParams Gargantua() {
        SmbhParams p;
        p.mass_solar = 1.0e8;
        p.spin_parameter = 0.998f;  // Thorne limit; 0.9999 was silently clamped here anyway
        p.inclination_deg = 85.0f;  // Nearly edge-on
        p.distance_mpc = 0.01f;     // Close for visualization
        p.ComputeDerived();
        return p;
    }
};

}  // namespace sirius::core
