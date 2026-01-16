// PHSM001A.h - Supermassive Black Hole Parameters
// Component ID: PHSM001A (Physics/Metric/SMBHParams)
//
// Astrophysically-scaled supermassive black hole parameters.
// Provides scaling between geometric units (G=c=1) and SI units.
//
// PHYSICAL FOUNDATION:
// ====================
// SMBH are found in galactic nuclei with masses M ~ 10^6 - 10^10 M_sun.
// Observable examples:
//   - M87*:    M ~ 6.5 × 10^9 M_sun, a* ~ 0.90, d = 16.8 Mpc
//   - Sgr A*:  M ~ 4.0 × 10^6 M_sun, a* ~ 0.50, d = 8.2 kpc
//
// Key scaling relations:
//   r_g = GM/c² (gravitational radius)
//   r_s = 2GM/c² = 2r_g (Schwarzschild radius)
//   r_isco = f(a*) × r_g (ISCO depends on spin)
//   L_Edd = 4πGMm_p c/σ_T (Eddington luminosity)
//
// REFERENCES:
// - Event Horizon Telescope Collaboration (2019) ApJL 875, L1-L6
// - Thorne (1974) "Disk-Accretion onto a Black Hole" for spin limit

#pragma once

#include <cmath>
#include <cstdint>

namespace Sirius {

//==============================================================================
// Physical Constants
//==============================================================================
namespace PhysicalConstants {
    constexpr double G = 6.67430e-11;           // Gravitational constant [m³/(kg·s²)]
    constexpr double c = 2.99792458e8;          // Speed of light [m/s]
    constexpr double M_sun = 1.98892e30;        // Solar mass [kg]
    constexpr double sigma_T = 6.6524587e-29;   // Thomson cross-section [m²]
    constexpr double m_p = 1.67262192e-27;      // Proton mass [kg]
    constexpr double pc = 3.08567758e16;        // Parsec [m]
    constexpr double Mpc = 3.08567758e22;       // Megaparsec [m]
    constexpr double year = 3.15576e7;          // Year [s]
    constexpr double k_B = 1.380649e-23;        // Boltzmann constant [J/K]
    constexpr double h = 6.62607015e-34;        // Planck constant [J·s]
    constexpr double sigma_SB = 5.670374e-8;    // Stefan-Boltzmann [W/(m²·K⁴)]
}

//==============================================================================
// Supermassive Black Hole Parameters
//==============================================================================
struct SMBHParams {
    // =========================================================================
    // Primary Parameters (User-Specified)
    // =========================================================================
    double mass_solar = 1.0e6;          ///< Mass in solar masses [M_sun]
    float spin_parameter = 0.0f;        ///< Dimensionless spin a* = a/M ∈ [-0.998, 0.998]
    float inclination_deg = 45.0f;      ///< Observer inclination angle [degrees]
    float distance_Mpc = 10.0f;         ///< Distance in megaparsecs

    // =========================================================================
    // Derived Quantities (Computed via computeDerived())
    // =========================================================================
    double M_kg = 0.0;                  ///< Mass in kilograms
    double M_geometric = 0.0;           ///< GM/c² in meters (gravitational radius)
    double r_g = 0.0;                   ///< Gravitational radius GM/c² [m]
    double r_s = 0.0;                   ///< Schwarzschild radius 2GM/c² [m]
    double r_isco = 0.0;                ///< ISCO radius in meters
    double r_isco_M = 0.0;              ///< ISCO radius in units of M
    double r_horizon = 0.0;             ///< Event horizon radius in meters
    double r_horizon_M = 0.0;           ///< Event horizon radius in units of M
    double distance_m = 0.0;            ///< Distance in meters
    double L_Edd = 0.0;                 ///< Eddington luminosity [W]
    double t_g = 0.0;                   ///< Gravitational time unit GM/c³ [s]

    // =========================================================================
    // Invariants
    // =========================================================================
    // mass_solar ∈ [10^5, 10^11] for SMBH class
    // |spin_parameter| ≤ 0.998 (Thorne limit)
    // distance_Mpc > 0
    // inclination_deg ∈ [0, 90]

    /// @brief Compute all derived quantities from primary parameters
    void computeDerived() {
        using namespace PhysicalConstants;

        // Clamp spin to Thorne limit
        if (spin_parameter > 0.998f) spin_parameter = 0.998f;
        if (spin_parameter < -0.998f) spin_parameter = -0.998f;

        // Mass conversions
        M_kg = mass_solar * M_sun;
        M_geometric = G * M_kg / (c * c);  // GM/c² in meters
        r_g = M_geometric;
        r_s = 2.0 * M_geometric;

        // Distance
        distance_m = distance_Mpc * Mpc;

        // ISCO radius (depends on spin)
        r_isco_M = computeISCO(spin_parameter);
        r_isco = r_isco_M * M_geometric;

        // Event horizon radius
        // r_+ = M + sqrt(M² - a²) where a = a* × M
        float a = std::abs(spin_parameter);
        r_horizon_M = 1.0 + std::sqrt(1.0 - a * a);
        r_horizon = r_horizon_M * M_geometric;

        // Eddington luminosity: L_Edd = 4πGMm_p c / σ_T
        L_Edd = 4.0 * M_PI * G * M_kg * m_p * c / sigma_T;

        // Gravitational time unit: GM/c³
        t_g = M_geometric / c;
    }

    /// @brief Compute ISCO radius in units of M
    /// r_ISCO = M {3 + Z₂ ∓ √[(3 - Z₁)(3 + Z₁ + 2Z₂)]}
    /// where ± is for prograde/retrograde orbits
    static double computeISCO(float a_star) {
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

    /// @brief Compute inner horizon radius in units of M
    /// r_- = M - sqrt(M² - a²)
    static double computeInnerHorizon(float a_star) {
        float a = std::abs(a_star);
        if (a > 1.0f) a = 1.0f;
        return 1.0 - std::sqrt(1.0 - a * a);
    }

    /// @brief Compute ergosphere radius at given theta
    /// r_ergo = M + sqrt(M² - a²cos²θ)
    double computeErgosphereRadius(double theta_rad) const {
        double a = std::abs(static_cast<double>(spin_parameter));
        double cos_theta = std::cos(theta_rad);
        return M_geometric * (1.0 + std::sqrt(1.0 - a * a * cos_theta * cos_theta));
    }

    /// @brief Convert radius from M units to meters
    double radiusToMeters(double r_M) const {
        return r_M * M_geometric;
    }

    /// @brief Convert radius from meters to M units
    double radiusFromMeters(double r_m) const {
        return r_m / M_geometric;
    }

    /// @brief Get angular size of r_g as seen from observer [radians]
    double angularSizeOfRg() const {
        return M_geometric / distance_m;
    }

    /// @brief Get angular size of the shadow (approximately 5.2 r_g for Schwarzschild)
    double shadowAngularSize() const {
        // Shadow radius ≈ sqrt(27) × M for Schwarzschild
        // For Kerr it's approximately similar but depends on spin and inclination
        return std::sqrt(27.0) * angularSizeOfRg();
    }

    // =========================================================================
    // Preset Configurations for Known SMBHs
    // =========================================================================

    /// @brief M87* configuration (EHT 2019)
    static SMBHParams M87Star() {
        SMBHParams p;
        p.mass_solar = 6.5e9;
        p.spin_parameter = 0.90f;
        p.inclination_deg = 17.0f;
        p.distance_Mpc = 16.8f;
        p.computeDerived();
        return p;
    }

    /// @brief Sagittarius A* configuration
    static SMBHParams SgrAStar() {
        SMBHParams p;
        p.mass_solar = 4.0e6;
        p.spin_parameter = 0.50f;
        p.inclination_deg = 45.0f;
        p.distance_Mpc = 0.0082f;  // 8.2 kpc
        p.computeDerived();
        return p;
    }

    /// @brief Gargantua-like configuration (Interstellar)
    static SMBHParams Gargantua() {
        SMBHParams p;
        p.mass_solar = 1.0e8;
        p.spin_parameter = 0.9999f;  // Near-extremal
        p.inclination_deg = 85.0f;   // Nearly edge-on
        p.distance_Mpc = 0.01f;      // Close for visualization
        p.computeDerived();
        return p;
    }
};

} // namespace Sirius
