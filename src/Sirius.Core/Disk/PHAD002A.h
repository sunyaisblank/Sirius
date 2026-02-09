// PHAD002A.h - Volumetric Accretion Disk Physics
// Component ID: PHAD002A
// Purpose: 3D disk model with vertical structure for Kerr black holes
//
// MATHEMATICAL BASIS:
// Extends the Novikov-Thorne thin disk (PHAD001A) to finite thickness using
// Shakura-Sunyaev α-disk vertical structure:
//
//   H(r) / r = (c_s / v_K) × (r / R_g)^ε
//
// Density profile (hydrostatic equilibrium):
//   ρ(r, z) = ρ_0(r) × exp(-z² / (2H(r)²))
//
// Radiative transfer along ray:
//   dI/ds = -κ ρ I + j = -κ ρ (I - S)
//
// where S = B(T) is the source function (LTE blackbody).
//
// REFERENCES:
// - Shakura & Sunyaev (1973) "Black Holes in Binary Systems"
// - Hubeny et al. (2001) "Vertical Structure of Accretion Disks"
// - Schnittman & Krolik (2009) "X-ray Polarization from Black Hole Accretion Disks"
//
// TESTS: TSDK002A.cpp

#pragma once

#include "PHAD001A.h"
#include <cmath>
#include <algorithm>

namespace Sirius {

// =============================================================================
// VolumetricDiskD: 3D disk with vertical structure
// =============================================================================

class VolumetricDiskD {
public:
    struct Config {
        // =====================================================================
        // Base parameters (inherited from thin disk model)
        // =====================================================================
        double M = 1.0;              // Black hole mass [M_sun]
        double a_star = 0.0;         // Dimensionless spin a/M ∈ [-1, 1]
        double Mdot = 1e-8;          // Accretion rate [M_sun/year]
        double r_inner = 0;          // Inner edge (0 = ISCO, auto-computed)
        double r_outer = 500;        // Outer edge [GM/c²]

        // =====================================================================
        // Vertical structure parameters
        // =====================================================================
        // Scale height at reference radius: H(r) = H_over_r × r × (r/r_ref)^H_power
        double H_over_r = 0.1;       // Scale height ratio H/r at r_ref
        double r_ref = 10.0;         // Reference radius for H/r [GM/c²]
        double H_power = 0.25;       // Flaring index: H/r ∝ r^H_power
                                     // H_power = 0: constant H/r (flat disk)
                                     // H_power > 0: flared disk (H/r increases outward)
                                     // H_power = 0.25: typical for radiation-supported disks

        // =====================================================================
        // Opacity parameters
        // =====================================================================
        // Midplane optical depth: τ_0(r) = tau_midplane × (r/r_ref)^tau_power
        double tau_midplane = 10.0;  // Vertical optical depth at r_ref
        double tau_power = -1.5;     // τ_0 ∝ r^tau_power (surface density falloff)
                                     // tau_power = -1.5: typical Shakura-Sunyaev

        // =====================================================================
        // Vertical temperature structure
        // =====================================================================
        // T(r, z) = T_mid(r) × [1 - (z/H)² × (1 - T_atm_ratio²)]^0.25
        // This gives T_mid at midplane and T_atm = T_atm_ratio × T_mid at surface
        double T_atm_ratio = 0.8;    // Atmosphere/midplane temperature ratio
                                     // 1.0 = isothermal, < 1.0 = cooler atmosphere

        // =====================================================================
        // Truncation
        // =====================================================================
        double z_truncation = 3.0;   // Disk extends to z_truncation × H(r)
                                     // 3.0 = 3σ truncation (contains 99.7% of mass)
    };

    // =========================================================================
    // Construction
    // =========================================================================

    VolumetricDiskD() : m_config() { init(); }

    explicit VolumetricDiskD(const Config& config)
        : m_config(config) { init(); }

private:
    void init() {
        // Compute ISCO radius
        m_r_isco = AccretionDiskD::computeISCO(m_config.a_star);

        // Set inner edge
        m_r_inner = (m_config.r_inner > 0) ? m_config.r_inner : m_r_isco;
        m_r_outer = m_config.r_outer;

        // Create underlying thin disk for temperature/flux calculations
        AccretionDiskD::Config thinConfig;
        thinConfig.M = m_config.M;
        thinConfig.a_star = m_config.a_star;
        thinConfig.Mdot = m_config.Mdot;
        thinConfig.r_inner = m_r_inner;
        thinConfig.r_outer = m_r_outer;
        m_thinDisk = AccretionDiskD(thinConfig);

        // Precompute normalization for density
        // ρ_0(r) such that vertical integral gives surface density Σ(r)
        // Σ = ∫ρ dz = ρ_0 × √(2π) × H
        // We normalize so that τ_vertical = κ × Σ = tau_midplane at r_ref
        // This gives κ × ρ_0 × √(2π) × H = tau_midplane
        // So: κ × ρ_0 = tau_midplane / (√(2π) × H)
        m_kappa_rho0_ref = m_config.tau_midplane /
                          (std::sqrt(2.0 * M_PI) * scaleHeight(m_config.r_ref));
    }

public:

    // =========================================================================
    // Scale Height
    // =========================================================================
    // H(r) = (H/r) × r = H_over_r × r × (r/r_ref)^H_power
    //
    // For α-disk: H/r ~ c_s/v_K ~ √(kT/μm_p) / √(GM/r)
    // With T ∝ r^(-3/4), this gives H/r ∝ r^(1/8) to r^(1/4)

    double scaleHeight(double r) const {
        if (r <= 0) return 0;

        double r_ratio = r / m_config.r_ref;
        double H_over_r = m_config.H_over_r * std::pow(r_ratio, m_config.H_power);

        // Clamp H/r to physical range [0.01, 0.5]
        H_over_r = std::clamp(H_over_r, 0.01, 0.5);

        return H_over_r * r;
    }

    // =========================================================================
    // Density Profile
    // =========================================================================
    // ρ(r, z) = ρ_0(r) × exp(-z² / (2H²))
    //
    // The midplane density ρ_0(r) is determined by the surface density Σ(r)
    // and scale height H(r): ρ_0 = Σ / (√(2π) × H)
    //
    // We parameterize via optical depth: κρ_0(r) ∝ τ_0(r) / H(r)

    double density(double r, double z) const {
        if (!isInDiskVolume(r, z)) return 0;

        double H = scaleHeight(r);
        if (H <= 0) return 0;

        // Gaussian vertical profile
        double z_over_H = z / H;
        double gaussian = std::exp(-0.5 * z_over_H * z_over_H);

        // Radial scaling of κρ_0
        double r_ratio = r / m_config.r_ref;
        double kappa_rho0 = m_kappa_rho0_ref * std::pow(r_ratio, m_config.tau_power) *
                           (scaleHeight(m_config.r_ref) / H);

        // Return ρ (we'll multiply by κ in opacity function)
        // For simplicity, return κρ directly since we don't separate κ and ρ
        return kappa_rho0 * gaussian;
    }

    // =========================================================================
    // Opacity × Density (combined for efficiency)
    // =========================================================================
    // Returns κ × ρ at position (r, z)
    // This is what we need for optical depth integration: dτ = κ ρ ds

    double opacityDensity(double r, double z) const {
        return density(r, z);  // density() already returns κρ
    }

    // =========================================================================
    // Temperature Profile
    // =========================================================================
    // T(r, z) = T_mid(r) × f(z/H)
    //
    // Vertical temperature structure:
    // - Midplane: T = T_mid (from thin disk model)
    // - Surface: T = T_atm = T_atm_ratio × T_mid
    //
    // Interpolation: T^4 ∝ 1 - (z/H)² × (1 - T_atm_ratio^4)
    // This approximates radiative equilibrium in plane-parallel atmosphere

    double temperature(double r, double z) const {
        if (!isInDiskVolume(r, z)) return 0;

        // Midplane temperature from thin disk model
        double T_mid = m_thinDisk.temperature(r);
        if (T_mid <= 0) return 0;

        // Vertical temperature profile
        double H = scaleHeight(r);
        if (H <= 0) return T_mid;

        double z_over_H = std::abs(z) / H;
        double z_over_H_sq = z_over_H * z_over_H;

        // Clamp to truncation boundary
        z_over_H_sq = std::min(z_over_H_sq, m_config.z_truncation * m_config.z_truncation);

        // T^4 interpolation
        double T_atm_ratio_4 = std::pow(m_config.T_atm_ratio, 4.0);
        double T4_factor = 1.0 - z_over_H_sq * (1.0 - T_atm_ratio_4) /
                          (m_config.z_truncation * m_config.z_truncation);

        // Clamp to avoid negative values
        T4_factor = std::max(T4_factor, T_atm_ratio_4);

        return T_mid * std::pow(T4_factor, 0.25);
    }

    // =========================================================================
    // Vertical Optical Depth
    // =========================================================================
    // τ(r, z) = ∫_z^∞ κρ dz' = τ_0(r) × [1 - erf(z / (√2 H))] / 2
    //
    // For Gaussian density: analytical solution via error function

    double verticalOpticalDepth(double r, double z) const {
        if (r < m_r_inner || r > m_r_outer) return 0;

        double H = scaleHeight(r);
        if (H <= 0) return 0;

        // Total vertical optical depth (midplane to infinity, one side)
        double r_ratio = r / m_config.r_ref;
        double tau_0_half = 0.5 * m_config.tau_midplane * std::pow(r_ratio, m_config.tau_power);

        // Complement of error function for z > 0
        double z_arg = std::abs(z) / (std::sqrt(2.0) * H);
        double erfc_val = std::erfc(z_arg);

        return tau_0_half * erfc_val;
    }

    // =========================================================================
    // Photosphere Height
    // =========================================================================
    // Find z where τ(r, z) = 2/3 (Eddington approximation)
    // Solve: erfc(z / (√2 H)) = (4/3) × τ_0^(-1)

    double photosphereHeight(double r) const {
        if (r < m_r_inner || r > m_r_outer) return 0;

        double H = scaleHeight(r);
        if (H <= 0) return 0;

        // Total vertical optical depth
        double r_ratio = r / m_config.r_ref;
        double tau_0_half = 0.5 * m_config.tau_midplane * std::pow(r_ratio, m_config.tau_power);

        if (tau_0_half < 2.0/3.0) {
            // Optically thin: no well-defined photosphere
            return 0;
        }

        // Solve erfc(x) = (2/3) / tau_0_half
        double target_erfc = (2.0/3.0) / tau_0_half;
        target_erfc = std::clamp(target_erfc, 0.0, 2.0);

        // Inverse erfc approximation (Newton's method would be more accurate)
        // For erfc(x) ≈ y, we have x ≈ erfcinv(y)
        // Use approximation: erfcinv(y) ≈ sqrt(-ln(y × sqrt(π)))  for small y
        double x;
        if (target_erfc < 0.1) {
            // High optical depth: use asymptotic
            x = std::sqrt(-std::log(target_erfc * std::sqrt(M_PI)));
        } else if (target_erfc > 1.9) {
            // Very low optical depth
            x = 0;
        } else {
            // Moderate: use rational approximation
            // erfcinv(1) = 0, erfcinv(0) = ∞
            // Linear interpolation as rough approximation
            x = 1.5 * (1.0 - target_erfc / 2.0);
        }

        return x * std::sqrt(2.0) * H;
    }

    // =========================================================================
    // Boundary Checks
    // =========================================================================

    /// @brief Check if point is within disk volume
    /// @param r Cylindrical radius [GM/c²]
    /// @param z Height above midplane [GM/c²]
    /// @return true if (r, z) is inside the disk
    bool isInDiskVolume(double r, double z) const {
        // Radial bounds
        if (r < m_r_inner || r > m_r_outer) return false;

        // Vertical bounds
        double H = scaleHeight(r);
        double z_max = m_config.z_truncation * H;

        return std::abs(z) <= z_max;
    }

    /// @brief Get disk surface height at radius r
    /// @param r Cylindrical radius [GM/c²]
    /// @return Maximum z coordinate of disk at this radius
    double surfaceHeight(double r) const {
        if (r < m_r_inner || r > m_r_outer) return 0;
        return m_config.z_truncation * scaleHeight(r);
    }

    // =========================================================================
    // Source Function (LTE)
    // =========================================================================
    // S(r, z) = B(T(r, z)) where B is Planck function
    // For gray opacity, this is just the blackbody intensity

    double sourceFunction(double r, double z) const {
        double T = temperature(r, z);
        if (T <= 0) return 0;

        // Stefan-Boltzmann: B = σT^4 / π (intensity, not flux)
        return Constants::Physical::SIGMA_SB * std::pow(T, 4.0) / M_PI;
    }

    // =========================================================================
    // Accessors
    // =========================================================================

    const Config& config() const { return m_config; }
    double iscoRadius() const { return m_r_isco; }
    double innerRadius() const { return m_r_inner; }
    double outerRadius() const { return m_r_outer; }

    /// @brief Get underlying thin disk model (for midplane properties)
    const AccretionDiskD& thinDisk() const { return m_thinDisk; }

private:
    Config m_config;

    // Derived quantities
    double m_r_isco;        // ISCO radius [GM/c²]
    double m_r_inner;       // Inner edge [GM/c²]
    double m_r_outer;       // Outer edge [GM/c²]
    double m_kappa_rho0_ref; // κρ_0 at reference radius

    // Underlying thin disk for temperature/flux
    AccretionDiskD m_thinDisk;
};

} // namespace Sirius
