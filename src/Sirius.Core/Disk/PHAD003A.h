// PHAD003A.h - Enhanced Volumetric Accretion Disk
// Component ID: PHAD003A (Physics/Accretion/VolumetricDisk)
//
// Combines Shakura-Sunyaev thin disk, Novikov-Thorne relativistic corrections,
// turbulent density perturbations, and corona emission into a unified
// volumetric accretion disk model.
//
// PHYSICAL FOUNDATION:
// ====================
// This model extends the standard thin disk with:
//   1. Vertical structure: ρ(r,z) = ρ_0(r) × exp(-z²/(2H²))
//   2. Turbulent perturbations: ρ' = ρ × (1 + δρ) from Kolmogorov cascade
//   3. Relativistic temperature: T(r) from Novikov-Thorne with Page & Thorne corrections
//   4. Corona emission: Inverse-Compton scattered X-rays
//
// The combined emissivity is:
//   j_ν = j_thermal(T(r)) × ρ'(r,z) + j_corona(r,z)
//
// REFERENCES:
// - Shakura & Sunyaev (1973) A&A 24, 337
// - Page & Thorne (1974) ApJ 191, 499
// - Novikov & Thorne (1973) in "Black Holes"
// - McKinney et al. (2012) MNRAS 423, 3083 (GRMHD)

#pragma once

#include "PHAD001A.h"  // Base AccretionDiskD
#include "PHTR001A.h"  // Turbulence
#include "PHCR001A.h"  // Corona
#include <cmath>
#include <algorithm>

namespace Sirius {

//==============================================================================
// Volumetric Disk Configuration
//==============================================================================
struct VolumetricDiskConfig {
    // =========================================================================
    // Disk Structure
    // =========================================================================
    float inner_radius_M = 0.0f;        ///< Inner edge (0 = use ISCO)
    float outer_radius_M = 500.0f;      ///< Outer edge in M
    float reference_radius_M = 10.0f;   ///< Reference radius for H/r
    float H_over_r = 0.1f;              ///< Scale height ratio at reference
    float H_power = 0.25f;              ///< Flaring index: H/r ∝ r^H_power
    float density_power = 1.5f;         ///< Midplane density ∝ r^(-density_power)

    // =========================================================================
    // Temperature Model
    // =========================================================================
    bool use_novikov_thorne = true;     ///< Use relativistic temperature profile
    float inner_temperature_K = 1e7f;   ///< Temperature at inner edge [K]
    float temperature_power = 0.75f;    ///< T ∝ r^(-temperature_power) for SS model

    // =========================================================================
    // Emission Properties
    // =========================================================================
    float emission_coefficient = 1.0f;  ///< Emission scaling factor
    float absorption_coefficient = 0.1f; ///< Absorption scaling factor
    float scattering_albedo = 0.3f;     ///< Single-scatter albedo ω
    float beaming_exponent = 3.0f;      ///< Doppler beaming power

    // =========================================================================
    // Sub-models
    // =========================================================================
    TurbulenceConfig turbulence;
    CoronaConfig corona;

    // =========================================================================
    // Ray Marching
    // =========================================================================
    int volumetric_samples = 64;        ///< Samples per ray through disk
    float optical_depth_max = 10.0f;    ///< Stop if τ exceeds this
    bool enabled = true;

    /// @brief Validate all parameters
    void validate() {
        inner_radius_M = std::max(inner_radius_M, 0.0f);
        outer_radius_M = std::max(outer_radius_M, inner_radius_M + 1.0f);
        reference_radius_M = std::clamp(reference_radius_M, inner_radius_M + 0.1f, outer_radius_M);
        H_over_r = std::clamp(H_over_r, 0.001f, 0.5f);
        H_power = std::clamp(H_power, 0.0f, 1.0f);
        density_power = std::clamp(density_power, 0.0f, 3.0f);
        temperature_power = std::clamp(temperature_power, 0.5f, 1.0f);
        volumetric_samples = std::clamp(volumetric_samples, 8, 256);
        turbulence.validate();
        corona.validate();
    }
};

//==============================================================================
// Enhanced Volumetric Disk Class
//==============================================================================
class VolumetricDisk {
public:
    explicit VolumetricDisk(const VolumetricDiskConfig& config = VolumetricDiskConfig())
        : m_Config(config) {
        m_Config.validate();
    }

    // =========================================================================
    // Geometry
    // =========================================================================

    /// @brief Compute scale height H at radius r
    /// H(r) = H_0 × (r/r_ref)^(1 + H_power) where H_0 = (H/r)_ref × r_ref
    float scaleHeight(float r) const {
        if (r <= 0.0f) return 0.0f;
        float H_ref = m_Config.H_over_r * m_Config.reference_radius_M;
        return H_ref * std::pow(r / m_Config.reference_radius_M, 1.0f + m_Config.H_power);
    }

    /// @brief Check if point is within disk volume
    bool isInDisk(float r, float theta, float isco) const {
        if (!m_Config.enabled) return false;

        float inner = (m_Config.inner_radius_M > 0) ? m_Config.inner_radius_M : isco;
        if (r < inner || r > m_Config.outer_radius_M) return false;

        // Vertical check: |z| < n × H for some n-sigma threshold
        float z = r * std::cos(theta);
        float rho_cyl = r * std::sin(theta);
        float H = scaleHeight(rho_cyl);

        return std::abs(z) < 4.0f * H;  // 4-sigma encompasses 99.99%
    }

    // =========================================================================
    // Density
    // =========================================================================

    /// @brief Midplane density at radius r (without turbulence)
    /// ρ_0(r) ∝ r^(-density_power)
    float midplaneDensity(float r, float isco) const {
        float inner = (m_Config.inner_radius_M > 0) ? m_Config.inner_radius_M : isco;
        if (r < inner) return 0.0f;

        return std::pow(inner / r, m_Config.density_power);
    }

    /// @brief Full density at position including vertical structure and turbulence
    /// ρ(r,θ,φ) = ρ_0(r) × exp(-z²/(2H²)) × (1 + δρ_turb)
    float density(float r, float theta, float phi, float isco) const {
        if (!isInDisk(r, theta, isco)) return 0.0f;

        float z = r * std::cos(theta);
        float rho_cyl = r * std::sin(theta);

        // Midplane density
        float rho_0 = midplaneDensity(rho_cyl, isco);

        // Vertical Gaussian falloff
        float H = scaleHeight(rho_cyl);
        float z_factor = std::exp(-0.5f * z * z / (H * H));

        // Turbulent perturbation
        float turb_factor = 1.0f;
        if (m_Config.turbulence.enabled) {
            turb_factor = TurbulenceNoise::sampleDensityPerturbation(
                r, theta, phi, m_Config.turbulence);
        }

        return rho_0 * z_factor * turb_factor;
    }

    // =========================================================================
    // Temperature
    // =========================================================================

    /// @brief Temperature at radius (midplane)
    /// Shakura-Sunyaev: T ∝ r^(-3/4)
    /// Novikov-Thorne: includes relativistic correction factor Q(r)
    float temperature(float r, float spin, float isco) const {
        float inner = (m_Config.inner_radius_M > 0) ? m_Config.inner_radius_M : isco;
        if (r < inner) return 0.0f;

        // Base temperature profile
        float T = m_Config.inner_temperature_K * std::pow(inner / r, m_Config.temperature_power);

        // Novikov-Thorne correction: T → 0 at ISCO
        if (m_Config.use_novikov_thorne) {
            float Q = novikovThorneCorrection(r, spin, isco);
            T *= std::pow(Q, 0.25f);  // T ∝ F^(1/4), Q ∝ F
        }

        return T;
    }

    /// @brief Novikov-Thorne correction factor Q(r)
    /// Q(r) = 1 - sqrt(r_isco/r) for simplified form
    /// Full form involves relativistic integrals
    float novikovThorneCorrection(float r, float spin, float isco) const {
        if (r <= isco) return 0.0f;

        // Simplified no-torque inner boundary condition
        float sqrt_ratio = std::sqrt(isco / r);
        float Q = 1.0f - sqrt_ratio;

        // Additional spin-dependent correction (approximate)
        float a = std::abs(spin);
        float spin_correction = 1.0f + 0.2f * a * std::pow(isco / r, 1.5f);

        return std::max(Q * spin_correction, 0.0f);
    }

    // =========================================================================
    // Emission
    // =========================================================================

    /// @brief Thermal emissivity at position
    /// j_ν ∝ ρ × B_ν(T) where B_ν is Planck function
    float thermalEmissivity(float r, float theta, float phi, float spin, float isco) const {
        float rho = density(r, theta, phi, isco);
        if (rho < 1e-10f) return 0.0f;

        float T = temperature(r, spin, isco);

        // Simplified emissivity proportional to ρ × T^4 (Stefan-Boltzmann scaling)
        float emission = rho * std::pow(T / m_Config.inner_temperature_K, 4.0f);

        return emission * m_Config.emission_coefficient;
    }

    /// @brief Absorption coefficient at position
    /// κ_ν ∝ ρ (gray approximation)
    float absorptionCoeff(float r, float theta, float phi, float isco) const {
        float rho = density(r, theta, phi, isco);
        return rho * m_Config.absorption_coefficient;
    }

    /// @brief Combined emissivity including corona
    float totalEmissivity(float r, float theta, float phi, float spin, float isco) const {
        float j_thermal = thermalEmissivity(r, theta, phi, spin, isco);

        // Add corona emission if enabled
        float j_corona = 0.0f;
        if (m_Config.corona.enabled &&
            CoronaPhysics::isInsideCorona(r, theta, phi, m_Config.corona, isco)) {
            j_corona = CoronaPhysics::emissivity(r, theta, m_Config.corona, isco);
        }

        return j_thermal + j_corona;
    }

    // =========================================================================
    // Ray Marching Integration
    // =========================================================================

    /// @brief Integrate emission along ray through disk
    /// Uses radiative transfer: dI/ds = j - κI
    struct RayMarchResult {
        float intensity;
        float optical_depth;
        float temperature_avg;
        bool hit_disk;
    };

    RayMarchResult integrateRay(float r_start, float theta_start, float phi_start,
                                 float r_end, float theta_end, float phi_end,
                                 float spin, float isco) const {
        RayMarchResult result = {0.0f, 0.0f, 0.0f, false};

        if (!m_Config.enabled) return result;

        int N = m_Config.volumetric_samples;
        float ds = 1.0f / static_cast<float>(N);

        float total_T = 0.0f;
        int T_samples = 0;

        // Ray marching with emission-absorption
        for (int i = 0; i < N; ++i) {
            float t = (static_cast<float>(i) + 0.5f) * ds;

            float r = r_start + t * (r_end - r_start);
            float theta = theta_start + t * (theta_end - theta_start);
            float phi = phi_start + t * (phi_end - phi_start);

            if (!isInDisk(r, theta, isco)) continue;
            result.hit_disk = true;

            // Path length element (approximate in spherical coords)
            float path_element = std::abs(r_end - r_start) * ds;

            // Emission and absorption
            float j = totalEmissivity(r, theta, phi, spin, isco);
            float kappa = absorptionCoeff(r, theta, phi, isco);

            // Optical depth increment
            float dtau = kappa * path_element;
            result.optical_depth += dtau;

            // Radiative transfer: I += j × ds × exp(-τ)
            float transmission = std::exp(-result.optical_depth);
            result.intensity += j * path_element * transmission;

            // Track average temperature
            float T = temperature(r, spin, isco);
            if (T > 0.0f) {
                total_T += T;
                T_samples++;
            }

            // Early termination if optically thick
            if (result.optical_depth > m_Config.optical_depth_max) break;
        }

        if (T_samples > 0) {
            result.temperature_avg = total_T / static_cast<float>(T_samples);
        }

        return result;
    }

    // =========================================================================
    // Configuration Access
    // =========================================================================
    const VolumetricDiskConfig& config() const { return m_Config; }
    void setConfig(const VolumetricDiskConfig& config) {
        m_Config = config;
        m_Config.validate();
    }

private:
    VolumetricDiskConfig m_Config;
};

} // namespace Sirius
