// PHCR001A.h - Corona Model for Accretion Disks
// Component ID: PHCR001A (Physics/Corona/CoronaModel)
//
// Inverse-Compton scattering corona model for X-ray emission regions
// surrounding black hole accretion disks.
//
// PHYSICAL FOUNDATION:
// ====================
// The corona is a hot, optically thin plasma above the accretion disk that
// produces hard X-rays via inverse-Compton scattering of soft disk photons.
//
// Key physics:
//   - Electron temperature: T_e ~ 10^8-10^9 K (10-100 keV)
//   - Thomson optical depth: τ ~ 0.1-2
//   - Geometry: Compact region near disk, often "lamppost" or "slab"
//   - Emission: Power-law X-ray spectrum with cutoff
//
// The Comptonization parameter y = 4kT_e/(m_e c²) × max(τ, τ²) determines
// the spectral hardening.
//
// REFERENCES:
// - Haardt & Maraschi (1991) ApJ 380, L51 (two-phase model)
// - Dove et al. (1997) ApJ 487, 759 (corona geometry)
// - Zdziarski & Gierliński (2004) Progress Theoretical Physics Supplement 155, 99

#pragma once

#include <cmath>
#include <cstdint>
#include <algorithm>

namespace Sirius {

//==============================================================================
// Corona Geometry Type
//==============================================================================
enum class CoronaGeometry : uint32_t {
    Slab = 0,       ///< Thin slab above disk (sandwich geometry)
    Lamppost = 1,   ///< Compact source on axis (jet base)
    Sphere = 2,     ///< Spherical cloud around black hole
    Extended = 3    ///< Extended corona following disk shape
};

//==============================================================================
// Corona Configuration
//==============================================================================
struct CoronaConfig {
    float temperature_keV = 100.0f;     ///< Electron temperature [keV]
    float optical_depth = 0.5f;         ///< Thomson optical depth τ
    float scale_height_M = 5.0f;        ///< Vertical extent in M
    float inner_radius_M = 0.0f;        ///< Corona inner boundary (0 = use ISCO)
    float outer_radius_M = 20.0f;       ///< Corona outer boundary in M
    CoronaGeometry geometry = CoronaGeometry::Extended;
    float lamppost_height_M = 10.0f;    ///< Height for lamppost geometry
    float emissivity_index = 3.0f;      ///< Emissivity ∝ r^(-q) index
    float intensity_scale = 1.0f;       ///< Overall intensity multiplier
    bool enabled = true;

    // =========================================================================
    // Physical Constants
    // =========================================================================
    static constexpr float m_e_keV = 511.0f;  // Electron rest mass in keV

    // =========================================================================
    // Invariants
    // =========================================================================
    // temperature_keV ∈ [10, 500]
    // optical_depth ∈ [0.1, 5.0]
    // inner_radius_M >= r_isco(a)
    // scale_height_M > 0

    /// @brief Validate and clamp parameters
    void validate() {
        temperature_keV = std::clamp(temperature_keV, 10.0f, 500.0f);
        optical_depth = std::clamp(optical_depth, 0.1f, 5.0f);
        scale_height_M = std::max(scale_height_M, 0.1f);
        inner_radius_M = std::max(inner_radius_M, 0.0f);
        outer_radius_M = std::max(outer_radius_M, inner_radius_M + 1.0f);
        emissivity_index = std::clamp(emissivity_index, 0.0f, 10.0f);
        intensity_scale = std::max(intensity_scale, 0.0f);
    }

    /// @brief Compute Comptonization parameter y
    /// y = 4kT_e/(m_e c²) × max(τ, τ²)
    float comptonizationParameter() const {
        float theta_e = temperature_keV / m_e_keV;  // Dimensionless temperature
        float tau_eff = std::max(optical_depth, optical_depth * optical_depth);
        return 4.0f * theta_e * tau_eff;
    }

    /// @brief Compute spectral index Γ from Comptonization
    /// For thermal Comptonization: Γ ≈ sqrt(9/4 + 3/(y × (1 + τ/3))) - 1/2
    float spectralIndex() const {
        float y = comptonizationParameter();
        float tau_term = 1.0f + optical_depth / 3.0f;
        float inner = 9.0f / 4.0f + 3.0f / (y * tau_term + 0.01f);
        return std::sqrt(inner) - 0.5f;
    }
};

//==============================================================================
// Corona Physics Functions (CPU Implementation)
//==============================================================================

namespace CoronaPhysics {

/// @brief Check if point is inside corona volume
/// @param r Radius in M
/// @param theta Polar angle [rad]
/// @param phi Azimuthal angle [rad]
/// @param config Corona configuration
/// @param isco ISCO radius in M (for inner boundary)
/// @return true if point is within corona
inline bool isInsideCorona(float r, float theta, float phi,
                           const CoronaConfig& config, float isco) {
    if (!config.enabled) return false;

    float inner = (config.inner_radius_M > 0) ? config.inner_radius_M : isco;
    if (r < inner || r > config.outer_radius_M) return false;

    // Vertical extent depends on geometry
    float z = r * std::cos(theta);
    float rho = r * std::sin(theta);  // Cylindrical radius

    switch (config.geometry) {
        case CoronaGeometry::Slab: {
            // Thin slab: |z| < scale_height
            return std::abs(z) < config.scale_height_M;
        }
        case CoronaGeometry::Lamppost: {
            // Compact source at height h on axis
            float dist = std::sqrt(rho * rho + (z - config.lamppost_height_M) * (z - config.lamppost_height_M));
            return dist < config.scale_height_M;
        }
        case CoronaGeometry::Sphere: {
            // Spherical shell around black hole
            return r < config.outer_radius_M;
        }
        case CoronaGeometry::Extended:
        default: {
            // Extended corona: |z| < H(r) where H scales with r
            float H_local = config.scale_height_M * std::sqrt(r / inner);
            return std::abs(z) < H_local;
        }
    }
}

/// @brief Compute corona emissivity at position
/// @param r Radius in M
/// @param theta Polar angle [rad]
/// @param config Corona configuration
/// @param isco ISCO radius in M
/// @return Emissivity (arbitrary units, normalized)
inline float emissivity(float r, float theta, const CoronaConfig& config, float isco) {
    if (!config.enabled) return 0.0f;

    float inner = (config.inner_radius_M > 0) ? config.inner_radius_M : isco;
    if (r < inner || r > config.outer_radius_M) return 0.0f;

    // Power-law emissivity ∝ r^(-q)
    float emiss = std::pow(inner / r, config.emissivity_index);

    // Vertical falloff (Gaussian)
    float z = r * std::cos(theta);
    float rho = r * std::sin(theta);

    float H_local = config.scale_height_M;
    if (config.geometry == CoronaGeometry::Extended) {
        H_local *= std::sqrt(r / inner);
    }

    float z_factor = std::exp(-0.5f * (z * z) / (H_local * H_local));

    return config.intensity_scale * emiss * z_factor;
}

/// @brief Compute corona optical depth along ray segment
/// @param r1, theta1, phi1 Start position
/// @param r2, theta2, phi2 End position
/// @param config Corona configuration
/// @param isco ISCO radius
/// @param num_samples Number of integration samples
/// @return Optical depth τ
inline float opticalDepthAlongRay(float r1, float theta1, float phi1,
                                   float r2, float theta2, float phi2,
                                   const CoronaConfig& config, float isco,
                                   int num_samples = 16) {
    if (!config.enabled) return 0.0f;

    // Simple trapezoidal integration
    float tau = 0.0f;
    float ds = 1.0f / static_cast<float>(num_samples);

    // Convert to Cartesian for path length
    auto toCart = [](float r, float th, float ph, float& x, float& y, float& z) {
        float sin_th = std::sin(th);
        x = r * sin_th * std::cos(ph);
        y = r * sin_th * std::sin(ph);
        z = r * std::cos(th);
    };

    float x1, y1, z1, x2, y2, z2;
    toCart(r1, theta1, phi1, x1, y1, z1);
    toCart(r2, theta2, phi2, x2, y2, z2);

    float path_length = std::sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1));

    for (int i = 0; i < num_samples; ++i) {
        float t = (static_cast<float>(i) + 0.5f) * ds;
        float r = r1 + t * (r2 - r1);
        float theta = theta1 + t * (theta2 - theta1);
        float phi = phi1 + t * (phi2 - phi1);

        if (isInsideCorona(r, theta, phi, config, isco)) {
            // Local opacity ∝ optical_depth / scale_height (normalized)
            float opacity = config.optical_depth / config.scale_height_M;
            tau += opacity * path_length * ds;
        }
    }

    return tau;
}

/// @brief Compute scattered intensity from corona
/// Uses simplified inverse-Compton model
/// @param incident_intensity Input photon intensity
/// @param tau Optical depth traversed
/// @param config Corona configuration
/// @return Scattered intensity
inline float scatteredIntensity(float incident_intensity, float tau,
                                 const CoronaConfig& config) {
    if (!config.enabled || tau < 1e-6f) return 0.0f;

    // Fraction scattered: 1 - exp(-τ)
    float scattered_fraction = 1.0f - std::exp(-tau);

    // Average energy boost from Comptonization
    float y = config.comptonizationParameter();
    float energy_boost = (y < 1.0f) ? (1.0f + y) : std::exp(y);

    return incident_intensity * scattered_fraction * energy_boost * config.intensity_scale;
}

} // namespace CoronaPhysics

} // namespace Sirius
