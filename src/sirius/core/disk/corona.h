#pragma once

// Inverse-Compton corona model for the hot, optically thin plasma above an
// accretion disk: geometry containment, power-law emissivity, ray optical
// depth, and Comptonised scattered intensity.
// Ported from PHCR001A.h.
// Reference: Haardt & Maraschi (1991) ApJ 380, L51; Dove et al. (1997) ApJ
// 487, 759; Zdziarski & Gierlinski (2004) PTPS 155, 99.

#include "sirius/base/contracts.h"

#include <algorithm>
#include <cmath>
#include <cstdint>

namespace sirius::core {

// Corona geometry family.
enum class CoronaGeometry : uint32_t {
    Slab = 0,      // Thin slab above the disk (sandwich geometry).
    Lamppost = 1,  // Compact on-axis source (jet base).
    Sphere = 2,    // Spherical cloud around the black hole.
    Extended = 3   // Extended corona following the disk shape.
};

// Corona configuration and derived spectral quantities.
struct CoronaConfig {
    float temperature_keV = 100.0f;  // Electron temperature [keV].
    float optical_depth = 0.5f;      // Thomson optical depth tau.
    float scale_height_M = 5.0f;     // Vertical extent [M].
    float inner_radius_M = 0.0f;     // Inner boundary [M] (0 = use ISCO).
    float outer_radius_M = 20.0f;    // Outer boundary [M].
    CoronaGeometry geometry = CoronaGeometry::Extended;
    float lamppost_height_M = 10.0f;  // Height for lamppost geometry [M].
    float emissivity_index = 3.0f;    // Emissivity ~ r^(-q) index.
    float intensity_scale = 1.0f;     // Overall intensity multiplier.
    bool enabled = true;

    static constexpr float kMeKeV = 511.0f;  // Electron rest mass in keV.

    // Physical validity ranges (enforced by clamping in Validate, not asserted):
    //   temperature_keV in [10, 500], optical_depth in [0.1, 5.0],
    //   inner_radius_M >= r_isco(a), scale_height_M > 0.

    // Clamp parameters into their valid ranges.
    void Validate() {
        temperature_keV = std::clamp(temperature_keV, 10.0f, 500.0f);
        optical_depth = std::clamp(optical_depth, 0.1f, 5.0f);
        scale_height_M = std::max(scale_height_M, 0.1f);
        inner_radius_M = std::max(inner_radius_M, 0.0f);
        outer_radius_M = std::max(outer_radius_M, inner_radius_M + 1.0f);
        emissivity_index = std::clamp(emissivity_index, 0.0f, 10.0f);
        intensity_scale = std::max(intensity_scale, 0.0f);
    }

    // Comptonisation parameter y = 4 k T_e / (m_e c^2) * max(tau, tau^2).
    float ComptonizationParameter() const {
        float theta_e = temperature_keV / kMeKeV;  // Dimensionless temperature.
        float tau_eff = std::max(optical_depth, optical_depth * optical_depth);
        return 4.0f * theta_e * tau_eff;
    }

    // Thermal Comptonisation spectral index
    // Gamma ~ sqrt(9/4 + 3/(y (1 + tau/3))) - 1/2.
    float SpectralIndex() const {
        float y = ComptonizationParameter();
        float tau_term = 1.0f + optical_depth / 3.0f;
        float inner = 9.0f / 4.0f + 3.0f / (y * tau_term + 0.01f);
        return std::sqrt(inner) - 0.5f;
    }
};

// Corona physics (CPU implementation).
namespace corona_physics {

// True when the point (r, theta, phi) lies within the corona volume;
// isco is the ISCO radius [M] used as the inner boundary when inner_radius_M is 0.
inline bool IsInsideCorona(float r, float theta, [[maybe_unused]] float phi,
                           const CoronaConfig& config, float isco) {
    if (!config.enabled) return false;

    float inner = (config.inner_radius_M > 0) ? config.inner_radius_M : isco;
    if (r < inner || r > config.outer_radius_M) return false;

    // Vertical extent depends on geometry.
    float z = r * std::cos(theta);
    float rho = r * std::sin(theta);  // Cylindrical radius.

    switch (config.geometry) {
        case CoronaGeometry::Slab: {
            // Thin slab: |z| < scale_height.
            return std::abs(z) < config.scale_height_M;
        }
        case CoronaGeometry::Lamppost: {
            // Compact source at height h on the axis.
            float dist = std::sqrt(rho * rho +
                                   (z - config.lamppost_height_M) * (z - config.lamppost_height_M));
            return dist < config.scale_height_M;
        }
        case CoronaGeometry::Sphere: {
            // Spherical shell around the black hole.
            return r < config.outer_radius_M;
        }
        case CoronaGeometry::Extended: {
            // Extended corona: |z| < H(r) where H scales with r.
            float H_local = config.scale_height_M * std::sqrt(r / inner);
            return std::abs(z) < H_local;
        }
        default:
            SIRIUS_PRE(false);
            return false;
    }
}

// Corona emissivity at (r, theta) in arbitrary normalised units.
inline float Emissivity(float r, float theta, const CoronaConfig& config, float isco) {
    if (!config.enabled) return 0.0f;

    float inner = (config.inner_radius_M > 0) ? config.inner_radius_M : isco;
    if (r < inner || r > config.outer_radius_M) return 0.0f;

    // Power-law emissivity ~ r^(-q).
    float emiss = std::pow(inner / r, config.emissivity_index);

    // Vertical (Gaussian) falloff.
    float z = r * std::cos(theta);
    [[maybe_unused]] float rho = r * std::sin(theta);

    float H_local = config.scale_height_M;
    if (config.geometry == CoronaGeometry::Extended) {
        H_local *= std::sqrt(r / inner);
    }

    float z_factor = std::exp(-0.5f * (z * z) / (H_local * H_local));

    return config.intensity_scale * emiss * z_factor;
}

// Corona optical depth tau along the ray segment from (r1, theta1, phi1) to
// (r2, theta2, phi2) by trapezoidal integration over num_samples samples.
inline float OpticalDepthAlongRay(float r1, float theta1, float phi1, float r2, float theta2,
                                  float phi2, const CoronaConfig& config, float isco,
                                  int num_samples = 16) {
    if (!config.enabled) return 0.0f;

    float tau = 0.0f;
    float ds = 1.0f / static_cast<float>(num_samples);

    // Convert to Cartesian for the path length.
    auto to_cart = [](float r, float th, float ph, float& x, float& y, float& z) {
        float sin_th = std::sin(th);
        x = r * sin_th * std::cos(ph);
        y = r * sin_th * std::sin(ph);
        z = r * std::cos(th);
    };

    float x1, y1, z1, x2, y2, z2;
    to_cart(r1, theta1, phi1, x1, y1, z1);
    to_cart(r2, theta2, phi2, x2, y2, z2);

    float path_length =
        std::sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1) + (z2 - z1) * (z2 - z1));

    for (int i = 0; i < num_samples; ++i) {
        float t = (static_cast<float>(i) + 0.5f) * ds;
        float r = r1 + t * (r2 - r1);
        float theta = theta1 + t * (theta2 - theta1);
        float phi = phi1 + t * (phi2 - phi1);

        if (IsInsideCorona(r, theta, phi, config, isco)) {
            // Local opacity ~ optical_depth / scale_height (normalised).
            float opacity = config.optical_depth / config.scale_height_M;
            tau += opacity * path_length * ds;
        }
    }

    return tau;
}

// Scattered intensity from the corona via a simplified inverse-Compton model.
inline float ScatteredIntensity(float incident_intensity, float tau, const CoronaConfig& config) {
    if (!config.enabled || tau < 1e-6f) return 0.0f;

    // Fraction scattered: 1 - exp(-tau).
    float scattered_fraction = 1.0f - std::exp(-tau);

    // Average energy boost from Comptonisation.
    float y = config.ComptonizationParameter();
    float energy_boost = (y < 1.0f) ? (1.0f + y) : std::exp(y);

    return incident_intensity * scattered_fraction * energy_boost * config.intensity_scale;
}

// Local source function used by live volumetric transfer before the optical
// depth recurrence is applied. Separating the source from the scattered
// fraction avoids applying (1-exp(-tau)) twice in a ray marcher.
inline float ComptonizedSource(float seed_intensity, const CoronaConfig& config) {
    if (!config.enabled || seed_intensity <= 0.0f) return 0.0f;
    const float y = config.ComptonizationParameter();
    const float energy_boost = (y < 1.0f) ? (1.0f + y) : std::exp(y);
    return seed_intensity * energy_boost;
}

}  // namespace corona_physics

}  // namespace sirius::core
