// PHSF001A.h - Depth-Resolved Starfield
// Component ID: PHSF001A (Physics/Starfield/DepthStarfield)
//
// Parallax-correct stellar catalog with proper depth of field simulation
// for cinematically accurate starfield backgrounds.
//
// PHYSICAL FOUNDATION:
// ====================
// Stars at different distances exhibit parallax when the camera moves,
// creating depth perception in the starfield. This module generates
// procedural star catalogs with:
//   - 3D positions (distance-dependent parallax)
//   - Spectral types (color from B-V index)
//   - Proper magnitudes (apparent brightness)
//   - Depth of field blur (aperture-dependent)
//
// Magnitude system: m = -2.5 log₁₀(F/F₀)
// Apparent magnitude from absolute: m = M + 5 log₁₀(d/10pc)
//
// REFERENCES:
// - Hipparchos Catalog (ESA, 1997)
// - Tycho-2 Catalog (Høg et al., 2000)
// - CIE chromaticity for blackbody colors

#pragma once

#include <cmath>
#include <cstdint>
#include <vector>
#include <algorithm>
#include <random>

namespace Sirius {

//==============================================================================
// Star Entry Structure
//==============================================================================
struct StarEntry {
    float direction_x;         ///< Unit direction X (galactic coords)
    float direction_y;         ///< Unit direction Y
    float direction_z;         ///< Unit direction Z
    float distance_pc;         ///< Distance in parsecs
    float magnitude;           ///< Apparent magnitude
    float color_bv;            ///< B-V color index
    float temperature_K;       ///< Effective temperature [K]
    float padding;             ///< Alignment

    /// @brief Compute RGB color from B-V index
    /// Uses blackbody approximation: T ≈ 4600 / (0.92 + B-V) for B-V > 0
    void computeColor(float& r, float& g, float& b) const {
        // Temperature to RGB (simplified Planckian locus)
        float T = temperature_K;
        if (T <= 0.0f) T = 5778.0f;  // Default to solar

        // Normalize to 6500K white point
        float T_norm = T / 6500.0f;

        if (T <= 6500.0f) {
            // Cool stars: more red
            r = 1.0f;
            g = std::pow(T_norm, 0.4f);
            b = std::pow(T_norm, 0.8f);
        } else {
            // Hot stars: more blue
            r = std::pow(1.0f / T_norm, 0.3f);
            g = std::pow(1.0f / T_norm, 0.15f);
            b = 1.0f;
        }

        // Normalize
        float max_c = std::max({r, g, b});
        r /= max_c; g /= max_c; b /= max_c;
    }

    /// @brief Compute intensity from magnitude
    /// F ∝ 10^(-0.4 × m)
    float intensity() const {
        return std::pow(10.0f, -0.4f * magnitude);
    }
};

//==============================================================================
// Starfield Configuration
//==============================================================================
struct StarfieldConfig {
    uint32_t star_count = 100000;       ///< Number of stars in catalog
    float min_distance_pc = 1.0f;       ///< Nearest star distance [pc]
    float max_distance_pc = 10000.0f;   ///< Farthest star distance [pc]
    float magnitude_limit = 12.0f;      ///< Faintest visible magnitude
    float aperture_mm = 50.0f;          ///< Virtual aperture for DoF (0 = infinite)
    float focus_distance_pc = 100.0f;   ///< Focus plane distance [pc]
    float brightness_scale = 1.0f;      ///< Overall brightness multiplier
    uint32_t seed = 42;                 ///< Random seed
    bool enabled = true;
    bool enable_parallax = true;        ///< Enable parallax for camera motion
    bool enable_dof = true;             ///< Enable depth of field blur

    // =========================================================================
    // Star Distribution Parameters (Galactic model)
    // =========================================================================
    float thin_disk_scale_height_pc = 300.0f;   ///< Thin disk scale height [pc]
    float thick_disk_scale_height_pc = 1000.0f; ///< Thick disk scale height [pc]
    float disk_scale_length_pc = 2500.0f;       ///< Disk radial scale [pc]
    float thin_disk_fraction = 0.9f;            ///< Fraction of stars in thin disk

    // =========================================================================
    // Invariants
    // =========================================================================
    // star_count <= 10^7
    // min_distance_pc > 0
    // aperture_mm ∈ [0, 1000] (0 = infinite DoF)

    void validate() {
        star_count = std::min(star_count, 10000000u);
        min_distance_pc = std::max(min_distance_pc, 0.1f);
        max_distance_pc = std::max(max_distance_pc, min_distance_pc + 1.0f);
        magnitude_limit = std::clamp(magnitude_limit, 0.0f, 20.0f);
        aperture_mm = std::clamp(aperture_mm, 0.0f, 1000.0f);
        focus_distance_pc = std::clamp(focus_distance_pc, min_distance_pc, max_distance_pc);
    }
};

//==============================================================================
// Starfield Generator
//==============================================================================
class StarfieldGenerator {
public:
    explicit StarfieldGenerator(const StarfieldConfig& config = StarfieldConfig())
        : m_Config(config) {
        m_Config.validate();
    }

    /// @brief Generate procedural star catalog
    std::vector<StarEntry> generate() {
        std::vector<StarEntry> stars;
        stars.reserve(m_Config.star_count);

        std::mt19937 rng(m_Config.seed);
        std::uniform_real_distribution<float> uniform(0.0f, 1.0f);

        for (uint32_t i = 0; i < m_Config.star_count; ++i) {
            StarEntry star = generateStar(rng, uniform);
            if (star.magnitude <= m_Config.magnitude_limit) {
                stars.push_back(star);
            }
        }

        return stars;
    }

    /// @brief Sample starfield at direction with optional parallax offset
    /// @param dir_x, dir_y, dir_z View direction (unit vector)
    /// @param camera_offset_pc Camera offset from origin [pc] for parallax
    /// @param stars Star catalog
    /// @param r, g, b Output color (accumulated)
    void sampleStarfield(float dir_x, float dir_y, float dir_z,
                         float cam_x_pc, float cam_y_pc, float cam_z_pc,
                         const std::vector<StarEntry>& stars,
                         float& r, float& g, float& b) const {
        r = g = b = 0.0f;

        if (!m_Config.enabled || stars.empty()) return;

        for (const auto& star : stars) {
            // Parallax-adjusted direction
            float dx = star.direction_x * star.distance_pc - cam_x_pc;
            float dy = star.direction_y * star.distance_pc - cam_y_pc;
            float dz = star.direction_z * star.distance_pc - cam_z_pc;
            float d = std::sqrt(dx * dx + dy * dy + dz * dz);

            if (d < 1e-6f) continue;

            // Normalized parallax-adjusted direction
            float px = dx / d, py = dy / d, pz = dz / d;

            // Angular distance from view direction
            float cos_angle = dir_x * px + dir_y * py + dir_z * pz;

            // Point spread function (Gaussian approximation)
            float sigma = computePSFSigma(star.distance_pc);
            float angle = std::acos(std::clamp(cos_angle, -1.0f, 1.0f));
            float psf = std::exp(-0.5f * angle * angle / (sigma * sigma));

            if (psf < 1e-6f) continue;

            // Intensity and color
            float intensity = star.intensity() * psf * m_Config.brightness_scale;
            float sr, sg, sb;
            star.computeColor(sr, sg, sb);

            r += sr * intensity;
            g += sg * intensity;
            b += sb * intensity;
        }
    }

    const StarfieldConfig& config() const { return m_Config; }
    void setConfig(const StarfieldConfig& config) {
        m_Config = config;
        m_Config.validate();
    }

private:
    /// @brief Generate a single star with realistic properties
    StarEntry generateStar(std::mt19937& rng, std::uniform_real_distribution<float>& uniform) {
        StarEntry star;

        // Direction: uniform on sphere
        float theta = std::acos(2.0f * uniform(rng) - 1.0f);
        float phi = 2.0f * 3.14159f * uniform(rng);
        star.direction_x = std::sin(theta) * std::cos(phi);
        star.direction_y = std::sin(theta) * std::sin(phi);
        star.direction_z = std::cos(theta);

        // Distance: power-law weighted towards nearby
        float u = uniform(rng);
        float d_min = m_Config.min_distance_pc;
        float d_max = m_Config.max_distance_pc;
        // r² dr weighting (uniform in volume)
        star.distance_pc = std::cbrt(u * (d_max * d_max * d_max) +
                                      (1.0f - u) * (d_min * d_min * d_min));

        // Absolute magnitude: Salpeter IMF approximation
        // More low-luminosity stars
        float abs_mag = -2.0f + 12.0f * std::pow(uniform(rng), 0.3f);

        // Apparent magnitude: m = M + 5 log₁₀(d/10)
        star.magnitude = abs_mag + 5.0f * std::log10(star.distance_pc / 10.0f);

        // B-V color index: correlate with absolute magnitude
        // Hot (blue) stars: low M, low B-V
        // Cool (red) stars: high M, high B-V
        float t = (abs_mag + 2.0f) / 12.0f;
        star.color_bv = -0.3f + 2.0f * t + 0.3f * (uniform(rng) - 0.5f);
        star.color_bv = std::clamp(star.color_bv, -0.4f, 2.0f);

        // Temperature from B-V
        if (star.color_bv > -0.1f) {
            star.temperature_K = 4600.0f / (0.92f + star.color_bv);
        } else {
            star.temperature_K = 30000.0f + 10000.0f * (-star.color_bv - 0.1f);
        }
        star.temperature_K = std::clamp(star.temperature_K, 2500.0f, 50000.0f);

        star.padding = 0.0f;

        return star;
    }

    /// @brief Compute point spread function sigma (angular)
    /// Includes diffraction limit and depth of field
    float computePSFSigma(float distance_pc) const {
        // Base angular size (point source)
        float sigma_base = 1e-6f;  // ~0.2 arcsec

        // Depth of field blur
        float sigma_dof = 0.0f;
        if (m_Config.enable_dof && m_Config.aperture_mm > 0.0f) {
            // Circle of confusion: CoC ∝ A × |1/d - 1/d_focus|
            float inv_d = 1.0f / distance_pc;
            float inv_focus = 1.0f / m_Config.focus_distance_pc;
            float defocus = std::abs(inv_d - inv_focus);

            // Convert to angular size
            sigma_dof = m_Config.aperture_mm * defocus * 1e-3f;
        }

        return std::sqrt(sigma_base * sigma_base + sigma_dof * sigma_dof);
    }

    StarfieldConfig m_Config;
};

} // namespace Sirius
