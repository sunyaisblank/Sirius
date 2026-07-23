#pragma once

// Depth-resolved procedural starfield: parallax-correct stellar catalog with
// aperture-dependent depth-of-field blur for cinematic backgrounds. Magnitude
// system m = -2.5 log10(F/F0), apparent m = M + 5 log10(d/10pc). References:
// Hipparcos (ESA 1997), Tycho-2 (Hoeg et al. 2000). Ported from PHSF001A.h.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <random>
#include <vector>

namespace sirius::core {

// A single catalog star: unit direction, distance, brightness, and colour.
struct StarEntry {
    float direction_x;    // Unit direction X (galactic coords)
    float direction_y;    // Unit direction Y
    float direction_z;    // Unit direction Z
    float distance_pc;    // Distance in parsecs
    float magnitude;      // Apparent magnitude
    float color_bv;       // B-V color index
    float temperature_K;  // Effective temperature [K]
    float padding;        // Alignment

    // RGB colour from a simplified Planckian locus (T ~ 4600 / (0.92 + B-V)
    // for B-V > 0), normalised so the brightest channel is 1.
    void ComputeColor(float& r, float& g, float& b) const {
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
        r /= max_c;
        g /= max_c;
        b /= max_c;
    }

    // Relative flux from magnitude: F propto 10^(-0.4 m).
    float Intensity() const { return std::pow(10.0f, -0.4f * magnitude); }
};

// Catalog generation and sampling parameters.
struct StarfieldConfig {
    uint32_t star_count = 100000;      // Number of stars in catalog
    float min_distance_pc = 1.0f;      // Nearest star distance [pc]
    float max_distance_pc = 10000.0f;  // Farthest star distance [pc]
    float magnitude_limit = 12.0f;     // Faintest visible magnitude
    float aperture_mm = 50.0f;         // Virtual aperture for DoF (0 = infinite)
    float focus_distance_pc = 100.0f;  // Focus plane distance [pc]
    float brightness_scale = 1.0f;     // Overall brightness multiplier
    uint32_t seed = 42;                // Random seed
    bool enabled = true;
    bool enable_parallax = true;  // Enable parallax for camera motion
    bool enable_dof = true;       // Enable depth of field blur

    // Star distribution parameters (galactic model)
    float thin_disk_scale_height_pc = 300.0f;    // Thin disk scale height [pc]
    float thick_disk_scale_height_pc = 1000.0f;  // Thick disk scale height [pc]
    float disk_scale_length_pc = 2500.0f;        // Disk radial scale [pc]
    float thin_disk_fraction = 0.9f;             // Fraction of stars in thin disk

    // Invariants (enforced by clamping, not assertion):
    //   star_count <= 10^7
    //   min_distance_pc > 0
    //   aperture_mm in [0, 1000] (0 = infinite DoF)
    void Validate() {
        star_count = std::min(star_count, 10000000u);
        min_distance_pc = std::max(min_distance_pc, 0.1f);
        max_distance_pc = std::max(max_distance_pc, min_distance_pc + 1.0f);
        magnitude_limit = std::clamp(magnitude_limit, 0.0f, 20.0f);
        aperture_mm = std::clamp(aperture_mm, 0.0f, 1000.0f);
        focus_distance_pc = std::clamp(focus_distance_pc, min_distance_pc, max_distance_pc);
    }
};

// Generates and samples a procedural star catalog.
class StarfieldGenerator {
  public:
    explicit StarfieldGenerator(const StarfieldConfig& config = StarfieldConfig())
        : config_(config) {
        config_.Validate();
    }

    // Generate the procedural star catalog (stars brighter than magnitude_limit).
    std::vector<StarEntry> Generate() {
        std::vector<StarEntry> stars;
        stars.reserve(config_.star_count);

        std::mt19937 rng(config_.seed);
        std::uniform_real_distribution<float> uniform(0.0f, 1.0f);

        for (uint32_t i = 0; i < config_.star_count; ++i) {
            StarEntry star = GenerateStar(rng, uniform);
            if (star.magnitude <= config_.magnitude_limit) {
                stars.push_back(star);
            }
        }

        return stars;
    }

    // Full deterministic catalogue for the filtered point-source star field (P3):
    // exactly star_count entries on the celestial sphere with magnitudes and
    // temperatures, no magnitude cull, so the catalogue size is guaranteed to meet
    // the >= 10^5 gate. The fixed seed makes it reproducible frame to frame, which
    // the anti-flicker measurement relies on.
    std::vector<StarEntry> GenerateCatalogue() const {
        std::vector<StarEntry> stars;
        stars.reserve(config_.star_count);
        std::mt19937 rng(config_.seed);
        std::uniform_real_distribution<float> uniform(0.0f, 1.0f);
        for (uint32_t i = 0; i < config_.star_count; ++i) {
            stars.push_back(GenerateStar(rng, uniform));
        }
        return stars;
    }

    // Accumulate catalogue stars as point sources filtered through a beam of
    // angular radius sigma (radians): each star's flux is weighted by a Gaussian
    // footprint exp(-theta^2 / (2 sigma^2)), DNGR's anti-flicker approach (James
    // et al. 2015, CQG 32 065001, section 3). A wide (lensed) beam gathers nearby
    // stars smoothly as the camera moves; a pinhole (tiny sigma) samples them
    // discontinuously, which is exactly the flicker the beam removes. dir_* is the
    // escape direction (Cartesian celestial-sphere unit vector).
    void AccumulateThroughBeam(float dir_x, float dir_y, float dir_z, float sigma,
                               const std::vector<StarEntry>& stars, float& r, float& g,
                               float& b) const {
        r = g = b = 0.0f;
        if (stars.empty() || sigma <= 0.0f) return;

        float dlen = std::sqrt(dir_x * dir_x + dir_y * dir_y + dir_z * dir_z);
        if (dlen < 1e-20f) return;
        dir_x /= dlen;
        dir_y /= dlen;
        dir_z /= dlen;

        // Stars beyond four sigma contribute nothing; the cosine test skips them
        // before the transcendental acos.
        float cos_cut = std::cos(std::min(4.0f * sigma, static_cast<float>(M_PI)));
        float inv_two_sigma2 = 1.0f / (2.0f * sigma * sigma);

        for (const auto& s : stars) {
            float ca = dir_x * s.direction_x + dir_y * s.direction_y + dir_z * s.direction_z;
            if (ca < cos_cut) continue;
            float angle = std::acos(std::clamp(ca, -1.0f, 1.0f));
            float w = std::exp(-angle * angle * inv_two_sigma2);
            float intensity = s.Intensity() * w * config_.brightness_scale;
            float sr, sg, sb;
            s.ComputeColor(sr, sg, sb);
            r += sr * intensity;
            g += sg * intensity;
            b += sb * intensity;
        }
    }

    // Accumulate starfield colour along view direction (dir_*) with a parallax
    // offset from the camera position (cam_*_pc) into (r, g, b).
    void SampleStarfield(float dir_x, float dir_y, float dir_z, float cam_x_pc, float cam_y_pc,
                         float cam_z_pc, const std::vector<StarEntry>& stars, float& r, float& g,
                         float& b) const {
        r = g = b = 0.0f;

        if (!config_.enabled || stars.empty()) return;

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
            float sigma = ComputePsfSigma(star.distance_pc);
            float angle = std::acos(std::clamp(cos_angle, -1.0f, 1.0f));
            float psf = std::exp(-0.5f * angle * angle / (sigma * sigma));

            if (psf < 1e-6f) continue;

            // Intensity and color
            float intensity = star.Intensity() * psf * config_.brightness_scale;
            float sr, sg, sb;
            star.ComputeColor(sr, sg, sb);

            r += sr * intensity;
            g += sg * intensity;
            b += sb * intensity;
        }
    }

    const StarfieldConfig& config() const { return config_; }
    void SetConfig(const StarfieldConfig& config) {
        config_ = config;
        config_.Validate();
    }

  private:
    // Generate a single star with realistic direction, distance, and colour.
    StarEntry GenerateStar(std::mt19937& rng,
                           std::uniform_real_distribution<float>& uniform) const {
        StarEntry star;

        // Direction: uniform on sphere
        float theta = std::acos(2.0f * uniform(rng) - 1.0f);
        float phi = 2.0f * static_cast<float>(M_PI) * uniform(rng);
        star.direction_x = std::sin(theta) * std::cos(phi);
        star.direction_y = std::sin(theta) * std::sin(phi);
        star.direction_z = std::cos(theta);

        // Distance: power-law weighted towards nearby
        float u = uniform(rng);
        float d_min = config_.min_distance_pc;
        float d_max = config_.max_distance_pc;
        // r^2 dr weighting (uniform in volume)
        star.distance_pc =
            std::cbrt(u * (d_max * d_max * d_max) + (1.0f - u) * (d_min * d_min * d_min));

        // Absolute magnitude: Salpeter IMF approximation
        // More low-luminosity stars
        float abs_mag = -2.0f + 12.0f * std::pow(uniform(rng), 0.3f);

        // Apparent magnitude: m = M + 5 log10(d/10)
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

    // Point spread function sigma (angular), combining the diffraction limit
    // with depth-of-field defocus (circle of confusion propto A * |1/d - 1/d_focus|).
    float ComputePsfSigma(float distance_pc) const {
        // Base angular size (point source)
        float sigma_base = 1e-6f;  // ~0.2 arcsec

        // Depth of field blur
        float sigma_dof = 0.0f;
        if (config_.enable_dof && config_.aperture_mm > 0.0f) {
            // Circle of confusion: CoC propto A * |1/d - 1/d_focus|
            float inv_d = 1.0f / distance_pc;
            float inv_focus = 1.0f / config_.focus_distance_pc;
            float defocus = std::abs(inv_d - inv_focus);

            // Convert to angular size
            sigma_dof = config_.aperture_mm * defocus * 1e-3f;
        }

        return std::sqrt(sigma_base * sigma_base + sigma_dof * sigma_dof);
    }

    StarfieldConfig config_;
};

}  // namespace sirius::core
