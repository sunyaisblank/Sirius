#pragma once

// Depth-resolved procedural starfield: a deterministic synthetic catalogue
// with parallax and aperture-dependent depth-of-field blur. It is not an
// imported astrometric catalogue. Magnitudes obey m = -2.5 log10(F/F0) and
// apparent m = M + 5 log10(d/10pc); effective temperature uses the Ballesteros
// (2012) B-V blackbody estimate. Ported from PHSF001A.h.

#include "sirius/core/spectral/blackbody.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numbers>
#include <random>
#include <span>
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

    // RGB colour from the same Planck/CIE/linear-sRGB authority used by disk
    // emission, normalised so the brightest channel is 1.
    void ComputeColor(float& r, float& g, float& b) const {
        float T = temperature_K;
        if (!std::isfinite(T) || T <= 0.0f) T = 5778.0f;  // Default to solar
        const spectral::Rgb color = spectral::BlackbodyToRgb(T);
        r = color.r;
        g = color.g;
        b = color.b;
    }

    // Relative flux from magnitude: F propto 10^(-0.4 m).
    float Intensity() const { return std::pow(10.0f, -0.4f * magnitude); }
};

static_assert(sizeof(StarEntry) == 8 * sizeof(float),
              "StarEntry must remain eight packed floats for the Slang upload ABI");

// Catalog generation and sampling parameters.
struct StarfieldConfig {
    std::uint32_t star_count = 100000;  // Number of stars in catalog
    float min_distance_pc = 1.0f;       // Nearest star distance [pc]
    float max_distance_pc = 10000.0f;   // Farthest star distance [pc]
    float magnitude_limit = 12.0f;      // Faintest visible magnitude
    float aperture_mm = 50.0f;          // Virtual aperture for DoF (0 = infinite)
    float focus_distance_pc = 100.0f;   // Focus plane distance [pc]
    // Catalogue magnitudes are relative fluxes. A factor of 100 maps their
    // zero point into the renderer's display-linear range; scale 1 survives in
    // HDR but quantises away in an ordinary P3 PNG.
    float brightness_scale = 100.0f;
    std::uint32_t seed = 42;  // Random seed
    bool enabled = true;
    bool enable_parallax = true;  // Enable parallax for camera motion
    bool enable_dof = true;       // Enable depth of field blur

    // Invariants (enforced by clamping, not assertion):
    //   star_count in [1, 10^7]
    //   min_distance_pc > 0
    //   aperture_mm in [0, 1000] (0 = infinite DoF)
    void Validate() {
        star_count = std::clamp(star_count, 1u, 10000000u);
        if (!std::isfinite(min_distance_pc)) min_distance_pc = 1.0f;
        if (!std::isfinite(max_distance_pc)) max_distance_pc = 10000.0f;
        if (!std::isfinite(magnitude_limit)) magnitude_limit = 12.0f;
        if (!std::isfinite(aperture_mm)) aperture_mm = 50.0f;
        if (!std::isfinite(focus_distance_pc)) focus_distance_pc = 100.0f;
        if (!std::isfinite(brightness_scale)) brightness_scale = 100.0f;
        min_distance_pc = std::max(min_distance_pc, 0.1f);
        if (!(max_distance_pc > min_distance_pc)) {
            max_distance_pc =
                std::nextafter(min_distance_pc, std::numeric_limits<float>::infinity());
            if (!std::isfinite(max_distance_pc)) {
                min_distance_pc = 1.0f;
                max_distance_pc = 10000.0f;
            }
        }
        magnitude_limit = std::clamp(magnitude_limit, 0.0f, 20.0f);
        aperture_mm = std::clamp(aperture_mm, 0.0f, 1000.0f);
        focus_distance_pc = std::clamp(focus_distance_pc, min_distance_pc, max_distance_pc);
        brightness_scale = std::clamp(brightness_scale, 0.0f, 1000000.0f);
    }
};

// Deterministic latitude/longitude acceleration structure for beam queries.
// The exact Gaussian and angular cutoff remain in StarfieldGenerator; this
// index only produces a conservative candidate superset, so it cannot change
// the represented star field. Its fixed topology is below 1 MiB for a 100k
// catalogue and avoids the previous O(stars * pixels) IMAX failure mode.
class StarfieldSpatialIndex {
  public:
    static constexpr int kThetaBins = 256;
    static constexpr int kPhiBins = 512;

    explicit StarfieldSpatialIndex(const std::vector<StarEntry>& stars) { Build(stars); }

    template <typename Callback>
    void ForEachCandidate(float dir_x, float dir_y, float dir_z, float sigma,
                          Callback&& callback) const {
        if (indices_.empty() || sigma <= 0.0f) return;
        const float length = std::sqrt(dir_x * dir_x + dir_y * dir_y + dir_z * dir_z);
        if (length < 1.0e-20f) return;
        dir_x /= length;
        dir_y /= length;
        dir_z /= length;

        constexpr float kStarfieldPi = 3.14159265358979323846f;
        constexpr float kStarfieldTwoPi = 2.0f * kStarfieldPi;
        const float theta = std::acos(std::clamp(dir_z, -1.0f, 1.0f));
        float phi = std::atan2(dir_y, dir_x);
        if (phi < 0.0f) phi += kStarfieldTwoPi;
        const float cutoff = std::min(4.0f * sigma, kStarfieldPi);

        const int first_theta =
            std::clamp(static_cast<int>(
                           std::floor(std::max(theta - cutoff, 0.0f) / kStarfieldPi * kThetaBins)),
                       0, kThetaBins - 1);
        const int last_theta =
            std::clamp(static_cast<int>(std::floor(std::min(theta + cutoff, kStarfieldPi) /
                                                   kStarfieldPi * kThetaBins)),
                       0, kThetaBins - 1);

        float phi_half_width = kStarfieldPi;
        const float sin_theta = std::sin(theta);
        if (theta > cutoff && kStarfieldPi - theta > cutoff && sin_theta > 1.0e-8f) {
            phi_half_width = std::asin(std::clamp(std::sin(cutoff) / sin_theta, 0.0f, 1.0f));
        }
        const float phi_bin_width = kStarfieldTwoPi / static_cast<float>(kPhiBins);
        const int centre_phi =
            std::min(static_cast<int>(std::floor(phi / kStarfieldTwoPi * kPhiBins)), kPhiBins - 1);
        const int phi_radius =
            std::min(static_cast<int>(std::ceil(phi_half_width / phi_bin_width)) + 1, kPhiBins);
        const bool all_phi = 2 * phi_radius + 1 >= kPhiBins;

        for (int theta_bin = first_theta; theta_bin <= last_theta; ++theta_bin) {
            if (all_phi) {
                for (int phi_bin = 0; phi_bin < kPhiBins; ++phi_bin) {
                    VisitCell(theta_bin, phi_bin, callback);
                }
                continue;
            }
            for (int delta = -phi_radius; delta <= phi_radius; ++delta) {
                int phi_bin = (centre_phi + delta) % kPhiBins;
                if (phi_bin < 0) phi_bin += kPhiBins;
                VisitCell(theta_bin, phi_bin, callback);
            }
        }
    }

    [[nodiscard]] std::size_t MemoryBytes() const {
        return offsets_.size() * sizeof(std::uint32_t) + indices_.size() * sizeof(std::uint32_t);
    }
    [[nodiscard]] std::span<const std::uint32_t> Offsets() const noexcept { return offsets_; }
    [[nodiscard]] std::span<const std::uint32_t> Indices() const noexcept { return indices_; }

  private:
    std::vector<std::uint32_t> offsets_;
    std::vector<std::uint32_t> indices_;

    static std::size_t Cell(int theta_bin, int phi_bin) {
        return static_cast<std::size_t>(theta_bin) * kPhiBins + phi_bin;
    }

    void Build(const std::vector<StarEntry>& stars) {
        constexpr float kStarfieldPi = 3.14159265358979323846f;
        constexpr float kStarfieldTwoPi = 2.0f * kStarfieldPi;
        constexpr std::size_t kCells = static_cast<std::size_t>(kThetaBins) * kPhiBins;
        offsets_.assign(kCells + 1, 0);
        std::vector<std::size_t> cells(stars.size());
        for (std::size_t i = 0; i < stars.size(); ++i) {
            const auto& star = stars[i];
            const float theta = std::acos(std::clamp(star.direction_z, -1.0f, 1.0f));
            float phi = std::atan2(star.direction_y, star.direction_x);
            if (phi < 0.0f) phi += kStarfieldTwoPi;
            const int theta_bin = std::clamp(
                static_cast<int>(std::floor(theta / kStarfieldPi * kThetaBins)), 0, kThetaBins - 1);
            const int phi_bin = std::clamp(
                static_cast<int>(std::floor(phi / kStarfieldTwoPi * kPhiBins)), 0, kPhiBins - 1);
            cells[i] = Cell(theta_bin, phi_bin);
            ++offsets_[cells[i] + 1];
        }
        for (std::size_t cell = 1; cell < offsets_.size(); ++cell) {
            offsets_[cell] += offsets_[cell - 1];
        }
        indices_.resize(stars.size());
        std::vector<std::uint32_t> cursor(offsets_.begin(), offsets_.end() - 1);
        for (std::size_t i = 0; i < cells.size(); ++i) {
            indices_[cursor[cells[i]]++] = static_cast<std::uint32_t>(i);
        }
    }

    template <typename Callback>
    void VisitCell(int theta_bin, int phi_bin, Callback& callback) const {
        const std::size_t cell = Cell(theta_bin, phi_bin);
        for (std::uint32_t at = offsets_[cell]; at < offsets_[cell + 1]; ++at) {
            callback(indices_[at]);
        }
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

        for (std::uint32_t i = 0; i < config_.star_count; ++i) {
            StarEntry star = GenerateStar(rng, uniform);
            if (star.magnitude <= config_.magnitude_limit) {
                stars.push_back(star);
            }
        }

        return stars;
    }

    // Full deterministic catalogue for the filtered point-source star field (P3):
    // exactly star_count entries on the celestial sphere with magnitudes and
    // temperatures and no magnitude cull. The function preserves the requested
    // count exactly; the operating profile requests at least 10^5. The fixed seed
    // makes it reproducible frame to frame, which the anti-flicker measurement
    // relies on.
    std::vector<StarEntry> GenerateCatalogue() const {
        std::vector<StarEntry> stars;
        stars.reserve(config_.star_count);
        std::mt19937 rng(config_.seed);
        std::uniform_real_distribution<float> uniform(0.0f, 1.0f);
        for (std::uint32_t i = 0; i < config_.star_count; ++i) {
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
        float cos_cut = std::cos(std::min(4.0f * sigma, static_cast<float>(std::numbers::pi)));
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

    // Indexed form of the exact beam accumulation above. The index returns a
    // conservative angular candidate superset; the same dot-product cutoff and
    // Gaussian decide every contribution.
    void AccumulateThroughBeam(float dir_x, float dir_y, float dir_z, float sigma,
                               const std::vector<StarEntry>& stars,
                               const StarfieldSpatialIndex& index, float& r, float& g,
                               float& b) const {
        r = g = b = 0.0f;
        if (stars.empty() || sigma <= 0.0f) return;
        const float dlen = std::sqrt(dir_x * dir_x + dir_y * dir_y + dir_z * dir_z);
        if (dlen < 1.0e-20f) return;
        dir_x /= dlen;
        dir_y /= dlen;
        dir_z /= dlen;
        const float cos_cut =
            std::cos(std::min(4.0f * sigma, static_cast<float>(std::numbers::pi)));
        const float inv_two_sigma2 = 1.0f / (2.0f * sigma * sigma);

        index.ForEachCandidate(dir_x, dir_y, dir_z, sigma, [&](std::uint32_t star_index) {
            const auto& star = stars[star_index];
            const float cosine =
                dir_x * star.direction_x + dir_y * star.direction_y + dir_z * star.direction_z;
            if (cosine < cos_cut) return;
            const float angle = std::acos(std::clamp(cosine, -1.0f, 1.0f));
            const float weight = std::exp(-angle * angle * inv_two_sigma2);
            const float intensity = star.Intensity() * weight * config_.brightness_scale;
            float sr = 0.0f;
            float sg = 0.0f;
            float sb = 0.0f;
            star.ComputeColor(sr, sg, sb);
            r += sr * intensity;
            g += sg * intensity;
            b += sb * intensity;
        });
    }

    // Elliptical DNGR footprint. `orientation` is the position angle of the
    // major axis in the deterministic tangent basis used by the ray-bundle
    // geometry extractor. The spatial index is queried with the major axis, a
    // conservative circular bound; the anisotropic Gaussian is the exact
    // contribution decision.
    void AccumulateThroughBeam(float dir_x, float dir_y, float dir_z, float sigma_major,
                               float sigma_minor, float orientation,
                               const std::vector<StarEntry>& stars,
                               const StarfieldSpatialIndex& index, float& r, float& g,
                               float& b) const {
        r = g = b = 0.0f;
        if (stars.empty() || sigma_major <= 0.0f || sigma_minor <= 0.0f) return;
        const float dlen = std::sqrt(dir_x * dir_x + dir_y * dir_y + dir_z * dir_z);
        if (dlen < 1.0e-20f) return;
        dir_x /= dlen;
        dir_y /= dlen;
        dir_z /= dlen;

        // Same transverse-basis convention as GeodesicTracer::FinaliseBundle:
        // project z unless the direction is near z, then project x.
        float ref_x = 0.0f;
        float ref_y = 0.0f;
        float ref_z = 1.0f;
        if (std::abs(dir_z) > 0.9f) {
            ref_x = 1.0f;
            ref_z = 0.0f;
        }
        const float ref_dot = ref_x * dir_x + ref_y * dir_y + ref_z * dir_z;
        float ex = ref_x - ref_dot * dir_x;
        float ey = ref_y - ref_dot * dir_y;
        float ez = ref_z - ref_dot * dir_z;
        const float elen = std::sqrt(ex * ex + ey * ey + ez * ez);
        ex /= elen;
        ey /= elen;
        ez /= elen;
        const float fx = dir_y * ez - dir_z * ey;
        const float fy = dir_z * ex - dir_x * ez;
        const float fz = dir_x * ey - dir_y * ex;

        const float major = std::max(sigma_major, sigma_minor);
        const float minor = std::min(sigma_major, sigma_minor);
        const float cos_orientation = std::cos(orientation);
        const float sin_orientation = std::sin(orientation);
        const float cos_cut =
            std::cos(std::min(4.0f * major, static_cast<float>(std::numbers::pi)));
        const float inv_major_squared = 1.0f / (major * major);
        const float inv_minor_squared = 1.0f / (minor * minor);

        index.ForEachCandidate(dir_x, dir_y, dir_z, major, [&](std::uint32_t star_index) {
            const auto& star = stars[star_index];
            const float cosine =
                dir_x * star.direction_x + dir_y * star.direction_y + dir_z * star.direction_z;
            if (cosine < cos_cut) return;
            const float angle = std::acos(std::clamp(cosine, -1.0f, 1.0f));
            float tangent_x = 0.0f;
            float tangent_y = 0.0f;
            const float sin_angle = std::sin(angle);
            if (sin_angle > 1.0e-8f) {
                const float tx = (star.direction_x - cosine * dir_x) / sin_angle;
                const float ty = (star.direction_y - cosine * dir_y) / sin_angle;
                const float tz = (star.direction_z - cosine * dir_z) / sin_angle;
                tangent_x = angle * (tx * ex + ty * ey + tz * ez);
                tangent_y = angle * (tx * fx + ty * fy + tz * fz);
            }
            const float along_major = cos_orientation * tangent_x + sin_orientation * tangent_y;
            const float along_minor = -sin_orientation * tangent_x + cos_orientation * tangent_y;
            const float exponent = -0.5f * (along_major * along_major * inv_major_squared +
                                            along_minor * along_minor * inv_minor_squared);
            const float weight = std::exp(exponent);
            const float intensity = star.Intensity() * weight * config_.brightness_scale;
            float sr = 0.0f;
            float sg = 0.0f;
            float sb = 0.0f;
            star.ComputeColor(sr, sg, sb);
            r += sr * intensity;
            g += sg * intensity;
            b += sb * intensity;
        });
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
            const float camera_x = config_.enable_parallax ? cam_x_pc : 0.0f;
            const float camera_y = config_.enable_parallax ? cam_y_pc : 0.0f;
            const float camera_z = config_.enable_parallax ? cam_z_pc : 0.0f;
            float dx = star.direction_x * star.distance_pc - camera_x;
            float dy = star.direction_y * star.distance_pc - camera_y;
            float dz = star.direction_z * star.distance_pc - camera_z;
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
        float phi = 2.0f * static_cast<float>(std::numbers::pi) * uniform(rng);
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

        // Ballesteros (2012), EPL 97 34008, blackbody estimate from Johnson B-V.
        const float colour_term = 0.92f * star.color_bv;
        star.temperature_K = 4600.0f * (1.0f / (colour_term + 1.7f) + 1.0f / (colour_term + 0.62f));
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
