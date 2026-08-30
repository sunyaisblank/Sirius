#pragma once

// Depth-resolved procedural starfield: a deterministic synthetic catalogue
// with parallax and aperture-dependent depth-of-field blur. It is not an
// imported astrometric catalogue. Magnitudes obey m = -2.5 log10(F/F0) and
// apparent m = M + 5 log10(d/10pc); effective temperature uses the Ballesteros
// (2012) B-V blackbody estimate. Ported from PHSF001A.h.

#include "sirius/base/contracts.h"
#include "sirius/core/celestial_tangent_basis.h"
#include "sirius/core/spectral/blackbody.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numbers>
#include <random>
#include <span>
#include <utility>
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
        SIRIUS_PRE(std::isfinite(temperature_K) && temperature_K > 0.0f);
        const spectral::Rgb color = spectral::BlackbodyToRgb(temperature_K);
        r = color.r;
        g = color.g;
        b = color.b;
    }

    // Relative flux from magnitude: F propto 10^(-0.4 m).
    float Intensity() const {
        // At m < -96 the exact relative flux exceeds binary32's finite range.
        SIRIUS_PRE(std::isfinite(magnitude) && magnitude >= -96.0f);
        return std::pow(10.0f, -0.4f * magnitude);
    }
};

[[nodiscard]] inline bool IsRepresentedStarEntry(const StarEntry& star) noexcept {
    const float direction_norm_squared = star.direction_x * star.direction_x +
                                         star.direction_y * star.direction_y +
                                         star.direction_z * star.direction_z;
    return std::isfinite(star.direction_x) && std::isfinite(star.direction_y) &&
           std::isfinite(star.direction_z) && std::isfinite(direction_norm_squared) &&
           std::abs(direction_norm_squared - 1.0f) <= 2.0e-4f && std::isfinite(star.distance_pc) &&
           star.distance_pc > 0.0f && std::isfinite(star.magnitude) && star.magnitude >= -96.0f &&
           std::isfinite(star.color_bv) && std::isfinite(star.temperature_K) &&
           star.temperature_K >= 100.0f && star.temperature_K <= 1.0e6f && star.padding == 0.0f;
}

static_assert(sizeof(StarEntry) == 8 * sizeof(float),
              "StarEntry must remain eight packed floats for the Slang upload ABI");

// Catalog generation and sampling parameters.
struct StarfieldConfig {
    friend bool operator==(const StarfieldConfig&, const StarfieldConfig&) = default;

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
};

// The catalogue and typed-session paths share one represented numerical
// domain. Invalid requests are rejected by their consuming boundary; this
// predicate never mutates a request into a different catalogue.
[[nodiscard]] inline bool IsRepresentedStarfieldConfig(const StarfieldConfig& config) noexcept {
    return config.star_count >= 1 && config.star_count <= 10000000u &&
           std::isfinite(config.min_distance_pc) && config.min_distance_pc >= 0.1f &&
           std::isfinite(config.max_distance_pc) &&
           config.max_distance_pc > config.min_distance_pc &&
           std::isfinite(config.magnitude_limit) && config.magnitude_limit >= 0.0f &&
           config.magnitude_limit <= 20.0f && std::isfinite(config.aperture_mm) &&
           config.aperture_mm >= 0.0f && config.aperture_mm <= 1000.0f &&
           std::isfinite(config.focus_distance_pc) &&
           config.focus_distance_pc >= config.min_distance_pc &&
           config.focus_distance_pc <= config.max_distance_pc &&
           std::isfinite(config.brightness_scale) && config.brightness_scale >= 0.0f &&
           config.brightness_scale <= 1000000.0f;
}

// Exact request consumed by the session's indexed DNGR point catalogue. The
// broader StarfieldConfig also carries magnitude-cull, parallax, and synthetic
// depth-of-field controls used by other sampling APIs; exposing them here would
// accept controls that this path cannot consume.
struct PointStarfieldConfig {
    friend bool operator==(const PointStarfieldConfig&, const PointStarfieldConfig&) = default;

    std::uint32_t star_count = 100000;
    float min_distance_pc = 1.0f;
    float max_distance_pc = 10000.0f;
    float brightness_scale = 100.0f;
    std::uint32_t seed = 42;
};

[[nodiscard]] inline bool IsRepresentedPointStarfieldConfig(
    const PointStarfieldConfig& config) noexcept {
    return config.star_count >= 1 && config.star_count <= 10000000u &&
           std::isfinite(config.min_distance_pc) && config.min_distance_pc >= 0.1f &&
           std::isfinite(config.max_distance_pc) &&
           config.max_distance_pc > config.min_distance_pc &&
           std::isfinite(config.brightness_scale) && config.brightness_scale >= 0.0f &&
           config.brightness_scale <= 1000000.0f;
}

[[nodiscard]] inline StarfieldConfig ExpandPointStarfieldConfig(
    const PointStarfieldConfig& config) {
    SIRIUS_PRE(IsRepresentedPointStarfieldConfig(config));
    StarfieldConfig expanded;
    expanded.star_count = config.star_count;
    expanded.min_distance_pc = config.min_distance_pc;
    expanded.max_distance_pc = config.max_distance_pc;
    // Generic depth-of-field sampling is not part of the point-catalogue
    // request, but the shared generator object still requires a represented
    // internal focus value. The nearest catalogue distance is exact and valid
    // for every admitted ordered domain.
    expanded.focus_distance_pc = config.min_distance_pc;
    expanded.brightness_scale = config.brightness_scale;
    expanded.seed = config.seed;
    return expanded;
}

// Deterministic latitude/longitude acceleration structure for beam queries.
// The exact Gaussian and angular cutoff remain in StarfieldGenerator; this
// index only produces a conservative candidate superset, so it cannot change
// the represented star field. Its fixed topology is below 1 MiB for a 100k
// catalogue and avoids the previous O(stars * pixels) IMAX failure mode.
class StarfieldSpatialIndex {
  public:
    static constexpr int kThetaBins = 256;
    static constexpr int kPhiBins = 512;

    explicit StarfieldSpatialIndex(std::vector<StarEntry> stars) : stars_(std::move(stars)) {
        SIRIUS_PRE(stars_.size() <= std::numeric_limits<std::uint32_t>::max());
        SIRIUS_PRE(std::all_of(stars_.begin(), stars_.end(), IsRepresentedStarEntry));
        Build();
    }

    [[nodiscard]] std::size_t Size() const noexcept { return stars_.size(); }
    [[nodiscard]] std::span<const StarEntry> Stars() const noexcept { return stars_; }

    template <typename Callback>
    void ForEachCandidate(float dir_x, float dir_y, float dir_z, float sigma,
                          Callback&& callback) const {
        SIRIUS_PRE(std::isfinite(dir_x) && std::isfinite(dir_y) && std::isfinite(dir_z));
        SIRIUS_PRE(std::isfinite(sigma) && sigma > 0.0f);
        if (indices_.empty()) return;
        const float length = std::sqrt(dir_x * dir_x + dir_y * dir_y + dir_z * dir_z);
        SIRIUS_PRE(std::isfinite(length) && length >= 1.0e-20f);
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
    std::vector<StarEntry> stars_;
    std::vector<std::uint32_t> offsets_;
    std::vector<std::uint32_t> indices_;

    static std::size_t Cell(int theta_bin, int phi_bin) {
        return static_cast<std::size_t>(theta_bin) * kPhiBins + phi_bin;
    }

    void Build() {
        constexpr float kStarfieldPi = 3.14159265358979323846f;
        constexpr float kStarfieldTwoPi = 2.0f * kStarfieldPi;
        constexpr std::size_t kCells = static_cast<std::size_t>(kThetaBins) * kPhiBins;
        offsets_.assign(kCells + 1, 0);
        std::vector<std::size_t> cells(stars_.size());
        for (std::size_t i = 0; i < stars_.size(); ++i) {
            const auto& star = stars_[i];
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
        indices_.resize(stars_.size());
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
        SIRIUS_PRE(IsRepresentedStarfieldConfig(config));
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
        SIRIUS_PRE(std::isfinite(dir_x) && std::isfinite(dir_y) && std::isfinite(dir_z));
        SIRIUS_PRE(std::isfinite(sigma) && sigma > 0.0f);
        SIRIUS_PRE(std::all_of(stars.begin(), stars.end(), IsRepresentedStarEntry));
        if (stars.empty()) return;

        float dlen = std::sqrt(dir_x * dir_x + dir_y * dir_y + dir_z * dir_z);
        SIRIUS_PRE(std::isfinite(dlen) && dlen >= 1.0e-20f);
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
                               const StarfieldSpatialIndex& index, float& r, float& g,
                               float& b) const {
        r = g = b = 0.0f;
        SIRIUS_PRE(std::isfinite(dir_x) && std::isfinite(dir_y) && std::isfinite(dir_z));
        SIRIUS_PRE(std::isfinite(sigma) && sigma > 0.0f);
        const std::span<const StarEntry> stars = index.Stars();
        if (stars.empty()) return;
        const float dlen = std::sqrt(dir_x * dir_x + dir_y * dir_y + dir_z * dir_z);
        SIRIUS_PRE(std::isfinite(dlen) && dlen >= 1.0e-20f);
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
                               const StarfieldSpatialIndex& index, float& r, float& g,
                               float& b) const {
        r = g = b = 0.0f;
        SIRIUS_PRE(std::isfinite(dir_x) && std::isfinite(dir_y) && std::isfinite(dir_z));
        SIRIUS_PRE(std::isfinite(sigma_major) && sigma_major > 0.0f);
        SIRIUS_PRE(std::isfinite(sigma_minor) && sigma_minor > 0.0f);
        SIRIUS_PRE(std::isfinite(orientation));
        const std::span<const StarEntry> stars = index.Stars();
        if (stars.empty()) return;
        const float dlen = std::sqrt(dir_x * dir_x + dir_y * dir_y + dir_z * dir_z);
        SIRIUS_PRE(std::isfinite(dlen) && dlen >= 1.0e-20f);
        dir_x /= dlen;
        dir_y /= dlen;
        dir_z /= dlen;

        // The ellipse angle was measured in this exact Sachs basis at the
        // terminal ray. Reusing it here preserves the footprint orientation.
        const auto tangent_basis =
            relativity::MakeCelestialTangentBasis(std::array{dir_x, dir_y, dir_z});
        SIRIUS_ASSERT(tangent_basis.has_value());
        if (!tangent_basis.has_value()) return;
        const float ex = tangent_basis->first[0];
        const float ey = tangent_basis->first[1];
        const float ez = tangent_basis->first[2];
        const float fx = tangent_basis->second[0];
        const float fy = tangent_basis->second[1];
        const float fz = tangent_basis->second[2];

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
        SIRIUS_PRE(IsRepresentedStarfieldConfig(config));
        config_ = config;
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
        const double u = uniform(rng);
        const double d_min = config_.min_distance_pc;
        const double d_max = config_.max_distance_pc;
        // r^2 dr weighting (uniform in volume)
        star.distance_pc = static_cast<float>(
            std::cbrt(u * (d_max * d_max * d_max) + (1.0 - u) * (d_min * d_min * d_min)));

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
