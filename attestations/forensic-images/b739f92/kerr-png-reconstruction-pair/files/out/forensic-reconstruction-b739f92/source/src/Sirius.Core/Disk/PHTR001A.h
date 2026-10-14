// PHTR001A.h - Turbulence Model for Accretion Disks
// Component ID: PHTR001A (Physics/Turbulence/TurbulenceModel)
//
// Kolmogorov cascade density perturbations for realistic GRMHD-informed
// volumetric plasma distribution in accretion disks.
//
// PHYSICAL FOUNDATION:
// ====================
// GRMHD simulations show accretion disks are highly turbulent, driven by
// the magnetorotational instability (MRI). The resulting density fluctuations
// follow Kolmogorov-type cascades in the inertial range.
//
// Power spectrum: P(k) ∝ k^(-5/3) for 3D Kolmogorov turbulence
// Density perturbations: δρ/ρ ~ 0.1-0.5 in typical AGN disks
//
// We use fractional Brownian motion (fBm) noise to generate density
// perturbations with the correct spectral slope.
//
// REFERENCES:
// - Kolmogorov (1941) "The local structure of turbulence"
// - Balbus & Hawley (1991) ApJ 376, 214 (MRI)
// - McKinney et al. (2012) MNRAS 423, 3083 (GRMHD disk simulations)

#pragma once

#include <cmath>
#include <cstdint>
#include <algorithm>

namespace Sirius {

//==============================================================================
// Turbulence Configuration
//==============================================================================
struct TurbulenceConfig {
    float kolmogorov_exponent = -5.0f / 3.0f;  ///< Spectral exponent (-5/3 for Kolmogorov)
    float outer_scale_M = 5.0f;                 ///< Outer turbulence scale in M (energy injection)
    float inner_scale_M = 0.1f;                 ///< Inner dissipation scale in M
    float amplitude = 0.3f;                     ///< Density fluctuation amplitude δρ/ρ [0, 1]
    uint32_t octaves = 6;                       ///< Noise octaves (fractal detail levels)
    uint32_t seed = 12345;                      ///< Random seed for reproducibility
    float lacunarity = 2.0f;                    ///< Frequency multiplier per octave
    float persistence = 0.5f;                   ///< Amplitude decay per octave
    bool enabled = true;

    // =========================================================================
    // Invariants
    // =========================================================================
    // kolmogorov_exponent ∈ [-2.0, -1.5]
    // outer_scale_M > inner_scale_M > 0
    // amplitude ∈ [0, 1]
    // octaves ∈ [1, 8]

    /// @brief Validate and clamp parameters to valid ranges
    void validate() {
        kolmogorov_exponent = std::clamp(kolmogorov_exponent, -2.0f, -1.5f);
        if (inner_scale_M >= outer_scale_M) {
            inner_scale_M = outer_scale_M * 0.01f;
        }
        inner_scale_M = std::max(inner_scale_M, 0.001f);
        outer_scale_M = std::max(outer_scale_M, 0.01f);
        amplitude = std::clamp(amplitude, 0.0f, 1.0f);
        octaves = std::clamp(octaves, 1u, 8u);
        lacunarity = std::max(lacunarity, 1.1f);
        persistence = std::clamp(persistence, 0.1f, 0.9f);
    }
};

//==============================================================================
// Perlin Noise Implementation (CPU)
// Gradient noise for turbulence generation
//==============================================================================

namespace TurbulenceNoise {

/// @brief Hash function for gradient generation
inline uint32_t hash(uint32_t x, uint32_t y, uint32_t z, uint32_t seed) {
    uint32_t h = seed;
    h ^= x * 374761393u;
    h ^= y * 668265263u;
    h ^= z * 1013904223u;
    h = ((h >> 16) ^ h) * 2654435761u;
    h = ((h >> 16) ^ h) * 2654435761u;
    return h;
}

/// @brief Generate gradient vector from hash
inline void gradient(uint32_t h, float& gx, float& gy, float& gz) {
    // 12 gradient directions (edges of cube + some face diagonals)
    const float grads[12][3] = {
        { 1, 1, 0}, {-1, 1, 0}, { 1,-1, 0}, {-1,-1, 0},
        { 1, 0, 1}, {-1, 0, 1}, { 1, 0,-1}, {-1, 0,-1},
        { 0, 1, 1}, { 0,-1, 1}, { 0, 1,-1}, { 0,-1,-1}
    };
    int idx = h % 12;
    gx = grads[idx][0];
    gy = grads[idx][1];
    gz = grads[idx][2];
}

/// @brief Smoothstep interpolation (Hermite quintic)
inline float fade(float t) {
    return t * t * t * (t * (t * 6.0f - 15.0f) + 10.0f);
}

/// @brief Linear interpolation
inline float lerp(float a, float b, float t) {
    return a + t * (b - a);
}

/// @brief 3D Perlin noise [-1, 1]
inline float perlin3D(float x, float y, float z, uint32_t seed) {
    // Integer cell coordinates
    int xi = static_cast<int>(std::floor(x));
    int yi = static_cast<int>(std::floor(y));
    int zi = static_cast<int>(std::floor(z));

    // Fractional coordinates within cell
    float xf = x - xi;
    float yf = y - yi;
    float zf = z - zi;

    // Fade curves
    float u = fade(xf);
    float v = fade(yf);
    float w = fade(zf);

    // Hash and gradient at each corner
    float n000, n001, n010, n011, n100, n101, n110, n111;

    auto cornerNoise = [&](int cx, int cy, int cz, float fx, float fy, float fz) {
        float gx, gy, gz;
        gradient(hash(xi + cx, yi + cy, zi + cz, seed), gx, gy, gz);
        return gx * fx + gy * fy + gz * fz;
    };

    n000 = cornerNoise(0, 0, 0, xf,       yf,       zf);
    n100 = cornerNoise(1, 0, 0, xf - 1.0f, yf,       zf);
    n010 = cornerNoise(0, 1, 0, xf,       yf - 1.0f, zf);
    n110 = cornerNoise(1, 1, 0, xf - 1.0f, yf - 1.0f, zf);
    n001 = cornerNoise(0, 0, 1, xf,       yf,       zf - 1.0f);
    n101 = cornerNoise(1, 0, 1, xf - 1.0f, yf,       zf - 1.0f);
    n011 = cornerNoise(0, 1, 1, xf,       yf - 1.0f, zf - 1.0f);
    n111 = cornerNoise(1, 1, 1, xf - 1.0f, yf - 1.0f, zf - 1.0f);

    // Trilinear interpolation
    float x1 = lerp(n000, n100, u);
    float x2 = lerp(n010, n110, u);
    float x3 = lerp(n001, n101, u);
    float x4 = lerp(n011, n111, u);

    float y1 = lerp(x1, x2, v);
    float y2 = lerp(x3, x4, v);

    return lerp(y1, y2, w);
}

/// @brief Fractional Brownian motion (fBm) noise
/// Multiple octaves of Perlin noise with decreasing amplitude and increasing frequency
inline float fBm3D(float x, float y, float z, const TurbulenceConfig& config) {
    float value = 0.0f;
    float amplitude = 1.0f;
    float frequency = 1.0f / config.outer_scale_M;
    float total_amplitude = 0.0f;

    for (uint32_t i = 0; i < config.octaves; ++i) {
        value += amplitude * perlin3D(x * frequency, y * frequency, z * frequency, config.seed + i);
        total_amplitude += amplitude;
        amplitude *= config.persistence;
        frequency *= config.lacunarity;
    }

    // Normalize to [-1, 1]
    return value / total_amplitude;
}

/// @brief Sample turbulent density perturbation at position
/// Returns multiplicative factor: ρ' = ρ × (1 + δ)
/// @param r Radius in M
/// @param theta Polar angle [rad]
/// @param phi Azimuthal angle [rad]
/// @param config Turbulence configuration
/// @return Density perturbation factor (1.0 = no perturbation)
inline float sampleDensityPerturbation(float r, float theta, float phi,
                                        const TurbulenceConfig& config) {
    if (!config.enabled || config.amplitude < 1e-6f) {
        return 1.0f;
    }

    // Convert spherical to Cartesian for noise sampling
    float sin_theta = std::sin(theta);
    float x = r * sin_theta * std::cos(phi);
    float y = r * sin_theta * std::sin(phi);
    float z = r * std::cos(theta);

    // Sample fBm noise
    float noise = fBm3D(x, y, z, config);

    // Scale to [1 - amplitude, 1 + amplitude]
    float perturbation = 1.0f + config.amplitude * noise;

    // Clamp to ensure ρ > 0
    return std::max(perturbation, 0.01f);
}

} // namespace TurbulenceNoise

} // namespace Sirius
