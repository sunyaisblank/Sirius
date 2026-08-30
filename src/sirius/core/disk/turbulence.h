#pragma once

// Deterministic procedural density perturbations for accretion-disk
// visualisation via fractional Brownian motion over 3D Perlin gradient noise.
// This is not a Kolmogorov, MRI, or GRMHD solution. Ported from PHTR001A.h.

#include "sirius/base/contracts.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>

namespace sirius::core {

// Procedural density-noise configuration.
struct TurbulenceConfig {
    friend bool operator==(const TurbulenceConfig&, const TurbulenceConfig&) = default;

    float outer_scale_M = 5.0f;  // Largest noise scale [M].
    float amplitude = 0.3f;      // Density fluctuation amplitude in [0, 1].
    std::uint32_t octaves = 6;   // Noise octaves (fractal detail levels).
    std::uint32_t seed = 12345;  // Random seed for reproducibility.
    float lacunarity = 2.0f;     // Frequency multiplier per octave.
    float persistence = 0.5f;    // Amplitude decay per octave.
};

[[nodiscard]] inline bool IsRepresentedTurbulenceConfig(const TurbulenceConfig& config) noexcept {
    return std::isfinite(config.outer_scale_M) && config.outer_scale_M >= 0.01f &&
           config.outer_scale_M <= 1.0e6f && std::isfinite(config.amplitude) &&
           config.amplitude >= 0.0f && config.amplitude <= 1.0f && config.octaves >= 1u &&
           config.octaves <= 8u && std::isfinite(config.lacunarity) && config.lacunarity >= 1.1f &&
           config.lacunarity <= 16.0f && std::isfinite(config.persistence) &&
           config.persistence >= 0.1f && config.persistence <= 0.9f;
}

// Perlin gradient noise for turbulence generation (CPU).
namespace turbulence_noise {

// Integer hash for gradient selection.
inline std::uint32_t Hash(std::uint32_t x, std::uint32_t y, std::uint32_t z, std::uint32_t seed) {
    std::uint32_t h = seed;
    h ^= x * 374761393u;
    h ^= y * 668265263u;
    h ^= z * 1013904223u;
    h = ((h >> 16) ^ h) * 2654435761u;
    h = ((h >> 16) ^ h) * 2654435761u;
    return h;
}

// Gradient vector from a hash (12 cube-edge directions).
inline void Gradient(std::uint32_t h, float& gx, float& gy, float& gz) {
    const float grads[12][3] = {{1, 1, 0}, {-1, 1, 0}, {1, -1, 0}, {-1, -1, 0},
                                {1, 0, 1}, {-1, 0, 1}, {1, 0, -1}, {-1, 0, -1},
                                {0, 1, 1}, {0, -1, 1}, {0, 1, -1}, {0, -1, -1}};
    int idx = h % 12;
    gx = grads[idx][0];
    gy = grads[idx][1];
    gz = grads[idx][2];
}

// Hermite quintic smoothstep.
inline float Fade(float t) { return t * t * t * (t * (t * 6.0f - 15.0f) + 10.0f); }

// Linear interpolation.
inline float Lerp(float a, float b, float t) { return a + t * (b - a); }

// 3D Perlin noise in [-1, 1].
inline float Perlin3D(float x, float y, float z, std::uint32_t seed) {
    constexpr float kMinimumCell = static_cast<float>(std::numeric_limits<int>::lowest() + 1);
    constexpr float kMaximumCell = static_cast<float>(std::numeric_limits<int>::max() - 1);
    SIRIUS_PRE(std::isfinite(x) && x >= kMinimumCell && x <= kMaximumCell);
    SIRIUS_PRE(std::isfinite(y) && y >= kMinimumCell && y <= kMaximumCell);
    SIRIUS_PRE(std::isfinite(z) && z >= kMinimumCell && z <= kMaximumCell);

    // Integer cell coordinates.
    int xi = static_cast<int>(std::floor(x));
    int yi = static_cast<int>(std::floor(y));
    int zi = static_cast<int>(std::floor(z));

    // Fractional coordinates within the cell.
    float xf = x - xi;
    float yf = y - yi;
    float zf = z - zi;

    // Fade curves.
    float u = Fade(xf);
    float v = Fade(yf);
    float w = Fade(zf);

    // Hash and gradient at each corner.
    float n000, n001, n010, n011, n100, n101, n110, n111;

    auto corner_noise = [&](int cx, int cy, int cz, float fx, float fy, float fz) {
        float gx, gy, gz;
        Gradient(Hash(xi + cx, yi + cy, zi + cz, seed), gx, gy, gz);
        return gx * fx + gy * fy + gz * fz;
    };

    n000 = corner_noise(0, 0, 0, xf, yf, zf);
    n100 = corner_noise(1, 0, 0, xf - 1.0f, yf, zf);
    n010 = corner_noise(0, 1, 0, xf, yf - 1.0f, zf);
    n110 = corner_noise(1, 1, 0, xf - 1.0f, yf - 1.0f, zf);
    n001 = corner_noise(0, 0, 1, xf, yf, zf - 1.0f);
    n101 = corner_noise(1, 0, 1, xf - 1.0f, yf, zf - 1.0f);
    n011 = corner_noise(0, 1, 1, xf, yf - 1.0f, zf - 1.0f);
    n111 = corner_noise(1, 1, 1, xf - 1.0f, yf - 1.0f, zf - 1.0f);

    // Trilinear interpolation.
    float x1 = Lerp(n000, n100, u);
    float x2 = Lerp(n010, n110, u);
    float x3 = Lerp(n001, n101, u);
    float x4 = Lerp(n011, n111, u);

    float y1 = Lerp(x1, x2, v);
    float y2 = Lerp(x3, x4, v);

    return Lerp(y1, y2, w);
}

// Fractional Brownian motion: octaves of Perlin noise with decreasing amplitude
// and increasing frequency, normalised to [-1, 1].
inline float Fbm3D(float x, float y, float z, const TurbulenceConfig& config) {
    SIRIUS_PRE(IsRepresentedTurbulenceConfig(config));
    SIRIUS_PRE(std::isfinite(x) && std::isfinite(y) && std::isfinite(z));
    float value = 0.0f;
    float amplitude = 1.0f;
    float frequency = 1.0f / config.outer_scale_M;
    float total_amplitude = 0.0f;

    for (std::uint32_t i = 0; i < config.octaves; ++i) {
        value += amplitude * Perlin3D(x * frequency, y * frequency, z * frequency, config.seed + i);
        total_amplitude += amplitude;
        amplitude *= config.persistence;
        frequency *= config.lacunarity;
    }

    return value / total_amplitude;
}

// Turbulent density perturbation factor at (r, theta, phi): rho' = rho (1 + delta),
// clamped to stay positive. Returns 1.0 for no perturbation.
inline float SampleDensityPerturbation(float r, float theta, float phi,
                                       const TurbulenceConfig& config) {
    SIRIUS_PRE(IsRepresentedTurbulenceConfig(config));
    SIRIUS_PRE(std::isfinite(r) && r >= 0.0f);
    SIRIUS_PRE(std::isfinite(theta) && theta >= 0.0f &&
               theta <= static_cast<float>(std::numbers::pi));
    SIRIUS_PRE(std::isfinite(phi));
    if (config.amplitude < 1e-6f) {
        return 1.0f;
    }

    // Spherical to Cartesian for noise sampling.
    float sin_theta = std::sin(theta);
    float x = r * sin_theta * std::cos(phi);
    float y = r * sin_theta * std::sin(phi);
    float z = r * std::cos(theta);

    float noise = Fbm3D(x, y, z, config);

    // Scale to [1 - amplitude, 1 + amplitude].
    float perturbation = 1.0f + config.amplitude * noise;

    // Clamp to keep rho > 0.
    return std::max(perturbation, 0.01f);
}

}  // namespace turbulence_noise

}  // namespace sirius::core
