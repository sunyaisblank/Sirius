#pragma once

// Post-processing operators for HDR image buffers: tonemapping, bloom, and
// colour grading, all on linear float RGBA. Reference: Academy Color Encoding
// System (ACES), Reinhard et al. (2002). Ported from PPOP001A.h.

#include "sirius/base/contracts.h"

#include <algorithm>
#include <cmath>
#include <vector>

namespace sirius::core {

// Tonemapping curve selector.
enum class TonemapType {
    None,      // No tonemapping (clamp only)
    Reinhard,  // Reinhard global operator
    Aces,      // ACES filmic curve (display-linear output; writer owns transfer)
    Filmic,    // Hable/Uncharted 2 filmic curve
    Exposure   // Simple exposure + gamma
};

// Tonemapping, bloom, and colour-grading parameters.
struct PostProcessConfig {
    // Tonemapping
    TonemapType tonemapper = TonemapType::Aces;
    float exposure = 1.0f;  // Exposure multiplier
    float gamma = 2.2f;     // Display gamma

    // Bloom
    bool enable_bloom = false;
    float bloom_intensity = 0.3f;
    float bloom_threshold = 0.8f;
    int bloom_radius = 8;

    // Colour grading
    float saturation = 1.0f;
    float contrast = 1.0f;
    float lift = 0.0f;  // Shadow lift
    float gain = 1.0f;  // Highlight gain
};

[[nodiscard]] inline bool IsRepresentedPostProcessConfig(const PostProcessConfig& config) noexcept {
    const bool known_tonemapper =
        config.tonemapper == TonemapType::None || config.tonemapper == TonemapType::Reinhard ||
        config.tonemapper == TonemapType::Aces || config.tonemapper == TonemapType::Filmic ||
        config.tonemapper == TonemapType::Exposure;
    if (!known_tonemapper || !std::isfinite(config.exposure) || config.exposure <= 0.0f ||
        config.exposure > 100.0f || !std::isfinite(config.gamma) || config.gamma <= 0.0f ||
        config.gamma > 10.0f || !std::isfinite(config.bloom_intensity) ||
        config.bloom_intensity < 0.0f || config.bloom_intensity > 5.0f ||
        !std::isfinite(config.bloom_threshold) || config.bloom_threshold < 0.0f ||
        config.bloom_threshold > 100.0f || config.bloom_radius < 0 || config.bloom_radius > 256 ||
        !std::isfinite(config.saturation) || config.saturation < 0.0f || config.saturation > 4.0f ||
        !std::isfinite(config.contrast) || config.contrast < 0.0f || config.contrast > 4.0f ||
        !std::isfinite(config.lift) || config.lift < -100.0f || config.lift > 100.0f ||
        !std::isfinite(config.gain) || config.gain < 0.0f || config.gain > 100.0f) {
        return false;
    }
    return config.enable_bloom || (config.bloom_intensity == PostProcessConfig{}.bloom_intensity &&
                                   config.bloom_threshold == PostProcessConfig{}.bloom_threshold &&
                                   config.bloom_radius == PostProcessConfig{}.bloom_radius);
}

// Scalar tonemapping curves and per-pixel dispatch.
namespace tonemap {

// ACES filmic curve (rational approximation), clamped to [0, 1].
inline float Aces(float x) {
    const float a = 2.51f;
    const float b = 0.03f;
    const float c = 2.43f;
    const float d = 0.59f;
    const float e = 0.14f;
    return std::clamp((x * (a * x + b)) / (x * (c * x + d) + e), 0.0f, 1.0f);
}

// Reinhard global operator x / (1 + x).
inline float Reinhard(float x) { return x / (1.0f + x); }

// Hable filmic curve (Uncharted 2), normalised by the white point.
inline float Filmic(float x) {
    const float A = 0.15f;
    const float B = 0.50f;
    const float C = 0.10f;
    const float D = 0.20f;
    const float E = 0.02f;
    const float F = 0.30f;
    auto curve = [=](float t) {
        return ((t * (A * t + C * B) + D * E) / (t * (A * t + B) + D * F)) - E / F;
    };
    const float white_point = 11.2f;
    return curve(x) / curve(white_point);
}

// Apply exposure then the selected tonemapping curve to an RGB triple.
inline void Apply(float& r, float& g, float& b, TonemapType type, float exposure) {
    SIRIUS_PRE(std::isfinite(r) && std::isfinite(g) && std::isfinite(b));
    SIRIUS_PRE(std::isfinite(exposure) && exposure > 0.0f && exposure <= 100.0f);
    SIRIUS_PRE(type == TonemapType::None || type == TonemapType::Reinhard ||
               type == TonemapType::Aces || type == TonemapType::Filmic ||
               type == TonemapType::Exposure);
    r *= exposure;
    g *= exposure;
    b *= exposure;

    switch (type) {
        case TonemapType::Aces:
            r = Aces(r);
            g = Aces(g);
            b = Aces(b);
            break;
        case TonemapType::Reinhard:
            r = Reinhard(r);
            g = Reinhard(g);
            b = Reinhard(b);
            break;
        case TonemapType::Filmic:
            r = Filmic(r);
            g = Filmic(g);
            b = Filmic(b);
            break;
        case TonemapType::Exposure:
            r = std::clamp(r, 0.0f, 1.0f);
            g = std::clamp(g, 0.0f, 1.0f);
            b = std::clamp(b, 0.0f, 1.0f);
            break;
        case TonemapType::None:
            r = std::clamp(r, 0.0f, 1.0f);
            g = std::clamp(g, 0.0f, 1.0f);
            b = std::clamp(b, 0.0f, 1.0f);
            break;
        default:
            break;
    }
}

}  // namespace tonemap

// Separable box-blur bloom over an RGBA float buffer.
class BloomFilter {
  public:
    // Add a bloom halo to an RGBA float buffer (4 floats per pixel) in place.
    static void Apply(std::vector<float>& buffer, int width, int height,
                      const PostProcessConfig& config) {
        SIRIUS_PRE(IsRepresentedPostProcessConfig(config));
        SIRIUS_PRE(width > 0 && height > 0);
        SIRIUS_PRE(buffer.size() ==
                   static_cast<std::size_t>(width) * static_cast<std::size_t>(height) * 4u);
        if (!config.enable_bloom || config.bloom_intensity <= 0.0f) {
            return;
        }

        // Extract bright pixels
        std::vector<float> bright(buffer.size(), 0.0f);
        for (int y = 0; y < height; ++y) {
            for (int x = 0; x < width; ++x) {
                int idx = (y * width + x) * 4;
                float luminance =
                    0.2126f * buffer[idx] + 0.7152f * buffer[idx + 1] + 0.0722f * buffer[idx + 2];

                if (luminance > config.bloom_threshold) {
                    float factor = (luminance - config.bloom_threshold) / luminance;
                    bright[idx] = buffer[idx] * factor;
                    bright[idx + 1] = buffer[idx + 1] * factor;
                    bright[idx + 2] = buffer[idx + 2] * factor;
                }
            }
        }

        // Simple box blur for bloom (separable)
        std::vector<float> blurred(buffer.size(), 0.0f);
        int radius = config.bloom_radius;

        // Horizontal pass
        for (int y = 0; y < height; ++y) {
            for (int x = 0; x < width; ++x) {
                float sum_r = 0, sum_g = 0, sum_b = 0;
                int count = 0;
                for (int dx = -radius; dx <= radius; ++dx) {
                    int nx = std::clamp(x + dx, 0, width - 1);
                    int idx = (y * width + nx) * 4;
                    sum_r += bright[idx];
                    sum_g += bright[idx + 1];
                    sum_b += bright[idx + 2];
                    count++;
                }
                int idx = (y * width + x) * 4;
                blurred[idx] = sum_r / count;
                blurred[idx + 1] = sum_g / count;
                blurred[idx + 2] = sum_b / count;
            }
        }

        // Vertical pass
        bright = blurred;
        for (int y = 0; y < height; ++y) {
            for (int x = 0; x < width; ++x) {
                float sum_r = 0, sum_g = 0, sum_b = 0;
                int count = 0;
                for (int dy = -radius; dy <= radius; ++dy) {
                    int ny = std::clamp(y + dy, 0, height - 1);
                    int idx = (ny * width + x) * 4;
                    sum_r += bright[idx];
                    sum_g += bright[idx + 1];
                    sum_b += bright[idx + 2];
                    count++;
                }
                int idx = (y * width + x) * 4;
                blurred[idx] = sum_r / count;
                blurred[idx + 1] = sum_g / count;
                blurred[idx + 2] = sum_b / count;
            }
        }

        // Add bloom to original
        for (std::size_t i = 0; i < buffer.size(); i += 4) {
            buffer[i] += blurred[i] * config.bloom_intensity;
            buffer[i + 1] += blurred[i + 1] * config.bloom_intensity;
            buffer[i + 2] += blurred[i + 2] * config.bloom_intensity;
        }
    }
};

// Lift/gain/contrast/saturation colour grading on an RGB triple.
class ColourGrading {
  public:
    // Grade an RGB triple (linear, 0-1) in place, clamped to [0, 1].
    static void Apply(float& r, float& g, float& b, const PostProcessConfig& config) {
        SIRIUS_PRE(IsRepresentedPostProcessConfig(config));
        SIRIUS_PRE(std::isfinite(r) && std::isfinite(g) && std::isfinite(b));
        // Lift (shadows)
        r += config.lift;
        g += config.lift;
        b += config.lift;

        // Gain (highlights)
        r *= config.gain;
        g *= config.gain;
        b *= config.gain;

        // Contrast (pivot at 0.5)
        r = (r - 0.5f) * config.contrast + 0.5f;
        g = (g - 0.5f) * config.contrast + 0.5f;
        b = (b - 0.5f) * config.contrast + 0.5f;

        // Saturation
        float lum = 0.2126f * r + 0.7152f * g + 0.0722f * b;
        r = lum + (r - lum) * config.saturation;
        g = lum + (g - lum) * config.saturation;
        b = lum + (b - lum) * config.saturation;

        // Clamp
        r = std::clamp(r, 0.0f, 1.0f);
        g = std::clamp(g, 0.0f, 1.0f);
        b = std::clamp(b, 0.0f, 1.0f);
    }
};

// Full post-processing pipeline: bloom, then per-pixel tonemap, grade, gamma.
class PostProcessor {
  public:
    // Run the pipeline over an RGBA HDR-linear buffer in place.
    static void Process(std::vector<float>& buffer, int width, int height,
                        const PostProcessConfig& config) {
        SIRIUS_PRE(IsRepresentedPostProcessConfig(config));
        SIRIUS_PRE(width > 0 && height > 0);
        SIRIUS_PRE(buffer.size() ==
                   static_cast<std::size_t>(width) * static_cast<std::size_t>(height) * 4u);
        // 1. Bloom (operates on HDR values)
        BloomFilter::Apply(buffer, width, height, config);

        // 2. Tonemapping and colour grading per pixel
        for (int y = 0; y < height; ++y) {
            for (int x = 0; x < width; ++x) {
                int idx = (y * width + x) * 4;

                float r = buffer[idx];
                float g = buffer[idx + 1];
                float b = buffer[idx + 2];

                // Tonemapping
                tonemap::Apply(r, g, b, config.tonemapper, config.exposure);

                // Colour grading
                ColourGrading::Apply(r, g, b, config);

                // Gamma correction
                r = std::pow(r, 1.0f / config.gamma);
                g = std::pow(g, 1.0f / config.gamma);
                b = std::pow(b, 1.0f / config.gamma);

                buffer[idx] = r;
                buffer[idx + 1] = g;
                buffer[idx + 2] = b;
            }
        }
    }
};

}  // namespace sirius::core
