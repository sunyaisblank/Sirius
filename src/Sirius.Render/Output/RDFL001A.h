// RDFL001A.h - IMAX 70mm Film Post-Processing Pipeline
// Component ID: RDFL001A (Render/Film/FilmPipeline)
//
// Post-processing pipeline implementing authentic film characteristics
// for IMAX 70mm 15-perf cinematic output.
//
// PIPELINE STAGES:
// ================
// 1. Grain → 2. Halation → 3. Color Grade → 4. Vignette → 5. Gate Weave
//
// Each stage is implemented as a GPU-accelerated pass that can be
// applied to the rendered framebuffer.
//
// REFERENCES:
// - Kodak Vision3 Technical Data
// - IMAX Corporation Film Standards
// - Interstellar VFX breakdown (Christopher Nolan, 2014)

#pragma once

#include "../Acceleration/OptiX/RDOP003A.h"
#include "../../Sirius.Infrastructure/Configuration/CRFM001A.h"
#include <cmath>
#include <cstdint>
#include <vector>
#include <algorithm>

namespace Sirius {

//==============================================================================
// Film Pipeline Interface
//==============================================================================
class FilmPipeline {
public:
    explicit FilmPipeline(const FilmConfig& config = FilmConfig::Interstellar())
        : m_Config(config) {}

    /// @brief Apply full film pipeline to framebuffer (CPU version)
    /// @param pixels RGBA float buffer (modified in place)
    /// @param width Image width
    /// @param height Image height
    /// @param frame_index Frame number (for temporal effects)
    void apply(float* pixels, int width, int height, uint32_t frame_index) {
        if (!m_Config.enabled) return;

        // Apply in order: grain, halation, color grade, vignette
        if (m_Config.grain_enabled) {
            applyGrain(pixels, width, height, frame_index);
        }

        if (m_Config.halation_enabled) {
            applyHalation(pixels, width, height);
        }

        applyColorGrade(pixels, width, height);

        if (m_Config.vignette_enabled) {
            applyVignette(pixels, width, height);
        }

        if (m_Config.bloom_enabled) {
            applyBloom(pixels, width, height);
        }
    }

    // =========================================================================
    // Individual Pipeline Stages (CPU implementations)
    // =========================================================================

    /// @brief Apply film grain
    /// Signal-dependent noise: σ² ∝ L (Poisson-like)
    void applyGrain(float* pixels, int width, int height, uint32_t frame_seed) {
        uint32_t seed = frame_seed * 1664525u + 1013904223u;

        for (int y = 0; y < height; ++y) {
            for (int x = 0; x < width; ++x) {
                int idx = (y * width + x) * 4;

                float r = pixels[idx + 0];
                float g = pixels[idx + 1];
                float b = pixels[idx + 2];

                // Luminance for signal-dependent noise
                float L = 0.299f * r + 0.587f * g + 0.114f * b;
                float sigma = m_Config.grain_intensity * std::sqrt(std::max(L, 0.01f));

                // Generate noise
                seed = seed * 1664525u + 1013904223u;
                float nr = gaussianNoise(seed) * sigma;
                seed = seed * 1664525u + 1013904223u;
                float ng = gaussianNoise(seed) * sigma;
                seed = seed * 1664525u + 1013904223u;
                float nb = gaussianNoise(seed) * sigma;

                // Correlation between channels
                float corr = m_Config.grain_uniformity;
                float common = (nr + ng + nb) / 3.0f;
                nr = nr * (1.0f - corr) + common * corr;
                ng = ng * (1.0f - corr) + common * corr;
                nb = nb * (1.0f - corr) + common * corr;

                pixels[idx + 0] = std::max(r + nr, 0.0f);
                pixels[idx + 1] = std::max(g + ng, 0.0f);
                pixels[idx + 2] = std::max(b + nb, 0.0f);
            }
        }
    }

    /// @brief Apply halation (light scatter from film base)
    void applyHalation(float* pixels, int width, int height) {
        // Create highlight mask
        std::vector<float> highlight(width * height * 3, 0.0f);

        for (int y = 0; y < height; ++y) {
            for (int x = 0; x < width; ++x) {
                int idx = (y * width + x) * 4;
                float L = 0.299f * pixels[idx] + 0.587f * pixels[idx + 1] + 0.114f * pixels[idx + 2];

                if (L > m_Config.halation_threshold) {
                    float excess = L - m_Config.halation_threshold;
                    int hidx = (y * width + x) * 3;
                    highlight[hidx + 0] = excess * m_Config.halation_color_r;
                    highlight[hidx + 1] = excess * m_Config.halation_color_g;
                    highlight[hidx + 2] = excess * m_Config.halation_color_b;
                }
            }
        }

        // Blur the highlight
        int radius = static_cast<int>(m_Config.halation_radius);
        std::vector<float> blurred(width * height * 3, 0.0f);
        boxBlur(highlight.data(), blurred.data(), width, height, 3, radius);

        // Add back to image
        for (int y = 0; y < height; ++y) {
            for (int x = 0; x < width; ++x) {
                int idx = (y * width + x) * 4;
                int hidx = (y * width + x) * 3;
                pixels[idx + 0] += blurred[hidx + 0] * m_Config.halation_strength;
                pixels[idx + 1] += blurred[hidx + 1] * m_Config.halation_strength;
                pixels[idx + 2] += blurred[hidx + 2] * m_Config.halation_strength;
            }
        }
    }

    /// @brief Apply color grading (exposure, contrast, saturation, S-curve)
    void applyColorGrade(float* pixels, int width, int height) {
        float exposure_mult = std::pow(2.0f, m_Config.exposure);

        for (int y = 0; y < height; ++y) {
            for (int x = 0; x < width; ++x) {
                int idx = (y * width + x) * 4;

                float r = pixels[idx + 0];
                float g = pixels[idx + 1];
                float b = pixels[idx + 2];

                // Exposure
                r *= exposure_mult;
                g *= exposure_mult;
                b *= exposure_mult;

                // Saturation (in HSL-like space)
                float L = 0.299f * r + 0.587f * g + 0.114f * b;
                r = L + (r - L) * m_Config.saturation;
                g = L + (g - L) * m_Config.saturation;
                b = L + (b - L) * m_Config.saturation;

                // S-curve (filmic tone mapping)
                r = filmicCurve(r);
                g = filmicCurve(g);
                b = filmicCurve(b);

                // Contrast around midpoint
                r = contrastCurve(r);
                g = contrastCurve(g);
                b = contrastCurve(b);

                pixels[idx + 0] = std::clamp(r, 0.0f, 1.0f);
                pixels[idx + 1] = std::clamp(g, 0.0f, 1.0f);
                pixels[idx + 2] = std::clamp(b, 0.0f, 1.0f);
            }
        }
    }

    /// @brief Apply vignette (corner darkening)
    void applyVignette(float* pixels, int width, int height) {
        float cx = width * 0.5f;
        float cy = height * 0.5f;
        float max_dist = std::sqrt(cx * cx + cy * cy);

        for (int y = 0; y < height; ++y) {
            for (int x = 0; x < width; ++x) {
                int idx = (y * width + x) * 4;

                float dx = (x - cx) / cx;
                float dy = (y - cy) / cy;
                float dist = std::sqrt(dx * dx + dy * dy);

                // Vignette falloff
                float v = 1.0f - m_Config.vignette_strength *
                          smoothstep(m_Config.vignette_radius - m_Config.vignette_softness,
                                    m_Config.vignette_radius + m_Config.vignette_softness,
                                    dist);

                pixels[idx + 0] *= v;
                pixels[idx + 1] *= v;
                pixels[idx + 2] *= v;
            }
        }
    }

    /// @brief Apply bloom (soft glow on bright areas)
    void applyBloom(float* pixels, int width, int height) {
        // Extract bright areas
        std::vector<float> bright(width * height * 3, 0.0f);

        for (int y = 0; y < height; ++y) {
            for (int x = 0; x < width; ++x) {
                int idx = (y * width + x) * 4;
                int bidx = (y * width + x) * 3;

                float L = std::max({pixels[idx], pixels[idx + 1], pixels[idx + 2]});
                if (L > m_Config.bloom_threshold) {
                    float excess = L - m_Config.bloom_threshold;
                    bright[bidx + 0] = pixels[idx + 0] * excess;
                    bright[bidx + 1] = pixels[idx + 1] * excess;
                    bright[bidx + 2] = pixels[idx + 2] * excess;
                }
            }
        }

        // Multi-pass blur
        std::vector<float> blurred(width * height * 3);
        std::copy(bright.begin(), bright.end(), blurred.begin());

        int radius = static_cast<int>(m_Config.bloom_radius);
        for (int i = 0; i < 3; ++i) {
            std::vector<float> temp(width * height * 3);
            boxBlur(blurred.data(), temp.data(), width, height, 3, radius);
            std::copy(temp.begin(), temp.end(), blurred.begin());
        }

        // Additive blend
        for (int y = 0; y < height; ++y) {
            for (int x = 0; x < width; ++x) {
                int idx = (y * width + x) * 4;
                int bidx = (y * width + x) * 3;

                pixels[idx + 0] += blurred[bidx + 0] * m_Config.bloom_intensity;
                pixels[idx + 1] += blurred[bidx + 1] * m_Config.bloom_intensity;
                pixels[idx + 2] += blurred[bidx + 2] * m_Config.bloom_intensity;
            }
        }
    }

    // =========================================================================
    // Configuration
    // =========================================================================
    const FilmConfig& config() const { return m_Config; }
    void setConfig(const FilmConfig& config) { m_Config = config; }

private:
    // =========================================================================
    // Helper Functions
    // =========================================================================

    /// @brief Gaussian random number from uniform seed (Box-Muller transform)
    float gaussianNoise(uint32_t seed) const {
        // Clamp u1 to (0, 1) exclusive to avoid log(0) = -inf and log(1+ε) causing NaN
        float u1 = (seed & 0xFFFF) / 65536.0f + 1e-6f;  // Use 65536 to keep in (0, 1)
        u1 = std::clamp(u1, 1e-6f, 1.0f - 1e-6f);
        float u2 = ((seed >> 16) & 0xFFFF) / 65535.0f;
        return std::sqrt(-2.0f * std::log(u1)) * std::cos(6.28318f * u2);
    }

    /// @brief Filmic S-curve (toe + shoulder)
    float filmicCurve(float x) const {
        // Handle negative and extreme values
        if (x <= 0.0f) return 0.0f;
        if (x >= 1.0f) return 1.0f;

        // Reinhard-style with shoulder
        float toe = m_Config.toe_strength;
        float shoulder = m_Config.shoulder_strength;
        float mid = std::clamp(m_Config.midtone_point, 0.01f, 0.99f);

        // Piecewise curve
        if (x < mid) {
            // Toe region
            float t = x / mid;
            return mid * std::pow(t, 1.0f + toe);
        } else {
            // Shoulder region
            float t = (x - mid) / (1.0f - mid);
            t = std::clamp(t, 0.0f, 1.0f);
            float y = 1.0f - std::pow(1.0f - t, 1.0f + shoulder);
            return mid + (1.0f - mid) * y;
        }
    }

    /// @brief Contrast curve around midpoint
    float contrastCurve(float x) const {
        float mid = 0.5f;
        return mid + (x - mid) * m_Config.contrast;
    }

    /// @brief Smoothstep interpolation
    float smoothstep(float edge0, float edge1, float x) const {
        float t = std::clamp((x - edge0) / (edge1 - edge0), 0.0f, 1.0f);
        return t * t * (3.0f - 2.0f * t);
    }

    /// @brief Simple box blur (separable)
    void boxBlur(const float* src, float* dst, int width, int height, int channels, int radius) {
        std::vector<float> temp(width * height * channels);

        // Horizontal pass
        for (int y = 0; y < height; ++y) {
            for (int x = 0; x < width; ++x) {
                for (int c = 0; c < channels; ++c) {
                    float sum = 0.0f;
                    int count = 0;
                    for (int dx = -radius; dx <= radius; ++dx) {
                        int sx = std::clamp(x + dx, 0, width - 1);
                        sum += src[(y * width + sx) * channels + c];
                        count++;
                    }
                    temp[(y * width + x) * channels + c] = sum / count;
                }
            }
        }

        // Vertical pass
        for (int y = 0; y < height; ++y) {
            for (int x = 0; x < width; ++x) {
                for (int c = 0; c < channels; ++c) {
                    float sum = 0.0f;
                    int count = 0;
                    for (int dy = -radius; dy <= radius; ++dy) {
                        int sy = std::clamp(y + dy, 0, height - 1);
                        sum += temp[(sy * width + x) * channels + c];
                        count++;
                    }
                    dst[(y * width + x) * channels + c] = sum / count;
                }
            }
        }
    }

    FilmConfig m_Config;
};

} // namespace Sirius
