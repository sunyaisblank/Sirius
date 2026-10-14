#pragma once

// IMAX 70mm film post-processing pipeline (CPU). Ported from RDFL001A.h.
//
// Stages, in order: grain -> halation -> colour grade -> vignette -> bloom.
// Each stage operates on an RGBA float framebuffer in place. FilmConfig now
// comes from render_config.h; the legacy header also pulled in the OptiX launch
// header solely so test code could reach FilmParamsGPU, and that include is
// dropped here (OptiX is retired; the GPU film path returns via the Vulkan
// backend seam later).
//
// References: Kodak Vision3 technical data; IMAX film standards.

#include "sirius/render/render_config.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

namespace sirius::render {

// Applies authentic film characteristics to a rendered RGBA framebuffer.
class FilmPipeline {
  public:
    explicit FilmPipeline(const FilmConfig& config = FilmConfig::Interstellar())
        : config_(config) {}

    // Apply the full pipeline. frame_index seeds temporal effects (grain).
    void Apply(float* pixels, int width, int height, uint32_t frame_index) {
        if (!config_.enabled) return;

        if (config_.grain_enabled) {
            ApplyGrain(pixels, width, height, frame_index);
        }

        if (config_.halation_enabled) {
            ApplyHalation(pixels, width, height);
        }

        ApplyColorGrade(pixels, width, height);

        if (config_.vignette_enabled) {
            ApplyVignette(pixels, width, height);
        }

        if (config_.bloom_enabled) {
            ApplyBloom(pixels, width, height);
        }
    }

    // Film grain: signal-dependent noise, sigma^2 proportional to luminance.
    void ApplyGrain(float* pixels, int width, int height, uint32_t frame_seed) {
        uint32_t seed = frame_seed * 1664525u + 1013904223u;

        for (int y = 0; y < height; ++y) {
            for (int x = 0; x < width; ++x) {
                int idx = (y * width + x) * 4;

                float r = pixels[idx + 0];
                float g = pixels[idx + 1];
                float b = pixels[idx + 2];

                float L = 0.299f * r + 0.587f * g + 0.114f * b;
                float sigma = config_.grain_intensity * std::sqrt(std::max(L, 0.01f));

                seed = seed * 1664525u + 1013904223u;
                float nr = GaussianNoise(seed) * sigma;
                seed = seed * 1664525u + 1013904223u;
                float ng = GaussianNoise(seed) * sigma;
                seed = seed * 1664525u + 1013904223u;
                float nb = GaussianNoise(seed) * sigma;

                float corr = config_.grain_uniformity;
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

    // Halation: light scatter from the film base, red-biased highlight bloom.
    void ApplyHalation(float* pixels, int width, int height) {
        std::vector<float> highlight(width * height * 3, 0.0f);

        for (int y = 0; y < height; ++y) {
            for (int x = 0; x < width; ++x) {
                int idx = (y * width + x) * 4;
                float L =
                    0.299f * pixels[idx] + 0.587f * pixels[idx + 1] + 0.114f * pixels[idx + 2];

                if (L > config_.halation_threshold) {
                    float excess = L - config_.halation_threshold;
                    int hidx = (y * width + x) * 3;
                    highlight[hidx + 0] = excess * config_.halation_color_r;
                    highlight[hidx + 1] = excess * config_.halation_color_g;
                    highlight[hidx + 2] = excess * config_.halation_color_b;
                }
            }
        }

        int radius = static_cast<int>(config_.halation_radius);
        std::vector<float> blurred(width * height * 3, 0.0f);
        BoxBlur(highlight.data(), blurred.data(), width, height, 3, radius);

        for (int y = 0; y < height; ++y) {
            for (int x = 0; x < width; ++x) {
                int idx = (y * width + x) * 4;
                int hidx = (y * width + x) * 3;
                pixels[idx + 0] += blurred[hidx + 0] * config_.halation_strength;
                pixels[idx + 1] += blurred[hidx + 1] * config_.halation_strength;
                pixels[idx + 2] += blurred[hidx + 2] * config_.halation_strength;
            }
        }
    }

    // Colour grade: exposure, saturation, filmic S-curve, contrast.
    void ApplyColorGrade(float* pixels, int width, int height) {
        float exposure_mult = std::pow(2.0f, config_.exposure);

        for (int y = 0; y < height; ++y) {
            for (int x = 0; x < width; ++x) {
                int idx = (y * width + x) * 4;

                float r = pixels[idx + 0];
                float g = pixels[idx + 1];
                float b = pixels[idx + 2];

                r *= exposure_mult;
                g *= exposure_mult;
                b *= exposure_mult;

                float L = 0.299f * r + 0.587f * g + 0.114f * b;
                r = L + (r - L) * config_.saturation;
                g = L + (g - L) * config_.saturation;
                b = L + (b - L) * config_.saturation;

                r = FilmicCurve(r);
                g = FilmicCurve(g);
                b = FilmicCurve(b);

                r = ContrastCurve(r);
                g = ContrastCurve(g);
                b = ContrastCurve(b);

                pixels[idx + 0] = std::clamp(r, 0.0f, 1.0f);
                pixels[idx + 1] = std::clamp(g, 0.0f, 1.0f);
                pixels[idx + 2] = std::clamp(b, 0.0f, 1.0f);
            }
        }
    }

    // Vignette: corner darkening via a smoothstep falloff.
    void ApplyVignette(float* pixels, int width, int height) {
        float cx = width * 0.5f;
        float cy = height * 0.5f;

        for (int y = 0; y < height; ++y) {
            for (int x = 0; x < width; ++x) {
                int idx = (y * width + x) * 4;

                float dx = (x - cx) / cx;
                float dy = (y - cy) / cy;
                float dist = std::sqrt(dx * dx + dy * dy);

                float v = 1.0f - config_.vignette_strength *
                                     Smoothstep(config_.vignette_radius - config_.vignette_softness,
                                                config_.vignette_radius + config_.vignette_softness,
                                                dist);

                pixels[idx + 0] *= v;
                pixels[idx + 1] *= v;
                pixels[idx + 2] *= v;
            }
        }
    }

    // Bloom: soft glow extracted from bright areas, multi-pass blurred.
    void ApplyBloom(float* pixels, int width, int height) {
        std::vector<float> bright(width * height * 3, 0.0f);

        for (int y = 0; y < height; ++y) {
            for (int x = 0; x < width; ++x) {
                int idx = (y * width + x) * 4;
                int bidx = (y * width + x) * 3;

                float L = std::max({pixels[idx], pixels[idx + 1], pixels[idx + 2]});
                if (L > config_.bloom_threshold) {
                    float excess = L - config_.bloom_threshold;
                    bright[bidx + 0] = pixels[idx + 0] * excess;
                    bright[bidx + 1] = pixels[idx + 1] * excess;
                    bright[bidx + 2] = pixels[idx + 2] * excess;
                }
            }
        }

        std::vector<float> blurred(width * height * 3);
        std::copy(bright.begin(), bright.end(), blurred.begin());

        int radius = static_cast<int>(config_.bloom_radius);
        for (int i = 0; i < 3; ++i) {
            std::vector<float> temp(width * height * 3);
            BoxBlur(blurred.data(), temp.data(), width, height, 3, radius);
            std::copy(temp.begin(), temp.end(), blurred.begin());
        }

        for (int y = 0; y < height; ++y) {
            for (int x = 0; x < width; ++x) {
                int idx = (y * width + x) * 4;
                int bidx = (y * width + x) * 3;

                pixels[idx + 0] += blurred[bidx + 0] * config_.bloom_intensity;
                pixels[idx + 1] += blurred[bidx + 1] * config_.bloom_intensity;
                pixels[idx + 2] += blurred[bidx + 2] * config_.bloom_intensity;
            }
        }
    }

    const FilmConfig& Config() const { return config_; }
    void SetConfig(const FilmConfig& config) { config_ = config; }

  private:
    // Gaussian sample from a uniform seed (Box-Muller); u1 floored off 0 and 1.
    float GaussianNoise(uint32_t seed) const {
        float u1 = (seed & 0xFFFF) / 65536.0f + 1e-6f;
        u1 = std::clamp(u1, 1e-6f, 1.0f - 1e-6f);
        float u2 = ((seed >> 16) & 0xFFFF) / 65535.0f;
        return std::sqrt(-2.0f * std::log(u1)) * std::cos(6.28318f * u2);
    }

    // Filmic S-curve with toe and shoulder about the midtone point.
    float FilmicCurve(float x) const {
        if (x <= 0.0f) return 0.0f;
        if (x >= 1.0f) return 1.0f;

        float toe = config_.toe_strength;
        float shoulder = config_.shoulder_strength;
        float mid = std::clamp(config_.midtone_point, 0.01f, 0.99f);

        if (x < mid) {
            float t = x / mid;
            return mid * std::pow(t, 1.0f + toe);
        } else {
            float t = (x - mid) / (1.0f - mid);
            t = std::clamp(t, 0.0f, 1.0f);
            float y = 1.0f - std::pow(1.0f - t, 1.0f + shoulder);
            return mid + (1.0f - mid) * y;
        }
    }

    // Contrast about the 0.5 pivot.
    float ContrastCurve(float x) const {
        float mid = 0.5f;
        return mid + (x - mid) * config_.contrast;
    }

    // Smoothstep interpolation.
    float Smoothstep(float edge0, float edge1, float x) const {
        float t = std::clamp((x - edge0) / (edge1 - edge0), 0.0f, 1.0f);
        return t * t * (3.0f - 2.0f * t);
    }

    // Separable box blur (horizontal then vertical pass).
    void BoxBlur(const float* src, float* dst, int width, int height, int channels, int radius) {
        std::vector<float> temp(width * height * channels);

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

    FilmConfig config_;
};

}  // namespace sirius::render
