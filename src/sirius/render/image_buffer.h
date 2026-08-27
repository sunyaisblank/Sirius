#pragma once

// Floating-point HDR image buffers (RGB and RGBA interleaved) plus a lightweight
// spectral-radiance helper for blackbody colour. Ported from OUIB001A.h; the
// arithmetic (sRGB transfer, CIE Gaussian blackbody fit) is bit-identical.
//
// This buffer holds linear radiance. The transfer encode is applied once by the
// owning writer (PNG applies sRGB, EXR stays linear), never here; ToSrgb8 exists
// for callers that ask for an 8-bit view explicitly.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

namespace sirius::render {

// Lightweight linear RGB radiance with a self-contained blackbody constructor;
// distinct from sirius::core::spectral::Rgb, kept local to the writer layer so
// the output path carries no core dependency.
struct SpectralRadiance {
    float r = 0.0f;
    float g = 0.0f;
    float b = 0.0f;

    SpectralRadiance() = default;
    SpectralRadiance(float r_, float g_, float b_) : r(r_), g(g_), b(b_) {}

    // Blackbody radiance for a temperature (K), integrated over the visible band
    // with a Gaussian CIE 1931 fit and normalised so the brightest channel is 1.
    static SpectralRadiance Blackbody(double temperature_k) {
        if (temperature_k <= 0) return SpectralRadiance(0, 0, 0);

        constexpr double h_planck = 6.62607015e-34;
        constexpr double c_light = 2.99792458e8;
        constexpr double k_boltzmann = 1.380649e-23;
        constexpr double planck_c1 = 2.0 * h_planck * c_light * c_light;
        constexpr double planck_c2 = h_planck * c_light / k_boltzmann;

        constexpr int n_samples = 32;
        constexpr double lambda_min = 380e-9;
        constexpr double lambda_max = 780e-9;
        constexpr double d_lambda = (lambda_max - lambda_min) / n_samples;

        double X = 0, Y = 0, Z = 0;

        for (int i = 0; i < n_samples; i++) {
            double lambda = lambda_min + (i + 0.5) * d_lambda;
            double lambda_nm = lambda * 1e9;

            // Planck radiance; the x > 700 guard avoids overflow in exp.
            double x = planck_c2 / (lambda * temperature_k);
            double radiance =
                (x > 700) ? 0.0 : planck_c1 / (std::pow(lambda, 5) * (std::exp(x) - 1.0));

            double t1, t2, t3;

            // CIE X.
            t1 = (lambda_nm - 442.0) * ((lambda_nm < 442.0) ? 0.0624 : 0.0374);
            t2 = (lambda_nm - 599.8) * ((lambda_nm < 599.8) ? 0.0264 : 0.0323);
            t3 = (lambda_nm - 501.1) * ((lambda_nm < 501.1) ? 0.0490 : 0.0382);
            double cie_x = 0.362 * std::exp(-0.5 * t1 * t1) + 1.056 * std::exp(-0.5 * t2 * t2) -
                           0.065 * std::exp(-0.5 * t3 * t3);

            // CIE Y.
            t1 = (lambda_nm - 568.8) * ((lambda_nm < 568.8) ? 0.0213 : 0.0247);
            t2 = (lambda_nm - 530.9) * ((lambda_nm < 530.9) ? 0.0613 : 0.0322);
            double cie_y = 0.821 * std::exp(-0.5 * t1 * t1) + 0.286 * std::exp(-0.5 * t2 * t2);

            // CIE Z.
            t1 = (lambda_nm - 437.0) * ((lambda_nm < 437.0) ? 0.0845 : 0.0278);
            t2 = (lambda_nm - 459.0) * ((lambda_nm < 459.0) ? 0.0385 : 0.0725);
            double cie_z = 1.217 * std::exp(-0.5 * t1 * t1) + 0.681 * std::exp(-0.5 * t2 * t2);

            X += cie_x * radiance * d_lambda;
            Y += cie_y * radiance * d_lambda;
            Z += cie_z * radiance * d_lambda;
        }

        // XYZ to linear sRGB.
        float rf = static_cast<float>(3.2404542 * X - 1.5371385 * Y - 0.4985314 * Z);
        float gf = static_cast<float>(-0.9692660 * X + 1.8760108 * Y + 0.0415560 * Z);
        float bf = static_cast<float>(0.0556434 * X - 0.2040259 * Y + 1.0572252 * Z);

        float max_val = std::max({rf, gf, bf, 0.001f});
        return SpectralRadiance(rf / max_val, gf / max_val, bf / max_val);
    }

    // Approximate redshift recolour by factor z (z > 0 reddens, z < 0 blueshifts).
    SpectralRadiance Redshifted(double z) const {
        if (std::abs(z) < 0.001) return *this;

        float factor = static_cast<float>(1.0 / (1.0 + z));
        if (z > 0) {
            return SpectralRadiance(r, g * factor, b * factor * factor);
        } else {
            float inv = static_cast<float>(1.0 + z);
            return SpectralRadiance(r * inv * inv, g * inv, b);
        }
    }
};

// Floating-point HDR image buffer, RGB interleaved (3 floats per pixel).
class ImageBuffer {
  public:
    int width = 0;
    int height = 0;
    std::vector<float> pixels;  // RGB interleaved (3 floats per pixel).

    ImageBuffer() = default;

    ImageBuffer(int w, int h) { Allocate(w, h); }

    void Allocate(int w, int h) {
        if (w <= 0 || h <= 0) {
            width = 0;
            height = 0;
            pixels.clear();
            return;
        }
        width = w;
        height = h;
        pixels.resize(static_cast<std::size_t>(w) * h * 3, 0.0f);
    }

    void Clear() { std::fill(pixels.begin(), pixels.end(), 0.0f); }

    std::size_t PixelCount() const {
        return width > 0 && height > 0 ? static_cast<std::size_t>(width) * height : 0;
    }
    std::size_t BufferSize() const { return pixels.size(); }
    bool HasValidShape() const { return PixelCount() > 0 && pixels.size() == PixelCount() * 3; }

    void SetPixel(int x, int y, float r, float g, float b) {
        if (!HasValidShape() || x < 0 || x >= width || y < 0 || y >= height) return;
        std::size_t idx = (static_cast<std::size_t>(y) * width + x) * 3;
        pixels[idx + 0] = r;
        pixels[idx + 1] = g;
        pixels[idx + 2] = b;
    }

    void GetPixel(int x, int y, float& r, float& g, float& b) const {
        if (!HasValidShape() || x < 0 || x >= width || y < 0 || y >= height) {
            r = g = b = 0.0f;
            return;
        }
        std::size_t idx = (static_cast<std::size_t>(y) * width + x) * 3;
        r = pixels[idx + 0];
        g = pixels[idx + 1];
        b = pixels[idx + 2];
    }

    void SetPixelFromSpectral(int x, int y, const SpectralRadiance& spec) {
        SetPixel(x, y, spec.r, spec.g, spec.b);
    }

    // Convert to 8-bit sRGB (gamma corrected); the display transfer applied once.
    std::vector<std::uint8_t> ToSrgb8() const {
        if (!HasValidShape()) return {};
        std::vector<std::uint8_t> result(PixelCount() * 3);

        for (std::size_t i = 0; i < PixelCount(); ++i) {
            for (int c = 0; c < 3; ++c) {
                float linear = pixels[i * 3 + c];
                float clamped = std::clamp(linear, 0.0f, 1.0f);

                float srgb;
                if (clamped <= 0.0031308f) {
                    srgb = 12.92f * clamped;
                } else {
                    srgb = 1.055f * std::pow(clamped, 1.0f / 2.4f) - 0.055f;
                }

                result[i * 3 + c] = static_cast<std::uint8_t>(srgb * 255.0f + 0.5f);
            }
        }

        return result;
    }

    // Convert to RGBA float with alpha = 1.0.
    std::vector<float> ToRgba() const {
        if (!HasValidShape()) return {};
        std::vector<float> result(PixelCount() * 4);

        for (std::size_t i = 0; i < PixelCount(); ++i) {
            result[i * 4 + 0] = pixels[i * 3 + 0];
            result[i * 4 + 1] = pixels[i * 3 + 1];
            result[i * 4 + 2] = pixels[i * 3 + 2];
            result[i * 4 + 3] = 1.0f;
        }

        return result;
    }
};

// Floating-point HDR image buffer, RGBA interleaved (4 floats per pixel).
class ImageBufferRGBA {
  public:
    int width = 0;
    int height = 0;
    std::vector<float> pixels;  // RGBA interleaved (4 floats per pixel).

    ImageBufferRGBA() = default;

    ImageBufferRGBA(int w, int h) { Allocate(w, h); }

    void Allocate(int w, int h) {
        if (w <= 0 || h <= 0) {
            width = 0;
            height = 0;
            pixels.clear();
            return;
        }
        width = w;
        height = h;
        pixels.resize(static_cast<std::size_t>(w) * h * 4, 0.0f);
    }

    void SetPixel(int x, int y, float r, float g, float b, float a = 1.0f) {
        if (!HasValidShape() || x < 0 || x >= width || y < 0 || y >= height) return;
        std::size_t idx = (static_cast<std::size_t>(y) * width + x) * 4;
        pixels[idx + 0] = r;
        pixels[idx + 1] = g;
        pixels[idx + 2] = b;
        pixels[idx + 3] = a;
    }

    std::size_t PixelCount() const {
        return width > 0 && height > 0 ? static_cast<std::size_t>(width) * height : 0;
    }
    bool HasValidShape() const { return PixelCount() > 0 && pixels.size() == PixelCount() * 4; }
};

}  // namespace sirius::render
