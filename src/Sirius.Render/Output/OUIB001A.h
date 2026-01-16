// OUIB001A.h - Image Buffer
// Component ID: OUIB001A
// Floating-point image buffer for HDR rendering
//
// Governance: docs/specification.md

#ifndef SIRIUS_RENDER_OUIB001A_H
#define SIRIUS_RENDER_OUIB001A_H

#include <vector>
#include <cstdint>
#include <algorithm>
#include <cmath>
#include <cstring>

// Forward declare spectral types if available
namespace Sirius { namespace Spectral { struct RGB; } }

namespace sirius::render {

//==============================================================================
// SpectralRadiance
// Lightweight spectral radiance representation
//==============================================================================

struct SpectralRadiance {
    float r = 0.0f;
    float g = 0.0f;
    float b = 0.0f;

    SpectralRadiance() = default;
    SpectralRadiance(float r_, float g_, float b_) : r(r_), g(g_), b(b_) {}

    //--------------------------------------------------------------------------
    // Blackbody Emission
    //--------------------------------------------------------------------------

    /// Create blackbody radiance for given temperature
    static SpectralRadiance blackbody(double temperatureK) {
        if (temperatureK <= 0) return SpectralRadiance(0, 0, 0);

        // Physical constants
        constexpr double h_PLANCK = 6.62607015e-34;
        constexpr double c_LIGHT = 2.99792458e8;
        constexpr double k_BOLTZMANN = 1.380649e-23;
        constexpr double PLANCK_C1 = 2.0 * h_PLANCK * c_LIGHT * c_LIGHT;
        constexpr double PLANCK_C2 = h_PLANCK * c_LIGHT / k_BOLTZMANN;

        // Integrate over visible spectrum
        constexpr int N_SAMPLES = 32;
        constexpr double LAMBDA_MIN = 380e-9;
        constexpr double LAMBDA_MAX = 780e-9;
        constexpr double D_LAMBDA = (LAMBDA_MAX - LAMBDA_MIN) / N_SAMPLES;

        double X = 0, Y = 0, Z = 0;

        for (int i = 0; i < N_SAMPLES; i++) {
            double lambda = LAMBDA_MIN + (i + 0.5) * D_LAMBDA;
            double lambda_nm = lambda * 1e9;

            // Planck radiance
            double x = PLANCK_C2 / (lambda * temperatureK);
            double radiance = (x > 700) ? 0.0 :
                PLANCK_C1 / (std::pow(lambda, 5) * (std::exp(x) - 1.0));

            // CIE color matching (Gaussian approximation)
            double t1, t2, t3;

            // X
            t1 = (lambda_nm - 442.0) * ((lambda_nm < 442.0) ? 0.0624 : 0.0374);
            t2 = (lambda_nm - 599.8) * ((lambda_nm < 599.8) ? 0.0264 : 0.0323);
            t3 = (lambda_nm - 501.1) * ((lambda_nm < 501.1) ? 0.0490 : 0.0382);
            double cie_x = 0.362 * std::exp(-0.5 * t1 * t1)
                         + 1.056 * std::exp(-0.5 * t2 * t2)
                         - 0.065 * std::exp(-0.5 * t3 * t3);

            // Y
            t1 = (lambda_nm - 568.8) * ((lambda_nm < 568.8) ? 0.0213 : 0.0247);
            t2 = (lambda_nm - 530.9) * ((lambda_nm < 530.9) ? 0.0613 : 0.0322);
            double cie_y = 0.821 * std::exp(-0.5 * t1 * t1)
                         + 0.286 * std::exp(-0.5 * t2 * t2);

            // Z
            t1 = (lambda_nm - 437.0) * ((lambda_nm < 437.0) ? 0.0845 : 0.0278);
            t2 = (lambda_nm - 459.0) * ((lambda_nm < 459.0) ? 0.0385 : 0.0725);
            double cie_z = 1.217 * std::exp(-0.5 * t1 * t1)
                         + 0.681 * std::exp(-0.5 * t2 * t2);

            // Accumulate
            X += cie_x * radiance * D_LAMBDA;
            Y += cie_y * radiance * D_LAMBDA;
            Z += cie_z * radiance * D_LAMBDA;
        }

        // XYZ to linear sRGB
        float rf = static_cast<float>( 3.2404542 * X - 1.5371385 * Y - 0.4985314 * Z);
        float gf = static_cast<float>(-0.9692660 * X + 1.8760108 * Y + 0.0415560 * Z);
        float bf = static_cast<float>( 0.0556434 * X - 0.2040259 * Y + 1.0572252 * Z);

        // Normalize to max = 1
        float maxVal = std::max({rf, gf, bf, 0.001f});
        return SpectralRadiance(rf / maxVal, gf / maxVal, bf / maxVal);
    }

    /// Apply redshift factor
    SpectralRadiance redshifted(double z) const {
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

//==============================================================================
// ImageBuffer
// Floating-point HDR image buffer (RGB)
//==============================================================================

class ImageBuffer {
public:
    int width = 0;
    int height = 0;
    std::vector<float> pixels;  // RGB interleaved (3 floats per pixel)

    ImageBuffer() = default;

    ImageBuffer(int w, int h) {
        allocate(w, h);
    }

    //--------------------------------------------------------------------------
    // Allocation
    //--------------------------------------------------------------------------

    void allocate(int w, int h) {
        width = w;
        height = h;
        pixels.resize(static_cast<size_t>(w) * h * 3, 0.0f);
    }

    void clear() {
        std::fill(pixels.begin(), pixels.end(), 0.0f);
    }

    size_t pixelCount() const { return static_cast<size_t>(width) * height; }
    size_t bufferSize() const { return pixels.size(); }

    //--------------------------------------------------------------------------
    // Pixel Access
    //--------------------------------------------------------------------------

    void setPixel(int x, int y, float r, float g, float b) {
        if (x < 0 || x >= width || y < 0 || y >= height) return;
        size_t idx = (static_cast<size_t>(y) * width + x) * 3;
        pixels[idx + 0] = r;
        pixels[idx + 1] = g;
        pixels[idx + 2] = b;
    }

    void getPixel(int x, int y, float& r, float& g, float& b) const {
        if (x < 0 || x >= width || y < 0 || y >= height) {
            r = g = b = 0.0f;
            return;
        }
        size_t idx = (static_cast<size_t>(y) * width + x) * 3;
        r = pixels[idx + 0];
        g = pixels[idx + 1];
        b = pixels[idx + 2];
    }

    void setPixelFromSpectral(int x, int y, const SpectralRadiance& spec) {
        setPixel(x, y, spec.r, spec.g, spec.b);
    }

    //--------------------------------------------------------------------------
    // Format Conversion
    //--------------------------------------------------------------------------

    /// Convert to 8-bit sRGB (gamma corrected)
    std::vector<uint8_t> toSRGB8() const {
        std::vector<uint8_t> result(pixelCount() * 3);

        for (size_t i = 0; i < pixelCount(); ++i) {
            for (int c = 0; c < 3; ++c) {
                float linear = pixels[i * 3 + c];
                float clamped = std::clamp(linear, 0.0f, 1.0f);

                // sRGB gamma
                float srgb;
                if (clamped <= 0.0031308f) {
                    srgb = 12.92f * clamped;
                } else {
                    srgb = 1.055f * std::pow(clamped, 1.0f / 2.4f) - 0.055f;
                }

                result[i * 3 + c] = static_cast<uint8_t>(srgb * 255.0f + 0.5f);
            }
        }

        return result;
    }

    /// Convert to RGBA with alpha = 1.0
    std::vector<float> toRGBA() const {
        std::vector<float> result(pixelCount() * 4);

        for (size_t i = 0; i < pixelCount(); ++i) {
            result[i * 4 + 0] = pixels[i * 3 + 0];
            result[i * 4 + 1] = pixels[i * 3 + 1];
            result[i * 4 + 2] = pixels[i * 3 + 2];
            result[i * 4 + 3] = 1.0f;
        }

        return result;
    }
};

//==============================================================================
// ImageBufferRGBA
// Floating-point HDR image buffer (RGBA)
//==============================================================================

class ImageBufferRGBA {
public:
    int width = 0;
    int height = 0;
    std::vector<float> pixels;  // RGBA interleaved (4 floats per pixel)

    ImageBufferRGBA() = default;

    ImageBufferRGBA(int w, int h) {
        allocate(w, h);
    }

    void allocate(int w, int h) {
        width = w;
        height = h;
        pixels.resize(static_cast<size_t>(w) * h * 4, 0.0f);
    }

    void setPixel(int x, int y, float r, float g, float b, float a = 1.0f) {
        if (x < 0 || x >= width || y < 0 || y >= height) return;
        size_t idx = (static_cast<size_t>(y) * width + x) * 4;
        pixels[idx + 0] = r;
        pixels[idx + 1] = g;
        pixels[idx + 2] = b;
        pixels[idx + 3] = a;
    }

    size_t pixelCount() const { return static_cast<size_t>(width) * height; }
};

} // namespace sirius::render

#endif // SIRIUS_RENDER_OUIB001A_H
