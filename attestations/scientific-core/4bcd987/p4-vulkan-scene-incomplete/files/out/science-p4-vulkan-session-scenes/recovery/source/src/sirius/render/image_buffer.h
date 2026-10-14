#pragma once

// Floating-point HDR image buffers (RGB and RGBA interleaved).
//
// This buffer holds linear radiance. The transfer encode is applied once by the
// owning writer (PNG applies sRGB, EXR stays linear), never here; ToSrgb8 exists
// for callers that ask for an 8-bit view explicitly.

#include "sirius/core/srgb_transfer.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

namespace sirius::render {

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

    // Convert to 8-bit sRGB (gamma corrected); the display transfer applied once.
    std::vector<std::uint8_t> ToSrgb8() const {
        if (!HasValidShape() || !std::all_of(pixels.begin(), pixels.end(),
                                             [](float value) { return std::isfinite(value); })) {
            return {};
        }
        std::vector<std::uint8_t> result(PixelCount() * 3);

        for (std::size_t i = 0; i < PixelCount(); ++i) {
            for (int c = 0; c < 3; ++c) {
                result[i * 3 + c] = *core::colour::TryEncodeSrgb8(pixels[i * 3 + c]);
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
