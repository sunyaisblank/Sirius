#pragma once

// PNG output via stb_image_write. Ported from OUPN001A.h.
//
// PNG is a display-target format: this writer owns the sRGB transfer encode and
// applies it exactly once (linear radiance in, sRGB 8-bit out). Alpha is passed
// through linearly (it carries no gamma).

#include "sirius/render/image_buffer.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <string>

namespace sirius::render {

// Static utility for writing PNG files.
class PNGWriter {
  public:
    // Write an RGB HDR buffer (sRGB-encoded on the way out).
    static bool Write(const std::string& path, const ImageBuffer& buffer);

    // Write an RGBA HDR buffer (sRGB-encoded; alpha linear).
    static bool Write(const std::string& path, const ImageBufferRGBA& buffer);

    // Write raw interleaved RGBA float data.
    static bool WriteRgba(const std::string& path, int width, int height, const float* pixels);

    // Write raw interleaved RGB float data.
    static bool WriteRgb(const std::string& path, int width, int height, const float* pixels);

  private:
    // Linear channel to 8-bit sRGB.
    static uint8_t ToSrgb8(float linear) {
        float clamped = std::clamp(linear, 0.0f, 1.0f);
        float srgb;
        if (clamped <= 0.0031308f) {
            srgb = 12.92f * clamped;
        } else {
            srgb = 1.055f * std::pow(clamped, 1.0f / 2.4f) - 0.055f;
        }
        return static_cast<uint8_t>(std::clamp(srgb * 255.0f + 0.5f, 0.0f, 255.0f));
    }
};

}  // namespace sirius::render
