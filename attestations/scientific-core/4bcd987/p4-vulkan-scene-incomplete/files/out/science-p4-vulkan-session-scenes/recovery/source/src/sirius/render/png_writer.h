#pragma once

// PNG output via stb_image_write. Ported from OUPN001A.h.
//
// PNG is a display-target format: this writer owns the sRGB transfer encode and
// applies it exactly once (linear radiance in, sRGB 8-bit out). Alpha is passed
// through linearly (it carries no gamma).

#include "sirius/core/srgb_transfer.h"
#include "sirius/render/image_buffer.h"

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
    static std::uint8_t ToSrgb8(float linear) { return *core::colour::TryEncodeSrgb8(linear); }
};

}  // namespace sirius::render
