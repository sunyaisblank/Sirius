// OUPN001A.h - PNG Writer
// Component ID: OUPN001A
// PNG file output using stb_image_write

#ifndef SIRIUS_RENDER_OUPN001A_H
#define SIRIUS_RENDER_OUPN001A_H

#include "OUIB001A.h"
#include <string>
#include <cmath>
#include <algorithm>

namespace sirius::render {

//==============================================================================
// PNGWriter
// Static utility for writing PNG files
//==============================================================================

class PNGWriter {
public:
    /// Write ImageBuffer (RGB) to PNG file
    static bool write(const std::string& path, const ImageBuffer& buffer);

    /// Write ImageBufferRGBA to PNG file
    static bool write(const std::string& path, const ImageBufferRGBA& buffer);

    /// Write raw RGBA float data to PNG file
    static bool writeRGBA(const std::string& path, int width, int height,
                          const float* pixels);

    /// Write raw RGB float data to PNG file
    static bool writeRGB(const std::string& path, int width, int height,
                         const float* pixels);

private:
    //--------------------------------------------------------------------------
    // sRGB Gamma Correction
    //--------------------------------------------------------------------------

    static uint8_t toSRGB8(float linear) {
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

} // namespace sirius::render

#endif // SIRIUS_RENDER_OUPN001A_H
