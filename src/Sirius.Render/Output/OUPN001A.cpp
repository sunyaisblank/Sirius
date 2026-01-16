// OUPN001A.cpp - PNG Output Driver Implementation
// Component ID: OUPN001A
// Actual stb_image_write integration
//
// Governance: docs/specification.md

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "../../../lib/stb/stb_image_write.h"

#include "OUPN001A.h"
#include <vector>

namespace sirius::render {

//==============================================================================
// PNGOutputDriver Implementation
//==============================================================================

bool PNGOutputDriver::writePNG() const {
    if (m_width <= 0 || m_height <= 0 || m_pixels.empty()) {
        return false;
    }

    // Convert float RGBA to sRGB uint8
    std::vector<uint8_t> rgb8(static_cast<size_t>(m_width) * m_height * 4);

    for (size_t i = 0; i < static_cast<size_t>(m_width) * m_height; ++i) {
        rgb8[i * 4 + 0] = toSRGB8(m_pixels[i * 4 + 0]);  // R
        rgb8[i * 4 + 1] = toSRGB8(m_pixels[i * 4 + 1]);  // G
        rgb8[i * 4 + 2] = toSRGB8(m_pixels[i * 4 + 2]);  // B
        rgb8[i * 4 + 3] = static_cast<uint8_t>(
            std::clamp(m_pixels[i * 4 + 3] * 255.0f + 0.5f, 0.0f, 255.0f));  // A
    }

    // Write PNG with stb
    int result = stbi_write_png(m_outputPath.c_str(), m_width, m_height,
                                4, rgb8.data(), m_width * 4);
    return result != 0;
}

//==============================================================================
// PNGWriter Implementation
//==============================================================================

bool PNGWriter::write(const std::string& path, const ImageBuffer& buffer) {
    return writeRGB(path, buffer.width, buffer.height, buffer.pixels.data());
}

bool PNGWriter::write(const std::string& path, const ImageBufferRGBA& buffer) {
    return writeRGBA(path, buffer.width, buffer.height, buffer.pixels.data());
}

bool PNGWriter::writeRGBA(const std::string& path, int width, int height,
                          const float* pixels) {
    if (width <= 0 || height <= 0 || !pixels) {
        return false;
    }

    // Convert float RGBA to sRGB uint8
    std::vector<uint8_t> rgb8(static_cast<size_t>(width) * height * 4);

    for (size_t i = 0; i < static_cast<size_t>(width) * height; ++i) {
        rgb8[i * 4 + 0] = toSRGB8(pixels[i * 4 + 0]);  // R
        rgb8[i * 4 + 1] = toSRGB8(pixels[i * 4 + 1]);  // G
        rgb8[i * 4 + 2] = toSRGB8(pixels[i * 4 + 2]);  // B
        rgb8[i * 4 + 3] = static_cast<uint8_t>(
            std::clamp(pixels[i * 4 + 3] * 255.0f + 0.5f, 0.0f, 255.0f));  // A
    }

    // Write PNG with stb
    int result = stbi_write_png(path.c_str(), width, height,
                                4, rgb8.data(), width * 4);
    return result != 0;
}

bool PNGWriter::writeRGB(const std::string& path, int width, int height,
                         const float* pixels) {
    if (width <= 0 || height <= 0 || !pixels) {
        return false;
    }

    // Convert float RGB to sRGB uint8
    std::vector<uint8_t> rgb8(static_cast<size_t>(width) * height * 3);

    for (size_t i = 0; i < static_cast<size_t>(width) * height; ++i) {
        rgb8[i * 3 + 0] = toSRGB8(pixels[i * 3 + 0]);  // R
        rgb8[i * 3 + 1] = toSRGB8(pixels[i * 3 + 1]);  // G
        rgb8[i * 3 + 2] = toSRGB8(pixels[i * 3 + 2]);  // B
    }

    // Write PNG with stb
    int result = stbi_write_png(path.c_str(), width, height,
                                3, rgb8.data(), width * 3);
    return result != 0;
}

} // namespace sirius::render
