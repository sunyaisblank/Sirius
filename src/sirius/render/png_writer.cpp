// PNG writer implementation over stb_image_write. Ported from OUPN001A.cpp.
//
#include "sirius/render/png_writer.h"

#include "stb_image_write.h"

#include <algorithm>
#include <cmath>
#include <vector>

namespace sirius::render {

bool PNGWriter::Write(const std::string& path, const ImageBuffer& buffer) {
    if (!buffer.HasValidShape()) return false;
    return WriteRgb(path, buffer.width, buffer.height, buffer.pixels.data());
}

bool PNGWriter::Write(const std::string& path, const ImageBufferRGBA& buffer) {
    if (!buffer.HasValidShape()) return false;
    return WriteRgba(path, buffer.width, buffer.height, buffer.pixels.data());
}

bool PNGWriter::WriteRgba(const std::string& path, int width, int height, const float* pixels) {
    if (width <= 0 || height <= 0 || !pixels) {
        return false;
    }
    const std::size_t sample_count = static_cast<std::size_t>(width) * height * 4;
    if (!std::all_of(pixels, pixels + sample_count,
                     [](float value) { return std::isfinite(value); })) {
        return false;
    }

    std::vector<uint8_t> rgb8(static_cast<size_t>(width) * height * 4);

    for (size_t i = 0; i < static_cast<size_t>(width) * height; ++i) {
        rgb8[i * 4 + 0] = ToSrgb8(pixels[i * 4 + 0]);  // R
        rgb8[i * 4 + 1] = ToSrgb8(pixels[i * 4 + 1]);  // G
        rgb8[i * 4 + 2] = ToSrgb8(pixels[i * 4 + 2]);  // B
        rgb8[i * 4 + 3] = static_cast<uint8_t>(
            std::clamp(pixels[i * 4 + 3] * 255.0f + 0.5f, 0.0f, 255.0f));  // A (linear)
    }

    int result = stbi_write_png(path.c_str(), width, height, 4, rgb8.data(), width * 4);
    return result != 0;
}

bool PNGWriter::WriteRgb(const std::string& path, int width, int height, const float* pixels) {
    if (width <= 0 || height <= 0 || !pixels) {
        return false;
    }
    const std::size_t sample_count = static_cast<std::size_t>(width) * height * 3;
    if (!std::all_of(pixels, pixels + sample_count,
                     [](float value) { return std::isfinite(value); })) {
        return false;
    }

    std::vector<uint8_t> rgb8(static_cast<size_t>(width) * height * 3);

    for (size_t i = 0; i < static_cast<size_t>(width) * height; ++i) {
        rgb8[i * 3 + 0] = ToSrgb8(pixels[i * 3 + 0]);  // R
        rgb8[i * 3 + 1] = ToSrgb8(pixels[i * 3 + 1]);  // G
        rgb8[i * 3 + 2] = ToSrgb8(pixels[i * 3 + 2]);  // B
    }

    int result = stbi_write_png(path.c_str(), width, height, 3, rgb8.data(), width * 3);
    return result != 0;
}

}  // namespace sirius::render
