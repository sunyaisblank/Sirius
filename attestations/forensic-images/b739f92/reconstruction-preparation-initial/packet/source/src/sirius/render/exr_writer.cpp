// EXR writer implementation over tinyexr. Ported from OUEW001A.cpp.
//
// This translation unit is the single home of TINYEXR_IMPLEMENTATION for the
// render library; miniz.c (compiled separately) provides the deflate codec.
// Values are stored as half float, so HDR radiance above 1.0 survives.

// miniz must precede tinyexr; it supplies the zlib-compatible codec.
#include "miniz.h"

#define TINYEXR_IMPLEMENTATION
#define TINYEXR_USE_MINIZ 1
#include "tinyexr.h"

#include "sirius/render/exr_writer.h"

#include <cstring>
#include <vector>

namespace sirius::render {

namespace {

// Write an interleaved RGB float buffer to EXR (channels stored B, G, R).
bool WriteRgbFloat(const std::string& path, int width, int height, const float* pixels) {
    if (width <= 0 || height <= 0 || !pixels) {
        return false;
    }

    std::vector<float> r(static_cast<size_t>(width) * height);
    std::vector<float> g(static_cast<size_t>(width) * height);
    std::vector<float> b(static_cast<size_t>(width) * height);

    for (size_t i = 0; i < static_cast<size_t>(width) * height; ++i) {
        r[i] = pixels[i * 3 + 0];
        g[i] = pixels[i * 3 + 1];
        b[i] = pixels[i * 3 + 2];
    }

    EXRHeader header;
    InitEXRHeader(&header);

    EXRImage image;
    InitEXRImage(&image);

    image.num_channels = 3;
    image.width = width;
    image.height = height;

    // EXR convention: channels in B, G, R order.
    float* images[3] = {b.data(), g.data(), r.data()};
    image.images = reinterpret_cast<unsigned char**>(images);

    EXRChannelInfo channelInfo[3];
    std::memset(channelInfo, 0, sizeof(channelInfo));

    std::strncpy(channelInfo[0].name, "B", 255);
    std::strncpy(channelInfo[1].name, "G", 255);
    std::strncpy(channelInfo[2].name, "R", 255);

    header.channels = channelInfo;
    header.num_channels = 3;

    // Stored as float in memory, requested as half for the file.
    int pixelTypes[3] = {TINYEXR_PIXELTYPE_FLOAT, TINYEXR_PIXELTYPE_FLOAT, TINYEXR_PIXELTYPE_FLOAT};
    int requestedPixelTypes[3] = {TINYEXR_PIXELTYPE_HALF, TINYEXR_PIXELTYPE_HALF,
                                  TINYEXR_PIXELTYPE_HALF};

    header.pixel_types = pixelTypes;
    header.requested_pixel_types = requestedPixelTypes;

    header.compression_type = TINYEXR_COMPRESSIONTYPE_ZIP;

    const char* err = nullptr;
    int ret = SaveEXRImageToFile(&image, &header, path.c_str(), &err);

    if (ret != TINYEXR_SUCCESS) {
        if (err) {
            FreeEXRErrorMessage(err);
        }
        return false;
    }

    return true;
}

// Write an interleaved RGBA float buffer to EXR (channels stored A, B, G, R).
bool WriteRgbaFloat(const std::string& path, int width, int height, const float* pixels) {
    if (width <= 0 || height <= 0 || !pixels) {
        return false;
    }

    std::vector<float> r(static_cast<size_t>(width) * height);
    std::vector<float> g(static_cast<size_t>(width) * height);
    std::vector<float> b(static_cast<size_t>(width) * height);
    std::vector<float> a(static_cast<size_t>(width) * height);

    for (size_t i = 0; i < static_cast<size_t>(width) * height; ++i) {
        r[i] = pixels[i * 4 + 0];
        g[i] = pixels[i * 4 + 1];
        b[i] = pixels[i * 4 + 2];
        a[i] = pixels[i * 4 + 3];
    }

    EXRHeader header;
    InitEXRHeader(&header);

    EXRImage image;
    InitEXRImage(&image);

    image.num_channels = 4;
    image.width = width;
    image.height = height;

    // EXR convention: channels in A, B, G, R order.
    float* images[4] = {a.data(), b.data(), g.data(), r.data()};
    image.images = reinterpret_cast<unsigned char**>(images);

    EXRChannelInfo channelInfo[4];
    std::memset(channelInfo, 0, sizeof(channelInfo));

    std::strncpy(channelInfo[0].name, "A", 255);
    std::strncpy(channelInfo[1].name, "B", 255);
    std::strncpy(channelInfo[2].name, "G", 255);
    std::strncpy(channelInfo[3].name, "R", 255);

    header.channels = channelInfo;
    header.num_channels = 4;

    int pixelTypes[4] = {TINYEXR_PIXELTYPE_FLOAT, TINYEXR_PIXELTYPE_FLOAT, TINYEXR_PIXELTYPE_FLOAT,
                         TINYEXR_PIXELTYPE_FLOAT};
    int requestedPixelTypes[4] = {TINYEXR_PIXELTYPE_HALF, TINYEXR_PIXELTYPE_HALF,
                                  TINYEXR_PIXELTYPE_HALF, TINYEXR_PIXELTYPE_HALF};

    header.pixel_types = pixelTypes;
    header.requested_pixel_types = requestedPixelTypes;

    header.compression_type = TINYEXR_COMPRESSIONTYPE_ZIP;

    const char* err = nullptr;
    int ret = SaveEXRImageToFile(&image, &header, path.c_str(), &err);

    if (ret != TINYEXR_SUCCESS) {
        if (err) {
            FreeEXRErrorMessage(err);
        }
        return false;
    }

    return true;
}

}  // namespace

bool EXRWriter::WriteExr(const std::string& path, const ImageBuffer& buffer,
                         [[maybe_unused]] const EXRMetadata& meta) {
    if (buffer.width <= 0 || buffer.height <= 0) return false;
    if (buffer.pixels.empty()) return false;

    return WriteRgbFloat(path, buffer.width, buffer.height, buffer.pixels.data());
}

bool EXRWriter::WriteExr(const std::string& path, const ImageBufferRGBA& buffer,
                         [[maybe_unused]] const EXRMetadata& meta) {
    if (buffer.width <= 0 || buffer.height <= 0) return false;
    if (buffer.pixels.empty()) return false;

    return WriteRgbaFloat(path, buffer.width, buffer.height, buffer.pixels.data());
}

}  // namespace sirius::render
