// OUEW001A.cpp - EXR Writer Implementation
// Component ID: OUEW001A
// Actual tinyexr integration
//
// Governance: docs/specification.md

// Include miniz implementation before tinyexr
// miniz.h includes both header and implementation when not in header-only mode
#include "miniz.h"

#define TINYEXR_IMPLEMENTATION
#define TINYEXR_USE_MINIZ 1
#include "tinyexr.h"

#include "OUEW001A.h"
#include <cstring>

namespace sirius::render {

//==============================================================================
// EXRWriterImpl
// Actual tinyexr-based implementation
//==============================================================================

class EXRWriterImpl {
public:
    /// Write RGB float buffer to EXR file
    static bool writeRGBFloat(const std::string& path,
                              int width, int height,
                              const float* pixels,
                              const EXRMetadata& meta) {
        if (width <= 0 || height <= 0 || !pixels) {
            return false;
        }

        // Separate RGB channels (tinyexr expects separate channel arrays)
        std::vector<float> r(static_cast<size_t>(width) * height);
        std::vector<float> g(static_cast<size_t>(width) * height);
        std::vector<float> b(static_cast<size_t>(width) * height);

        for (size_t i = 0; i < static_cast<size_t>(width) * height; ++i) {
            r[i] = pixels[i * 3 + 0];
            g[i] = pixels[i * 3 + 1];
            b[i] = pixels[i * 3 + 2];
        }

        // Set up tinyexr structures
        EXRHeader header;
        InitEXRHeader(&header);

        EXRImage image;
        InitEXRImage(&image);

        image.num_channels = 3;
        image.width = width;
        image.height = height;

        // Channel pointers (B, G, R order for EXR convention)
        float* images[3] = {b.data(), g.data(), r.data()};
        image.images = reinterpret_cast<unsigned char**>(images);

        // Channel info
        EXRChannelInfo channelInfo[3];
        std::memset(channelInfo, 0, sizeof(channelInfo));

        std::strncpy(channelInfo[0].name, "B", 255);
        std::strncpy(channelInfo[1].name, "G", 255);
        std::strncpy(channelInfo[2].name, "R", 255);

        header.channels = channelInfo;
        header.num_channels = 3;

        // Pixel types
        int pixelTypes[3] = {TINYEXR_PIXELTYPE_FLOAT, TINYEXR_PIXELTYPE_FLOAT, TINYEXR_PIXELTYPE_FLOAT};
        int requestedPixelTypes[3] = {TINYEXR_PIXELTYPE_HALF, TINYEXR_PIXELTYPE_HALF, TINYEXR_PIXELTYPE_HALF};

        header.pixel_types = pixelTypes;
        header.requested_pixel_types = requestedPixelTypes;

        // Compression
        header.compression_type = TINYEXR_COMPRESSIONTYPE_ZIP;

        // Write file
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

    /// Write RGBA float buffer to EXR file
    static bool writeRGBAFloat(const std::string& path,
                               int width, int height,
                               const float* pixels,
                               const EXRMetadata& meta) {
        if (width <= 0 || height <= 0 || !pixels) {
            return false;
        }

        // Separate RGBA channels
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

        // Set up tinyexr structures
        EXRHeader header;
        InitEXRHeader(&header);

        EXRImage image;
        InitEXRImage(&image);

        image.num_channels = 4;
        image.width = width;
        image.height = height;

        // Channel pointers (A, B, G, R order for EXR convention)
        float* images[4] = {a.data(), b.data(), g.data(), r.data()};
        image.images = reinterpret_cast<unsigned char**>(images);

        // Channel info
        EXRChannelInfo channelInfo[4];
        std::memset(channelInfo, 0, sizeof(channelInfo));

        std::strncpy(channelInfo[0].name, "A", 255);
        std::strncpy(channelInfo[1].name, "B", 255);
        std::strncpy(channelInfo[2].name, "G", 255);
        std::strncpy(channelInfo[3].name, "R", 255);

        header.channels = channelInfo;
        header.num_channels = 4;

        // Pixel types - store as float, request half for file
        int pixelTypes[4] = {TINYEXR_PIXELTYPE_FLOAT, TINYEXR_PIXELTYPE_FLOAT,
                             TINYEXR_PIXELTYPE_FLOAT, TINYEXR_PIXELTYPE_FLOAT};
        int requestedPixelTypes[4] = {TINYEXR_PIXELTYPE_HALF, TINYEXR_PIXELTYPE_HALF,
                                      TINYEXR_PIXELTYPE_HALF, TINYEXR_PIXELTYPE_HALF};

        header.pixel_types = pixelTypes;
        header.requested_pixel_types = requestedPixelTypes;

        // Compression
        header.compression_type = TINYEXR_COMPRESSIONTYPE_ZIP;

        // Write file
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
};

//==============================================================================
// EXRWriter static method implementations
//==============================================================================

bool EXRWriter::writeEXR(const std::string& path, const ImageBuffer& buffer,
                         const EXRMetadata& meta) {
    if (buffer.width <= 0 || buffer.height <= 0) return false;
    if (buffer.pixels.empty()) return false;

    return EXRWriterImpl::writeRGBFloat(path, buffer.width, buffer.height,
                                        buffer.pixels.data(), meta);
}

bool EXRWriter::writeEXR(const std::string& path, const ImageBufferRGBA& buffer,
                         const EXRMetadata& meta) {
    if (buffer.width <= 0 || buffer.height <= 0) return false;
    if (buffer.pixels.empty()) return false;

    return EXRWriterImpl::writeRGBAFloat(path, buffer.width, buffer.height,
                                         buffer.pixels.data(), meta);
}

} // namespace sirius::render
