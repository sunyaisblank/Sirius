// EXR writer implementation over tinyexr. Ported from OUEW001A.cpp.
//
// This translation unit is the single home of TINYEXR_IMPLEMENTATION for the
// render library; miniz.c (compiled separately) provides the deflate codec.
// Values are stored as float so every finite fp32 radiance value survives
// without half-float overflow.

// miniz must precede tinyexr; it supplies the zlib-compatible codec.
#include "miniz.h"

#define TINYEXR_IMPLEMENTATION
#define TINYEXR_USE_MINIZ 1
#include "sirius/render/exr_writer.h"

#include "tinyexr.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <string>
#include <vector>

namespace sirius::render {

namespace {

// Write an interleaved RGB float buffer to EXR (channels stored B, G, R).
struct AttributeStorage {
    std::vector<std::string> values;
    std::vector<EXRAttribute> attributes;
};

AttributeStorage BuildAttributes(const EXRMetadata& meta) {
    AttributeStorage storage;
    storage.values = {
        meta.software_version,
        meta.camera_model,
        meta.metric_type,
        std::to_string(meta.black_hole_mass),
        std::to_string(meta.black_hole_spin),
        std::to_string(meta.observer_distance),
        std::to_string(meta.observer_inclination),
        std::to_string(meta.samples_per_pixel),
        std::to_string(meta.render_time_seconds),
        meta.color_space,
        meta.chromaticities,
    };
    constexpr std::array<const char*, 11> names = {
        "siriusSoftware",
        "siriusCamera",
        "siriusMetric",
        "siriusBlackHoleMass",
        "siriusBlackHoleSpin",
        "siriusObserverDistance",
        "siriusObserverInclination",
        "siriusSamplesPerPixel",
        "siriusRenderTimeSeconds",
        "siriusColorSpace",
        "siriusChromaticities",
    };
    storage.attributes.resize(names.size());
    for (std::size_t i = 0; i < names.size(); ++i) {
        EXRAttribute& attribute = storage.attributes[i];
        std::memset(&attribute, 0, sizeof(attribute));
        std::strncpy(attribute.name, names[i], sizeof(attribute.name) - 1);
        std::strncpy(attribute.type, "string", sizeof(attribute.type) - 1);
        attribute.value = reinterpret_cast<unsigned char*>(storage.values[i].data());
        attribute.size = static_cast<int>(storage.values[i].size());
    }
    return storage;
}

bool AllFinite(const float* pixels, std::size_t count) {
    return std::all_of(pixels, pixels + count, [](float value) { return std::isfinite(value); });
}

bool WriteRgbFloat(const std::string& path, int width, int height, const float* pixels,
                   const EXRMetadata& meta) {
    if (width <= 0 || height <= 0 || !pixels) {
        return false;
    }
    const std::size_t pixel_count = static_cast<std::size_t>(width) * height;
    if (!AllFinite(pixels, pixel_count * 3)) return false;

    std::vector<float> r(pixel_count);
    std::vector<float> g(pixel_count);
    std::vector<float> b(pixel_count);

    for (std::size_t i = 0; i < pixel_count; ++i) {
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

    // Store fp32 in the file: converting a finite value above 65504 to half
    // would create infinity and violate the output-integrity contract.
    int pixelTypes[3] = {TINYEXR_PIXELTYPE_FLOAT, TINYEXR_PIXELTYPE_FLOAT, TINYEXR_PIXELTYPE_FLOAT};
    int requestedPixelTypes[3] = {TINYEXR_PIXELTYPE_FLOAT, TINYEXR_PIXELTYPE_FLOAT,
                                  TINYEXR_PIXELTYPE_FLOAT};

    header.pixel_types = pixelTypes;
    header.requested_pixel_types = requestedPixelTypes;

    AttributeStorage attributes = BuildAttributes(meta);
    header.num_custom_attributes = static_cast<int>(attributes.attributes.size());
    header.custom_attributes = attributes.attributes.data();
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
bool WriteRgbaFloat(const std::string& path, int width, int height, const float* pixels,
                    const EXRMetadata& meta) {
    if (width <= 0 || height <= 0 || !pixels) {
        return false;
    }
    const std::size_t pixel_count = static_cast<std::size_t>(width) * height;
    if (!AllFinite(pixels, pixel_count * 4)) return false;

    std::vector<float> r(pixel_count);
    std::vector<float> g(pixel_count);
    std::vector<float> b(pixel_count);
    std::vector<float> a(pixel_count);

    for (std::size_t i = 0; i < pixel_count; ++i) {
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
    int requestedPixelTypes[4] = {TINYEXR_PIXELTYPE_FLOAT, TINYEXR_PIXELTYPE_FLOAT,
                                  TINYEXR_PIXELTYPE_FLOAT, TINYEXR_PIXELTYPE_FLOAT};

    header.pixel_types = pixelTypes;
    header.requested_pixel_types = requestedPixelTypes;

    AttributeStorage attributes = BuildAttributes(meta);
    header.num_custom_attributes = static_cast<int>(attributes.attributes.size());
    header.custom_attributes = attributes.attributes.data();
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
                         const EXRMetadata& meta) {
    if (!buffer.HasValidShape()) return false;

    return WriteRgbFloat(path, buffer.width, buffer.height, buffer.pixels.data(), meta);
}

bool EXRWriter::WriteExr(const std::string& path, const ImageBufferRGBA& buffer,
                         const EXRMetadata& meta) {
    if (!buffer.HasValidShape()) return false;

    return WriteRgbaFloat(path, buffer.width, buffer.height, buffer.pixels.data(), meta);
}

}  // namespace sirius::render
