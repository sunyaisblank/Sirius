#pragma once

// OpenEXR output via tinyexr. Ported from OUEW001A.h.
//
// EXR is the HDR interchange path: the writer receives the linear radiance
// buffer untouched and stores it as fp32. It applies no transfer encode,
// tonemap, or grade; those are display concerns owned by the PNG/PPM writers.
// This is the transfer-encode single authority the specification preserves.

#include "sirius/core/srgb_transfer.h"
#include "sirius/render/image_buffer.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <sstream>
#include <string>

namespace sirius::render {

// Metadata embedded in (or accompanying) an EXR file.
struct EXRMetadata {
    std::string software_version = "Sirius 1.0";
    std::string camera_model;
    std::string metric_type;

    // Geometric scene parameters. Sirius does not infer a solar-mass scale.
    double black_hole_mass = 1.0;
    double black_hole_spin = 0.0;
    double observer_distance = 50.0;
    double observer_inclination = 90.0;

    // Render parameters.
    int samples_per_pixel = 1;
    double render_time_seconds = 0.0;

    // Sirius currently emits its linear working RGB values without an ACES
    // primary conversion. State that exactly instead of attaching false ACES
    // chromaticity metadata.
    std::string color_space = "Sirius linear RGB";
    std::string chromaticities = "unspecified";
};

// Static utility for writing EXR and explicit PPM/PFM files.
class EXRWriter {
  public:
    // Human-readable form of the same metadata embedded in EXR attributes.
    static std::string GenerateMetadataHeader(const EXRMetadata& meta) {
        std::ostringstream oss;

        oss << "# " << meta.color_space << "\n";
        oss << "# Sirius Ray Tracer\n";
        oss << "# Software: " << meta.software_version << "\n";
        oss << "#\n";
        oss << "# Geometric Scene Parameters:\n";
        oss << "#   Metric Mass Parameter: " << meta.black_hole_mass << " coordinate units\n";
        oss << "#   Dimensionless Spin a/M: " << meta.black_hole_spin << "\n";
        oss << "#   Observer Coordinate Radius: " << meta.observer_distance
            << " coordinate units\n";
        oss << "#   Observer Inclination: " << meta.observer_inclination << " deg\n";
        oss << "#\n";
        oss << "# Render Settings:\n";
        oss << "#   Samples Per Pixel: " << meta.samples_per_pixel << "\n";
        oss << "#   Render Time: " << meta.render_time_seconds << " s\n";
        oss << "#\n";
        oss << "# Color Space: " << meta.color_space << "\n";

        return oss.str();
    }

    // Write an RGB HDR buffer to EXR. Returns true on success.
    static bool WriteExr(const std::string& path, const ImageBuffer& buffer,
                         const EXRMetadata& meta = EXRMetadata());

    // Write an RGBA HDR buffer to EXR. Returns true on success.
    static bool WriteExr(const std::string& path, const ImageBufferRGBA& buffer,
                         const EXRMetadata& meta = EXRMetadata());

    // Explicit PPM writer (portable pixmap, 8-bit sRGB-encoded).
    static bool WritePpm(const std::string& path, const ImageBuffer& buffer) {
        if (!buffer.HasValidShape() ||
            !std::all_of(buffer.pixels.begin(), buffer.pixels.end(),
                         [](float value) { return std::isfinite(value); })) {
            return false;
        }
        FILE* fp = std::fopen(path.c_str(), "wb");
        if (!fp) return false;

        if (std::fprintf(fp, "P6\n%d %d\n255\n", buffer.width, buffer.height) < 0) {
            std::fclose(fp);
            return false;
        }

        std::vector<std::uint8_t> srgb(buffer.pixels.size());
        std::transform(buffer.pixels.begin(), buffer.pixels.end(), srgb.begin(), EncodeSrgb8);
        const std::size_t written = std::fwrite(srgb.data(), 1, srgb.size(), fp);
        const int close_result = std::fclose(fp);
        return written == srgb.size() && close_result == 0;
    }

    // Write display-linear RGBA pixels as an RGB PPM. This is the render
    // session's PPM boundary: alpha is discarded and this writer applies the
    // one and only sRGB transfer encode.
    static bool WritePpmRgba(const std::string& path, int width, int height, const float* pixels) {
        if (width <= 0 || height <= 0 || pixels == nullptr) return false;
        const std::size_t pixel_count = static_cast<std::size_t>(width) * height;
        if (!std::all_of(pixels, pixels + pixel_count * 4,
                         [](float value) { return std::isfinite(value); })) {
            return false;
        }

        FILE* fp = std::fopen(path.c_str(), "wb");
        if (!fp) return false;
        if (std::fprintf(fp, "P6\n%d %d\n255\n", width, height) < 0) {
            std::fclose(fp);
            return false;
        }

        bool write_ok = true;
        for (std::size_t i = 0; i < pixel_count; ++i) {
            const std::uint8_t rgb[3] = {EncodeSrgb8(pixels[i * 4]), EncodeSrgb8(pixels[i * 4 + 1]),
                                         EncodeSrgb8(pixels[i * 4 + 2])};
            if (std::fwrite(rgb, 1, sizeof(rgb), fp) != sizeof(rgb)) {
                write_ok = false;
                break;
            }
        }
        const int close_result = std::fclose(fp);
        return write_ok && close_result == 0;
    }

    // Explicit PFM writer (portable float map, HDR, little-endian, bottom-to-top).
    static bool WritePfm(const std::string& path, const ImageBuffer& buffer) {
        if (!buffer.HasValidShape() ||
            !std::all_of(buffer.pixels.begin(), buffer.pixels.end(),
                         [](float value) { return std::isfinite(value); })) {
            return false;
        }
        FILE* fp = std::fopen(path.c_str(), "wb");
        if (!fp) return false;

        if (std::fprintf(fp, "PF\n%d %d\n-1.0\n", buffer.width, buffer.height) < 0) {
            std::fclose(fp);
            return false;
        }

        bool write_ok = true;
        for (int y = buffer.height - 1; y >= 0; --y) {
            const std::size_t idx = static_cast<std::size_t>(y) * buffer.width * 3;
            const std::size_t written = std::fwrite(&buffer.pixels[idx], sizeof(float),
                                                    static_cast<std::size_t>(buffer.width) * 3, fp);
            if (written != static_cast<std::size_t>(buffer.width) * 3) {
                write_ok = false;
                break;
            }
        }

        const int close_result = std::fclose(fp);
        return write_ok && close_result == 0;
    }

  private:
    static std::uint8_t EncodeSrgb8(float linear) { return *core::colour::TryEncodeSrgb8(linear); }
};

}  // namespace sirius::render
