#pragma once

// OpenEXR output via tinyexr. Ported from OUEW001A.h.
//
// EXR is the HDR interchange path: the writer receives the linear radiance
// buffer untouched and stores it as fp32. It applies no transfer encode,
// tonemap, or grade; those are display concerns owned by the PNG/PPM writers.
// This is the transfer-encode single authority the specification preserves.

#include "sirius/render/image_buffer.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <sstream>
#include <string>

namespace sirius::render {

// Metadata embedded in (or accompanying) an EXR file.
struct EXRMetadata {
    std::string softwareVersion = "Sirius 1.0";
    std::string cameraModel;
    std::string metricType;

    // Physical parameters.
    double blackHoleMass = 1.0;
    double blackHoleSpin = 0.0;
    double observerDistance = 50.0;
    double observerInclination = 90.0;

    // Render parameters.
    int samplesPerPixel = 1;
    double renderTimeSeconds = 0.0;

    // Sirius currently emits its linear working RGB values without an ACES
    // primary conversion. State that exactly instead of attaching false ACES
    // chromaticity metadata.
    std::string colorSpace = "Sirius linear RGB";
    std::string chromaticities = "unspecified";
};

// Static utility for writing EXR and explicit PPM/PFM files.
class EXRWriter {
  public:
    // Human-readable form of the same metadata embedded in EXR attributes.
    static std::string GenerateMetadataHeader(const EXRMetadata& meta) {
        std::ostringstream oss;

        oss << "# " << meta.colorSpace << "\n";
        oss << "# Sirius Ray Tracer\n";
        oss << "# Software: " << meta.softwareVersion << "\n";
        oss << "#\n";
        oss << "# Physical Parameters:\n";
        oss << "#   Black Hole Mass: " << meta.blackHoleMass << " M_sun\n";
        oss << "#   Black Hole Spin: " << meta.blackHoleSpin << "\n";
        oss << "#   Observer Distance: " << meta.observerDistance << " M\n";
        oss << "#   Observer Inclination: " << meta.observerInclination << " deg\n";
        oss << "#\n";
        oss << "# Render Settings:\n";
        oss << "#   Samples Per Pixel: " << meta.samplesPerPixel << "\n";
        oss << "#   Render Time: " << meta.renderTimeSeconds << " s\n";
        oss << "#\n";
        oss << "# Color Space: " << meta.colorSpace << "\n";

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

        auto srgb = buffer.ToSrgb8();
        const std::size_t written = std::fwrite(srgb.data(), 1, srgb.size(), fp);
        const int close_result = std::fclose(fp);
        return written == srgb.size() && close_result == 0;
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
            const size_t idx = static_cast<size_t>(y) * buffer.width * 3;
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
};

}  // namespace sirius::render
