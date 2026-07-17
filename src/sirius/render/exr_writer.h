#pragma once

// OpenEXR output via tinyexr. Ported from OUEW001A.h.
//
// EXR is the HDR interchange path: the writer receives the linear radiance
// buffer untouched and stores it (as half float). It applies no transfer encode,
// tonemap, or grade; those are display concerns owned by the PNG/PPM writers.
// This is the transfer-encode single authority the specification preserves.

#include "sirius/render/image_buffer.h"

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

    // ACES colour space metadata.
    std::string colorSpace = "ACES AP0";
    std::string chromaticities = "ACES";
};

// Static utility for writing EXR (and PPM/PFM fallback) files.
class EXRWriter {
  public:
    // ACES metadata header string (informational; not part of the EXR binary).
    static std::string GenerateAcesHeader(const EXRMetadata& meta) {
        std::ostringstream oss;

        oss << "# ACES AP0 Linear\n";
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

    // PPM fallback (portable pixmap, 8-bit sRGB-encoded).
    static bool WritePpm(const std::string& path, const ImageBuffer& buffer) {
        FILE* fp = std::fopen(path.c_str(), "wb");
        if (!fp) return false;

        std::fprintf(fp, "P6\n%d %d\n255\n", buffer.width, buffer.height);

        auto srgb = buffer.ToSrgb8();
        std::fwrite(srgb.data(), 1, srgb.size(), fp);

        std::fclose(fp);
        return true;
    }

    // PFM fallback (portable float map, HDR, little-endian, bottom-to-top).
    static bool WritePfm(const std::string& path, const ImageBuffer& buffer) {
        FILE* fp = std::fopen(path.c_str(), "wb");
        if (!fp) return false;

        std::fprintf(fp, "PF\n%d %d\n-1.0\n", buffer.width, buffer.height);

        for (int y = buffer.height - 1; y >= 0; --y) {
            for (int x = 0; x < buffer.width; ++x) {
                size_t idx = (static_cast<size_t>(y) * buffer.width + x) * 3;
                std::fwrite(&buffer.pixels[idx], sizeof(float), 3, fp);
            }
        }

        std::fclose(fp);
        return true;
    }
};

}  // namespace sirius::render
