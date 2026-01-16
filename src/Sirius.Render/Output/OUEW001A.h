// OUEW001A.h - EXR Writer
// Component ID: OUEW001A
// OpenEXR file output using tinyexr
//
// Governance: docs/specification.md

#ifndef SIRIUS_RENDER_OUEW001A_H
#define SIRIUS_RENDER_OUEW001A_H

#include "OUIB001A.h"
#include "../Session/SRRS001A.h"
#include <string>
#include <vector>
#include <sstream>
#include <ctime>

namespace sirius::render {

//==============================================================================
// EXRMetadata
// Metadata to be embedded in EXR file
//==============================================================================

struct EXRMetadata {
    std::string softwareVersion = "Sirius 1.0";
    std::string cameraModel;
    std::string metricType;

    // Physical parameters
    double blackHoleMass = 1.0;
    double blackHoleSpin = 0.0;
    double observerDistance = 50.0;
    double observerInclination = 90.0;

    // Render parameters
    int samplesPerPixel = 1;
    double renderTimeSeconds = 0.0;

    // ACES color space metadata
    std::string colorSpace = "ACES AP0";
    std::string chromaticities = "ACES";
};

//==============================================================================
// EXRWriter
// Static utility for writing EXR files
//==============================================================================

class EXRWriter {
public:
    //--------------------------------------------------------------------------
    // Buffer Conversion
    //--------------------------------------------------------------------------

    /// Convert RenderSession tiles to ImageBuffer
    static ImageBuffer sessionToBuffer(const RenderSession& session) {
        const auto& tiles = session.tiles();
        if (tiles.empty()) return ImageBuffer();

        // Get dimensions from session config (need to compute from tiles)
        int maxX = 0, maxY = 0;
        for (const auto& tile : tiles) {
            maxX = std::max(maxX, tile.x + tile.width);
            maxY = std::max(maxY, tile.y + tile.height);
        }

        ImageBuffer buffer;
        buffer.allocate(maxX, maxY);

        // In real implementation, would copy tile data
        // For now, return empty buffer of correct size

        return buffer;
    }

    //--------------------------------------------------------------------------
    // Header Generation
    //--------------------------------------------------------------------------

    /// Generate ACES metadata header string
    static std::string generateACESHeader(const EXRMetadata& meta) {
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

    //--------------------------------------------------------------------------
    // File I/O
    //--------------------------------------------------------------------------

    /// Write ImageBuffer to EXR file
    /// Returns true on success, false on failure
    /// Uses tinyexr for actual EXR output
    static bool writeEXR(const std::string& path, const ImageBuffer& buffer,
                         const EXRMetadata& meta = EXRMetadata());

    /// Write ImageBufferRGBA to EXR file
    static bool writeEXR(const std::string& path, const ImageBufferRGBA& buffer,
                         const EXRMetadata& meta = EXRMetadata());

    //--------------------------------------------------------------------------
    // Fallback Formats
    //--------------------------------------------------------------------------

    /// Write as PPM (portable pixmap)
    static bool writePPM(const std::string& path, const ImageBuffer& buffer) {
        FILE* fp = fopen(path.c_str(), "wb");
        if (!fp) return false;

        fprintf(fp, "P6\n%d %d\n255\n", buffer.width, buffer.height);

        auto srgb = buffer.toSRGB8();
        fwrite(srgb.data(), 1, srgb.size(), fp);

        fclose(fp);
        return true;
    }

    /// Write as PFM (portable float map - HDR)
    static bool writePFM(const std::string& path, const ImageBuffer& buffer) {
        FILE* fp = fopen(path.c_str(), "wb");
        if (!fp) return false;

        // PFM header (little-endian)
        fprintf(fp, "PF\n%d %d\n-1.0\n", buffer.width, buffer.height);

        // PFM stores bottom-to-top, so flip Y
        for (int y = buffer.height - 1; y >= 0; --y) {
            for (int x = 0; x < buffer.width; ++x) {
                size_t idx = (static_cast<size_t>(y) * buffer.width + x) * 3;
                fwrite(&buffer.pixels[idx], sizeof(float), 3, fp);
            }
        }

        fclose(fp);
        return true;
    }
};

} // namespace sirius::render

#endif // SIRIUS_RENDER_OUEW001A_H
