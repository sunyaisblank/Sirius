// OUMD001A.h - Output Driver Interface
// Component ID: OUMD001A
// Abstract interface for render output destinations
//
// Governance: docs/specification.md

#ifndef SIRIUS_RENDER_OUMD001A_H
#define SIRIUS_RENDER_OUMD001A_H

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace sirius::render {

//==============================================================================
// OutputPass Enumeration
//==============================================================================

enum class OutputPass {
    COMBINED,    // Final composited RGBA image
    DEPTH,       // Z-depth (coordinate distance)
    REDSHIFT,    // Gravitational redshift factor
    TEMPERATURE, // Effective temperature (K)
    ALPHA        // Alpha/mask channel
};

/// Get human-readable name for output pass
inline const char* passName(OutputPass pass) {
    switch (pass) {
        case OutputPass::COMBINED:    return "Combined";
        case OutputPass::DEPTH:       return "Depth";
        case OutputPass::REDSHIFT:    return "Redshift";
        case OutputPass::TEMPERATURE: return "Temperature";
        case OutputPass::ALPHA:       return "Alpha";
        default:                      return "Unknown";
    }
}

//==============================================================================
// OutputTile
// A rectangular region of rendered pixels
//==============================================================================

struct OutputTile {
    int x = 0;          // Tile X position in image
    int y = 0;          // Tile Y position in image
    int width = 0;      // Tile width in pixels
    int height = 0;     // Tile height in pixels
    OutputPass pass = OutputPass::COMBINED;
    std::vector<float> pixels;  // RGBA float data (4 floats per pixel)

    size_t pixelCount() const { return static_cast<size_t>(width) * height; }
    size_t bufferSize() const { return pixelCount() * 4; }
};

//==============================================================================
// OutputMetadata
// Information about the completed render
//==============================================================================

struct OutputMetadata {
    int width = 0;
    int height = 0;
    int tilesRendered = 0;
    int samplesPerPixel = 1;        // Samples per pixel (for adaptive sampling)
    double renderTimeSeconds = 0.0;
    std::string metricName;
    double blackHoleSpin = 0.0;
    double massParameter = 1.0;
    std::string comment;            // User comment/notes
};

//==============================================================================
// IOutputDriver - Abstract Interface
//==============================================================================

class IOutputDriver {
public:
    virtual ~IOutputDriver() = default;

    /// Configure driver for a new render session
    virtual bool configure(int width, int height,
                           const std::vector<OutputPass>& passes) = 0;

    /// Begin receiving tiles
    virtual bool start() = 0;

    /// Write a completed tile
    virtual bool writeTile(const OutputTile& tile) = 0;

    /// Finalize output and write metadata
    virtual bool finish(const OutputMetadata& metadata) = 0;

    /// Cancel and clean up partial output
    virtual void cancel() = 0;

    /// Check if driver is currently active
    virtual bool isActive() const = 0;
};

} // namespace sirius::render

#endif // SIRIUS_RENDER_OUMD001A_H
