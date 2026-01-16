// SRIO001A.h - Session I/O Utilities
// Component ID: SRIO001A
// File I/O for render sessions (checkpoint saving/loading, image export)
//
// Governance: docs/specification.md

#ifndef SIRIUS_RENDER_SRIO001A_H
#define SIRIUS_RENDER_SRIO001A_H

#include "SRCK001A.h"
#include <cstdint>
#include <fstream>
#include <string>

namespace sirius::render {

//==============================================================================
// Checkpoint I/O
//==============================================================================

/// Save checkpoint to file
inline bool saveCheckpoint(const Checkpoint& cp, const std::string& path) {
    std::ofstream file(path, std::ios::binary);
    if (!file) return false;

    file.write(reinterpret_cast<const char*>(&cp.magic), sizeof(cp.magic));
    file.write(reinterpret_cast<const char*>(&cp.version), sizeof(cp.version));
    file.write(reinterpret_cast<const char*>(&cp.tilesX), sizeof(cp.tilesX));
    file.write(reinterpret_cast<const char*>(&cp.tilesY), sizeof(cp.tilesY));
    file.write(reinterpret_cast<const char*>(&cp.tilesCompleted), sizeof(cp.tilesCompleted));
    file.write(reinterpret_cast<const char*>(&cp.imageWidth), sizeof(cp.imageWidth));
    file.write(reinterpret_cast<const char*>(&cp.imageHeight), sizeof(cp.imageHeight));
    file.write(reinterpret_cast<const char*>(&cp.tileSize), sizeof(cp.tileSize));

    uint32_t maskSize = static_cast<uint32_t>(cp.tileCompletionMask.size());
    file.write(reinterpret_cast<const char*>(&maskSize), sizeof(maskSize));
    if (maskSize > 0) {
        file.write(reinterpret_cast<const char*>(cp.tileCompletionMask.data()), maskSize);
    }

    return file.good();
}

/// Load checkpoint from file
inline bool loadCheckpoint(Checkpoint& cp, const std::string& path) {
    std::ifstream file(path, std::ios::binary);
    if (!file) return false;

    file.read(reinterpret_cast<char*>(&cp.magic), sizeof(cp.magic));
    file.read(reinterpret_cast<char*>(&cp.version), sizeof(cp.version));

    // Validate header
    if (cp.magic != CHECKPOINT_MAGIC || cp.version != CHECKPOINT_VERSION) {
        return false;
    }

    file.read(reinterpret_cast<char*>(&cp.tilesX), sizeof(cp.tilesX));
    file.read(reinterpret_cast<char*>(&cp.tilesY), sizeof(cp.tilesY));
    file.read(reinterpret_cast<char*>(&cp.tilesCompleted), sizeof(cp.tilesCompleted));
    file.read(reinterpret_cast<char*>(&cp.imageWidth), sizeof(cp.imageWidth));
    file.read(reinterpret_cast<char*>(&cp.imageHeight), sizeof(cp.imageHeight));
    file.read(reinterpret_cast<char*>(&cp.tileSize), sizeof(cp.tileSize));

    uint32_t maskSize = 0;
    file.read(reinterpret_cast<char*>(&maskSize), sizeof(maskSize));
    if (maskSize > 0 && maskSize < 1000000) {  // Sanity limit
        cp.tileCompletionMask.resize(maskSize);
        file.read(reinterpret_cast<char*>(cp.tileCompletionMask.data()), maskSize);
    }

    return file.good() && cp.isValid();
}

//==============================================================================
// Image I/O Utilities
//==============================================================================

/// Write raw float RGBA image to file
inline bool writeRawImage(const float* pixels, int width, int height,
                          const std::string& path) {
    std::ofstream file(path, std::ios::binary);
    if (!file) return false;

    // Simple header: width, height
    file.write(reinterpret_cast<const char*>(&width), sizeof(width));
    file.write(reinterpret_cast<const char*>(&height), sizeof(height));

    // Pixel data (RGBA float)
    size_t dataSize = static_cast<size_t>(width) * height * 4 * sizeof(float);
    file.write(reinterpret_cast<const char*>(pixels), dataSize);

    return file.good();
}

} // namespace sirius::render

#endif // SIRIUS_RENDER_SRIO001A_H
