// SRTL001A.h - Tile Layout and Render Configuration
// Component ID: SRTL001A
// Tile decomposition for batch rendering, render configuration, and progress.
// Replaces the deprecated RenderSession (SRRS001A).
//
// Governance: docs/specification.md

#ifndef SIRIUS_RENDER_SRTL001A_H
#define SIRIUS_RENDER_SRTL001A_H

#include "../Camera/CMBL001A.h"
#include <algorithm>
#include <string>
#include <vector>

namespace sirius::render {

//==============================================================================
// RenderConfig
// Configuration for a render session
//==============================================================================

struct RenderConfig {
    // Image dimensions
    int imageWidth = 1920;
    int imageHeight = 1080;
    int tileSize = 256;

    // Quality settings
    int samplesPerPixel = 1;
    int maxBounces = 1;
    double tolerance = 1e-6;

    // Camera
    CameraModel camera;

    // Metric parameters
    std::string metricType = "Schwarzschild";
    double blackHoleMass = 1.0;
    double blackHoleSpin = 0.0;

    // Output settings
    std::string outputPath;
    bool saveCheckpoints = true;
    int checkpointInterval = 10;  // Save every N tiles
};

//==============================================================================
// TileInfo
// Information about a single render tile
//==============================================================================

struct TileInfo {
    int index = 0;
    int x = 0;
    int y = 0;
    int width = 0;
    int height = 0;
    bool completed = false;
};

//==============================================================================
// RenderProgress
//==============================================================================

struct RenderProgress {
    int tilesCompleted = 0;
    int tilesTotal = 0;
    double percentComplete = 0.0;
    double elapsedSeconds = 0.0;
    double estimatedRemaining = 0.0;
};

//==============================================================================
// TileLayout
// Decomposes an image into rectangular tiles for batch rendering.
//==============================================================================

class TileLayout {
public:
    TileLayout(int imageWidth, int imageHeight, int tileSize) {
        int tilesX = (imageWidth + tileSize - 1) / tileSize;
        int tilesY = (imageHeight + tileSize - 1) / tileSize;

        int idx = 0;
        for (int ty = 0; ty < tilesY; ++ty) {
            for (int tx = 0; tx < tilesX; ++tx) {
                TileInfo tile;
                tile.index = idx++;
                tile.x = tx * tileSize;
                tile.y = ty * tileSize;
                tile.width = std::min(tileSize, imageWidth - tile.x);
                tile.height = std::min(tileSize, imageHeight - tile.y);
                tile.completed = false;
                m_tiles.push_back(tile);
            }
        }
    }

    const std::vector<TileInfo>& tiles() const { return m_tiles; }

    void markTileCompleted(int tileIndex) {
        if (tileIndex >= 0 && tileIndex < static_cast<int>(m_tiles.size())) {
            m_tiles[tileIndex].completed = true;
        }
    }

private:
    std::vector<TileInfo> m_tiles;
};

} // namespace sirius::render

#endif // SIRIUS_RENDER_SRTL001A_H
