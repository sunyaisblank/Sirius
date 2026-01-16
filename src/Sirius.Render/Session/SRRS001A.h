// SRRS001A.h - Render Session
// Component ID: SRRS001A
// Manages tile-based rendering sessions
//
// Governance: docs/specification.md

#ifndef SIRIUS_RENDER_SRRS001A_H
#define SIRIUS_RENDER_SRRS001A_H

#include "../Camera/CMBL001A.h"
#include <functional>
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
// SessionState
//==============================================================================

enum class SessionState {
    IDLE,       // Not started
    RUNNING,    // Actively rendering
    PAUSED,     // Paused by user
    COMPLETED,  // Finished successfully
    CANCELLED,  // Stopped by user
    ERROR       // Failed with error
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
// RenderSession
// Manages tile-based rendering
//==============================================================================

class RenderSession {
public:
    using ProgressCallback = std::function<void(const RenderProgress&)>;

    explicit RenderSession(const RenderConfig& config)
        : m_config(config)
        , m_state(SessionState::IDLE)
    {
        initializeTiles();
    }

    //--------------------------------------------------------------------------
    // Tile Access
    //--------------------------------------------------------------------------

    const std::vector<TileInfo>& tiles() const { return m_tiles; }

    //--------------------------------------------------------------------------
    // State
    //--------------------------------------------------------------------------

    SessionState state() const { return m_state; }

    RenderProgress getProgress() const {
        RenderProgress prog;
        prog.tilesTotal = static_cast<int>(m_tiles.size());
        for (const auto& t : m_tiles) {
            if (t.completed) ++prog.tilesCompleted;
        }
        if (prog.tilesTotal > 0) {
            prog.percentComplete = 100.0 * prog.tilesCompleted / prog.tilesTotal;
        }
        return prog;
    }

    //--------------------------------------------------------------------------
    // Callbacks
    //--------------------------------------------------------------------------

    void setProgressCallback(ProgressCallback callback) {
        m_progressCallback = std::move(callback);
    }

    //--------------------------------------------------------------------------
    // Control
    //--------------------------------------------------------------------------

    void start() {
        m_state = SessionState::RUNNING;
    }

    void pause() {
        if (m_state == SessionState::RUNNING) {
            m_state = SessionState::PAUSED;
        }
    }

    void resume() {
        if (m_state == SessionState::PAUSED) {
            m_state = SessionState::RUNNING;
        }
    }

    void cancel() {
        m_state = SessionState::CANCELLED;
    }

    void markTileCompleted(int tileIndex) {
        if (tileIndex >= 0 && tileIndex < static_cast<int>(m_tiles.size())) {
            m_tiles[tileIndex].completed = true;
            if (m_progressCallback) {
                m_progressCallback(getProgress());
            }
        }
    }

private:
    void initializeTiles() {
        m_tiles.clear();

        int tilesX = (m_config.imageWidth + m_config.tileSize - 1) / m_config.tileSize;
        int tilesY = (m_config.imageHeight + m_config.tileSize - 1) / m_config.tileSize;

        int idx = 0;
        for (int ty = 0; ty < tilesY; ++ty) {
            for (int tx = 0; tx < tilesX; ++tx) {
                TileInfo tile;
                tile.index = idx++;
                tile.x = tx * m_config.tileSize;
                tile.y = ty * m_config.tileSize;
                tile.width = std::min(m_config.tileSize, m_config.imageWidth - tile.x);
                tile.height = std::min(m_config.tileSize, m_config.imageHeight - tile.y);
                tile.completed = false;
                m_tiles.push_back(tile);
            }
        }
    }

    RenderConfig m_config;
    SessionState m_state;
    std::vector<TileInfo> m_tiles;
    ProgressCallback m_progressCallback;
};

} // namespace sirius::render

#endif // SIRIUS_RENDER_SRRS001A_H
