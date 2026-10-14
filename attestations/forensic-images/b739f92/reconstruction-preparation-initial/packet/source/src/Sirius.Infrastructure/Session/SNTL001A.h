// SNTL001A.h - Tile Scheduler
// Component ID: SNTL001A (Session/Tile)
//
// Manages tile decomposition and scheduling order.
// Uses spiral ordering from centre outward (Cycles pattern).

#pragma once

#include <vector>
#include <atomic>
#include <mutex>
#include <cstdint>
#include <algorithm>
#include <cmath>

namespace Sirius {

//==============================================================================
// Tile State
//==============================================================================
enum class TileState : uint8_t {
    Pending,    ///< Not yet started
    Active,     ///< Currently rendering
    Complete,   ///< Finished successfully
    Failed      ///< Error occurred
};

//==============================================================================
// Tile Definition
//==============================================================================
struct Tile {
    int id;             ///< Unique tile ID
    int x, y;           ///< Top-left corner (pixels)
    int width, height;  ///< Tile dimensions
    TileState state = TileState::Pending;
    int priority = 0;   ///< Lower = higher priority
    
    // Calculated properties
    int centreX() const { return x + width / 2; }
    int centreY() const { return y + height / 2; }
    int pixelCount() const { return width * height; }
};

//==============================================================================
// Tile Scheduler
//==============================================================================
class TileScheduler {
public:
    /// @brief Create tiles for image with spiral ordering
    /// @param imageWidth Full image width
    /// @param imageHeight Full image height
    /// @param tileSize Target tile size (tiles at edges may be smaller)
    void initialise(int imageWidth, int imageHeight, int tileSize = 64) {
        std::lock_guard<std::mutex> lock(m_Mutex);
        
        m_ImageWidth = imageWidth;
        m_ImageHeight = imageHeight;
        m_TileSize = tileSize;
        m_Tiles.clear();
        m_NextTileIndex = 0;
        
        // Generate tiles
        int tilesX = (imageWidth + tileSize - 1) / tileSize;
        int tilesY = (imageHeight + tileSize - 1) / tileSize;
        int tileId = 0;
        
        for (int ty = 0; ty < tilesY; ++ty) {
            for (int tx = 0; tx < tilesX; ++tx) {
                Tile tile;
                tile.id = tileId++;
                tile.x = tx * tileSize;
                tile.y = ty * tileSize;
                tile.width = std::min(tileSize, imageWidth - tile.x);
                tile.height = std::min(tileSize, imageHeight - tile.y);
                tile.state = TileState::Pending;
                m_Tiles.push_back(tile);
            }
        }
        
        // Sort by spiral order from centre
        sortSpiralOrder();
    }
    
    /// @brief Get next pending tile (thread-safe)
    /// @return Pointer to tile, or nullptr if none available
    Tile* getNextTile() {
        std::lock_guard<std::mutex> lock(m_Mutex);
        
        for (size_t i = m_NextTileIndex; i < m_Tiles.size(); ++i) {
            if (m_Tiles[i].state == TileState::Pending) {
                m_Tiles[i].state = TileState::Active;
                m_NextTileIndex = i + 1;
                return &m_Tiles[i];
            }
        }
        return nullptr;
    }
    
    /// @brief Mark tile as complete
    void completeTile(int tileId) {
        std::lock_guard<std::mutex> lock(m_Mutex);
        for (auto& tile : m_Tiles) {
            if (tile.id == tileId) {
                tile.state = TileState::Complete;
                m_CompletedCount++;
                return;
            }
        }
    }
    
    /// @brief Mark tile as failed
    void failTile(int tileId) {
        std::lock_guard<std::mutex> lock(m_Mutex);
        for (auto& tile : m_Tiles) {
            if (tile.id == tileId) {
                tile.state = TileState::Failed;
                return;
            }
        }
    }
    
    /// @brief Reset all tiles to pending
    void reset() {
        std::lock_guard<std::mutex> lock(m_Mutex);
        for (auto& tile : m_Tiles) {
            tile.state = TileState::Pending;
        }
        m_NextTileIndex = 0;
        m_CompletedCount = 0;
    }
    
    // Queries
    int getTileCount() const { return static_cast<int>(m_Tiles.size()); }
    int getCompletedCount() const { return m_CompletedCount.load(); }
    int getPendingCount() const { return getTileCount() - getCompletedCount(); }
    float getProgress() const {
        int total = getTileCount();
        return total > 0 ? static_cast<float>(getCompletedCount()) / total : 0.0f;
    }
    bool allComplete() const { return getCompletedCount() == getTileCount(); }
    
    // Tile access (for display)
    const std::vector<Tile>& getTiles() const { return m_Tiles; }
    
private:
    /// @brief Sort tiles in spiral order from centre
    void sortSpiralOrder() {
        float centreX = m_ImageWidth / 2.0f;
        float centreY = m_ImageHeight / 2.0f;
        
        // Calculate distance from centre for each tile
        for (auto& tile : m_Tiles) {
            float dx = tile.centreX() - centreX;
            float dy = tile.centreY() - centreY;
            // Prioritise by distance + small angle component for spiral
            float angle = std::atan2(dy, dx);
            tile.priority = static_cast<int>(std::sqrt(dx*dx + dy*dy) * 100 + angle * 10);
        }
        
        // Sort by priority (lower = closer to centre)
        std::sort(m_Tiles.begin(), m_Tiles.end(),
                  [](const Tile& a, const Tile& b) { return a.priority < b.priority; });
        
        // Reassign IDs after sorting
        for (size_t i = 0; i < m_Tiles.size(); ++i) {
            m_Tiles[i].id = static_cast<int>(i);
        }
    }
    
    int m_ImageWidth = 0;
    int m_ImageHeight = 0;
    int m_TileSize = 64;
    std::vector<Tile> m_Tiles;
    size_t m_NextTileIndex = 0;
    std::atomic<int> m_CompletedCount{0};
    mutable std::mutex m_Mutex;
};

} // namespace Sirius
