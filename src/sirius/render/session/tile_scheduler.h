#pragma once

// Tile decomposition and spiral (centre-outward) scheduling. Ported from
// SNTL001A.h. The spiral order is preserved exactly: tiles are prioritised by
// distance from the image centre plus a small angular term, so the reference
// path visits tiles in the identical sequence.

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <mutex>
#include <vector>

namespace sirius::render {

// Per-tile lifecycle state.
enum class TileState : uint8_t {
    Pending,   // Not yet started.
    Active,    // Currently rendering.
    Complete,  // Finished successfully.
    Failed     // Error occurred.
};

// One rectangular tile of the image.
struct Tile {
    int id;
    int x, y;           // Top-left corner (pixels).
    int width, height;  // Tile dimensions.
    TileState state = TileState::Pending;
    int priority = 0;  // Lower = higher priority (rendered earlier).

    int CentreX() const { return x + width / 2; }
    int CentreY() const { return y + height / 2; }
    int PixelCount() const { return width * height; }
};

// Generates tiles and hands them out in spiral order (thread-safe).
class TileScheduler {
  public:
    // Build the tile grid and sort it into spiral order.
    void Initialise(int image_width, int image_height, int tile_size = 64) {
        std::lock_guard<std::mutex> lock(mutex_);

        image_width_ = image_width;
        image_height_ = image_height;
        tile_size_ = tile_size;
        tiles_.clear();
        next_tile_index_ = 0;

        int tiles_x = (image_width + tile_size - 1) / tile_size;
        int tiles_y = (image_height + tile_size - 1) / tile_size;
        int tile_id = 0;

        for (int ty = 0; ty < tiles_y; ++ty) {
            for (int tx = 0; tx < tiles_x; ++tx) {
                Tile tile;
                tile.id = tile_id++;
                tile.x = tx * tile_size;
                tile.y = ty * tile_size;
                tile.width = std::min(tile_size, image_width - tile.x);
                tile.height = std::min(tile_size, image_height - tile.y);
                tile.state = TileState::Pending;
                tiles_.push_back(tile);
            }
        }

        SortSpiralOrder();
    }

    // Next pending tile (marked Active), or nullptr when none remain.
    Tile* GetNextTile() {
        std::lock_guard<std::mutex> lock(mutex_);

        for (size_t i = next_tile_index_; i < tiles_.size(); ++i) {
            if (tiles_[i].state == TileState::Pending) {
                tiles_[i].state = TileState::Active;
                next_tile_index_ = i + 1;
                return &tiles_[i];
            }
        }
        return nullptr;
    }

    void CompleteTile(int tile_id) {
        std::lock_guard<std::mutex> lock(mutex_);
        for (auto& tile : tiles_) {
            if (tile.id == tile_id) {
                tile.state = TileState::Complete;
                completed_count_++;
                return;
            }
        }
    }

    void FailTile(int tile_id) {
        std::lock_guard<std::mutex> lock(mutex_);
        for (auto& tile : tiles_) {
            if (tile.id == tile_id) {
                tile.state = TileState::Failed;
                return;
            }
        }
    }

    void Reset() {
        std::lock_guard<std::mutex> lock(mutex_);
        for (auto& tile : tiles_) {
            tile.state = TileState::Pending;
        }
        next_tile_index_ = 0;
        completed_count_ = 0;
    }

    int GetTileCount() const { return static_cast<int>(tiles_.size()); }
    int GetCompletedCount() const { return completed_count_.load(); }
    int GetPendingCount() const { return GetTileCount() - GetCompletedCount(); }
    float GetProgress() const {
        int total = GetTileCount();
        return total > 0 ? static_cast<float>(GetCompletedCount()) / total : 0.0f;
    }
    bool AllComplete() const { return GetCompletedCount() == GetTileCount(); }

    const std::vector<Tile>& GetTiles() const { return tiles_; }

  private:
    // Sort tiles by distance-from-centre plus a small angular term (spiral),
    // then reassign ids so id matches render order.
    void SortSpiralOrder() {
        float centre_x = image_width_ / 2.0f;
        float centre_y = image_height_ / 2.0f;

        for (auto& tile : tiles_) {
            float dx = tile.CentreX() - centre_x;
            float dy = tile.CentreY() - centre_y;
            float angle = std::atan2(dy, dx);
            tile.priority = static_cast<int>(std::sqrt(dx * dx + dy * dy) * 100 + angle * 10);
        }

        std::sort(tiles_.begin(), tiles_.end(),
                  [](const Tile& a, const Tile& b) { return a.priority < b.priority; });

        for (size_t i = 0; i < tiles_.size(); ++i) {
            tiles_[i].id = static_cast<int>(i);
        }
    }

    int image_width_ = 0;
    int image_height_ = 0;
    int tile_size_ = 64;
    std::vector<Tile> tiles_;
    size_t next_tile_index_ = 0;
    std::atomic<int> completed_count_{0};
    mutable std::mutex mutex_;
};

}  // namespace sirius::render
