#pragma once

// Tile decomposition and spiral (centre-outward) scheduling. Ported from
// SNTL001A.h. Tiles are prioritised by distance from the image centre plus a
// small angular term. Optional work
// groups keep neighboring tiles with one worker without changing their bounds
// or publication ledger; the default retains the original spiral sequence.

#include "sirius/base/contracts.h"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <mutex>
#include <span>
#include <vector>

namespace sirius::render {

// Per-tile lifecycle state.
enum class TileState : std::uint8_t {
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
    void Initialise(int image_width, int image_height, int tile_size = 64,
                    int minimum_work_edge = 0) {
        SIRIUS_PRE(image_width > 0 && image_height > 0 && tile_size > 0);
        SIRIUS_PRE(minimum_work_edge >= 0 &&
                   (minimum_work_edge <= tile_size || minimum_work_edge % tile_size == 0));
        std::lock_guard<std::mutex> lock(mutex_);

        image_width_ = image_width;
        image_height_ = image_height;
        tile_size_ = tile_size;
        work_edge_ = std::max(tile_size, minimum_work_edge);
        tiles_.clear();
        next_tile_index_ = 0;
        completed_count_ = 0;

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
        const auto tiles = AcquireTiles(false);
        return tiles.empty() ? nullptr : tiles.data();
    }

    // One worker exclusively owns this group until its tiles finish. Storage
    // belongs to the scheduler and remains stable until Initialise/Reset.
    std::span<Tile> GetNextTileGroup() { return AcquireTiles(true); }

    void CompleteTile(int tile_id) {
        std::lock_guard<std::mutex> lock(mutex_);
        if (tile_id >= 0 && static_cast<std::size_t>(tile_id) < tiles_.size()) {
            auto& tile = tiles_[static_cast<std::size_t>(tile_id)];
            if (tile.state == TileState::Complete) return;
            tile.state = TileState::Complete;
            completed_count_++;
        }
    }

    void FailTile(int tile_id) {
        std::lock_guard<std::mutex> lock(mutex_);
        if (tile_id >= 0 && static_cast<std::size_t>(tile_id) < tiles_.size()) {
            auto& tile = tiles_[static_cast<std::size_t>(tile_id)];
            if (tile.state == TileState::Complete) --completed_count_;
            tile.state = TileState::Failed;
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
    std::span<Tile> AcquireTiles(bool grouped) {
        std::lock_guard<std::mutex> lock(mutex_);
        while (next_tile_index_ < tiles_.size() &&
               tiles_[next_tile_index_].state != TileState::Pending)
            ++next_tile_index_;
        if (next_tile_index_ == tiles_.size()) return {};
        const auto first = next_tile_index_++;
        const auto& anchor = tiles_[first];
        if (grouped && work_edge_ > tile_size_)
            while (next_tile_index_ < tiles_.size()) {
                const auto& next = tiles_[next_tile_index_];
                if (next.state != TileState::Pending ||
                    next.x / work_edge_ != anchor.x / work_edge_ ||
                    next.y / work_edge_ != anchor.y / work_edge_)
                    break;
                ++next_tile_index_;
            }
        auto group = std::span(tiles_).subspan(first, next_tile_index_ - first);
        for (auto& tile : group) tile.state = TileState::Active;
        return group;
    }

    // Sort tiles by distance-from-centre plus a small angular term (spiral),
    // then reassign ids so id matches render order.
    void SortSpiralOrder() {
        float centre_x = image_width_ / 2.0f;
        float centre_y = image_height_ / 2.0f;

        for (auto& tile : tiles_) {
            float x = static_cast<float>(tile.CentreX());
            float y = static_cast<float>(tile.CentreY());
            if (work_edge_ > tile_size_) {
                const int left = (tile.x / work_edge_) * work_edge_;
                const int top = (tile.y / work_edge_) * work_edge_;
                x = left + std::min(work_edge_, image_width_ - left) * .5f;
                y = top + std::min(work_edge_, image_height_ - top) * .5f;
            }
            const float dx = x - centre_x;
            const float dy = y - centre_y;
            float angle = std::atan2(dy, dx);
            tile.priority = static_cast<int>(std::sqrt(dx * dx + dy * dy) * 100 + angle * 10);
        }

        std::sort(tiles_.begin(), tiles_.end(), [&](const Tile& a, const Tile& b) {
            if (work_edge_ <= tile_size_ || a.priority != b.priority)
                return a.priority < b.priority;
            // Equal integer priorities must not interleave different groups.
            if (a.y / work_edge_ != b.y / work_edge_) return a.y / work_edge_ < b.y / work_edge_;
            if (a.x / work_edge_ != b.x / work_edge_) return a.x / work_edge_ < b.x / work_edge_;
            return a.y == b.y ? a.x < b.x : a.y < b.y;
        });

        for (std::size_t i = 0; i < tiles_.size(); ++i) {
            tiles_[i].id = static_cast<int>(i);
        }
    }

    int image_width_ = 0;
    int image_height_ = 0;
    int tile_size_ = 64;
    int work_edge_ = 64;
    std::vector<Tile> tiles_;
    std::size_t next_tile_index_ = 0;
    std::atomic<int> completed_count_{0};
    mutable std::mutex mutex_;
};

}  // namespace sirius::render
