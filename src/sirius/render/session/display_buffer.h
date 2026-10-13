#pragma once

// Thread-safe RGBA float framebuffer for live preview during rendering.

#include "sirius/base/contracts.h"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <functional>
#include <mutex>
#include <optional>
#include <span>
#include <vector>

namespace sirius::render {

// Advisory completed physical radiance. The span is borrowed only for the call;
// consumers copy it before returning. Delivery can occur on a tile worker.
using CompletedTileCallback =
    std::function<void(int x, int y, int width, int height, std::span<const float> rgba)>;

// Holds the accumulating render as linear RGBA float; workers write tiles and
// readers take stable snapshots.
class DisplayBuffer {
  public:
    // Allocate to the given dimensions and clear to black.
    void Initialise(int width, int height) {
        SIRIUS_PRE(width >= 1 && width <= 8192 && height >= 1 && height <= 8192);
        std::lock_guard<std::mutex> lock(mutex_);
        width_ = width;
        height_ = height;
        const std::size_t channels =
            static_cast<std::size_t>(width) * static_cast<std::size_t>(height) * 4;
        pixel_data_.resize(channels, 0.0f);  // RGBA float.
        update_counter_ = 0;
    }

    // Copy a tile region into the buffer (called from a render thread).
    void UpdateTile(int tile_x, int tile_y, int tile_width, int tile_height,
                    const float* tile_data) {
        SIRIUS_PRE(tile_data != nullptr);
        SIRIUS_PRE(tile_x >= 0 && tile_y >= 0 && tile_width > 0 && tile_height > 0);
        std::lock_guard<std::mutex> lock(mutex_);
        SIRIUS_PRE(width_ > 0 && height_ > 0);
        SIRIUS_PRE(tile_x < width_ && tile_y < height_);
        SIRIUS_PRE(tile_width <= width_ - tile_x && tile_height <= height_ - tile_y);

        for (int y = 0; y < tile_height; ++y) {
            const int dest_y = tile_y + y;
            for (int x = 0; x < tile_width; ++x) {
                const int dest_x = tile_x + x;

                const std::size_t src_idx = (static_cast<std::size_t>(y) * tile_width + x) * 4;
                const std::size_t dst_idx =
                    (static_cast<std::size_t>(dest_y) * width_ + dest_x) * 4;

                pixel_data_[dst_idx + 0] = tile_data[src_idx + 0];
                pixel_data_[dst_idx + 1] = tile_data[src_idx + 1];
                pixel_data_[dst_idx + 2] = tile_data[src_idx + 2];
                pixel_data_[dst_idx + 3] = tile_data[src_idx + 3];
            }
        }

        update_counter_++;
    }

    // Apply an owning-session transformation while retaining the buffer lock.
    // Raw references are never published across the thread-safety boundary.
    template <typename Callback>
    void MutateFloatData(Callback&& callback) {
        std::lock_guard<std::mutex> lock(mutex_);
        const std::size_t expected_size = pixel_data_.size();
        callback(pixel_data_);
        SIRIUS_POST(pixel_data_.size() == expected_size);
    }

    // Copy one completed rectangle without exposing storage across the lock.
    [[nodiscard]] std::vector<float> SnapshotTileData(int x, int y, int width, int height) const {
        std::lock_guard<std::mutex> lock(mutex_);
        SIRIUS_PRE(x >= 0 && y >= 0 && width > 0 && height > 0);
        SIRIUS_PRE(x < width_ && y < height_ && width <= width_ - x && height <= height_ - y);
        std::vector<float> tile(static_cast<std::size_t>(width) * height * 4);
        for (int row = 0; row < height; ++row) {
            const auto source =
                pixel_data_.begin() + (static_cast<std::size_t>(y + row) * width_ + x) * 4;
            std::copy_n(source, static_cast<std::size_t>(width) * 4,
                        tile.begin() + static_cast<std::size_t>(row) * width * 4);
        }
        return tile;
    }

    // Stable copy for consumers that may overlap a render-thread update.
    [[nodiscard]] std::vector<float> SnapshotFloatData() const {
        std::lock_guard<std::mutex> lock(mutex_);
        return pixel_data_;
    }

    [[nodiscard]] std::optional<std::size_t> FirstNonFiniteIndex() const {
        std::lock_guard<std::mutex> lock(mutex_);
        for (std::size_t i = 0; i < pixel_data_.size(); ++i) {
            if (!std::isfinite(pixel_data_[i])) {
                return i;
            }
        }
        return std::nullopt;
    }

    void Clear() {
        std::lock_guard<std::mutex> lock(mutex_);
        std::fill(pixel_data_.begin(), pixel_data_.end(), 0.0f);
    }

    int GetWidth() const {
        std::lock_guard<std::mutex> lock(mutex_);
        return width_;
    }
    int GetHeight() const {
        std::lock_guard<std::mutex> lock(mutex_);
        return height_;
    }
    std::uint64_t GetUpdateCounter() const { return update_counter_.load(); }

  private:
    int width_ = 0;
    int height_ = 0;
    std::vector<float> pixel_data_;  // RGBA float [0, inf).
    std::atomic<std::uint64_t> update_counter_{0};
    mutable std::mutex mutex_;
};

}  // namespace sirius::render
