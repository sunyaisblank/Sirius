#pragma once

// Thread-safe RGBA float framebuffer for live preview during rendering, with a
// gamma-corrected 8-bit view for texture upload. Ported from SNDP001A.h.

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <mutex>
#include <optional>
#include <vector>

namespace sirius::render {

// Holds the accumulating render as linear RGBA float; workers write tiles, the
// UI reads an 8-bit gamma-encoded copy.
class DisplayBuffer {
  public:
    // Allocate to the given dimensions and clear to black.
    void Initialise(int width, int height) {
        std::lock_guard<std::mutex> lock(mutex_);
        if (width < 1 || width > 8192 || height < 1 || height > 8192) {
            width_ = 0;
            height_ = 0;
            pixel_data_.clear();
            byte_data_.clear();
            update_counter_ = 0;
            dirty_ = false;
            return;
        }
        width_ = width;
        height_ = height;
        const std::size_t channels =
            static_cast<std::size_t>(width) * static_cast<std::size_t>(height) * 4;
        pixel_data_.resize(channels, 0.0f);  // RGBA float.
        byte_data_.resize(channels, 0);      // RGBA uint8.
        update_counter_ = 0;
        dirty_ = true;
    }

    // Copy a tile region into the buffer (called from a render thread).
    void UpdateTile(int tile_x, int tile_y, int tile_width, int tile_height,
                    const float* tile_data) {
        if (tile_data == nullptr || tile_width <= 0 || tile_height <= 0) return;
        std::lock_guard<std::mutex> lock(mutex_);

        for (int y = 0; y < tile_height; ++y) {
            int dest_y = tile_y + y;
            if (dest_y < 0 || dest_y >= height_) continue;

            for (int x = 0; x < tile_width; ++x) {
                int dest_x = tile_x + x;
                if (dest_x < 0 || dest_x >= width_) continue;

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
        dirty_ = true;
    }

    // Gamma-corrected 8-bit RGBA for GL upload (recomputed when dirty).
    const std::uint8_t* GetByteData(float gamma = 2.2f) {
        if (!std::isfinite(gamma) || gamma <= 0.0f) return nullptr;
        std::lock_guard<std::mutex> lock(mutex_);

        if (dirty_) {
            float inv_gamma = 1.0f / gamma;
            for (std::size_t i = 0; i < pixel_data_.size(); i += 4) {
                byte_data_[i + 0] = FloatToByte(std::pow(pixel_data_[i + 0], inv_gamma));
                byte_data_[i + 1] = FloatToByte(std::pow(pixel_data_[i + 1], inv_gamma));
                byte_data_[i + 2] = FloatToByte(std::pow(pixel_data_[i + 2], inv_gamma));
                byte_data_[i + 3] = FloatToByte(pixel_data_[i + 3]);
            }
            dirty_ = false;
        }

        return byte_data_.data();
    }

    // Raw linear float data (for file output).
    const float* GetFloatData() const { return pixel_data_.data(); }
    std::vector<float>& GetFloatBuffer() { return pixel_data_; }

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
        dirty_ = true;
    }

    int GetWidth() const { return width_; }
    int GetHeight() const { return height_; }
    std::uint64_t GetUpdateCounter() const { return update_counter_.load(); }

  private:
    static std::uint8_t FloatToByte(float v) {
        if (v <= 0.0f) return 0;
        if (v >= 1.0f) return 255;
        return static_cast<std::uint8_t>(v * 255.0f + 0.5f);
    }

    int width_ = 0;
    int height_ = 0;
    std::vector<float> pixel_data_;        // RGBA float [0, inf).
    std::vector<std::uint8_t> byte_data_;  // RGBA uint8 [0, 255].
    std::atomic<std::uint64_t> update_counter_{0};
    bool dirty_ = false;
    mutable std::mutex mutex_;
};

}  // namespace sirius::render
