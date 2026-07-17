#pragma once

// Render progress, ETA estimation, and a cancellation token. Ported from
// SNPR001A.h. ETA uses an exponentially weighted moving average of the tile
// rate so the estimate stays stable under bursty tile completion.

#include <atomic>
#include <chrono>
#include <cmath>
#include <functional>
#include <string>

namespace sirius::render {

// Cooperative cancellation flag shared with the render workers.
class CancellationToken {
  public:
    bool IsCancelled() const { return cancelled_.load(); }
    void Cancel() { cancelled_.store(true); }
    void Reset() { cancelled_.store(false); }

  private:
    std::atomic<bool> cancelled_{false};
};

using ProgressCallback =
    std::function<void(float progress, int tilesComplete, int tilesTotal, double etaSeconds)>;

// Tracks tile completion and derives a smoothed ETA.
class ProgressTracker {
  public:
    // Begin timing and reset counters.
    void Start() {
        start_time_ = std::chrono::high_resolution_clock::now();
        last_update_time_ = start_time_;
        tiles_complete_ = 0;
        tiles_total_ = 0;
        samples_complete_ = 0;
        samples_total_ = 0;
        smoothed_tiles_per_second_ = 0.0;
        last_tile_count_ = 0;
    }

    // Set the total tile and sample counts.
    void SetTotals(int tiles, int samples_per_tile) {
        tiles_total_ = tiles;
        samples_total_ = tiles * samples_per_tile;
    }

    // Record one completed tile; updates the smoothed rate and fires the callback.
    void CompleteTile(int samples_in_tile = 1) {
        tiles_complete_++;
        samples_complete_ += samples_in_tile;

        auto now = std::chrono::high_resolution_clock::now();
        double dt = std::chrono::duration<double>(now - last_update_time_).count();

        if (dt > 0.1) {  // Update at most every 100 ms to suppress noise.
            int tiles_processed = tiles_complete_.load() - last_tile_count_;
            double current_rate = tiles_processed / dt;

            // alpha = 0.3 balances responsiveness against stability.
            constexpr double kSmoothingAlpha = 0.3;
            if (smoothed_tiles_per_second_ <= 0) {
                smoothed_tiles_per_second_ = current_rate;
            } else {
                smoothed_tiles_per_second_ =
                    kSmoothingAlpha * current_rate + (1.0 - kSmoothingAlpha) * smoothed_tiles_per_second_;
            }

            last_update_time_ = now;
            last_tile_count_ = tiles_complete_.load();
        }

        if (callback_) {
            callback_(GetProgress(), tiles_complete_, tiles_total_, GetEta());
        }
    }

    // Progress fraction in [0, 1].
    float GetProgress() const {
        if (tiles_total_ == 0) return 0.0f;
        return static_cast<float>(tiles_complete_) / tiles_total_;
    }

    double GetElapsedSeconds() const {
        auto now = std::chrono::high_resolution_clock::now();
        return std::chrono::duration<double>(now - start_time_).count();
    }

    // Estimated seconds remaining (-1 when not yet estimable).
    double GetEta() const {
        int remaining = tiles_total_ - tiles_complete_.load();
        if (remaining <= 0) return 0.0;

        double rate = smoothed_tiles_per_second_;
        if (rate <= 0) {
            double elapsed = GetElapsedSeconds();
            int completed = tiles_complete_.load();
            if (completed > 0 && elapsed > 0) {
                rate = completed / elapsed;
            }
        }

        if (rate <= 0) return -1.0;
        return remaining / rate;
    }

    // ETA formatted "Xh Ym Zs".
    std::string GetEtaString() const {
        double eta = GetEta();
        if (eta < 0) return "calculating...";

        int total_seconds = static_cast<int>(eta);
        int hours = total_seconds / 3600;
        int minutes = (total_seconds % 3600) / 60;
        int seconds = total_seconds % 60;

        std::string result;
        if (hours > 0) result += std::to_string(hours) + "h ";
        if (minutes > 0 || hours > 0) result += std::to_string(minutes) + "m ";
        result += std::to_string(seconds) + "s";
        return result;
    }

    double GetTilesPerSecond() const {
        if (smoothed_tiles_per_second_ > 0) {
            return smoothed_tiles_per_second_;
        }
        double elapsed = GetElapsedSeconds();
        if (elapsed <= 0) return 0.0;
        return tiles_complete_ / elapsed;
    }

    int GetTilesComplete() const { return tiles_complete_; }
    int GetTilesTotal() const { return tiles_total_; }
    CancellationToken& GetCancellationToken() { return cancel_token_; }

    void SetCallback(ProgressCallback callback) { callback_ = callback; }

  private:
    std::chrono::high_resolution_clock::time_point start_time_;
    std::chrono::high_resolution_clock::time_point last_update_time_;
    std::atomic<int> tiles_complete_{0};
    int tiles_total_ = 0;
    std::atomic<int> samples_complete_{0};
    int samples_total_ = 0;
    CancellationToken cancel_token_;
    ProgressCallback callback_;

    double smoothed_tiles_per_second_ = 0.0;
    int last_tile_count_ = 0;
};

}  // namespace sirius::render
