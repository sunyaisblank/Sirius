#pragma once

// Render progress, ETA estimation, and a cancellation token. Ported from
// SNPR001A.h. ETA uses complete elapsed render time so a burst of finished
// tiles cannot discard the time spent tracing that worker wave.

#include <atomic>
#include <chrono>
#include <cstdint>
#include <functional>
#include <mutex>
#include <string>
#include <utility>

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

// Tracks tile completion and estimates remaining time from overall throughput.
class ProgressTracker {
  public:
    using Clock = std::chrono::steady_clock;
    using NowFunction = std::function<Clock::time_point()>;

    explicit ProgressTracker(NowFunction now = Clock::now) : now_(std::move(now)) {}

    // Begin timing and reset counters.
    void Start() {
        std::lock_guard<std::mutex> lock(mutex_);
        start_time_ = now_();
        tiles_complete_ = 0;
        tiles_total_ = 0;
        samples_complete_ = 0;
        samples_total_ = 0;
    }

    // Set the total tile and sample counts.
    void SetTotals(int tiles, int samples_per_tile) {
        std::lock_guard<std::mutex> lock(mutex_);
        tiles_total_ = tiles;
        samples_total_ =
            static_cast<std::int64_t>(tiles) * static_cast<std::int64_t>(samples_per_tile);
    }

    // Record one completed tile and fire the callback.
    void CompleteTile(int samples_in_tile = 1) {
        ProgressCallback callback;
        float progress = 0.0f;
        int complete = 0;
        int total = 0;
        double eta = -1.0;
        {
            std::lock_guard<std::mutex> lock(mutex_);
            ++tiles_complete_;
            samples_complete_ += samples_in_tile;

            const auto now = now_();
            complete = tiles_complete_;
            total = tiles_total_;
            progress = ProgressLocked();
            eta = EtaLocked(now);
            callback = callback_;
        }

        if (callback) {
            try {
                callback(progress, complete, total, eta);
            } catch (...) {
                // Observer callbacks cannot invalidate render progress or
                // terminate a worker thread.
            }
        }
    }

    // Progress fraction in [0, 1].
    float GetProgress() const {
        std::lock_guard<std::mutex> lock(mutex_);
        return ProgressLocked();
    }

    double GetElapsedSeconds() const {
        std::lock_guard<std::mutex> lock(mutex_);
        const auto now = now_();
        return std::chrono::duration<double>(now - start_time_).count();
    }

    // Estimated seconds remaining (-1 when not yet estimable).
    double GetEta() const {
        std::lock_guard<std::mutex> lock(mutex_);
        return EtaLocked(now_());
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
        std::lock_guard<std::mutex> lock(mutex_);
        return RateLocked(now_());
    }

    int GetTilesComplete() const {
        std::lock_guard<std::mutex> lock(mutex_);
        return tiles_complete_;
    }
    int GetTilesTotal() const {
        std::lock_guard<std::mutex> lock(mutex_);
        return tiles_total_;
    }
    CancellationToken& GetCancellationToken() { return cancel_token_; }
    const CancellationToken& GetCancellationToken() const { return cancel_token_; }

    void SetCallback(ProgressCallback callback) {
        std::lock_guard<std::mutex> lock(mutex_);
        callback_ = std::move(callback);
    }

  private:
    [[nodiscard]] float ProgressLocked() const {
        if (tiles_total_ == 0) return 0.0f;
        return static_cast<float>(tiles_complete_) / static_cast<float>(tiles_total_);
    }

    [[nodiscard]] double EtaLocked(Clock::time_point now) const {
        const int remaining = tiles_total_ - tiles_complete_;
        if (remaining <= 0) return 0.0;

        const double rate = RateLocked(now);
        if (rate <= 0) return -1.0;
        return static_cast<double>(remaining) / rate;
    }

    [[nodiscard]] double RateLocked(Clock::time_point now) const {
        const double elapsed = std::chrono::duration<double>(now - start_time_).count();
        return elapsed > 0 ? static_cast<double>(tiles_complete_) / elapsed : 0.0;
    }

    NowFunction now_;
    Clock::time_point start_time_;
    int tiles_complete_ = 0;
    int tiles_total_ = 0;
    std::int64_t samples_complete_ = 0;
    std::int64_t samples_total_ = 0;
    CancellationToken cancel_token_;
    ProgressCallback callback_;

    mutable std::mutex mutex_;
};

}  // namespace sirius::render
