// SNPR001A.h - Progress Tracker
// Component ID: SNPR001A (Session/Progress)
//
// Tracks render progress, calculates ETA, and provides cancellation token.
// Uses exponential smoothing for more stable ETA estimates.

#pragma once

#include <atomic>
#include <chrono>
#include <functional>
#include <string>
#include <cmath>

namespace Sirius {

//==============================================================================
// Cancellation Token
//==============================================================================
class CancellationToken {
public:
    bool isCancelled() const { return m_Cancelled.load(); }
    void cancel() { m_Cancelled.store(true); }
    void reset() { m_Cancelled.store(false); }
    
private:
    std::atomic<bool> m_Cancelled{false};
};

//==============================================================================
// Progress Callback
//==============================================================================
using ProgressCallback = std::function<void(float progress, int tilesComplete, int tilesTotal, double etaSeconds)>;

//==============================================================================
// Progress Tracker
//==============================================================================
/// @brief Tracks render progress with exponential smoothing for stable ETA.
///
/// Uses exponentially weighted moving average (EWMA) for rate estimation:
///   rate_smooth = α × rate_current + (1-α) × rate_smooth
/// where α is the smoothing factor (higher = more responsive, lower = more stable)
class ProgressTracker {
public:
    /// @brief Start timing
    void start() {
        m_StartTime = std::chrono::high_resolution_clock::now();
        m_LastUpdateTime = m_StartTime;
        m_TilesComplete = 0;
        m_TilesTotal = 0;
        m_SamplesComplete = 0;
        m_SamplesTotal = 0;
        m_SmoothedTilesPerSecond = 0.0;
        m_LastTileCount = 0;
    }

    /// @brief Set total counts
    void setTotals(int tiles, int samplesPerTile) {
        m_TilesTotal = tiles;
        m_SamplesTotal = tiles * samplesPerTile;
    }

    /// @brief Update tile completion
    void completeTile(int samplesInTile = 1) {
        int prevTiles = m_TilesComplete.load();
        m_TilesComplete++;
        m_SamplesComplete += samplesInTile;

        // Update exponential smoothing for ETA calculation
        auto now = std::chrono::high_resolution_clock::now();
        double dt = std::chrono::duration<double>(now - m_LastUpdateTime).count();

        if (dt > 0.1) {  // Update every 100ms minimum to avoid noise
            int tilesProcessed = m_TilesComplete.load() - m_LastTileCount;
            double currentRate = tilesProcessed / dt;

            // Exponential smoothing: α = 0.3 balances responsiveness and stability
            constexpr double SMOOTHING_ALPHA = 0.3;
            if (m_SmoothedTilesPerSecond <= 0) {
                m_SmoothedTilesPerSecond = currentRate;
            } else {
                m_SmoothedTilesPerSecond = SMOOTHING_ALPHA * currentRate +
                                           (1.0 - SMOOTHING_ALPHA) * m_SmoothedTilesPerSecond;
            }

            m_LastUpdateTime = now;
            m_LastTileCount = m_TilesComplete.load();
        }

        if (m_Callback) {
            m_Callback(getProgress(), m_TilesComplete, m_TilesTotal, getETA());
        }
    }

    /// @brief Get progress [0, 1]
    float getProgress() const {
        if (m_TilesTotal == 0) return 0.0f;
        return static_cast<float>(m_TilesComplete) / m_TilesTotal;
    }

    /// @brief Get elapsed time in seconds
    double getElapsedSeconds() const {
        auto now = std::chrono::high_resolution_clock::now();
        return std::chrono::duration<double>(now - m_StartTime).count();
    }

    /// @brief Get estimated time remaining in seconds (exponentially smoothed)
    double getETA() const {
        int remaining = m_TilesTotal - m_TilesComplete.load();
        if (remaining <= 0) return 0.0;

        // Use smoothed rate if available, otherwise fall back to overall average
        double rate = m_SmoothedTilesPerSecond;
        if (rate <= 0) {
            double elapsed = getElapsedSeconds();
            int completed = m_TilesComplete.load();
            if (completed > 0 && elapsed > 0) {
                rate = completed / elapsed;
            }
        }

        if (rate <= 0) return -1.0;  // Cannot estimate yet
        return remaining / rate;
    }
    
    /// @brief Format ETA as string "Xm Ys"
    std::string getETAString() const {
        double eta = getETA();
        if (eta < 0) return "calculating...";
        
        int totalSeconds = static_cast<int>(eta);
        int hours = totalSeconds / 3600;
        int minutes = (totalSeconds % 3600) / 60;
        int seconds = totalSeconds % 60;
        
        std::string result;
        if (hours > 0) result += std::to_string(hours) + "h ";
        if (minutes > 0 || hours > 0) result += std::to_string(minutes) + "m ";
        result += std::to_string(seconds) + "s";
        return result;
    }
    
    /// @brief Get tiles per second rate (exponentially smoothed)
    double getTilesPerSecond() const {
        if (m_SmoothedTilesPerSecond > 0) {
            return m_SmoothedTilesPerSecond;
        }
        // Fallback to overall average
        double elapsed = getElapsedSeconds();
        if (elapsed <= 0) return 0.0;
        return m_TilesComplete / elapsed;
    }
    
    // Accessors
    int getTilesComplete() const { return m_TilesComplete; }
    int getTilesTotal() const { return m_TilesTotal; }
    CancellationToken& getCancellationToken() { return m_CancelToken; }
    
    void setCallback(ProgressCallback callback) { m_Callback = callback; }
    
private:
    std::chrono::high_resolution_clock::time_point m_StartTime;
    std::chrono::high_resolution_clock::time_point m_LastUpdateTime;
    std::atomic<int> m_TilesComplete{0};
    int m_TilesTotal = 0;
    std::atomic<int> m_SamplesComplete{0};
    int m_SamplesTotal = 0;
    CancellationToken m_CancelToken;
    ProgressCallback m_Callback;

    // Exponential smoothing state
    double m_SmoothedTilesPerSecond = 0.0;
    int m_LastTileCount = 0;
};

} // namespace Sirius
