// SRPG001A.h - ETA Calculator (Progress)
// Component ID: SRPG001A
// Estimates remaining render time based on tile completion times
//
// Governance: docs/specification.md

#ifndef SIRIUS_RENDER_SRPG001A_H
#define SIRIUS_RENDER_SRPG001A_H

#include <algorithm>
#include <cmath>
#include <vector>

namespace sirius::render {

//==============================================================================
// ETACalculator
// Tracks tile completion times and estimates remaining render time
// Uses exponentially weighted moving average to favor recent samples
//==============================================================================

class ETACalculator {
public:
    static constexpr size_t MAX_SAMPLES = 64;

    ETACalculator() = default;

    //--------------------------------------------------------------------------
    // Sample Recording
    //--------------------------------------------------------------------------

    /// Record a tile completion event
    /// @param timestamp Time in seconds since render start
    void recordTileCompletion(double timestamp) {
        m_samples.push_back(timestamp);
        // Keep only most recent samples
        if (m_samples.size() > MAX_SAMPLES) {
            m_samples.erase(m_samples.begin());
        }
    }

    //--------------------------------------------------------------------------
    // Query
    //--------------------------------------------------------------------------

    /// Number of recorded samples
    size_t sampleCount() const {
        return m_samples.size();
    }

    /// Elapsed time from first to last sample
    double elapsedSeconds() const {
        if (m_samples.size() < 2) return 0.0;
        return m_samples.back() - m_samples.front();
    }

    /// Current tiles per second rate (weighted toward recent)
    double tilesPerSecond() const {
        if (m_samples.size() < 2) return 0.0;

        // Use recent window for rate calculation (favor recent speed)
        // Smaller window = more responsive to recent changes
        size_t windowSize = std::min(m_samples.size(), static_cast<size_t>(10));
        size_t startIdx = m_samples.size() - windowSize;

        double windowDuration = m_samples.back() - m_samples[startIdx];
        if (windowDuration <= 0.0) return 0.0;

        return static_cast<double>(windowSize - 1) / windowDuration;
    }

    /// Estimate remaining time in seconds
    /// @param completedTiles Number of tiles already rendered
    /// @param totalTiles Total tiles to render
    /// @return Estimated seconds remaining (0 if insufficient data)
    double estimateRemaining(int completedTiles, int totalTiles) const {
        double rate = tilesPerSecond();
        if (rate <= 0.0) return 0.0;

        int remaining = totalTiles - completedTiles;
        if (remaining <= 0) return 0.0;

        return static_cast<double>(remaining) / rate;
    }

    //--------------------------------------------------------------------------
    // Control
    //--------------------------------------------------------------------------

    /// Reset all recorded data
    void reset() {
        m_samples.clear();
    }

private:
    std::vector<double> m_samples;
};

} // namespace sirius::render

#endif // SIRIUS_RENDER_SRPG001A_H
