// KNTA001A.h - Temporal Sampling / Motion Blur
// Component ID: KNTA001A
// Temporal sampling for motion blur rendering
//
// Governance: docs/specification.md

#ifndef SIRIUS_RENDER_KNTA001A_H
#define SIRIUS_RENDER_KNTA001A_H

#include <algorithm>
#include <cmath>
#include <map>
#include <random>
#include <vector>

#include <MTSB001A.h>  // For SpectralRadiance

namespace sirius::kernel {

//==============================================================================
// ShutterConfig
// Configuration for temporal sampling / motion blur
//==============================================================================

struct ShutterConfig {
    double shutterAngle = 180.0;  // Shutter angle in degrees
    double frameRate = 24.0;      // Frame rate in fps
    int temporalSamples = 8;      // Number of temporal samples

    /// Compute shutter open time in seconds
    double shutterTime() const {
        return (shutterAngle / 360.0) / frameRate;
    }
};

//==============================================================================
// TemporalSample
// A single temporal sample within the shutter interval
//==============================================================================

struct TemporalSample {
    double t = 0.0;       // Time offset from frame start
    double weight = 1.0;  // Sample weight
    int stratumIndex = 0; // Index of the stratum this sample belongs to
};

//==============================================================================
// TemporalSampler
// Generates temporal samples for motion blur
//==============================================================================

class TemporalSampler {
public:
    explicit TemporalSampler(const ShutterConfig& config)
        : m_config(config), m_rng(42) {}

    /// Generate stratified jittered samples
    std::vector<TemporalSample> generateSamples() {
        std::vector<TemporalSample> samples;
        samples.reserve(m_config.temporalSamples);

        double shutterTime = m_config.shutterTime();
        double stratum = shutterTime / m_config.temporalSamples;
        std::uniform_real_distribution<double> dist(0.0, stratum);

        for (int i = 0; i < m_config.temporalSamples; ++i) {
            TemporalSample s;
            s.t = i * stratum + dist(m_rng);
            s.weight = 1.0 / m_config.temporalSamples;
            s.stratumIndex = i;
            samples.push_back(s);
        }

        return samples;
    }

    /// Generate uniform (non-jittered) samples
    std::vector<TemporalSample> generateUniformSamples() const {
        std::vector<TemporalSample> samples;
        samples.reserve(m_config.temporalSamples);

        double shutterTime = m_config.shutterTime();
        double spacing = shutterTime / m_config.temporalSamples;

        for (int i = 0; i < m_config.temporalSamples; ++i) {
            TemporalSample s;
            s.t = (i + 0.5) * spacing;  // Centre of each stratum
            s.weight = 1.0 / m_config.temporalSamples;
            s.stratumIndex = i;
            samples.push_back(s);
        }

        return samples;
    }

private:
    ShutterConfig m_config;
    std::mt19937 m_rng;
};

//==============================================================================
// CameraState
// Represents camera position/orientation at a specific time
//==============================================================================

struct CameraState {
    double t = 0.0;       // Time parameter
    double r = 100.0;     // Radial coordinate (Boyer-Lindquist)
    double theta = M_PI / 2.0;  // Polar angle
    double phi = 0.0;     // Azimuthal angle

    /// Linear interpolation between two states
    static CameraState lerp(const CameraState& a, const CameraState& b, double t) {
        CameraState result;
        result.t = a.t + t * (b.t - a.t);
        result.r = a.r + t * (b.r - a.r);
        result.theta = a.theta + t * (b.theta - a.theta);
        result.phi = a.phi + t * (b.phi - a.phi);
        return result;
    }
};

//==============================================================================
// CameraWorldline
// Camera path through spacetime for motion blur
//==============================================================================

class CameraWorldline {
public:
    /// Add a keyframe at specified time
    void addKeyframe(double t, const CameraState& state) {
        m_keyframes[t] = state;
        m_keyframes[t].t = t;
    }

    /// Interpolate camera state at time t
    CameraState interpolate(double t) const {
        if (m_keyframes.empty()) {
            return CameraState();
        }

        // Find bounding keyframes
        auto upper = m_keyframes.upper_bound(t);
        if (upper == m_keyframes.begin()) {
            return upper->second;
        }
        if (upper == m_keyframes.end()) {
            return std::prev(upper)->second;
        }

        auto lower = std::prev(upper);
        double t0 = lower->first;
        double t1 = upper->first;
        double u = (t - t0) / (t1 - t0);

        return CameraState::lerp(lower->second, upper->second, u);
    }

    /// Get number of keyframes
    size_t keyframeCount() const { return m_keyframes.size(); }

private:
    std::map<double, CameraState> m_keyframes;
};

//==============================================================================
// DiskRotation
// Computes Keplerian orbital motion for disk material
//==============================================================================

class DiskRotation {
public:
    DiskRotation() : m_M(1.0), m_a(0.0) {}
    DiskRotation(double M, double a) : m_M(M), m_a(a) {}

    /// Angular velocity at radius r (Keplerian)
    /// Ω = sqrt(M) / (r^(3/2) + a*sqrt(M))
    double omega(double r) const {
        double sqrtM = std::sqrt(m_M);
        return sqrtM / (std::pow(r, 1.5) + m_a * sqrtM);
    }

    /// Orbital period at radius r
    double period(double r) const {
        return 2.0 * M_PI / omega(r);
    }

    /// Compute phi evolution over time dt
    double phiEvolution(double r, double dt) const {
        return omega(r) * dt;
    }

    /// Compute disk phi at time t (absolute position)
    double diskPhi(double r, double t) const {
        return omega(r) * t;
    }

private:
    double m_M;  // Black hole mass parameter
    double m_a;  // Black hole spin parameter
};

//==============================================================================
// TemporalAccumulator
// Accumulates weighted spectral samples across temporal samples
//==============================================================================

class TemporalAccumulator {
public:
    /// Add a weighted sample
    void addSample(const sirius::spectral::SpectralRadiance& sample, double weight = 1.0) {
        if (m_sampleCount == 0) {
            m_accumulated = sample * weight;
        } else {
            m_accumulated = m_accumulated + sample * weight;
        }
        m_totalWeight += weight;
        ++m_sampleCount;
    }

    /// Get final accumulated radiance (weighted average)
    sirius::spectral::SpectralRadiance result() const {
        if (m_totalWeight > 0) {
            return m_accumulated * (1.0 / m_totalWeight);
        }
        return sirius::spectral::SpectralRadiance();
    }

    /// Alias for result() for API compatibility
    sirius::spectral::SpectralRadiance getResult() const { return result(); }

    /// Reset accumulator
    void reset() {
        m_accumulated = sirius::spectral::SpectralRadiance();
        m_totalWeight = 0.0;
        m_sampleCount = 0;
    }

    /// Get sample count
    int sampleCount() const { return m_sampleCount; }

    /// Get total weight
    double totalWeight() const { return m_totalWeight; }

private:
    sirius::spectral::SpectralRadiance m_accumulated;
    double m_totalWeight = 0.0;
    int m_sampleCount = 0;
};

} // namespace sirius::kernel

namespace sirius::spectral {
    // Allow tests to reference spectral namespace
    using sirius::kernel::ShutterConfig;
    using sirius::kernel::TemporalSample;
    using sirius::kernel::TemporalSampler;
    using sirius::kernel::CameraState;
    using sirius::kernel::CameraWorldline;
    using sirius::kernel::DiskRotation;
    using sirius::kernel::TemporalAccumulator;
}

#endif // SIRIUS_RENDER_KNTA001A_H
