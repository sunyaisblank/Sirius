// KNTA001A.h - Temporal Accumulation and Motion Blur
// Component ID: KNTA001A
// Purpose: Stratified temporal sampling for motion blur
//
// MATHEMATICAL BASIS:
// Motion blur integrates radiance over the shutter interval [t₀, t₀+Δt]:
//   L_blur = (1/Δt) ∫_{t₀}^{t₀+Δt} L(t) dt
//
// We approximate this with stratified sampling:
//   L_blur ≈ (1/N) Σᵢ L(tᵢ) where tᵢ ∈ [t₀ + iΔt/N, t₀ + (i+1)Δt/N]
//
// REFERENCE: James et al. (2015) "DNGR" Section 5.2

#pragma once

#include "../Sirius.Core/Symplectic/MTSB001A.h"
#include <vector>
#include <random>
#include <cmath>

namespace sirius::kernel {

using namespace sirius::spectral;

//==============================================================================
// Shutter Configuration
//==============================================================================

struct ShutterConfig {
    double shutterAngle = 180.0;  // Degrees (180° = 50% duty cycle)
    double frameRate = 24.0;      // Frames per second
    int temporalSamples = 8;      // Samples per pixel for motion blur
    
    // Shutter open time in seconds
    double shutterTime() const {
        return (shutterAngle / 360.0) / frameRate;
    }
    
    // Frame duration in seconds
    double frameDuration() const {
        return 1.0 / frameRate;
    }
};

//==============================================================================
// Temporal Sample
//==============================================================================

struct TemporalSample {
    double t;           // Time offset from frame start [0, shutterTime]
    double weight;      // Sample weight (usually 1/N)
    int stratumIndex;   // Which stratum this sample belongs to
};

//==============================================================================
// TemporalSampler: Generate stratified time samples
//==============================================================================

class TemporalSampler {
public:
    explicit TemporalSampler(const ShutterConfig& config, uint32_t seed = 42)
        : m_config(config), m_rng(seed) {}
    
    //--------------------------------------------------------------------------
    // Generate stratified samples for a pixel
    //--------------------------------------------------------------------------
    
    std::vector<TemporalSample> generateSamples() {
        std::vector<TemporalSample> samples;
        samples.reserve(m_config.temporalSamples);
        
        double shutterTime = m_config.shutterTime();
        double stratumSize = shutterTime / m_config.temporalSamples;
        double weight = 1.0 / m_config.temporalSamples;
        
        std::uniform_real_distribution<double> dist(0.0, 1.0);
        
        for (int i = 0; i < m_config.temporalSamples; ++i) {
            TemporalSample sample;
            sample.stratumIndex = i;
            
            // Stratified jitter within stratum
            double jitter = dist(m_rng);
            sample.t = (i + jitter) * stratumSize;
            sample.weight = weight;
            
            samples.push_back(sample);
        }
        
        return samples;
    }
    
    //--------------------------------------------------------------------------
    // Generate uniform samples (no jitter, for testing)
    //--------------------------------------------------------------------------
    
    std::vector<TemporalSample> generateUniformSamples() {
        std::vector<TemporalSample> samples;
        samples.reserve(m_config.temporalSamples);
        
        double shutterTime = m_config.shutterTime();
        double weight = 1.0 / m_config.temporalSamples;
        
        for (int i = 0; i < m_config.temporalSamples; ++i) {
            TemporalSample sample;
            sample.stratumIndex = i;
            sample.t = (i + 0.5) * shutterTime / m_config.temporalSamples;
            sample.weight = weight;
            samples.push_back(sample);
        }
        
        return samples;
    }
    
private:
    ShutterConfig m_config;
    std::mt19937 m_rng;
};

//==============================================================================
// TemporalAccumulator: Accumulate samples over time
//==============================================================================

class TemporalAccumulator {
public:
    explicit TemporalAccumulator(int numSamples = 8)
        : m_numSamples(numSamples), m_accumulated(SpectralRadiance::zero()), 
          m_sampleCount(0), m_totalWeight(0) {}
    
    //--------------------------------------------------------------------------
    // Add a sample
    //--------------------------------------------------------------------------
    
    void addSample(const SpectralRadiance& radiance, double weight = 1.0) {
        m_accumulated = m_accumulated + radiance * weight;
        m_totalWeight += weight;
        m_sampleCount++;
    }
    
    //--------------------------------------------------------------------------
    // Get normalised result
    //--------------------------------------------------------------------------
    
    SpectralRadiance getResult() const {
        if (m_totalWeight > 0) {
            return m_accumulated * (1.0 / m_totalWeight);
        }
        return SpectralRadiance::zero();
    }
    
    //--------------------------------------------------------------------------
    // Reset accumulator
    //--------------------------------------------------------------------------
    
    void reset() {
        m_accumulated = SpectralRadiance::zero();
        m_sampleCount = 0;
        m_totalWeight = 0;
    }
    
    int sampleCount() const { return m_sampleCount; }
    double totalWeight() const { return m_totalWeight; }
    
private:
    int m_numSamples;
    SpectralRadiance m_accumulated;
    int m_sampleCount;
    double m_totalWeight;
};

//==============================================================================
// Camera Worldline: Interpolate camera state over time
//==============================================================================

struct CameraState {
    double r;       // Radial distance
    double theta;   // Inclination
    double phi;     // Azimuth
    double v_r;     // Radial velocity (dr/dt)
    double v_theta; // Angular velocity in theta
    double v_phi;   // Angular velocity in phi (orbital)
};

class CameraWorldline {
public:
    //--------------------------------------------------------------------------
    // Add keyframe
    //--------------------------------------------------------------------------
    
    void addKeyframe(double t, const CameraState& state) {
        m_keyframes.push_back({t, state});
        // Keep sorted by time
        std::sort(m_keyframes.begin(), m_keyframes.end(),
            [](const auto& a, const auto& b) { return a.first < b.first; });
    }
    
    //--------------------------------------------------------------------------
    // Interpolate camera state at time t
    //--------------------------------------------------------------------------
    
    CameraState interpolate(double t) const {
        if (m_keyframes.empty()) {
            return CameraState{100, M_PI/4, 0, 0, 0, 0};  // Default
        }
        
        if (m_keyframes.size() == 1 || t <= m_keyframes.front().first) {
            return m_keyframes.front().second;
        }
        
        if (t >= m_keyframes.back().first) {
            return m_keyframes.back().second;
        }
        
        // Find surrounding keyframes
        for (size_t i = 0; i < m_keyframes.size() - 1; ++i) {
            if (t >= m_keyframes[i].first && t < m_keyframes[i+1].first) {
                double t0 = m_keyframes[i].first;
                double t1 = m_keyframes[i+1].first;
                double u = (t - t0) / (t1 - t0);  // Normalised [0,1]
                
                const auto& s0 = m_keyframes[i].second;
                const auto& s1 = m_keyframes[i+1].second;
                
                // Linear interpolation (could use Hermite for smoothness)
                CameraState result;
                result.r = s0.r + u * (s1.r - s0.r);
                result.theta = s0.theta + u * (s1.theta - s0.theta);
                result.phi = s0.phi + u * (s1.phi - s0.phi);
                result.v_r = s0.v_r + u * (s1.v_r - s0.v_r);
                result.v_theta = s0.v_theta + u * (s1.v_theta - s0.v_theta);
                result.v_phi = s0.v_phi + u * (s1.v_phi - s0.v_phi);
                
                return result;
            }
        }
        
        return m_keyframes.back().second;
    }
    
    size_t keyframeCount() const { return m_keyframes.size(); }
    
private:
    std::vector<std::pair<double, CameraState>> m_keyframes;
};

//==============================================================================
// Disk Rotation: Keplerian angular velocity
//==============================================================================

class DiskRotation {
public:
    explicit DiskRotation(double M = 1.0, double a = 0.0) 
        : m_M(M), m_a(a) {}
    
    //--------------------------------------------------------------------------
    // Keplerian angular velocity at radius r
    // For Schwarzschild: Ω = √(M/r³)
    // For Kerr: Ω = √(M) / (r^(3/2) + a√M)
    //--------------------------------------------------------------------------
    
    double omega(double r) const {
        if (r <= 0) return 0;
        
        if (std::abs(m_a) < 1e-10) {
            // Schwarzschild
            return std::sqrt(m_M / (r * r * r));
        } else {
            // Kerr (Boyer-Lindquist)
            double sqrtM = std::sqrt(m_M);
            return sqrtM / (std::pow(r, 1.5) + m_a * sqrtM);
        }
    }
    
    //--------------------------------------------------------------------------
    // Disk azimuthal angle at radius r and time t
    //--------------------------------------------------------------------------
    
    double diskPhi(double r, double t, double phi0 = 0) const {
        return phi0 + omega(r) * t;
    }
    
    //--------------------------------------------------------------------------
    // Orbital period at radius r
    //--------------------------------------------------------------------------
    
    double period(double r) const {
        double w = omega(r);
        return (w > 0) ? 2 * M_PI / w : 0;
    }
    
private:
    double m_M;  // Black hole mass
    double m_a;  // Spin parameter
};

} // namespace sirius::kernel
