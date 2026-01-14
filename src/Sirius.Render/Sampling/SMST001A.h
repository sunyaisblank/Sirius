// SMST001A.h - Stratified Sampling
// Component ID: SMST001A (Sampling/Stratified)
//
// Stratified sampling divides the sample domain into strata
// and samples each stratum, reducing variance vs pure random.
//
// MATHEMATICAL BASIS:
// For 1D stratified sampling with N strata:
//   x_i = (i + ξ_i) / N,  where ξ_i ∈ [0,1) is uniform random
//
// Variance reduction factor ≈ 1/N compared to uniform.
//
// For 2D: Divide into N×M grid, jitter within each cell.

#pragma once

#include <vector>
#include <random>
#include <cmath>

namespace Sirius::Sampling {

//==============================================================================
// Random Number Generator (thread-local)
//==============================================================================
class RNG {
public:
    explicit RNG(uint64_t seed = 0) : m_Gen(seed) {}
    
    /// @brief Uniform [0, 1)
    double uniform() {
        return m_Dist(m_Gen);
    }
    
    /// @brief Uniform [a, b)
    double uniform(double a, double b) {
        return a + (b - a) * uniform();
    }
    
    /// @brief Integer [0, n)
    int uniformInt(int n) {
        return static_cast<int>(uniform() * n);
    }
    
    /// @brief Reseed
    void seed(uint64_t s) { m_Gen.seed(s); }
    
private:
    std::mt19937_64 m_Gen;
    std::uniform_real_distribution<double> m_Dist{0.0, 1.0};
};

//==============================================================================
// 2D Sample Point
//==============================================================================
struct Sample2D {
    double u, v;  ///< Coordinates in [0,1)²
};

//==============================================================================
// Stratified Sampler
//==============================================================================
class StratifiedSampler {
public:
    /// @brief Create sampler with N×M strata
    StratifiedSampler(int strataX, int strataY, uint64_t seed = 0)
        : m_StrataX(strataX)
        , m_StrataY(strataY)
        , m_RNG(seed)
    {}
    
    /// @brief Generate all stratified samples
    std::vector<Sample2D> generate() {
        std::vector<Sample2D> samples;
        samples.reserve(m_StrataX * m_StrataY);
        
        double invX = 1.0 / m_StrataX;
        double invY = 1.0 / m_StrataY;
        
        for (int y = 0; y < m_StrataY; ++y) {
            for (int x = 0; x < m_StrataX; ++x) {
                Sample2D s;
                s.u = (x + m_RNG.uniform()) * invX;
                s.v = (y + m_RNG.uniform()) * invY;
                samples.push_back(s);
            }
        }
        
        return samples;
    }
    
    /// @brief Generate single sample for stratum (x, y)
    Sample2D sampleStratum(int x, int y) {
        Sample2D s;
        s.u = (x + m_RNG.uniform()) / m_StrataX;
        s.v = (y + m_RNG.uniform()) / m_StrataY;
        return s;
    }
    
    /// @brief Total sample count
    int totalSamples() const { return m_StrataX * m_StrataY; }
    
    /// @brief Reset RNG
    void reset(uint64_t seed) { m_RNG.seed(seed); }
    
private:
    int m_StrataX, m_StrataY;
    RNG m_RNG;
};

//==============================================================================
// Latin Hypercube Sampling
//==============================================================================
class LatinHypercubeSampler {
public:
    explicit LatinHypercubeSampler(int n, uint64_t seed = 0)
        : m_N(n), m_RNG(seed) {}
    
    /// @brief Generate N samples with Latin hypercube property
    std::vector<Sample2D> generate() {
        std::vector<Sample2D> samples(m_N);
        std::vector<int> permU(m_N), permV(m_N);
        
        // Initialize permutations
        for (int i = 0; i < m_N; ++i) {
            permU[i] = i;
            permV[i] = i;
        }
        
        // Shuffle permutations
        for (int i = m_N - 1; i > 0; --i) {
            std::swap(permU[i], permU[m_RNG.uniformInt(i + 1)]);
            std::swap(permV[i], permV[m_RNG.uniformInt(i + 1)]);
        }
        
        // Generate samples
        double inv = 1.0 / m_N;
        for (int i = 0; i < m_N; ++i) {
            samples[i].u = (permU[i] + m_RNG.uniform()) * inv;
            samples[i].v = (permV[i] + m_RNG.uniform()) * inv;
        }
        
        return samples;
    }
    
private:
    int m_N;
    RNG m_RNG;
};

} // namespace Sirius::Sampling
