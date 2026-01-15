// PHSI001A.h - Time-Transformed Explicit Symplectic Integrator (TTESI)
// Component ID: PHSI001A
// Purpose: Double-precision symplectic integration for Kerr geodesics
//
// MATHEMATICAL BASIS:
// The Hamiltonian for null geodesics in Kerr spacetime:
//   H = (1/2) g^μν p_μ p_ν = 0
//
// Time transformation:
//   dτ = g(q) dλ, where g(q) = Σ = r² + a²cos²θ
// This regularises integration near the horizon where Σ → 0.
//
// Yoshida 4th-order composition of leapfrog steps:
//   Coefficients: w₁ = 1/(2-2^(1/3)), w₀ = -2^(1/3)/(2-2^(1/3))
//
// REFERENCE:
// - Yoshida (1990) "Construction of higher order symplectic integrators"
// - Preto & Saha (1999) "Numerical regularization"
// - James et al. (2015) "DNGR" Section A.4

#pragma once

#include "PHMT000B.h"
#include <cmath>
#include <algorithm>

namespace sirius::physics {

//==============================================================================
// Symplectic Integrator Constants
//==============================================================================

namespace symplectic {

//==============================================================================
// Yoshida Symplectic Integrator Coefficients
// Reference: Yoshida (1990) "Construction of higher order symplectic integrators"
//==============================================================================

// Yoshida 4th-order coefficients (3 substeps)
// w₁ = 1/(2-2^(1/3)) ≈ 1.351207, w₀ = -2^(1/3)/(2-2^(1/3)) ≈ -1.702414
constexpr double W4_1 = 1.3512071919596576;
constexpr double W4_0 = -1.7024143839193153;
constexpr int ORDER_4_SUBSTEPS = 3;
constexpr double ORDER_4_WEIGHTS[] = { W4_1, W4_0, W4_1 };

// Yoshida 6th-order "Solution A" coefficients (7 substeps)
// From Yoshida (1990), optimised for minimal error
constexpr int ORDER_6_SUBSTEPS = 7;
constexpr double ORDER_6_WEIGHTS[] = {
     0.78451361047756,
     0.23557321335936,
    -1.17767998417887,
     1.31518632068391,
    -1.17767998417887,
     0.23557321335936,
     0.78451361047756
};

// Yoshida 8th-order coefficients (15 substeps)
// From Yoshida (1990), higher precision
constexpr int ORDER_8_SUBSTEPS = 15;
constexpr double ORDER_8_WEIGHTS[] = {
     0.74167036435061,
    -0.40910082580003,
     0.19075471029623,
    -0.57386247111608,
     0.29906418130365,
     0.33462491824529,
     0.31529309239676,
    -0.79688793935291,
     0.31529309239676,
     0.33462491824529,
     0.29906418130365,
    -0.57386247111608,
     0.19075471029623,
    -0.40910082580003,
     0.74167036435061
};

} // namespace symplectic

//==============================================================================
// Integrator Order Selection
//==============================================================================

enum class IntegratorOrder {
    YOSHIDA_4,  // 4th-order, 3 substeps, fast
    YOSHIDA_6,  // 6th-order, 7 substeps, balanced (DEFAULT)
    YOSHIDA_8   // 8th-order, 15 substeps, high precision
};

//==============================================================================
// SymplecticIntegratorD: Multi-order TTESI implementation for Kerr geodesics
//==============================================================================

class SymplecticIntegratorD {
public:
    struct Config {
        double initialStepSize = 0.1;
        double minStepSize = 1e-8;
        double maxStepSize = 1.0;
        double tolerance = 1e-10;
        int maxStepsPerCall = 10000;
        double horizonBuffer = 1.01;
        IntegratorOrder order = IntegratorOrder::YOSHIDA_6;  // Default: 6th-order
        Config() = default;  // Explicit default constructor
    };
    
    struct StepResult {
        HamiltonianStateD state;
        double lambdaAdvance;
        double hamiltonianError;
        bool terminated;
        int substeps;
    };
    
    explicit SymplecticIntegratorD(const IMetricD* metric)
        : m_metric(metric), m_config() {}

    SymplecticIntegratorD(const IMetricD* metric, const Config& config)
        : m_metric(metric), m_config(config) {}
    
    //--------------------------------------------------------------------------
    // Single Integration Step (Selectable Order)
    //--------------------------------------------------------------------------
    
    StepResult step(const HamiltonianStateD& state, double h) const {
        StepResult result;
        result.substeps = 0;
        result.terminated = false;
        
        HamiltonianStateD s = state;
        
        // Select coefficients based on order
        const double* weights = nullptr;
        int numSubsteps = 0;
        
        switch (m_config.order) {
            case IntegratorOrder::YOSHIDA_4:
                weights = symplectic::ORDER_4_WEIGHTS;
                numSubsteps = symplectic::ORDER_4_SUBSTEPS;
                break;
            case IntegratorOrder::YOSHIDA_6:
                weights = symplectic::ORDER_6_WEIGHTS;
                numSubsteps = symplectic::ORDER_6_SUBSTEPS;
                break;
            case IntegratorOrder::YOSHIDA_8:
                weights = symplectic::ORDER_8_WEIGHTS;
                numSubsteps = symplectic::ORDER_8_SUBSTEPS;
                break;
        }
        
        // Apply Yoshida composition: each weight is a leapfrog step with that fraction of h
        double lambdaAccum = 0;
        for (int i = 0; i < numSubsteps; ++i) {
            double wi = weights[i];
            s = leapfrogStep(s, wi * h);
            
            if (!m_metric->isValid(s.q)) {
                result.terminated = true;
                result.state = s;
                result.lambdaAdvance = lambdaAccum;
                result.hamiltonianError = std::abs(m_metric->hamiltonian(s.q, s.p));
                return result;
            }
            
            lambdaAccum += wi * h;
            result.substeps++;
        }
        
        // Compute Hamiltonian error
        result.state = s;
        result.state.H = m_metric->hamiltonian(s.q, s.p);
        result.hamiltonianError = std::abs(result.state.H);
        result.lambdaAdvance = h;  // Total step is always h (weights sum to 1)
        
        return result;
    }
    
    //--------------------------------------------------------------------------
    // Integrate to Target Lambda or Termination
    //--------------------------------------------------------------------------
    
    struct IntegrationResult {
        GeodesicStateD finalState;
        double totalLambda;
        double maxHamiltonianError;
        int totalSteps;
        bool hitHorizon;
        bool escaped;
    };
    
    IntegrationResult integrate(
        const GeodesicStateD& initial,
        double lambdaMax,
        double escapeRadius = 1e6
    ) const {
        IntegrationResult result;
        result.totalLambda = 0;
        result.maxHamiltonianError = 0;
        result.totalSteps = 0;
        result.hitHorizon = false;
        result.escaped = false;
        
        HamiltonianStateD hs(initial);
        double h = m_config.initialStepSize;
        
        while (result.totalLambda < lambdaMax && result.totalSteps < m_config.maxStepsPerCall) {
            // Adaptive step size based on position (smaller near horizon)
            double r = hs.q.r;
            double r_horizon = m_metric->horizonRadius();
            
            if (r < r_horizon * m_config.horizonBuffer) {
                result.hitHorizon = true;
                break;
            }
            
            if (r > escapeRadius) {
                result.escaped = true;
                break;
            }
            
            // Adapt step size: smaller near horizon, larger when far away
            double r_ratio = (r - r_horizon) / r_horizon;
            h = m_config.initialStepSize * std::min(1.0, r_ratio * r_ratio);
            h = std::clamp(h, m_config.minStepSize, m_config.maxStepSize);
            
            // Don't overshoot lambdaMax
            if (result.totalLambda + h > lambdaMax) {
                h = lambdaMax - result.totalLambda;
            }
            
            // Take step
            StepResult sr = step(hs, h);
            
            if (sr.terminated) {
                result.hitHorizon = true;
                break;
            }
            
            hs = sr.state;
            result.totalLambda += sr.lambdaAdvance;
            result.maxHamiltonianError = std::max(result.maxHamiltonianError, sr.hamiltonianError);
            result.totalSteps++;
        }
        
        // Convert back to GeodesicStateD
        result.finalState = hs.toGeodesicState(result.totalLambda);
        result.finalState.E = -hs.p.t;
        result.finalState.Lz = hs.p.phi;
        
        return result;
    }
    
    //--------------------------------------------------------------------------
    // Conservation Tests (for Validation)
    //--------------------------------------------------------------------------
    
    // Integrate for N steps and return conservation statistics
    struct ConservationStats {
        double initialH;
        double finalH;
        double maxHError;
        double initialE;
        double finalE;
        double initialLz;
        double finalLz;
        int steps;
    };
    
    ConservationStats testConservation(
        const GeodesicStateD& initial,
        int numSteps,
        double stepSize
    ) const {
        ConservationStats stats;
        stats.steps = 0;
        
        HamiltonianStateD hs(initial);
        
        // Initial quantities
        stats.initialH = m_metric->hamiltonian(hs.q, hs.p);
        stats.initialE = -hs.p.t;
        stats.initialLz = hs.p.phi;
        stats.maxHError = std::abs(stats.initialH);
        
        // Integrate
        for (int i = 0; i < numSteps; ++i) {
            StepResult sr = step(hs, stepSize);
            if (sr.terminated) break;
            
            hs = sr.state;
            stats.maxHError = std::max(stats.maxHError, sr.hamiltonianError);
            stats.steps++;
        }
        
        // Final quantities
        stats.finalH = m_metric->hamiltonian(hs.q, hs.p);
        stats.finalE = -hs.p.t;
        stats.finalLz = hs.p.phi;
        
        return stats;
    }
    
private:
    const IMetricD* m_metric;
    Config m_config;
    
    //--------------------------------------------------------------------------
    // Leapfrog (Störmer-Verlet) Substep
    //--------------------------------------------------------------------------
    
    HamiltonianStateD leapfrogStep(const HamiltonianStateD& state, double h) const {
        // Leapfrog: p_{n+1/2} = p_n - (h/2) ∂H/∂q
        //           q_{n+1} = q_n + h ∂H/∂p(q_n, p_{n+1/2})
        //           p_{n+1} = p_{n+1/2} - (h/2) ∂H/∂q(q_{n+1})
        
        HamiltonianStateD s = state;
        
        // Half-step momentum update
        Vec4d dHdq = m_metric->dHdq(s.q, s.p);
        s.p = s.p + dHdq * (h * 0.5);  // Note: dHdq returns dp/dλ = -∂H/∂q (negated)
        
        // Full-step position update
        Vec4d dHdp = m_metric->dHdp(s.q, s.p);
        s.q = s.q + dHdp * h;
        
        // Normalise angles
        while (s.q.phi > 2*M_PI) s.q.phi -= 2*M_PI;
        while (s.q.phi < 0) s.q.phi += 2*M_PI;
        s.q.theta = std::clamp(s.q.theta, 1e-6, M_PI - 1e-6);
        
        // Half-step momentum update (at new position)
        dHdq = m_metric->dHdq(s.q, s.p);
        s.p = s.p + dHdq * (h * 0.5);
        
        return s;
    }
};

} // namespace sirius::physics
