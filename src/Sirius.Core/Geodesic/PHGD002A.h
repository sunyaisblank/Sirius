// =============================================================================
// Sirius.Physics/Core/PHGD002A.h - Numerical Geodesic Integrator
// =============================================================================
// Geodesic integrator for numerical spacetimes using RK4.
// 
// Mathematical Foundation:
//   Geodesic equation: d²x^μ/dλ² + Γ^μ_αβ (dx^α/dλ)(dx^β/dλ) = 0
//   Rewritten as first-order system:
//     dx^μ/dλ = u^μ
//     du^μ/dλ = -Γ^μ_αβ u^α u^β
// =============================================================================
#pragma once

#include <MTTN001A.h>
#include "PHMT000A.h"
#include <functional>

namespace Sirius {
namespace Physics {

// Geodesic state: position and 4-velocity
struct GeodesicState {
    Vec4 x;  // Position (t, r, θ, φ) or (t, x, y, z)
    Vec4 u;  // 4-velocity (dx^μ/dλ)
};

// Integration result
struct IntegrationResult {
    GeodesicState finalState;
    double affineDist;      // Total affine parameter distance
    int steps;              // Number of integration steps
    bool hitHorizon;        // Ray fell into horizon
    bool escaped;           // Ray escaped to infinity
    bool error;             // Numerical error occurred
};

// =============================================================================
// NumericalGeodesicIntegrator - RK4 integration for numerical spacetimes
// =============================================================================
class NumericalGeodesicIntegrator {
public:
    NumericalGeodesicIntegrator();
    ~NumericalGeodesicIntegrator();
    
    // Set the metric to use
    void setMetric(IMetric* metric);
    
    // Integration parameters
    void setMaxSteps(int steps) { m_MaxSteps = steps; }
    void setStepSize(double h) { m_StepSize = h; }
    void setTolerance(double tol) { m_Tolerance = tol; }
    void setMaxDistance(double d) { m_MaxDistance = d; }
    void setHorizonRadius(double r) { m_HorizonRadius = r; }
    
    // Integrate geodesic from initial state
    IntegrationResult integrate(const GeodesicState& initial);
    
    // Single RK4 step
    GeodesicState rk4Step(const GeodesicState& state, double h);
    
    // Compute derivative (geodesic equation RHS)
    GeodesicState computeDerivative(const GeodesicState& state);
    
private:
    IMetric* m_Metric = nullptr;
    
    // Integration parameters
    int m_MaxSteps = 1000;
    double m_StepSize = 0.1;
    double m_Tolerance = 1e-6;
    double m_MaxDistance = 100.0;
    double m_HorizonRadius = 2.0;
    
    // Get radial coordinate from state
    double getRadius(const GeodesicState& state) const;
    
    // Check for NaN/Inf
    bool isValidState(const GeodesicState& state) const;
};

} // namespace Physics
} // namespace Sirius
