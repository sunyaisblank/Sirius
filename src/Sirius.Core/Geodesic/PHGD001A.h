// PHGD001A.h - Geodesic Integration for Relativistic Ray Tracing
//
// Integrates null geodesics: d²xᵘ/dλ² + Γᵘₐᵦ (dxᵅ/dλ)(dxᵝ/dλ) = 0
// Constraint: gₘᵥ kᵘ kᵛ = 0 (null condition)
//
// Integration: RK45/Dormand-Prince with adaptive step control
// Termination: horizon (r < r_s), escape (r > r_max), or NaN/Inf
//
// Tests: TSPH003A.cpp (unit), TSDG003A.cpp (conservation), TSIN001A.cpp (paths)

#pragma once
#include <MTTN001A.h>
#include "PHMT000A.h"

// =============================================================================
// Lightray - State vector for null geodesic integration
// =============================================================================
struct Lightray {
    Vec4 position;           ///< 4-position (t, r, θ, φ) or (t, x, y, z)
    Vec4 velocity;           ///< 4-velocity kᵘ = dxᵘ/dλ
    Vec4 acceleration;       ///< 4-acceleration d²xᵘ/dλ²
    float ku_uobsu;          ///< Initial k·u (for redshift calculation)
    float proper_time;       ///< Accumulated affine parameter λ
    float coordinate_time;   ///< Accumulated coordinate time t
    float running_dlambda_dnew;  ///< Jacobian for intensity calculation
    int terminated;          ///< Termination status (0 = active)
    int sx, sy;              ///< Screen pixel coordinates
    int bounce_count;        ///< Number of reflections/interactions
    float step_size;         ///< Current adaptive step size
    int padding;             ///< Alignment padding
};

// =============================================================================
// ObserverState - Local observer reference frame
// =============================================================================
struct ObserverState {
    Vec4 position;   ///< Observer 4-position
    Vec4 velocity;   ///< Observer 4-velocity (normalized: g(u,u) = -1)
    Vec4 e0, e1, e2, e3;  ///< Orthonormal tetrad (e0 = u, e1-e3 spatial)
    bool is_timelike;     ///< True if velocity is timelike
};

// =============================================================================
// IntegratorConfig - Configuration for geodesic integration
// =============================================================================
struct IntegratorConfig {
    float abs_tolerance = 1e-6f;   ///< Absolute error tolerance
    float rel_tolerance = 1e-6f;   ///< Relative error tolerance  
    float min_step = 1e-6f;        ///< Minimum step size
    float max_step = 0.1f;         ///< Maximum step size
    float initial_step = 0.01f;    ///< Initial step size
    float safety_factor = 0.9f;    ///< Safety factor for step adaptation (0.8-0.95)
    float step_grow_max = 2.0f;    ///< Maximum step growth factor
    float step_shrink_min = 0.1f;  ///< Minimum step shrink factor
    bool use_rk45 = true;          ///< Use RK45 (true) or RK4 (false)
};

// =============================================================================
// RK45State - Internal state for RK45 integrator
// =============================================================================
struct RK45State {
    Vec4 x;          ///< Position
    Vec4 p;          ///< Covariant momentum
    Vec4 k;          ///< Contravariant velocity (derived from p)
    float error;     ///< Estimated local truncation error
    bool accepted;   ///< Whether step was accepted
};

// =============================================================================
// Geodesic - Static methods for geodesic integration
// =============================================================================
class Geodesic {
public:
    /// @brief Integrate one step of geodesic equation
    /// @param ray [in/out] Light ray state
    /// @param metric Spacetime metric
    /// @param min_step Minimum allowed step size (default 1e-6)
    /// @param max_step Maximum allowed step size (default 0.1)
    /// @return true if step succeeded, false if failed/terminated
    /// PRECONDITION: ray.position in valid domain, metric != nullptr
    /// POSTCONDITION: ray.velocity satisfies null condition
    static bool integrateStep(Lightray& ray, IMetric* metric, float min_step = 1e-6f, float max_step = 0.1f);
    
    /// @brief Compute geodesic acceleration from Christoffel symbols
    /// @return aᵘ = -Γᵘₐᵦ vᵅ vᵝ
    /// PRECONDITION: position in valid domain
    /// POSTCONDITION: finite acceleration vector
    static Vec4 calculateAcceleration(const Vec4& velocity, const Vec4& position, IMetric* metric);
    
    /// @brief Check if ray should be terminated
    static bool checkTermination(const Lightray& ray, IMetric* metric);
    
    /// @brief Calculate gravitational redshift z = (λ_obs - λ_emit)/λ_emit
    /// @return Redshift factor (>0 for redshift, <0 for blueshift)
    static float calculateRedshift(const Lightray& ray, const ObserverState& observer, IMetric* metric);
    
    /// @brief Create and initialize observer state
    /// @param position Observer 4-position
    /// @param velocity Observer 4-velocity (will be normalized)
    /// POSTCONDITION: observer.velocity satisfies g(u,u) = -1
    static ObserverState createObserver(const Vec4& position, const Vec4& velocity, IMetric* metric);
    
    /// @brief Construct orthonormal tetrad for observer via Gram-Schmidt
    /// POSTCONDITION: g(eₐ, eᵦ) = ηₐᵦ (Minkowski metric in tetrad frame)
    static void calculateTetrads(ObserverState& observer, IMetric* metric);
    
    // =========================================================================
    // RK45 Adaptive Integration (Phase 2 Performance Optimization)
    // =========================================================================
    
    /// @brief Integrate one step using RK45 (Dormand-Prince) with error control
    /// @param ray [in/out] Light ray state
    /// @param metric Spacetime metric
    /// @param config Integration configuration (tolerances, step limits)
    /// @return true if step succeeded, false if failed/terminated
    /// MATHEMATICAL BASIS: Embedded 4th/5th order Runge-Kutta pair
    /// ERROR ESTIMATE: ||y5 - y4|| used for step adaptation
    static bool integrateStepRK45(Lightray& ray, IMetric* metric, const IntegratorConfig& config);
    
    /// @brief Compute optimal step size from error estimate
    /// FORMULA: h_new = h * safety * (tol/err)^(1/5)
    static float computeOptimalStep(float h, float error, float tolerance, 
                                     const IntegratorConfig& config);
    
    /// @brief Get default integrator configuration
    static IntegratorConfig getDefaultConfig();
};
