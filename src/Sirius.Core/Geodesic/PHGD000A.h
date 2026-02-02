// =============================================================================
// PHGD000A.h - Geodesic Integrator Interface
// Component ID: PHGD000A (Physics/Geodesic/Interface)
// =============================================================================
//
// PURPOSE:
// Defines the common interface for all geodesic integrators in Sirius.
// This allows the rendering pipeline to work with different integration
// schemes (RK4, RK45, symplectic) through a unified interface.
//
// IMPLEMENTATIONS:
// - PHGD001A: RK45 Dormand-Prince (adaptive, high accuracy)
// - PHGD002A: RK4 Fixed-step (simple, predictable)
// - GPU: RDOP002A uses symplectic Verlet (energy-conserving)
//
// MATHEMATICAL FOUNDATION:
// Geodesic equation: d²x^μ/dλ² + Γ^μ_αβ (dx^α/dλ)(dx^β/dλ) = 0
// Constraint: g_μν k^μ k^ν = 0 (null condition for light rays)
//
// TESTS: TSPH003A.cpp, TSDG003A.cpp (conservation)
// =============================================================================

#ifndef PHGD000A_H
#define PHGD000A_H

#include "MTTN001A.h"
#include "PHMT000A.h"
#include <optional>

namespace Sirius {

// =============================================================================
// Integration Configuration
// =============================================================================

/// @brief Configuration for geodesic integration
struct IntegratorConfig {
    // Error tolerances (for adaptive methods)
    float abs_tolerance = 1e-6f;   ///< Absolute error tolerance
    float rel_tolerance = 1e-6f;   ///< Relative error tolerance

    // Step size limits
    float min_step = 1e-6f;        ///< Minimum step size
    float max_step = 0.1f;         ///< Maximum step size
    float initial_step = 0.01f;    ///< Initial step size

    // Adaptive step control
    float safety_factor = 0.9f;    ///< Safety factor (0.8-0.95)
    float step_grow_max = 2.0f;    ///< Maximum step growth factor
    float step_shrink_min = 0.1f;  ///< Minimum step shrink factor

    // Termination
    int max_steps = 10000;         ///< Maximum integration steps
    float max_distance = 1000.0f;  ///< Maximum affine parameter
    float horizon_buffer = 0.01f;  ///< Buffer distance from horizon

    // Method selection
    bool use_adaptive = true;      ///< Use adaptive stepping
    bool preserve_null = true;     ///< Re-normalise to null condition
};

// =============================================================================
// Integration State
// =============================================================================

/// @brief State of geodesic integration
struct GeodesicState {
    Vec4 position;                 ///< 4-position x^μ
    Vec4 velocity;                 ///< 4-velocity k^μ = dx^μ/dλ
    Vec4 covariant_momentum;       ///< Covariant momentum p_μ = g_μν k^ν
    float affine_param = 0.0f;     ///< Affine parameter λ
    float step_size = 0.01f;       ///< Current step size
    float error_estimate = 0.0f;   ///< Local truncation error estimate
};

// =============================================================================
// Integration Result
// =============================================================================

/// @brief Result of geodesic integration
struct IntegrationResult {
    enum class Status {
        Active,        ///< Still integrating
        Escaped,       ///< Reached r_max (escaped to infinity)
        Horizon,       ///< Crossed event horizon
        MaxSteps,      ///< Reached maximum steps
        Error,         ///< Numerical error (NaN, etc.)
        DiskHit        ///< Hit accretion disk
    };

    GeodesicState final_state;
    Status status = Status::Active;
    int steps_taken = 0;
    float total_distance = 0.0f;
    float final_null_error = 0.0f;  ///< Final |g_μν k^μ k^ν|
};

// =============================================================================
// IGeodesicIntegrator - Geodesic Integration Interface
// =============================================================================

class IGeodesicIntegrator {
public:
    virtual ~IGeodesicIntegrator() = default;

    // =========================================================================
    // Identification
    // =========================================================================

    /// @brief Get integrator name (e.g., "RK45", "Verlet")
    virtual const char* name() const = 0;

    /// @brief Get integrator order (accuracy order)
    virtual int order() const = 0;

    /// @brief Is this an adaptive method?
    virtual bool isAdaptive() const = 0;

    /// @brief Is this a symplectic method?
    virtual bool isSymplectic() const = 0;

    // =========================================================================
    // Configuration
    // =========================================================================

    /// @brief Get current configuration
    virtual IntegratorConfig getConfig() const = 0;

    /// @brief Set configuration
    virtual void setConfig(const IntegratorConfig& config) = 0;

    /// @brief Set metric for integration
    virtual void setMetric(IMetric* metric) = 0;

    // =========================================================================
    // Integration
    // =========================================================================

    /// @brief Integrate single step
    /// @param state [in/out] Current state (updated in place)
    /// @return true if step succeeded
    /// PRECONDITION: state.position in valid domain
    /// POSTCONDITION: state.velocity satisfies null condition (if preserve_null)
    virtual bool step(GeodesicState& state) = 0;

    /// @brief Integrate until termination
    /// @param initial Initial state
    /// @return Final state and termination reason
    virtual IntegrationResult integrate(const GeodesicState& initial) = 0;

    // =========================================================================
    // Diagnostics
    // =========================================================================

    /// @brief Compute null condition error |g_μν k^μ k^ν|
    virtual float nullConditionError(const GeodesicState& state) const = 0;

    /// @brief Compute Killing energy E = -g_tμ k^μ
    virtual float killingEnergy(const GeodesicState& state) const = 0;

    /// @brief Compute Killing angular momentum L
    virtual float killingAngularMomentum(const GeodesicState& state) const = 0;
};

// =============================================================================
// Factory Function (to be implemented by concrete integrators)
// =============================================================================

/// @brief Create default integrator (RK45)
/// @param metric Metric to use
/// @param config Optional configuration
/// @return Unique pointer to integrator
// std::unique_ptr<IGeodesicIntegrator> createDefaultIntegrator(
//     IMetric* metric,
//     const IntegratorConfig& config = IntegratorConfig{});

} // namespace Sirius

#endif // PHGD000A_H
