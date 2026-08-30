#pragma once

// Null geodesic integration for relativistic ray tracing: the CPU reference
// tracer's Dormand-Prince RK45 Hamiltonian integrator. Ported from PHGD001A.h.
//
// Integrates d^2 x^mu/dlambda^2 + Gamma^mu_alpha_beta (dx^alpha/dlambda)
// (dx^beta/dlambda) = 0 under the null constraint g_mu_nu k^mu k^nu = 0, with
// adaptive step control and termination on horizon capture, escape, or NaN/Inf.

#include "sirius/core/metrics/metric.h"
#include "sirius/core/tensor.h"

#include <cmath>

namespace sirius::core {

// State vector for a null geodesic.
struct Lightray {
    Vec4 position;          // 4-position (t, r, theta, phi) or (t, x, y, z).
    Vec4 velocity;          // 4-velocity k^mu = dx^mu/dlambda.
    Vec4 acceleration;      // 4-acceleration d^2 x^mu/dlambda^2.
    float ku_uobsu;         // Positive launch-observer frequency -k.u_obs.
    float proper_time;      // Accumulated affine parameter lambda.
    float coordinate_time;  // Accumulated coordinate time t.
    int terminated;         // Termination status (0 = active).
    int sx, sy;             // Screen pixel coordinates.
    float step_size;        // Current adaptive step size.
    int padding;            // Alignment padding.
};

// Local observer reference frame.
struct ObserverState {
    Vec4 position;        // Observer 4-position.
    Vec4 velocity;        // Observer 4-velocity, normalised g(u, u) = -1.
    Vec4 e0, e1, e2, e3;  // Orthonormal tetrad (e0 = u, e1-e3 spatial).
    bool is_timelike;     // True when velocity is timelike.
};

// Configuration for adaptive integration.
struct IntegratorConfig {
    float abs_tolerance = 1e-6f;
    float rel_tolerance = 1e-6f;
    float min_step = 1e-6f;
    float max_step = 0.1f;
    float initial_step = 0.01f;
    float safety_factor = 0.9f;  // Step adaptation safety factor (0.8-0.95).
    float step_grow_max = 2.0f;
    float step_shrink_min = 0.1f;
};

// Parameters consumed by every RK45 step. initial_step is deliberately absent:
// the caller consumes it once when constructing the ray's current step.
[[nodiscard]] inline bool IsRepresentedIntegratorStepControl(
    const IntegratorConfig& config) noexcept {
    return std::isfinite(config.abs_tolerance) && config.abs_tolerance > 0.0f &&
           std::isfinite(config.rel_tolerance) && config.rel_tolerance > 0.0f &&
           std::isfinite(config.min_step) && config.min_step > 0.0f &&
           std::isfinite(config.max_step) && config.max_step >= config.min_step &&
           std::isfinite(config.safety_factor) && config.safety_factor >= 0.8f &&
           config.safety_factor <= 0.95f && std::isfinite(config.step_grow_max) &&
           config.step_grow_max >= 1.0f && std::isfinite(config.step_shrink_min) &&
           config.step_shrink_min > 0.0f && config.step_shrink_min <= 1.0f &&
           config.step_shrink_min <= config.step_grow_max;
}

[[nodiscard]] inline bool IsRepresentedIntegratorConfig(const IntegratorConfig& config) noexcept {
    return IsRepresentedIntegratorStepControl(config) && std::isfinite(config.initial_step) &&
           config.initial_step >= config.min_step && config.initial_step <= config.max_step;
}

// Internal state for the RK45 integrator.
struct Rk45State {
    Vec4 x;         // Position.
    Vec4 p;         // Covariant momentum.
    Vec4 k;         // Contravariant velocity (derived from p).
    float error;    // Estimated local truncation error.
    bool accepted;  // Whether the step was accepted.
};

// Static methods for geodesic integration.
class Geodesic {
  public:
    // Attempt one adaptive Dormand-Prince RK45 step. A false result with
    // terminated == 0 is a
    // recoverable rejection: step_size is reduced and the state is unchanged.
    static bool IntegrateStep(Lightray& ray, IMetric* metric, float min_step = 1e-6f,
                              float max_step = 0.1f);

    // Geodesic acceleration a^mu = -Gamma^mu_alpha_beta v^alpha v^beta.
    static Vec4 CalculateAcceleration(const Vec4& velocity, const Vec4& position, IMetric* metric);

    // Whether the ray should terminate.
    static bool CheckTermination(const Lightray& ray, IMetric* metric);

    // Static-observer frequency shift z = nu_emit/nu_obs - 1 where both static
    // worldlines are timelike. Live disk transfer uses invariant emitter and
    // observer contractions instead.
    static float CalculateRedshift(const Lightray& ray, const ObserverState& observer,
                                   IMetric* metric);

    // Create and initialise an observer. The supplied coordinate velocity must
    // be finite and timelike; invalid worldlines contract-fail rather than
    // silently becoming a different static observer. The velocity is normalised to
    // g(u, u) = -1.
    static ObserverState CreateObserver(const Vec4& position, const Vec4& velocity,
                                        IMetric* metric);

    // Construct an orthonormal tetrad for the observer via Gram-Schmidt, so
    // g(e_a, e_b) = eta_a_b.
    static void CalculateTetrads(ObserverState& observer, IMetric* metric);

    // Integrate one step with the RK45 (Dormand-Prince) embedded 4th/5th order
    // pair; the step is adapted from ||y5 - y4||. A false result with
    // terminated == 0 is a recoverable rejection and must be retried.
    static bool IntegrateStepRk45(Lightray& ray, IMetric* metric, const IntegratorConfig& config);

    // Optimal step from the error estimate: h_new = h safety (tol/err)^(1/5).
    static float ComputeOptimalStep(float h, float error, float tolerance,
                                    const IntegratorConfig& config);

    // Default integrator configuration.
    static IntegratorConfig GetDefaultConfig();
};

}  // namespace sirius::core
