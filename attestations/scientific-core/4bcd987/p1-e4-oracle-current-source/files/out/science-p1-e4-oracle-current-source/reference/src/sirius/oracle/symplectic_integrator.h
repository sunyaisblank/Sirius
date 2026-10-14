// Canonical symplectic integrator for the double-precision Kerr Hamiltonian:
// Yoshida composition of symmetric implicit-midpoint maps. Implicit midpoint
// is required because H(q,p) = (1/2) g^uv(q) p_u p_v is non-separable; an
// explicit kick-drift-kick map is not symplectic for this Hamiltonian.
// Optional state-dependent step selection and null projection are operational
// stabilisers outside the strict fixed-step symplectic claim. Off the render
// path.
// Reference: Yoshida (1990); implicit midpoint as a symplectic Runge-Kutta map.

#pragma once

#include "sirius/base/contracts.h"
#include "sirius/oracle/metric_interface.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <optional>

namespace sirius::oracle {

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
constexpr double kW41 = 1.3512071919596576;
constexpr double kW40 = -1.7024143839193153;
constexpr int kOrder4Substeps = 3;
constexpr double kOrder4Weights[] = {kW41, kW40, kW41};

// Yoshida 6th-order "Solution A" coefficients (7 substeps)
// From Yoshida (1990), optimised for minimal error
constexpr int kOrder6Substeps = 7;
constexpr double kOrder6Weights[] = {0.78451361047756, 0.23557321335936,  -1.17767998417887,
                                     1.31518632068391, -1.17767998417887, 0.23557321335936,
                                     0.78451361047756};

// Yoshida 8th-order coefficients (15 substeps)
// From Yoshida (1990), higher precision
constexpr int kOrder8Substeps = 15;
constexpr double kOrder8Weights[] = {
    0.74167036435061, -0.40910082580003, 0.19075471029623,  -0.57386247111608, 0.29906418130365,
    0.33462491824529, 0.31529309239676,  -0.79688793935291, 0.31529309239676,  0.33462491824529,
    0.29906418130365, -0.57386247111608, 0.19075471029623,  -0.40910082580003, 0.74167036435061};

}  // namespace symplectic

//==============================================================================
// Integrator Order Selection
//==============================================================================

enum class IntegratorOrder {
    kYoshida4,  // 4th-order, 3 substeps, fast
    kYoshida6,  // 6th-order, 7 substeps, balanced (DEFAULT)
    kYoshida8   // 8th-order, 15 substeps, high precision
};

//==============================================================================
// SymplecticIntegratorD: multi-order fixed-step canonical composition
//==============================================================================

class SymplecticIntegratorD {
  public:
    enum class StepOutcome { kAccepted, kChartDomainExit, kConstraintFailure };

    struct Config {
        double initial_step_size = 0.1;
        double min_step_size = 1e-8;
        double max_step_size = 1.0;
        int max_steps_per_call = 10000;
        double horizon_buffer = 1.01;
        IntegratorOrder order = IntegratorOrder::kYoshida6;  // Default: 6th-order

        // Null condition renormalization settings
        bool enforce_null_condition = true;      // Renormalize to g_μν k^μ k^ν = 0
        double null_condition_tolerance = 1e-6;  // Tolerance for null condition drift
        int renormalize_every_n_steps = 10;      // Frequency of renormalization (0 = every Step)

        Config() = default;  // Explicit default constructor
    };

    [[nodiscard]] static bool IsRepresentedConfig(const Config& config) noexcept {
        const bool known_order = config.order == IntegratorOrder::kYoshida4 ||
                                 config.order == IntegratorOrder::kYoshida6 ||
                                 config.order == IntegratorOrder::kYoshida8;
        return std::isfinite(config.initial_step_size) && std::isfinite(config.min_step_size) &&
               config.min_step_size > 0.0 && std::isfinite(config.max_step_size) &&
               config.max_step_size >= config.min_step_size &&
               config.initial_step_size >= config.min_step_size &&
               config.initial_step_size <= config.max_step_size && config.max_steps_per_call > 0 &&
               std::isfinite(config.horizon_buffer) && config.horizon_buffer > 1.0 && known_order &&
               std::isfinite(config.null_condition_tolerance) &&
               config.null_condition_tolerance > 0.0 && config.renormalize_every_n_steps >= 0;
    }

    struct StepResult {
        HamiltonianStateD state;
        double lambda_advance;
        double hamiltonian_error;
        double null_condition_error;  // |g_μν k^μ k^ν|
        StepOutcome outcome;
        int substeps;
    };

    explicit SymplecticIntegratorD(const IMetricD* metric)
        : SymplecticIntegratorD(metric, Config()) {}

    SymplecticIntegratorD(const IMetricD* metric, const Config& config)
        : metric_(metric), config_(config) {
        SIRIUS_PRE(metric != nullptr);
        SIRIUS_PRE(IsRepresentedConfig(config));
    }

    //--------------------------------------------------------------------------
    // Single Integration Step (Selectable Order)
    //--------------------------------------------------------------------------

    // accepted_steps_before makes the optional periodic null projection a
    // property of the represented trajectory, not mutable integrator history.
    StepResult Step(const HamiltonianStateD& state, double h, int accepted_steps_before = 0) const {
        SIRIUS_PRE(std::isfinite(h) && h != 0.0);
        SIRIUS_PRE(accepted_steps_before >= 0);
        const auto failure = [](const HamiltonianStateD& failed_state, StepOutcome outcome,
                                double lambda_advance, int substeps) {
            const double infinity = std::numeric_limits<double>::infinity();
            return StepResult{failed_state, lambda_advance, infinity, infinity, outcome, substeps};
        };
        if (!IsFiniteVector(state.q) || !metric_->IsValid(state.q)) {
            return failure(state, StepOutcome::kChartDomainExit, 0.0, 0);
        }
        if (!IsFiniteVector(state.p)) {
            return failure(state, StepOutcome::kConstraintFailure, 0.0, 0);
        }
        StepResult result;
        result.substeps = 0;
        result.outcome = StepOutcome::kAccepted;
        result.null_condition_error = 0;

        HamiltonianStateD s = state;

        // Select coefficients based on order
        const double* weights = nullptr;
        int substep_count = 0;

        switch (config_.order) {
            case IntegratorOrder::kYoshida4:
                weights = symplectic::kOrder4Weights;
                substep_count = symplectic::kOrder4Substeps;
                break;
            case IntegratorOrder::kYoshida6:
                weights = symplectic::kOrder6Weights;
                substep_count = symplectic::kOrder6Substeps;
                break;
            case IntegratorOrder::kYoshida8:
                weights = symplectic::kOrder8Weights;
                substep_count = symplectic::kOrder8Substeps;
                break;
        }

        // Apply Yoshida composition: each weight is one symmetric second-order
        // implicit-midpoint map with that fraction of h.
        double lambda_accumulator = 0;
        for (int i = 0; i < substep_count; ++i) {
            double wi = weights[i];
            const MidpointResult midpoint = TryImplicitMidpointStep(s, wi * h);
            if (midpoint.outcome != StepOutcome::kAccepted) {
                return failure(midpoint.state, midpoint.outcome, lambda_accumulator,
                               result.substeps);
            }
            s = midpoint.state;

            if (!metric_->IsValid(s.q) || !IsFiniteVector(s.p)) {
                // Boyer-Lindquist is outside its chart domain here. Do not
                // evaluate the metric merely to decorate a termination result.
                return failure(s, StepOutcome::kChartDomainExit, lambda_accumulator,
                               result.substeps);
            }

            lambda_accumulator += wi * h;
            result.substeps++;
        }

        // Compute null condition error before potential renormalization
        result.null_condition_error = ComputeNullError(s);
        if (!std::isfinite(result.null_condition_error)) {
            return failure(s, StepOutcome::kConstraintFailure, lambda_accumulator, result.substeps);
        }

        // Apply null condition renormalization if enabled and needed
        if (config_.enforce_null_condition) {
            const bool scheduled = config_.renormalize_every_n_steps == 0 ||
                                   accepted_steps_before % config_.renormalize_every_n_steps ==
                                       config_.renormalize_every_n_steps - 1;
            const bool should_renormalize =
                scheduled || result.null_condition_error > config_.null_condition_tolerance;

            if (should_renormalize) {
                const std::optional<HamiltonianStateD> renormalized = TryRenormalizeNull(s);
                if (!renormalized) {
                    return failure(s, StepOutcome::kConstraintFailure, lambda_accumulator,
                                   result.substeps);
                }
                s = *renormalized;
                // Recompute error after renormalization
                result.null_condition_error = ComputeNullError(s);
                if (!std::isfinite(result.null_condition_error) ||
                    result.null_condition_error > config_.null_condition_tolerance) {
                    return failure(s, StepOutcome::kConstraintFailure, lambda_accumulator,
                                   result.substeps);
                }
            }
        }

        // Compute Hamiltonian error
        result.state = s;
        result.state.H = metric_->Hamiltonian(s.q, s.p);
        if (!std::isfinite(result.state.H)) {
            return failure(s, StepOutcome::kConstraintFailure, lambda_accumulator, result.substeps);
        }
        result.hamiltonian_error = std::abs(result.state.H);
        result.lambda_advance = h;  // Total Step is always h (weights sum to 1)

        return result;
    }

    //--------------------------------------------------------------------------
    // Integrate to Target Lambda or Termination
    //--------------------------------------------------------------------------

    enum class IntegrationOutcome {
        kAffineLimit,
        kStepLimit,
        kCaptured,
        kEscaped,
        kChartDomainExit,
        kConstraintFailure
    };

    struct IntegrationResult {
        GeodesicStateD final_state;
        double total_lambda;
        double max_hamiltonian_error;
        int total_steps;
        IntegrationOutcome outcome;
    };

    IntegrationResult Integrate(const GeodesicStateD& initial, double lambdaMax,
                                double escape_radius = 1e6) const {
        SIRIUS_PRE(std::isfinite(lambdaMax) && lambdaMax >= 0.0);
        SIRIUS_PRE(std::isfinite(escape_radius) && escape_radius > 0.0);
        IntegrationResult result;
        result.total_lambda = 0;
        result.max_hamiltonian_error = 0;
        result.total_steps = 0;
        result.outcome = IntegrationOutcome::kStepLimit;

        HamiltonianStateD hs(initial);
        double h = config_.initial_step_size;

        while (result.total_lambda < lambdaMax && result.total_steps < config_.max_steps_per_call) {
            // Adaptive Step size based on position (smaller near horizon)
            double r = hs.q.r;
            double r_horizon = metric_->HorizonRadius();

            if (!IsFiniteVector(hs.q)) {
                result.outcome = IntegrationOutcome::kChartDomainExit;
                break;
            }
            if (!IsFiniteVector(hs.p)) {
                result.outcome = IntegrationOutcome::kConstraintFailure;
                break;
            }
            if (!metric_->IsValid(hs.q)) {
                result.outcome = IntegrationOutcome::kChartDomainExit;
                break;
            }

            if (r <= r_horizon * config_.horizon_buffer) {
                result.outcome = IntegrationOutcome::kCaptured;
                break;
            }

            if (r > escape_radius) {
                result.outcome = IntegrationOutcome::kEscaped;
                break;
            }

            // Adapt Step size: smaller near horizon, larger when far away
            if (r_horizon > 0.0) {
                const double r_ratio = (r - r_horizon) / r_horizon;
                h = config_.initial_step_size * std::min(1.0, r_ratio * r_ratio);
            } else {
                h = config_.initial_step_size;
            }
            h = std::clamp(h, config_.min_step_size, config_.max_step_size);

            // Don't overshoot lambdaMax
            if (result.total_lambda + h > lambdaMax) {
                h = lambdaMax - result.total_lambda;
            }

            // Take Step
            StepResult sr = Step(hs, h, result.total_steps);

            if (sr.outcome != StepOutcome::kAccepted) {
                if (sr.outcome == StepOutcome::kConstraintFailure) {
                    result.outcome = IntegrationOutcome::kConstraintFailure;
                } else {
                    // Capture is admitted only by the represented pre-step
                    // radial event above. An invalid substep is not localised,
                    // so its radius alone cannot prove which chart boundary
                    // was crossed first.
                    result.outcome = IntegrationOutcome::kChartDomainExit;
                }
                break;
            }

            hs = sr.state;
            result.total_lambda += sr.lambda_advance;
            result.max_hamiltonian_error =
                std::max(result.max_hamiltonian_error, sr.hamiltonian_error);
            result.total_steps++;
        }

        if (result.total_lambda >= lambdaMax) result.outcome = IntegrationOutcome::kAffineLimit;

        // Convert back to GeodesicStateD
        result.final_state = hs.ToGeodesicState(result.total_lambda);
        result.final_state.E = -hs.p.t;
        result.final_state.Lz = hs.p.phi;

        return result;
    }

    //--------------------------------------------------------------------------
    // Conservation Tests (for Validation)
    //--------------------------------------------------------------------------

    // Integrate for N steps and return conservation statistics
    struct ConservationStats {
        double initial_h;
        double final_h;
        double max_h_error;
        double initial_e;
        double final_e;
        double initial_lz;
        double final_lz;
        int steps;
    };

    ConservationStats TestConservation(const GeodesicStateD& initial, int numSteps,
                                       double step_size) const {
        ConservationStats stats;
        stats.steps = 0;

        HamiltonianStateD hs(initial);

        // Initial quantities
        stats.initial_h = metric_->Hamiltonian(hs.q, hs.p);
        stats.initial_e = -hs.p.t;
        stats.initial_lz = hs.p.phi;
        stats.max_h_error = std::abs(stats.initial_h);

        // Integrate
        for (int i = 0; i < numSteps; ++i) {
            StepResult sr = Step(hs, step_size, i);
            if (sr.outcome != StepOutcome::kAccepted) break;

            hs = sr.state;
            stats.max_h_error = std::max(stats.max_h_error, sr.hamiltonian_error);
            stats.steps++;
        }

        // Final quantities
        stats.final_h = metric_->Hamiltonian(hs.q, hs.p);
        stats.final_e = -hs.p.t;
        stats.final_lz = hs.p.phi;

        return stats;
    }

  private:
    const IMetricD* metric_;
    Config config_;
    [[nodiscard]] static bool IsFiniteVector(const Vec4d& vector) noexcept {
        return std::isfinite(vector.t) && std::isfinite(vector.r) && std::isfinite(vector.theta) &&
               std::isfinite(vector.phi);
    }

    //--------------------------------------------------------------------------
    // Null Condition Renormalization (Energy-Preserving)
    //--------------------------------------------------------------------------
    //
    // For null geodesics (photons), the constraint g^μν p_μ p_ν = 0 must hold.
    // Numerical integration causes this to drift.
    //
    // IMPORTANT: We must preserve the Killing conserved quantities:
    //   E  = -p_t  (from ∂/∂t Killing vector)
    //   Lz = p_φ   (from ∂/∂φ Killing vector)
    //
    // Therefore we renormalize by solving for p_r (NOT p_t) from null condition:
    //   g^rr p_r² + (other terms with fixed p_t, p_θ, p_φ) = 0
    //
    // This preserves E and Lz while restoring the null condition.
    //
    // Reference: core/constants.h "Null Constraint Preservation"
    //--------------------------------------------------------------------------

    [[nodiscard]] std::optional<HamiltonianStateD> TryRenormalizeNull(
        const HamiltonianStateD& state) const {
        HamiltonianStateD s = state;

        // Get inverse metric at current position
        double g[4][4], g_inv[4][4];
        metric_->Evaluate(s.q, g, g_inv);

        // Extract momenta (p_t and p_φ are CONSERVED and must not change)
        double p_t = s.p.t;  // Conserved: E = -p_t
        double p_th = s.p.theta;
        double p_ph = s.p.phi;  // Conserved: Lz = p_φ

        // Null condition: g^μν p_μ p_ν = 0
        // Solve for p_r:
        //   g^rr p_r² + 2g^rθ p_r p_θ + (terms without p_r) = 0
        //
        // In Boyer-Lindquist/Kerr, g^rθ = 0, simplifying to:
        //   g^rr p_r² = -[g^tt p_t² + 2g^tφ p_t p_φ + g^θθ p_θ² + g^φφ p_φ²]

        // Coefficient for p_r²
        double A = g_inv[1][1];  // g^rr

        // Coefficient for p_r (from cross terms g^rθ, g^rt, g^rφ)
        double B = 2.0 * (g_inv[0][1] * p_t + g_inv[1][2] * p_th + g_inv[1][3] * p_ph);

        // Constant term (all terms without p_r)
        double C = g_inv[0][0] * p_t * p_t + g_inv[2][2] * p_th * p_th + g_inv[3][3] * p_ph * p_ph +
                   2.0 * g_inv[0][2] * p_t * p_th + 2.0 * g_inv[0][3] * p_t * p_ph +
                   2.0 * g_inv[2][3] * p_th * p_ph;

        // Solve quadratic: A*p_r² + B*p_r + C = 0
        double disc = B * B - 4 * A * C;

        if (!std::isfinite(A) || !std::isfinite(B) || !std::isfinite(C) || !std::isfinite(disc) ||
            disc < 0.0 || A == 0.0) {
            return std::nullopt;
        }

        double sqrtDisc = std::sqrt(disc);
        double p_r_plus = (-B + sqrtDisc) / (2.0 * A);
        double p_r_minus = (-B - sqrtDisc) / (2.0 * A);

        // Choose the root closest to current p_r to preserve direction of motion
        double p_r_new;
        if (std::abs(p_r_plus - s.p.r) < std::abs(p_r_minus - s.p.r)) {
            p_r_new = p_r_plus;
        } else {
            p_r_new = p_r_minus;
        }

        if (!std::isfinite(p_r_new)) return std::nullopt;
        s.p.r = p_r_new;
        return s;
    }

    /// @brief Compute null condition error |g^μν p_μ p_ν|
    double ComputeNullError(const HamiltonianStateD& state) const {
        // This is equivalent to 2|H| for the Hamiltonian H = (1/2)g^μν p_μ p_ν
        return std::abs(2.0 * metric_->Hamiltonian(state.q, state.p));
    }

    //--------------------------------------------------------------------------
    // Symmetric implicit-midpoint substep
    //--------------------------------------------------------------------------

    struct MidpointResult {
        HamiltonianStateD state;
        StepOutcome outcome;
    };

    MidpointResult TryImplicitMidpointStep(const HamiltonianStateD& state, double h) const {
        // z_(n+1) = z_n + h J grad H((z_n + z_(n+1))/2). Fixed-point
        // iteration solves the nonlinear midpoint equations to near machine
        // precision. Every iterate is derived from the original state so the
        // converged map remains the implicit Runge-Kutta map, not a chain of
        // explicit Euler updates.
        HamiltonianStateD next = state;
        next.q = state.q + metric_->dHdp(state.q, state.p) * h;
        next.p = state.p + metric_->dHdq(state.q, state.p) * h;

        constexpr int kMaximumIterations = 64;
        constexpr double kRelativeTolerance = 2.0e-14;
        bool converged = false;
        for (int iteration = 0; iteration < kMaximumIterations; ++iteration) {
            HamiltonianStateD midpoint;
            midpoint.q = (state.q + next.q) * 0.5;
            midpoint.p = (state.p + next.p) * 0.5;
            if (!IsFiniteVector(midpoint.q) || !metric_->IsValid(midpoint.q)) {
                return MidpointResult{next, StepOutcome::kChartDomainExit};
            }
            if (!IsFiniteVector(midpoint.p)) {
                return MidpointResult{next, StepOutcome::kConstraintFailure};
            }

            HamiltonianStateD candidate = state;
            candidate.q = state.q + metric_->dHdp(midpoint.q, midpoint.p) * h;
            candidate.p = state.p + metric_->dHdq(midpoint.q, midpoint.p) * h;
            if (!IsFiniteVector(candidate.q) || !IsFiniteVector(candidate.p)) {
                return MidpointResult{candidate, StepOutcome::kConstraintFailure};
            }

            double error = 0.0;
            double scale = 1.0;
            for (int component = 0; component < 4; ++component) {
                error = std::max(error, std::abs(candidate.q[component] - next.q[component]));
                error = std::max(error, std::abs(candidate.p[component] - next.p[component]));
                scale = std::max(scale, std::abs(candidate.q[component]));
                scale = std::max(scale, std::abs(candidate.p[component]));
            }
            next = candidate;
            if (error <= kRelativeTolerance * scale) {
                converged = true;
                break;
            }
        }
        if (!converged) return MidpointResult{next, StepOutcome::kConstraintFailure};

        next.q.phi = std::fmod(next.q.phi, 2 * std::numbers::pi);
        if (next.q.phi < 0) next.q.phi += 2 * std::numbers::pi;
        return MidpointResult{next, StepOutcome::kAccepted};
    }
};

}  // namespace sirius::oracle
