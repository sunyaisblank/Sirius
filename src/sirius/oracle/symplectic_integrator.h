// Canonical symplectic integrator for the double-precision Kerr Hamiltonian:
// Yoshida composition of symmetric implicit-midpoint maps. Implicit midpoint
// is required because H(q,p) = (1/2) g^uv(q) p_u p_v is non-separable; an
// explicit kick-drift-kick map is not symplectic for this Hamiltonian.
// Optional state-dependent step selection and null projection are operational
// stabilisers outside the strict fixed-step symplectic claim. Off the render
// path.
// Reference: Yoshida (1990); implicit midpoint as a symplectic Runge-Kutta map.
// The single
// core-constants reference (ANGULAR_CLAMP_EPS) is inlined as its literal
// 1e-6 to keep the oracle self-contained; see core/constants.h
// (kAngularClampEps).

#pragma once

#include "sirius/base/contracts.h"
#include "sirius/oracle/metric_interface.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>

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
    struct Config {
        double initial_step_size = 0.1;
        double min_step_size = 1e-8;
        double max_step_size = 1.0;
        double tolerance = 1e-10;
        int max_steps_per_call = 10000;
        double horizon_buffer = 1.01;
        IntegratorOrder order = IntegratorOrder::kYoshida6;  // Default: 6th-order

        // Null condition renormalization settings
        bool enforce_null_condition = true;      // Renormalize to g_μν k^μ k^ν = 0
        double null_condition_tolerance = 1e-6;  // Tolerance for null condition drift
        int renormalize_every_n_steps = 10;      // Frequency of renormalization (0 = every Step)

        Config() = default;  // Explicit default constructor
    };

    struct StepResult {
        HamiltonianStateD state;
        double lambda_advance;
        double hamiltonian_error;
        double null_condition_error;  // |g_μν k^μ k^ν|
        bool terminated;
        int substeps;
    };

    explicit SymplecticIntegratorD(const IMetricD* metric) : metric_(metric), config_() {
        SIRIUS_PRE(metric != nullptr);
    }

    SymplecticIntegratorD(const IMetricD* metric, const Config& config)
        : metric_(metric), config_(config) {
        SIRIUS_PRE(metric != nullptr);
    }

    //--------------------------------------------------------------------------
    // Single Integration Step (Selectable Order)
    //--------------------------------------------------------------------------

    StepResult Step(const HamiltonianStateD& state, double h) const {
        StepResult result;
        result.substeps = 0;
        result.terminated = false;
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
            s = ImplicitMidpointStep(s, wi * h);

            if (!metric_->IsValid(s.q)) {
                result.terminated = true;
                result.state = s;
                result.lambda_advance = lambda_accumulator;
                // Boyer-Lindquist is outside its chart domain here. Do not
                // evaluate the metric merely to decorate a termination result.
                result.hamiltonian_error = std::numeric_limits<double>::infinity();
                result.null_condition_error = std::numeric_limits<double>::infinity();
                return result;
            }

            lambda_accumulator += wi * h;
            result.substeps++;
        }

        // Compute null condition error before potential renormalization
        result.null_condition_error = ComputeNullError(s);

        // Apply null condition renormalization if enabled and needed
        if (config_.enforce_null_condition) {
            steps_since_renorm_++;

            bool should_renormalize = false;
            if (config_.renormalize_every_n_steps == 0) {
                // Renormalize every Step
                should_renormalize = true;
            } else if (steps_since_renorm_ >= config_.renormalize_every_n_steps) {
                should_renormalize = true;
            } else if (result.null_condition_error > config_.null_condition_tolerance) {
                // Renormalize if error exceeds tolerance
                should_renormalize = true;
            }

            if (should_renormalize) {
                s = RenormalizeNull(s);
                steps_since_renorm_ = 0;
                // Recompute error after renormalization
                result.null_condition_error = ComputeNullError(s);
            }
        }

        // Compute Hamiltonian error
        result.state = s;
        result.state.H = metric_->Hamiltonian(s.q, s.p);
        result.hamiltonian_error = std::abs(result.state.H);
        result.lambda_advance = h;  // Total Step is always h (weights sum to 1)

        return result;
    }

    //--------------------------------------------------------------------------
    // Integrate to Target Lambda or Termination
    //--------------------------------------------------------------------------

    struct IntegrationResult {
        GeodesicStateD final_state;
        double total_lambda;
        double max_hamiltonian_error;
        int total_steps;
        bool hit_horizon;
        bool escaped;
    };

    IntegrationResult Integrate(const GeodesicStateD& initial, double lambdaMax,
                                double escape_radius = 1e6) const {
        IntegrationResult result;
        result.total_lambda = 0;
        result.max_hamiltonian_error = 0;
        result.total_steps = 0;
        result.hit_horizon = false;
        result.escaped = false;

        HamiltonianStateD hs(initial);
        double h = config_.initial_step_size;

        while (result.total_lambda < lambdaMax && result.total_steps < config_.max_steps_per_call) {
            // Adaptive Step size based on position (smaller near horizon)
            double r = hs.q.r;
            double r_horizon = metric_->HorizonRadius();

            if (r < r_horizon * config_.horizon_buffer) {
                result.hit_horizon = true;
                break;
            }

            if (r > escape_radius) {
                result.escaped = true;
                break;
            }

            // Adapt Step size: smaller near horizon, larger when far away
            double r_ratio = (r - r_horizon) / r_horizon;
            h = config_.initial_step_size * std::min(1.0, r_ratio * r_ratio);
            h = std::clamp(h, config_.min_step_size, config_.max_step_size);

            // Don't overshoot lambdaMax
            if (result.total_lambda + h > lambdaMax) {
                h = lambdaMax - result.total_lambda;
            }

            // Take Step
            StepResult sr = Step(hs, h);

            if (sr.terminated) {
                result.hit_horizon = true;
                break;
            }

            hs = sr.state;
            result.total_lambda += sr.lambda_advance;
            result.max_hamiltonian_error =
                std::max(result.max_hamiltonian_error, sr.hamiltonian_error);
            result.total_steps++;
        }

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
            StepResult sr = Step(hs, step_size);
            if (sr.terminated) break;

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
    mutable int steps_since_renorm_ = 0;  // Track steps since last renormalization

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

    HamiltonianStateD RenormalizeNull(const HamiltonianStateD& state) const {
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

        if (disc < 0 || std::abs(A) < 1e-15) {
            // No real solution or degenerate - return unchanged
            // This can happen near coordinate singularities
            return s;
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

    HamiltonianStateD ImplicitMidpointStep(const HamiltonianStateD& state, double h) const {
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
            if (!metric_->IsValid(midpoint.q)) return next;

            HamiltonianStateD candidate = state;
            candidate.q = state.q + metric_->dHdp(midpoint.q, midpoint.p) * h;
            candidate.p = state.p + metric_->dHdq(midpoint.q, midpoint.p) * h;

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
        SIRIUS_ASSERT(converged);

        next.q.phi = std::fmod(next.q.phi, 2 * std::numbers::pi);
        if (next.q.phi < 0) next.q.phi += 2 * std::numbers::pi;
        // 1e-6 == core/constants.h kAngularClampEps; inlined to keep the
        // oracle self-contained (STYLE.md self-containment over shared const).
        next.q.theta = std::clamp(next.q.theta, 1e-6, std::numbers::pi - 1e-6);
        return next;
    }
};

}  // namespace sirius::oracle
