// Time-transformed explicit symplectic integrator (Yoshida composition of
// Stormer-Verlet leapfrog steps) for double-precision Kerr geodesics: the
// oracle's energy and angular-momentum conservation reference for the live
// RK45 path. Off the render path.
//   H = (1/2) g^uv p_u p_v = 0 (null),  time transform d tau = Sigma d lambda
// Reference: Yoshida (1990); Preto and Saha (1999); James et al. (2015) A.4.
// Ported from PHSI001A.h; numeric content bit-identical. The single
// core-constants reference (ANGULAR_CLAMP_EPS) is inlined as its literal
// 1e-6 to keep the oracle self-contained; see core/constants.h
// (kAngularClampEps).

#pragma once

#include "sirius/base/contracts.h"
#include "sirius/oracle/metric_interface.h"

#include <algorithm>
#include <cmath>
#include <limits>

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
        IntegratorOrder order = IntegratorOrder::kYoshida6;  // Default: 6th-order

        // Null condition renormalization settings
        bool enforceNullCondition = true;  // Renormalize to g_μν k^μ k^ν = 0
        double nullConditionTol = 1e-6;    // Tolerance for null condition drift
        int renormalizeEveryNSteps = 10;   // Frequency of renormalization (0 = every Step)

        Config() = default;  // Explicit default constructor
    };

    struct StepResult {
        HamiltonianStateD state;
        double lambdaAdvance;
        double hamiltonianError;
        double nullConditionError;  // |g_μν k^μ k^ν|
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
        result.nullConditionError = 0;

        HamiltonianStateD s = state;

        // Select coefficients based on order
        const double* weights = nullptr;
        int numSubsteps = 0;

        switch (config_.order) {
            case IntegratorOrder::kYoshida4:
                weights = symplectic::kOrder4Weights;
                numSubsteps = symplectic::kOrder4Substeps;
                break;
            case IntegratorOrder::kYoshida6:
                weights = symplectic::kOrder6Weights;
                numSubsteps = symplectic::kOrder6Substeps;
                break;
            case IntegratorOrder::kYoshida8:
                weights = symplectic::kOrder8Weights;
                numSubsteps = symplectic::kOrder8Substeps;
                break;
        }

        // Apply Yoshida composition: each weight is a leapfrog Step with that fraction of h
        double lambdaAccum = 0;
        for (int i = 0; i < numSubsteps; ++i) {
            double wi = weights[i];
            s = LeapfrogStep(s, wi * h);

            if (!metric_->IsValid(s.q)) {
                result.terminated = true;
                result.state = s;
                result.lambdaAdvance = lambdaAccum;
                // Boyer-Lindquist is outside its chart domain here. Do not
                // evaluate the metric merely to decorate a termination result.
                result.hamiltonianError = std::numeric_limits<double>::infinity();
                result.nullConditionError = std::numeric_limits<double>::infinity();
                return result;
            }

            lambdaAccum += wi * h;
            result.substeps++;
        }

        // Compute null condition error before potential renormalization
        result.nullConditionError = ComputeNullError(s);

        // Apply null condition renormalization if enabled and needed
        if (config_.enforceNullCondition) {
            steps_since_renorm_++;

            bool shouldRenorm = false;
            if (config_.renormalizeEveryNSteps == 0) {
                // Renormalize every Step
                shouldRenorm = true;
            } else if (steps_since_renorm_ >= config_.renormalizeEveryNSteps) {
                shouldRenorm = true;
            } else if (result.nullConditionError > config_.nullConditionTol) {
                // Renormalize if error exceeds tolerance
                shouldRenorm = true;
            }

            if (shouldRenorm) {
                s = RenormalizeNull(s);
                steps_since_renorm_ = 0;
                // Recompute error after renormalization
                result.nullConditionError = ComputeNullError(s);
            }
        }

        // Compute Hamiltonian error
        result.state = s;
        result.state.H = metric_->Hamiltonian(s.q, s.p);
        result.hamiltonianError = std::abs(result.state.H);
        result.lambdaAdvance = h;  // Total Step is always h (weights sum to 1)

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

    IntegrationResult Integrate(const GeodesicStateD& initial, double lambdaMax,
                                double escapeRadius = 1e6) const {
        IntegrationResult result;
        result.totalLambda = 0;
        result.maxHamiltonianError = 0;
        result.totalSteps = 0;
        result.hitHorizon = false;
        result.escaped = false;

        HamiltonianStateD hs(initial);
        double h = config_.initialStepSize;

        while (result.totalLambda < lambdaMax && result.totalSteps < config_.maxStepsPerCall) {
            // Adaptive Step size based on position (smaller near horizon)
            double r = hs.q.r;
            double r_horizon = metric_->HorizonRadius();

            if (r < r_horizon * config_.horizonBuffer) {
                result.hitHorizon = true;
                break;
            }

            if (r > escapeRadius) {
                result.escaped = true;
                break;
            }

            // Adapt Step size: smaller near horizon, larger when far away
            double r_ratio = (r - r_horizon) / r_horizon;
            h = config_.initialStepSize * std::min(1.0, r_ratio * r_ratio);
            h = std::clamp(h, config_.minStepSize, config_.maxStepSize);

            // Don't overshoot lambdaMax
            if (result.totalLambda + h > lambdaMax) {
                h = lambdaMax - result.totalLambda;
            }

            // Take Step
            StepResult sr = Step(hs, h);

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
        result.finalState = hs.ToGeodesicState(result.totalLambda);
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

    ConservationStats TestConservation(const GeodesicStateD& initial, int numSteps,
                                       double stepSize) const {
        ConservationStats stats;
        stats.steps = 0;

        HamiltonianStateD hs(initial);

        // Initial quantities
        stats.initialH = metric_->Hamiltonian(hs.q, hs.p);
        stats.initialE = -hs.p.t;
        stats.initialLz = hs.p.phi;
        stats.maxHError = std::abs(stats.initialH);

        // Integrate
        for (int i = 0; i < numSteps; ++i) {
            StepResult sr = Step(hs, stepSize);
            if (sr.terminated) break;

            hs = sr.state;
            stats.maxHError = std::max(stats.maxHError, sr.hamiltonianError);
            stats.steps++;
        }

        // Final quantities
        stats.finalH = metric_->Hamiltonian(hs.q, hs.p);
        stats.finalE = -hs.p.t;
        stats.finalLz = hs.p.phi;

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
    // Leapfrog (Störmer-Verlet) Substep
    //--------------------------------------------------------------------------

    HamiltonianStateD LeapfrogStep(const HamiltonianStateD& state, double h) const {
        // Leapfrog: p_{n+1/2} = p_n - (h/2) ∂H/∂q
        //           q_{n+1} = q_n + h ∂H/∂p(q_n, p_{n+1/2})
        //           p_{n+1} = p_{n+1/2} - (h/2) ∂H/∂q(q_{n+1})

        HamiltonianStateD s = state;

        // Half-Step momentum update
        Vec4d dHdq = metric_->dHdq(s.q, s.p);
        s.p = s.p + dHdq * (h * 0.5);  // Note: dHdq returns dp/dλ = -∂H/∂q (negated)

        // Full-Step position update
        Vec4d dHdp = metric_->dHdp(s.q, s.p);
        s.q = s.q + dHdp * h;

        // Normalise angles
        s.q.phi = std::fmod(s.q.phi, 2 * M_PI);
        if (s.q.phi < 0) s.q.phi += 2 * M_PI;
        // 1e-6 == core/constants.h kAngularClampEps; inlined to keep the
        // oracle self-contained (STYLE.md self-containment over shared const).
        s.q.theta = std::clamp(s.q.theta, 1e-6, M_PI - 1e-6);

        // A composed substep can cross the Boyer-Lindquist horizon even when
        // the enclosing step began outside it. Return the position for the
        // caller's termination check before evaluating dH/dq outside the
        // metric chart.
        if (!metric_->IsValid(s.q)) {
            return s;
        }

        // Half-Step momentum update (at new position)
        dHdq = metric_->dHdq(s.q, s.p);
        s.p = s.p + dHdq * (h * 0.5);

        return s;
    }
};

}  // namespace sirius::oracle
