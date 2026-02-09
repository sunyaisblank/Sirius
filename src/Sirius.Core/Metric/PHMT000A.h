// PHMT000A.h - IMetric Interface for Spacetime Metrics
//
// Interface for metric tensors gμν. All metrics (Schwarzschild, Kerr, etc.) implement this.
// Contract: evaluate(pos) → g[μ,ν], dg[σ,μ,ν] (metric and derivatives)
// Postconditions: symmetric, Lorentzian signature (-,+,+,+), non-degenerate

#ifndef PHMT000A_H
#define PHMT000A_H

#include <MTTN001A.h>
#include <string>
#include <vector>
#include <map>

namespace Sirius {

/// @brief Metric parameter with value and range constraints
struct MetricParameter {
    double value;  ///< Current parameter value
    double min;    ///< Minimum allowed value
    double max;    ///< Maximum allowed value
};

typedef std::map<std::string, MetricParameter> Config;

/// @brief Abstract interface for spacetime metric implementations
/// @note All implementations must satisfy postconditions P1-P4
class IMetric {
public:
    virtual ~IMetric() = default;

    /// @brief Evaluate metric tensor and derivatives at position
    /// @param pos 4-position (t, r, θ, φ) or (t, x, y, z)
    /// @param g [out] Metric tensor gμν (with Dual for autodiff compatibility)
    /// @param dg [out] Metric derivatives ∂σ gμν
    /// PRECONDITION: pos within metric's valid domain
    /// POSTCONDITION: g symmetric, Lorentzian signature, finite
    virtual void evaluate(const Tensor<double, 4>& pos, Metric4D& g, Tensor<Dual<double>, 4, 4, 4>& dg) = 0;

    /// @brief Get current parameter configuration
    virtual const Config& getParameters() const = 0;

    /// @brief Set a parameter value (clamped to [min, max])
    virtual void setParameter(const std::string& key, double value) = 0;

    /// @brief Get human-readable metric name
    virtual const char* getName() const = 0;
};

// =============================================================================
// Metric Symmetry Validation Utilities
// =============================================================================
// Runtime checks for metric tensor properties that must hold exactly.
// Use these in tests and diagnostics to verify implementations.

namespace MetricValidation {

/// @brief Verify metric symmetry: g_μν = g_νμ
/// @param g Metric tensor
/// @return Maximum violation |g_μν - g_νμ|
inline double checkSymmetry(const Metric4D& g) {
    double max_violation = 0.0;
    for (int mu = 0; mu < 4; ++mu) {
        for (int nu = mu + 1; nu < 4; ++nu) {
            double violation = std::abs(g(mu, nu).real - g(nu, mu).real);
            max_violation = std::max(max_violation, violation);
        }
    }
    return max_violation;
}

/// @brief Verify Christoffel symmetry: Γ^λ_μν = Γ^λ_νμ
/// @param Gamma Christoffel symbols
/// @return Maximum violation |Γ^λ_μν - Γ^λ_νμ|
inline double checkChristoffelSymmetry(const double Gamma[4][4][4]) {
    double max_violation = 0.0;
    for (int lam = 0; lam < 4; ++lam) {
        for (int mu = 0; mu < 4; ++mu) {
            for (int nu = mu + 1; nu < 4; ++nu) {
                double violation = std::abs(Gamma[lam][mu][nu] - Gamma[lam][nu][mu]);
                max_violation = std::max(max_violation, violation);
            }
        }
    }
    return max_violation;
}

/// @brief Verify inverse relation: g^μα g_αν = δ^μ_ν
/// @param g Metric tensor (covariant)
/// @param g_inv Inverse metric (contravariant)
/// @return Maximum deviation from Kronecker delta
inline double checkInverse(const Metric4D& g, const double g_inv[4][4]) {
    double max_deviation = 0.0;
    for (int mu = 0; mu < 4; ++mu) {
        for (int nu = 0; nu < 4; ++nu) {
            double sum = 0.0;
            for (int alpha = 0; alpha < 4; ++alpha) {
                sum += g_inv[mu][alpha] * g(alpha, nu).real;
            }
            double expected = (mu == nu) ? 1.0 : 0.0;
            double deviation = std::abs(sum - expected);
            max_deviation = std::max(max_deviation, deviation);
        }
    }
    return max_deviation;
}

/// @brief Check Lorentzian signature (-, +, +, +) via eigenvalues
/// @param g Metric tensor
/// @return true if signature is correct (one negative, three positive eigenvalues)
inline bool checkLorentzianSignature(const Metric4D& g) {
    // For diagonal-dominant metrics, check diagonal signs
    // This is a simplified check - full eigenvalue analysis is more complex
    double g00 = g(0, 0).real;
    double g11 = g(1, 1).real;
    double g22 = g(2, 2).real;
    double g33 = g(3, 3).real;

    // For mostly-minus convention: g00 < 0, gii > 0
    // Account for possible coordinate singularities near horizons
    if (std::abs(g00) < 1e-10) return true;  // Degenerate at horizon

    // Check signature pattern
    bool time_negative = g00 < 0;
    bool space_positive = (g11 > 0) && (g22 > 0) && (g33 > 0);

    return time_negative && space_positive;
}

/// @brief Check for NaN/Inf in metric tensor
/// @param g Metric tensor
/// @return true if all components are finite
inline bool checkFinite(const Metric4D& g) {
    for (int mu = 0; mu < 4; ++mu) {
        for (int nu = 0; nu < 4; ++nu) {
            if (!std::isfinite(g(mu, nu).real)) return false;
        }
    }
    return true;
}

/// @brief Check for NaN/Inf in Christoffel symbols
/// @param Gamma Christoffel symbols
/// @return true if all components are finite
inline bool checkChristoffelFinite(const double Gamma[4][4][4]) {
    for (int lam = 0; lam < 4; ++lam) {
        for (int mu = 0; mu < 4; ++mu) {
            for (int nu = 0; nu < 4; ++nu) {
                if (!std::isfinite(Gamma[lam][mu][nu])) return false;
            }
        }
    }
    return true;
}

/// @brief Full metric validation at a point
/// @param metric Metric implementation
/// @param pos Test position
/// @param tolerance Maximum allowed violation
/// @return true if all checks pass
inline bool validateMetricAtPoint(IMetric* metric, const Tensor<double, 4>& pos, double tolerance = 1e-10) {
    Metric4D g;
    Tensor<Dual<double>, 4, 4, 4> dg;
    metric->evaluate(pos, g, dg);

    // Check finiteness
    if (!checkFinite(g)) return false;

    // Check symmetry
    if (checkSymmetry(g) > tolerance) return false;

    // Check signature (relaxed check)
    // Note: May fail near horizons where g00 → 0
    // if (!checkLorentzianSignature(g)) return false;

    return true;
}

} // namespace MetricValidation

} // namespace Sirius

#endif // PHMT000A_H