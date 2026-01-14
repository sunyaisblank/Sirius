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

#endif // PHMT000A_H