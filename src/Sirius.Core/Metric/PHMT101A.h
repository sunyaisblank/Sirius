// =============================================================================
// PHMT101A.h - Unified Morris-Thorne Wormhole Family
// Component ID: PHMT101A (Physics/Metric/Morris-Thorne Family)
// =============================================================================
//
// UNIFIED WORMHOLE FAMILY
// =======================
// This implementation covers traversable wormholes through parameterization:
//
// | b(r)        | Φ(r) | Spacetime           |
// |-------------|------|---------------------|
// | b₀²/r       | 0    | Ellis drainhole     |
// | b₀          | 0    | Zero-tidal wormhole |
// | b₀(2-b₀/r)  | 0    | Absurdly benign     |
//
// Mathematical Foundation
// =======================
// Morris-Thorne metric in Schwarzschild-like coordinates (t, r, θ, φ):
//
//   ds² = -e^{2Φ(r)} dt² + dr²/(1 - b(r)/r) + r²(dθ² + sin²θ dφ²)
//
// where:
//   Φ(r) = redshift function (determines gravitational time dilation)
//   b(r) = shape function (determines embedding geometry)
//   b₀   = throat radius (minimum circumference location)
//
// For traversability: b(b₀) = b₀, Φ finite everywhere
//
// Cartesian Embedding
// ===================
// To avoid coordinate singularity at throat, we use proper radial distance:
//   l = ±∫[b₀ to r] dr / √(1 - b(r)/r)
//
// In Cartesian coords (t, x, y, z):
//   r(l) = obtained by inverting l(r)
//   The metric becomes regular at the throat.
//
// Reference: Morris & Thorne, Am. J. Phys. 56, 395 (1988)
// =============================================================================

#ifndef PHMT101A_H
#define PHMT101A_H

#include "PHMT000A.h"
#include <cmath>
#include <algorithm>
#include <functional>

namespace Sirius {

// =============================================================================
// Shape Function Types
// =============================================================================
enum class WormholeShapeType {
    Ellis,          // b(r) = b₀²/r
    ZeroTidal,      // b(r) = b₀ (constant)
    AbsurdlyBenign, // b(r) = b₀(2 - b₀/r)
    Custom          // User-defined
};

// =============================================================================
// Custom Shape Function Type
// =============================================================================
/// @brief Function signature for custom shape function b(r)
/// @param r Radial coordinate
/// @return Shape function value b(r)
using ShapeFunction = std::function<double(double)>;

/// @brief Function signature for shape function derivative db/dr
using ShapeDerivative = std::function<double(double)>;

// =============================================================================
// Morris-Thorne Family Parameters
// =============================================================================
struct MorrisThorneParams {
    double b0 = 1.0;               // Throat radius
    double Phi0 = 0.0;             // Redshift at throat (0 = zero-tidal)
    WormholeShapeType shapeType = WormholeShapeType::Ellis;

    // Custom shape function callbacks (used when shapeType == Custom)
    ShapeFunction customShapeFunc = nullptr;
    ShapeDerivative customShapeDerivFunc = nullptr;

    // Convenience constructors
    static MorrisThorneParams Ellis(double b0) {
        return {b0, 0.0, WormholeShapeType::Ellis, nullptr, nullptr};
    }
    static MorrisThorneParams ZeroTidal(double b0) {
        return {b0, 0.0, WormholeShapeType::ZeroTidal, nullptr, nullptr};
    }

    /// @brief Create custom wormhole with user-defined shape function
    /// @param b0 Throat radius
    /// @param shapeFunc Shape function b(r), must satisfy b(b0) = b0
    /// @param derivFunc Optional derivative db/dr, computed numerically if null
    static MorrisThorneParams Custom(double b0, ShapeFunction shapeFunc,
                                      ShapeDerivative derivFunc = nullptr) {
        MorrisThorneParams p;
        p.b0 = b0;
        p.Phi0 = 0.0;
        p.shapeType = WormholeShapeType::Custom;
        p.customShapeFunc = std::move(shapeFunc);
        p.customShapeDerivFunc = std::move(derivFunc);
        return p;
    }
};

// =============================================================================
// Morris-Thorne Family Metric Class
// =============================================================================
class MorrisThorneFamily : public IMetric {
public:
    MorrisThorneFamily();
    explicit MorrisThorneFamily(const MorrisThorneParams& params);
    
    // IMetric interface
    void evaluate(const Tensor<double, 4>& pos, Metric4D& g, 
                  Tensor<Dual<double>, 4, 4, 4>& dg) override;
    
    const Config& getParameters() const override { return m_Config; }
    void setParameter(const std::string& key, double value) override;
    const char* getName() const override;
    
    // Morris-Thorne specific methods
    void setParams(const MorrisThorneParams& params);
    MorrisThorneParams getParams() const;
    
    // Shape function b(r)
    double shapeFunction(double r) const;
    double shapeFunctionDerivative(double r) const;

    // Redshift function Φ(r)
    double redshiftFunction(double r) const;
    double redshiftFunctionDerivative(double r) const;

    /// @brief Validate flare-out condition at throat
    /// The flare-out condition ensures the wormhole geometry is traversable:
    ///   (b - r·db/dr) / (2b²) > 0 at r = b₀
    /// Equivalently: b'(b₀) < 1
    /// @return true if flare-out condition is satisfied
    bool validateFlareOutCondition() const;

    /// @brief Get flare-out parameter value at given radius
    /// @param r Radial coordinate
    /// @return (b - r·db/dr) / (2b²)
    double flareOutParameter(double r) const;

private:
    Config m_Config;
    MorrisThorneParams m_params;
};

// =============================================================================
// Inline Implementation
// =============================================================================

inline MorrisThorneFamily::MorrisThorneFamily() {
    m_Config["throat_radius"] = {1.0, 0.1, 100.0};
    m_Config["redshift"] = {0.0, -10.0, 10.0};
    m_params = MorrisThorneParams::Ellis(1.0);
}

inline MorrisThorneFamily::MorrisThorneFamily(const MorrisThorneParams& params) 
    : MorrisThorneFamily() {
    setParams(params);
}

inline void MorrisThorneFamily::setParams(const MorrisThorneParams& params) {
    m_params = params;
    m_Config["throat_radius"].value = params.b0;
    m_Config["redshift"].value = params.Phi0;
}

inline MorrisThorneParams MorrisThorneFamily::getParams() const {
    return m_params;
}

inline void MorrisThorneFamily::setParameter(const std::string& key, double value) {
    if (m_Config.find(key) != m_Config.end()) {
        m_Config[key].value = std::clamp(value, m_Config[key].min, m_Config[key].max);
    }
    m_params.b0 = m_Config["throat_radius"].value;
    m_params.Phi0 = m_Config["redshift"].value;
}

inline const char* MorrisThorneFamily::getName() const {
    switch (m_params.shapeType) {
        case WormholeShapeType::Ellis: return "Ellis Drainhole";
        case WormholeShapeType::ZeroTidal: return "Zero-Tidal Wormhole";
        case WormholeShapeType::AbsurdlyBenign: return "Absurdly Benign Wormhole";
        default: return "Morris-Thorne Wormhole";
    }
}

inline double MorrisThorneFamily::shapeFunction(double r) const {
    double b0 = m_params.b0;
    r = std::max(r, b0);  // Ensure r >= b0

    switch (m_params.shapeType) {
        case WormholeShapeType::Ellis:
            return b0 * b0 / r;
        case WormholeShapeType::ZeroTidal:
            return b0;
        case WormholeShapeType::AbsurdlyBenign:
            return b0 * (2.0 - b0 / r);
        case WormholeShapeType::Custom:
            if (m_params.customShapeFunc) {
                return m_params.customShapeFunc(r);
            }
            return b0 * b0 / r;  // Fallback to Ellis
        default:
            return b0 * b0 / r;  // Default to Ellis
    }
}

inline double MorrisThorneFamily::shapeFunctionDerivative(double r) const {
    double b0 = m_params.b0;
    r = std::max(r, b0);

    switch (m_params.shapeType) {
        case WormholeShapeType::Ellis:
            return -b0 * b0 / (r * r);
        case WormholeShapeType::ZeroTidal:
            return 0.0;
        case WormholeShapeType::AbsurdlyBenign:
            return b0 * b0 / (r * r);
        case WormholeShapeType::Custom:
            // Use user-provided derivative if available
            if (m_params.customShapeDerivFunc) {
                return m_params.customShapeDerivFunc(r);
            }
            // Otherwise compute numerically via central differences
            if (m_params.customShapeFunc) {
                const double h = 1e-6 * r;
                double b_plus = m_params.customShapeFunc(r + h);
                double b_minus = m_params.customShapeFunc(r - h);
                return (b_plus - b_minus) / (2.0 * h);
            }
            return -b0 * b0 / (r * r);
        default:
            return -b0 * b0 / (r * r);
    }
}

inline double MorrisThorneFamily::flareOutParameter(double r) const {
    double b = shapeFunction(r);
    double db_dr = shapeFunctionDerivative(r);

    // Flare-out parameter: (b - r·db/dr) / (2b²)
    // This must be > 0 for traversable wormhole geometry
    double b2 = b * b;
    if (b2 < 1e-20) return 0;

    return (b - r * db_dr) / (2.0 * b2);
}

inline bool MorrisThorneFamily::validateFlareOutCondition() const {
    // Check flare-out condition at throat: b'(b0) < 1
    // Equivalently: flareOutParameter(b0) > 0

    double b0 = m_params.b0;
    double db_dr_at_throat = shapeFunctionDerivative(b0);

    // Condition: b'(b0) < 1
    // For b(b0) = b0 (throat condition), this becomes db/dr < 1
    if (db_dr_at_throat >= 1.0) {
        return false;  // Violates flare-out
    }

    // Also check that flare-out parameter is positive at throat
    double fop = flareOutParameter(b0);
    return fop > 0;
}

inline double MorrisThorneFamily::redshiftFunction([[maybe_unused]] double r) const {
    // For zero-tidal wormhole, Φ = 0 everywhere (r-independent)
    return m_params.Phi0;
}

inline double MorrisThorneFamily::redshiftFunctionDerivative([[maybe_unused]] double r) const {
    // For constant redshift, derivative is 0 (r-independent)
    return 0.0;
}

inline void MorrisThorneFamily::evaluate(const Tensor<double, 4>& pos, Metric4D& g,
                                          Tensor<Dual<double>, 4, 4, 4>& dg) {
    [[maybe_unused]] double t = pos(0);  // Time coordinate (unused in static metric)
    double r = pos(1);
    double theta = pos(2);
    // phi coordinate not needed for spherically symmetric metric evaluation

    double b0 = m_params.b0;
    r = std::max(r, b0 * 1.001);  // Stay just outside throat
    
    double sin_theta = std::sin(theta);
    double sin2 = sin_theta * sin_theta;
    double r2 = r * r;
    
    // Shape and redshift functions
    double b = shapeFunction(r);
    double db_dr = shapeFunctionDerivative(r);
    double Phi = redshiftFunction(r);
    double dPhi_dr = redshiftFunctionDerivative(r);
    
    double exp2Phi = std::exp(2.0 * Phi);
    double one_minus_b_over_r = 1.0 - b / r;
    one_minus_b_over_r = std::max(one_minus_b_over_r, 1e-10);  // Regularize at throat
    
    // Initialize metric
    g.zero();
    
    // g_tt = -e^{2Φ}
    g(0, 0) = Dual<double>(-exp2Phi);
    
    // g_rr = 1/(1 - b/r)
    g(1, 1) = Dual<double>(1.0 / one_minus_b_over_r);
    
    // g_θθ = r²
    g(2, 2) = Dual<double>(r2);
    
    // g_φφ = r² sin²θ
    g(3, 3) = Dual<double>(r2 * sin2);
    
    // Compute derivatives
    dg.zero();
    
    // ∂g_tt/∂r = -2 dΦ/dr e^{2Φ}
    dg(1, 0, 0) = Dual<double>(-2.0 * dPhi_dr * exp2Phi);
    
    // ∂g_rr/∂r = -(1-b/r)^{-2} * d(1-b/r)/dr
    //          = -(1-b/r)^{-2} * (b/r² - db/dr/r)
    //          = (1-b/r)^{-2} * (db/dr/r - b/r²)
    double d_one_minus = db_dr / r - b / r2;
    dg(1, 1, 1) = Dual<double>(-d_one_minus / (one_minus_b_over_r * one_minus_b_over_r));
    
    // ∂g_θθ/∂r = 2r
    dg(1, 2, 2) = Dual<double>(2.0 * r);
    
    // ∂g_φφ/∂r = 2r sin²θ
    dg(1, 3, 3) = Dual<double>(2.0 * r * sin2);
    
    // ∂g_φφ/∂θ = 2r² sinθ cosθ
    double cos_theta = std::cos(theta);
    dg(2, 3, 3) = Dual<double>(2.0 * r2 * sin_theta * cos_theta);
}

} // namespace Sirius

#endif // PHMT101A_H
