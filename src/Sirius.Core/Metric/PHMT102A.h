// =============================================================================
// PHMT102A.h - Unified Warp Drive Family
// Component ID: PHMT102A (Physics/Metric/Warp Drive Family)
// =============================================================================
//
// UNIFIED WARP DRIVE FAMILY
// =========================
// This implementation covers Alcubierre-class warp geometries:
//
// | Shape f(rs) | Shift | Spacetime                |
// |-------------|-------|--------------------------|
// | Top-hat     | v_s   | Classic Alcubierre       |
// | Smooth      | v_s   | Van den Broeck           |
// | Zero-exp    | N^i   | Natário                  |
//
// Mathematical Foundation
// =======================
// Warp drive metric in Cartesian coordinates (t, x, y, z):
//
//   ds² = -dt² + (dx - v_s(t) f(r_s) dt)² + dy² + dz²
//
// where:
//   r_s = √((x - x_s(t))² + y² + z²)  (distance from bubble center)
//   x_s(t) = trajectory of bubble center
//   v_s = dx_s/dt = bubble velocity
//   f(r_s) = shape function (1 inside, 0 outside, smooth transition)
//
// Standard Alcubierre shape function:
//   f(r_s) = [tanh(σ(r_s + R)) - tanh(σ(r_s - R))] / [2 tanh(σR)]
//
// KEY PROPERTY: Already in Cartesian coordinates - no poles!
//
// Reference: Alcubierre, Class. Quantum Grav. 11, L73 (1994)
// =============================================================================

#ifndef PHMT102A_H
#define PHMT102A_H

#include "PHMT000A.h"
#include <cmath>
#include <algorithm>

namespace Sirius {

// =============================================================================
// Warp Drive Family Parameters
// =============================================================================
struct WarpDriveParams {
    double vs = 1.0;      // Bubble velocity (can be > 1 for superluminal)
    double sigma = 8.0;   // Wall thickness parameter (larger = sharper)
    double R = 1.0;       // Bubble radius
    double xs = 0.0;      // Current x position of bubble center
    double ys = 0.0;      // Current y position
    double zs = 0.0;      // Current z position
    
    // Convenience constructors
    static WarpDriveParams Alcubierre(double vs, double R, double sigma = 8.0) {
        return {vs, sigma, R, 0.0, 0.0, 0.0};
    }
    static WarpDriveParams Subluminal(double R) {
        return {0.5, 8.0, R, 0.0, 0.0, 0.0};
    }
    static WarpDriveParams Superluminal(double vs, double R) {
        return {vs, 8.0, R, 0.0, 0.0, 0.0};
    }
};

// =============================================================================
// Warp Drive Family Metric Class
// =============================================================================
class WarpDriveFamily : public IMetric {
public:
    WarpDriveFamily();
    explicit WarpDriveFamily(const WarpDriveParams& params);
    
    // IMetric interface
    void evaluate(const Tensor<double, 4>& pos, Metric4D& g, 
                  Tensor<Dual<double>, 4, 4, 4>& dg) override;
    
    const Config& getParameters() const override { return m_Config; }
    void setParameter(const std::string& key, double value) override;
    const char* getName() const override { return "Alcubierre Warp Drive"; }
    
    // Warp drive specific methods
    void setParams(const WarpDriveParams& params);
    WarpDriveParams getParams() const;
    
    // Shape function f(r_s)
    double shapeFunction(double rs) const;
    double shapeFunctionDerivative(double rs) const;
    
    // Update bubble position (for animation)
    void updateBubblePosition(double t);
    
private:
    Config m_Config;
    WarpDriveParams m_params;
};

// =============================================================================
// Inline Implementation
// =============================================================================

inline WarpDriveFamily::WarpDriveFamily() {
    m_Config["velocity"] = {1.0, 0.0, 10.0};
    m_Config["sigma"] = {8.0, 1.0, 50.0};
    m_Config["radius"] = {1.0, 0.1, 100.0};
    m_params = WarpDriveParams::Alcubierre(1.0, 1.0);
}

inline WarpDriveFamily::WarpDriveFamily(const WarpDriveParams& params) 
    : WarpDriveFamily() {
    setParams(params);
}

inline void WarpDriveFamily::setParams(const WarpDriveParams& params) {
    m_params = params;
    m_Config["velocity"].value = params.vs;
    m_Config["sigma"].value = params.sigma;
    m_Config["radius"].value = params.R;
}

inline WarpDriveParams WarpDriveFamily::getParams() const {
    return m_params;
}

inline void WarpDriveFamily::setParameter(const std::string& key, double value) {
    if (m_Config.find(key) != m_Config.end()) {
        m_Config[key].value = std::clamp(value, m_Config[key].min, m_Config[key].max);
    }
    m_params.vs = m_Config["velocity"].value;
    m_params.sigma = m_Config["sigma"].value;
    m_params.R = m_Config["radius"].value;
}

inline double WarpDriveFamily::shapeFunction(double rs) const {
    double sigma = m_params.sigma;
    double R = m_params.R;
    
    // Alcubierre shape: smooth top-hat
    double tanh_plus = std::tanh(sigma * (rs + R));
    double tanh_minus = std::tanh(sigma * (rs - R));
    double tanh_R = std::tanh(sigma * R);
    
    return (tanh_plus - tanh_minus) / (2.0 * tanh_R);
}

inline double WarpDriveFamily::shapeFunctionDerivative(double rs) const {
    double sigma = m_params.sigma;
    double R = m_params.R;
    
    double sech2_plus = 1.0 / std::cosh(sigma * (rs + R));
    sech2_plus *= sech2_plus;
    double sech2_minus = 1.0 / std::cosh(sigma * (rs - R));
    sech2_minus *= sech2_minus;
    double tanh_R = std::tanh(sigma * R);
    
    return sigma * (sech2_plus - sech2_minus) / (2.0 * tanh_R);
}

inline void WarpDriveFamily::updateBubblePosition(double t) {
    // Simple linear motion along x-axis
    m_params.xs = m_params.vs * t;
}

inline void WarpDriveFamily::evaluate(const Tensor<double, 4>& pos, Metric4D& g,
                                       Tensor<Dual<double>, 4, 4, 4>& dg) {
    // Note: t is unused directly; bubble position updated separately via updateBubblePosition()
    [[maybe_unused]] double t = pos(0);
    double x = pos(1);
    double y = pos(2);
    double z = pos(3);
    
    double vs = m_params.vs;
    double xs = m_params.xs;
    double ys = m_params.ys;
    double zs = m_params.zs;
    
    // Distance from bubble center
    double dx = x - xs;
    double dy = y - ys;
    double dz = z - zs;
    double rs = std::sqrt(dx*dx + dy*dy + dz*dz);
    rs = std::max(rs, 1e-10);
    
    // Shape function and derivative
    double f = shapeFunction(rs);
    double df_drs = shapeFunctionDerivative(rs);
    
    // Derivatives of rs with respect to position
    double drs_dx = dx / rs;
    double drs_dy = dy / rs;
    double drs_dz = dz / rs;
    
    // Derivatives of f with respect to position
    double df_dx = df_drs * drs_dx;
    double df_dy = df_drs * drs_dy;
    double df_dz = df_drs * drs_dz;
    
    // Initialize metric (Minkowski + warp perturbation)
    g.zero();
    
    // Metric: ds² = -(1 - v_s²f²)dt² - 2v_sf dx dt + dx² + dy² + dz²
    double vsf = vs * f;
    double vsf2 = vsf * vsf;
    
    g(0, 0) = Dual<double>(-(1.0 - vsf2));    // g_tt
    g(0, 1) = Dual<double>(-vsf);              // g_tx
    g(1, 0) = Dual<double>(-vsf);              // g_xt
    g(1, 1) = Dual<double>(1.0);               // g_xx
    g(2, 2) = Dual<double>(1.0);               // g_yy
    g(3, 3) = Dual<double>(1.0);               // g_zz
    
    // Metric derivatives
    dg.zero();
    
    // ∂g_tt/∂x = 2 v_s² f df/dx
    dg(1, 0, 0) = Dual<double>(2.0 * vs * vs * f * df_dx);
    dg(2, 0, 0) = Dual<double>(2.0 * vs * vs * f * df_dy);
    dg(3, 0, 0) = Dual<double>(2.0 * vs * vs * f * df_dz);
    
    // ∂g_tx/∂x = -v_s df/dx
    dg(1, 0, 1) = Dual<double>(-vs * df_dx);
    dg(1, 1, 0) = Dual<double>(-vs * df_dx);
    dg(2, 0, 1) = Dual<double>(-vs * df_dy);
    dg(2, 1, 0) = Dual<double>(-vs * df_dy);
    dg(3, 0, 1) = Dual<double>(-vs * df_dz);
    dg(3, 1, 0) = Dual<double>(-vs * df_dz);
}

} // namespace Sirius

#endif // PHMT102A_H
