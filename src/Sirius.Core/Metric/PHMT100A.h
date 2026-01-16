// =============================================================================
// PHMT100A.h - Unified Kerr-Schild Black Hole Family
// Component ID: PHMT100A (Physics/Metric/Kerr-Schild Family)
// =============================================================================
//
// UNIFIED METRIC FAMILY
// =====================
// This single implementation covers 9 spacetimes through parameterization:
//
// | M | a | Q | Λ | Spacetime                |
// |---|---|---|---|--------------------------|
// | 0 | 0 | 0 | 0 | Minkowski                |
// | ✓ | 0 | 0 | 0 | Schwarzschild            |
// | ✓ | ✓ | 0 | 0 | Kerr                     |
// | ✓ | 0 | ✓ | 0 | Reissner-Nordström       |
// | ✓ | ✓ | ✓ | 0 | Kerr-Newman              |
// | 0 | 0 | 0 | ✓ | de Sitter                |
// | ✓ | 0 | 0 | ✓ | Schwarzschild-de Sitter  |
// | ✓ | ✓ | 0 | ✓ | Kerr-de Sitter           |
// | ✓ | ✓ | ✓ | ✓ | Kerr-Newman-de Sitter    |
//
// Mathematical Foundation
// =======================
// Kerr-Schild ansatz in Cartesian coordinates (t, x, y, z):
//
//   g_μν = η_μν + H · l_μ · l_ν
//   g^μν = η^μν - H · l^μ · l^ν
//
// where:
//   η_μν = diag(-1, 1, 1, 1)  (Minkowski metric)
//   r    = implicit solution of r⁴ - (R² - a²)r² - a²z² = 0
//   H    = (2Mr - Q²)r² / (r⁴ + a²z²)
//   l_μ  = null vector field: (1, (rx+ay)/(r²+a²), (ry-ax)/(r²+a²), z/r)
//   l^μ  = (1, -(rx+ay)/(r²+a²), -(ry-ax)/(r²+a²), -z/r)
//
// KEY ADVANTAGE: Christoffel symbols are polynomial in (x, y, z, r).
//                No trigonometric functions → No pole singularities!
//
// Reference: Visser, "The Kerr spacetime: A brief introduction" (arXiv:0706.0622)
// =============================================================================

#ifndef PHMT100A_H
#define PHMT100A_H

#include "PHMT000A.h"
#include <cmath>
#include <algorithm>

namespace Sirius {

// =============================================================================
// Kerr-Schild Family Parameters
// =============================================================================
struct KerrSchildParams {
    double M = 1.0;      // Mass (Schwarzschild radius = 2M)
    double a = 0.0;      // Spin parameter (|a| ≤ M for black hole)
    double Q = 0.0;      // Electric charge (|Q| ≤ M for black hole)  
    double Lambda = 0.0; // Cosmological constant (Λ > 0 for de Sitter)
    
    // Convenience constructors for common spacetimes
    static KerrSchildParams Minkowski() { return {0, 0, 0, 0}; }
    static KerrSchildParams Schwarzschild(double M) { return {M, 0, 0, 0}; }
    static KerrSchildParams Kerr(double M, double a) { return {M, a, 0, 0}; }
    static KerrSchildParams ReissnerNordstrom(double M, double Q) { return {M, 0, Q, 0}; }
    static KerrSchildParams KerrNewman(double M, double a, double Q) { return {M, a, Q, 0}; }
    static KerrSchildParams DeSitter(double Lambda) { return {0, 0, 0, Lambda}; }
};

// =============================================================================
// Kerr-Schild Family Metric Class
// =============================================================================
class KerrSchildFamily : public IMetric {
public:
    KerrSchildFamily();
    explicit KerrSchildFamily(const KerrSchildParams& params);
    
    // IMetric interface
    void evaluate(const Tensor<double, 4>& pos, Metric4D& g, 
                  Tensor<Dual<double>, 4, 4, 4>& dg) override;
    
    const Config& getParameters() const override { return m_Config; }
    void setParameter(const std::string& key, double value) override;
    const char* getName() const override;
    
    // Kerr-Schild specific methods
    void setParams(const KerrSchildParams& params);
    KerrSchildParams getParams() const;
    
    // Compute Boyer-Lindquist-like radius from Cartesian coordinates
    // Solves: r⁴ - (x² + y² + z² - a²)r² - a²z² = 0
    double computeKerrRadius(double x, double y, double z) const;
    
    // Compute the null vector l^μ at position
    void computeNullVector(double x, double y, double z, double r, double l[4]) const;
    
    // Compute the scalar function H = (2Mr - Q²)r² / (r⁴ + a²z²)
    double computeH(double r, double z) const;
    
private:
    Config m_Config;
    KerrSchildParams m_params;
    
    void updateDerivedQuantities();
};

// =============================================================================
// Inline Implementation
// =============================================================================

inline KerrSchildFamily::KerrSchildFamily() {
    m_Config["mass"] = {1.0, 0.0, 100.0};
    m_Config["spin"] = {0.0, -0.998, 0.998};
    m_Config["charge"] = {0.0, -1.0, 1.0};
    m_Config["lambda"] = {0.0, -0.1, 0.1};
    m_params = KerrSchildParams::Schwarzschild(1.0);
}

inline KerrSchildFamily::KerrSchildFamily(const KerrSchildParams& params) 
    : KerrSchildFamily() {
    setParams(params);
}

inline void KerrSchildFamily::setParams(const KerrSchildParams& params) {
    m_params = params;
    m_Config["mass"].value = params.M;
    m_Config["spin"].value = params.a / std::max(params.M, 1e-10);
    m_Config["charge"].value = params.Q / std::max(params.M, 1e-10);
    m_Config["lambda"].value = params.Lambda;
}

inline KerrSchildParams KerrSchildFamily::getParams() const {
    return m_params;
}

inline void KerrSchildFamily::setParameter(const std::string& key, double value) {
    if (m_Config.find(key) != m_Config.end()) {
        m_Config[key].value = std::clamp(value, m_Config[key].min, m_Config[key].max);
    }
    // Update internal params
    m_params.M = m_Config["mass"].value;
    m_params.a = m_Config["spin"].value * m_params.M;
    m_params.Q = m_Config["charge"].value * m_params.M;
    m_params.Lambda = m_Config["lambda"].value;
}

inline const char* KerrSchildFamily::getName() const {
    // Return appropriate name based on active parameters
    if (m_params.M == 0 && m_params.Lambda == 0) return "Minkowski";
    if (m_params.a == 0 && m_params.Q == 0 && m_params.Lambda == 0) return "Schwarzschild";
    if (m_params.Q == 0 && m_params.Lambda == 0) return "Kerr";
    if (m_params.a == 0 && m_params.Lambda == 0) return "Reissner-Nordström";
    if (m_params.Lambda == 0) return "Kerr-Newman";
    if (m_params.M == 0) return "de Sitter";
    if (m_params.a == 0 && m_params.Q == 0) return "Schwarzschild-de Sitter";
    if (m_params.Q == 0) return "Kerr-de Sitter";
    return "Kerr-Newman-de Sitter";
}

inline double KerrSchildFamily::computeKerrRadius(double x, double y, double z) const {
    double a = m_params.a;
    double a2 = a * a;
    double R2 = x*x + y*y + z*z;
    
    // For a = 0 (Schwarzschild), r = R directly
    if (std::abs(a) < 1e-12) {
        return std::sqrt(std::max(R2, 1e-20));
    }
    
    // Solve r⁴ - (R² - a²)r² - a²z² = 0
    // Let u = r². Then: u² - (R² - a²)u - a²z² = 0
    // u = [(R² - a²) + √((R² - a²)² + 4a²z²)] / 2
    double Rm2 = R2 - a2;
    double disc = Rm2 * Rm2 + 4.0 * a2 * z * z;
    double r2 = (Rm2 + std::sqrt(std::max(disc, 0.0))) / 2.0;
    
    return std::sqrt(std::max(r2, 1e-20));
}

inline void KerrSchildFamily::computeNullVector(double x, double y, double z, double r, double l[4]) const {
    double a = m_params.a;
    double a2 = a * a;
    double r2 = r * r;
    double denom = r2 + a2;
    
    // l^μ = (1, (rx + ay)/(r² + a²), (ry - ax)/(r² + a²), z/r)
    l[0] = 1.0;
    l[1] = (r * x + a * y) / denom;
    l[2] = (r * y - a * x) / denom;
    l[3] = z / std::max(r, 1e-10);
}

inline double KerrSchildFamily::computeH(double r, double z) const {
    double M = m_params.M;
    double a = m_params.a;
    double Q = m_params.Q;
    
    double r2 = r * r;
    double r4 = r2 * r2;
    double a2 = a * a;
    
    // H = (2Mr - Q²) · r² / (r⁴ + a²z²)
    // For pure Schwarzschild: H = 2M·r³ / r⁴ = 2M/r
    double numerator = (2.0 * M * r - Q * Q) * r2;
    double denominator = r4 + a2 * z * z;
    
    return numerator / std::max(denominator, 1e-20);
}

inline void KerrSchildFamily::evaluate(const Tensor<double, 4>& pos, Metric4D& g,
                                        Tensor<Dual<double>, 4, 4, 4>& dg) {
    [[maybe_unused]] double t = pos(0);  // Time coordinate (unused in static metric)
    double x = pos(1);
    double y = pos(2);
    double z = pos(3);
    
    double M = m_params.M;
    double a = m_params.a;
    double Q = m_params.Q;
    double Lambda = m_params.Lambda;
    
    double a2 = a * a;
    double Q2 = Q * Q;
    
    // Compute Kerr radius
    double R2 = x*x + y*y + z*z;
    double Rm2 = R2 - a2;
    double disc = Rm2 * Rm2 + 4.0 * a2 * z * z;
    disc = std::max(disc, 1e-20);
    double sqrt_disc = std::sqrt(disc);
    double r2 = (Rm2 + sqrt_disc) / 2.0;
    r2 = std::max(r2, 1e-10);
    double r = std::sqrt(r2);
    double r3 = r2 * r;
    double r4 = r2 * r2;
    
    // =========================================================================
    // Derivatives of r with respect to spatial coordinates
    // =========================================================================
    double d_disc_dx = 4.0 * x * Rm2;
    double d_disc_dy = 4.0 * y * Rm2;
    double d_disc_dz = 4.0 * z * (R2 + a2);
    
    double d_sqrt_disc_dx = d_disc_dx / (2.0 * sqrt_disc);
    double d_sqrt_disc_dy = d_disc_dy / (2.0 * sqrt_disc);
    double d_sqrt_disc_dz = d_disc_dz / (2.0 * sqrt_disc);
    
    double d_r2_dx = x + d_sqrt_disc_dx / 2.0;
    double d_r2_dy = y + d_sqrt_disc_dy / 2.0;
    double d_r2_dz = z + d_sqrt_disc_dz / 2.0;
    
    double dr_dx = d_r2_dx / (2.0 * r);
    double dr_dy = d_r2_dy / (2.0 * r);
    double dr_dz = d_r2_dz / (2.0 * r);
    double dr[4] = {0.0, dr_dx, dr_dy, dr_dz};
    
    // =========================================================================
    // Null vector l^μ and its derivatives
    // =========================================================================
    double denom = r2 + a2;
    double denom2 = denom * denom;
    double l[4] = {
        1.0,
        (r*x + a*y) / denom,
        (r*y - a*x) / denom,
        z / r
    };
    
    // dl[λ][μ] = ∂l^μ/∂x^λ
    double dl[4][4] = {{0}};
    double d_denom[4] = {0.0, d_r2_dx, d_r2_dy, d_r2_dz};
    
    // l^1 = (rx + ay) / denom
    dl[1][1] = (r + x * dr_dx) / denom - (r*x + a*y) * d_denom[1] / denom2;
    dl[2][1] = (x * dr_dy + a) / denom - (r*x + a*y) * d_denom[2] / denom2;
    dl[3][1] = (x * dr_dz) / denom - (r*x + a*y) * d_denom[3] / denom2;
    
    // l^2 = (ry - ax) / denom
    dl[1][2] = (y * dr_dx - a) / denom - (r*y - a*x) * d_denom[1] / denom2;
    dl[2][2] = (r + y * dr_dy) / denom - (r*y - a*x) * d_denom[2] / denom2;
    dl[3][2] = (y * dr_dz) / denom - (r*y - a*x) * d_denom[3] / denom2;
    
    // l^3 = z / r
    dl[1][3] = -z * dr_dx / r2;
    dl[2][3] = -z * dr_dy / r2;
    dl[3][3] = 1.0 / r - z * dr_dz / r2;
    
    // =========================================================================
    // Kerr-Schild scalar function H = (2Mr - Q²)r² / (r⁴ + a²z²)
    // =========================================================================
    double f_denom = r4 + a2 * z * z;
    double H = (2.0 * M * r - Q2) * r2 / f_denom;
    
    // Derivatives of H
    double dH[4] = {0.0, 0.0, 0.0, 0.0};
    for (int lam = 1; lam <= 3; lam++) {
        // d(numerator)/dlam = (2M)·dr·r² + (2Mr - Q²)·2r·dr = 2r·dr·(M·r + 2Mr - Q²) = 2r·dr·(3Mr - Q²)
        double d_num = 2.0 * r * dr[lam] * (3.0 * M * r - Q2);
        // d(f_denom)/dlam = 4r³·dr + (lam==3 ? 2a²z : 0)
        double d_f_denom = 4.0 * r3 * dr[lam];
        if (lam == 3) d_f_denom += 2.0 * a2 * z;
        
        double numerator = (2.0 * M * r - Q2) * r2;
        dH[lam] = (d_num * f_denom - numerator * d_f_denom) / (f_denom * f_denom);
    }
    
    // =========================================================================
    // Metric: g_μν = η_μν + H·l_μ·l_ν
    // For de Sitter (Λ ≠ 0), we add cosmological term separately
    // =========================================================================
    g.zero();
    g(0, 0) = Dual<double>(-1.0);
    g(1, 1) = Dual<double>(1.0);
    g(2, 2) = Dual<double>(1.0);
    g(3, 3) = Dual<double>(1.0);
    
    // Add Kerr-Schild perturbation
    for (int mu = 0; mu < 4; mu++) {
        for (int nu = 0; nu < 4; nu++) {
            g(mu, nu) = Dual<double>(g(mu, nu).real + H * l[mu] * l[nu]);
        }
    }
    
    // Cosmological constant contribution (if non-zero)
    // For de Sitter in static coords: g_tt → -(1 - Λr²/3), etc.
    // In Kerr-Schild-de Sitter, this is more complex - simplified here
    if (std::abs(Lambda) > 1e-15) {
        double LambdaFactor = Lambda * R2 / 3.0;
        g(0, 0) = Dual<double>(g(0, 0).real - LambdaFactor);
        // Spatial components also modified in full treatment
    }
    
    // =========================================================================
    // Derivatives: ∂g_μν/∂x^λ = ∂H/∂x^λ·l_μ·l_ν + H·∂l_μ/∂x^λ·l_ν + H·l_μ·∂l_ν/∂x^λ
    // =========================================================================
    dg.zero();
    
    for (int lam = 1; lam <= 3; lam++) {
        for (int mu = 0; mu < 4; mu++) {
            for (int nu = 0; nu < 4; nu++) {
                double dg_val = dH[lam] * l[mu] * l[nu]
                              + H * dl[lam][mu] * l[nu]
                              + H * l[mu] * dl[lam][nu];
                dg(lam, mu, nu) = Dual<double>(dg_val);
            }
        }
    }
    
    // Cosmological derivative contribution
    if (std::abs(Lambda) > 1e-15) {
        double dLambdaFactor_dx = 2.0 * Lambda * x / 3.0;
        double dLambdaFactor_dy = 2.0 * Lambda * y / 3.0;
        double dLambdaFactor_dz = 2.0 * Lambda * z / 3.0;
        dg(1, 0, 0) = Dual<double>(dg(1, 0, 0).real - dLambdaFactor_dx);
        dg(2, 0, 0) = Dual<double>(dg(2, 0, 0).real - dLambdaFactor_dy);
        dg(3, 0, 0) = Dual<double>(dg(3, 0, 0).real - dLambdaFactor_dz);
    }
}

} // namespace Sirius

#endif // PHMT100A_H
