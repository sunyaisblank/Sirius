// PHMT100B.h - Kerr-Schild Family Double-Precision Implementation
// Component ID: PHMT100B
// Purpose: Double-precision Kerr metric for offline rendering
//
// MATHEMATICAL BASIS:
// Boyer-Lindquist coordinates (t, r, θ, φ) with metric:
//   ds² = -(1 - 2Mr/Σ)dt² - (4aMr sin²θ/Σ)dt dφ + (Σ/Δ)dr² + Σ dθ²
//        + sin²θ[(r² + a²)² - a²Δ sin²θ]/Σ dφ²
//
// where: Σ = r² + a²cos²θ,  Δ = r² - 2Mr + a²
//
// REFERENCE: Misner, Thorne & Wheeler (1973), James et al. (2015)

#pragma once

#include "PHMT000B.h"
#include <cmath>
#include <algorithm>

namespace sirius::physics {

//==============================================================================
// KerrMetricD: Double-precision Kerr metric implementation
//==============================================================================

class KerrMetricD : public IMetricD {
public:
    explicit KerrMetricD(const MetricParamsD& params)
        : m_M(params.M), m_a(params.a), m_Q(params.Q)
        , m_rplus(params.rplus), m_rminus(params.rminus) {
        // Clamp spin to avoid naked singularity
        if (std::abs(m_a) > 0.9999 * m_M) {
            m_a = std::copysign(0.9999 * m_M, m_a);
        }
    }
    
    KerrMetricD(double M = 1.0, double a = 0.0, double Q = 0.0) 
        : KerrMetricD(MetricParamsD(M, a, Q)) {}
    
    //--------------------------------------------------------------------------
    // IMetricD Implementation
    //--------------------------------------------------------------------------
    
    void evaluate(const Vec4d& x, double g[4][4], double g_inv[4][4]) const override {
        double r = x.r;
        double theta = x.theta;
        
        // Handle pole singularities
        double sinth = std::sin(theta);
        double costh = std::cos(theta);
        if (std::abs(sinth) < 1e-10) {
            sinth = std::copysign(1e-10, sinth);
        }
        double sin2th = sinth * sinth;
        double cos2th = costh * costh;
        
        // Kerr metric functions
        double Sigma = r*r + m_a*m_a * cos2th;
        double Delta = r*r - 2*m_M*r + m_a*m_a;
        double A = (r*r + m_a*m_a)*(r*r + m_a*m_a) - m_a*m_a * Delta * sin2th;
        
        // Initialize to zero
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
                g[i][j] = g_inv[i][j] = 0.0;
        
        // Covariant metric g_μν
        g[0][0] = -(1 - 2*m_M*r/Sigma);              // g_tt
        g[0][3] = g[3][0] = -2*m_M*m_a*r*sin2th/Sigma;  // g_tφ
        g[1][1] = Sigma/Delta;                        // g_rr
        g[2][2] = Sigma;                              // g_θθ
        g[3][3] = A*sin2th/Sigma;                     // g_φφ
        
        // Contravariant metric g^μν
        // For block-diagonal structure, invert 2x2 (t,φ) and 1x1 (r), (θ) blocks separately
        // det(2x2 t-φ block) = g_tt * g_φφ - g_tφ²
        double det_tphi = g[0][0] * g[3][3] - g[0][3] * g[0][3];
        if (std::abs(det_tphi) < 1e-20) det_tphi = std::copysign(1e-20, det_tphi);
        
        // 2x2 inverse: [[a,b],[b,c]]^-1 = (1/det) * [[c,-b],[-b,a]]
        g_inv[0][0] = g[3][3] / det_tphi;
        g_inv[0][3] = g_inv[3][0] = -g[0][3] / det_tphi;
        g_inv[3][3] = g[0][0] / det_tphi;
        g_inv[1][1] = 1.0 / g[1][1];
        g_inv[2][2] = 1.0 / g[2][2];
    }
    
    void christoffel(const Vec4d& x, double Gamma[4][4][4]) const override {
        // Initialize to zero
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
                for (int k = 0; k < 4; ++k)
                    Gamma[i][j][k] = 0.0;
        
        double r = x.r;
        double theta = x.theta;
        
        double sinth = std::sin(theta);
        double costh = std::cos(theta);
        if (std::abs(sinth) < 1e-10) sinth = std::copysign(1e-10, sinth);
        double sin2th = sinth * sinth;
        double cos2th = costh * costh;
        
        double a = m_a;
        double M = m_M;
        double a2 = a*a;
        double r2 = r*r;
        
        double Sigma = r2 + a2*cos2th;
        double Delta = r2 - 2*M*r + a2;
        double Sigma2 = Sigma*Sigma;
        double Sigma3 = Sigma2*Sigma;
        
        // Partial derivatives of Sigma
        double dSigma_dr = 2*r;
        double dSigma_dth = -2*a2*costh*sinth;
        
        // Partial derivatives of Delta
        double dDelta_dr = 2*r - 2*M;
        
        // Key metric components
        double A = (r2 + a2)*(r2 + a2) - a2*Delta*sin2th;
        
        // Christoffel symbols (analytic formulas for Kerr)
        // Using MTW conventions
        
        // Γ^t components
        double Mr_term = M*(r2 - a2*cos2th);
        Gamma[0][0][1] = Mr_term * (r2 + a2) / (Sigma2 * Delta);
        Gamma[0][0][2] = -2*M*a2*r*sinth*costh / Sigma2;
        Gamma[0][1][3] = -a*sin2th * (Mr_term - Sigma*r) / (Sigma2 * Delta) 
                         + a*M*(r2 + a2)*sin2th / (Sigma2 * Delta);
        Gamma[0][2][3] = 2*M*a*r*(r2 + a2)*sinth*costh / Sigma2;
        
        // Symmetrize
        Gamma[0][1][0] = Gamma[0][0][1];
        Gamma[0][2][0] = Gamma[0][0][2];
        Gamma[0][3][1] = Gamma[0][1][3];
        Gamma[0][3][2] = Gamma[0][2][3];
        
        // Γ^r components
        Gamma[1][0][0] = M*Delta*(r2 - a2*cos2th) / Sigma3;
        Gamma[1][0][3] = -M*a*Delta*sin2th*(r2 - a2*cos2th) / Sigma3;
        Gamma[1][1][1] = (r*(a2 - r2) + M*(r2 - a2*cos2th)) / (Sigma * Delta);
        Gamma[1][1][2] = -a2*sinth*costh / Sigma;
        Gamma[1][2][2] = -r*Delta / Sigma;
        Gamma[1][3][3] = -Delta*sin2th*(r*Sigma2 + M*a2*sin2th*(a2*cos2th - r2)) / Sigma3;
        
        Gamma[1][3][0] = Gamma[1][0][3];
        Gamma[1][2][1] = Gamma[1][1][2];
        
        // Γ^θ components
        Gamma[2][0][0] = -2*M*a2*r*sinth*costh / Sigma3;
        Gamma[2][0][3] = 2*M*a*r*(r2 + a2)*sinth*costh / Sigma3;
        Gamma[2][1][1] = a2*sinth*costh / (Sigma * Delta);
        Gamma[2][1][2] = r / Sigma;
        Gamma[2][2][2] = -a2*sinth*costh / Sigma;
        Gamma[2][3][3] = -sinth*costh*(A*Sigma + 2*M*a2*r*sin2th*(r2 + a2)) / Sigma3;
        
        Gamma[2][3][0] = Gamma[2][0][3];
        Gamma[2][2][1] = Gamma[2][1][2];
        
        // Γ^φ components
        Gamma[3][0][1] = M*a*(a2*cos2th - r2) / (Sigma2 * Delta);
        Gamma[3][0][2] = -2*M*a*r*costh / (Sigma2 * sinth);
        Gamma[3][1][3] = (r*Sigma2 - M*(r2 - a2*cos2th)*(r2 + a2)/Delta) / (Sigma2 * Sigma);
        Gamma[3][2][3] = (a2*sinth*costh*Sigma + 2*M*a2*r*sin2th*costh/sinth) / Sigma2;
        
        Gamma[3][1][0] = Gamma[3][0][1];
        Gamma[3][2][0] = Gamma[3][0][2];
        Gamma[3][3][1] = Gamma[3][1][3];
        Gamma[3][3][2] = Gamma[3][2][3];
    }
    
    void riemann(const Vec4d& x, double R[4][4][4][4]) const override {
        // Initialize to zero
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
                for (int k = 0; k < 4; ++k)
                    for (int l = 0; l < 4; ++l)
                        R[i][j][k][l] = 0.0;
        
        // Compute Riemann tensor via finite differences of Christoffel symbols
        // R^μ_νρσ = ∂_ρ Γ^μ_νσ - ∂_σ Γ^μ_νρ + Γ^μ_λρ Γ^λ_νσ - Γ^μ_λσ Γ^λ_νρ
        
        const double h = 1e-6;  // Step size for finite differences
        
        // Get Christoffel at current point
        double Gamma[4][4][4];
        christoffel(x, Gamma);
        
        // Compute derivatives via central differences
        double dGamma_dr[4][4][4], dGamma_dth[4][4][4];
        Vec4d x_plus, x_minus;
        double Gamma_plus[4][4][4], Gamma_minus[4][4][4];
        
        // ∂Γ/∂r
        x_plus = x; x_plus.r += h;
        x_minus = x; x_minus.r -= h;
        christoffel(x_plus, Gamma_plus);
        christoffel(x_minus, Gamma_minus);
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
                for (int k = 0; k < 4; ++k)
                    dGamma_dr[i][j][k] = (Gamma_plus[i][j][k] - Gamma_minus[i][j][k]) / (2*h);
        
        // ∂Γ/∂θ
        x_plus = x; x_plus.theta += h;
        x_minus = x; x_minus.theta -= h;
        christoffel(x_plus, Gamma_plus);
        christoffel(x_minus, Gamma_minus);
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
                for (int k = 0; k < 4; ++k)
                    dGamma_dth[i][j][k] = (Gamma_plus[i][j][k] - Gamma_minus[i][j][k]) / (2*h);
        
        // Compute Riemann components
        for (int mu = 0; mu < 4; ++mu) {
            for (int nu = 0; nu < 4; ++nu) {
                for (int rho = 0; rho < 4; ++rho) {
                    for (int sigma = 0; sigma < 4; ++sigma) {
                        double val = 0;
                        
                        // ∂_ρ Γ^μ_νσ - ∂_σ Γ^μ_νρ
                        if (rho == 1) val += dGamma_dr[mu][nu][sigma];
                        if (rho == 2) val += dGamma_dth[mu][nu][sigma];
                        if (sigma == 1) val -= dGamma_dr[mu][nu][rho];
                        if (sigma == 2) val -= dGamma_dth[mu][nu][rho];
                        
                        // Γ^μ_λρ Γ^λ_νσ - Γ^μ_λσ Γ^λ_νρ
                        for (int lambda = 0; lambda < 4; ++lambda) {
                            val += Gamma[mu][lambda][rho] * Gamma[lambda][nu][sigma];
                            val -= Gamma[mu][lambda][sigma] * Gamma[lambda][nu][rho];
                        }
                        
                        R[mu][nu][rho][sigma] = val;
                    }
                }
            }
        }
    }
    
    double hamiltonian(const Vec4d& q, const Vec4d& p) const override {
        double g_inv[4][4], g[4][4];
        evaluate(q, g, g_inv);
        
        // H = (1/2) g^μν p_μ p_ν
        double H = 0;
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
                H += g_inv[i][j] * p[i] * p[j];
        return 0.5 * H;
    }
    
    Vec4d dHdp(const Vec4d& q, const Vec4d& p) const override {
        // ∂H/∂p_μ = g^μν p_ν
        double g_inv[4][4], g[4][4];
        evaluate(q, g, g_inv);
        
        Vec4d result;
        for (int i = 0; i < 4; ++i) {
            result[i] = 0;
            for (int j = 0; j < 4; ++j)
                result[i] += g_inv[i][j] * p[j];
        }
        return result;
    }
    
    Vec4d dHdq(const Vec4d& q, const Vec4d& p) const override {
        // Compute dp/dλ FULLY ANALYTICALLY using geodesic equation form:
        // dp_μ/dλ = (1/2)(∂g_αβ/∂q^μ) v^α v^β
        // where v^α = g^αβ p_β is the contravariant 4-velocity
        //
        // This is equivalent to dp/dλ = -∂H/∂q but uses simpler covariant
        // metric derivatives instead of inverse metric derivatives.
        //
        // Reference: Chandrasekhar "Mathematical Theory of Black Holes"
        
        double r = q.r;
        double theta = q.theta;
        double M = m_M;
        double a = m_a;
        double a2 = a*a;
        double r2 = r*r;
        
        double sinth = std::sin(theta);
        double costh = std::cos(theta);
        if (std::abs(sinth) < 1e-10) sinth = std::copysign(1e-10, sinth);
        double sin2th = sinth * sinth;
        double cos2th = costh * costh;
        double sin2theta = 2*sinth*costh;  // sin(2θ)
        
        // Kerr metric auxiliary functions
        double Sigma = r2 + a2*cos2th;
        double Delta = r2 - 2*M*r + a2;
        double A = (r2 + a2)*(r2 + a2) - a2*Delta*sin2th;
        
        double Sigma2 = Sigma * Sigma;
        double Delta2 = Delta * Delta;
        
        // Derivatives of auxiliary functions
        double dSigma_dr = 2*r;
        double dSigma_dth = -a2*sin2theta;
        double dDelta_dr = 2*(r - M);
        double dA_dr = 4*r*(r2 + a2) - a2*dDelta_dr*sin2th;
        double dA_dth = -a2*Delta*sin2theta;
        
        // ============================================================
        // Covariant metric derivatives ∂g_αβ/∂r
        // ============================================================
        
        // g_tt = -(1 - 2Mr/Σ) = -1 + 2Mr/Σ
        // ∂g_tt/∂r = 2M(Σ - r·dΣ/dr)/Σ² = 2M(Σ - 2r²)/Σ²
        double dg_tt_dr = 2*M*(Sigma - 2*r2)/Sigma2;
        
        // g_tφ = -2Mar·sin²θ/Σ
        // ∂g_tφ/∂r = -2Ma·sin²θ·(Σ - r·dΣ/dr)/Σ² = -2Ma·sin²θ·(Σ - 2r²)/Σ²
        double dg_tphi_dr = -2*M*a*sin2th*(Sigma - 2*r2)/Sigma2;
        
        // g_rr = Σ/Δ
        // ∂g_rr/∂r = (dΣ/dr·Δ - Σ·dΔ/dr)/Δ² = (2r·Δ - Σ·2(r-M))/Δ²
        double dg_rr_dr = (2*r*Delta - Sigma*dDelta_dr)/Delta2;
        
        // g_θθ = Σ
        // ∂g_θθ/∂r = dΣ/dr = 2r
        double dg_thth_dr = dSigma_dr;
        
        // g_φφ = A·sin²θ/Σ
        // ∂g_φφ/∂r = sin²θ·(dA/dr·Σ - A·dΣ/dr)/Σ²
        double dg_phiphi_dr = sin2th*(dA_dr*Sigma - A*dSigma_dr)/Sigma2;
        
        // ============================================================
        // Covariant metric derivatives ∂g_αβ/∂θ
        // ============================================================
        
        // ∂g_tt/∂θ = -2Mr·dΣ/dθ/Σ² = 2Mra²sin2θ/Σ²
        double dg_tt_dth = -2*M*r*dSigma_dth/Sigma2;
        
        // ∂g_tφ/∂θ = -2Mar·[2sinθcosθ/Σ - sin²θ·dΣ/dθ/Σ²]
        //          = -2Mar·[sin2θ/Σ + sin²θ·a²sin2θ/Σ²]
        //          = -2Mar·sin2θ·(Σ + a²sin²θ)/Σ²
        double dg_tphi_dth = -2*M*a*r*sin2theta*(Sigma + a2*sin2th)/Sigma2;
        
        // ∂g_rr/∂θ = (dΣ/dθ)/Δ = (-a²sin2θ)/Δ
        double dg_rr_dth = dSigma_dth/Delta;
        
        // ∂g_θθ/∂θ = dΣ/dθ = -a²sin2θ
        double dg_thth_dth = dSigma_dth;
        
        // ∂g_φφ/∂θ = ∂(A·sin²θ/Σ)/∂θ
        //          = (dA/dθ·sin²θ + A·sin2θ)/Σ - A·sin²θ·dΣ/dθ/Σ²
        double dg_phiphi_dth = (dA_dth*sin2th + A*sin2theta)/Sigma 
                             - A*sin2th*dSigma_dth/Sigma2;
        
        // ============================================================
        // Get inverse metric and compute contravariant velocity v^α = g^αβ p_β
        // ============================================================
        
        double g[4][4], g_inv[4][4];
        evaluate(q, g, g_inv);
        
        double vt = g_inv[0][0]*p.t + g_inv[0][3]*p.phi;
        double vr = g_inv[1][1]*p.r;
        double vth = g_inv[2][2]*p.theta;
        double vphi = g_inv[3][0]*p.t + g_inv[3][3]*p.phi;
        
        // ============================================================
        // Compute dp_μ/dλ = (1/2)(∂g_αβ/∂q^μ) v^α v^β
        // ============================================================
        
        // For μ = r:
        double dp_r = 0.5 * (
            dg_tt_dr * vt * vt +
            2 * dg_tphi_dr * vt * vphi +
            dg_rr_dr * vr * vr +
            dg_thth_dr * vth * vth +
            dg_phiphi_dr * vphi * vphi
        );
        
        // For μ = θ:
        double dp_th = 0.5 * (
            dg_tt_dth * vt * vt +
            2 * dg_tphi_dth * vt * vphi +
            dg_rr_dth * vr * vr +
            dg_thth_dth * vth * vth +
            dg_phiphi_dth * vphi * vphi
        );
        
        // Return dp/dλ
        Vec4d result;
        result.t = 0;       // dp_t/dλ = 0 (∂g/∂t = 0, stationary)
        result.r = dp_r;
        result.theta = dp_th;
        result.phi = 0;     // dp_φ/dλ = 0 (∂g/∂φ = 0, axisymmetric)
        
        return result;
    }
    
    bool isValid(const Vec4d& x) const override {
        // Check radial coordinate
        if (x.r <= m_rplus * 1.001) return false;  // Inside or at horizon
        if (x.r > 1e6 * m_M) return false;  // Too far (numerical issues)
        
        // Check angular coordinate
        if (x.theta <= 1e-6 || x.theta >= M_PI - 1e-6) return false;  // At poles
        
        return true;
    }
    
    double horizonRadius() const override { return m_rplus; }
    double innerHorizonRadius() const override { return m_rminus; }
    
    double ergosphereRadius(double theta) const override {
        // r_ergo = M + √(M² - a²cos²θ)
        double cos2th = std::cos(theta) * std::cos(theta);
        return m_M + std::sqrt(m_M*m_M - m_a*m_a*cos2th);
    }
    
    double iscoRadius() const override {
        // ISCO for prograde orbit
        if (std::abs(m_a) < 1e-10) return 6.0 * m_M;  // Schwarzschild
        
        double a_star = m_a / m_M;
        double Z1 = 1 + std::cbrt(1 - a_star*a_star) * 
                    (std::cbrt(1 + a_star) + std::cbrt(1 - a_star));
        double Z2 = std::sqrt(3*a_star*a_star + Z1*Z1);
        return m_M * (3 + Z2 - std::sqrt((3 - Z1)*(3 + Z1 + 2*Z2)));
    }
    
    double photonSphereRadius() const override {
        // Prograde photon sphere for Kerr
        if (std::abs(m_a) < 1e-10) return 3.0 * m_M;  // Schwarzschild
        
        double a_star = m_a / m_M;
        return 2 * m_M * (1 + std::cos(2.0/3.0 * std::acos(-a_star)));
    }
    
    double mass() const override { return m_M; }
    double spin() const override { return m_a; }
    double charge() const override { return m_Q; }
    
    double timeTransformationFunction(const Vec4d& q) const override {
        // g(q) = Σ for TTESI regularisation
        double cos2th = std::cos(q.theta) * std::cos(q.theta);
        return q.r*q.r + m_a*m_a*cos2th;
    }
    
    double gFactor(const Vec4d& x, const Vec4d& k, const Vec4d& u_emitter) const override {
        // g = (k · u_obs) / (k · u_emit)
        // For static observer: u_obs = (1/α, 0, 0, ω/α) where α = lapse, ω = frame-dragging
        
        double g[4][4], g_inv[4][4];
        evaluate(x, g, g_inv);
        
        // Static observer 4-velocity
        double alpha2 = -1.0 / g_inv[0][0];
        double alpha = std::sqrt(std::max(alpha2, 1e-10));
        double omega = -g_inv[0][3] / g_inv[0][0];
        
        Vec4d u_obs;
        u_obs.t = 1.0 / alpha;
        u_obs.r = 0;
        u_obs.theta = 0;
        u_obs.phi = omega / alpha;
        
        // Contract with covariant k
        double k_dot_u_obs = k.t * u_obs.t + k.phi * u_obs.phi;  // Only t, φ non-zero
        double k_dot_u_emit = 0;
        for (int i = 0; i < 4; ++i)
            k_dot_u_emit += k[i] * u_emitter[i];
        
        if (std::abs(k_dot_u_emit) < 1e-20) return 1.0;
        return k_dot_u_obs / k_dot_u_emit;
    }
    
private:
    double m_M;       // Mass
    double m_a;       // Spin
    double m_Q;       // Charge (for future Kerr-Newman extension)
    double m_rplus;   // Outer horizon
    double m_rminus;  // Inner horizon
};

} // namespace sirius::physics
