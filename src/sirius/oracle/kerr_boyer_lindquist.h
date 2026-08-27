// Double-precision Kerr metric in Boyer-Lindquist coordinates with analytic
// Christoffel and Riemann tensors: the oracle's analytic reference for the
// live single-precision Kerr-Schild path. Off the render path.
//   ds^2 = -(1 - 2Mr/Sigma) dt^2 - (4aMr sin^2 th / Sigma) dt dphi
//          + (Sigma/Delta) dr^2 + Sigma dth^2
//          + sin^2 th [(r^2+a^2)^2 - a^2 Delta sin^2 th]/Sigma dphi^2
//   Sigma = r^2 + a^2 cos^2 th,  Delta = r^2 - 2Mr + a^2
// Reference: Misner, Thorne & Wheeler (1973); James et al. (2015).
// Ported from PHMT100B.h; numeric content bit-identical.

#pragma once

#include "sirius/base/contracts.h"
#include "sirius/oracle/metric_interface.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>

namespace sirius::oracle {

//==============================================================================
// Analytic Boyer-Lindquist metric derivatives
//==============================================================================

// dg[sigma][mu][nu] = d_sigma g_munu for Kerr; stationary and axisymmetric, so
// only sigma = r, theta are non-zero. These expressions are the single source
// consumed by dHdq and by the connection below; their finite-difference gates
// are independent of either consumer. Reference: Misner, Thorne & Wheeler
// (1973), chapter 33.
inline void KerrMetricDerivatives(double M, double a, const Vec4d& x, double dg[4][4][4]) {
    for (int s = 0; s < 4; ++s)
        for (int mu = 0; mu < 4; ++mu)
            for (int nu = 0; nu < 4; ++nu) dg[s][mu][nu] = 0.0;

    double r = x.r;
    double theta = x.theta;
    double sinth = std::sin(theta);
    double costh = std::cos(theta);
    if (std::abs(sinth) < 1e-10) sinth = std::copysign(1e-10, sinth);
    double sin2th = sinth * sinth;
    double cos2th = costh * costh;
    double sin2theta = 2.0 * sinth * costh;  // sin(2 theta)

    double a2 = a * a;
    double r2 = r * r;
    double Sigma = r2 + a2 * cos2th;
    double Delta = r2 - 2.0 * M * r + a2;
    double A = (r2 + a2) * (r2 + a2) - a2 * Delta * sin2th;
    double Sigma2 = Sigma * Sigma;
    double Delta2 = Delta * Delta;

    double dSigma_dr = 2.0 * r;
    double dSigma_dth = -a2 * sin2theta;
    double dDelta_dr = 2.0 * (r - M);
    double dA_dr = 4.0 * r * (r2 + a2) - a2 * dDelta_dr * sin2th;
    double dA_dth = -a2 * Delta * sin2theta;

    // d/dr g_munu.
    dg[1][0][0] = 2.0 * M * (Sigma - 2.0 * r2) / Sigma2;
    dg[1][0][3] = dg[1][3][0] = -2.0 * M * a * sin2th * (Sigma - 2.0 * r2) / Sigma2;
    dg[1][1][1] = (2.0 * r * Delta - Sigma * dDelta_dr) / Delta2;
    dg[1][2][2] = dSigma_dr;
    dg[1][3][3] = sin2th * (dA_dr * Sigma - A * dSigma_dr) / Sigma2;

    // d/dtheta g_munu.
    dg[2][0][0] = -2.0 * M * r * dSigma_dth / Sigma2;
    dg[2][0][3] = dg[2][3][0] = -2.0 * M * a * r * sin2theta * (Sigma + a2 * sin2th) / Sigma2;
    dg[2][1][1] = dSigma_dth / Delta;
    dg[2][2][2] = dSigma_dth;
    dg[2][3][3] = (dA_dth * sin2th + A * sin2theta) / Sigma - A * sin2th * dSigma_dth / Sigma2;
}

// ddg[s1][s2][mu][nu] = d_s1 d_s2 g_munu for Kerr; non-zero only for
// s1, s2 in {r, theta}, symmetric in (s1, s2). Differentiated by hand from
// the first derivatives above and pinned against central differences of
// KerrMetricDerivatives by KerrMetricDTest.SecondDerivativesMatchFiniteDifference;
// together with the first derivatives this is the complete analytic input the
// Riemann assembly below contracts. Pole clamping matches KerrMetricDerivatives
// so the two stay consistent term by term.
inline void KerrMetricSecondDerivatives(double M, double a, const Vec4d& x,
                                        double ddg[4][4][4][4]) {
    for (int s1 = 0; s1 < 4; ++s1)
        for (int s2 = 0; s2 < 4; ++s2)
            for (int mu = 0; mu < 4; ++mu)
                for (int nu = 0; nu < 4; ++nu) ddg[s1][s2][mu][nu] = 0.0;

    double r = x.r;
    double theta = x.theta;
    double sinth = std::sin(theta);
    double costh = std::cos(theta);
    if (std::abs(sinth) < 1e-10) sinth = std::copysign(1e-10, sinth);
    double sin2th = sinth * sinth;
    double cos2th = costh * costh;
    double sin2theta = 2.0 * sinth * costh;  // sin(2 theta)
    double cos2theta = cos2th - sin2th;      // cos(2 theta)

    double a2 = a * a;
    double r2 = r * r;
    double Sigma = r2 + a2 * cos2th;
    double Delta = r2 - 2.0 * M * r + a2;
    double A = (r2 + a2) * (r2 + a2) - a2 * Delta * sin2th;
    double Sigma2 = Sigma * Sigma;
    double Sigma3 = Sigma2 * Sigma;
    double Delta2 = Delta * Delta;
    double Delta3 = Delta2 * Delta;

    double dSigma_dr = 2.0 * r;
    double ddSigma_drr = 2.0;
    double dSigma_dth = -a2 * sin2theta;
    double ddSigma_dthth = -2.0 * a2 * cos2theta;
    double dDelta_dr = 2.0 * (r - M);
    double ddDelta_drr = 2.0;
    double dA_dr = 4.0 * r * (r2 + a2) - a2 * dDelta_dr * sin2th;
    double ddA_drr = 4.0 * (3.0 * r2 + a2) - 2.0 * a2 * sin2th;
    double dA_dth = -a2 * Delta * sin2theta;
    double ddA_drth = -a2 * dDelta_dr * sin2theta;
    double ddA_dthth = -2.0 * a2 * Delta * cos2theta;

    // Second derivatives of the shared profile f = 2Mr/Sigma, from which
    // g_tt = -1 + f and g_tphi = -a sin^2(theta) f both descend.
    double ddf_drr = 4.0 * M * r * (4.0 * r2 - 3.0 * Sigma) / Sigma3;
    double ddf_drth = 2.0 * M * dSigma_dth * (4.0 * r2 - Sigma) / Sigma3;
    double ddf_dthth =
        -2.0 * M * r * (ddSigma_dthth * Sigma - 2.0 * dSigma_dth * dSigma_dth) / Sigma3;
    double df_dr = 2.0 * M * (Sigma - 2.0 * r2) / Sigma2;

    // d^2/dr^2 g_munu.
    ddg[1][1][0][0] = ddf_drr;
    ddg[1][1][0][3] = ddg[1][1][3][0] = -a * sin2th * ddf_drr;
    ddg[1][1][1][1] = ((ddSigma_drr * Delta - Sigma * ddDelta_drr) * Delta -
                       2.0 * dDelta_dr * (dSigma_dr * Delta - Sigma * dDelta_dr)) /
                      Delta3;
    ddg[1][1][2][2] = ddSigma_drr;
    ddg[1][1][3][3] = sin2th *
                      ((ddA_drr * Sigma - A * ddSigma_drr) * Sigma -
                       2.0 * dSigma_dr * (dA_dr * Sigma - A * dSigma_dr)) /
                      Sigma3;

    // d^2/dr dtheta g_munu (symmetric in the derivative pair).
    ddg[1][2][0][0] = ddf_drth;
    ddg[1][2][0][3] = ddg[1][2][3][0] = -a * (sin2theta * df_dr + sin2th * ddf_drth);
    ddg[1][2][1][1] = -dSigma_dth * dDelta_dr / Delta2;
    ddg[1][2][2][2] = 0.0;
    ddg[1][2][3][3] = sin2theta * (dA_dr * Sigma - A * dSigma_dr) / Sigma2 +
                      sin2th *
                          ((ddA_drth * Sigma + dA_dr * dSigma_dth - dA_dth * dSigma_dr) * Sigma -
                           2.0 * dSigma_dth * (dA_dr * Sigma - A * dSigma_dr)) /
                          Sigma3;
    for (int mu = 0; mu < 4; ++mu)
        for (int nu = 0; nu < 4; ++nu) ddg[2][1][mu][nu] = ddg[1][2][mu][nu];

    // d^2/dtheta^2 g_munu.
    double N_tphi = sin2theta * Sigma - sin2th * dSigma_dth;
    double dN_tphi_dth = 2.0 * cos2theta * Sigma - sin2th * ddSigma_dthth;
    ddg[2][2][0][0] = ddf_dthth;
    ddg[2][2][0][3] = ddg[2][2][3][0] =
        -2.0 * M * a * r * (dN_tphi_dth * Sigma - 2.0 * dSigma_dth * N_tphi) / Sigma3;
    ddg[2][2][1][1] = ddSigma_dthth / Delta;
    ddg[2][2][2][2] = ddSigma_dthth;
    ddg[2][2][3][3] =
        (ddA_dthth * sin2th + 2.0 * dA_dth * sin2theta + 2.0 * A * cos2theta) / Sigma -
        2.0 * (dA_dth * sin2th + A * sin2theta) * dSigma_dth / Sigma2 -
        A * sin2th * ddSigma_dthth / Sigma2 + 2.0 * A * sin2th * dSigma_dth * dSigma_dth / Sigma3;
}

//==============================================================================
// KerrMetricD: Double-precision Kerr metric implementation
//==============================================================================

class KerrMetricD : public IMetricD {
  public:
    explicit KerrMetricD(const MetricParamsD& params)
        : mass_(params.M),
          spin_(params.a),
          charge_(params.Q),
          r_plus_(params.rplus),
          r_minus_(params.rminus) {
        // Clamp spin to avoid naked singularity
        if (std::abs(spin_) > 0.9999 * mass_) {
            spin_ = std::copysign(0.9999 * mass_, spin_);
        }
    }

    KerrMetricD(double M = 1.0, double a = 0.0, double Q = 0.0)
        : KerrMetricD(MetricParamsD(M, a, Q)) {}

    //--------------------------------------------------------------------------
    // IMetricD Implementation
    //--------------------------------------------------------------------------

    void Evaluate(const Vec4d& x, double g[4][4], double g_inv[4][4]) const override {
        double r = x.r;
        double theta = x.theta;

        // Boyer-Lindquist is singular at the horizon. The oracle caller must
        // remain in its valid domain; it may not acquire flat-space evidence.
        SIRIUS_PRE(r > r_plus_);

        // Handle pole singularities
        double sinth = std::sin(theta);
        double costh = std::cos(theta);
        if (std::abs(sinth) < 1e-10) {
            sinth = std::copysign(1e-10, sinth);
        }
        double sin2th = sinth * sinth;
        double cos2th = costh * costh;

        // Kerr metric functions
        double Sigma = r * r + spin_ * spin_ * cos2th;
        double Delta = r * r - 2 * mass_ * r + spin_ * spin_;
        double A =
            (r * r + spin_ * spin_) * (r * r + spin_ * spin_) - spin_ * spin_ * Delta * sin2th;

        // Initialize to Zero
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j) g[i][j] = g_inv[i][j] = 0.0;

        // Covariant metric g_μν
        g[0][0] = -(1 - 2 * mass_ * r / Sigma);                       // g_tt
        g[0][3] = g[3][0] = -2 * mass_ * spin_ * r * sin2th / Sigma;  // g_tφ
        g[1][1] = Sigma / Delta;                                      // g_rr
        g[2][2] = Sigma;                                              // g_θθ
        g[3][3] = A * sin2th / Sigma;                                 // g_φφ

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

    void Christoffel(const Vec4d& x, double Gamma[4][4][4]) const override {
        // Gamma^mu_nu_rho = (1/2) g^mu_sigma (d_nu g_sigma_rho + d_rho g_sigma_nu
        // - d_sigma g_nu_rho), contracted from the exact inverse metric and the
        // analytic derivatives above. Replaces hand-expanded per-component
        // formulas that disagreed with finite differences of this class's own
        // metric at O(1) (worst: Gamma^phi_theta_phi lost its cot(theta) term);
        // the defect was latent because the integration path uses dHdq. Pinned
        // by KerrMetricDTest.ChristoffelMatchesFiniteDifferencesOfMetric.
        double g[4][4], g_inv[4][4];
        Evaluate(x, g, g_inv);
        double dg[4][4][4];
        KerrMetricDerivatives(mass_, spin_, x, dg);

        for (int mu = 0; mu < 4; ++mu) {
            for (int nu = 0; nu < 4; ++nu) {
                for (int rho = 0; rho < 4; ++rho) {
                    double s = 0.0;
                    for (int sigma = 0; sigma < 4; ++sigma) {
                        s += g_inv[mu][sigma] *
                             (dg[nu][sigma][rho] + dg[rho][sigma][nu] - dg[sigma][nu][rho]);
                    }
                    Gamma[mu][nu][rho] = 0.5 * s;
                }
            }
        }
    }

    void Riemann(const Vec4d& x, double R[4][4][4][4]) const override {
        // R^rho_smn = d_m Gamma^rho_ns - d_n Gamma^rho_ms
        //           + Gamma^rho_ml Gamma^l_ns - Gamma^rho_nl Gamma^l_ms,
        // with the connection derivative expanded through the product rule:
        //   d_s Gamma^m_nr = 1/2 (d_s g^ml)(d_n g_lr + d_r g_ln - d_l g_nr)
        //                  + 1/2 g^ml (d_s d_n g_lr + d_s d_r g_ln - d_s d_l g_nr),
        //   d_s g^ml = -g^ma g^lb d_s g_ab.
        // Contracted from the same analytic derivative source as Christoffel,
        // extended one order deeper by KerrMetricSecondDerivatives. Replaces a
        // hand-expanded component table that was incomplete off the diagonal
        // blocks (the frame-dragging t-phi sector was partial and several
        // components disagreed with finite differences of this class's own
        // connection); the beam integrator was the only consumer and its gates
        // were structural, so the defect stayed latent. Pinned by
        // KerrMetricD_RiemannMatchesFiniteDifferenceChristoffel and the
        // vacuum-Ricci, lowered-symmetry, and Kretschmann-contraction gates.
        double g[4][4], g_inv[4][4];
        Evaluate(x, g, g_inv);
        double dg[4][4][4];
        KerrMetricDerivatives(mass_, spin_, x, dg);
        double ddg[4][4][4][4];
        KerrMetricSecondDerivatives(mass_, spin_, x, ddg);

        double Gamma[4][4][4];
        Christoffel(x, Gamma);

        // d_s g^ab; the metric is stationary and axisymmetric, so only
        // s = r, theta contribute.
        double dginv[4][4][4];
        for (int s = 0; s < 4; ++s) {
            for (int m = 0; m < 4; ++m) {
                for (int n = 0; n < 4; ++n) {
                    double v = 0.0;
                    if (s == 1 || s == 2) {
                        for (int p = 0; p < 4; ++p)
                            for (int q = 0; q < 4; ++q)
                                v -= g_inv[m][p] * g_inv[n][q] * dg[s][p][q];
                    }
                    dginv[s][m][n] = v;
                }
            }
        }

        // dGamma[s][m][n][r] = d_s Gamma^m_nr.
        double dGamma[4][4][4][4];
        for (int s = 0; s < 4; ++s) {
            for (int m = 0; m < 4; ++m) {
                for (int n = 0; n < 4; ++n) {
                    for (int r = 0; r < 4; ++r) {
                        double v = 0.0;
                        if (s == 1 || s == 2) {
                            for (int l = 0; l < 4; ++l) {
                                v += dginv[s][m][l] * (dg[n][l][r] + dg[r][l][n] - dg[l][n][r]) +
                                     g_inv[m][l] *
                                         (ddg[s][n][l][r] + ddg[s][r][l][n] - ddg[s][l][n][r]);
                            }
                            v *= 0.5;
                        }
                        dGamma[s][m][n][r] = v;
                    }
                }
            }
        }

        for (int rho = 0; rho < 4; ++rho) {
            for (int sig = 0; sig < 4; ++sig) {
                for (int mu = 0; mu < 4; ++mu) {
                    for (int nu = 0; nu < 4; ++nu) {
                        double v = dGamma[mu][rho][nu][sig] - dGamma[nu][rho][mu][sig];
                        for (int lam = 0; lam < 4; ++lam) {
                            v += Gamma[rho][mu][lam] * Gamma[lam][nu][sig] -
                                 Gamma[rho][nu][lam] * Gamma[lam][mu][sig];
                        }
                        R[rho][sig][mu][nu] = v;
                    }
                }
            }
        }
    }

    /// @brief Compute Kretschmann scalar K = R_αβγδ R^αβγδ (analytic)
    /// For validation of Riemann tensor
    ///
    /// FORMULA:
    /// For Kerr spacetime in Boyer-Lindquist coordinates:
    ///   K = 48M² × (r² - a²cos²θ) × [(r² + a²cos²θ)² - 16a²r²cos²θ] / Σ⁶
    ///
    /// where Σ = r² + a²cos²θ. The bracket is Σ² - 16a²r²cos²θ; an earlier
    /// version squared (r² - a²cos²θ) there instead, which coincides with the
    /// correct value only where a·cosθ = 0 (Schwarzschild, or the Kerr
    /// equator) — exactly the two regimes the historical gates sampled. The
    /// contraction gate KerrMetricD_KretschmannMatchesRiemannContraction now
    /// pins this form against the assembled tensor off-equator.
    ///
    /// Special cases:
    /// - Schwarzschild (a=0): K = 48M²/r⁶
    /// - Extremal Kerr (a=M): K finite everywhere except ring singularity
    ///
    /// Reference: Henry, R.C. (2000), "Kretschmann Scalar for a Kerr-Newman Black Hole"
    ///            Astrophys. J. 535:350-353, equation (18) at Q = 0.
    double Kretschmann(const Vec4d& x) const {
        double r = x.r;
        double theta = x.theta;
        double M = mass_;
        double a = spin_;

        double cos2th = std::cos(theta) * std::cos(theta);
        double r2 = r * r;
        double a2 = a * a;

        double Sigma = r2 + a2 * cos2th;
        double Sigma6 = std::pow(Sigma, 6);

        // Protect against singularity
        if (Sigma6 < 1e-60) {
            return std::numeric_limits<double>::infinity();
        }

        // r² - a²cos²θ (appears in numerator)
        double r_term = r2 - a2 * cos2th;

        // K = 48M² × r_term × [Σ² - 16a²r²cos²θ] / Σ⁶ (Henry 2000, eq. 18).
        double bracket = Sigma * Sigma - 16.0 * a2 * r2 * cos2th;

        return 48.0 * M * M * r_term * bracket / Sigma6;
    }

    /// @brief Verify Riemann tensor symmetries (diagnostic)
    /// Returns maximum violation of antisymmetry R^α_βγδ = -R^α_βδγ
    double VerifyRiemannSymmetries(const Vec4d& x) const {
        double R[4][4][4][4];
        Riemann(x, R);

        double max_violation = 0;
        for (int a = 0; a < 4; ++a) {
            for (int b = 0; b < 4; ++b) {
                for (int c = 0; c < 4; ++c) {
                    for (int d = 0; d < 4; ++d) {
                        // Antisymmetry in last two indices
                        double violation = std::abs(R[a][b][c][d] + R[a][b][d][c]);
                        max_violation = std::max(max_violation, violation);
                    }
                }
            }
        }
        return max_violation;
    }

    double Hamiltonian(const Vec4d& q, const Vec4d& p) const override {
        double g_inv[4][4], g[4][4];
        Evaluate(q, g, g_inv);

        // H = (1/2) g^μν p_μ p_ν
        double H = 0;
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j) H += g_inv[i][j] * p[i] * p[j];
        return 0.5 * H;
    }

    Vec4d dHdp(const Vec4d& q, const Vec4d& p) const override {
        // ∂H/∂p_μ = g^μν p_ν
        double g_inv[4][4], g[4][4];
        Evaluate(q, g, g_inv);

        Vec4d result;
        for (int i = 0; i < 4; ++i) {
            result[i] = 0;
            for (int j = 0; j < 4; ++j) result[i] += g_inv[i][j] * p[j];
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
        double M = mass_;
        double a = spin_;
        double a2 = a * a;
        double r2 = r * r;

        double sinth = std::sin(theta);
        double costh = std::cos(theta);
        if (std::abs(sinth) < 1e-10) sinth = std::copysign(1e-10, sinth);
        double sin2th = sinth * sinth;
        double cos2th = costh * costh;
        double sin2theta = 2 * sinth * costh;  // sin(2θ)

        // Kerr metric auxiliary functions
        double Sigma = r2 + a2 * cos2th;
        double Delta = r2 - 2 * M * r + a2;
        double A = (r2 + a2) * (r2 + a2) - a2 * Delta * sin2th;

        double Sigma2 = Sigma * Sigma;
        double Delta2 = Delta * Delta;

        // Derivatives of auxiliary functions
        double dSigma_dr = 2 * r;
        double dSigma_dth = -a2 * sin2theta;
        double dDelta_dr = 2 * (r - M);
        double dA_dr = 4 * r * (r2 + a2) - a2 * dDelta_dr * sin2th;
        double dA_dth = -a2 * Delta * sin2theta;

        // ============================================================
        // Covariant metric derivatives ∂g_αβ/∂r
        // ============================================================

        // g_tt = -(1 - 2Mr/Σ) = -1 + 2Mr/Σ
        // ∂g_tt/∂r = 2M(Σ - r·dΣ/dr)/Σ² = 2M(Σ - 2r²)/Σ²
        double dg_tt_dr = 2 * M * (Sigma - 2 * r2) / Sigma2;

        // g_tφ = -2Mar·sin²θ/Σ
        // ∂g_tφ/∂r = -2Ma·sin²θ·(Σ - r·dΣ/dr)/Σ² = -2Ma·sin²θ·(Σ - 2r²)/Σ²
        double dg_tphi_dr = -2 * M * a * sin2th * (Sigma - 2 * r2) / Sigma2;

        // g_rr = Σ/Δ
        // ∂g_rr/∂r = (dΣ/dr·Δ - Σ·dΔ/dr)/Δ² = (2r·Δ - Σ·2(r-M))/Δ²
        double dg_rr_dr = (2 * r * Delta - Sigma * dDelta_dr) / Delta2;

        // g_θθ = Σ
        // ∂g_θθ/∂r = dΣ/dr = 2r
        double dg_thth_dr = dSigma_dr;

        // g_φφ = A·sin²θ/Σ
        // ∂g_φφ/∂r = sin²θ·(dA/dr·Σ - A·dΣ/dr)/Σ²
        double dg_phiphi_dr = sin2th * (dA_dr * Sigma - A * dSigma_dr) / Sigma2;

        // ============================================================
        // Covariant metric derivatives ∂g_αβ/∂θ
        // ============================================================

        // ∂g_tt/∂θ = -2Mr·dΣ/dθ/Σ² = 2Mra²sin2θ/Σ²
        double dg_tt_dth = -2 * M * r * dSigma_dth / Sigma2;

        // ∂g_tφ/∂θ = -2Mar·[2sinθcosθ/Σ - sin²θ·dΣ/dθ/Σ²]
        //          = -2Mar·[sin2θ/Σ + sin²θ·a²sin2θ/Σ²]
        //          = -2Mar·sin2θ·(Σ + a²sin²θ)/Σ²
        double dg_tphi_dth = -2 * M * a * r * sin2theta * (Sigma + a2 * sin2th) / Sigma2;

        // ∂g_rr/∂θ = (dΣ/dθ)/Δ = (-a²sin2θ)/Δ
        double dg_rr_dth = dSigma_dth / Delta;

        // ∂g_θθ/∂θ = dΣ/dθ = -a²sin2θ
        double dg_thth_dth = dSigma_dth;

        // ∂g_φφ/∂θ = ∂(A·sin²θ/Σ)/∂θ
        //          = (dA/dθ·sin²θ + A·sin2θ)/Σ - A·sin²θ·dΣ/dθ/Σ²
        double dg_phiphi_dth =
            (dA_dth * sin2th + A * sin2theta) / Sigma - A * sin2th * dSigma_dth / Sigma2;

        // ============================================================
        // Get inverse metric and compute contravariant velocity v^α = g^αβ p_β
        // ============================================================

        double g[4][4], g_inv[4][4];
        Evaluate(q, g, g_inv);

        double vt = g_inv[0][0] * p.t + g_inv[0][3] * p.phi;
        double vr = g_inv[1][1] * p.r;
        double vth = g_inv[2][2] * p.theta;
        double vphi = g_inv[3][0] * p.t + g_inv[3][3] * p.phi;

        // ============================================================
        // Compute dp_μ/dλ = (1/2)(∂g_αβ/∂q^μ) v^α v^β
        // ============================================================

        // For μ = r:
        double dp_r = 0.5 * (dg_tt_dr * vt * vt + 2 * dg_tphi_dr * vt * vphi + dg_rr_dr * vr * vr +
                             dg_thth_dr * vth * vth + dg_phiphi_dr * vphi * vphi);

        // For μ = θ:
        double dp_th =
            0.5 * (dg_tt_dth * vt * vt + 2 * dg_tphi_dth * vt * vphi + dg_rr_dth * vr * vr +
                   dg_thth_dth * vth * vth + dg_phiphi_dth * vphi * vphi);

        // Return dp/dλ
        Vec4d result;
        result.t = 0;  // dp_t/dλ = 0 (∂g/∂t = 0, stationary)
        result.r = dp_r;
        result.theta = dp_th;
        result.phi = 0;  // dp_φ/dλ = 0 (∂g/∂φ = 0, axisymmetric)

        return result;
    }

    bool IsValid(const Vec4d& x) const override {
        // Check radial coordinate
        if (x.r <= r_plus_ * 1.001) return false;  // Inside or at horizon
        if (x.r > 1e6 * mass_) return false;       // Too far (numerical issues)

        // Check angular coordinate
        if (x.theta <= 1e-6 || x.theta >= std::numbers::pi - 1e-6) return false;  // At poles

        return true;
    }

    double HorizonRadius() const override { return r_plus_; }
    double InnerHorizonRadius() const override { return r_minus_; }

    double ErgosphereRadius(double theta) const override {
        // r_ergo = M + √(M² - a²cos²θ)
        double cos2th = std::cos(theta) * std::cos(theta);
        return mass_ + std::sqrt(mass_ * mass_ - spin_ * spin_ * cos2th);
    }

    double IscoRadius() const override {
        // ISCO for prograde orbit
        if (std::abs(spin_) < 1e-10) return 6.0 * mass_;  // Schwarzschild

        double a_star = spin_ / mass_;
        double Z1 =
            1 + std::cbrt(1 - a_star * a_star) * (std::cbrt(1 + a_star) + std::cbrt(1 - a_star));
        double Z2 = std::sqrt(3 * a_star * a_star + Z1 * Z1);
        return mass_ * (3 + Z2 - std::sqrt((3 - Z1) * (3 + Z1 + 2 * Z2)));
    }

    double PhotonSphereRadius() const override {
        // Prograde photon sphere for Kerr
        if (std::abs(spin_) < 1e-10) return 3.0 * mass_;  // Schwarzschild

        double a_star = spin_ / mass_;
        return 2 * mass_ * (1 + std::cos(2.0 / 3.0 * std::acos(-a_star)));
    }

    double mass() const override { return mass_; }
    double spin() const override { return spin_; }
    double charge() const override { return charge_; }

    double TimeTransformationFunction(const Vec4d& q) const override {
        // g(q) = Σ for TTESI regularisation
        double cos2th = std::cos(q.theta) * std::cos(q.theta);
        return q.r * q.r + spin_ * spin_ * cos2th;
    }

    double GFactor(const Vec4d& x, const Vec4d& k, const Vec4d& u_emitter) const override {
        // g = (k · u_obs) / (k · u_emit)
        // For static observer: u_obs = (1/α, 0, 0, ω/α) where α = lapse, ω = frame-dragging

        double g[4][4], g_inv[4][4];
        Evaluate(x, g, g_inv);

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
        double k_dot_u_obs = k.t * u_obs.t + k.phi * u_obs.phi;  // Only t, φ non-Zero
        double k_dot_u_emit = 0;
        for (int i = 0; i < 4; ++i) k_dot_u_emit += k[i] * u_emitter[i];

        if (std::abs(k_dot_u_emit) < 1e-20) return 1.0;
        return k_dot_u_obs / k_dot_u_emit;
    }

  private:
    double mass_;     // Mass
    double spin_;     // Spin
    double charge_;   // Charge (for future Kerr-Newman extension)
    double r_plus_;   // Outer horizon
    double r_minus_;  // Inner horizon
};

}  // namespace sirius::oracle
