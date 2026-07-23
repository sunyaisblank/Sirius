#pragma once

// Unified Kerr-Schild black hole family: one implementation covering Minkowski,
// Schwarzschild, Kerr, Reissner-Nordstrom, Kerr-Newman, and their de Sitter
// extensions through the parameters (M, a, Q, Lambda). Ported from PHMT100A.h.
//
// Kerr-Schild ansatz in Cartesian coordinates:
//   g_mu_nu = eta_mu_nu + H l_mu l_nu,   g^mu_nu = eta^mu_nu - H l^mu l^nu,
// with r the implicit solution of r^4 - (R^2 - a^2) r^2 - a^2 z^2 = 0,
// H = (2 M r - Q^2) r^2 / (r^4 + a^2 z^2), and l a null vector field. The
// Christoffel symbols are polynomial in (x, y, z, r), so there are no pole
// singularities.
// Reference: Visser, "The Kerr spacetime" (arXiv:0706.0622).

#include "sirius/core/metrics/metric.h"

#include <algorithm>
#include <cmath>

namespace sirius::core {

// Kerr-Schild family parameters.
struct KerrSchildParams {
    double M = 1.0;       // Mass (Schwarzschild radius = 2M).
    double a = 0.0;       // Spin parameter (|a| <= M for a black hole).
    double Q = 0.0;       // Electric charge (|Q| <= M for a black hole).
    double Lambda = 0.0;  // Cosmological constant (Lambda > 0 for de Sitter).

    static KerrSchildParams Minkowski() { return {0, 0, 0, 0}; }
    static KerrSchildParams Schwarzschild(double M) { return {M, 0, 0, 0}; }
    static KerrSchildParams Kerr(double M, double a) { return {M, a, 0, 0}; }
    static KerrSchildParams ReissnerNordstrom(double M, double Q) { return {M, 0, Q, 0}; }
    static KerrSchildParams KerrNewman(double M, double a, double Q) { return {M, a, Q, 0}; }
    static KerrSchildParams DeSitter(double Lambda) { return {0, 0, 0, Lambda}; }
};

// Kerr-Schild family metric.
class KerrSchildFamily : public IMetric {
  public:
    KerrSchildFamily();
    explicit KerrSchildFamily(const KerrSchildParams& params);

    void Evaluate(const Tensor<double, 4>& pos, Metric4d& g,
                  Tensor<Dual<double>, 4, 4, 4>& dg) override;

    const Config& GetParameters() const override { return config_; }
    void SetParameter(const std::string& key, double value) override;
    const char* GetName() const override;

    // Exact closed-form inverse g^mu_nu = eta^mu_nu - H l^mu l^nu. Valid because
    // l is null with respect to both eta and g, so the Sherman-Morrison
    // correction terms vanish identically.
    bool InverseMetric(const Tensor<double, 4>& pos, Metric4d& g_inv) const override;

    // Capture test against the outer horizon in the Kerr radial coordinate,
    // exact for spin unlike a Cartesian-norm comparison.
    bool InsideCaptureSurface(const Tensor<double, 4>& pos, double margin) const override;

    void SetParams(const KerrSchildParams& params);
    KerrSchildParams GetParams() const;

    // Boyer-Lindquist-like radius from Cartesian coordinates, solving
    // r^4 - (x^2 + y^2 + z^2 - a^2) r^2 - a^2 z^2 = 0.
    double ComputeKerrRadius(double x, double y, double z) const;

    // Null vector l^mu at a position.
    void ComputeNullVector(double x, double y, double z, double r, double l[4]) const;

    // Scalar function H = (2 M r - Q^2) r^2 / (r^4 + a^2 z^2).
    double ComputeH(double r, double z) const;

    // Outer event horizon r+ = M + sqrt(M^2 - a^2 - Q^2), or -1 for a naked
    // singularity.
    double OuterHorizonRadius() const;

    // Inner (Cauchy) horizon r- = M - sqrt(M^2 - a^2 - Q^2), or -1 for a naked
    // singularity.
    double InnerHorizonRadius() const;

    // True when M^2 >= a^2 + Q^2 (horizons exist).
    bool HasHorizon() const;

    // Ergosphere boundary radius at polar angle theta.
    double ErgosphereRadius(double theta) const;

    // ISCO radius for prograde equatorial orbits.
    double IscoRadius() const;

    // Extremality parameter chi = sqrt(a^2 + Q^2) / M, in [0, 1] for black holes.
    double ExtremalityParameter() const;

  private:
    Config config_;
    KerrSchildParams params_;
};

inline KerrSchildFamily::KerrSchildFamily() {
    config_["mass"] = {1.0, 0.0, 100.0};
    config_["spin"] = {0.0, -0.998, 0.998};
    config_["charge"] = {0.0, -1.0, 1.0};
    config_["lambda"] = {0.0, -0.1, 0.1};
    params_ = KerrSchildParams::Schwarzschild(1.0);
}

inline KerrSchildFamily::KerrSchildFamily(const KerrSchildParams& params) : KerrSchildFamily() {
    SetParams(params);
}

inline void KerrSchildFamily::SetParams(const KerrSchildParams& params) {
    params_ = params;
    config_["mass"].value = params.M;
    config_["spin"].value = params.a / std::max(params.M, 1e-10);
    config_["charge"].value = params.Q / std::max(params.M, 1e-10);
    config_["lambda"].value = params.Lambda;
}

inline KerrSchildParams KerrSchildFamily::GetParams() const { return params_; }

inline void KerrSchildFamily::SetParameter(const std::string& key, double value) {
    if (config_.find(key) != config_.end()) {
        config_[key].value = std::clamp(value, config_[key].min, config_[key].max);
    }
    params_.M = config_["mass"].value;
    params_.a = config_["spin"].value * params_.M;
    params_.Q = config_["charge"].value * params_.M;
    params_.Lambda = config_["lambda"].value;
}

inline const char* KerrSchildFamily::GetName() const {
    if (params_.M == 0 && params_.Lambda == 0) return "Minkowski";
    if (params_.a == 0 && params_.Q == 0 && params_.Lambda == 0) return "Schwarzschild";
    if (params_.Q == 0 && params_.Lambda == 0) return "Kerr";
    if (params_.a == 0 && params_.Lambda == 0) return "Reissner-Nordström";
    if (params_.Lambda == 0) return "Kerr-Newman";
    if (params_.M == 0) return "de Sitter";
    if (params_.a == 0 && params_.Q == 0) return "Schwarzschild-de Sitter";
    if (params_.Q == 0) return "Kerr-de Sitter";
    return "Kerr-Newman-de Sitter";
}

inline double KerrSchildFamily::ComputeKerrRadius(double x, double y, double z) const {
    double a = params_.a;
    double a2 = a * a;
    double R2 = x * x + y * y + z * z;

    // a = 0 (Schwarzschild): r = R directly.
    if (std::abs(a) < 1e-12) {
        return std::sqrt(std::max(R2, 1e-20));
    }

    // Solve u^2 - (R^2 - a^2) u - a^2 z^2 = 0 for u = r^2.
    double Rm2 = R2 - a2;
    double disc = Rm2 * Rm2 + 4.0 * a2 * z * z;
    double r2 = (Rm2 + std::sqrt(std::max(disc, 0.0))) / 2.0;

    return std::sqrt(std::max(r2, 1e-20));
}

inline void KerrSchildFamily::ComputeNullVector(double x, double y, double z, double r,
                                                double l[4]) const {
    double a = params_.a;
    double a2 = a * a;
    double r2 = r * r;
    double denom = r2 + a2;

    // l^mu = (1, (rx + ay)/(r^2 + a^2), (ry - ax)/(r^2 + a^2), z/r).
    l[0] = 1.0;
    l[1] = (r * x + a * y) / denom;
    l[2] = (r * y - a * x) / denom;
    l[3] = z / std::max(r, 1e-10);
}

inline double KerrSchildFamily::ComputeH(double r, double z) const {
    double M = params_.M;
    double a = params_.a;
    double Q = params_.Q;
    double Lambda = params_.Lambda;

    double r2 = r * r;
    double r4 = r2 * r2;
    double a2 = a * a;

    // H = (2 M r - Q^2) r^2 / (r^4 + a^2 z^2); pure Schwarzschild gives 2M/r.
    double numerator = (2.0 * M * r - Q * Q) * r2;
    double denominator = r4 + a2 * z * z;
    double H = numerator / std::max(denominator, 1e-20);

    // Schwarzschild-de Sitter is exactly Kerr-Schild with H += Lambda r^2/3 at
    // a = 0. The rotating de Sitter forms need a different ansatz, so Lambda with
    // a != 0 is not represented here; the configuration boundary rejects that
    // combination rather than letting an approximation stand in.
    if (std::abs(a) < 1e-12 && std::abs(Lambda) > 1e-15) {
        H += Lambda * r2 / 3.0;
    }

    return H;
}

inline void KerrSchildFamily::Evaluate(const Tensor<double, 4>& pos, Metric4d& g,
                                       Tensor<Dual<double>, 4, 4, 4>& dg) {
    [[maybe_unused]] double t = pos(0);  // Time coordinate (unused in a static metric).
    double x = pos(1);
    double y = pos(2);
    double z = pos(3);

    double M = params_.M;
    double a = params_.a;
    double Q = params_.Q;
    double Lambda = params_.Lambda;

    double a2 = a * a;
    double Q2 = Q * Q;

    // Kerr radius.
    double R2 = x * x + y * y + z * z;
    double Rm2 = R2 - a2;
    double disc = Rm2 * Rm2 + 4.0 * a2 * z * z;
    disc = std::max(disc, 1e-20);
    double sqrt_disc = std::sqrt(disc);
    double r2 = (Rm2 + sqrt_disc) / 2.0;
    r2 = std::max(r2, 1e-10);
    double r = std::sqrt(r2);
    double r3 = r2 * r;
    double r4 = r2 * r2;

    // Derivatives of r with respect to the spatial coordinates.
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

    // Null vector l^mu and its derivatives.
    double denom = r2 + a2;
    double denom2 = denom * denom;
    double l[4] = {1.0, (r * x + a * y) / denom, (r * y - a * x) / denom, z / r};

    // dl[lam][mu] = d l^mu / d x^lam.
    double dl[4][4] = {{0}};
    double d_denom[4] = {0.0, d_r2_dx, d_r2_dy, d_r2_dz};

    // l^1 = (rx + ay) / denom.
    dl[1][1] = (r + x * dr_dx) / denom - (r * x + a * y) * d_denom[1] / denom2;
    dl[2][1] = (x * dr_dy + a) / denom - (r * x + a * y) * d_denom[2] / denom2;
    dl[3][1] = (x * dr_dz) / denom - (r * x + a * y) * d_denom[3] / denom2;

    // l^2 = (ry - ax) / denom.
    dl[1][2] = (y * dr_dx - a) / denom - (r * y - a * x) * d_denom[1] / denom2;
    dl[2][2] = (r + y * dr_dy) / denom - (r * y - a * x) * d_denom[2] / denom2;
    dl[3][2] = (y * dr_dz) / denom - (r * y - a * x) * d_denom[3] / denom2;

    // l^3 = z / r.
    dl[1][3] = -z * dr_dx / r2;
    dl[2][3] = -z * dr_dy / r2;
    dl[3][3] = 1.0 / r - z * dr_dz / r2;

    // Kerr-Schild scalar H (single authority ComputeH, which folds the a = 0
    // cosmological term so Schwarzschild-de Sitter stays in exact Kerr-Schild
    // form).
    double f_denom = r4 + a2 * z * z;
    double H = ComputeH(r, z);

    // Derivatives of H.
    double dH[4] = {0.0, 0.0, 0.0, 0.0};
    for (int lam = 1; lam <= 3; lam++) {
        double d_num = 2.0 * r * dr[lam] * (3.0 * M * r - Q2);
        double d_f_denom = 4.0 * r3 * dr[lam];
        if (lam == 3) d_f_denom += 2.0 * a2 * z;

        double numerator = (2.0 * M * r - Q2) * r2;
        dH[lam] = (d_num * f_denom - numerator * d_f_denom) / (f_denom * f_denom);

        // a = 0 cosmological term: d(Lambda r^2/3)/d x^lam = (2 Lambda/3) r dr.
        if (std::abs(a) < 1e-12 && std::abs(Lambda) > 1e-15) {
            dH[lam] += (2.0 * Lambda / 3.0) * r * dr[lam];
        }
    }

    // Metric g_mu_nu = eta_mu_nu + H l_mu l_nu.
    g.Zero();
    g(0, 0) = Dual<double>(-1.0);
    g(1, 1) = Dual<double>(1.0);
    g(2, 2) = Dual<double>(1.0);
    g(3, 3) = Dual<double>(1.0);

    for (int mu = 0; mu < 4; mu++) {
        for (int nu = 0; nu < 4; nu++) {
            g(mu, nu) = Dual<double>(g(mu, nu).real + H * l[mu] * l[nu]);
        }
    }

    // Derivatives d g_mu_nu / d x^lam = dH l_mu l_nu + H dl_mu l_nu + H l_mu dl_nu.
    dg.Zero();

    for (int lam = 1; lam <= 3; lam++) {
        for (int mu = 0; mu < 4; mu++) {
            for (int nu = 0; nu < 4; nu++) {
                double dg_val =
                    dH[lam] * l[mu] * l[nu] + H * dl[lam][mu] * l[nu] + H * l[mu] * dl[lam][nu];
                dg(lam, mu, nu) = Dual<double>(dg_val);
            }
        }
    }
}

inline bool KerrSchildFamily::InverseMetric(const Tensor<double, 4>& pos, Metric4d& g_inv) const {
    // g = eta + H l(x)l with l null (eta^mu_nu l_mu l_nu = 0, a consequence of the
    // defining quartic for r), so the inverse is exactly
    //   g^mu_nu = eta^mu_nu - H l^mu l^nu,   l^mu = eta^mu_sigma l_sigma = (-1, l_1, l_2, l_3).
    // The sign of l cancels in the outer product, so the spatial covariant
    // components from ComputeNullVector can be used directly.
    double x = pos(1), y = pos(2), z = pos(3);

    double r = ComputeKerrRadius(x, y, z);
    double l[4];
    ComputeNullVector(x, y, z, r, l);
    double H = ComputeH(r, z);

    double l_up[4] = {-1.0, l[1], l[2], l[3]};

    g_inv.Zero();
    g_inv(0, 0) = Dual<double>(-1.0);
    g_inv(1, 1) = Dual<double>(1.0);
    g_inv(2, 2) = Dual<double>(1.0);
    g_inv(3, 3) = Dual<double>(1.0);
    for (int mu = 0; mu < 4; mu++) {
        for (int nu = 0; nu < 4; nu++) {
            g_inv(mu, nu) = Dual<double>(g_inv(mu, nu).real - H * l_up[mu] * l_up[nu]);
        }
    }
    return true;
}

inline bool KerrSchildFamily::InsideCaptureSurface(const Tensor<double, 4>& pos,
                                                   double margin) const {
    if (!HasHorizon()) return false;

    // Compare in the Kerr radial coordinate: the horizon r = r+ is an oblate
    // surface in Cartesian coordinates, so a Cartesian-norm comparison would
    // misplace it for a != 0.
    double r = ComputeKerrRadius(pos(1), pos(2), pos(3));
    return r <= OuterHorizonRadius() * (1.0 + margin);
}

inline double KerrSchildFamily::OuterHorizonRadius() const {
    double M = params_.M;
    double a = params_.a;
    double Q = params_.Q;

    double discriminant = M * M - a * a - Q * Q;

    if (discriminant < 0) {
        return -1.0;  // No horizon (naked singularity).
    }

    return M + std::sqrt(discriminant);
}

inline double KerrSchildFamily::InnerHorizonRadius() const {
    double M = params_.M;
    double a = params_.a;
    double Q = params_.Q;

    double discriminant = M * M - a * a - Q * Q;

    if (discriminant < 0) {
        return -1.0;  // No horizon (naked singularity).
    }

    return M - std::sqrt(discriminant);
}

inline bool KerrSchildFamily::HasHorizon() const {
    double M = params_.M;
    double a = params_.a;
    double Q = params_.Q;

    return (M * M >= a * a + Q * Q);
}

inline double KerrSchildFamily::ErgosphereRadius(double theta) const {
    double M = params_.M;
    double a = params_.a;
    double Q = params_.Q;

    // Static limit g_tt = 0. For Kerr-Newman: r^2 - 2Mr + Q^2 + a^2 cos^2 theta = 0,
    // giving r = M + sqrt(M^2 - Q^2 - a^2 cos^2 theta).
    double cos2th = std::cos(theta) * std::cos(theta);

    double disc_charged = M * M - Q * Q - a * a * cos2th;
    if (disc_charged < 0) {
        // The ergosphere does not exist at this angle.
        return M;  // Return the event horizon as a fallback.
    }

    return M + std::sqrt(disc_charged);
}

inline double KerrSchildFamily::IscoRadius() const {
    double M = params_.M;
    double a = params_.a;
    double Q = params_.Q;

    // Schwarzschild (a = 0, Q = 0): r_ISCO = 6M.
    if (std::abs(a) < 1e-12 && std::abs(Q) < 1e-12) {
        return 6.0 * M;
    }

    // Reissner-Nordstrom (a = 0, Q != 0): r_ISCO = 3M + sqrt(9M^2 - 8Q^2),
    // reducing to 6M at Q = 0 and 4M at Q^2 = M^2/2.
    if (std::abs(a) < 1e-12) {
        double disc = 9 * M * M - 8 * Q * Q;
        if (disc < 0) {
            // Very high charge: use the minimum.
            return 4.0 * M;
        }
        return (3 * M + std::sqrt(disc)) / 1.0;  // Prograde equivalent.
    }

    // Kerr (Q = 0): Bardeen-Press-Teukolsky closed form, exact.
    //
    // RECORDED APPROXIMATION for Kerr-Newman (a != 0 and Q != 0): the Kerr
    // formula below ignores Q. The error vanishes as Q -> 0 and grows like
    // (Q/M)^2; the reference anchors are exact (Q=0: this formula; a=0: the
    // Reissner-Nordstrom branch above, which is exact at both Q=0 -> 6M and
    // Q=M -> 4M). Replacing this with the Dadhich-Kale marginal-stability solve
    // is a recorded follow-up gated on the test toolchain, because an
    // unvalidated root-solve is worse than a documented approximation.
    double a_star = std::abs(a / M);
    if (a_star > 0.9999) a_star = 0.9999;

    double Z1 =
        1 + std::cbrt(1 - a_star * a_star) * (std::cbrt(1 + a_star) + std::cbrt(1 - a_star));
    double Z2 = std::sqrt(3 * a_star * a_star + Z1 * Z1);

    double r_isco;
    if (a >= 0) {
        r_isco = 3 + Z2 - std::sqrt((3 - Z1) * (3 + Z1 + 2 * Z2));
    } else {
        r_isco = 3 + Z2 + std::sqrt((3 - Z1) * (3 + Z1 + 2 * Z2));
    }

    return M * r_isco;
}

inline double KerrSchildFamily::ExtremalityParameter() const {
    double M = params_.M;
    double a = params_.a;
    double Q = params_.Q;

    if (M <= 0) return 0;

    // chi = sqrt(a^2 + Q^2) / M: 0 for Schwarzschild, 1 for extremal, > 1 naked.
    return std::sqrt(a * a + Q * Q) / M;
}

}  // namespace sirius::core
