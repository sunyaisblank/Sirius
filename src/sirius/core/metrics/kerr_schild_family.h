#pragma once

// Unified Kerr-Schild metric family: one implementation covering Minkowski,
// Schwarzschild, Kerr, Reissner-Nordstrom, Kerr-Newman, de Sitter, and the
// spherical Schwarzschild-de Sitter (Kottler) sector through (M, a, Q, Lambda).
// Rotating or charged cosmological sectors are not represented. Ported from
// PHMT100A.h.
//
// Kerr-Schild ansatz in Cartesian coordinates:
//   g_mu_nu = eta_mu_nu + H l_mu l_nu,   g^mu_nu = eta^mu_nu - H l^mu l^nu,
// with r the implicit solution of r^4 - (R^2 - a^2) r^2 - a^2 z^2 = 0,
// H = (2 M r - Q^2) r^2 / (r^4 + a^2 z^2), and l a null vector field. The
// Christoffel symbols are polynomial in (x, y, z, r), so there are no pole
// singularities.
// Reference: Visser, "The Kerr spacetime" (arXiv:0706.0622).

#include "sirius/base/contracts.h"
#include "sirius/core/coordinates.h"
#include "sirius/core/metrics/metric.h"
#include "sirius/core/metrics/registry.h"

#include <algorithm>
#include <cmath>
#include <limits>

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
    bool IsValidEvent(const Tensor<double, 4>& pos) const override;

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

    // Outer black-hole horizon. This is r+ for the asymptotically flat
    // Kerr-Newman sector and the smaller positive Kottler root for Lambda > 0;
    // returns -1 when no black-hole horizon exists.
    double OuterHorizonRadius() const;

    // Inner (Cauchy) horizon r- in the asymptotically flat Kerr-Newman sector.
    // The spherical uncharged sector retains the Schwarzschild limit r-=0;
    // returns -1 when no black-hole horizon exists.
    double InnerHorizonRadius() const;

    // True when the represented parameters contain a black-hole horizon.
    bool HasHorizon() const;

    // Larger positive Kottler root for Lambda > 0, including sqrt(3/Lambda)
    // for pure de Sitter; -1 outside the spherical positive-Lambda sector or
    // when the Kottler horizons do not exist.
    double CosmologicalHorizonRadius() const;

    // Static-coordinate lapse function f(r)=1-2M/r-Lambda*r^2/3 for the
    // spherical uncharged sector. Horizon-penetrating renderer observers are
    // ADM/Eulerian and are not silently treated as static observers.
    double KottlerStaticLapse(double radius) const;

    // Kerr-Newman ergosphere boundary radius at polar angle theta. The
    // cosmological sector has two static-limit horizons and must use the
    // explicit Kottler horizon authorities instead.
    double ErgosphereRadius(double theta) const;

    // ISCO radius for prograde equatorial orbits in the uncharged,
    // asymptotically flat Schwarzschild/Kerr sector.
    double IscoRadius() const;

    // Asymptotically flat Kerr-Newman extremality parameter
    // chi=sqrt(a^2+Q^2)/M, in [0,1] for black holes.
    double ExtremalityParameter() const;

  private:
    Config config_;
    KerrSchildParams params_;
};

inline KerrSchildFamily::KerrSchildFamily() {
    config_["mass"] = {1.0, 0.0, 100.0};
    config_["spin"] = {0.0, -0.998, 0.998};
    config_["charge"] = {0.0, -1.0, 1.0};
    config_["lambda"] = {0.0, 0.0, 0.1};
    params_ = KerrSchildParams::Schwarzschild(1.0);
}

inline KerrSchildFamily::KerrSchildFamily(const KerrSchildParams& params) : KerrSchildFamily() {
    SetParams(params);
}

inline bool IsRepresentedKerrSchildParameters(const KerrSchildParams& params) {
    if (!std::isfinite(params.M) || !std::isfinite(params.a) || !std::isfinite(params.Q) ||
        !std::isfinite(params.Lambda) || params.M < 0.0 || params.Lambda < 0.0) {
        return false;
    }
    if (params.M == 0.0 && (params.a != 0.0 || params.Q != 0.0)) return false;
    // This Cartesian ansatz represents the cosmological sector exactly only
    // for de Sitter and Schwarzschild-de Sitter.
    return params.Lambda == 0.0 || (params.a == 0.0 && params.Q == 0.0);
}

inline void KerrSchildFamily::SetParams(const KerrSchildParams& params) {
    SIRIUS_PRE(IsRepresentedKerrSchildParameters(params));
    if (!IsRepresentedKerrSchildParameters(params)) return;

    params_ = params;
    config_["mass"].value = params.M;
    config_["spin"].value = params.M == 0.0 ? 0.0 : params.a / params.M;
    config_["charge"].value = params.M == 0.0 ? 0.0 : params.Q / params.M;
    config_["lambda"].value = params.Lambda;
}

inline KerrSchildParams KerrSchildFamily::GetParams() const { return params_; }

inline bool KerrSchildFamily::IsValidEvent(const Tensor<double, 4>& pos) const {
    for (int component = 0; component < 4; ++component) {
        if (!std::isfinite(pos(component))) return false;
    }

    // The exact flat limit is defined everywhere.  Pure de Sitter is also
    // regular at its Cartesian origin even though the auxiliary radial null
    // direction is not unique there; Evaluate supplies the analytic limit.
    if (params_.M == 0.0 && params_.a == 0.0 && params_.Q == 0.0 && params_.Lambda == 0.0) {
        return true;
    }
    if (params_.M == 0.0 && params_.a == 0.0 && params_.Q == 0.0 && pos(1) == 0.0 &&
        pos(2) == 0.0 && pos(3) == 0.0) {
        return true;
    }

    const coordinates::Vec4Cart event{pos(0), pos(1), pos(2), pos(3)};
    return coordinates::TryKerrSchildRadiusDifferential(event, params_.a).has_value();
}

inline void KerrSchildFamily::SetParameter(const std::string& key, double value) {
    const auto found = config_.find(key);
    SIRIUS_PRE(found != config_.end());
    if (found == config_.end()) return;

    const bool in_range =
        std::isfinite(value) && value >= found->second.min && value <= found->second.max;
    SIRIUS_PRE(in_range);
    if (!in_range) return;

    const double previous = found->second.value;
    found->second.value = value;
    KerrSchildParams next;
    next.M = config_["mass"].value;
    next.a = config_["spin"].value * next.M;
    next.Q = config_["charge"].value * next.M;
    next.Lambda = config_["lambda"].value;
    const bool represented = IsRepresentedKerrSchildParameters(next);
    SIRIUS_PRE(represented);
    if (!represented) {
        found->second.value = previous;
        return;
    }
    params_ = next;
}

inline const char* KerrSchildFamily::GetName() const {
    if (params_.M == 0 && params_.Lambda == 0) return "Minkowski";
    if (params_.M == 0 && params_.Lambda != 0) return "de Sitter";
    if (params_.Lambda != 0) return "Schwarzschild-de Sitter";
    if (params_.a == 0 && params_.Q == 0 && params_.Lambda == 0) return "Schwarzschild";
    if (params_.Q == 0 && params_.Lambda == 0) return "Kerr";
    if (params_.a == 0 && params_.Lambda == 0) return "Reissner-Nordström";
    if (params_.Lambda == 0) return "Kerr-Newman";
    SIRIUS_ASSERT(false);
    return "Invalid Kerr-Schild metric";
}

inline double KerrSchildFamily::ComputeKerrRadius(double x, double y, double z) const {
    return coordinates::KerrSchildRadius(coordinates::Vec4Cart{0.0, x, y, z}, params_.a);
}

inline void KerrSchildFamily::ComputeNullVector(double x, double y, double z, double r,
                                                double l[4]) const {
    const bool represented = std::isfinite(x) && std::isfinite(y) && std::isfinite(z) &&
                             std::isfinite(r) && r > 0.0 && l != nullptr;
    SIRIUS_PRE(represented);
    if (!represented) return;
    double a = params_.a;
    double a2 = a * a;
    double r2 = r * r;
    double denom = r2 + a2;

    // l^mu = (1, (rx + ay)/(r^2 + a^2), (ry - ax)/(r^2 + a^2), z/r).
    l[0] = 1.0;
    l[1] = (r * x + a * y) / denom;
    l[2] = (r * y - a * x) / denom;
    l[3] = z / r;
}

inline double KerrSchildFamily::ComputeH(double r, double z) const {
    double M = params_.M;
    double a = params_.a;
    double Q = params_.Q;
    double Lambda = params_.Lambda;

    const bool finite = std::isfinite(r) && r >= 0.0 && std::isfinite(z);
    SIRIUS_PRE(finite);
    if (!finite) return std::numeric_limits<double>::quiet_NaN();

    const double r2 = r * r;
    const double a2 = a * a;

    // Divide numerator and denominator by r^2:
    //   H = (2 M r - Q^2) / (r^2 + a^2 (z/r)^2).
    // This is algebraically identical to the Kerr-Schild form but does not
    // square an O(r^2) denominator when the whole represented scene is small.
    double H = 0.0;
    if (M != 0.0 || Q != 0.0) {
        const bool positive_radius = r > 0.0;
        SIRIUS_PRE(positive_radius);
        if (!positive_radius) return std::numeric_limits<double>::quiet_NaN();
        const double cosine = z / r;
        const double sigma = r2 + a2 * cosine * cosine;
        const bool represented = std::isfinite(sigma) && sigma > 0.0;
        SIRIUS_PRE(represented);
        if (!represented) return std::numeric_limits<double>::quiet_NaN();
        H = (2.0 * M * r - Q * Q) / sigma;
    }

    // Schwarzschild-de Sitter is exactly Kerr-Schild with H += Lambda r^2/3 at
    // a = 0. The rotating de Sitter forms need a different ansatz, so Lambda with
    // a != 0 is not represented here; the configuration boundary rejects that
    // combination rather than letting an approximation stand in.
    if (a == 0.0 && Lambda != 0.0) {
        H += Lambda * r2 / 3.0;
    }

    return H;
}

inline void KerrSchildFamily::Evaluate(const Tensor<double, 4>& pos, Metric4d& g,
                                       Tensor<Dual<double>, 4, 4, 4>& dg) {
    const bool represented = IsValidEvent(pos);
    SIRIUS_PRE(represented);
    if (!represented) {
        const Dual<double> nan(std::numeric_limits<double>::quiet_NaN());
        for (int mu = 0; mu < 4; ++mu)
            for (int nu = 0; nu < 4; ++nu) g(mu, nu) = nan;
        for (int axis = 0; axis < 4; ++axis)
            for (int mu = 0; mu < 4; ++mu)
                for (int nu = 0; nu < 4; ++nu) dg(axis, mu, nu) = nan;
        return;
    }

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

    const auto initialise_minkowski = [&] {
        g.Zero();
        g(0, 0) = Dual<double>(-1.0);
        g(1, 1) = Dual<double>(1.0);
        g(2, 2) = Dual<double>(1.0);
        g(3, 3) = Dual<double>(1.0);
        dg.Zero();
    };

    // Exact analytic limits that do not need the non-unique radial null vector.
    if (M == 0.0 && a == 0.0 && Q == 0.0 && Lambda == 0.0) {
        initialise_minkowski();
        return;
    }
    if (M == 0.0 && a == 0.0 && Q == 0.0 && x == 0.0 && y == 0.0 && z == 0.0) {
        initialise_minkowski();
        return;
    }

    // One scale-safe radial authority supplies both r and its exact Cartesian
    // gradient.  In particular, no interval of small non-zero a is identified
    // with Schwarzschild and no r=0 event is moved off the singular sheet.
    const coordinates::Vec4Cart event{pos(0), x, y, z};
    const auto radial = coordinates::TryKerrSchildRadiusDifferential(event, a);
    SIRIUS_ASSERT(radial.has_value());
    if (!radial) {
        const Dual<double> nan(std::numeric_limits<double>::quiet_NaN());
        for (int mu = 0; mu < 4; ++mu)
            for (int nu = 0; nu < 4; ++nu) g(mu, nu) = nan;
        for (int axis = 0; axis < 4; ++axis)
            for (int mu = 0; mu < 4; ++mu)
                for (int nu = 0; nu < 4; ++nu) dg(axis, mu, nu) = nan;
        return;
    }
    double r = radial->radius;
    double r2 = r * r;

    // Derivatives of r with respect to the spatial coordinates.
    double dr_dx = radial->dx;
    double dr_dy = radial->dy;
    double dr_dz = radial->dz;
    double d_r2_dx = 2.0 * r * dr_dx;
    double d_r2_dy = 2.0 * r * dr_dy;
    double d_r2_dz = 2.0 * r * dr_dz;
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
    double H = ComputeH(r, z);

    // Derivatives of H.
    double dH[4] = {0.0, 0.0, 0.0, 0.0};
    const double cosine = z / r;
    const double sigma = r2 + a2 * cosine * cosine;
    const double asymptotic_h = (M != 0.0 || Q != 0.0) ? (2.0 * M * r - Q2) / sigma : 0.0;
    for (int lam = 1; lam <= 3; lam++) {
        if (M != 0.0 || Q != 0.0) {
            const double d_cosine = ((lam == 3) ? 1.0 / r : 0.0) - z * dr[lam] / r2;
            const double d_sigma = 2.0 * r * dr[lam] + 2.0 * a2 * cosine * d_cosine;
            const double d_numerator = 2.0 * M * dr[lam];
            dH[lam] = (d_numerator - asymptotic_h * d_sigma) / sigma;
        }

        // a = 0 cosmological term: d(Lambda r^2/3)/d x^lam = (2 Lambda/3) r dr.
        if (a == 0.0 && Lambda != 0.0) {
            dH[lam] += (2.0 * Lambda / 3.0) * r * dr[lam];
        }
    }

    // Metric g_mu_nu = eta_mu_nu + H l_mu l_nu.
    initialise_minkowski();

    for (int mu = 0; mu < 4; mu++) {
        for (int nu = 0; nu < 4; nu++) {
            g(mu, nu) = Dual<double>(g(mu, nu).real + H * l[mu] * l[nu]);
        }
    }

    // Derivatives d g_mu_nu / d x^lam = dH l_mu l_nu + H dl_mu l_nu + H l_mu dl_nu.
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
    const bool represented = IsValidEvent(pos);
    SIRIUS_PRE(represented);
    if (!represented) {
        const Dual<double> nan(std::numeric_limits<double>::quiet_NaN());
        for (int mu = 0; mu < 4; ++mu)
            for (int nu = 0; nu < 4; ++nu) g_inv(mu, nu) = nan;
        return false;
    }

    double x = pos(1), y = pos(2), z = pos(3);

    g_inv.Zero();
    g_inv(0, 0) = Dual<double>(-1.0);
    g_inv(1, 1) = Dual<double>(1.0);
    g_inv(2, 2) = Dual<double>(1.0);
    g_inv(3, 3) = Dual<double>(1.0);

    const bool flat =
        params_.M == 0.0 && params_.a == 0.0 && params_.Q == 0.0 && params_.Lambda == 0.0;
    const bool de_sitter_origin = params_.M == 0.0 && params_.a == 0.0 && params_.Q == 0.0 &&
                                  x == 0.0 && y == 0.0 && z == 0.0;
    if (flat || de_sitter_origin) return true;

    double r = ComputeKerrRadius(x, y, z);
    double l[4];
    ComputeNullVector(x, y, z, r, l);
    double H = ComputeH(r, z);

    double l_up[4] = {-1.0, l[1], l[2], l[3]};

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
    if (!std::isfinite(margin) || margin < 0.0) return false;
    for (int component = 0; component < 4; ++component) {
        if (!std::isfinite(pos(component))) return false;
    }

    // Compare in the Kerr radial coordinate: the horizon r = r+ is an oblate
    // surface in Cartesian coordinates, so a Cartesian-norm comparison would
    // misplace it for a != 0.
    double r = ComputeKerrRadius(pos(1), pos(2), pos(3));
    return r <= OuterHorizonRadius() * (1.0 + margin);
}

inline double KerrSchildFamily::KottlerStaticLapse(double radius) const {
    SIRIUS_PRE(params_.a == 0.0 && params_.Q == 0.0 && radius > 0.0 && std::isfinite(radius));
    if (params_.a != 0.0 || params_.Q != 0.0 || !(radius > 0.0) || !std::isfinite(radius)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return sirius::core::KottlerStaticLapse(params_.M, params_.Lambda, radius);
}

inline double KerrSchildFamily::OuterHorizonRadius() const {
    if (!HasHorizon()) return -1.0;

    double M = params_.M;
    double a = params_.a;
    double Q = params_.Q;

    if (params_.Lambda > 0.0) {
        return KottlerBlackHoleHorizonRadius(M, params_.Lambda);
    }

    double discriminant = M * M - a * a - Q * Q;
    return M + std::sqrt(discriminant);
}

inline double KerrSchildFamily::InnerHorizonRadius() const {
    if (!HasHorizon()) return -1.0;

    double M = params_.M;
    double a = params_.a;
    double Q = params_.Q;

    if (params_.Lambda > 0.0) return 0.0;

    double discriminant = M * M - a * a - Q * Q;
    return M - std::sqrt(discriminant);
}

inline bool KerrSchildFamily::HasHorizon() const {
    double M = params_.M;
    double a = params_.a;
    double Q = params_.Q;

    if (!(M > 0.0)) return false;
    if (params_.Lambda > 0.0) {
        SIRIUS_ASSERT(a == 0.0 && Q == 0.0);
        return 9.0 * params_.Lambda * M * M <= 1.0;
    }
    return (M * M >= a * a + Q * Q);
}

inline double KerrSchildFamily::CosmologicalHorizonRadius() const {
    if (!(params_.Lambda > 0.0) || params_.a != 0.0 || params_.Q != 0.0) return -1.0;
    return KottlerCosmologicalHorizonRadius(params_.M, params_.Lambda);
}

inline double KerrSchildFamily::ErgosphereRadius(double theta) const {
    SIRIUS_PRE(params_.Lambda == 0.0);
    if (params_.Lambda != 0.0) return std::numeric_limits<double>::quiet_NaN();

    double M = params_.M;
    double a = params_.a;
    double Q = params_.Q;

    // Static limit g_tt = 0. For Kerr-Newman: r^2 - 2Mr + Q^2 + a^2 cos^2 theta = 0,
    // giving r = M + sqrt(M^2 - Q^2 - a^2 cos^2 theta).
    double cos2th = std::cos(theta) * std::cos(theta);

    double disc_charged = M * M - Q * Q - a * a * cos2th;
    if (disc_charged < 0) {
        // The static-limit surface does not exist at this angle.
        return std::numeric_limits<double>::quiet_NaN();
    }

    return M + std::sqrt(disc_charged);
}

inline double KerrSchildFamily::IscoRadius() const {
    double M = params_.M;
    double a = params_.a;
    double Q = params_.Q;

    const bool represented = params_.Lambda == 0.0 && Q == 0.0 && M > 0.0 && std::isfinite(M) &&
                             std::isfinite(a) && std::abs(a) <= M;
    SIRIUS_PRE(represented);
    if (!represented) return std::numeric_limits<double>::quiet_NaN();

    // Schwarzschild (a = 0, Q = 0): r_ISCO = 6M.
    if (a == 0.0) {
        return 6.0 * M;
    }

    // Kerr (Q = 0): Bardeen-Press-Teukolsky closed form, exact.
    double a_star = std::abs(a / M);

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
    SIRIUS_PRE(params_.Lambda == 0.0);
    if (params_.Lambda != 0.0) return std::numeric_limits<double>::quiet_NaN();

    double M = params_.M;
    double a = params_.a;
    double Q = params_.Q;

    if (M <= 0) return 0;

    // chi = sqrt(a^2 + Q^2) / M: 0 for Schwarzschild, 1 for extremal, > 1 naked.
    return std::sqrt(a * a + Q * Q) / M;
}

}  // namespace sirius::core
