#pragma once

// Novikov-Thorne geometrically-thin, optically-thick accretion disk around a
// Kerr black hole, with the full Page & Thorne (1974) relativistic flux. Ported
// from PHAD001A.h.
//
// Temperature profile T(r) = [3 G M Mdot/(8 pi sigma r^3) f(r)]^(1/4), with f(r)
// the Page-Thorne correction. References: Novikov & Thorne (1973); Page & Thorne
// (1974); James et al. (2015), "DNGR" Section 3.2.
//
// Oracle-decoupling note: the legacy header pulled in the double-precision Kerr
// Boyer-Lindquist oracle (PHMT000B/PHMT100B) solely to name the BL 4-vector type
// used by IntersectRay. Core must not depend on the oracle layer, so that type is
// replaced by the core coordinates::Vec4Bl (same {t, r, theta, phi} layout and
// member access); behaviour is unchanged.

#include "sirius/base/contracts.h"
#include "sirius/core/constants.h"
#include "sirius/core/coordinates.h"
#include "sirius/core/disk/disk_model.h"
#include "sirius/core/spectral/spectral_radiance.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>

namespace sirius::core {

// Novikov-Thorne thin disk model, implementing IDiskModel.
class AccretionDiskD : public IDiskModel {
  public:
    struct Config {
        double M = 1.0;                             // Black hole mass [M_sun].
        double a_star = 0.0;                        // Dimensionless spin a/M in [-1, 1].
        double Mdot = 1e-8;                         // Accretion rate [M_sun/year].
        double r_inner = 0;                         // Inner edge (0 = ISCO, auto-computed).
        double r_outer = 500;                       // Outer edge [GM/c^2].
        double inclination = std::numbers::pi / 4;  // Disk inclination to the observer [rad].

        void Validate() {
            if (!std::isfinite(M) || M <= 0) M = 1.0;
            if (!std::isfinite(a_star)) a_star = 0.0;
            a_star = std::clamp(a_star, -1.0, 1.0);
            if (!std::isfinite(Mdot) || Mdot <= 0) Mdot = 1e-8;
            if (!std::isfinite(r_inner) || r_inner < 0) r_inner = 0.0;
            if (!std::isfinite(r_outer)) r_outer = 500.0;
            if (!(r_outer > r_inner)) {
                r_outer = std::nextafter(r_inner, std::numeric_limits<double>::infinity());
                if (!std::isfinite(r_outer)) {
                    r_inner = 0.0;
                    r_outer = 500.0;
                }
            }
            if (!std::isfinite(inclination)) inclination = std::numbers::pi / 4;
            inclination = std::clamp(inclination, 0.0, std::numbers::pi);
        }
    };

    AccretionDiskD() : config_() { Init(); }

    explicit AccretionDiskD(const Config& config) : config_(config) { Init(); }

  private:
    void Init() {
        config_.Validate();
        // Derived quantities.
        m_kg_ = config_.M * constants::physical::kSolarMass;
        gm_ = constants::physical::kGravitation * m_kg_;
        rs_ = 2 * gm_ / (constants::physical::kSpeedOfLight * constants::physical::kSpeedOfLight);
        a_ = config_.a_star;  // a/M in geometric units.

        r_isco_ = ComputeIsco(a_);

        r_inner_ = (config_.r_inner > 0) ? config_.r_inner : r_isco_;
        r_outer_ = config_.r_outer;

        // Accretion rate M_sun/yr -> kg/s.
        mdot_si_ = config_.Mdot * constants::physical::kSolarMass / (365.25 * 24 * 3600);
    }

  public:
    // ISCO radius in units of M:
    // r_ISCO = M {3 + Z2 -/+ sqrt[(3 - Z1)(3 + Z1 + 2 Z2)]}, -/+ prograde/retrograde.
    static double ComputeIsco(double a_star) {
        double a = std::abs(a_star);
        if (a > 1) a = 1;  // Clamp to extremal.

        double Z1 = 1 + std::cbrt(1 - a * a) * (std::cbrt(1 + a) + std::cbrt(1 - a));
        double Z2 = std::sqrt(3 * a * a + Z1 * Z1);

        double r_isco;
        if (a_star >= 0) {
            r_isco = 3 + Z2 - std::sqrt((3 - Z1) * (3 + Z1 + 2 * Z2));
        } else {
            r_isco = 3 + Z2 + std::sqrt((3 - Z1) * (3 + Z1 + 2 * Z2));
        }

        return r_isco;  // In units of M.
    }

    // Specific energy for a circular orbit (Page & Thorne 1974, Eq. 15.11):
    // E(r) = (r^2 - 2Mr + a sqrt(Mr)) / (r sqrt(r^2 - 3Mr + 2a sqrt(Mr))).
    static double SpecificEnergy(double r, double M, double a) {
        double sqrtMr = std::sqrt(M * r);
        double numerator = r * r - 2 * M * r + a * sqrtMr;
        double denominator = r * std::sqrt(r * r - 3 * M * r + 2 * a * sqrtMr);
        if (std::abs(denominator) < constants::tolerances::kDivisionSafeEps) return 1.0;
        return numerator / denominator;
    }

    // Specific angular momentum for a circular orbit (Page & Thorne 1974, Eq. 15.12):
    // L(r) = sqrt(M) (r^2 - 2a sqrt(Mr) + a^2) /
    //        (sqrt(r) sqrt(r^2 - 3Mr + 2a sqrt(Mr))).
    static double SpecificAngularMomentum(double r, double M, double a) {
        double sqrtM = std::sqrt(M);
        double sqrtMr = std::sqrt(M * r);
        double numerator = sqrtM * (r * r - 2 * a * sqrtMr + a * a);
        double denominator = std::sqrt(r) * std::sqrt(r * r - 3 * M * r + 2 * a * sqrtMr);
        if (std::abs(denominator) < constants::tolerances::kDivisionSafeEps) return 0.0;
        return numerator / denominator;
    }

    // Angular velocity for a circular orbit (Bardeen 1970):
    // Omega(r) = sqrt(M) / (r^(3/2) + a sqrt(M)).
    static double AngularVelocity(double r, double M, double a) {
        double sqrtM = std::sqrt(M);
        return sqrtM / (std::pow(r, 1.5) + a * sqrtM);
    }

    // dL/dr, the analytic derivative of the specific angular momentum required by
    // the Page-Thorne flux integral (derived from Page & Thorne 1974, Eq. 15.12
    // by the quotient rule; the derivative name is kept in domain notation).
    static double dLdr(double r, double M, double a) {
        if (r <= 0 || M <= 0) return 0;

        const double sqrtM = std::sqrt(M);
        const double sqrtMr = std::sqrt(M * r);
        const double sqrtR = std::sqrt(r);
        const double N = r * r - 2 * a * sqrtMr + a * a;
        const double radicand = r * r - 3 * M * r + 2 * a * sqrtMr;
        if (radicand <= constants::tolerances::kDivisionSafeEps) return 0;

        const double denominator = sqrtR * std::sqrt(radicand);
        const double dN = 2 * r - a * sqrtM / sqrtR;
        const double dRadicand = 2 * r - 3 * M + a * sqrtM / sqrtR;
        const double dLogDenominator = 0.5 / r + 0.5 * dRadicand / radicand;
        return sqrtM * (dN - N * dLogDenominator) / denominator;
    }

    // dOmega/dr, the analytic derivative of the angular velocity (from Bardeen
    // 1970): dOmega/dr = -3 sqrt(Mr) / [2 (r^(3/2) + a sqrt(M))^2].
    static double dOmegadr(double r, double M, double a) {
        if (r <= 0 || M <= 0) return 0;

        double sqrtM = std::sqrt(M);
        double sqrtMr = std::sqrt(M * r);
        double r32 = std::pow(r, 1.5);   // r^(3/2)
        double denom = r32 + a * sqrtM;  // r^(3/2) + a sqrt(M)

        if (std::abs(denom) < constants::tolerances::kFluxIntegralEps) return 0;

        return -3.0 * sqrtMr / (2.0 * denom * denom);
    }

    // Full Page & Thorne (1974) relativistic flux, using 16-point Gauss-Legendre
    // quadrature of the flux integral over [r_ISCO, r].
    double FullPageThorneFlux(double r) const {
        if (r <= r_isco_ || r > r_outer_) return 0;

        double M = 1.0;  // Working in units of M = 1.
        double a = a_;

        double E_r = SpecificEnergy(r, M, a);
        double L_r = SpecificAngularMomentum(r, M, a);
        double Omega_r = AngularVelocity(r, M, a);

        // 16-point Gauss-Legendre nodes and weights on [-1, 1]; the interior
        // nodes avoid the ISCO boundary where the integrand is singular
        // (Abramowitz & Stegun, Table 25.4).
        static constexpr double kGlNodes[16] = {
            -0.9894009349916499, -0.9445750230732326, -0.8656312023878318, -0.7554044083550030,
            -0.6178762444026438, -0.4580167776572274, -0.2816035507792589, -0.0950125098376374,
            0.0950125098376374,  0.2816035507792589,  0.4580167776572274,  0.6178762444026438,
            0.7554044083550030,  0.8656312023878318,  0.9445750230732326,  0.9894009349916499};
        static constexpr double kGlWeights[16] = {
            0.0271524594117541, 0.0622535239386479, 0.0951585116824928, 0.1246289712555339,
            0.1495959888165767, 0.1691565193950025, 0.1826034150449236, 0.1894506104550685,
            0.1894506104550685, 0.1826034150449236, 0.1691565193950025, 0.1495959888165767,
            0.1246289712555339, 0.0951585116824928, 0.0622535239386479, 0.0271524594117541};

        // Map [-1, 1] -> [r_ISCO, r].
        double half_width = 0.5 * (r - r_isco_);
        double midpoint = 0.5 * (r + r_isco_);

        double integral = 0;
        for (int i = 0; i < 16; ++i) {
            double r_i = midpoint + half_width * kGlNodes[i];

            double E_i = SpecificEnergy(r_i, M, a);
            double L_i = SpecificAngularMomentum(r_i, M, a);
            double Omega_i = AngularVelocity(r_i, M, a);
            double dL_i = dLdr(r_i, M, a);

            double integrand = (E_i - Omega_i * L_i) * dL_i;
            integral += kGlWeights[i] * integrand;
        }
        integral *= half_width;

        double E_minus_OmegaL = E_r - Omega_r * L_r;
        if (std::abs(E_minus_OmegaL) < constants::tolerances::kDivisionSafeEps) return 0;

        double dOmega_dr = dOmegadr(r, M, a);

        // Page-Thorne Q-factor (relativistic correction):
        // Q = (-Omega,r / (E - Omega L)^2) integral (E - Omega L) L,r dr.
        double Q = (-dOmega_dr / (E_minus_OmegaL * E_minus_OmegaL)) * integral;

        // Classical Newtonian flux prefactor in SI.
        double r_physical = r * rs_ / 2;  // r in metres.
        double F_newt =
            (3 * gm_ * mdot_si_) / (8 * std::numbers::pi * r_physical * r_physical * r_physical);

        double F = F_newt * Q;

        if (!std::isfinite(F) || F < 0) return 0;

        return F;
    }

    // Radiative flux F(r) = (3 G M Mdot)/(8 pi r^3) Q(r) (IDiskModel::Flux).
    double Flux(double r) const override {
        if (r < r_inner_ || r > r_outer_) return 0;

        double r_M = r;  // r is in units of M.

        // Must be outside the ISCO for emission.
        if (r_M <= r_isco_) return 0;

        // Full Page-Thorne everywhere outside the inner-edge buffer. A numerical
        // failure must not substitute a different relativistic correction.
        if (r_M > r_isco_ * constants::disk::kInnerEdgeBuffer) {
            const double F_PT = FullPageThorneFlux(r);
            SIRIUS_ASSERT(std::isfinite(F_PT) && F_PT > 0.0);
            return std::isfinite(F_PT) && F_PT > 0.0 ? F_PT : 0.0;
        }

        // The zero-torque boundary is exactly dark inside the guarded buffer.
        return 0.0;
    }

    // Temperature T(r) = [F(r) / sigma_SB]^(1/4) (IDiskModel::Temperature).
    double Temperature(double r, [[maybe_unused]] double z = 0) const override {
        double F = Flux(r);
        if (F <= 0) return 0;

        return std::pow(F / constants::physical::kStefanBoltzmann, 0.25);
    }

    // Peak temperature, near 1.5 r_ISCO.
    double PeakTemperature() const {
        double r_peak = 1.5 * r_isco_;
        return Temperature(r_peak);
    }

    // Blackbody emission spectrum at the local temperature.
    SpectralRadiance EmissionSpectrum(double r) const {
        double T = Temperature(r);
        if (T <= 0) return SpectralRadiance::Zero();

        return SpectralRadiance::Blackbody(T);
    }

    // Limb darkening I(theta) = I0 (1 + u cos theta) / (1 + u).
    static double LimbDarkening(double cos_theta, double u = 0.6) {
        if (cos_theta <= 0) return 0;
        return (1 + u * cos_theta) / (1 + u);
    }

    // Whether a point is within the disk (IDiskModel::IsInDisk).
    bool IsInDisk(double r, double theta) const override {
        // Equatorial plane: |theta - pi/2| < epsilon.
        double equatorial_distance = std::abs(theta - std::numbers::pi / 2);
        if (equatorial_distance > 0.01) return false;  // ~0.5 degrees from equator.

        return (r >= r_inner_ && r <= r_outer_);
    }

    // Ray/disk intersection: returns the affine parameter at intersection, or -1
    // for no intersection. The BL 4-vector type is the core coordinates::Vec4Bl
    // (see the oracle-decoupling note above).
    double IntersectRay(const coordinates::Vec4Bl& x, const coordinates::Vec4Bl& k,
                        double lambda_max = 1000) const {
        // Check when theta crosses pi/2.
        double theta0 = x.theta;
        double dtheta = k.theta;  // Approximate: k_theta ~ dtheta/dlambda.

        if (std::abs(dtheta) < 1e-10) return -1;  // Parallel to the equator.

        double lambda_cross = (std::numbers::pi / 2 - theta0) / dtheta;

        if (lambda_cross < 0 || lambda_cross > lambda_max) return -1;

        // Estimate r at crossing.
        double r_cross = x.r + k.r * lambda_cross;

        if (r_cross >= r_inner_ && r_cross <= r_outer_) {
            return lambda_cross;
        }

        return -1;
    }

    // Accessors.
    const Config& GetConfig() const { return config_; }
    double IscoRadius() const { return r_isco_; }
    double SpinParameter() const { return a_; }
    double MassKg() const { return m_kg_; }
    double SchwarzschildRadius() const { return rs_; }

    // Alias for Temperature() for API compatibility.
    double EffectiveTemperature(double r) const { return Temperature(r); }

    // --- IDiskModel interface ---

    const char* ModelName() const override { return "Novikov-Thorne"; }

    double InnerRadius() const override { return r_inner_; }

    double OuterRadius() const override { return r_outer_; }

    // Thin disk: H = 0.
    double HalfThickness([[maybe_unused]] double r) const override { return 0.0; }

    // Thin disk: normalised value for optical-depth calculations.
    double Density([[maybe_unused]] double r, [[maybe_unused]] double z = 0) const override {
        return 1.0;
    }

    double AngularVelocity(double r) const override {
        return AngularVelocity(r, 1.0, a_);  // Static version with M = 1.
    }

    void SetBlackHoleParameters(double mass, double spin) override {
        config_.M = mass;
        config_.a_star = spin;
        Init();  // Recompute derived quantities.
    }

    double BlackHoleMass() const override { return config_.M; }

    double BlackHoleSpin() const override { return config_.a_star; }

  private:
    Config config_;

    // Derived quantities.
    double m_kg_;     // Mass in kg.
    double gm_;       // GM in SI.
    double rs_;       // Schwarzschild radius in metres.
    double a_;        // Dimensionless spin.
    double r_isco_;   // ISCO radius in units of M.
    double r_inner_;  // Inner edge in units of M.
    double r_outer_;  // Outer edge in units of M.
    double mdot_si_;  // Accretion rate in kg/s.
};

}  // namespace sirius::core
