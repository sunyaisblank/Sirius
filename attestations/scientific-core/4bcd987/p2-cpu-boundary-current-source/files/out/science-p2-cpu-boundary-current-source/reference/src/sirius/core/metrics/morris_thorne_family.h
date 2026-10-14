#pragma once

// Unified Morris-Thorne traversable wormhole family in Schwarzschild-like
// coordinates. Ported from PHMT101A.h.
//
//   ds^2 = -e^{2 Phi(r)} dt^2 + dr^2/(1 - b(r)/r) + r^2 (dtheta^2 + sin^2 theta dphi^2)
// with Phi(r) the redshift function and b(r) the shape function; traversability
// needs b(b0) = b0 and Phi finite everywhere.
// Reference: Morris & Thorne, Am. J. Phys. 56, 395 (1988).

#include "sirius/base/contracts.h"
#include "sirius/core/metrics/metric.h"
#include "sirius/core/metrics/registry.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <numbers>

namespace sirius::core {

enum class WormholeShapeType {
    Ellis,           // b(r) = b0^2 / r.
    ZeroTidal,       // b(r) = b0 (constant).
    AbsurdlyBenign,  // b(r) = b0 (2 - b0/r).
    Custom           // User-defined.
};

// Callback for a custom shape function b(r).
using ShapeFunctionCallback = std::function<double(double)>;

// Callback for a custom shape derivative db/dr.
using ShapeDerivativeCallback = std::function<double(double)>;

// Morris-Thorne family parameters.
struct MorrisThorneParams {
    double b0 = 1.0;    // Throat radius.
    double Phi0 = 0.0;  // Redshift at the throat (0 = zero-tidal).
    WormholeShapeType shape_type = WormholeShapeType::Ellis;

    // Custom shape callbacks, used when shape_type == Custom.
    ShapeFunctionCallback custom_shape_func = nullptr;
    ShapeDerivativeCallback custom_shape_deriv_func = nullptr;

    static MorrisThorneParams Ellis(double b0) {
        return {b0, 0.0, WormholeShapeType::Ellis, nullptr, nullptr};
    }
    static MorrisThorneParams ZeroTidal(double b0) {
        return {b0, 0.0, WormholeShapeType::ZeroTidal, nullptr, nullptr};
    }

    // Custom wormhole with a user-defined shape function b(r), which must satisfy
    // b(b0) = b0; the derivative is computed numerically when null.
    static MorrisThorneParams Custom(double b0, ShapeFunctionCallback shape_func,
                                     ShapeDerivativeCallback deriv_func = nullptr) {
        SIRIUS_PRE(static_cast<bool>(shape_func));
        MorrisThorneParams p;
        p.b0 = b0;
        p.Phi0 = 0.0;
        p.shape_type = WormholeShapeType::Custom;
        p.custom_shape_func = std::move(shape_func);
        p.custom_shape_deriv_func = std::move(deriv_func);
        return p;
    }
};

// Morris-Thorne family metric.
class MorrisThorneFamily : public IMetric {
  public:
    MorrisThorneFamily();
    explicit MorrisThorneFamily(const MorrisThorneParams& params);

    void Evaluate(const Tensor<double, 4>& pos, Metric4d& g,
                  Tensor<Dual<double>, 4, 4, 4>& dg) override;

    const Config& GetParameters() const override { return config_; }
    void SetParameter(const std::string& key, double value) override;
    const char* GetName() const override;
    bool IsValidEvent(const Tensor<double, 4>& pos) const override;

    void SetParams(const MorrisThorneParams& params);
    MorrisThorneParams GetParams() const;

    // Shape function b(r) and its derivative.
    double ShapeFunction(double r) const;
    double ShapeFunctionDerivative(double r) const;

    // Redshift function Phi(r) and its derivative.
    double RedshiftFunction(double r) const;
    double RedshiftFunctionDerivative(double r) const;

    // Flare-out condition for traversability: (b - r db/dr)/(2 b^2) > 0 at r = b0,
    // equivalently b'(b0) < 1.
    bool ValidateFlareOutCondition() const;

    // Flare-out parameter (b - r db/dr)/(2 b^2) at a radius.
    double FlareOutParameter(double r) const;

  private:
    Config config_;
    MorrisThorneParams params_;
};

inline MorrisThorneFamily::MorrisThorneFamily() {
    config_["throat_radius"] = {kDefaultMorrisThorneThroatRadius, kMinMorrisThorneThroatRadius,
                                kMaxMorrisThorneThroatRadius};
    config_["redshift"] = {0.0, -10.0, 10.0};
    params_ = MorrisThorneParams::Ellis(1.0);
}

inline MorrisThorneFamily::MorrisThorneFamily(const MorrisThorneParams& params)
    : MorrisThorneFamily() {
    SetParams(params);
}

inline void MorrisThorneFamily::SetParams(const MorrisThorneParams& params) {
    const bool known_shape = params.shape_type == WormholeShapeType::Ellis ||
                             params.shape_type == WormholeShapeType::ZeroTidal ||
                             params.shape_type == WormholeShapeType::AbsurdlyBenign ||
                             params.shape_type == WormholeShapeType::Custom;
    const bool valid = known_shape && std::isfinite(params.b0) &&
                       params.b0 >= kMinMorrisThorneThroatRadius &&
                       params.b0 <= kMaxMorrisThorneThroatRadius && std::isfinite(params.Phi0) &&
                       params.Phi0 >= -10.0 && params.Phi0 <= 10.0 &&
                       (params.shape_type != WormholeShapeType::Custom ||
                        static_cast<bool>(params.custom_shape_func));
    SIRIUS_PRE(valid);
    if (!valid) return;

    params_ = params;
    config_["throat_radius"].value = params.b0;
    config_["redshift"].value = params.Phi0;
}

inline MorrisThorneParams MorrisThorneFamily::GetParams() const { return params_; }

inline bool MorrisThorneFamily::IsValidEvent(const Tensor<double, 4>& pos) const {
    for (int component = 0; component < 4; ++component) {
        if (!std::isfinite(pos(component))) return false;
    }
    const double r = pos(1);
    const double theta = pos(2);
    if (!(r > params_.b0 && theta > 0.0 && theta < std::numbers::pi)) return false;
    const double b = ShapeFunction(r);
    return std::isfinite(b) && 1.0 - b / r > 0.0;
}

inline void MorrisThorneFamily::SetParameter(const std::string& key, double value) {
    const auto found = config_.find(key);
    SIRIUS_PRE(found != config_.end());
    if (found == config_.end()) return;
    const bool in_range =
        std::isfinite(value) && value >= found->second.min && value <= found->second.max;
    SIRIUS_PRE(in_range);
    if (!in_range) return;
    found->second.value = value;
    params_.b0 = config_["throat_radius"].value;
    params_.Phi0 = config_["redshift"].value;
}

inline const char* MorrisThorneFamily::GetName() const {
    switch (params_.shape_type) {
        case WormholeShapeType::Ellis:
            return "Ellis Drainhole";
        case WormholeShapeType::ZeroTidal:
            return "Zero-Tidal Wormhole";
        case WormholeShapeType::AbsurdlyBenign:
            return "Absurdly Benign Wormhole";
        case WormholeShapeType::Custom:
            return "Custom Morris-Thorne Wormhole";
        default:
            SIRIUS_PRE(false);
            return "Invalid Morris-Thorne Wormhole";
    }
}

inline double MorrisThorneFamily::ShapeFunction(double r) const {
    const double b0 = params_.b0;
    const bool represented = std::isfinite(r) && r >= b0;
    SIRIUS_PRE(represented);
    if (!represented) return std::numeric_limits<double>::quiet_NaN();

    switch (params_.shape_type) {
        case WormholeShapeType::Ellis:
            return b0 * b0 / r;
        case WormholeShapeType::ZeroTidal:
            return b0;
        case WormholeShapeType::AbsurdlyBenign:
            return b0 * (2.0 - b0 / r);
        case WormholeShapeType::Custom:
            if (params_.custom_shape_func) {
                return params_.custom_shape_func(r);
            }
            SIRIUS_ASSERT(false);
            return std::numeric_limits<double>::quiet_NaN();
        default:
            SIRIUS_ASSERT(false);
            return std::numeric_limits<double>::quiet_NaN();
    }
}

inline double MorrisThorneFamily::ShapeFunctionDerivative(double r) const {
    const double b0 = params_.b0;
    const bool represented = std::isfinite(r) && r >= b0;
    SIRIUS_PRE(represented);
    if (!represented) return std::numeric_limits<double>::quiet_NaN();

    switch (params_.shape_type) {
        case WormholeShapeType::Ellis:
            return -b0 * b0 / (r * r);
        case WormholeShapeType::ZeroTidal:
            return 0.0;
        case WormholeShapeType::AbsurdlyBenign:
            return b0 * b0 / (r * r);
        case WormholeShapeType::Custom:
            if (params_.custom_shape_deriv_func) {
                return params_.custom_shape_deriv_func(r);
            }
            // Otherwise compute numerically via central differences.
            if (params_.custom_shape_func) {
                const double h = 1e-6 * r;
                if (r - h < b0) {
                    const double b_here = params_.custom_shape_func(r);
                    const double b_plus = params_.custom_shape_func(r + h);
                    const double b_twice_plus = params_.custom_shape_func(r + 2.0 * h);
                    return (-3.0 * b_here + 4.0 * b_plus - b_twice_plus) / (2.0 * h);
                }
                const double b_plus = params_.custom_shape_func(r + h);
                const double b_minus = params_.custom_shape_func(r - h);
                return (b_plus - b_minus) / (2.0 * h);
            }
            SIRIUS_ASSERT(false);
            return std::numeric_limits<double>::quiet_NaN();
        default:
            SIRIUS_ASSERT(false);
            return std::numeric_limits<double>::quiet_NaN();
    }
}

inline double MorrisThorneFamily::FlareOutParameter(double r) const {
    double b = ShapeFunction(r);
    double db_dr = ShapeFunctionDerivative(r);

    // (b - r db/dr) / (2 b^2), which must be > 0 for a traversable geometry.
    double b2 = b * b;
    if (b2 < 1e-20) return 0;

    return (b - r * db_dr) / (2.0 * b2);
}

inline bool MorrisThorneFamily::ValidateFlareOutCondition() const {
    // Flare-out at the throat requires b'(b0) < 1, equivalently
    // FlareOutParameter(b0) > 0.
    double b0 = params_.b0;
    double db_dr_at_throat = ShapeFunctionDerivative(b0);

    if (db_dr_at_throat >= 1.0) {
        return false;  // Violates flare-out.
    }

    double fop = FlareOutParameter(b0);
    return fop > 0;
}

inline double MorrisThorneFamily::RedshiftFunction([[maybe_unused]] double r) const {
    // RECORDED DECISION: this family implements the zero-tidal-force subclass of
    // Morris-Thorne wormholes, Phi(r) = Phi0 (constant). That is a physically
    // legitimate subclass (Morris & Thorne 1988, Section III.F), not an
    // approximation of a more general wormhole; variable Phi(r) is out of scope
    // until a use case needs it. The derivative below is exactly zero as a
    // consequence, not as a placeholder.
    return params_.Phi0;
}

inline double MorrisThorneFamily::RedshiftFunctionDerivative([[maybe_unused]] double r) const {
    // Constant redshift, so the derivative is r-independent zero.
    return 0.0;
}

// Exact isotropic Cartesian chart of the zero-tidal Ellis member rendered by
// Sirius.  The spherical family above remains the areal-radius authority:
//
//   ds^2 = -e^{2 Phi0} dt^2 + dr^2/(1 - b0^2/r^2) + r^2 dOmega^2.
//
// With proper radial distance l and isotropic radius rho,
//
//   r = rho + b0^2/(4 rho),       l = rho - b0^2/(4 rho),
//
// the exact spatial metric is conformally flat:
//
//   g_ij = A(rho)^2 delta_ij,     A = 1 + b0^2/(4 rho^2).
//
// Unlike areal radius, rho is a regular chart coordinate at the throat
// rho=b0/2.  The punctured isotropic chart is global across both asymptotic
// ends: 0<rho<b0/2 is the second sheet and rho->0 is its spatial infinity.
// Geometry exposes the throat identity but does not classify it as a horizon;
// OneSheetCapture versus TwoSheet is a tracer boundary condition.
[[nodiscard]] inline double EllisArealRadiusFromIsotropic(double b0, double rho) {
    const bool represented = std::isfinite(b0) && b0 > 0.0 && std::isfinite(rho) && rho > 0.0;
    SIRIUS_PRE(represented);
    if (!represented) return std::numeric_limits<double>::quiet_NaN();
    return rho + b0 * b0 / (4.0 * rho);
}

[[nodiscard]] inline double EllisProperRadialDistanceFromIsotropic(double b0, double rho) {
    const bool represented = std::isfinite(b0) && b0 > 0.0 && std::isfinite(rho) && rho > 0.0;
    SIRIUS_PRE(represented);
    if (!represented) return std::numeric_limits<double>::quiet_NaN();
    return rho - b0 * b0 / (4.0 * rho);
}

// The inversion rho'=(b0/2)^2/rho exchanges the two asymptotic ends while
// preserving areal radius and reversing signed proper radial distance.
[[nodiscard]] inline double EllisInvertedIsotropicRadius(double b0, double rho) {
    const bool represented = std::isfinite(b0) && b0 > 0.0 && std::isfinite(rho) && rho > 0.0;
    SIRIUS_PRE(represented);
    if (!represented) return std::numeric_limits<double>::quiet_NaN();
    return b0 * b0 / (4.0 * rho);
}

// Convert a local sky direction expressed in the isotropic x-chart at the
// second end into the inversion-related asymptotically Cartesian chart.  The
// inversion Jacobian is a radial Householder reflection; its conformal scale
// cancels from orthonormal direction components.
[[nodiscard]] inline std::optional<Vec4> MapEllisSecondSheetSkyDirection(
    const Vec4& position, const Vec4& chart_direction) {
    const double radius = std::hypot(position(1), position(2), position(3));
    const double direction_norm =
        std::hypot(chart_direction(1), chart_direction(2), chart_direction(3));
    if (!(std::isfinite(radius) && radius > 0.0 && std::isfinite(direction_norm) &&
          direction_norm > 0.0)) {
        return std::nullopt;
    }
    for (int component = 0; component < 4; ++component) {
        if (!std::isfinite(position(component)) || !std::isfinite(chart_direction(component))) {
            return std::nullopt;
        }
    }
    const double nx = position(1) / radius;
    const double ny = position(2) / radius;
    const double nz = position(3) / radius;
    const double radial =
        nx * chart_direction(1) + ny * chart_direction(2) + nz * chart_direction(3);
    Vec4 mapped = chart_direction;
    mapped(1) -= 2.0 * radial * nx;
    mapped(2) -= 2.0 * radial * ny;
    mapped(3) -= 2.0 * radial * nz;
    return mapped;
}

class MorrisThorneCartesian : public IMetric {
  public:
    MorrisThorneCartesian() : family_() {
        config_["throat_radius"] = family_.GetParameters().at("throat_radius");
    }
    explicit MorrisThorneCartesian(const MorrisThorneParams& params) : family_(params) {
        const bool exact_live_ellis = params.shape_type == WormholeShapeType::Ellis &&
                                      params.Phi0 == 0.0 && !params.custom_shape_func &&
                                      !params.custom_shape_deriv_func;
        SIRIUS_PRE(exact_live_ellis);
        if (!exact_live_ellis) return;
        config_["throat_radius"] = family_.GetParameters().at("throat_radius");
    }

    void Evaluate(const Tensor<double, 4>& pos, Metric4d& g,
                  Tensor<Dual<double>, 4, 4, 4>& dg) override;

    const Config& GetParameters() const override { return config_; }
    void SetParameter(const std::string& key, double value) override {
        const bool represented_parameter = key == "throat_radius";
        SIRIUS_PRE(represented_parameter);
        if (!represented_parameter) return;
        family_.SetParameter(key, value);
        config_["throat_radius"] = family_.GetParameters().at("throat_radius");
    }
    const char* GetName() const override { return family_.GetName(); }

    bool IsValidEvent(const Tensor<double, 4>& pos) const override;
    bool InverseMetric(const Tensor<double, 4>& pos, Metric4d& g_inv) const override;
    std::optional<double> IsotropicEllisThroatRadius() const override {
        return IsotropicThroatRadius();
    }
    [[nodiscard]] double IsotropicThroatRadius() const { return 0.5 * family_.GetParams().b0; }

    const MorrisThorneFamily& SphericalFamily() const { return family_; }

  private:
    // The renderer represents the asymptotically normalised zero-tidal Ellis
    // member, not the wider spherical comparison family.  Keeping a separate
    // configuration surface prevents the spherical family's coordinate-time
    // rescaling and non-Ellis shape parameters from becoming live capabilities.
    Config config_;
    MorrisThorneFamily family_;
};

inline bool MorrisThorneCartesian::IsValidEvent(const Tensor<double, 4>& pos) const {
    if (family_.GetParams().shape_type != WormholeShapeType::Ellis) return false;
    for (int component = 0; component < 4; ++component) {
        if (!std::isfinite(pos(component))) return false;
    }
    const double x = pos(1), y = pos(2), z = pos(3);
    // rho=0 is the other asymptotic infinity of the isotropic chart, not a
    // finite event.  The regular throat and all positive-rho second-sheet
    // stages are represented exactly.
    return x * x + y * y + z * z > 0.0;
}

inline void MorrisThorneCartesian::Evaluate(const Tensor<double, 4>& pos, Metric4d& g,
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

    const double x = pos(1), y = pos(2), z = pos(3);
    const double rho2 = x * x + y * y + z * z;
    const double rho = std::sqrt(rho2);
    const double b0 = family_.GetParams().b0;
    const double q = (b0 * b0) / (4.0 * rho2);
    const double conformal_base = 1.0 + q;
    const double conformal = conformal_base * conformal_base;
    g.Zero();
    g(0, 0) = Dual<double>(-1.0);
    for (int i = 0; i < 3; ++i) {
        g(i + 1, i + 1) = Dual<double>(conformal);
    }

    dg.Zero();
    // dC/drho = -4 q (1+q)/rho for C=(1+q)^2 and
    // d_k C = (dC/drho) x_k/rho.
    const double d_conformal_drho = -4.0 * q * conformal_base / rho;
    const double coordinates[3] = {x, y, z};
    for (int k = 0; k < 3; ++k) {
        for (int i = 0; i < 3; ++i) {
            dg(k + 1, i + 1, i + 1) = Dual<double>(d_conformal_drho * coordinates[k] / rho);
        }
    }
}

inline bool MorrisThorneCartesian::InverseMetric(const Tensor<double, 4>& pos,
                                                 Metric4d& g_inv) const {
    const bool represented = IsValidEvent(pos);
    SIRIUS_PRE(represented);
    if (!represented) {
        const Dual<double> nan(std::numeric_limits<double>::quiet_NaN());
        for (int mu = 0; mu < 4; ++mu)
            for (int nu = 0; nu < 4; ++nu) g_inv(mu, nu) = nan;
        return false;
    }

    const double x = pos(1), y = pos(2), z = pos(3);
    const double rho2 = x * x + y * y + z * z;
    const double b0 = family_.GetParams().b0;
    const double conformal_base = 1.0 + (b0 * b0) / (4.0 * rho2);
    const double inverse_conformal = 1.0 / (conformal_base * conformal_base);
    g_inv.Zero();
    g_inv(0, 0) = Dual<double>(-1.0);
    for (int i = 0; i < 3; ++i) {
        g_inv(i + 1, i + 1) = Dual<double>(inverse_conformal);
    }
    return true;
}

inline void MorrisThorneFamily::Evaluate(const Tensor<double, 4>& pos, Metric4d& g,
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
    double r = pos(1);
    double theta = pos(2);
    // The phi coordinate is not needed for a spherically symmetric evaluation.

    double sin_theta = std::sin(theta);
    double sin2 = sin_theta * sin_theta;
    double r2 = r * r;

    double b = ShapeFunction(r);
    double db_dr = ShapeFunctionDerivative(r);
    double Phi = RedshiftFunction(r);
    double dPhi_dr = RedshiftFunctionDerivative(r);

    double exp2Phi = std::exp(2.0 * Phi);
    double one_minus_b_over_r = 1.0 - b / r;

    g.Zero();

    // g_tt = -e^{2 Phi}.
    g(0, 0) = Dual<double>(-exp2Phi);

    // g_rr = 1/(1 - b/r).
    g(1, 1) = Dual<double>(1.0 / one_minus_b_over_r);

    // g_theta_theta = r^2.
    g(2, 2) = Dual<double>(r2);

    // g_phi_phi = r^2 sin^2 theta.
    g(3, 3) = Dual<double>(r2 * sin2);

    dg.Zero();

    // d g_tt / dr = -2 dPhi/dr e^{2 Phi}.
    dg(1, 0, 0) = Dual<double>(-2.0 * dPhi_dr * exp2Phi);

    // d g_rr / dr = (1 - b/r)^{-2} (db/dr / r - b/r^2).
    double d_one_minus = db_dr / r - b / r2;
    dg(1, 1, 1) = Dual<double>(-d_one_minus / (one_minus_b_over_r * one_minus_b_over_r));

    // d g_theta_theta / dr = 2r.
    dg(1, 2, 2) = Dual<double>(2.0 * r);

    // d g_phi_phi / dr = 2r sin^2 theta.
    dg(1, 3, 3) = Dual<double>(2.0 * r * sin2);

    // d g_phi_phi / dtheta = 2 r^2 sin theta cos theta.
    double cos_theta = std::cos(theta);
    dg(2, 3, 3) = Dual<double>(2.0 * r2 * sin_theta * cos_theta);
}

}  // namespace sirius::core
