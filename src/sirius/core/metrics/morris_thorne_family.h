#pragma once

// Unified Morris-Thorne traversable wormhole family in Schwarzschild-like
// coordinates. Ported from PHMT101A.h.
//
//   ds^2 = -e^{2 Phi(r)} dt^2 + dr^2/(1 - b(r)/r) + r^2 (dtheta^2 + sin^2 theta dphi^2)
// with Phi(r) the redshift function and b(r) the shape function; traversability
// needs b(b0) = b0 and Phi finite everywhere.
// Reference: Morris & Thorne, Am. J. Phys. 56, 395 (1988).

#include "sirius/core/metrics/metric.h"

#include <algorithm>
#include <cmath>
#include <functional>

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
    config_["throat_radius"] = {1.0, 0.1, 100.0};
    config_["redshift"] = {0.0, -10.0, 10.0};
    params_ = MorrisThorneParams::Ellis(1.0);
}

inline MorrisThorneFamily::MorrisThorneFamily(const MorrisThorneParams& params)
    : MorrisThorneFamily() {
    SetParams(params);
}

inline void MorrisThorneFamily::SetParams(const MorrisThorneParams& params) {
    params_ = params;
    config_["throat_radius"].value = params.b0;
    config_["redshift"].value = params.Phi0;
}

inline MorrisThorneParams MorrisThorneFamily::GetParams() const {
    return params_;
}

inline void MorrisThorneFamily::SetParameter(const std::string& key, double value) {
    if (config_.find(key) != config_.end()) {
        config_[key].value = std::clamp(value, config_[key].min, config_[key].max);
    }
    params_.b0 = config_["throat_radius"].value;
    params_.Phi0 = config_["redshift"].value;
}

inline const char* MorrisThorneFamily::GetName() const {
    switch (params_.shape_type) {
        case WormholeShapeType::Ellis: return "Ellis Drainhole";
        case WormholeShapeType::ZeroTidal: return "Zero-Tidal Wormhole";
        case WormholeShapeType::AbsurdlyBenign: return "Absurdly Benign Wormhole";
        default: return "Morris-Thorne Wormhole";
    }
}

inline double MorrisThorneFamily::ShapeFunction(double r) const {
    double b0 = params_.b0;
    r = std::max(r, b0);  // Ensure r >= b0.

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
            return b0 * b0 / r;  // Fallback to Ellis.
        default:
            return b0 * b0 / r;  // Default to Ellis.
    }
}

inline double MorrisThorneFamily::ShapeFunctionDerivative(double r) const {
    double b0 = params_.b0;
    r = std::max(r, b0);

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
                double b_plus = params_.custom_shape_func(r + h);
                double b_minus = params_.custom_shape_func(r - h);
                return (b_plus - b_minus) / (2.0 * h);
            }
            return -b0 * b0 / (r * r);
        default:
            return -b0 * b0 / (r * r);
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

inline void MorrisThorneFamily::Evaluate(const Tensor<double, 4>& pos, Metric4d& g,
                                         Tensor<Dual<double>, 4, 4, 4>& dg) {
    [[maybe_unused]] double t = pos(0);  // Time coordinate (unused in a static metric).
    double r = pos(1);
    double theta = pos(2);
    // The phi coordinate is not needed for a spherically symmetric evaluation.

    double b0 = params_.b0;
    r = std::max(r, b0 * 1.001);  // Stay just outside the throat.

    double sin_theta = std::sin(theta);
    double sin2 = sin_theta * sin_theta;
    double r2 = r * r;

    double b = ShapeFunction(r);
    double db_dr = ShapeFunctionDerivative(r);
    double Phi = RedshiftFunction(r);
    double dPhi_dr = RedshiftFunctionDerivative(r);

    double exp2Phi = std::exp(2.0 * Phi);
    double one_minus_b_over_r = 1.0 - b / r;
    one_minus_b_over_r = std::max(one_minus_b_over_r, 1e-10);  // Regularise at the throat.

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
