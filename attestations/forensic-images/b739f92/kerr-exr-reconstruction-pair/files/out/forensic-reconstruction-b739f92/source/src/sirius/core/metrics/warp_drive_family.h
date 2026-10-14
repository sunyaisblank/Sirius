#pragma once

// Unified Alcubierre-class warp drive family in Cartesian coordinates. Ported
// from PHMT102A.h.
//
//   ds^2 = -dt^2 + (dx - v_s(t) f(r_s) dt)^2 + dy^2 + dz^2
// with r_s the distance from the bubble centre and f the shape function
//   f(r_s) = [tanh(sigma (r_s + R)) - tanh(sigma (r_s - R))] / [2 tanh(sigma R)].
// Already Cartesian, so there are no poles.
// Reference: Alcubierre, Class. Quantum Grav. 11, L73 (1994).

#include "sirius/core/metrics/metric.h"

#include <algorithm>
#include <cmath>

namespace sirius::core {

// Warp drive family parameters.
struct WarpDriveParams {
    double vs = 1.0;     // Bubble velocity (may exceed 1 for superluminal).
    double sigma = 8.0;  // Wall thickness parameter (larger = sharper).
    double R = 1.0;      // Bubble radius.
    double xs = 0.0;     // Current x position of the bubble centre.
    double ys = 0.0;     // Current y position.
    double zs = 0.0;     // Current z position.

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

// Warp drive family metric.
class WarpDriveFamily : public IMetric {
  public:
    WarpDriveFamily();
    explicit WarpDriveFamily(const WarpDriveParams& params);

    void Evaluate(const Tensor<double, 4>& pos, Metric4d& g,
                  Tensor<Dual<double>, 4, 4, 4>& dg) override;

    const Config& GetParameters() const override { return config_; }
    void SetParameter(const std::string& key, double value) override;
    const char* GetName() const override { return "Alcubierre Warp Drive"; }

    void SetParams(const WarpDriveParams& params);
    WarpDriveParams GetParams() const;

    // Shape function f(r_s) and its derivative.
    double ShapeFunction(double rs) const;
    double ShapeFunctionDerivative(double rs) const;

    // Reposition the bubble centre (for animation between renders).
    void UpdateBubblePosition(double t);

  private:
    Config config_;
    WarpDriveParams params_;
};

inline WarpDriveFamily::WarpDriveFamily() {
    config_["velocity"] = {1.0, 0.0, 10.0};
    config_["sigma"] = {8.0, 1.0, 50.0};
    config_["radius"] = {1.0, 0.1, 100.0};
    params_ = WarpDriveParams::Alcubierre(1.0, 1.0);
}

inline WarpDriveFamily::WarpDriveFamily(const WarpDriveParams& params) : WarpDriveFamily() {
    SetParams(params);
}

inline void WarpDriveFamily::SetParams(const WarpDriveParams& params) {
    params_ = params;
    config_["velocity"].value = params.vs;
    config_["sigma"].value = params.sigma;
    config_["radius"].value = params.R;
}

inline WarpDriveParams WarpDriveFamily::GetParams() const {
    return params_;
}

inline void WarpDriveFamily::SetParameter(const std::string& key, double value) {
    if (config_.find(key) != config_.end()) {
        config_[key].value = std::clamp(value, config_[key].min, config_[key].max);
    }
    params_.vs = config_["velocity"].value;
    params_.sigma = config_["sigma"].value;
    params_.R = config_["radius"].value;
}

inline double WarpDriveFamily::ShapeFunction(double rs) const {
    double sigma = params_.sigma;
    double R = params_.R;

    // Alcubierre smooth top-hat.
    double tanh_plus = std::tanh(sigma * (rs + R));
    double tanh_minus = std::tanh(sigma * (rs - R));
    double tanh_R = std::tanh(sigma * R);

    return (tanh_plus - tanh_minus) / (2.0 * tanh_R);
}

inline double WarpDriveFamily::ShapeFunctionDerivative(double rs) const {
    double sigma = params_.sigma;
    double R = params_.R;

    double sech2_plus = 1.0 / std::cosh(sigma * (rs + R));
    sech2_plus *= sech2_plus;
    double sech2_minus = 1.0 / std::cosh(sigma * (rs - R));
    sech2_minus *= sech2_minus;
    double tanh_R = std::tanh(sigma * R);

    return sigma * (sech2_plus - sech2_minus) / (2.0 * tanh_R);
}

inline void WarpDriveFamily::UpdateBubblePosition(double t) {
    // Sets the reference (t = 0) bubble position. Evaluate() adds vs*t on top of
    // this reference, so calling it per integration step would double-count the
    // motion; it exists for repositioning the bubble between renders.
    params_.xs = params_.vs * t;
}

inline void WarpDriveFamily::Evaluate(const Tensor<double, 4>& pos, Metric4d& g,
                                      Tensor<Dual<double>, 4, 4, 4>& dg) {
    double t = pos(0);
    double x = pos(1);
    double y = pos(2);
    double z = pos(3);

    double vs = params_.vs;
    // The bubble moves along x at coordinate speed vs: xs(t) = xs0 + vs t.
    // Previously Evaluate() ignored t and UpdateBubblePosition() was never called
    // during integration, so the metric was silently static and its time
    // derivative silently zero: the geodesic equation was solved for a different
    // spacetime than the one advertised.
    double xs = params_.xs + vs * t;
    double ys = params_.ys;
    double zs = params_.zs;

    // Distance from the bubble centre.
    double dx = x - xs;
    double dy = y - ys;
    double dz = z - zs;
    double rs = std::sqrt(dx * dx + dy * dy + dz * dz);
    rs = std::max(rs, 1e-10);

    double f = ShapeFunction(rs);
    double df_drs = ShapeFunctionDerivative(rs);

    // Derivatives of rs with respect to position.
    double drs_dx = dx / rs;
    double drs_dy = dy / rs;
    double drs_dz = dz / rs;

    // Derivatives of f with respect to position.
    double df_dx = df_drs * drs_dx;
    double df_dy = df_drs * drs_dy;
    double df_dz = df_drs * drs_dz;

    g.Zero();

    // ds^2 = -(1 - v_s^2 f^2) dt^2 - 2 v_s f dx dt + dx^2 + dy^2 + dz^2.
    double vsf = vs * f;
    double vsf2 = vsf * vsf;

    g(0, 0) = Dual<double>(-(1.0 - vsf2));  // g_tt
    g(0, 1) = Dual<double>(-vsf);           // g_tx
    g(1, 0) = Dual<double>(-vsf);           // g_xt
    g(1, 1) = Dual<double>(1.0);            // g_xx
    g(2, 2) = Dual<double>(1.0);            // g_yy
    g(3, 3) = Dual<double>(1.0);            // g_zz

    dg.Zero();

    // d g_tt / dx = 2 v_s^2 f df/dx.
    dg(1, 0, 0) = Dual<double>(2.0 * vs * vs * f * df_dx);
    dg(2, 0, 0) = Dual<double>(2.0 * vs * vs * f * df_dy);
    dg(3, 0, 0) = Dual<double>(2.0 * vs * vs * f * df_dz);

    // d g_tx / dx = -v_s df/dx.
    dg(1, 0, 1) = Dual<double>(-vs * df_dx);
    dg(1, 1, 0) = Dual<double>(-vs * df_dx);
    dg(2, 0, 1) = Dual<double>(-vs * df_dy);
    dg(2, 1, 0) = Dual<double>(-vs * df_dy);
    dg(3, 0, 1) = Dual<double>(-vs * df_dz);
    dg(3, 1, 0) = Dual<double>(-vs * df_dz);

    // Time derivatives, exactly, by the chain rule: the metric depends on t only
    // through dx = x - xs(t) with dxs/dt = vs, so d_t g = -vs d_x g.
    for (int mu = 0; mu < 4; mu++) {
        for (int nu = 0; nu < 4; nu++) {
            dg(0, mu, nu) = Dual<double>(-vs * dg(1, mu, nu).real);
        }
    }
}

}  // namespace sirius::core
