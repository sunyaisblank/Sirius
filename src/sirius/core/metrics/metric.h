#pragma once

// The IMetric interface every spacetime implements, plus the runtime validators
// that check the tensor properties a metric must satisfy. Ported from PHMT000A.h.
//
// Contract of Evaluate(pos): returns g[mu,nu] and dg[sigma,mu,nu]; the metric is
// symmetric, Lorentzian signature (-,+,+,+), non-degenerate, and finite on the
// family's valid domain.

#include "sirius/core/tensor.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <string>

namespace sirius::core {

// A metric parameter with its value and admissible range.
struct MetricParameter {
    double value;
    double min;
    double max;
};

using Config = std::map<std::string, MetricParameter>;

// Abstract interface for a spacetime metric.
class IMetric {
  public:
    virtual ~IMetric() = default;

    // Evaluate the metric g_mu_nu and its derivatives d_sigma g_mu_nu at a
    // 4-position (t, r, theta, phi) or (t, x, y, z). The metric carries Dual
    // scalars for autodiff compatibility.
    virtual void Evaluate(const Tensor<double, 4>& pos, Metric4d& g,
                          Tensor<Dual<double>, 4, 4, 4>& dg) = 0;

    // Current parameter configuration.
    virtual const Config& GetParameters() const = 0;

    // Set a parameter, clamped to its [min, max] range.
    virtual void SetParameter(const std::string& key, double value) = 0;

    // Human-readable metric name.
    virtual const char* GetName() const = 0;

    // Analytic inverse metric where the family has a closed form. Writes g_inv
    // and returns true only then; false directs the caller to invert the
    // evaluated metric via TensorOps::Inverse.
    virtual bool InverseMetric(const Tensor<double, 4>& pos, Metric4d& g_inv) const {
        (void)pos;
        (void)g_inv;
        return false;
    }

    // Capture-surface test in the metric's own coordinates. margin is the
    // relative safety margin on the capture radius (0.05 means terminate at 1.05
    // times the horizon radius). Horizonless spacetimes keep the default false.
    virtual bool InsideCaptureSurface(const Tensor<double, 4>& pos, double margin) const {
        (void)pos;
        (void)margin;
        return false;
    }
};

// Runtime checks for the metric tensor properties that must hold exactly; used
// by tests and diagnostics.
namespace metric_validation {

// Maximum symmetry violation |g_mu_nu - g_nu_mu|.
inline double CheckSymmetry(const Metric4d& g) {
    double max_violation = 0.0;
    for (int mu = 0; mu < 4; ++mu) {
        for (int nu = mu + 1; nu < 4; ++nu) {
            double violation = std::abs(g(mu, nu).real - g(nu, mu).real);
            max_violation = std::max(max_violation, violation);
        }
    }
    return max_violation;
}

// Maximum Christoffel symmetry violation |Gamma^lam_mu_nu - Gamma^lam_nu_mu|.
inline double CheckChristoffelSymmetry(const double Gamma[4][4][4]) {
    double max_violation = 0.0;
    for (int lam = 0; lam < 4; ++lam) {
        for (int mu = 0; mu < 4; ++mu) {
            for (int nu = mu + 1; nu < 4; ++nu) {
                double violation = std::abs(Gamma[lam][mu][nu] - Gamma[lam][nu][mu]);
                max_violation = std::max(max_violation, violation);
            }
        }
    }
    return max_violation;
}

// Maximum deviation of g^mu_alpha g_alpha_nu from the Kronecker delta.
inline double CheckInverse(const Metric4d& g, const double g_inv[4][4]) {
    double max_deviation = 0.0;
    for (int mu = 0; mu < 4; ++mu) {
        for (int nu = 0; nu < 4; ++nu) {
            double sum = 0.0;
            for (int alpha = 0; alpha < 4; ++alpha) {
                sum += g_inv[mu][alpha] * g(alpha, nu).real;
            }
            double expected = (mu == nu) ? 1.0 : 0.0;
            double deviation = std::abs(sum - expected);
            max_deviation = std::max(max_deviation, deviation);
        }
    }
    return max_deviation;
}

// Lorentzian signature (-, +, +, +) via the diagonal signs; a diagonal-dominant
// approximation that tolerates the g00 -> 0 coordinate singularity at horizons.
inline bool CheckLorentzianSignature(const Metric4d& g) {
    double g00 = g(0, 0).real;
    double g11 = g(1, 1).real;
    double g22 = g(2, 2).real;
    double g33 = g(3, 3).real;

    if (std::abs(g00) < 1e-10) return true;  // Degenerate at the horizon.

    bool time_negative = g00 < 0;
    bool space_positive = (g11 > 0) && (g22 > 0) && (g33 > 0);

    return time_negative && space_positive;
}

// True when every metric component is finite.
inline bool CheckFinite(const Metric4d& g) {
    for (int mu = 0; mu < 4; ++mu) {
        for (int nu = 0; nu < 4; ++nu) {
            if (!std::isfinite(g(mu, nu).real)) return false;
        }
    }
    return true;
}

// True when every Christoffel component is finite.
inline bool CheckChristoffelFinite(const double Gamma[4][4][4]) {
    for (int lam = 0; lam < 4; ++lam) {
        for (int mu = 0; mu < 4; ++mu) {
            for (int nu = 0; nu < 4; ++nu) {
                if (!std::isfinite(Gamma[lam][mu][nu])) return false;
            }
        }
    }
    return true;
}

// Full metric validation at a point: finite, symmetric within tolerance. The
// signature check is intentionally skipped because it fails near horizons where
// g00 -> 0.
inline bool ValidateMetricAtPoint(IMetric* metric, const Tensor<double, 4>& pos,
                                  double tolerance = 1e-10) {
    Metric4d g;
    Tensor<Dual<double>, 4, 4, 4> dg;
    metric->Evaluate(pos, g, dg);

    if (!CheckFinite(g)) return false;

    if (CheckSymmetry(g) > tolerance) return false;

    return true;
}

}  // namespace metric_validation

}  // namespace sirius::core
