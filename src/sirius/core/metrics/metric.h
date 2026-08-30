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
#include <limits>
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

    // Set a known finite parameter in its declared [min, max] range. Unknown or
    // out-of-domain requests violate the concrete family's contract rather than
    // being ignored or rewritten as a different spacetime.
    virtual void SetParameter(const std::string& key, double value) = 0;

    // Human-readable metric name.
    virtual const char* GetName() const = 0;

    // Whether a finite chart event belongs to the concrete metric's numerical
    // domain.  All-chart families inherit the finite-coordinate predicate;
    // families with singular sheets or chart boundaries tighten it.  Adaptive
    // integrators consult this before evaluating any Runge-Kutta stage so an
    // unrepresented event is rejected, never replaced by a nearby metric.
    virtual bool IsValidEvent(const Tensor<double, 4>& pos) const {
        for (int component = 0; component < 4; ++component) {
            if (!std::isfinite(pos(component))) return false;
        }
        return true;
    }

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

// Lorentzian signature (-,+,+,+) from the inertia of the full symmetric
// tensor. Diagonal signs are not a signature test when time-space cross terms
// are present (for example inside a Kerr ergoregion), and accepting g00 == 0
// would allow genuinely degenerate tensors. A Jacobi eigensolve is small and
// deterministic for this fixed 4x4 diagnostic.
inline bool CheckLorentzianSignature(const Metric4d& g) {
    double a[4][4] = {};
    double scale = 0.0;
    for (int row = 0; row < 4; ++row) {
        for (int col = 0; col < 4; ++col) {
            const double value = g(row, col).real;
            if (!std::isfinite(value)) return false;
            a[row][col] = value;
            scale = std::max(scale, std::abs(value));
        }
    }
    if (scale == 0.0) return false;

    const double tolerance = 64.0 * std::numeric_limits<double>::epsilon() * scale;
    for (int row = 0; row < 4; ++row) {
        for (int col = row + 1; col < 4; ++col) {
            if (std::abs(a[row][col] - a[col][row]) > tolerance) return false;
            const double symmetric = 0.5 * (a[row][col] + a[col][row]);
            a[row][col] = symmetric;
            a[col][row] = symmetric;
        }
    }

    for (int iteration = 0; iteration < 64; ++iteration) {
        int p = 0;
        int q = 1;
        double largest = 0.0;
        for (int row = 0; row < 4; ++row) {
            for (int col = row + 1; col < 4; ++col) {
                const double candidate = std::abs(a[row][col]);
                if (candidate > largest) {
                    largest = candidate;
                    p = row;
                    q = col;
                }
            }
        }
        if (largest <= tolerance) break;

        const double apq = a[p][q];
        const double tau = (a[q][q] - a[p][p]) / (2.0 * apq);
        const double t = std::copysign(1.0 / (std::abs(tau) + std::sqrt(1.0 + tau * tau)), tau);
        const double c = 1.0 / std::sqrt(1.0 + t * t);
        const double s = t * c;
        const double app = a[p][p];
        const double aqq = a[q][q];

        for (int index = 0; index < 4; ++index) {
            if (index == p || index == q) continue;
            const double aip = a[index][p];
            const double aiq = a[index][q];
            a[index][p] = c * aip - s * aiq;
            a[p][index] = a[index][p];
            a[index][q] = s * aip + c * aiq;
            a[q][index] = a[index][q];
        }
        a[p][p] = app - t * apq;
        a[q][q] = aqq + t * apq;
        a[p][q] = 0.0;
        a[q][p] = 0.0;
    }

    int negative = 0;
    int positive = 0;
    for (int index = 0; index < 4; ++index) {
        if (a[index][index] < -tolerance) {
            ++negative;
        } else if (a[index][index] > tolerance) {
            ++positive;
        }
    }
    return negative == 1 && positive == 3;
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

// Full metric validation at a point: finite, symmetric within tolerance, and
// non-degenerate with exactly one timelike direction.
inline bool ValidateMetricAtPoint(IMetric* metric, const Tensor<double, 4>& pos,
                                  double tolerance = 1e-10) {
    Metric4d g;
    Tensor<Dual<double>, 4, 4, 4> dg;
    metric->Evaluate(pos, g, dg);

    if (!CheckFinite(g)) return false;

    if (CheckSymmetry(g) > tolerance) return false;

    return CheckLorentzianSignature(g);
}

}  // namespace metric_validation

}  // namespace sirius::core
