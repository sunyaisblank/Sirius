// Double-precision polarisation transport for the validation oracle
// (sirius::oracle): parallel transport of a polarisation vector along a null
// geodesic in the analytic Boyer-Lindquist Kerr metric, and the conserved
// Walker-Penrose complex constant that certifies it. Off the render path; the
// reference the live Kerr-Schild transport (core/polarisation/walker_penrose.h)
// is measured against.
//
// Transport equations (contravariant tangent k^mu and polarisation f^mu,
// second-order geodesic form so k and f share one connection evaluation):
//   dx^mu/dl = k^mu,  dk^mu/dl = -Gamma^mu_ab k^a k^b,
//   df^mu/dl = -Gamma^mu_ab k^a f^b,   with f.k = 0 and f spacelike.
// Parallel transport preserves every inner product, so g_uv f^u f^v and
// g_uv f^u k^v are integration-quality monitors, and the Walker-Penrose
// constant below is the polarisation-specific invariant.
//
// Connection source: KerrMetricD::Christoffel, the single connection authority.
// This workstream originally found that method disagreeing with its own metric
// (finite differences matched only to ~1e-1; the defect was latent because the
// integration path uses dHdq); it has since been rebuilt from the analytic
// metric derivatives and is pinned against finite differences by
// KerrMetricDTest.ChristoffelMatchesFiniteDifferencesOfMetric. KerrChristoffel
// below remains as the transport-local name and simply delegates.
//
// WALKER-PENROSE CONSTANT (Kerr, Petrov type D). For a null geodesic with
// parallel-propagated f (f.k = 0),
//   kappa_WP = (A - i B)(r - i a cos theta)
// is conserved, with the Boyer-Lindquist projections
//   A = (k^t f^r - k^r f^t) + a sin^2(theta) (k^r f^phi - k^phi f^r),
//   B = [(r^2 + a^2)(k^phi f^theta - k^theta f^phi)
//        - a (k^t f^theta - k^theta f^t)] sin(theta).
// Derivation and citation: the imaginary part equals the Killing-Yano
// contraction Im(kappa_WP) = Y_uv k^u f^v, where the Kerr Killing-Yano 2-form
// (Penrose & Floyd 1973; Chandrasekhar 1983, section 61) has non-zero
// Boyer-Lindquist components
//   Y_{tr} = -a cos theta,        Y_{t theta} = a r sin theta,
//   Y_{r phi} = -a^2 cos theta sin^2 theta,  Y_{theta phi} = r (r^2+a^2) sin theta.
// Because a Killing-Yano tensor obeys grad_(lam Y_mu)nu = 0, the derivative
// grad_lam Y_munu is totally antisymmetric; contracted with the symmetric
// k^lam k^mu of a geodesic it vanishes, and the geodesic/transport conditions
// kill the remaining terms, so d/dl (Y_uv k^u f^v) = 0. Expanding Y_uv k^u f^v
// reproduces exactly Im[(A - iB)(r - i a cos theta)] = -(a A cos theta + r B),
// and the real part is the dual Killing-Yano (closed conformal Killing-Yano)
// contraction, likewise conserved; hence the full complex kappa_WP is
// invariant. References: Walker & Penrose (1970), Commun. Math. Phys. 18, 265;
// Chandrasekhar, "The Mathematical Theory of Black Holes" (1983), section 60;
// component form as used for black-hole polarimetry in Connors, Piran & Stark
// (1980), ApJ 235, 224.

#pragma once

#include "sirius/base/contracts.h"
#include "sirius/oracle/kerr_boyer_lindquist.h"
#include "sirius/oracle/transport_types.h"

#include <algorithm>
#include <cmath>
#include <complex>

namespace sirius::oracle {

//==============================================================================
// Tolerances and step defaults (self-contained, mirroring core/constants.h)
//==============================================================================

// Relative conservation the double-precision RK4 transport reaches for the
// Walker-Penrose constant over a full geodesic. The invariant is exact in the
// continuum; a fourth-order integrator on a well-conditioned Kerr geodesic with
// adaptive step control accumulates error to the double round-off floor scaled
// by the path length, and 1e-10 is the oracle-tier conservation bar the
// specification sets for E, L_z and the null constraint (core/constants.h
// geodesic::kNullConditionTolCpu, kHamiltonianTol). Kept local so the oracle
// header depends only on sirius_base, exactly as symplectic_integrator.h does.
inline constexpr double kWalkerPenroseConservationTol = 1e-10;

// Base affine step for the polarised integrator. Small enough that RK4 clears
// the conservation bar above near the photon sphere; the integrator shrinks it
// further where the geodesic bends most (see AdaptiveStepSize).
inline constexpr double kPolarisedTransportStep = 0.01;

// Theta clamp away from the polar axis singularity of the Boyer-Lindquist
// chart; 1e-6 == core/constants.h coordinates::kPoleEpsilon, inlined to keep
// the oracle self-contained (STYLE.md self-containment over shared const).
inline constexpr double kPolarisedPoleClamp = 1e-6;

//==============================================================================
// PolarisedStateD: geodesic tangent plus a transported polarisation vector
//==============================================================================
//
// Unlike GeodesicStateD, which stores the covariant momentum p_mu, this state
// stores the contravariant tangent k^mu = dx^mu/dl directly, because parallel
// transport is stated in contravariant components. Raise with g^uv p_v to
// convert an existing GeodesicStateD tangent.

struct PolarisedStateD {
    Vec4d x;              // Position (t, r, theta, phi) in Boyer-Lindquist.
    Vec4d k;              // Contravariant tangent k^mu = dx^mu/dl.
    Vec4d f;              // Contravariant polarisation f^mu, f.k = 0, spacelike.
    double lambda = 0.0;  // Affine parameter.

    PolarisedStateD() = default;
    PolarisedStateD(const Vec4d& x_, const Vec4d& k_, const Vec4d& f_) : x(x_), k(k_), f(f_) {}
};

//==============================================================================
// Metric-contraction helpers (double, index order t, r, theta, phi)
//==============================================================================

// Inner product g_uv u^u v^v from a covariant metric block.
[[nodiscard]] inline double InnerProductD(const double g[4][4], const Vec4d& u, const Vec4d& v) {
    double s = 0.0;
    for (int mu = 0; mu < 4; ++mu)
        for (int nu = 0; nu < 4; ++nu) s += g[mu][nu] * u[mu] * v[nu];
    return s;
}

// Raise an index: v^mu = g^uv v_nu.
[[nodiscard]] inline Vec4d RaiseIndexD(const double g_inv[4][4], const Vec4d& v_cov) {
    Vec4d out;
    for (int mu = 0; mu < 4; ++mu) {
        double s = 0.0;
        for (int nu = 0; nu < 4; ++nu) s += g_inv[mu][nu] * v_cov[nu];
        out[mu] = s;
    }
    return out;
}

// -Gamma^mu_ab u^a v^b: the connection term shared by the geodesic (u = v = k)
// and the transport (u = k, v = f) right-hand sides.
[[nodiscard]] inline Vec4d NegConnectionContraction(const double Gamma[4][4][4], const Vec4d& u,
                                                    const Vec4d& v) {
    Vec4d out;
    for (int mu = 0; mu < 4; ++mu) {
        double s = 0.0;
        for (int a = 0; a < 4; ++a)
            for (int b = 0; b < 4; ++b) s += Gamma[mu][a][b] * u[a] * v[b];
        out[mu] = -s;
    }
    return out;
}

//==============================================================================
// Analytic Kerr Boyer-Lindquist connection (correct, self-contained)
//==============================================================================

// KerrMetricDerivatives now lives beside the metric it differentiates
// (kerr_boyer_lindquist.h), the single connection source.

// Transport-local name for the connection; delegates to the single authority
// (KerrMetricD::Christoffel, rebuilt from the analytic derivatives and pinned
// against finite differences; see the file header).
inline void KerrChristoffel(const KerrMetricD& metric, const Vec4d& x, double Gamma[4][4][4]) {
    metric.Christoffel(x, Gamma);
}

//==============================================================================
// Walker-Penrose constant (Kerr). Header block above carries the derivation.
//==============================================================================

// kappa_WP = (A - iB)(r - i a cos theta) from the contravariant tangent and
// polarisation of state s at its current Boyer-Lindquist point, spin a.
// Reference: Walker & Penrose (1970); Chandrasekhar (1983), section 60.
[[nodiscard]] inline std::complex<double> WalkerPenroseConstant(const PolarisedStateD& s,
                                                                double a) {
    const double r = s.x.r;
    const double theta = s.x.theta;
    const double sinth = std::sin(theta);
    const double costh = std::cos(theta);
    const double sin2th = sinth * sinth;

    const double kt = s.k.t, kr = s.k.r, kth = s.k.theta, kph = s.k.phi;
    const double ft = s.f.t, fr = s.f.r, fth = s.f.theta, fph = s.f.phi;

    const double A = (kt * fr - kr * ft) + a * sin2th * (kr * fph - kph * fr);
    const double B =
        ((r * r + a * a) * (kph * fth - kth * fph) - a * (kt * fth - kth * ft)) * sinth;

    return std::complex<double>(A, -B) * std::complex<double>(r, -a * costh);
}

// Overload pulling the spin from the metric.
[[nodiscard]] inline std::complex<double> WalkerPenroseConstant(const PolarisedStateD& s,
                                                                const KerrMetricD& metric) {
    return WalkerPenroseConstant(s, metric.spin());
}

//==============================================================================
// Initial-data construction: a null tangent and an orthonormal polarisation
//==============================================================================

// Solve g_uv k^u k^v = 0 for the future-directed k^t given spatial contravariant
// components (k^r, k^theta, k^phi). In Boyer-Lindquist the only t off-diagonal is
// g_t phi, so the quadratic is g_tt (k^t)^2 + 2 g_t phi k^phi k^t + (spatial) = 0.
[[nodiscard]] inline Vec4d MakeNullTangent(const KerrMetricD& metric, const Vec4d& x, double kr,
                                           double kth, double kph) {
    double g[4][4], g_inv[4][4];
    metric.Evaluate(x, g, g_inv);

    const double a_coef = g[0][0];
    const double b_coef = 2.0 * g[0][3] * kph;  // g_tr = g_ttheta = 0 in Kerr BL.
    const double c_coef = g[1][1] * kr * kr + g[2][2] * kth * kth + g[3][3] * kph * kph;

    const double disc = b_coef * b_coef - 4.0 * a_coef * c_coef;
    SIRIUS_PRE(disc >= 0.0);  // A physical null direction exists.
    const double sqrt_disc = std::sqrt(std::max(disc, 0.0));

    // g_tt < 0 outside the ergosphere: the future-directed (k^t > 0) root is
    // (-b - sqrt_disc)/(2 a_coef); pick whichever root is positive to stay
    // future-directed inside the ergosphere as well.
    const double root_plus = (-b_coef + sqrt_disc) / (2.0 * a_coef);
    const double root_minus = (-b_coef - sqrt_disc) / (2.0 * a_coef);
    const double kt = (root_plus > 0.0) ? root_plus : root_minus;

    return Vec4d(kt, kr, kth, kph);
}

// Orthonormalise a trial vector against a null tangent: f = trial - (trial.k /
// e_t.k) e_t with e_t = (1,0,0,0), which gives f.k = 0 because e_t.k = -k_t != 0,
// then scale to g_uv f^u f^v = 1. Precondition: the projected f is spacelike.
// Adding a multiple of k leaves f.k and the Walker-Penrose constant unchanged
// (kappa_WP vanishes when f = k), so this fixes the residual polarisation gauge.
[[nodiscard]] inline Vec4d MakeOrthonormalPolarisation(const KerrMetricD& metric, const Vec4d& x,
                                                       const Vec4d& k, const Vec4d& trial) {
    double g[4][4], g_inv[4][4];
    metric.Evaluate(x, g, g_inv);

    const Vec4d e_t(1.0, 0.0, 0.0, 0.0);
    const double trial_dot_k = InnerProductD(g, trial, k);
    const double et_dot_k = InnerProductD(g, e_t, k);
    SIRIUS_PRE(std::abs(et_dot_k) > 0.0);

    Vec4d f = trial - e_t * (trial_dot_k / et_dot_k);
    const double norm2 = InnerProductD(g, f, f);
    SIRIUS_PRE(norm2 > 0.0);  // Spacelike polarisation.
    return f / std::sqrt(norm2);
}

//==============================================================================
// ParallelTransportStep: one RK4 step of (x, k, f)
//==============================================================================

namespace detail {

struct PolarisedDerivD {
    Vec4d dx, dk, df;
};

// Right-hand side of the coupled transport system at state s. One connection
// evaluation drives both the geodesic (dk) and the polarisation (df), so k and
// f see an identical connection and their inner products transport together.
[[nodiscard]] inline PolarisedDerivD PolarisedRhs(const KerrMetricD& metric,
                                                  const PolarisedStateD& s) {
    double Gamma[4][4][4];
    KerrChristoffel(metric, s.x, Gamma);
    PolarisedDerivD d;
    d.dx = s.k;
    d.dk = NegConnectionContraction(Gamma, s.k, s.k);
    d.df = NegConnectionContraction(Gamma, s.k, s.f);
    return d;
}

// s + factor * derivative, position/tangent/polarisation componentwise.
[[nodiscard]] inline PolarisedStateD Advance(const PolarisedStateD& s, const PolarisedDerivD& d,
                                             double factor) {
    PolarisedStateD out;
    out.x = s.x + d.dx * factor;
    out.k = s.k + d.dk * factor;
    out.f = s.f + d.df * factor;
    out.lambda = s.lambda;
    return out;
}

}  // namespace detail

// Classical fourth-order Runge-Kutta step of size h on the coupled system.
// RK4 rather than the oracle's symplectic leapfrog because the transported f
// has no Hamiltonian structure to preserve; the shared connection evaluation
// keeps k and f consistent, which is what conserves the Walker-Penrose
// constant. Reference: geodesic and Jacobi transport of Kerr null rays, James
// et al. (2015), Appendix; matches the RK4 style of beam_integrator.h.
[[nodiscard]] inline PolarisedStateD ParallelTransportStep(const KerrMetricD& metric,
                                                           const PolarisedStateD& s, double h) {
    const detail::PolarisedDerivD k1 = detail::PolarisedRhs(metric, s);
    const detail::PolarisedDerivD k2 =
        detail::PolarisedRhs(metric, detail::Advance(s, k1, 0.5 * h));
    const detail::PolarisedDerivD k3 =
        detail::PolarisedRhs(metric, detail::Advance(s, k2, 0.5 * h));
    const detail::PolarisedDerivD k4 = detail::PolarisedRhs(metric, detail::Advance(s, k3, h));

    PolarisedStateD out;
    out.x = s.x + (k1.dx + k2.dx * 2.0 + k3.dx * 2.0 + k4.dx) * (h / 6.0);
    out.k = s.k + (k1.dk + k2.dk * 2.0 + k3.dk * 2.0 + k4.dk) * (h / 6.0);
    out.f = s.f + (k1.df + k2.df * 2.0 + k3.df * 2.0 + k4.df) * (h / 6.0);
    out.lambda = s.lambda + h;

    // Keep the next Christoffel evaluation away from the axis singularity; rays
    // integrated here stay well clear of it, so the clamp is inert in practice.
    out.x.theta = std::clamp(out.x.theta, kPolarisedPoleClamp, M_PI - kPolarisedPoleClamp);
    return out;
}

//==============================================================================
// PolarisedGeodesicIntegratorD: advance (x, k, f) to escape or capture
//==============================================================================

class PolarisedGeodesicIntegratorD {
  public:
    struct Config {
        double baseStep = kPolarisedTransportStep;
        double minStep = 1e-5;
        double maxStep = 0.05;
        double escapeRadius = 60.0;   // Units of M; a full deflected geodesic.
        double horizonBuffer = 1.02;  // Terminate at r <= r_+ * this.
        int maxSteps = 400000;
    };

    struct Result {
        PolarisedStateD state;
        int steps = 0;
        bool escaped = false;
        bool captured = false;
    };

    explicit PolarisedGeodesicIntegratorD(const KerrMetricD* metric) : metric_(metric), config_() {
        SIRIUS_PRE(metric != nullptr);
    }

    PolarisedGeodesicIntegratorD(const KerrMetricD* metric, const Config& config)
        : metric_(metric), config_(config) {
        SIRIUS_PRE(metric != nullptr);
    }

    // Step shrinks quadratically as the geodesic approaches the horizon, exactly
    // as the oracle symplectic integrator adapts, concentrating resolution where
    // the connection is strongest (near the photon sphere).
    [[nodiscard]] double AdaptiveStepSize(double r) const {
        const double r_h = metric_->HorizonRadius();
        const double ratio = (r - r_h) / r_h;
        const double h = config_.baseStep * std::min(1.0, ratio * ratio);
        return std::clamp(h, config_.minStep, config_.maxStep);
    }

    [[nodiscard]] Result Integrate(const PolarisedStateD& initial) const {
        Result result;
        PolarisedStateD s = initial;
        const double r_h = metric_->HorizonRadius();

        for (int i = 0; i < config_.maxSteps; ++i) {
            if (s.x.r <= r_h * config_.horizonBuffer) {
                result.captured = true;
                break;
            }
            if (s.x.r > config_.escapeRadius) {
                result.escaped = true;
                break;
            }
            s = ParallelTransportStep(*metric_, s, AdaptiveStepSize(s.x.r));
            ++result.steps;
        }
        result.state = s;
        return result;
    }

  private:
    const KerrMetricD* metric_;
    Config config_;
};

}  // namespace sirius::oracle
