// Implementation of the geodesic integrator declared in geodesic_integrator.h.
// Ported from PHGD001A.cpp; arithmetic order, casts, and literals bit-identical.
//
// Hamiltonian formulation with covariant momenta makes the continuum null
// constraint H = (1/2) g^mu_nu p_mu p_nu = 0 a conserved quantity. Numerical
// steps are accepted only when both truncation error and relative constraint
// residual are within tolerance; admitted roundoff is projected back onto H=0.
// Hamilton's equations:
//   dx^mu/dlambda = g^mu_nu p_nu,   dp_mu/dlambda = (1/2)(d g_alpha_beta/dx^mu) k^alpha k^beta.

#include "sirius/core/geodesic_integrator.h"

#include "sirius/base/contracts.h"
#include "sirius/core/constants.h"
#include "sirius/core/trace_boundary.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace sirius::core {

// Dormand-Prince RK45 Butcher tableau (nodes c_i, matrix a_ij, weights b_i, b*_i,
// and error coefficients e_i = b_i - b*_i), kept in domain notation.
namespace dp45 {
[[maybe_unused]] constexpr double c2 = 1.0 / 5.0;
[[maybe_unused]] constexpr double c3 = 3.0 / 10.0;
[[maybe_unused]] constexpr double c4 = 4.0 / 5.0;
[[maybe_unused]] constexpr double c5 = 8.0 / 9.0;
[[maybe_unused]] constexpr double c6 = 1.0;
[[maybe_unused]] constexpr double c7 = 1.0;

constexpr double a21 = 1.0 / 5.0;

constexpr double a31 = 3.0 / 40.0;
constexpr double a32 = 9.0 / 40.0;

constexpr double a41 = 44.0 / 45.0;
constexpr double a42 = -56.0 / 15.0;
constexpr double a43 = 32.0 / 9.0;

constexpr double a51 = 19372.0 / 6561.0;
constexpr double a52 = -25360.0 / 2187.0;
constexpr double a53 = 64448.0 / 6561.0;
constexpr double a54 = -212.0 / 729.0;

constexpr double a61 = 9017.0 / 3168.0;
constexpr double a62 = -355.0 / 33.0;
constexpr double a63 = 46732.0 / 5247.0;
constexpr double a64 = 49.0 / 176.0;
constexpr double a65 = -5103.0 / 18656.0;

[[maybe_unused]] constexpr double a71 = 35.0 / 384.0;
// a72 = 0
[[maybe_unused]] constexpr double a73 = 500.0 / 1113.0;
[[maybe_unused]] constexpr double a74 = 125.0 / 192.0;
[[maybe_unused]] constexpr double a75 = -2187.0 / 6784.0;
[[maybe_unused]] constexpr double a76 = 11.0 / 84.0;

// 5th order weights (solution).
constexpr double b1 = 35.0 / 384.0;
// b2 = 0
constexpr double b3 = 500.0 / 1113.0;
constexpr double b4 = 125.0 / 192.0;
constexpr double b5 = -2187.0 / 6784.0;
constexpr double b6 = 11.0 / 84.0;
// b7 = 0

// 4th order weights (error estimate).
constexpr double bs1 = 5179.0 / 57600.0;
// bs2 = 0
constexpr double bs3 = 7571.0 / 16695.0;
constexpr double bs4 = 393.0 / 640.0;
constexpr double bs5 = -92097.0 / 339200.0;
constexpr double bs6 = 187.0 / 2100.0;
constexpr double bs7 = 1.0 / 40.0;

// Error coefficients e_i = b_i - b*_i.
constexpr double e1 = b1 - bs1;  // 71/57600
// e2 = 0
constexpr double e3 = b3 - bs3;  // -71/16695
constexpr double e4 = b4 - bs4;  // 71/1920
constexpr double e5 = b5 - bs5;  // -17253/339200
constexpr double e6 = b6 - bs6;  // 22/525
constexpr double e7 = -bs7;      // -1/40
}  // namespace dp45

// Inverse metric at a position, preferring the family's closed form. Kerr-Schild
// metrics supply g^mu_nu = eta^mu_nu - H l^mu l^nu exactly; anything else falls
// back to the full Cramer inverse. A degenerate metric yields non-finite entries
// there, which HasInvalidState converts into ray termination; no flat-space
// stand-in is fabricated.
static Metric4d InverseAt(IMetric* metric, const Vec4& pos, const Metric4d& g) {
    Metric4d g_inv;
    if (!metric->InverseMetric(pos, g_inv)) {
        g_inv = TensorOps::Inverse(g);
    }
    return g_inv;
}

// dp_mu/dlambda = (1/2)(d g_rho_sigma/dx^mu) k^rho k^sigma, using the identity
// d g^alpha_beta/dx^mu = -g^alpha_rho g^beta_sigma (d g_rho_sigma/dx^mu) with
// k^mu = g^mu_nu p_nu the contravariant velocity.
static Vec4 MomentumDerivative([[maybe_unused]] const Vec4& p, const Vec4& k,
                               const Tensor<Dual<double>, 4, 4, 4>& dg) {
    Vec4 dp;

    for (int mu = 0; mu < 4; mu++) {
        double sum = 0.0;
        for (int rho = 0; rho < 4; rho++) {
            for (int sigma = 0; sigma < 4; sigma++) {
                // Metric derivatives live in the real part of the dual.
                sum += dg(mu, rho, sigma).real * k(rho) * k(sigma);
            }
        }
        dp(mu) = 0.5 * sum;
    }
    return dp;
}

// Forward declarations for helpers defined later.
static bool EvaluateRk45Stage(const Vec4& x, const Vec4& p, IMetric* metric, Vec4& k_x, Vec4& k_p);
static Vec4 ComputeMomentum(const Vec4& velocity, const Metric4d& g);

static bool RejectUnrepresentedStage(Lightray& ray, float minimum_step) {
    if (ray.step_size <= minimum_step) {
        ray.terminated = 3;
    } else {
        ray.step_size = std::max(minimum_step, ray.step_size * 0.5f);
    }
    return false;
}

bool Geodesic::IntegrateStep(Lightray& ray, IMetric* metric, float min_step, float max_step) {
    SIRIUS_PRE(metric != nullptr);

    Vec4 x0 = ray.position;
    float h = ray.step_size;

    if (!metric->IsValidEvent(x0)) {
        ray.terminated = 3;
        return false;
    }

    // Metric and initial momentum.
    Metric4d g0;
    Tensor<Dual<double>, 4, 4, 4> dg0;
    metric->Evaluate(x0, g0, dg0);
    Vec4 p0 = ComputeMomentum(ray.velocity, g0);
    Vec4 k0 = ray.velocity;

    // RK4 stages.
    Vec4 k1_x = k0, k1_p = MomentumDerivative(p0, k0, dg0);

    Vec4 x1 = x0 + k1_x * (0.5f * h), p1 = p0 + k1_p * (0.5f * h);
    Vec4 k2_x, k2_p;
    if (!EvaluateRk45Stage(x1, p1, metric, k2_x, k2_p)) {
        return RejectUnrepresentedStage(ray, min_step);
    }

    Vec4 x2 = x0 + k2_x * (0.5f * h), p2 = p0 + k2_p * (0.5f * h);
    Vec4 k3_x, k3_p;
    if (!EvaluateRk45Stage(x2, p2, metric, k3_x, k3_p)) {
        return RejectUnrepresentedStage(ray, min_step);
    }

    Vec4 x3 = x0 + k3_x * h, p3 = p0 + k3_p * h;
    Vec4 k4_x, k4_p;
    if (!EvaluateRk45Stage(x3, p3, metric, k4_x, k4_p)) {
        return RejectUnrepresentedStage(ray, min_step);
    }

    // Combine RK4.
    Vec4 new_position = x0 + (k1_x + k2_x * 2.0f + k3_x * 2.0f + k4_x) * (h / 6.0f);
    Vec4 new_momentum = p0 + (k1_p + k2_p * 2.0f + k3_p * 2.0f + k4_p) * (h / 6.0f);

    // New velocity.
    Metric4d g_new;
    Tensor<Dual<double>, 4, 4, 4> dg_new;
    if (!metric->IsValidEvent(new_position)) {
        return RejectUnrepresentedStage(ray, min_step);
    }
    metric->Evaluate(new_position, g_new, dg_new);
    Vec4 new_velocity = TensorOps::RaiseIndex(new_momentum, InverseAt(metric, new_position, g_new));

    // Adaptive step control.
    const float velocity_change = static_cast<float>((new_velocity - k0).Length());
    const float position_change = static_cast<float>((new_position - x0).Length());
    const float target_velocity_change = 0.01f, max_position_change = 0.1f;

    if (velocity_change > target_velocity_change * 2.0f || position_change > max_position_change) {
        ray.step_size = std::max(ray.step_size * 0.5f, min_step);
        if (ray.step_size <= min_step) {
            ray.terminated = 5;
            return false;
        }
        return false;
    }
    if (velocity_change < target_velocity_change * 0.5f && ray.step_size < max_step) {
        ray.step_size = std::min(ray.step_size * 1.2f, max_step);
    }

    // Update ray state.
    ray.position = new_position;
    ray.velocity = new_velocity;
    ray.acceleration = CalculateAcceleration(new_velocity, new_position, metric);
    ray.proper_time += h;
    ray.coordinate_time += static_cast<float>(h * std::abs(new_velocity(0)));
    return true;
}

Vec4 Geodesic::CalculateAcceleration(const Vec4& velocity, const Vec4& position, IMetric* metric) {
    Metric4d g;
    Tensor<Dual<double>, 4, 4, 4> dg;
    metric->Evaluate(position, g, dg);

    // Direct acceleration bypasses Christoffel construction (~1.5x faster).
    return TensorOps::GeodesicAccelerationDirect(velocity, InverseAt(metric, position, g), dg);
}

bool Geodesic::CheckTermination(const Lightray& ray, IMetric* metric) {
    // Positions are Cartesian Kerr-Schild (t, x, y, z); the earlier version read
    // position(1) as a spherical radius and compared it to a hardcoded
    // Schwarzschild horizon, which was wrong on both counts.
    using namespace constants::termination;

    const double x = ray.position(1);
    const double y = ray.position(2);
    const double z = ray.position(3);
    const double R = std::sqrt(x * x + y * y + z * z);

    // Outward radial rate d(R)/dlambda = (x vx + y vy + z vz) / R.
    const double dR_dlambda =
        (x * ray.velocity(1) + y * ray.velocity(2) + z * ray.velocity(3)) / std::max(R, 1e-12);

    // Escape: far away and moving outward, the geodesic is essentially straight.
    if (R > kEscapeRadius && dR_dlambda > 0.0) return true;

    // Unconditional background hit.
    if (R > kBackgroundRadius) return true;

    // Capture: the metric family decides in its own coordinates (exact horizon
    // for Kerr-Schild; horizonless spacetimes never report capture).
    if (metric->InsideCaptureSurface(ray.position, kCaptureMargin)) return true;

    // Affine-parameter budget against runaway integration.
    if (ray.proper_time > kMaxAffineParameter) return true;

    // A vanishing velocity means the ray is numerically stuck.
    if (ray.velocity.Length() < kStalledVelocity) return true;

    return false;
}

float Geodesic::CalculateRedshift(const Lightray& ray, const ObserverState& observer,
                                  IMetric* metric) {
    SIRIUS_PRE(metric != nullptr);
    SIRIUS_PRE(observer.is_timelike);
    SIRIUS_PRE(std::isfinite(ray.ku_uobsu) && ray.ku_uobsu > 0.0f);
    Metric4d g;
    Tensor<Dual<double>, 4, 4, 4> dg;
    metric->Evaluate(ray.position, g, dg);

    Vec4 observer_lower = TensorOps::LowerIndex(observer.velocity, g);

    double dot_product = 0.0;
    for (int mu = 0; mu < 4; mu++) {
        dot_product += ray.velocity(mu) * observer_lower(mu);
    }

    // Frequency is positive for either affine orientation. Render rays are
    // past-directed camera-to-source tangents; forward transport rays may be
    // future-directed. Their physical frequency differs only by the sign of k.
    const double measured_frequency = std::abs(dot_product);
    SIRIUS_ASSERT(std::isfinite(measured_frequency) && measured_frequency > 0.0);
    // 1 + z = lambda_obs/lambda_emit = nu_emit/nu_obs.
    return static_cast<float>(measured_frequency / ray.ku_uobsu - 1.0);
}

ObserverState Geodesic::CreateObserver(const Vec4& position, const Vec4& velocity,
                                       IMetric* metric) {
    SIRIUS_PRE(metric != nullptr);
    ObserverState observer;
    observer.position = position;
    observer.velocity = velocity;

    Metric4d g;
    Tensor<Dual<double>, 4, 4, 4> dg;
    metric->Evaluate(position, g, dg);
    double velocity_norm = TensorOps::InnerProduct(velocity, velocity, g);
    observer.is_timelike = std::isfinite(velocity_norm) && velocity_norm < 0.0;
    SIRIUS_PRE(observer.is_timelike);

    double normalization = 1.0 / std::sqrt(-velocity_norm);
    observer.velocity *= normalization;

    CalculateTetrads(observer, metric);

    return observer;
}

void Geodesic::CalculateTetrads(ObserverState& observer, IMetric* metric) {
    SIRIUS_PRE(metric != nullptr);
    SIRIUS_PRE(observer.is_timelike);
    Metric4d g;
    Tensor<Dual<double>, 4, 4, 4> dg;
    metric->Evaluate(observer.position, g, dg);

    observer.e0 = observer.velocity;

    const auto dot = [&g](const Vec4& lhs, const Vec4& rhs) {
        return TensorOps::InnerProduct(lhs, rhs, g);
    };
    const auto normalise_spacelike = [&dot](Vec4& vector) {
        const double norm_squared = dot(vector, vector);
        SIRIUS_ASSERT(std::isfinite(norm_squared) && norm_squared > 1.0e-20);
        vector /= std::sqrt(norm_squared);
    };

    // Lorentzian Gram-Schmidt.  Since g(e0,e0)=-1, removing the
    // timelike projection is v + g(v,e0)e0 (the previous subtraction had
    // the wrong sign).  Each spacelike vector is normalised before it is used
    // to project the next seed.
    observer.e1 = Vec4();
    observer.e1(1) = 1.0;
    observer.e1 += observer.e0 * dot(observer.e1, observer.e0);
    normalise_spacelike(observer.e1);

    observer.e2 = Vec4();
    observer.e2(2) = 1.0;
    observer.e2 += observer.e0 * dot(observer.e2, observer.e0);
    observer.e2 -= observer.e1 * dot(observer.e2, observer.e1);
    normalise_spacelike(observer.e2);

    observer.e3 = Vec4();
    observer.e3(3) = 1.0;
    observer.e3 += observer.e0 * dot(observer.e3, observer.e0);
    observer.e3 -= observer.e1 * dot(observer.e3, observer.e1);
    observer.e3 -= observer.e2 * dot(observer.e3, observer.e2);
    normalise_spacelike(observer.e3);
}

IntegratorConfig Geodesic::GetDefaultConfig() {
    IntegratorConfig config;
    config.abs_tolerance = 1e-6f;
    config.rel_tolerance = 1e-6f;
    config.min_step = 1e-6f;
    config.max_step = 0.1f;
    config.initial_step = 0.01f;
    config.safety_factor = 0.9f;
    config.step_grow_max = 2.0f;
    config.step_shrink_min = 0.1f;
    return config;
}

float Geodesic::ComputeOptimalStep(float h, float error, float tolerance,
                                   const IntegratorConfig& config) {
    SIRIUS_PRE(IsRepresentedIntegratorStepControl(config));
    SIRIUS_PRE(std::isfinite(h) && h >= config.min_step && h <= config.max_step);
    SIRIUS_PRE(std::isfinite(error) && error >= 0.0f);
    SIRIUS_PRE(std::isfinite(tolerance) && tolerance > 0.0f);
    if (error < 1e-15f) {
        return std::min(h * config.step_grow_max, config.max_step);
    }

    float ratio = tolerance / error;
    float factor = config.safety_factor * std::pow(ratio, 0.2f);
    factor = std::max(config.step_shrink_min, std::min(config.step_grow_max, factor));

    float new_step = h * factor;
    return std::max(config.min_step, std::min(config.max_step, new_step));
}

// Covariant momentum p_mu = g_mu_nu k^nu.
static Vec4 ComputeMomentum(const Vec4& velocity, const Metric4d& g) {
    Vec4 p;
    for (int mu = 0; mu < 4; mu++) {
        p(mu) = 0;
        for (int nu = 0; nu < 4; nu++) {
            p(mu) += g(mu, nu).real * velocity(nu);
        }
    }
    return p;
}

// Evaluate one RK45 stage, returning the position and momentum derivatives.
static bool EvaluateRk45Stage(const Vec4& x, const Vec4& p, IMetric* metric, Vec4& k_x, Vec4& k_p) {
    if (!metric->IsValidEvent(x)) return false;
    Metric4d g;
    Tensor<Dual<double>, 4, 4, 4> dg;
    metric->Evaluate(x, g, dg);
    k_x = TensorOps::RaiseIndex(p, InverseAt(metric, x, g));
    k_p = MomentumDerivative(p, k_x, dg);
    return true;
}

// RK45 error norm.
static float ComputeRk45ErrorNorm(const Vec4& error_x, const Vec4& error_p,
                                  const Vec4& old_position, const Vec4& old_momentum,
                                  const Vec4& new_position, const Vec4& new_momentum,
                                  const IntegratorConfig& config) {
    double error_norm = 0.0;
    for (int i = 0; i < 4; i++) {
        const double position_scale =
            config.abs_tolerance +
            config.rel_tolerance * std::max(std::abs(old_position(i)), std::abs(new_position(i)));
        const double momentum_scale =
            config.abs_tolerance +
            config.rel_tolerance * std::max(std::abs(old_momentum(i)), std::abs(new_momentum(i)));
        const double position_error = std::abs(error_x(i)) / position_scale;
        const double momentum_error = std::abs(error_p(i)) / momentum_scale;
        error_norm += position_error * position_error + momentum_error * momentum_error;
    }
    return static_cast<float>(std::sqrt(error_norm / 8.0));
}

// NaN/Inf check on the ray state.
static bool HasInvalidState(const Vec4& position, const Vec4& velocity) {
    for (int i = 0; i < 4; i++) {
        if (std::isnan(position(i)) || std::isinf(position(i)) || std::isnan(velocity(i)) ||
            std::isinf(velocity(i))) {
            return true;
        }
    }
    return false;
}

// Scale-free null residual. Affine reparameterisation k -> Ck multiplies both
// numerator and denominator by C^2, so admission does not depend on the chosen
// photon-frequency normalisation or large coordinate components near a chart
// boundary.
static double RelativeNullResidual(const Vec4& velocity, const Metric4d& metric) {
    double contraction = 0.0;
    double absolute_scale = 0.0;
    for (int mu = 0; mu < 4; ++mu) {
        for (int nu = 0; nu < 4; ++nu) {
            const double term = metric(mu, nu).real * velocity(mu) * velocity(nu);
            contraction += term;
            absolute_scale += std::abs(term);
        }
    }
    if (!std::isfinite(contraction) || !std::isfinite(absolute_scale) ||
        absolute_scale <= std::numeric_limits<double>::min()) {
        return std::numeric_limits<double>::infinity();
    }
    return std::abs(contraction) / absolute_scale;
}

namespace {

bool FiniteVector(const Vec4& value) {
    for (int component = 0; component < 4; ++component)
        if (!std::isfinite(value(component))) return false;
    return true;
}

bool ValidCoupledControl(const Rk45CoupledState& state) {
    if (!std::isfinite(state.length_scale) || !(state.length_scale > 0.0) ||
        !std::isfinite(state.frequency_scale) || !(state.frequency_scale > 0.0) ||
        !std::isfinite(state.tolerance) || !(state.tolerance > 0.0))
        return false;
    for (const auto& column : state.variations)
        if (!FiniteVector(column.displacement) || !FiniteVector(column.derivative)) return false;
    return true;
}

struct PhaseVariation {
    Vec4 x;
    Vec4 p;
};
using PhaseVariations = std::array<PhaseVariation, 4>;

Vec4 CoordinateVariation(const GeodesicVariation& variation, const Vec4& tangent,
                         const ChristoffelSymbols& connection) {
    Vec4 coordinate = variation.derivative;
    for (int mu = 0; mu < 4; ++mu)
        for (int a = 0; a < 4; ++a)
            for (int b = 0; b < 4; ++b)
                coordinate(mu) -=
                    connection.gamma(mu, a, b).real * tangent(a) * variation.displacement(b);
    return coordinate;
}

// Differentiate the same Hamiltonian stage as the central integrator. The
// Hessian is a fourth-order derivative of the metric's own first derivatives.
// Samples must be finite, distinct and in the concrete metric chart.
bool EvaluateVariationStage(IMetric& metric, const Vec4& position, const Vec4& momentum,
                            const PhaseVariations& columns, PhaseVariations& rhs,
                            Rk45CoupledState& control) {
    ++control.variation_stages;
    if (!metric.IsValidEvent(position)) return false;
    Metric4d g;
    Tensor<Dual<double>, 4, 4, 4> dg;
    metric.Evaluate(position, g, dg);
    ++control.variation_metric_evaluations;
    const auto inverse = InverseAt(&metric, position, g);
    const Vec4 tangent = TensorOps::RaiseIndex(momentum, inverse);
    if (!FiniteVector(tangent)) return false;
    double second[4][4][4][4]{};
    constexpr std::array<int, 4> offsets{-2, -1, 1, 2};
    const double radius = std::hypot(position(1), position(2), position(3));
    const double local_scale = radius > 0.0 ? radius : control.length_scale;
    for (int axis = control.stationary ? 1 : 0; axis < 4; ++axis) {
        double spacing = 2.5e-4 * local_scale;
        std::array<Vec4, 4> nodes;
        bool represented = false;
        for (int retry = 0; retry < std::numeric_limits<double>::digits; ++retry) {
            represented = std::isfinite(spacing) && spacing > 0.0;
            for (std::size_t node = 0; node < nodes.size(); ++node) {
                nodes[node] = position;
                nodes[node](axis) += offsets[node] * spacing;
                if (nodes[node](axis) == position(axis) || !FiniteVector(nodes[node]) ||
                    !metric.IsValidEvent(nodes[node]))
                    represented = false;
                for (std::size_t previous = 0; previous < node; ++previous)
                    if (nodes[node](axis) == nodes[previous](axis)) represented = false;
            }
            if (represented) break;
            spacing *= 0.5;
        }
        if (!represented) return false;
        std::array<Tensor<Dual<double>, 4, 4, 4>, 4> samples;
        for (std::size_t node = 0; node < nodes.size(); ++node) {
            Metric4d unused;
            metric.Evaluate(nodes[node], unused, samples[node]);
            ++control.variation_metric_evaluations;
        }
        for (int mu = 0; mu < 4; ++mu)
            for (int a = 0; a < 4; ++a)
                for (int b = 0; b < 4; ++b) {
                    second[axis][mu][a][b] =
                        ((samples[0](mu, a, b).real - samples[3](mu, a, b).real) +
                         8.0 * (samples[2](mu, a, b).real - samples[1](mu, a, b).real)) /
                        (12.0 * spacing);
                    if (!std::isfinite(second[axis][mu][a][b])) return false;
                }
    }
    for (std::size_t column = 0; column < columns.size(); ++column) {
        Vec4 covector = columns[column].p;
        for (int a = 0; a < 4; ++a)
            for (int b = 0; b < 4; ++b)
                for (int axis = 0; axis < 4; ++axis)
                    covector(a) -= dg(axis, a, b).real * columns[column].x(axis) * tangent(b);
        rhs[column].x = TensorOps::RaiseIndex(covector, inverse);
        for (int mu = 0; mu < 4; ++mu) {
            double value = 0.0;
            for (int a = 0; a < 4; ++a)
                for (int b = 0; b < 4; ++b) {
                    value += dg(mu, a, b).real * tangent(a) * rhs[column].x(b);
                    for (int axis = 0; axis < 4; ++axis)
                        value += 0.5 * second[axis][mu][a][b] * columns[column].x(axis) *
                                 tangent(a) * tangent(b);
                }
            rhs[column].p(mu) = value;
        }
        if (!FiniteVector(rhs[column].x) || !FiniteVector(rhs[column].p)) return false;
    }
    return true;
}

bool VariationPair(IMetric& metric, const std::array<Vec4, 7>& positions,
                   const std::array<Vec4, 7>& momenta, const Vec4& initial_tangent, double interval,
                   Rk45CoupledState& control, PhaseVariations& fifth, PhaseVariations& fourth,
                   CoupledSegmentIncrement& fifth_increment,
                   CoupledSegmentIncrement& fourth_increment) {
    using namespace dp45;
    constexpr double weights[7][7] = {{},
                                      {a21},
                                      {a31, a32},
                                      {a41, a42, a43},
                                      {a51, a52, a53, a54},
                                      {a61, a62, a63, a64, a65},
                                      {b1, 0.0, b3, b4, b5, b6}};
    constexpr std::array<double, 7> errors{e1, 0.0, e3, e4, e5, e6, e7};
    Metric4d g;
    Tensor<Dual<double>, 4, 4, 4> dg;
    metric.Evaluate(positions[0], g, dg);
    ++control.variation_metric_evaluations;
    const auto connection = TensorOps::Christoffel(g, dg);
    PhaseVariations initial;
    for (std::size_t column = 0; column < initial.size(); ++column) {
        initial[column].x = control.variations[column].displacement;
        const Vec4 coordinate =
            CoordinateVariation(control.variations[column], initial_tangent, connection);
        for (int a = 0; a < 4; ++a) {
            initial[column].p(a) = 0.0;
            for (int b = 0; b < 4; ++b) {
                initial[column].p(a) += g(a, b).real * coordinate(b);
                for (int axis = 0; axis < 4; ++axis)
                    initial[column].p(a) +=
                        dg(axis, a, b).real * initial[column].x(axis) * initial_tangent(b);
            }
        }
    }
    std::array<PhaseVariations, 7> derivatives;
    PhaseVariations stage;
    for (std::size_t index = 0; index < derivatives.size(); ++index) {
        stage = initial;
        for (std::size_t column = 0; column < stage.size(); ++column)
            for (std::size_t previous = 0; previous < index; ++previous) {
                stage[column].x +=
                    derivatives[previous][column].x * (interval * weights[index][previous]);
                stage[column].p +=
                    derivatives[previous][column].p * (interval * weights[index][previous]);
            }
        if (!EvaluateVariationStage(metric, positions[index], momenta[index], stage,
                                    derivatives[index], control))
            return false;
    }
    fifth = stage;
    fourth = fifth;
    for (std::size_t column = 0; column < fourth.size(); ++column) {
        for (std::size_t index = 0; index < derivatives.size(); ++index) {
            fifth_increment.displacement[column] +=
                derivatives[index][column].x * (interval * weights[6][index]);
            fourth[column].x -= derivatives[index][column].x * (interval * errors[index]);
            fourth[column].p -= derivatives[index][column].p * (interval * errors[index]);
        }
        fourth_increment.displacement[column] = fifth_increment.displacement[column];
        for (std::size_t index = 0; index < derivatives.size(); ++index)
            fourth_increment.displacement[column] -=
                derivatives[index][column].x * (interval * errors[index]);
    }
    return true;
}

struct NullProjection {
    Vec4 tangent;
    int component;
};

}  // namespace

static std::optional<NullProjection> ProjectNullTangentWithBranch(const Vec4& tangent,
                                                                  const Metric4d& metric) {
    for (int component = 0; component < 4; ++component) {
        if (!std::isfinite(tangent(component))) return std::nullopt;
    }

    double metric_scale = 0.0;
    for (int mu = 0; mu < 4; ++mu) {
        for (int nu = 0; nu < 4; ++nu) {
            const double value = metric(mu, nu).real;
            if (!std::isfinite(value)) return std::nullopt;
            metric_scale = std::max(metric_scale, std::abs(value));
        }
    }
    if (!(metric_scale > 0.0)) return std::nullopt;

    std::optional<NullProjection> best;
    double best_correction = std::numeric_limits<double>::infinity();
    const double roundoff = 64.0 * std::numeric_limits<double>::epsilon();
    bool temporal_represented = false;

    // A coordinate-time root need not exist inside a stationary limit. Solve
    // it first to retain the established spatial ray direction; only when it
    // has no represented root, solve the three spatial quadratics and retain
    // their smallest normalized correction. That fallback leaves k^0 unchanged
    // and therefore cannot switch the time orientation of the cone.
    for (int component = 0; component < 4; ++component) {
        if (component > 0 && temporal_represented) break;
        const double a = metric(component, component).real;
        double b = 0.0;
        double c = 0.0;
        for (int i = 0; i < 4; ++i) {
            if (i == component) continue;
            b += (metric(component, i).real + metric(i, component).real) * tangent(i);
            for (int j = 0; j < 4; ++j) {
                if (j == component) continue;
                c += metric(i, j).real * tangent(i) * tangent(j);
            }
        }
        if (!std::isfinite(a) || !std::isfinite(b) || !std::isfinite(c)) continue;

        double roots[2]{};
        int root_count = 0;
        if (a == 0.0) {
            if (b == 0.0) continue;
            roots[root_count++] = -c / b;
        } else {
            const double four_ac = 4.0 * a * c;
            double discriminant = b * b - four_ac;
            const double discriminant_scale = b * b + std::abs(four_ac);
            if (!std::isfinite(discriminant) || discriminant < -roundoff * discriminant_scale) {
                continue;
            }
            discriminant = std::max(discriminant, 0.0);
            const double square_root = std::sqrt(discriminant);
            const double q = -0.5 * (b + (b >= 0.0 ? square_root : -square_root));
            roots[root_count++] = q == 0.0 ? -b / (2.0 * a) : q / a;
            roots[root_count++] = q == 0.0 ? roots[0] : c / q;
        }

        for (int root_index = 0; root_index < root_count; ++root_index) {
            const double root = roots[root_index];
            if (!std::isfinite(root)) continue;
            const double correction =
                std::abs(root - tangent(component)) / (1.0 + std::abs(tangent(component)));
            if (correction >= best_correction) continue;
            Vec4 candidate = tangent;
            candidate(component) = root;
            best = NullProjection{candidate, component};
            best_correction = correction;
            if (component == 0) temporal_represented = true;
        }
    }
    return best;
}

namespace {

std::optional<GeodesicVariations> ProjectVariations(
    const Metric4d& metric, const Metric4d& inverse,
    const Tensor<Dual<double>, 4, 4, 4>& derivatives, const Vec4& unprojected,
    const NullProjection& projection, const PhaseVariations& phase) {
    const auto connection = TensorOps::Christoffel(metric, derivatives);
    const int component = projection.component;
    double denominator = 0.0;
    double denominator_scale = 0.0;
    for (int a = 0; a < 4; ++a) {
        const double term =
            (metric(component, a).real + metric(a, component).real) * projection.tangent(a);
        denominator += term;
        denominator_scale += std::abs(term);
    }
    if (!std::isfinite(denominator) || !std::isfinite(denominator_scale) ||
        std::abs(denominator) <= 256.0 * std::numeric_limits<double>::epsilon() * denominator_scale)
        return std::nullopt;
    GeodesicVariations result;
    for (std::size_t column = 0; column < result.size(); ++column) {
        Vec4 covector = phase[column].p;
        for (int a = 0; a < 4; ++a)
            for (int b = 0; b < 4; ++b)
                for (int axis = 0; axis < 4; ++axis)
                    covector(a) -=
                        derivatives(axis, a, b).real * phase[column].x(axis) * unprojected(b);
        Vec4 coordinate = TensorOps::RaiseIndex(covector, inverse);
        double numerator = 0.0;
        for (int a = 0; a < 4; ++a)
            for (int b = 0; b < 4; ++b) {
                for (int axis = 0; axis < 4; ++axis)
                    numerator += derivatives(axis, a, b).real * phase[column].x(axis) *
                                 projection.tangent(a) * projection.tangent(b);
                if (b != component)
                    numerator += (metric(a, b).real + metric(b, a).real) * projection.tangent(a) *
                                 coordinate(b);
            }
        coordinate(component) = -numerator / denominator;
        result[column].displacement = phase[column].x;
        result[column].derivative = coordinate;
        for (int mu = 0; mu < 4; ++mu)
            for (int a = 0; a < 4; ++a)
                for (int b = 0; b < 4; ++b)
                    result[column].derivative(mu) += connection.gamma(mu, a, b).real *
                                                     projection.tangent(a) * phase[column].x(b);
        if (!FiniteVector(result[column].displacement) || !FiniteVector(result[column].derivative))
            return std::nullopt;
    }
    return result;
}

}  // namespace

std::optional<Vec4> Geodesic::ProjectNullTangentPreservingBranch(const Vec4& tangent,
                                                                 const Metric4d& metric) {
    const auto projected = ProjectNullTangentWithBranch(tangent, metric);
    if (!projected) return std::nullopt;
    return projected->tangent;
}

// Compare the actual represented outputs, including projection/event
// amplification. Every component is checked before max/RMS aggregation so
// NaNs cannot disappear through std::max ordering.
double Geodesic::CoupledStateError(const Lightray& first,
                                   const GeodesicVariations& first_variations,
                                   const Lightray& second,
                                   const GeodesicVariations& second_variations,
                                   const IntegratorConfig& config,
                                   const Rk45CoupledState& control) {
    const double invalid = std::numeric_limits<double>::infinity();
    if (!ValidCoupledControl(control)) return invalid;
    double central_sum = 0.0;
    for (int component = 0; component < 4; ++component) {
        const double position_scale =
            config.abs_tolerance * control.length_scale +
            config.rel_tolerance *
                std::max(std::abs(first.position(component)), std::abs(second.position(component)));
        const double tangent_scale =
            config.abs_tolerance * control.frequency_scale +
            config.rel_tolerance *
                std::max(std::abs(first.velocity(component)), std::abs(second.velocity(component)));
        const double x_error =
            std::abs(first.position(component) - second.position(component)) / position_scale;
        const double k_error =
            std::abs(first.velocity(component) - second.velocity(component)) / tangent_scale;
        if (!std::isfinite(x_error) || !std::isfinite(k_error) || !std::isfinite(position_scale) ||
            !std::isfinite(tangent_scale) || !(position_scale > 0.0) || !(tangent_scale > 0.0))
            return invalid;
        central_sum += x_error * x_error + k_error * k_error;
    }
    double ratio = std::sqrt(central_sum / 8.0);
    for (std::size_t column = 0; column < first_variations.size(); ++column)
        for (int component = 0; component < 4; ++component) {
            const auto& first_column = first_variations[column];
            const auto& second_column = second_variations[column];
            const double x_scale =
                control.tolerance *
                (control.length_scale + std::max(std::abs(first_column.displacement(component)),
                                                 std::abs(second_column.displacement(component))));
            const double v_scale =
                control.tolerance *
                (control.frequency_scale + std::max(std::abs(first_column.derivative(component)),
                                                    std::abs(second_column.derivative(component))));
            const double x_error = std::abs(first_column.displacement(component) -
                                            second_column.displacement(component)) /
                                   x_scale;
            const double v_error =
                std::abs(first_column.derivative(component) - second_column.derivative(component)) /
                v_scale;
            if (!std::isfinite(x_scale) || !std::isfinite(v_scale) || !(x_scale > 0.0) ||
                !(v_scale > 0.0) || !std::isfinite(x_error) || !std::isfinite(v_error))
                return invalid;
            ratio = std::max(ratio, std::max(x_error, v_error));
        }
    return std::isfinite(ratio) ? ratio : invalid;
}

std::optional<CoupledSegmentSample> Geodesic::SampleCoupledSegment(
    IMetric* metric, const Lightray& start, const GeodesicVariations& start_variations,
    const Lightray& end, const GeodesicVariations& end_variations, double interval, double fraction,
    const Vec4* event_normal, std::uint64_t* metric_evaluations,
    const CoupledSegmentIncrement* increment) {
    if (!metric || !std::isfinite(interval) || !(interval > 0.0) || !std::isfinite(fraction) ||
        fraction < 0.0 || fraction > 1.0 || !metric->IsValidEvent(start.position) ||
        !metric->IsValidEvent(end.position))
        return std::nullopt;
    const double s = fraction;
    CoupledSegmentSample result;
    result.ray = end;
    const Vec4 displacement = increment ? increment->position : end.position - start.position;
    const auto central =
        SampleAcceptedTraceSegment(start.position, start.velocity, end.position, end.velocity,
                                   interval, s, &displacement, &result.ray.acceleration);
    result.ray.position = central.position;
    result.ray.velocity = central.tangent;
    if (!FiniteVector(result.ray.position) || !FiniteVector(result.ray.velocity) ||
        !FiniteVector(result.ray.acceleration) || !metric->IsValidEvent(result.ray.position))
        return std::nullopt;
    Metric4d g;
    Tensor<Dual<double>, 4, 4, 4> dg;
    metric->Evaluate(start.position, g, dg);
    if (metric_evaluations) ++*metric_evaluations;
    const auto initial_connection = TensorOps::Christoffel(g, dg);
    metric->Evaluate(end.position, g, dg);
    if (metric_evaluations) ++*metric_evaluations;
    const auto final_connection = TensorOps::Christoffel(g, dg);
    metric->Evaluate(result.ray.position, g, dg);
    if (metric_evaluations) ++*metric_evaluations;
    const auto event_connection = TensorOps::Christoffel(g, dg);
    double denominator = 0.0, denominator_scale = 0.0;
    if (event_normal) {
        if (!FiniteVector(*event_normal)) return std::nullopt;
        for (int component = 0; component < 4; ++component) {
            const double term = (*event_normal)(component)*result.ray.velocity(component);
            denominator += term;
            denominator_scale += std::abs(term);
        }
        if (!std::isfinite(denominator) || !std::isfinite(denominator_scale) ||
            std::abs(denominator) <=
                256.0 * std::numeric_limits<double>::epsilon() * denominator_scale)
            return std::nullopt;
    }
    // A moving event differentiates the geodesic flow: X=xi+k*eta and
    // coordinate K=dki+a*eta. Use its ODE acceleration here; roundoff in
    // the cubic location interpolant's second derivative is amplified by
    // eta and would introduce a fictitious covariant acceleration.
    Vec4 event_acceleration;
    if (event_normal) {
        for (int mu = 0; mu < 4; ++mu)
            for (int a = 0; a < 4; ++a)
                for (int b = 0; b < 4; ++b)
                    event_acceleration(mu) -= event_connection.gamma(mu, a, b).real *
                                              result.ray.velocity(a) * result.ray.velocity(b);
        if (!FiniteVector(event_acceleration)) return std::nullopt;
    }
    for (std::size_t column = 0; column < result.variations.size(); ++column) {
        const Vec4 initial_coordinate =
            CoordinateVariation(start_variations[column], start.velocity, initial_connection);
        const Vec4 final_coordinate =
            CoordinateVariation(end_variations[column], end.velocity, final_connection);
        const Vec4 variation_increment =
            increment ? increment->displacement[column]
                      : end_variations[column].displacement - start_variations[column].displacement;
        const auto variation_polynomial =
            detail::MakeIncrementHermite(start_variations[column].displacement, initial_coordinate,
                                         final_coordinate, interval, variation_increment);
        const auto variation_sample = variation_polynomial.Sample(s);
        Vec4 displacement = variation_sample.position;
        Vec4 coordinate = variation_sample.tangent;
        if (event_normal) {
            double numerator = 0.0;
            for (int component = 0; component < 4; ++component)
                numerator += (*event_normal)(component)*displacement(component);
            const double shift = -numerator / denominator;
            displacement += result.ray.velocity * shift;
            coordinate += event_acceleration * shift;
        }
        result.variations[column] = {displacement, coordinate};
        for (int mu = 0; mu < 4; ++mu)
            for (int a = 0; a < 4; ++a)
                for (int b = 0; b < 4; ++b)
                    result.variations[column].derivative(mu) +=
                        event_connection.gamma(mu, a, b).real * result.ray.velocity(a) *
                        displacement(b);
        if (!FiniteVector(displacement) || !FiniteVector(result.variations[column].derivative))
            return std::nullopt;
    }
    return result;
}

static bool IntegrateStepRk45Candidate(Lightray& ray, IMetric* metric,
                                       const IntegratorConfig& config, Rk45CoupledState* coupled,
                                       Rk45CoupledComparison* comparison) {
    using namespace dp45;

    SIRIUS_PRE(metric != nullptr);
    SIRIUS_PRE(IsRepresentedIntegratorStepControl(config));
    SIRIUS_PRE(std::isfinite(ray.step_size) && ray.step_size >= config.min_step &&
               ray.step_size <= config.max_step);

    Vec4 x0 = ray.position;
    float h = ray.step_size;

    if (!metric->IsValidEvent(x0)) {
        ray.terminated = 3;
        return false;
    }

    // Metric and initial momentum.
    Metric4d g0;
    Tensor<Dual<double>, 4, 4, 4> dg0;
    metric->Evaluate(x0, g0, dg0);
    Vec4 p0 = ComputeMomentum(ray.velocity, g0);
    Vec4 k0 = ray.velocity;
    constexpr double kMaximumRelativeNullResidual = 1e-6;
    if (RelativeNullResidual(k0, g0) > kMaximumRelativeNullResidual) {
        ray.terminated = 6;
        return false;
    }

    if (coupled) ++coupled->central_stages;

    // Stage 1.
    Vec4 k1_x = k0, k1_p = MomentumDerivative(p0, k0, dg0);

    // Stages 2-6.
    Vec4 k2_x, k2_p, k3_x, k3_p, k4_x, k4_p, k5_x, k5_p, k6_x, k6_p;

    Vec4 x2 = x0 + k1_x * (a21 * h);
    Vec4 p2 = p0 + k1_p * (a21 * h);
    if (coupled) ++coupled->central_stages;
    if (!EvaluateRk45Stage(x2, p2, metric, k2_x, k2_p)) {
        return RejectUnrepresentedStage(ray, config.min_step);
    }

    Vec4 x3 = x0 + k1_x * (a31 * h) + k2_x * (a32 * h);
    Vec4 p3 = p0 + k1_p * (a31 * h) + k2_p * (a32 * h);
    if (coupled) ++coupled->central_stages;
    if (!EvaluateRk45Stage(x3, p3, metric, k3_x, k3_p)) {
        return RejectUnrepresentedStage(ray, config.min_step);
    }

    Vec4 x4 = x0 + k1_x * (a41 * h) + k2_x * (a42 * h) + k3_x * (a43 * h);
    Vec4 p4 = p0 + k1_p * (a41 * h) + k2_p * (a42 * h) + k3_p * (a43 * h);
    if (coupled) ++coupled->central_stages;
    if (!EvaluateRk45Stage(x4, p4, metric, k4_x, k4_p)) {
        return RejectUnrepresentedStage(ray, config.min_step);
    }

    Vec4 x5 = x0 + k1_x * (a51 * h) + k2_x * (a52 * h) + k3_x * (a53 * h) + k4_x * (a54 * h);
    Vec4 p5 = p0 + k1_p * (a51 * h) + k2_p * (a52 * h) + k3_p * (a53 * h) + k4_p * (a54 * h);
    if (coupled) ++coupled->central_stages;
    if (!EvaluateRk45Stage(x5, p5, metric, k5_x, k5_p)) {
        return RejectUnrepresentedStage(ray, config.min_step);
    }

    Vec4 x6 = x0 + k1_x * (a61 * h) + k2_x * (a62 * h) + k3_x * (a63 * h) + k4_x * (a64 * h) +
              k5_x * (a65 * h);
    Vec4 p6 = p0 + k1_p * (a61 * h) + k2_p * (a62 * h) + k3_p * (a63 * h) + k4_p * (a64 * h) +
              k5_p * (a65 * h);
    if (coupled) ++coupled->central_stages;
    if (!EvaluateRk45Stage(x6, p6, metric, k6_x, k6_p)) {
        return RejectUnrepresentedStage(ray, config.min_step);
    }

    // 5th order solution.
    const Vec4 position_increment = (k1_x * b1 + k3_x * b3 + k4_x * b4 + k5_x * b5 + k6_x * b6) * h;
    Vec4 new_position = x0 + position_increment;
    Vec4 new_momentum = p0 + (k1_p * b1 + k3_p * b3 + k4_p * b4 + k5_p * b5 + k6_p * b6) * h;

    // Stage 7 (FSAL).
    Vec4 k7_x, k7_p;
    if (coupled) ++coupled->central_stages;
    if (!EvaluateRk45Stage(new_position, new_momentum, metric, k7_x, k7_p)) {
        return RejectUnrepresentedStage(ray, config.min_step);
    }
    Vec4 new_velocity = k7_x;

    PhaseVariations fifth_variations, fourth_variations;
    CoupledSegmentIncrement fifth_increment, fourth_increment;
    fifth_increment.position = position_increment;
    if (coupled &&
        !VariationPair(*metric, {x0, x2, x3, x4, x5, x6, new_position},
                       {p0, p2, p3, p4, p5, p6, new_momentum}, k0, h, *coupled, fifth_variations,
                       fourth_variations, fifth_increment, fourth_increment)) {
        coupled->failure = CoupledStepFailure::DerivativeDomain;
        return RejectUnrepresentedStage(ray, config.min_step);
    }

    // Error estimation.
    Vec4 error_x = (k1_x * e1 + k3_x * e3 + k4_x * e4 + k5_x * e5 + k6_x * e6 + k7_x * e7) * h;
    Vec4 error_p = (k1_p * e1 + k3_p * e3 + k4_p * e4 + k5_p * e5 + k6_p * e6 + k7_p * e7) * h;
    float error_norm =
        ComputeRk45ErrorNorm(error_x, error_p, x0, p0, new_position, new_momentum, config);

    // Step acceptance.
    if (error_norm > 1.0f) {
        // Reaching the configured floor is not itself a failure: the caller
        // must be allowed to retry once at that floor.  Only a candidate that
        // was already evaluated at min_step may terminate.  The former order
        // classified any rejection that *selected* min_step as terminal and
        // created a large false shadow around near-critical Kerr rays.
        if (h <= config.min_step) {
            ray.terminated = 5;
        } else {
            ray.step_size = Geodesic::ComputeOptimalStep(h, error_norm, 1.0f, config);
        }
        return false;
    }

    if (HasInvalidState(new_position, new_velocity)) {
        ray.terminated = 3;
        return false;
    }

    // The Hamiltonian flow must remain on the null constraint surface. Bound
    // the unprojected defect first, then remove its accumulated roundoff with
    // the nearest represented temporal root. When no temporal root exists, the
    // smallest normalized spatial correction leaves k^0 unchanged. An
    // unrepresented projection still fails closed.
    {
        Metric4d g_check;
        Tensor<Dual<double>, 4, 4, 4> dg_check;
        metric->Evaluate(new_position, g_check, dg_check);
        const double unprojected_residual = RelativeNullResidual(new_velocity, g_check);
        const auto projected = ProjectNullTangentWithBranch(new_velocity, g_check);
        const bool represented_projection =
            projected.has_value() && RelativeNullResidual(projected->tangent, g_check) <=
                                         256.0 * std::numeric_limits<double>::epsilon();
        if (unprojected_residual > kMaximumRelativeNullResidual || !represented_projection) {
            if (h <= config.min_step) {
                ray.terminated = 6;
            } else {
                ray.step_size = std::max(config.min_step, h * 0.5f);
            }
            return false;
        }
        if (coupled) {
            const auto inverse = InverseAt(metric, new_position, g_check);
            const auto fifth = ProjectVariations(g_check, inverse, dg_check, new_velocity,
                                                 *projected, fifth_variations);
            Lightray lower = ray;
            lower.position = new_position - error_x;
            const Vec4 lower_momentum = new_momentum - error_p;
            if (!metric->IsValidEvent(lower.position)) {
                coupled->failure = CoupledStepFailure::DerivativeDomain;
                return RejectUnrepresentedStage(ray, config.min_step);
            }
            Metric4d lower_metric;
            Tensor<Dual<double>, 4, 4, 4> lower_derivatives;
            metric->Evaluate(lower.position, lower_metric, lower_derivatives);
            ++coupled->variation_metric_evaluations;
            const auto lower_inverse = InverseAt(metric, lower.position, lower_metric);
            const Vec4 lower_unprojected = TensorOps::RaiseIndex(lower_momentum, lower_inverse);
            const auto lower_projection =
                ProjectNullTangentWithBranch(lower_unprojected, lower_metric);
            const auto fourth =
                lower_projection
                    ? ProjectVariations(lower_metric, lower_inverse, lower_derivatives,
                                        lower_unprojected, *lower_projection, fourth_variations)
                    : std::nullopt;
            if (!fifth || !fourth || projected->component != lower_projection->component) {
                coupled->failure = CoupledStepFailure::Projection;
                if (h <= config.min_step)
                    ray.terminated = 6;
                else
                    ray.step_size = std::max(config.min_step, h * 0.5f);
                return false;
            }
            lower.velocity = lower_projection->tangent;
            Lightray upper = ray;
            upper.position = new_position;
            upper.velocity = projected->tangent;
            const double projected_error =
                Geodesic::CoupledStateError(upper, *fifth, lower, *fourth, config, *coupled);
            error_norm = std::max(error_norm, static_cast<float>(projected_error));
            if (!std::isfinite(projected_error) || error_norm > 1.0f) {
                coupled->failure = CoupledStepFailure::Projection;
                if (h <= config.min_step)
                    ray.terminated = 5;
                else
                    ray.step_size = std::max(config.min_step, h * 0.5f);
                return false;
            }
            coupled->variations = *fifth;
            fourth_increment.position = position_increment - error_x;
            comparison->full_increment = fifth_increment;
            comparison->lower_increment = fourth_increment;
            comparison->lower_order = lower;
            comparison->lower_variations = *fourth;
            comparison->error_ratio = projected_error;
        }
        new_velocity = projected->tangent;
    }

    // Commit only after both embedded error and the physical constraint admit
    // the candidate state.
    ray.position = new_position;
    ray.velocity = new_velocity;
    ray.acceleration = Geodesic::CalculateAcceleration(new_velocity, new_position, metric);
    ray.proper_time += h;
    ray.coordinate_time += static_cast<float>(h * std::abs(new_velocity(0)));
    ray.step_size = Geodesic::ComputeOptimalStep(h, error_norm, 1.0f, config);

    return true;
}

bool Geodesic::IntegrateStepRk45(Lightray& ray, IMetric* metric, const IntegratorConfig& config,
                                 Rk45CoupledState* coupled, Rk45CoupledComparison* comparison) {
    SIRIUS_PRE((coupled == nullptr) == (comparison == nullptr));
    if (!coupled) return IntegrateStepRk45Candidate(ray, metric, config, nullptr, nullptr);
    coupled->failure = CoupledStepFailure::None;
    if (!ValidCoupledControl(*coupled)) {
        coupled->failure = CoupledStepFailure::InvalidState;
        ray.terminated = 3;
        return false;
    }
    const Lightray previous = ray;
    const Rk45CoupledState previous_coupled = *coupled;
    if (!IntegrateStepRk45Candidate(ray, metric, config, coupled, comparison)) return false;

    // An embedded endpoint pair cannot expose error shared by its Hermite
    // interpolants. Compare against a separately integrated midpoint before
    // allowing the tracer to consume any part of this private interval.
    Lightray midpoint = previous;
    midpoint.step_size = previous.step_size * 0.5f;
    Rk45CoupledState midpoint_coupled = previous_coupled;
    midpoint_coupled.central_stages = coupled->central_stages;
    midpoint_coupled.variation_stages = coupled->variation_stages;
    midpoint_coupled.variation_metric_evaluations = coupled->variation_metric_evaluations;
    Rk45CoupledComparison midpoint_comparison;
    // The minimum bounds accepted physical intervals, not private diagnostic
    // stages. As with the fractional DP stages, these half-steps never commit
    // affine distance or source effects. An unrepresentable half still fails.
    auto diagnostic_config = config;
    diagnostic_config.min_step = config.min_step * 0.5f;
    const bool represented = IsRepresentedIntegratorStepControl(diagnostic_config) &&
                             midpoint.step_size >= diagnostic_config.min_step &&
                             IntegrateStepRk45Candidate(midpoint, metric, diagnostic_config,
                                                        &midpoint_coupled, &midpoint_comparison);
    coupled->central_stages = midpoint_coupled.central_stages;
    coupled->variation_stages = midpoint_coupled.variation_stages;
    coupled->variation_metric_evaluations = midpoint_coupled.variation_metric_evaluations;
    const auto interpolated =
        represented ? SampleCoupledSegment(metric, previous, previous_coupled.variations, ray,
                                           coupled->variations, previous.step_size, 0.5, nullptr,
                                           &coupled->variation_metric_evaluations,
                                           &comparison->full_increment)
                    : std::nullopt;
    const double interior_error =
        interpolated ? CoupledStateError(midpoint, midpoint_coupled.variations, interpolated->ray,
                                         interpolated->variations, config, *coupled)
                     : std::numeric_limits<double>::infinity();
    if (!std::isfinite(interior_error) || interior_error > 1.0) {
        coupled->variations = previous_coupled.variations;
        coupled->failure = CoupledStepFailure::Interpolation;
        ray = previous;
        if (previous.step_size <= config.min_step)
            ray.terminated = 5;
        else
            ray.step_size = std::max(config.min_step, previous.step_size * 0.5f);
        return false;
    }
    // Complete the independent subdivided trajectory. Its two dense pieces
    // let the tracer test the same terminal event without sharing the full
    // interval's Hermite error, including amplification by the event root.
    Lightray refined = midpoint;
    refined.step_size = previous.step_size * 0.5f;
    auto refined_coupled = midpoint_coupled;
    refined_coupled.variation_metric_evaluations = coupled->variation_metric_evaluations;
    Rk45CoupledComparison refined_comparison;
    const bool refined_accepted = IntegrateStepRk45Candidate(refined, metric, diagnostic_config,
                                                             &refined_coupled, &refined_comparison);
    coupled->central_stages = refined_coupled.central_stages;
    coupled->variation_stages = refined_coupled.variation_stages;
    coupled->variation_metric_evaluations = refined_coupled.variation_metric_evaluations;
    const double refined_error =
        refined_accepted ? CoupledStateError(ray, coupled->variations, refined,
                                             refined_coupled.variations, config, *coupled)
                         : std::numeric_limits<double>::infinity();
    if (!std::isfinite(refined_error) || refined_error > 1.0) {
        coupled->variations = previous_coupled.variations;
        coupled->failure = CoupledStepFailure::Interpolation;
        ray = previous;
        if (previous.step_size <= config.min_step)
            ray.terminated = 5;
        else
            ray.step_size = std::max(config.min_step, previous.step_size * 0.5f);
        return false;
    }
    comparison->midpoint_increment = midpoint_comparison.full_increment;
    comparison->refined_increment = refined_comparison.full_increment;
    comparison->midpoint = midpoint;
    comparison->midpoint_variations = midpoint_coupled.variations;
    comparison->refined_endpoint = refined;
    comparison->refined_variations = refined_coupled.variations;
    comparison->error_ratio = std::max({comparison->error_ratio, interior_error, refined_error});
    return true;
}

}  // namespace sirius::core
