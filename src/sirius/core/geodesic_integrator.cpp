// Implementation of the geodesic integrator declared in geodesic_integrator.h.
// Ported from PHGD001A.cpp; arithmetic order, casts, and literals bit-identical.
//
// Hamiltonian formulation with covariant momenta preserves the null constraint
// H = (1/2) g^mu_nu p_mu p_nu = 0 automatically. Hamilton's equations:
//   dx^mu/dlambda = g^mu_nu p_nu,   dp_mu/dlambda = (1/2)(d g_alpha_beta/dx^mu) k^alpha k^beta.

#include "sirius/core/geodesic_integrator.h"

#include "sirius/base/contracts.h"
#include "sirius/core/constants.h"

#include <algorithm>
#include <cmath>

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
static void EvaluateRk45Stage(const Vec4& x, const Vec4& p, IMetric* metric, Vec4& k_x, Vec4& k_p);
static Vec4 ComputeMomentum(const Vec4& velocity, const Metric4d& g);

bool Geodesic::IntegrateStep(Lightray& ray, IMetric* metric, float min_step, float max_step) {
    SIRIUS_PRE(metric != nullptr);

    Vec4 x0 = ray.position;
    float h = ray.step_size;

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
    EvaluateRk45Stage(x1, p1, metric, k2_x, k2_p);

    Vec4 x2 = x0 + k2_x * (0.5f * h), p2 = p0 + k2_p * (0.5f * h);
    Vec4 k3_x, k3_p;
    EvaluateRk45Stage(x2, p2, metric, k3_x, k3_p);

    Vec4 x3 = x0 + k3_x * h, p3 = p0 + k3_p * h;
    Vec4 k4_x, k4_p;
    EvaluateRk45Stage(x3, p3, metric, k4_x, k4_p);

    // Combine RK4.
    Vec4 new_position = x0 + (k1_x + k2_x * 2.0f + k3_x * 2.0f + k4_x) * (h / 6.0f);
    Vec4 new_momentum = p0 + (k1_p + k2_p * 2.0f + k3_p * 2.0f + k4_p) * (h / 6.0f);

    // New velocity.
    Metric4d g_new;
    Tensor<Dual<double>, 4, 4, 4> dg_new;
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
    ray.running_dlambda_dnew *= (1.0f + velocity_change * 0.1f);
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
    Metric4d g;
    Tensor<Dual<double>, 4, 4, 4> dg;
    metric->Evaluate(ray.position, g, dg);

    Vec4 observer_lower = TensorOps::LowerIndex(observer.velocity, g);

    double dot_product = 0.0;
    for (int mu = 0; mu < 4; mu++) {
        dot_product += ray.velocity(mu) * observer_lower(mu);
    }

    return static_cast<float>(dot_product / ray.ku_uobsu);
}

ObserverState Geodesic::CreateObserver(const Vec4& position, const Vec4& velocity,
                                       IMetric* metric) {
    ObserverState observer;
    observer.position = position;
    observer.velocity = velocity;

    Metric4d g;
    Tensor<Dual<double>, 4, 4, 4> dg;
    metric->Evaluate(position, g, dg);
    double velocity_norm = TensorOps::InnerProduct(velocity, velocity, g);
    observer.is_timelike = velocity_norm < 0.0;

    if (!observer.is_timelike) {
        observer.velocity = Vec4();  // Default constructor initialises to zeros.
        observer.velocity(0) = 1.0;  // Set the time component.
        velocity_norm = TensorOps::InnerProduct(observer.velocity, observer.velocity, g);
    }

    double normalization = 1.0 / std::sqrt(-velocity_norm);
    observer.velocity *= normalization;

    CalculateTetrads(observer, metric);

    return observer;
}

void Geodesic::CalculateTetrads(ObserverState& observer, IMetric* metric) {
    Metric4d g;
    Tensor<Dual<double>, 4, 4, 4> dg;
    metric->Evaluate(observer.position, g, dg);

    observer.e0 = observer.velocity;

    observer.e1 = Vec4();
    observer.e1(1) = 1.0f;
    Vec4 e1_parallel = TensorOps::LowerIndex(observer.e1, g);
    double e1_dot_e0 = 0.0;
    for (int mu = 0; mu < 4; mu++) {
        e1_dot_e0 += observer.e0(mu) * e1_parallel(mu);
    }
    observer.e1 = observer.e1 - observer.e0 * e1_dot_e0;

    observer.e2 = Vec4();
    observer.e2(2) = 1.0f;
    Vec4 e2_parallel = TensorOps::LowerIndex(observer.e2, g);
    double e2_dot_e0 = 0.0;
    double e2_dot_e1 = 0.0;
    for (int mu = 0; mu < 4; mu++) {
        e2_dot_e0 += observer.e0(mu) * e2_parallel(mu);
        e2_dot_e1 += observer.e1(mu) * e2_parallel(mu);
    }
    observer.e2 = observer.e2 - observer.e0 * e2_dot_e0 - observer.e1 * e2_dot_e1;

    observer.e3 = Vec4();
    observer.e3(3) = 1.0f;
    Vec4 e3_parallel = TensorOps::LowerIndex(observer.e3, g);
    double e3_dot_e0 = 0.0;
    double e3_dot_e1 = 0.0;
    double e3_dot_e2 = 0.0;
    for (int mu = 0; mu < 4; mu++) {
        e3_dot_e0 += observer.e0(mu) * e3_parallel(mu);
        e3_dot_e1 += observer.e1(mu) * e3_parallel(mu);
        e3_dot_e2 += observer.e2(mu) * e3_parallel(mu);
    }
    observer.e3 =
        observer.e3 - observer.e0 * e3_dot_e0 - observer.e1 * e3_dot_e1 - observer.e2 * e3_dot_e2;

    double e1_norm = std::sqrt(std::abs(TensorOps::InnerProduct(observer.e1, observer.e1, g)));
    double e2_norm = std::sqrt(std::abs(TensorOps::InnerProduct(observer.e2, observer.e2, g)));
    double e3_norm = std::sqrt(std::abs(TensorOps::InnerProduct(observer.e3, observer.e3, g)));

    if (e1_norm > 1e-10) observer.e1 /= e1_norm;
    if (e2_norm > 1e-10) observer.e2 /= e2_norm;
    if (e3_norm > 1e-10) observer.e3 /= e3_norm;
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
    config.use_rk45 = true;
    return config;
}

float Geodesic::ComputeOptimalStep(float h, float error, float tolerance,
                                   const IntegratorConfig& config) {
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
static void EvaluateRk45Stage(const Vec4& x, const Vec4& p, IMetric* metric, Vec4& k_x, Vec4& k_p) {
    Metric4d g;
    Tensor<Dual<double>, 4, 4, 4> dg;
    metric->Evaluate(x, g, dg);
    k_x = TensorOps::RaiseIndex(p, InverseAt(metric, x, g));
    k_p = MomentumDerivative(p, k_x, dg);
}

// RK45 error norm.
static float ComputeRk45ErrorNorm(const Vec4& error_x, const Vec4& new_position,
                                  const IntegratorConfig& config) {
    float error_norm = 0.0f;
    for (int i = 0; i < 4; i++) {
        const float scale = static_cast<float>(config.abs_tolerance +
                                               config.rel_tolerance * std::abs(new_position(i)));
        const float err_i = static_cast<float>(std::abs(error_x(i)) / scale);
        error_norm += err_i * err_i;
    }
    return std::sqrt(error_norm / 4.0f);
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

bool Geodesic::IntegrateStepRk45(Lightray& ray, IMetric* metric, const IntegratorConfig& config) {
    using namespace dp45;

    Vec4 x0 = ray.position;
    float h = ray.step_size;

    // Metric and initial momentum.
    Metric4d g0;
    Tensor<Dual<double>, 4, 4, 4> dg0;
    metric->Evaluate(x0, g0, dg0);
    Vec4 p0 = ComputeMomentum(ray.velocity, g0);
    Vec4 k0 = ray.velocity;

    // Stage 1.
    Vec4 k1_x = k0, k1_p = MomentumDerivative(p0, k0, dg0);

    // Stages 2-6.
    Vec4 k2_x, k2_p, k3_x, k3_p, k4_x, k4_p, k5_x, k5_p, k6_x, k6_p;

    Vec4 x2 = x0 + k1_x * (a21 * h);
    Vec4 p2 = p0 + k1_p * (a21 * h);
    EvaluateRk45Stage(x2, p2, metric, k2_x, k2_p);

    Vec4 x3 = x0 + k1_x * static_cast<float>(a31 * h) + k2_x * static_cast<float>(a32 * h);
    Vec4 p3 = p0 + k1_p * static_cast<float>(a31 * h) + k2_p * static_cast<float>(a32 * h);
    EvaluateRk45Stage(x3, p3, metric, k3_x, k3_p);

    Vec4 x4 = x0 + k1_x * static_cast<float>(a41 * h) + k2_x * static_cast<float>(a42 * h) +
              k3_x * static_cast<float>(a43 * h);
    Vec4 p4 = p0 + k1_p * static_cast<float>(a41 * h) + k2_p * static_cast<float>(a42 * h) +
              k3_p * static_cast<float>(a43 * h);
    EvaluateRk45Stage(x4, p4, metric, k4_x, k4_p);

    Vec4 x5 = x0 + k1_x * static_cast<float>(a51 * h) + k2_x * static_cast<float>(a52 * h) +
              k3_x * static_cast<float>(a53 * h) + k4_x * static_cast<float>(a54 * h);
    Vec4 p5 = p0 + k1_p * static_cast<float>(a51 * h) + k2_p * static_cast<float>(a52 * h) +
              k3_p * static_cast<float>(a53 * h) + k4_p * static_cast<float>(a54 * h);
    EvaluateRk45Stage(x5, p5, metric, k5_x, k5_p);

    Vec4 x6 = x0 + k1_x * static_cast<float>(a61 * h) + k2_x * static_cast<float>(a62 * h) +
              k3_x * static_cast<float>(a63 * h) + k4_x * static_cast<float>(a64 * h) +
              k5_x * static_cast<float>(a65 * h);
    Vec4 p6 = p0 + k1_p * static_cast<float>(a61 * h) + k2_p * static_cast<float>(a62 * h) +
              k3_p * static_cast<float>(a63 * h) + k4_p * static_cast<float>(a64 * h) +
              k5_p * static_cast<float>(a65 * h);
    EvaluateRk45Stage(x6, p6, metric, k6_x, k6_p);

    // 5th order solution.
    Vec4 new_position = x0 + (k1_x * b1 + k3_x * b3 + k4_x * b4 + k5_x * b5 + k6_x * b6) * h;
    Vec4 new_momentum = p0 + (k1_p * b1 + k3_p * b3 + k4_p * b4 + k5_p * b5 + k6_p * b6) * h;

    // Stage 7 (FSAL).
    Vec4 k7_x, k7_p;
    EvaluateRk45Stage(new_position, new_momentum, metric, k7_x, k7_p);
    Vec4 new_velocity = k7_x;

    // Error estimation.
    Vec4 error_x = (k1_x * e1 + k3_x * e3 + k4_x * e4 + k5_x * e5 + k6_x * e6 + k7_x * e7) * h;
    float error_norm = ComputeRk45ErrorNorm(error_x, new_position, config);

    // Step acceptance.
    if (error_norm > 1.0f) {
        ray.step_size = ComputeOptimalStep(h, error_norm, 1.0f, config);
        if (ray.step_size <= config.min_step) {
            ray.terminated = 5;
            return false;
        }
        return false;
    }

    // Update ray state.
    ray.position = new_position;
    ray.velocity = new_velocity;
    ray.acceleration = CalculateAcceleration(new_velocity, new_position, metric);
    ray.proper_time += h;
    ray.coordinate_time += static_cast<float>(h * std::abs(new_velocity(0)));
    ray.step_size = ComputeOptimalStep(h, error_norm, 1.0f, config);

    if (HasInvalidState(new_position, new_velocity)) {
        ray.terminated = 3;
        return false;
    }

    // Periodic null re-normalisation: g_mu_nu k^mu k^nu = 0 drifts during
    // integration from numerical error, so the velocity is re-normalised when the
    // drift exceeds the threshold. Checked every step for tight control.
    {
        Metric4d g_check;
        Tensor<Dual<double>, 4, 4, 4> dg_check;
        metric->Evaluate(new_position, g_check, dg_check);

        // Null condition violation H = g_mu_nu k^mu k^nu (should be 0).
        double null_violation = TensorOps::InnerProduct(new_velocity, new_velocity, g_check);

        // Phase 3 P3: tightened from 1e-4 to 1e-6 for better geodesic accuracy.
        constexpr double kNullRenormThreshold = 1e-6;
        if (std::abs(null_violation) > kNullRenormThreshold) {
            Vec4 renormalized = TensorOps::NormalizeNull(new_velocity, g_check);

            // Reject a numerically failed normalisation.
            bool valid = true;
            for (int i = 0; i < 4; i++) {
                if (std::isnan(renormalized(i)) || std::isinf(renormalized(i))) {
                    valid = false;
                    break;
                }
            }

            if (valid) {
                ray.velocity = renormalized;
            }
        }
    }

    return true;
}

}  // namespace sirius::core
