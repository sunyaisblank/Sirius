// PHGD001A.cpp - Geodesic Integration Implementation
//
// Hamiltonian formulation with covariant momenta for constraint-preserving
// null geodesic integration. H = (1/2) g^μν p_μ p_ν = 0 is preserved automatically.
//
// Hamilton's equations: dx^μ/dλ = g^μν p_ν, dp_μ/dλ = (1/2)(∂g_αβ/∂x^μ) k^α k^β

#include "PHGD001A.h"
#include <algorithm>
#include <cmath>

// =============================================================================
// Dormand-Prince RK45 Coefficients
// =============================================================================
namespace DP45 {
    // Node coefficients (c_i)
    constexpr double c2 = 1.0 / 5.0;
    constexpr double c3 = 3.0 / 10.0;
    constexpr double c4 = 4.0 / 5.0;
    constexpr double c5 = 8.0 / 9.0;
    constexpr double c6 = 1.0;
    constexpr double c7 = 1.0;
    
    // Runge-Kutta matrix (a_ij) - row by row
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
    
    constexpr double a71 = 35.0 / 384.0;
    // a72 = 0
    constexpr double a73 = 500.0 / 1113.0;
    constexpr double a74 = 125.0 / 192.0;
    constexpr double a75 = -2187.0 / 6784.0;
    constexpr double a76 = 11.0 / 84.0;
    
    // 5th order weights (b_i) - for the solution
    constexpr double b1 = 35.0 / 384.0;
    // b2 = 0
    constexpr double b3 = 500.0 / 1113.0;
    constexpr double b4 = 125.0 / 192.0;
    constexpr double b5 = -2187.0 / 6784.0;
    constexpr double b6 = 11.0 / 84.0;
    // b7 = 0
    
    // 4th order weights (b*_i) - for the error estimate
    constexpr double bs1 = 5179.0 / 57600.0;
    // bs2 = 0
    constexpr double bs3 = 7571.0 / 16695.0;
    constexpr double bs4 = 393.0 / 640.0;
    constexpr double bs5 = -92097.0 / 339200.0;
    constexpr double bs6 = 187.0 / 2100.0;
    constexpr double bs7 = 1.0 / 40.0;
    
    // Error coefficients (e_i = b_i - b*_i) for direct error computation
    constexpr double e1 = b1 - bs1;    // 71/57600
    // e2 = 0
    constexpr double e3 = b3 - bs3;    // -71/16695
    constexpr double e4 = b4 - bs4;    // 71/1920
    constexpr double e5 = b5 - bs5;    // -17253/339200
    constexpr double e6 = b6 - bs6;    // 22/525
    constexpr double e7 = -bs7;        // -1/40
}

// Helper: Compute inverse of 4x4 metric using Cramer's rule
static void invertMetric4x4(const double m[4][4], double g_inv[4][4]) {
    // Compute 2x2 minors for cofactor expansion
    double A2323 = m[2][2] * m[3][3] - m[2][3] * m[3][2];
    double A1323 = m[2][1] * m[3][3] - m[2][3] * m[3][1];
    double A1223 = m[2][1] * m[3][2] - m[2][2] * m[3][1];
    double A0323 = m[2][0] * m[3][3] - m[2][3] * m[3][0];
    double A0223 = m[2][0] * m[3][2] - m[2][2] * m[3][0];
    double A0123 = m[2][0] * m[3][1] - m[2][1] * m[3][0];
    double A2313 = m[1][2] * m[3][3] - m[1][3] * m[3][2];
    double A1313 = m[1][1] * m[3][3] - m[1][3] * m[3][1];
    double A1213 = m[1][1] * m[3][2] - m[1][2] * m[3][1];
    double A2312 = m[1][2] * m[2][3] - m[1][3] * m[2][2];
    double A1312 = m[1][1] * m[2][3] - m[1][3] * m[2][1];
    double A1212 = m[1][1] * m[2][2] - m[1][2] * m[2][1];
    double A0313 = m[1][0] * m[3][3] - m[1][3] * m[3][0];
    double A0213 = m[1][0] * m[3][2] - m[1][2] * m[3][0];
    double A0312 = m[1][0] * m[2][3] - m[1][3] * m[2][0];
    double A0212 = m[1][0] * m[2][2] - m[1][2] * m[2][0];
    double A0113 = m[1][0] * m[3][1] - m[1][1] * m[3][0];
    double A0112 = m[1][0] * m[2][1] - m[1][1] * m[2][0];

    double det = m[0][0] * (m[1][1] * A2323 - m[1][2] * A1323 + m[1][3] * A1223)
               - m[0][1] * (m[1][0] * A2323 - m[1][2] * A0323 + m[1][3] * A0223)
               + m[0][2] * (m[1][0] * A1323 - m[1][1] * A0323 + m[1][3] * A0123)
               - m[0][3] * (m[1][0] * A1223 - m[1][1] * A0223 + m[1][2] * A0123);
    
    if (std::abs(det) < 1e-10) {
        // Degenerate metric - fall back to Minkowski
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                g_inv[i][j] = (i == j) ? (i == 0 ? -1.0 : 1.0) : 0.0;
        return;
    }
    
    double invDet = 1.0 / det;
    g_inv[0][0] =  invDet * (m[1][1] * A2323 - m[1][2] * A1323 + m[1][3] * A1223);
    g_inv[0][1] = -invDet * (m[0][1] * A2323 - m[0][2] * A1323 + m[0][3] * A1223);
    g_inv[0][2] =  invDet * (m[0][1] * A2313 - m[0][2] * A1313 + m[0][3] * A1213);
    g_inv[0][3] = -invDet * (m[0][1] * A2312 - m[0][2] * A1312 + m[0][3] * A1212);
    g_inv[1][0] = -invDet * (m[1][0] * A2323 - m[1][2] * A0323 + m[1][3] * A0223);
    g_inv[1][1] =  invDet * (m[0][0] * A2323 - m[0][2] * A0323 + m[0][3] * A0223);
    g_inv[1][2] = -invDet * (m[0][0] * A2313 - m[0][2] * A0313 + m[0][3] * A0213);
    g_inv[1][3] =  invDet * (m[0][0] * A2312 - m[0][2] * A0312 + m[0][3] * A0212);
    g_inv[2][0] =  invDet * (m[1][0] * A1323 - m[1][1] * A0323 + m[1][3] * A0123);
    g_inv[2][1] = -invDet * (m[0][0] * A1323 - m[0][1] * A0323 + m[0][3] * A0123);
    g_inv[2][2] =  invDet * (m[0][0] * A1313 - m[0][1] * A0313 + m[0][3] * A0113);
    g_inv[2][3] = -invDet * (m[0][0] * A1312 - m[0][1] * A0312 + m[0][3] * A0112);
    g_inv[3][0] = -invDet * (m[1][0] * A1223 - m[1][1] * A0223 + m[1][2] * A0123);
    g_inv[3][1] =  invDet * (m[0][0] * A1223 - m[0][1] * A0223 + m[0][2] * A0123);
    g_inv[3][2] = -invDet * (m[0][0] * A1213 - m[0][1] * A0213 + m[0][2] * A0113);
    g_inv[3][3] =  invDet * (m[0][0] * A1212 - m[0][1] * A0212 + m[0][2] * A0112);
}

// Helper: Compute velocity from covariant momentum k^μ = g^μν p_ν
static Vec4 velocityFromMomentum(const Vec4& p, const Metric4D& g) {
    Vec4 k;
    
    // Extract metric components
    double m[4][4];
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            m[i][j] = g(i, j).real;
    
    // Compute inverse metric
    double g_inv[4][4];
    invertMetric4x4(m, g_inv);
    
    // Compute k^μ = g^μν p_ν
    for (int mu = 0; mu < 4; mu++) {
        k(mu) = 0;
        for (int nu = 0; nu < 4; nu++) {
            k(mu) += static_cast<float>(g_inv[mu][nu] * p(nu));
        }
    }
    return k;
}

// Helper: compute dp_μ/dλ = -(1/2)(∂g^αβ/∂x^μ) p_α p_β
// Using the identity: (∂g^αβ/∂x^μ) = -g^αρ g^βσ (∂g_ρσ/∂x^μ)
// So: dp_μ/dλ = (1/2) g^αρ g^βσ (∂g_ρσ/∂x^μ) p_α p_β
//             = (1/2) (∂g_ρσ/∂x^μ) k^ρ k^σ 
// where k^μ = g^μν p_ν is the contravariant velocity
static Vec4 momentumDerivative(const Vec4& p, const Vec4& k, const Tensor<Dual<double>, 4, 4, 4>& dg) {
    Vec4 dp;
    
    for (int mu = 0; mu < 4; mu++) {
        double sum = 0.0;
        for (int rho = 0; rho < 4; rho++) {
            for (int sigma = 0; sigma < 4; sigma++) {
                // Metric derivatives are stored in the REAL part of Dual<double>
                // (confirmed by debug output: dg(1,0,0).real=0.02 .dual=0)
                sum += dg(mu, rho, sigma).real * k(rho) * k(sigma);
            }
        }
        dp(mu) = static_cast<float>(0.5 * sum);
    }
    return dp;
}

// Forward declarations for helpers defined later
static void evaluateRK45Stage(const Vec4& x, const Vec4& p, IMetric* metric, Vec4& k_x, Vec4& k_p);
static Vec4 computeMomentum(const Vec4& velocity, const Metric4D& g);

bool Geodesic::integrateStep(Lightray& ray, IMetric* metric, float min_step, float max_step) {
    Vec4 x0 = ray.position;
    float h = ray.step_size;
    
    // Get metric and compute initial momentum
    Metric4D g0;
    Tensor<Dual<double>, 4, 4, 4> dg0;
    metric->evaluate(x0, g0, dg0);
    Vec4 p0 = computeMomentum(ray.velocity, g0);
    Vec4 k0 = ray.velocity;
    
    // RK4 stages using helper
    Vec4 k1_x = k0, k1_p = momentumDerivative(p0, k0, dg0);
    
    Vec4 x1 = x0 + k1_x * (0.5f * h), p1 = p0 + k1_p * (0.5f * h);
    Vec4 k2_x, k2_p;
    evaluateRK45Stage(x1, p1, metric, k2_x, k2_p);
    
    Vec4 x2 = x0 + k2_x * (0.5f * h), p2 = p0 + k2_p * (0.5f * h);
    Vec4 k3_x, k3_p;
    evaluateRK45Stage(x2, p2, metric, k3_x, k3_p);
    
    Vec4 x3 = x0 + k3_x * h, p3 = p0 + k3_p * h;
    Vec4 k4_x, k4_p;
    evaluateRK45Stage(x3, p3, metric, k4_x, k4_p);
    
    // Combine RK4
    Vec4 new_position = x0 + (k1_x + k2_x * 2.0f + k3_x * 2.0f + k4_x) * (h / 6.0f);
    Vec4 new_momentum = p0 + (k1_p + k2_p * 2.0f + k3_p * 2.0f + k4_p) * (h / 6.0f);
    
    // Get new velocity
    Metric4D g_new;
    Tensor<Dual<double>, 4, 4, 4> dg_new;
    metric->evaluate(new_position, g_new, dg_new);
    Vec4 new_velocity = velocityFromMomentum(new_momentum, g_new);
    
    // Adaptive step control
    float velocity_change = (new_velocity - k0).length();
    float position_change = (new_position - x0).length();
    const float target_velocity_change = 0.01f, max_position_change = 0.1f;
    
    if (velocity_change > target_velocity_change * 2.0f || position_change > max_position_change) {
        ray.step_size = std::max(ray.step_size * 0.5f, min_step);
        if (ray.step_size <= min_step) { ray.terminated = 5; return false; }
        return false;
    }
    if (velocity_change < target_velocity_change * 0.5f && ray.step_size < max_step) {
        ray.step_size = std::min(ray.step_size * 1.2f, max_step);
    }
    
    // Update ray state
    ray.position = new_position;
    ray.velocity = new_velocity;
    ray.acceleration = calculateAcceleration(new_velocity, new_position, metric);
    ray.proper_time += h;
    ray.coordinate_time += h * static_cast<float>(std::abs(new_velocity(0)));
    ray.running_dlambda_dnew *= (1.0f + velocity_change * 0.1f);
    return true;
}

Vec4 Geodesic::calculateAcceleration(const Vec4& velocity, const Vec4& position, IMetric* metric) {
    Metric4D g;
    Tensor<Dual<double>, 4, 4, 4> dg;
    metric->evaluate(position, g, dg);
    
    ChristoffelSymbols gamma = TensorOps::christoffel(g, dg);
    
    return TensorOps::geodesicAcceleration(velocity, gamma);
}

bool Geodesic::checkTermination(const Lightray& ray, IMetric* metric) {
    // =========================================================================
    // EARLY RAY TERMINATION (Task 2.2 - Performance Optimization)
    // =========================================================================
    // Terminates rays early when they are clearly:
    // - Escaping to infinity (r > r_escape AND dr/dλ > 0)
    // - Captured by horizon (r < r_horizon + ε)
    // - Hit background (for flat/asymptotic regions)
    // =========================================================================
    
    const float r = ray.position(1);
    const float dr_dlambda = ray.velocity(1);  // Radial velocity component
    
    // -------------------------------------------------------------------------
    // 1. ESCAPE DETECTION
    // -------------------------------------------------------------------------
    // A ray is escaping if:
    // - It's far from the source (r > r_escape)
    // - AND moving outward (dr/dλ > 0)
    // 
    // At large r, the geodesic is essentially straight - no need to continue
    // -------------------------------------------------------------------------
    constexpr float R_ESCAPE = 50.0f;  // Escape radius threshold
    constexpr float R_BACKGROUND = 100.0f;  // Background hit radius
    
    if (r > R_ESCAPE && dr_dlambda > 0.0f) {
        // Ray is escaping - terminate early
        return true;
    }
    
    if (r > R_BACKGROUND) {
        // Definitely hit background regardless of direction
        return true;
    }
    
    // -------------------------------------------------------------------------
    // 2. HORIZON CAPTURE DETECTION
    // -------------------------------------------------------------------------
    // For black hole metrics, rays approaching the horizon will never escape.
    // We can terminate when r < r_horizon + ε
    //
    // For Schwarzschild: r_s = 2M (we use M=1, so r_s = 2)
    // For Kerr: r_+ = M + sqrt(M² - a²) (outer horizon)
    // -------------------------------------------------------------------------
    constexpr float R_HORIZON_SCHWARZSCHILD = 2.0f;  // For M=1
    constexpr float HORIZON_EPSILON = 0.1f;  // Safety margin
    
    if (r < R_HORIZON_SCHWARZSCHILD + HORIZON_EPSILON) {
        // Near or inside horizon - ray captured
        return true;
    }
    
    // -------------------------------------------------------------------------
    // 3. INTEGRATION LIMITS
    // -------------------------------------------------------------------------
    // Prevent runaway integration
    // -------------------------------------------------------------------------
    constexpr float MAX_PROPER_TIME = 100.0f;
    
    if (ray.proper_time > MAX_PROPER_TIME) {
        return true;
    }
    
    // -------------------------------------------------------------------------
    // 4. VELOCITY SANITY CHECK
    // -------------------------------------------------------------------------
    // If velocity becomes too small, ray is stuck
    // -------------------------------------------------------------------------
    if (ray.velocity.length() < 1e-8f) {
        return true;
    }
    
    // -------------------------------------------------------------------------
    // 5. THETA BOUNDARY CHECK  
    // -------------------------------------------------------------------------
    // Rays at theta ≈ 0 or theta ≈ π can have coordinate singularities
    // -------------------------------------------------------------------------
    const float theta = ray.position(2);
    constexpr float THETA_MIN = 0.01f;
    constexpr float THETA_MAX = 3.13f;  // π - 0.01
    
    if (theta < THETA_MIN || theta > THETA_MAX) {
        // Near coordinate singularity at poles
        return true;
    }
    
    return false;
}

float Geodesic::calculateRedshift(const Lightray& ray, const ObserverState& observer, IMetric* metric) {
    Metric4D g;
    Tensor<Dual<double>, 4, 4, 4> dg;
    metric->evaluate(ray.position, g, dg);
    
    Vec4 observer_lower = TensorOps::lowerIndex(observer.velocity, g);
    
    double dot_product = 0.0;
    for (int mu = 0; mu < 4; mu++) {
        dot_product += ray.velocity(mu) * observer_lower(mu);
    }
    
    return static_cast<float>(dot_product / ray.ku_uobsu);
}

ObserverState Geodesic::createObserver(const Vec4& position, const Vec4& velocity, IMetric* metric) {
    ObserverState observer;
    observer.position = position;
    observer.velocity = velocity;
    
    Metric4D g;
    Tensor<Dual<double>, 4, 4, 4> dg;
    metric->evaluate(position, g, dg);
    double velocity_norm = TensorOps::innerProduct(velocity, velocity, g);
    observer.is_timelike = velocity_norm < 0.0;
    
    if (!observer.is_timelike) {
        observer.velocity = Vec4(); // Default constructor initializes to zeros
        observer.velocity(0) = 1.0f; // Set time component
        velocity_norm = TensorOps::innerProduct(observer.velocity, observer.velocity, g);
    }
    
    float normalization = static_cast<float>(1.0 / sqrt(-velocity_norm));
    observer.velocity *= normalization;
    
    calculateTetrads(observer, metric);
    
    return observer;
}

void Geodesic::calculateTetrads(ObserverState& observer, IMetric* metric) {
    Metric4D g;
    Tensor<Dual<double>, 4, 4, 4> dg;
    metric->evaluate(observer.position, g, dg);
    
    observer.e0 = observer.velocity;
    
    observer.e1 = Vec4(); observer.e1(1) = 1.0f;
    Vec4 e1_parallel = TensorOps::lowerIndex(observer.e1, g);
    float e1_dot_e0 = 0.0f;
    for (int mu = 0; mu < 4; mu++) {
        e1_dot_e0 += observer.e0(mu) * e1_parallel(mu);
    }
    observer.e1 = observer.e1 - observer.e0 * e1_dot_e0;
    
    observer.e2 = Vec4(); observer.e2(2) = 1.0f;
    Vec4 e2_parallel = TensorOps::lowerIndex(observer.e2, g);
    float e2_dot_e0 = 0.0f;
    float e2_dot_e1 = 0.0f;
    for (int mu = 0; mu < 4; mu++) {
        e2_dot_e0 += observer.e0(mu) * e2_parallel(mu);
        e2_dot_e1 += observer.e1(mu) * e2_parallel(mu);
    }
    observer.e2 = observer.e2 - observer.e0 * e2_dot_e0 - observer.e1 * e2_dot_e1;
    
    observer.e3 = Vec4(); observer.e3(3) = 1.0f;
    Vec4 e3_parallel = TensorOps::lowerIndex(observer.e3, g);
    float e3_dot_e0 = 0.0f;
    float e3_dot_e1 = 0.0f;
    float e3_dot_e2 = 0.0f;
    for (int mu = 0; mu < 4; mu++) {
        e3_dot_e0 += observer.e0(mu) * e3_parallel(mu);
        e3_dot_e1 += observer.e1(mu) * e3_parallel(mu);
        e3_dot_e2 += observer.e2(mu) * e3_parallel(mu);
    }
    observer.e3 = observer.e3 - observer.e0 * e3_dot_e0 - observer.e1 * e3_dot_e1 - observer.e2 * e3_dot_e2;
    
    float e1_norm = static_cast<float>(sqrt(std::abs(TensorOps::innerProduct(observer.e1, observer.e1, g))));
    float e2_norm = static_cast<float>(sqrt(std::abs(TensorOps::innerProduct(observer.e2, observer.e2, g))));
    float e3_norm = static_cast<float>(sqrt(std::abs(TensorOps::innerProduct(observer.e3, observer.e3, g))));
    
    if (e1_norm > 1e-10f) observer.e1 /= e1_norm;
    if (e2_norm > 1e-10f) observer.e2 /= e2_norm;
    if (e3_norm > 1e-10f) observer.e3 /= e3_norm;
}

// =============================================================================
// RK45 Adaptive Integration Implementation (Phase 2)
// =============================================================================

IntegratorConfig Geodesic::getDefaultConfig() {
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

float Geodesic::computeOptimalStep(float h, float error, float tolerance,
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

// Helper: Compute covariant momentum from velocity p_μ = g_μν k^ν
static Vec4 computeMomentum(const Vec4& velocity, const Metric4D& g) {
    Vec4 p;
    for (int mu = 0; mu < 4; mu++) {
        p(mu) = 0;
        for (int nu = 0; nu < 4; nu++) {
            p(mu) += g(mu, nu).real * velocity(nu);
        }
    }
    return p;
}

// Helper: Evaluate RK45 stage and return derivatives
static void evaluateRK45Stage(const Vec4& x, const Vec4& p, IMetric* metric,
                               Vec4& k_x, Vec4& k_p) {
    Metric4D g;
    Tensor<Dual<double>, 4, 4, 4> dg;
    metric->evaluate(x, g, dg);
    k_x = velocityFromMomentum(p, g);
    k_p = momentumDerivative(p, k_x, dg);
}

// Helper: Compute RK45 error norm
static float computeRK45ErrorNorm(const Vec4& error_x, const Vec4& new_position,
                                   const IntegratorConfig& config) {
    float error_norm = 0.0f;
    for (int i = 0; i < 4; i++) {
        float scale = config.abs_tolerance + config.rel_tolerance * std::abs(new_position(i));
        float err_i = std::abs(error_x(i)) / scale;
        error_norm += err_i * err_i;
    }
    return std::sqrt(error_norm / 4.0f);
}

// Helper: Check for NaN/Inf in ray state
static bool hasInvalidState(const Vec4& position, const Vec4& velocity) {
    for (int i = 0; i < 4; i++) {
        if (std::isnan(position(i)) || std::isinf(position(i)) ||
            std::isnan(velocity(i)) || std::isinf(velocity(i))) {
            return true;
        }
    }
    return false;
}

bool Geodesic::integrateStepRK45(Lightray& ray, IMetric* metric, const IntegratorConfig& config) {
    using namespace DP45;
    
    Vec4 x0 = ray.position;
    float h = ray.step_size;
    
    // Get metric and compute initial momentum
    Metric4D g0;
    Tensor<Dual<double>, 4, 4, 4> dg0;
    metric->evaluate(x0, g0, dg0);
    Vec4 p0 = computeMomentum(ray.velocity, g0);
    Vec4 k0 = ray.velocity;
    
    // Stage 1
    Vec4 k1_x = k0, k1_p = momentumDerivative(p0, k0, dg0);
    
    // Stage 2-6: Use helper for each stage
    Vec4 k2_x, k2_p, k3_x, k3_p, k4_x, k4_p, k5_x, k5_p, k6_x, k6_p;
    
    Vec4 x2 = x0 + k1_x * static_cast<float>(a21 * h);
    Vec4 p2 = p0 + k1_p * static_cast<float>(a21 * h);
    evaluateRK45Stage(x2, p2, metric, k2_x, k2_p);
    
    Vec4 x3 = x0 + k1_x * static_cast<float>(a31 * h) + k2_x * static_cast<float>(a32 * h);
    Vec4 p3 = p0 + k1_p * static_cast<float>(a31 * h) + k2_p * static_cast<float>(a32 * h);
    evaluateRK45Stage(x3, p3, metric, k3_x, k3_p);
    
    Vec4 x4 = x0 + k1_x * static_cast<float>(a41 * h) + k2_x * static_cast<float>(a42 * h) + k3_x * static_cast<float>(a43 * h);
    Vec4 p4 = p0 + k1_p * static_cast<float>(a41 * h) + k2_p * static_cast<float>(a42 * h) + k3_p * static_cast<float>(a43 * h);
    evaluateRK45Stage(x4, p4, metric, k4_x, k4_p);
    
    Vec4 x5 = x0 + k1_x * static_cast<float>(a51 * h) + k2_x * static_cast<float>(a52 * h) + k3_x * static_cast<float>(a53 * h) + k4_x * static_cast<float>(a54 * h);
    Vec4 p5 = p0 + k1_p * static_cast<float>(a51 * h) + k2_p * static_cast<float>(a52 * h) + k3_p * static_cast<float>(a53 * h) + k4_p * static_cast<float>(a54 * h);
    evaluateRK45Stage(x5, p5, metric, k5_x, k5_p);
    
    Vec4 x6 = x0 + k1_x * static_cast<float>(a61 * h) + k2_x * static_cast<float>(a62 * h) + k3_x * static_cast<float>(a63 * h) + k4_x * static_cast<float>(a64 * h) + k5_x * static_cast<float>(a65 * h);
    Vec4 p6 = p0 + k1_p * static_cast<float>(a61 * h) + k2_p * static_cast<float>(a62 * h) + k3_p * static_cast<float>(a63 * h) + k4_p * static_cast<float>(a64 * h) + k5_p * static_cast<float>(a65 * h);
    evaluateRK45Stage(x6, p6, metric, k6_x, k6_p);
    
    // 5th order solution
    Vec4 new_position = x0 + (k1_x * static_cast<float>(b1) + k3_x * static_cast<float>(b3) + k4_x * static_cast<float>(b4) + k5_x * static_cast<float>(b5) + k6_x * static_cast<float>(b6)) * h;
    Vec4 new_momentum = p0 + (k1_p * static_cast<float>(b1) + k3_p * static_cast<float>(b3) + k4_p * static_cast<float>(b4) + k5_p * static_cast<float>(b5) + k6_p * static_cast<float>(b6)) * h;
    
    // Stage 7 (FSAL)
    Vec4 k7_x, k7_p;
    evaluateRK45Stage(new_position, new_momentum, metric, k7_x, k7_p);
    Vec4 new_velocity = k7_x;
    
    // Error estimation
    Vec4 error_x = (k1_x * static_cast<float>(e1) + k3_x * static_cast<float>(e3) + k4_x * static_cast<float>(e4) + k5_x * static_cast<float>(e5) + k6_x * static_cast<float>(e6) + k7_x * static_cast<float>(e7)) * h;
    float error_norm = computeRK45ErrorNorm(error_x, new_position, config);
    
    // Step acceptance
    if (error_norm > 1.0f) {
        ray.step_size = computeOptimalStep(h, error_norm, 1.0f, config);
        if (ray.step_size <= config.min_step) { ray.terminated = 5; return false; }
        return false;
    }
    
    // Update ray state
    ray.position = new_position;
    ray.velocity = new_velocity;
    ray.acceleration = calculateAcceleration(new_velocity, new_position, metric);
    ray.proper_time += h;
    ray.coordinate_time += h * static_cast<float>(std::abs(new_velocity(0)));
    ray.step_size = computeOptimalStep(h, error_norm, 1.0f, config);
    
    if (hasInvalidState(new_position, new_velocity)) { ray.terminated = 3; return false; }
    return true;
}