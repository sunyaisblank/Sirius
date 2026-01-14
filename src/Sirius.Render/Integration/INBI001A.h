// KNBI001A.h - Beam State and Integration Structures
// Component ID: KNBI001A
// Purpose: Beam propagation with Jacobian tracking for offline rendering
//
// MATHEMATICAL BASIS:
// A beam is a bundle of neighbouring geodesics. The Jacobian J^i_j = ∂x^i/∂x^j_0
// tracks how the beam spreads as it propagates through curved spacetime.
//
// The geodesic deviation equation:
//   D²ξ^μ/dλ² + R^μ_νρσ k^ν ξ^ρ k^σ = 0
//
// Related work: James et al. (2015) "DNGR" Appendix B
//
// MEMORY LAYOUT: 320 bytes per beam (optimised for GPU cache lines)

#pragma once

#include "../Sirius.Math/MTTP001A.h"
#include "../Sirius.Physics/Core/PHMT000B.h"
#include <cmath>

// CPU/GPU compatibility macro
#ifdef __CUDACC__
#define SIRIUS_HOST_DEVICE __host__ __device__
#else
#define SIRIUS_HOST_DEVICE
#endif

namespace sirius::kernel {

using namespace sirius::math;
using namespace sirius::physics;

//==============================================================================
// BeamStateD: Complete state for beam propagation
//==============================================================================

struct BeamStateD {
    // Central geodesic state (64 bytes)
    Vec4d x;        // Position: (t, r, θ, φ)
    Vec4d k;        // Wave 4-vector (covariant): (k_t, k_r, k_θ, k_φ)
    
    // Jacobian matrix J^i_j = ∂x^i/∂x^j_0 (128 bytes)
    // Maps initial conditions to current position
    double J[4][4];
    
    // Jacobian velocity dJ^i_j/dλ (128 bytes)
    // Needed for second-order geodesic deviation equation
    double dJ[4][4];
    
    // Affine parameter
    double lambda;
    
    // Conserved quantities (computed at initialisation)
    double E;       // Energy: E = -k_t
    double Lz;      // Angular momentum: L_z = k_φ
    double Q;       // Carter constant (Kerr only)
    
    // Derived beam properties (cached)
    double solidAngle;      // Beam solid angle on sky
    double magnification;   // |det(J_angular)|^{-1}
    double majorAxis;       // Ellipse semi-major axis [rad]
    double minorAxis;       // Ellipse semi-minor axis [rad]
    double orientation;     // PA of major axis [rad]
    
    // Beam status
    bool terminated;        // Hit horizon or escaped
    bool atCaustic;         // det(J) ≈ 0 (magnification → ∞)
    
    // Initial pixel solid angle (set at creation)
    double initialPixelSolidAngle;
    
    //--------------------------------------------------------------------------
    // Initialisation
    //--------------------------------------------------------------------------
    
    SIRIUS_HOST_DEVICE void initialise() {
        // Identity Jacobian (beam starts as a point)
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                J[i][j] = (i == j) ? 1.0 : 0.0;
                dJ[i][j] = 0.0;
            }
        }
        
        lambda = 0;
        solidAngle = 0;
        magnification = 1.0;
        majorAxis = 0;
        minorAxis = 0;
        orientation = 0;
        terminated = false;
        atCaustic = false;
        initialPixelSolidAngle = 0;
    }
    
    //--------------------------------------------------------------------------
    // Extract beam geometry from Jacobian
    //--------------------------------------------------------------------------
    
    SIRIUS_HOST_DEVICE void updateGeometry() {
        // Extract 2×2 angular submatrix (θ,φ) → (θ₀,φ₀) mapping
        // This gives the beam ellipse on the sky
        
        double a = J[2][2];  // ∂θ/∂θ₀
        double b = J[2][3];  // ∂θ/∂φ₀
        double c = J[3][2];  // ∂φ/∂θ₀
        double d = J[3][3];  // ∂φ/∂φ₀
        
        // Determinant of angular submatrix
        double det = a*d - b*c;
        
        // Check for caustic (det → 0)
        if (std::abs(det) < 1e-12) {
            atCaustic = true;
            magnification = 1e12;  // Cap at large value
        } else {
            atCaustic = false;
            magnification = 1.0 / std::abs(det);
        }
        
        // Singular value decomposition for ellipse axes
        // For matrix [[a,b],[c,d]], singular values are:
        // σ = √[(p ± √(p² - 4q²))/2] where p = a² + b² + c² + d², q = det
        
        double p = a*a + b*b + c*c + d*d;
        double q = det;
        double disc = p*p - 4*q*q;
        
        if (disc < 0) disc = 0;  // Numerical protection
        double s = std::sqrt(disc);
        
        majorAxis = std::sqrt(std::max(0.0, (p + s) / 2));
        minorAxis = std::sqrt(std::max(0.0, (p - s) / 2));
        
        // Orientation: angle of major axis from θ direction
        // From SVD: tan(2φ) = 2(ab + cd) / (a² + b² - c² - d²)
        double num = 2*(a*b + c*d);
        double den = a*a + b*b - c*c - d*d;
        orientation = 0.5 * std::atan2(num, den);
        
        // Solid angle = π × σ_max × σ_min × (initial pixel solid angle)
        solidAngle = M_PI * majorAxis * minorAxis * initialPixelSolidAngle;
    }
    
    //--------------------------------------------------------------------------
    // Convert from GeodesicStateD
    //--------------------------------------------------------------------------
    
    SIRIUS_HOST_DEVICE void fromGeodesic(const GeodesicStateD& geo) {
        x = geo.x;
        k = geo.k;
        lambda = geo.lambda;
        E = geo.E;
        Lz = geo.Lz;
        Q = geo.Q;
    }
    
    //--------------------------------------------------------------------------
    // Convert to GeodesicStateD
    //--------------------------------------------------------------------------
    
    SIRIUS_HOST_DEVICE GeodesicStateD toGeodesic() const {
        GeodesicStateD geo;
        geo.x = x;
        geo.k = k;
        geo.lambda = lambda;
        geo.E = E;
        geo.Lz = Lz;
        geo.Q = Q;
        return geo;
    }
};

//==============================================================================
// BeamIntegratorD: Simultaneous geodesic and Jacobian integration
//==============================================================================

class BeamIntegratorD {
public:
    struct Config {
        double stepSize = 0.1;
        double minStepSize = 1e-6;
        double maxStepSize = 1.0;
        int maxSteps = 100000;
        double escapeRadius = 1e6;
        double causticThreshold = 1e-12;
    };
    
    explicit BeamIntegratorD(const IMetricD* metric, const Config& config = Config())
        : m_metric(metric), m_config(config) {}
    
    //--------------------------------------------------------------------------
    // Single integration step
    //--------------------------------------------------------------------------
    
    SIRIUS_HOST_DEVICE bool step(BeamStateD& beam, double h) const {
        if (beam.terminated) return false;
        
        // Check if outside valid region
        if (!m_metric->isValid(beam.x)) {
            beam.terminated = true;
            return false;
        }
        
        // Check if escaped
        if (beam.x.r > m_config.escapeRadius) {
            beam.terminated = true;
            return false;
        }
        
        // ============================================================
        // GEODESIC INTEGRATION: Symplectic leapfrog (Störmer-Verlet)
        // Uses the metric's dHdq for momentum evolution
        // ============================================================
        
        // Get metric and inverse metric
        double g[4][4], g_inv[4][4];
        m_metric->evaluate(beam.x, g, g_inv);
        
        // Contravariant wave vector k^μ = g^μν k_ν (for position update)
        auto computeKup = [&](const Vec4d& k_cov, double* k_up) {
            for (int mu = 0; mu < 4; ++mu) {
                k_up[mu] = 0;
                for (int nu = 0; nu < 4; ++nu) {
                    k_up[mu] += g_inv[mu][nu] * k_cov[nu];
                }
            }
        };
        
        double k_up[4];
        computeKup(beam.k, k_up);
        
        // Leapfrog step (kick-drift-kick):
        // 1. Half-step momentum: k' = k + (h/2) dk/dλ
        Vec4d dk_dlambda = m_metric->dHdq(beam.x, beam.k);
        Vec4d k_half = beam.k;
        for (int mu = 0; mu < 4; ++mu) {
            k_half[mu] += 0.5 * h * dk_dlambda[mu];
        }
        
        // Update g_inv for midpoint
        m_metric->evaluate(beam.x, g, g_inv);
        double k_up_half[4];
        computeKup(k_half, k_up_half);
        
        // 2. Full-step position: x' = x + h * k^μ
        Vec4d x_new = beam.x;
        for (int mu = 0; mu < 4; ++mu) {
            x_new[mu] += h * k_up_half[mu];
        }
        
        // 3. Half-step momentum at new position: k'' = k' + (h/2) dk/dλ|_x'
        Vec4d dk_dlambda_new = m_metric->dHdq(x_new, k_half);
        Vec4d k_new = k_half;
        for (int mu = 0; mu < 4; ++mu) {
            k_new[mu] += 0.5 * h * dk_dlambda_new[mu];
        }
        
        // ============================================================
        // JACOBIAN INTEGRATION: RK4 for geodesic deviation equation
        // d²J/dλ² = -R·k·k·J (Jacobi equation)
        // Converted to first-order: dY/dλ = f(Y) where Y = (J, dJ)
        // ============================================================
        
        // Get Riemann tensor at current position
        double R[4][4][4][4];
        m_metric->riemann(beam.x, R);
        
        // RK4 coefficients for Jacobian
        double J_k1[4][4], dJ_k1[4][4];
        double J_k2[4][4], dJ_k2[4][4];
        double J_k3[4][4], dJ_k3[4][4];
        double J_k4[4][4], dJ_k4[4][4];
        
        // Lambda to compute dJ/dλ and d²J/dλ² given current state
        auto computeJacobianRHS = [&](const double J_in[4][4], const double dJ_in[4][4],
                                      const double* kup, double J_out[4][4], double dJ_out[4][4]) {
            for (int i = 0; i < 4; ++i) {
                for (int j = 0; j < 4; ++j) {
                    // dJ/dλ = dJ (trivial)
                    J_out[i][j] = dJ_in[i][j];
                    
                    // d²J/dλ² = -R^i_νρσ k^ν k^σ J^ρ_j
                    double ddJ = 0;
                    for (int nu = 0; nu < 4; ++nu) {
                        for (int rho = 0; rho < 4; ++rho) {
                            for (int sigma = 0; sigma < 4; ++sigma) {
                                ddJ -= R[i][nu][rho][sigma] * kup[nu] * kup[sigma] * J_in[rho][j];
                            }
                        }
                    }
                    dJ_out[i][j] = ddJ;
                }
            }
        };
        
        // RK4 stage 1: k1 = f(y_n)
        computeJacobianRHS(beam.J, beam.dJ, k_up, J_k1, dJ_k1);
        
        // RK4 stage 2: k2 = f(y_n + h/2 * k1)
        double J_tmp[4][4], dJ_tmp[4][4];
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                J_tmp[i][j] = beam.J[i][j] + 0.5 * h * J_k1[i][j];
                dJ_tmp[i][j] = beam.dJ[i][j] + 0.5 * h * dJ_k1[i][j];
            }
        }
        computeJacobianRHS(J_tmp, dJ_tmp, k_up, J_k2, dJ_k2);
        
        // RK4 stage 3: k3 = f(y_n + h/2 * k2)
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                J_tmp[i][j] = beam.J[i][j] + 0.5 * h * J_k2[i][j];
                dJ_tmp[i][j] = beam.dJ[i][j] + 0.5 * h * dJ_k2[i][j];
            }
        }
        computeJacobianRHS(J_tmp, dJ_tmp, k_up, J_k3, dJ_k3);
        
        // RK4 stage 4: k4 = f(y_n + h * k3)
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                J_tmp[i][j] = beam.J[i][j] + h * J_k3[i][j];
                dJ_tmp[i][j] = beam.dJ[i][j] + h * dJ_k3[i][j];
            }
        }
        computeJacobianRHS(J_tmp, dJ_tmp, k_up, J_k4, dJ_k4);
        
        // Final RK4 update: y_{n+1} = y_n + h/6 * (k1 + 2k2 + 2k3 + k4)
        double J_new[4][4], dJ_new[4][4];
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                J_new[i][j] = beam.J[i][j] + (h/6.0) * (J_k1[i][j] + 2*J_k2[i][j] + 2*J_k3[i][j] + J_k4[i][j]);
                dJ_new[i][j] = beam.dJ[i][j] + (h/6.0) * (dJ_k1[i][j] + 2*dJ_k2[i][j] + 2*dJ_k3[i][j] + dJ_k4[i][j]);
            }
        }
        
        // ============================================================
        // Update beam state
        // ============================================================
        
        beam.x = x_new;
        beam.k = k_new;
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                beam.J[i][j] = J_new[i][j];
                beam.dJ[i][j] = dJ_new[i][j];
            }
        }
        beam.lambda += h;
        
        // Update beam geometry
        beam.updateGeometry();
        
        return true;
    }
    
private:
    const IMetricD* m_metric;
    Config m_config;
};

} // namespace sirius::kernel
