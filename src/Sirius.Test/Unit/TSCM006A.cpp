// TSCM006A.cpp - G-Factor Validation Tests (DNGR Phase 6.2)
// Tests: DNGR Eq A.16 blue shift formula, SR/GR factors, frame dragging.

#include <gtest/gtest.h>
#include <cmath>

namespace sirius::test {

constexpr double kEpsilon = 1e-6;
constexpr double kPi = 3.14159265358979323846;

// =============================================================================
// Test Fixture
// =============================================================================

class GFactorTests : public ::testing::Test {
protected:
    // DNGR blue shift formula (Eq A.16) reference implementation
    double computeDNGRBlueShift(
        double beta,     // Camera speed relative to FIDO
        double N_y,      // Ray direction component along camera velocity
        double b,        // Ray angular momentum
        double omega,    // Frame dragging angular velocity
        double alpha)    // Lapse function
    {
        // Special relativistic factor
        double sr_factor = std::sqrt(1.0 - beta * beta) / (1.0 - beta * N_y);
        
        // Gravitational factor
        double gr_factor = (1.0 - b * omega) / alpha;
        
        return sr_factor * gr_factor;
    }
    
    // Invariant g-factor: g = (k·u_obs) / (k·u_emit)
    // For observer at infinity: u_obs = (1, 0, 0, 0)
    double computeInvariantG(
        const double k[4],      // Photon 4-momentum
        const double u_emit[4], // Emitter 4-velocity
        const double g_metric[4][4]) // Metric tensor
    {
        // k·u_obs = g_tt * k^t + g_tφ * k^φ (for Kerr)
        // Simplified: for static observer at infinity, k·u_obs ≈ -k^t
        double u_obs[4] = {1.0, 0.0, 0.0, 0.0};
        
        double k_dot_u_obs = 0.0;
        double k_dot_u_emit = 0.0;
        
        for (int mu = 0; mu < 4; ++mu) {
            for (int nu = 0; nu < 4; ++nu) {
                k_dot_u_obs += g_metric[mu][nu] * k[mu] * u_obs[nu];
                k_dot_u_emit += g_metric[mu][nu] * k[mu] * u_emit[nu];
            }
        }
        
        return k_dot_u_obs / k_dot_u_emit;
    }
    
    // Compute Kerr metric at given position
    void computeKerrMetric(
        double r, double theta, double M, double a,
        double g[4][4])
    {
        double r2 = r * r;
        double a2 = a * a;
        double cos_theta = std::cos(theta);
        double sin_theta = std::sin(theta);
        double cos2 = cos_theta * cos_theta;
        double sin2 = sin_theta * sin_theta;
        
        double Sigma = r2 + a2 * cos2;
        double Delta = r2 - 2.0 * M * r + a2;
        
        // Initialize to zero
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
                g[i][j] = 0.0;
        
        // g_tt
        g[0][0] = -(1.0 - 2.0 * M * r / Sigma);
        
        // g_tφ = g_φt (frame dragging)
        g[0][3] = -2.0 * M * a * r * sin2 / Sigma;
        g[3][0] = g[0][3];
        
        // g_rr
        g[1][1] = Sigma / Delta;
        
        // g_θθ
        g[2][2] = Sigma;
        
        // g_φφ
        g[3][3] = (r2 + a2 + 2.0 * M * a2 * r * sin2 / Sigma) * sin2;
    }
    
    void SetUp() override {}
    void TearDown() override {}
};

// =============================================================================
// DNGR Formula Component Tests
// =============================================================================

// Test: Static camera (β=0) gives pure gravitational redshift
TEST_F(GFactorTests, StaticCameraGravitationalRedshift) {
    double beta = 0.0;
    double N_y = 0.0;
    double b = 0.0;  // Radial ray
    double omega = 0.0;  // Schwarzschild
    double alpha = 0.8;  // Typical lapse at r = 6M
    
    double g = computeDNGRBlueShift(beta, N_y, b, omega, alpha);
    
    // For static camera, SR factor = 1, so g = 1/α
    EXPECT_NEAR(g, 1.0 / alpha, kEpsilon)
        << "Static camera should give g = 1/α";
}

// Test: SR factor is 1 when β = 0
TEST_F(GFactorTests, SRFactorOneWhenStationary) {
    double beta = 0.0;
    double N_y = 0.5;  // Any direction
    
    double sr_factor = std::sqrt(1.0 - beta * beta) / (1.0 - beta * N_y);
    
    EXPECT_NEAR(sr_factor, 1.0, kEpsilon)
        << "SR factor should be 1 for stationary camera";
}

// Test: Approaching ray (N_y = β direction) gives blueshift
TEST_F(GFactorTests, ApproachingRayBlueshift) {
    double beta = 0.3;
    double N_y = 1.0;  // Ray coming from direction camera is moving toward
    double alpha = 1.0;  // Flat space for isolation
    double omega = 0.0;
    double b = 0.0;
    
    double g = computeDNGRBlueShift(beta, N_y, b, omega, alpha);
    
    // Classical Doppler: approaching gives g > 1
    EXPECT_GT(g, 1.0) << "Approaching ray should be blueshifted (g > 1)";
}

// Test: Receding ray (N_y = -β direction) gives redshift
TEST_F(GFactorTests, RecedingRayRedshift) {
    double beta = 0.3;
    double N_y = -1.0;  // Ray coming from direction camera is moving away from
    double alpha = 1.0;
    double omega = 0.0;
    double b = 0.0;
    
    double g = computeDNGRBlueShift(beta, N_y, b, omega, alpha);
    
    // Classical Doppler: receding gives g < 1
    EXPECT_LT(g, 1.0) << "Receding ray should be redshifted (g < 1)";
}

// Test: Transverse ray (N_y = 0) gives transverse Doppler
TEST_F(GFactorTests, TransverseRayDoppler) {
    double beta = 0.5;
    double N_y = 0.0;  // Perpendicular to camera motion
    double alpha = 1.0;
    double omega = 0.0;
    double b = 0.0;
    
    double g = computeDNGRBlueShift(beta, N_y, b, omega, alpha);
    
    // Transverse Doppler: g = √(1 - β²) = 1/γ < 1 (redshift)
    double expected = std::sqrt(1.0 - beta * beta);
    EXPECT_NEAR(g, expected, kEpsilon)
        << "Transverse ray should show time dilation redshift";
}

// =============================================================================
// Frame Dragging Effects Tests
// =============================================================================

// Test: Frame dragging modifies g-factor for nonzero b and ω
TEST_F(GFactorTests, FrameDraggingEffect) {
    double beta = 0.0;  // Static camera
    double N_y = 0.0;
    double alpha = 0.9;
    
    // Case 1: No angular momentum
    double b1 = 0.0;
    double omega = 0.1;  // Moderate frame dragging
    double g1 = computeDNGRBlueShift(beta, N_y, b1, omega, alpha);
    
    // Case 2: Ray with positive angular momentum (corotating)
    double b2 = 2.0;
    double g2 = computeDNGRBlueShift(beta, N_y, b2, omega, alpha);
    
    // Corotating ray should have different g
    EXPECT_NE(g1, g2) << "Frame dragging should affect rays with angular momentum";
    
    // Specifically: (1 - bω) factor reduces g for corotating rays
    EXPECT_LT(g2, g1) << "Corotating ray should have lower g (more redshifted)";
}

// Test: Frame dragging has no effect when a = 0 (ω = 0)
TEST_F(GFactorTests, NoFrameDraggingInSchwarzschild) {
    double beta = 0.2;
    double N_y = 0.3;
    double alpha = 0.8;
    double omega = 0.0;  // Schwarzschild
    
    double b1 = 0.0;
    double b2 = 5.0;  // Large angular momentum
    
    double g1 = computeDNGRBlueShift(beta, N_y, b1, omega, alpha);
    double g2 = computeDNGRBlueShift(beta, N_y, b2, omega, alpha);
    
    EXPECT_NEAR(g1, g2, kEpsilon)
        << "Angular momentum shouldn't matter when ω = 0";
}

// =============================================================================
// Known Value Tests
// =============================================================================

// Test: g-factor at r = 6M, a = 0 (Schwarzschild ISCO)
TEST_F(GFactorTests, SchwarzschildISCOGravitationalRedshift) {
    double M = 1.0;
    double r = 6.0 * M;  // ISCO
    double a = 0.0;
    
    // Lapse at ISCO: α = √(1 - 2M/r) = √(1 - 1/3) = √(2/3)
    double alpha = std::sqrt(1.0 - 2.0 * M / r);
    double expected_alpha = std::sqrt(2.0 / 3.0);
    
    EXPECT_NEAR(alpha, expected_alpha, kEpsilon);
    
    // For static observer, g = 1/α = √(3/2) ≈ 1.225
    double g = computeDNGRBlueShift(0.0, 0.0, 0.0, 0.0, alpha);
    double expected_g = std::sqrt(3.0 / 2.0);
    
    EXPECT_NEAR(g, expected_g, kEpsilon)
        << "Gravitational redshift at ISCO should give g = √(3/2)";
}

// Test: g-factor approaches 1 at large r
TEST_F(GFactorTests, AsymptoticGFactorApproachesOne) {
    double M = 1.0;
    double a = 0.5;
    
    std::vector<double> radii = {100.0, 500.0, 1000.0};
    
    for (double r : radii) {
        // Approximate lapse at large r
        double alpha = std::sqrt(1.0 - 2.0 * M / r);
        double omega = 2.0 * a * M * r / std::pow(r * r + a * a, 2);
        
        // Static camera, radial ray
        double g = computeDNGRBlueShift(0.0, 0.0, 0.0, omega, alpha);
        
        // Should approach 1 at large r  
        // At r=100M: g ≈ 1/α = 1/0.9899 ≈ 1.0101
        EXPECT_GT(g, 0.98) << "g should approach 1 at r = " << r;
        EXPECT_LT(g, 1.02) << "g should approach 1 at r = " << r;
    }
}

// =============================================================================
// Consistency Tests
// =============================================================================

// Test: G-factor is always positive
TEST_F(GFactorTests, GFactorAlwaysPositive) {
    std::vector<double> betas = {0.0, 0.1, 0.3, 0.5, 0.8};
    std::vector<double> N_ys = {-1.0, -0.5, 0.0, 0.5, 1.0};
    std::vector<double> alphas = {0.5, 0.8, 1.0};
    
    for (double beta : betas) {
        for (double N_y : N_ys) {
            // Skip unphysical: |β N_y| >= 1 would mean superluminal
            if (std::abs(beta * N_y) >= 0.999) continue;
            
            for (double alpha : alphas) {
                double g = computeDNGRBlueShift(beta, N_y, 0.0, 0.0, alpha);
                
                EXPECT_GT(g, 0.0) 
                    << "g must be positive for β=" << beta 
                    << ", N_y=" << N_y << ", α=" << alpha;
            }
        }
    }
}

// Test: G-factor is finite for valid inputs
TEST_F(GFactorTests, GFactorIsFinite) {
    // Normal case
    double g = computeDNGRBlueShift(0.3, 0.5, 1.0, 0.05, 0.8);
    
    EXPECT_TRUE(std::isfinite(g)) << "g must be finite for valid inputs";
    EXPECT_GT(g, 0.0) << "g must be positive";
    EXPECT_LT(g, 100.0) << "g should be reasonable magnitude";
}

// =============================================================================
// Relativistic Limit Tests
// =============================================================================

// Test: High speed approaching gives large blueshift
TEST_F(GFactorTests, HighSpeedApproachingBlueshift) {
    double beta = 0.9;
    double N_y = 1.0;  // Head-on
    double alpha = 1.0;
    double omega = 0.0;
    double b = 0.0;
    
    double g = computeDNGRBlueShift(beta, N_y, b, omega, alpha);
    
    // Doppler factor for head-on: √((1+β)/(1-β)) ≈ 4.36 for β=0.9
    double expected = std::sqrt((1.0 + beta) / (1.0 - beta));
    
    EXPECT_NEAR(g, expected, 0.01)
        << "Head-on relativistic Doppler should give large blueshift";
}

// Test: High speed receding gives large redshift
TEST_F(GFactorTests, HighSpeedRecedingRedshift) {
    double beta = 0.9;
    double N_y = -1.0;  // Tail-on
    double alpha = 1.0;
    double omega = 0.0;
    double b = 0.0;
    
    double g = computeDNGRBlueShift(beta, N_y, b, omega, alpha);
    
    // Doppler factor for tail-on: √((1-β)/(1+β)) ≈ 0.23 for β=0.9
    double expected = std::sqrt((1.0 - beta) / (1.0 + beta));
    
    EXPECT_NEAR(g, expected, 0.01)
        << "Tail-on relativistic Doppler should give large redshift";
}

} // namespace sirius::test
