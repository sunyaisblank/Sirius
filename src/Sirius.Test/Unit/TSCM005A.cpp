// TSCM005A.cpp - FIDO Frame Consistency Tests (DNGR Phase 6.1)
// Tests: ρ, Δ, Σ, α, ω, ϖ quantities, FIDO basis, circular geodesic speed.

#include <gtest/gtest.h>
#include <cmath>

namespace sirius::test {

constexpr double kEpsilon = 1e-6;
constexpr double kPi = 3.14159265358979323846;

// =============================================================================
// Test Fixture
// =============================================================================

class FIDOFrameTests : public ::testing::Test {
protected:
    // Kerr metric quantities from DNGR Eq A.2
    struct KerrQuantities {
        double rho;       // ρ = √(r² + a²cos²θ)
        double Delta;     // Δ = r² - 2Mr + a²
        double Sigma;     // Σ = √((r² + a²)² - a²Δsin²θ)
        double alpha;     // α = ρ√Δ/Σ (lapse)
        double omega;     // ω = 2ar/Σ² (frame dragging)
        double varpi;     // ϖ = Σsinθ/ρ (circumference factor)
    };
    
    // Compute Kerr metric quantities analytically (reference implementation)
    KerrQuantities computeKerrQuantities(double r, double theta, double M, double a) {
        KerrQuantities kq;
        
        double r2 = r * r;
        double a2 = a * a;
        double cos_theta = std::cos(theta);
        double sin_theta = std::sin(theta);
        double cos2 = cos_theta * cos_theta;
        double sin2 = sin_theta * sin_theta;
        
        // ρ = √(r² + a²cos²θ)
        kq.rho = std::sqrt(r2 + a2 * cos2);
        
        // Δ = r² - 2Mr + a²
        kq.Delta = r2 - 2.0 * M * r + a2;
        
        // Σ² = (r² + a²)² - a²Δsin²θ
        double r2_plus_a2 = r2 + a2;
        double Sigma_sq = r2_plus_a2 * r2_plus_a2 - a2 * kq.Delta * sin2;
        kq.Sigma = std::sqrt(Sigma_sq);
        
        // α = ρ√Δ/Σ
        kq.alpha = kq.rho * std::sqrt(kq.Delta) / kq.Sigma;
        
        // ω = 2ar/Σ²
        kq.omega = 2.0 * a * r / (kq.Sigma * kq.Sigma);
        
        // ϖ = Σsinθ/ρ
        kq.varpi = kq.Sigma * sin_theta / kq.rho;
        
        return kq;
    }
    
    // Compute camera speed for circular geodesic (DNGR Eq A.7)
    double computeCircularGeodesicSpeed(double r, double M, double a, const KerrQuantities& kq) {
        // Ω = √M / (r^(3/2) + a√M)
        double sqrtM = std::sqrt(M);
        double r32 = std::pow(r, 1.5);
        double Omega = sqrtM / (r32 + a * sqrtM);
        
        // β = (ϖ/α)(Ω - ω)
        double beta = (kq.varpi / kq.alpha) * (Omega - kq.omega);
        
        return beta;
    }
    
    void SetUp() override {}
    void TearDown() override {}
};

// =============================================================================
// Kerr Metric Quantity Tests (DNGR Eq A.2)
// =============================================================================

// Test: ρ computation matches analytical formula
TEST_F(FIDOFrameTests, RhoComputationCorrect) {
    double M = 1.0;
    double a = 0.5;  // Moderate spin
    double r = 6.0;  // Outside ISCO
    double theta = kPi / 2.0;  // Equatorial
    
    KerrQuantities kq = computeKerrQuantities(r, theta, M, a);
    
    // At equator (θ = π/2), cos²θ = 0, so ρ = r
    EXPECT_NEAR(kq.rho, r, kEpsilon)
        << "At equator, ρ should equal r";
}

// Test: Δ computation matches analytical formula
TEST_F(FIDOFrameTests, DeltaComputationCorrect) {
    double M = 1.0;
    double a = 0.5;
    double r = 6.0;
    double theta = kPi / 2.0;
    
    KerrQuantities kq = computeKerrQuantities(r, theta, M, a);
    
    // Δ = r² - 2Mr + a² = 36 - 12 + 0.25 = 24.25
    double expected_Delta = r*r - 2.0*M*r + a*a;
    EXPECT_NEAR(kq.Delta, expected_Delta, kEpsilon)
        << "Delta should match analytical formula";
}

// Test: Frame dragging ω is zero for Schwarzschild (a=0)
TEST_F(FIDOFrameTests, FrameDraggingZeroForSchwarzschild) {
    double M = 1.0;
    double a = 0.0;  // Schwarzschild
    double r = 6.0;
    double theta = kPi / 2.0;
    
    KerrQuantities kq = computeKerrQuantities(r, theta, M, a);
    
    EXPECT_NEAR(kq.omega, 0.0, kEpsilon)
        << "Frame dragging ω should be zero for Schwarzschild";
}

// Test: Frame dragging increases with spin
TEST_F(FIDOFrameTests, FrameDraggingIncreasesWithSpin) {
    double M = 1.0;
    double r = 6.0;
    double theta = kPi / 2.0;
    
    double a_low = 0.1;
    double a_high = 0.9;
    
    KerrQuantities kq_low = computeKerrQuantities(r, theta, M, a_low);
    KerrQuantities kq_high = computeKerrQuantities(r, theta, M, a_high);
    
    EXPECT_GT(kq_high.omega, kq_low.omega)
        << "Higher spin should produce more frame dragging";
}

// Test: Lapse function α < 1 outside horizon (gravitational redshift)
TEST_F(FIDOFrameTests, LapseFunctionLessThanOne) {
    double M = 1.0;
    double a = 0.5;
    double r = 6.0;
    double theta = kPi / 2.0;
    
    KerrQuantities kq = computeKerrQuantities(r, theta, M, a);
    
    EXPECT_GT(kq.alpha, 0.0) << "Lapse should be positive outside horizon";
    EXPECT_LT(kq.alpha, 1.0) << "Lapse should be < 1 near black hole";
}

// Test: Lapse function approaches 1 at large r (asymptotic flatness)
TEST_F(FIDOFrameTests, LapseFunctionApproachesOneAtInfinity) {
    double M = 1.0;
    double a = 0.5;
    double theta = kPi / 2.0;
    
    std::vector<double> radii = {100.0, 500.0, 1000.0};
    
    for (double r : radii) {
        KerrQuantities kq = computeKerrQuantities(r, theta, M, a);
        
        // α should approach 1 as r → ∞
        // At r=100M: α ≈ √(1 - 2/100) = √0.98 ≈ 0.9899
        EXPECT_GT(kq.alpha, 0.98) 
            << "Lapse should approach 1 at r = " << r;
    }
}

// =============================================================================
// FIDO Basis Vector Tests (DNGR Eq A.3)
// =============================================================================

// Test: FIDO basis vectors are orthonormal at equator
TEST_F(FIDOFrameTests, FIDOBasisOrthonormalAtEquator) {
    // At equator (θ = π/2), φ = 0:
    // e_r = (1, 0, 0)   [sinθcosφ, cosθ, sinθsinφ] = [1, 0, 0]
    // e_θ = (0, -1, 0)  [cosθcosφ, -sinθ, cosθsinφ] = [0, -1, 0]
    // e_φ = (0, 0, 1)   [-sinφ, 0, cosφ] = [0, 0, 1]
    
    double theta = kPi / 2.0;
    double phi = 0.0;
    
    double sin_theta = std::sin(theta);
    double cos_theta = std::cos(theta);
    double sin_phi = std::sin(phi);
    double cos_phi = std::cos(phi);
    
    // e_r components
    double e_r_x = sin_theta * cos_phi;
    double e_r_y = cos_theta;
    double e_r_z = sin_theta * sin_phi;
    
    // e_θ components
    double e_theta_x = cos_theta * cos_phi;
    double e_theta_y = -sin_theta;
    double e_theta_z = cos_theta * sin_phi;
    
    // e_φ components
    double e_phi_x = -sin_phi;
    double e_phi_y = 0.0;
    double e_phi_z = cos_phi;
    
    // Check orthogonality: e_r · e_θ = 0
    double dot_r_theta = e_r_x * e_theta_x + e_r_y * e_theta_y + e_r_z * e_theta_z;
    EXPECT_NEAR(dot_r_theta, 0.0, kEpsilon) << "e_r and e_θ should be orthogonal";
    
    // Check orthogonality: e_r · e_φ = 0
    double dot_r_phi = e_r_x * e_phi_x + e_r_y * e_phi_y + e_r_z * e_phi_z;
    EXPECT_NEAR(dot_r_phi, 0.0, kEpsilon) << "e_r and e_φ should be orthogonal";
    
    // Check orthogonality: e_θ · e_φ = 0
    double dot_theta_phi = e_theta_x * e_phi_x + e_theta_y * e_phi_y + e_theta_z * e_phi_z;
    EXPECT_NEAR(dot_theta_phi, 0.0, kEpsilon) << "e_θ and e_φ should be orthogonal";
    
    // Check normalization: |e_r| = 1
    double norm_r = std::sqrt(e_r_x*e_r_x + e_r_y*e_r_y + e_r_z*e_r_z);
    EXPECT_NEAR(norm_r, 1.0, kEpsilon) << "e_r should be unit length";
    
    // Check normalization: |e_θ| = 1
    double norm_theta = std::sqrt(e_theta_x*e_theta_x + e_theta_y*e_theta_y + e_theta_z*e_theta_z);
    EXPECT_NEAR(norm_theta, 1.0, kEpsilon) << "e_θ should be unit length";
    
    // Check normalization: |e_φ| = 1
    double norm_phi = std::sqrt(e_phi_x*e_phi_x + e_phi_y*e_phi_y + e_phi_z*e_phi_z);
    EXPECT_NEAR(norm_phi, 1.0, kEpsilon) << "e_φ should be unit length";
}

// =============================================================================
// Camera Speed Tests (DNGR Eq A.7)
// =============================================================================

// Test: Camera speed is subluminal for circular geodesic
TEST_F(FIDOFrameTests, CircularGeodesicSpeedSubluminal) {
    double M = 1.0;
    double a = 0.5;
    double r = 6.0;  // Just outside ISCO for a=0.5
    double theta = kPi / 2.0;
    
    KerrQuantities kq = computeKerrQuantities(r, theta, M, a);
    double beta = computeCircularGeodesicSpeed(r, M, a, kq);
    
    EXPECT_GT(beta, -1.0) << "Speed should be > -1";
    EXPECT_LT(beta, 1.0) << "Speed should be < 1 (subluminal)";
}

// Test: Camera speed increases at smaller radii
TEST_F(FIDOFrameTests, CircularGeodesicSpeedIncreasesNearHorizon) {
    double M = 1.0;
    double a = 0.5;
    double theta = kPi / 2.0;
    
    double r_far = 20.0;
    double r_near = 6.0;
    
    KerrQuantities kq_far = computeKerrQuantities(r_far, theta, M, a);
    KerrQuantities kq_near = computeKerrQuantities(r_near, theta, M, a);
    
    double beta_far = computeCircularGeodesicSpeed(r_far, M, a, kq_far);
    double beta_near = computeCircularGeodesicSpeed(r_near, M, a, kq_near);
    
    EXPECT_GT(std::abs(beta_near), std::abs(beta_far))
        << "Speed should be higher at smaller radii";
}

// Test: Camera speed for Schwarzschild matches Keplerian prediction
TEST_F(FIDOFrameTests, SchwarzschildCircularSpeed) {
    double M = 1.0;
    double a = 0.0;  // Schwarzschild
    double r = 10.0;
    double theta = kPi / 2.0;
    
    KerrQuantities kq = computeKerrQuantities(r, theta, M, a);
    double beta = computeCircularGeodesicSpeed(r, M, a, kq);
    
    // For Schwarzschild, the orbital velocity at infinity normalization is
    // v = sqrt(M/r), but measured locally it's different due to redshift
    // The FIDO-measured speed should be positive (prograde) and reasonable
    EXPECT_GT(beta, 0.0) << "Prograde orbit should have positive speed";
    EXPECT_LT(beta, 0.5) << "Speed at r=10M should be moderate";
}

// =============================================================================
// Schwarzschild Limit Tests
// =============================================================================

// Test: Kerr quantities reduce to Schwarzschild when a=0
TEST_F(FIDOFrameTests, KerrReducesToSchwarzschildAtZeroSpin) {
    double M = 1.0;
    double a = 0.0;
    double r = 6.0;
    double theta = kPi / 2.0;
    
    KerrQuantities kq = computeKerrQuantities(r, theta, M, a);
    
    // For Schwarzschild at equator:
    // ρ = r (since a=0)
    EXPECT_NEAR(kq.rho, r, kEpsilon);
    
    // Δ = r² - 2Mr = r(r - 2M)
    double expected_Delta = r * (r - 2.0 * M);
    EXPECT_NEAR(kq.Delta, expected_Delta, kEpsilon);
    
    // ω = 0 (no frame dragging)
    EXPECT_NEAR(kq.omega, 0.0, kEpsilon);
    
    // α = √(1 - 2M/r) = √(1 - 1/3) ≈ 0.816
    double expected_alpha = std::sqrt(1.0 - 2.0 * M / r);
    EXPECT_NEAR(kq.alpha, expected_alpha, kEpsilon);
}

} // namespace sirius::test
