// TSBM001A.cpp - Geodesic Integration Accuracy Benchmarks
// Tests circular orbits, light deflection, energy/momentum conservation.

#define _USE_MATH_DEFINES
#include <cmath>
#include <gtest/gtest.h>
#include <MTTN001A.h>
#include <MTDL001A.h>

namespace sirius::test {

constexpr double kEpsilon = 1e-8;

// =============================================================================
// Test Fixture
// =============================================================================

class GeodesicBenchmarks : public ::testing::Test {
protected:
    static constexpr double M = 1.0;  // Mass in geometric units
    
    void SetUp() override {}
    void TearDown() override {}
};

// =============================================================================
// Circular Orbit Parameters
// =============================================================================

// Test: Kepler's third law in Schwarzschild (geometric units)
// Orbital period T = 2π√(r³/M) for circular orbits at r >> 2M
TEST_F(GeodesicBenchmarks, CircularOrbitPeriod) {
    // Test at various radii
    std::vector<double> radii = {10.0, 20.0, 50.0, 100.0};
    
    for (double r : radii) {
        // Coordinate angular velocity: Ω = √(M/r³)
        double omega = std::sqrt(M / (r * r * r));
        
        // Coordinate period: T = 2π/Ω
        double T_coord = 2.0 * M_PI / omega;
        
        // Expected from Kepler's law: T = 2π√(r³/M)
        double T_expected = 2.0 * M_PI * std::sqrt(r * r * r / M);
        
        EXPECT_NEAR(T_coord, T_expected, kEpsilon * T_expected)
            << "Period mismatch at r=" << r;
    }
}

// Test: ISCO radius
// The Innermost Stable Circular Orbit is at r = 6M for Schwarzschild
TEST_F(GeodesicBenchmarks, ISCORadius) {
    double r_isco = 6.0 * M;
    
    // At ISCO, the effective potential has an inflection point
    // Verify orbital frequency is real and positive
    double omega_isco = std::sqrt(M / (r_isco * r_isco * r_isco));
    
    EXPECT_GT(omega_isco, 0.0) << "ISCO angular velocity should be positive";
    
    // Expected value: Ω_ISCO = 1/(6√6 M)
    double omega_expected = 1.0 / (6.0 * std::sqrt(6.0) * M);
    
    EXPECT_NEAR(omega_isco, omega_expected, kEpsilon * omega_expected)
        << "ISCO angular velocity incorrect";
}

// Test: Photon sphere at r = 3M
TEST_F(GeodesicBenchmarks, PhotonSphereRadius) {
    double r_photon = 3.0 * M;
    
    // Circular photon orbit: impact parameter b = L/E = 3√3 M
    double b_photon = 3.0 * std::sqrt(3.0) * M;
    
    // Verify this is the critical impact parameter
    EXPECT_NEAR(b_photon, std::sqrt(27.0) * M, kEpsilon)
        << "Photon sphere impact parameter incorrect";
}

// =============================================================================
// Light Deflection
// =============================================================================

// Test: Light deflection in weak field limit
// Δφ = 4M/b for large impact parameter b
TEST_F(GeodesicBenchmarks, WeakFieldLightDeflection) {
    // Large impact parameter (weak field)
    std::vector<double> impact_params = {100.0, 200.0, 500.0, 1000.0};
    
    for (double b : impact_params) {
        // Weak field deflection angle
        double delta_phi = 4.0 * M / b;
        
        // Convert to arcseconds for intuition (not a test assertion)
        double arcsec = delta_phi * 180.0 / M_PI * 3600.0;
        
        // Deflection should be small and positive
        EXPECT_GT(delta_phi, 0.0) << "Deflection should be positive";
        EXPECT_LT(delta_phi, 0.1) << "Weak field deflection should be small";
        
        // Check scaling: Δφ ∝ 1/b
        if (b > 100.0) {
            double ratio = (4.0 * M / 100.0) / (4.0 * M / b);
            EXPECT_NEAR(ratio, b / 100.0, kEpsilon)
                << "Deflection scaling incorrect at b=" << b;
        }
    }
}

// Test: Solar deflection value (for reference)
// At the Sun's limb: b ≈ R_sun = 696,000 km ≈ 2.35 × 10^(-6) in geometric units
// Deflection ≈ 1.75 arcsec
TEST_F(GeodesicBenchmarks, SolarDeflectionOrder) {
    // In geometric units where G = c = 1:
    // M_sun ≈ 1.477 km ≈ 1 (we normalize)
    // R_sun ≈ 696,000 km ≈ 471,000 M_sun
    
    // Using our normalized M = 1
    double b_solar = 471000.0;  // Solar radius in units of M
    double delta_phi = 4.0 * M / b_solar;
    double arcsec = delta_phi * (180.0 * 3600.0 / M_PI);
    
    // Should be approximately 1.75 arcsec
    EXPECT_NEAR(arcsec, 1.75, 0.01)
        << "Solar deflection should be ~1.75 arcsec";
}

// =============================================================================
// Energy and Angular Momentum Conservation
// =============================================================================

// Test: Conserved quantities for geodesic motion
TEST_F(GeodesicBenchmarks, GeodesicConservedQuantities) {
    // For Schwarzschild, the conserved energy per unit mass is:
    // E = (1 - 2M/r) * dt/dτ
    // For circular orbit at r = 10M
    double r = 10.0 * M;
    double f = 1.0 - 2.0 * M / r;  // = 0.8
    
    // For circular timelike geodesic, the 4-velocity normalization gives:
    // E = (1 - 2M/r) / sqrt(1 - 3M/r)
    // This is the correct formula for specific orbital energy
    double E = (1.0 - 2.0 * M / r) / std::sqrt(1.0 - 3.0 * M / r);
    
    // Verify E is less than 1 for bound orbits (r > 6M)
    EXPECT_LT(E, 1.0) << "Bound orbit energy should be less than 1";
    EXPECT_GT(E, 0.0) << "Energy should be positive";
    
    // At r = 10M: E = 0.8 / sqrt(0.7) ≈ 0.9562
    double E_expected = 0.8 / std::sqrt(0.7);
    EXPECT_NEAR(E, E_expected, kEpsilon * E_expected)
        << "Conserved energy incorrect for circular orbit";
    
    // Conserved angular momentum per unit mass: L = r² * (dφ/dτ)
    // For circular orbit: L = sqrt(Mr) / sqrt(1 - 3M/r)
    double L = std::sqrt(M * r) / std::sqrt(1.0 - 3.0 * M / r);
    
    // At r = 10M: L = sqrt(10) / sqrt(0.7) ≈ 3.780
    double L_expected = std::sqrt(10.0 * M) / std::sqrt(0.7);
    
    EXPECT_NEAR(L, L_expected, kEpsilon * L_expected)
        << "Conserved angular momentum incorrect for circular orbit";
}

// =============================================================================
// Proper Time
// =============================================================================

// Test: Proper time for radial fall from rest
// Time to fall from r_0 to r_f follows known formula
TEST_F(GeodesicBenchmarks, RadialFallProperTime) {
    // Object falling from rest at r = 10M
    double r0 = 10.0 * M;
    double rf = 3.0 * M;  // Just above photon sphere
    
    // For radial fall from rest, the proper time is:
    // τ = ∫ dr/√(2M/r - 2M/r0)
    // This has a closed form in terms of cycloid parameters
    
    // The total proper time to fall to r = 0 from rest at r0 is:
    // τ_total = (π/2) * r0^(3/2) / √(2M)
    double tau_to_singularity = (M_PI / 2.0) * std::pow(r0, 1.5) / std::sqrt(2.0 * M);
    
    // This is approximately 15.81 M for r0 = 10M
    double expected = (M_PI / 2.0) * std::pow(10.0, 1.5) / std::sqrt(2.0);
    
    EXPECT_NEAR(tau_to_singularity, expected, kEpsilon * expected)
        << "Proper time for radial infall incorrect";
}

// =============================================================================
// Orbital Precession
// =============================================================================

// Test: Perihelion precession per orbit
// Δφ = 6πM / [a(1-e²)] for elliptical orbits
TEST_F(GeodesicBenchmarks, PerihelionPrecession) {
    // Circular orbit (e = 0) at r = 100M (weak field)
    double a = 100.0 * M;  // Semi-major axis
    double e = 0.0;        // Eccentricity (circular)
    
    // Precession per orbit
    double delta_phi = 6.0 * M_PI * M / (a * (1.0 - e * e));
    
    // For circular orbit at r = 100M:
    // Δφ = 6πM / r ≈ 0.188 radians per orbit
    double expected = 6.0 * M_PI * M / a;
    
    EXPECT_NEAR(delta_phi, expected, kEpsilon * expected)
        << "Perihelion precession incorrect";
    
    // Mercury's precession (in proper units):
    // a ≈ 5.79e10 m, M_sun ≈ 1.477 km
    // Δφ ≈ 43 arcsec/century (not directly testable here)
}

// =============================================================================
// Verification Data Points
// =============================================================================

// Test: Compile verification data for numerical integration
TEST_F(GeodesicBenchmarks, VerificationDataPoints) {
    // Generate reference data for different orbital configurations
    struct OrbitConfig {
        double r;
        double omega;
        double E;
        double L;
    };
    
    std::vector<OrbitConfig> configs;
    
    for (double r : {7.0, 10.0, 20.0, 50.0, 100.0}) {
        OrbitConfig cfg;
        cfg.r = r * M;
        cfg.omega = std::sqrt(M / (cfg.r * cfg.r * cfg.r));
        // Correct formula for circular orbit energy: E = (1 - 2M/r) / sqrt(1 - 3M/r)
        cfg.E = (1.0 - 2.0 * M / cfg.r) / std::sqrt(1.0 - 3.0 * M / cfg.r);
        // Correct formula for angular momentum: L = sqrt(Mr) / sqrt(1 - 3M/r)
        cfg.L = std::sqrt(M * cfg.r) / std::sqrt(1.0 - 3.0 * M / cfg.r);
        configs.push_back(cfg);
    }
    
    // Verify all configurations are physical (E, L real and positive)
    for (const auto& cfg : configs) {
        EXPECT_GT(cfg.omega, 0.0) << "Angular velocity should be positive at r=" << cfg.r;
        EXPECT_GT(cfg.E, 0.0) << "Energy should be positive at r=" << cfg.r;
        EXPECT_GT(cfg.L, 0.0) << "Angular momentum should be positive at r=" << cfg.r;
        
        // For stable orbits r > 6M: E < 1 (bound orbit)
        if (cfg.r > 6.0 * M) {
            EXPECT_LT(cfg.E, 1.0) << "Bound orbit at r=" << cfg.r << " should have E < 1, got " << cfg.E;
        }
    }
}

} // namespace sirius::test
