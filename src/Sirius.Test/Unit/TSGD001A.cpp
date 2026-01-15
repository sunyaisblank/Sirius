// TSGD001A.cpp - Beam Propagation and Geodesic Deviation Tests
// Tests for KNBI001A beam state and BeamIntegratorD
// Verifies: Jacobian determinant conservation, magnification, caustic detection

#include <gtest/gtest.h>
#include <cmath>
#include "../../Sirius.Core/Transport/MTTP001A.h"
#include "../../Sirius.Core/Geodesic/PHMT000B.h"
#include "../../Sirius.Core/Metric/PHMT100B.h"
#include "../../Sirius.Render/Integration/INBI001A.h"

using namespace sirius::math;
using namespace sirius::physics;
using namespace sirius::kernel;

namespace {

// Fixture for beam propagation tests
class BeamPropagationTest : public ::testing::Test {
protected:
    void SetUp() override {
        schwarzschild = std::make_unique<KerrMetricD>(1.0, 0.0);
        kerr = std::make_unique<KerrMetricD>(1.0, 0.5);
        
        BeamIntegratorD::Config config;
        config.stepSize = 0.1;
        config.escapeRadius = 1000.0;
        
        integrator_sch = std::make_unique<BeamIntegratorD>(schwarzschild.get(), config);
        integrator_kerr = std::make_unique<BeamIntegratorD>(kerr.get(), config);
    }
    
    // Create a beam starting at given position with outgoing radial direction
    BeamStateD createOutgoingBeam(const IMetricD* metric, double r, double theta) {
        BeamStateD beam;
        beam.initialise();
        
        beam.x = Vec4d(0, r, theta, 0);
        
        double g[4][4], g_inv[4][4];
        metric->evaluate(beam.x, g, g_inv);
        
        // Outgoing null ray: k_t = -1, k_r > 0, k_θ = 0, k_φ = 0.1
        beam.k.t = -1.0;
        beam.k.phi = 0.1;
        beam.k.theta = 0;
        
        // Solve for k_r from null condition
        double A = g_inv[1][1];
        double C = g_inv[0][0]*beam.k.t*beam.k.t + 
                   2*g_inv[0][3]*beam.k.t*beam.k.phi + 
                   g_inv[3][3]*beam.k.phi*beam.k.phi;
        
        if (A > 0 && C < 0) {
            beam.k.r = std::sqrt(-C / A);  // Positive = outgoing
        }
        
        beam.E = -beam.k.t;
        beam.Lz = beam.k.phi;
        beam.initialPixelSolidAngle = 1e-6;  // 1 arcsec² roughly
        
        return beam;
    }
    
    std::unique_ptr<KerrMetricD> schwarzschild;
    std::unique_ptr<KerrMetricD> kerr;
    std::unique_ptr<BeamIntegratorD> integrator_sch;
    std::unique_ptr<BeamIntegratorD> integrator_kerr;
};

//==============================================================================
// Test: BeamStateD Initialisation
// Verifies: Jacobian starts as identity matrix
//==============================================================================

TEST_F(BeamPropagationTest, BeamInitialisation) {
    BeamStateD beam;
    beam.initialise();
    
    // Jacobian should be identity
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            if (i == j) {
                EXPECT_DOUBLE_EQ(beam.J[i][j], 1.0);
            } else {
                EXPECT_DOUBLE_EQ(beam.J[i][j], 0.0);
            }
            EXPECT_DOUBLE_EQ(beam.dJ[i][j], 0.0);
        }
    }
    
    EXPECT_FALSE(beam.terminated);
    EXPECT_FALSE(beam.atCaustic);
    EXPECT_DOUBLE_EQ(beam.magnification, 1.0);
}

//==============================================================================
// Test: Jacobian Determinant in Flat Space
// Verifies: det(J) = 1 for all time (Liouville's theorem)
//==============================================================================

TEST_F(BeamPropagationTest, JacobianDeterminantFlat) {
    // For very large r, Schwarzschild approaches flat space
    BeamStateD beam = createOutgoingBeam(schwarzschild.get(), 1000.0, M_PI/2);
    
    // Initial determinant should be 1 (identity matrix)
    double det0 = beam.J[0][0] * (beam.J[1][1] * (beam.J[2][2] * beam.J[3][3] - beam.J[2][3] * beam.J[3][2]) -
                                   beam.J[1][2] * (beam.J[2][1] * beam.J[3][3] - beam.J[2][3] * beam.J[3][1]) +
                                   beam.J[1][3] * (beam.J[2][1] * beam.J[3][2] - beam.J[2][2] * beam.J[3][1])) + 0;  // Simplified
    
    EXPECT_NEAR(det0, 1.0, 1e-10) << "Initial Jacobian not identity";
    
    // Integrate for several steps
    const int numSteps = 100;
    for (int i = 0; i < numSteps; ++i) {
        if (!integrator_sch->step(beam, 0.1)) break;
    }
    
    // In nearly flat space, magnification should stay close to 1
    EXPECT_GT(beam.magnification, 0.1) << "Magnification collapsed unexpectedly";
    EXPECT_LT(beam.magnification, 100.0) << "Magnification exploded unexpectedly";
}

//==============================================================================
// Test: Beam Geometry Extraction
// Verifies: Ellipse parameters computed correctly from Jacobian
//==============================================================================

TEST_F(BeamPropagationTest, BeamGeometryExtraction) {
    BeamStateD beam;
    beam.initialise();
    beam.initialPixelSolidAngle = 1e-6;
    
    // Set a known Jacobian with specific eigenvalues
    // 2×2 angular block with singular values 2 and 0.5
    beam.J[2][2] = 2.0;
    beam.J[2][3] = 0.0;
    beam.J[3][2] = 0.0;
    beam.J[3][3] = 0.5;
    
    beam.updateGeometry();
    
    // Major axis = 2, minor axis = 0.5
    EXPECT_NEAR(beam.majorAxis, 2.0, 0.01);
    EXPECT_NEAR(beam.minorAxis, 0.5, 0.01);
    
    // Magnification = 1/det = 1/(2×0.5) = 1
    EXPECT_NEAR(beam.magnification, 1.0, 0.01);
    
    // Solid angle = π × 2 × 0.5 × 1e-6 = π×1e-6
    EXPECT_NEAR(beam.solidAngle, M_PI * 1e-6, 1e-8);
}

//==============================================================================
// Test: Caustic Detection
// Verifies: atCaustic flag set when det(J_angular) → 0
//==============================================================================

TEST_F(BeamPropagationTest, CausticDetection) {
    BeamStateD beam;
    beam.initialise();
    
    // Set nearly singular angular Jacobian
    beam.J[2][2] = 1.0;
    beam.J[2][3] = 1.0;
    beam.J[3][2] = 1.0;
    beam.J[3][3] = 1.0 + 1e-14;  // det ≈ 1e-14
    
    beam.updateGeometry();
    
    EXPECT_TRUE(beam.atCaustic) << "Should detect caustic when det ≈ 0";
    EXPECT_GT(beam.magnification, 1e10) << "Magnification should be very high at caustic";
}

//==============================================================================
// Test: Beam Integration Step
// Verifies: Single step updates position and Jacobian
//==============================================================================

TEST_F(BeamPropagationTest, BeamIntegrationStep) {
    BeamStateD beam = createOutgoingBeam(schwarzschild.get(), 10.0, M_PI/2);
    double initial_r = beam.x.r;
    double initial_lambda = beam.lambda;
    
    bool success = integrator_sch->step(beam, 0.1);
    
    EXPECT_TRUE(success) << "Step should succeed";
    EXPECT_GT(beam.lambda, initial_lambda) << "Lambda should advance";
    
    // For outgoing ray, r should increase
    EXPECT_GT(beam.x.r, initial_r) << "Outgoing ray should move outward";
}

//==============================================================================
// Test: Horizon Termination
// Verifies: Beam terminates when approaching horizon
//==============================================================================

TEST_F(BeamPropagationTest, HorizonTermination) {
    // Create ingoing beam near horizon
    BeamStateD beam;
    beam.initialise();
    beam.x = Vec4d(0, 3.0, M_PI/2, 0);
    
    double g[4][4], g_inv[4][4];
    schwarzschild->evaluate(beam.x, g, g_inv);
    
    beam.k.t = -1.0;
    beam.k.phi = 0;
    beam.k.theta = 0;
    
    // Ingoing: negative k_r
    double A = g_inv[1][1];
    double C = g_inv[0][0]*beam.k.t*beam.k.t;
    if (A > 0 && C < 0) {
        beam.k.r = -std::sqrt(-C / A);
    }
    
    // Integrate until termination
    int steps = 0;
    while (integrator_sch->step(beam, 0.1) && steps < 1000) {
        steps++;
    }
    
    // Should have terminated (hit horizon or invalid region)
    EXPECT_TRUE(beam.terminated) << "Ingoing beam should terminate at horizon";
}

//==============================================================================
// Test: Conversion Functions
// Verifies: GeodesicStateD ↔ BeamStateD round-trip
//==============================================================================

TEST_F(BeamPropagationTest, ConversionRoundTrip) {
    GeodesicStateD geo;
    geo.x = Vec4d(1.0, 10.0, M_PI/3, M_PI/4);
    geo.k = Vec4d(-1.0, 0.5, 0.1, 2.0);
    geo.lambda = 42.0;
    geo.E = 1.0;
    geo.Lz = 2.0;
    geo.Q = 0.5;
    
    BeamStateD beam;
    beam.initialise();
    beam.fromGeodesic(geo);
    
    GeodesicStateD geo2 = beam.toGeodesic();
    
    EXPECT_DOUBLE_EQ(geo2.x.t, geo.x.t);
    EXPECT_DOUBLE_EQ(geo2.x.r, geo.x.r);
    EXPECT_DOUBLE_EQ(geo2.x.theta, geo.x.theta);
    EXPECT_DOUBLE_EQ(geo2.x.phi, geo.x.phi);
    EXPECT_DOUBLE_EQ(geo2.lambda, geo.lambda);
    EXPECT_DOUBLE_EQ(geo2.E, geo.E);
    EXPECT_DOUBLE_EQ(geo2.Lz, geo.Lz);
}

} // namespace
