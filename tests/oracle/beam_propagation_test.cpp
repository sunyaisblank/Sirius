// Beam-propagation oracle: Jacobian determinant, magnification, caustic
// detection, geodesic <-> beam round-trip.
// Ported from TSGD001A.cpp; assertions and tolerances unchanged.

#include "sirius/oracle/beam_integrator.h"
#include "sirius/oracle/kerr_boyer_lindquist.h"
#include "sirius/oracle/metric_interface.h"
#include "sirius/oracle/transport_types.h"

#include <gtest/gtest.h>

#include <cmath>
#include <numbers>

using namespace sirius::oracle;

namespace {

// Fixture for beam propagation tests
class BeamPropagationTest : public ::testing::Test {
  protected:
    void SetUp() override {
        schwarzschild = std::make_unique<KerrMetricD>(1.0, 0.0);
        kerr = std::make_unique<KerrMetricD>(1.0, 0.5);

        BeamIntegratorD::Config config;
        config.step_size = 0.1;
        config.escape_radius = 1000.0;

        integrator_sch = std::make_unique<BeamIntegratorD>(schwarzschild.get(), config);
        integrator_kerr = std::make_unique<BeamIntegratorD>(kerr.get(), config);
    }

    // Create a beam starting at given position with outgoing radial direction
    BeamStateD createOutgoingBeam(const IMetricD* metric, double r, double theta) {
        BeamStateD beam;
        beam.Initialise();

        beam.x = Vec4d(0, r, theta, 0);

        double g[4][4], g_inv[4][4];
        metric->Evaluate(beam.x, g, g_inv);

        // Outgoing null ray: k_t = -1, k_r > 0, k_θ = 0, k_φ = 0.1
        beam.k.t = -1.0;
        beam.k.phi = 0.1;
        beam.k.theta = 0;

        // Solve for k_r from null condition
        double A = g_inv[1][1];
        double C = g_inv[0][0] * beam.k.t * beam.k.t + 2 * g_inv[0][3] * beam.k.t * beam.k.phi +
                   g_inv[3][3] * beam.k.phi * beam.k.phi;

        if (A > 0 && C < 0) {
            beam.k.r = std::sqrt(-C / A);  // Positive = outgoing
        }

        beam.E = -beam.k.t;
        beam.Lz = beam.k.phi;
        beam.initial_pixel_solid_angle = 1e-6;  // 1 arcsec² roughly

        return beam;
    }

    std::unique_ptr<KerrMetricD> schwarzschild;
    std::unique_ptr<KerrMetricD> kerr;
    std::unique_ptr<BeamIntegratorD> integrator_sch;
    std::unique_ptr<BeamIntegratorD> integrator_kerr;
};

void ClearJacobiState(BeamStateD& beam) {
    for (int mu = 0; mu < 4; ++mu) {
        for (int column = 0; column < 4; ++column) {
            beam.J[mu][column] = 0.0;
            beam.dJ[mu][column] = 0.0;
        }
    }
}

// BeamStateD stores V = Dxi/dlambda. Analytic Schwarzschild solutions are
// conventionally written as coordinate xi and dxi/dlambda, so this performs the
// non-optional connection conversion at the initial event.
void SetJacobiColumn(const IMetricD& metric, BeamStateD& beam, int column, const Vec4d& xi,
                     const Vec4d& coordinate_derivative) {
    double g[4][4], g_inv[4][4];
    metric.Evaluate(beam.x, g, g_inv);
    const Vec4d k = metric.dHdp(beam.x, beam.k);

    double Gamma[4][4][4];
    metric.Christoffel(beam.x, Gamma);
    for (int mu = 0; mu < 4; ++mu) {
        beam.J[mu][column] = xi[mu];
        double covariant_derivative = coordinate_derivative[mu];
        for (int nu = 0; nu < 4; ++nu) {
            for (int rho = 0; rho < 4; ++rho) {
                covariant_derivative += Gamma[mu][nu][rho] * k[nu] * xi[rho];
            }
        }
        beam.dJ[mu][column] = covariant_derivative;
    }
}

double MetricInnerProduct(const IMetricD& metric, const Vec4d& x, const Vec4d& a, const Vec4d& b) {
    double g[4][4], g_inv[4][4];
    metric.Evaluate(x, g, g_inv);
    double result = 0.0;
    for (int mu = 0; mu < 4; ++mu) {
        for (int nu = 0; nu < 4; ++nu) {
            result += g[mu][nu] * a[mu] * b[nu];
        }
    }
    return result;
}

Vec4d JacobiColumn(const BeamStateD& beam, int column) {
    return Vec4d(beam.J[0][column], beam.J[1][column], beam.J[2][column], beam.J[3][column]);
}

//==============================================================================
// Test: BeamStateD Initialisation
// Verifies: Jacobian starts as Identity matrix
//==============================================================================

TEST_F(BeamPropagationTest, BeamInitialisation) {
    BeamStateD beam;
    beam.Initialise();

    // Jacobian should be Identity
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
    EXPECT_FALSE(beam.at_caustic);
    EXPECT_DOUBLE_EQ(beam.magnification, 1.0);
}

//==============================================================================
// Test: Jacobian Determinant in Flat Space
// Verifies: det(J) = 1 for all time (Liouville's theorem)
//==============================================================================

TEST_F(BeamPropagationTest, JacobianDeterminantFlat) {
    // For very large r, Schwarzschild approaches flat space
    BeamStateD beam = createOutgoingBeam(schwarzschild.get(), 1000.0, std::numbers::pi / 2);

    // Initial Determinant should be 1 (Identity matrix)
    double det0 = beam.J[0][0] *
                      (beam.J[1][1] * (beam.J[2][2] * beam.J[3][3] - beam.J[2][3] * beam.J[3][2]) -
                       beam.J[1][2] * (beam.J[2][1] * beam.J[3][3] - beam.J[2][3] * beam.J[3][1]) +
                       beam.J[1][3] * (beam.J[2][1] * beam.J[3][2] - beam.J[2][2] * beam.J[3][1])) +
                  0;  // Simplified

    EXPECT_NEAR(det0, 1.0, 1e-10) << "Initial Jacobian not Identity";

    // Integrate for several steps
    const int numSteps = 100;
    for (int i = 0; i < numSteps; ++i) {
        if (!integrator_sch->Step(beam, 0.1)) break;
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
    beam.Initialise();
    beam.initial_pixel_solid_angle = 1e-6;

    // Set a known Jacobian with specific eigenvalues
    // 2×2 angular block with singular values 2 and 0.5
    beam.J[2][2] = 2.0;
    beam.J[2][3] = 0.0;
    beam.J[3][2] = 0.0;
    beam.J[3][3] = 0.5;

    beam.UpdateGeometry();

    // Major axis = 2, minor axis = 0.5
    EXPECT_NEAR(beam.major_axis, 2.0, 0.01);
    EXPECT_NEAR(beam.minor_axis, 0.5, 0.01);

    // Magnification = 1/det = 1/(2×0.5) = 1
    EXPECT_NEAR(beam.magnification, 1.0, 0.01);

    // Solid angle = π × 2 × 0.5 × 1e-6 = π×1e-6
    EXPECT_NEAR(beam.solid_angle, std::numbers::pi * 1e-6, 1e-8);
}

TEST_F(BeamPropagationTest, OrientationDescribesOutputEllipseRatherThanInputBasis) {
    BeamStateD beam;
    beam.Initialise();

    // J = R(0.4) diag(3, 1) R(-0.2). The right rotation changes only the input
    // basis; the rendered output ellipse must remain oriented at +0.4 radians.
    constexpr double output_angle = 0.4;
    constexpr double input_angle = -0.2;
    const double co = std::cos(output_angle);
    const double so = std::sin(output_angle);
    const double ci = std::cos(input_angle);
    const double si = std::sin(input_angle);
    beam.J[2][2] = 3.0 * co * ci - so * si;
    beam.J[2][3] = -3.0 * co * si - so * ci;
    beam.J[3][2] = 3.0 * so * ci + co * si;
    beam.J[3][3] = -3.0 * so * si + co * ci;

    beam.UpdateGeometry();

    EXPECT_NEAR(beam.major_axis, 3.0, 1.0e-12);
    EXPECT_NEAR(beam.minor_axis, 1.0, 1.0e-12);
    EXPECT_NEAR(beam.orientation, output_angle, 1.0e-12);
}

//==============================================================================
// Test: Caustic Detection
// Verifies: at_caustic flag set when det(J_angular) → 0
//==============================================================================

TEST_F(BeamPropagationTest, CausticDetection) {
    BeamStateD beam;
    beam.Initialise();

    // Set nearly singular angular Jacobian
    beam.J[2][2] = 1.0;
    beam.J[2][3] = 1.0;
    beam.J[3][2] = 1.0;
    beam.J[3][3] = 1.0 + 1e-14;  // det ≈ 1e-14

    beam.UpdateGeometry();

    EXPECT_TRUE(beam.at_caustic) << "Should detect caustic when det ≈ 0";
    EXPECT_GT(beam.magnification, 1e10) << "Magnification should be very high at caustic";
}

//==============================================================================
// Test: Beam Integration Step
// Verifies: Single Step updates position and Jacobian
//==============================================================================

TEST_F(BeamPropagationTest, BeamIntegrationStep) {
    BeamStateD beam = createOutgoingBeam(schwarzschild.get(), 10.0, std::numbers::pi / 2);
    double initial_r = beam.x.r;
    double initial_lambda = beam.lambda;

    bool success = integrator_sch->Step(beam, 0.1);

    EXPECT_TRUE(success) << "Step should succeed";
    EXPECT_GT(beam.lambda, initial_lambda) << "Lambda should advance";

    // For outgoing ray, r should increase
    EXPECT_GT(beam.x.r, initial_r) << "Outgoing ray should move outward";
}

// Closed-form null-deviation solutions are from Morales-Ruiz and Raposo,
// arXiv:2308.07098, equations (19)-(24). For an outgoing radial ray,
// r(lambda)=r0+E lambda. Choosing point-source screen data xi(0)=0 and unit
// physical derivative gives r*xi_theta = r*sin(theta)*xi_phi = lambda exactly.
TEST_F(BeamPropagationTest, SchwarzschildRadialCongruenceMatchesClosedFormToOnePartPerMillion) {
    constexpr double mass = 1.0;
    constexpr double energy = 1.0;
    constexpr double initial_radius = 10.0;
    constexpr double final_lambda = 2.0;
    constexpr double step = 1.0e-3;

    BeamStateD beam;
    beam.Initialise();
    beam.x = Vec4d(0.0, initial_radius, std::numbers::pi / 2.0, 0.0);
    beam.k = Vec4d(-energy, energy * initial_radius / (initial_radius - 2.0 * mass), 0.0, 0.0);
    beam.E = energy;
    ClearJacobiState(beam);

    SetJacobiColumn(*schwarzschild, beam, 0, Vec4d(), Vec4d(0.0, 0.0, 1.0 / initial_radius, 0.0));
    SetJacobiColumn(*schwarzschild, beam, 1, Vec4d(), Vec4d(0.0, 0.0, 0.0, 1.0 / initial_radius));

    const int steps = static_cast<int>(final_lambda / step);
    for (int i = 0; i < steps; ++i) {
        ASSERT_TRUE(integrator_sch->Step(beam, step));
    }

    const Vec4d polar = JacobiColumn(beam, 0);
    const Vec4d azimuthal = JacobiColumn(beam, 1);
    const double polar_axis =
        std::sqrt(std::abs(MetricInnerProduct(*schwarzschild, beam.x, polar, polar)));
    const double azimuthal_axis =
        std::sqrt(std::abs(MetricInnerProduct(*schwarzschild, beam.x, azimuthal, azimuthal)));
    const double cross = MetricInnerProduct(*schwarzschild, beam.x, polar, azimuthal);

    EXPECT_NEAR(beam.x.r, initial_radius + energy * final_lambda, 1.0e-10);
    EXPECT_NEAR(polar_axis / final_lambda, 1.0, 1.0e-6);
    EXPECT_NEAR(azimuthal_axis / final_lambda, 1.0, 1.0e-6);
    EXPECT_NEAR(cross, 0.0, 1.0e-12);
}

// Equations (16)-(18) of the same independent source give the photon-sphere
// orbit and its constant-coefficient variational system. Equal physical screen
// axes at lambda=0 evolve as cosh(E lambda/(sqrt(3) M)) in the radial direction
// and cos(E lambda/(sqrt(3) M)) normal to the orbital plane.
TEST_F(BeamPropagationTest,
       SchwarzschildCircularPhotonCongruenceMatchesClosedFormToOnePartPerMillion) {
    constexpr double mass = 1.0;
    constexpr double energy = 1.0;
    constexpr double initial_axis = 1.0e-3;
    constexpr double final_lambda = 1.0;
    constexpr double step = 5.0e-4;
    constexpr double sqrt_three = 1.7320508075688772935;

    BeamStateD beam;
    beam.Initialise();
    beam.x = Vec4d(0.0, 3.0 * mass, std::numbers::pi / 2.0, 0.0);
    beam.k = Vec4d(-energy, 0.0, 0.0, 3.0 * sqrt_three * mass * energy);
    beam.E = energy;
    beam.Lz = beam.k.phi;
    ClearJacobiState(beam);

    const double radial_coordinate_axis = initial_axis / sqrt_three;
    const Vec4d radial_xi(0.0, radial_coordinate_axis, 0.0, 0.0);
    const Vec4d radial_coordinate_derivative(
        -2.0 * energy * radial_coordinate_axis / mass, 0.0, 0.0,
        -2.0 * energy * radial_coordinate_axis / (3.0 * sqrt_three * mass * mass));
    SetJacobiColumn(*schwarzschild, beam, 0, radial_xi, radial_coordinate_derivative);

    const Vec4d polar_xi(0.0, 0.0, initial_axis / (3.0 * mass), 0.0);
    SetJacobiColumn(*schwarzschild, beam, 1, polar_xi, Vec4d());

    const int steps = static_cast<int>(final_lambda / step);
    for (int i = 0; i < steps; ++i) {
        ASSERT_TRUE(integrator_sch->Step(beam, step));
    }

    const Vec4d radial = JacobiColumn(beam, 0);
    const Vec4d polar = JacobiColumn(beam, 1);
    const double radial_axis =
        std::sqrt(std::abs(MetricInnerProduct(*schwarzschild, beam.x, radial, radial)));
    const double polar_axis =
        std::sqrt(std::abs(MetricInnerProduct(*schwarzschild, beam.x, polar, polar)));
    const double omega_lambda = energy * final_lambda / (sqrt_three * mass);
    const double expected_radial = initial_axis * std::cosh(omega_lambda);
    const double expected_polar = initial_axis * std::abs(std::cos(omega_lambda));
    const double cross = MetricInnerProduct(*schwarzschild, beam.x, radial, polar);

    EXPECT_NEAR(beam.x.r, 3.0 * mass, 1.0e-10);
    EXPECT_NEAR(radial_axis / expected_radial, 1.0, 1.0e-6);
    EXPECT_NEAR(polar_axis / expected_polar, 1.0, 1.0e-6);
    EXPECT_NEAR(cross, 0.0, 1.0e-12);
}

//==============================================================================
// Test: Horizon Termination
// Verifies: Beam terminates when approaching horizon
//==============================================================================

TEST_F(BeamPropagationTest, HorizonTermination) {
    // Create ingoing beam near horizon
    BeamStateD beam;
    beam.Initialise();
    beam.x = Vec4d(0, 3.0, std::numbers::pi / 2, 0);

    double g[4][4], g_inv[4][4];
    schwarzschild->Evaluate(beam.x, g, g_inv);

    beam.k.t = -1.0;
    beam.k.phi = 0;
    beam.k.theta = 0;

    // Ingoing: negative k_r
    double A = g_inv[1][1];
    double C = g_inv[0][0] * beam.k.t * beam.k.t;
    if (A > 0 && C < 0) {
        beam.k.r = -std::sqrt(-C / A);
    }

    // Integrate until termination
    int steps = 0;
    while (integrator_sch->Step(beam, 0.1) && steps < 1000) {
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
    geo.x = Vec4d(1.0, 10.0, std::numbers::pi / 3, std::numbers::pi / 4);
    geo.k = Vec4d(-1.0, 0.5, 0.1, 2.0);
    geo.lambda = 42.0;
    geo.E = 1.0;
    geo.Lz = 2.0;
    geo.Q = 0.5;

    BeamStateD beam;
    beam.Initialise();
    beam.FromGeodesic(geo);

    GeodesicStateD geo2 = beam.ToGeodesic();

    EXPECT_DOUBLE_EQ(geo2.x.t, geo.x.t);
    EXPECT_DOUBLE_EQ(geo2.x.r, geo.x.r);
    EXPECT_DOUBLE_EQ(geo2.x.theta, geo.x.theta);
    EXPECT_DOUBLE_EQ(geo2.x.phi, geo.x.phi);
    EXPECT_DOUBLE_EQ(geo2.lambda, geo.lambda);
    EXPECT_DOUBLE_EQ(geo2.E, geo.E);
    EXPECT_DOUBLE_EQ(geo2.Lz, geo.Lz);
}

}  // namespace
