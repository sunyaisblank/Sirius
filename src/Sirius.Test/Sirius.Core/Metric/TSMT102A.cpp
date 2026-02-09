// TSMT102A.cpp - Alcubierre Warp Drive Metric Property Tests
// Component ID: TSMT102A (Test/Metric/AlcubierreProperties)
// Tests: PHMT102A.h (WarpDriveFamily)

#include <gtest/gtest.h>
#include <cmath>
#include <MTTN001A.h>
#include <MTDL001A.h>
#include <PHMT102A.h>

namespace sirius::test {

constexpr double kEps = 1e-8;

// =============================================================================
// Test Fixture
// =============================================================================

class AlcubierreMetricTests : public ::testing::Test {
protected:
    Sirius::WarpDriveFamily metric;
    Sirius::WarpDriveParams defaultParams;

    void SetUp() override {
        defaultParams = Sirius::WarpDriveParams::Alcubierre(1.0, 1.0);
        metric.setParams(defaultParams);
    }

    void evaluateAt(double x, double y, double z, Metric4D& g) {
        Tensor<double, 4> pos;
        pos(0) = 0.0; pos(1) = x; pos(2) = y; pos(3) = z;
        Tensor<Dual<double>, 4, 4, 4> dg;
        metric.evaluate(pos, g, dg);
    }
};

// =============================================================================
// Metric Symmetry Tests
// =============================================================================

TEST_F(AlcubierreMetricTests, MetricSymmetry) {
    // Test symmetry g_mu_nu = g_nu_mu at multiple positions
    double positions[][3] = {{0.0, 0.0, 0.0}, {0.5, 0.0, 0.0}, {2.0, 1.0, 0.5}, {10.0, 0.0, 0.0}};
    for (auto& p : positions) {
        Metric4D g;
        evaluateAt(p[0], p[1], p[2], g);
        for (int mu = 0; mu < 4; ++mu) {
            for (int nu = mu + 1; nu < 4; ++nu) {
                EXPECT_NEAR(g(mu, nu).real, g(nu, mu).real, kEps)
                    << "Asymmetry at (" << p[0] << "," << p[1] << "," << p[2]
                    << ") indices (" << mu << "," << nu << ")";
            }
        }
    }
}

TEST_F(AlcubierreMetricTests, LorentzianSignature) {
    // At any point, the metric should have Lorentzian signature (-,+,+,+)
    double positions[][3] = {{0.0, 0.0, 0.0}, {0.5, 0.0, 0.0}, {5.0, 0.0, 0.0}};
    for (auto& p : positions) {
        Metric4D g;
        evaluateAt(p[0], p[1], p[2], g);

        // For the Alcubierre metric in these coordinates:
        // g_yy = 1, g_zz = 1 (positive)
        EXPECT_NEAR(g(2, 2).real, 1.0, kEps);
        EXPECT_NEAR(g(3, 3).real, 1.0, kEps);

        // g_xx = 1 (positive)
        EXPECT_NEAR(g(1, 1).real, 1.0, kEps);

        // det(2x2 block of g_tt, g_tx) should be negative for Lorentzian signature
        // det = g_tt * g_xx - g_tx^2 = -(1 - vs^2 f^2) * 1 - (-vs*f)^2
        //     = -(1 - vs^2 f^2) - vs^2 f^2 = -1
        double det2 = g(0, 0).real * g(1, 1).real - g(0, 1).real * g(1, 0).real;
        EXPECT_NEAR(det2, -1.0, kEps)
            << "2x2 determinant should be -1 at (" << p[0] << "," << p[1] << "," << p[2] << ")";
    }
}

TEST_F(AlcubierreMetricTests, ReducesToMinkowskiOutside) {
    // Far from bubble centre, f(rs) -> 0, metric -> Minkowski
    Metric4D g;
    evaluateAt(100.0, 0.0, 0.0, g);

    EXPECT_NEAR(g(0, 0).real, -1.0, 1e-6);
    EXPECT_NEAR(g(0, 1).real, 0.0, 1e-6);
    EXPECT_NEAR(g(1, 0).real, 0.0, 1e-6);
    EXPECT_NEAR(g(1, 1).real, 1.0, kEps);
    EXPECT_NEAR(g(2, 2).real, 1.0, kEps);
    EXPECT_NEAR(g(3, 3).real, 1.0, kEps);
}

TEST_F(AlcubierreMetricTests, ShiftVectorAtCentre) {
    // At bubble centre, f(0) ≈ 1, so g_tx = -vs * f ≈ -vs
    Metric4D g;
    evaluateAt(0.0, 0.0, 0.0, g);

    double vs = defaultParams.vs;
    double f = metric.shapeFunction(0.0);
    EXPECT_NEAR(f, 1.0, 1e-4); // shape function at centre ≈ 1

    EXPECT_NEAR(g(0, 1).real, -vs * f, kEps);
}

TEST_F(AlcubierreMetricTests, GttAtCentre) {
    // g_tt = -(1 - vs^2 * f^2)
    Metric4D g;
    evaluateAt(0.0, 0.0, 0.0, g);

    double vs = defaultParams.vs;
    double f = metric.shapeFunction(0.0);
    double expected_gtt = -(1.0 - vs * vs * f * f);
    EXPECT_NEAR(g(0, 0).real, expected_gtt, kEps);
}

TEST_F(AlcubierreMetricTests, SpatialComponentsFlat) {
    // g_xx = 1, g_yy = 1, g_zz = 1 everywhere (no spatial curvature)
    double positions[][3] = {{0.0, 0.0, 0.0}, {0.5, 0.3, 0.1}, {10.0, 5.0, 3.0}};
    for (auto& p : positions) {
        Metric4D g;
        evaluateAt(p[0], p[1], p[2], g);
        EXPECT_NEAR(g(1, 1).real, 1.0, kEps);
        EXPECT_NEAR(g(2, 2).real, 1.0, kEps);
        EXPECT_NEAR(g(3, 3).real, 1.0, kEps);

        // Off-diagonal spatial: g_xy = g_xz = g_yz = 0
        EXPECT_NEAR(g(1, 2).real, 0.0, kEps);
        EXPECT_NEAR(g(1, 3).real, 0.0, kEps);
        EXPECT_NEAR(g(2, 3).real, 0.0, kEps);
    }
}

// =============================================================================
// Shape Function Tests
// =============================================================================

TEST_F(AlcubierreMetricTests, ShapeFunctionAtCentre) {
    EXPECT_NEAR(metric.shapeFunction(0.0), 1.0, 1e-4);
}

TEST_F(AlcubierreMetricTests, ShapeFunctionFarField) {
    EXPECT_NEAR(metric.shapeFunction(10.0 * defaultParams.R), 0.0, 1e-6);
}

TEST_F(AlcubierreMetricTests, ShapeFunctionMonotoneDecrease) {
    double prev = metric.shapeFunction(0.0);
    for (double rs = 0.1; rs <= 5.0; rs += 0.1) {
        double val = metric.shapeFunction(rs);
        EXPECT_LE(val, prev + kEps)
            << "Shape function not monotone at rs = " << rs;
        prev = val;
    }
}

TEST_F(AlcubierreMetricTests, ShapeFunctionDerivativeFinite) {
    // No NaN in derivative, even near bubble wall
    for (double rs = 0.0; rs <= 5.0; rs += 0.05) {
        double df = metric.shapeFunctionDerivative(rs);
        EXPECT_FALSE(std::isnan(df)) << "NaN in derivative at rs = " << rs;
        EXPECT_FALSE(std::isinf(df)) << "Inf in derivative at rs = " << rs;
    }
}

// =============================================================================
// Parameter Configuration Tests
// =============================================================================

TEST_F(AlcubierreMetricTests, SubluminalConstruction) {
    auto params = Sirius::WarpDriveParams::Subluminal(2.0);
    EXPECT_LT(params.vs, 1.0);
    EXPECT_DOUBLE_EQ(params.R, 2.0);
}

TEST_F(AlcubierreMetricTests, SuperluminalConstruction) {
    auto params = Sirius::WarpDriveParams::Superluminal(3.0, 1.5);
    EXPECT_GT(params.vs, 1.0);
}

TEST_F(AlcubierreMetricTests, SetParameterClamping) {
    metric.setParameter("velocity", -5.0);
    auto params = metric.getParams();
    EXPECT_GE(params.vs, 0.0); // should be clamped to valid range
}

// =============================================================================
// Numerical Stability
// =============================================================================

TEST_F(AlcubierreMetricTests, NoNaNAtBubbleWall) {
    // Check metric at r_s = R (the wall)
    double R = defaultParams.R;
    Metric4D g;
    evaluateAt(R, 0.0, 0.0, g);

    for (int mu = 0; mu < 4; ++mu) {
        for (int nu = 0; nu < 4; ++nu) {
            EXPECT_FALSE(std::isnan(g(mu, nu).real))
                << "NaN at bubble wall, indices (" << mu << "," << nu << ")";
        }
    }
}

TEST_F(AlcubierreMetricTests, NoNaNAtOrigin) {
    Metric4D g;
    evaluateAt(0.0, 0.0, 0.0, g);
    for (int mu = 0; mu < 4; ++mu) {
        for (int nu = 0; nu < 4; ++nu) {
            EXPECT_FALSE(std::isnan(g(mu, nu).real))
                << "NaN at origin, indices (" << mu << "," << nu << ")";
        }
    }
}

TEST_F(AlcubierreMetricTests, BubblePositionUpdate) {
    metric.updateBubblePosition(5.0);
    auto params = metric.getParams();
    EXPECT_NEAR(params.xs, defaultParams.vs * 5.0, kEps);
}

} // namespace sirius::test
