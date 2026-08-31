// TSMT102A.cpp - Alcubierre Warp Drive Metric Property Tests
// Component ID: TSMT102A (Test/Metric/AlcubierreProperties)
// Tests: PHMT102A.h (WarpDriveFamily)

#include "sirius/core/dual_number.h"
#include "sirius/core/metrics/warp_drive_family.h"
#include "sirius/core/tensor.h"

#include <gtest/gtest.h>

#include "metric_curvature_oracle.h"

#include <array>
#include <cmath>
#include <numbers>

namespace sirius::test {
using namespace sirius::core;

constexpr double kEps = 1e-8;

// =============================================================================
// Test Fixture
// =============================================================================

class AlcubierreMetricTests : public ::testing::Test {
  protected:
    sirius::core::WarpDriveFamily metric;
    sirius::core::WarpDriveParams defaultParams;

    void SetUp() override {
        defaultParams = sirius::core::WarpDriveParams::Alcubierre(1.0, 1.0);
        metric.SetParams(defaultParams);
    }

    void evaluateAt(double x, double y, double z, Metric4d& g) {
        Tensor<double, 4> pos;
        pos(0) = 0.0;
        pos(1) = x;
        pos(2) = y;
        pos(3) = z;
        Tensor<Dual<double>, 4, 4, 4> dg;
        metric.Evaluate(pos, g, dg);
    }
};

// =============================================================================
// Metric Symmetry Tests
// =============================================================================

TEST_F(AlcubierreMetricTests, MetricSymmetry) {
    // Test symmetry g_mu_nu = g_nu_mu at multiple positions
    double positions[][3] = {{0.0, 0.0, 0.0}, {0.5, 0.0, 0.0}, {2.0, 1.0, 0.5}, {10.0, 0.0, 0.0}};
    for (auto& p : positions) {
        Metric4d g;
        evaluateAt(p[0], p[1], p[2], g);
        for (int mu = 0; mu < 4; ++mu) {
            for (int nu = mu + 1; nu < 4; ++nu) {
                EXPECT_NEAR(g(mu, nu).real, g(nu, mu).real, kEps)
                    << "Asymmetry at (" << p[0] << "," << p[1] << "," << p[2] << ") indices (" << mu
                    << "," << nu << ")";
            }
        }
    }
}

TEST_F(AlcubierreMetricTests, LorentzianSignature) {
    // At any point, the metric should have Lorentzian signature (-,+,+,+)
    double positions[][3] = {{0.0, 0.0, 0.0}, {0.5, 0.0, 0.0}, {5.0, 0.0, 0.0}};
    for (auto& p : positions) {
        Metric4d g;
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
    Metric4d g;
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
    Metric4d g;
    evaluateAt(0.0, 0.0, 0.0, g);

    double vs = defaultParams.vs;
    double f = metric.ShapeFunction(0.0);
    EXPECT_NEAR(f, 1.0, 1e-4);  // shape function at centre ≈ 1

    EXPECT_NEAR(g(0, 1).real, -vs * f, kEps);
}

TEST_F(AlcubierreMetricTests, GttAtCentre) {
    // g_tt = -(1 - vs^2 * f^2)
    Metric4d g;
    evaluateAt(0.0, 0.0, 0.0, g);

    double vs = defaultParams.vs;
    double f = metric.ShapeFunction(0.0);
    double expected_gtt = -(1.0 - vs * vs * f * f);
    EXPECT_NEAR(g(0, 0).real, expected_gtt, kEps);
}

TEST_F(AlcubierreMetricTests, SpatialComponentsFlat) {
    // g_xx = 1, g_yy = 1, g_zz = 1 everywhere (no spatial curvature)
    double positions[][3] = {{0.0, 0.0, 0.0}, {0.5, 0.3, 0.1}, {10.0, 5.0, 3.0}};
    for (auto& p : positions) {
        Metric4d g;
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

TEST_F(AlcubierreMetricTests, AnalyticInverseMatchesExactUnitDeterminantBlock) {
    for (const auto& params : {
             WarpDriveParams::Alcubierre(0.5, 1.0, 0.5),
             WarpDriveParams::Alcubierre(2.0, 1.0, 8.0),
             WarpDriveParams::Alcubierre(-1.25, 4.0, 0.5),
         }) {
        metric.SetParams(params);
        for (const auto& point : {
                 std::array<double, 4>{0.0, 0.0, 0.0, 0.0},
                 std::array<double, 4>{0.25, 1.0, -0.3, 0.7},
                 std::array<double, 4>{-0.5, 8.0, 2.0, -1.0},
             }) {
            Vec4 position;
            for (int component = 0; component < 4; ++component) {
                position(component) = point[component];
            }
            Metric4d g;
            Tensor<Dual<double>, 4, 4, 4> dg;
            metric.Evaluate(position, g, dg);
            Metric4d inverse;
            ASSERT_TRUE(metric.InverseMetric(position, inverse));
            for (int mu = 0; mu < 4; ++mu) {
                for (int nu = 0; nu < 4; ++nu) {
                    double product = 0.0;
                    for (int alpha = 0; alpha < 4; ++alpha) {
                        product += inverse(mu, alpha).real * g(alpha, nu).real;
                    }
                    EXPECT_NEAR(product, mu == nu ? 1.0 : 0.0, 3.0e-14);
                }
            }
        }
    }
}

// With unit lapse, flat spatial slices, and shift beta^x=-v_s f, the ADM
// Hamiltonian constraint gives the Eulerian energy density (G=c=1)
//
//   T_mu_nu n^mu n^nu = -v_s^2 [(d_y f)^2 + (d_z f)^2] / (32 pi).
//
// A fourth-order numerical Ricci oracle reconstructs G_mu_nu from the live
// metric/connection and compares it with that independent closed form. This
// proves the defining negative-energy warp wall rather than only its metric
// component layout.
TEST_F(AlcubierreMetricTests, EulerianWarpWallEnergyMatchesEinsteinConstraint) {
    constexpr std::array<std::array<double, 3>, 3> directions = {{
        {0.0, 1.0, 0.0},
        {0.707106781186548, 0.707106781186548, 0.0},
        {0.267261241912424, -0.534522483824849, 0.801783725737273},
    }};
    for (const auto& params : {
             WarpDriveParams::Alcubierre(0.5, 1.0, 1.0),
             WarpDriveParams::Alcubierre(2.0, 1.0, 8.0),
             WarpDriveParams::Alcubierre(-1.25, 4.0, 0.5),
         }) {
        metric.SetParams(params);
        const double feature_scale = std::min(params.R, 1.0 / params.sigma);
        for (const auto& direction : directions) {
            Vec4 position;
            position(1) = params.xs + params.R * direction[0];
            position(2) = params.ys + params.R * direction[1];
            position(3) = params.zs + params.R * direction[2];

            Metric4d g;
            Tensor<Dual<double>, 4, 4, 4> dg;
            metric.Evaluate(position, g, dg);
            const double f = metric.ShapeFunction(params.R);
            const double normal[4] = {1.0, params.vs * f, 0.0, 0.0};
            double normal_norm = 0.0;
            for (int mu = 0; mu < 4; ++mu) {
                for (int nu = 0; nu < 4; ++nu) {
                    normal_norm += g(mu, nu).real * normal[mu] * normal[nu];
                }
            }
            EXPECT_NEAR(normal_norm, -1.0, 3.0e-14);

            const sirius::test_support::RicciSample ricci =
                sirius::test_support::RicciFromConnectionFiniteDifference(metric, position,
                                                                          2.0e-4 * feature_scale);
            double ricci_normal = 0.0;
            for (int mu = 0; mu < 4; ++mu) {
                for (int nu = 0; nu < 4; ++nu) {
                    ricci_normal += ricci.covariant[mu][nu] * normal[mu] * normal[nu];
                }
            }
            const double measured_density =
                (ricci_normal + 0.5 * ricci.scalar) / (8.0 * std::numbers::pi);
            const double df_drs = metric.ShapeFunctionDerivative(params.R);
            const double transverse_fraction =
                direction[1] * direction[1] + direction[2] * direction[2];
            const double expected_density = -params.vs * params.vs * df_drs * df_drs *
                                            transverse_fraction / (32.0 * std::numbers::pi);
            const double tolerance = std::max(2.0e-8 / (feature_scale * feature_scale),
                                              2.0e-4 * std::abs(expected_density));
            EXPECT_NEAR(measured_density, expected_density, tolerance)
                << "vs=" << params.vs << " R=" << params.R << " sigma=" << params.sigma;
            EXPECT_LT(measured_density, 0.0)
                << "every sampled off-axis warp-wall observer must measure negative energy";
        }
    }
}

// =============================================================================
// Shape Function Tests
// =============================================================================

TEST_F(AlcubierreMetricTests, ShapeFunctionAtCentre) {
    EXPECT_NEAR(metric.ShapeFunction(0.0), 1.0, 1e-4);
}

TEST_F(AlcubierreMetricTests, ShapeFunctionFarField) {
    const double exterior_scale = std::max(defaultParams.R, 1.0 / defaultParams.sigma);
    EXPECT_NEAR(metric.ShapeFunction(10.0 * exterior_scale), 0.0, 1e-6);
}

TEST_F(AlcubierreMetricTests, ShapeFunctionMonotoneDecrease) {
    double prev = metric.ShapeFunction(0.0);
    for (double rs = 0.1; rs <= 5.0; rs += 0.1) {
        double val = metric.ShapeFunction(rs);
        EXPECT_LE(val, prev + kEps) << "Shape function not monotone at rs = " << rs;
        prev = val;
    }
}

TEST_F(AlcubierreMetricTests, ShapeFunctionDerivativeFinite) {
    // No NaN in derivative, even near bubble wall
    for (double rs = 0.0; rs <= 5.0; rs += 0.05) {
        double df = metric.ShapeFunctionDerivative(rs);
        EXPECT_FALSE(std::isnan(df)) << "NaN in derivative at rs = " << rs;
        EXPECT_FALSE(std::isinf(df)) << "Inf in derivative at rs = " << rs;
    }
}

// =============================================================================
// Parameter Configuration Tests
// =============================================================================

TEST_F(AlcubierreMetricTests, SubluminalConstruction) {
    auto params = sirius::core::WarpDriveParams::Subluminal(2.0);
    EXPECT_LT(params.vs, 1.0);
    EXPECT_DOUBLE_EQ(params.R, 2.0);
}

TEST_F(AlcubierreMetricTests, SuperluminalConstruction) {
    auto params = sirius::core::WarpDriveParams::Superluminal(3.0, 1.5);
    EXPECT_GT(params.vs, 1.0);
}

TEST_F(AlcubierreMetricTests, SetParameterMatchesTheOperatorVelocityDomain) {
    metric.SetParameter("velocity", -5.0);
    auto params = metric.GetParams();
    EXPECT_DOUBLE_EQ(params.vs, -5.0);
}

TEST_F(AlcubierreMetricTests, DirectMetricRejectsUnresolvedProfiles) {
    EXPECT_DEATH(
        {
            WarpDriveFamily invalid(WarpDriveParams::Alcubierre(1.0, 1.0, 0.05));
            static_cast<void>(invalid);
        },
        "violated");
    EXPECT_DEATH(
        {
            WarpDriveFamily invalid(WarpDriveParams::Alcubierre(1.0, 1.0, 101.0));
            static_cast<void>(invalid);
        },
        "violated");
    EXPECT_DEATH(metric.SetParameter("sigma", 0.05), "violated");

    WarpDriveFamily lower(WarpDriveParams::Alcubierre(1.0, 1.0, 0.1));
    WarpDriveFamily upper(WarpDriveParams::Alcubierre(1.0, 1.0, 100.0));
    EXPECT_DOUBLE_EQ(lower.GetParams().sigma * lower.GetParams().R, kMinAlcubierreSigmaRadius);
    EXPECT_DOUBLE_EQ(upper.GetParams().sigma * upper.GetParams().R, kMaxAlcubierreSigmaRadius);
}

// =============================================================================
// Numerical Stability
// =============================================================================

TEST_F(AlcubierreMetricTests, NoNaNAtBubbleWall) {
    // Check metric at r_s = R (the wall)
    double R = defaultParams.R;
    Metric4d g;
    evaluateAt(R, 0.0, 0.0, g);

    for (int mu = 0; mu < 4; ++mu) {
        for (int nu = 0; nu < 4; ++nu) {
            EXPECT_FALSE(std::isnan(g(mu, nu).real))
                << "NaN at bubble wall, indices (" << mu << "," << nu << ")";
        }
    }
}

TEST_F(AlcubierreMetricTests, NoNaNAtOrigin) {
    Metric4d g;
    evaluateAt(0.0, 0.0, 0.0, g);
    for (int mu = 0; mu < 4; ++mu) {
        for (int nu = 0; nu < 4; ++nu) {
            EXPECT_FALSE(std::isnan(g(mu, nu).real))
                << "NaN at origin, indices (" << mu << "," << nu << ")";
        }
    }
}

TEST_F(AlcubierreMetricTests, BubblePositionUpdate) {
    metric.UpdateBubblePosition(5.0);
    auto params = metric.GetParams();
    EXPECT_NEAR(params.xs, defaultParams.vs * 5.0, kEps);
}

}  // namespace sirius::test
