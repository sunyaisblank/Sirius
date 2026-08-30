// Production Schwarzschild gates in Sirius's Cartesian Kerr-Schild chart.

#include "sirius/core/metrics/kerr_schild_family.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <numbers>

namespace sirius::test {
using namespace sirius::core;

namespace {

constexpr double kMass = 1.0;
constexpr double kTolerance = 2.0e-12;

Vec4 Position(double x, double y, double z) {
    Vec4 position;
    position(0) = 0.0;
    position(1) = x;
    position(2) = y;
    position(3) = z;
    return position;
}

Metric4d Evaluate(KerrSchildFamily& metric, const Vec4& position) {
    Metric4d value;
    Tensor<Dual<double>, 4, 4, 4> derivatives;
    metric.Evaluate(position, value, derivatives);
    return value;
}

double MaxMinkowskiDifference(const Metric4d& metric) {
    double maximum = 0.0;
    for (int mu = 0; mu < 4; ++mu) {
        for (int nu = 0; nu < 4; ++nu) {
            const double eta = mu == nu ? (mu == 0 ? -1.0 : 1.0) : 0.0;
            maximum = std::max(maximum, std::abs(metric(mu, nu).real - eta));
        }
    }
    return maximum;
}

}  // namespace

TEST(SchwarzschildTests, MetricMatchesIndependentCartesianKerrSchildForm) {
    KerrSchildFamily metric(KerrSchildParams::Schwarzschild(kMass));
    const std::array<Vec4, 3> positions = {
        Position(10.0, 0.0, 0.0),
        Position(4.0, -3.0, 2.0),
        Position(-7.0, 5.0, 1.0),
    };
    for (const Vec4& position : positions) {
        const Metric4d actual = Evaluate(metric, position);
        const double x = position(1), y = position(2), z = position(3);
        const double radius = std::sqrt(x * x + y * y + z * z);
        const double h = 2.0 * kMass / radius;
        const double l[4] = {1.0, x / radius, y / radius, z / radius};
        for (int mu = 0; mu < 4; ++mu) {
            for (int nu = 0; nu < 4; ++nu) {
                const double eta = mu == nu ? (mu == 0 ? -1.0 : 1.0) : 0.0;
                EXPECT_NEAR(actual(mu, nu).real, eta + h * l[mu] * l[nu], kTolerance)
                    << "component (" << mu << "," << nu << ")";
            }
        }
    }
}

TEST(SchwarzschildTests, AnalyticDerivativesMatchIndependentFiniteDifferences) {
    KerrSchildFamily metric(KerrSchildParams::Schwarzschild(kMass));
    const Vec4 position = Position(8.0, -2.0, 3.0);
    Metric4d value;
    Tensor<Dual<double>, 4, 4, 4> derivatives;
    metric.Evaluate(position, value, derivatives);

    constexpr double epsilon = 1.0e-5;
    for (int axis = 1; axis < 4; ++axis) {
        Vec4 plus = position;
        Vec4 minus = position;
        plus(axis) += epsilon;
        minus(axis) -= epsilon;
        const Metric4d g_plus = Evaluate(metric, plus);
        const Metric4d g_minus = Evaluate(metric, minus);
        for (int mu = 0; mu < 4; ++mu) {
            for (int nu = 0; nu < 4; ++nu) {
                const double finite_difference =
                    (g_plus(mu, nu).real - g_minus(mu, nu).real) / (2.0 * epsilon);
                EXPECT_NEAR(derivatives(axis, mu, nu).real, finite_difference, 5.0e-9)
                    << "derivative (" << axis << "," << mu << "," << nu << ")";
            }
        }
    }
}

TEST(SchwarzschildTests, HorizonAndCaptureUseTheExactArealRadius) {
    KerrSchildFamily metric(KerrSchildParams::Schwarzschild(kMass));
    ASSERT_TRUE(metric.HasHorizon());
    EXPECT_DOUBLE_EQ(metric.OuterHorizonRadius(), 2.0 * kMass);
    EXPECT_DOUBLE_EQ(metric.InnerHorizonRadius(), 0.0);
    EXPECT_TRUE(metric.InsideCaptureSurface(Position(1.999, 0.0, 0.0), 0.0));
    EXPECT_FALSE(metric.InsideCaptureSurface(Position(2.001, 0.0, 0.0), 0.0));
}

TEST(KottlerHorizonTests, ExactRootsSeparateBlackHoleCaptureFromCosmologicalHorizon) {
    constexpr double lambda = 0.01;
    KerrSchildParams parameters = KerrSchildParams::Schwarzschild(kMass);
    parameters.Lambda = lambda;
    KerrSchildFamily metric(parameters);

    const double root_angle = std::acos(-3.0 * kMass * std::sqrt(lambda));
    const double expected_black_hole =
        2.0 / std::sqrt(lambda) * std::cos((root_angle + 4.0 * std::numbers::pi) / 3.0);
    const double expected_cosmological = 2.0 / std::sqrt(lambda) * std::cos(root_angle / 3.0);
    const double black_hole = metric.OuterHorizonRadius();
    const double cosmological = metric.CosmologicalHorizonRadius();

    ASSERT_TRUE(metric.HasHorizon());
    EXPECT_NEAR(black_hole, expected_black_hole, 2.0e-14);
    EXPECT_NEAR(cosmological, expected_cosmological, 2.0e-14);
    EXPECT_GT(black_hole, 2.0 * kMass);
    EXPECT_LT(cosmological, std::sqrt(3.0 / lambda));
    EXPECT_NEAR(metric.KottlerStaticLapse(black_hole), 0.0, 2.0e-16);
    EXPECT_NEAR(metric.KottlerStaticLapse(cosmological), 0.0, 2.0e-15);
    EXPECT_GT(metric.KottlerStaticLapse(0.5 * (black_hole + cosmological)), 0.0);
    EXPECT_TRUE(metric.InsideCaptureSurface(Position(0.999 * black_hole, 0.0, 0.0), 0.0));
    EXPECT_FALSE(metric.InsideCaptureSurface(Position(1.001 * black_hole, 0.0, 0.0), 0.0));
    EXPECT_DEATH((void)metric.ErgosphereRadius(0.0), "precondition.*enforced, terminating");
    EXPECT_DEATH((void)metric.IscoRadius(), "precondition.*enforced, terminating");
    EXPECT_DEATH((void)metric.ExtremalityParameter(), "precondition.*enforced, terminating");
}

TEST(KottlerHorizonTests, PureDeSitterHasOnlyACosmologicalHorizon) {
    constexpr double lambda = 0.01;
    KerrSchildFamily metric(KerrSchildParams::DeSitter(lambda));
    EXPECT_FALSE(metric.HasHorizon());
    EXPECT_EQ(metric.OuterHorizonRadius(), -1.0);
    EXPECT_EQ(metric.InnerHorizonRadius(), -1.0);
    EXPECT_DOUBLE_EQ(metric.CosmologicalHorizonRadius(), std::sqrt(3.0 / lambda));
    EXPECT_NEAR(metric.KottlerStaticLapse(metric.CosmologicalHorizonRadius()), 0.0, 3.0e-16);
    EXPECT_FALSE(metric.InsideCaptureSurface(Position(0.001, 0.0, 0.0), 0.0));
}

TEST(KottlerHorizonTests, RegularDeSitterOriginUsesTheExactAnalyticLimit) {
    KerrSchildFamily de_sitter(KerrSchildParams::DeSitter(0.01));
    KerrSchildFamily minkowski(KerrSchildParams::Minkowski());
    KerrSchildFamily schwarzschild(KerrSchildParams::Schwarzschild(kMass));
    const Vec4 origin = Position(0.0, 0.0, 0.0);

    EXPECT_TRUE(de_sitter.IsValidEvent(origin));
    EXPECT_TRUE(minkowski.IsValidEvent(origin));
    EXPECT_FALSE(schwarzschild.IsValidEvent(origin));

    Metric4d value, inverse;
    Tensor<Dual<double>, 4, 4, 4> derivatives;
    de_sitter.Evaluate(origin, value, derivatives);
    ASSERT_TRUE(de_sitter.InverseMetric(origin, inverse));
    for (int mu = 0; mu < 4; ++mu) {
        for (int nu = 0; nu < 4; ++nu) {
            const double eta = mu == nu ? (mu == 0 ? -1.0 : 1.0) : 0.0;
            EXPECT_DOUBLE_EQ(value(mu, nu).real, eta);
            EXPECT_DOUBLE_EQ(inverse(mu, nu).real, eta);
            for (int axis = 0; axis < 4; ++axis) {
                EXPECT_DOUBLE_EQ(derivatives(axis, mu, nu).real, 0.0);
            }
        }
    }

    EXPECT_DEATH((void)Evaluate(schwarzschild, origin), "precondition.*enforced, terminating");
}

TEST(SchwarzschildTests, FarFieldPerturbationHasExactInverseRadiusScaling) {
    KerrSchildFamily metric(KerrSchildParams::Schwarzschild(kMass));
    const double near_difference = MaxMinkowskiDifference(Evaluate(metric, Position(100.0, 0, 0)));
    const double far_difference = MaxMinkowskiDifference(Evaluate(metric, Position(200.0, 0, 0)));
    ASSERT_GT(far_difference, 0.0);
    EXPECT_NEAR(near_difference / far_difference, 2.0, kTolerance);
}

TEST(SchwarzschildTests, StationaryClockRateMatchesTheExactStaticPotential) {
    KerrSchildFamily metric(KerrSchildParams::Schwarzschild(kMass));
    for (const double radius : {3.0, 5.0, 10.0, 100.0}) {
        const Metric4d value = Evaluate(metric, Position(radius, 0.0, 0.0));
        ASSERT_LT(value(0, 0).real, 0.0);
        EXPECT_NEAR(std::sqrt(-value(0, 0).real), std::sqrt(1.0 - 2.0 * kMass / radius),
                    kTolerance);
    }
}

}  // namespace sirius::test
