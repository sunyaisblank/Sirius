// Production Reissner-Nordstrom gates.  Expected values are formed from the
// independent Cartesian Kerr-Schild definition, never from a test-local metric.

#include "sirius/core/metrics/kerr_schild_family.h"

#include <gtest/gtest.h>

#include <array>
#include <cmath>

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

}  // namespace

TEST(ReissnerNordstromTests, MetricMatchesIndependentCartesianKerrSchildForm) {
    constexpr double charge = 0.6;
    KerrSchildFamily metric(KerrSchildParams::ReissnerNordstrom(kMass, charge));
    const std::array<Vec4, 3> positions = {
        Position(8.0, 0.0, 0.0),
        Position(4.0, -3.0, 2.0),
        Position(-6.0, 2.0, 5.0),
    };

    for (const Vec4& position : positions) {
        Metric4d actual;
        Tensor<Dual<double>, 4, 4, 4> derivatives;
        metric.Evaluate(position, actual, derivatives);

        const double x = position(1), y = position(2), z = position(3);
        const double radius = std::sqrt(x * x + y * y + z * z);
        const double h = 2.0 * kMass / radius - charge * charge / (radius * radius);
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

TEST(ReissnerNordstromTests, AnalyticDerivativesMatchIndependentFiniteDifferences) {
    KerrSchildFamily metric(KerrSchildParams::ReissnerNordstrom(kMass, 0.7));
    const Vec4 position = Position(7.0, -2.0, 3.0);
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

TEST(ReissnerNordstromTests, ZeroChargeEqualsProductionSchwarzschild) {
    KerrSchildFamily charged(KerrSchildParams::ReissnerNordstrom(kMass, 0.0));
    KerrSchildFamily schwarzschild(KerrSchildParams::Schwarzschild(kMass));
    for (const Vec4& position : {Position(8.0, 0.0, 0.0), Position(4.0, 3.0, -2.0)}) {
        Metric4d charged_value, schwarzschild_value;
        Tensor<Dual<double>, 4, 4, 4> charged_derivatives, schwarzschild_derivatives;
        charged.Evaluate(position, charged_value, charged_derivatives);
        schwarzschild.Evaluate(position, schwarzschild_value, schwarzschild_derivatives);
        for (int axis = 0; axis < 4; ++axis) {
            for (int mu = 0; mu < 4; ++mu) {
                EXPECT_DOUBLE_EQ(charged_value(axis, mu).real, schwarzschild_value(axis, mu).real);
                for (int nu = 0; nu < 4; ++nu) {
                    EXPECT_DOUBLE_EQ(charged_derivatives(axis, mu, nu).real,
                                     schwarzschild_derivatives(axis, mu, nu).real);
                }
            }
        }
    }
}

TEST(ReissnerNordstromTests, HorizonAuthoritiesSatisfyTheIndependentPolynomial) {
    for (const double charge : {0.0, 0.3, 0.7, 0.99}) {
        KerrSchildFamily metric(KerrSchildParams::ReissnerNordstrom(kMass, charge));
        ASSERT_TRUE(metric.HasHorizon());
        const double outer = metric.OuterHorizonRadius();
        const double inner = metric.InnerHorizonRadius();
        const auto horizon_polynomial = [charge](double radius) {
            return radius * radius - 2.0 * kMass * radius + charge * charge;
        };
        EXPECT_NEAR(horizon_polynomial(outer), 0.0, kTolerance);
        EXPECT_NEAR(horizon_polynomial(inner), 0.0, kTolerance);
        EXPECT_NEAR(outer + inner, 2.0 * kMass, kTolerance);
        EXPECT_NEAR(outer * inner, charge * charge, kTolerance);

        EXPECT_TRUE(metric.InsideCaptureSurface(Position(outer * 0.999, 0.0, 0.0), 0.0));
        EXPECT_FALSE(metric.InsideCaptureSurface(Position(outer * 1.001, 0.0, 0.0), 0.0));
    }
}

TEST(ReissnerNordstromTests, SuperExtremalDomainHasNoHorizonOrCaptureSurface) {
    KerrSchildFamily metric(KerrSchildParams::ReissnerNordstrom(kMass, 1.01));
    EXPECT_FALSE(metric.HasHorizon());
    EXPECT_EQ(metric.OuterHorizonRadius(), -1.0);
    EXPECT_EQ(metric.InnerHorizonRadius(), -1.0);
    EXPECT_GT(metric.ExtremalityParameter(), 1.0);
    EXPECT_FALSE(metric.InsideCaptureSurface(Position(0.5, 0.0, 0.0), 0.0));
}

TEST(ReissnerNordstromTests, ClosedFormInverseMultipliesToIdentity) {
    KerrSchildFamily metric(KerrSchildParams::ReissnerNordstrom(kMass, 0.8));
    const Vec4 position = Position(5.0, -2.0, 1.0);
    const Metric4d covariant = Evaluate(metric, position);
    Metric4d inverse;
    ASSERT_TRUE(metric.InverseMetric(position, inverse));

    for (int row = 0; row < 4; ++row) {
        for (int column = 0; column < 4; ++column) {
            double product = 0.0;
            for (int inner = 0; inner < 4; ++inner) {
                product += covariant(row, inner).real * inverse(inner, column).real;
            }
            EXPECT_NEAR(product, row == column ? 1.0 : 0.0, kTolerance)
                << "matrix product (" << row << "," << column << ")";
        }
    }
}

}  // namespace sirius::test
