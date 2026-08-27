// Production Kerr gates in Sirius's Cartesian Kerr-Schild chart.

#include "sirius/core/metrics/kerr_schild_family.h"

#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <numbers>

namespace sirius::test {
using namespace sirius::core;

namespace {

constexpr double kMass = 1.0;
constexpr double kSpin = 0.9;
constexpr double kTolerance = 3.0e-12;

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

TEST(KerrTests, MetricMatchesIndependentCartesianKerrSchildForm) {
    KerrSchildFamily metric(KerrSchildParams::Kerr(kMass, kSpin));
    const std::array<Vec4, 3> positions = {
        Position(8.0, 0.0, 0.0),
        Position(4.0, -3.0, 2.0),
        Position(-7.0, 5.0, 1.0),
    };
    for (const Vec4& position : positions) {
        const Metric4d actual = Evaluate(metric, position);
        const double x = position(1), y = position(2), z = position(3);
        const double radial2 = x * x + y * y + z * z;
        const double reduced = radial2 - kSpin * kSpin;
        const double r2 =
            0.5 * (reduced + std::sqrt(reduced * reduced + 4.0 * kSpin * kSpin * z * z));
        const double radius = std::sqrt(r2);
        const double h = 2.0 * kMass * radius * r2 / (r2 * r2 + kSpin * kSpin * z * z);
        const double l[4] = {
            1.0,
            (radius * x + kSpin * y) / (r2 + kSpin * kSpin),
            (radius * y - kSpin * x) / (r2 + kSpin * kSpin),
            z / radius,
        };
        for (int mu = 0; mu < 4; ++mu) {
            for (int nu = 0; nu < 4; ++nu) {
                const double eta = mu == nu ? (mu == 0 ? -1.0 : 1.0) : 0.0;
                EXPECT_NEAR(actual(mu, nu).real, eta + h * l[mu] * l[nu], kTolerance)
                    << "component (" << mu << "," << nu << ")";
            }
        }
    }
}

TEST(KerrTests, ComputedRadiusSatisfiesTheDefiningOblateQuartic) {
    KerrSchildFamily metric(KerrSchildParams::Kerr(kMass, kSpin));
    for (const Vec4& position :
         {Position(8.0, 0.0, 0.0), Position(4.0, -3.0, 2.0), Position(-2.0, 1.0, 5.0)}) {
        const double x = position(1), y = position(2), z = position(3);
        const double radius = metric.ComputeKerrRadius(x, y, z);
        const double r2 = radius * radius;
        const double cartesian2 = x * x + y * y + z * z;
        const double residual = r2 * r2 - (cartesian2 - kSpin * kSpin) * r2 - kSpin * kSpin * z * z;
        EXPECT_NEAR(residual, 0.0, 2.0e-10) << "position " << x << "," << y << "," << z;
    }
}

TEST(KerrTests, HorizonsAndCaptureFollowTheIndependentKerrPolynomial) {
    KerrSchildFamily metric(KerrSchildParams::Kerr(kMass, kSpin));
    ASSERT_TRUE(metric.HasHorizon());
    const double outer = metric.OuterHorizonRadius();
    const double inner = metric.InnerHorizonRadius();
    const auto polynomial = [](double radius) {
        return radius * radius - 2.0 * kMass * radius + kSpin * kSpin;
    };
    EXPECT_NEAR(polynomial(outer), 0.0, kTolerance);
    EXPECT_NEAR(polynomial(inner), 0.0, kTolerance);

    const auto equatorial_position = [](double radius) {
        return Position(std::sqrt(radius * radius + kSpin * kSpin), 0.0, 0.0);
    };
    EXPECT_TRUE(metric.InsideCaptureSurface(equatorial_position(outer * 0.999), 0.0));
    EXPECT_FALSE(metric.InsideCaptureSurface(equatorial_position(outer * 1.001), 0.0));
}

TEST(KerrTests, StaticLimitAuthorityIsAZeroOfProductionGtt) {
    KerrSchildFamily metric(KerrSchildParams::Kerr(kMass, kSpin));
    for (const double theta : {0.2, std::numbers::pi / 4.0, std::numbers::pi / 2.0}) {
        const double radius = metric.ErgosphereRadius(theta);
        ASSERT_TRUE(std::isfinite(radius));
        const double x = std::sqrt(radius * radius + kSpin * kSpin) * std::sin(theta);
        const double z = radius * std::cos(theta);
        const Metric4d value = Evaluate(metric, Position(x, 0.0, z));
        EXPECT_NEAR(value(0, 0).real, 0.0, 2.0e-12) << "theta=" << theta;
    }
}

TEST(KerrTests, ZeroSpinEqualsProductionSchwarzschild) {
    KerrSchildFamily kerr(KerrSchildParams::Kerr(kMass, 0.0));
    KerrSchildFamily schwarzschild(KerrSchildParams::Schwarzschild(kMass));
    for (const Vec4& position : {Position(8.0, 0.0, 0.0), Position(4.0, 3.0, -2.0)}) {
        Metric4d kerr_value, schwarzschild_value;
        Tensor<Dual<double>, 4, 4, 4> kerr_derivatives, schwarzschild_derivatives;
        kerr.Evaluate(position, kerr_value, kerr_derivatives);
        schwarzschild.Evaluate(position, schwarzschild_value, schwarzschild_derivatives);
        for (int axis = 0; axis < 4; ++axis) {
            for (int mu = 0; mu < 4; ++mu) {
                EXPECT_DOUBLE_EQ(kerr_value(axis, mu).real, schwarzschild_value(axis, mu).real);
                for (int nu = 0; nu < 4; ++nu) {
                    EXPECT_DOUBLE_EQ(kerr_derivatives(axis, mu, nu).real,
                                     schwarzschild_derivatives(axis, mu, nu).real);
                }
            }
        }
    }
}

TEST(KerrTests, ClosedFormInverseMultipliesToIdentity) {
    KerrSchildFamily metric(KerrSchildParams::Kerr(kMass, kSpin));
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
