// MorrisThorneCartesian: the exact isotropic Cartesian chart of the zero-tidal
// Ellis member, the chart the CPU tracer drives.
//
// Gates: coordinate-transformed agreement with the spherical areal-radius
// authority on both sheets, a finite exact throat, analytic derivatives,
// inverse identity, chart domain, and complete accepted-segment capture
// events.  The latter includes equal-sign endpoints and a tangent contact.

#include "sirius/core/metrics/morris_thorne_family.h"
#include "sirius/core/tensor.h"
#include "sirius/core/trace_boundary.h"

#include <gtest/gtest.h>

#include <cmath>
#include <limits>
#include <numbers>
#include <vector>

namespace {

using namespace sirius::core;

struct Direction {
    double x, y, z;
};

// Unit-normalised sample directions, none axis-aligned except the last, so
// the derivative and inverse gates see generic points of the chart.
std::vector<Direction> SampleDirections() {
    std::vector<Direction> dirs = {
        {1.0, 0.0, 0.0},
        {0.6, 0.8, 0.0},
        {0.577350269189626, 0.577350269189626, 0.577350269189626},
        {-0.267261241912424, 0.534522483824849, 0.801783725737273},
        {0.0, 0.0, 1.0},
    };
    return dirs;
}

Vec4 CartPos(double r, const Direction& n) {
    Vec4 p;
    p(1) = r * n.x;
    p(2) = r * n.y;
    p(3) = r * n.z;
    return p;
}

// A unit vector orthogonal to n (any one; the metric is isotropic about n).
Direction Orthogonal(const Direction& n) {
    // Cross with the axis least aligned with n.
    double ax = std::abs(n.x), ay = std::abs(n.y), az = std::abs(n.z);
    double ex = 0, ey = 0, ez = 0;
    if (ax <= ay && ax <= az)
        ex = 1;
    else if (ay <= az)
        ey = 1;
    else
        ez = 1;
    double cx = n.y * ez - n.z * ey;
    double cy = n.z * ex - n.x * ez;
    double cz = n.x * ey - n.y * ex;
    double len = std::sqrt(cx * cx + cy * cy + cz * cz);
    return {cx / len, cy / len, cz / len};
}

std::vector<MorrisThorneParams> SampleShapes() {
    return {
        MorrisThorneParams::Ellis(1.0),
        MorrisThorneParams::Ellis(2.5),
    };
}

// -----------------------------------------------------------------------------
// Under r=rho+b0^2/(4rho), the areal-radius radial coefficient transforms with
// (dr/drho)^2 and an isotropic-coordinate unit tangent spans angle 1/rho. Both
// must equal A^2 for A=1+b0^2/(4rho^2), including the mathematically continued
// second sheet rho<b0/2.
// -----------------------------------------------------------------------------
TEST(MorrisThorneCartesianTests, ChartAgreementWithSphericalFamily) {
    for (const auto& params : SampleShapes()) {
        MorrisThorneCartesian cart(params);
        MorrisThorneFamily sph(params);
        double b0 = params.b0;

        for (double rho_factor : {0.25, 0.55, 0.75, 1.0, 2.0, 10.0}) {
            const double rho = rho_factor * b0;
            const double q = (b0 * b0) / (4.0 * rho * rho);
            const double areal_radius = rho * (1.0 + q);
            const double dr_drho = 1.0 - q;
            const double expected_conformal = (1.0 + q) * (1.0 + q);
            for (const auto& n : SampleDirections()) {
                Metric4d g;
                Tensor<Dual<double>, 4, 4, 4> dg;
                cart.Evaluate(CartPos(rho, n), g, dg);

                Metric4d gs;
                Tensor<Dual<double>, 4, 4, 4> dgs;
                Vec4 spos;
                spos(1) = areal_radius;
                spos(2) = std::numbers::pi / 2.0;
                sph.Evaluate(spos, gs, dgs);

                double nn[3] = {n.x, n.y, n.z};
                Direction u = Orthogonal(n);
                double uu[3] = {u.x, u.y, u.z};

                double radial = 0.0, tangential = 0.0, cross = 0.0;
                for (int i = 0; i < 3; ++i) {
                    for (int j = 0; j < 3; ++j) {
                        double gij = g(i + 1, j + 1).real;
                        radial += nn[i] * nn[j] * gij;
                        tangential += uu[i] * uu[j] * gij;
                        cross += nn[i] * uu[j] * gij;
                    }
                }

                const double transformed_radial = gs(1, 1).real * dr_drho * dr_drho;
                const double transformed_tangential = gs(2, 2).real / (rho * rho);
                EXPECT_NEAR(radial, transformed_radial, 2e-11 * std::abs(transformed_radial))
                    << "radial component, b0=" << b0 << " rho=" << rho;
                EXPECT_NEAR(tangential, transformed_tangential,
                            2e-12 * std::abs(transformed_tangential))
                    << "tangential component, b0=" << b0 << " rho=" << rho;
                EXPECT_NEAR(radial, expected_conformal, 2e-11 * expected_conformal);
                EXPECT_NEAR(tangential, expected_conformal, 2e-12 * expected_conformal);
                EXPECT_NEAR(cross, 0.0, 1e-12) << "radial-tangential cross term";
                EXPECT_NEAR(g(0, 0).real, gs(0, 0).real, 1e-14) << "g_tt";
                // Time-space block is exactly zero (static metric).
                for (int i = 1; i < 4; ++i) {
                    EXPECT_EQ(g(0, i).real, 0.0);
                    EXPECT_EQ(g(i, 0).real, 0.0);
                }
            }
        }
    }
}

// -----------------------------------------------------------------------------
// The analytic derivative block matches central finite differences of
// Evaluate at generic points. This is the completeness gate for the chain
// rule through the isotropic conformal factor.
// -----------------------------------------------------------------------------
TEST(MorrisThorneCartesianTests, DerivativesMatchFiniteDifferencesOfMetric) {
    const double h = 1e-6;
    for (const auto& params : SampleShapes()) {
        MorrisThorneCartesian cart(params);
        double b0 = params.b0;

        for (double rf : {0.25, 0.5, 0.75, 1.5, 3.0, 10.0}) {
            for (const auto& n : SampleDirections()) {
                Vec4 x = CartPos(rf * b0, n);

                Metric4d g;
                Tensor<Dual<double>, 4, 4, 4> dg;
                cart.Evaluate(x, g, dg);

                for (int k = 0; k < 4; ++k) {
                    Vec4 xp = x, xm = x;
                    xp(k) += h;
                    xm(k) -= h;
                    Metric4d gp, gm;
                    Tensor<Dual<double>, 4, 4, 4> scratch;
                    cart.Evaluate(xp, gp, scratch);
                    cart.Evaluate(xm, gm, scratch);

                    for (int i = 0; i < 4; ++i) {
                        for (int j = 0; j < 4; ++j) {
                            double fd = (gp(i, j).real - gm(i, j).real) / (2.0 * h);
                            EXPECT_NEAR(dg(k, i, j).real, fd, 5e-7 * std::max(1.0, std::abs(fd)))
                                << "d_" << k << " g_" << i << j << " at b0=" << b0
                                << " r=" << rf * b0;
                        }
                    }
                }
            }
        }
    }
}

// -----------------------------------------------------------------------------
// The closed-form diagonal inverse satisfies g^mu_alpha g_alpha_nu = delta^mu_nu
// to near machine precision, so the bar below is rounding.
// -----------------------------------------------------------------------------
TEST(MorrisThorneCartesianTests, AnalyticInverseIsExact) {
    for (const auto& params : SampleShapes()) {
        MorrisThorneCartesian cart(params);
        double b0 = params.b0;

        for (double rf : {0.25, 0.5, 0.75, 1.5, 3.0, 50.0}) {
            for (const auto& n : SampleDirections()) {
                Vec4 x = CartPos(rf * b0, n);

                Metric4d g;
                Tensor<Dual<double>, 4, 4, 4> dg;
                cart.Evaluate(x, g, dg);

                Metric4d gi;
                ASSERT_TRUE(cart.InverseMetric(x, gi));

                for (int mu = 0; mu < 4; ++mu) {
                    for (int nu = 0; nu < 4; ++nu) {
                        double s = 0.0;
                        for (int al = 0; al < 4; ++al) s += gi(mu, al).real * g(al, nu).real;
                        double expected = (mu == nu) ? 1.0 : 0.0;
                        EXPECT_NEAR(s, expected, 1e-13)
                            << "inverse identity at b0=" << b0 << " r=" << rf * b0;
                    }
                }
            }
        }
    }
}

// -----------------------------------------------------------------------------
// The physical areal throat r=b0 is rho=b0/2 in the isotropic chart. It is a
// regular topology event, not an intrinsic horizon or capture surface.
// -----------------------------------------------------------------------------
TEST(MorrisThorneCartesianTests, ThroatIsARegularTopologyBoundaryNotACaptureSurface) {
    for (double b0 : {1.0, 2.5}) {
        MorrisThorneCartesian cart(MorrisThorneParams::Ellis(b0));

        EXPECT_DOUBLE_EQ(cart.IsotropicThroatRadius(), 0.5 * b0);
        ASSERT_TRUE(cart.IsotropicEllisThroatRadius().has_value());
        EXPECT_DOUBLE_EQ(*cart.IsotropicEllisThroatRadius(), 0.5 * b0);
        EXPECT_FALSE(cart.InsideCaptureSurface(CartPos(0.5 * b0, {1, 0, 0}), 0.0));
        EXPECT_FALSE(cart.InsideCaptureSurface(CartPos(0.25 * b0, {0, 0, 1}), 0.05));
    }
}

TEST(MorrisThorneCartesianTests, InversionExchangesEndsAndPreservesArealRadius) {
    for (double b0 : {0.1, 1.0, 1000.0}) {
        for (double rho : {0.125 * b0, 0.5 * b0, 8.0 * b0, 200.0 * b0}) {
            const double inverted = EllisInvertedIsotropicRadius(b0, rho);
            EXPECT_NEAR(EllisInvertedIsotropicRadius(b0, inverted), rho, 2e-13 * rho);
            EXPECT_NEAR(EllisArealRadiusFromIsotropic(b0, inverted),
                        EllisArealRadiusFromIsotropic(b0, rho),
                        2e-13 * EllisArealRadiusFromIsotropic(b0, rho));
            EXPECT_NEAR(EllisProperRadialDistanceFromIsotropic(b0, inverted),
                        -EllisProperRadialDistanceFromIsotropic(b0, rho), 2e-13 * b0);
        }
    }
}

TEST(MorrisThorneCartesianTests, SecondSheetSkyUsesTheExactInversionJacobian) {
    Vec4 position = CartPos(0.25, {1.0, 0.0, 0.0});
    Vec4 radial;
    radial(1) = -1.0;
    const auto mapped_radial = MapEllisSecondSheetSkyDirection(position, radial);
    ASSERT_TRUE(mapped_radial.has_value());
    EXPECT_DOUBLE_EQ((*mapped_radial)(1), 1.0);

    Vec4 tangent;
    tangent(2) = 1.0;
    const auto mapped_tangent = MapEllisSecondSheetSkyDirection(position, tangent);
    ASSERT_TRUE(mapped_tangent.has_value());
    EXPECT_DOUBLE_EQ((*mapped_tangent)(1), 0.0);
    EXPECT_DOUBLE_EQ((*mapped_tangent)(2), 1.0);

    Vec4 origin;
    EXPECT_FALSE(MapEllisSecondSheetSkyDirection(origin, tangent).has_value());
}

TEST(MorrisThorneCartesianTests, ThroatAndPositiveRadiusSecondSheetAreFiniteAndUnclamped) {
    MorrisThorneCartesian cart(MorrisThorneParams::Ellis(2.0));
    for (double rho : {0.25, 0.5, 1.0, 2.0}) {
        Vec4 position = CartPos(rho, {0.6, 0.8, 0.0});
        EXPECT_TRUE(cart.IsValidEvent(position));
        Metric4d metric;
        Tensor<Dual<double>, 4, 4, 4> derivative;
        cart.Evaluate(position, metric, derivative);
        const double q = 1.0 / (rho * rho);  // b0^2/(4 rho^2), b0=2.
        const double expected = (1.0 + q) * (1.0 + q);
        for (int axis = 1; axis < 4; ++axis) {
            EXPECT_NEAR(metric(axis, axis).real, expected, 1e-13 * expected);
        }
        EXPECT_TRUE(metric_validation::CheckLorentzianSignature(metric));
    }

    Vec4 origin;
    EXPECT_FALSE(cart.IsValidEvent(origin));
    Vec4 nonfinite = CartPos(1.0, {1, 0, 0});
    nonfinite(2) = std::numeric_limits<double>::infinity();
    EXPECT_FALSE(cart.IsValidEvent(nonfinite));
}

TEST(MorrisThorneCartesianTests, SphericalArealChartDeclinesItsCoordinateSingularities) {
    MorrisThorneFamily spherical(MorrisThorneParams::Ellis(1.0));
    Vec4 exterior;
    exterior(1) = 1.01;
    exterior(2) = std::numbers::pi / 2.0;
    EXPECT_TRUE(spherical.IsValidEvent(exterior));

    Vec4 throat = exterior;
    throat(1) = 1.0;
    EXPECT_FALSE(spherical.IsValidEvent(throat));
    Vec4 below = exterior;
    below(1) = 0.99;
    EXPECT_FALSE(spherical.IsValidEvent(below));
    Vec4 pole = exterior;
    pole(2) = 0.0;
    EXPECT_FALSE(spherical.IsValidEvent(pole));
}

TEST(MorrisThorneCartesianTests, AcceptedSegmentFindsHiddenAndTangentThroatContacts) {
    Vec4 start, finish, start_tangent, finish_tangent;
    start(1) = 1.1;
    finish(0) = 1.0;
    finish(1) = 1.1;
    start_tangent(0) = 1.0;
    finish_tangent(0) = 1.0;
    start_tangent(1) = -0.8;
    finish_tangent(1) = 0.8;

    const auto hidden =
        FindSphericalCaptureEvent(start, start_tangent, finish, finish_tangent, 1.0, 1.0);
    ASSERT_TRUE(hidden.has_value());
    const double expected_fraction = (0.8 - std::sqrt(0.32)) / 1.6;
    EXPECT_NEAR(hidden->fraction, expected_fraction, 2e-14);
    EXPECT_NEAR(hidden->position(1), 1.0, 2e-14);
    EXPECT_NEAR(hidden->tangent(1), -std::sqrt(0.32), 2e-13);
    EXPECT_GT(std::abs(start(1)), 1.0);
    EXPECT_GT(std::abs(finish(1)), 1.0);
    const auto midpoint =
        SampleAcceptedTraceSegment(start, start_tangent, finish, finish_tangent, 1.0, 0.5);
    EXPECT_NEAR(midpoint.position(1), 0.9, 1e-15);

    start_tangent(1) = -0.4;
    finish_tangent(1) = 0.4;
    const auto tangent =
        FindSphericalCaptureEvent(start, start_tangent, finish, finish_tangent, 1.0, 1.0);
    ASSERT_TRUE(tangent.has_value());
    EXPECT_NEAR(tangent->fraction, 0.5, 2e-12);
    EXPECT_NEAR(tangent->position(1), 1.0, 2e-13);
    EXPECT_NEAR(tangent->tangent(1), 0.0, 2e-12);

    start_tangent(1) = -0.3;
    finish_tangent(1) = 0.3;
    EXPECT_FALSE(FindSphericalCaptureEvent(start, start_tangent, finish, finish_tangent, 1.0, 1.0)
                     .has_value());
}

TEST(MorrisThorneCartesianTests, DirectionalBoundaryIgnoresTangentAndSelectsInwardCrossing) {
    Vec4 start, finish, start_tangent, finish_tangent;
    start(1) = 1.1;
    finish(0) = 1.0;
    finish(1) = 1.1;
    start_tangent(0) = 1.0;
    finish_tangent(0) = 1.0;
    start_tangent(1) = -0.8;
    finish_tangent(1) = 0.8;
    const auto inward =
        FindSphericalBoundaryEvent(start, start_tangent, finish, finish_tangent, 1.0, 1.0,
                                   SphericalBoundarySense::DecreasingRadius);
    ASSERT_TRUE(inward.has_value());
    EXPECT_LT(inward->tangent(1), 0.0);

    start_tangent(1) = -0.4;
    finish_tangent(1) = 0.4;
    EXPECT_FALSE(FindSphericalBoundaryEvent(start, start_tangent, finish, finish_tangent, 1.0, 1.0,
                                            SphericalBoundarySense::DecreasingRadius)
                     .has_value());
    EXPECT_TRUE(FindSphericalCaptureEvent(start, start_tangent, finish, finish_tangent, 1.0, 1.0)
                    .has_value());
}

TEST(MorrisThorneCartesianTests, NonEllisCartesianRequestsFailClosed) {
    EXPECT_DEATH(
        {
            MorrisThorneCartesian rejected(MorrisThorneParams::ZeroTidal(1.0));
            (void)rejected;
        },
        "violated");
}

TEST(MorrisThorneCartesianTests, ParameterBoundsMatchConfigAuthority) {
    auto boundary = MorrisThorneParams::Ellis(kMinMorrisThorneThroatRadius);
    boundary.Phi0 = 10.0;
    MorrisThorneFamily accepted(boundary);
    EXPECT_DOUBLE_EQ(accepted.GetParams().b0, kMinMorrisThorneThroatRadius);
    EXPECT_DOUBLE_EQ(accepted.GetParams().Phi0, 10.0);

    EXPECT_DEATH(
        {
            auto outside = MorrisThorneParams::Ellis(1.0);
            outside.Phi0 = 10.01;
            MorrisThorneFamily rejected(outside);
            (void)rejected;
        },
        "violated");
    EXPECT_DEATH(
        {
            MorrisThorneFamily rejected(
                MorrisThorneParams::Ellis(kMinMorrisThorneThroatRadius - 0.01));
            (void)rejected;
        },
        "violated");
}

}  // namespace
