// MorrisThorneCartesian: the exact isotropic Cartesian chart of the zero-tidal
// Ellis member, the chart the CPU tracer drives.
//
// Gates: coordinate-transformed agreement with the spherical areal-radius
// authority on both sheets, a finite exact throat, analytic derivatives,
// inverse identity, independent Ricci/energy-condition invariants, chart
// domain, and complete accepted-segment capture events.  The latter includes
// equal-sign endpoints and a tangent contact.

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

struct RicciSample {
    double covariant[4][4] = {};
    double scalar = 0.0;
};

ChristoffelSymbols ConnectionAt(MorrisThorneCartesian& metric, const Vec4& position) {
    Metric4d g;
    Tensor<Dual<double>, 4, 4, 4> dg;
    metric.Evaluate(position, g, dg);
    return TensorOps::Christoffel(g, dg);
}

// Independent numerical curvature oracle.  It differentiates the production
// analytic connection with a fourth-order centred stencil and assembles
//
//   R^rho_{ sigma mu nu} = d_mu Gamma^rho_{nu sigma}
//                         - d_nu Gamma^rho_{mu sigma}
//                         + Gamma^rho_{mu lambda} Gamma^lambda_{nu sigma}
//                         - Gamma^rho_{nu lambda} Gamma^lambda_{mu sigma}.
//
// This does not reuse the closed-form Ellis Ricci scalar or stress tensor that
// the test compares against, so a mutually wrong chart/connection pair cannot
// satisfy the oracle by construction.
RicciSample RicciFromConnectionFiniteDifference(MorrisThorneCartesian& metric, const Vec4& position,
                                                double h) {
    const ChristoffelSymbols gamma = ConnectionAt(metric, position);
    double d_gamma[4][4][4][4] = {};
    for (int derivative = 0; derivative < 4; ++derivative) {
        Vec4 plus_one = position;
        Vec4 minus_one = position;
        Vec4 plus_two = position;
        Vec4 minus_two = position;
        plus_one(derivative) += h;
        minus_one(derivative) -= h;
        plus_two(derivative) += 2.0 * h;
        minus_two(derivative) -= 2.0 * h;
        const ChristoffelSymbols gp1 = ConnectionAt(metric, plus_one);
        const ChristoffelSymbols gm1 = ConnectionAt(metric, minus_one);
        const ChristoffelSymbols gp2 = ConnectionAt(metric, plus_two);
        const ChristoffelSymbols gm2 = ConnectionAt(metric, minus_two);
        for (int upper = 0; upper < 4; ++upper) {
            for (int lower_one = 0; lower_one < 4; ++lower_one) {
                for (int lower_two = 0; lower_two < 4; ++lower_two) {
                    d_gamma[derivative][upper][lower_one][lower_two] =
                        (-gp2.gamma(upper, lower_one, lower_two).real +
                         8.0 * gp1.gamma(upper, lower_one, lower_two).real -
                         8.0 * gm1.gamma(upper, lower_one, lower_two).real +
                         gm2.gamma(upper, lower_one, lower_two).real) /
                        (12.0 * h);
                }
            }
        }
    }

    RicciSample result;
    for (int sigma = 0; sigma < 4; ++sigma) {
        for (int nu = 0; nu < 4; ++nu) {
            for (int rho = 0; rho < 4; ++rho) {
                double component = d_gamma[rho][rho][nu][sigma] - d_gamma[nu][rho][rho][sigma];
                for (int lambda = 0; lambda < 4; ++lambda) {
                    component +=
                        gamma.gamma(rho, rho, lambda).real * gamma.gamma(lambda, nu, sigma).real -
                        gamma.gamma(rho, nu, lambda).real * gamma.gamma(lambda, rho, sigma).real;
                }
                result.covariant[sigma][nu] += component;
            }
        }
    }

    Metric4d inverse;
    EXPECT_TRUE(metric.InverseMetric(position, inverse));
    for (int mu = 0; mu < 4; ++mu) {
        for (int nu = 0; nu < 4; ++nu) {
            result.scalar += inverse(mu, nu).real * result.covariant[mu][nu];
        }
    }
    return result;
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
// In proper radial distance l the exact Ellis metric is
//
//   ds^2 = -dt^2 + dl^2 + (l^2+b0^2) dOmega^2.
//
// Therefore R = -2 b0^2/r_areal^4.  The same value is R_mu_nu k^mu k^nu for
// either future radial null vector normalised to k^t=1, proving strict radial
// null-energy-condition violation through Einstein's equation.  The identity
// is invariant under the isotropic chart and under exchange of its two ends.
// -----------------------------------------------------------------------------
TEST(MorrisThorneCartesianTests, CurvatureAndRadialNullEnergyConditionMatchExactEllis) {
    constexpr Direction radial{0.267261241912424, -0.534522483824849, 0.801783725737273};
    const Direction angular = Orthogonal(radial);

    for (const double b0 : {0.25, 1.0, 10.0}) {
        MorrisThorneCartesian metric(MorrisThorneParams::Ellis(b0));
        for (const double rho_factor : {0.25, 0.5, 1.0, 3.0, 10.0}) {
            const double rho = rho_factor * b0;
            const Vec4 position = CartPos(rho, radial);
            const RicciSample ricci =
                RicciFromConnectionFiniteDifference(metric, position, 5.0e-5 * b0);
            const double areal_radius = EllisArealRadiusFromIsotropic(b0, rho);
            const double expected = -2.0 * b0 * b0 / std::pow(areal_radius, 4.0);
            const double absolute_floor = 2.0e-8 / (b0 * b0);
            const double tolerance = std::max(absolute_floor, 2.0e-5 * std::abs(expected));

            EXPECT_NEAR(ricci.scalar, expected, tolerance)
                << "Ricci scalar, b0=" << b0 << " rho=" << rho;

            const double q = b0 * b0 / (4.0 * rho * rho);
            const double inverse_conformal_root = 1.0 / (1.0 + q);
            const double radial_null[4] = {
                1.0,
                radial.x * inverse_conformal_root,
                radial.y * inverse_conformal_root,
                radial.z * inverse_conformal_root,
            };
            Metric4d metric_tensor;
            Tensor<Dual<double>, 4, 4, 4> metric_derivative;
            metric.Evaluate(position, metric_tensor, metric_derivative);
            double null_norm = 0.0;
            for (int mu = 0; mu < 4; ++mu) {
                for (int nu = 0; nu < 4; ++nu) {
                    null_norm += metric_tensor(mu, nu).real * radial_null[mu] * radial_null[nu];
                }
            }
            EXPECT_NEAR(null_norm, 0.0, 2.0e-14)
                << "radial vector must be null, b0=" << b0 << " rho=" << rho;

            double radial_null_ricci = 0.0;
            for (int mu = 0; mu < 4; ++mu) {
                for (int nu = 0; nu < 4; ++nu) {
                    radial_null_ricci +=
                        ricci.covariant[mu][nu] * radial_null[mu] * radial_null[nu];
                }
            }
            EXPECT_NEAR(radial_null_ricci, expected, tolerance)
                << "radial null Ricci contraction, b0=" << b0 << " rho=" << rho;
            EXPECT_LT(radial_null_ricci, 0.0)
                << "the represented traversable throat must require radial NEC violation";

            double angular_ricci = 0.0;
            const double angular_unit[3] = {
                angular.x * inverse_conformal_root,
                angular.y * inverse_conformal_root,
                angular.z * inverse_conformal_root,
            };
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    angular_ricci +=
                        ricci.covariant[i + 1][j + 1] * angular_unit[i] * angular_unit[j];
                }
            }
            EXPECT_NEAR(angular_ricci, 0.0, tolerance)
                << "angular Ricci eigenvalue, b0=" << b0 << " rho=" << rho;
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

TEST(MorrisThorneCartesianTests, LiveEllisRejectsUnnormalisedLapseAndIrrelevantShapeData) {
    EXPECT_DEATH(
        {
            auto unnormalised = MorrisThorneParams::Ellis(1.0);
            unnormalised.Phi0 = 0.25;
            MorrisThorneCartesian rejected(unnormalised);
            (void)rejected;
        },
        "violated");

    EXPECT_DEATH(
        {
            auto irrelevant_callback = MorrisThorneParams::Ellis(1.0);
            irrelevant_callback.custom_shape_func = [](double radius) { return radius; };
            MorrisThorneCartesian rejected(irrelevant_callback);
            (void)rejected;
        },
        "violated");

    MorrisThorneCartesian represented;
    EXPECT_EQ(represented.GetParameters().size(), 1U);
    EXPECT_TRUE(represented.GetParameters().contains("throat_radius"));
    EXPECT_FALSE(represented.GetParameters().contains("redshift"));
    EXPECT_DEATH(represented.SetParameter("redshift", 0.0), "violated");
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
