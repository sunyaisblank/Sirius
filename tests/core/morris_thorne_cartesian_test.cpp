// MorrisThorneCartesian: the flat-plus-rank-one Cartesian embedding of the
// spherical Morris-Thorne family, the chart the CPU tracer drives.
//
// Gates: agreement with the spherical authority component by component
// (radial, tangential, cross, time), analytic derivatives against finite
// differences of the metric itself, exactness of the Sherman-Morrison inverse,
// and the throat capture surface. Together these are the preconditions the
// tracer's geodesic step relies on; the spherical family stays the
// shape-function authority so a defect here cannot silently fork the physics.

#include "sirius/core/metrics/morris_thorne_family.h"
#include "sirius/core/tensor.h"

#include <gtest/gtest.h>

#include <cmath>
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
        MorrisThorneParams::ZeroTidal(1.0),
        {1.0, 0.0, WormholeShapeType::AbsurdlyBenign, nullptr, nullptr},
    };
}

// -----------------------------------------------------------------------------
// The embedding reproduces the spherical chart's line element: the radial
// sectional component is 1/(1 - b/r), any tangential unit direction has unit
// length, radial and tangential directions are orthogonal, and g_tt matches.
// -----------------------------------------------------------------------------
TEST(MorrisThorneCartesianTests, ChartAgreementWithSphericalFamily) {
    for (const auto& params : SampleShapes()) {
        MorrisThorneCartesian cart(params);
        MorrisThorneFamily sph(params);
        double b0 = params.b0;

        for (double rf : {1.5, 2.0, 5.0, 20.0}) {
            double r = rf * b0;
            for (const auto& n : SampleDirections()) {
                Metric4d g;
                Tensor<Dual<double>, 4, 4, 4> dg;
                cart.Evaluate(CartPos(r, n), g, dg);

                Metric4d gs;
                Tensor<Dual<double>, 4, 4, 4> dgs;
                Vec4 spos;
                spos(1) = r;
                spos(2) = M_PI / 2.0;
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

                EXPECT_NEAR(radial, gs(1, 1).real, 1e-12 * std::abs(gs(1, 1).real))
                    << "radial component, b0=" << b0 << " r=" << r;
                EXPECT_NEAR(tangential, 1.0, 1e-12)
                    << "tangential component, b0=" << b0 << " r=" << r;
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
// rule through f(r) and the radial frame.
// -----------------------------------------------------------------------------
TEST(MorrisThorneCartesianTests, DerivativesMatchFiniteDifferencesOfMetric) {
    const double h = 1e-6;
    for (const auto& params : SampleShapes()) {
        MorrisThorneCartesian cart(params);
        double b0 = params.b0;

        for (double rf : {1.5, 3.0, 10.0}) {
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
// The closed-form inverse satisfies g^mu_alpha g_alpha_nu = delta^mu_nu to
// near machine precision; Sherman-Morrison is exact, so the bar is rounding.
// -----------------------------------------------------------------------------
TEST(MorrisThorneCartesianTests, AnalyticInverseIsExact) {
    for (const auto& params : SampleShapes()) {
        MorrisThorneCartesian cart(params);
        double b0 = params.b0;

        for (double rf : {1.2, 1.5, 3.0, 50.0}) {
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
// The throat is the capture surface: inside at r <= b0 (1 + margin), outside
// beyond it. The tracer's horizon_factor drives margin exactly as for black
// holes, so the same termination path serves the one-sheet wormhole.
// -----------------------------------------------------------------------------
TEST(MorrisThorneCartesianTests, ThroatIsTheCaptureSurface) {
    for (double b0 : {1.0, 2.5}) {
        MorrisThorneCartesian cart(MorrisThorneParams::Ellis(b0));
        const double margin = 0.05;

        EXPECT_TRUE(cart.InsideCaptureSurface(CartPos(1.049 * b0, {1, 0, 0}), margin));
        EXPECT_TRUE(cart.InsideCaptureSurface(CartPos(0.5 * b0, {0, 0, 1}), margin));
        EXPECT_FALSE(cart.InsideCaptureSurface(CartPos(1.051 * b0, {0.6, 0.8, 0}), margin));
        EXPECT_FALSE(cart.InsideCaptureSurface(CartPos(10.0 * b0, {1, 0, 0}), margin));
    }
}

}  // namespace
