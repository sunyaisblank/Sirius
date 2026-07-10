// TSGD003A.cpp - Live-Path Conservation and Capture-Surface Tests
// Component ID: TSGD003A
// Tests for: PHGD001A (integrateStepRK45, checkTermination),
//            PHMT100A (insideCaptureSurface), PHMT102A
//
// The conservation suite historically validated the double-precision
// Boyer-Lindquist oracle stack, which is not the code that renders. These
// tests drive the LIVE path: the single-precision-configured Cartesian
// Kerr-Schild family through Geodesic::integrateStepRK45, monitoring the
// conserved quantities of the stationary axisymmetric geometry:
//   energy            E  = -g_{t mu} k^mu          (Killing vector d/dt)
//   angular momentum  L_z = g_{mu nu} xi^mu k^nu,  xi = (0, -y, x, 0)
// and the null condition g_{mu nu} k^mu k^nu = 0. Tolerances follow
// docs-era spec values: conservation drift < 1e-4 relative, null < 1e-6.

#include <gtest/gtest.h>
#include <PHGD001A.h>
#include <PHMT100A.h>
#include <PHMT102A.h>
#include <PHCN001A.h>
#include <cmath>

using namespace Sirius;

namespace {

double energyOf(const Lightray& ray, const Metric4D& g) {
    double E = 0.0;
    for (int nu = 0; nu < 4; nu++) {
        E -= g(0, nu).real * ray.velocity(nu);
    }
    return E;
}

double angularMomentumOf(const Lightray& ray, const Metric4D& g) {
    // xi = d/dphi in Cartesian components: (0, -y, x, 0)
    Vec4 xi;
    xi(1) = -ray.position(2);
    xi(2) = ray.position(1);
    return TensorOps::innerProduct(xi, ray.velocity, g);
}

Lightray makeRay(IMetric& metric, double x, double y, double z,
                 double vx, double vy, double vz) {
    Lightray ray{};
    ray.position(0) = 0.0;
    ray.position(1) = x;
    ray.position(2) = y;
    ray.position(3) = z;
    ray.velocity(1) = vx;
    ray.velocity(2) = vy;
    ray.velocity(3) = vz;

    Metric4D g;
    Tensor<Dual<double>, 4, 4, 4> dg;
    metric.evaluate(ray.position, g, dg);
    ray.velocity = TensorOps::normalizeNull(ray.velocity, g);

    ray.step_size = 0.01f;
    ray.proper_time = 0.0f;
    ray.terminated = 0;
    return ray;
}

} // namespace

class LivePathConservationTests : public ::testing::Test {
protected:
    IntegratorConfig config = Geodesic::getDefaultConfig();

    // Integrate up to maxSteps accepted steps or until termination; returns
    // the worst relative drift of E and L_z and the worst null violation.
    void run(IMetric& metric, Lightray ray, int maxSteps,
             double& worstEDrift, double& worstLDrift, double& worstNull) {
        Metric4D g;
        Tensor<Dual<double>, 4, 4, 4> dg;
        metric.evaluate(ray.position, g, dg);
        const double E0 = energyOf(ray, g);
        const double L0 = angularMomentumOf(ray, g);
        ASSERT_GT(std::abs(E0), 1e-12) << "degenerate initial energy";

        worstEDrift = worstLDrift = worstNull = 0.0;
        int accepted = 0, attempts = 0;
        while (accepted < maxSteps && attempts < maxSteps * 10) {
            ++attempts;
            if (!Geodesic::integrateStepRK45(ray, &metric, config)) {
                if (ray.terminated != 0) break;
                continue;  // step rejected, retry with adapted step
            }
            ++accepted;
            if (Geodesic::checkTermination(ray, &metric)) break;

            metric.evaluate(ray.position, g, dg);
            worstEDrift = std::max(worstEDrift,
                                   std::abs(energyOf(ray, g) - E0) / std::abs(E0));
            if (std::abs(L0) > 1e-6) {
                worstLDrift = std::max(worstLDrift,
                                       std::abs(angularMomentumOf(ray, g) - L0) / std::abs(L0));
            }
            worstNull = std::max(worstNull,
                                 std::abs(TensorOps::innerProduct(ray.velocity,
                                                                  ray.velocity, g)));
        }
        ASSERT_GT(accepted, 100) << "integrator made too little progress";
    }
};

TEST_F(LivePathConservationTests, SchwarzschildEnergyAndAngularMomentum) {
    KerrSchildFamily metric(KerrSchildParams::Schwarzschild(1.0));
    // Tangential launch at r = 10M in the equatorial plane
    Lightray ray = makeRay(metric, 10.0, 0.0, 0.0, 0.0, 1.0, 0.0);

    double eDrift, lDrift, nullViolation;
    run(metric, ray, 1500, eDrift, lDrift, nullViolation);

    EXPECT_LT(eDrift, Constants::Geodesic::CONSERVATION_TOL);
    EXPECT_LT(lDrift, Constants::Geodesic::CONSERVATION_TOL);
    EXPECT_LT(nullViolation, 1e-6);
}

TEST_F(LivePathConservationTests, KerrEnergyAndAngularMomentum) {
    KerrSchildFamily metric(KerrSchildParams::Kerr(1.0, 0.9));
    // Off-plane launch so every metric component participates
    Lightray ray = makeRay(metric, 8.0, 0.0, 2.0, 0.1, 1.0, 0.0);

    double eDrift, lDrift, nullViolation;
    run(metric, ray, 1500, eDrift, lDrift, nullViolation);

    EXPECT_LT(eDrift, Constants::Geodesic::CONSERVATION_TOL);
    EXPECT_LT(lDrift, Constants::Geodesic::CONSERVATION_TOL);
    EXPECT_LT(nullViolation, 1e-6);
}

//==============================================================================
// Capture surface and termination
//==============================================================================

TEST(CaptureSurfaceTests, SchwarzschildHorizonIsCaptured) {
    KerrSchildFamily metric(KerrSchildParams::Schwarzschild(1.0));
    Tensor<double, 4> inside;
    inside(1) = 1.5;  // r = 1.5M < 2M
    EXPECT_TRUE(metric.insideCaptureSurface(inside, 0.05));

    Tensor<double, 4> outside;
    outside(1) = 6.0;
    EXPECT_FALSE(metric.insideCaptureSurface(outside, 0.05));
}

TEST(CaptureSurfaceTests, KerrOblateHorizonUsesKerrRadius) {
    // At a = 0.9 the horizon r+ = 1 + sqrt(1 - 0.81) is an oblate surface:
    // an equatorial point can sit inside it (Kerr r < r+) while its
    // Cartesian norm exceeds r+; a Cartesian-norm comparison misses this.
    KerrSchildFamily metric(KerrSchildParams::Kerr(1.0, 0.9));
    const double rPlus = metric.outerHorizonRadius();
    ASSERT_NEAR(rPlus, 1.0 + std::sqrt(1.0 - 0.81), 1e-12);

    const double rInside = 0.98 * rPlus;
    const double a = 0.9;
    // Equatorial (z = 0): Cartesian norm R = sqrt(r^2 + a^2) > r+
    const double R = std::sqrt(rInside * rInside + a * a);
    ASSERT_GT(R, rPlus * 1.05);

    Tensor<double, 4> pos;
    pos(1) = R;  // place along x-axis: Kerr radius there is rInside
    EXPECT_NEAR(metric.computeKerrRadius(R, 0.0, 0.0), rInside, 1e-9);
    EXPECT_TRUE(metric.insideCaptureSurface(pos, 0.05));
}

TEST(CaptureSurfaceTests, HorizonlessSpacetimesNeverCapture) {
    WarpDriveFamily warp(WarpDriveParams::Alcubierre(1.0, 1.0));
    Tensor<double, 4> pos;
    pos(1) = 0.1;  // even deep inside the bubble
    EXPECT_FALSE(warp.insideCaptureSurface(pos, 0.05));

    KerrSchildFamily naked(KerrSchildParams::KerrNewman(1.0, 0.9, 0.9));
    ASSERT_FALSE(naked.hasHorizon());
    Tensor<double, 4> origin;
    origin(1) = 0.5;
    EXPECT_FALSE(naked.insideCaptureSurface(origin, 0.05));
}

TEST(CaptureSurfaceTests, CheckTerminationUsesCartesianNormAndCapture) {
    KerrSchildFamily metric(KerrSchildParams::Schwarzschild(1.0));

    Lightray escaping = makeRay(metric, 60.0, 0.0, 0.0, 1.0, 0.0, 0.0);
    EXPECT_TRUE(Geodesic::checkTermination(escaping, &metric))
        << "outward ray beyond the escape radius must terminate";

    Lightray background = makeRay(metric, 0.0, 150.0, 0.0, 0.0, -1.0, 0.0);
    EXPECT_TRUE(Geodesic::checkTermination(background, &metric))
        << "any ray beyond the background radius must terminate";

    Lightray captured = makeRay(metric, 1.2, 0.0, 0.0, 0.0, 1.0, 0.0);
    EXPECT_TRUE(Geodesic::checkTermination(captured, &metric))
        << "ray inside the horizon must terminate as captured";

    // A ray with |y| below the old THETA_MIN = 0.01 must NOT terminate: the
    // deleted pole check misread Cartesian y as a polar angle and killed
    // half of every image.
    Lightray equatorial = makeRay(metric, 10.0, 0.001, 0.0, 0.0, 1.0, 0.0);
    EXPECT_FALSE(Geodesic::checkTermination(equatorial, &metric));
}
