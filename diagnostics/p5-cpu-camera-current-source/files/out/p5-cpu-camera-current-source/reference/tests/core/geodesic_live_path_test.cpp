// Live-Path Conservation and Capture-Surface Tests
// Ported from TSGD003A.cpp; assertions and tolerances unchanged.
// Tests for: geodesic_integrator (IntegrateStepRk45, CheckTermination),
//            kerr_schild_family (InsideCaptureSurface), warp_drive_family
//
// The conservation suite historically validated the double-precision
// Boyer-Lindquist oracle stack, which is not the code that renders. These
// tests drive the LIVE path: the single-precision-configured Cartesian
// Kerr-Schild family through Geodesic::IntegrateStepRk45, monitoring the
// conserved quantities of the stationary axisymmetric geometry:
//   energy            E  = -g_{t mu} k^mu          (Killing vector d/dt)
//   angular momentum  L_z = g_{mu nu} xi^mu k^nu,  xi = (0, -y, x, 0)
// and the null condition g_{mu nu} k^mu k^nu = 0. Tolerances follow
// docs-era spec values: conservation drift < 1e-4 relative, null < 1e-6.

#include "sirius/core/constants.h"
#include "sirius/core/coordinates.h"
#include "sirius/core/geodesic_integrator.h"
#include "sirius/core/metrics/kerr_schild_family.h"
#include "sirius/core/metrics/warp_drive_family.h"
#include "sirius/core/observer_frame.h"

#include <gtest/gtest.h>

#include <array>
#include <cmath>

using namespace sirius::core;

namespace {

double energyOf(const Lightray& ray, const Metric4d& g) {
    double E = 0.0;
    for (int nu = 0; nu < 4; nu++) {
        E -= g(0, nu).real * ray.velocity(nu);
    }
    return E;
}

double angularMomentumOf(const Lightray& ray, const Metric4d& g) {
    // xi = d/dphi in Cartesian components: (0, -y, x, 0)
    Vec4 xi;
    xi(1) = -ray.position(2);
    xi(2) = ray.position(1);
    return TensorOps::InnerProduct(xi, ray.velocity, g);
}

double carterConstantOf(const Lightray& ray, const Metric4d& g, double mass, double spin) {
    const coordinates::Vec4Cart position{ray.position(0), ray.position(1), ray.position(2),
                                         ray.position(3)};
    const coordinates::Vec4Cart velocity{ray.velocity(0), ray.velocity(1), ray.velocity(2),
                                         ray.velocity(3)};
    const coordinates::Vec4Bl position_bl = coordinates::KerrSchildCartToBl(position, spin);
    const coordinates::Vec4Bl velocity_bl =
        coordinates::TransformVectorKerrSchildCartToBl(velocity, position, mass, spin);
    const double energy = energyOf(ray, g);
    const double angular_momentum = angularMomentumOf(ray, g);
    const double cosine = std::cos(position_bl.theta);
    const double sine = std::sin(position_bl.theta);
    const double sigma = position_bl.r * position_bl.r + spin * spin * cosine * cosine;
    const double polar_momentum = sigma * velocity_bl.theta;
    return polar_momentum * polar_momentum +
           cosine * cosine *
               (-spin * spin * energy * energy +
                angular_momentum * angular_momentum / (sine * sine));
}

bool reachedPhysicalBoundary(const Lightray& ray, const IMetric& metric) {
    using namespace constants::termination;
    const double x = ray.position(1);
    const double y = ray.position(2);
    const double z = ray.position(3);
    const double vx = ray.velocity(1);
    const double vy = ray.velocity(2);
    const double vz = ray.velocity(3);
    if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z) || !std::isfinite(vx) ||
        !std::isfinite(vy) || !std::isfinite(vz)) {
        return false;
    }
    const double radius = std::sqrt(x * x + y * y + z * z);
    const double radial_rate = (x * vx + y * vy + z * vz) / std::max(radius, 1.0e-12);
    return (radius > kEscapeRadius && radial_rate > 0.0) || radius > kBackgroundRadius ||
           metric.InsideCaptureSurface(ray.position, kCaptureMargin);
}

Lightray makeRay(IMetric& metric, double x, double y, double z, double vx, double vy, double vz) {
    Lightray ray{};
    ray.position(0) = 0.0;
    ray.position(1) = x;
    ray.position(2) = y;
    ray.position(3) = z;
    ray.velocity(1) = vx;
    ray.velocity(2) = vy;
    ray.velocity(3) = vz;

    Metric4d g;
    Tensor<Dual<double>, 4, 4, 4> dg;
    metric.Evaluate(ray.position, g, dg);
    ray.velocity = TensorOps::NormalizeNull(ray.velocity, g);

    ray.step_size = 0.01f;
    ray.proper_time = 0.0f;
    ray.terminated = 0;
    return ray;
}

}  // namespace

class LivePathConservationTests : public ::testing::Test {
  protected:
    IntegratorConfig config = Geodesic::GetDefaultConfig();

    // Integrate up to max_steps accepted steps or until termination; returns
    // the worst relative drift of E and L_z and the worst null violation. Kerr
    // callers may also monitor the Carter constant through the exact
    // Kerr-Schild-to-Boyer-Lindquist vector transform.
    void run(IMetric& metric, Lightray ray, int max_steps, double& worstEDrift, double& worstLDrift,
             double& worstNull, double* worstCarterDrift = nullptr, double mass = 0.0,
             double spin = 0.0, bool* reached_physical_boundary = nullptr) {
        Metric4d g;
        Tensor<Dual<double>, 4, 4, 4> dg;
        metric.Evaluate(ray.position, g, dg);
        const double E0 = energyOf(ray, g);
        const double L0 = angularMomentumOf(ray, g);
        const double carter0 =
            worstCarterDrift != nullptr ? carterConstantOf(ray, g, mass, spin) : 0.0;
        ASSERT_GT(std::abs(E0), 1e-12) << "degenerate initial energy";
        if (worstCarterDrift != nullptr) {
            ASSERT_GT(std::abs(carter0), 1e-12) << "degenerate initial Carter constant";
            *worstCarterDrift = 0.0;
        }

        worstEDrift = worstLDrift = worstNull = 0.0;
        if (reached_physical_boundary != nullptr) *reached_physical_boundary = false;
        int accepted = 0, attempts = 0;
        while (accepted < max_steps && attempts < max_steps * 10) {
            ++attempts;
            if (!Geodesic::IntegrateStepRk45(ray, &metric, config)) {
                if (ray.terminated != 0) {
                    if (reached_physical_boundary != nullptr) {
                        *reached_physical_boundary = reachedPhysicalBoundary(ray, metric);
                    }
                    break;
                }
                continue;  // step rejected, retry with adapted step
            }
            ++accepted;
            const bool terminated = Geodesic::CheckTermination(ray, &metric);
            const bool captured =
                metric.InsideCaptureSurface(ray.position, constants::termination::kCaptureMargin);

            // The escaping terminal state is part of the complete-ray witness.
            // A captured state may already be inside the family's safe metric
            // domain, so its immediately preceding exterior sample remains the
            // last conservation sample.
            if (!captured) {
                metric.Evaluate(ray.position, g, dg);
                worstEDrift = std::max(worstEDrift, std::abs(energyOf(ray, g) - E0) / std::abs(E0));
                if (std::abs(L0) > 1e-6) {
                    worstLDrift = std::max(worstLDrift,
                                           std::abs(angularMomentumOf(ray, g) - L0) / std::abs(L0));
                }
                worstNull = std::max(
                    worstNull, std::abs(TensorOps::InnerProduct(ray.velocity, ray.velocity, g)));
                if (worstCarterDrift != nullptr) {
                    *worstCarterDrift =
                        std::max(*worstCarterDrift,
                                 std::abs(carterConstantOf(ray, g, mass, spin) - carter0) /
                                     std::abs(carter0));
                }
            }
            if (terminated) {
                if (reached_physical_boundary != nullptr) {
                    *reached_physical_boundary = reachedPhysicalBoundary(ray, metric);
                }
                break;
            }
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

    EXPECT_LT(eDrift, constants::geodesic::kConservationTol);
    EXPECT_LT(lDrift, constants::geodesic::kConservationTol);
    EXPECT_LT(nullViolation, 1e-6);
}

TEST_F(LivePathConservationTests, KerrEnergyAndAngularMomentum) {
    KerrSchildFamily metric(KerrSchildParams::Kerr(1.0, 0.9));
    // Off-plane launch so every metric component participates
    Lightray ray = makeRay(metric, 8.0, 0.0, 2.0, 0.1, 1.0, 0.0);

    double eDrift, lDrift, nullViolation;
    run(metric, ray, 1500, eDrift, lDrift, nullViolation);

    EXPECT_LT(eDrift, constants::geodesic::kConservationTol);
    EXPECT_LT(lDrift, constants::geodesic::kConservationTol);
    EXPECT_LT(nullViolation, 1e-6);
}

TEST_F(LivePathConservationTests, NearExtremalKerrEnergyAngularMomentumAndCarter) {
    KerrSchildFamily metric(KerrSchildParams::Kerr(1.0, 0.998));
    Lightray ray = makeRay(metric, 12.0, 0.0, 3.0, 0.08, 1.0, -0.03);

    double eDrift, lDrift, nullViolation, carterDrift;
    bool reached_physical_boundary = false;
    run(metric, ray, 2000, eDrift, lDrift, nullViolation, &carterDrift, 1.0, 0.998,
        &reached_physical_boundary);

    EXPECT_TRUE(reached_physical_boundary)
        << "the live witness ended numerically instead of at escape/capture";
    EXPECT_LT(eDrift, constants::geodesic::kConservationTol);
    EXPECT_LT(lDrift, constants::geodesic::kConservationTol);
    EXPECT_LT(carterDrift, constants::geodesic::kConservationTol);
    EXPECT_LT(nullViolation, 1e-6);
}

//==============================================================================
// Capture surface and termination
//==============================================================================

TEST(CaptureSurfaceTests, SchwarzschildHorizonIsCaptured) {
    KerrSchildFamily metric(KerrSchildParams::Schwarzschild(1.0));
    Tensor<double, 4> inside;
    inside(1) = 1.5;  // r = 1.5M < 2M
    EXPECT_TRUE(metric.InsideCaptureSurface(inside, 0.05));

    Tensor<double, 4> outside;
    outside(1) = 6.0;
    EXPECT_FALSE(metric.InsideCaptureSurface(outside, 0.05));
}

TEST(CaptureSurfaceTests, KerrOblateHorizonUsesKerrRadius) {
    // At a = 0.9 the horizon r+ = 1 + sqrt(1 - 0.81) is an oblate surface:
    // an equatorial point can sit inside it (Kerr r < r+) while its
    // Cartesian norm exceeds r+; a Cartesian-norm comparison misses this.
    KerrSchildFamily metric(KerrSchildParams::Kerr(1.0, 0.9));
    const double rPlus = metric.OuterHorizonRadius();
    ASSERT_NEAR(rPlus, 1.0 + std::sqrt(1.0 - 0.81), 1e-12);

    const double rInside = 0.98 * rPlus;
    const double a = 0.9;
    // Equatorial (z = 0): Cartesian norm R = sqrt(r^2 + a^2) > r+
    const double R = std::sqrt(rInside * rInside + a * a);
    ASSERT_GT(R, rPlus * 1.05);

    Tensor<double, 4> pos;
    pos(1) = R;  // place along x-axis: Kerr radius there is rInside
    EXPECT_NEAR(metric.ComputeKerrRadius(R, 0.0, 0.0), rInside, 1e-9);
    EXPECT_TRUE(metric.InsideCaptureSurface(pos, 0.05));
}

TEST(CaptureSurfaceTests, HorizonlessSpacetimesNeverCapture) {
    WarpDriveFamily warp(WarpDriveParams::Alcubierre(1.0, 1.0));
    Tensor<double, 4> pos;
    pos(1) = 0.1;  // even deep inside the bubble
    EXPECT_FALSE(warp.InsideCaptureSurface(pos, 0.05));

    KerrSchildFamily naked(KerrSchildParams::KerrNewman(1.0, 0.9, 0.9));
    ASSERT_FALSE(naked.HasHorizon());
    Tensor<double, 4> origin;
    origin(1) = 0.5;
    EXPECT_FALSE(naked.InsideCaptureSurface(origin, 0.05));
}

TEST(CaptureSurfaceTests, CheckTerminationUsesCartesianNormAndCapture) {
    KerrSchildFamily metric(KerrSchildParams::Schwarzschild(1.0));

    Lightray escaping = makeRay(metric, 60.0, 0.0, 0.0, 1.0, 0.0, 0.0);
    EXPECT_TRUE(Geodesic::CheckTermination(escaping, &metric))
        << "outward ray beyond the escape radius must terminate";

    Lightray background = makeRay(metric, 0.0, 150.0, 0.0, 0.0, -1.0, 0.0);
    EXPECT_TRUE(Geodesic::CheckTermination(background, &metric))
        << "any ray beyond the background radius must terminate";

    Lightray captured = makeRay(metric, 1.2, 0.0, 0.0, 0.0, 1.0, 0.0);
    EXPECT_TRUE(Geodesic::CheckTermination(captured, &metric))
        << "ray inside the horizon must terminate as captured";

    // A ray with |y| below the old THETA_MIN = 0.01 must NOT terminate: the
    // deleted pole check misread Cartesian y as a polar angle and killed
    // half of every image.
    Lightray equatorial = makeRay(metric, 10.0, 0.001, 0.0, 0.0, 1.0, 0.0);
    EXPECT_FALSE(Geodesic::CheckTermination(equatorial, &metric));
}

TEST(ObserverFrameTests, MovingKerrSchildObserverProducesAnOrthonormalTetrad) {
    KerrSchildFamily metric(KerrSchildParams::Kerr(1.0, 0.7));
    Vec4 position;
    position(1) = 12.0;
    position(2) = 1.5;
    position(3) = 0.8;
    Vec4 coordinate_velocity;
    coordinate_velocity(0) = 1.0;
    coordinate_velocity(1) = 0.08;
    coordinate_velocity(2) = -0.03;
    coordinate_velocity(3) = 0.02;

    const ObserverState observer = Geodesic::CreateObserver(position, coordinate_velocity, &metric);
    Metric4d g;
    Tensor<Dual<double>, 4, 4, 4> dg;
    metric.Evaluate(position, g, dg);
    const std::array<Vec4, 4> basis = {observer.e0, observer.e1, observer.e2, observer.e3};
    for (int row = 0; row < 4; ++row) {
        for (int column = 0; column < 4; ++column) {
            const double expected = (row == column) ? (row == 0 ? -1.0 : 1.0) : 0.0;
            EXPECT_NEAR(TensorOps::InnerProduct(basis[row], basis[column], g), expected, 2.0e-13)
                << "row=" << row << " column=" << column;
        }
    }
}

TEST(ObserverFrameTests, InvalidWorldlineCannotBecomeASilentStaticObserver) {
    KerrSchildFamily metric(KerrSchildParams::Minkowski());
    Vec4 position;
    Vec4 spacelike_velocity;
    spacelike_velocity(1) = 1.0;
    EXPECT_DEATH((void)Geodesic::CreateObserver(position, spacelike_velocity, &metric),
                 "precondition.*enforced, terminating");
}

TEST(ObserverFrameTests, MovingCameraLaunchMatchesIndependentLorentzTransform) {
    KerrSchildFamily metric(KerrSchildParams::Minkowski());
    Vec4 position;
    Metric4d g;
    Tensor<Dual<double>, 4, 4, 4> derivatives;
    metric.Evaluate(position, g, derivatives);

    std::array<Vec4, 3> seeds;
    seeds[0](1) = 1.0;
    seeds[1](2) = 1.0;
    seeds[2](3) = 1.0;
    const Metric4d inverse = TensorOps::Inverse(g);
    const auto reference = relativity::EulerianObserverFrame(g, inverse, seeds);
    ASSERT_TRUE(reference.has_value());
    const std::array<double, 3> beta{0.25, -0.4, 0.1};
    const auto moving = relativity::BoostObserverFrame(*reference, beta);
    ASSERT_TRUE(moving.has_value());
    const std::array<double, 3> screen_direction{-0.8, 0.3, std::sqrt(0.27)};
    const auto past_ray = relativity::PastDirectedCameraRay(*moving, screen_direction);
    ASSERT_TRUE(past_ray.has_value());

    const double beta_squared = beta[0] * beta[0] + beta[1] * beta[1] + beta[2] * beta[2];
    const double gamma = 1.0 / std::sqrt(1.0 - beta_squared);
    const double beta_dot_screen = beta[0] * screen_direction[0] + beta[1] * screen_direction[1] +
                                   beta[2] * screen_direction[2];
    // Independent boost of the physical photon q'=(1,-n) from the moving
    // camera rest frame into the static Minkowski frame.
    const double expected_time = gamma * (1.0 - beta_dot_screen);
    const double coefficient = gamma - (gamma - 1.0) * beta_dot_screen / beta_squared;
    EXPECT_NEAR(-(*past_ray)(0), expected_time, 2.0e-14);
    for (int component = 0; component < 3; ++component) {
        const double expected_spatial =
            -screen_direction[component] + coefficient * beta[component];
        EXPECT_NEAR(-(*past_ray)(component + 1), expected_spatial, 2.0e-14);
    }
}

TEST(ObserverFrameTests, KerrCameraRayIsPastNullAndLaunchFrequencyIsOne) {
    KerrSchildFamily metric(KerrSchildParams::Kerr(1.0, 0.9));
    Vec4 position;
    position(1) = 30.0;
    position(2) = 2.0;
    position(3) = 4.0;
    Metric4d g;
    Tensor<Dual<double>, 4, 4, 4> derivatives;
    metric.Evaluate(position, g, derivatives);

    std::array<Vec4, 3> seeds;
    seeds[0](1) = 1.0;
    seeds[1](2) = 1.0;
    seeds[2](3) = 1.0;
    Metric4d inverse;
    ASSERT_TRUE(metric.InverseMetric(position, inverse));
    const auto reference = relativity::EulerianObserverFrame(g, inverse, seeds);
    ASSERT_TRUE(reference.has_value());
    const auto moving = relativity::BoostObserverFrame(*reference, {0.7, -0.2, 0.1});
    ASSERT_TRUE(moving.has_value());
    const auto past_ray = relativity::PastDirectedCameraRay(*moving, {-0.9, 0.3, std::sqrt(0.1)});
    ASSERT_TRUE(past_ray.has_value());
    const auto screen = relativity::ObserverScreenBasis(*moving, {-0.9, 0.3, std::sqrt(0.1)});
    ASSERT_TRUE(screen.has_value());

    EXPECT_LT((*past_ray)(0), 0.0);
    EXPECT_NEAR(TensorOps::InnerProduct(*past_ray, *past_ray, g), 0.0, 3.0e-13);
    const Vec4 physical_photon = *past_ray * -1.0;
    EXPECT_NEAR(-TensorOps::InnerProduct(physical_photon, moving->time, g), 1.0, 3.0e-13);
    const std::array<Vec4, 4> tetrad = {moving->time, moving->spatial[0], moving->spatial[1],
                                        moving->spatial[2]};
    for (int row = 0; row < 4; ++row) {
        for (int column = 0; column < 4; ++column) {
            const double expected = row == column ? (row == 0 ? -1.0 : 1.0) : 0.0;
            EXPECT_NEAR(TensorOps::InnerProduct(tetrad[row], tetrad[column], g), expected, 3.0e-13);
        }
    }
    for (int axis = 0; axis < 2; ++axis) {
        EXPECT_NEAR(TensorOps::InnerProduct((*screen)[axis], moving->time, g), 0.0, 3.0e-13);
        EXPECT_NEAR(TensorOps::InnerProduct((*screen)[axis], *past_ray, g), 0.0, 3.0e-13);
        EXPECT_NEAR(TensorOps::InnerProduct((*screen)[axis], (*screen)[axis], g), 1.0, 3.0e-13);
    }
    EXPECT_NEAR(TensorOps::InnerProduct((*screen)[0], (*screen)[1], g), 0.0, 3.0e-13);
}

TEST(ObserverFrameTests, EulerianReferenceRemainsTimelikeInsideTheKerrErgosphere) {
    KerrSchildFamily metric(KerrSchildParams::Kerr(1.0, 0.998));
    const double radius = 1.5;
    Vec4 position;
    position(1) = std::sqrt(radius * radius + 0.998 * 0.998);
    Metric4d g;
    Metric4d inverse;
    Tensor<Dual<double>, 4, 4, 4> derivatives;
    metric.Evaluate(position, g, derivatives);
    ASSERT_TRUE(metric.InverseMetric(position, inverse));
    ASSERT_GT(g(0, 0).real, 0.0) << "probe must lie inside the static limit";

    std::array<Vec4, 3> seeds;
    seeds[0](1) = 1.0;
    seeds[1](2) = 1.0;
    seeds[2](3) = 1.0;
    const auto frame = relativity::EulerianObserverFrame(g, inverse, seeds);
    ASSERT_TRUE(frame.has_value());
    EXPECT_GT(frame->time(0), 0.0);
    EXPECT_NEAR(TensorOps::InnerProduct(frame->time, frame->time, g), -1.0, 3.0e-13);

    const auto ray = relativity::PastDirectedCameraRay(*frame, {-0.8, 0.3, 0.2});
    ASSERT_TRUE(ray.has_value());
    Vec4 trial;
    trial(2) = 1.0;
    const auto screen = relativity::ProjectToObserverScreen(g, frame->time, *ray, trial);
    ASSERT_TRUE(screen.has_value());
    EXPECT_NEAR(TensorOps::InnerProduct(*screen, frame->time, g), 0.0, 3.0e-13);
    EXPECT_NEAR(TensorOps::InnerProduct(*screen, *ray, g), 0.0, 3.0e-13);
    EXPECT_NEAR(TensorOps::InnerProduct(*screen, *screen, g), 1.0, 3.0e-13);
}
