// Camera-worldline aberration tests (specification P5).
//
// An arbitrary camera four-velocity applies special-relativistic aberration to
// every generated ray in the camera's local frame (core/camera.h). For motion
// along the view axis the aberrated angle satisfies the light-cone Lorentz form
// cos theta' = (cos theta - beta) / (1 - beta cos theta) (MTW eq 2.29), which
// this suite pins to 1e-12; the perpendicular scaling keeps the direction a unit
// vector. beta = 0 is exactly the un-aberrated ray, so every lens model composes
// and the pinned render is unmoved (the byte-pin cmp confirms that separately).

#include "sirius/core/camera.h"

#include <gtest/gtest.h>

#include <array>
#include <cmath>

namespace {

using namespace sirius::core;

CameraRay RayAlong(double dir1, double dir2, double dir3) {
    CameraRay r;
    r.direction(1) = dir1;
    r.direction(2) = dir2;
    r.direction(3) = dir3;
    return r;
}

// -----------------------------------------------------------------------------
// Analytic aberration along the view axis (direction index 1) at several angles
// and speeds, to double precision.
// -----------------------------------------------------------------------------
TEST(CameraAberrationTest, MatchesAnalyticFormulaAlongViewAxis) {
    for (double beta : {0.1, 0.3, 0.6, 0.9}) {
        double gamma = 1.0 / std::sqrt(1.0 - beta * beta);
        for (double theta = 0.1; theta < 3.1; theta += 0.37) {
            double c = std::cos(theta), s = std::sin(theta);
            // Direction at angle theta from axis 1, in the (1, 2) plane.
            CameraRay ray = RayAlong(c, s, 0.0);
            AberrateRay(ray, beta, 0.0, 0.0);

            double denom = 1.0 - beta * c;
            double cos_prime = (c - beta) / denom;
            double sin_prime = s / (gamma * denom);
            // The result must be unit and match the analytic components.
            double len = std::sqrt(ray.direction(1) * ray.direction(1) +
                                   ray.direction(2) * ray.direction(2) +
                                   ray.direction(3) * ray.direction(3));
            EXPECT_NEAR(len, 1.0, 1e-12);
            EXPECT_NEAR(ray.direction(1), cos_prime, 1e-12)
                << "beta=" << beta << " theta=" << theta;
            EXPECT_NEAR(ray.direction(2), sin_prime, 1e-12)
                << "beta=" << beta << " theta=" << theta;
            EXPECT_NEAR(ray.direction(3), 0.0, 1e-12);
        }
    }
}

TEST(CameraAberrationTest, MatchesGeneralLorentzBoostForThreeAxisVelocities) {
    const std::array<std::array<double, 3>, 5> velocities = {
        std::array<double, 3>{0.1, 0.2, -0.3},   std::array<double, 3>{-0.4, 0.05, 0.2},
        std::array<double, 3>{0.0, -0.7, 0.1},   std::array<double, 3>{0.55, 0.55, 0.55},
        std::array<double, 3>{-0.97, 0.1, 0.05},
    };
    const std::array<std::array<double, 3>, 4> directions = {
        std::array<double, 3>{1.0, 0.0, 0.0},
        std::array<double, 3>{0.0, 1.0, 0.0},
        std::array<double, 3>{0.3, -0.4, std::sqrt(0.75)},
        std::array<double, 3>{-0.2, 0.9, std::sqrt(0.15)},
    };

    for (const auto& beta : velocities) {
        const double beta_squared = beta[0] * beta[0] + beta[1] * beta[1] + beta[2] * beta[2];
        ASSERT_LT(beta_squared, 1.0);
        const double gamma = 1.0 / std::sqrt(1.0 - beta_squared);
        for (const auto& direction : directions) {
            CameraRay actual = RayAlong(direction[0], direction[1], direction[2]);
            AberrateRay(actual, beta[0], beta[1], beta[2]);

            // Independent Lorentz transform of photon four-momentum k=(1,n):
            // k'^0=gamma(1-beta.n), k'=n+[((gamma-1)beta.n/beta^2)-gamma]beta.
            const double dot =
                beta[0] * direction[0] + beta[1] * direction[1] + beta[2] * direction[2];
            const double time = gamma * (1.0 - dot);
            const double coefficient = (gamma - 1.0) * dot / beta_squared - gamma;
            for (int component = 0; component < 3; ++component) {
                const double expected =
                    (direction[component] + coefficient * beta[component]) / time;
                EXPECT_NEAR(actual.direction(component + 1), expected, 2.0e-12);
                EXPECT_TRUE(std::isfinite(actual.direction(component + 1)));
            }
            const double length = std::sqrt(actual.direction(1) * actual.direction(1) +
                                            actual.direction(2) * actual.direction(2) +
                                            actual.direction(3) * actual.direction(3));
            EXPECT_NEAR(length, 1.0, 2.0e-12);
        }
    }
}

// -----------------------------------------------------------------------------
// beta = 0 is an exact no-op: the aberrated ray is bit-identical to the input.
// -----------------------------------------------------------------------------
TEST(CameraAberrationTest, ZeroBetaIsExactNoOp) {
    CameraRay ray = RayAlong(0.3, -0.4, 0.86602540378);
    CameraRay original = ray;
    AberrateRay(ray, 0.0, 0.0, 0.0);
    EXPECT_EQ(ray.direction(1), original.direction(1));
    EXPECT_EQ(ray.direction(2), original.direction(2));
    EXPECT_EQ(ray.direction(3), original.direction(3));

    EXPECT_DEATH(
        {
            CameraRay invalid = RayAlong(1.0, 0.0, 0.0);
            AberrateRay(invalid, 1.0, 0.0, 0.0);
        },
        "precondition.*enforced, terminating");
}

// -----------------------------------------------------------------------------
// Forward motion beams the aberrated ray toward the boost axis (the headlight
// effect): the forward component increases for an off-axis ray.
// -----------------------------------------------------------------------------
TEST(CameraAberrationTest, ForwardMotionBeamsTowardAxis) {
    CameraRay ray = RayAlong(std::cos(1.0), std::sin(1.0), 0.0);
    double before = ray.direction(1);
    AberrateRay(ray, -0.5, 0.0, 0.0);  // Motion toward +axis-1 (see the formula sign).
    EXPECT_GT(ray.direction(1), before);
}

// -----------------------------------------------------------------------------
// Aberration composes over every lens model through GenerateRayAberrated: the
// boosted ray differs from the static one, and beta = 0 leaves it identical.
// -----------------------------------------------------------------------------
TEST(CameraAberrationTest, ComposesOverLensModels) {
    CameraConfig cfg;
    cfg.width = 64;
    cfg.height = 64;
    cfg.beta_x = 0.3;  // Along the view axis.

    for (LensType lens : {LensType::Pinhole, LensType::ThinLens, LensType::Fisheye}) {
        auto cam = CreateCamera(lens, cfg);
        CameraRay plain = cam->GenerateRay(20, 24, 0.5f, 0.5f);
        CameraRay aberrated = cam->GenerateRayAberrated(20, 24, 0.5f, 0.5f);
        double diff = std::abs(plain.direction(1) - aberrated.direction(1)) +
                      std::abs(plain.direction(2) - aberrated.direction(2)) +
                      std::abs(plain.direction(3) - aberrated.direction(3));
        EXPECT_GT(diff, 1e-6) << "aberration had no effect for lens " << static_cast<int>(lens);
        EXPECT_TRUE(std::isfinite(aberrated.direction(1)));
    }

    // beta = 0 -> GenerateRayAberrated equals GenerateRay for every model.
    CameraConfig cfg0 = cfg;
    cfg0.beta_x = 0.0;
    auto cam0 = CreateCamera(LensType::Pinhole, cfg0);
    CameraRay a = cam0->GenerateRay(11, 40, 0.5f, 0.5f);
    CameraRay b = cam0->GenerateRayAberrated(11, 40, 0.5f, 0.5f);
    EXPECT_EQ(a.direction(1), b.direction(1));
    EXPECT_EQ(a.direction(2), b.direction(2));
    EXPECT_EQ(a.direction(3), b.direction(3));
}

}  // namespace
