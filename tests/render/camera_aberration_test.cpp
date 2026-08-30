// Typed camera-worldline composition tests (specification P5).
//
// Lens projection produces a direction in the camera's instantaneous rest
// frame. GenerateRayForObserver binds the finite timelike worldline without
// applying a Euclidean direction edit; the metric-aware CPU/GPU launch paths
// perform the tetrad boost. The exact Lorentz and null-frame witnesses live in
// the non-render ObserverFrameTests suite.

#include "sirius/core/camera.h"

#include <gtest/gtest.h>

#include <cmath>
#include <limits>

namespace {

using namespace sirius::core;

TEST(CameraWorldlineTest, RestScreenRayAndWorldlineComposeOverLensModels) {
    CameraConfig config;
    config.width = 64;
    config.height = 64;
    config.beta_x = 0.3;
    config.beta_y = -0.2;
    config.beta_z = 0.1;

    for (LensType lens : {LensType::Pinhole, LensType::ThinLens, LensType::Fisheye}) {
        auto camera = CreateCamera(lens, config);
        const CameraRay rest_ray = camera->GenerateRay(20, 24, 0.5f, 0.5f);
        const CameraRay bound_ray = camera->GenerateRayForObserver(20, 24, 0.5f, 0.5f);
        for (int component = 0; component < 4; ++component) {
            EXPECT_EQ(bound_ray.direction(component), rest_ray.direction(component));
        }
        EXPECT_DOUBLE_EQ(bound_ray.beta_forward, config.beta_x);
        EXPECT_DOUBLE_EQ(bound_ray.beta_up, config.beta_y);
        EXPECT_DOUBLE_EQ(bound_ray.beta_right, config.beta_z);
    }
}

TEST(CameraWorldlineTest, ZeroVelocityIsExactlyRepresented) {
    PinholeCamera camera;
    const CameraRay rest_ray = camera.GenerateRay(11, 40, 0.5f, 0.5f);
    const CameraRay bound_ray = camera.GenerateRayForObserver(11, 40, 0.5f, 0.5f);
    for (int component = 0; component < 4; ++component) {
        EXPECT_EQ(bound_ray.direction(component), rest_ray.direction(component));
    }
    EXPECT_DOUBLE_EQ(bound_ray.beta_forward, 0.0);
    EXPECT_DOUBLE_EQ(bound_ray.beta_up, 0.0);
    EXPECT_DOUBLE_EQ(bound_ray.beta_right, 0.0);
}

TEST(CameraWorldlineTest, InvalidInternalWorldlineFailsClosed) {
    CameraConfig invalid;
    invalid.beta_x = 1.0;
    EXPECT_TRUE(CameraConfigIssue(LensType::Pinhole, invalid).has_value());
    EXPECT_DEATH((void)PinholeCamera(invalid), "precondition.*enforced, terminating");

    invalid.beta_x = std::numeric_limits<double>::quiet_NaN();
    EXPECT_TRUE(CameraConfigIssue(LensType::Pinhole, invalid).has_value());
    EXPECT_DEATH((void)PinholeCamera(invalid), "precondition.*enforced, terminating");
}

}  // namespace
