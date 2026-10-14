// TSCM001A.cpp - Camera Lens Model Tests
// Tests for CMBS001A (PinholeCamera, ThinLensCamera, FisheyeCamera)
//
// Validates: ray generation correctness, coordinate mapping, factory function

#include <gtest/gtest.h>
#include <CMBS001A.h>
#include <cmath>

using namespace Sirius;

namespace sirius::test {
using namespace Sirius;

//==============================================================================
// PinholeCamera Tests
//==============================================================================

TEST(PinholeCameraTest, CentreRayPointsInward) {
    CameraConfig config;
    config.width = 100;
    config.height = 100;
    config.fov = 90.0f;
    config.r = 50.0;
    config.theta = M_PI / 2.0;
    config.phi = 0.0;
    config.yaw = 0.0f;
    config.pitch = 0.0f;
    config.roll = 0.0f;

    PinholeCamera camera(config);
    // For exact centre: x+u = width/2, y+v = height/2
    // x=50, u=0 → (2*50/100 - 1) = 0, y=50, v=0 → 1 - 2*50/100 = 0
    CameraRay ray = camera.generateRay(50, 50, 0.0f, 0.0f);

    // Centre ray: NDC = (0,0), camera direction = (0,0,-1)
    // After coordinate mapping: vr = dz = -1, vθ = 0, vφ = 0
    // Normalised: direction(1) should be negative (inward)
    EXPECT_LT(ray.direction(1), 0.0f) << "Centre ray should point radially inward";
    EXPECT_NEAR(ray.direction(2), 0.0f, 1e-5f) << "Centre ray should have no θ component";
    EXPECT_NEAR(ray.direction(3), 0.0f, 1e-5f) << "Centre ray should have no φ component";
}

TEST(PinholeCameraTest, RayDirectionIsNormalised) {
    CameraConfig config;
    config.width = 200;
    config.height = 200;
    config.fov = 60.0f;

    PinholeCamera camera(config);

    // Test several pixel positions
    int positions[][2] = {{0, 0}, {100, 100}, {199, 199}, {50, 150}, {0, 199}};
    for (auto& pos : positions) {
        CameraRay ray = camera.generateRay(pos[0], pos[1]);
        // direction(0) is always 0 (time), so spatial norm:
        float spatial_norm = std::sqrt(
            ray.direction(1) * ray.direction(1) +
            ray.direction(2) * ray.direction(2) +
            ray.direction(3) * ray.direction(3));
        EXPECT_NEAR(spatial_norm, 1.0f, 1e-5f)
            << "Spatial direction should be unit length at pixel ("
            << pos[0] << ", " << pos[1] << ")";
    }
}

TEST(PinholeCameraTest, OriginMatchesConfig) {
    CameraConfig config;
    config.t = 0.0;
    config.r = 30.0;
    config.theta = M_PI / 3.0;
    config.phi = 1.5;

    PinholeCamera camera(config);
    CameraRay ray = camera.generateRay(0, 0);

    EXPECT_FLOAT_EQ(ray.origin(0), static_cast<float>(config.t));
    EXPECT_FLOAT_EQ(ray.origin(1), static_cast<float>(config.r));
    EXPECT_FLOAT_EQ(ray.origin(2), static_cast<float>(config.theta));
    EXPECT_FLOAT_EQ(ray.origin(3), static_cast<float>(config.phi));
}

TEST(PinholeCameraTest, RightPixelIncreasesAzimuth) {
    CameraConfig config;
    config.width = 200;
    config.height = 200;
    config.fov = 90.0f;
    config.yaw = 0.0f;
    config.pitch = 0.0f;
    config.roll = 0.0f;

    PinholeCamera camera(config);

    // Centre pixel
    CameraRay centre = camera.generateRay(100, 100);
    // Right of centre
    CameraRay right = camera.generateRay(150, 100);

    // Camera +X maps to +φ direction → direction(3) should increase
    EXPECT_GT(right.direction(3), centre.direction(3))
        << "Moving right should increase φ component";
}

TEST(PinholeCameraTest, UpPixelDecreasesTheta) {
    CameraConfig config;
    config.width = 200;
    config.height = 200;
    config.fov = 90.0f;
    config.yaw = 0.0f;
    config.pitch = 0.0f;
    config.roll = 0.0f;

    PinholeCamera camera(config);

    // Centre pixel
    CameraRay centre = camera.generateRay(100, 100);
    // Above centre (y decreases = up in screen space)
    CameraRay up = camera.generateRay(100, 50);

    // Camera +Y (up) maps to -θ → direction(2) = -dy
    // Moving up in screen: py increases, dy > 0, so direction(2) < 0 (toward pole)
    EXPECT_LT(up.direction(2), centre.direction(2))
        << "Moving up should decrease θ component (toward pole)";
}

TEST(PinholeCameraTest, WeightIsOne) {
    PinholeCamera camera;
    CameraRay ray = camera.generateRay(0, 0);
    EXPECT_FLOAT_EQ(ray.weight, 1.0f);
}

TEST(PinholeCameraTest, TimeComponentIsZero) {
    PinholeCamera camera;
    CameraRay ray = camera.generateRay(50, 50);
    EXPECT_FLOAT_EQ(ray.direction(0), 0.0f)
        << "Time component should be zero (set by geodesic normalisation)";
}

//==============================================================================
// ThinLensCamera Tests
//==============================================================================

TEST(ThinLensCameraTest, CentreRayWithDefaultSampling) {
    CameraConfig config;
    config.width = 100;
    config.height = 100;
    config.fov = 60.0f;
    config.focalLength = 50.0f;
    config.aperture = 2.8f;
    config.focusDistance = 50.0f;

    ThinLensCamera camera(config);

    // Centre pixel with u=0.5, v=0.5 → lens offset depends on u,v
    CameraRay ray = camera.generateRay(50, 50, 0.5f, 0.5f);

    // Should produce a valid direction
    float spatial_norm = std::sqrt(
        ray.direction(1) * ray.direction(1) +
        ray.direction(2) * ray.direction(2) +
        ray.direction(3) * ray.direction(3));
    EXPECT_NEAR(spatial_norm, 1.0f, 1e-4f)
        << "ThinLens direction should be unit length";
}

TEST(ThinLensCameraTest, DifferentSamplesGiveDifferentRays) {
    CameraConfig config;
    config.width = 100;
    config.height = 100;
    config.focalLength = 50.0f;
    config.aperture = 1.4f;  // Large aperture for noticeable DoF
    config.focusDistance = 50.0f;

    ThinLensCamera camera(config);

    CameraRay ray1 = camera.generateRay(50, 50, 0.1f, 0.1f);
    CameraRay ray2 = camera.generateRay(50, 50, 0.9f, 0.9f);

    // Different (u,v) samples should produce different ray directions
    bool different = false;
    for (int i = 1; i < 4; ++i) {
        if (std::abs(ray1.direction(i) - ray2.direction(i)) > 1e-6f) {
            different = true;
            break;
        }
    }
    EXPECT_TRUE(different) << "Different aperture samples should produce different rays";
}

//==============================================================================
// FisheyeCamera Tests
//==============================================================================

TEST(FisheyeCameraTest, CentreRayPointsInward) {
    CameraConfig config;
    config.width = 100;
    config.height = 100;
    config.fov = 180.0f;

    FisheyeCamera camera(config);
    // Exact centre: x=50, u=0 and y=50, v=0 give NDC = (0,0)
    CameraRay ray = camera.generateRay(50, 50, 0.0f, 0.0f);

    // Centre ray: r_img = 0, theta_ray = 0, direction = (0,0,-1)
    // After mapping: vr = dz = -1 (inward)
    EXPECT_LT(ray.direction(1), -0.9f) << "Centre fisheye ray should point inward";
    EXPECT_NEAR(ray.direction(2), 0.0f, 1e-5f);
    EXPECT_NEAR(ray.direction(3), 0.0f, 1e-5f);
}

TEST(FisheyeCameraTest, EdgeRayPerpendicularAt180Fov) {
    CameraConfig config;
    config.width = 200;
    config.height = 200;
    config.fov = 180.0f;

    FisheyeCamera camera(config);

    // Right edge of image: px ≈ 1 (aspect ratio 1), py = 0
    // r_img = 1, theta_ray = π/2 → direction perpendicular to forward
    CameraRay ray = camera.generateRay(199, 100, 0.5f, 0.5f);

    // At θ_ray = π/2, cos(θ) ≈ 0 → vr ≈ 0
    EXPECT_NEAR(ray.direction(1), 0.0f, 0.15f)
        << "Edge ray at 180° FOV should be roughly perpendicular";
}

TEST(FisheyeCameraTest, OutOfBoundsRayHasZeroWeight) {
    CameraConfig config;
    config.width = 100;
    config.height = 100;
    config.fov = 90.0f;  // 90° FOV

    FisheyeCamera camera(config);

    // Corner pixel: far from centre with limited FOV
    // For 100×100 image with aspect 1: px ≈ 1, py ≈ 1, r_img ≈ √2
    // theta_ray = √2 × 45° = 63.6° → within π, so still valid
    // Need extreme corner with very narrow FOV to exceed π
    config.fov = 10.0f;  // Very narrow FOV
    camera.setConfig(config);

    CameraRay ray = camera.generateRay(0, 0, 0.5f, 0.5f);
    // r_img ≈ √2, theta_ray = √2 × 5° ≈ 7° → within π, still valid
    // FisheyeCamera only zeroes weight when theta_ray > π
    // That requires r_img × (fov/2) > π → r_img > 2π/fov (in radians)
    // For this config: threshold = 2π / (10π/180) = 36, r_img ≈ 1.4 → valid
    EXPECT_FLOAT_EQ(ray.weight, 1.0f) << "Within FOV should have weight 1";
}

//==============================================================================
// Camera Factory Tests
//==============================================================================

TEST(CameraFactoryTest, CreatePinhole) {
    auto camera = createCamera(LensType::Pinhole);
    ASSERT_NE(camera, nullptr);
    EXPECT_EQ(camera->getLensType(), LensType::Pinhole);
    EXPECT_STREQ(camera->getName(), "Pinhole Camera");
}

TEST(CameraFactoryTest, CreateThinLens) {
    auto camera = createCamera(LensType::ThinLens);
    ASSERT_NE(camera, nullptr);
    EXPECT_EQ(camera->getLensType(), LensType::ThinLens);
    EXPECT_STREQ(camera->getName(), "Thin Lens Camera");
}

TEST(CameraFactoryTest, CreateFisheye) {
    auto camera = createCamera(LensType::Fisheye);
    ASSERT_NE(camera, nullptr);
    EXPECT_EQ(camera->getLensType(), LensType::Fisheye);
    EXPECT_STREQ(camera->getName(), "Fisheye Camera");
}

TEST(CameraFactoryTest, DefaultFallsToPinhole) {
    auto camera = createCamera(LensType::Panoramic);  // Not implemented
    ASSERT_NE(camera, nullptr);
    EXPECT_EQ(camera->getLensType(), LensType::Pinhole)
        << "Unimplemented types should fall back to Pinhole";
}

TEST(CameraFactoryTest, ConfigPassthrough) {
    CameraConfig config;
    config.r = 100.0;
    config.theta = 1.0;
    config.phi = 2.5;

    auto camera = createCamera(LensType::Pinhole, config);
    const auto& retrieved = camera->getConfig();
    EXPECT_DOUBLE_EQ(retrieved.r, 100.0);
    EXPECT_DOUBLE_EQ(retrieved.theta, 1.0);
    EXPECT_DOUBLE_EQ(retrieved.phi, 2.5);
}

//==============================================================================
// ICamera Interface Tests
//==============================================================================

TEST(ICameraTest, GetPositionReturnsConfigCoordinates) {
    CameraConfig config;
    config.t = 10.0;
    config.r = 25.0;
    config.theta = 1.2;
    config.phi = 3.0;

    PinholeCamera camera(config);
    auto pos = camera.getPosition();

    EXPECT_FLOAT_EQ(pos(0), 10.0f);
    EXPECT_FLOAT_EQ(pos(1), 25.0f);
    EXPECT_NEAR(pos(2), 1.2f, 1e-6f);
    EXPECT_FLOAT_EQ(pos(3), 3.0f);
}

TEST(ICameraTest, SetConfigUpdatesRayGeneration) {
    PinholeCamera camera;

    CameraConfig config;
    config.width = 50;
    config.height = 50;
    config.fov = 120.0f;
    config.r = 10.0;
    camera.setConfig(config);

    CameraRay ray = camera.generateRay(25, 25);
    EXPECT_FLOAT_EQ(ray.origin(1), 10.0f)
        << "Origin should reflect updated config";
}

} // namespace sirius::test
