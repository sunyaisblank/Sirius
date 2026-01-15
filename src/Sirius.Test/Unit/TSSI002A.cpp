// TSSI002A.cpp - Scene Intersection and Anti-Aliasing Tests
// Tests for KNAA001A (beam sampling) and KNSI001A (scene intersection)

#include <gtest/gtest.h>
#include <cmath>
#include "../../Sirius.Render/Integration/KNAA001A.h"
#include "../../Sirius.Render/Integration/KNSI001A.h"
#include "../../Sirius.Core/Disk/PHAD001A.h"

using namespace sirius::kernel;
using namespace sirius::physics;

namespace {

//==============================================================================
// Beam Sampler Tests
//==============================================================================

TEST(BeamSamplerTest, GaussianWeightDecay) {
    BeamSampler sampler;
    
    // Create a beam with known ellipse parameters
    BeamStateD beam;
    beam.majorAxis = 0.01;   // 10 mrad
    beam.minorAxis = 0.005;  // 5 mrad
    beam.orientation = 0.0;
    beam.initialPixelSolidAngle = 1e-6;
    
    // Create catalog with stars at varying distances
    std::vector<StarEntry> catalog;
    
    // Star at beam centre
    catalog.push_back({0.0, 0.0, 5.0, 6000.0});
    
    // Star at 1σ along major axis
    catalog.push_back({0.01, 0.0, 5.0, 6000.0});
    
    // Star at 3σ (should have low weight)
    catalog.push_back({0.03, 0.0, 5.0, 6000.0});
    
    auto samples = sampler.sampleStarfield(beam, catalog, 0.0, 0.0);
    
    EXPECT_GE(samples.size(), 2) << "Should include at least 2 stars";
    
    // Centre star should have weight ~1
    if (samples.size() > 0) {
        EXPECT_NEAR(samples[0].weight, 1.0, 0.01)
            << "Centre star should have weight ~1";
    }
    
    // 1σ star should have weight ~0.61 (exp(-0.5))
    if (samples.size() > 1) {
        EXPECT_NEAR(samples[1].weight, std::exp(-0.5), 0.01)
            << "1σ star should have Gaussian weight";
    }
}

TEST(BeamSamplerTest, EllipseOrientation) {
    BeamSampler sampler;
    
    BeamStateD beam;
    beam.majorAxis = 0.02;
    beam.minorAxis = 0.01;
    beam.orientation = M_PI/4;  // 45° rotation
    beam.initialPixelSolidAngle = 1e-6;
    
    // Star along rotated major axis
    double dalpha = 0.01 * std::cos(M_PI/4);
    double ddelta = 0.01 * std::sin(M_PI/4);
    
    std::vector<StarEntry> catalog;
    catalog.push_back({dalpha, ddelta, 5.0, 6000.0});
    
    auto samples = sampler.sampleStarfield(beam, catalog, 0.0, 0.0);
    
    EXPECT_EQ(samples.size(), 1) << "Should find 1 star";
    if (samples.size() > 0) {
        // At 0.5σ along major axis
        double expected_weight = std::exp(-0.5 * 0.25);  // (0.5)^2 / 2
        EXPECT_GT(samples[0].weight, 0.7)
            << "Star along major axis should have high weight";
    }
}

TEST(BeamSamplerTest, SolidAngleComputation) {
    BeamStateD beam;
    beam.majorAxis = 0.01;
    beam.minorAxis = 0.005;
    
    double omega = BeamSampler::beamSolidAngle(beam);
    double expected = M_PI * 0.01 * 0.005;  // π × a × b
    
    EXPECT_NEAR(omega, expected, 1e-10)
        << "Solid angle should be π×a×b";
}

//==============================================================================
// Scene Intersection Tests
//==============================================================================

TEST(SceneIntersectorTest, DiskJacobianTransform) {
    AccretionDiskD::Config diskConfig;
    diskConfig.M = 10.0;
    AccretionDiskD disk(diskConfig);
    
    SceneIntersector intersector(&disk);
    
    // Create beam at disk
    BeamStateD beam;
    beam.x = Vec4d(0, 10.0, M_PI/2, 0);  // At equator, r=10M
    beam.initialPixelSolidAngle = 1e-8;
    
    // Set up Jacobian with specific test values
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            beam.J[i][j] = 0.0;
        }
    }
    // Set the relevant elements for disk projection
    beam.J[1][2] = 2.0;  // ∂r/∂θ₀
    beam.J[1][3] = 0.5;  // ∂r/∂φ₀
    beam.J[3][2] = 0.1;  // ∂φ/∂θ₀
    beam.J[3][3] = 1.5;  // ∂φ/∂φ₀
    
    auto diskJ = intersector.transformJacobianToDisk(beam);
    
    // Check extracted elements
    EXPECT_NEAR(diskJ.J[0][0], 2.0, 1e-10);   // ∂r/∂θ₀
    EXPECT_NEAR(diskJ.J[0][1], 0.5, 1e-10);   // ∂r/∂φ₀
    EXPECT_NEAR(diskJ.J[1][0], 0.1, 1e-10);   // ∂φ/∂θ₀
    EXPECT_NEAR(diskJ.J[1][1], 1.5, 1e-10);   // ∂φ/∂φ₀
    
    // Determinant = 2.0 × 1.5 - 0.5 × 0.1 = 2.95
    EXPECT_NEAR(diskJ.determinant, 2.95, 1e-10);
}

TEST(SceneIntersectorTest, DiskIntersectionDetection) {
    AccretionDiskD::Config diskConfig;
    diskConfig.M = 10.0;
    diskConfig.a_star = 0.0;  // Schwarzschild
    AccretionDiskD disk(diskConfig);
    
    SceneIntersector intersector(&disk);
    
    // Beam at equator within disk
    BeamStateD beamInDisk;
    beamInDisk.x = Vec4d(0, 15.0, M_PI/2, 0);
    beamInDisk.terminated = false;
    beamInDisk.initialPixelSolidAngle = 1e-8;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            beamInDisk.J[i][j] = (i == j) ? 1.0 : 0.0;
    
    auto result = intersector.intersectScene(beamInDisk);
    
    EXPECT_EQ(result.type, SceneIntersector::IntersectionResult::DISK)
        << "Should intersect disk at equator";
    EXPECT_NEAR(result.r, 15.0, 1e-10);
}

TEST(SceneIntersectorTest, CelestialSphereDetection) {
    AccretionDiskD::Config diskConfig;
    AccretionDiskD disk(diskConfig);
    
    SceneIntersector::Config config;
    config.escapeRadius = 100;
    SceneIntersector intersector(&disk, config);
    
    // Beam escaped to celestial sphere
    BeamStateD beamEscaped;
    beamEscaped.x = Vec4d(0, 200.0, M_PI/3, 0);
    beamEscaped.terminated = false;
    beamEscaped.solidAngle = 1e-8;
    beamEscaped.initialPixelSolidAngle = 1e-8;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            beamEscaped.J[i][j] = (i == j) ? 1.0 : 0.0;
    
    auto result = intersector.intersectScene(beamEscaped);
    
    EXPECT_EQ(result.type, SceneIntersector::IntersectionResult::CELESTIAL_SPHERE)
        << "Should detect celestial sphere";
}

TEST(SceneIntersectorTest, HorizonDetection) {
    AccretionDiskD::Config diskConfig;
    AccretionDiskD disk(diskConfig);
    
    SceneIntersector intersector(&disk);
    
    // Beam terminated inside ISCO
    BeamStateD beamHorizon;
    beamHorizon.x = Vec4d(0, 1.5, M_PI/2, 0);
    beamHorizon.terminated = true;
    
    auto result = intersector.intersectScene(beamHorizon);
    
    EXPECT_EQ(result.type, SceneIntersector::IntersectionResult::HORIZON)
        << "Should detect horizon";
}

TEST(SceneIntersectorTest, DiskEmissionHasRedshift) {
    AccretionDiskD::Config diskConfig;
    diskConfig.M = 10.0;
    diskConfig.Mdot = 1e-8;
    AccretionDiskD disk(diskConfig);
    
    SceneIntersector intersector(&disk);
    
    // Beam at disk
    BeamStateD beam;
    beam.x = Vec4d(0, 20.0, M_PI/2, 0);  // r = 20M
    beam.terminated = false;
    beam.initialPixelSolidAngle = 1e-8;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            beam.J[i][j] = (i == j) ? 1.0 : 0.0;
    
    auto result = intersector.intersectScene(beam);
    
    EXPECT_EQ(result.type, SceneIntersector::IntersectionResult::DISK);
    EXPECT_GT(result.g_factor, 0) << "g-factor should be positive";
    EXPECT_LT(result.g_factor, 1) << "g-factor should be < 1 (redshifted)";
    EXPECT_GT(result.emission.totalEnergy(), 0)
        << "Disk should have non-zero emission";
}

} // namespace
