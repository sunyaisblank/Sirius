// TSIN003A.cpp - Render Pipeline Integration Tests

#include <gtest/gtest.h>
#include <cmath>
#include <memory>

#include "PHMT100A.h"
#include "PHGD001A.h"

namespace sirius::test {
using namespace Sirius;

// Constants
constexpr double M = 1.0;  // Black hole mass
constexpr double PI = 3.14159265358979323846;

// Coordinate conversion helpers
inline double getRadius(const Vec4& pos) {
    return std::sqrt(pos(1)*pos(1) + pos(2)*pos(2) + pos(3)*pos(3));
}

inline double getTheta(const Vec4& pos) {
    double r = getRadius(pos);
    if (r < 1e-10) return 0.0;
    return std::acos(pos(3) / r);
}

inline double getPhi(const Vec4& pos) {
    return std::atan2(pos(2), pos(1));
}

// Test fixture providing Schwarzschild, Kerr, and Minkowski metrics
class RenderPipelineTests : public ::testing::Test {
protected:
    Sirius::KerrSchildFamily schwarzschild{Sirius::KerrSchildParams::Schwarzschild(M)};
    Sirius::KerrSchildFamily kerr{Sirius::KerrSchildParams::Kerr(M, 0.5)};
    Sirius::KerrSchildFamily minkowski{Sirius::KerrSchildParams::Minkowski()};

    void SetUp() override {
        schwarzschild.setParameter("mass", M);
        kerr.setParameter("mass", M);
        kerr.setParameter("spin", 0.5);
    }

    // Create lightray with spherical position/velocity converted to Cartesian
    Lightray createRay(double t, double r, double theta, double phi,
                       double ut, double ur, double utheta, double uphi,
                       IMetric* metric) {
        Lightray ray;
        
        // Position: spherical (t,r,θ,φ) → Cartesian (t,x,y,z)
        ray.position(0) = t;
        ray.position(1) = r * std::sin(theta) * std::cos(phi);
        ray.position(2) = r * std::sin(theta) * std::sin(phi);
        ray.position(3) = r * std::cos(theta);
        
        // Velocity: spherical basis → Cartesian basis
        double st = std::sin(theta);
        double ct = std::cos(theta);
        double sp = std::sin(phi);
        double cp = std::cos(phi);
        
        ray.velocity(0) = ut;
        ray.velocity(1) = ur * st * cp + r * utheta * ct * cp - r * uphi * st * sp;
        ray.velocity(2) = ur * st * sp + r * utheta * ct * sp + r * uphi * st * cp;
        ray.velocity(3) = ur * ct - r * utheta * st;
        
        ray.proper_time = 0.0f;
        ray.coordinate_time = 0.0f;
        ray.terminated = 0;
        ray.step_size = 0.05f;
        ray.bounce_count = 0;
        ray.padding = 0;
        ray.ku_uobsu = 1.0f;
        ray.running_dlambda_dnew = 1.0f;
        
        ray.acceleration = Geodesic::calculateAcceleration(ray.velocity, ray.position, metric);
        
        return ray;
    }
};

// --- Lightray Structure Tests ---

// Verify lightray initialization with correct position and velocity
TEST_F(RenderPipelineTests, LightrayInitialization) {
    Lightray ray = createRay(0.0, 10.0, PI / 2.0, 0.0, 1.0, 0.5, 0.0, 0.1, &schwarzschild);
    
    EXPECT_DOUBLE_EQ(getRadius(ray.position), 10.0);
    EXPECT_DOUBLE_EQ(ray.velocity(0), 1.0);
    EXPECT_EQ(ray.terminated, 0);
}

// Verify all position and velocity components are finite
TEST_F(RenderPipelineTests, LightrayPositionFinite) {
    Lightray ray = createRay(0.0, 20.0, PI/2.0, PI, 1.0, 0.5, 0.0, 0.1, &schwarzschild);
    
    for (int mu = 0; mu < 4; mu++) {
        EXPECT_TRUE(std::isfinite(ray.position(mu)));
        EXPECT_TRUE(std::isfinite(ray.velocity(mu)));
    }
}

// --- Geodesic Acceleration Tests ---

// Verify acceleration calculation returns finite values
TEST_F(RenderPipelineTests, AccelerationIsFinite) {
    Lightray ray = createRay(0.0, 10.0, PI / 2.0, PI / 4.0, 1.0, 0.5, 0.1, 0.2, &schwarzschild);
    
    for (int mu = 0; mu < 4; mu++) {
        EXPECT_TRUE(std::isfinite(ray.acceleration(mu))) 
            << "Acceleration component " << mu << " should be finite";
    }
}

// Verify gravitational acceleration increases closer to horizon
TEST_F(RenderPipelineTests, SchwarzschildAccelerationNearHorizon) {
    Lightray ray_far = createRay(0.0, 10.0, PI / 2.0, 0.0, 1.0, -0.5, 0.0, 0.0, &schwarzschild);
    Lightray ray_near = createRay(0.0, 3.0, PI / 2.0, 0.0, 1.0, -0.5, 0.0, 0.0, &schwarzschild);
    
    double mag_far = std::abs(ray_far.acceleration(1));
    double mag_near = std::abs(ray_near.acceleration(1));
    
    EXPECT_GT(mag_near, mag_far) 
        << "Acceleration should be stronger closer to black hole";
}

// --- Integration Step Tests ---

// Verify integration step produces valid state with finite coordinates
TEST_F(RenderPipelineTests, IntegrateStepProducesValidState) {
    Lightray ray = createRay(0.0, 15.0, PI / 2.0, 0.0, 1.0, 0.5, 0.0, 0.1, &schwarzschild);
    
    bool success = Geodesic::integrateStep(ray, &schwarzschild);
    
    EXPECT_TRUE(success) << "Integration step should succeed";
    EXPECT_TRUE(std::isfinite(getRadius(ray.position))) << "r should be finite";
    EXPECT_TRUE(std::isfinite(getTheta(ray.position))) << "theta should be finite";
}

// Verify multiple integration steps progress ray position
TEST_F(RenderPipelineTests, MultipleStepsProgressRay) {
    Lightray ray = createRay(0.0, 20.0, PI / 2.0, 0.0, 1.0, 1.0, 0.0, 0.0, &schwarzschild);
    
    double initial_r = getRadius(ray.position);
    
    for (int i = 0; i < 10; i++) {
        bool success = Geodesic::integrateStep(ray, &schwarzschild);
        if (!success || ray.terminated) break;
    }
    
    EXPECT_GT(getRadius(ray.position), initial_r) 
        << "Outward ray should increase r";
}

// --- RK45 Integration Tests ---

// Verify RK45 integrator completes step successfully
TEST_F(RenderPipelineTests, RK45StepSucceeds) {
    Lightray ray = createRay(0.0, 15.0, PI / 2.0, PI / 4.0, 1.0, 0.3, 0.1, -0.1, &schwarzschild);
    
    IntegratorConfig config = Geodesic::getDefaultConfig();
    config.use_rk45 = true;
    
    bool success = Geodesic::integrateStepRK45(ray, &schwarzschild, config);
    
    EXPECT_TRUE(success) << "RK45 step should succeed";
}

// Verify RK45 integration maintains finite values throughout trajectory
TEST_F(RenderPipelineTests, RK45ProducesFiniteValues) {
    Lightray ray = createRay(0.0, 10.0, PI / 3.0, PI / 6.0, 1.0, -0.3, 0.2, 0.1, &schwarzschild);
    
    IntegratorConfig config = Geodesic::getDefaultConfig();
    
    for (int i = 0; i < 20; i++) {
        bool success = Geodesic::integrateStepRK45(ray, &schwarzschild, config);
        if (!success || ray.terminated) break;
        
        for (int mu = 0; mu < 4; mu++) {
            EXPECT_TRUE(std::isfinite(ray.position(mu)));
            EXPECT_TRUE(std::isfinite(ray.velocity(mu)));
        }
    }
}

// --- Cross-Metric Consistency Tests ---

// Verify different metrics produce different trajectories from same initial conditions
TEST_F(RenderPipelineTests, DifferentMetricsGiveDifferentPaths) {
    Lightray ray_mink = createRay(0.0, 20.0, PI / 2.0, 0.0, 1.0, -0.5, 0.0, 0.1, &minkowski);
    Lightray ray_schw = createRay(0.0, 20.0, PI / 2.0, 0.0, 1.0, -0.5, 0.0, 0.1, &schwarzschild);
    
    for (int i = 0; i < 50; i++) {
        Geodesic::integrateStep(ray_mink, &minkowski);
        Geodesic::integrateStep(ray_schw, &schwarzschild);
        if (ray_mink.terminated || ray_schw.terminated) break;
    }
    
    double diff = std::abs(getRadius(ray_mink.position) - getRadius(ray_schw.position));
    EXPECT_GT(diff, 0.01) 
        << "Different metrics should produce different paths";
}

// --- Edge Case Tests ---

// Verify integration handles θ near poles without singularities
TEST_F(RenderPipelineTests, HandlesNearPoleTheta) {
    Lightray ray = createRay(0.0, 10.0, 0.01, 0.0, 1.0, 0.5, 0.1, 0.0, &schwarzschild);
    
    for (int i = 0; i < 10; i++) {
        bool success = Geodesic::integrateStep(ray, &schwarzschild);
        if (!success || ray.terminated) break;
        
        EXPECT_TRUE(std::isfinite(getTheta(ray.position))) << "Theta should be finite";
    }
}

// Verify integration handles large radius (weak-field limit)
TEST_F(RenderPipelineTests, HandlesLargeRadius) {
    Lightray ray = createRay(0.0, 100.0, PI / 2.0, 0.0, 1.0, 1.0, 0.0, 0.0, &schwarzschild);
    ray.step_size = 0.5f;
    
    bool success = Geodesic::integrateStep(ray, &schwarzschild);
    
    if (success) {
        EXPECT_TRUE(std::isfinite(getRadius(ray.position))) << "r should be finite if step succeeds";
    }
    EXPECT_TRUE(std::isfinite(getRadius(ray.position))) << "Position should be finite";
}

// Verify integration is deterministic with identical initial conditions
TEST_F(RenderPipelineTests, IntegrationIsDeterministic) {
    Lightray ray1 = createRay(0.0, 15.0, PI / 2.0, 0.0, 1.0, 0.5, 0.1, 0.2, &schwarzschild);
    Lightray ray2 = createRay(0.0, 15.0, PI / 2.0, 0.0, 1.0, 0.5, 0.1, 0.2, &schwarzschild);
    
    for (int i = 0; i < 20; i++) {
        Geodesic::integrateStep(ray1, &schwarzschild);
        Geodesic::integrateStep(ray2, &schwarzschild);
        if (ray1.terminated || ray2.terminated) break;
    }
    
    for (int mu = 0; mu < 4; mu++) {
        EXPECT_DOUBLE_EQ(ray1.position(mu), ray2.position(mu))
            << "Same inputs should give same outputs";
    }
}

// Verify Kerr frame-dragging causes φ motion for radially infalling rays
TEST_F(RenderPipelineTests, KerrFrameDragging) {
    Lightray ray = createRay(0.0, 5.0, PI / 2.0, 0.0, 1.0, -0.5, 0.0, 0.0, &kerr);
    
    for (int i = 0; i < 30; i++) {
        Geodesic::integrateStep(ray, &kerr);
        if (ray.terminated) break;
    }
    
    EXPECT_NE(getPhi(ray.position), 0.0) 
        << "Frame dragging should cause phi to change";
}

} // namespace sirius::test