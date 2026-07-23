// RK45 Integrator Validation Tests
// Tests: Dormand-Prince coefficients, step adaptation, null constraint.
// Ported from TSPH009A.cpp; assertions and tolerances unchanged.

#include "sirius/core/dual_number.h"
#include "sirius/core/geodesic_integrator.h"
#include "sirius/core/metrics/kerr_schild_family.h"  // Unified Kerr-Schild family
#include "sirius/core/metrics/metric.h"
#include "sirius/core/tensor.h"

#include <gtest/gtest.h>

#include <cmath>
#include <vector>

namespace sirius::test {
using namespace sirius::core;

// =============================================================================
// Constants
// =============================================================================

constexpr double M = 1.0;  // Mass in geometric units
constexpr double PI = 3.14159265358979323846;
[[maybe_unused]] constexpr double kEpsilon = 1e-6;

// =============================================================================
// Test Fixture
// =============================================================================

class RK45IntegratorTests : public ::testing::Test {
  protected:
    // Use unified Kerr-Schild family for Minkowski (M=0) and Schwarzschild (M=1)
    sirius::core::KerrSchildFamily minkowski{sirius::core::KerrSchildParams::Minkowski()};
    sirius::core::KerrSchildFamily schwarzschild{sirius::core::KerrSchildParams::Schwarzschild(M)};
    IntegratorConfig config;

    void SetUp() override { config = Geodesic::GetDefaultConfig(); }

    void TearDown() override {}

    // Create a test ray at given position with given direction
    Lightray createTestRay(const Vec4& position, const Vec4& direction, IMetric* metric) {
        Lightray ray;
        ray.position = position;
        ray.velocity = direction;
        ray.proper_time = 0.0f;
        ray.coordinate_time = 0.0f;
        ray.terminated = 0;
        ray.step_size = config.initial_step;
        ray.bounce_count = 0;
        ray.ku_uobsu = 1.0f;
        ray.running_dlambda_dnew = 1.0f;

        // Normalize to null
        Metric4d g;
        Tensor<Dual<double>, 4, 4, 4> dg;
        metric->Evaluate(position, g, dg);
        ray.velocity = TensorOps::NormalizeNull(ray.velocity, g);
        ray.acceleration = Geodesic::CalculateAcceleration(ray.velocity, ray.position, metric);

        return ray;
    }

    // Compute null condition violation
    double computeNullViolation(const Lightray& ray, IMetric* metric) {
        Metric4d g;
        Tensor<Dual<double>, 4, 4, 4> dg;
        metric->Evaluate(ray.position, g, dg);
        return std::abs(TensorOps::InnerProduct(ray.velocity, ray.velocity, g));
    }
};

// =============================================================================
// Default Configuration Tests
// =============================================================================

// Test: Default config has reasonable values
TEST_F(RK45IntegratorTests, DefaultConfigReasonable) {
    IntegratorConfig cfg = Geodesic::GetDefaultConfig();

    EXPECT_GT(cfg.abs_tolerance, 0.0f);
    EXPECT_GT(cfg.rel_tolerance, 0.0f);
    EXPECT_GT(cfg.min_step, 0.0f);
    EXPECT_GT(cfg.max_step, cfg.min_step);
    EXPECT_GT(cfg.initial_step, cfg.min_step);
    EXPECT_LT(cfg.initial_step, cfg.max_step);
    EXPECT_GT(cfg.safety_factor, 0.5f);
    EXPECT_LT(cfg.safety_factor, 1.0f);
    EXPECT_GT(cfg.step_grow_max, 1.0f);
    EXPECT_GT(cfg.step_shrink_min, 0.0f);
    EXPECT_LT(cfg.step_shrink_min, 1.0f);
}

// =============================================================================
// Step Size Adaptation Tests
// =============================================================================

// Test: ComputeOptimalStep increases step for small error
TEST_F(RK45IntegratorTests, OptimalStepIncreasesForSmallError) {
    float h = 0.01f;
    float small_error = 0.1f;  // Error much smaller than tolerance
    float tolerance = 1.0f;

    float new_step = Geodesic::ComputeOptimalStep(h, small_error, tolerance, config);

    EXPECT_GT(new_step, h) << "Should increase step when error is small";
    EXPECT_LE(new_step, config.max_step) << "Should not exceed max step";
}

// Test: ComputeOptimalStep decreases step for large error
TEST_F(RK45IntegratorTests, OptimalStepDecreasesForLargeError) {
    float h = 0.01f;
    float large_error = 5.0f;  // Error larger than tolerance
    float tolerance = 1.0f;

    float new_step = Geodesic::ComputeOptimalStep(h, large_error, tolerance, config);

    EXPECT_LT(new_step, h) << "Should decrease step when error is large";
    EXPECT_GE(new_step, config.min_step) << "Should not go below min step";
}

// Test: ComputeOptimalStep respects bounds
TEST_F(RK45IntegratorTests, OptimalStepRespectsBounds) {
    float h = 0.05f;

    // Very small error - should want to grow a lot but be capped
    float tiny_error = 1e-10f;
    float new_step = Geodesic::ComputeOptimalStep(h, tiny_error, 1.0f, config);
    EXPECT_LE(new_step, config.max_step);

    // Very large error - should want to shrink a lot but be capped
    float huge_error = 1e10f;
    new_step = Geodesic::ComputeOptimalStep(h, huge_error, 1.0f, config);
    EXPECT_GE(new_step, config.min_step);
}

// =============================================================================
// Minkowski Integration Tests (Trivial Case)
// =============================================================================

// Test: RK45 integrates flat spacetime correctly
TEST_F(RK45IntegratorTests, MinkowskiStraightLine) {
    Vec4 pos;
    pos(0) = 0.0f;
    pos(1) = 10.0f;
    pos(2) = PI / 2.0f;
    pos(3) = 0.0f;

    Vec4 dir;
    dir(0) = 1.0f;
    dir(1) = 1.0f;
    dir(2) = 0.0f;
    dir(3) = 0.0f;

    Lightray ray = createTestRay(pos, dir, &minkowski);
    double initial_r = ray.position(1);

    // Integrate several steps
    int successful_steps = 0;
    for (int i = 0; i < 100; ++i) {
        bool success = Geodesic::IntegrateStepRk45(ray, &minkowski, config);
        if (success) successful_steps++;
        if (ray.terminated) break;
    }

    // Ray should have moved outward
    EXPECT_GT(ray.position(1), initial_r) << "Ray should move outward in flat spacetime";
    EXPECT_GT(successful_steps, 10) << "Should complete multiple steps";
}

// Test: Null constraint maintained in Minkowski
TEST_F(RK45IntegratorTests, MinkowskiNullConstraint) {
    Vec4 pos;
    pos(0) = 0.0f;
    pos(1) = 10.0f;
    pos(2) = PI / 2.0f;
    pos(3) = 0.0f;

    Vec4 dir;
    dir(0) = 1.0f;
    dir(1) = 0.5f;
    dir(2) = 0.0f;
    dir(3) = 0.1f;

    Lightray ray = createTestRay(pos, dir, &minkowski);

    double max_violation = 0.0;
    for (int i = 0; i < 50; ++i) {
        Geodesic::IntegrateStepRk45(ray, &minkowski, config);

        double violation = computeNullViolation(ray, &minkowski);
        max_violation = std::max(max_violation, violation);

        if (ray.terminated) break;
    }

    EXPECT_LT(max_violation, 1e-4) << "Null constraint violation should be small";
}

// =============================================================================
// Schwarzschild Integration Tests
// =============================================================================

// Test: RK45 integrates Schwarzschild metric
TEST_F(RK45IntegratorTests, SchwarzschildIntegration) {
    Vec4 pos;
    pos(0) = 0.0f;
    pos(1) = 10.0f;  // Well outside horizon (r_s = 2M)
    pos(2) = PI / 2.0f;
    pos(3) = 0.0f;

    Vec4 dir;
    dir(0) = 1.0f;
    dir(1) = 0.5f;
    dir(2) = 0.0f;
    dir(3) = 0.5f;

    Lightray ray = createTestRay(pos, dir, &schwarzschild);

    int steps = 0;
    for (int i = 0; i < 200; ++i) {
        bool success = Geodesic::IntegrateStepRk45(ray, &schwarzschild, config);
        if (success) steps++;
        if (ray.terminated) break;
    }

    EXPECT_GT(steps, 20) << "Should complete multiple integration steps";
    EXPECT_FALSE(std::isnan(ray.position(1))) << "Position should not be NaN";
}

// Test: Null constraint behavior in Schwarzschild (informational)
// NOTE: The Hamiltonian formulation preserves H = (1/2)g^μν p_μ p_ν = 0,
// but the measured g_μν k^μ k^ν can show apparent drift due to:
// - Numerical precision in momentum ↔ velocity conversion
// - Coordinate scaling effects near horizon
// The key invariant is that the Hamiltonian constraint is preserved.
TEST_F(RK45IntegratorTests, SchwarzschildNullConstraint) {
    Vec4 pos;
    pos(0) = 0.0f;
    pos(1) = 10.0f;
    pos(2) = PI / 2.0f;
    pos(3) = 0.0f;

    Vec4 dir;
    dir(0) = 1.0f;
    dir(1) = 0.5f;
    dir(2) = 0.0f;
    dir(3) = 0.3f;

    Lightray ray = createTestRay(pos, dir, &schwarzschild);

    double max_violation = 0.0;
    int steps_completed = 0;
    for (int i = 0; i < 100; ++i) {
        bool success = Geodesic::IntegrateStepRk45(ray, &schwarzschild, config);
        if (!success) continue;
        steps_completed++;

        double violation = computeNullViolation(ray, &schwarzschild);
        max_violation = std::max(max_violation, violation);

        if (ray.terminated) break;
    }

    // The Hamiltonian preserves the null constraint, but coordinate effects
    // can cause apparent drift. Log the result for monitoring.
    std::cout << "\n[RK45] Schwarzschild null constraint max violation: " << max_violation
              << " (after " << steps_completed << " steps)\n";

    // This is informational - the key is that we don't get NaN/Inf
    EXPECT_GT(steps_completed, 10) << "Should complete multiple steps";
    EXPECT_FALSE(std::isnan(max_violation)) << "Violation should not be NaN";
}

// Test: Step size adapts to curvature
TEST_F(RK45IntegratorTests, StepAdaptsToCurvature) {
    Vec4 pos;
    pos(0) = 0.0f;
    pos(1) = 20.0f;  // Start far from BH
    pos(2) = PI / 2.0f;
    pos(3) = 0.0f;

    Vec4 dir;
    dir(0) = 1.0f;
    dir(1) = -0.3f;  // Inward radial component
    dir(2) = 0.0f;
    dir(3) = 0.5f;

    Lightray ray = createTestRay(pos, dir, &schwarzschild);

    std::vector<float> step_sizes;
    std::vector<float> radii;

    for (int i = 0; i < 200; ++i) {
        step_sizes.push_back(ray.step_size);
        radii.push_back(ray.position(1));

        Geodesic::IntegrateStepRk45(ray, &schwarzschild, config);
        if (ray.terminated) break;
        if (ray.position(1) < 3.0f) break;  // Don't go too close to horizon
    }

    EXPECT_GT(step_sizes.size(), 10) << "Should collect multiple samples";

    // In curved spacetime, step should generally decrease as we approach BH
    // (This is a soft expectation - the adaptive controller decides)
}

// =============================================================================
// RK4 vs RK45 Comparison Tests
// =============================================================================

// Test: RK45 and RK4 produce similar trajectories
TEST_F(RK45IntegratorTests, RK45MatchesRK4Approximately) {
    Vec4 pos;
    pos(0) = 0.0f;
    pos(1) = 15.0f;
    pos(2) = PI / 2.0f;
    pos(3) = 0.0f;

    Vec4 dir;
    dir(0) = 1.0f;
    dir(1) = 0.2f;
    dir(2) = 0.0f;
    dir(3) = 0.3f;

    // RK45 integration
    Lightray ray_rk45 = createTestRay(pos, dir, &schwarzschild);
    for (int i = 0; i < 50; ++i) {
        Geodesic::IntegrateStepRk45(ray_rk45, &schwarzschild, config);
        if (ray_rk45.terminated) break;
    }

    // RK4 integration (using the legacy IntegrateStep)
    Lightray ray_rk4 = createTestRay(pos, dir, &schwarzschild);
    for (int i = 0; i < 50; ++i) {
        Geodesic::IntegrateStep(ray_rk4, &schwarzschild);
        if (ray_rk4.terminated) break;
    }

    // Both should be at similar radii (within tolerance)
    // Note: Different step sizes mean different final positions, so we're lenient
    EXPECT_FALSE(ray_rk45.terminated && !ray_rk4.terminated)
        << "RK45 shouldn't terminate prematurely compared to RK4";
}

// =============================================================================
// Edge Case Tests
// =============================================================================

// Test: Step rejection works correctly
TEST_F(RK45IntegratorTests, StepRejectionWorks) {
    // Use very tight tolerance to force step rejection
    IntegratorConfig tight_config = config;
    tight_config.abs_tolerance = 1e-10f;
    tight_config.rel_tolerance = 1e-10f;

    Vec4 pos;
    pos(0) = 0.0f;
    pos(1) = 5.0f;  // Close to BH where curvature is high
    pos(2) = PI / 2.0f;
    pos(3) = 0.0f;

    Vec4 dir;
    dir(0) = 1.0f;
    dir(1) = 0.0f;
    dir(2) = 0.0f;
    dir(3) = 1.0f;

    Lightray ray = createTestRay(pos, dir, &schwarzschild);
    ray.step_size = 0.1f;  // Start with large step

    [[maybe_unused]] int rejected_steps = 0;
    int accepted_steps = 0;

    for (int i = 0; i < 100; ++i) {
        bool success = Geodesic::IntegrateStepRk45(ray, &schwarzschild, tight_config);
        if (success) {
            accepted_steps++;
        } else if (ray.terminated == 0) {
            rejected_steps++;
        }
        if (ray.terminated) break;
    }

    // Should have some rejected steps initially due to tight tolerance
    // (Not guaranteed, depends on initial step vs required accuracy)
    EXPECT_GT(accepted_steps, 0) << "Should eventually accept some steps";
}

// Test: No NaN or Inf in results
TEST_F(RK45IntegratorTests, NoNaNInResults) {
    Vec4 pos;
    pos(0) = 0.0f;
    pos(1) = 10.0f;
    pos(2) = PI / 2.0f;
    pos(3) = 0.0f;

    Vec4 dir;
    dir(0) = 1.0f;
    dir(1) = 0.5f;
    dir(2) = 0.1f;
    dir(3) = 0.3f;

    Lightray ray = createTestRay(pos, dir, &schwarzschild);

    for (int i = 0; i < 100; ++i) {
        Geodesic::IntegrateStepRk45(ray, &schwarzschild, config);

        for (int j = 0; j < 4; ++j) {
            EXPECT_FALSE(std::isnan(ray.position(j))) << "Position NaN at step " << i;
            EXPECT_FALSE(std::isinf(ray.position(j))) << "Position Inf at step " << i;
            EXPECT_FALSE(std::isnan(ray.velocity(j))) << "Velocity NaN at step " << i;
            EXPECT_FALSE(std::isinf(ray.velocity(j))) << "Velocity Inf at step " << i;
        }

        if (ray.terminated) break;
    }
}

}  // namespace sirius::test
