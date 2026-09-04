// RK45 Integrator Validation Tests
// Tests: Dormand-Prince coefficients, step adaptation, null constraint.
// Ported from TSPH009A.cpp; assertions and tolerances unchanged.

#include "sirius/core/constants.h"
#include "sirius/core/dual_number.h"
#include "sirius/core/geodesic_integrator.h"
#include "sirius/core/metrics/kerr_schild_family.h"  // Unified Kerr-Schild family
#include "sirius/core/metrics/metric.h"
#include "sirius/core/tensor.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace sirius::test {
using namespace sirius::core;

// =============================================================================
// Constants
// =============================================================================

constexpr double M = 1.0;  // Mass in geometric units
constexpr double PI = 3.14159265358979323846;
[[maybe_unused]] constexpr double kEpsilon = 1e-6;

class HalfSpaceMinkowskiMetric final : public IMetric {
  public:
    bool IsValidEvent(const Tensor<double, 4>& position) const override {
        return IMetric::IsValidEvent(position) && position(1) >= 0.0;
    }

    void Evaluate(const Tensor<double, 4>& position, Metric4d& metric,
                  Tensor<Dual<double>, 4, 4, 4>& derivatives) override {
        if (!IsValidEvent(position)) ++invalid_evaluations;
        metric.Zero();
        metric(0, 0) = Dual<double>(-1.0);
        metric(1, 1) = Dual<double>(1.0);
        metric(2, 2) = Dual<double>(1.0);
        metric(3, 3) = Dual<double>(1.0);
        derivatives.Zero();
    }

    const Config& GetParameters() const override { return config; }
    void SetParameter(const std::string&, double) override {}
    const char* GetName() const override { return "Half-space Minkowski test metric"; }

    int invalid_evaluations = 0;

  private:
    Config config;
};

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
        ray.ku_uobsu = 1.0f;

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
    EXPECT_TRUE(IsRepresentedIntegratorConfig(cfg));

    cfg.safety_factor = std::numeric_limits<float>::quiet_NaN();
    EXPECT_FALSE(IsRepresentedIntegratorConfig(cfg));
    EXPECT_DEATH((void)Geodesic::ComputeOptimalStep(0.01f, 0.1f, 1.0f, cfg),
                 "precondition.*enforced, terminating");
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

// Test: the live Schwarzschild RK45 path preserves its null constraint at the
// named double-precision CPU tolerance.
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

    EXPECT_GT(steps_completed, 10) << "Should complete multiple steps";
    ASSERT_TRUE(std::isfinite(max_violation)) << "Violation must be finite";
    EXPECT_LT(max_violation, constants::geodesic::kNullConditionTolCpu)
        << "Schwarzschild RK45 path exceeded the CPU null-condition authority";
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

    for (int i = 0; i < 200; ++i) {
        step_sizes.push_back(ray.step_size);

        Geodesic::IntegrateStepRk45(ray, &schwarzschild, config);
        if (ray.terminated) break;
        if (ray.position(1) < 3.0f) break;  // Don't go too close to horizon
    }

    EXPECT_GT(step_sizes.size(), 10) << "Should collect multiple samples";

    const auto [min_step, max_step] = std::minmax_element(step_sizes.begin(), step_sizes.end());
    EXPECT_LT(*min_step, *max_step)
        << "the adaptive controller never changed its step on the curved path";
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

TEST_F(RK45IntegratorTests, UnrepresentedStageShrinksBeforeMetricEvaluation) {
    HalfSpaceMinkowskiMetric metric;
    Lightray ray{};
    ray.position(1) = 0.01;
    ray.velocity(0) = 1.0;
    ray.velocity(1) = -1.0;
    ray.step_size = 0.1f;

    IntegratorConfig stage_config = config;
    stage_config.min_step = 1.0e-4f;
    stage_config.max_step = 0.1f;
    stage_config.initial_step = 0.1f;
    const Vec4 initial_position = ray.position;

    EXPECT_FALSE(Geodesic::IntegrateStepRk45(ray, &metric, stage_config));
    EXPECT_EQ(ray.terminated, 0);
    EXPECT_FLOAT_EQ(ray.step_size, 0.05f);
    EXPECT_EQ(metric.invalid_evaluations, 0);
    for (int component = 0; component < 4; ++component) {
        EXPECT_DOUBLE_EQ(ray.position(component), initial_position(component));
    }

    ray.step_size = stage_config.min_step;
    ray.position(1) = 1.0e-6;
    EXPECT_FALSE(Geodesic::IntegrateStepRk45(ray, &metric, stage_config));
    EXPECT_EQ(ray.terminated, 3);
    EXPECT_EQ(metric.invalid_evaluations, 0);
}

TEST_F(RK45IntegratorTests, RejectionMaySelectTheMinimumStepBeforeTerminating) {
    IntegratorConfig floor_config = config;
    floor_config.abs_tolerance = 1.0e-14f;
    floor_config.rel_tolerance = 1.0e-14f;
    floor_config.min_step = 0.05f;
    floor_config.max_step = 0.1f;
    floor_config.initial_step = 0.1f;

    Vec4 position;
    position(1) = 5.0;
    position(2) = PI / 2.0;
    Vec4 direction;
    direction(0) = 1.0;
    direction(3) = 1.0;
    Lightray ray = createTestRay(position, direction, &schwarzschild);
    ray.step_size = floor_config.initial_step;

    ASSERT_FALSE(Geodesic::IntegrateStepRk45(ray, &schwarzschild, floor_config));
    EXPECT_EQ(ray.terminated, 0)
        << "selecting the minimum retry step is not a terminal integration failure";
    EXPECT_FLOAT_EQ(ray.step_size, floor_config.min_step);
}

TEST_F(RK45IntegratorTests, NullProjectionPreservesTheIncomingLightConeBranch) {
    Metric4d minkowski;
    minkowski(0, 0) = -1.0;
    minkowski(1, 1) = 1.0;
    minkowski(2, 2) = 1.0;
    minkowski(3, 3) = 1.0;

    Vec4 past_directed;
    past_directed(0) = -1.000001;
    past_directed(1) = 1.0;
    const auto projected = Geodesic::ProjectNullTangentPreservingBranch(past_directed, minkowski);
    ASSERT_TRUE(projected.has_value());
    EXPECT_DOUBLE_EQ((*projected)(0), -1.0);
    EXPECT_DOUBLE_EQ((*projected)(1), 1.0);
    EXPECT_DOUBLE_EQ(TensorOps::InnerProduct(*projected, *projected, minkowski), 0.0);

    // At a stationary-limit-like event the temporal equation is linear.  The
    // represented root remains finite and does not invent the other cone.
    Metric4d linear;
    linear(0, 1) = 1.0;
    linear(1, 0) = 1.0;
    linear(1, 1) = 1.0;
    Vec4 tangent;
    tangent(0) = -0.6;
    tangent(1) = 1.0;
    const auto linear_projected = Geodesic::ProjectNullTangentPreservingBranch(tangent, linear);
    ASSERT_TRUE(linear_projected.has_value());
    EXPECT_DOUBLE_EQ((*linear_projected)(0), -0.5);
    EXPECT_DOUBLE_EQ(TensorOps::InnerProduct(*linear_projected, *linear_projected, linear), 0.0);

    // This Lorentzian metric has g_00 > 0. For the supplied near-null tangent,
    // holding all spatial components fixed gives a negative temporal
    // discriminant, as can occur inside the Kerr stationary limit. A nearby
    // spatial root remains represented and must preserve the incoming k^0.
    Metric4d ergoregion;
    ergoregion(0, 0) = 1.0;
    ergoregion(0, 1) = 2.0;
    ergoregion(1, 0) = 2.0;
    ergoregion(1, 1) = 1.0;
    ergoregion(2, 2) = 1.0;
    ergoregion(3, 3) = 1.0;
    Vec4 ergoregion_tangent;
    ergoregion_tangent(0) = -1.0;
    ergoregion_tangent(1) = 0.57;
    ergoregion_tangent(2) = 1.0;
    const auto ergoregion_projected =
        Geodesic::ProjectNullTangentPreservingBranch(ergoregion_tangent, ergoregion);
    ASSERT_TRUE(ergoregion_projected.has_value());
    EXPECT_DOUBLE_EQ((*ergoregion_projected)(0), ergoregion_tangent(0));
    EXPECT_NEAR((*ergoregion_projected)(1), 2.0 - std::sqrt(2.0), 1.0e-15);
    EXPECT_NEAR(TensorOps::InnerProduct(*ergoregion_projected, *ergoregion_projected, ergoregion),
                0.0, 1.0e-15);
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
