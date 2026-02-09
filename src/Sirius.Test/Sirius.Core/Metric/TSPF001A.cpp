// TSPF001A.cpp - FPS Threshold Performance Tests
// Tests: frame budget validation, metric evaluation time, per-ray budgets.

#define _USE_MATH_DEFINES
#include <cmath>
#include <gtest/gtest.h>
#include <chrono>
#include <vector>
#include <numeric>
#include "MTTN001A.h"
#include "MTDL001A.h"
#include "PHMT000A.h"
#include "PHMT100A.h"

namespace sirius::test {

// =============================================================================
// Performance Targets (from DevelopmentRoadmap.md)
// =============================================================================

namespace FPSThresholds {
    // Resolution
    constexpr int WIDTH_1080P = 1920;
    constexpr int HEIGHT_1080P = 1080;
    constexpr int TOTAL_PIXELS_1080P = WIDTH_1080P * HEIGHT_1080P;  // 2,073,600
    
    // FPS targets and bounds
    constexpr double MINKOWSKI_TARGET_FPS = 120.0;
    constexpr double MINKOWSKI_BOUND_FPS = 60.0;
    constexpr double SCHWARZSCHILD_TARGET_FPS = 60.0;
    constexpr double SCHWARZSCHILD_BOUND_FPS = 30.0;
    constexpr double KERR_TARGET_FPS = 30.0;
    constexpr double KERR_BOUND_FPS = 15.0;
    constexpr double NUMERICAL_TARGET_FPS = 10.0;
    constexpr double NUMERICAL_BOUND_FPS = 5.0;
    
    // Frame time budgets (ms)
    constexpr double MINKOWSKI_FRAME_BUDGET_MS = 1000.0 / MINKOWSKI_TARGET_FPS;      // 8.33ms
    constexpr double SCHWARZSCHILD_FRAME_BUDGET_MS = 1000.0 / SCHWARZSCHILD_TARGET_FPS; // 16.67ms
    constexpr double KERR_FRAME_BUDGET_MS = 1000.0 / KERR_TARGET_FPS;                // 33.33ms
    constexpr double NUMERICAL_FRAME_BUDGET_MS = 1000.0 / NUMERICAL_TARGET_FPS;      // 100ms
    
    // Integration parameters
    constexpr int TYPICAL_STEPS_PER_RAY = 500;
    constexpr double STEP_SIZE = 0.01;
}

// =============================================================================
// Test Fixture
// =============================================================================

class FPSThresholdTests : public ::testing::Test {
protected:
    static constexpr double M = 1.0;  // Mass in geometric units
    
    // Create standard test position (Cartesian x=10, y=0, z=0)
    Vec4 getStandardPosition() {
        Vec4 pos;
        pos(0) = 0.0;
        pos(1) = 10.0;
        pos(2) = 0.0;
        pos(3) = 0.0;
        return pos;
    }
    
    // Measure time to evaluate metric N times (in microseconds)
    double measureMetricEvaluationTime(IMetric* metric, int iterations) {
        using namespace std::chrono;
        
        Vec4 pos = getStandardPosition();
        Metric4D g;
        Tensor<Dual<double>, 4, 4, 4> dg;
        
        auto start = high_resolution_clock::now();
        
        for (int i = 0; i < iterations; ++i) {
            metric->evaluate(pos, g, dg);
        }
        
        auto end = high_resolution_clock::now();
        return duration_cast<nanoseconds>(end - start).count() / 1000.0;  // µs
    }
    
    // Calculate estimated FPS based on per-evaluation time
    struct FPSEstimate {
        double eval_time_us;        // Time per evaluation (µs)
        double evals_per_pixel;     // Evaluations per pixel (typically steps_per_ray)
        double frame_time_ms;       // Total frame time (ms)
        double estimated_fps;       // Estimated FPS
        bool meets_bound;           // Meets minimum bound?
        bool meets_target;          // Meets target?
    };
    
    FPSEstimate estimateFPS(double eval_time_us, int evals_per_pixel, 
                            double target_fps, double bound_fps) {
        FPSEstimate est;
        est.eval_time_us = eval_time_us;
        est.evals_per_pixel = evals_per_pixel;
        
        // Total evaluations per frame
        double total_evals = FPSThresholds::TOTAL_PIXELS_1080P * evals_per_pixel;
        
        // Total frame time in ms (assuming no parallelism for worst case)
        // Note: GPU parallelization reduces this by ~1000x or more
        est.frame_time_ms = (total_evals * eval_time_us) / 1000.0;
        
        // Estimated FPS (theoretical sequential)
        est.estimated_fps = 1000.0 / est.frame_time_ms;
        
        // Check thresholds
        est.meets_bound = est.estimated_fps >= bound_fps;
        est.meets_target = est.estimated_fps >= target_fps;
        
        return est;
    }
    
    void SetUp() override {}
    void TearDown() override {}
};

// =============================================================================
// Frame Budget Validation Tests
// =============================================================================

// Test: Verify frame budget calculations are consistent
TEST_F(FPSThresholdTests, FrameBudgetCalculations) {
    // Verify frame time = 1000 / FPS
    EXPECT_NEAR(FPSThresholds::MINKOWSKI_FRAME_BUDGET_MS, 
                1000.0 / FPSThresholds::MINKOWSKI_TARGET_FPS, 0.01);
    EXPECT_NEAR(FPSThresholds::SCHWARZSCHILD_FRAME_BUDGET_MS, 
                1000.0 / FPSThresholds::SCHWARZSCHILD_TARGET_FPS, 0.01);
    EXPECT_NEAR(FPSThresholds::KERR_FRAME_BUDGET_MS, 
                1000.0 / FPSThresholds::KERR_TARGET_FPS, 0.01);
    EXPECT_NEAR(FPSThresholds::NUMERICAL_FRAME_BUDGET_MS, 
                1000.0 / FPSThresholds::NUMERICAL_TARGET_FPS, 0.01);
}

// Test: Verify resolution calculations
TEST_F(FPSThresholdTests, ResolutionCalculations) {
    EXPECT_EQ(FPSThresholds::TOTAL_PIXELS_1080P, 1920 * 1080);
    EXPECT_EQ(FPSThresholds::TOTAL_PIXELS_1080P, 2073600);
}

// Test: Verify target > bound for all metrics
TEST_F(FPSThresholdTests, TargetsExceedBounds) {
    EXPECT_GT(FPSThresholds::MINKOWSKI_TARGET_FPS, FPSThresholds::MINKOWSKI_BOUND_FPS);
    EXPECT_GT(FPSThresholds::SCHWARZSCHILD_TARGET_FPS, FPSThresholds::SCHWARZSCHILD_BOUND_FPS);
    EXPECT_GT(FPSThresholds::KERR_TARGET_FPS, FPSThresholds::KERR_BOUND_FPS);
    EXPECT_GT(FPSThresholds::NUMERICAL_TARGET_FPS, FPSThresholds::NUMERICAL_BOUND_FPS);
}

// =============================================================================
// Metric Evaluation Time Tests
// =============================================================================

// Test: Minkowski metric evaluation time (should be fastest)
TEST_F(FPSThresholdTests, MinkowskiEvaluationTime) {
    Sirius::KerrSchildFamily mink{Sirius::KerrSchildParams::Minkowski()};
    
    int iterations = 10000;
    double total_time_us = measureMetricEvaluationTime(&mink, iterations);
    double per_eval_us = total_time_us / iterations;
    
    std::cout << "\n[PERF] Minkowski metric evaluation: " << per_eval_us << " µs\n";
    
    // Minkowski should be very fast (constant metric, no computation)
    EXPECT_LT(per_eval_us, 10.0) << "Minkowski should evaluate in < 10µs";
    
    // Report estimated FPS
    auto est = estimateFPS(per_eval_us, 1, 
                           FPSThresholds::MINKOWSKI_TARGET_FPS,
                           FPSThresholds::MINKOWSKI_BOUND_FPS);
    std::cout << "  Per-eval time: " << est.eval_time_us << " µs\n";
    std::cout << "  Frame time (sequential): " << est.frame_time_ms << " ms\n";
    std::cout << "  Note: GPU parallelism reduces this by ~1000x+\n";
}

// Test: Schwarzschild metric evaluation time
TEST_F(FPSThresholdTests, SchwarzschildEvaluationTime) {
    Sirius::KerrSchildFamily schw{Sirius::KerrSchildParams::Schwarzschild(M)};
    schw.setParameter("mass", M);
    
    int iterations = 10000;
    double total_time_us = measureMetricEvaluationTime(&schw, iterations);
    double per_eval_us = total_time_us / iterations;
    
    std::cout << "\n[PERF] Schwarzschild metric evaluation: " << per_eval_us << " µs\n";
    
    // Schwarzschild should be reasonably fast
    EXPECT_LT(per_eval_us, 50.0) << "Schwarzschild should evaluate in < 50µs";
}

// Test: Kerr metric evaluation time
TEST_F(FPSThresholdTests, KerrEvaluationTime) {
    Sirius::KerrSchildFamily kerr{Sirius::KerrSchildParams::Kerr(M, 0.9)};
    kerr.setParameter("mass", M);
    kerr.setParameter("spin", 0.9);
    
    int iterations = 10000;
    double total_time_us = measureMetricEvaluationTime(&kerr, iterations);
    double per_eval_us = total_time_us / iterations;
    
    std::cout << "\n[PERF] Kerr metric evaluation: " << per_eval_us << " µs\n";
    
    // Kerr is more complex but should still be reasonable
    EXPECT_LT(per_eval_us, 100.0) << "Kerr should evaluate in < 100µs";
}

// Test: Metric complexity ordering (Minkowski < Schwarzschild < Kerr)
TEST_F(FPSThresholdTests, MetricComplexityOrdering) {
    Sirius::KerrSchildFamily mink{Sirius::KerrSchildParams::Minkowski()};
    Sirius::KerrSchildFamily schw{Sirius::KerrSchildParams::Schwarzschild(M)};
    Sirius::KerrSchildFamily kerr{Sirius::KerrSchildParams::Kerr(M, 0.9)};
    
    schw.setParameter("mass", M);
    kerr.setParameter("mass", M);
    kerr.setParameter("spin", 0.9);
    
    int iterations = 5000;
    
    double mink_time = measureMetricEvaluationTime(&mink, iterations) / iterations;
    double schw_time = measureMetricEvaluationTime(&schw, iterations) / iterations;
    double kerr_time = measureMetricEvaluationTime(&kerr, iterations) / iterations;
    
    std::cout << "\n[PERF] Metric complexity ordering:\n";
    std::cout << "  Minkowski:     " << mink_time << " µs\n";
    std::cout << "  Schwarzschild: " << schw_time << " µs\n";
    std::cout << "  Kerr:          " << kerr_time << " µs\n";
    
    // Minkowski should be fastest (or at least as fast)
    EXPECT_LE(mink_time, schw_time + 1.0)  // Allow 1µs tolerance
        << "Minkowski should be at least as fast as Schwarzschild";
}

// =============================================================================
// Per-Ray Budget Tests
// =============================================================================

// Test: Budget per ray for Minkowski
TEST_F(FPSThresholdTests, MinkowskiPerRayBudget) {
    // Given target FPS and resolution, calculate allowed time per ray
    double frame_budget_ms = FPSThresholds::MINKOWSKI_FRAME_BUDGET_MS;
    double rays_per_frame = FPSThresholds::TOTAL_PIXELS_1080P;
    double budget_per_ray_us = (frame_budget_ms * 1000.0) / rays_per_frame;
    
    std::cout << "\n[BUDGET] Minkowski @ " << FPSThresholds::MINKOWSKI_TARGET_FPS << " FPS:\n";
    std::cout << "  Frame budget: " << frame_budget_ms << " ms\n";
    std::cout << "  Rays per frame: " << rays_per_frame << "\n";
    std::cout << "  Budget per ray: " << budget_per_ray_us << " µs\n";
    
    // At 120 FPS, 1080p, each ray gets ~4µs on a single thread
    // GPU parallelism makes this achievable
    EXPECT_GT(budget_per_ray_us, 0.0);
}

// Test: Budget per ray for Schwarzschild
TEST_F(FPSThresholdTests, SchwarzschildPerRayBudget) {
    double frame_budget_ms = FPSThresholds::SCHWARZSCHILD_FRAME_BUDGET_MS;
    double rays_per_frame = FPSThresholds::TOTAL_PIXELS_1080P;
    double budget_per_ray_us = (frame_budget_ms * 1000.0) / rays_per_frame;
    
    std::cout << "\n[BUDGET] Schwarzschild @ " << FPSThresholds::SCHWARZSCHILD_TARGET_FPS << " FPS:\n";
    std::cout << "  Frame budget: " << frame_budget_ms << " ms\n";
    std::cout << "  Budget per ray: " << budget_per_ray_us << " µs\n";
    
    // At 60 FPS, each ray gets ~8µs
    EXPECT_GT(budget_per_ray_us, 0.0);
}

// Test: Budget per ray for Kerr
TEST_F(FPSThresholdTests, KerrPerRayBudget) {
    double frame_budget_ms = FPSThresholds::KERR_FRAME_BUDGET_MS;
    double rays_per_frame = FPSThresholds::TOTAL_PIXELS_1080P;
    double budget_per_ray_us = (frame_budget_ms * 1000.0) / rays_per_frame;
    
    std::cout << "\n[BUDGET] Kerr @ " << FPSThresholds::KERR_TARGET_FPS << " FPS:\n";
    std::cout << "  Frame budget: " << frame_budget_ms << " ms\n";
    std::cout << "  Budget per ray: " << budget_per_ray_us << " µs\n";
    
    // At 30 FPS, each ray gets ~16µs
    EXPECT_GT(budget_per_ray_us, 0.0);
}

// =============================================================================
// Integration Step Budget Tests  
// =============================================================================

// Test: Per-step budget with typical ray length
TEST_F(FPSThresholdTests, PerStepBudgetAnalysis) {
    int steps_per_ray = FPSThresholds::TYPICAL_STEPS_PER_RAY;
    
    std::vector<std::pair<std::string, double>> metrics = {
        {"Minkowski", FPSThresholds::MINKOWSKI_TARGET_FPS},
        {"Schwarzschild", FPSThresholds::SCHWARZSCHILD_TARGET_FPS},
        {"Kerr", FPSThresholds::KERR_TARGET_FPS},
        {"Numerical", FPSThresholds::NUMERICAL_TARGET_FPS}
    };
    
    std::cout << "\n[BUDGET] Per-step budget @ " << steps_per_ray << " steps/ray:\n";
    
    for (const auto& [name, target_fps] : metrics) {
        double frame_budget_ms = 1000.0 / target_fps;
        double rays_per_frame = FPSThresholds::TOTAL_PIXELS_1080P;
        double budget_per_ray_us = (frame_budget_ms * 1000.0) / rays_per_frame;
        double budget_per_step_ns = (budget_per_ray_us * 1000.0) / steps_per_ray;
        
        std::cout << "  " << name << " @ " << target_fps << " FPS: "
                  << budget_per_step_ns << " ns/step\n";
        
        EXPECT_GT(budget_per_step_ns, 0.0);
    }
}

// =============================================================================
// GPU Parallelism Assumptions
// =============================================================================

// Test: Document GPU parallelism assumptions
TEST_F(FPSThresholdTests, GPUParallelismAssumptions) {
    // Modern GPUs have thousands of CUDA cores
    // RTX 3080 has 8704 CUDA cores
    // RTX 4090 has 16384 CUDA cores
    
    constexpr int TYPICAL_CUDA_CORES = 8000;
    constexpr double PARALLEL_EFFICIENCY = 0.7;  // Occupancy and overhead
    
    double effective_parallelism = TYPICAL_CUDA_CORES * PARALLEL_EFFICIENCY;
    
    // With parallelism, frame time is reduced
    double sequential_time_ms = 1000.0;  // Hypothetical 1 second sequential
    double parallel_time_ms = sequential_time_ms / effective_parallelism;
    
    std::cout << "\n[GPU] Parallelism assumptions:\n";
    std::cout << "  CUDA cores (typical): " << TYPICAL_CUDA_CORES << "\n";
    std::cout << "  Effective parallelism: " << effective_parallelism << "\n";
    std::cout << "  Speedup factor: " << effective_parallelism << "x\n";
    std::cout << "  1s sequential → " << parallel_time_ms << " ms parallel\n";
    
    EXPECT_GT(effective_parallelism, 1000.0) 
        << "GPU should provide >1000x parallelism";
}

// =============================================================================
// Consistency Checks
// =============================================================================

// Test: Evaluate same position multiple times gives same result
TEST_F(FPSThresholdTests, EvaluationConsistency) {
    Sirius::KerrSchildFamily schw{Sirius::KerrSchildParams::Schwarzschild(M)};
    schw.setParameter("mass", M);
    
    Vec4 pos = getStandardPosition();
    Metric4D g1, g2;
    Tensor<Dual<double>, 4, 4, 4> dg;
    
    schw.evaluate(pos, g1, dg);
    
    // Evaluate 100 more times
    for (int i = 0; i < 100; ++i) {
        schw.evaluate(pos, g2, dg);
        
        // Results should be identical
        for (int mu = 0; mu < 4; ++mu) {
            for (int nu = 0; nu < 4; ++nu) {
                EXPECT_DOUBLE_EQ(g1(mu, nu).real, g2(mu, nu).real)
                    << "Inconsistent evaluation at iteration " << i;
            }
        }
    }
}

// Test: No NaN or Inf in performance-critical paths
TEST_F(FPSThresholdTests, NoNaNInCriticalPath) {
    std::vector<IMetric*> metrics;
    Sirius::KerrSchildFamily mink{Sirius::KerrSchildParams::Minkowski()};
    Sirius::KerrSchildFamily schw{Sirius::KerrSchildParams::Schwarzschild(M)};
    Sirius::KerrSchildFamily kerr{Sirius::KerrSchildParams::Kerr(M, 0.9)};
    
    schw.setParameter("mass", M);
    kerr.setParameter("mass", M);
    kerr.setParameter("spin", 0.9);
    
    metrics = {&mink, &schw, &kerr};
    
    Vec4 pos = getStandardPosition();
    
    for (auto* metric : metrics) {
        Metric4D g;
        Tensor<Dual<double>, 4, 4, 4> dg;
        
        for (int i = 0; i < 100; ++i) {
            metric->evaluate(pos, g, dg);
            
            // Check for NaN/Inf
            for (int mu = 0; mu < 4; ++mu) {
                for (int nu = 0; nu < 4; ++nu) {
                    EXPECT_FALSE(std::isnan(g(mu, nu).real))
                        << "NaN detected for " << metric->getName();
                    EXPECT_FALSE(std::isinf(g(mu, nu).real))
                        << "Inf detected for " << metric->getName();
                }
            }
        }
    }
}

} // namespace sirius::test
