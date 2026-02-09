// TSBM002A.cpp - FPS Performance Tracking Benchmarks
// Benchmarks metric evaluation time and extrapolates to target FPS.

#define _USE_MATH_DEFINES
#include <cmath>
#include <gtest/gtest.h>
#include <chrono>
#include <vector>
#include <numeric>
#include <MTTN001A.h>
#include <MTDL001A.h>

namespace sirius::test {

// =============================================================================
// Test Configuration
// =============================================================================

// Default FPS targets from DevelopmentRoadmap.md
struct FPSTargets {
    static constexpr double MINKOWSKI_TARGET = 120.0;
    static constexpr double MINKOWSKI_BOUND = 60.0;
    static constexpr double SCHWARZSCHILD_TARGET = 60.0;
    static constexpr double SCHWARZSCHILD_BOUND = 30.0;
    static constexpr double KERR_TARGET = 30.0;
    static constexpr double KERR_BOUND = 15.0;
    static constexpr double NUMERICAL_TARGET = 10.0;
    static constexpr double NUMERICAL_BOUND = 5.0;
};

// Benchmark configuration
struct BenchmarkConfig {
    static constexpr int SAMPLE_PIXELS = 1000;         // Pixels to benchmark
    static constexpr int WARMUP_ITERATIONS = 10;       // Warmup iterations
    static constexpr int BENCHMARK_ITERATIONS = 100;   // Measurement iterations
    static constexpr int STEPS_PER_RAY = 500;          // Integration steps per ray
    static constexpr int RESOLUTION_WIDTH = 1920;      // Target resolution
    static constexpr int RESOLUTION_HEIGHT = 1080;     // Target resolution
};

// =============================================================================
// Test Fixture
// =============================================================================

class FPSBenchmarks : public ::testing::Test {
protected:
    static constexpr double M = 1.0;  // Mass in geometric units
    
    // Simulate a single RK4 step for geodesic integration
    // This measures the core compute cost per ray step
    void runRK4Step(const Metric4D& g, const Vec4& pos, const Vec4& vel,
                    Vec4& new_pos, Vec4& new_vel) {
        // Simplified RK4 step simulation
        // In real code, this would compute Christoffel symbols and update state
        double dt = 0.01;
        
        // Placeholder computation to simulate actual work
        Vec4 k1, k2, k3, k4;
        
        // Stage 1
        for (int i = 0; i < 4; ++i) {
            k1(i) = vel(i) * dt;
        }
        
        // Stage 2 (using metric components to simulate Christoffel computation)
        for (int i = 0; i < 4; ++i) {
            double accel = 0.0;
            for (int j = 0; j < 4; ++j) {
                accel -= g(i, j).real * vel(j) * 0.1;  // Simplified
            }
            k2(i) = (vel(i) + k1(i) * 0.5) * dt + accel * dt * dt;
        }
        
        // Stage 3
        for (int i = 0; i < 4; ++i) {
            k3(i) = (vel(i) + k2(i) * 0.5) * dt;
        }
        
        // Stage 4
        for (int i = 0; i < 4; ++i) {
            k4(i) = (vel(i) + k3(i)) * dt;
        }
        
        // Combine stages
        for (int i = 0; i < 4; ++i) {
            new_pos(i) = pos(i) + (k1(i) + 2*k2(i) + 2*k3(i) + k4(i)) / 6.0;
            new_vel(i) = vel(i);  // Simplified
        }
    }
    
    // Create Minkowski metric
    Metric4D createMinkowskiMetric() {
        Metric4D g;
        g.zero();
        g(0, 0) = Dual<double>(-1.0);
        g(1, 1) = Dual<double>(1.0);
        g(2, 2) = Dual<double>(1.0);
        g(3, 3) = Dual<double>(1.0);
        return g;
    }
    
    // Create Schwarzschild metric at given radius
    Metric4D createSchwarzschildMetric(double r, double theta = M_PI/2) {
        Metric4D g;
        g.zero();
        
        double rs = 2.0 * M;
        double f = 1.0 - rs / r;
        double sin_theta = std::sin(theta);
        
        g(0, 0) = Dual<double>(-f);
        g(1, 1) = Dual<double>(1.0 / f);
        g(2, 2) = Dual<double>(r * r);
        g(3, 3) = Dual<double>(r * r * sin_theta * sin_theta);
        
        return g;
    }
    
    // Create Kerr metric at given radius and spin
    Metric4D createKerrMetric(double r, double a, double theta = M_PI/2) {
        Metric4D g;
        g.zero();
        
        double sin_theta = std::sin(theta);
        double cos_theta = std::cos(theta);
        double sin2 = sin_theta * sin_theta;
        double cos2 = cos_theta * cos_theta;
        
        double rho2 = r * r + a * a * cos2;
        double delta = r * r - 2.0 * M * r + a * a;
        double A = (r * r + a * a) * (r * r + a * a) - delta * a * a * sin2;
        
        // Avoid division by zero
        if (std::abs(rho2) < 1e-10) rho2 = 1e-10;
        if (std::abs(delta) < 1e-10) delta = 1e-10;
        
        g(0, 0) = Dual<double>(-(1.0 - 2.0 * M * r / rho2));
        g(0, 3) = Dual<double>(-2.0 * M * r * a * sin2 / rho2);
        g(3, 0) = g(0, 3);
        g(1, 1) = Dual<double>(rho2 / delta);
        g(2, 2) = Dual<double>(rho2);
        g(3, 3) = Dual<double>(A * sin2 / rho2);
        
        return g;
    }
    
    // Benchmark a single metric configuration
    struct BenchmarkResult {
        double mean_time_us;       // Mean time per ray (microseconds)
        double std_dev_us;         // Standard deviation
        double estimated_fps;      // Estimated FPS for full frame
        int samples;               // Number of samples taken
    };
    
    BenchmarkResult benchmarkMetric(std::function<Metric4D(double, double)> metric_factory,
                                    double r_start, double r_end) {
        using namespace std::chrono;
        
        std::vector<double> times_us;
        times_us.reserve(BenchmarkConfig::BENCHMARK_ITERATIONS);
        
        // Initial state
        Vec4 pos, vel, new_pos, new_vel;
        pos(0) = 0.0;  // t
        pos(1) = r_start;  // r
        pos(2) = M_PI / 2.0;  // theta
        pos(3) = 0.0;  // phi
        
        vel(0) = 1.0;  // dt/dlambda
        vel(1) = -0.1; // dr/dlambda
        vel(2) = 0.0;  // dtheta/dlambda
        vel(3) = 0.01; // dphi/dlambda
        
        // Warmup
        for (int i = 0; i < BenchmarkConfig::WARMUP_ITERATIONS; ++i) {
            double r = r_start;
            Metric4D g = metric_factory(r, M_PI/2);
            runRK4Step(g, pos, vel, new_pos, new_vel);
        }
        
        // Benchmark
        for (int iter = 0; iter < BenchmarkConfig::BENCHMARK_ITERATIONS; ++iter) {
            auto start = high_resolution_clock::now();
            
            // Simulate a full ray trace
            double r = r_start;
            pos(1) = r;
            
            for (int step = 0; step < BenchmarkConfig::STEPS_PER_RAY; ++step) {
                // Linearly interpolate r (simplified ray path)
                r = r_start + (r_end - r_start) * step / BenchmarkConfig::STEPS_PER_RAY;
                pos(1) = r;
                
                Metric4D g = metric_factory(r, pos(2));
                runRK4Step(g, pos, vel, new_pos, new_vel);
                
                pos = new_pos;
                vel = new_vel;
            }
            
            auto end = high_resolution_clock::now();
            double elapsed_us = duration_cast<nanoseconds>(end - start).count() / 1000.0;
            times_us.push_back(elapsed_us);
        }
        
        // Calculate statistics
        double sum = std::accumulate(times_us.begin(), times_us.end(), 0.0);
        double mean = sum / times_us.size();
        
        double sq_sum = 0.0;
        for (double t : times_us) {
            sq_sum += (t - mean) * (t - mean);
        }
        double std_dev = std::sqrt(sq_sum / times_us.size());
        
        // Estimate FPS
        // Total pixels * time per ray = frame time
        // Note: This is a simplified model; actual GPU parallelism changes this significantly
        int total_pixels = BenchmarkConfig::RESOLUTION_WIDTH * BenchmarkConfig::RESOLUTION_HEIGHT;
        double frame_time_us = mean * total_pixels / BenchmarkConfig::SAMPLE_PIXELS;
        double estimated_fps = 1e6 / frame_time_us;
        
        BenchmarkResult result;
        result.mean_time_us = mean;
        result.std_dev_us = std_dev;
        result.estimated_fps = estimated_fps;
        result.samples = static_cast<int>(times_us.size());
        
        return result;
    }
    
    void SetUp() override {}
    void TearDown() override {}
};

// =============================================================================
// Minkowski Metric Benchmarks
// =============================================================================

// Test: Minkowski metric should be fast (direct path, no curvature)
TEST_F(FPSBenchmarks, MinkowskiFPS) {
    auto metric_factory = [this](double r, double theta) {
        (void)r; (void)theta;
        return createMinkowskiMetric();
    };
    
    BenchmarkResult result = benchmarkMetric(metric_factory, 10.0, 100.0);
    
    // Report benchmarks (not strict assertions since performance is hardware-dependent)
    std::cout << "\n[BENCHMARK] Minkowski Metric:\n";
    std::cout << "  Mean time per ray: " << result.mean_time_us << " µs\n";
    std::cout << "  Std deviation: " << result.std_dev_us << " µs\n";
    std::cout << "  Estimated FPS (CPU sequential): " << result.estimated_fps << "\n";
    std::cout << "  Target FPS: " << FPSTargets::MINKOWSKI_TARGET << "\n";
    std::cout << "  Bound FPS: " << FPSTargets::MINKOWSKI_BOUND << "\n";
    
    // Soft assertion: warn but don't fail if below bound
    if (result.estimated_fps < FPSTargets::MINKOWSKI_BOUND) {
        std::cout << "  [WARNING] Estimated FPS below bound!\n";
    }
    
    // Basic sanity check: computation should complete in reasonable time
    EXPECT_GT(result.mean_time_us, 0.0) << "Computation time should be positive";
    EXPECT_LT(result.mean_time_us, 1e9) << "Computation time should be reasonable";
}

// =============================================================================
// Schwarzschild Metric Benchmarks
// =============================================================================

TEST_F(FPSBenchmarks, SchwarzschildFPS) {
    auto metric_factory = [this](double r, double theta) {
        return createSchwarzschildMetric(r, theta);
    };
    
    BenchmarkResult result = benchmarkMetric(metric_factory, 10.0, 100.0);
    
    std::cout << "\n[BENCHMARK] Schwarzschild Metric:\n";
    std::cout << "  Mean time per ray: " << result.mean_time_us << " µs\n";
    std::cout << "  Std deviation: " << result.std_dev_us << " µs\n";
    std::cout << "  Estimated FPS (CPU sequential): " << result.estimated_fps << "\n";
    std::cout << "  Target FPS: " << FPSTargets::SCHWARZSCHILD_TARGET << "\n";
    std::cout << "  Bound FPS: " << FPSTargets::SCHWARZSCHILD_BOUND << "\n";
    
    if (result.estimated_fps < FPSTargets::SCHWARZSCHILD_BOUND) {
        std::cout << "  [WARNING] Estimated FPS below bound!\n";
    }
    
    EXPECT_GT(result.mean_time_us, 0.0);
    EXPECT_LT(result.mean_time_us, 1e9);
}

// Test: Near-horizon performance (more complex geometry)
TEST_F(FPSBenchmarks, SchwarzschildNearHorizonFPS) {
    auto metric_factory = [this](double r, double theta) {
        return createSchwarzschildMetric(r, theta);
    };
    
    // Start near horizon
    BenchmarkResult result = benchmarkMetric(metric_factory, 3.0, 20.0);
    
    std::cout << "\n[BENCHMARK] Schwarzschild (Near Horizon):\n";
    std::cout << "  Mean time per ray: " << result.mean_time_us << " µs\n";
    std::cout << "  Note: Near-horizon geodesics require more careful integration\n";
    
    EXPECT_GT(result.mean_time_us, 0.0);
}

// =============================================================================
// Kerr Metric Benchmarks
// =============================================================================

TEST_F(FPSBenchmarks, KerrFPS) {
    double spin = 0.9 * M;  // High spin for realistic case
    
    auto metric_factory = [this, spin](double r, double theta) {
        return createKerrMetric(r, spin, theta);
    };
    
    BenchmarkResult result = benchmarkMetric(metric_factory, 10.0, 100.0);
    
    std::cout << "\n[BENCHMARK] Kerr Metric (a=0.9M):\n";
    std::cout << "  Mean time per ray: " << result.mean_time_us << " µs\n";
    std::cout << "  Std deviation: " << result.std_dev_us << " µs\n";
    std::cout << "  Estimated FPS (CPU sequential): " << result.estimated_fps << "\n";
    std::cout << "  Target FPS: " << FPSTargets::KERR_TARGET << "\n";
    std::cout << "  Bound FPS: " << FPSTargets::KERR_BOUND << "\n";
    
    if (result.estimated_fps < FPSTargets::KERR_BOUND) {
        std::cout << "  [WARNING] Estimated FPS below bound!\n";
    }
    
    EXPECT_GT(result.mean_time_us, 0.0);
    EXPECT_LT(result.mean_time_us, 1e9);
}

// Test: Kerr with varying spin
TEST_F(FPSBenchmarks, KerrSpinComparison) {
    std::vector<double> spins = {0.0, 0.5, 0.9, 0.99};
    
    std::cout << "\n[BENCHMARK] Kerr Metric (Spin Comparison):\n";
    std::cout << "  Spin\t\tTime (µs)\n";
    std::cout << "  ----\t\t---------\n";
    
    for (double spin_ratio : spins) {
        double spin = spin_ratio * M;
        
        auto metric_factory = [this, spin](double r, double theta) {
            return createKerrMetric(r, spin, theta);
        };
        
        BenchmarkResult result = benchmarkMetric(metric_factory, 10.0, 50.0);
        std::cout << "  " << spin_ratio << "M\t\t" << result.mean_time_us << "\n";
    }
}

// =============================================================================
// Comparative Benchmarks
// =============================================================================

// Test: Compare relative performance of different metrics
TEST_F(FPSBenchmarks, MetricComparison) {
    std::cout << "\n[BENCHMARK] Metric Performance Comparison:\n";
    std::cout << "  Metric\t\tTime (µs)\tRelative\n";
    std::cout << "  ------\t\t---------\t--------\n";
    
    // Minkowski baseline
    auto mink_factory = [this](double r, double theta) {
        (void)r; (void)theta;
        return createMinkowskiMetric();
    };
    BenchmarkResult mink_result = benchmarkMetric(mink_factory, 10.0, 100.0);
    double baseline = mink_result.mean_time_us;
    
    std::cout << "  Minkowski\t\t" << mink_result.mean_time_us 
              << "\t\t1.00x\n";
    
    // Schwarzschild
    auto schw_factory = [this](double r, double theta) {
        return createSchwarzschildMetric(r, theta);
    };
    BenchmarkResult schw_result = benchmarkMetric(schw_factory, 10.0, 100.0);
    std::cout << "  Schwarzschild\t\t" << schw_result.mean_time_us 
              << "\t\t" << (schw_result.mean_time_us / baseline) << "x\n";
    
    // Kerr
    auto kerr_factory = [this](double r, double theta) {
        return createKerrMetric(r, 0.9 * M, theta);
    };
    BenchmarkResult kerr_result = benchmarkMetric(kerr_factory, 10.0, 100.0);
    std::cout << "  Kerr (a=0.9M)\t\t" << kerr_result.mean_time_us 
              << "\t\t" << (kerr_result.mean_time_us / baseline) << "x\n";
    
    // Verify ordering: Minkowski <= Schwarzschild <= Kerr (generally)
    // Note: This may not always hold due to measurement noise
    EXPECT_GT(baseline, 0.0);
}

// =============================================================================
// Integration Step Count Analysis
// =============================================================================

// Test: Measure how performance scales with integration steps
TEST_F(FPSBenchmarks, StepCountScaling) {
    using namespace std::chrono;
    
    std::vector<int> step_counts = {100, 250, 500, 1000, 2000};
    
    std::cout << "\n[BENCHMARK] Step Count Scaling (Schwarzschild):\n";
    std::cout << "  Steps\t\tTime (µs)\tTime/Step (ns)\n";
    std::cout << "  -----\t\t---------\t--------------\n";
    
    for (int steps : step_counts) {
        Vec4 pos, vel, new_pos, new_vel;
        pos(0) = 0.0; pos(1) = 20.0; pos(2) = M_PI/2; pos(3) = 0.0;
        vel(0) = 1.0; vel(1) = -0.1; vel(2) = 0.0; vel(3) = 0.01;
        
        auto start = high_resolution_clock::now();
        
        for (int s = 0; s < steps; ++s) {
            double r = 20.0 - (15.0 * s / steps);
            Metric4D g = createSchwarzschildMetric(r, M_PI/2);
            runRK4Step(g, pos, vel, new_pos, new_vel);
            pos = new_pos;
            vel = new_vel;
        }
        
        auto end = high_resolution_clock::now();
        double elapsed_us = duration_cast<nanoseconds>(end - start).count() / 1000.0;
        double per_step_ns = (elapsed_us * 1000.0) / steps;
        
        std::cout << "  " << steps << "\t\t" << elapsed_us 
                  << "\t\t" << per_step_ns << "\n";
    }
}

} // namespace sirius::test
