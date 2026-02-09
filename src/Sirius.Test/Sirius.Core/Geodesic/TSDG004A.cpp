// TSDG004A.cpp - Determinism Verification Tests
// Tests: bit-identical results across iterations for all computations.

#define _USE_MATH_DEFINES
#include <cmath>
#include <gtest/gtest.h>
#include <vector>
#include <cstring>
#include "MTDL001A.h"
#include "MTTN001A.h"
#include "PHMT000A.h"
#include "PHMT100A.h"  // Unified Kerr-Schild Family
#include "PHGD001A.h"

namespace sirius::test {

// =============================================================================
// Constants
// =============================================================================

constexpr int DETERMINISM_ITERATIONS = 10;

// =============================================================================
// Test Fixture
// =============================================================================

class DeterminismTests : public ::testing::Test {
protected:
    Sirius::KerrSchildFamily schwarzschild;
    Sirius::KerrSchildFamily kerr;
    
    void SetUp() override {
        schwarzschild.setParameter("mass", 1.0);
        kerr.setParameter("mass", 1.0);
        kerr.setParameter("spin", 0.5);
    }
    
    // Bit-level comparison of doubles
    bool bitIdentical(double a, double b) {
        return std::memcmp(&a, &b, sizeof(double)) == 0;
    }
    
    // Hash a Vec4 for comparison
    uint64_t hashVec4(const Vec4& v) {
        uint64_t hash = 0;
        for (int i = 0; i < 4; ++i) {
            uint64_t bits;
            std::memcpy(&bits, &v(i), sizeof(double));
            hash ^= bits + 0x9e3779b97f4a7c15 + (hash << 6) + (hash >> 2);
        }
        return hash;
    }
    
    // Hash a Metric4D for comparison
    uint64_t hashMetric(const Metric4D& g) {
        uint64_t hash = 0;
        for (int mu = 0; mu < 4; ++mu) {
            for (int nu = 0; nu < 4; ++nu) {
                uint64_t bits;
                double val = g(mu, nu).real;
                std::memcpy(&bits, &val, sizeof(double));
                hash ^= bits + 0x9e3779b97f4a7c15 + (hash << 6) + (hash >> 2);
            }
        }
        return hash;
    }
};

// =============================================================================
// Dual Number Determinism
// =============================================================================

TEST_F(DeterminismTests, DualArithmeticDeterminism)
{
    // FORMAL SPECIFICATION D1:
    // Dual number arithmetic must be bit-identical across runs
    
    Dual<double> a(3.14159265358979, 1.0);
    Dual<double> b(2.71828182845904, 0.5);
    
    std::vector<Dual<double>> results;
    
    for (int iter = 0; iter < DETERMINISM_ITERATIONS; ++iter) {
        Dual<double> sum = a + b;
        Dual<double> prod = a * b;
        Dual<double> quot = a / b;
        Dual<double> sq = sqrt(a);
        Dual<double> s = sin(a);
        Dual<double> c = cos(a);
        
        if (iter == 0) {
            results = {sum, prod, quot, sq, s, c};
        } else {
            // All iterations must match first iteration bit-identically
            EXPECT_TRUE(bitIdentical(sum.real, results[0].real)) 
                << "Addition not deterministic at iteration " << iter;
            EXPECT_TRUE(bitIdentical(prod.real, results[1].real))
                << "Multiplication not deterministic at iteration " << iter;
            EXPECT_TRUE(bitIdentical(quot.real, results[2].real))
                << "Division not deterministic at iteration " << iter;
            EXPECT_TRUE(bitIdentical(sq.real, results[3].real))
                << "Square root not deterministic at iteration " << iter;
            EXPECT_TRUE(bitIdentical(s.real, results[4].real))
                << "Sine not deterministic at iteration " << iter;
            EXPECT_TRUE(bitIdentical(c.real, results[5].real))
                << "Cosine not deterministic at iteration " << iter;
        }
    }
}

// =============================================================================
// Metric Evaluation Determinism
// =============================================================================

TEST_F(DeterminismTests, SchwarzschildMetricDeterminism)
{
    // FORMAL SPECIFICATION D1:
    // Metric evaluation must produce identical results across calls
    
    Vec4 pos;
    // Spherical (10.0, PI/2, 0.0) -> Cartesian (10.0, 0.0, 0.0)
    pos(0) = 0.0;
    pos(1) = 10.0; // x
    pos(2) = 0.0;  // y
    pos(3) = 0.0;  // z
    
    uint64_t reference_hash = 0;
    
    for (int iter = 0; iter < DETERMINISM_ITERATIONS; ++iter) {
        Metric4D g;
        Tensor<Dual<double>, 4, 4, 4> dg;
        schwarzschild.evaluate(pos, g, dg);
        
        uint64_t current_hash = hashMetric(g);
        
        if (iter == 0) {
            reference_hash = current_hash;
        } else {
            EXPECT_EQ(current_hash, reference_hash)
                << "Schwarzschild metric not deterministic at iteration " << iter;
        }
    }
}

TEST_F(DeterminismTests, KerrMetricDeterminism)
{
    // FORMAL SPECIFICATION D1:
    // Kerr metric evaluation must be deterministic
    
    Vec4 pos;
    // Spherical (5.0, PI/3, 1.0) -> Cartesian
    double r = 5.0;
    double th = M_PI / 3.0;
    double ph = 1.0;
    pos(0) = 0.0;
    pos(1) = r * std::sin(th) * std::cos(ph); // x
    pos(2) = r * std::sin(th) * std::sin(ph); // y
    pos(3) = r * std::cos(th);                // z
    
    uint64_t reference_hash = 0;
    
    for (int iter = 0; iter < DETERMINISM_ITERATIONS; ++iter) {
        Metric4D g;
        Tensor<Dual<double>, 4, 4, 4> dg;
        kerr.evaluate(pos, g, dg);
        
        uint64_t current_hash = hashMetric(g);
        
        if (iter == 0) {
            reference_hash = current_hash;
        } else {
            EXPECT_EQ(current_hash, reference_hash)
                << "Kerr metric not deterministic at iteration " << iter;
        }
    }
}

// =============================================================================
// Christoffel Symbol Determinism
// =============================================================================

TEST_F(DeterminismTests, ChristoffelDeterminism)
{
    // FORMAL SPECIFICATION D1:
    // Christoffel symbol computation must be bit-identical
    
    Vec4 pos;
    // Spherical (10.0, PI/2, 0.0) -> Cartesian (10.0, 0.0, 0.0)
    pos(0) = 0.0;
    pos(1) = 10.0;
    pos(2) = 0.0;
    pos(3) = 0.0;
    
    Metric4D g;
    Tensor<Dual<double>, 4, 4, 4> dg;
    schwarzschild.evaluate(pos, g, dg);
    
    uint64_t reference_hash = 0;
    
    for (int iter = 0; iter < DETERMINISM_ITERATIONS; ++iter) {
        ChristoffelSymbols christoffel = TensorOps::christoffel(g, dg);
        
        // Hash all Christoffel components
        uint64_t hash = 0;
        for (int l = 0; l < 4; ++l) {
            for (int m = 0; m < 4; ++m) {
                for (int n = 0; n < 4; ++n) {
                    uint64_t bits;
                    double val = christoffel.gamma(l, m, n).real;
                    std::memcpy(&bits, &val, sizeof(double));
                    hash ^= bits + 0x9e3779b97f4a7c15 + (hash << 6) + (hash >> 2);
                }
            }
        }
        
        if (iter == 0) {
            reference_hash = hash;
        } else {
            EXPECT_EQ(hash, reference_hash)
                << "Christoffel symbols not deterministic at iteration " << iter;
        }
    }
}

// =============================================================================
// Geodesic Integration Determinism
// =============================================================================

TEST_F(DeterminismTests, GeodesicIntegrationDeterminism)
{
    // FORMAL SPECIFICATION D1:
    // Geodesic integration with same initial conditions must be deterministic
    
    Vec4 initial_pos;
    // Spherical (15.0, PI/2, 0.0) -> Cartesian (15.0, 0.0, 0.0)
    initial_pos(0) = 0.0;
    initial_pos(1) = 15.0;
    initial_pos(2) = 0.0;
    initial_pos(3) = 0.0;
    
    // Cartesian velocity (arbitrary but valid)
    Vec4 initial_vel;
    initial_vel(0) = 1.0;
    initial_vel(1) = -0.5; // vx
    initial_vel(2) = 0.0;  // vy
    initial_vel(3) = 0.1;  // vz
    
    uint64_t reference_hash = 0;
    
    for (int iter = 0; iter < DETERMINISM_ITERATIONS; ++iter) {
        // Create identical rays
        Lightray ray;
        ray.position = initial_pos;
        ray.velocity = initial_vel;
        ray.proper_time = 0.0f;
        ray.coordinate_time = 0.0f;
        ray.terminated = 0;
        ray.step_size = 0.01f;
        ray.bounce_count = 0;
        
        // Normalize to null
        Metric4D g;
        Tensor<Dual<double>, 4, 4, 4> dg;
        schwarzschild.evaluate(initial_pos, g, dg);
        ray.velocity = TensorOps::normalizeNull(ray.velocity, g);
        ray.acceleration = Geodesic::calculateAcceleration(ray.velocity, ray.position, &schwarzschild);
        
        // Integrate 100 steps
        for (int step = 0; step < 100; ++step) {
            if (ray.terminated) break;
            Geodesic::integrateStep(ray, &schwarzschild);
        }
        
        // Hash final state
        uint64_t hash = hashVec4(ray.position) ^ hashVec4(ray.velocity);
        
        if (iter == 0) {
            reference_hash = hash;
        } else {
            EXPECT_EQ(hash, reference_hash)
                << "Geodesic integration not deterministic at iteration " << iter;
        }
    }
}

// =============================================================================
// Inner Product Determinism
// =============================================================================

TEST_F(DeterminismTests, InnerProductDeterminism)
{
    // FORMAL SPECIFICATION D1:
    // Inner product computation must be deterministic
    
    Vec4 u, v;
    u(0) = 1.0; u(1) = 0.5; u(2) = 0.25; u(3) = 0.1;
    v(0) = 0.8; v(1) = 0.6; v(2) = 0.4; v(3) = 0.2;
    
    Vec4 pos;
    // Spherical (5.0, PI/3, 0.0) -> Cartesian
    double r = 5.0;
    double th = M_PI / 3.0;
    double ph = 0.0;
    pos(0) = 0.0; 
    pos(1) = r * std::sin(th) * std::cos(ph); // x
    pos(2) = r * std::sin(th) * std::sin(ph); // y
    pos(3) = r * std::cos(th);                // z
    
    Metric4D g;
    Tensor<Dual<double>, 4, 4, 4> dg;
    schwarzschild.evaluate(pos, g, dg);
    
    double reference = 0.0;
    bool first = true;
    
    for (int iter = 0; iter < DETERMINISM_ITERATIONS; ++iter) {
        double result = TensorOps::innerProduct(u, v, g);
        
        if (first) {
            reference = result;
            first = false;
        } else {
            EXPECT_TRUE(bitIdentical(result, reference))
                << "Inner product not deterministic at iteration " << iter;
        }
    }
}

// =============================================================================
// Null Normalization Determinism
// =============================================================================

TEST_F(DeterminismTests, NullNormalizationDeterminism)
{
    // FORMAL SPECIFICATION D1:
    // Null vector normalization must be deterministic
    
    Vec4 pos;
    // Cartesian (on x-axis)
    pos(0) = 0.0; pos(1) = 10.0; pos(2) = 0.0; pos(3) = 0.0;
    
    Vec4 spatial;
    spatial(0) = 0.0; spatial(1) = -1.0; spatial(2) = 0.0; spatial(3) = 0.5;
    
    Metric4D g;
    Tensor<Dual<double>, 4, 4, 4> dg;
    schwarzschild.evaluate(pos, g, dg);
    
    uint64_t reference_hash = 0;
    
    for (int iter = 0; iter < DETERMINISM_ITERATIONS; ++iter) {
        Vec4 k = TensorOps::normalizeNull(spatial, g);
        uint64_t hash = hashVec4(k);
        
        if (iter == 0) {
            reference_hash = hash;
        } else {
            EXPECT_EQ(hash, reference_hash)
                << "Null normalization not deterministic at iteration " << iter;
        }
    }
}

} // namespace sirius::test
