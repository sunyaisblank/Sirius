// TSMT002A.cpp - First-Principles Validation for Tensor Operations
// Tests: vector ops, Minkowski metric, inner products, Christoffel symbols.

#include <gtest/gtest.h>
#include <cmath>
#include <MTTN001A.h>
#include <MTDL001A.h>

namespace sirius::test {

constexpr double kEpsilon = 1e-10;

// =============================================================================
// Test Fixture
// =============================================================================

class TensorTests : public ::testing::Test {
protected:
    void SetUp() override {}
    void TearDown() override {}
};

// =============================================================================
// Vector Operations (Tensor<T, 4>)
// =============================================================================

// Test: Zero initialization
TEST_F(TensorTests, VectorZeroInitialization) {
    Vec4 v;
    v.zero();
    
    for (int i = 0; i < 4; ++i) {
        EXPECT_NEAR(v(i), 0.0, kEpsilon);
    }
}

// Test: Vector addition
TEST_F(TensorTests, VectorAddition) {
    Vec4 u, v;
    u(0) = 1.0; u(1) = 2.0; u(2) = 3.0; u(3) = 4.0;
    v(0) = 5.0; v(1) = 6.0; v(2) = 7.0; v(3) = 8.0;
    
    Vec4 w = u + v;
    
    EXPECT_NEAR(w(0), 6.0, kEpsilon);
    EXPECT_NEAR(w(1), 8.0, kEpsilon);
    EXPECT_NEAR(w(2), 10.0, kEpsilon);
    EXPECT_NEAR(w(3), 12.0, kEpsilon);
}

// Test: Vector subtraction
TEST_F(TensorTests, VectorSubtraction) {
    Vec4 u, v;
    u(0) = 5.0; u(1) = 6.0; u(2) = 7.0; u(3) = 8.0;
    v(0) = 1.0; v(1) = 2.0; v(2) = 3.0; v(3) = 4.0;
    
    Vec4 w = u - v;
    
    EXPECT_NEAR(w(0), 4.0, kEpsilon);
    EXPECT_NEAR(w(1), 4.0, kEpsilon);
    EXPECT_NEAR(w(2), 4.0, kEpsilon);
    EXPECT_NEAR(w(3), 4.0, kEpsilon);
}

// Test: Scalar multiplication
TEST_F(TensorTests, VectorScalarMultiplication) {
    Vec4 u;
    u(0) = 1.0; u(1) = 2.0; u(2) = 3.0; u(3) = 4.0;
    
    Vec4 v = u * 2.0;
    Vec4 w = 3.0 * u;
    
    EXPECT_NEAR(v(0), 2.0, kEpsilon);
    EXPECT_NEAR(v(3), 8.0, kEpsilon);
    EXPECT_NEAR(w(0), 3.0, kEpsilon);
    EXPECT_NEAR(w(3), 12.0, kEpsilon);
}

// Test: Vector length (Euclidean)
TEST_F(TensorTests, VectorLength) {
    Vec4 u;
    u(0) = 1.0; u(1) = 2.0; u(2) = 2.0; u(3) = 0.0;
    
    double len = u.length();
    
    EXPECT_NEAR(len, 3.0, kEpsilon);  // sqrt(1 + 4 + 4 + 0) = 3
}

// Test: Unary minus
TEST_F(TensorTests, VectorUnaryMinus) {
    Vec4 u;
    u(0) = 1.0; u(1) = -2.0; u(2) = 3.0; u(3) = -4.0;
    
    Vec4 v = -u;
    
    EXPECT_NEAR(v(0), -1.0, kEpsilon);
    EXPECT_NEAR(v(1), 2.0, kEpsilon);
    EXPECT_NEAR(v(2), -3.0, kEpsilon);
    EXPECT_NEAR(v(3), 4.0, kEpsilon);
}

// =============================================================================
// Matrix Operations (Tensor<T, 4, 4>)
// =============================================================================

// Test: Matrix zero initialization
TEST_F(TensorTests, MatrixZeroInitialization) {
    Metric4D g;
    g.zero();
    
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            EXPECT_NEAR(g(i, j).real, 0.0, kEpsilon);
            EXPECT_NEAR(g(i, j).dual, 0.0, kEpsilon);
        }
    }
}

// =============================================================================
// Minkowski Metric Properties
// =============================================================================

// Helper: Create Minkowski metric η_μν = diag(-1, 1, 1, 1)
Metric4D createMinkowskiMetric() {
    Metric4D eta;
    eta.zero();
    eta(0, 0) = Dual<double>(-1.0, 0.0);
    eta(1, 1) = Dual<double>(1.0, 0.0);
    eta(2, 2) = Dual<double>(1.0, 0.0);
    eta(3, 3) = Dual<double>(1.0, 0.0);
    return eta;
}

// Test: Minkowski metric is symmetric
TEST_F(TensorTests, MinkowskiSymmetry) {
    Metric4D eta = createMinkowskiMetric();
    
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            EXPECT_NEAR(eta(i, j).real, eta(j, i).real, kEpsilon)
                << "Metric not symmetric at (" << i << "," << j << ")";
        }
    }
}

// Test: Minkowski metric is diagonal
TEST_F(TensorTests, MinkowskiDiagonal) {
    Metric4D eta = createMinkowskiMetric();
    
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            if (i != j) {
                EXPECT_NEAR(eta(i, j).real, 0.0, kEpsilon)
                    << "Off-diagonal element non-zero at (" << i << "," << j << ")";
            }
        }
    }
}

// Test: Minkowski metric signature (-,+,+,+)
TEST_F(TensorTests, MinkowskiSignature) {
    Metric4D eta = createMinkowskiMetric();
    
    EXPECT_LT(eta(0, 0).real, 0.0) << "Temporal component should be negative";
    EXPECT_GT(eta(1, 1).real, 0.0) << "Spatial component should be positive";
    EXPECT_GT(eta(2, 2).real, 0.0) << "Spatial component should be positive";
    EXPECT_GT(eta(3, 3).real, 0.0) << "Spatial component should be positive";
}

// =============================================================================
// TensorOps Tests
// =============================================================================

// Test: Inner product with Minkowski metric
TEST_F(TensorTests, InnerProductMinkowski) {
    Metric4D eta = createMinkowskiMetric();
    
    Vec4 u, v;
    u(0) = 1.0; u(1) = 0.0; u(2) = 0.0; u(3) = 0.0;  // Temporal unit vector
    v(0) = 1.0; v(1) = 0.0; v(2) = 0.0; v(3) = 0.0;
    
    double inner = TensorOps::innerProduct(u, v, eta);
    
    // <(1,0,0,0), (1,0,0,0)>_η = -1*1 + 0 + 0 + 0 = -1
    EXPECT_NEAR(inner, -1.0, kEpsilon);
}

// Test: Spatial vectors have positive inner product
TEST_F(TensorTests, InnerProductSpatial) {
    Metric4D eta = createMinkowskiMetric();
    
    Vec4 u;
    u(0) = 0.0; u(1) = 1.0; u(2) = 0.0; u(3) = 0.0;  // Spatial unit vector
    
    double inner = TensorOps::innerProduct(u, u, eta);
    
    EXPECT_NEAR(inner, 1.0, kEpsilon);
}

// Test: Null vector has zero inner product
TEST_F(TensorTests, NullVectorInnerProduct) {
    Metric4D eta = createMinkowskiMetric();
    
    // Light-like vector: (1, 1, 0, 0)
    Vec4 k;
    k(0) = 1.0; k(1) = 1.0; k(2) = 0.0; k(3) = 0.0;
    
    double inner = TensorOps::innerProduct(k, k, eta);
    
    // <(1,1,0,0), (1,1,0,0)>_η = -1 + 1 + 0 + 0 = 0
    EXPECT_NEAR(inner, 0.0, kEpsilon);
}

// =============================================================================
// 3-Tensor Operations (Tensor<T, 4, 4, 4>)
// =============================================================================

// Test: 3-tensor zero initialization
TEST_F(TensorTests, Tensor3ZeroInitialization) {
    Tensor<double, 4, 4, 4> T;
    T.zero();
    
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            for (int k = 0; k < 4; ++k) {
                EXPECT_NEAR(T(i, j, k), 0.0, kEpsilon);
            }
        }
    }
}

// =============================================================================
// Christoffel Symbol Properties
// =============================================================================

// Test: Christoffel symbols of flat space are zero
TEST_F(TensorTests, FlatSpaceChristoffelZero) {
    Metric4D eta = createMinkowskiMetric();
    
    // Create zero derivative tensor (flat space has constant metric)
    Tensor<Dual<double>, 4, 4, 4> dg;
    dg.zero();
    
    ChristoffelSymbols gamma = TensorOps::christoffel(eta, dg);
    
    // All Christoffel symbols should be zero for flat space
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            for (int k = 0; k < 4; ++k) {
                EXPECT_NEAR(gamma.gamma(i, j, k).real, 0.0, kEpsilon)
                    << "Non-zero Christoffel at Γ^" << i << "_" << j << k;
            }
        }
    }
}

// Test: Christoffel symbols are symmetric in lower indices
// Γ^λ_μν = Γ^λ_νμ (torsion-free connection)
TEST_F(TensorTests, ChristoffelSymmetry) {
    // Use a non-trivial metric (simplified Schwarzschild-like)
    Metric4D g;
    g.zero();
    double r = 10.0;  // Some radius
    double rs = 1.0;  // Schwarzschild radius
    double f = 1.0 - rs/r;
    
    g(0, 0) = Dual<double>(-f, 0.0);
    g(1, 1) = Dual<double>(1.0/f, 0.0);
    g(2, 2) = Dual<double>(r*r, 0.0);
    g(3, 3) = Dual<double>(r*r, 0.0);  // Simplified (ignoring sin²θ)
    
    // Create metric derivatives (simplified)
    Tensor<Dual<double>, 4, 4, 4> dg;
    dg.zero();
    // ∂g_00/∂r = -∂f/∂r = -rs/r²
    dg(1, 0, 0) = Dual<double>(rs/(r*r), 0.0);
    // ∂g_11/∂r = ∂(1/f)/∂r = rs/(r²f²)
    dg(1, 1, 1) = Dual<double>(-rs/(r*r*f*f), 0.0);
    
    ChristoffelSymbols gamma = TensorOps::christoffel(g, dg);
    
    // Check symmetry: Γ^λ_μν = Γ^λ_νμ
    for (int lam = 0; lam < 4; ++lam) {
        for (int mu = 0; mu < 4; ++mu) {
            for (int nu = 0; nu < 4; ++nu) {
                EXPECT_NEAR(gamma.gamma(lam, mu, nu).real, 
                           gamma.gamma(lam, nu, mu).real, kEpsilon)
                    << "Torsion-free symmetry violated at Γ^" << lam << "_" << mu << nu;
            }
        }
    }
}

} // namespace sirius::test
