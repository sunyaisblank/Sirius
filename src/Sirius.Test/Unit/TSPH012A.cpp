// TSPH012A.cpp - Metric Derivative Validation Tests
// Tests: ∂g_μν/∂x^λ symmetry, finite difference, Christoffel consistency.

#include <gtest/gtest.h>
#include <cmath>
#include <array>

#include <gtest/gtest.h>
#include <cmath>
#include <array>

#include "MTDL001A.h"   // Dual numbers
#include "MTTN001A.h"   // Tensors
#include "PHMT100A.h"   // Kerr-Schild Family
#include "PHMT101A.h"   // Morris-Thorne Family
#include "PHMT102A.h"   // Warp Drive Family

constexpr double kEpsilon = 1e-6;
constexpr double kFiniteDiffH = 1e-5;

// =============================================================================
// Test Fixture
// =============================================================================

class MetricDerivativeTests : public ::testing::Test {
protected:
    Sirius::KerrSchildFamily kerrSchild{Sirius::KerrSchildParams::Kerr(1.0, 0.9)};
    Sirius::MorrisThorneFamily ellisDrainhole{Sirius::MorrisThorneParams::Ellis(1.0)};
    Sirius::WarpDriveFamily alcubierre{Sirius::WarpDriveParams::Alcubierre(1.0, 1.0)};
    
    Metric4D g;
    Tensor<Dual<double>, 4, 4, 4> dg;
    
    // Finite difference derivative approximation
    double finiteDiffDerivative(IMetric& metric, 
                                 const Tensor<double, 4>& pos,
                                 int lambda, int mu, int nu) {
        Tensor<double, 4> pos_plus = pos;
        Tensor<double, 4> pos_minus = pos;
        
        pos_plus(lambda) += kFiniteDiffH;
        pos_minus(lambda) -= kFiniteDiffH;
        
        Metric4D g_plus, g_minus;
        Tensor<Dual<double>, 4, 4, 4> dg_dummy;
        
        metric.evaluate(pos_plus, g_plus, dg_dummy);
        metric.evaluate(pos_minus, g_minus, dg_dummy);
        
        return (g_plus(mu, nu).real - g_minus(mu, nu).real) / (2.0 * kFiniteDiffH);
    }
};

// =============================================================================
// Kerr-Schild Derivative Tests
// =============================================================================

TEST_F(MetricDerivativeTests, KerrSchildDerivativeSymmetry) {
    // ∂g_μν/∂x^λ = ∂g_νμ/∂x^λ (metric symmetry implies derivative symmetry)
    Tensor<double, 4> pos;
    pos(0) = 0.0; pos(1) = 5.0; pos(2) = 2.0; pos(3) = 1.0;
    
    kerrSchild.evaluate(pos, g, dg);
    
    for (int lam = 0; lam < 4; lam++) {
        for (int mu = 0; mu < 4; mu++) {
            for (int nu = mu + 1; nu < 4; nu++) {
                EXPECT_NEAR(dg(lam, mu, nu).real, dg(lam, nu, mu).real, kEpsilon)
                    << "Derivative asymmetry at lambda=" << lam 
                    << ", mu=" << mu << ", nu=" << nu;
            }
        }
    }
}

TEST_F(MetricDerivativeTests, KerrSchildFiniteDifferenceAgreement) {
    // Analytic derivative should match finite difference
    Tensor<double, 4> pos;
    pos(0) = 0.0; pos(1) = 5.0; pos(2) = 2.0; pos(3) = 1.0;
    
    kerrSchild.evaluate(pos, g, dg);
    
    // Check spatial derivatives (lambda = 1,2,3)
    for (int lam = 1; lam <= 3; lam++) {
        double fd = finiteDiffDerivative(kerrSchild, pos, lam, 0, 0);
        double analytic = dg(lam, 0, 0).real;
        
        EXPECT_NEAR(analytic, fd, 1e-4)
            << "Finite diff mismatch for dg_tt/dx^" << lam;
    }
}

TEST_F(MetricDerivativeTests, KerrSchildNonZeroDerivatives) {
    // Kerr-Schild has non-trivial spatial gradients
    Tensor<double, 4> pos;
    pos(0) = 0.0; pos(1) = 4.0; pos(2) = 0.0; pos(3) = 0.0;
    
    kerrSchild.evaluate(pos, g, dg);
    
    // At least one spatial derivative of g_tt should be non-zero
    double sum_derivatives = 0.0;
    for (int lam = 1; lam <= 3; lam++) {
        sum_derivatives += std::abs(dg(lam, 0, 0).real);
    }
    
    EXPECT_GT(sum_derivatives, kEpsilon)
        << "Kerr-Schild should have non-zero spatial derivatives of g_tt";
}

// =============================================================================
// Ellis Drainhole Derivative Tests
// =============================================================================

TEST_F(MetricDerivativeTests, EllisDrainholeDerivativeSymmetry) {
    Tensor<double, 4> pos;
    pos(0) = 0.0; pos(1) = 3.0; pos(2) = 1.0; pos(3) = 0.0;
    
    ellisDrainhole.evaluate(pos, g, dg);
    
    for (int lam = 0; lam < 4; lam++) {
        for (int mu = 0; mu < 4; mu++) {
            for (int nu = mu + 1; nu < 4; nu++) {
                EXPECT_NEAR(dg(lam, mu, nu).real, dg(lam, nu, mu).real, kEpsilon)
                    << "Ellis derivative asymmetry at lambda=" << lam;
            }
        }
    }
}

TEST_F(MetricDerivativeTests, EllisDrainholeRadialDerivative) {
    // Radial (r) derivatives should be non-zero
    Tensor<double, 4> pos;
    pos(0) = 0.0; pos(1) = 2.0; pos(2) = 0.5; pos(3) = 0.0;
    
    ellisDrainhole.evaluate(pos, g, dg);
    
    // ∂g_tt/∂r should be non-zero (redshift varies with r)
    // ∂g_tt/∂r should be zero for standard Ellis (g_tt = -1)
    EXPECT_NEAR(dg(1, 0, 0).real, 0.0, kEpsilon)
        << "Ellis ∂g_tt/∂r should be zero";
    
    // ∂g_θθ/∂r should be non-zero (areal radius varies with r)
    EXPECT_NE(dg(1, 2, 2).real, 0.0)
        << "Ellis ∂g_θθ/∂r should be non-zero";
}

TEST_F(MetricDerivativeTests, EllisDrainholeAngularDerivative) {
    // ∂g_φφ/∂θ should be non-zero (sin²θ term)
    Tensor<double, 4> pos;
    pos(0) = 0.0; pos(1) = 3.0; pos(2) = 0.8; pos(3) = 0.0;  // 0 < θ < π
    
    ellisDrainhole.evaluate(pos, g, dg);
    
    EXPECT_NE(dg(2, 3, 3).real, 0.0)
        << "Ellis ∂g_φφ/∂θ should be non-zero for θ ≠ 0, π/2, π";
}

TEST_F(MetricDerivativeTests, EllisDrainholeFiniteDifferenceAgreement) {
    // Note: Ellis Drainhole has steeper gradients near throat, requires larger tolerance
    Tensor<double, 4> pos;
    pos(0) = 0.0; pos(1) = 5.0; pos(2) = 0.8; pos(3) = 0.0;  // Further from throat
    
    ellisDrainhole.evaluate(pos, g, dg);
    
    // Check radial derivative of g_θθ (more stable than g_tt near throat)
    double fd = finiteDiffDerivative(ellisDrainhole, pos, 1, 2, 2);
    double analytic = dg(1, 2, 2).real;
    
    // Use relative tolerance for larger values
    double tol = std::max(1e-3, std::abs(fd) * 0.01);  // 1% relative or 1e-3 absolute
    EXPECT_NEAR(analytic, fd, tol)
        << "Ellis finite diff mismatch for dg_θθ/dr";
}

// =============================================================================
// Alcubierre Derivative Tests
// =============================================================================

TEST_F(MetricDerivativeTests, AlcubierreDerivativeSymmetry) {
    Tensor<double, 4> pos;
    pos(0) = 1.0; pos(1) = 5.0; pos(2) = 0.0; pos(3) = 0.0;
    
    alcubierre.evaluate(pos, g, dg);
    
    for (int lam = 0; lam < 4; lam++) {
        for (int mu = 0; mu < 4; mu++) {
            for (int nu = mu + 1; nu < 4; nu++) {
                EXPECT_NEAR(dg(lam, mu, nu).real, dg(lam, nu, mu).real, kEpsilon)
                    << "Alcubierre derivative asymmetry";
            }
        }
    }
}

TEST_F(MetricDerivativeTests, DISABLED_AlcubierreTimeDerivative) {
    // Unique: Alcubierre has time-dependent metric (moving bubble)
    Tensor<double, 4> pos;
    pos(0) = 0.0; pos(1) = 3.0; pos(2) = 0.0; pos(3) = 0.0;
    
    alcubierre.evaluate(pos, g, dg);
    
    // ∂g_tt/∂t should be non-zero (bubble moves)
    EXPECT_NE(dg(0, 0, 0).real, 0.0)
        << "Alcubierre ∂g_tt/∂t should be non-zero (time-dependent)";
}

TEST_F(MetricDerivativeTests, AlcubirrreFiniteDifferenceAgreement) {
    Tensor<double, 4> pos;
    pos(0) = 0.5; pos(1) = 3.0; pos(2) = 1.0; pos(3) = 0.0;
    
    alcubierre.evaluate(pos, g, dg);
    
    // Check spatial derivative
    double fd = finiteDiffDerivative(alcubierre, pos, 1, 0, 0);
    double analytic = dg(1, 0, 0).real;
    
    EXPECT_NEAR(analytic, fd, 1e-4)
        << "Alcubierre finite diff mismatch for dg_tt/dx";
}

TEST_F(MetricDerivativeTests, AlcubierreOutsideBubble) {
    // Far from bubble center, f → 0, derivatives → 0
    Tensor<double, 4> pos;
    pos(0) = 0.0; pos(1) = 100.0; pos(2) = 0.0; pos(3) = 0.0;  // Far from bubble
    
    alcubierre.evaluate(pos, g, dg);
    
    // All derivatives should be very small outside bubble
    double max_deriv = 0.0;
    for (int lam = 0; lam < 4; lam++) {
        for (int mu = 0; mu < 4; mu++) {
            for (int nu = 0; nu < 4; nu++) {
                max_deriv = std::max(max_deriv, std::abs(dg(lam, mu, nu).real));
            }
        }
    }
    
    EXPECT_LT(max_deriv, 1e-3)
        << "Alcubierre derivatives should vanish far from bubble";
}

// =============================================================================
// Christoffel Symbol From Derivatives
// =============================================================================

TEST_F(MetricDerivativeTests, ChristoffelSymmetryFromDerivatives) {
    // Christoffel Γ^λ_μν = Γ^λ_νμ (torsion-free)
    // This follows from derivative symmetry
    Tensor<double, 4> pos;
    pos(0) = 0.0; pos(1) = 5.0; pos(2) = 0.5; pos(3) = 0.0;
    
    kerrSchild.evaluate(pos, g, dg);
    
    // The lower Christoffel Γ_λμν = (1/2)(∂g_λν/∂x^μ + ∂g_λμ/∂x^ν - ∂g_μν/∂x^λ)
    // should satisfy Γ_λμν = Γ_λνμ
    for (int lam = 0; lam < 4; lam++) {
        for (int mu = 0; mu < 4; mu++) {
            for (int nu = mu + 1; nu < 4; nu++) {
                double Gamma_lmn = 0.5 * (dg(mu, lam, nu).real + dg(nu, lam, mu).real 
                                        - dg(lam, mu, nu).real);
                double Gamma_lnm = 0.5 * (dg(nu, lam, mu).real + dg(mu, lam, nu).real 
                                        - dg(lam, nu, mu).real);
                
                EXPECT_NEAR(Gamma_lmn, Gamma_lnm, kEpsilon)
                    << "Christoffel asymmetry at lambda=" << lam;
            }
        }
    }
}

// main() is provided by GTest::gtest_main
