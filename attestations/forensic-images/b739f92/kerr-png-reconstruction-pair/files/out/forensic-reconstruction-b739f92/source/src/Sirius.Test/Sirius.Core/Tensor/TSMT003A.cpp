// TSMT003A.cpp - Tensor Inverse and Determinant Tests
// Component ID: TSMT003A
// Tests for: MTTN001A (inverse, determinant), PHMT100A (analytic inverse)
//
// Regression suite for the corrected general 4x4 inverse. The previous
// implementation zeroed every off-diagonal component of the inverse and
// mis-computed the diagonal for any metric with off-diagonal terms, which
// includes every Kerr-Schild metric (eta + H l otimes l has off-diagonal
// terms throughout the domain). No earlier test could catch this because the
// conservation suite exercised the separate double-precision oracle stack.
//
// Two exact identities pin correctness for the Kerr-Schild family:
//   det g = det eta = -1        (l is null, so the rank-one update is unimodular)
//   g^{mu nu} = eta^{mu nu} - H l^mu l^nu
// and the generic contraction g^{mu alpha} g_{alpha nu} = delta^mu_nu holds
// for any non-degenerate metric.

#include <gtest/gtest.h>
#include <MTTN001A.h>
#include <PHMT100A.h>
#include <cmath>

using namespace Sirius;

namespace {

// Contract g_inv with g and return the maximum deviation from the identity.
double inverseDeviation(const Metric4D& g, const Metric4D& g_inv) {
    double worst = 0.0;
    for (int mu = 0; mu < 4; mu++) {
        for (int nu = 0; nu < 4; nu++) {
            double sum = 0.0;
            for (int alpha = 0; alpha < 4; alpha++) {
                sum += g_inv(mu, alpha).real * g(alpha, nu).real;
            }
            double expected = (mu == nu) ? 1.0 : 0.0;
            worst = std::max(worst, std::abs(sum - expected));
        }
    }
    return worst;
}

Metric4D evaluateAt(KerrSchildFamily& metric, double x, double y, double z) {
    Metric4D g;
    Tensor<Dual<double>, 4, 4, 4> dg;
    Tensor<double, 4> pos;
    pos(0) = 0.0; pos(1) = x; pos(2) = y; pos(3) = z;
    metric.evaluate(pos, g, dg);
    return g;
}

} // namespace

TEST(TensorInverseTests, DiagonalMetricInverseIsExact) {
    Metric4D g;
    g(0, 0) = Dual<double>(-2.0);
    g(1, 1) = Dual<double>(4.0);
    g(2, 2) = Dual<double>(0.5);
    g(3, 3) = Dual<double>(8.0);

    Metric4D g_inv = TensorOps::inverse(g);
    EXPECT_LT(inverseDeviation(g, g_inv), 1e-14);
    EXPECT_NEAR(TensorOps::determinant(g), -32.0, 1e-12);
}

TEST(TensorInverseTests, OffDiagonalMetricInverseIsExact) {
    // A symmetric Lorentzian-like matrix with every off-diagonal populated.
    // The defective implementation returned zero off-diagonals here.
    Metric4D g;
    const double m[4][4] = {
        {-1.2,  0.3,  0.1, -0.2},
        { 0.3,  1.5, -0.4,  0.2},
        { 0.1, -0.4,  2.0,  0.5},
        {-0.2,  0.2,  0.5,  1.1},
    };
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            g(i, j) = Dual<double>(m[i][j]);

    Metric4D g_inv = TensorOps::inverse(g);
    EXPECT_LT(inverseDeviation(g, g_inv), 1e-12);

    // The inverse of a symmetric matrix is symmetric, and for this matrix the
    // off-diagonals are genuinely non-zero.
    EXPECT_NEAR(g_inv(0, 1).real, g_inv(1, 0).real, 1e-14);
    EXPECT_GT(std::abs(g_inv(0, 1).real), 1e-3);
}

TEST(TensorInverseTests, KerrSchildDeterminantIsMinusOne) {
    // det(eta + H l l^T) = det(eta) (1 + H l^T eta^-1 l) = -1 exactly,
    // because the Kerr-Schild null vector satisfies eta^{mu nu} l_mu l_nu = 0.
    KerrSchildFamily kerr(KerrSchildParams::Kerr(1.0, 0.9));
    const double points[][3] = {
        {6.0, 0.0, 0.0}, {3.0, 2.0, 1.0}, {0.5, 0.5, 4.0}, {-2.0, 5.0, -1.5},
    };
    for (const auto& p : points) {
        Metric4D g = evaluateAt(kerr, p[0], p[1], p[2]);
        EXPECT_NEAR(TensorOps::determinant(g), -1.0, 1e-10)
            << "at (" << p[0] << ", " << p[1] << ", " << p[2] << ")";
    }
}

TEST(TensorInverseTests, KerrSchildNumericInverseSatisfiesIdentity) {
    KerrSchildFamily kn(KerrSchildParams::KerrNewman(1.0, 0.7, 0.3));
    Metric4D g = evaluateAt(kn, 4.0, 1.0, 2.0);
    Metric4D g_inv = TensorOps::inverse(g);
    EXPECT_LT(inverseDeviation(g, g_inv), 1e-10);
}

TEST(TensorInverseTests, AnalyticInverseMatchesNumericInverse) {
    KerrSchildFamily kerr(KerrSchildParams::Kerr(1.0, 0.9));
    const double points[][3] = {
        {6.0, 0.0, 0.0}, {3.0, 2.0, 1.0}, {2.0, -1.0, 0.5},
    };
    for (const auto& p : points) {
        Tensor<double, 4> pos;
        pos(0) = 0.0; pos(1) = p[0]; pos(2) = p[1]; pos(3) = p[2];

        Metric4D g = evaluateAt(kerr, p[0], p[1], p[2]);
        Metric4D analytic;
        ASSERT_TRUE(kerr.inverseMetric(pos, analytic));
        EXPECT_LT(inverseDeviation(g, analytic), 1e-10)
            << "analytic inverse fails the identity at ("
            << p[0] << ", " << p[1] << ", " << p[2] << ")";

        Metric4D numeric = TensorOps::inverse(g);
        for (int mu = 0; mu < 4; mu++) {
            for (int nu = 0; nu < 4; nu++) {
                EXPECT_NEAR(analytic(mu, nu).real, numeric(mu, nu).real, 1e-9);
            }
        }
    }
}

TEST(TensorInverseTests, ChristoffelSymbolsAreSymmetricInLowerIndices) {
    // christoffel() routes through the corrected inverse; the torsion-free
    // symmetry must be exact for the off-axis Kerr point.
    KerrSchildFamily kerr(KerrSchildParams::Kerr(1.0, 0.9));
    Metric4D g;
    Tensor<Dual<double>, 4, 4, 4> dg;
    Tensor<double, 4> pos;
    pos(0) = 0.0; pos(1) = 3.0; pos(2) = 2.0; pos(3) = 1.0;
    kerr.evaluate(pos, g, dg);

    ChristoffelSymbols gamma = TensorOps::christoffel(g, dg);
    for (int lam = 0; lam < 4; lam++) {
        for (int mu = 0; mu < 4; mu++) {
            for (int nu = mu + 1; nu < 4; nu++) {
                EXPECT_NEAR(gamma.gamma(lam, mu, nu).real,
                            gamma.gamma(lam, nu, mu).real, 1e-12);
            }
        }
    }
}
