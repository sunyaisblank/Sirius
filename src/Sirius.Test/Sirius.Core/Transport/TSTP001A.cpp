// TSTP001A.cpp - Double-Precision Transport Type Tests
// Component ID: TSTP001A (Test/Transport/Vec4dTypes)
// Tests: MTTP001A.h (Vec4d, GeodesicStateD, HamiltonianStateD, Mat4d)

#include <gtest/gtest.h>
#include <cmath>
#include <MTTP001A.h>

namespace sirius::test {

constexpr double kEps = 1e-14;

// =============================================================================
// Vec4d Tests
// =============================================================================

class Vec4dTests : public ::testing::Test {};

TEST_F(Vec4dTests, DefaultConstructionAllZeros) {
    sirius::math::Vec4d v;
    EXPECT_DOUBLE_EQ(v.t, 0.0);
    EXPECT_DOUBLE_EQ(v.r, 0.0);
    EXPECT_DOUBLE_EQ(v.theta, 0.0);
    EXPECT_DOUBLE_EQ(v.phi, 0.0);
}

TEST_F(Vec4dTests, ParameterisedConstruction) {
    sirius::math::Vec4d v(1.0, 2.0, 3.0, 4.0);
    EXPECT_DOUBLE_EQ(v.t, 1.0);
    EXPECT_DOUBLE_EQ(v.r, 2.0);
    EXPECT_DOUBLE_EQ(v.theta, 3.0);
    EXPECT_DOUBLE_EQ(v.phi, 4.0);
}

TEST_F(Vec4dTests, IndexedAccess) {
    sirius::math::Vec4d v(10.0, 20.0, 30.0, 40.0);
    EXPECT_DOUBLE_EQ(v[0], 10.0);
    EXPECT_DOUBLE_EQ(v[1], 20.0);
    EXPECT_DOUBLE_EQ(v[2], 30.0);
    EXPECT_DOUBLE_EQ(v[3], 40.0);
}

TEST_F(Vec4dTests, IndexedMutation) {
    sirius::math::Vec4d v;
    v[0] = 5.0;
    v[3] = 8.0;
    EXPECT_DOUBLE_EQ(v.t, 5.0);
    EXPECT_DOUBLE_EQ(v.phi, 8.0);
}

TEST_F(Vec4dTests, Addition) {
    sirius::math::Vec4d a(1.0, 2.0, 3.0, 4.0);
    sirius::math::Vec4d b(0.5, 1.5, 2.5, 3.5);
    auto c = a + b;
    EXPECT_DOUBLE_EQ(c.t, 1.5);
    EXPECT_DOUBLE_EQ(c.r, 3.5);
    EXPECT_DOUBLE_EQ(c.theta, 5.5);
    EXPECT_DOUBLE_EQ(c.phi, 7.5);
}

TEST_F(Vec4dTests, Subtraction) {
    sirius::math::Vec4d a(5.0, 4.0, 3.0, 2.0);
    sirius::math::Vec4d b(1.0, 1.0, 1.0, 1.0);
    auto c = a - b;
    EXPECT_DOUBLE_EQ(c.t, 4.0);
    EXPECT_DOUBLE_EQ(c.r, 3.0);
    EXPECT_DOUBLE_EQ(c.theta, 2.0);
    EXPECT_DOUBLE_EQ(c.phi, 1.0);
}

TEST_F(Vec4dTests, ScalarMultiplication) {
    sirius::math::Vec4d v(1.0, 2.0, 3.0, 4.0);
    auto w = v * 3.0;
    EXPECT_DOUBLE_EQ(w.t, 3.0);
    EXPECT_DOUBLE_EQ(w.r, 6.0);
    EXPECT_DOUBLE_EQ(w.theta, 9.0);
    EXPECT_DOUBLE_EQ(w.phi, 12.0);
}

TEST_F(Vec4dTests, ScalarMultiplicationCommutative) {
    sirius::math::Vec4d v(1.0, 2.0, 3.0, 4.0);
    auto a = v * 2.0;
    auto b = 2.0 * v;
    EXPECT_DOUBLE_EQ(a.t, b.t);
    EXPECT_DOUBLE_EQ(a.r, b.r);
    EXPECT_DOUBLE_EQ(a.theta, b.theta);
    EXPECT_DOUBLE_EQ(a.phi, b.phi);
}

TEST_F(Vec4dTests, ScalarDivision) {
    sirius::math::Vec4d v(10.0, 20.0, 30.0, 40.0);
    auto w = v / 2.0;
    EXPECT_DOUBLE_EQ(w.t, 5.0);
    EXPECT_DOUBLE_EQ(w.r, 10.0);
}

TEST_F(Vec4dTests, CompoundAssignment) {
    sirius::math::Vec4d v(1.0, 2.0, 3.0, 4.0);
    sirius::math::Vec4d w(0.5, 0.5, 0.5, 0.5);
    v += w;
    EXPECT_DOUBLE_EQ(v.t, 1.5);

    v -= w;
    EXPECT_DOUBLE_EQ(v.t, 1.0);

    v *= 2.0;
    EXPECT_DOUBLE_EQ(v.t, 2.0);
}

TEST_F(Vec4dTests, Negation) {
    sirius::math::Vec4d v(1.0, -2.0, 3.0, -4.0);
    auto w = -v;
    EXPECT_DOUBLE_EQ(w.t, -1.0);
    EXPECT_DOUBLE_EQ(w.r, 2.0);
    EXPECT_DOUBLE_EQ(w.theta, -3.0);
    EXPECT_DOUBLE_EQ(w.phi, 4.0);
}

TEST_F(Vec4dTests, Norm2) {
    sirius::math::Vec4d v(1.0, 2.0, 3.0, 4.0);
    EXPECT_DOUBLE_EQ(v.norm2(), 1.0 + 4.0 + 9.0 + 16.0);
}

TEST_F(Vec4dTests, IsZero) {
    sirius::math::Vec4d zero;
    EXPECT_TRUE(zero.isZero());

    sirius::math::Vec4d nonzero(0.0, 0.0, 0.0, 1e-20);
    EXPECT_FALSE(nonzero.isZero());
}

// =============================================================================
// GeodesicStateD Tests
// =============================================================================

class GeodesicStateDTests : public ::testing::Test {};

TEST_F(GeodesicStateDTests, DefaultConstruction) {
    sirius::math::GeodesicStateD g;
    EXPECT_TRUE(g.x.isZero());
    EXPECT_TRUE(g.k.isZero());
    EXPECT_DOUBLE_EQ(g.lambda, 0.0);
    EXPECT_DOUBLE_EQ(g.E, 0.0);
    EXPECT_DOUBLE_EQ(g.Lz, 0.0);
}

TEST_F(GeodesicStateDTests, ConstructionFromPosAndMom) {
    sirius::math::Vec4d pos(0.0, 10.0, 1.5708, 0.0);
    sirius::math::Vec4d mom(-1.0, 0.5, 0.0, 2.0);
    sirius::math::GeodesicStateD g(pos, mom);

    // E = -k_t = -(-1.0) = 1.0
    EXPECT_DOUBLE_EQ(g.E, 1.0);
    // Lz = k_phi = 2.0
    EXPECT_DOUBLE_EQ(g.Lz, 2.0);
}

TEST_F(GeodesicStateDTests, ConservedQuantitiesSchwarzschildLike) {
    sirius::math::Vec4d pos(0.0, 10.0, M_PI / 2.0, 0.0);
    sirius::math::Vec4d mom(-0.9, 0.3, 0.0, 3.0);
    sirius::math::GeodesicStateD g(pos, mom);
    g.computeConservedQuantities(0.0); // a=0 (Schwarzschild)

    EXPECT_DOUBLE_EQ(g.E, 0.9);
    EXPECT_DOUBLE_EQ(g.Lz, 3.0);
}

// =============================================================================
// HamiltonianStateD Tests
// =============================================================================

class HamiltonianStateDTests : public ::testing::Test {};

TEST_F(HamiltonianStateDTests, DefaultConstruction) {
    sirius::math::HamiltonianStateD h;
    EXPECT_TRUE(h.q.isZero());
    EXPECT_TRUE(h.p.isZero());
    EXPECT_DOUBLE_EQ(h.H, 0.0);
}

TEST_F(HamiltonianStateDTests, RoundTripViaGeodesicState) {
    sirius::math::Vec4d q(0.0, 10.0, 1.5, 0.5);
    sirius::math::Vec4d p(-1.0, 0.3, 0.1, 2.0);
    sirius::math::HamiltonianStateD h(q, p);

    auto geo = h.toGeodesicState(5.0);
    EXPECT_NEAR(geo.x.r, 10.0, kEps);
    EXPECT_NEAR(geo.k.phi, 2.0, kEps);
    EXPECT_DOUBLE_EQ(geo.lambda, 5.0);
    EXPECT_DOUBLE_EQ(geo.E, 1.0);
    EXPECT_DOUBLE_EQ(geo.Lz, 2.0);
}

TEST_F(HamiltonianStateDTests, ConstructFromGeodesicState) {
    sirius::math::Vec4d pos(0.0, 8.0, 1.2, 0.8);
    sirius::math::Vec4d mom(-0.5, 0.2, 0.0, 1.5);
    sirius::math::GeodesicStateD geo(pos, mom);
    sirius::math::HamiltonianStateD h(geo);

    EXPECT_NEAR(h.q.r, 8.0, kEps);
    EXPECT_NEAR(h.p.phi, 1.5, kEps);
}

// =============================================================================
// Mat4d Tests
// =============================================================================

class Mat4dTests : public ::testing::Test {};

TEST_F(Mat4dTests, DefaultIsZero) {
    sirius::math::Mat4d m;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            EXPECT_DOUBLE_EQ(m(i, j), 0.0);
}

TEST_F(Mat4dTests, IdentityDiagonal) {
    auto I = sirius::math::Mat4d::identity();
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            EXPECT_DOUBLE_EQ(I(i, j), (i == j) ? 1.0 : 0.0);
        }
    }
}

TEST_F(Mat4dTests, Trace) {
    auto I = sirius::math::Mat4d::identity();
    EXPECT_DOUBLE_EQ(I.trace(), 4.0);

    sirius::math::Mat4d m;
    m(0, 0) = 2.0; m(1, 1) = 3.0; m(2, 2) = 4.0; m(3, 3) = 5.0;
    EXPECT_DOUBLE_EQ(m.trace(), 14.0);
}

TEST_F(Mat4dTests, IdentityDeterminant) {
    auto I = sirius::math::Mat4d::identity();
    EXPECT_NEAR(I.determinant(), 1.0, kEps);
}

TEST_F(Mat4dTests, ZeroDeterminant) {
    auto Z = sirius::math::Mat4d::zero();
    EXPECT_NEAR(Z.determinant(), 0.0, kEps);
}

TEST_F(Mat4dTests, KnownDeterminant) {
    // Minkowski metric: diag(-1, 1, 1, 1), det = -1
    auto eta = sirius::math::Mat4d::identity();
    eta(0, 0) = -1.0;
    EXPECT_NEAR(eta.determinant(), -1.0, kEps);
}

TEST_F(Mat4dTests, MatrixVectorMultiplication) {
    auto I = sirius::math::Mat4d::identity();
    sirius::math::Vec4d v(1.0, 2.0, 3.0, 4.0);
    auto w = I * v;
    EXPECT_DOUBLE_EQ(w.t, 1.0);
    EXPECT_DOUBLE_EQ(w.r, 2.0);
    EXPECT_DOUBLE_EQ(w.theta, 3.0);
    EXPECT_DOUBLE_EQ(w.phi, 4.0);
}

TEST_F(Mat4dTests, ScalarMultiplication) {
    auto I = sirius::math::Mat4d::identity();
    auto M = I * 3.0;
    EXPECT_DOUBLE_EQ(M(0, 0), 3.0);
    EXPECT_DOUBLE_EQ(M(1, 1), 3.0);
    EXPECT_DOUBLE_EQ(M(0, 1), 0.0);
}

TEST_F(Mat4dTests, Addition) {
    auto I = sirius::math::Mat4d::identity();
    auto M = I + I;
    EXPECT_DOUBLE_EQ(M(0, 0), 2.0);
    EXPECT_DOUBLE_EQ(M(1, 1), 2.0);
    EXPECT_DOUBLE_EQ(M(0, 1), 0.0);
}

// =============================================================================
// Conversion Utilities
// =============================================================================

class ConversionTests : public ::testing::Test {};

TEST_F(ConversionTests, ToVec4dFromFloats) {
    auto v = sirius::math::toVec4d(1.0f, 2.0f, 3.0f, 4.0f);
    EXPECT_DOUBLE_EQ(v.t, 1.0);
    EXPECT_DOUBLE_EQ(v.r, 2.0);
}

TEST_F(ConversionTests, ToFloat4RoundTrip) {
    sirius::math::Vec4d v(1.5, 2.5, 3.5, 4.5);
    float t, r, theta, phi;
    sirius::math::toFloat4(v, t, r, theta, phi);
    EXPECT_FLOAT_EQ(t, 1.5f);
    EXPECT_FLOAT_EQ(r, 2.5f);
    EXPECT_FLOAT_EQ(theta, 3.5f);
    EXPECT_FLOAT_EQ(phi, 4.5f);
}

} // namespace sirius::test
