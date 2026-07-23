// Stokes vector and Mueller matrix tests. Ported from TSPL001A.cpp; assertions
// and tolerances unchanged. Tests stokes.h (StokesVector, MuellerMatrix,
// polarised_emission, parallel_transport).

#include "sirius/core/polarisation/stokes.h"

#include <gtest/gtest.h>

#include <cmath>

namespace sirius::test {
using namespace sirius::core;

constexpr float kEps = 1e-5f;
constexpr float kPi = 3.14159265358979f;

// =============================================================================
// StokesVector Tests
// =============================================================================

class StokesVectorTests : public ::testing::Test {};

TEST_F(StokesVectorTests, UnpolarisedConstruction) {
    auto s = sirius::core::StokesVector::Unpolarised(2.0f);
    EXPECT_FLOAT_EQ(s.I, 2.0f);
    EXPECT_FLOAT_EQ(s.Q, 0.0f);
    EXPECT_FLOAT_EQ(s.U, 0.0f);
    EXPECT_FLOAT_EQ(s.V, 0.0f);
}

TEST_F(StokesVectorTests, HorizontalPolarisation) {
    auto s = sirius::core::StokesVector::Horizontal(1.0f);
    EXPECT_FLOAT_EQ(s.I, 1.0f);
    EXPECT_FLOAT_EQ(s.Q, 1.0f);
    EXPECT_FLOAT_EQ(s.U, 0.0f);
    EXPECT_FLOAT_EQ(s.V, 0.0f);
    EXPECT_TRUE(s.IsPhysical());
}

TEST_F(StokesVectorTests, VerticalPolarisation) {
    auto s = sirius::core::StokesVector::Vertical(1.0f);
    EXPECT_FLOAT_EQ(s.Q, -1.0f);
    EXPECT_TRUE(s.IsPhysical());
}

TEST_F(StokesVectorTests, CircularPolarisation) {
    auto r = sirius::core::StokesVector::RightCircular(1.0f);
    auto l = sirius::core::StokesVector::LeftCircular(1.0f);
    EXPECT_FLOAT_EQ(r.V, 1.0f);
    EXPECT_FLOAT_EQ(l.V, -1.0f);
    EXPECT_TRUE(r.IsPhysical());
    EXPECT_TRUE(l.IsPhysical());
}

TEST_F(StokesVectorTests, PhysicalConstraint) {
    // Q^2 + U^2 + V^2 <= I^2 for all factory methods
    auto check = [](const sirius::core::StokesVector& s) {
        float p2 = s.Q * s.Q + s.U * s.U + s.V * s.V;
        EXPECT_LE(p2, s.I * s.I * 1.001f);
    };
    check(sirius::core::StokesVector::Unpolarised(1.0f));
    check(sirius::core::StokesVector::Horizontal(1.0f));
    check(sirius::core::StokesVector::Vertical(1.0f));
    check(sirius::core::StokesVector::Diagonal45(1.0f));
    check(sirius::core::StokesVector::RightCircular(1.0f));
    check(sirius::core::StokesVector::LeftCircular(1.0f));
}

TEST_F(StokesVectorTests, PolarisationDegree) {
    EXPECT_NEAR(sirius::core::StokesVector::Unpolarised(1.0f).PolarisationDegree(), 0.0f, kEps);
    EXPECT_NEAR(sirius::core::StokesVector::Horizontal(1.0f).PolarisationDegree(), 1.0f, kEps);
    EXPECT_NEAR(sirius::core::StokesVector::RightCircular(1.0f).PolarisationDegree(), 1.0f, kEps);
}

TEST_F(StokesVectorTests, LinearPolarisationDegree) {
    EXPECT_NEAR(sirius::core::StokesVector::Horizontal(1.0f).LinearPolarisationDegree(), 1.0f,
                kEps);
    EXPECT_NEAR(sirius::core::StokesVector::RightCircular(1.0f).LinearPolarisationDegree(), 0.0f,
                kEps);
}

TEST_F(StokesVectorTests, CircularPolarisationDegree) {
    EXPECT_NEAR(sirius::core::StokesVector::RightCircular(1.0f).CircularPolarisationDegree(), 1.0f,
                kEps);
    EXPECT_NEAR(sirius::core::StokesVector::Horizontal(1.0f).CircularPolarisationDegree(), 0.0f,
                kEps);
    EXPECT_NEAR(sirius::core::StokesVector::LeftCircular(1.0f).CircularPolarisationDegree(), -1.0f,
                kEps);
}

TEST_F(StokesVectorTests, EVPAComputation) {
    // Horizontal: EVPA = 0.5 * atan2(0, 1) = 0
    EXPECT_NEAR(sirius::core::StokesVector::Horizontal(1.0f).Evpa(), 0.0f, kEps);
    // +45 degree: EVPA = 0.5 * atan2(1, 0) = pi/4
    EXPECT_NEAR(sirius::core::StokesVector::Diagonal45(1.0f).Evpa(), kPi / 4.0f, kEps);
    // Vertical: EVPA = 0.5 * atan2(0, -1) = pi/2
    EXPECT_NEAR(sirius::core::StokesVector::Vertical(1.0f).Evpa(), kPi / 2.0f, kEps);
}

TEST_F(StokesVectorTests, IsPhysicalRejectsUnphysical) {
    sirius::core::StokesVector bad(1.0f, 2.0f, 0.0f, 0.0f);  // Q^2 > I^2
    EXPECT_FALSE(bad.IsPhysical());
}

TEST_F(StokesVectorTests, NormalisationProjection) {
    sirius::core::StokesVector s(1.0f, 2.0f, 0.0f, 0.0f);
    s.Normalise();
    EXPECT_TRUE(s.IsPhysical());
    EXPECT_LE(s.PolarisationDegree(), 1.0f + kEps);
}

TEST_F(StokesVectorTests, ZeroIntensityHandling) {
    sirius::core::StokesVector s(0.0f, 0.0f, 0.0f, 0.0f);
    EXPECT_FLOAT_EQ(s.PolarisationDegree(), 0.0f);
    EXPECT_TRUE(s.IsPhysical());
    s.Normalise();
    EXPECT_FLOAT_EQ(s.I, 0.0f);
}

TEST_F(StokesVectorTests, ArithmeticOperators) {
    auto a = sirius::core::StokesVector::Horizontal(1.0f);
    auto b = sirius::core::StokesVector::Vertical(1.0f);
    auto sum = a + b;
    EXPECT_FLOAT_EQ(sum.I, 2.0f);
    EXPECT_FLOAT_EQ(sum.Q, 0.0f);  // +1 + -1

    auto scaled = a * 3.0f;
    EXPECT_FLOAT_EQ(scaled.I, 3.0f);
    EXPECT_FLOAT_EQ(scaled.Q, 3.0f);
}

// =============================================================================
// MuellerMatrix Tests
// =============================================================================

class MuellerMatrixTests : public ::testing::Test {};

TEST_F(MuellerMatrixTests, IdentityPreservesStokes) {
    sirius::core::MuellerMatrix M;  // default = identity
    auto s = sirius::core::StokesVector(1.0f, 0.5f, 0.3f, 0.1f);
    auto out = M.Apply(s);
    EXPECT_NEAR(out.I, s.I, kEps);
    EXPECT_NEAR(out.Q, s.Q, kEps);
    EXPECT_NEAR(out.U, s.U, kEps);
    EXPECT_NEAR(out.V, s.V, kEps);
}

TEST_F(MuellerMatrixTests, HorizontalPolariserOnUnpolarised) {
    auto M = sirius::core::MuellerMatrix::HorizontalPolariser();
    auto s = sirius::core::StokesVector::Unpolarised(1.0f);
    auto out = M.Apply(s);
    EXPECT_NEAR(out.I, 0.5f, kEps);
    EXPECT_NEAR(out.Q, 0.5f, kEps);
    EXPECT_NEAR(out.U, 0.0f, kEps);
    EXPECT_NEAR(out.V, 0.0f, kEps);
}

TEST_F(MuellerMatrixTests, VerticalPolariserOnUnpolarised) {
    auto M = sirius::core::MuellerMatrix::VerticalPolariser();
    auto s = sirius::core::StokesVector::Unpolarised(1.0f);
    auto out = M.Apply(s);
    EXPECT_NEAR(out.I, 0.5f, kEps);
    EXPECT_NEAR(out.Q, -0.5f, kEps);
}

TEST_F(MuellerMatrixTests, CrossedPolariersExtinguish) {
    // Horizontal then vertical = zero output
    auto H = sirius::core::MuellerMatrix::HorizontalPolariser();
    auto V = sirius::core::MuellerMatrix::VerticalPolariser();
    auto M = V * H;
    auto s = sirius::core::StokesVector::Unpolarised(1.0f);
    auto out = M.Apply(s);
    EXPECT_NEAR(out.I, 0.0f, kEps);
}

TEST_F(MuellerMatrixTests, QuarterWavePlateConvertsToCircular) {
    auto QWP = sirius::core::MuellerMatrix::QuarterWavePlate();
    auto s = sirius::core::StokesVector::Diagonal45(1.0f);  // (1, 0, 1, 0)
    auto out = QWP.Apply(s);
    // +45 linear through QWP(horizontal fast axis) → right circular
    EXPECT_NEAR(out.I, 1.0f, kEps);
    EXPECT_NEAR(out.Q, 0.0f, kEps);
    EXPECT_NEAR(out.U, 0.0f, kEps);
    EXPECT_NEAR(out.V, -1.0f, kEps);
}

TEST_F(MuellerMatrixTests, CompositionAssociativity) {
    auto M1 = sirius::core::MuellerMatrix::Rotation(0.3f);
    auto M2 = sirius::core::MuellerMatrix::HorizontalPolariser();
    auto s = sirius::core::StokesVector(1.0f, 0.5f, 0.3f, 0.1f);

    // (M1 * M2).Apply(s) should equal M1.Apply(M2.Apply(s))
    auto composed = (M1 * M2).Apply(s);
    auto sequential = M1.Apply(M2.Apply(s));

    EXPECT_NEAR(composed.I, sequential.I, kEps);
    EXPECT_NEAR(composed.Q, sequential.Q, kEps);
    EXPECT_NEAR(composed.U, sequential.U, kEps);
    EXPECT_NEAR(composed.V, sequential.V, kEps);
}

TEST_F(MuellerMatrixTests, DepolariserReducesPolarisation) {
    auto M = sirius::core::MuellerMatrix::Depolariser(0.0f);  // full depolarisation
    auto s = sirius::core::StokesVector::Horizontal(1.0f);
    auto out = M.Apply(s);
    EXPECT_NEAR(out.I, 1.0f, kEps);
    EXPECT_NEAR(out.Q, 0.0f, kEps);
    EXPECT_NEAR(out.U, 0.0f, kEps);
    EXPECT_NEAR(out.V, 0.0f, kEps);
}

TEST_F(MuellerMatrixTests, HalfWavePlateFlipsHandedness) {
    auto HWP = sirius::core::MuellerMatrix::HalfWavePlate();
    auto s = sirius::core::StokesVector::RightCircular(1.0f);  // (1, 0, 0, 1)
    auto out = HWP.Apply(s);
    EXPECT_NEAR(out.V, -1.0f, kEps);  // right → left
}

// =============================================================================
// polarised_emission Tests
// =============================================================================

class PolarisedEmissionTests : public ::testing::Test {};

TEST_F(PolarisedEmissionTests, SynchrotronPolarisationDegree) {
    // p=2: (2+1)/(2+7/3) = 3/(13/3) = 9/13 ≈ 0.6923
    float pi_L = sirius::core::polarised_emission::SynchrotronPolarisationDegree(2.0f);
    EXPECT_NEAR(pi_L, 9.0f / 13.0f, kEps);

    // p=3: (3+1)/(3+7/3) = 4/(16/3) = 12/16 = 0.75
    pi_L = sirius::core::polarised_emission::SynchrotronPolarisationDegree(3.0f);
    EXPECT_NEAR(pi_L, 0.75f, kEps);
}

TEST_F(PolarisedEmissionTests, SynchrotronEmissionIsPhysical) {
    auto s = sirius::core::polarised_emission::SynchrotronEmission(1.0f, 0.7f, 0.0f);
    EXPECT_TRUE(s.IsPhysical());
    EXPECT_NEAR(s.V, 0.0f, kEps);  // no circular from synchrotron
}

TEST_F(PolarisedEmissionTests, ThomsonScatteringAt90Degrees) {
    // At 90 degrees (cos_theta = 0): pi = sin^2(90)/(1+cos^2(90)) = 1/1 = 1
    float pi = sirius::core::polarised_emission::ThomsonPolarisationDegree(0.0f);
    EXPECT_NEAR(pi, 1.0f, kEps);
}

TEST_F(PolarisedEmissionTests, ThomsonScatteringForward) {
    // At 0 degrees (cos_theta = 1): pi = 0/(1+1) = 0
    float pi = sirius::core::polarised_emission::ThomsonPolarisationDegree(1.0f);
    EXPECT_NEAR(pi, 0.0f, kEps);
}

// =============================================================================
// parallel_transport Tests
// =============================================================================

class ParallelTransportTests : public ::testing::Test {};

TEST_F(ParallelTransportTests, ZeroSpinNoRotation) {
    float angle =
        sirius::core::parallel_transport::GravitationalFaradayRotation(0.0f, 10.0f, kPi / 2.0f);
    EXPECT_NEAR(angle, 0.0f, kEps);
}

TEST_F(ParallelTransportTests, RotationIncreasesWithSpin) {
    float a1 =
        sirius::core::parallel_transport::GravitationalFaradayRotation(0.5f, 10.0f, kPi / 2.0f);
    float a2 =
        sirius::core::parallel_transport::GravitationalFaradayRotation(0.9f, 10.0f, kPi / 2.0f);
    EXPECT_GT(std::abs(a2), std::abs(a1));
}

TEST_F(ParallelTransportTests, ApplyPreservesIntensity) {
    auto s = sirius::core::StokesVector::Horizontal(1.0f);
    auto out = sirius::core::parallel_transport::ApplyParallelTransport(s, 0.5f);
    EXPECT_NEAR(out.I, s.I, kEps);  // rotation preserves I
}

}  // namespace sirius::test
