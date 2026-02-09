// TSPL001A.cpp - Stokes Vector and Mueller Matrix Tests
// Component ID: TSPL001A (Test/Polarisation/StokesVector)
// Tests: PHPL001A.h (StokesVector, MuellerMatrix, PolarisedEmission, ParallelTransport)

#include <gtest/gtest.h>
#include <cmath>
#include <PHPL001A.h>

namespace sirius::test {

constexpr float kEps = 1e-5f;
constexpr float kPi = 3.14159265358979f;

// =============================================================================
// StokesVector Tests
// =============================================================================

class StokesVectorTests : public ::testing::Test {};

TEST_F(StokesVectorTests, UnpolarisedConstruction) {
    auto s = Sirius::StokesVector::unpolarised(2.0f);
    EXPECT_FLOAT_EQ(s.I, 2.0f);
    EXPECT_FLOAT_EQ(s.Q, 0.0f);
    EXPECT_FLOAT_EQ(s.U, 0.0f);
    EXPECT_FLOAT_EQ(s.V, 0.0f);
}

TEST_F(StokesVectorTests, HorizontalPolarisation) {
    auto s = Sirius::StokesVector::horizontal(1.0f);
    EXPECT_FLOAT_EQ(s.I, 1.0f);
    EXPECT_FLOAT_EQ(s.Q, 1.0f);
    EXPECT_FLOAT_EQ(s.U, 0.0f);
    EXPECT_FLOAT_EQ(s.V, 0.0f);
    EXPECT_TRUE(s.isPhysical());
}

TEST_F(StokesVectorTests, VerticalPolarisation) {
    auto s = Sirius::StokesVector::vertical(1.0f);
    EXPECT_FLOAT_EQ(s.Q, -1.0f);
    EXPECT_TRUE(s.isPhysical());
}

TEST_F(StokesVectorTests, CircularPolarisation) {
    auto r = Sirius::StokesVector::rightCircular(1.0f);
    auto l = Sirius::StokesVector::leftCircular(1.0f);
    EXPECT_FLOAT_EQ(r.V, 1.0f);
    EXPECT_FLOAT_EQ(l.V, -1.0f);
    EXPECT_TRUE(r.isPhysical());
    EXPECT_TRUE(l.isPhysical());
}

TEST_F(StokesVectorTests, PhysicalConstraint) {
    // Q^2 + U^2 + V^2 <= I^2 for all factory methods
    auto check = [](const Sirius::StokesVector& s) {
        float p2 = s.Q * s.Q + s.U * s.U + s.V * s.V;
        EXPECT_LE(p2, s.I * s.I * 1.001f);
    };
    check(Sirius::StokesVector::unpolarised(1.0f));
    check(Sirius::StokesVector::horizontal(1.0f));
    check(Sirius::StokesVector::vertical(1.0f));
    check(Sirius::StokesVector::diagonal45(1.0f));
    check(Sirius::StokesVector::rightCircular(1.0f));
    check(Sirius::StokesVector::leftCircular(1.0f));
}

TEST_F(StokesVectorTests, PolarisationDegree) {
    EXPECT_NEAR(Sirius::StokesVector::unpolarised(1.0f).polarisationDegree(), 0.0f, kEps);
    EXPECT_NEAR(Sirius::StokesVector::horizontal(1.0f).polarisationDegree(), 1.0f, kEps);
    EXPECT_NEAR(Sirius::StokesVector::rightCircular(1.0f).polarisationDegree(), 1.0f, kEps);
}

TEST_F(StokesVectorTests, LinearPolarisationDegree) {
    EXPECT_NEAR(Sirius::StokesVector::horizontal(1.0f).linearPolarisationDegree(), 1.0f, kEps);
    EXPECT_NEAR(Sirius::StokesVector::rightCircular(1.0f).linearPolarisationDegree(), 0.0f, kEps);
}

TEST_F(StokesVectorTests, CircularPolarisationDegree) {
    EXPECT_NEAR(Sirius::StokesVector::rightCircular(1.0f).circularPolarisationDegree(), 1.0f, kEps);
    EXPECT_NEAR(Sirius::StokesVector::horizontal(1.0f).circularPolarisationDegree(), 0.0f, kEps);
    EXPECT_NEAR(Sirius::StokesVector::leftCircular(1.0f).circularPolarisationDegree(), -1.0f, kEps);
}

TEST_F(StokesVectorTests, EVPAComputation) {
    // Horizontal: EVPA = 0.5 * atan2(0, 1) = 0
    EXPECT_NEAR(Sirius::StokesVector::horizontal(1.0f).EVPA(), 0.0f, kEps);
    // +45 degree: EVPA = 0.5 * atan2(1, 0) = pi/4
    EXPECT_NEAR(Sirius::StokesVector::diagonal45(1.0f).EVPA(), kPi / 4.0f, kEps);
    // Vertical: EVPA = 0.5 * atan2(0, -1) = pi/2
    EXPECT_NEAR(Sirius::StokesVector::vertical(1.0f).EVPA(), kPi / 2.0f, kEps);
}

TEST_F(StokesVectorTests, IsPhysicalRejectsUnphysical) {
    Sirius::StokesVector bad(1.0f, 2.0f, 0.0f, 0.0f); // Q^2 > I^2
    EXPECT_FALSE(bad.isPhysical());
}

TEST_F(StokesVectorTests, NormalisationProjection) {
    Sirius::StokesVector s(1.0f, 2.0f, 0.0f, 0.0f);
    s.normalise();
    EXPECT_TRUE(s.isPhysical());
    EXPECT_LE(s.polarisationDegree(), 1.0f + kEps);
}

TEST_F(StokesVectorTests, ZeroIntensityHandling) {
    Sirius::StokesVector s(0.0f, 0.0f, 0.0f, 0.0f);
    EXPECT_FLOAT_EQ(s.polarisationDegree(), 0.0f);
    EXPECT_TRUE(s.isPhysical());
    s.normalise();
    EXPECT_FLOAT_EQ(s.I, 0.0f);
}

TEST_F(StokesVectorTests, ArithmeticOperators) {
    auto a = Sirius::StokesVector::horizontal(1.0f);
    auto b = Sirius::StokesVector::vertical(1.0f);
    auto sum = a + b;
    EXPECT_FLOAT_EQ(sum.I, 2.0f);
    EXPECT_FLOAT_EQ(sum.Q, 0.0f); // +1 + -1

    auto scaled = a * 3.0f;
    EXPECT_FLOAT_EQ(scaled.I, 3.0f);
    EXPECT_FLOAT_EQ(scaled.Q, 3.0f);
}

// =============================================================================
// MuellerMatrix Tests
// =============================================================================

class MuellerMatrixTests : public ::testing::Test {};

TEST_F(MuellerMatrixTests, IdentityPreservesStokes) {
    Sirius::MuellerMatrix M; // default = identity
    auto s = Sirius::StokesVector(1.0f, 0.5f, 0.3f, 0.1f);
    auto out = M.apply(s);
    EXPECT_NEAR(out.I, s.I, kEps);
    EXPECT_NEAR(out.Q, s.Q, kEps);
    EXPECT_NEAR(out.U, s.U, kEps);
    EXPECT_NEAR(out.V, s.V, kEps);
}

TEST_F(MuellerMatrixTests, HorizontalPolariserOnUnpolarised) {
    auto M = Sirius::MuellerMatrix::horizontalPolariser();
    auto s = Sirius::StokesVector::unpolarised(1.0f);
    auto out = M.apply(s);
    EXPECT_NEAR(out.I, 0.5f, kEps);
    EXPECT_NEAR(out.Q, 0.5f, kEps);
    EXPECT_NEAR(out.U, 0.0f, kEps);
    EXPECT_NEAR(out.V, 0.0f, kEps);
}

TEST_F(MuellerMatrixTests, VerticalPolariserOnUnpolarised) {
    auto M = Sirius::MuellerMatrix::verticalPolariser();
    auto s = Sirius::StokesVector::unpolarised(1.0f);
    auto out = M.apply(s);
    EXPECT_NEAR(out.I, 0.5f, kEps);
    EXPECT_NEAR(out.Q, -0.5f, kEps);
}

TEST_F(MuellerMatrixTests, CrossedPolariersExtinguish) {
    // Horizontal then vertical = zero output
    auto H = Sirius::MuellerMatrix::horizontalPolariser();
    auto V = Sirius::MuellerMatrix::verticalPolariser();
    auto M = V * H;
    auto s = Sirius::StokesVector::unpolarised(1.0f);
    auto out = M.apply(s);
    EXPECT_NEAR(out.I, 0.0f, kEps);
}

TEST_F(MuellerMatrixTests, QuarterWavePlateConvertsToCircular) {
    auto QWP = Sirius::MuellerMatrix::quarterWavePlate();
    auto s = Sirius::StokesVector::diagonal45(1.0f); // (1, 0, 1, 0)
    auto out = QWP.apply(s);
    // +45 linear through QWP(horizontal fast axis) → right circular
    EXPECT_NEAR(out.I, 1.0f, kEps);
    EXPECT_NEAR(out.Q, 0.0f, kEps);
    EXPECT_NEAR(out.U, 0.0f, kEps);
    EXPECT_NEAR(out.V, -1.0f, kEps);
}

TEST_F(MuellerMatrixTests, CompositionAssociativity) {
    auto M1 = Sirius::MuellerMatrix::rotation(0.3f);
    auto M2 = Sirius::MuellerMatrix::horizontalPolariser();
    auto s = Sirius::StokesVector(1.0f, 0.5f, 0.3f, 0.1f);

    // (M1 * M2).apply(s) should equal M1.apply(M2.apply(s))
    auto composed = (M1 * M2).apply(s);
    auto sequential = M1.apply(M2.apply(s));

    EXPECT_NEAR(composed.I, sequential.I, kEps);
    EXPECT_NEAR(composed.Q, sequential.Q, kEps);
    EXPECT_NEAR(composed.U, sequential.U, kEps);
    EXPECT_NEAR(composed.V, sequential.V, kEps);
}

TEST_F(MuellerMatrixTests, DepolariserReducesPolarisation) {
    auto M = Sirius::MuellerMatrix::depolariser(0.0f); // full depolarisation
    auto s = Sirius::StokesVector::horizontal(1.0f);
    auto out = M.apply(s);
    EXPECT_NEAR(out.I, 1.0f, kEps);
    EXPECT_NEAR(out.Q, 0.0f, kEps);
    EXPECT_NEAR(out.U, 0.0f, kEps);
    EXPECT_NEAR(out.V, 0.0f, kEps);
}

TEST_F(MuellerMatrixTests, HalfWavePlateFlipsHandedness) {
    auto HWP = Sirius::MuellerMatrix::halfWavePlate();
    auto s = Sirius::StokesVector::rightCircular(1.0f); // (1, 0, 0, 1)
    auto out = HWP.apply(s);
    EXPECT_NEAR(out.V, -1.0f, kEps); // right → left
}

// =============================================================================
// PolarisedEmission Tests
// =============================================================================

class PolarisedEmissionTests : public ::testing::Test {};

TEST_F(PolarisedEmissionTests, SynchrotronPolarisationDegree) {
    // p=2: (2+1)/(2+7/3) = 3/(13/3) = 9/13 ≈ 0.6923
    float pi_L = Sirius::PolarisedEmission::synchrotronPolarisationDegree(2.0f);
    EXPECT_NEAR(pi_L, 9.0f / 13.0f, kEps);

    // p=3: (3+1)/(3+7/3) = 4/(16/3) = 12/16 = 0.75
    pi_L = Sirius::PolarisedEmission::synchrotronPolarisationDegree(3.0f);
    EXPECT_NEAR(pi_L, 0.75f, kEps);
}

TEST_F(PolarisedEmissionTests, SynchrotronEmissionIsPhysical) {
    auto s = Sirius::PolarisedEmission::synchrotronEmission(1.0f, 0.7f, 0.0f);
    EXPECT_TRUE(s.isPhysical());
    EXPECT_NEAR(s.V, 0.0f, kEps); // no circular from synchrotron
}

TEST_F(PolarisedEmissionTests, ThomsonScatteringAt90Degrees) {
    // At 90 degrees (cos_theta = 0): pi = sin^2(90)/(1+cos^2(90)) = 1/1 = 1
    float pi = Sirius::PolarisedEmission::thomsonPolarisationDegree(0.0f);
    EXPECT_NEAR(pi, 1.0f, kEps);
}

TEST_F(PolarisedEmissionTests, ThomsonScatteringForward) {
    // At 0 degrees (cos_theta = 1): pi = 0/(1+1) = 0
    float pi = Sirius::PolarisedEmission::thomsonPolarisationDegree(1.0f);
    EXPECT_NEAR(pi, 0.0f, kEps);
}

// =============================================================================
// ParallelTransport Tests
// =============================================================================

class ParallelTransportTests : public ::testing::Test {};

TEST_F(ParallelTransportTests, ZeroSpinNoRotation) {
    float angle = Sirius::ParallelTransport::gravitationalFaradayRotation(0.0f, 10.0f, kPi / 2.0f);
    EXPECT_NEAR(angle, 0.0f, kEps);
}

TEST_F(ParallelTransportTests, RotationIncreasesWithSpin) {
    float a1 = Sirius::ParallelTransport::gravitationalFaradayRotation(0.5f, 10.0f, kPi / 2.0f);
    float a2 = Sirius::ParallelTransport::gravitationalFaradayRotation(0.9f, 10.0f, kPi / 2.0f);
    EXPECT_GT(std::abs(a2), std::abs(a1));
}

TEST_F(ParallelTransportTests, ApplyPreservesIntensity) {
    auto s = Sirius::StokesVector::horizontal(1.0f);
    auto out = Sirius::ParallelTransport::applyParallelTransport(s, 0.5f);
    EXPECT_NEAR(out.I, s.I, kEps); // rotation preserves I
}

} // namespace sirius::test
