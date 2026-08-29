// Stokes vector and Mueller matrix tests. Ported from TSPL001A.cpp; assertions
// and tolerances unchanged. Tests stokes.h (StokesVector, MuellerMatrix,
// polarised_emission, parallel_transport).

#include "sirius/core/polarisation/stokes.h"

#include <gtest/gtest.h>

#include <cmath>
#include <limits>

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

TEST_F(PolarisedEmissionTests, ChandrasekharAtmosphereHasPhysicalEndpointPolarisation) {
    const auto edge = polarised_emission::ChandrasekharElectronScatteringAtmosphere(0.0f);
    const auto face = polarised_emission::ChandrasekharElectronScatteringAtmosphere(1.0f);
    ASSERT_TRUE(edge.has_value());
    ASSERT_TRUE(face.has_value());
    EXPECT_NEAR(edge->linear_polarisation_degree, 0.1171f, kEps);
    EXPECT_NEAR(face->linear_polarisation_degree, 0.0f, kEps);
    EXPECT_GT(edge->intensity_scale, 0.0f);
    EXPECT_GT(face->intensity_scale, edge->intensity_scale);
}

TEST_F(PolarisedEmissionTests, ChandrasekharAtmospherePreservesHemisphericFlux) {
    // The analytic integral of 2 mu (1+a mu)/(1+2a/3) from zero to one is one.
    constexpr double a = 2.06;
    const double integrated_flux = 2.0 * (0.5 + a / 3.0) / (1.0 + 2.0 * a / 3.0);
    EXPECT_NEAR(integrated_flux, 1.0, 2.0e-15);
}

TEST_F(PolarisedEmissionTests, ChandrasekharAtmosphereRejectsInvalidDirectionCosines) {
    EXPECT_FALSE(polarised_emission::ChandrasekharElectronScatteringAtmosphere(-0.01f).has_value());
    EXPECT_FALSE(polarised_emission::ChandrasekharElectronScatteringAtmosphere(1.01f).has_value());
    EXPECT_FALSE(polarised_emission::ChandrasekharElectronScatteringAtmosphere(
                     std::numeric_limits<float>::quiet_NaN())
                     .has_value());
}

// =============================================================================
// parallel_transport Tests
// =============================================================================

class ParallelTransportTests : public ::testing::Test {};

TEST_F(ParallelTransportTests, ApplyPreservesIntensity) {
    auto s = sirius::core::StokesVector::Horizontal(1.0f);
    auto out = sirius::core::parallel_transport::ApplyParallelTransport(s, 0.5f);
    EXPECT_NEAR(out.I, s.I, kEps);  // rotation preserves I
}

}  // namespace sirius::test
