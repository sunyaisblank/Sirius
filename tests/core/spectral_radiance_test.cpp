// Spectral radiance tests: Planck blackbody spectrum, redshift energy behaviour,
// and CIE colour conversion. Ported from TSSP001A.cpp.

#include "sirius/core/spectral/spectral_radiance.h"

#include "sirius/core/spectral/blackbody.h"

#include <gtest/gtest.h>

#include <cmath>
#include <limits>

using namespace sirius::core;

namespace {

// Wien's law: a 6500K blackbody peaks in the blue-green.
TEST(SpectralRadianceTest, BlackbodyPeakWavelength) {
    SpectralRadiance bb = SpectralRadiance::Blackbody(6500);

    int peakBin = 0;
    double peakValue = 0;
    for (int i = 0; i < kNumWavelengthBins; ++i) {
        if (bb.L[i] > peakValue) {
            peakValue = bb.L[i];
            peakBin = i;
        }
    }

    double peakWavelength = SpectralRadiance::Wavelength(peakBin);

    EXPECT_GT(peakWavelength, 420) << "Peak wavelength too low";
    EXPECT_LT(peakWavelength, 520) << "Peak wavelength too high";
}

// A 6500K blackbody produces approximately neutral XYZ (D65 white point).
TEST(SpectralRadianceTest, BlackbodyWhitePoint) {
    SpectralRadiance bb = SpectralRadiance::Blackbody(6504);  // D65 approximation.

    auto xyz = bb.ToXyz();

    double xRatio = xyz.X / xyz.Y;
    double zRatio = xyz.Z / xyz.Y;

    EXPECT_GT(xRatio, 0.85) << "X/Y ratio too low for white";
    EXPECT_LT(xRatio, 1.05) << "X/Y ratio too high for white";
    EXPECT_GT(zRatio, 0.90) << "Z/Y ratio too low for white";
    EXPECT_LT(zRatio, 1.20) << "Z/Y ratio too high for white";
}

TEST(SpectralRadianceTest, BlackbodyBinsMatchPlanckAuthorityAndRejectInvalidTemperature) {
    const SpectralRadiance spectrum = SpectralRadiance::Blackbody(5778.0);
    for (int bin = 0; bin < kNumWavelengthBins; ++bin) {
        const double wavelength = SpectralRadiance::Wavelength(bin) * 1.0e-9;
        const double expected_per_nanometre = spectral::PlanckRadiance(wavelength, 5778.0) * 1.0e-9;
        EXPECT_LT(std::abs(spectrum.L[bin] - expected_per_nanometre) / spectrum.L[bin], 2.0e-14);
    }
    for (const double invalid : {0.0, -1.0, std::numeric_limits<double>::infinity(),
                                 std::numeric_limits<double>::quiet_NaN()}) {
        const SpectralRadiance rejected = SpectralRadiance::Blackbody(invalid);
        EXPECT_DOUBLE_EQ(rejected.TotalEnergy(), 0.0);
    }
}

TEST(SpectralRadianceTest, RedshiftDepositsExactlyOneGFourthWeightedBin) {
    SpectralRadiance original = SpectralRadiance::Zero();
    constexpr int kEmitBin = 9;
    constexpr int kObservedBin = 19;
    original.L[kEmitBin] = 3.25;
    const double g =
        SpectralRadiance::Wavelength(kEmitBin) / SpectralRadiance::Wavelength(kObservedBin);
    const double expected = original.L[kEmitBin] * std::pow(g, 4);
    const SpectralRadiance shifted = original.ApplyRedshift(g);

    EXPECT_DOUBLE_EQ(shifted.L[kObservedBin], expected);
    EXPECT_DOUBLE_EQ(shifted.TotalEnergy(), expected * kLambdaStep);
    EXPECT_DOUBLE_EQ(original.ApplyRedshift(0.0).TotalEnergy(), 0.0);
    EXPECT_DOUBLE_EQ(original.ApplyRedshift(-1.0).TotalEnergy(), 0.0);
    EXPECT_DOUBLE_EQ(original.ApplyRedshift(std::numeric_limits<double>::quiet_NaN()).TotalEnergy(),
                     0.0);
}

// Redshift moves spectral content to longer wavelengths.
TEST(SpectralRadianceTest, RedshiftWavelengthShift) {
    SpectralRadiance original = SpectralRadiance::Zero();
    int centreBin500 = SpectralRadiance::BinIndex(500);
    original.L[centreBin500] = 1.0;

    // g = 0.8 -> lambda_obs = 500 / 0.8 = 625 nm.
    double g = 0.8;
    SpectralRadiance shifted = original.ApplyRedshift(g);

    int expectedBin = SpectralRadiance::BinIndex(625);

    EXPECT_GT(shifted.L[expectedBin], 0) << "Redshifted energy should appear at longer wavelength";
}

// sRGB output is clamped to [0, 1].
TEST(SpectralRadianceTest, SRGBConversionRange) {
    SpectralRadiance bb = SpectralRadiance::Blackbody(5500);

    bb *= 1e-12;  // Scale down to a reasonable display range.

    auto rgb = bb.ToSrgb();

    EXPECT_GE(rgb.r, 0) << "sRGB red should be >= 0";
    EXPECT_LE(rgb.r, 1) << "sRGB red should be <= 1";
    EXPECT_GE(rgb.g, 0) << "sRGB green should be >= 0";
    EXPECT_LE(rgb.g, 1) << "sRGB green should be <= 1";
    EXPECT_GE(rgb.b, 0) << "sRGB blue should be >= 0";
    EXPECT_LE(rgb.b, 1) << "sRGB blue should be <= 1";
}

// ACES output is positive for a natural (blackbody) source.
TEST(SpectralRadianceTest, ACESConversion) {
    SpectralRadiance bb = SpectralRadiance::Blackbody(5500);

    auto aces = bb.ToAces();

    EXPECT_GT(aces.r, 0) << "ACES red should be positive for blackbody";
    EXPECT_GT(aces.g, 0) << "ACES green should be positive for blackbody";
    EXPECT_GT(aces.b, 0) << "ACES blue should be positive for blackbody";
}

// Addition and scalar multiplication behave componentwise.
TEST(SpectralRadianceTest, SpectralArithmetic) {
    SpectralRadiance a = SpectralRadiance::Zero();
    SpectralRadiance b = SpectralRadiance::Zero();

    a.L[10] = 1.0;
    b.L[10] = 2.0;
    b.L[20] = 3.0;

    SpectralRadiance c = a + b;
    EXPECT_DOUBLE_EQ(c.L[10], 3.0);
    EXPECT_DOUBLE_EQ(c.L[20], 3.0);

    SpectralRadiance d = a * 5.0;
    EXPECT_DOUBLE_EQ(d.L[10], 5.0);

    a += b;
    EXPECT_DOUBLE_EQ(a.L[10], 3.0);
    EXPECT_DOUBLE_EQ(a.L[20], 3.0);
}

// Wavelength and bin index map consistently.
TEST(SpectralRadianceTest, WavelengthBinIndexing) {
    double w0 = SpectralRadiance::Wavelength(0);
    EXPECT_NEAR(w0, 380 + kLambdaStep / 2, 1e-6);

    double w31 = SpectralRadiance::Wavelength(31);
    EXPECT_NEAR(w31, 780 - kLambdaStep / 2, 1e-6);

    int bin500 = SpectralRadiance::BinIndex(500);
    double w500 = SpectralRadiance::Wavelength(bin500);
    EXPECT_LT(std::abs(w500 - 500), kLambdaStep);
    EXPECT_EQ(SpectralRadiance::BinIndex(kLambdaMin - 0.01), -1);
    EXPECT_EQ(SpectralRadiance::BinIndex(kLambdaMin), 0);
    EXPECT_EQ(SpectralRadiance::BinIndex(kLambdaMax - 0.01), kNumWavelengthBins - 1);
    EXPECT_EQ(SpectralRadiance::BinIndex(kLambdaMax), kNumWavelengthBins);
}

}  // namespace
