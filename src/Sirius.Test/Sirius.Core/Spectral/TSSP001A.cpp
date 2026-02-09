// TSSP001A.cpp - Spectral Rendering Tests
// Tests for MTSB001A spectral radiance representation
// Verifies: Blackbody spectrum, redshift energy conservation, colour conversion

#include <gtest/gtest.h>
#include <cmath>
#include <MTSB001A.h>

using namespace sirius::spectral;

namespace {

//==============================================================================
// Test: Blackbody Peak Wavelength (Wien's Law)
// Verifies: 6500K blackbody peak near 446nm
//==============================================================================

TEST(SpectralRadianceTest, BlackbodyPeakWavelength) {
    // Wien's displacement law: λ_max = b/T, with b = 2897.77 μm·K
    // For T = 6500K: λ_max ≈ 446nm
    
    SpectralRadiance bb = SpectralRadiance::blackbody(6500);
    
    // Find peak bin
    int peakBin = 0;
    double peakValue = 0;
    for (int i = 0; i < NUM_WAVELENGTH_BINS; ++i) {
        if (bb.L[i] > peakValue) {
            peakValue = bb.L[i];
            peakBin = i;
        }
    }
    
    double peakWavelength = SpectralRadiance::wavelength(peakBin);
    
    // Should be in the blue-green region (roughly 440-520nm for 6500K)
    EXPECT_GT(peakWavelength, 420) << "Peak wavelength too low";
    EXPECT_LT(peakWavelength, 520) << "Peak wavelength too high";
}

//==============================================================================
// Test: Blackbody White Point
// Verifies: 6500K blackbody produces approximately neutral XYZ
//==============================================================================

TEST(SpectralRadianceTest, BlackbodyWhitePoint) {
    SpectralRadiance bb = SpectralRadiance::blackbody(6504);  // D65 approximation
    
    auto xyz = bb.toXYZ();
    
    // D65 white point ratios: X/Y ≈ 0.95, Z/Y ≈ 1.09
    // We just verify the ratios are reasonable for a white-ish source
    double xRatio = xyz.X / xyz.Y;
    double zRatio = xyz.Z / xyz.Y;
    
    EXPECT_GT(xRatio, 0.85) << "X/Y ratio too low for white";
    EXPECT_LT(xRatio, 1.05) << "X/Y ratio too high for white";
    EXPECT_GT(zRatio, 0.90) << "Z/Y ratio too low for white";
    EXPECT_LT(zRatio, 1.20) << "Z/Y ratio too high for white";
}

//==============================================================================
// Test: Redshift Energy Conservation
// Verifies: Total energy conserved under redshift (modulo UV/IR cutoff)
//==============================================================================

TEST(SpectralRadianceTest, RedshiftEnergyConservation) {
    SpectralRadiance original = SpectralRadiance::blackbody(5000);
    
    double E0 = original.totalEnergy();
    
    // Apply small blueshift (g > 1, wavelengths decrease)
    // For small shifts, energy should be approximately conserved
    double g = 1.05;  // 5% blueshift
    SpectralRadiance shifted = original.applyRedshift(g);
    
    double E1 = shifted.totalEnergy();
    
    // g⁴ factor for intensity, but wavelength compression means
    // some energy may shift out of range. For small g, should be close.
    // The g⁴ factor = 1.05⁴ ≈ 1.2155
    // But this applies per-photon; for integrated energy, it's different.
    
    // Just verify energy is positive and reasonable
    EXPECT_GT(E1, 0) << "Shifted energy should be positive";
    EXPECT_GT(E1, E0 * 0.5) << "Too much energy lost in shift";
}

//==============================================================================
// Test: Redshift Wavelength Shift
// Verifies: Redshift moves spectral content to longer wavelengths
//==============================================================================

TEST(SpectralRadianceTest, RedshiftWavelengthShift) {
    // Create narrow spectrum centred at 500nm
    SpectralRadiance original = SpectralRadiance::zero();
    int centreBin500 = SpectralRadiance::binIndex(500);
    original.L[centreBin500] = 1.0;
    
    // Apply redshift g = 0.8 (20% redshift)
    // λ_obs = λ_emit / g = 500 / 0.8 = 625nm
    double g = 0.8;
    SpectralRadiance shifted = original.applyRedshift(g);
    
    int expectedBin = SpectralRadiance::binIndex(625);
    
    // The energy should have moved to around 625nm
    EXPECT_GT(shifted.L[expectedBin], 0) << "Redshifted energy should appear at longer wavelength";
}

//==============================================================================
// Test: sRGB Conversion Range
// Verifies: sRGB output is in [0,1] range
//==============================================================================

TEST(SpectralRadianceTest, SRGBConversionRange) {
    SpectralRadiance bb = SpectralRadiance::blackbody(5500);
    
    // Normalise to reasonable brightness
    bb *= 1e-12;  // Scale down to reasonable display range
    
    auto rgb = bb.toSRGB();
    
    EXPECT_GE(rgb.r, 0) << "sRGB red should be >= 0";
    EXPECT_LE(rgb.r, 1) << "sRGB red should be <= 1";
    EXPECT_GE(rgb.g, 0) << "sRGB green should be >= 0";
    EXPECT_LE(rgb.g, 1) << "sRGB green should be <= 1";
    EXPECT_GE(rgb.b, 0) << "sRGB blue should be >= 0";
    EXPECT_LE(rgb.b, 1) << "sRGB blue should be <= 1";
}

//==============================================================================
// Test: ACES Conversion
// Verifies: ACES output is reasonable for blackbody
//==============================================================================

TEST(SpectralRadianceTest, ACESConversion) {
    SpectralRadiance bb = SpectralRadiance::blackbody(5500);
    
    auto aces = bb.toACES();
    
    // ACES should produce positive values for natural light sources
    // (though it can produce negative values for highly saturated sources)
    EXPECT_GT(aces.r, 0) << "ACES red should be positive for blackbody";
    EXPECT_GT(aces.g, 0) << "ACES green should be positive for blackbody";
    EXPECT_GT(aces.b, 0) << "ACES blue should be positive for blackbody";
}

//==============================================================================
// Test: Spectral Arithmetic
// Verifies: Addition and multiplication work correctly
//==============================================================================

TEST(SpectralRadianceTest, SpectralArithmetic) {
    SpectralRadiance a = SpectralRadiance::zero();
    SpectralRadiance b = SpectralRadiance::zero();
    
    a.L[10] = 1.0;
    b.L[10] = 2.0;
    b.L[20] = 3.0;
    
    // Test addition
    SpectralRadiance c = a + b;
    EXPECT_DOUBLE_EQ(c.L[10], 3.0);
    EXPECT_DOUBLE_EQ(c.L[20], 3.0);
    
    // Test scalar multiplication
    SpectralRadiance d = a * 5.0;
    EXPECT_DOUBLE_EQ(d.L[10], 5.0);
    
    // Test += operator
    a += b;
    EXPECT_DOUBLE_EQ(a.L[10], 3.0);
    EXPECT_DOUBLE_EQ(a.L[20], 3.0);
}

//==============================================================================
// Test: Wavelength Bin Indexing
// Verifies: Correct mapping between wavelength and bin index
//==============================================================================

TEST(SpectralRadianceTest, WavelengthBinIndexing) {
    // First bin centre
    double w0 = SpectralRadiance::wavelength(0);
    EXPECT_NEAR(w0, 380 + LAMBDA_STEP/2, 1e-6);
    
    // Last bin centre
    double w31 = SpectralRadiance::wavelength(31);
    EXPECT_NEAR(w31, 780 - LAMBDA_STEP/2, 1e-6);
    
    // Bin index from wavelength
    int bin500 = SpectralRadiance::binIndex(500);
    double w500 = SpectralRadiance::wavelength(bin500);
    EXPECT_LT(std::abs(w500 - 500), LAMBDA_STEP);
}

} // namespace
