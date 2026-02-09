// TSPH011A.cpp - Spectral Rendering Utility Tests
// Tests: blackbody radiation, CIE color matching, XYZ→sRGB, redshift.

#include <gtest/gtest.h>
#include <cmath>

#include "PHSP001A.h"

namespace sirius::test {
using namespace Sirius;

using namespace Sirius::Spectral;
using namespace Sirius::Constants::Physical;

// =============================================================================
// Physical Constants Tests
// =============================================================================

// Test: Physical constants are positive and reasonable
TEST(SpectralUtilsTests, PhysicalConstantsValid) {
    EXPECT_GT(h_PLANCK, 0) << "Planck constant should be positive";
    EXPECT_GT(c_LIGHT, 0) << "Speed of light should be positive";
    EXPECT_GT(k_BOLTZMANN, 0) << "Boltzmann constant should be positive";

    // Order of magnitude checks
    EXPECT_GT(h_PLANCK, 1e-35);
    EXPECT_LT(h_PLANCK, 1e-33);
    EXPECT_GT(c_LIGHT, 2e8);
    EXPECT_LT(c_LIGHT, 4e8);
}

// =============================================================================
// Blackbody Radiation Tests
// =============================================================================

// Test: Planck radiance is positive for valid inputs
TEST(SpectralUtilsTests, PlanckRadiancePositive) {
    double lambda = 550e-9;  // Green light
    double T = 5778;  // Sun's surface temperature
    
    double radiance = planckRadiance(lambda, T);
    
    EXPECT_GT(radiance, 0) << "Blackbody radiance should be positive";
}

// Test: Planck radiance increases with temperature (at fixed wavelength)
TEST(SpectralUtilsTests, PlanckRadianceIncreasesWithTemp) {
    double lambda = 500e-9;
    
    double rad_3000K = planckRadiance(lambda, 3000);
    double rad_6000K = planckRadiance(lambda, 6000);
    
    EXPECT_GT(rad_6000K, rad_3000K) 
        << "Higher temperature should give higher radiance";
}

// Test: Planck radiance returns 0 for invalid inputs
TEST(SpectralUtilsTests, PlanckRadianceHandlesInvalid) {
    EXPECT_EQ(planckRadiance(0, 5000), 0);
    EXPECT_EQ(planckRadiance(-1e-9, 5000), 0);
    EXPECT_EQ(planckRadiance(500e-9, 0), 0);
    EXPECT_EQ(planckRadiance(500e-9, -100), 0);
}

// Test: Wien's law gives reasonable peak wavelengths
TEST(SpectralUtilsTests, WienPeakWavelengthReasonable) {
    // Sun (~5778K) peaks around 500nm
    double sun_peak = wienPeakWavelength(5778);
    EXPECT_GT(sun_peak, 400e-9);
    EXPECT_LT(sun_peak, 600e-9);
    
    // Hot disk (~10000K) peaks in UV
    double hot_peak = wienPeakWavelength(10000);
    EXPECT_LT(hot_peak, 400e-9);
    
    // Cool disk (~3000K) peaks in IR
    double cool_peak = wienPeakWavelength(3000);
    EXPECT_GT(cool_peak, 700e-9);
}

// =============================================================================
// CIE Color Matching Tests
// =============================================================================

// Test: Color matching functions are positive in visible range
TEST(SpectralUtilsTests, ColorMatchingPositiveInVisible) {
    // Sample wavelengths in visible range
    double wavelengths[] = {450, 500, 550, 600, 650};
    
    for (double lambda : wavelengths) {
        double y = cie_y(lambda);
        EXPECT_GE(y, 0) << "Y should be non-negative at " << lambda << "nm";
    }
}

// Test: Y peaks around 555nm (human eye sensitivity)
TEST(SpectralUtilsTests, YPeaksNearGreen) {
    double y_500 = cie_y(500);
    double y_555 = cie_y(555);
    double y_600 = cie_y(600);
    
    EXPECT_GT(y_555, y_500) << "Y at 555nm should exceed Y at 500nm";
    EXPECT_GT(y_555, y_600) << "Y at 555nm should exceed Y at 600nm";
}

// Test: Wavelength to XYZ gives zero outside visible range
TEST(SpectralUtilsTests, WavelengthToXYZOutOfRange) {
    XYZ uv = wavelengthToXYZ(300);  // UV
    XYZ ir = wavelengthToXYZ(900);  // IR
    
    EXPECT_EQ(uv.X, 0);
    EXPECT_EQ(uv.Y, 0);
    EXPECT_EQ(uv.Z, 0);
    EXPECT_EQ(ir.X, 0);
    EXPECT_EQ(ir.Y, 0);
    EXPECT_EQ(ir.Z, 0);
}

// =============================================================================
// Color Conversion Tests
// =============================================================================

// Test: XYZ to RGB conversion produces valid output
TEST(SpectralUtilsTests, XYZToRGBValid) {
    XYZ white(0.95047f, 1.0f, 1.08883f);  // D65 white point
    
    RGB rgb = xyzToLinearRGB(white);
    
    // Should be close to (1, 1, 1)
    EXPECT_NEAR(rgb.r, 1.0f, 0.1f);
    EXPECT_NEAR(rgb.g, 1.0f, 0.1f);
    EXPECT_NEAR(rgb.b, 1.0f, 0.1f);
}

// Test: sRGB gamma correction works correctly
TEST(SpectralUtilsTests, SRGBGammaCorrection) {
    // Linear 0 -> sRGB 0
    EXPECT_NEAR(srgbGamma(0.0f), 0.0f, 0.001f);
    
    // Linear 1 -> sRGB 1
    EXPECT_NEAR(srgbGamma(1.0f), 1.0f, 0.001f);
    
    // Linear 0.5 -> sRGB ~0.735 (brighter)
    float srgb_half = srgbGamma(0.5f);
    EXPECT_GT(srgb_half, 0.5f);
    EXPECT_LT(srgb_half, 0.8f);
}

// Test: Linear to sRGB clamps values
TEST(SpectralUtilsTests, LinearToSRGBClamps) {
    RGB over(1.5f, 2.0f, -0.5f);
    RGB clamped = linearToSRGB(over);
    
    EXPECT_LE(clamped.r, 1.0f);
    EXPECT_LE(clamped.g, 1.0f);
    EXPECT_GE(clamped.b, 0.0f);
}

// =============================================================================
// Blackbody to RGB Tests
// =============================================================================

// Test: Blackbody color is warmer at lower temperatures
TEST(SpectralUtilsTests, BlackbodyColorTemperature) {
    RGB cool = blackbodyToRGB(3000);  // Warm/red
    RGB hot = blackbodyToRGB(10000); // Cool/blue
    
    // Cooler temp should be more red
    EXPECT_GT(cool.r, cool.b);
    
    // Hotter temp should be more blue
    EXPECT_GT(hot.b, hot.r);
}

// Test: Sun temperature gives yellowish-white
TEST(SpectralUtilsTests, SunColorReasonable) {
    RGB sun = blackbodyToRGB(5778);
    
    // Sun should be fairly neutral (all channels close to 1)
    EXPECT_GT(sun.r, 0.5f);
    EXPECT_GT(sun.g, 0.5f);
    EXPECT_GT(sun.b, 0.5f);
}

// =============================================================================
// Redshift Tests
// =============================================================================

// Test: Zero velocity gives Doppler factor of 1
TEST(SpectralUtilsTests, DopplerFactorZeroVelocity) {
    double factor = dopplerFactor(0.0);
    EXPECT_NEAR(factor, 1.0, 1e-10);
}

// Test: Positive velocity (away) gives redshift
TEST(SpectralUtilsTests, DopplerFactorReceding) {
    double factor = dopplerFactor(0.5);  // v = 0.5c away
    EXPECT_LT(factor, 1.0) << "Receding source should have factor < 1";
}

// Test: Negative velocity (approach) gives blueshift
TEST(SpectralUtilsTests, DopplerFactorApproaching) {
    double factor = dopplerFactor(-0.5);  // v = 0.5c towards
    EXPECT_GT(factor, 1.0) << "Approaching source should have factor > 1";
}

// Test: Apply redshift shifts color towards red
TEST(SpectralUtilsTests, ApplyRedshiftEffect) {
    RGB white(1.0f, 1.0f, 1.0f);
    
    RGB redshifted = applyRedshift(white, 0.5);  // z = 0.5
    
    // Blue should be reduced more than red
    EXPECT_LT(redshifted.b, redshifted.r);
}

// Test: Apply blueshift shifts color towards blue
TEST(SpectralUtilsTests, ApplyBlueshiftEffect) {
    RGB white(1.0f, 1.0f, 1.0f);
    
    RGB blueshifted = applyRedshift(white, -0.3);  // z = -0.3
    
    // Red should be reduced more than blue
    EXPECT_LT(blueshifted.r, blueshifted.b);
}

// Test: Total redshift combines gravitational and Doppler
TEST(SpectralUtilsTests, TotalRedshiftCombined) {
    // Schwarzschild metric: g_tt = -(1 - rs/r)
    double g_tt_emit = -0.9;   // Closer to BH
    double g_tt_obs = -1.0;    // Far from BH
    double velocity = 0.1;     // Slightly receding
    
    double z = totalRedshift(g_tt_emit, g_tt_obs, velocity);
    
    // Should have net redshift (gravitational + Doppler both positive)
    EXPECT_GT(z, 0) << "Combined effect should give net redshift";
}

} // namespace sirius::test

// Note: main() provided by GTest::gtest_main
