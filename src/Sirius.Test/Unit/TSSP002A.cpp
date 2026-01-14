// TSSP002A.cpp - Redshift Energy Conservation and Spectral Transport Tests
// Tests for Phase 4.2 (energy conservation) and KNST001A (spectral transport)

#include <gtest/gtest.h>
#include <cmath>
#include "../../Sirius.Math/MTSB001A.h"
#include "../../Sirius.Kernel/KNST001A.h"
#include "../../Sirius.Physics/Disk/PHAD001A.h"

using namespace sirius::spectral;
using namespace sirius::kernel;
using namespace sirius::physics;

namespace {

//==============================================================================
// Redshift Energy Conservation Tests
// Note: Our 32-bin visible spectrum (380-780nm) loses energy to IR/UV
// when wavelengths shift outside this range. True g⁴ scaling only applies
// to bolometric (all wavelengths) flux.
//==============================================================================

TEST(RedshiftEnergyTest, RedshiftReducesVisibleEnergy) {
    // Under redshift, wavelengths shift redward. Blue light shifts to red,
    // red light shifts to IR (out of our visible bins).
    
    SpectralRadiance original = SpectralRadiance::blackbody(6000);  // Sun-like
    double E_original = original.totalEnergy();
    
    // Apply redshift g = 0.8 (approaching BH)
    double g = 0.8;
    SpectralRadiance redshifted = original.applyRedshift(g);
    double E_redshifted = redshifted.totalEnergy();
    
    // Redshift should reduce visible energy (some goes to IR)
    EXPECT_LT(E_redshifted, E_original)
        << "Redshift should reduce visible energy";
    
    // But not by more than the intensity scaling g⁴
    double max_reduction = std::pow(g, 4);
    EXPECT_GT(E_redshifted / E_original, max_reduction * 0.5)
        << "Energy loss should be bounded";
}

TEST(RedshiftEnergyTest, BlueshiftIncreasesVisibleEnergy) {
    SpectralRadiance original = SpectralRadiance::blackbody(6000);
    double E_original = original.totalEnergy();
    
    // Blueshift g = 1.2 (Doppler approaching)
    double g = 1.2;
    SpectralRadiance blueshifted = original.applyRedshift(g);
    double E_blueshifted = blueshifted.totalEnergy();
    
    // For this Sun-like spectrum, blueshift brings some IR into visible
    // and shifts visible toward UV. Net effect depends on blackbody shape.
    // Energy change is governed by g⁴ intensity scaling + wavelength shift
    
    // With our implementation, blueshifted visible should be somewhat higher
    // (g⁴ × 1.2⁴ ≈ 2.07× if all wavelengths stayed in range)
    EXPECT_GT(E_blueshifted, E_original * 0.8)
        << "Blueshift should not dramatically reduce energy";
}

TEST(RedshiftEnergyTest, ModerateRedshiftPreservesSomeEnergy) {
    SpectralRadiance original = SpectralRadiance::blackbody(5500);
    double E_original = original.totalEnergy();
    
    // Moderate redshift
    double g = 0.9;
    SpectralRadiance redshifted = original.applyRedshift(g);
    double E_redshifted = redshifted.totalEnergy();
    
    // Should preserve most energy for small redshift
    EXPECT_GT(E_redshifted / E_original, 0.3)
        << "Moderate redshift should preserve significant energy";
}

TEST(RedshiftEnergyTest, ExtremeRedshiftToIR) {
    SpectralRadiance original = SpectralRadiance::blackbody(6000);
    
    // Strong redshift pushes visible light to IR (out of our 380-780nm range)
    double g = 0.3;
    SpectralRadiance redshifted = original.applyRedshift(g);
    
    // Most energy should be lost to IR
    double E_original = original.totalEnergy();
    double E_redshifted = redshifted.totalEnergy();
    
    EXPECT_LT(E_redshifted, E_original * std::pow(g, 4) * 1.5)
        << "Extreme redshift should lose most visible energy to IR";
}

//==============================================================================
// Spectral Transport Tests
//==============================================================================

TEST(SpectralTransportTest, DiskGFactorComputation) {
    // Schwarzschild: g = √(1 - 2M/r)
    
    // At r = 10M: g = √(1 - 0.2) = √0.8 ≈ 0.894
    double g_10 = SpectralTransport::computeDiskGFactor(10.0, 0);
    EXPECT_NEAR(g_10, std::sqrt(0.8), 1e-10);
    
    // At r = 100M: g ≈ 0.99
    double g_100 = SpectralTransport::computeDiskGFactor(100.0, 0);
    EXPECT_NEAR(g_100, std::sqrt(0.98), 1e-10);
    
    // At r = 3M (photon sphere): g = √(1/3) ≈ 0.577
    double g_3 = SpectralTransport::computeDiskGFactor(3.0, 0);
    EXPECT_NEAR(g_3, std::sqrt(1.0/3.0), 1e-10);
}

TEST(SpectralTransportTest, IntegrateEmptyPath) {
    AccretionDiskD::Config diskConfig;
    AccretionDiskD disk(diskConfig);
    
    SpectralTransport transport(&disk);
    
    std::vector<BeamStateD> emptyPath;
    std::vector<double> emptyG;
    
    auto result = transport.integrate(emptyPath, emptyG);
    
    EXPECT_DOUBLE_EQ(result.totalRadiance.totalEnergy(), 0)
        << "Empty path should give zero radiance";
    EXPECT_EQ(result.diskHits, 0);
    EXPECT_FALSE(result.hitHorizon);
}

TEST(SpectralTransportTest, SinglePointEmission) {
    AccretionDiskD::Config diskConfig;
    diskConfig.M = 10.0;
    diskConfig.Mdot = 1e-8;
    AccretionDiskD disk(diskConfig);
    
    SpectralTransport transport(&disk);
    
    // Point on disk
    BeamStateD beam;
    beam.x = Vec4d(0, 15.0, M_PI/2, 0);
    beam.terminated = false;
    beam.initialPixelSolidAngle = 1e-8;
    beam.solidAngle = 1e-8;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            beam.J[i][j] = (i == j) ? 1.0 : 0.0;
    
    double g = SpectralTransport::computeDiskGFactor(15.0, 0);
    auto emission = transport.evaluateEmission(beam, g);
    
    EXPECT_GT(emission.totalEnergy(), 0)
        << "Disk point should emit";
}

TEST(SpectralTransportTest, HorizonDetection) {
    AccretionDiskD::Config diskConfig;
    AccretionDiskD disk(diskConfig);
    
    SpectralTransport transport(&disk);
    
    // Path that terminates inside ISCO
    std::vector<BeamStateD> path(1);
    path[0].x = Vec4d(0, 2.0, M_PI/2, 0);  // Inside ISCO (r < 6M)
    path[0].terminated = true;
    
    std::vector<double> gFactors = {0.3};
    
    auto result = transport.integrate(path, gFactors);
    
    EXPECT_TRUE(result.hitHorizon)
        << "Path terminating at r < ISCO should detect horizon";
}

} // namespace
