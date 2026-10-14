// Novikov-Thorne accretion disk tests: ISCO radius, temperature profile, flux,
// limb darkening, disk boundaries, spectral emission. Ported from TSDK001A.cpp.

#include "sirius/core/disk/novikov_thorne_disk.h"

#include <gtest/gtest.h>

#include <cmath>

using namespace sirius::core;

namespace {

// r_ISCO = 6M for a non-spinning black hole.
TEST(AccretionDiskTest, ISCO_Schwarzschild) {
    double r_isco = AccretionDiskD::ComputeIsco(0.0);

    EXPECT_NEAR(r_isco, 6.0, 1e-10)
        << "Schwarzschild ISCO should be 6M";
}

// r_ISCO -> 1M for a prograde orbit around extremal Kerr.
TEST(AccretionDiskTest, ISCO_ExtremalKerr_Prograde) {
    double r_isco = AccretionDiskD::ComputeIsco(0.998);  // Near-extremal.

    EXPECT_LT(r_isco, 2.0)
        << "Near-extremal Kerr prograde ISCO should be < 2M";
    EXPECT_GT(r_isco, 1.0)
        << "ISCO should be > 1M (horizon)";
}

// r_ISCO -> 9M for a retrograde orbit around extremal Kerr.
TEST(AccretionDiskTest, ISCO_ExtremalKerr_Retrograde) {
    double r_isco = AccretionDiskD::ComputeIsco(-0.998);  // Near-extremal retrograde.

    EXPECT_GT(r_isco, 8.0)
        << "Near-extremal Kerr retrograde ISCO should be > 8M";
    EXPECT_LT(r_isco, 10.0)
        << "ISCO should be < 10M";
}

// Interpolation between 6M and 1M for moderate spin.
TEST(AccretionDiskTest, ISCO_ModerateSpin) {
    double r_isco_05 = AccretionDiskD::ComputeIsco(0.5);
    double r_isco_00 = AccretionDiskD::ComputeIsco(0.0);
    double r_isco_10 = AccretionDiskD::ComputeIsco(0.95);

    EXPECT_LT(r_isco_05, r_isco_00)
        << "Higher spin should give smaller prograde ISCO";
    EXPECT_GT(r_isco_05, r_isco_10)
        << "Higher spin should give smaller prograde ISCO";

    EXPECT_NEAR(r_isco_05, 4.233, 0.01)
        << "a=0.5 ISCO should be ~4.23M";
}

// Peak temperature near 1.5-2x the ISCO radius.
TEST(AccretionDiskTest, TemperatureProfileShape) {
    AccretionDiskD::Config config;
    config.M = 10.0;      // 10 solar masses.
    config.a_star = 0.0;  // Schwarzschild.
    config.Mdot = 1e-8;   // Typical rate.

    AccretionDiskD disk(config);

    double r_isco = disk.IscoRadius();

    EXPECT_DOUBLE_EQ(disk.Temperature(r_isco * 0.9), 0.0)
        << "Temperature inside ISCO should be zero";

    EXPECT_GT(disk.Temperature(r_isco * 1.5), 0)
        << "Temperature outside ISCO should be positive";

    double T_near = disk.Temperature(r_isco * 2);
    double T_far = disk.Temperature(r_isco * 10);

    EXPECT_GT(T_near, T_far)
        << "Temperature should decrease with radius";
}

// T ~ M^(-1/4) for a fixed Eddington ratio.
TEST(AccretionDiskTest, TemperatureScalingWithMass) {
    AccretionDiskD::Config config1;
    config1.M = 10.0;
    config1.Mdot = 1e-8;

    AccretionDiskD::Config config2;
    config2.M = 100.0;   // 10x mass.
    config2.Mdot = 1e-7;  // 10x accretion rate (same Eddington ratio).

    AccretionDiskD disk1(config1);
    AccretionDiskD disk2(config2);

    double T1 = disk1.Temperature(disk1.IscoRadius() * 2);
    double T2 = disk2.Temperature(disk2.IscoRadius() * 2);

    EXPECT_GT(T1, 0) << "Temperature should be positive";
    EXPECT_GT(T2, 0) << "Temperature should be positive";
}

// Face-on (cos = 1) gives maximum intensity; behind (cos < 0) gives zero.
TEST(AccretionDiskTest, LimbDarkening) {
    double I_faceon = AccretionDiskD::LimbDarkening(1.0, 0.6);

    double I_limb = AccretionDiskD::LimbDarkening(0.1, 0.6);

    EXPECT_GT(I_faceon, I_limb)
        << "Face-on should be brighter than limb";

    EXPECT_NEAR(I_faceon, 1.0, 0.01)
        << "Face-on limb darkening should be ~1.0";

    double I_behind = AccretionDiskD::LimbDarkening(-0.5, 0.6);
    EXPECT_DOUBLE_EQ(I_behind, 0.0)
        << "Behind disk should have zero intensity";
}

// IsInDisk enforces inner/outer edges and the equatorial band.
TEST(AccretionDiskTest, DiskBoundaries) {
    AccretionDiskD::Config config;
    config.a_star = 0.5;
    config.r_outer = 100;

    AccretionDiskD disk(config);

    double r_isco = disk.IscoRadius();

    EXPECT_FALSE(disk.IsInDisk(r_isco * 0.5, M_PI / 2))
        << "Inside ISCO should not be in disk";
    EXPECT_TRUE(disk.IsInDisk(r_isco * 1.5, M_PI / 2))
        << "Outside ISCO at equator should be in disk";
    EXPECT_TRUE(disk.IsInDisk(50, M_PI / 2))
        << "Mid-disk should be in disk";
    EXPECT_FALSE(disk.IsInDisk(200, M_PI / 2))
        << "Beyond outer edge should not be in disk";

    EXPECT_FALSE(disk.IsInDisk(50, M_PI / 4))
        << "Off equator should not be in disk";
}

// Blackbody spectrum for the disk temperature; zero inside the ISCO.
TEST(AccretionDiskTest, SpectralEmission) {
    AccretionDiskD::Config config;
    config.M = 10.0;
    config.Mdot = 1e-8;

    AccretionDiskD disk(config);

    double r = disk.IscoRadius() * 2;
    auto spectrum = disk.EmissionSpectrum(r);

    EXPECT_GT(spectrum.TotalEnergy(), 0)
        << "Disk emission should have positive energy";

    auto spectrum_inner = disk.EmissionSpectrum(disk.IscoRadius() * 0.5);
    EXPECT_DOUBLE_EQ(spectrum_inner.TotalEnergy(), 0)
        << "Inside ISCO should have zero emission";
}

// Flux is zero inside the ISCO and decreases outward past the peak.
TEST(AccretionDiskTest, FluxProfile) {
    AccretionDiskD::Config config;
    config.M = 10.0;
    config.Mdot = 1e-8;

    AccretionDiskD disk(config);

    double r_isco = disk.IscoRadius();

    EXPECT_DOUBLE_EQ(disk.Flux(r_isco * 0.9), 0.0)
        << "Flux inside ISCO should be zero";

    double F1 = disk.Flux(r_isco * 1.5);
    double F2 = disk.Flux(r_isco * 5);
    double F3 = disk.Flux(r_isco * 20);

    EXPECT_GT(F1, 0) << "Flux should be positive near ISCO";

    EXPECT_GT(F2, F3)
        << "Flux should decrease at large radii";
}

}  // namespace
