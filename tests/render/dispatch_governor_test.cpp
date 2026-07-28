// Dispatch governor gate (dispatch_governor.h): band heights honour their
// [1, remaining] bound, growth is damped, overshoot shrinks proportionally,
// disabling restores one-dispatch-per-tile, and the environment override is
// loud on garbage. Pure-arithmetic suite; no Vulkan device is touched.

#include "sirius/render/dispatch_governor.h"

#include <gtest/gtest.h>

#include "support/scoped_environment.h"

namespace {

using sirius::render::BandController;
using sirius::render::kBandGrowthCap;
using sirius::render::kDefaultDispatchTargetMs;
using sirius::render::kInitialBandRows;
using sirius::render::ResolveDispatchTargetMs;
using sirius::test::ScopedEnvironmentVariable;

TEST(DispatchGovernor, FirstBandIsTheInitialHeightClampedToTheTile) {
    BandController wide(4096, 250.0);
    EXPECT_EQ(wide.NextRows(4096), kInitialBandRows);

    BandController narrow(4, 250.0);
    EXPECT_EQ(narrow.NextRows(4), 4);
}

TEST(DispatchGovernor, BandsNeverExceedRemainingRowsNorDropBelowOne) {
    BandController bands(1024, 250.0);
    // Fast measurements grow the band; the remaining-rows clamp still binds.
    for (int i = 0; i < 12; ++i) {
        bands.Record(1.0);
    }
    EXPECT_EQ(bands.NextRows(3), 3);
    EXPECT_EQ(bands.NextRows(1), 1);

    // Massive overshoot shrinks to the floor of one row, never zero.
    bands.Record(1e9);
    EXPECT_EQ(bands.NextRows(1024), 1);
}

TEST(DispatchGovernor, GrowthPerStepIsBoundedByTheCap) {
    BandController bands(4096, 250.0);
    const int before = bands.NextRows(4096);
    bands.Record(0.001);  // absurdly fast: uncapped growth would explode
    const int after = bands.NextRows(4096);
    EXPECT_LE(after, static_cast<int>(before * kBandGrowthCap));
    EXPECT_GT(after, before);
}

TEST(DispatchGovernor, ZeroMeasurementTakesTheCappedGrowthStep) {
    BandController bands(4096, 250.0);
    const int before = bands.NextRows(4096);
    bands.Record(0.0);  // below clock resolution: no rate information
    EXPECT_EQ(bands.NextRows(4096), static_cast<int>(before * kBandGrowthCap));
}

TEST(DispatchGovernor, OvershootShrinksProportionallyInOneStep) {
    BandController bands(4096, 250.0);
    bands.Record(1000.0);  // 4x the target: the next band must be ~4x shorter
    EXPECT_EQ(bands.NextRows(4096), kInitialBandRows / 4);
}

TEST(DispatchGovernor, DisabledControllerDispatchesWholeTilesAndIgnoresFeedback) {
    BandController bands(4096, 0.0);
    EXPECT_FALSE(bands.Enabled());
    EXPECT_EQ(bands.NextRows(2864), 2864);
    bands.Record(1e9);
    EXPECT_EQ(bands.NextRows(2864), 2864);
}

TEST(DispatchGovernor, TargetDefaultsWhenTheEnvironmentIsUnset) {
    ScopedEnvironmentVariable clear("SIRIUS_DISPATCH_TARGET_MS", nullptr);
    const auto target = ResolveDispatchTargetMs();
    ASSERT_TRUE(target.has_value());
    EXPECT_DOUBLE_EQ(*target, kDefaultDispatchTargetMs);
}

TEST(DispatchGovernor, TargetHonoursTheOverrideIncludingZero) {
    {
        ScopedEnvironmentVariable set("SIRIUS_DISPATCH_TARGET_MS", "125.5");
        const auto target = ResolveDispatchTargetMs();
        ASSERT_TRUE(target.has_value());
        EXPECT_DOUBLE_EQ(*target, 125.5);
    }
    {
        ScopedEnvironmentVariable zero("SIRIUS_DISPATCH_TARGET_MS", "0");
        const auto target = ResolveDispatchTargetMs();
        ASSERT_TRUE(target.has_value());
        EXPECT_DOUBLE_EQ(*target, 0.0);
        EXPECT_FALSE(BandController(64, *target).Enabled());
    }
}

TEST(DispatchGovernor, TargetFailsLoudOnGarbageAndNegatives) {
    for (const char* bad : {"fast", "-1", "250ms", "nan"}) {
        ScopedEnvironmentVariable set("SIRIUS_DISPATCH_TARGET_MS", bad);
        EXPECT_FALSE(ResolveDispatchTargetMs().has_value()) << bad;
    }
}

}  // namespace
