// Dispatch governor gate (dispatch_governor.h): band heights honour their
// [1, remaining] bound, growth is damped, overshoot shrinks proportionally,
// feedback is attributed to the work actually dispatched, learned area
// normalises across band widths, disabling restores one-dispatch-per-tile,
// and the environment override is loud on garbage. The two regression tests
// pin the confirmed findings of the 2026-07-28 adversarial review.
// Pure-arithmetic suite; no Vulkan device is touched.

#include "sirius/render/dispatch_governor.h"

#include <gtest/gtest.h>

#include "support/scoped_environment.h"

#include <cstdint>

namespace {

using sirius::render::BandController;
using sirius::render::kBandGrowthCap;
using sirius::render::kDefaultDispatchTargetMs;
using sirius::render::kInitialBandRows;
using sirius::render::ResolveDispatchTargetMs;
using sirius::test::ScopedEnvironmentVariable;

// Area of a full-width band `rows` high on a `width`-wide tile.
[[nodiscard]] constexpr std::int64_t Area(int rows, int width) {
    return static_cast<std::int64_t>(rows) * width;
}

TEST(DispatchGovernor, FirstBandIsTheInitialHeightClampedToTheTile) {
    BandController wide(4096, 250.0);
    EXPECT_EQ(wide.NextRows(4096, 4096), kInitialBandRows);

    BandController narrow(4, 250.0);
    EXPECT_EQ(narrow.NextRows(4, 4), 4);
}

TEST(DispatchGovernor, BandsNeverExceedRemainingRowsNorDropBelowOne) {
    BandController bands(1024, 250.0);
    // Fast measurements grow the band; the remaining-rows clamp still binds.
    for (int i = 0; i < 12; ++i) {
        const int rows = bands.NextRows(1024, 1024);
        bands.Record(Area(rows, 1024), 1.0);
    }
    EXPECT_EQ(bands.NextRows(3, 1024), 3);
    EXPECT_EQ(bands.NextRows(1, 1024), 1);

    // Massive overshoot shrinks to the floor, never zero.
    bands.Record(Area(1024, 1024), 1e9);
    EXPECT_EQ(bands.NextRows(1024, 1024), 1);
}

TEST(DispatchGovernor, GrowthPerStepIsBoundedByTheCap) {
    BandController bands(4096, 250.0);
    const int before = bands.NextRows(4096, 4096);
    bands.Record(Area(before, 4096), 0.001);  // absurdly fast
    const int after = bands.NextRows(4096, 4096);
    EXPECT_LE(after, static_cast<int>(before * kBandGrowthCap));
    EXPECT_GT(after, before);
}

TEST(DispatchGovernor, ZeroMeasurementTakesTheCappedGrowthStep) {
    BandController bands(4096, 250.0);
    const int before = bands.NextRows(4096, 4096);
    bands.Record(Area(before, 4096), 0.0);  // below clock resolution
    EXPECT_EQ(bands.NextRows(4096, 4096), static_cast<int>(before * kBandGrowthCap));
}

TEST(DispatchGovernor, OvershootShrinksProportionallyInOneStep) {
    BandController bands(4096, 250.0);
    bands.Record(Area(kInitialBandRows, 4096), 1000.0);  // 4x the target
    EXPECT_EQ(bands.NextRows(4096, 4096), kInitialBandRows / 4);
}

// Regression (adversarial review 2026-07-28, confirmed major): a tile's
// truncated tail band must feed back only the work it actually covered. The
// row-count controller read an 864-row tail measured at 108 ms as "the full
// 2000-row band ran fast", took the capped growth step to the whole tile
// edge, and submitted the next tile as a single whole-tile dispatch.
TEST(DispatchGovernor, TruncatedTailBandFeedsBackOnlyItsOwnWork) {
    constexpr int kEdge = 2864;
    BandController bands(kEdge, 250.0);
    // Steady state: 2000 rows on the clock at exactly the target.
    bands.Record(Area(2000, kEdge), 250.0);
    ASSERT_EQ(bands.NextRows(kEdge, kEdge), 2000);

    // Tile tail: only 864 rows remain; they run proportionally fast (108 ms).
    const int tail = bands.NextRows(864, kEdge);
    ASSERT_EQ(tail, 864);
    bands.Record(Area(tail, kEdge), 108.0);

    // The next band may grow from the TAIL's area (864 -> 1728 at the cap),
    // never from the desired area; a whole-tile band would be 358 ms.
    const int next = bands.NextRows(kEdge, kEdge);
    EXPECT_EQ(next, 1728);
    EXPECT_LT(next, kEdge);
}

// Regression (adversarial review 2026-07-28, confirmed critical): the learned
// state must be width-safe. A height learned on the narrow last-column tile
// of the tile grid costs proportionally more at full width; carrying the row
// count across produced a full-tile dispatch on the next full-width tile.
TEST(DispatchGovernor, LearnedAreaNormalisesAcrossBandWidths) {
    constexpr int kEdge = 2864;
    constexpr int kNarrow = 1232;  // 4096-wide image: second-column tile width
    BandController bands(kEdge, 250.0);
    // Learn on the narrow tile: 2000 rows at exactly the target.
    bands.Record(Area(2000, kNarrow), 250.0);
    ASSERT_EQ(bands.NextRows(kEdge, kNarrow), 2000);

    // Back at full width the same area is fewer rows, in proportion.
    EXPECT_EQ(bands.NextRows(kEdge, kEdge), (Area(2000, kNarrow)) / kEdge);
}

TEST(DispatchGovernor, DisabledControllerDispatchesWholeTilesAndIgnoresFeedback) {
    BandController bands(4096, 0.0);
    EXPECT_FALSE(bands.Enabled());
    EXPECT_EQ(bands.NextRows(2864, 2864), 2864);
    bands.Record(Area(2864, 2864), 1e9);
    EXPECT_EQ(bands.NextRows(2864, 2864), 2864);
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

// "inf" and overflowing literals parse as non-finite doubles; accepting them
// silently disabled banding while the log claimed the governor was active
// (adversarial review 2026-07-28, confirmed minor).
TEST(DispatchGovernor, TargetFailsLoudOnGarbageNegativesAndNonFinite) {
    for (const char* bad : {"fast", "-1", "250ms", "nan", "inf", "INFINITY", "1e400", "-1e400"}) {
        ScopedEnvironmentVariable set("SIRIUS_DISPATCH_TARGET_MS", bad);
        EXPECT_FALSE(ResolveDispatchTargetMs().has_value()) << bad;
    }
}

}  // namespace
