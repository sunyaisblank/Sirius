// Dispatch governor gate (dispatch_governor.h): band heights honour their
// [1, remaining] bound, growth doubles and overshoot halves actual work,
// feedback is attributed to the work actually dispatched, learned area
// normalises across band widths, hard caps survive disabled adaptation,
// and the environment override is loud on garbage.
// Pure-arithmetic suite; no Vulkan device is touched.

#include "sirius/render/dispatch_governor.h"

#include "sirius/render/vulkan_renderer.h"

#include <gtest/gtest.h>

#include "support/scoped_environment.h"

#include <cstdint>
#include <limits>

namespace {

using sirius::render::BandController;
using sirius::render::kBandGrowthCap;
using sirius::render::kDefaultDispatchTargetMs;
using sirius::render::kHeavyFp32DispatchTargetMs;
using sirius::render::kInitialBandRows;
using sirius::render::kMaxTileEdge;
using sirius::render::kWatchdogSafeMaxBandRows;
using sirius::render::kWatchdogSafeMaxTileEdge;
using sirius::render::PrecisionRung;
using sirius::render::ResolveDispatchTargetMs;
using sirius::render::ResolveVulkanDispatchLimits;
using sirius::test::ScopedEnvironmentVariable;

// Area of a full-width band `rows` high on a `width`-wide tile.
[[nodiscard]] constexpr std::int64_t Area(int rows, int width) {
    return static_cast<std::int64_t>(rows) * width;
}

TEST(DispatchGovernor, FirstBandUsesTheMinimumFullWidthRowBeforeMeasurement) {
    static_assert(kInitialBandRows == 1);
    BandController wide(4096, 250.0);
    EXPECT_EQ(wide.NextRows(4096, 4096), kInitialBandRows);

    BandController narrow(4, 250.0);
    EXPECT_EQ(narrow.NextRows(4, 4), kInitialBandRows);

    BandController fp64_limited(64, 500.0, 1, 64);
    for (int i = 0; i < 8; ++i) fp64_limited.Record(64, 0.001);
    EXPECT_EQ(fp64_limited.NextRows(64, 64), 1);
}

TEST(DispatchGovernor, ExpensivePrecisionAndBundleWorkloadsUseTheStrictPhysicalFootprint) {
    const auto ordinary = ResolveVulkanDispatchLimits(PrecisionRung::Fp32, false, false);
    EXPECT_EQ(ordinary.tile_edge_cap, kMaxTileEdge);
    EXPECT_EQ(ordinary.max_band_rows, kMaxTileEdge);
    EXPECT_DOUBLE_EQ(ordinary.default_target_ms, 250.0);

    const auto ordinary_comp = ResolveVulkanDispatchLimits(PrecisionRung::Fp32Comp, false, false);
    EXPECT_EQ(ordinary_comp.tile_edge_cap, ordinary.tile_edge_cap);
    EXPECT_EQ(ordinary_comp.max_band_width, ordinary.max_band_width);
    EXPECT_EQ(ordinary_comp.max_band_rows, ordinary.max_band_rows);
    EXPECT_EQ(ordinary_comp.max_pixels, ordinary.max_pixels);
    EXPECT_DOUBLE_EQ(ordinary_comp.default_target_ms, 250.0);
    const auto ordinary_fp64 = ResolveVulkanDispatchLimits(PrecisionRung::Fp64, false, false);
    EXPECT_EQ(ordinary_fp64.max_band_rows, 1);
    EXPECT_EQ(ordinary_fp64.max_pixels, 64);
    EXPECT_DOUBLE_EQ(ordinary_fp64.default_target_ms, 250.0);

    for (const auto precision : {PrecisionRung::Fp64, PrecisionRung::Fp32Comp}) {
        for (const auto limits : {
                 ResolveVulkanDispatchLimits(precision, true, false),
                 ResolveVulkanDispatchLimits(precision, false, true),
                 ResolveVulkanDispatchLimits(precision, true, true),
             }) {
            EXPECT_EQ(limits.tile_edge_cap, kWatchdogSafeMaxTileEdge);
            EXPECT_EQ(limits.max_band_width, 64);
            EXPECT_EQ(limits.max_band_rows, 1);
            EXPECT_EQ(limits.max_pixels, 64);
            EXPECT_DOUBLE_EQ(limits.default_target_ms, 250.0);
        }
    }
    for (const auto limits : {
             ResolveVulkanDispatchLimits(PrecisionRung::Fp32, true, false),
             ResolveVulkanDispatchLimits(PrecisionRung::Fp32, false, true),
             ResolveVulkanDispatchLimits(PrecisionRung::Fp32, true, true),
         }) {
        EXPECT_EQ(limits.tile_edge_cap, kWatchdogSafeMaxTileEdge);
        EXPECT_EQ(limits.max_band_width, 64);
        EXPECT_EQ(limits.max_band_rows, kWatchdogSafeMaxBandRows);
        EXPECT_EQ(limits.max_pixels, 256);
        EXPECT_DOUBLE_EQ(limits.default_target_ms, 750.0);
    }
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

TEST(DispatchGovernor, ZeroMeasurementFallsBackToMinimumWork) {
    BandController bands(4096, 250.0);
    const int before = bands.NextRows(4096, 4096);
    bands.Record(Area(before, 4096), 0.0);  // no trustworthy throughput estimate
    EXPECT_EQ(bands.NextRows(4096, 4096), 1);
}

TEST(DispatchGovernor, OvershootHalvesTheObservedWork) {
    BandController bands(4096, 250.0);
    ASSERT_TRUE(bands.Record(Area(8, 4096), 500.0));
    EXPECT_EQ(bands.NextRows(4096, 4096), 4);
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

TEST(DispatchGovernor, DisabledAdaptationPreservesCapsAndIgnoresOrdinaryFeedback) {
    BandController bands(4096, 0.0);
    EXPECT_FALSE(bands.Enabled());
    EXPECT_EQ(bands.NextRows(2864, 2864), 2864);
    bands.Record(Area(2864, 2864), 500.0);
    EXPECT_EQ(bands.NextRows(2864, 2864), 2864);
    EXPECT_TRUE(bands.Record(Area(2864, 2864), 1100.0));
    EXPECT_EQ(bands.NextRows(2864, 2864), 1)
        << "safety fallback remains active when adaptation is disabled";
}

TEST(DispatchGovernor, StrictCapsApplyToEveryTailEvenWithoutAdaptation) {
    for (const double target : {0.0, kDefaultDispatchTargetMs, kHeavyFp32DispatchTargetMs}) {
        for (int max_rows : {1, 4}) {
            for (int width = 1; width <= 64; ++width) {
                for (int remaining = 1; remaining <= 64; ++remaining) {
                    BandController bands(width, target, max_rows, 64 * max_rows);
                    for (int observation = 0; observation < 8; ++observation) {
                        const int rows = bands.NextRows(remaining, width);
                        ASSERT_GE(rows, 1);
                        ASSERT_LE(rows, remaining);
                        ASSERT_LE(rows, max_rows);
                        ASSERT_LE(rows * width, 64 * max_rows);
                        ASSERT_TRUE(bands.Record(rows * width, 0.001));
                    }
                }
            }
        }
    }
    // A row limit remains independent of an area cap on narrow tails.
    BandController one_row(64, 0.0, 1);
    EXPECT_EQ(one_row.NextRows(64, 1), 1);
}

TEST(DispatchGovernor, NearTargetLatencyStillGrowsWholeBands) {
    BandController bands(64, kHeavyFp32DispatchTargetMs, 4, 256);
    EXPECT_EQ(bands.NextRows(64, 64), 1);
    ASSERT_TRUE(bands.Record(64, 700.0));
    EXPECT_EQ(bands.NextRows(64, 64), 2)
        << "fractional proportional growth would remain trapped at one row";
    ASSERT_TRUE(bands.Record(128, 720.0));
    EXPECT_EQ(bands.NextRows(64, 64), 4);
    ASSERT_TRUE(bands.Record(256, 800.0));
    EXPECT_EQ(bands.NextRows(64, 64), 2);
    ASSERT_TRUE(bands.Record(128, 900.0));
    EXPECT_EQ(bands.NextRows(64, 64), 1);
}

TEST(DispatchGovernor, InvalidTimingAndIrreducibleOvershootDecline) {
    for (double bad : {-1.0, std::numeric_limits<double>::infinity(),
                       std::numeric_limits<double>::quiet_NaN()}) {
        BandController bands(8, 250.0, 8, 64);
        EXPECT_FALSE(bands.Record(8, bad));
        EXPECT_EQ(bands.NextRows(64, 8), 1);
    }
    for (double target : {0.0, kDefaultDispatchTargetMs, kHeavyFp32DispatchTargetMs}) {
        BandController bands(8, target, 8, 64);
        EXPECT_TRUE(bands.Record(64, 1100.0));
        EXPECT_EQ(bands.NextRows(64, 8), 1)
            << "an excessive observation forces minimum work even with adaptation disabled";
        EXPECT_FALSE(bands.Record(8, 1100.0));
    }
}

TEST(DispatchGovernor, TargetDefaultsWhenTheEnvironmentIsUnset) {
    for (const char* unset : {static_cast<const char*>(nullptr), ""}) {
        ScopedEnvironmentVariable clear("SIRIUS_DISPATCH_TARGET_MS", unset);
        const auto target = ResolveDispatchTargetMs();
        ASSERT_TRUE(target.has_value());
        EXPECT_DOUBLE_EQ(*target, kDefaultDispatchTargetMs);
        const auto heavy_target = ResolveDispatchTargetMs(kHeavyFp32DispatchTargetMs);
        ASSERT_TRUE(heavy_target.has_value());
        EXPECT_DOUBLE_EQ(*heavy_target, kHeavyFp32DispatchTargetMs);
    }
}

TEST(DispatchGovernor, TargetHonoursTheOverrideIncludingZero) {
    for (double profile_default : {kDefaultDispatchTargetMs, kHeavyFp32DispatchTargetMs}) {
        ScopedEnvironmentVariable set("SIRIUS_DISPATCH_TARGET_MS", "125.5");
        const auto target = ResolveDispatchTargetMs(profile_default);
        ASSERT_TRUE(target.has_value());
        EXPECT_DOUBLE_EQ(*target, 125.5);
    }
    for (double profile_default : {kDefaultDispatchTargetMs, kHeavyFp32DispatchTargetMs}) {
        ScopedEnvironmentVariable zero("SIRIUS_DISPATCH_TARGET_MS", "0");
        const auto target = ResolveDispatchTargetMs(profile_default);
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
        EXPECT_FALSE(ResolveDispatchTargetMs(kHeavyFp32DispatchTargetMs).has_value()) << bad;
    }
}

}  // namespace
