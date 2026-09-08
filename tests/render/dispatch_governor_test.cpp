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

#include <array>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

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
    EXPECT_EQ(ordinary.max_band_width, 512);
    EXPECT_EQ(ordinary.max_band_rows, 16);
    EXPECT_EQ(ordinary.max_pixels, 8192);
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

TEST(DispatchGovernor, OrdinaryPrecisionBoundsSurviveGrowthAndDisabledAdaptation) {
    for (const auto precision : {PrecisionRung::Fp32, PrecisionRung::Fp32Comp}) {
        const auto limits = ResolveVulkanDispatchLimits(precision, false, false);
        // A larger memory tile must not enlarge either an initial submission or
        // a band grown after cheap sky rays into a more expensive disk region.
        for (const int tile_width : {1, 17, 512, 513, 4096}) {
            for (const double target : {0.0, 250.0, 750.0}) {
                for (int x = 0; x < tile_width; x += limits.max_band_width) {
                    const int width = std::min(limits.max_band_width, tile_width - x);
                    EXPECT_LE(width, 512);
                    BandController bands(width, target, limits.max_band_rows, limits.max_pixels);
                    for (int observation = 0; observation < 20; ++observation) {
                        const int rows = bands.NextRows(4096, width);
                        EXPECT_GE(rows, 1);
                        EXPECT_LE(rows, 16);
                        EXPECT_LE(Area(rows, width), 8192);
                        ASSERT_TRUE(bands.Record(Area(rows, width), 0.001));
                    }
                    EXPECT_EQ(bands.NextRows(3, width), 3);
                }
            }
        }
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
        EXPECT_TRUE(bands.Record(8, 1100.0)) << "a full row is still reducible";
        EXPECT_FALSE(bands.Record(1, 1100.0));
    }
}

TEST(DispatchGovernor, SafetyFallbackRestartsAllSamplesWithExactPixelCoverage) {
    using sirius::base::Expected;
    using sirius::render::CameraSample;
    using sirius::render::DispatchRegion;
    using sirius::render::ExecuteDispatchRegions;
    std::vector<CameraSample> samples;
    sirius::render::ForEachCameraSample(3, [&](const auto& sample) { samples.push_back(sample); });
    ASSERT_EQ(samples.size(), 3U);

    for (const auto shape : {std::array{1, 3}, std::array{3, 1}, std::array{7, 3}, std::array{9, 5},
                             std::array{17, 11}, std::array{64, 1}}) {
        for (const double target : {0.0, 250.0, 750.0}) {
            for (const double fallback_ms : {0.0, 1013.575332}) {
                for (const int late_sample : {1, 2}) {
                    SCOPED_TRACE(::testing::Message()
                                 << shape[0] << 'x' << shape[1] << " target=" << target
                                 << " timing=" << fallback_ms << " sample=" << late_sample);
                    const DispatchRegion original{5, 7, shape[0], shape[1]};
                    BandController bands(original.width, target, original.height,
                                         original.Pixels());
                    const auto size = static_cast<std::size_t>(original.Pixels());
                    std::vector<std::array<float, 4>> packed(size, {-999, -999, -999, -999});
                    auto frame = packed;
                    std::vector<int> coverage(size, 0);
                    std::vector<std::int64_t> expected_areas;
                    if (fallback_ms == 0.0 || size == 3) {
                        expected_areas.assign(size, 1);
                    } else if (shape[0] == 7) {
                        expected_areas = {7, 7, 7};
                    } else if (shape[0] == 9) {
                        expected_areas = {18, 9, 18};
                    } else if (shape[0] == 17) {
                        expected_areas = {85, 51, 51};
                    } else {
                        expected_areas = {32, 32};
                    }
                    const auto expected_cap = fallback_ms == 0.0 ? 1 : original.Pixels() / 2;
                    std::vector<std::int64_t> completed_areas;
                    bool injected = false;
                    int attempts = 0;
                    auto submitted = [&](const DispatchRegion& region, const CameraSample& sample,
                                         int index) -> Expected<double> {
                        ++attempts;
                        EXPECT_GE(index, 0);
                        EXPECT_LT(index, 3);
                        if (index < 0 || index >= 3) {
                            return sirius::base::Fail(sirius::base::ErrorDomain::kInternal,
                                                      "test dispatch", "invalid sample index");
                        }
                        EXPECT_EQ(sample.image_u, samples[index].image_u);
                        EXPECT_EQ(sample.image_v, samples[index].image_v);
                        EXPECT_EQ(sample.pupil_u, samples[index].pupil_u);
                        EXPECT_EQ(sample.pupil_v, samples[index].pupil_v);
                        if (injected) {
                            EXPECT_LE(region.Pixels(), expected_cap);
                            EXPECT_EQ(index, (attempts - late_sample - 2) % 3);
                        }
                        for (int row = 0; row < region.height; ++row) {
                            for (int column = 0; column < region.width; ++column) {
                                const int absolute = (region.y + row) * 128 + region.x + column;
                                auto& pixel =
                                    packed[static_cast<std::size_t>(row * region.width + column)];
                                for (int channel = 0; channel < 4; ++channel) {
                                    const float value =
                                        static_cast<float>(absolute * 32 + channel * 4 + index);
                                    pixel[channel] =
                                        index == 0 ? value
                                                   : (pixel[channel] * index + value) / (index + 1);
                                }
                            }
                        }
                        if (!injected && index == late_sample) {
                            injected = true;
                            return fallback_ms;
                        }
                        return 1.0;
                    };
                    auto completed = [&](const DispatchRegion& region) -> Expected<void> {
                        EXPECT_TRUE(injected);
                        completed_areas.push_back(region.Pixels());
                        EXPECT_LE(region.Pixels(), expected_cap);
                        EXPECT_GE(region.x, original.x);
                        EXPECT_GE(region.y, original.y);
                        EXPECT_LE(region.x + region.width, original.x + original.width);
                        EXPECT_LE(region.y + region.height, original.y + original.height);
                        if (region.x < original.x || region.y < original.y ||
                            region.x + region.width > original.x + original.width ||
                            region.y + region.height > original.y + original.height) {
                            return sirius::base::Fail(sirius::base::ErrorDomain::kInternal,
                                                      "test completion", "outside region");
                        }
                        for (int row = 0; row < region.height; ++row) {
                            for (int column = 0; column < region.width; ++column) {
                                const auto destination = static_cast<std::size_t>(
                                    (region.y + row - original.y) * original.width + region.x +
                                    column - original.x);
                                ++coverage[destination];
                                frame[destination] =
                                    packed[static_cast<std::size_t>(row * region.width + column)];
                            }
                        }
                        return {};
                    };
                    const auto result =
                        ExecuteDispatchRegions(original, 3, bands, submitted, completed);
                    ASSERT_TRUE(result.has_value()) << result.error().Description();
                    EXPECT_EQ(attempts, late_sample + 1 + expected_areas.size() * 3);
                    EXPECT_EQ(completed_areas, expected_areas);
                    EXPECT_TRUE(bands.SafetyFallback());
                    EXPECT_EQ(bands.SafetyPixelCap(), expected_cap);
                    for (int row = 0; row < original.height; ++row) {
                        for (int column = 0; column < original.width; ++column) {
                            const auto position =
                                static_cast<std::size_t>(row * original.width + column);
                            EXPECT_EQ(coverage[position], 1);
                            const int absolute = (original.y + row) * 128 + original.x + column;
                            for (int channel = 0; channel < 4; ++channel) {
                                // The independent mean of sample ordinals {0,1,2} is exactly 1.
                                EXPECT_EQ(frame[position][channel],
                                          static_cast<float>(absolute * 32 + channel * 4 + 1));
                            }
                        }
                    }
                }
            }
        }
    }
}

TEST(DispatchGovernor, MeasuredSafetyCapRetainsUsefulWorkAndBoundsRepeatedOvershoots) {
    using sirius::base::Expected;
    using sirius::render::CameraSample;
    using sirius::render::DispatchRegion;
    using sirius::render::ExecuteDispatchRegions;
    for (const double target : {0.0, 250.0, 750.0}) {
        for (const bool repeated : {false, true}) {
            SCOPED_TRACE(::testing::Message() << target << " repeated=" << repeated);
            BandController bands(64, target, 4, 256);
            std::vector<std::array<int, 4>> actual;
            std::vector<std::array<int, 3>> completed;
            auto submit = [&](const DispatchRegion& region, const CameraSample&,
                              int index) -> Expected<double> {
                actual.push_back({region.x, region.y, region.width, index});
                EXPECT_EQ(region.height, 1);
                EXPECT_LE(region.Pixels(), bands.SafetyPixelCap());
                if (region.width == 64 && index == 1) return 1013.575332;
                if (repeated && region.width == 32 && index == 2) return 1000.01;
                return 1.0;
            };
            auto complete = [&](const DispatchRegion& region) -> Expected<void> {
                completed.push_back({region.x, region.y, region.width});
                return {};
            };
            const auto first = ExecuteDispatchRegions({2, 3, 64, 1}, 3, bands, submit, complete);
            ASSERT_TRUE(first.has_value()) << first.error().Description();
            EXPECT_EQ(actual.size(), repeated ? 17U : 8U);
            EXPECT_EQ(bands.SafetyPixelCap(), repeated ? 16 : 32);
            // A new logical row must retain the measured cap after fast samples.
            const auto next = ExecuteDispatchRegions({2, 4, 64, 1}, 3, bands, submit, complete);
            ASSERT_TRUE(next.has_value()) << next.error().Description();
            EXPECT_EQ(bands.SafetyPixelCap(), repeated ? 16 : 32);
            EXPECT_EQ(bands.NextRows(64, 64), 1);
            std::vector<std::array<int, 4>> expected;
            auto samples = [&](int x, int y, int width, int count) {
                for (int index = 0; index < count; ++index) {
                    expected.push_back({x, y, width, index});
                }
            };
            samples(2, 3, 64, 2);
            if (repeated) samples(2, 3, 32, 3);
            std::vector<std::array<int, 3>> expected_completed;
            const int width = repeated ? 16 : 32;
            for (const int y : {3, 4}) {
                for (int x = 2; x < 66; x += width) {
                    samples(x, y, width, 3);
                    expected_completed.push_back({x, y, width});
                }
            }
            EXPECT_EQ(actual, expected);
            EXPECT_EQ(completed, expected_completed);
        }
        BandController bands(64, target, 4, 256);
        std::vector<std::int64_t> areas;
        int completions = 0;
        const auto refused = ExecuteDispatchRegions(
            {2, 3, 64, 1}, 3, bands,
            [&](const DispatchRegion& region, const CameraSample&, int index) -> Expected<double> {
                areas.push_back(region.Pixels());
                EXPECT_EQ(index, 0);
                return 1000.01;
            },
            [&](const DispatchRegion&) -> Expected<void> {
                ++completions;
                return {};
            });
        ASSERT_FALSE(refused.has_value());
        EXPECT_EQ(areas, (std::vector<std::int64_t>{64, 32, 16, 8, 4, 2, 1}));
        EXPECT_EQ(completions, 0);
        EXPECT_EQ(bands.SafetyPixelCap(), 1);
        BandController boundary(1, target, 1, 1);
        EXPECT_TRUE(boundary.Record(1, 1000.0));
        EXPECT_FALSE(boundary.SafetyFallback());
    }
}

TEST(DispatchGovernor, OnePixelSafetyBoundaryRemainsFatalWithStickyFallback) {
    using sirius::base::Expected;
    using sirius::render::CameraSample;
    using sirius::render::DispatchRegion;
    using sirius::render::ExecuteDispatchRegions;
    for (const double target : {0.0, 250.0, 750.0}) {
        for (const double bad : {1000.01, -1.0, std::numeric_limits<double>::infinity(),
                                 std::numeric_limits<double>::quiet_NaN()}) {
            BandController bands(64, target, 4, 256);
            int attempts = 0;
            int completions = 0;
            const auto result = ExecuteDispatchRegions(
                {2, 3, 64, 1}, 3, bands,
                [&](const DispatchRegion& region, const CameraSample&,
                    int index) -> Expected<double> {
                    ++attempts;
                    EXPECT_EQ(index, 0);
                    if (attempts == 1) return 0.0;
                    EXPECT_EQ(region.Pixels(), 1);
                    return bad;
                },
                [&](const DispatchRegion&) -> Expected<void> {
                    ++completions;
                    return {};
                });
            EXPECT_FALSE(result.has_value());
            EXPECT_EQ(attempts, 2);
            EXPECT_EQ(completions, 0);
        }
        BandController bands(1, target, 1, 1);
        int attempts = 0;
        int completions = 0;
        const auto zero = ExecuteDispatchRegions(
            {0, 0, 1, 1}, 3, bands,
            [&](const DispatchRegion&, const CameraSample&, int) -> Expected<double> {
                ++attempts;
                return 0.0;
            },
            [&](const DispatchRegion&) -> Expected<void> {
                ++completions;
                return {};
            });
        EXPECT_TRUE(zero.has_value()) << "zero timing at minimum work must not retry forever";
        EXPECT_EQ(attempts, 3);
        EXPECT_EQ(completions, 1);
    }
}

TEST(DispatchGovernor, RegionExecutionPropagatesCancellationAndBoundaryErrors) {
    using sirius::base::Expected;
    using sirius::render::CameraSample;
    using sirius::render::DispatchRegion;
    using sirius::render::ExecuteDispatchRegions;
    for (int failure_stage : {0, 1, 2}) {
        SCOPED_TRACE(failure_stage);
        BandController bands(9, 750.0, 3, 27);
        int attempts = 0;
        int completions = 0;
        const auto result = ExecuteDispatchRegions(
            {5, 7, 9, 3}, 3, bands,
            [&](const DispatchRegion&, const CameraSample&, int) -> Expected<double> {
                ++attempts;
                if (failure_stage == 1 && attempts == 2) {
                    return sirius::base::Fail(sirius::base::ErrorDomain::kDevice, "test submit",
                                              "submission failed");
                }
                return 1.0;
            },
            [&](const DispatchRegion&) -> Expected<void> {
                ++completions;
                return sirius::base::Fail(sirius::base::ErrorDomain::kDevice, "test readback",
                                          "readback failed");
            },
            [&] { return failure_stage == 0 && attempts == 1; });
        ASSERT_FALSE(result.has_value());
        EXPECT_EQ(attempts, failure_stage + 1);
        EXPECT_EQ(completions, failure_stage == 2 ? 1 : 0);
        const std::string expected = failure_stage == 0   ? "cancelled"
                                     : failure_stage == 1 ? "submission failed"
                                                          : "readback failed";
        EXPECT_NE(result.error().Description().find(expected), std::string::npos);
    }
    for (int failure_stage : {0, 1, 2}) {
        SCOPED_TRACE(::testing::Message() << "pending siblings, stage " << failure_stage);
        BandController bands(3, 750.0, 3, 9);
        int attempts = 0;
        int completions = 0;
        const auto result = ExecuteDispatchRegions(
            {5, 7, 3, 3}, 3, bands,
            [&](const DispatchRegion& region, const CameraSample&, int) -> Expected<double> {
                ++attempts;
                if (attempts == 1) return 1013.575332;
                EXPECT_EQ(region.Pixels(), 3);
                if (failure_stage == 1 && completions == 1) {
                    return sirius::base::Fail(sirius::base::ErrorDomain::kDevice, "test submit",
                                              "pending submission failed");
                }
                return 1.0;
            },
            [&](const DispatchRegion&) -> Expected<void> {
                ++completions;
                if (failure_stage == 2 && completions == 2) {
                    return sirius::base::Fail(sirius::base::ErrorDomain::kDevice, "test readback",
                                              "pending readback failed");
                }
                return {};
            },
            [&] { return failure_stage == 0 && completions == 1; });
        ASSERT_FALSE(result.has_value());
        EXPECT_EQ(attempts, (std::array{4, 5, 7})[failure_stage]);
        EXPECT_EQ(completions, failure_stage == 2 ? 2 : 1);
        const std::string expected = failure_stage == 0   ? "cancelled"
                                     : failure_stage == 1 ? "pending submission failed"
                                                          : "pending readback failed";
        EXPECT_NE(result.error().Description().find(expected), std::string::npos);
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
