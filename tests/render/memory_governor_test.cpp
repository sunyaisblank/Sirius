// Memory-governor unit gates (specification programme 4, docs/ARCHITECTURE.md
// section 6). Pure arithmetic over a reported budget, so these run with no
// Vulkan device: they pin that the governor derives a sane tile across budgets,
// shrinks the tile as the budget shrinks, never over-commits, and declines
// loudly when the budget cannot seat a minimal tile.

#include "sirius/render/memory_governor.h"

#include <gtest/gtest.h>

#include "support/scoped_environment.h"

#include <cstdint>

namespace {

using sirius::render::DeriveTilePlan;
using sirius::render::kMaxTileEdge;
using sirius::render::kMinTileEdge;
using sirius::render::kResidencyFraction;
using sirius::render::kTileWorkingSetBytesPerPixel;
using sirius::render::ResolveBudgetBytes;
using sirius::test::ScopedEnvironmentVariable;

constexpr std::uint64_t kMiB = 1024ULL * 1024ULL;
constexpr std::uint64_t kGiB = 1024ULL * kMiB;

// An IMAX-class target frame, large enough that the tile is budget-bound rather
// than image-bound, so budget differences show up as tile-size differences.
constexpr int kImaxWidth = 5616;
constexpr int kImaxHeight = 4096;

TEST(MemoryGovernor, TwoGigabyteBudgetSeatsAWorkableTile) {
    const auto plan = DeriveTilePlan(2 * kGiB, kImaxWidth, kImaxHeight, 0);
    ASSERT_TRUE(plan.has_value()) << plan.error().Description();
    EXPECT_GE(plan->tile_edge, kMinTileEdge);
    EXPECT_LE(plan->tile_edge, kMaxTileEdge);
    // The working set plus overhead never exceeds the usable fraction.
    EXPECT_LE(plan->tile_working_set_bytes + plan->fixed_overhead_bytes, plan->usable_bytes);
}

TEST(MemoryGovernor, SmallerBudgetYieldsSmallerTile) {
    const auto big = DeriveTilePlan(2 * kGiB, kImaxWidth, kImaxHeight, 0);
    const auto small = DeriveTilePlan(256 * kMiB, kImaxWidth, kImaxHeight, 0);
    ASSERT_TRUE(big.has_value() && small.has_value());
    EXPECT_LT(small->tile_edge, big->tile_edge)
        << "a 256 MiB budget must give a strictly smaller tile than 2 GiB";
    EXPECT_GE(small->tile_edge, kMinTileEdge);
}

TEST(MemoryGovernor, TinyBudgetDeclinesLoudly) {
    // 1 KiB cannot seat even an 8x8 tile (8*8*16 = 1024 bytes) after the fraction.
    const auto plan = DeriveTilePlan(1024, 1920, 1080, 0);
    ASSERT_FALSE(plan.has_value());
    EXPECT_EQ(plan.error().domain(), sirius::base::ErrorDomain::kDevice);
}

TEST(MemoryGovernor, OverheadLargerThanUsableBudgetDeclines) {
    // A 100 MiB budget (50 MiB usable) cannot seat a 200 MiB fixed residency.
    const auto plan = DeriveTilePlan(100 * kMiB, 1920, 1080, 200 * kMiB);
    ASSERT_FALSE(plan.has_value());
    EXPECT_EQ(plan.error().domain(), sirius::base::ErrorDomain::kDevice);
}

TEST(MemoryGovernor, TileNeverExceedsImageExtent) {
    // A generous budget against a small frame clamps the tile to the image.
    const auto plan = DeriveTilePlan(2 * kGiB, 128, 96, 0);
    ASSERT_TRUE(plan.has_value());
    EXPECT_EQ(plan->tile_edge, 96) << "tile edge is min(image dims) when budget is ample";
}

TEST(MemoryGovernor, WorkingSetMatchesTheDerivedTile) {
    const auto plan = DeriveTilePlan(512 * kMiB, kImaxWidth, kImaxHeight, 8 * kMiB);
    ASSERT_TRUE(plan.has_value());
    const std::uint64_t expected = static_cast<std::uint64_t>(plan->tile_edge) * plan->tile_edge *
                                   kTileWorkingSetBytesPerPixel;
    EXPECT_EQ(plan->tile_working_set_bytes, expected);
    EXPECT_LE(plan->tile_working_set_bytes + plan->fixed_overhead_bytes, plan->usable_bytes);
}

TEST(MemoryGovernor, EnvironmentOverrideResolvesBudget) {
    ScopedEnvironmentVariable clean_environment("SIRIUS_MEMORY_BUDGET_MB", nullptr);
    {
        ScopedEnvironmentVariable budget("SIRIUS_MEMORY_BUDGET_MB", "256");
        EXPECT_EQ(ResolveBudgetBytes(9999 * kMiB), 256 * kMiB)
            << "override wins over device budget";
    }
    EXPECT_EQ(ResolveBudgetBytes(2 * kGiB), 2 * kGiB) << "device budget used when no override";
}

}  // namespace
