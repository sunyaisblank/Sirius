#include "sirius/render/session/tile_scheduler.h"

#include <gtest/gtest.h>

#include <array>
#include <atomic>
#include <latch>
#include <set>
#include <thread>
#include <utility>
#include <vector>

namespace sirius::test {
namespace {
using render::TileScheduler;
using render::TileState;

TEST(TileScheduler, ReinitialiseResetsCompletionLedger) {
    sirius::render::TileScheduler scheduler;
    scheduler.Initialise(128, 64, 64);
    ASSERT_EQ(scheduler.GetTileCount(), 2);
    const sirius::render::Tile* tile = scheduler.GetNextTile();
    ASSERT_NE(tile, nullptr);
    scheduler.CompleteTile(tile->id);
    ASSERT_EQ(scheduler.GetCompletedCount(), 1);

    scheduler.Initialise(64, 64, 64);
    EXPECT_EQ(scheduler.GetTileCount(), 1);
    EXPECT_EQ(scheduler.GetCompletedCount(), 0);
    EXPECT_FALSE(scheduler.AllComplete());
}

TEST(TileScheduler, WorkGroupsPreserveEveryOriginalTileAndOwnEachRegionOnce) {
    constexpr int width = 65, height = 33;
    for (const int edge : {1, 8, 16, 32, 64}) {
        SCOPED_TRACE(edge);
        TileScheduler scheduler;
        scheduler.Initialise(width, height, edge, 32);
        const int work_edge = std::max(32, edge);
        const int expected_tiles = ((width + edge - 1) / edge) * ((height + edge - 1) / edge);
        ASSERT_EQ(scheduler.GetTileCount(), expected_tiles);
        std::set<std::pair<int, int>> groups;
        std::set<int> ids;
        std::array<int, width * height> covered{};
        for (;;) {
            const auto group = scheduler.GetNextTileGroup();
            if (group.empty()) break;
            const std::pair anchor{group.front().x / work_edge, group.front().y / work_edge};
            EXPECT_TRUE(groups.insert(anchor).second);
            EXPECT_LE(group.size(),
                      static_cast<std::size_t>((work_edge / edge) * (work_edge / edge)));
            for (const auto& tile : group) {
                EXPECT_TRUE(ids.insert(tile.id).second);
                EXPECT_EQ(tile.state, TileState::Active);
                EXPECT_EQ((std::pair{tile.x / work_edge, tile.y / work_edge}), anchor);
                EXPECT_EQ(tile.x % edge, 0);
                EXPECT_EQ(tile.y % edge, 0);
                EXPECT_EQ(tile.width, std::min(edge, width - tile.x));
                EXPECT_EQ(tile.height, std::min(edge, height - tile.y));
                for (int y = tile.y; y < tile.y + tile.height; ++y)
                    for (int x = tile.x; x < tile.x + tile.width; ++x) ++covered[y * width + x];
                scheduler.CompleteTile(tile.id);
            }
        }
        EXPECT_EQ(groups.size(), static_cast<std::size_t>(((width + work_edge - 1) / work_edge) *
                                                          ((height + work_edge - 1) / work_edge)));
        EXPECT_EQ(ids.size(), static_cast<std::size_t>(expected_tiles));
        for (const int count : covered) EXPECT_EQ(count, 1);
        EXPECT_TRUE(scheduler.AllComplete());
        EXPECT_EQ(scheduler.GetCompletedCount(), expected_tiles);
    }
}

TEST(TileScheduler, ParallelAcquisitionNeverDuplicatesAGroupOrTile) {
    constexpr int width = 257, height = 193, edge = 8, work_edge = 32;
    constexpr int groups_x = (width + work_edge - 1) / work_edge;
    constexpr int groups_y = (height + work_edge - 1) / work_edge;
    TileScheduler scheduler;
    scheduler.Initialise(width, height, edge, work_edge);
    std::vector<std::atomic<int>> tile_visits(static_cast<std::size_t>(scheduler.GetTileCount()));
    std::array<std::atomic<int>, groups_x * groups_y> group_visits{};
    std::latch start(1);
    std::array<std::thread, 4> workers;
    for (auto& worker : workers)
        worker = std::thread([&] {
            start.wait();
            for (;;) {
                const auto group = scheduler.GetNextTileGroup();
                if (group.empty()) break;
                ++group_visits[(group.front().y / work_edge) * groups_x +
                               group.front().x / work_edge];
                std::this_thread::yield();
                for (const auto& tile : group) {
                    ++tile_visits[static_cast<std::size_t>(tile.id)];
                    scheduler.CompleteTile(tile.id);
                }
            }
        });
    start.count_down();
    for (auto& worker : workers) worker.join();
    for (const auto& count : group_visits) EXPECT_EQ(count.load(), 1);
    for (const auto& count : tile_visits) EXPECT_EQ(count.load(), 1);
    EXPECT_TRUE(scheduler.AllComplete());
    EXPECT_EQ(scheduler.GetCompletedCount(), scheduler.GetTileCount());
}

TEST(TileScheduler, CompletionLedgerIsIdempotentAndResetRetainsGroups) {
    TileScheduler scheduler;
    scheduler.Initialise(65, 33, 8, 32);
    const auto first = scheduler.GetNextTileGroup();
    ASSERT_FALSE(first.empty());
    const int id = first.front().id;
    scheduler.CompleteTile(id);
    scheduler.CompleteTile(id);
    EXPECT_EQ(scheduler.GetCompletedCount(), 1);
    scheduler.FailTile(id);
    EXPECT_EQ(scheduler.GetCompletedCount(), 0);
    EXPECT_EQ(first.front().state, TileState::Failed);
    scheduler.CompleteTile(id);
    scheduler.CompleteTile(-1);
    scheduler.CompleteTile(scheduler.GetTileCount());
    scheduler.FailTile(-1);
    EXPECT_EQ(scheduler.GetCompletedCount(), 1);
    scheduler.Reset();
    EXPECT_EQ(scheduler.GetCompletedCount(), 0);
    std::set<std::pair<int, int>> groups;
    for (;;) {
        const auto group = scheduler.GetNextTileGroup();
        if (group.empty()) break;
        EXPECT_TRUE(groups.emplace(group.front().x / 32, group.front().y / 32).second);
        for (const auto& tile : group) scheduler.CompleteTile(tile.id);
    }
    EXPECT_EQ(groups.size(), 6u);
    EXPECT_TRUE(scheduler.AllComplete());
}

TEST(TileScheduler, DefaultGroupsRetainTheOriginalSpiralSequence) {
    constexpr std::array<std::pair<int, int>, 6> expected{
        {{64, 0}, {64, 64}, {0, 0}, {128, 0}, {128, 64}, {0, 64}}};
    TileScheduler scheduler;
    scheduler.Initialise(192, 128, 64);
    for (const auto& coordinate : expected) {
        const auto group = scheduler.GetNextTileGroup();
        ASSERT_EQ(group.size(), 1u);
        EXPECT_EQ((std::pair{group.front().x, group.front().y}), coordinate);
    }
    EXPECT_TRUE(scheduler.GetNextTileGroup().empty());
    scheduler.Reset();
    for (const auto& coordinate : expected) {
        const auto* tile = scheduler.GetNextTile();
        ASSERT_NE(tile, nullptr);
        EXPECT_EQ((std::pair{tile->x, tile->y}), coordinate);
    }
    EXPECT_EQ(scheduler.GetNextTile(), nullptr);
}

}  // namespace
}  // namespace sirius::test
