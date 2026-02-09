// TSOF004A.cpp - Tile Scheduler Tests
// Component ID: TSOF004A
// Tests for: OFTM001A.h

#include <gtest/gtest.h>
#include "Sirius.Render/Scheduling/SCTM001A.h"
#include <set>
#include <algorithm>

using namespace sirius::render;

//==============================================================================
// Coverage Verification (All orders must include all tiles exactly once)
//==============================================================================

class TileOrderCoverageTest : public ::testing::TestWithParam<TileOrder> {};

TEST_P(TileOrderCoverageTest, AllTilesCoveredExactlyOnce) {
    TileOrder order = GetParam();
    constexpr int tilesX = 8;
    constexpr int tilesY = 6;
    constexpr int totalTiles = tilesX * tilesY;
    
    auto indices = TileScheduler::generateOrder(tilesX, tilesY, order);
    
    // Correct count
    EXPECT_EQ(indices.size(), static_cast<size_t>(totalTiles));
    
    // All unique
    std::set<int> unique(indices.begin(), indices.end());
    EXPECT_EQ(unique.size(), static_cast<size_t>(totalTiles));
    
    // All in valid range
    for (int idx : indices) {
        EXPECT_GE(idx, 0);
        EXPECT_LT(idx, totalTiles);
    }
}

INSTANTIATE_TEST_SUITE_P(
    AllOrders,
    TileOrderCoverageTest,
    ::testing::Values(
        TileOrder::SCANLINE,
        TileOrder::SPIRAL_CENTER,
        TileOrder::HILBERT,
        TileOrder::RANDOM
    )
);

//==============================================================================
// Scanline Order
//==============================================================================

TEST(TileSchedulerTest, ScanlineOrderIsSequential) {
    auto order = TileScheduler::generateOrder(4, 3, TileOrder::SCANLINE);
    
    ASSERT_EQ(order.size(), 12);
    for (int i = 0; i < 12; ++i) {
        EXPECT_EQ(order[i], i);
    }
}

//==============================================================================
// Spiral Order
//==============================================================================

TEST(TileSchedulerTest, SpiralStartsNearCenter) {
    constexpr int tilesX = 8;
    constexpr int tilesY = 6;
    
    auto order = TileScheduler::generateOrder(tilesX, tilesY, TileOrder::SPIRAL_CENTER);
    
    ASSERT_FALSE(order.empty());
    
    // First tile should be near centre
    int firstIdx = order[0];
    int firstX = firstIdx % tilesX;
    int firstY = firstIdx / tilesX;
    
    float centreX = (tilesX - 1) / 2.0f;
    float centreY = (tilesY - 1) / 2.0f;
    
    float distFromCentre = std::sqrt(
        (firstX - centreX) * (firstX - centreX) +
        (firstY - centreY) * (firstY - centreY)
    );
    
    // First tile should be within 1 unit of centre
    EXPECT_LT(distFromCentre, 1.5f);
}

TEST(TileSchedulerTest, SpiralDistancesMonotonicallyIncrease) {
    constexpr int tilesX = 8;
    constexpr int tilesY = 6;
    
    auto order = TileScheduler::generateOrder(tilesX, tilesY, TileOrder::SPIRAL_CENTER);
    
    float centreX = (tilesX - 1) / 2.0f;
    float centreY = (tilesY - 1) / 2.0f;
    
    float prevDist = -1.0f;
    for (int idx : order) {
        int x = idx % tilesX;
        int y = idx / tilesX;
        float dist = std::sqrt(
            (x - centreX) * (x - centreX) +
            (y - centreY) * (y - centreY)
        );
        // Distance should be >= previous (may be equal for same ring)
        EXPECT_GE(dist, prevDist - 0.001f);
        prevDist = dist;
    }
}

//==============================================================================
// Hilbert Order
//==============================================================================

TEST(TileSchedulerTest, HilbertCoversAllTiles) {
    constexpr int tilesX = 5;  // Non-power-of-2
    constexpr int tilesY = 7;
    
    auto order = TileScheduler::generateOrder(tilesX, tilesY, TileOrder::HILBERT);
    
    EXPECT_EQ(order.size(), static_cast<size_t>(tilesX * tilesY));
    
    std::set<int> unique(order.begin(), order.end());
    EXPECT_EQ(unique.size(), static_cast<size_t>(tilesX * tilesY));
}

TEST(TileSchedulerTest, HilbertAdjacencyProperty) {
    constexpr int tilesX = 4;  // Power of 2 for clean Hilbert
    constexpr int tilesY = 4;
    
    auto order = TileScheduler::generateOrder(tilesX, tilesY, TileOrder::HILBERT);
    
    // Adjacent indices in Hilbert order should be spatially close
    int adjacentCount = 0;
    for (size_t i = 1; i < order.size(); ++i) {
        int idx1 = order[i-1];
        int idx2 = order[i];
        
        int x1 = idx1 % tilesX, y1 = idx1 / tilesX;
        int x2 = idx2 % tilesX, y2 = idx2 / tilesX;
        
        int manhattan = std::abs(x2 - x1) + std::abs(y2 - y1);
        if (manhattan == 1) {
            adjacentCount++;
        }
    }
    
    // Most consecutive pairs should be adjacent (Hilbert property)
    double adjacencyRatio = static_cast<double>(adjacentCount) / (order.size() - 1);
    EXPECT_GT(adjacencyRatio, 0.8) << "Hilbert curve should have high adjacency";
}

//==============================================================================
// Random Order
//==============================================================================

TEST(TileSchedulerTest, RandomOrderDeterministic) {
    constexpr int tilesX = 8;
    constexpr int tilesY = 6;
    constexpr uint32_t seed = 12345;
    
    auto order1 = TileScheduler::generateOrder(tilesX, tilesY, TileOrder::RANDOM, seed);
    auto order2 = TileScheduler::generateOrder(tilesX, tilesY, TileOrder::RANDOM, seed);
    
    EXPECT_EQ(order1, order2) << "Same seed should produce same order";
}

TEST(TileSchedulerTest, RandomOrderDifferentSeeds) {
    constexpr int tilesX = 8;
    constexpr int tilesY = 6;
    
    auto order1 = TileScheduler::generateOrder(tilesX, tilesY, TileOrder::RANDOM, 1);
    auto order2 = TileScheduler::generateOrder(tilesX, tilesY, TileOrder::RANDOM, 2);
    
    EXPECT_NE(order1, order2) << "Different seeds should produce different orders";
}

TEST(TileSchedulerTest, RandomOrderNotScanline) {
    constexpr int tilesX = 8;
    constexpr int tilesY = 6;
    
    auto order = TileScheduler::generateOrder(tilesX, tilesY, TileOrder::RANDOM);
    auto scanline = TileScheduler::generateOrder(tilesX, tilesY, TileOrder::SCANLINE);
    
    EXPECT_NE(order, scanline) << "Random order should differ from scanline";
}

//==============================================================================
// Edge Cases
//==============================================================================

TEST(TileSchedulerTest, SingleTile) {
    auto order = TileScheduler::generateOrder(1, 1, TileOrder::SPIRAL_CENTER);
    ASSERT_EQ(order.size(), 1);
    EXPECT_EQ(order[0], 0);
}

TEST(TileSchedulerTest, SingleRow) {
    auto order = TileScheduler::generateOrder(10, 1, TileOrder::SPIRAL_CENTER);
    EXPECT_EQ(order.size(), 10);
    
    // Centre tiles should be first
    int firstX = order[0];
    EXPECT_GE(firstX, 4);
    EXPECT_LE(firstX, 5);
}

TEST(TileSchedulerTest, SingleColumn) {
    auto order = TileScheduler::generateOrder(1, 10, TileOrder::SPIRAL_CENTER);
    EXPECT_EQ(order.size(), 10);
}

//==============================================================================
// String Conversion
//==============================================================================

TEST(TileSchedulerTest, OrderNames) {
    EXPECT_STREQ(TileScheduler::orderName(TileOrder::SCANLINE), "scanline");
    EXPECT_STREQ(TileScheduler::orderName(TileOrder::SPIRAL_CENTER), "spiral");
    EXPECT_STREQ(TileScheduler::orderName(TileOrder::HILBERT), "hilbert");
    EXPECT_STREQ(TileScheduler::orderName(TileOrder::RANDOM), "random");
}

TEST(TileSchedulerTest, ParseOrder) {
    EXPECT_EQ(TileScheduler::parseOrder("scanline"), TileOrder::SCANLINE);
    EXPECT_EQ(TileScheduler::parseOrder("spiral"), TileOrder::SPIRAL_CENTER);
    EXPECT_EQ(TileScheduler::parseOrder("hilbert"), TileOrder::HILBERT);
    EXPECT_EQ(TileScheduler::parseOrder("random"), TileOrder::RANDOM);
    EXPECT_EQ(TileScheduler::parseOrder("unknown"), TileOrder::SCANLINE);  // Default
}
