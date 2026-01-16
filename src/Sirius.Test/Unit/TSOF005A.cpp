// TSOF005A.cpp - ETA Calculator Tests
// Component ID: TSOF005A
// Tests for: OFPG001A.h

#include <gtest/gtest.h>
#include "Sirius.Render/Session/SRPG001A.h"
#include <cmath>

using namespace sirius::render;

//==============================================================================
// Basic Functionality
//==============================================================================

TEST(ETACalculatorTest, InitiallyZero) {
    ETACalculator eta;
    EXPECT_EQ(eta.sampleCount(), 0);
    EXPECT_DOUBLE_EQ(eta.tilesPerSecond(), 0.0);
    EXPECT_DOUBLE_EQ(eta.elapsedSeconds(), 0.0);
    EXPECT_DOUBLE_EQ(eta.estimateRemaining(0, 100), 0.0);
}

TEST(ETACalculatorTest, SingleSampleNotEnough) {
    ETACalculator eta;
    eta.recordTileCompletion(1.0);
    
    EXPECT_EQ(eta.sampleCount(), 1);
    EXPECT_DOUBLE_EQ(eta.estimateRemaining(1, 100), 0.0);
    EXPECT_DOUBLE_EQ(eta.tilesPerSecond(), 0.0);  // Need at least 2 samples
}

TEST(ETACalculatorTest, TwoSamplesGivesEstimate) {
    ETACalculator eta;
    eta.recordTileCompletion(0.0);
    eta.recordTileCompletion(1.0);
    
    EXPECT_EQ(eta.sampleCount(), 2);
    EXPECT_GT(eta.tilesPerSecond(), 0.0);
    EXPECT_GT(eta.estimateRemaining(2, 100), 0.0);
}

TEST(ETACalculatorTest, Reset) {
    ETACalculator eta;
    eta.recordTileCompletion(0.0);
    eta.recordTileCompletion(1.0);
    
    eta.reset();
    
    EXPECT_EQ(eta.sampleCount(), 0);
    EXPECT_DOUBLE_EQ(eta.tilesPerSecond(), 0.0);
}

//==============================================================================
// ETA Estimation Accuracy
//==============================================================================

TEST(ETACalculatorTest, LinearProgressEstimate) {
    ETACalculator eta;
    
    // Simulate constant-speed render: 1 tile per second
    for (int i = 0; i <= 50; ++i) {
        eta.recordTileCompletion(static_cast<double>(i));
    }
    
    // 50 tiles done, 50 remaining at 1 tile/sec = ~50 seconds remaining
    double remaining = eta.estimateRemaining(50, 100);
    
    // Should be close to 50 seconds (allow some weighted average variation)
    EXPECT_GT(remaining, 40.0);
    EXPECT_LT(remaining, 60.0);
}

TEST(ETACalculatorTest, AcceleratingRenderFavorsRecent) {
    ETACalculator eta;
    
    // First 10 tiles: slow (1 tile per 2 seconds)
    for (int i = 0; i < 10; ++i) {
        eta.recordTileCompletion(i * 2.0);
    }
    
    // Next 10 tiles: fast (1 tile per 0.5 seconds)
    for (int i = 0; i < 10; ++i) {
        eta.recordTileCompletion(20.0 + i * 0.5);
    }
    
    // 20 tiles done at time 25.0
    double remaining = eta.estimateRemaining(20, 100);
    
    // Should favor recent speed (0.5 sec/tile) over old (2 sec/tile)
    // 80 tiles * 0.5 sec = 40 seconds if fully favoring recent
    // 80 tiles * 2 sec = 160 seconds if fully favoring old
    // Weighted average should be closer to recent
    EXPECT_LT(remaining, 100.0) << "Should favor recent faster speed";
}

TEST(ETACalculatorTest, DeceleratingRenderFavorsRecent) {
    ETACalculator eta;
    
    // First 10 tiles: fast (1 tile per 0.5 seconds)
    for (int i = 0; i < 10; ++i) {
        eta.recordTileCompletion(i * 0.5);
    }
    
    // Next 10 tiles: slow (1 tile per 2 seconds)
    for (int i = 0; i < 10; ++i) {
        eta.recordTileCompletion(5.0 + i * 2.0);
    }
    
    // 20 tiles done at time 25.0
    double remaining = eta.estimateRemaining(20, 100);
    
    // Should favor recent slower speed
    // 80 tiles * 2 sec = 160 seconds if fully favoring recent
    EXPECT_GT(remaining, 80.0) << "Should favor recent slower speed";
}

//==============================================================================
// Tiles Per Second
//==============================================================================

TEST(ETACalculatorTest, TilesPerSecondAccuracy) {
    ETACalculator eta;
    
    // 10 tiles in 5 seconds = 2 tiles/sec
    for (int i = 0; i <= 10; ++i) {
        eta.recordTileCompletion(i * 0.5);
    }
    
    double rate = eta.tilesPerSecond();
    EXPECT_NEAR(rate, 2.0, 0.1);
}

TEST(ETACalculatorTest, ElapsedSeconds) {
    ETACalculator eta;
    
    eta.recordTileCompletion(0.0);
    eta.recordTileCompletion(1.5);
    eta.recordTileCompletion(3.0);
    
    EXPECT_DOUBLE_EQ(eta.elapsedSeconds(), 3.0);
}

//==============================================================================
// Edge Cases
//==============================================================================

TEST(ETACalculatorTest, ZeroRemainingTiles) {
    ETACalculator eta;
    eta.recordTileCompletion(0.0);
    eta.recordTileCompletion(1.0);
    
    // All tiles complete
    EXPECT_DOUBLE_EQ(eta.estimateRemaining(100, 100), 0.0);
}

TEST(ETACalculatorTest, MoreThanMaxSamples) {
    ETACalculator eta;
    
    // Record more than MAX_SAMPLES (64)
    for (int i = 0; i < 100; ++i) {
        eta.recordTileCompletion(static_cast<double>(i));
    }
    
    // Should cap at 64 samples
    EXPECT_LE(eta.sampleCount(), 64u);
    
    // Should still give valid estimate
    double remaining = eta.estimateRemaining(100, 200);
    EXPECT_GT(remaining, 0.0);
}

TEST(ETACalculatorTest, VerySmallIntervals) {
    ETACalculator eta;
    
    // Very fast tiles (microseconds apart)
    for (int i = 0; i < 10; ++i) {
        eta.recordTileCompletion(i * 0.0001);
    }
    
    // Should handle without division issues
    double rate = eta.tilesPerSecond();
    EXPECT_GT(rate, 0.0);
    EXPECT_FALSE(std::isnan(rate));
    EXPECT_FALSE(std::isinf(rate));
}

TEST(ETACalculatorTest, NegativeTimeNotExpected) {
    ETACalculator eta;
    
    // If times go backwards (shouldn't happen), should still not crash
    eta.recordTileCompletion(10.0);
    eta.recordTileCompletion(9.0);  // Time went backwards
    eta.recordTileCompletion(11.0);
    
    // Should not crash, estimate might be weird but no NaN/Inf
    double remaining = eta.estimateRemaining(3, 100);
    EXPECT_FALSE(std::isnan(remaining));
    EXPECT_FALSE(std::isinf(remaining));
}
