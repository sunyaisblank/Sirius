// TSIN004A.cpp - Offline Renderer Integration Tests
// Component ID: TSIN004A
// Tests full integration of offline renderer components

#include <gtest/gtest.h>
#include "OFCT001A.h"  // Cancellation token
#include "OFER001A.h"  // Error accumulator
#include "OFPG001A.h"  // ETA calculator
#include "OFTM001A.h"  // Tile scheduler
#include "OFOD001A.h"  // Output driver interface
#include "OFOD002A.h"  // Output driver implementations
#include "OFCK001A.h"  // Checkpointing
#include "OFDD001A.h"  // Display driver
#include <thread>
#include <chrono>
#include <filesystem>

using namespace sirius::offline;

//==============================================================================
// Component Integration Tests
//==============================================================================

/// Test that components work together in a simulated render loop
TEST(OfflineIntegrationTest, SimulatedRenderLoop) {
    // Setup components
    CancellationToken cancelToken;
    ErrorAccumulator errors;
    ETACalculator eta;
    MemoryOutputDriver output;
    NullDisplayDriver display;
    
    // Configure
    constexpr int tilesX = 4;
    constexpr int tilesY = 4;
    constexpr int tileSize = 64;
    constexpr int imageW = tilesX * tileSize;
    constexpr int imageH = tilesY * tileSize;
    
    ASSERT_TRUE(output.configure(imageW, imageH, {OutputPass::COMBINED}));
    ASSERT_TRUE(output.start());
    ASSERT_TRUE(display.initialize(imageW, imageH));
    
    // Generate spiral tile order
    auto tileOrder = TileScheduler::generateOrder(tilesX, tilesY, TileOrder::SPIRAL_CENTER);
    EXPECT_EQ(tileOrder.size(), tilesX * tilesY);
    
    // Simulate render
    int tilesCompleted = 0;
    auto startTime = std::chrono::steady_clock::now();
    
    for (int idx : tileOrder) {
        if (cancelToken.isCancelled()) break;
        
        int tx = idx % tilesX;
        int ty = idx / tilesX;
        
        // Record elapsed time for ETA
        auto now = std::chrono::steady_clock::now();
        double elapsed = std::chrono::duration<double>(now - startTime).count();
        eta.recordTileCompletion(elapsed);
        
        // Create tile data
        OutputTile tile;
        tile.x = tx * tileSize;
        tile.y = ty * tileSize;
        tile.width = tileSize;
        tile.height = tileSize;
        tile.pass = OutputPass::COMBINED;
        tile.pixels.resize(tile.bufferSize(), 0.5f);  // Grey
        
        // Write to output
        EXPECT_TRUE(output.writeTile(tile));
        
        // Increment progress
        tilesCompleted++;
        
        // Check ETA is reasonable after a few tiles
        if (tilesCompleted > 2) {
            double remaining = eta.estimateRemaining(tilesCompleted, tilesX * tilesY);
            EXPECT_GE(remaining, 0.0);
        }
    }
    
    // Finish
    OutputMetadata meta;
    meta.width = imageW;
    meta.height = imageH;
    meta.samplesPerPixel = 1;
    EXPECT_TRUE(output.finish(meta));
    
    display.shutdown();
    
    // Verify
    EXPECT_EQ(tilesCompleted, tilesX * tilesY);
    EXPECT_FALSE(errors.hasFatalError());
    EXPECT_EQ(output.tilesWritten(), tilesX * tilesY);
}

/// Test cancellation during render
TEST(OfflineIntegrationTest, CancellationDuringRender) {
    CancellationToken cancelToken;
    NullOutputDriver output;
    
    constexpr int tilesX = 8;
    constexpr int tilesY = 8;
    
    output.configure(512, 512, {});
    output.start();
    
    auto tileOrder = TileScheduler::generateOrder(tilesX, tilesY, TileOrder::SCANLINE);
    int tilesCompleted = 0;
    
    // Cancel after 25% completion
    for (int idx : tileOrder) {
        if (cancelToken.isCancelled()) break;
        
        OutputTile tile{0, 0, 64, 64, OutputPass::COMBINED, {}};
        tile.pixels.resize(tile.bufferSize());
        output.writeTile(tile);
        
        tilesCompleted++;
        
        if (tilesCompleted == 16) {  // 25% of 64
            cancelToken.request();
        }
    }
    
    output.cancel();
    
    EXPECT_EQ(tilesCompleted, 16);
    EXPECT_TRUE(cancelToken.isCancelled());
}

/// Test checkpoint save and restore
TEST(OfflineIntegrationTest, CheckpointRoundTrip) {
    // Setup initial render state
    Checkpoint cp;
    cp.configHash = hashConfig(512, 512, 64, 32, 0.9, 500.0, 75.0);
    cp.initializeMask(8, 8);
    cp.imageWidth = 512;
    cp.imageHeight = 512;
    
    // Simulate partial completion (32 of 64 tiles)
    for (int i = 0; i < 32; ++i) {
        cp.setTileCompleted(i);
    }
    cp.elapsedSeconds = 120.5;
    
    // Add spectral data
    cp.spectralData.resize(512 * 512 * 4, 0.0f);
    for (size_t i = 0; i < cp.spectralData.size(); i += 4) {
        cp.spectralData[i] = 0.5f;  // R
    }
    
    // Save
    std::string path = std::filesystem::temp_directory_path().string() + "/integration_test.ckpt";
    ASSERT_TRUE(writeCheckpoint(path, cp));
    
    // Restore
    Checkpoint restored;
    ASSERT_TRUE(readCheckpoint(path, restored));
    
    // Validate compatibility
    uint64_t hash = hashConfig(512, 512, 64, 32, 0.9, 500.0, 75.0);
    EXPECT_TRUE(validateCheckpoint(restored, hash));
    
    // Verify state
    EXPECT_EQ(restored.tilesCompleted, 32);
    EXPECT_DOUBLE_EQ(restored.elapsedSeconds, 120.5);
    
    // Find remaining tiles
    int remaining = 0;
    for (int i = 0; i < 64; ++i) {
        if (!restored.isTileCompleted(i)) {
            remaining++;
        }
    }
    EXPECT_EQ(remaining, 32);
    
    // Cleanup
    std::filesystem::remove(path);
}

/// Test error accumulation during render
TEST(OfflineIntegrationTest, ErrorHandlingGracefulDegradation) {
    ErrorAccumulator errors;
    MemoryOutputDriver output;
    
    output.configure(256, 256, {});
    output.start();
    
    // Simulate render with some warnings
    for (int i = 0; i < 16; ++i) {
        if (i % 5 == 0) {
            errors.add(RenderError::warning("Tile " + std::to_string(i) + " had numerical issues"));
        }
        
        OutputTile tile{(i % 4) * 64, (i / 4) * 64, 64, 64, OutputPass::COMBINED, {}};
        tile.pixels.resize(tile.bufferSize(), 0.3f);
        output.writeTile(tile);
    }
    
    output.finish({});
    
    // Render completed despite warnings
    EXPECT_EQ(output.tilesWritten(), 16);
    EXPECT_GT(errors.warningCount(), 0);
    EXPECT_FALSE(errors.hasFatalError());
    EXPECT_EQ(errors.warningCount(), 4);  // Tiles 0, 5, 10, 15
}

//==============================================================================
// Performance Regression Test
//==============================================================================

TEST(OfflineIntegrationTest, TileSchedulerPerformance) {
    // Ensure tile scheduling is fast even for large grids
    constexpr int tilesX = 64;
    constexpr int tilesY = 64;
    
    auto start = std::chrono::high_resolution_clock::now();
    
    auto spiral = TileScheduler::generateOrder(tilesX, tilesY, TileOrder::SPIRAL_CENTER);
    auto hilbert = TileScheduler::generateOrder(tilesX, tilesY, TileOrder::HILBERT);
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    
    EXPECT_EQ(spiral.size(), tilesX * tilesY);
    EXPECT_EQ(hilbert.size(), tilesX * tilesY);
    
    // Should complete in under 100ms
    EXPECT_LT(duration.count(), 100) << "Tile scheduling took too long: " << duration.count() << "ms";
}
