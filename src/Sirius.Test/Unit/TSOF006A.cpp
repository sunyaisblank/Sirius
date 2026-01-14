// TSOF006A.cpp - Output Driver Tests
// Component ID: TSOF006A
// Tests for: OFOD001A.h, OFOD002A.h

#include <gtest/gtest.h>
#include "OFOD001A.h"
#include "OFOD002A.h"
#include <cmath>

using namespace sirius::offline;

//==============================================================================
// OutputPass Name Tests
//==============================================================================

TEST(OutputPassTest, PassNames) {
    EXPECT_STREQ(passName(OutputPass::COMBINED), "Combined");
    EXPECT_STREQ(passName(OutputPass::DEPTH), "Depth");
    EXPECT_STREQ(passName(OutputPass::REDSHIFT), "Redshift");
    EXPECT_STREQ(passName(OutputPass::TEMPERATURE), "Temperature");
    EXPECT_STREQ(passName(OutputPass::ALPHA), "Alpha");
}

//==============================================================================
// OutputTile Tests
//==============================================================================

TEST(OutputTileTest, PixelCount) {
    OutputTile tile;
    tile.width = 64;
    tile.height = 48;
    
    EXPECT_EQ(tile.pixelCount(), 64 * 48);
    EXPECT_EQ(tile.bufferSize(), 64 * 48 * 4);
}

//==============================================================================
// NullOutputDriver Tests
//==============================================================================

TEST(NullOutputDriverTest, BasicLifecycle) {
    NullOutputDriver driver;
    
    EXPECT_FALSE(driver.isActive());
    
    EXPECT_TRUE(driver.configure(1920, 1080, {OutputPass::COMBINED}));
    EXPECT_TRUE(driver.start());
    EXPECT_TRUE(driver.isActive());
    
    OutputTile tile{0, 0, 64, 64, OutputPass::COMBINED, {}};
    tile.pixels.resize(tile.bufferSize());
    EXPECT_TRUE(driver.writeTile(tile));
    
    OutputMetadata metadata;
    EXPECT_TRUE(driver.finish(metadata));
    EXPECT_FALSE(driver.isActive());
}

TEST(NullOutputDriverTest, TileCount) {
    NullOutputDriver driver;
    driver.configure(256, 256, {});
    driver.start();
    
    OutputTile tile{0, 0, 64, 64, OutputPass::COMBINED, {}};
    tile.pixels.resize(tile.bufferSize());
    
    for (int i = 0; i < 16; ++i) {
        driver.writeTile(tile);
    }
    
    EXPECT_EQ(driver.tilesWritten(), 16);
}

TEST(NullOutputDriverTest, Cancel) {
    NullOutputDriver driver;
    driver.configure(256, 256, {});
    driver.start();
    EXPECT_TRUE(driver.isActive());
    
    driver.cancel();
    EXPECT_FALSE(driver.isActive());
}

//==============================================================================
// MemoryOutputDriver Tests
//==============================================================================

TEST(MemoryOutputDriverTest, StoresPixels) {
    MemoryOutputDriver driver;
    driver.configure(4, 4, {OutputPass::COMBINED});
    driver.start();
    
    // Create a 2x2 tile at position (1,1) with known values
    OutputTile tile{1, 1, 2, 2, OutputPass::COMBINED, {}};
    tile.pixels.resize(tile.bufferSize());
    
    // Fill with test pattern
    for (size_t i = 0; i < tile.bufferSize(); i += 4) {
        tile.pixels[i + 0] = 1.0f;  // R
        tile.pixels[i + 1] = 0.5f;  // G
        tile.pixels[i + 2] = 0.25f; // B
        tile.pixels[i + 3] = 1.0f;  // A
    }
    
    driver.writeTile(tile);
    driver.finish({});
    
    // Check pixels were stored at correct positions
    EXPECT_FLOAT_EQ(driver.pixel(1, 1, 0), 1.0f);   // R
    EXPECT_FLOAT_EQ(driver.pixel(1, 1, 1), 0.5f);   // G
    EXPECT_FLOAT_EQ(driver.pixel(1, 1, 2), 0.25f);  // B
    
    // Pixels outside tile should be 0
    EXPECT_FLOAT_EQ(driver.pixel(0, 0, 0), 0.0f);
    EXPECT_FLOAT_EQ(driver.pixel(3, 3, 0), 0.0f);
}

TEST(MemoryOutputDriverTest, MultipleTiles) {
    MemoryOutputDriver driver;
    driver.configure(4, 4, {OutputPass::COMBINED});
    driver.start();
    
    // Write 4 tiles covering the whole image
    for (int ty = 0; ty < 2; ++ty) {
        for (int tx = 0; tx < 2; ++tx) {
            OutputTile tile{tx * 2, ty * 2, 2, 2, OutputPass::COMBINED, {}};
            tile.pixels.resize(tile.bufferSize());
            
            float value = static_cast<float>(ty * 2 + tx) / 4.0f;
            for (size_t i = 0; i < tile.bufferSize(); i += 4) {
                tile.pixels[i] = value;
            }
            
            driver.writeTile(tile);
        }
    }
    
    driver.finish({});
    
    EXPECT_EQ(driver.tilesWritten(), 4);
    
    // Each quadrant should have different value
    EXPECT_FLOAT_EQ(driver.pixel(0, 0, 0), 0.0f);    // Top-left
    EXPECT_FLOAT_EQ(driver.pixel(2, 0, 0), 0.25f);   // Top-right
    EXPECT_FLOAT_EQ(driver.pixel(0, 2, 0), 0.5f);    // Bottom-left
    EXPECT_FLOAT_EQ(driver.pixel(2, 2, 0), 0.75f);   // Bottom-right
}

TEST(MemoryOutputDriverTest, StoresMetadata) {
    MemoryOutputDriver driver;
    driver.configure(256, 256, {});
    driver.start();
    
    OutputMetadata meta;
    meta.width = 256;
    meta.height = 256;
    meta.samplesPerPixel = 128;
    meta.blackHoleSpin = 0.99;
    meta.comment = "Test render";
    
    driver.finish(meta);
    
    EXPECT_EQ(driver.metadata().width, 256);
    EXPECT_EQ(driver.metadata().samplesPerPixel, 128);
    EXPECT_DOUBLE_EQ(driver.metadata().blackHoleSpin, 0.99);
    EXPECT_EQ(driver.metadata().comment, "Test render");
}

//==============================================================================
// OpenEXR/PPM Output Driver Tests
//==============================================================================

TEST(OpenEXROutputDriverTest, ConfigureAndStart) {
    // Write to temp directory
    OpenEXROutputDriver driver("test_output");
    
    EXPECT_TRUE(driver.configure(64, 64, {OutputPass::COMBINED}));
    EXPECT_TRUE(driver.start());
    EXPECT_TRUE(driver.isActive());
    
    driver.cancel();
    EXPECT_FALSE(driver.isActive());
}

TEST(OpenEXROutputDriverTest, WriteTile) {
    OpenEXROutputDriver driver("test_output");
    driver.configure(64, 64, {OutputPass::COMBINED});
    driver.start();
    
    OutputTile tile{0, 0, 32, 32, OutputPass::COMBINED, {}};
    tile.pixels.resize(tile.bufferSize(), 0.5f);
    
    EXPECT_TRUE(driver.writeTile(tile));
    
    driver.cancel();
}

// Skip full write test as it requires filesystem access
// Integration tests will verify file output

//==============================================================================
// Edge Cases
//==============================================================================

TEST(OutputDriverTest, EmptyTile) {
    MemoryOutputDriver driver;
    driver.configure(64, 64, {});
    driver.start();
    
    OutputTile tile{0, 0, 0, 0, OutputPass::COMBINED, {}};
    
    // Empty tile should not crash
    EXPECT_TRUE(driver.writeTile(tile));
}

TEST(OutputDriverTest, TileOutsideBounds) {
    MemoryOutputDriver driver;
    driver.configure(64, 64, {});
    driver.start();
    
    // Tile completely outside image bounds
    OutputTile tile{100, 100, 32, 32, OutputPass::COMBINED, {}};
    tile.pixels.resize(tile.bufferSize(), 1.0f);
    
    // Should not crash (pixels outside bounds are ignored)
    EXPECT_TRUE(driver.writeTile(tile));
}

TEST(OutputDriverTest, PartialOverlapTile) {
    MemoryOutputDriver driver;
    driver.configure(64, 64, {});
    driver.start();
    
    // Tile partially overlapping image edge
    OutputTile tile{56, 56, 32, 32, OutputPass::COMBINED, {}};
    tile.pixels.resize(tile.bufferSize(), 1.0f);
    
    // Should write only the overlapping portion
    EXPECT_TRUE(driver.writeTile(tile));
}
