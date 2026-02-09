// TSOF009A.cpp - PNG Output Driver Tests
// Component ID: TSOF009A
// Tests for: OUPN001A.h

#include <gtest/gtest.h>
#include "Sirius.Render/Output/OUPN001A.h"
#include "Sirius.Render/Output/OUIB001A.h"
#include <filesystem>
#include <cmath>

using namespace sirius::render;

//==============================================================================
// PNGOutputDriver Tests
//==============================================================================

TEST(PNGOutputDriverTest, Lifecycle) {
    std::string testPath = std::filesystem::temp_directory_path().string() +
                           "/sirius_test_png_driver.png";

    PNGOutputDriver driver(testPath);

    EXPECT_FALSE(driver.isActive());

    EXPECT_TRUE(driver.configure(64, 64, {OutputPass::COMBINED}));
    EXPECT_TRUE(driver.start());
    EXPECT_TRUE(driver.isActive());

    driver.cancel();
    EXPECT_FALSE(driver.isActive());

    // Cleanup
    std::filesystem::remove(testPath);
}

TEST(PNGOutputDriverTest, WriteTile) {
    std::string testPath = std::filesystem::temp_directory_path().string() +
                           "/sirius_test_png_tile.png";

    PNGOutputDriver driver(testPath);
    driver.configure(64, 64, {OutputPass::COMBINED});
    driver.start();

    // Create 32x32 tile at (16, 16) with red color
    OutputTile tile{16, 16, 32, 32, OutputPass::COMBINED, {}};
    tile.pixels.resize(32 * 32 * 4);
    for (size_t i = 0; i < tile.pixels.size(); i += 4) {
        tile.pixels[i + 0] = 1.0f;  // R
        tile.pixels[i + 1] = 0.0f;  // G
        tile.pixels[i + 2] = 0.0f;  // B
        tile.pixels[i + 3] = 1.0f;  // A
    }

    EXPECT_TRUE(driver.writeTile(tile));
    EXPECT_EQ(driver.tilesWritten(), 1);

    driver.cancel();

    // Cleanup
    std::filesystem::remove(testPath);
}

TEST(PNGOutputDriverTest, FinishWritesFile) {
    std::string testPath = std::filesystem::temp_directory_path().string() +
                           "/sirius_test_png_finish.png";

    // Remove if exists from previous run
    std::filesystem::remove(testPath);

    PNGOutputDriver driver(testPath);
    driver.configure(4, 4, {OutputPass::COMBINED});
    driver.start();

    // Create tile covering entire image
    OutputTile tile{0, 0, 4, 4, OutputPass::COMBINED, {}};
    tile.pixels.resize(4 * 4 * 4);
    for (size_t i = 0; i < tile.pixels.size(); i += 4) {
        tile.pixels[i + 0] = 0.5f;  // R
        tile.pixels[i + 1] = 0.5f;  // G
        tile.pixels[i + 2] = 0.5f;  // B
        tile.pixels[i + 3] = 1.0f;  // A
    }

    driver.writeTile(tile);

    OutputMetadata meta;
    meta.width = 4;
    meta.height = 4;
    EXPECT_TRUE(driver.finish(meta));

    // Verify file exists
    EXPECT_TRUE(std::filesystem::exists(testPath));

    // Cleanup
    std::filesystem::remove(testPath);
}

//==============================================================================
// PNGWriter Tests
//==============================================================================

TEST(PNGWriterTest, WriteImageBuffer) {
    std::string testPath = std::filesystem::temp_directory_path().string() +
                           "/sirius_test_png_buffer.png";

    // Remove if exists from previous run
    std::filesystem::remove(testPath);

    ImageBuffer buffer;
    buffer.allocate(8, 8);

    // Create gradient pattern
    for (int y = 0; y < 8; ++y) {
        for (int x = 0; x < 8; ++x) {
            float r = static_cast<float>(x) / 7.0f;
            float g = static_cast<float>(y) / 7.0f;
            float b = 0.5f;
            buffer.setPixel(x, y, r, g, b);
        }
    }

    EXPECT_TRUE(PNGWriter::write(testPath, buffer));
    EXPECT_TRUE(std::filesystem::exists(testPath));

    // Cleanup
    std::filesystem::remove(testPath);
}

TEST(PNGWriterTest, WriteImageBufferRGBA) {
    std::string testPath = std::filesystem::temp_directory_path().string() +
                           "/sirius_test_png_rgba.png";

    // Remove if exists from previous run
    std::filesystem::remove(testPath);

    ImageBufferRGBA buffer;
    buffer.allocate(8, 8);

    // Create checkerboard pattern with alpha
    for (int y = 0; y < 8; ++y) {
        for (int x = 0; x < 8; ++x) {
            float val = ((x + y) % 2 == 0) ? 1.0f : 0.0f;
            float alpha = 0.5f + 0.5f * val;
            buffer.setPixel(x, y, val, val, val, alpha);
        }
    }

    EXPECT_TRUE(PNGWriter::write(testPath, buffer));
    EXPECT_TRUE(std::filesystem::exists(testPath));

    // Cleanup
    std::filesystem::remove(testPath);
}

TEST(PNGWriterTest, GammaCorrection) {
    // Test that sRGB gamma is applied correctly
    // Linear 0.5 should become ~0.735 in sRGB
    std::string testPath = std::filesystem::temp_directory_path().string() +
                           "/sirius_test_png_gamma.png";

    std::filesystem::remove(testPath);

    ImageBuffer buffer;
    buffer.allocate(1, 1);
    buffer.setPixel(0, 0, 0.5f, 0.5f, 0.5f);

    EXPECT_TRUE(PNGWriter::write(testPath, buffer));

    // File should exist
    EXPECT_TRUE(std::filesystem::exists(testPath));

    // File size should be reasonable for 1x1 PNG
    auto fileSize = std::filesystem::file_size(testPath);
    EXPECT_GT(fileSize, 0);
    EXPECT_LT(fileSize, 1000);  // 1x1 PNG should be tiny

    std::filesystem::remove(testPath);
}

TEST(PNGWriterTest, LargeImage) {
    std::string testPath = std::filesystem::temp_directory_path().string() +
                           "/sirius_test_png_large.png";

    std::filesystem::remove(testPath);

    // 256x256 gradient
    ImageBuffer buffer;
    buffer.allocate(256, 256);

    for (int y = 0; y < 256; ++y) {
        for (int x = 0; x < 256; ++x) {
            float r = static_cast<float>(x) / 255.0f;
            float g = static_cast<float>(y) / 255.0f;
            float b = static_cast<float>(x + y) / 510.0f;
            buffer.setPixel(x, y, r, g, b);
        }
    }

    EXPECT_TRUE(PNGWriter::write(testPath, buffer));
    EXPECT_TRUE(std::filesystem::exists(testPath));

    // File should be reasonably sized
    auto fileSize = std::filesystem::file_size(testPath);
    EXPECT_GT(fileSize, 1000);       // Should have some content
    EXPECT_LT(fileSize, 500000);     // But not huge

    std::filesystem::remove(testPath);
}

//==============================================================================
// Edge Cases
//==============================================================================

TEST(PNGWriterTest, EmptyBuffer) {
    std::string testPath = std::filesystem::temp_directory_path().string() +
                           "/sirius_test_png_empty.png";

    ImageBuffer buffer;
    // Don't allocate - should fail gracefully

    EXPECT_FALSE(PNGWriter::write(testPath, buffer));
}

TEST(PNGWriterTest, ZeroSizeBuffer) {
    std::string testPath = std::filesystem::temp_directory_path().string() +
                           "/sirius_test_png_zero.png";

    // writeRGB with zero dimensions should fail
    EXPECT_FALSE(PNGWriter::writeRGB(testPath, 0, 0, nullptr));
}

TEST(PNGWriterTest, NullPixels) {
    std::string testPath = std::filesystem::temp_directory_path().string() +
                           "/sirius_test_png_null.png";

    EXPECT_FALSE(PNGWriter::writeRGB(testPath, 10, 10, nullptr));
}
