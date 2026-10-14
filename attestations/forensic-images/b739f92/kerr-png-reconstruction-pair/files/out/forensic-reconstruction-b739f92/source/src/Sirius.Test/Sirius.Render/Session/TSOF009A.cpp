// TSOF009A.cpp - PNG Writer Tests
// Component ID: TSOF009A
// Tests for: OUPN001A.h

#include <gtest/gtest.h>
#include "Sirius.Render/Output/OUPN001A.h"
#include "Sirius.Render/Output/OUIB001A.h"

// Decode side of the round-trip; nothing else in the test target defines the
// stb_image implementation (OUPN001A.cpp defines only the writer's).
#define STB_IMAGE_IMPLEMENTATION
#include <stb_image.h>

#include <filesystem>
#include <cmath>

using namespace sirius::render;

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

//==============================================================================
// Decode Round-Trip
// The writer's earlier tests only asserted that a file appeared on disk; this
// decodes the PNG and verifies the sRGB-encoded pixel values, closing the loop
// on the transfer function (linear in, sRGB out, one encode exactly).
//==============================================================================

namespace {

uint8_t expectedSRGB8(float linear) {
    float clamped = std::clamp(linear, 0.0f, 1.0f);
    float srgb = (clamped <= 0.0031308f)
        ? 12.92f * clamped
        : 1.055f * std::pow(clamped, 1.0f / 2.4f) - 0.055f;
    return static_cast<uint8_t>(std::clamp(srgb * 255.0f + 0.5f, 0.0f, 255.0f));
}

} // namespace

TEST(PNGWriterTest, DecodeRoundTripMatchesSRGBEncoding) {
    std::string testPath = std::filesystem::temp_directory_path().string() +
                           "/sirius_test_png_roundtrip.png";
    std::filesystem::remove(testPath);

    const int W = 8, H = 4;
    ImageBuffer buffer;
    buffer.allocate(W, H);
    for (int y = 0; y < H; ++y) {
        for (int x = 0; x < W; ++x) {
            buffer.setPixel(x, y,
                            static_cast<float>(x) / (W - 1),
                            static_cast<float>(y) / (H - 1),
                            0.5f);
        }
    }

    ASSERT_TRUE(PNGWriter::write(testPath, buffer));

    int width = 0, height = 0, channels = 0;
    unsigned char* data = stbi_load(testPath.c_str(), &width, &height, &channels, 3);
    ASSERT_NE(data, nullptr) << stbi_failure_reason();
    ASSERT_EQ(width, W);
    ASSERT_EQ(height, H);

    int worst = 0;
    for (int y = 0; y < H; ++y) {
        for (int x = 0; x < W; ++x) {
            float r, g, b;
            buffer.getPixel(x, y, r, g, b);
            const uint8_t expected[3] = {expectedSRGB8(r), expectedSRGB8(g), expectedSRGB8(b)};
            for (int c = 0; c < 3; ++c) {
                int actual = data[(y * W + x) * 3 + c];
                worst = std::max(worst, std::abs(actual - expected[c]));
            }
        }
    }
    // PNG is lossless; only the quantisation rounding itself may differ.
    EXPECT_LE(worst, 1);

    stbi_image_free(data);
    std::filesystem::remove(testPath);
}
