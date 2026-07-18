// PNG writer tests. Ported from TSOF009A.cpp.
//
// Exercises the sRGB encode on write and decodes the result back (the transfer
// function is applied exactly once: linear in, sRGB out). Assertions and
// tolerances are identical to the legacy suite.

#include "sirius/render/image_buffer.h"
#include "sirius/render/png_writer.h"

#include <gtest/gtest.h>

// Decode side of the round-trip. The stb_image decoder implementation is
// provided once by the render library (stb_impl.cpp); this TU only needs the
// declarations and links against those symbols, so it must NOT define the
// implementation macro itself.
#include <stb_image.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <string>

using namespace sirius::render;

TEST(PNGWriterTest, WriteImageBuffer) {
    std::string testPath =
        std::filesystem::temp_directory_path().string() + "/sirius_test_png_buffer.png";

    std::filesystem::remove(testPath);

    ImageBuffer buffer;
    buffer.Allocate(8, 8);

    for (int y = 0; y < 8; ++y) {
        for (int x = 0; x < 8; ++x) {
            float r = static_cast<float>(x) / 7.0f;
            float g = static_cast<float>(y) / 7.0f;
            float b = 0.5f;
            buffer.SetPixel(x, y, r, g, b);
        }
    }

    EXPECT_TRUE(PNGWriter::Write(testPath, buffer));
    EXPECT_TRUE(std::filesystem::exists(testPath));

    std::filesystem::remove(testPath);
}

TEST(PNGWriterTest, WriteImageBufferRGBA) {
    std::string testPath =
        std::filesystem::temp_directory_path().string() + "/sirius_test_png_rgba.png";

    std::filesystem::remove(testPath);

    ImageBufferRGBA buffer;
    buffer.Allocate(8, 8);

    for (int y = 0; y < 8; ++y) {
        for (int x = 0; x < 8; ++x) {
            float val = ((x + y) % 2 == 0) ? 1.0f : 0.0f;
            float alpha = 0.5f + 0.5f * val;
            buffer.SetPixel(x, y, val, val, val, alpha);
        }
    }

    EXPECT_TRUE(PNGWriter::Write(testPath, buffer));
    EXPECT_TRUE(std::filesystem::exists(testPath));

    std::filesystem::remove(testPath);
}

TEST(PNGWriterTest, GammaCorrection) {
    // Linear 0.5 should become ~0.735 in sRGB.
    std::string testPath =
        std::filesystem::temp_directory_path().string() + "/sirius_test_png_gamma.png";

    std::filesystem::remove(testPath);

    ImageBuffer buffer;
    buffer.Allocate(1, 1);
    buffer.SetPixel(0, 0, 0.5f, 0.5f, 0.5f);

    EXPECT_TRUE(PNGWriter::Write(testPath, buffer));

    EXPECT_TRUE(std::filesystem::exists(testPath));

    auto fileSize = std::filesystem::file_size(testPath);
    EXPECT_GT(fileSize, 0);
    EXPECT_LT(fileSize, 1000);  // 1x1 PNG should be tiny.

    std::filesystem::remove(testPath);
}

TEST(PNGWriterTest, LargeImage) {
    std::string testPath =
        std::filesystem::temp_directory_path().string() + "/sirius_test_png_large.png";

    std::filesystem::remove(testPath);

    ImageBuffer buffer;
    buffer.Allocate(256, 256);

    for (int y = 0; y < 256; ++y) {
        for (int x = 0; x < 256; ++x) {
            float r = static_cast<float>(x) / 255.0f;
            float g = static_cast<float>(y) / 255.0f;
            float b = static_cast<float>(x + y) / 510.0f;
            buffer.SetPixel(x, y, r, g, b);
        }
    }

    EXPECT_TRUE(PNGWriter::Write(testPath, buffer));
    EXPECT_TRUE(std::filesystem::exists(testPath));

    auto fileSize = std::filesystem::file_size(testPath);
    EXPECT_GT(fileSize, 1000);
    EXPECT_LT(fileSize, 500000);

    std::filesystem::remove(testPath);
}

TEST(PNGWriterTest, EmptyBuffer) {
    std::string testPath =
        std::filesystem::temp_directory_path().string() + "/sirius_test_png_empty.png";

    ImageBuffer buffer;
    // Not allocated - should fail gracefully.

    EXPECT_FALSE(PNGWriter::Write(testPath, buffer));
}

TEST(PNGWriterTest, ZeroSizeBuffer) {
    std::string testPath =
        std::filesystem::temp_directory_path().string() + "/sirius_test_png_zero.png";

    EXPECT_FALSE(PNGWriter::WriteRgb(testPath, 0, 0, nullptr));
}

TEST(PNGWriterTest, NullPixels) {
    std::string testPath =
        std::filesystem::temp_directory_path().string() + "/sirius_test_png_null.png";

    EXPECT_FALSE(PNGWriter::WriteRgb(testPath, 10, 10, nullptr));
}

namespace {

uint8_t expectedSRGB8(float linear) {
    float clamped = std::clamp(linear, 0.0f, 1.0f);
    float srgb = (clamped <= 0.0031308f) ? 12.92f * clamped
                                         : 1.055f * std::pow(clamped, 1.0f / 2.4f) - 0.055f;
    return static_cast<uint8_t>(std::clamp(srgb * 255.0f + 0.5f, 0.0f, 255.0f));
}

}  // namespace

TEST(PNGWriterTest, DecodeRoundTripMatchesSRGBEncoding) {
    std::string testPath =
        std::filesystem::temp_directory_path().string() + "/sirius_test_png_roundtrip.png";
    std::filesystem::remove(testPath);

    const int W = 8, H = 4;
    ImageBuffer buffer;
    buffer.Allocate(W, H);
    for (int y = 0; y < H; ++y) {
        for (int x = 0; x < W; ++x) {
            buffer.SetPixel(x, y, static_cast<float>(x) / (W - 1), static_cast<float>(y) / (H - 1),
                            0.5f);
        }
    }

    ASSERT_TRUE(PNGWriter::Write(testPath, buffer));

    int width = 0, height = 0, channels = 0;
    unsigned char* data = stbi_load(testPath.c_str(), &width, &height, &channels, 3);
    ASSERT_NE(data, nullptr) << stbi_failure_reason();
    ASSERT_EQ(width, W);
    ASSERT_EQ(height, H);

    int worst = 0;
    for (int y = 0; y < H; ++y) {
        for (int x = 0; x < W; ++x) {
            float r, g, b;
            buffer.GetPixel(x, y, r, g, b);
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
