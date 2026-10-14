// EXR round-trip tests. Ported from TSOF011A.cpp.
//
// Writes a known HDR gradient and reads it back with tinyexr's loader (the same
// library the writer uses), validating the encode-decode pair end to end. EXR is
// the HDR interchange format: values above 1.0 must survive.

#include "sirius/render/exr_writer.h"
#include "sirius/render/image_buffer.h"

#include <gtest/gtest.h>

#include <tinyexr.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <string>

using namespace sirius::render;

namespace {

std::string tempPath(const char* name) {
    return (std::filesystem::temp_directory_path() / name).string();
}

}  // namespace

TEST(EXRRoundTripTests, HDRGradientSurvivesWriteAndRead) {
    const std::string path = tempPath("sirius_test_roundtrip.exr");
    std::filesystem::remove(path);

    const int W = 16, H = 8;
    ImageBufferRGBA buffer(W, H);
    for (int y = 0; y < H; ++y) {
        for (int x = 0; x < W; ++x) {
            // HDR values beyond 1.0 on purpose.
            buffer.SetPixel(x, y, 0.25f * x,  // R: 0 .. 3.75
                            2.0f * y,          // G: 0 .. 14
                            0.5f,              // B: constant
                            1.0f);
        }
    }

    EXRMetadata meta;
    meta.metricType = "Kerr";
    ASSERT_TRUE(EXRWriter::WriteExr(path, buffer, meta));
    ASSERT_TRUE(std::filesystem::exists(path));

    float* rgba = nullptr;
    int width = 0, height = 0;
    const char* err = nullptr;
    int ret = LoadEXR(&rgba, &width, &height, path.c_str(), &err);
    ASSERT_EQ(ret, TINYEXR_SUCCESS) << (err ? err : "unknown tinyexr error");
    ASSERT_EQ(width, W);
    ASSERT_EQ(height, H);

    // Half-float storage carries ~11 bits of significand; tolerate the
    // corresponding relative error.
    double worstRel = 0.0;
    for (int y = 0; y < H; ++y) {
        for (int x = 0; x < W; ++x) {
            size_t idx = (static_cast<size_t>(y) * W + x) * 4;
            const float expected[3] = {0.25f * x, 2.0f * y, 0.5f};
            for (int c = 0; c < 3; ++c) {
                double denom = std::max(1.0, static_cast<double>(expected[c]));
                worstRel = std::max(worstRel, std::abs(rgba[idx + c] - expected[c]) / denom);
            }
        }
    }
    EXPECT_LT(worstRel, 2e-3) << "pixel data degraded beyond half-float precision";

    std::free(rgba);
    std::filesystem::remove(path);
}

TEST(EXRRoundTripTests, WriteFailsCleanlyOnBadPath) {
    ImageBufferRGBA buffer(4, 4);
    EXPECT_FALSE(EXRWriter::WriteExr("/nonexistent-dir-sirius/out.exr", buffer));
}
