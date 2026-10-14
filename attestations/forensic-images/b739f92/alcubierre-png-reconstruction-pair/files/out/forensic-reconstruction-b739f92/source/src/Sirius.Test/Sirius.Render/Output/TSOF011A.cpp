// TSOF011A.cpp - EXR Round-Trip Tests
// Component ID: TSOF011A
// Tests for: OUEW001A (EXRWriter)
//
// The EXR writer previously had a single test asserting that a metadata
// string contained two substrings; the pixel path was never exercised. This
// suite writes a known HDR gradient and reads it back with tinyexr's loader
// (the same library the writer uses, compiled once via OUEW001A.cpp), so the
// encode-decode pair is validated end to end. EXR is the HDR interchange
// format: values above 1.0 must survive, which no 8-bit path can show.

#include <gtest/gtest.h>
#include "Sirius.Render/Output/OUEW001A.h"
#include "Sirius.Render/Output/OUIB001A.h"

#include <tinyexr.h>

#include <cstdlib>
#include <filesystem>
#include <string>

using namespace sirius::render;

namespace {

std::string tempPath(const char* name) {
    return (std::filesystem::temp_directory_path() / name).string();
}

} // namespace

TEST(EXRRoundTripTests, HDRGradientSurvivesWriteAndRead) {
    const std::string path = tempPath("sirius_test_roundtrip.exr");
    std::filesystem::remove(path);

    const int W = 16, H = 8;
    ImageBufferRGBA buffer(W, H);
    for (int y = 0; y < H; ++y) {
        for (int x = 0; x < W; ++x) {
            // HDR values beyond 1.0 on purpose
            buffer.setPixel(x, y,
                            0.25f * x,           // R: 0 .. 3.75
                            2.0f * y,            // G: 0 .. 14
                            0.5f,                // B: constant
                            1.0f);
        }
    }

    EXRMetadata meta;
    meta.metricType = "Kerr";
    ASSERT_TRUE(EXRWriter::writeEXR(path, buffer, meta));
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
                worstRel = std::max(worstRel,
                                    std::abs(rgba[idx + c] - expected[c]) / denom);
            }
        }
    }
    EXPECT_LT(worstRel, 2e-3) << "pixel data degraded beyond half-float precision";

    std::free(rgba);
    std::filesystem::remove(path);
}

TEST(EXRRoundTripTests, WriteFailsCleanlyOnBadPath) {
    ImageBufferRGBA buffer(4, 4);
    EXPECT_FALSE(EXRWriter::writeEXR("/nonexistent-dir-sirius/out.exr", buffer));
}
