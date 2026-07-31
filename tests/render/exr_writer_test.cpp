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
#include <limits>
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
                            2.0f * y,         // G: 0 .. 14
                            0.5f,             // B: constant
                            1.0f);
        }
    }

    EXRMetadata meta;
    meta.metricType = "Kerr";
    ASSERT_TRUE(EXRWriter::WriteExr(path, buffer, meta));
    ASSERT_TRUE(std::filesystem::exists(path));

    EXRVersion version{};
    ASSERT_EQ(ParseEXRVersionFromFile(&version, path.c_str()), TINYEXR_SUCCESS);
    EXRHeader header{};
    InitEXRHeader(&header);
    const char* header_error = nullptr;
    ASSERT_EQ(ParseEXRHeaderFromFile(&header, &version, path.c_str(), &header_error),
              TINYEXR_SUCCESS)
        << (header_error ? header_error : "unknown EXR header error");
    std::string metric_attribute;
    std::string colour_attribute;
    for (int i = 0; i < header.num_custom_attributes; ++i) {
        const EXRAttribute& attribute = header.custom_attributes[i];
        const std::string value(reinterpret_cast<const char*>(attribute.value),
                                static_cast<std::size_t>(attribute.size));
        if (std::string(attribute.name) == "siriusMetric") metric_attribute = value;
        if (std::string(attribute.name) == "siriusColorSpace") colour_attribute = value;
    }
    EXPECT_EQ(metric_attribute, "Kerr");
    EXPECT_EQ(colour_attribute, "Sirius linear RGB");
    FreeEXRHeader(&header);

    float* rgba = nullptr;
    int width = 0, height = 0;
    const char* err = nullptr;
    int ret = LoadEXR(&rgba, &width, &height, path.c_str(), &err);
    ASSERT_EQ(ret, TINYEXR_SUCCESS) << (err ? err : "unknown tinyexr error");
    ASSERT_EQ(width, W);
    ASSERT_EQ(height, H);

    // The file stores fp32 so the round trip should be substantially tighter
    // than half-float precision.
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
    EXPECT_LT(worstRel, 1e-6) << "pixel data degraded beyond fp32 precision";

    std::free(rgba);
    std::filesystem::remove(path);
}

TEST(EXRRoundTripTests, WriteFailsCleanlyOnBadPath) {
    ImageBufferRGBA buffer(4, 4);
    EXPECT_FALSE(EXRWriter::WriteExr("/nonexistent-dir-sirius/out.exr", buffer));
}

TEST(EXRRoundTripTests, NonFiniteRadianceIsRejected) {
    ImageBufferRGBA buffer(1, 1);
    buffer.pixels[0] = std::numeric_limits<float>::infinity();
    const std::string path = tempPath("sirius_nonfinite_must_not_write.exr");
    std::filesystem::remove(path);
    EXPECT_FALSE(EXRWriter::WriteExr(path, buffer));
    EXPECT_FALSE(std::filesystem::exists(path));
}

TEST(EXRRoundTripTests, MalformedBufferShapesAreRejectedByEveryPublicWriter) {
    ImageBuffer buffer(2, 2);
    buffer.pixels.pop_back();
    const std::string exr_path = tempPath("sirius_malformed.exr");
    const std::string ppm_path = tempPath("sirius_malformed.ppm");
    const std::string pfm_path = tempPath("sirius_malformed.pfm");
    std::filesystem::remove(exr_path);
    std::filesystem::remove(ppm_path);
    std::filesystem::remove(pfm_path);

    EXPECT_FALSE(EXRWriter::WriteExr(exr_path, buffer));
    EXPECT_FALSE(EXRWriter::WritePpm(ppm_path, buffer));
    EXPECT_FALSE(EXRWriter::WritePfm(pfm_path, buffer));
    EXPECT_FALSE(std::filesystem::exists(exr_path));
    EXPECT_FALSE(std::filesystem::exists(ppm_path));
    EXPECT_FALSE(std::filesystem::exists(pfm_path));

    ImageBufferRGBA rgba(2, 2);
    rgba.pixels.pop_back();
    EXPECT_FALSE(EXRWriter::WriteExr(exr_path, rgba));
}
