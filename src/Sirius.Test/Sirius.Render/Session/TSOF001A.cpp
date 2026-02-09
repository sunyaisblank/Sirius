// TSOF001A.cpp - Render Session Tests
// Component ID: TSOF001A
// Tests for Phase 6: Offline render session management

#include <gtest/gtest.h>
#include <cmath>

#include "Sirius.Render/Session/SRTL001A.h"
#include "Sirius.Render/Camera/CMFM001A.h"
#include "Sirius.Render/Camera/CMBL001A.h"
#include "Sirius.Render/Output/OUIB001A.h"
#include "Sirius.Render/Output/OUEW001A.h"

using namespace sirius::render;

namespace {

//==============================================================================
// RenderSession Tests
//==============================================================================

TEST(TileLayoutTest, TileDecomposition) {
    TileLayout layout(1920, 1080, 256);

    // Expected tiles: ceil(1920/256) × ceil(1080/256) = 8 × 5 = 40
    int expectedTilesX = (1920 + 255) / 256;
    int expectedTilesY = (1080 + 255) / 256;
    int expectedTotal = expectedTilesX * expectedTilesY;

    EXPECT_EQ(layout.tiles().size(), expectedTotal)
        << "Should have correct number of tiles";

    EXPECT_EQ(expectedTilesX, 8);
    EXPECT_EQ(expectedTilesY, 5);
}

TEST(TileLayoutTest, TileEdgeSizes) {
    TileLayout layout(1000, 500, 256);

    // Last tile in X should be 1000 - 3*256 = 232 pixels wide
    // Last tile in Y should be 500 - 256 = 244 pixels tall

    const auto& tiles = layout.tiles();

    // Find edge tiles
    const auto& rightEdgeTile = tiles[3];  // Tile at x=3
    const auto& bottomEdgeTile = tiles[4]; // Tile at y=1

    EXPECT_EQ(rightEdgeTile.width, 1000 - 3*256)
        << "Right edge tile should have correct width";
    EXPECT_EQ(bottomEdgeTile.height, 500 - 256)
        << "Bottom edge tile should have correct height";
}

TEST(TileLayoutTest, InitialState) {
    TileLayout layout(256, 256, 256);

    EXPECT_EQ(layout.tiles().size(), 1)
        << "Single-tile image should have exactly one tile";
    EXPECT_FALSE(layout.tiles()[0].completed)
        << "Tile should not be completed initially";
}

TEST(TileLayoutTest, MarkTileCompleted) {
    TileLayout layout(256, 256, 256);

    layout.markTileCompleted(0);

    EXPECT_TRUE(layout.tiles()[0].completed)
        << "Tile should be marked completed";
}

//==============================================================================
// Camera Model Tests
//==============================================================================

TEST(CameraModelTest, IMAX70mmFormat) {
    auto format = FilmFormat::IMAX70mm();

    EXPECT_NEAR(format.gateWidth, 69.6, 0.1)
        << "IMAX 70mm gate width should be ~70mm";
    EXPECT_NEAR(format.aspectRatio, 1.43, 0.01)
        << "IMAX 70mm aspect ratio should be ~1.43:1";
}

TEST(CameraModelTest, FOVCalculation) {
    auto camera = presets::IMAX_70mm_Standard();

    double hfov = camera.horizontalFOV();
    double vfov = camera.verticalFOV();

    // With 50mm lens on 69.6mm gate: FOV ≈ 70°
    EXPECT_GT(hfov, 60) << "IMAX horizontal FOV should be wide";
    EXPECT_LT(hfov, 80) << "IMAX horizontal FOV should be < 80°";
    EXPECT_LT(vfov, hfov) << "Vertical FOV should be less than horizontal";
}

TEST(CameraModelTest, ResolutionCalculation) {
    auto camera = presets::IMAX_70mm_Standard();

    // At 50 pixels/mm (high quality print)
    auto res = camera.calculateResolution(50.0);

    EXPECT_GT(res.width, 3000) << "IMAX should produce high resolution";
    EXPECT_GT(res.megapixels, 5) << "Should be > 5 megapixels";
}

TEST(CameraModelTest, RayGeneration) {
    CameraModel camera;

    // Center ray
    auto centerRay = camera.rayDirection(0, 0);
    EXPECT_NEAR(centerRay[2], -1.0, 0.1)
        << "Center ray should point along -Z";
    EXPECT_NEAR(centerRay[0], 0, 0.01)
        << "Center ray X should be ~0";
    EXPECT_NEAR(centerRay[1], 0, 0.01)
        << "Center ray Y should be ~0";

    // Corner ray
    auto cornerRay = camera.rayDirection(1, 1);
    EXPECT_GT(cornerRay[0], 0) << "Right corner should have positive X";
    EXPECT_GT(cornerRay[1], 0) << "Top corner should have positive Y";
}

//==============================================================================
// Image Buffer Tests
//==============================================================================

TEST(ImageBufferTest, Allocation) {
    ImageBuffer buffer;
    buffer.allocate(100, 50);

    EXPECT_EQ(buffer.width, 100);
    EXPECT_EQ(buffer.height, 50);
    EXPECT_EQ(buffer.pixels.size(), 100 * 50 * 3);
}

TEST(ImageBufferTest, PixelAccess) {
    ImageBuffer buffer;
    buffer.allocate(10, 10);

    buffer.setPixel(5, 3, 1.0f, 0.5f, 0.25f);

    int idx = (3 * 10 + 5) * 3;
    EXPECT_FLOAT_EQ(buffer.pixels[idx + 0], 1.0f);
    EXPECT_FLOAT_EQ(buffer.pixels[idx + 1], 0.5f);
    EXPECT_FLOAT_EQ(buffer.pixels[idx + 2], 0.25f);
}

TEST(ImageBufferTest, SpectralConversion) {
    ImageBuffer buffer;
    buffer.allocate(1, 1);

    // Daylight white point approximation (6500K)
    auto white = SpectralRadiance::blackbody(6500);
    buffer.setPixelFromSpectral(0, 0, white);

    // Should have reasonable RGB values
    EXPECT_GT(buffer.pixels[0], 0) << "Red should be > 0";
    EXPECT_GT(buffer.pixels[1], 0) << "Green should be > 0";
    EXPECT_GT(buffer.pixels[2], 0) << "Blue should be > 0";
}

//==============================================================================
// EXR Writer Tests
//==============================================================================

TEST(EXRWriterTest, MetadataGeneration) {
    EXRMetadata meta;
    meta.softwareVersion = "Sirius 1.0";
    meta.blackHoleMass = 4.3e6;  // Sgr A*
    meta.blackHoleSpin = 0.9;

    std::string header = EXRWriter::generateACESHeader(meta);

    EXPECT_TRUE(header.find("ACES AP0") != std::string::npos)
        << "Header should mention ACES";
    EXPECT_TRUE(header.find("Sirius") != std::string::npos)
        << "Header should mention software";
}

} // namespace
