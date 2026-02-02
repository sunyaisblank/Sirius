// TSOF001A.cpp - Render Session Tests
// Component ID: TSOF001A
// Tests for Phase 6: Offline render session management

#include <gtest/gtest.h>
#include <cmath>

// MIGRATION NOTE: These tests use the deprecated sirius::render::RenderSession.
// The deprecation warnings are suppressed for backwards compatibility testing.
// New tests should use Sirius::RenderSession from SNRS001A.h.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable: 4996)
#endif

#include "Sirius.Render/Session/SRRS001A.h"
#include "Sirius.Render/Camera/CMFM001A.h"
#include "Sirius.Render/Camera/CMBL001A.h"
#include "Sirius.Render/Output/OUIB001A.h"
#include "Sirius.Render/Output/OUEW001A.h"

using namespace sirius::render;

namespace {

//==============================================================================
// RenderSession Tests
//==============================================================================

TEST(RenderSessionTest, TileDecomposition) {
    RenderConfig config;
    config.imageWidth = 1920;
    config.imageHeight = 1080;
    config.tileSize = 256;

    RenderSession session(config);

    // Expected tiles: ceil(1920/256) × ceil(1080/256) = 8 × 5 = 40
    int expectedTilesX = (1920 + 255) / 256;
    int expectedTilesY = (1080 + 255) / 256;
    int expectedTotal = expectedTilesX * expectedTilesY;

    EXPECT_EQ(session.tiles().size(), expectedTotal)
        << "Should have correct number of tiles";

    EXPECT_EQ(expectedTilesX, 8);
    EXPECT_EQ(expectedTilesY, 5);
}

TEST(RenderSessionTest, TileEdgeSizes) {
    RenderConfig config;
    config.imageWidth = 1000;  // 1000 / 256 = 3.9 → 4 tiles
    config.imageHeight = 500;  // 500 / 256 = 1.9 → 2 tiles
    config.tileSize = 256;

    RenderSession session(config);

    // Last tile in X should be 1000 - 3*256 = 232 pixels wide
    // Last tile in Y should be 500 - 256 = 244 pixels tall

    const auto& tiles = session.tiles();

    // Find edge tiles
    const auto& rightEdgeTile = tiles[3];  // Tile at x=3
    const auto& bottomEdgeTile = tiles[4]; // Tile at y=1

    EXPECT_EQ(rightEdgeTile.width, 1000 - 3*256)
        << "Right edge tile should have correct width";
    EXPECT_EQ(bottomEdgeTile.height, 500 - 256)
        << "Bottom edge tile should have correct height";
}

TEST(RenderSessionTest, InitialState) {
    RenderConfig config;
    RenderSession session(config);

    EXPECT_EQ(session.state(), SessionState::IDLE)
        << "Initial state should be IDLE";

    auto progress = session.getProgress();
    EXPECT_EQ(progress.tilesCompleted, 0)
        << "No tiles should be completed initially";
    EXPECT_EQ(progress.percentComplete, 0)
        << "Should be 0% complete";
}

TEST(RenderSessionTest, ProgressCallback) {
    RenderConfig config;
    config.imageWidth = 256;
    config.imageHeight = 256;
    config.tileSize = 256;  // Single tile
    config.saveCheckpoints = false;

    RenderSession session(config);

    int callbackCount = 0;
    session.setProgressCallback([&](const RenderProgress& /* p */) {
        callbackCount++;
    });

    // Verify callback was set
    EXPECT_EQ(callbackCount, 0)
        << "Callback should not be called until tile is marked complete";

    // Mark tile completed - callback should fire
    session.start();
    session.markTileCompleted(0);

    EXPECT_EQ(callbackCount, 1)
        << "Callback should be called once after marking tile complete";
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

TEST(EXRWriterTest, SessionToBuffer) {
    RenderConfig config;
    config.imageWidth = 64;
    config.imageHeight = 64;
    config.tileSize = 64;
    config.saveCheckpoints = false;

    RenderSession session(config);

    auto buffer = EXRWriter::sessionToBuffer(session);

    EXPECT_EQ(buffer.width, 64);
    EXPECT_EQ(buffer.height, 64);
    EXPECT_EQ(buffer.pixels.size(), 64 * 64 * 3);
}

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

#pragma GCC diagnostic pop
#ifdef _MSC_VER
#pragma warning(pop)
#endif
