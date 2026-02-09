// TSOF010A.cpp - Offline Renderer Integration Tests
// Component ID: TSOF010A
// Tests for: SROR001A.h (OfflineRenderer)

#include <gtest/gtest.h>
#include "Sirius.Render/Session/SROR001A.h"
#include <atomic>
#include <cmath>

using namespace sirius::render;

//==============================================================================
// Basic Lifecycle Tests
//==============================================================================

TEST(OfflineRendererTest, DefaultConstruction) {
    OfflineRenderer renderer;
    EXPECT_FALSE(renderer.isCancelled());
    EXPECT_FALSE(renderer.hasErrors());
}

TEST(OfflineRendererTest, Configuration) {
    OfflineRenderer renderer;

    OfflineRenderConfig config;
    config.imageWidth = 256;
    config.imageHeight = 256;
    config.tileSize = 64;
    config.tileOrder = TileOrder::SPIRAL_CENTER;
    config.outputDriver = "memory";

    renderer.configure(config);

    // Tile layout should be created
    ASSERT_NE(renderer.tileLayout(), nullptr);
    EXPECT_EQ(renderer.tileLayout()->tiles().size(), 16);  // 4x4 tiles
}

//==============================================================================
// Rendering Tests
//==============================================================================

TEST(OfflineRendererTest, SimpleRender) {
    OfflineRenderer renderer;

    OfflineRenderConfig config;
    config.imageWidth = 64;
    config.imageHeight = 64;
    config.tileSize = 32;  // 2x2 = 4 tiles
    config.outputDriver = "memory";

    renderer.configure(config);

    // Set up a simple render callback that fills with solid color
    renderer.setTileRenderCallback([](const TileInfo& tile,
                                      std::vector<float>& pixels,
                                      const RenderConfig& /* config */) {
        // Fill with red
        for (int y = 0; y < tile.height; ++y) {
            for (int x = 0; x < tile.width; ++x) {
                size_t idx = (y * tile.width + x) * 4;
                pixels[idx + 0] = 1.0f;  // R
                pixels[idx + 1] = 0.0f;  // G
                pixels[idx + 2] = 0.0f;  // B
                pixels[idx + 3] = 1.0f;  // A
            }
        }
        return true;
    });

    auto result = renderer.render();

    EXPECT_TRUE(result.success);
    EXPECT_EQ(result.tilesCompleted, 4);
    EXPECT_EQ(result.tilesTotal, 4);
    EXPECT_GT(result.elapsedSeconds, 0.0);
    EXPECT_TRUE(result.errors.empty());
}

TEST(OfflineRendererTest, ProgressCallback) {
    OfflineRenderer renderer;

    OfflineRenderConfig config;
    config.imageWidth = 64;
    config.imageHeight = 64;
    config.tileSize = 32;  // 4 tiles
    config.outputDriver = "memory";

    renderer.configure(config);

    std::atomic<int> progressCount{0};
    std::atomic<int> lastTilesCompleted{0};

    renderer.setProgressCallback([&](const RenderProgress& progress) {
        progressCount++;
        lastTilesCompleted = progress.tilesCompleted;
        EXPECT_GE(progress.percentComplete, 0.0);
        EXPECT_LE(progress.percentComplete, 100.0);
    });

    renderer.setTileRenderCallback([](const TileInfo&,
                                      std::vector<float>& pixels,
                                      const RenderConfig&) {
        std::fill(pixels.begin(), pixels.end(), 0.5f);
        return true;
    });

    auto result = renderer.render();

    EXPECT_TRUE(result.success);
    EXPECT_EQ(progressCount.load(), 4);  // One per tile
    EXPECT_EQ(lastTilesCompleted.load(), 4);
}

TEST(OfflineRendererTest, TileOrdering) {
    OfflineRenderer renderer;

    OfflineRenderConfig config;
    config.imageWidth = 128;
    config.imageHeight = 128;
    config.tileSize = 64;  // 2x2 = 4 tiles
    config.tileOrder = TileOrder::SPIRAL_CENTER;
    config.outputDriver = "memory";

    renderer.configure(config);

    std::vector<int> renderOrder;
    renderer.setTileRenderCallback([&](const TileInfo& tile,
                                       std::vector<float>& pixels,
                                       const RenderConfig&) {
        renderOrder.push_back(tile.index);
        std::fill(pixels.begin(), pixels.end(), 0.5f);
        return true;
    });

    auto result = renderer.render();

    EXPECT_TRUE(result.success);
    EXPECT_EQ(renderOrder.size(), 4);

    // Spiral order for 2x2 should start from center-ish
    // Verify all tiles were rendered
    std::sort(renderOrder.begin(), renderOrder.end());
    EXPECT_EQ(renderOrder[0], 0);
    EXPECT_EQ(renderOrder[1], 1);
    EXPECT_EQ(renderOrder[2], 2);
    EXPECT_EQ(renderOrder[3], 3);
}

//==============================================================================
// Error Handling Tests
//==============================================================================

TEST(OfflineRendererTest, NoCallbackError) {
    OfflineRenderer renderer;

    OfflineRenderConfig config;
    config.imageWidth = 64;
    config.imageHeight = 64;
    config.outputDriver = "memory";

    renderer.configure(config);
    // Don't set callback

    auto result = renderer.render();

    EXPECT_FALSE(result.success);
    EXPECT_FALSE(result.errors.empty());
}

TEST(OfflineRendererTest, TileRenderFailure) {
    OfflineRenderer renderer;

    OfflineRenderConfig config;
    config.imageWidth = 64;
    config.imageHeight = 64;
    config.tileSize = 32;  // 4 tiles
    config.outputDriver = "memory";

    renderer.configure(config);

    int tileCount = 0;
    renderer.setTileRenderCallback([&](const TileInfo&,
                                       std::vector<float>& pixels,
                                       const RenderConfig&) {
        tileCount++;
        if (tileCount == 2) {
            return false;  // Fail second tile
        }
        std::fill(pixels.begin(), pixels.end(), 0.5f);
        return true;
    });

    auto result = renderer.render();

    // Should still complete (just skip failed tile)
    EXPECT_TRUE(result.success);
    EXPECT_EQ(result.tilesCompleted, 3);  // 4 - 1 failed
}

TEST(OfflineRendererTest, TileRenderException) {
    OfflineRenderer renderer;

    OfflineRenderConfig config;
    config.imageWidth = 64;
    config.imageHeight = 64;
    config.tileSize = 64;  // 1 tile
    config.outputDriver = "memory";

    renderer.configure(config);

    std::atomic<bool> errorCallbackCalled{false};
    renderer.setErrorCallback([&](const RenderError& error) {
        errorCallbackCalled = true;
        EXPECT_EQ(error.severity, ErrorSeverity::ERROR);
    });

    renderer.setTileRenderCallback([](const TileInfo&,
                                      std::vector<float>&,
                                      const RenderConfig&) -> bool {
        throw std::runtime_error("Test exception");
    });

    auto result = renderer.render();

    EXPECT_TRUE(result.success);  // Still succeeds overall
    EXPECT_EQ(result.tilesCompleted, 0);
    EXPECT_TRUE(errorCallbackCalled);
}

//==============================================================================
// Output Driver Tests
//==============================================================================

TEST(OfflineRendererTest, NullOutputDriver) {
    OfflineRenderer renderer;

    OfflineRenderConfig config;
    config.imageWidth = 64;
    config.imageHeight = 64;
    config.tileSize = 64;
    config.outputDriver = "null";

    renderer.configure(config);
    renderer.setTileRenderCallback([](const TileInfo&,
                                      std::vector<float>& pixels,
                                      const RenderConfig&) {
        std::fill(pixels.begin(), pixels.end(), 1.0f);
        return true;
    });

    auto result = renderer.render();

    EXPECT_TRUE(result.success);
    EXPECT_EQ(result.imageData, nullptr);  // Null driver has no image
}

TEST(OfflineRendererTest, MemoryOutputDriver) {
    OfflineRenderer renderer;

    OfflineRenderConfig config;
    config.imageWidth = 4;
    config.imageHeight = 4;
    config.tileSize = 4;  // Single tile
    config.outputDriver = "memory";

    renderer.configure(config);
    renderer.setTileRenderCallback([](const TileInfo& tile,
                                      std::vector<float>& pixels,
                                      const RenderConfig&) {
        // Fill with gradient
        for (int y = 0; y < tile.height; ++y) {
            for (int x = 0; x < tile.width; ++x) {
                size_t idx = (y * tile.width + x) * 4;
                pixels[idx + 0] = static_cast<float>(x) / 3.0f;  // R
                pixels[idx + 1] = static_cast<float>(y) / 3.0f;  // G
                pixels[idx + 2] = 0.5f;                          // B
                pixels[idx + 3] = 1.0f;                          // A
            }
        }
        return true;
    });

    auto result = renderer.render();

    EXPECT_TRUE(result.success);
    ASSERT_NE(result.imageData, nullptr);

    // Check pixel values
    EXPECT_FLOAT_EQ(result.imageData->pixel(0, 0, 0), 0.0f);    // R at (0,0)
    EXPECT_FLOAT_EQ(result.imageData->pixel(3, 0, 0), 1.0f);    // R at (3,0)
    EXPECT_FLOAT_EQ(result.imageData->pixel(0, 3, 1), 1.0f);    // G at (0,3)
    EXPECT_FLOAT_EQ(result.imageData->pixel(0, 0, 2), 0.5f);    // B at (0,0)
}

//==============================================================================
// Cancellation Tests
//==============================================================================

TEST(OfflineRendererTest, Cancellation) {
    OfflineRenderer renderer;

    OfflineRenderConfig config;
    config.imageWidth = 256;
    config.imageHeight = 256;
    config.tileSize = 32;  // 8x8 = 64 tiles
    config.outputDriver = "memory";

    renderer.configure(config);

    std::atomic<int> tilesRendered{0};
    renderer.setTileRenderCallback([&](const TileInfo&,
                                       std::vector<float>& pixels,
                                       const RenderConfig&) {
        tilesRendered++;
        if (tilesRendered >= 10) {
            renderer.cancel();  // Cancel after 10 tiles
        }
        std::fill(pixels.begin(), pixels.end(), 0.5f);
        return true;
    });

    auto result = renderer.render();

    EXPECT_FALSE(result.success);  // Cancelled = not successful
    EXPECT_LE(result.tilesCompleted, 12);  // Should stop shortly after cancel
    EXPECT_LT(result.tilesCompleted, 64);  // Should not complete all tiles
}

//==============================================================================
// Large Render Test
//==============================================================================

TEST(OfflineRendererTest, LargeRender) {
    OfflineRenderer renderer;

    OfflineRenderConfig config;
    config.imageWidth = 512;
    config.imageHeight = 512;
    config.tileSize = 64;  // 8x8 = 64 tiles
    config.tileOrder = TileOrder::HILBERT;
    config.outputDriver = "memory";

    renderer.configure(config);
    renderer.setTileRenderCallback([](const TileInfo& tile,
                                      std::vector<float>& pixels,
                                      const RenderConfig&) {
        // Simple distance-based pattern
        int cx = tile.x + tile.width / 2;
        int cy = tile.y + tile.height / 2;

        for (int y = 0; y < tile.height; ++y) {
            for (int x = 0; x < tile.width; ++x) {
                int gx = tile.x + x;
                int gy = tile.y + y;

                float dist = std::sqrt(static_cast<float>((gx - 256) * (gx - 256) +
                                                          (gy - 256) * (gy - 256)));
                float val = 1.0f - std::min(dist / 256.0f, 1.0f);

                size_t idx = (y * tile.width + x) * 4;
                pixels[idx + 0] = val;
                pixels[idx + 1] = val * 0.5f;
                pixels[idx + 2] = val * 0.25f;
                pixels[idx + 3] = 1.0f;
            }
        }
        return true;
    });

    auto result = renderer.render();

    EXPECT_TRUE(result.success);
    EXPECT_EQ(result.tilesCompleted, 64);
    EXPECT_GT(result.elapsedSeconds, 0.0);
    ASSERT_NE(result.imageData, nullptr);

    // Center should be bright
    float centerR = result.imageData->pixel(256, 256, 0);
    EXPECT_GT(centerR, 0.9f);

    // Corner should be dark
    float cornerR = result.imageData->pixel(0, 0, 0);
    EXPECT_LT(cornerR, 0.3f);
}
