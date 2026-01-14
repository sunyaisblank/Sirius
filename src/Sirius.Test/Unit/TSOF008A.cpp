// TSOF008A.cpp - Display Driver Tests
// Component ID: TSOF008A
// Tests for: OFDD001A.h

#include <gtest/gtest.h>
#include "OFDD001A.h"
#include <cmath>
#include <cstring>

using namespace sirius::offline;

//==============================================================================
// DisplayFormat Tests
//==============================================================================

TEST(DisplayFormatTest, BytesPerPixel) {
    EXPECT_EQ(bytesPerPixel(DisplayFormat::RGBA_FLOAT), 16);
    EXPECT_EQ(bytesPerPixel(DisplayFormat::RGBA_HALF), 8);
    EXPECT_EQ(bytesPerPixel(DisplayFormat::RGBA_UINT8), 4);
    EXPECT_EQ(bytesPerPixel(DisplayFormat::RGB_UINT8), 3);
}

//==============================================================================
// DisplayRegion Tests
//==============================================================================

TEST(DisplayRegionTest, BufferSize) {
    DisplayRegion region;
    region.width = 64;
    region.height = 48;
    
    EXPECT_EQ(region.bufferSize(DisplayFormat::RGBA_FLOAT), 64 * 48 * 16);
    EXPECT_EQ(region.bufferSize(DisplayFormat::RGBA_UINT8), 64 * 48 * 4);
}

//==============================================================================
// NullDisplayDriver Tests
//==============================================================================

TEST(NullDisplayDriverTest, Lifecycle) {
    NullDisplayDriver driver;
    
    EXPECT_FALSE(driver.isReady());
    
    EXPECT_TRUE(driver.initialize(1920, 1080));
    EXPECT_TRUE(driver.isReady());
    EXPECT_EQ(driver.width(), 1920);
    EXPECT_EQ(driver.height(), 1080);
    
    driver.shutdown();
    EXPECT_FALSE(driver.isReady());
}

TEST(NullDisplayDriverTest, UpdateAndPresent) {
    NullDisplayDriver driver;
    driver.initialize(256, 256);
    
    DisplayRegion region{0, 0, 64, 64, {}};
    region.pixels.resize(64 * 64 * 16);
    
    EXPECT_TRUE(driver.updateRegion(region));
    EXPECT_EQ(driver.updateCount(), 1);
    
    driver.present();
    EXPECT_EQ(driver.presentCount(), 1);
}

TEST(NullDisplayDriverTest, Format) {
    NullDisplayDriver driver;
    driver.initialize(256, 256, DisplayFormat::RGBA_UINT8);
    
    EXPECT_EQ(driver.format(), DisplayFormat::RGBA_UINT8);
}

//==============================================================================
// MemoryDisplayDriver Tests
//==============================================================================

TEST(MemoryDisplayDriverTest, Lifecycle) {
    MemoryDisplayDriver driver;
    
    EXPECT_FALSE(driver.isReady());
    
    EXPECT_TRUE(driver.initialize(64, 64, DisplayFormat::RGBA_UINT8));
    EXPECT_TRUE(driver.isReady());
    EXPECT_EQ(driver.width(), 64);
    EXPECT_EQ(driver.height(), 64);
    
    // Framebuffer should be allocated
    EXPECT_EQ(driver.framebuffer().size(), 64 * 64 * 4);
    
    driver.shutdown();
    EXPECT_FALSE(driver.isReady());
    EXPECT_TRUE(driver.framebuffer().empty());
}

TEST(MemoryDisplayDriverTest, UpdateRegion) {
    MemoryDisplayDriver driver;
    driver.initialize(8, 8, DisplayFormat::RGBA_UINT8);
    
    // Create 2x2 region at (2, 2) with known pattern
    DisplayRegion region{2, 2, 2, 2, {}};
    region.pixels.resize(2 * 2 * 4);
    
    // Fill with 255,0,0,255 (red)
    for (size_t i = 0; i < region.pixels.size(); i += 4) {
        region.pixels[i + 0] = 255;  // R
        region.pixels[i + 1] = 0;    // G
        region.pixels[i + 2] = 0;    // B
        region.pixels[i + 3] = 255;  // A
    }
    
    driver.updateRegion(region);
    
    // Check pixel at (2,2) is red
    EXPECT_EQ(driver.getPixel<uint8_t>(2, 2, 0), 255);  // R
    EXPECT_EQ(driver.getPixel<uint8_t>(2, 2, 1), 0);    // G
    
    // Check pixel at (0,0) is still black
    EXPECT_EQ(driver.getPixel<uint8_t>(0, 0, 0), 0);
}

TEST(MemoryDisplayDriverTest, UpdateOutOfBounds) {
    MemoryDisplayDriver driver;
    driver.initialize(8, 8, DisplayFormat::RGBA_UINT8);
    
    // Region partially outside bounds
    DisplayRegion region{6, 6, 4, 4, {}};
    region.pixels.resize(4 * 4 * 4, 128);
    
    // Should not crash, only copy valid portion
    EXPECT_TRUE(driver.updateRegion(region));
    
    // Pixel at (7,7) should be updated
    EXPECT_EQ(driver.getPixel<uint8_t>(7, 7, 0), 128);
}

TEST(MemoryDisplayDriverTest, Clear) {
    MemoryDisplayDriver driver;
    driver.initialize(4, 4, DisplayFormat::RGBA_UINT8);
    
    // First set some pixels
    DisplayRegion region{0, 0, 4, 4, {}};
    region.pixels.resize(4 * 4 * 4, 100);
    driver.updateRegion(region);
    
    // Clear to white
    driver.clear(1.0f, 1.0f, 1.0f, 1.0f);
    
    // All pixels should be ~255 (sRGB white)
    EXPECT_GT(driver.getPixel<uint8_t>(0, 0, 0), 200);
    EXPECT_GT(driver.getPixel<uint8_t>(2, 2, 1), 200);
}

//==============================================================================
// Factory Tests
//==============================================================================

TEST(DisplayDriverFactoryTest, CreateNull) {
    auto driver = createDisplayDriver("null");
    ASSERT_NE(driver, nullptr);
    
    driver->initialize(100, 100);
    EXPECT_TRUE(driver->isReady());
}

TEST(DisplayDriverFactoryTest, CreateMemory) {
    auto driver = createDisplayDriver("memory");
    ASSERT_NE(driver, nullptr);
    
    driver->initialize(100, 100);
    EXPECT_TRUE(driver->isReady());
}

TEST(DisplayDriverFactoryTest, UnknownDefaultsToNull) {
    auto driver = createDisplayDriver("opengl");  // Not implemented
    ASSERT_NE(driver, nullptr);
    
    // Should still work (defaults to null)
    driver->initialize(100, 100);
    EXPECT_TRUE(driver->isReady());
}

//==============================================================================
// Edge Cases
//==============================================================================

TEST(MemoryDisplayDriverTest, ZeroSizeRegion) {
    MemoryDisplayDriver driver;
    driver.initialize(64, 64, DisplayFormat::RGBA_UINT8);
    
    DisplayRegion region{0, 0, 0, 0, {}};
    
    // Should not crash
    EXPECT_TRUE(driver.updateRegion(region));
}

TEST(MemoryDisplayDriverTest, FloatFormat) {
    MemoryDisplayDriver driver;
    driver.initialize(4, 4, DisplayFormat::RGBA_FLOAT);
    
    EXPECT_EQ(driver.framebuffer().size(), 4 * 4 * 16);
    
    DisplayRegion region{0, 0, 2, 2, {}};
    region.pixels.resize(2 * 2 * 16);
    
    // Set first pixel to (1.0, 0.5, 0.25, 1.0)
    float vals[] = {1.0f, 0.5f, 0.25f, 1.0f};
    std::memcpy(region.pixels.data(), vals, 16);
    
    driver.updateRegion(region);
    
    EXPECT_FLOAT_EQ(driver.getPixel<float>(0, 0, 0), 1.0f);
    EXPECT_FLOAT_EQ(driver.getPixel<float>(0, 0, 1), 0.5f);
}
