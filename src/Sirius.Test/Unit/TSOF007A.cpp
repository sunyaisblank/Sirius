// TSOF007A.cpp - Checkpoint Tests
// Component ID: TSOF007A
// Tests for: OFCK001A.h

#include <gtest/gtest.h>
#include "OFCK001A.h"
#include <filesystem>
#include <cmath>

using namespace sirius::offline;

//==============================================================================
// Checkpoint Structure Tests
//==============================================================================

TEST(CheckpointTest, DefaultValues) {
    Checkpoint cp;
    
    EXPECT_EQ(cp.magic, CHECKPOINT_MAGIC);
    EXPECT_EQ(cp.version, CHECKPOINT_VERSION);
    EXPECT_EQ(cp.tilesCompleted, 0);
    EXPECT_EQ(cp.tilesX, 0);
    EXPECT_EQ(cp.tilesY, 0);
}

TEST(CheckpointTest, InitializeMask) {
    Checkpoint cp;
    cp.initializeMask(8, 6);
    
    EXPECT_EQ(cp.tilesX, 8);
    EXPECT_EQ(cp.tilesY, 6);
    EXPECT_EQ(cp.totalTiles(), 48);
    
    // Mask should be ceil(48/8) = 6 bytes
    EXPECT_EQ(cp.tileCompletionMask.size(), 6);
}

TEST(CheckpointTest, SetAndCheckTileCompleted) {
    Checkpoint cp;
    cp.initializeMask(4, 4);
    
    // Initially no tiles completed
    for (int i = 0; i < 16; ++i) {
        EXPECT_FALSE(cp.isTileCompleted(i));
    }
    
    // Mark some tiles
    cp.setTileCompleted(0);
    cp.setTileCompleted(5);
    cp.setTileCompleted(15);
    
    EXPECT_TRUE(cp.isTileCompleted(0));
    EXPECT_FALSE(cp.isTileCompleted(1));
    EXPECT_TRUE(cp.isTileCompleted(5));
    EXPECT_FALSE(cp.isTileCompleted(6));
    EXPECT_TRUE(cp.isTileCompleted(15));
    
    EXPECT_EQ(cp.tilesCompleted, 3);
}

TEST(CheckpointTest, CompletionPercent) {
    Checkpoint cp;
    cp.initializeMask(10, 10);
    
    EXPECT_DOUBLE_EQ(cp.completionPercent(), 0.0);
    
    for (int i = 0; i < 50; ++i) {
        cp.setTileCompleted(i);
    }
    
    EXPECT_DOUBLE_EQ(cp.completionPercent(), 50.0);
}

TEST(CheckpointTest, OutOfBoundsTile) {
    Checkpoint cp;
    cp.initializeMask(4, 4);
    
    // Out of bounds should not crash and return false
    EXPECT_FALSE(cp.isTileCompleted(-1));
    EXPECT_FALSE(cp.isTileCompleted(16));
    EXPECT_FALSE(cp.isTileCompleted(1000));
    
    // Setting out of bounds should not crash or change count
    int before = cp.tilesCompleted;
    cp.setTileCompleted(-1);
    cp.setTileCompleted(1000);
    EXPECT_EQ(cp.tilesCompleted, before);
}

//==============================================================================
// Config Hash Tests
//==============================================================================

TEST(ConfigHashTest, SameConfigSameHash) {
    uint64_t hash1 = hashConfig(1920, 1080, 256, 64, 0.9, 500.0, 75.0);
    uint64_t hash2 = hashConfig(1920, 1080, 256, 64, 0.9, 500.0, 75.0);
    
    EXPECT_EQ(hash1, hash2);
}

TEST(ConfigHashTest, DifferentConfigDifferentHash) {
    uint64_t base = hashConfig(1920, 1080, 256, 64, 0.9, 500.0, 75.0);
    
    // Each parameter change should produce different hash
    EXPECT_NE(base, hashConfig(1921, 1080, 256, 64, 0.9, 500.0, 75.0));  // width
    EXPECT_NE(base, hashConfig(1920, 1081, 256, 64, 0.9, 500.0, 75.0));  // height
    EXPECT_NE(base, hashConfig(1920, 1080, 128, 64, 0.9, 500.0, 75.0));  // tileSize
    EXPECT_NE(base, hashConfig(1920, 1080, 256, 128, 0.9, 500.0, 75.0)); // samples
    EXPECT_NE(base, hashConfig(1920, 1080, 256, 64, 0.8, 500.0, 75.0));  // spin
    EXPECT_NE(base, hashConfig(1920, 1080, 256, 64, 0.9, 600.0, 75.0));  // distance
    EXPECT_NE(base, hashConfig(1920, 1080, 256, 64, 0.9, 500.0, 80.0));  // inclination
}

TEST(ConfigHashTest, HashIsConsistent) {
    // Hash should be deterministic
    for (int i = 0; i < 100; ++i) {
        uint64_t h = hashConfig(1920, 1080, 256, 64, 0.9, 500.0, 75.0);
        EXPECT_EQ(h, hashConfig(1920, 1080, 256, 64, 0.9, 500.0, 75.0));
    }
}

//==============================================================================
// File I/O Tests
//==============================================================================

class CheckpointFileTest : public ::testing::Test {
protected:
    std::string testPath;
    
    void SetUp() override {
        testPath = std::filesystem::temp_directory_path().string() + "/sirius_test_checkpoint.bin";
    }
    
    void TearDown() override {
        std::filesystem::remove(testPath);
    }
};

TEST_F(CheckpointFileTest, WriteAndRead) {
    Checkpoint original;
    original.configHash = 12345;
    original.initializeMask(4, 4);
    original.setTileCompleted(0);
    original.setTileCompleted(5);
    original.setTileCompleted(10);
    original.elapsedSeconds = 42.5;
    original.imageWidth = 256;
    original.imageHeight = 256;
    original.spectralData = {1.0f, 2.0f, 3.0f, 4.0f};
    
    // Write
    ASSERT_TRUE(writeCheckpoint(testPath, original));
    
    // Read
    Checkpoint loaded;
    ASSERT_TRUE(readCheckpoint(testPath, loaded));
    
    // Verify
    EXPECT_EQ(loaded.magic, CHECKPOINT_MAGIC);
    EXPECT_EQ(loaded.version, CHECKPOINT_VERSION);
    EXPECT_EQ(loaded.configHash, 12345u);
    EXPECT_EQ(loaded.tilesX, 4);
    EXPECT_EQ(loaded.tilesY, 4);
    EXPECT_EQ(loaded.tilesCompleted, 3);
    EXPECT_TRUE(loaded.isTileCompleted(0));
    EXPECT_TRUE(loaded.isTileCompleted(5));
    EXPECT_TRUE(loaded.isTileCompleted(10));
    EXPECT_FALSE(loaded.isTileCompleted(1));
    EXPECT_DOUBLE_EQ(loaded.elapsedSeconds, 42.5);
    EXPECT_EQ(loaded.imageWidth, 256);
    EXPECT_EQ(loaded.imageHeight, 256);
    ASSERT_EQ(loaded.spectralData.size(), 4);
    EXPECT_FLOAT_EQ(loaded.spectralData[0], 1.0f);
    EXPECT_FLOAT_EQ(loaded.spectralData[3], 4.0f);
}

TEST_F(CheckpointFileTest, ReadInvalidFile) {
    // Write garbage
    std::ofstream file(testPath, std::ios::binary);
    file << "not a checkpoint";
    file.close();
    
    Checkpoint cp;
    EXPECT_FALSE(readCheckpoint(testPath, cp));
}

TEST_F(CheckpointFileTest, ReadNonexistent) {
    Checkpoint cp;
    EXPECT_FALSE(readCheckpoint("/nonexistent/path/checkpoint.bin", cp));
}

TEST_F(CheckpointFileTest, LargeSpectralData) {
    Checkpoint original;
    original.initializeMask(8, 8);
    original.imageWidth = 256;
    original.imageHeight = 256;
    original.spectralData.resize(256 * 256 * 4);
    
    // Fill with test pattern
    for (size_t i = 0; i < original.spectralData.size(); ++i) {
        original.spectralData[i] = static_cast<float>(i % 256) / 255.0f;
    }
    
    ASSERT_TRUE(writeCheckpoint(testPath, original));
    
    Checkpoint loaded;
    ASSERT_TRUE(readCheckpoint(testPath, loaded));
    
    ASSERT_EQ(loaded.spectralData.size(), original.spectralData.size());
    for (size_t i = 0; i < original.spectralData.size(); ++i) {
        EXPECT_FLOAT_EQ(loaded.spectralData[i], original.spectralData[i]);
    }
}

//==============================================================================
// Validation Tests
//==============================================================================

TEST(CheckpointValidationTest, ValidCheckpoint) {
    Checkpoint cp;
    uint64_t hash = hashConfig(1920, 1080, 256, 64, 0.9, 500.0, 75.0);
    cp.configHash = hash;
    
    EXPECT_TRUE(validateCheckpoint(cp, hash));
}

TEST(CheckpointValidationTest, InvalidMagic) {
    Checkpoint cp;
    cp.magic = 0xDEADBEEF;
    
    EXPECT_FALSE(validateCheckpoint(cp, 0));
}

TEST(CheckpointValidationTest, InvalidVersion) {
    Checkpoint cp;
    cp.version = 999;
    
    EXPECT_FALSE(validateCheckpoint(cp, 0));
}

TEST(CheckpointValidationTest, ConfigMismatch) {
    Checkpoint cp;
    cp.configHash = 12345;
    
    EXPECT_FALSE(validateCheckpoint(cp, 54321));
}

TEST(CheckpointValidationTest, ErrorMessages) {
    Checkpoint cp;
    
    cp.magic = 0xBAD;
    EXPECT_NE(checkpointValidationError(cp, 0).find("magic"), std::string::npos);
    
    cp.magic = CHECKPOINT_MAGIC;
    cp.version = 999;
    EXPECT_NE(checkpointValidationError(cp, 0).find("version"), std::string::npos);
    
    cp.version = CHECKPOINT_VERSION;
    cp.configHash = 123;
    EXPECT_NE(checkpointValidationError(cp, 456).find("Configuration"), std::string::npos);
    
    cp.configHash = 456;
    EXPECT_EQ(checkpointValidationError(cp, 456), "");
}
