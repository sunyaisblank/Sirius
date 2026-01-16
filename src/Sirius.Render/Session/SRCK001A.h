// SRCK001A.h - Checkpoint System
// Component ID: SRCK001A
// Enables render session persistence and resumption
//
// Governance: docs/specification.md

#ifndef SIRIUS_RENDER_SRCK001A_H
#define SIRIUS_RENDER_SRCK001A_H

#include <cstdint>
#include <fstream>
#include <string>
#include <vector>

namespace sirius::render {

//==============================================================================
// Magic Number and Version
//==============================================================================

constexpr uint32_t CHECKPOINT_MAGIC = 0x53495243;  // "SIRC" (Sirius Checkpoint)
constexpr uint32_t CHECKPOINT_VERSION = 1;

//==============================================================================
// Checkpoint Structure
// Serializable render session state for resumption
//==============================================================================

struct Checkpoint {
    uint32_t magic = CHECKPOINT_MAGIC;
    uint32_t version = CHECKPOINT_VERSION;

    // Tile grid configuration
    int tilesX = 0;
    int tilesY = 0;
    int tilesCompleted = 0;

    // Render parameters (for validation on resume)
    int imageWidth = 0;
    int imageHeight = 0;
    int tileSize = 64;

    // Configuration hash for resume validation
    uint64_t configHash = 0;

    // Elapsed render time (for progress reporting)
    double elapsedSeconds = 0.0;

    // Completion bitmask (1 bit per tile)
    std::vector<uint8_t> tileCompletionMask;

    // Spectral radiance data (optional accumulated spectra)
    std::vector<float> spectralData;

    //--------------------------------------------------------------------------
    // Mask Operations
    //--------------------------------------------------------------------------

    /// Initialize completion mask for tile grid
    void initializeMask(int tx, int ty) {
        tilesX = tx;
        tilesY = ty;
        tilesCompleted = 0;
        int totalTiles = tx * ty;
        int maskBytes = (totalTiles + 7) / 8;
        tileCompletionMask.assign(maskBytes, 0);
    }

    /// Total number of tiles in grid
    int totalTiles() const {
        return tilesX * tilesY;
    }

    /// Check if a specific tile is completed
    bool isTileCompleted(int tileIdx) const {
        if (tileIdx < 0 || tileIdx >= totalTiles()) return false;
        int byteIdx = tileIdx / 8;
        int bitIdx = tileIdx % 8;
        if (byteIdx >= static_cast<int>(tileCompletionMask.size())) return false;
        return (tileCompletionMask[byteIdx] & (1 << bitIdx)) != 0;
    }

    /// Mark a tile as completed
    void setTileCompleted(int tileIdx) {
        if (tileIdx < 0 || tileIdx >= totalTiles()) return;
        int byteIdx = tileIdx / 8;
        int bitIdx = tileIdx % 8;
        if (byteIdx >= static_cast<int>(tileCompletionMask.size())) return;

        uint8_t mask = 1 << bitIdx;
        if ((tileCompletionMask[byteIdx] & mask) == 0) {
            tileCompletionMask[byteIdx] |= mask;
            ++tilesCompleted;
        }
    }

    /// Get completion percentage
    double completionPercent() const {
        int total = totalTiles();
        if (total == 0) return 0.0;
        return 100.0 * static_cast<double>(tilesCompleted) / total;
    }

    //--------------------------------------------------------------------------
    // Validation
    //--------------------------------------------------------------------------

    /// Check if checkpoint is valid
    bool isValid() const {
        return magic == CHECKPOINT_MAGIC &&
               version == CHECKPOINT_VERSION &&
               tilesX > 0 && tilesY > 0;
    }
};

//==============================================================================
// Configuration Hash
// Creates hash from render parameters for resumption validation
//==============================================================================

inline uint64_t hashConfig(int width, int height, int spp, int tileSize,
                           double spin, double rOuter, double inclination) {
    // Simple FNV-1a-like hash
    uint64_t hash = 14695981039346656037ULL;
    auto mix = [&hash](uint64_t val) {
        hash ^= val;
        hash *= 1099511628211ULL;
    };
    mix(static_cast<uint64_t>(width));
    mix(static_cast<uint64_t>(height));
    mix(static_cast<uint64_t>(spp));
    mix(static_cast<uint64_t>(tileSize));
    mix(*reinterpret_cast<const uint64_t*>(&spin));
    mix(*reinterpret_cast<const uint64_t*>(&rOuter));
    mix(*reinterpret_cast<const uint64_t*>(&inclination));
    return hash;
}

//==============================================================================
// Checkpoint I/O Functions
//==============================================================================

/// Write checkpoint to file
inline bool writeCheckpoint(const std::string& path, const Checkpoint& cp) {
    std::ofstream file(path, std::ios::binary);
    if (!file) return false;

    // Write header
    file.write(reinterpret_cast<const char*>(&cp.magic), sizeof(cp.magic));
    file.write(reinterpret_cast<const char*>(&cp.version), sizeof(cp.version));
    file.write(reinterpret_cast<const char*>(&cp.tilesX), sizeof(cp.tilesX));
    file.write(reinterpret_cast<const char*>(&cp.tilesY), sizeof(cp.tilesY));
    file.write(reinterpret_cast<const char*>(&cp.tilesCompleted), sizeof(cp.tilesCompleted));
    file.write(reinterpret_cast<const char*>(&cp.imageWidth), sizeof(cp.imageWidth));
    file.write(reinterpret_cast<const char*>(&cp.imageHeight), sizeof(cp.imageHeight));
    file.write(reinterpret_cast<const char*>(&cp.tileSize), sizeof(cp.tileSize));
    file.write(reinterpret_cast<const char*>(&cp.configHash), sizeof(cp.configHash));
    file.write(reinterpret_cast<const char*>(&cp.elapsedSeconds), sizeof(cp.elapsedSeconds));

    // Write tile completion mask
    size_t maskSize = cp.tileCompletionMask.size();
    file.write(reinterpret_cast<const char*>(&maskSize), sizeof(maskSize));
    if (maskSize > 0) {
        file.write(reinterpret_cast<const char*>(cp.tileCompletionMask.data()), maskSize);
    }

    // Write spectral data
    size_t specSize = cp.spectralData.size();
    file.write(reinterpret_cast<const char*>(&specSize), sizeof(specSize));
    if (specSize > 0) {
        file.write(reinterpret_cast<const char*>(cp.spectralData.data()),
                   specSize * sizeof(float));
    }

    return file.good();
}

/// Read checkpoint from file
inline bool readCheckpoint(const std::string& path, Checkpoint& cp) {
    std::ifstream file(path, std::ios::binary);
    if (!file) return false;

    // Read header
    file.read(reinterpret_cast<char*>(&cp.magic), sizeof(cp.magic));
    file.read(reinterpret_cast<char*>(&cp.version), sizeof(cp.version));
    file.read(reinterpret_cast<char*>(&cp.tilesX), sizeof(cp.tilesX));
    file.read(reinterpret_cast<char*>(&cp.tilesY), sizeof(cp.tilesY));
    file.read(reinterpret_cast<char*>(&cp.tilesCompleted), sizeof(cp.tilesCompleted));
    file.read(reinterpret_cast<char*>(&cp.imageWidth), sizeof(cp.imageWidth));
    file.read(reinterpret_cast<char*>(&cp.imageHeight), sizeof(cp.imageHeight));
    file.read(reinterpret_cast<char*>(&cp.tileSize), sizeof(cp.tileSize));
    file.read(reinterpret_cast<char*>(&cp.configHash), sizeof(cp.configHash));
    file.read(reinterpret_cast<char*>(&cp.elapsedSeconds), sizeof(cp.elapsedSeconds));

    // Read tile completion mask
    size_t maskSize = 0;
    file.read(reinterpret_cast<char*>(&maskSize), sizeof(maskSize));
    if (maskSize > 0 && maskSize < 1000000) {  // Sanity check
        cp.tileCompletionMask.resize(maskSize);
        file.read(reinterpret_cast<char*>(cp.tileCompletionMask.data()), maskSize);
    }

    // Read spectral data
    size_t specSize = 0;
    file.read(reinterpret_cast<char*>(&specSize), sizeof(specSize));
    if (specSize > 0 && specSize < 100000000) {  // Sanity check
        cp.spectralData.resize(specSize);
        file.read(reinterpret_cast<char*>(cp.spectralData.data()),
                  specSize * sizeof(float));
    }

    return file.good();
}

//==============================================================================
// Checkpoint Validation
//==============================================================================

/// Validate checkpoint against expected configuration hash
inline bool validateCheckpoint(const Checkpoint& cp, uint64_t expectedHash) {
    if (cp.magic != CHECKPOINT_MAGIC) return false;
    if (cp.version != CHECKPOINT_VERSION) return false;
    if (cp.configHash != expectedHash) return false;
    return true;
}

/// Get validation error message (empty if valid)
inline std::string checkpointValidationError(const Checkpoint& cp, uint64_t expectedHash) {
    if (cp.magic != CHECKPOINT_MAGIC) {
        return "Invalid checkpoint magic number";
    }
    if (cp.version != CHECKPOINT_VERSION) {
        return "Incompatible checkpoint version";
    }
    if (cp.configHash != expectedHash) {
        return "Configuration mismatch - cannot resume";
    }
    return "";
}

} // namespace sirius::render

#endif // SIRIUS_RENDER_SRCK001A_H
