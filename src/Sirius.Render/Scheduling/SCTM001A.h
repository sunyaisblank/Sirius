// SCTM001A.h - Tile Scheduler
// Component ID: SCTM001A
// Generates tile rendering orders for batch rendering
//
// Governance: docs/specification.md

#ifndef SIRIUS_RENDER_SCTM001A_H
#define SIRIUS_RENDER_SCTM001A_H

#include <algorithm>
#include <cmath>
#include <cstring>
#include <numeric>
#include <random>
#include <vector>

namespace sirius::render {

//==============================================================================
// TileOrder Enumeration
//==============================================================================

enum class TileOrder {
    SCANLINE,       // Left-to-right, top-to-bottom (default)
    SPIRAL_CENTER,  // Start from center, spiral outward
    HILBERT,        // Hilbert curve (locality-preserving)
    RANDOM          // Randomized order
};

//==============================================================================
// TileScheduler
// Generates tile indices in various orders for batch rendering
//==============================================================================

class TileScheduler {
public:
    /// Get string name for tile order
    static const char* orderName(TileOrder order) {
        switch (order) {
            case TileOrder::SCANLINE: return "scanline";
            case TileOrder::SPIRAL_CENTER: return "spiral";
            case TileOrder::HILBERT: return "hilbert";
            case TileOrder::RANDOM: return "random";
            default: return "scanline";
        }
    }

    /// Parse tile order from string (defaults to SCANLINE if unknown)
    static TileOrder parseOrder(const char* name) {
        if (std::strcmp(name, "scanline") == 0) return TileOrder::SCANLINE;
        if (std::strcmp(name, "spiral") == 0) return TileOrder::SPIRAL_CENTER;
        if (std::strcmp(name, "hilbert") == 0) return TileOrder::HILBERT;
        if (std::strcmp(name, "random") == 0) return TileOrder::RANDOM;
        return TileOrder::SCANLINE;  // Default
    }

    /// Generate tile indices in specified order
    /// @param tilesX Number of tiles horizontally
    /// @param tilesY Number of tiles vertically
    /// @param order Desired tile order
    /// @param seed Random seed (only used for RANDOM order)
    /// @return Vector of tile indices in rendering order
    static std::vector<int> generateOrder(int tilesX, int tilesY, TileOrder order,
                                          unsigned int seed = 42) {
        int totalTiles = tilesX * tilesY;
        std::vector<int> indices;
        indices.reserve(totalTiles);

        switch (order) {
            case TileOrder::SCANLINE:
                return generateScanline(tilesX, tilesY);

            case TileOrder::SPIRAL_CENTER:
                return generateSpiralCenter(tilesX, tilesY);

            case TileOrder::HILBERT:
                return generateHilbert(tilesX, tilesY);

            case TileOrder::RANDOM:
                return generateRandom(tilesX, tilesY, seed);

            default:
                return generateScanline(tilesX, tilesY);
        }
    }

private:
    //--------------------------------------------------------------------------
    // Scanline Order
    //--------------------------------------------------------------------------

    static std::vector<int> generateScanline(int tilesX, int tilesY) {
        int totalTiles = tilesX * tilesY;
        std::vector<int> indices(totalTiles);
        std::iota(indices.begin(), indices.end(), 0);
        return indices;
    }

    //--------------------------------------------------------------------------
    // Spiral Center Order
    //--------------------------------------------------------------------------

    static std::vector<int> generateSpiralCenter(int tilesX, int tilesY) {
        int totalTiles = tilesX * tilesY;

        // Calculate distances from center for each tile
        float centerX = (tilesX - 1) / 2.0f;
        float centerY = (tilesY - 1) / 2.0f;

        std::vector<std::pair<float, int>> tileDistances;
        tileDistances.reserve(totalTiles);

        for (int y = 0; y < tilesY; ++y) {
            for (int x = 0; x < tilesX; ++x) {
                int idx = y * tilesX + x;
                float dx = x - centerX;
                float dy = y - centerY;
                float dist = std::sqrt(dx * dx + dy * dy);
                // Add angle component to create spiral rather than rings
                float angle = std::atan2(dy, dx);
                float sortKey = dist + angle * 0.01f;
                tileDistances.emplace_back(sortKey, idx);
            }
        }

        // Sort by distance from center
        std::sort(tileDistances.begin(), tileDistances.end());

        std::vector<int> indices;
        indices.reserve(totalTiles);
        for (const auto& p : tileDistances) {
            indices.push_back(p.second);
        }
        return indices;
    }

    //--------------------------------------------------------------------------
    // Hilbert Curve Order
    //--------------------------------------------------------------------------

    static std::vector<int> generateHilbert(int tilesX, int tilesY) {
        // Find smallest power of 2 that contains the grid
        int size = 1;
        while (size < tilesX || size < tilesY) {
            size *= 2;
        }

        std::vector<std::pair<int, int>> hilbertCoords;
        hilbertCoords.reserve(size * size);

        // Generate Hilbert curve coordinates
        for (int i = 0; i < size * size; ++i) {
            auto [x, y] = hilbertD2xy(size, i);
            if (x < tilesX && y < tilesY) {
                hilbertCoords.emplace_back(x, y);
            }
        }

        // Convert to tile indices
        std::vector<int> indices;
        indices.reserve(tilesX * tilesY);
        for (const auto& [x, y] : hilbertCoords) {
            indices.push_back(y * tilesX + x);
        }
        return indices;
    }

    // Convert Hilbert curve index to (x, y) coordinates
    static std::pair<int, int> hilbertD2xy(int n, int d) {
        int x = 0, y = 0;
        for (int s = 1; s < n; s *= 2) {
            int rx = 1 & (d / 2);
            int ry = 1 & (d ^ rx);
            if (ry == 0) {
                if (rx == 1) {
                    x = s - 1 - x;
                    y = s - 1 - y;
                }
                std::swap(x, y);
            }
            x += s * rx;
            y += s * ry;
            d /= 4;
        }
        return {x, y};
    }

    //--------------------------------------------------------------------------
    // Random Order
    //--------------------------------------------------------------------------

    static std::vector<int> generateRandom(int tilesX, int tilesY, unsigned int seed) {
        int totalTiles = tilesX * tilesY;
        std::vector<int> indices(totalTiles);
        std::iota(indices.begin(), indices.end(), 0);

        std::mt19937 rng(seed);
        std::shuffle(indices.begin(), indices.end(), rng);
        return indices;
    }
};

} // namespace sirius::render

#endif // SIRIUS_RENDER_SCTM001A_H
