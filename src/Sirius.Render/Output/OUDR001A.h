// OUDR001A.h - Output Driver Implementations
// Component ID: OUDR001A
// Concrete output drivers: Null, Memory, File
//
// Governance: docs/specification.md

#ifndef SIRIUS_RENDER_OUDR001A_H
#define SIRIUS_RENDER_OUDR001A_H

#include "OUMD001A.h"
#include "OUIB001A.h"
#include "OUEW001A.h"
#include <atomic>
#include <mutex>

namespace sirius::render {

//==============================================================================
// NullOutputDriver
// Discards all output (useful for benchmarking)
//==============================================================================

class NullOutputDriver : public IOutputDriver {
public:
    bool configure(int width, int height,
                   const std::vector<OutputPass>& passes) override {
        m_width = width;
        m_height = height;
        m_passes = passes;
        m_tilesWritten = 0;
        return true;
    }

    bool start() override {
        m_active.store(true, std::memory_order_release);
        return true;
    }

    bool writeTile(const OutputTile& /* tile */) override {
        if (!m_active.load(std::memory_order_acquire)) return false;
        ++m_tilesWritten;
        return true;
    }

    bool finish(const OutputMetadata& /* metadata */) override {
        m_active.store(false, std::memory_order_release);
        return true;
    }

    void cancel() override {
        m_active.store(false, std::memory_order_release);
    }

    bool isActive() const override {
        return m_active.load(std::memory_order_acquire);
    }

    // Statistics
    int tilesWritten() const { return m_tilesWritten; }

private:
    int m_width = 0;
    int m_height = 0;
    std::vector<OutputPass> m_passes;
    std::atomic<bool> m_active{false};
    int m_tilesWritten = 0;
};

//==============================================================================
// MemoryOutputDriver
// Stores rendered image in memory
//==============================================================================

class MemoryOutputDriver : public IOutputDriver {
public:
    bool configure(int width, int height,
                   const std::vector<OutputPass>& passes) override {
        m_width = width;
        m_height = height;
        m_passes = passes;
        m_pixels.resize(static_cast<size_t>(width) * height * 4, 0.0f);
        m_tilesWritten = 0;
        return true;
    }

    bool start() override {
        m_active.store(true, std::memory_order_release);
        return true;
    }

    bool writeTile(const OutputTile& tile) override {
        if (!m_active.load(std::memory_order_acquire)) return false;

        std::lock_guard<std::mutex> lock(m_mutex);

        // Copy tile pixels to framebuffer
        for (int ty = 0; ty < tile.height; ++ty) {
            int destY = tile.y + ty;
            if (destY < 0 || destY >= m_height) continue;

            for (int tx = 0; tx < tile.width; ++tx) {
                int destX = tile.x + tx;
                if (destX < 0 || destX >= m_width) continue;

                size_t srcIdx = (ty * tile.width + tx) * 4;
                size_t destIdx = (destY * m_width + destX) * 4;

                if (srcIdx + 3 < tile.pixels.size() &&
                    destIdx + 3 < m_pixels.size()) {
                    m_pixels[destIdx + 0] = tile.pixels[srcIdx + 0];
                    m_pixels[destIdx + 1] = tile.pixels[srcIdx + 1];
                    m_pixels[destIdx + 2] = tile.pixels[srcIdx + 2];
                    m_pixels[destIdx + 3] = tile.pixels[srcIdx + 3];
                }
            }
        }

        ++m_tilesWritten;
        return true;
    }

    bool finish(const OutputMetadata& metadata) override {
        m_metadata = metadata;
        m_active.store(false, std::memory_order_release);
        return true;
    }

    void cancel() override {
        m_active.store(false, std::memory_order_release);
    }

    bool isActive() const override {
        return m_active.load(std::memory_order_acquire);
    }

    // Access framebuffer
    const std::vector<float>& pixels() const { return m_pixels; }
    int width() const { return m_width; }
    int height() const { return m_height; }
    int tilesWritten() const { return m_tilesWritten; }

    /// Get stored metadata
    const OutputMetadata& metadata() const { return m_metadata; }

    /// Get single pixel channel value at (x, y, channel)
    /// channel: 0=R, 1=G, 2=B, 3=A
    float pixel(int x, int y, int channel) const {
        if (x < 0 || x >= m_width || y < 0 || y >= m_height ||
            channel < 0 || channel > 3) {
            return 0.0f;
        }
        size_t idx = (y * m_width + x) * 4 + channel;
        return m_pixels[idx];
    }

    /// Get pixel RGBA at (x, y)
    std::array<float, 4> getPixel(int x, int y) const {
        if (x < 0 || x >= m_width || y < 0 || y >= m_height) {
            return {0.0f, 0.0f, 0.0f, 0.0f};
        }
        size_t idx = (y * m_width + x) * 4;
        return {
            m_pixels[idx + 0],
            m_pixels[idx + 1],
            m_pixels[idx + 2],
            m_pixels[idx + 3]
        };
    }

private:
    int m_width = 0;
    int m_height = 0;
    std::vector<OutputPass> m_passes;
    std::vector<float> m_pixels;
    std::atomic<bool> m_active{false};
    mutable std::mutex m_mutex;
    int m_tilesWritten = 0;
    OutputMetadata m_metadata;
};

//==============================================================================
// OpenEXROutputDriver
// Writes rendered image to OpenEXR file via MemoryOutputDriver + EXRWriter
//==============================================================================

class OpenEXROutputDriver : public IOutputDriver {
public:
    explicit OpenEXROutputDriver(const std::string& basePath)
        : m_basePath(basePath) {}

    bool configure(int width, int height,
                   const std::vector<OutputPass>& passes) override {
        return m_memDriver.configure(width, height, passes);
    }

    bool start() override {
        return m_memDriver.start();
    }

    bool writeTile(const OutputTile& tile) override {
        return m_memDriver.writeTile(tile);
    }

    bool finish(const OutputMetadata& metadata) override {
        // Convert accumulated pixels to ImageBufferRGBA
        ImageBufferRGBA buffer;
        buffer.allocate(m_memDriver.width(), m_memDriver.height());
        const auto& px = m_memDriver.pixels();
        if (buffer.pixels.size() == px.size()) {
            std::copy(px.begin(), px.end(), buffer.pixels.begin());
        }

        // Build EXR metadata from OutputMetadata
        EXRMetadata meta;
        meta.metricType = metadata.metricName;
        meta.blackHoleSpin = metadata.blackHoleSpin;
        meta.blackHoleMass = metadata.massParameter;
        meta.samplesPerPixel = metadata.samplesPerPixel;
        meta.renderTimeSeconds = metadata.renderTimeSeconds;

        EXRWriter::writeEXR(m_basePath, buffer, meta);

        return m_memDriver.finish(metadata);
    }

    void cancel() override {
        m_memDriver.cancel();
    }

    bool isActive() const override {
        return m_memDriver.isActive();
    }

    // Statistics
    int tilesWritten() const { return m_memDriver.tilesWritten(); }
    const std::string& basePath() const { return m_basePath; }

private:
    std::string m_basePath;
    MemoryOutputDriver m_memDriver;
};

} // namespace sirius::render

#endif // SIRIUS_RENDER_OUDR001A_H
