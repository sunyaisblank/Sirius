// OUPN001A.h - PNG Output Driver
// Component ID: OUPN001A
// PNG file output using stb_image_write
//
// Governance: docs/specification.md

#ifndef SIRIUS_RENDER_OUPN001A_H
#define SIRIUS_RENDER_OUPN001A_H

#include "OUMD001A.h"
#include "OUIB001A.h"
#include <atomic>
#include <mutex>
#include <string>
#include <cmath>
#include <algorithm>

namespace sirius::render {

//==============================================================================
// PNGOutputDriver
// Writes rendered image to PNG file with sRGB gamma correction
//==============================================================================

class PNGOutputDriver : public IOutputDriver {
public:
    explicit PNGOutputDriver(const std::string& outputPath)
        : m_outputPath(outputPath) {}

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
        // Clear buffer
        std::fill(m_pixels.begin(), m_pixels.end(), 0.0f);
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
        return writePNG();
    }

    void cancel() override {
        m_active.store(false, std::memory_order_release);
    }

    bool isActive() const override {
        return m_active.load(std::memory_order_acquire);
    }

    int tilesWritten() const { return m_tilesWritten; }
    const std::string& outputPath() const { return m_outputPath; }

private:
    //--------------------------------------------------------------------------
    // sRGB Gamma Correction
    //--------------------------------------------------------------------------

    static uint8_t toSRGB8(float linear) {
        float clamped = std::clamp(linear, 0.0f, 1.0f);
        float srgb;
        if (clamped <= 0.0031308f) {
            srgb = 12.92f * clamped;
        } else {
            srgb = 1.055f * std::pow(clamped, 1.0f / 2.4f) - 0.055f;
        }
        return static_cast<uint8_t>(std::clamp(srgb * 255.0f + 0.5f, 0.0f, 255.0f));
    }

    //--------------------------------------------------------------------------
    // PNG Writing
    //--------------------------------------------------------------------------

    bool writePNG() const;

    std::string m_outputPath;
    int m_width = 0;
    int m_height = 0;
    std::vector<OutputPass> m_passes;
    std::vector<float> m_pixels;  // RGBA float
    std::atomic<bool> m_active{false};
    mutable std::mutex m_mutex;
    int m_tilesWritten = 0;
    OutputMetadata m_metadata;
};

//==============================================================================
// PNGWriter
// Static utility for writing PNG files
//==============================================================================

class PNGWriter {
public:
    /// Write ImageBuffer (RGB) to PNG file
    static bool write(const std::string& path, const ImageBuffer& buffer);

    /// Write ImageBufferRGBA to PNG file
    static bool write(const std::string& path, const ImageBufferRGBA& buffer);

    /// Write raw RGBA float data to PNG file
    static bool writeRGBA(const std::string& path, int width, int height,
                          const float* pixels);

    /// Write raw RGB float data to PNG file
    static bool writeRGB(const std::string& path, int width, int height,
                         const float* pixels);

private:
    //--------------------------------------------------------------------------
    // sRGB Gamma Correction
    //--------------------------------------------------------------------------

    static uint8_t toSRGB8(float linear) {
        float clamped = std::clamp(linear, 0.0f, 1.0f);
        float srgb;
        if (clamped <= 0.0031308f) {
            srgb = 12.92f * clamped;
        } else {
            srgb = 1.055f * std::pow(clamped, 1.0f / 2.4f) - 0.055f;
        }
        return static_cast<uint8_t>(std::clamp(srgb * 255.0f + 0.5f, 0.0f, 255.0f));
    }
};

} // namespace sirius::render

#endif // SIRIUS_RENDER_OUPN001A_H
