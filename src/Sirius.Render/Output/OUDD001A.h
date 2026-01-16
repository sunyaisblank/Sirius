// OUDD001A.h - Display Driver
// Component ID: OUDD001A
// Interface for progressive render preview display
//
// Governance: docs/specification.md

#ifndef SIRIUS_RENDER_OUDD001A_H
#define SIRIUS_RENDER_OUDD001A_H

#include <cstddef>
#include <cstdint>
#include <cstring>
#include <memory>
#include <mutex>
#include <string>
#include <type_traits>
#include <vector>

namespace sirius::render {

//==============================================================================
// DisplayFormat Enumeration
//==============================================================================

enum class DisplayFormat {
    RGBA_FLOAT,  // 4x32-bit float (16 bytes/pixel)
    RGBA_HALF,   // 4x16-bit half-float (8 bytes/pixel)
    RGBA_UINT8,  // 4x8-bit unsigned (4 bytes/pixel)
    RGB_UINT8    // 3x8-bit unsigned (3 bytes/pixel)
};

/// Get bytes per pixel for a format
inline size_t bytesPerPixel(DisplayFormat format) {
    switch (format) {
        case DisplayFormat::RGBA_FLOAT: return 16;
        case DisplayFormat::RGBA_HALF:  return 8;
        case DisplayFormat::RGBA_UINT8: return 4;
        case DisplayFormat::RGB_UINT8:  return 3;
        default:                        return 4;
    }
}

//==============================================================================
// DisplayRegion
// A rectangular region to update on the display
//==============================================================================

struct DisplayRegion {
    int x = 0;
    int y = 0;
    int width = 0;
    int height = 0;
    std::vector<uint8_t> pixels;  // Raw pixel data

    size_t bufferSize(DisplayFormat format) const {
        return static_cast<size_t>(width) * height * bytesPerPixel(format);
    }
};

//==============================================================================
// IDisplayDriver - Abstract Interface
//==============================================================================

class IDisplayDriver {
public:
    virtual ~IDisplayDriver() = default;

    /// Initialize display surface
    virtual bool initialize(int width, int height,
                            DisplayFormat format = DisplayFormat::RGBA_FLOAT) = 0;

    /// Shutdown display
    virtual void shutdown() = 0;

    /// Check if ready to receive updates
    virtual bool isReady() const = 0;

    /// Get current dimensions
    virtual int width() const = 0;
    virtual int height() const = 0;

    /// Get current format
    virtual DisplayFormat format() const = 0;

    /// Update a region of the display
    virtual bool updateRegion(const DisplayRegion& region) = 0;

    /// Present the updated framebuffer
    virtual void present() = 0;
};

//==============================================================================
// NullDisplayDriver
// Discards all updates (for headless rendering)
//==============================================================================

class NullDisplayDriver : public IDisplayDriver {
public:
    bool initialize(int w, int h, DisplayFormat fmt = DisplayFormat::RGBA_FLOAT) override {
        m_width = w;
        m_height = h;
        m_format = fmt;
        m_ready = true;
        m_updateCount = 0;
        m_presentCount = 0;
        return true;
    }

    void shutdown() override {
        m_ready = false;
    }

    bool isReady() const override { return m_ready; }
    int width() const override { return m_width; }
    int height() const override { return m_height; }
    DisplayFormat format() const override { return m_format; }

    bool updateRegion(const DisplayRegion& /* region */) override {
        if (!m_ready) return false;
        ++m_updateCount;
        return true;
    }

    void present() override {
        if (m_ready) ++m_presentCount;
    }

    // Statistics
    int updateCount() const { return m_updateCount; }
    int presentCount() const { return m_presentCount; }

private:
    int m_width = 0;
    int m_height = 0;
    DisplayFormat m_format = DisplayFormat::RGBA_FLOAT;
    bool m_ready = false;
    int m_updateCount = 0;
    int m_presentCount = 0;
};

//==============================================================================
// MemoryDisplayDriver
// Stores display in memory buffer
//==============================================================================

class MemoryDisplayDriver : public IDisplayDriver {
public:
    bool initialize(int w, int h, DisplayFormat fmt = DisplayFormat::RGBA_FLOAT) override {
        m_width = w;
        m_height = h;
        m_format = fmt;
        m_buffer.resize(w * h * bytesPerPixel(fmt), 0);
        m_ready = true;
        return true;
    }

    void shutdown() override {
        m_ready = false;
        m_buffer.clear();
    }

    bool isReady() const override { return m_ready; }
    int width() const override { return m_width; }
    int height() const override { return m_height; }
    DisplayFormat format() const override { return m_format; }

    bool updateRegion(const DisplayRegion& region) override {
        if (!m_ready) return false;

        std::lock_guard<std::mutex> lock(m_mutex);
        size_t bpp = bytesPerPixel(m_format);

        for (int ry = 0; ry < region.height; ++ry) {
            int destY = region.y + ry;
            if (destY < 0 || destY >= m_height) continue;

            for (int rx = 0; rx < region.width; ++rx) {
                int destX = region.x + rx;
                if (destX < 0 || destX >= m_width) continue;

                size_t srcIdx = (ry * region.width + rx) * bpp;
                size_t destIdx = (destY * m_width + destX) * bpp;

                if (srcIdx + bpp <= region.pixels.size() &&
                    destIdx + bpp <= m_buffer.size()) {
                    for (size_t i = 0; i < bpp; ++i) {
                        m_buffer[destIdx + i] = region.pixels[srcIdx + i];
                    }
                }
            }
        }
        return true;
    }

    void present() override {
        // Memory driver doesn't need to present
    }

    // Access buffer
    const std::vector<uint8_t>& buffer() const { return m_buffer; }
    const std::vector<uint8_t>& framebuffer() const { return m_buffer; }

    /// Clear the display to a solid color
    void clear(float r, float g, float b, float a) {
        std::lock_guard<std::mutex> lock(m_mutex);
        size_t bpp = bytesPerPixel(m_format);

        if (m_format == DisplayFormat::RGBA_UINT8 || m_format == DisplayFormat::RGB_UINT8) {
            uint8_t ru = static_cast<uint8_t>(r * 255.0f);
            uint8_t gu = static_cast<uint8_t>(g * 255.0f);
            uint8_t bu = static_cast<uint8_t>(b * 255.0f);
            uint8_t au = static_cast<uint8_t>(a * 255.0f);

            for (size_t i = 0; i < m_buffer.size(); i += bpp) {
                m_buffer[i] = ru;
                m_buffer[i + 1] = gu;
                m_buffer[i + 2] = bu;
                if (bpp > 3) m_buffer[i + 3] = au;
            }
        } else if (m_format == DisplayFormat::RGBA_FLOAT) {
            float rgba[4] = {r, g, b, a};
            for (size_t i = 0; i < m_buffer.size(); i += 16) {
                std::memcpy(&m_buffer[i], rgba, 16);
            }
        }
    }

    /// Get a single pixel value (templated for type safety)
    template<typename T>
    T getPixel(int x, int y, int channel) const {
        if (x < 0 || x >= m_width || y < 0 || y >= m_height) {
            return T{};
        }

        size_t bpp = bytesPerPixel(m_format);
        size_t idx = (y * m_width + x) * bpp;

        if constexpr (std::is_same_v<T, uint8_t>) {
            if (m_format == DisplayFormat::RGBA_UINT8 || m_format == DisplayFormat::RGB_UINT8) {
                return m_buffer[idx + channel];
            }
            return 0;
        } else if constexpr (std::is_same_v<T, float>) {
            if (m_format == DisplayFormat::RGBA_FLOAT) {
                float val;
                std::memcpy(&val, &m_buffer[idx + channel * sizeof(float)], sizeof(float));
                return val;
            }
            return 0.0f;
        }
        return T{};
    }

private:
    int m_width = 0;
    int m_height = 0;
    DisplayFormat m_format = DisplayFormat::RGBA_FLOAT;
    std::vector<uint8_t> m_buffer;
    bool m_ready = false;
    mutable std::mutex m_mutex;
};

//==============================================================================
// Factory Function
//==============================================================================

/// Create a display driver by name
/// Supported names: "null", "memory"
/// Unknown names default to NullDisplayDriver
inline std::unique_ptr<IDisplayDriver> createDisplayDriver(const std::string& name) {
    if (name == "memory") {
        return std::make_unique<MemoryDisplayDriver>();
    }
    // Default to null driver for unknown types
    return std::make_unique<NullDisplayDriver>();
}

} // namespace sirius::render

#endif // SIRIUS_RENDER_OUDD001A_H
