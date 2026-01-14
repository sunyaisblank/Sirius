// SNDP001A.h - Display Buffer
// Component ID: SNDP001A (Session/Display)
//
// Thread-safe RGBA buffer for live preview during rendering.
// Supports tile-by-tile updates and OpenGL texture upload.

#pragma once

#include <vector>
#include <mutex>
#include <atomic>
#include <cstdint>
#include <cstring>

namespace Sirius {

//==============================================================================
// Display Buffer
//==============================================================================
class DisplayBuffer {
public:
    /// @brief Initialise buffer with dimensions
    void initialise(int width, int height) {
        std::lock_guard<std::mutex> lock(m_Mutex);
        m_Width = width;
        m_Height = height;
        m_PixelData.resize(width * height * 4, 0.0f);  // RGBA float
        m_ByteData.resize(width * height * 4, 0);       // RGBA uint8
        m_UpdateCounter = 0;
        m_Dirty = true;
    }
    
    /// @brief Update a tile region (from render thread)
    void updateTile(int tileX, int tileY, int tileWidth, int tileHeight,
                    const float* tileData) {
        std::lock_guard<std::mutex> lock(m_Mutex);
        
        for (int y = 0; y < tileHeight; ++y) {
            int destY = tileY + y;
            if (destY < 0 || destY >= m_Height) continue;
            
            for (int x = 0; x < tileWidth; ++x) {
                int destX = tileX + x;
                if (destX < 0 || destX >= m_Width) continue;
                
                int srcIdx = (y * tileWidth + x) * 4;
                int dstIdx = (destY * m_Width + destX) * 4;
                
                m_PixelData[dstIdx + 0] = tileData[srcIdx + 0];
                m_PixelData[dstIdx + 1] = tileData[srcIdx + 1];
                m_PixelData[dstIdx + 2] = tileData[srcIdx + 2];
                m_PixelData[dstIdx + 3] = tileData[srcIdx + 3];
            }
        }
        
        m_UpdateCounter++;
        m_Dirty = true;
    }
    
    /// @brief Get 8-bit RGBA data for OpenGL upload (from UI thread)
    /// @param gamma Apply gamma correction (default 2.2)
    const uint8_t* getByteData(float gamma = 2.2f) {
        std::lock_guard<std::mutex> lock(m_Mutex);
        
        if (m_Dirty) {
            float invGamma = 1.0f / gamma;
            for (size_t i = 0; i < m_PixelData.size(); i += 4) {
                m_ByteData[i + 0] = floatToByte(std::pow(m_PixelData[i + 0], invGamma));
                m_ByteData[i + 1] = floatToByte(std::pow(m_PixelData[i + 1], invGamma));
                m_ByteData[i + 2] = floatToByte(std::pow(m_PixelData[i + 2], invGamma));
                m_ByteData[i + 3] = floatToByte(m_PixelData[i + 3]);
            }
            m_Dirty = false;
        }
        
        return m_ByteData.data();
    }
    
    /// @brief Get raw float data (for file output)
    const float* getFloatData() const { return m_PixelData.data(); }
    std::vector<float>& getFloatBuffer() { return m_PixelData; }
    
    /// @brief Clear buffer to black
    void clear() {
        std::lock_guard<std::mutex> lock(m_Mutex);
        std::fill(m_PixelData.begin(), m_PixelData.end(), 0.0f);
        m_Dirty = true;
    }
    
    // Accessors
    int getWidth() const { return m_Width; }
    int getHeight() const { return m_Height; }
    uint64_t getUpdateCounter() const { return m_UpdateCounter.load(); }
    
private:
    static uint8_t floatToByte(float v) {
        if (v <= 0.0f) return 0;
        if (v >= 1.0f) return 255;
        return static_cast<uint8_t>(v * 255.0f + 0.5f);
    }
    
    int m_Width = 0;
    int m_Height = 0;
    std::vector<float> m_PixelData;     // RGBA float [0, inf)
    std::vector<uint8_t> m_ByteData;    // RGBA uint8 [0, 255]
    std::atomic<uint64_t> m_UpdateCounter{0};
    bool m_Dirty = false;
    mutable std::mutex m_Mutex;
};

} // namespace Sirius
