// BFIO001A.h - Buffer I/O Operations
// Component ID: BFIO001A (Buffer/IO)
//
// File output for rendered framebuffers.
// Supports PPM (LDR) and EXR (HDR) formats.

#pragma once

#include <string>
#include <vector>
#include <fstream>
#include <cmath>
#include <algorithm>

namespace Sirius::Buffer {

//==============================================================================
// Image Format Enum
//==============================================================================
enum class ImageFormat : uint8_t {
    PPM,    ///< Portable Pixmap (LDR, 8-bit)
    EXR,    ///< OpenEXR (HDR, 16/32-bit float)
    PNG     ///< PNG (LDR, 8-bit, compressed)
};

//==============================================================================
// Tonemapping for LDR Output
//==============================================================================
inline float reinhard(float x) {
    return x / (1.0f + x);
}

inline float aces(float x) {
    float a = 2.51f, b = 0.03f, c = 2.43f, d = 0.59f, e = 0.14f;
    return std::clamp((x*(a*x+b))/(x*(c*x+d)+e), 0.0f, 1.0f);
}

//==============================================================================
// Buffer Writer
//==============================================================================
class BufferWriter {
public:
    /// @brief Write PPM (LDR, tonemapped)
    static bool writePPM(const std::string& path, const float* rgba, 
                         int width, int height, float exposure = 1.0f) {
        std::ofstream file(path, std::ios::binary);
        if (!file) return false;
        
        file << "P6\n" << width << " " << height << "\n255\n";
        
        for (int i = 0; i < width * height; ++i) {
            float r = aces(rgba[i*4 + 0] * exposure);
            float g = aces(rgba[i*4 + 1] * exposure);
            float b = aces(rgba[i*4 + 2] * exposure);
            
            // Gamma correction
            r = std::pow(r, 1.0f/2.2f);
            g = std::pow(g, 1.0f/2.2f);
            b = std::pow(b, 1.0f/2.2f);
            
            uint8_t rgb[3] = {
                static_cast<uint8_t>(std::clamp(r * 255.0f, 0.0f, 255.0f)),
                static_cast<uint8_t>(std::clamp(g * 255.0f, 0.0f, 255.0f)),
                static_cast<uint8_t>(std::clamp(b * 255.0f, 0.0f, 255.0f))
            };
            file.write(reinterpret_cast<char*>(rgb), 3);
        }
        
        return file.good();
    }
    
    /// @brief Write raw float buffer (for debugging)
    static bool writeRaw(const std::string& path, const float* data,
                         int width, int height, int channels = 4) {
        std::ofstream file(path, std::ios::binary);
        if (!file) return false;
        
        // Header: width, height, channels
        file.write(reinterpret_cast<const char*>(&width), sizeof(int));
        file.write(reinterpret_cast<const char*>(&height), sizeof(int));
        file.write(reinterpret_cast<const char*>(&channels), sizeof(int));
        
        // Data
        size_t byteCount = width * height * channels * sizeof(float);
        file.write(reinterpret_cast<const char*>(data), byteCount);
        
        return file.good();
    }
    
    /// @brief Read raw float buffer
    static std::vector<float> readRaw(const std::string& path, 
                                       int& width, int& height, int& channels) {
        std::ifstream file(path, std::ios::binary);
        if (!file) return {};
        
        file.read(reinterpret_cast<char*>(&width), sizeof(int));
        file.read(reinterpret_cast<char*>(&height), sizeof(int));
        file.read(reinterpret_cast<char*>(&channels), sizeof(int));
        
        std::vector<float> data(width * height * channels);
        file.read(reinterpret_cast<char*>(data.data()), data.size() * sizeof(float));
        
        return data;
    }
};

} // namespace Sirius::Buffer
