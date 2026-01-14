// PPOP001A.h - Post-Processing Operators
// Component ID: PPOP001A (Post-Processing/Operators)
//
// Tonemapping, bloom, and colour grading operators for HDR image processing.
// All operators work on linear HDR float buffers.
//
// Reference: Academy Color Encoding System (ACES), Reinhard et al. (2002)

#pragma once

#include <vector>
#include <cmath>
#include <algorithm>
#include <functional>

namespace Sirius {

//==============================================================================
// Tonemapping Operator Type
//==============================================================================
enum class TonemapType {
    None,           ///< No tonemapping (clamp only)
    Reinhard,       ///< Reinhard global operator
    ACES,           ///< ACES filmic curve (sRGB output)
    Filmic,         ///< Hable/Uncharted 2 filmic curve
    AgX,            ///< AgX (Blender 3.6+)
    Exposure        ///< Simple exposure + gamma
};

//==============================================================================
// Post-Processing Configuration
//==============================================================================
struct PostProcessConfig {
    // Tonemapping
    TonemapType tonemapper = TonemapType::ACES;
    float exposure = 1.0f;      ///< Exposure multiplier
    float gamma = 2.2f;         ///< Display gamma
    
    // Bloom
    bool enableBloom = false;
    float bloomIntensity = 0.3f;
    float bloomThreshold = 0.8f;
    int bloomRadius = 8;
    
    // Colour grading
    float saturation = 1.0f;
    float contrast = 1.0f;
    float lift = 0.0f;          ///< Shadow lift
    float gain = 1.0f;          ///< Highlight gain
};

//==============================================================================
// Tonemapping Functions
//==============================================================================
namespace Tonemap {

/// @brief ACES filmic curve (approximation)
/// @param x Linear HDR value
/// @return Tonemapped value [0, 1]
inline float ACES(float x) {
    const float a = 2.51f;
    const float b = 0.03f;
    const float c = 2.43f;
    const float d = 0.59f;
    const float e = 0.14f;
    return std::clamp((x * (a * x + b)) / (x * (c * x + d) + e), 0.0f, 1.0f);
}

/// @brief Reinhard global operator
inline float Reinhard(float x) {
    return x / (1.0f + x);
}

/// @brief Hable filmic curve (Uncharted 2)
inline float Filmic(float x) {
    const float A = 0.15f;
    const float B = 0.50f;
    const float C = 0.10f;
    const float D = 0.20f;
    const float E = 0.02f;
    const float F = 0.30f;
    auto curve = [=](float t) {
        return ((t * (A * t + C * B) + D * E) / (t * (A * t + B) + D * F)) - E / F;
    };
    const float whitePoint = 11.2f;
    return curve(x) / curve(whitePoint);
}

/// @brief AgX base contrast curve
inline float AgX(float x) {

    // AgX sigmoid with punch
    // const float offset = 0.0f;
    // const float slope = 1.0f;
    const float power = 1.3f;
    // const float saturation = 1.0f;
    
    float v = x;
    v = std::max(0.0f, v);
    v = std::pow(v, power);

    v = v / (v + 1.0f);
    return std::clamp(v, 0.0f, 1.0f);
}

/// @brief Apply tonemapping to RGB values
inline void apply(float& r, float& g, float& b, TonemapType type, float exposure) {
    r *= exposure;
    g *= exposure;
    b *= exposure;
    
    switch (type) {
        case TonemapType::ACES:
            r = ACES(r);
            g = ACES(g);
            b = ACES(b);
            break;
        case TonemapType::Reinhard:
            r = Reinhard(r);
            g = Reinhard(g);
            b = Reinhard(b);
            break;
        case TonemapType::Filmic:
            r = Filmic(r);
            g = Filmic(g);
            b = Filmic(b);
            break;
        case TonemapType::AgX:
            r = AgX(r);
            g = AgX(g);
            b = AgX(b);
            break;
        case TonemapType::Exposure:
            r = std::clamp(r, 0.0f, 1.0f);
            g = std::clamp(g, 0.0f, 1.0f);
            b = std::clamp(b, 0.0f, 1.0f);
            break;
        case TonemapType::None:
        default:
            r = std::clamp(r, 0.0f, 1.0f);
            g = std::clamp(g, 0.0f, 1.0f);
            b = std::clamp(b, 0.0f, 1.0f);
            break;
    }
}

} // namespace Tonemap

//==============================================================================
// Bloom Filter
//==============================================================================
class BloomFilter {
public:
    /// @brief Apply bloom effect to RGBA float buffer
    /// @param buffer [in/out] RGBA float buffer (4 floats per pixel)
    /// @param width Image width
    /// @param height Image height
    /// @param config Post-processing configuration
    static void apply(std::vector<float>& buffer, int width, int height, 
                      const PostProcessConfig& config) {
        if (!config.enableBloom || config.bloomIntensity <= 0.0f) {
            return;
        }
        
        // Extract bright pixels
        std::vector<float> bright(buffer.size(), 0.0f);
        for (int y = 0; y < height; ++y) {
            for (int x = 0; x < width; ++x) {
                int idx = (y * width + x) * 4;
                float luminance = 0.2126f * buffer[idx] + 0.7152f * buffer[idx+1] + 0.0722f * buffer[idx+2];
                
                if (luminance > config.bloomThreshold) {
                    float factor = (luminance - config.bloomThreshold) / luminance;
                    bright[idx] = buffer[idx] * factor;
                    bright[idx+1] = buffer[idx+1] * factor;
                    bright[idx+2] = buffer[idx+2] * factor;
                }
            }
        }
        
        // Simple box blur for bloom (separable)
        std::vector<float> blurred(buffer.size(), 0.0f);
        int radius = config.bloomRadius;
        
        // Horizontal pass
        for (int y = 0; y < height; ++y) {
            for (int x = 0; x < width; ++x) {
                float sumR = 0, sumG = 0, sumB = 0;
                int count = 0;
                for (int dx = -radius; dx <= radius; ++dx) {
                    int nx = std::clamp(x + dx, 0, width - 1);
                    int idx = (y * width + nx) * 4;
                    sumR += bright[idx];
                    sumG += bright[idx+1];
                    sumB += bright[idx+2];
                    count++;
                }
                int idx = (y * width + x) * 4;
                blurred[idx] = sumR / count;
                blurred[idx+1] = sumG / count;
                blurred[idx+2] = sumB / count;
            }
        }
        
        // Vertical pass
        bright = blurred;
        for (int y = 0; y < height; ++y) {
            for (int x = 0; x < width; ++x) {
                float sumR = 0, sumG = 0, sumB = 0;
                int count = 0;
                for (int dy = -radius; dy <= radius; ++dy) {
                    int ny = std::clamp(y + dy, 0, height - 1);
                    int idx = (ny * width + x) * 4;
                    sumR += bright[idx];
                    sumG += bright[idx+1];
                    sumB += bright[idx+2];
                    count++;
                }
                int idx = (y * width + x) * 4;
                blurred[idx] = sumR / count;
                blurred[idx+1] = sumG / count;
                blurred[idx+2] = sumB / count;
            }
        }
        
        // Add bloom to original
        for (size_t i = 0; i < buffer.size(); i += 4) {
            buffer[i] += blurred[i] * config.bloomIntensity;
            buffer[i+1] += blurred[i+1] * config.bloomIntensity;
            buffer[i+2] += blurred[i+2] * config.bloomIntensity;
        }
    }
};

//==============================================================================
// Colour Grading
//==============================================================================
class ColourGrading {
public:
    /// @brief Apply colour grading to RGB values
    /// @param r, g, b [in/out] Colour values (linear, 0-1)
    /// @param config Post-processing configuration
    static void apply(float& r, float& g, float& b, const PostProcessConfig& config) {
        // Lift (shadows)
        r += config.lift;
        g += config.lift;
        b += config.lift;
        
        // Gain (highlights)
        r *= config.gain;
        g *= config.gain;
        b *= config.gain;
        
        // Contrast (pivot at 0.5)
        r = (r - 0.5f) * config.contrast + 0.5f;
        g = (g - 0.5f) * config.contrast + 0.5f;
        b = (b - 0.5f) * config.contrast + 0.5f;
        
        // Saturation
        float lum = 0.2126f * r + 0.7152f * g + 0.0722f * b;
        r = lum + (r - lum) * config.saturation;
        g = lum + (g - lum) * config.saturation;
        b = lum + (b - lum) * config.saturation;
        
        // Clamp
        r = std::clamp(r, 0.0f, 1.0f);
        g = std::clamp(g, 0.0f, 1.0f);
        b = std::clamp(b, 0.0f, 1.0f);
    }
};

//==============================================================================
// Post-Processing Pipeline
//==============================================================================
class PostProcessor {
public:
    /// @brief Apply full post-processing pipeline to RGBA buffer
    /// @param buffer [in/out] RGBA float buffer (HDR linear)
    /// @param width Image width
    /// @param height Image height
    /// @param config Configuration
    static void process(std::vector<float>& buffer, int width, int height,
                        const PostProcessConfig& config) {
        // 1. Bloom (operates on HDR values)
        BloomFilter::apply(buffer, width, height, config);
        
        // 2. Tonemapping and colour grading per pixel
        for (int y = 0; y < height; ++y) {
            for (int x = 0; x < width; ++x) {
                int idx = (y * width + x) * 4;
                
                float r = buffer[idx];
                float g = buffer[idx+1];
                float b = buffer[idx+2];
                
                // Tonemapping
                Tonemap::apply(r, g, b, config.tonemapper, config.exposure);
                
                // Colour grading
                ColourGrading::apply(r, g, b, config);
                
                // Gamma correction
                r = std::pow(r, 1.0f / config.gamma);
                g = std::pow(g, 1.0f / config.gamma);
                b = std::pow(b, 1.0f / config.gamma);
                
                buffer[idx] = r;
                buffer[idx+1] = g;
                buffer[idx+2] = b;
            }
        }
    }
};

} // namespace Sirius
