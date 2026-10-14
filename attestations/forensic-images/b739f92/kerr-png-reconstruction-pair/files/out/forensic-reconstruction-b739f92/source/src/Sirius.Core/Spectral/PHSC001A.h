// PHSC001A.h - Astronomical Spectral Coloring Modes
// Component ID: PHSC001A (Physics/Spectral/ColorModes)
//
// Implements multiple coloring strategies for scientific visualization,
// inspired by the Hubble Space Telescope imaging pipeline.
//
// COLOR MODES:
// ============
// 1. TrueColor:      Physical blackbody → XYZ → sRGB (human vision simulation)
// 2. TemperatureMap: False color temperature gradient (cold→blue, hot→red)
// 3. RedshiftMap:    g-factor visualization (redshift→red, blueshift→blue)
// 4. Narrowband:     Hubble palette mapping emission bands to RGB channels
// 5. Polarisation:   Polarisation degree and EVPA as false color
//
// SCIENTIFIC BACKGROUND:
// ======================
// Real astronomical images are often false color representations that map
// different data channels to RGB for visualization:
// - Broadband filters (visible light) → True color reconstruction
// - Narrowband filters (emission lines) → Hubble palette (SII→R, Hα→G, OIII→B)
// - Infrared/UV data → Mapped to visible spectrum
//
// Reference: "How scientists colorize photos of space" (Vox)
//            NASA Hubble Space Telescope imaging guidelines
//
// TESTS: TSSP002A.cpp

#pragma once

#include "PHSP001A.h"
#include "PHPL001A.h"
#include <cmath>
#include <algorithm>

namespace Sirius::ColorModes {

//==============================================================================
// Color Mode Enumeration (mirrors SessionConfig::ColorMode)
//==============================================================================
enum class Mode {
    TrueColor,       ///< Physical blackbody colors
    TemperatureMap,  ///< False color temperature visualization
    RedshiftMap,     ///< g-factor visualization
    Narrowband,      ///< Hubble palette (emission line mapping)
    Polarisation     ///< Polarisation degree and EVPA
};

//==============================================================================
// Temperature Map Coloring
//==============================================================================
// Maps temperature to a perceptually uniform color gradient.
// Uses a modified plasma/inferno colormap for scientific visualization.
//==============================================================================
namespace TemperatureMap {

/// @brief Convert normalized temperature (0-1) to RGB
/// @param t Normalized temperature [0, 1] where 0=cold, 1=hot
/// @return RGB color in linear space
inline Spectral::RGB temperatureToRGB(float t) {
    t = std::clamp(t, 0.0f, 1.0f);

    // Modified plasma colormap (perceptually uniform)
    // Cold (0) → Dark blue → Purple → Orange → Hot (1) → White
    Spectral::RGB color;

    if (t < 0.25f) {
        // Dark blue to purple
        float s = t / 0.25f;
        color.r = 0.05f + 0.2f * s;
        color.g = 0.0f + 0.05f * s;
        color.b = 0.2f + 0.4f * s;
    } else if (t < 0.5f) {
        // Purple to red
        float s = (t - 0.25f) / 0.25f;
        color.r = 0.25f + 0.55f * s;
        color.g = 0.05f;
        color.b = 0.6f - 0.45f * s;
    } else if (t < 0.75f) {
        // Red to orange
        float s = (t - 0.5f) / 0.25f;
        color.r = 0.8f + 0.15f * s;
        color.g = 0.05f + 0.5f * s;
        color.b = 0.15f - 0.1f * s;
    } else {
        // Orange to yellow/white
        float s = (t - 0.75f) / 0.25f;
        color.r = 0.95f + 0.05f * s;
        color.g = 0.55f + 0.45f * s;
        color.b = 0.05f + 0.5f * s;
    }

    return color;
}

/// @brief Map disk temperature to color
/// @param T_emit Emitted temperature (normalized units, 0-1 typical range)
/// @param T_min Minimum temperature for mapping
/// @param T_max Maximum temperature for mapping
/// @return RGB color in linear space
inline Spectral::RGB mapTemperature(float T_emit, float T_min = 0.1f, float T_max = 1.5f) {
    // Logarithmic mapping for better dynamic range
    float log_T = std::log10(std::max(T_emit, 0.001f));
    float log_min = std::log10(T_min);
    float log_max = std::log10(T_max);
    float t = (log_T - log_min) / (log_max - log_min);
    return temperatureToRGB(t);
}

} // namespace TemperatureMap

//==============================================================================
// Redshift Map Coloring
//==============================================================================
// Visualizes the g-factor (combined gravitational + Doppler redshift).
// g < 1: Redshifted → Red tint
// g = 1: No shift → White/neutral
// g > 1: Blueshifted → Blue tint
//==============================================================================
namespace RedshiftMap {

/// @brief Convert g-factor to RGB color
/// @param g Redshift factor (g < 1 = redshifted, g > 1 = blueshifted)
/// @return RGB color in linear space
inline Spectral::RGB gfactorToRGB(float g) {
    Spectral::RGB color;

    // Clamp to reasonable range
    g = std::clamp(g, 0.1f, 3.0f);

    if (g < 1.0f) {
        // Redshifted: white → red
        float t = 1.0f - g;  // 0 at g=1, 0.9 at g=0.1
        color.r = 1.0f;
        color.g = 1.0f - t * 0.9f;
        color.b = 1.0f - t * 0.95f;
    } else {
        // Blueshifted: white → blue
        float t = (g - 1.0f) / 2.0f;  // 0 at g=1, 1 at g=3
        t = std::min(t, 1.0f);
        color.r = 1.0f - t * 0.9f;
        color.g = 1.0f - t * 0.5f;
        color.b = 1.0f;
    }

    return color;
}

/// @brief Map g-factor with intensity modulation
/// @param g Redshift factor
/// @param intensity Base intensity (e.g., from T^4)
/// @return RGB color with intensity applied
inline Spectral::RGB mapRedshift(float g, float intensity) {
    Spectral::RGB color = gfactorToRGB(g);
    color.r *= intensity;
    color.g *= intensity;
    color.b *= intensity;
    return color;
}

} // namespace RedshiftMap

//==============================================================================
// Narrowband (Hubble Palette) Coloring
//==============================================================================
// Maps emission line strengths to RGB channels.
//
// Hubble's classic palette (SHO):
//   S-II (672nm) → Red channel (ionized sulfur)
//   H-α (656nm) → Green channel (hydrogen alpha)
//   O-III (501nm) → Blue channel (doubly ionized oxygen)
//
// For accretion disks, we can approximate emission line strength from
// temperature regions:
//   Hot inner disk → High ionization (O-III analog) → Blue
//   Mid disk → Hydrogen emission → Green
//   Cool outer disk → Low ionization (S-II analog) → Red
//==============================================================================
namespace Narrowband {

/// @brief Simulate narrowband emission based on temperature
/// This approximates emission line strength from disk temperature.
/// @param T Temperature in normalized units (0-1)
/// @return Emission strengths for SII, Hα, OIII
struct EmissionLines {
    float SII = 0.0f;   ///< Sulfur II (cool regions)
    float Halpha = 0.0f; ///< Hydrogen alpha (warm regions)
    float OIII = 0.0f;  ///< Oxygen III (hot regions)
};

inline EmissionLines estimateEmission(float T) {
    EmissionLines lines;

    // Temperature-dependent emission strength (simplified model)
    // Real emission depends on ionization equilibrium and density

    // S-II: peaks at cooler temperatures (outer disk)
    // Modeled as Gaussian centered at T=0.3
    float T_SII = 0.3f;
    float sigma_SII = 0.2f;
    lines.SII = std::exp(-0.5f * std::pow((T - T_SII) / sigma_SII, 2.0f));

    // H-alpha: broad emission across temperature range
    // Peaks at moderate temperature
    float T_Ha = 0.5f;
    float sigma_Ha = 0.3f;
    lines.Halpha = std::exp(-0.5f * std::pow((T - T_Ha) / sigma_Ha, 2.0f));

    // O-III: peaks at hot temperatures (inner disk)
    float T_OIII = 0.8f;
    float sigma_OIII = 0.25f;
    lines.OIII = std::exp(-0.5f * std::pow((T - T_OIII) / sigma_OIII, 2.0f));

    return lines;
}

/// @brief Convert emission lines to Hubble palette RGB
/// @param lines Emission line strengths
/// @param intensity Base intensity multiplier
/// @return RGB color (SII→R, Hα→G, OIII→B)
inline Spectral::RGB hubblePalette(const EmissionLines& lines, float intensity = 1.0f) {
    Spectral::RGB color;
    color.r = lines.SII * intensity;
    color.g = lines.Halpha * intensity;
    color.b = lines.OIII * intensity;
    return color;
}

/// @brief Map temperature to Hubble palette
/// @param T Temperature in normalized units
/// @param intensity Base intensity
/// @return RGB in Hubble palette
inline Spectral::RGB mapNarrowband(float T, float intensity = 1.0f) {
    EmissionLines lines = estimateEmission(T);
    return hubblePalette(lines, intensity);
}

} // namespace Narrowband

//==============================================================================
// Polarisation Visualization
//==============================================================================
// Maps polarisation degree and EVPA to color.
// Multiple visualization schemes:
//   1. Degree as brightness, EVPA as hue (HSV-based)
//   2. False color with degree saturation
//==============================================================================
namespace PolarisationVis {

/// @brief Convert polarisation to RGB using HSV scheme
/// @param stokes Stokes vector (I, Q, U, V)
/// @return RGB color where hue=EVPA, saturation=degree, value=intensity
inline Spectral::RGB stokesToRGB_HSV(const StokesVector& stokes) {
    float I = stokes.I;
    float p = stokes.polarisationDegree();  // [0, 1]
    float chi = stokes.EVPA();  // [-π/2, π/2]

    // Map EVPA to hue [0, 1]
    float hue = (chi + static_cast<float>(M_PI / 2.0)) / static_cast<float>(M_PI);

    // Polarisation degree as saturation
    float sat = p;

    // Intensity as value (normalized)
    float val = std::min(I, 1.0f);

    // HSV to RGB conversion
    float h = hue * 6.0f;
    int i = static_cast<int>(h);
    float f = h - i;
    float q = val * (1.0f - sat);
    float t = val * (1.0f - sat * (1.0f - f));
    float p_val = val * (1.0f - sat * f);

    Spectral::RGB color;
    switch (i % 6) {
        case 0: color = {val, t, q}; break;
        case 1: color = {p_val, val, q}; break;
        case 2: color = {q, val, t}; break;
        case 3: color = {q, p_val, val}; break;
        case 4: color = {t, q, val}; break;
        case 5: color = {val, q, p_val}; break;
    }

    return color;
}

/// @brief Create false-color polarisation map
/// Uses a diverging colormap for EVPA and intensity modulation for degree
/// @param stokes Stokes vector
/// @return RGB color
inline Spectral::RGB stokesToRGB_FalseColor(const StokesVector& stokes) {
    float I = stokes.I;
    float p = stokes.polarisationDegree();
    float chi = stokes.EVPA();

    // Map EVPA to blue-white-red diverging colormap
    float t = (chi + static_cast<float>(M_PI / 2.0)) / static_cast<float>(M_PI);  // [0, 1]

    Spectral::RGB color;
    if (t < 0.5f) {
        // Blue to white
        float s = t / 0.5f;
        color.r = s;
        color.g = s;
        color.b = 1.0f;
    } else {
        // White to red
        float s = (t - 0.5f) / 0.5f;
        color.r = 1.0f;
        color.g = 1.0f - s;
        color.b = 1.0f - s;
    }

    // Modulate by polarisation degree and intensity
    float scale = p * std::sqrt(I);
    color.r *= scale;
    color.g *= scale;
    color.b *= scale;

    return color;
}

/// @brief Draw polarisation vectors (returns line orientation angle)
/// @param stokes Stokes vector
/// @return EVPA angle in radians for drawing polarisation tick
inline float getEVPA(const StokesVector& stokes) {
    return stokes.EVPA();
}

} // namespace PolarisationVis

//==============================================================================
// Unified Color Mapping Interface
//==============================================================================

/// @brief Apply color mode to disk emission
/// @param mode Color mode to use
/// @param T_emit Emitted temperature (normalized)
/// @param g Redshift g-factor
/// @param intensity Stefan-Boltzmann intensity (T^4)
/// @param stokes Optional Stokes vector for polarisation mode
/// @return RGB color in linear space
inline Spectral::RGB applyColorMode(
    Mode mode,
    float T_emit,
    float g,
    float intensity,
    const StokesVector* stokes = nullptr)
{
    switch (mode) {
        case Mode::TrueColor: {
            // Physical blackbody → sRGB
            float T_obs = T_emit * g;
            constexpr float T_INNER_KELVIN = 30000.0f;
            float T_kelvin = T_obs * T_INNER_KELVIN;
            T_kelvin = std::clamp(T_kelvin, 1000.0f, 100000.0f);
            Spectral::RGB color = Spectral::blackbodyToRGB(static_cast<double>(T_kelvin));

            // Apply relativistic beaming
            float g4 = g * g * g * g;
            color.r *= intensity * g4;
            color.g *= intensity * g4;
            color.b *= intensity * g4;
            return color;
        }

        case Mode::TemperatureMap:
            return TemperatureMap::mapTemperature(T_emit);

        case Mode::RedshiftMap:
            return RedshiftMap::mapRedshift(g, intensity);

        case Mode::Narrowband:
            return Narrowband::mapNarrowband(T_emit, intensity);

        case Mode::Polarisation:
            if (stokes) {
                return PolarisationVis::stokesToRGB_HSV(*stokes);
            } else {
                // Fallback to true color if no Stokes data
                return applyColorMode(Mode::TrueColor, T_emit, g, intensity, nullptr);
            }

        default:
            return Spectral::RGB{1.0f, 0.0f, 1.0f};  // Magenta for unknown mode
    }
}

} // namespace Sirius::ColorModes
