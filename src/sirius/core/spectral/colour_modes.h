#pragma once

// Scientific-visualisation colouring strategies for disk emission: physical
// true colour, false-colour temperature and g-factor maps, the Hubble
// narrowband palette, and polarisation visualisation.
// Ported from PHSC001A.h.
// Reference: NASA Hubble Space Telescope imaging guidelines.

#include "sirius/core/polarisation/stokes.h"
#include "sirius/core/spectral/blackbody.h"

#include <algorithm>
#include <cmath>

namespace sirius::core::color_modes {

// Colour mode selection (mirrors SessionConfig::ColorMode).
enum class Mode {
    TrueColor,       // Physical blackbody colours.
    TemperatureMap,  // False-colour temperature visualisation.
    RedshiftMap,     // g-factor visualisation.
    Narrowband,      // Hubble palette (emission-line mapping).
    Polarisation     // Polarisation degree and EVPA.
};

// Maps temperature to a perceptually uniform (modified plasma) colour gradient.
namespace temperature_map {

// Normalised temperature t in [0, 1] (0 = cold, 1 = hot) to linear RGB.
inline spectral::Rgb TemperatureToRgb(float t) {
    t = std::clamp(t, 0.0f, 1.0f);

    // Modified plasma colormap: dark blue -> purple -> orange -> white.
    spectral::Rgb color;

    if (t < 0.25f) {
        // Dark blue to purple.
        float s = t / 0.25f;
        color.r = 0.05f + 0.2f * s;
        color.g = 0.0f + 0.05f * s;
        color.b = 0.2f + 0.4f * s;
    } else if (t < 0.5f) {
        // Purple to red.
        float s = (t - 0.25f) / 0.25f;
        color.r = 0.25f + 0.55f * s;
        color.g = 0.05f;
        color.b = 0.6f - 0.45f * s;
    } else if (t < 0.75f) {
        // Red to orange.
        float s = (t - 0.5f) / 0.25f;
        color.r = 0.8f + 0.15f * s;
        color.g = 0.05f + 0.5f * s;
        color.b = 0.15f - 0.1f * s;
    } else {
        // Orange to yellow/white.
        float s = (t - 0.75f) / 0.25f;
        color.r = 0.95f + 0.05f * s;
        color.g = 0.55f + 0.45f * s;
        color.b = 0.05f + 0.5f * s;
    }

    return color;
}

// Map disk temperature to colour with a logarithmic scale for dynamic range.
inline spectral::Rgb MapTemperature(float T_emit, float T_min = 0.1f, float T_max = 1.5f) {
    float log_T = std::log10(std::max(T_emit, 0.001f));
    float log_min = std::log10(T_min);
    float log_max = std::log10(T_max);
    float t = (log_T - log_min) / (log_max - log_min);
    return TemperatureToRgb(t);
}

}  // namespace temperature_map

// Visualises the g-factor (combined gravitational + Doppler redshift):
// g < 1 red-tinted, g = 1 neutral, g > 1 blue-tinted.
namespace redshift_map {

// Redshift factor g (g < 1 redshifted, g > 1 blueshifted) to linear RGB.
inline spectral::Rgb GfactorToRgb(float g) {
    spectral::Rgb color;

    g = std::clamp(g, 0.1f, 3.0f);

    if (g < 1.0f) {
        // Redshifted: white -> red.
        float t = 1.0f - g;  // 0 at g=1, 0.9 at g=0.1.
        color.r = 1.0f;
        color.g = 1.0f - t * 0.9f;
        color.b = 1.0f - t * 0.95f;
    } else {
        // Blueshifted: white -> blue.
        float t = (g - 1.0f) / 2.0f;  // 0 at g=1, 1 at g=3.
        t = std::min(t, 1.0f);
        color.r = 1.0f - t * 0.9f;
        color.g = 1.0f - t * 0.5f;
        color.b = 1.0f;
    }

    return color;
}

// Map g-factor with intensity modulation (e.g. from T^4).
inline spectral::Rgb MapRedshift(float g, float intensity) {
    spectral::Rgb color = GfactorToRgb(g);
    color.r *= intensity;
    color.g *= intensity;
    color.b *= intensity;
    return color;
}

}  // namespace redshift_map

// Maps emission-line strengths to RGB channels (Hubble SHO palette):
// S-II -> R, H-alpha -> G, O-III -> B. Line strengths approximated from disk
// temperature: hot inner disk high-ionisation (blue), cool outer disk (red).
namespace narrowband {

// Emission-line strengths for a temperature region.
struct EmissionLines {
    float SII = 0.0f;     // Sulfur II (cool regions).
    float Halpha = 0.0f;  // Hydrogen alpha (warm regions).
    float OIII = 0.0f;    // Oxygen III (hot regions).
};

// Approximate emission-line strength from normalised temperature T in [0, 1].
inline EmissionLines EstimateEmission(float T) {
    EmissionLines lines;

    // Temperature-dependent emission strength (simplified Gaussian model; real
    // emission depends on ionisation equilibrium and density).

    // S-II peaks at cooler temperatures (outer disk).
    float T_SII = 0.3f;
    float sigma_SII = 0.2f;
    lines.SII = std::exp(-0.5f * std::pow((T - T_SII) / sigma_SII, 2.0f));

    // H-alpha: broad emission peaking at moderate temperature.
    float T_Ha = 0.5f;
    float sigma_Ha = 0.3f;
    lines.Halpha = std::exp(-0.5f * std::pow((T - T_Ha) / sigma_Ha, 2.0f));

    // O-III peaks at hot temperatures (inner disk).
    float T_OIII = 0.8f;
    float sigma_OIII = 0.25f;
    lines.OIII = std::exp(-0.5f * std::pow((T - T_OIII) / sigma_OIII, 2.0f));

    return lines;
}

// Emission lines to Hubble palette RGB (SII -> R, Halpha -> G, OIII -> B).
inline spectral::Rgb HubblePalette(const EmissionLines& lines, float intensity = 1.0f) {
    spectral::Rgb color;
    color.r = lines.SII * intensity;
    color.g = lines.Halpha * intensity;
    color.b = lines.OIII * intensity;
    return color;
}

// Map normalised temperature to the Hubble palette.
inline spectral::Rgb MapNarrowband(float T, float intensity = 1.0f) {
    EmissionLines lines = EstimateEmission(T);
    return HubblePalette(lines, intensity);
}

}  // namespace narrowband

// Maps polarisation degree and EVPA to colour.
namespace polarisation_vis {

// HSV scheme: hue = EVPA, saturation = degree, value = intensity.
inline spectral::Rgb StokesToRgbHsv(const StokesVector& stokes) {
    float I = stokes.I;
    float p = stokes.PolarisationDegree();  // [0, 1].
    float chi = stokes.Evpa();              // [-pi/2, pi/2].

    // Map EVPA to hue [0, 1].
    float hue = (chi + static_cast<float>(M_PI / 2.0)) / static_cast<float>(M_PI);

    // Polarisation degree as saturation.
    float sat = p;

    // Intensity as value (normalised).
    float val = std::min(I, 1.0f);

    // HSV to RGB conversion.
    float h = hue * 6.0f;
    int i = static_cast<int>(h);
    float f = h - i;
    float q = val * (1.0f - sat);
    float t = val * (1.0f - sat * (1.0f - f));
    float p_val = val * (1.0f - sat * f);

    spectral::Rgb color;
    switch (i % 6) {
        case 0:
            color = {val, t, q};
            break;
        case 1:
            color = {p_val, val, q};
            break;
        case 2:
            color = {q, val, t};
            break;
        case 3:
            color = {q, p_val, val};
            break;
        case 4:
            color = {t, q, val};
            break;
        case 5:
            color = {val, q, p_val};
            break;
    }

    return color;
}

// Diverging false-colour map for EVPA, modulated by degree and intensity.
inline spectral::Rgb StokesToRgbFalseColor(const StokesVector& stokes) {
    float I = stokes.I;
    float p = stokes.PolarisationDegree();
    float chi = stokes.Evpa();

    // Map EVPA to a blue-white-red diverging colormap.
    float t = (chi + static_cast<float>(M_PI / 2.0)) / static_cast<float>(M_PI);  // [0, 1].

    spectral::Rgb color;
    if (t < 0.5f) {
        // Blue to white.
        float s = t / 0.5f;
        color.r = s;
        color.g = s;
        color.b = 1.0f;
    } else {
        // White to red.
        float s = (t - 0.5f) / 0.5f;
        color.r = 1.0f;
        color.g = 1.0f - s;
        color.b = 1.0f - s;
    }

    // Modulate by polarisation degree and intensity.
    float scale = p * std::sqrt(I);
    color.r *= scale;
    color.g *= scale;
    color.b *= scale;

    return color;
}

// EVPA angle (radians) for drawing a polarisation tick.
inline float GetEvpa(const StokesVector& stokes) { return stokes.Evpa(); }

}  // namespace polarisation_vis

// Apply a colour mode to disk emission. T_emit is the emitted temperature
// (normalised), g the redshift factor, intensity the Stefan-Boltzmann term,
// and stokes optional data for the polarisation mode.
inline spectral::Rgb ApplyColorMode(Mode mode, float T_emit, float g, float intensity,
                                    const StokesVector* stokes = nullptr) {
    switch (mode) {
        case Mode::TrueColor: {
            // Physical blackbody -> sRGB.
            float T_obs = T_emit * g;
            constexpr float kTInnerKelvin = 30000.0f;
            float T_kelvin = T_obs * kTInnerKelvin;
            T_kelvin = std::clamp(T_kelvin, 1000.0f, 100000.0f);
            spectral::Rgb color = spectral::BlackbodyToRgb(static_cast<double>(T_kelvin));

            // Apply relativistic beaming.
            float g4 = g * g * g * g;
            color.r *= intensity * g4;
            color.g *= intensity * g4;
            color.b *= intensity * g4;
            return color;
        }

        case Mode::TemperatureMap:
            return temperature_map::MapTemperature(T_emit);

        case Mode::RedshiftMap:
            return redshift_map::MapRedshift(g, intensity);

        case Mode::Narrowband:
            return narrowband::MapNarrowband(T_emit, intensity);

        case Mode::Polarisation:
            if (stokes) {
                return polarisation_vis::StokesToRgbHsv(*stokes);
            } else {
                // Fall back to true colour when no Stokes data is present.
                return ApplyColorMode(Mode::TrueColor, T_emit, g, intensity, nullptr);
            }

        default:
            return spectral::Rgb{1.0f, 0.0f, 1.0f};  // Magenta for unknown mode.
    }
}

}  // namespace sirius::core::color_modes
