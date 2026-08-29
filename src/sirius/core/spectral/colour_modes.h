#pragma once

// Scientific-visualisation colouring strategies for disk emission: physical
// true colour, false-colour temperature and g-factor maps, and polarisation
// visualisation.
// Ported from PHSC001A.h.
// Reference: NASA Hubble Space Telescope imaging guidelines.

#include "sirius/base/contracts.h"
#include "sirius/core/polarisation/stokes.h"
#include "sirius/core/spectral/blackbody.h"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <span>

namespace sirius::core::color_modes {

// Bolometric Liouville scaling for disk emission. This is the CPU authority for
// the one g^4 factor applied after emission; thin and volumetric paths both call
// it rather than re-deriving the exponent independently.
inline float ObservedBolometricIntensity(float emitted_intensity, float redshift) {
    const float g2 = redshift * redshift;
    return emitted_intensity * g2 * g2;
}

// Colour mode selection (mirrors SessionConfig::ColorMode).
enum class Mode {
    TrueColor,       // Physical blackbody colours.
    TemperatureMap,  // False-colour temperature visualisation.
    RedshiftMap,     // g-factor visualisation.
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

// Maps polarisation degree and EVPA to colour.
namespace polarisation_vis {

// HSV scheme: hue = EVPA, saturation = degree, value = intensity.
inline spectral::Rgb StokesToRgbHsv(const StokesVector& stokes) {
    float I = stokes.I;
    float p = stokes.PolarisationDegree();  // [0, 1].
    float chi = stokes.Evpa();              // [-pi/2, pi/2].

    // Map EVPA to hue [0, 1].
    float hue =
        (chi + static_cast<float>(std::numbers::pi / 2.0)) / static_cast<float>(std::numbers::pi);

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
    float t = (chi + static_cast<float>(std::numbers::pi / 2.0)) /
              static_cast<float>(std::numbers::pi);  // [0, 1].

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
                                    const StokesVector* stokes = nullptr,
                                    float temperature_scale_kelvin = 30000.0f) {
    switch (mode) {
        case Mode::TrueColor: {
            // Physical blackbody -> linear RGB in the sRGB/D65 primaries.
            float T_obs = T_emit * g;
            float T_kelvin = T_obs * temperature_scale_kelvin;
            spectral::Rgb color = spectral::BlackbodyToRgb(static_cast<double>(T_kelvin));

            // Apply relativistic beaming.
            const float observed_intensity = ObservedBolometricIntensity(intensity, g);
            color.r *= observed_intensity;
            color.g *= observed_intensity;
            color.b *= observed_intensity;
            return color;
        }

        case Mode::TemperatureMap:
            return temperature_map::MapTemperature(T_emit);

        case Mode::RedshiftMap:
            return redshift_map::MapRedshift(g, intensity);

        case Mode::Polarisation:
            SIRIUS_PRE(stokes != nullptr);
            return stokes != nullptr ? polarisation_vis::StokesToRgbHsv(*stokes)
                                     : spectral::Rgb{1.0f, 0.0f, 1.0f};

        default:
            SIRIUS_PRE(false);
            return spectral::Rgb{1.0f, 0.0f, 1.0f};
    }
}

// Average radiance/diagnostic colour over temporal redshift samples. Averaging
// g before applying the nonlinear g^4 and blackbody transformations is not a
// temporal radiance integral and systematically biases motion-blurred output.
inline spectral::Rgb AverageTemporalColorMode(Mode mode, float T_emit,
                                              std::span<const float> redshifts,
                                              float emitted_intensity,
                                              float temperature_scale_kelvin = 30000.0f) {
    SIRIUS_PRE(mode != Mode::Polarisation);
    SIRIUS_PRE(!redshifts.empty());
    if (mode == Mode::Polarisation || redshifts.empty()) {
        return spectral::Rgb{1.0f, 0.0f, 1.0f};
    }

    spectral::Rgb total;
    for (const float g : redshifts) {
        const float observed_intensity = ObservedBolometricIntensity(emitted_intensity, g);
        const float mode_intensity =
            mode == Mode::TrueColor ? emitted_intensity : observed_intensity;
        const spectral::Rgb sample =
            ApplyColorMode(mode, T_emit, g, mode_intensity, nullptr, temperature_scale_kelvin);
        total.r += sample.r;
        total.g += sample.g;
        total.b += sample.b;
    }

    const float inverse_count = 1.0f / static_cast<float>(redshifts.size());
    total.r *= inverse_count;
    total.g *= inverse_count;
    total.b *= inverse_count;
    return total;
}

}  // namespace sirius::core::color_modes
