// CRFM001A.h - Film Simulation Configuration
// Component ID: CRFM001A (Core/Configuration/FilmConfig)
//
// IMAX 70mm 15-perf film simulation for cinematic rendering.
// Implements authentic film characteristics including grain, halation,
// and color science based on Kodak Vision3 film stock.
//
// PHYSICAL FOUNDATION:
// ====================
// IMAX 70mm 15-perf format:
//   - Frame size: 70.41mm × 52.63mm (15 perforations wide)
//   - Aspect ratio: 1.43:1 (native IMAX)
//   - Resolution: ~18K equivalent at full scan
//
// Film characteristics modeled:
//   - Grain: Signal-dependent noise following Poisson statistics
//   - Halation: Light scatter from film base reflection
//   - Color response: Sensitometric curves based on Kodak Vision3 500T
//   - Gate weave: Sub-pixel frame registration instability
//
// REFERENCES:
// - Kodak Vision3 500T Technical Data Sheet (5219)
// - IMAX Corporation Technical Specifications (2019)
// - Interstellar Technical Press Kit (2014)

#pragma once

#include <cstdint>
#include <cmath>

namespace Sirius {

//==============================================================================
// Film Format Presets
//==============================================================================
enum class FilmFormat : uint32_t {
    IMAX70mm_15perf = 0,    ///< 70mm 15-perf (1.43:1) - Interstellar
    IMAX70mm_5perf = 1,     ///< 70mm 5-perf (2.20:1) - standard IMAX
    Panavision70mm = 2,     ///< 65mm 5-perf (2.20:1) - 2001: A Space Odyssey
    VistaVision = 3,        ///< 35mm 8-perf horizontal (1.66:1)
    Academy35mm = 4,        ///< 35mm 4-perf (1.37:1) - classic
    Anamorphic235 = 5,      ///< 35mm anamorphic (2.35:1)
    Digital = 6             ///< No film simulation (pure digital)
};

//==============================================================================
// Film Stock Presets
//==============================================================================
enum class FilmStock : uint32_t {
    KodakVision3_500T = 0,  ///< Tungsten balanced, ISO 500 (Interstellar)
    KodakVision3_250D = 1,  ///< Daylight balanced, ISO 250
    KodakVision3_50D = 2,   ///< Fine grain daylight, ISO 50
    KodakEktachrome = 3,    ///< Slide film look (high saturation)
    FujiEterna = 4,         ///< Fuji cinema stock
    Custom = 5              ///< User-defined parameters
};

//==============================================================================
// Film Simulation Configuration
//==============================================================================
struct FilmConfig {
    // =========================================================================
    // Format Settings
    // =========================================================================
    FilmFormat format = FilmFormat::IMAX70mm_15perf;
    float aspect_ratio = 1.43f;         ///< Width / Height (1.43 for IMAX 15-perf)
    uint32_t width = 4096;              ///< Target render width
    uint32_t height = 2864;             ///< Height (computed from aspect if 0)
    bool enabled = true;                ///< Enable film simulation

    // =========================================================================
    // Film Stock Characteristics
    // =========================================================================
    FilmStock stock = FilmStock::KodakVision3_500T;
    float iso = 500.0f;                 ///< Film speed (50-3200)

    // Grain model: σ² = A × L + B × L² where L is luminance
    float grain_intensity = 0.02f;      ///< Base grain strength [0, 1]
    float grain_size = 1.5f;            ///< Grain particle size [pixels]
    float grain_uniformity = 0.7f;      ///< Color channel correlation [0, 1]
    bool grain_enabled = true;

    // =========================================================================
    // Halation (Light Scatter)
    // =========================================================================
    // Light reflecting off film base creates red-orange halo around highlights
    float halation_radius = 8.0f;       ///< Scatter radius [pixels]
    float halation_strength = 0.15f;    ///< Halation intensity [0, 1]
    float halation_threshold = 0.8f;    ///< Brightness threshold for halation
    float halation_color_r = 1.0f;      ///< Red channel contribution
    float halation_color_g = 0.5f;      ///< Green channel (typically lower)
    float halation_color_b = 0.2f;      ///< Blue channel (minimal)
    bool halation_enabled = true;

    // =========================================================================
    // Color Science
    // =========================================================================
    float saturation = 1.0f;            ///< Color saturation multiplier
    float contrast = 1.0f;              ///< Gamma curve steepness
    float exposure = 0.0f;              ///< Exposure compensation [stops]
    float color_temperature_K = 5600.0f; ///< White point temperature
    float tint = 0.0f;                  ///< Green-magenta tint [-1, 1]

    // Tone curve (S-curve for film characteristic)
    float toe_strength = 0.5f;          ///< Shadow compression
    float shoulder_strength = 0.5f;     ///< Highlight rolloff
    float midtone_point = 0.18f;        ///< Middle gray (18% gray card)

    // =========================================================================
    // Gate Weave (Frame Instability)
    // =========================================================================
    float weave_amplitude_x = 0.5f;     ///< Horizontal weave [pixels]
    float weave_amplitude_y = 0.3f;     ///< Vertical weave [pixels]
    float weave_frequency = 1.0f;       ///< Oscillation rate [per frame]
    bool weave_enabled = false;         ///< Off by default for stills

    // =========================================================================
    // Vignette (Corner Darkening)
    // =========================================================================
    float vignette_strength = 0.3f;     ///< Darkening amount [0, 1]
    float vignette_radius = 1.2f;       ///< Falloff radius (1.0 = frame edge)
    float vignette_softness = 0.5f;     ///< Transition smoothness
    bool vignette_enabled = true;

    // =========================================================================
    // Bloom (Lens Glow)
    // =========================================================================
    float bloom_intensity = 0.15f;      ///< Bloom strength [0, 1]
    float bloom_threshold = 0.9f;       ///< Brightness threshold
    float bloom_radius = 16.0f;         ///< Blur radius [pixels]
    bool bloom_enabled = true;

    // =========================================================================
    // Invariants
    // =========================================================================
    // aspect_ratio ∈ [1.0, 3.0]
    // width ∈ [1920, 8192]
    // iso ∈ [25, 3200]
    // grain_intensity ∈ [0, 1]
    // halation_strength ∈ [0, 1]

    /// @brief Compute height from width and aspect ratio
    void computeHeight() {
        if (aspect_ratio > 0.0f && width > 0) {
            height = static_cast<uint32_t>(std::round(static_cast<float>(width) / aspect_ratio));
            // Ensure even dimensions for video encoding compatibility
            height = (height / 2) * 2;
        }
    }

    /// @brief Apply format preset
    void applyFormat(FilmFormat fmt) {
        format = fmt;
        switch (fmt) {
            case FilmFormat::IMAX70mm_15perf:
                aspect_ratio = 1.43f;
                break;
            case FilmFormat::IMAX70mm_5perf:
            case FilmFormat::Panavision70mm:
                aspect_ratio = 2.20f;
                break;
            case FilmFormat::VistaVision:
                aspect_ratio = 1.66f;
                break;
            case FilmFormat::Academy35mm:
                aspect_ratio = 1.37f;
                break;
            case FilmFormat::Anamorphic235:
                aspect_ratio = 2.35f;
                break;
            case FilmFormat::Digital:
                aspect_ratio = 16.0f / 9.0f;  // Default to 16:9
                break;
        }
        computeHeight();
    }

    /// @brief Apply film stock preset
    void applyStock(FilmStock stk) {
        stock = stk;
        switch (stk) {
            case FilmStock::KodakVision3_500T:
                iso = 500.0f;
                grain_intensity = 0.025f;
                saturation = 0.95f;
                contrast = 1.05f;
                color_temperature_K = 3200.0f;  // Tungsten balanced
                halation_strength = 0.15f;
                break;
            case FilmStock::KodakVision3_250D:
                iso = 250.0f;
                grain_intensity = 0.015f;
                saturation = 1.0f;
                contrast = 1.0f;
                color_temperature_K = 5600.0f;  // Daylight balanced
                halation_strength = 0.12f;
                break;
            case FilmStock::KodakVision3_50D:
                iso = 50.0f;
                grain_intensity = 0.005f;
                saturation = 1.0f;
                contrast = 1.02f;
                color_temperature_K = 5600.0f;
                halation_strength = 0.08f;
                break;
            case FilmStock::KodakEktachrome:
                iso = 100.0f;
                grain_intensity = 0.01f;
                saturation = 1.3f;
                contrast = 1.15f;
                color_temperature_K = 5600.0f;
                halation_strength = 0.05f;
                break;
            case FilmStock::FujiEterna:
                iso = 500.0f;
                grain_intensity = 0.02f;
                saturation = 0.9f;
                contrast = 0.95f;
                color_temperature_K = 3200.0f;
                halation_strength = 0.1f;
                break;
            case FilmStock::Custom:
                // Keep current values
                break;
        }
    }

    // =========================================================================
    // Preset Configurations
    // =========================================================================

    /// @brief Interstellar-style IMAX configuration
    static FilmConfig Interstellar() {
        FilmConfig cfg;
        cfg.applyFormat(FilmFormat::IMAX70mm_15perf);
        cfg.applyStock(FilmStock::KodakVision3_500T);
        cfg.width = 4096;
        cfg.computeHeight();
        cfg.bloom_enabled = true;
        cfg.bloom_intensity = 0.2f;
        cfg.vignette_enabled = true;
        cfg.vignette_strength = 0.25f;
        return cfg;
    }

    /// @brief 2001: A Space Odyssey style
    static FilmConfig SpaceOdyssey2001() {
        FilmConfig cfg;
        cfg.applyFormat(FilmFormat::Panavision70mm);
        cfg.applyStock(FilmStock::KodakVision3_50D);  // Fine grain
        cfg.width = 4096;
        cfg.computeHeight();
        cfg.saturation = 0.8f;  // Desaturated look
        cfg.contrast = 0.9f;
        cfg.grain_intensity = 0.008f;
        return cfg;
    }

    /// @brief Clean digital output (no film simulation)
    static FilmConfig DigitalClean() {
        FilmConfig cfg;
        cfg.format = FilmFormat::Digital;
        cfg.aspect_ratio = 16.0f / 9.0f;
        cfg.width = 3840;  // 4K UHD
        cfg.computeHeight();
        cfg.enabled = false;
        cfg.grain_enabled = false;
        cfg.halation_enabled = false;
        cfg.weave_enabled = false;
        cfg.vignette_enabled = false;
        return cfg;
    }
};

} // namespace Sirius
