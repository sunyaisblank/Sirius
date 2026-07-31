#pragma once

// Configuration structs the render session and film pipeline consume, ported
// into the render layer to keep this workstream self-contained.
//
// Origins and re-home note:
//   - FilmConfig / FilmFormat / FilmStock  <- CRFM001A.h (destined for
//     app/config/film_config.h)
//   - SiriusConfig and its sub-structs      <- CRCF003A.h (destined for
//     app/config/config_schema.h)
// Both headers are APP-workstream scope. The APP layer must ALIAS these types
// rather than redefine them; the JSON (nlohmann) serialisation macros the
// legacy CRCF003A carried are an app/config concern and are intentionally
// absent here (the render library must not depend on the JSON codec). Field
// names and defaults are preserved exactly so the alias is a no-op.

#include "sirius/base/contracts.h"
#include "sirius/core/constants.h"

#include <cmath>
#include <cstdint>
#include <string>

namespace sirius::render {

// ============================================================================
// Film simulation configuration (origin: CRFM001A.h)
// ============================================================================

// Film format presets (frame geometry).
enum class FilmFormat : uint32_t {
    IMAX70mm_15perf = 0,  // 70mm 15-perf (1.43:1) - Interstellar.
    IMAX70mm_5perf = 1,   // 70mm 5-perf (2.20:1) - standard IMAX.
    Panavision70mm = 2,   // 65mm 5-perf (2.20:1) - 2001: A Space Odyssey.
    VistaVision = 3,      // 35mm 8-perf horizontal (1.66:1).
    Academy35mm = 4,      // 35mm 4-perf (1.37:1) - classic.
    Anamorphic235 = 5,    // 35mm anamorphic (2.35:1).
    Digital = 6           // No film simulation (pure digital).
};

// Film stock presets (grain and colour response).
enum class FilmStock : uint32_t {
    KodakVision3_500T = 0,  // Tungsten balanced, ISO 500 (Interstellar).
    KodakVision3_250D = 1,  // Daylight balanced, ISO 250.
    KodakVision3_50D = 2,   // Fine grain daylight, ISO 50.
    KodakEktachrome = 3,    // Slide film look (high saturation).
    FujiEterna = 4,         // Fuji cinema stock.
    Custom = 5              // User-defined parameters.
};

// IMAX 70mm film-simulation parameters (grain, halation, grade, vignette, bloom).
struct FilmConfig {
    // Format.
    FilmFormat format = FilmFormat::IMAX70mm_15perf;
    float aspect_ratio = 1.43f;
    uint32_t width = 4096;
    uint32_t height = 2864;
    bool enabled = true;

    // Film stock.
    FilmStock stock = FilmStock::KodakVision3_500T;
    float iso = 500.0f;

    // Grain (signal-dependent noise).
    float grain_intensity = 0.02f;
    float grain_size = 1.5f;
    float grain_uniformity = 0.7f;
    bool grain_enabled = true;

    // Halation (light scatter off the film base).
    float halation_radius = 8.0f;
    float halation_strength = 0.15f;
    float halation_threshold = 0.8f;
    float halation_color_r = 1.0f;
    float halation_color_g = 0.5f;
    float halation_color_b = 0.2f;
    bool halation_enabled = true;

    // Colour science.
    float saturation = 1.0f;
    float contrast = 1.0f;
    float exposure = 0.0f;
    float color_temperature_K = 5600.0f;
    float tint = 0.0f;

    // Tone curve (S-curve).
    float toe_strength = 0.5f;
    float shoulder_strength = 0.5f;
    float midtone_point = 0.18f;

    // Gate weave (frame instability).
    float weave_amplitude_x = 0.5f;
    float weave_amplitude_y = 0.3f;
    float weave_frequency = 1.0f;
    bool weave_enabled = false;

    // Vignette.
    float vignette_strength = 0.3f;
    float vignette_radius = 1.2f;
    float vignette_softness = 0.5f;
    bool vignette_enabled = true;

    // Bloom.
    float bloom_intensity = 0.15f;
    float bloom_threshold = 0.9f;
    float bloom_radius = 16.0f;
    bool bloom_enabled = true;

    // Height from width and aspect ratio; kept even for encoder compatibility.
    void ComputeHeight() {
        if (aspect_ratio > 0.0f && width > 0) {
            height = static_cast<uint32_t>(std::round(static_cast<float>(width) / aspect_ratio));
            height = (height / 2) * 2;
        }
    }

    // Apply a format preset (sets aspect_ratio, recomputes height).
    void ApplyFormat(FilmFormat fmt) {
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
                aspect_ratio = 16.0f / 9.0f;
                break;
            default:
                SIRIUS_PRE(false);
                return;
        }
        format = fmt;
        ComputeHeight();
    }

    // Apply a film-stock preset (ISO, grain, colour response).
    void ApplyStock(FilmStock stk) {
        switch (stk) {
            case FilmStock::KodakVision3_500T:
                iso = 500.0f;
                grain_intensity = 0.025f;
                saturation = 0.95f;
                contrast = 1.05f;
                color_temperature_K = 3200.0f;  // Tungsten balanced.
                halation_strength = 0.15f;
                break;
            case FilmStock::KodakVision3_250D:
                iso = 250.0f;
                grain_intensity = 0.015f;
                saturation = 1.0f;
                contrast = 1.0f;
                color_temperature_K = 5600.0f;  // Daylight balanced.
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
                break;  // Keep current values.
            default:
                SIRIUS_PRE(false);
                return;
        }
        stock = stk;
    }

    // Interstellar-style IMAX configuration.
    static FilmConfig Interstellar() {
        FilmConfig cfg;
        cfg.ApplyFormat(FilmFormat::IMAX70mm_15perf);
        cfg.ApplyStock(FilmStock::KodakVision3_500T);
        cfg.width = 4096;
        cfg.ComputeHeight();
        cfg.bloom_enabled = true;
        cfg.bloom_intensity = 0.2f;
        cfg.vignette_enabled = true;
        cfg.vignette_strength = 0.25f;
        return cfg;
    }

    // 2001: A Space Odyssey style.
    static FilmConfig SpaceOdyssey2001() {
        FilmConfig cfg;
        cfg.ApplyFormat(FilmFormat::Panavision70mm);
        cfg.ApplyStock(FilmStock::KodakVision3_50D);  // Fine grain.
        cfg.width = 4096;
        cfg.ComputeHeight();
        cfg.saturation = 0.8f;  // Desaturated look.
        cfg.contrast = 0.9f;
        cfg.grain_intensity = 0.008f;
        return cfg;
    }

    // Clean digital output (no film simulation).
    static FilmConfig DigitalClean() {
        FilmConfig cfg;
        cfg.format = FilmFormat::Digital;
        cfg.aspect_ratio = 16.0f / 9.0f;
        cfg.width = 3840;  // 4K UHD.
        cfg.ComputeHeight();
        cfg.enabled = false;
        cfg.grain_enabled = false;
        cfg.halation_enabled = false;
        cfg.weave_enabled = false;
        cfg.vignette_enabled = false;
        return cfg;
    }
};

// ============================================================================
// Unified render configuration schema (origin: CRCF003A.h, JSON macros dropped)
// ============================================================================

// Render (resolution, sampling, tiling, output path) settings.
struct RenderConfig {
    int width = 1920;
    int height = 1080;
    int samplesPerPixel = 64;
    int tileSize = 64;
    int threadCount = 0;  // 0 = auto-detect.
    std::string outputPath = "render.ppm";
};

// Compile-time identity used after the external string configuration has crossed
// its validation/conversion boundary. Render code must not infer one model from
// an unrecognised string.
enum class DiskTemperatureModel : uint8_t {
    NovikovThorne,
    ShakuraSunyaev,
};

// Metric (spacetime) settings; name is parsed via the core registry.
struct MetricConfig {
    std::string name = "Schwarzschild";
    double mass = 1.0;
    double spin = 0.0;
    double charge = 0.0;
    double lambda = 0.0;  // Cosmological constant (de Sitter family; requires spin = 0).
    std::string temperatureModel = "NovikovThorne";  // NovikovThorne or ShakuraSunyaev.
    float diskTemperature = 50000.0f;                // T_scale (Kelvin).

    // Morris-Thorne wormhole.
    double throatRadius = 1.0;  // b0 in M.
    // The represented chart has one exterior sheet and a dark throat capture
    // surface. TwoSheet remains a named request so it can decline precisely.
    std::string wormholeTopology = "OneSheetCapture";

    // Alcubierre warp drive.
    double warpVelocity = 0.5;  // vs (sub-luminal).
    double bubbleRadius = 1.0;  // R.
    double bubbleSigma = 0.5;   // sigma (wall thickness).
};

// Observer/camera settings.
struct ObserverConfig {
    double distance = 50.0;     // Distance from origin in M.
    double inclination = 90.0;  // Inclination in degrees.
    double azimuth = 0.0;       // Azimuthal angle in degrees.
    double fov = 60.0;          // Field of view in degrees.

    // Camera four-velocity (P5): spatial beta in the local ray-component frame
    // (index 1 = forward/view axis). |beta| < 1. Zero is a static camera and
    // leaves generated rays unaberrated.
    double cameraBetaForward = 0.0;  // beta along the view axis.
    double cameraBetaUp = 0.0;       // beta along the image-up axis.
    double cameraBetaRight = 0.0;    // beta along the image-right axis.

    // Live lens selection (P5). ThinLens uses the deterministic per-sample
    // aperture construction in core/camera.h.
    std::string lensModel = "Pinhole";  // Pinhole or ThinLens.
    float focalLength = 50.0f;          // Millimetres.
    float aperture = 2.8f;              // f-number.
    float focusDistance = 50.0f;        // Geometric units M.
};

// Post-processing settings (display path).
struct PostProcessConfig {
    bool enableBloom = true;
    float bloomIntensity = 0.3f;
    float bloomThreshold = 0.3f;
    float exposure = 1.0f;
    float contrast = 1.0f;
    float saturation = 1.0f;
    std::string tonemapper = "ACES";  // ACES, Reinhard, Uncharted2, Filmic, None, Linear.
};

// Volumetric disk settings.
struct VolumetricConfig {
    bool enabled = false;
    float hOverR = 0.1f;
    float hPower = 0.25f;
    float tauMidplane = 10.0f;
    int samples = 64;
    bool enableTurbulence = false;
    bool enableCorona = false;
};

// Temporal integration of the rotating thin-disk emitter.
struct MotionBlurConfig {
    bool enabled = false;
    float shutterTime = 0.1f;
    int samples = 3;
};

// Film simulation settings (a subset toggled from config; presets fill the rest).
struct FilmSimConfig {
    bool enabled = false;
    std::string preset = "Interstellar";  // Interstellar, SpaceOdyssey2001, DigitalClean.
    float grainIntensity = 0.15f;
    float halationStrength = 0.15f;
    float vignetteStrength = 0.3f;
};

// Backend (renderer) settings.
struct BackendConfig {
    std::string preferred = "auto";  // auto, cpu, vulkan (optix retired).
    bool enableDenoiser = false;
    int cudaDevice = 0;
};

// Root configuration structure.
struct SiriusConfig {
    RenderConfig render;
    MetricConfig metric;
    ObserverConfig observer;
    PostProcessConfig postprocess;
    BackendConfig backend;
    VolumetricConfig volumetric;
    MotionBlurConfig motionBlur;
    FilmSimConfig film;

    // Thin/volumetric accretion disk. Charged and cosmological black holes
    // currently require this false because their circular-orbit/flux model is
    // not represented by the Kerr-only Page-Thorne implementation.
    bool diskEnabled = true;

    // Disk Doppler beaming toggle (P4): true is full physics and the pinned
    // render; false suppresses the disk brightness asymmetry.
    bool dopplerBeaming = true;

    // Filtered point-source star field (P3): true renders stars as beam-filtered
    // points instead of the equirectangular texture.
    bool pointStarfield = false;

    // Ray bundles (P2): propagate geodesic deviation and derive per-pixel beam
    // footprints. Default off keeps the pinned render.
    bool rayBundles = false;

    // Astronomical output representation. Polarisation selects the CPU
    // parallel-transport/Stokes path; diagnostic modes select CPU false-colour
    // paths. The Vulkan capability boundary declines every non-TrueColor mode.
    std::string colorMode = "TrueColor";

    static SiriusConfig defaults() { return SiriusConfig{}; }

    double inclinationRadians() const {
        return observer.inclination * core::constants::math::kPi / 180.0;
    }

    double fovRadians() const { return observer.fov * core::constants::math::kPi / 180.0; }
};

}  // namespace sirius::render
