#pragma once

// Typed film-pipeline configuration. Operator-facing JSON configuration belongs
// to sirius::app and is projected into this render-owned type at the app/render
// boundary; the render layer therefore has no dependency on the JSON codec.

#include <cmath>

namespace sirius::render {

// Bounded display finishing parameters. These are deliberately named for the
// operations the implementation actually performs; stock sensitivity, colour
// temperature, gate weave, and film-gate geometry are absent until represented.
struct FilmConfig {
    friend bool operator==(const FilmConfig&, const FilmConfig&) = default;

    // Grain (signal-dependent noise).
    float grain_intensity = 0.02f;
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

    // Tone curve (S-curve).
    float toe_strength = 0.5f;
    float shoulder_strength = 0.5f;
    float midtone_point = 0.18f;

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

    // Interstellar-inspired display finish. This names a look, not a stock or
    // film-gate simulation.
    static FilmConfig Interstellar() {
        FilmConfig cfg;
        cfg.grain_intensity = 0.025f;
        cfg.saturation = 0.95f;
        cfg.contrast = 1.05f;
        cfg.halation_strength = 0.15f;
        cfg.bloom_enabled = true;
        cfg.bloom_intensity = 0.2f;
        cfg.vignette_enabled = true;
        cfg.vignette_strength = 0.25f;
        return cfg;
    }

    // 2001-inspired lower-grain, desaturated display finish.
    static FilmConfig SpaceOdyssey2001() {
        FilmConfig cfg;
        cfg.grain_intensity = 0.008f;
        cfg.halation_strength = 0.08f;
        cfg.saturation = 0.8f;
        cfg.contrast = 0.9f;
        return cfg;
    }
};

[[nodiscard]] inline bool IsRepresentedFilmConfig(const FilmConfig& config) noexcept {
    const auto in_range = [](float value, float minimum, float maximum) {
        return std::isfinite(value) && value >= minimum && value <= maximum;
    };
    if (!in_range(config.grain_intensity, 0.0f, 1.0f) ||
        !in_range(config.grain_uniformity, 0.0f, 1.0f) ||
        !in_range(config.halation_radius, 0.0f, 256.0f) ||
        !in_range(config.halation_strength, 0.0f, 5.0f) ||
        !in_range(config.halation_threshold, 0.0f, 100.0f) ||
        !in_range(config.halation_color_r, 0.0f, 10.0f) ||
        !in_range(config.halation_color_g, 0.0f, 10.0f) ||
        !in_range(config.halation_color_b, 0.0f, 10.0f) ||
        !in_range(config.saturation, 0.0f, 4.0f) || !in_range(config.contrast, 0.0f, 4.0f) ||
        !in_range(config.exposure, -100.0f, 100.0f) ||
        !in_range(config.toe_strength, 0.0f, 10.0f) ||
        !in_range(config.shoulder_strength, 0.0f, 10.0f) ||
        !in_range(config.midtone_point, 0.01f, 0.99f) ||
        !in_range(config.vignette_strength, 0.0f, 2.0f) ||
        !in_range(config.vignette_radius, 0.0f, 10.0f) ||
        !in_range(config.vignette_softness, 0.0f, 10.0f) ||
        !in_range(config.bloom_intensity, 0.0f, 5.0f) ||
        !in_range(config.bloom_threshold, 0.0f, 100.0f) ||
        !in_range(config.bloom_radius, 0.0f, 256.0f)) {
        return false;
    }
    if (!config.grain_enabled &&
        (config.grain_intensity != 0.0f || config.grain_uniformity != 0.0f)) {
        return false;
    }
    if (!config.halation_enabled &&
        (config.halation_radius != 0.0f || config.halation_strength != 0.0f ||
         config.halation_threshold != 0.0f || config.halation_color_r != 0.0f ||
         config.halation_color_g != 0.0f || config.halation_color_b != 0.0f)) {
        return false;
    }
    if (!config.vignette_enabled &&
        (config.vignette_strength != 0.0f || config.vignette_radius != 0.0f ||
         config.vignette_softness != 0.0f)) {
        return false;
    }
    if (config.vignette_enabled && config.vignette_softness <= 0.0f) return false;
    return config.bloom_enabled || (config.bloom_intensity == 0.0f &&
                                    config.bloom_threshold == 0.0f && config.bloom_radius == 0.0f);
}

}  // namespace sirius::render
