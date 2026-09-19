#pragma once

#include <cmath>
#include <cstdint>

namespace sirius::core {

// Exact request consumed by the session's indexed DNGR point catalogue. The
// broader StarfieldConfig also carries magnitude-cull, parallax, and synthetic
// depth-of-field controls used by other sampling APIs; exposing them here would
// accept controls that this path cannot consume.
struct PointStarfieldConfig {
    friend bool operator==(const PointStarfieldConfig&, const PointStarfieldConfig&) = default;

    std::uint32_t star_count = 100000;
    float min_distance_pc = 1.0f;
    float max_distance_pc = 10000.0f;
    float brightness_scale = 100.0f;
    std::uint32_t seed = 42;
};

[[nodiscard]] inline bool IsRepresentedPointStarfieldConfig(
    const PointStarfieldConfig& config) noexcept {
    return config.star_count >= 1 && config.star_count <= 10000000u &&
           std::isfinite(config.min_distance_pc) && config.min_distance_pc >= 0.1f &&
           std::isfinite(config.max_distance_pc) &&
           config.max_distance_pc > config.min_distance_pc &&
           std::isfinite(config.brightness_scale) && config.brightness_scale >= 0.0f &&
           config.brightness_scale <= 1000000.0f;
}

}  // namespace sirius::core
