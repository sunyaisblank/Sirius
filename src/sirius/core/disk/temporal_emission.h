#pragma once

// Exact temporal sampling of the primary thin-disk Doppler factor. The render
// session and its convergence/cardinality gate share this authority so an
// operator-requested sample count cannot be silently reduced or reinterpreted.

#include "sirius/base/contracts.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

namespace sirius::core::disk_emission {

inline std::vector<float> SampleTemporalRedshifts(float gravitational_factor, float gamma,
                                                  float orbital_velocity, float cosine_coefficient,
                                                  float sine_coefficient, float delta_phi_max,
                                                  int sample_count) {
    SIRIUS_PRE(std::isfinite(gravitational_factor));
    SIRIUS_PRE(std::isfinite(gamma) && gamma > 0.0f);
    SIRIUS_PRE(std::isfinite(orbital_velocity));
    SIRIUS_PRE(std::isfinite(cosine_coefficient));
    SIRIUS_PRE(std::isfinite(sine_coefficient));
    SIRIUS_PRE(std::isfinite(delta_phi_max));
    SIRIUS_PRE(sample_count >= 1 && sample_count <= 4096);

    constexpr float kDopplerDenominatorMin = 0.1f;
    constexpr float kDopplerDenominatorMax = 10.0f;
    constexpr float kGFactorMin = 0.1f;
    constexpr float kGFactorMax = 5.0f;

    std::vector<float> redshifts;
    redshifts.reserve(static_cast<std::size_t>(sample_count));
    for (int sample = 0; sample < sample_count; ++sample) {
        const float time =
            sample_count > 1 ? static_cast<float>(sample) / (sample_count - 1) : 0.5f;
        const float delta_phi = delta_phi_max * (time - 0.5f);
        const float projected_velocity =
            orbital_velocity *
            (cosine_coefficient * std::cos(delta_phi) + sine_coefficient * std::sin(delta_phi));
        const float denominator = std::clamp(gamma * (1.0f - projected_velocity),
                                             kDopplerDenominatorMin, kDopplerDenominatorMax);
        redshifts.push_back(
            std::clamp(gravitational_factor / denominator, kGFactorMin, kGFactorMax));
    }
    SIRIUS_POST(redshifts.size() == static_cast<std::size_t>(sample_count));
    return redshifts;
}

}  // namespace sirius::core::disk_emission
