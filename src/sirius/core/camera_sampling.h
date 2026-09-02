#pragma once

#include "sirius/base/contracts.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <numeric>

namespace sirius::core {

struct CameraSample {
    float image_u = 0.5f;
    float image_v = 0.5f;
    float pupil_u = 0.5f;
    float pupil_v = 0.5f;
};

// A separate low-discrepancy dimension for the finite pupil. Using the image
// coordinates again would collapse the four-dimensional film/pupil integral
// onto a two-dimensional diagonal and can bias depth-of-field visibility.
[[nodiscard]] inline float RadicalInverse(std::uint32_t ordinal, std::uint32_t base) noexcept {
    float value = 0.0f;
    float place = 1.0f / static_cast<float>(base);
    while (ordinal != 0) {
        value += static_cast<float>(ordinal % base) * place;
        ordinal /= base;
        place /= static_cast<float>(base);
    }
    return value;
}

// Deterministic stratified pixel sampling. Exactly `samples_per_pixel` samples
// are emitted even when the requested count is not a perfect square.
template <typename Callback>
int ForEachPixelSample(int samples_per_pixel, Callback&& callback) {
    SIRIUS_PRE(samples_per_pixel >= 1 && samples_per_pixel <= 4096);
    const int count = samples_per_pixel;
    const int grid = static_cast<int>(std::ceil(std::sqrt(static_cast<double>(count))));
    if (grid * grid != count) {
        // A truncated square grid clusters the remainder in its first rows.
        // Use a deterministic two-axis Latin hypercube for exact non-square
        // counts: each marginal receives one sample in every stratum.
        int stride = std::max(1, static_cast<int>(0.6180339887498948 * count));
        while (std::gcd(stride, count) != 1) {
            ++stride;
        }
        for (int sample = 0; sample < count; ++sample) {
            const float u = (static_cast<float>(sample) + 0.5f) / static_cast<float>(count);
            const float v =
                (static_cast<float>((static_cast<std::int64_t>(sample) * stride) % count) + 0.5f) /
                static_cast<float>(count);
            callback(u, v);
        }
        return count;
    }

    int emitted = 0;
    for (int y = 0; y < grid && emitted < count; ++y) {
        for (int x = 0; x < grid && emitted < count; ++x) {
            const float u = (static_cast<float>(x) + 0.5f) / static_cast<float>(grid);
            const float v = (static_cast<float>(y) + 0.5f) / static_cast<float>(grid);
            callback(u, v);
            ++emitted;
        }
    }
    return emitted;
}

// Deterministic four-dimensional camera sampling. Film coordinates retain the
// exact-count stratification above; bases 5 and 7 provide distinct pupil
// dimensions without pretending that two film coordinates can also integrate
// a finite aperture.
template <typename Callback>
int ForEachCameraSample(int samples_per_pixel, Callback&& callback) {
    std::uint32_t ordinal = 1;
    return ForEachPixelSample(samples_per_pixel, [&](float image_u, float image_v) {
        callback(CameraSample{
            .image_u = image_u,
            .image_v = image_v,
            .pupil_u = RadicalInverse(ordinal, 5),
            .pupil_v = RadicalInverse(ordinal, 7),
        });
        ++ordinal;
    });
}

}  // namespace sirius::core
