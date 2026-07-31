#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <numeric>

namespace sirius::render {

// Deterministic stratified pixel sampling. Exactly `samples_per_pixel` samples
// are emitted even when the requested count is not a perfect square.
template <typename Callback>
int ForEachPixelSample(int samples_per_pixel, Callback&& callback) {
    const int count = std::max(1, samples_per_pixel);
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

}  // namespace sirius::render
