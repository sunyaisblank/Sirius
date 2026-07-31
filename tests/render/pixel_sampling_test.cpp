#include "sirius/render/pixel_sampling.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <utility>
#include <vector>

namespace sirius::render::test {

TEST(PixelSampling, EmitsExactlyEveryRequestedNonSquareCount) {
    for (const int requested : {1, 2, 3, 5, 7, 10, 63, 65}) {
        std::vector<std::pair<float, float>> samples;
        const int emitted =
            ForEachPixelSample(requested, [&](float u, float v) { samples.emplace_back(u, v); });
        EXPECT_EQ(emitted, requested);
        ASSERT_EQ(samples.size(), static_cast<std::size_t>(requested));
        for (const auto& [u, v] : samples) {
            EXPECT_GT(u, 0.0f);
            EXPECT_LT(u, 1.0f);
            EXPECT_GT(v, 0.0f);
            EXPECT_LT(v, 1.0f);
        }
    }
}

TEST(PixelSampling, NonPositiveInputStillEmitsOneDefensiveSample) {
    int callbacks = 0;
    const int emitted = ForEachPixelSample(0, [&](float, float) { ++callbacks; });
    EXPECT_EQ(emitted, 1);
    EXPECT_EQ(callbacks, 1);
}

TEST(PixelSampling, NonSquareCountsCoverBothAxesWithoutRemainderBias) {
    for (const int requested : {2, 3, 5, 7, 10, 63, 65}) {
        std::vector<std::pair<float, float>> samples;
        ForEachPixelSample(requested, [&](float u, float v) { samples.emplace_back(u, v); });

        double u_sum = 0.0;
        double v_sum = 0.0;
        std::vector<int> u_strata;
        std::vector<int> v_strata;
        for (const auto& [u, v] : samples) {
            u_sum += u;
            v_sum += v;
            u_strata.push_back(
                std::min(requested - 1, static_cast<int>(std::floor(u * requested))));
            v_strata.push_back(
                std::min(requested - 1, static_cast<int>(std::floor(v * requested))));
        }
        std::sort(u_strata.begin(), u_strata.end());
        std::sort(v_strata.begin(), v_strata.end());

        EXPECT_NEAR(u_sum / requested, 0.5, 1e-6);
        EXPECT_NEAR(v_sum / requested, 0.5, 1e-6);
        for (int stratum = 0; stratum < requested; ++stratum) {
            EXPECT_EQ(u_strata[static_cast<std::size_t>(stratum)], stratum);
            EXPECT_EQ(v_strata[static_cast<std::size_t>(stratum)], stratum);
        }
    }
}

TEST(PixelSampling, PatternIsDeterministic) {
    std::vector<std::pair<float, float>> first;
    std::vector<std::pair<float, float>> second;
    ForEachPixelSample(5, [&](float u, float v) { first.emplace_back(u, v); });
    ForEachPixelSample(5, [&](float u, float v) { second.emplace_back(u, v); });
    EXPECT_EQ(first, second);
}

}  // namespace sirius::render::test
