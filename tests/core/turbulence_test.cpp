#include "sirius/core/disk/turbulence.h"

#include <gtest/gtest.h>

#include <cmath>
#include <limits>

namespace sirius::test {

TEST(TurbulenceTest, ValidateRestoresFiniteOrderedDomain) {
    core::TurbulenceConfig config;
    const float nan = std::numeric_limits<float>::quiet_NaN();
    const float infinity = std::numeric_limits<float>::infinity();
    config.kolmogorov_exponent = nan;
    config.outer_scale_M = infinity;
    config.inner_scale_M = nan;
    config.amplitude = infinity;
    config.lacunarity = nan;
    config.persistence = infinity;
    config.octaves = 0;
    config.Validate();

    EXPECT_TRUE(std::isfinite(config.kolmogorov_exponent));
    EXPECT_TRUE(std::isfinite(config.outer_scale_M));
    EXPECT_TRUE(std::isfinite(config.inner_scale_M));
    EXPECT_TRUE(std::isfinite(config.amplitude));
    EXPECT_TRUE(std::isfinite(config.lacunarity));
    EXPECT_TRUE(std::isfinite(config.persistence));
    EXPECT_GT(config.outer_scale_M, config.inner_scale_M);
    EXPECT_GE(config.octaves, 1u);
    EXPECT_LE(config.octaves, 8u);

    const float sample =
        core::turbulence_noise::SampleDensityPerturbation(10.0f, 1.0f, 0.5f, config);
    EXPECT_TRUE(std::isfinite(sample));
    EXPECT_GT(sample, 0.0f);
}

}  // namespace sirius::test
