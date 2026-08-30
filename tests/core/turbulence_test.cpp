#include "sirius/core/disk/turbulence.h"

#include <gtest/gtest.h>

#include <cmath>
#include <limits>

namespace sirius::test {

TEST(TurbulenceTest, MalformedProceduralDomainFailsClosedWithoutRewritingTheRequest) {
    core::TurbulenceConfig config;
    ASSERT_TRUE(core::IsRepresentedTurbulenceConfig(config));
    const float represented_sample =
        core::turbulence_noise::SampleDensityPerturbation(10.0f, 1.0f, 0.5f, config);
    EXPECT_TRUE(std::isfinite(represented_sample));
    EXPECT_GT(represented_sample, 0.0f);

    const float nan = std::numeric_limits<float>::quiet_NaN();
    const float infinity = std::numeric_limits<float>::infinity();
    config.outer_scale_M = infinity;
    config.amplitude = infinity;
    config.lacunarity = nan;
    config.persistence = infinity;
    config.octaves = 0;
    EXPECT_FALSE(core::IsRepresentedTurbulenceConfig(config));
    EXPECT_DEATH((void)core::turbulence_noise::SampleDensityPerturbation(10.0f, 1.0f, 0.5f, config),
                 "precondition.*enforced, terminating");
    EXPECT_TRUE(std::isinf(config.outer_scale_M));
    EXPECT_TRUE(std::isinf(config.amplitude));
    EXPECT_TRUE(std::isnan(config.lacunarity));
    EXPECT_TRUE(std::isinf(config.persistence));
    EXPECT_EQ(config.octaves, 0u);

    config = {};
    EXPECT_DEATH((void)core::turbulence_noise::SampleDensityPerturbation(
                     std::numeric_limits<float>::max(), 1.0f, 0.5f, config),
                 "precondition.*enforced, terminating");
    EXPECT_DEATH((void)core::turbulence_noise::SampleDensityPerturbation(-1.0f, 1.0f, 0.5f, config),
                 "precondition.*enforced, terminating");
}

}  // namespace sirius::test
