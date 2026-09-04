#include "sirius/core/camera_sampling.h"

#include <gtest/gtest.h>

namespace sirius::core::test {

TEST(CameraSampling, SequenceDoesNotCollapseFilmAndPupilDimensions) {
    constexpr int kSamples = 4096;
    double joint_square_difference = 0.0;
    int observed = 0;
    const int emitted = ForEachCameraSample(kSamples, [&](const CameraSample& sample) {
        for (const float coordinate :
             {sample.image_u, sample.image_v, sample.pupil_u, sample.pupil_v}) {
            EXPECT_GT(coordinate, 0.0f);
            EXPECT_LT(coordinate, 1.0f);
        }
        const double du = static_cast<double>(sample.image_u - sample.pupil_u);
        const double dv = static_cast<double>(sample.image_v - sample.pupil_v);
        joint_square_difference += du * du + dv * dv;
        ++observed;
    });

    ASSERT_EQ(emitted, kSamples);
    ASSERT_EQ(observed, kSamples);
    // For four independent unit-uniform dimensions,
    // E[(film_u-pupil_u)^2 + (film_v-pupil_v)^2] = 1/3. Reusing the two film
    // coordinates for the pupil is an adversarial negative control that gives
    // exactly zero and therefore cannot pass this bounded quadrature oracle.
    EXPECT_NEAR(joint_square_difference / kSamples, 1.0 / 3.0, 0.01);
}

}  // namespace sirius::core::test
