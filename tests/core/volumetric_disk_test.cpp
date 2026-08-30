#include "sirius/core/disk/volumetric_disk.h"

#include <gtest/gtest.h>

#include <cmath>
#include <numbers>

namespace {

using sirius::core::volumetric_disk::kTruncatedGaussianColumnFraction;
using sirius::core::volumetric_disk::kTruncatedGaussianColumnIntegral;
using sirius::core::volumetric_disk::kVerticalTruncationSigma;
using sirius::core::volumetric_disk::TruncatedGaussianOpacityDensity;

TEST(VolumetricDiskClosure, TruncatedGaussianColumnEqualsDeclaredOpticalDepth) {
    constexpr double tau_ref = 7.0;
    constexpr double reference_radius = 6.0;
    constexpr double radius = 15.0;
    constexpr double scale_height = 1.75;
    constexpr int panels = 16384;  // Even: composite Simpson quadrature.
    constexpr double lower = -kVerticalTruncationSigma * scale_height;
    constexpr double upper = kVerticalTruncationSigma * scale_height;
    constexpr double step = (upper - lower) / panels;

    double weighted_sum =
        TruncatedGaussianOpacityDensity(tau_ref, radius, reference_radius, scale_height, lower) +
        TruncatedGaussianOpacityDensity(tau_ref, radius, reference_radius, scale_height, upper);
    for (int panel = 1; panel < panels; ++panel) {
        const double height = lower + panel * step;
        const double weight = panel % 2 == 0 ? 2.0 : 4.0;
        weighted_sum += weight * TruncatedGaussianOpacityDensity(tau_ref, radius, reference_radius,
                                                                 scale_height, height);
    }
    const double integrated_column = weighted_sum * step / 3.0;
    const double declared_column = tau_ref * std::pow(radius / reference_radius, -1.5);

    EXPECT_NEAR(kTruncatedGaussianColumnFraction,
                std::erf(kVerticalTruncationSigma / std::numbers::sqrt2), 1.0e-16);
    EXPECT_NEAR(kTruncatedGaussianColumnIntegral,
                std::sqrt(2.0 * std::numbers::pi) * kTruncatedGaussianColumnFraction, 5.0e-16);
    EXPECT_NEAR(integrated_column, declared_column, 2.0e-12 * declared_column);
    const double former_infinite_support_normalisation =
        declared_column * kTruncatedGaussianColumnFraction;
    EXPECT_LT(former_infinite_support_normalisation, declared_column);
    EXPECT_NEAR(former_infinite_support_normalisation / declared_column,
                kTruncatedGaussianColumnFraction, 1.0e-16);
    EXPECT_DOUBLE_EQ(
        TruncatedGaussianOpacityDensity(tau_ref, radius, reference_radius, scale_height,
                                        (kVerticalTruncationSigma + 1.0e-6) * scale_height),
        0.0);
}

}  // namespace
