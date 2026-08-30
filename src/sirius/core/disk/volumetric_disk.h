#pragma once

// Declared stationary grey closure for the phenomenological volumetric disk.
// The atmosphere has Gaussian vertical shape but finite support |z| <= 3 H.
// Therefore tau_ref denotes the column through that represented finite support,
// not through an unrepresented Gaussian tail extending to infinity.

#include <cmath>

namespace sirius::core::volumetric_disk {

inline constexpr double kVerticalTruncationSigma = 3.0;
inline constexpr double kTruncatedGaussianColumnFraction = 0.99730020393673979;
inline constexpr double kTruncatedGaussianColumnIntegral = 2.4998608894830947;

// Extinction coefficient kappa*rho at cylindrical (radius, height). The radial
// target column is tau_ref * (radius/reference_radius)^(-3/2). Returning zero
// for malformed or out-of-support samples keeps this bounded closure fail-safe.
[[nodiscard]] inline double TruncatedGaussianOpacityDensity(double tau_ref, double radius,
                                                            double reference_radius,
                                                            double scale_height,
                                                            double height) noexcept {
    if (!std::isfinite(tau_ref) || tau_ref < 0.0 || !std::isfinite(radius) || !(radius > 0.0) ||
        !std::isfinite(reference_radius) || !(reference_radius > 0.0) ||
        !std::isfinite(scale_height) || !(scale_height > 0.0) || !std::isfinite(height)) {
        return 0.0;
    }
    const double z_over_h = height / scale_height;
    if (!std::isfinite(z_over_h) || std::abs(z_over_h) > kVerticalTruncationSigma) {
        return 0.0;
    }
    const double radial_column = tau_ref * std::pow(radius / reference_radius, -1.5);
    const double density = radial_column * std::exp(-0.5 * z_over_h * z_over_h) /
                           (kTruncatedGaussianColumnIntegral * scale_height);
    return std::isfinite(density) && density >= 0.0 ? density : 0.0;
}

}  // namespace sirius::core::volumetric_disk
