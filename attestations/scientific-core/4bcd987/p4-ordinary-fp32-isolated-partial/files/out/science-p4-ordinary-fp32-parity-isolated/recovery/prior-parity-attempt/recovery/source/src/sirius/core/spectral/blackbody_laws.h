#pragma once

// Single host authority for the represented blackbody laws. Quantities are
// named with their physical measure and units:
//
//   B_lambda(lambda, T)  spectral radiance per metre [W sr^-1 m^-3]
//   lambda_max(T)        Wien peak wavelength [m]
//   M(T)                 hemispheric radiant exitance [W m^-2]
//
// Malformed inputs and finite inputs whose result cannot be represented as a
// finite double decline. A physically valid radiance below the double range is
// represented by zero rather than by a substituted wavelength or temperature.

#include "sirius/core/constants.h"

#include <cmath>
#include <limits>
#include <optional>

namespace sirius::core::spectral {

[[nodiscard]] inline std::optional<double> TryPlanckSpectralRadiancePerMetre(
    double wavelength_metres, double temperature_kelvin) {
    if (!std::isfinite(wavelength_metres) || !(wavelength_metres > 0.0) ||
        !std::isfinite(temperature_kelvin) || !(temperature_kelvin > 0.0)) {
        return std::nullopt;
    }

    const double log_wavelength = std::log(wavelength_metres);
    const double log_temperature = std::log(temperature_kelvin);
    const double log_exponent =
        std::log(constants::physical::kPlanckC2) - log_wavelength - log_temperature;
    const double log_maximum = std::log(std::numeric_limits<double>::max());

    // If the dimensionless Planck exponent itself exceeds the finite range,
    // its exponential dominates every representable lambda^-5 numerator.
    if (log_exponent > log_maximum) return 0.0;

    const double wavelength_temperature = wavelength_metres * temperature_kelvin;
    const double exponent = std::isfinite(wavelength_temperature) && wavelength_temperature > 0.0
                                ? constants::physical::kPlanckC2 / wavelength_temperature
                                : std::exp(log_exponent);

    double exponential_minus_one = std::numeric_limits<double>::infinity();
    if (exponent <= log_maximum) exponential_minus_one = std::expm1(exponent);

    // Preserve the direct expression in its stable range. This also keeps the
    // common physical domain as close as possible to the defining equation.
    const double denominator = std::pow(wavelength_metres, 5.0) * exponential_minus_one;
    if (std::isfinite(denominator) && denominator > 0.0) {
        const double radiance = constants::physical::kPlanckC1 / denominator;
        if (std::isfinite(radiance) && radiance >= 0.0) return radiance;
    }

    // Resolve intermediate overflow/underflow without changing the law. For
    // large x, log(expm1(x)) = x + log(1-exp(-x)); for an exponent that rounded
    // to zero, log(x) remains available from the input logarithms.
    double log_exponential_minus_one = log_exponent;
    if (exponent > 50.0) {
        log_exponential_minus_one = exponent + std::log1p(-std::exp(-exponent));
    } else if (exponential_minus_one > 0.0) {
        log_exponential_minus_one = std::log(exponential_minus_one);
    }
    const double log_radiance =
        std::log(constants::physical::kPlanckC1) - 5.0 * log_wavelength - log_exponential_minus_one;
    if (!std::isfinite(log_radiance) || log_radiance > log_maximum) {
        return std::nullopt;
    }
    const double radiance = std::exp(log_radiance);
    if (!std::isfinite(radiance) || radiance < 0.0) return std::nullopt;
    return radiance;
}

[[nodiscard]] inline std::optional<double> TryWienPeakWavelength(double temperature_kelvin) {
    if (!std::isfinite(temperature_kelvin) || !(temperature_kelvin > 0.0)) {
        return std::nullopt;
    }
    const double wavelength = constants::physical::kWienB / temperature_kelvin;
    if (!std::isfinite(wavelength) || !(wavelength > 0.0)) return std::nullopt;
    return wavelength;
}

[[nodiscard]] inline std::optional<double> TryStefanBoltzmannExitance(double temperature_kelvin) {
    if (!std::isfinite(temperature_kelvin) || !(temperature_kelvin > 0.0)) {
        return std::nullopt;
    }
    const double temperature_squared = temperature_kelvin * temperature_kelvin;
    const double exitance =
        constants::physical::kStefanBoltzmann * temperature_squared * temperature_squared;
    if (!std::isfinite(exitance) || exitance < 0.0) return std::nullopt;
    return exitance;
}

}  // namespace sirius::core::spectral
