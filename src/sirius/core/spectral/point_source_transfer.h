#pragma once

#include "sirius/core/spectral/blackbody.h"

#include <array>
#include <cmath>
#include <expected>

namespace sirius::core::spectral {

enum class PointSourceTransferFailure { InvalidInput, Arithmetic };

// Relative reference-normalized Planck colour radiance, not calibrated SI or
// bolometric flux. g = nu_camera / nu_source. The catalogue's original flux F
// and unshifted N(T) are retained: F * response_density * C(g*T) / N(T).
// C uses the existing 32-midpoint visible colour response, including its gamut
// clipping. I_lambda/nu^5 transport is already contained in C(g*T): multiplying
// by another g^4, determinant, or pixel area would count that factor twice.
//
// No lookup-temperature clamp is imposed. Finite positive T, g and finite
// nonnegative F, density are valid inputs, subject to representable arithmetic
// in the existing float source normalizer and double shifted integral/output. Planck underflow
// retains the law helper's zero convention. Errors never carry partial RGB.
// This local transfer assumes the caller has established that g and the source
// map are sufficiently constant over the detector support, or refined it.
[[nodiscard]] inline std::expected<std::array<double, 3>, PointSourceTransferFailure>
TransferPointSourceBand(double temperature_kelvin, double camera_over_source_frequency,
                        double relative_flux, double response_density) {
    if (!std::isfinite(temperature_kelvin) || !(temperature_kelvin > 0.0) ||
        !std::isfinite(camera_over_source_frequency) || !(camera_over_source_frequency > 0.0) ||
        !std::isfinite(relative_flux) || relative_flux < 0.0 || !std::isfinite(response_density) ||
        response_density < 0.0) {
        return std::unexpected(PointSourceTransferFailure::InvalidInput);
    }
    const double shifted_temperature = temperature_kelvin * camera_over_source_frequency;
    if (!std::isfinite(shifted_temperature) || !(shifted_temperature > 0.0))
        return std::unexpected(PointSourceTransferFailure::Arithmetic);
    const auto source = BlackbodyBandRgb(temperature_kelvin);
    if (!source || !std::isfinite(source->r) || !std::isfinite(source->g) ||
        !std::isfinite(source->b))
        return std::unexpected(PointSourceTransferFailure::Arithmetic);
    const float normalizer = std::max({source->r, source->g, source->b, 0.001f});
    std::array<double, 3> channels{};
    {
        // Same 32-bin response, evaluated before narrowing: a faint shifted
        // band can become representable after multiplication by source flux.
        // The SOURCE normalizer above deliberately retains catalogue rounding.
        constexpr int samples = 32;
        constexpr double step =
            (constants::spectral::kLambdaMax - constants::spectral::kLambdaMin) / samples;
        std::array<double, 3> xyz{};
        for (int i = 0; i < samples; ++i) {
            const double wavelength = constants::spectral::kLambdaMin + (i + .5) * step;
            const auto radiance =
                TryPlanckSpectralRadiancePerMetre(wavelength, shifted_temperature);
            if (!radiance) return std::unexpected(PointSourceTransferFailure::Arithmetic);
            const Xyz matching = WavelengthToXyz(wavelength * 1e9);
            const double weight = *radiance * step;
            xyz[0] += matching.X * weight;
            xyz[1] += matching.Y * weight;
            xyz[2] += matching.Z * weight;
        }
        for (double value : xyz)
            if (!std::isfinite(value))
                return std::unexpected(PointSourceTransferFailure::Arithmetic);
        const auto rgb = colour::XyzD65ToLinearSrgb(xyz[0], xyz[1], xyz[2]);
        channels = {rgb.r, rgb.g, rgb.b};
        for (auto& channel : channels) {
            if (!std::isfinite(channel))
                return std::unexpected(PointSourceTransferFailure::Arithmetic);
            channel = std::max(channel, 0.0);
        }
    }
    std::array<double, 3> result{};
    for (std::size_t channel = 0; channel < result.size(); ++channel) {
        const double relative_colour = channels[channel] / normalizer;
        result[channel] = relative_colour * relative_flux * response_density;
        if (channels[channel] > 0.0 && relative_flux > 0.0 && response_density > 0.0 &&
            (!std::isfinite(result[channel]) || result[channel] == 0.0)) {
            // Recover a representable product when its chosen intermediate
            // multiplication overflowed or underflowed. The ordinary path
            // above avoids unnecessary exponent manipulation.
            int colour_exponent, normalization_exponent, flux_exponent, density_exponent;
            const double mantissa =
                (std::frexp(channels[channel], &colour_exponent) /
                 std::frexp(static_cast<double>(normalizer), &normalization_exponent)) *
                std::frexp(relative_flux, &flux_exponent) *
                std::frexp(response_density, &density_exponent);
            result[channel] = std::scalbn(mantissa, colour_exponent - normalization_exponent +
                                                        flux_exponent + density_exponent);
        }
        if (!std::isfinite(result[channel]))
            return std::unexpected(PointSourceTransferFailure::Arithmetic);
    }
    return result;
}

}  // namespace sirius::core::spectral
