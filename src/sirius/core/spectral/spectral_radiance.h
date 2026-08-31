#pragma once

// 32-bin spectral radiance across the visible band (380-780 nm, 12.5 nm/bin) for
// colour-accurate rendering, with Planck emission, invariant redshift rebinning, and CIE
// tristimulus integration. Ported from MTSB001A.h.
// Reference: James et al. (2015), "DNGR", Section 3.4.

#include "sirius/core/cie1931_observer.h"
#include "sirius/core/constants.h"
#include "sirius/core/spectral/blackbody_laws.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace sirius::core {

constexpr int kNumWavelengthBins = 32;
constexpr double kLambdaMin = 380.0;                                            // nm
constexpr double kLambdaMax = 780.0;                                            // nm
constexpr double kLambdaStep = (kLambdaMax - kLambdaMin) / kNumWavelengthBins;  // 12.5 nm

// 32-bin visible spectrum representation.
struct SpectralRadiance {
    double L[kNumWavelengthBins];  // Radiance per bin [W/(m^2 sr nm)].

    SpectralRadiance() {
        for (int i = 0; i < kNumWavelengthBins; ++i) L[i] = 0;
    }

    static SpectralRadiance Zero() { return SpectralRadiance(); }

    static double Wavelength(int bin) { return kLambdaMin + (bin + 0.5) * kLambdaStep; }

    static int BinIndex(double wavelength) {
        if (!std::isfinite(wavelength)) return -1;
        return static_cast<int>(std::floor((wavelength - kLambdaMin) / kLambdaStep));
    }

    SpectralRadiance operator+(const SpectralRadiance& o) const {
        SpectralRadiance result;
        for (int i = 0; i < kNumWavelengthBins; ++i) {
            result.L[i] = L[i] + o.L[i];
        }
        return result;
    }

    SpectralRadiance operator*(double s) const {
        SpectralRadiance result;
        for (int i = 0; i < kNumWavelengthBins; ++i) {
            result.L[i] = L[i] * s;
        }
        return result;
    }

    SpectralRadiance& operator+=(const SpectralRadiance& o) {
        for (int i = 0; i < kNumWavelengthBins; ++i) {
            L[i] += o.L[i];
        }
        return *this;
    }

    SpectralRadiance& operator*=(double s) {
        for (int i = 0; i < kNumWavelengthBins; ++i) {
            L[i] *= s;
        }
        return *this;
    }

    // Planck blackbody spectrum B(lambda, T) = (2hc^2/lambda^5)/(exp(hc/lambda k T) - 1).
    static SpectralRadiance Blackbody(double temperature) {
        if (!std::isfinite(temperature) || temperature <= 0.0) return Zero();

        SpectralRadiance result;

        for (int i = 0; i < kNumWavelengthBins; ++i) {
            const double wavelength_metres = Wavelength(i) * 1.0e-9;
            const std::optional<double> radiance =
                spectral::TryPlanckSpectralRadiancePerMetre(wavelength_metres, temperature);
            if (!radiance.has_value()) return Zero();
            result.L[i] = *radiance * 1.0e-9;  // Convert density from per-metre to per-nm.
        }

        return result;
    }

    // Gravitational/Doppler transfer for wavelength-specific radiance.
    // Liouville invariance I_nu/nu^3 gives
    //
    //   lambda_obs = lambda_emit/g,
    //   I_lambda,obs(lambda_obs) = g^5 I_lambda,emit(g lambda_obs).
    //
    // Each stored value is a density per nanometre over its bin. Rebinning by
    // interval overlap applies the g^5 density and the 1/g wavelength stretch,
    // so the integrated bolometric radiance is exactly g^4 when the shifted
    // support remains in the represented band.
    SpectralRadiance ApplyRedshift(double g) const {
        SpectralRadiance result;
        if (!std::isfinite(g) || g <= 0.0) return result;
        const double g2 = g * g;
        const double g5 = g2 * g2 * g;

        for (int i = 0; i < kNumWavelengthBins; ++i) {
            if (!std::isfinite(L[i]) || L[i] == 0.0) continue;
            const double observed_lower = (kLambdaMin + i * kLambdaStep) / g;
            const double observed_upper = (kLambdaMin + (i + 1) * kLambdaStep) / g;
            const int first_bin = std::max(0, BinIndex(observed_lower));
            const int last_bin = std::min(
                kNumWavelengthBins - 1,
                BinIndex(std::nextafter(observed_upper, -std::numeric_limits<double>::infinity())));
            for (int j = first_bin; j <= last_bin; ++j) {
                const double destination_lower = kLambdaMin + j * kLambdaStep;
                const double destination_upper = destination_lower + kLambdaStep;
                const double overlap =
                    std::max(0.0, std::min(observed_upper, destination_upper) -
                                      std::max(observed_lower, destination_lower));
                result.L[j] += L[i] * g5 * overlap / kLambdaStep;
            }
        }

        return result;
    }

    // Total energy, the integral over the spectrum.
    double TotalEnergy() const {
        double total = 0;
        for (int i = 0; i < kNumWavelengthBins; ++i) {
            total += L[i] * kLambdaStep;
        }
        return total;
    }

    struct Xyz {
        double X, Y, Z;
    };

    // Spectrum to unnormalised CIE XYZ tristimulus radiance. L is a density per
    // nanometre and kLambdaStep is in nanometres, so the result retains the
    // radiometric scale W/(m^2 sr). No illuminant, white, or relative-exposure
    // normalisation is introduced at this boundary.
    Xyz ToXyz() const {
        Xyz xyz = {0, 0, 0};

        for (int i = 0; i < kNumWavelengthBins; ++i) {
            const colour::CieXyzMatching matching = colour::Cie1931TwoDegreeFit(Wavelength(i));
            xyz.X += L[i] * matching.x_bar * kLambdaStep;
            xyz.Y += L[i] * matching.y_bar * kLambdaStep;
            xyz.Z += L[i] * matching.z_bar * kLambdaStep;
        }

        return xyz;
    }
};

}  // namespace sirius::core
