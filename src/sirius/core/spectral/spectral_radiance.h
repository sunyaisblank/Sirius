#pragma once

// 32-bin spectral radiance across the visible band (380-780 nm, 12.5 nm/bin) for
// colour-accurate rendering, with Planck emission, g^4 redshift scaling, and CIE
// colour conversion. Ported from MTSB001A.h.
// Reference: James et al. (2015), "DNGR", Section 3.4.

#include "sirius/core/constants.h"

#include <algorithm>
#include <cmath>

namespace sirius::core {

constexpr int kNumWavelengthBins = 32;
constexpr double kLambdaMin = 380.0;                                            // nm
constexpr double kLambdaMax = 780.0;                                            // nm
constexpr double kLambdaStep = (kLambdaMax - kLambdaMin) / kNumWavelengthBins;  // 12.5 nm

// CIE 1931 2-degree observer colour matching functions, sampled at bin centres.
namespace cie1931 {

// x-bar(lambda).
constexpr double kXBar[32] = {0.014310, 0.043510, 0.134380, 0.283900, 0.348280, 0.336200, 0.290800,
                              0.195360, 0.095640, 0.032010, 0.004900, 0.009300, 0.063270, 0.165500,
                              0.290400, 0.433450, 0.594500, 0.762100, 0.916300, 1.026300, 1.062200,
                              1.002600, 0.854450, 0.642400, 0.447900, 0.283500, 0.164900, 0.087400,
                              0.046770, 0.022700, 0.011359, 0.005790};

// y-bar(lambda) (luminance).
constexpr double kYBar[32] = {0.000396, 0.001210, 0.004000, 0.011600, 0.023000, 0.038000, 0.060000,
                              0.091000, 0.139020, 0.208020, 0.323000, 0.503000, 0.710000, 0.862000,
                              0.954000, 0.995000, 0.995000, 0.952000, 0.870000, 0.757000, 0.631000,
                              0.503000, 0.381000, 0.265000, 0.175000, 0.107000, 0.061000, 0.032000,
                              0.017000, 0.008210, 0.004102, 0.002091};

// z-bar(lambda).
constexpr double kZBar[32] = {0.067850, 0.207400, 0.645600, 1.385600, 1.747060, 1.772110, 1.669200,
                              1.287640, 0.813000, 0.465180, 0.272000, 0.158200, 0.078200, 0.042200,
                              0.020300, 0.008700, 0.003900, 0.002100, 0.001650, 0.001100, 0.000800,
                              0.000340, 0.000190, 0.000050, 0.000020, 0.000000, 0.000000, 0.000000,
                              0.000000, 0.000000, 0.000000, 0.000000};

}  // namespace cie1931

// 32-bin visible spectrum representation.
struct SpectralRadiance {
    double L[kNumWavelengthBins];  // Radiance per bin [W/(m^2 sr nm)].

    SpectralRadiance() {
        for (int i = 0; i < kNumWavelengthBins; ++i) L[i] = 0;
    }

    static SpectralRadiance Zero() { return SpectralRadiance(); }

    static double Wavelength(int bin) { return kLambdaMin + (bin + 0.5) * kLambdaStep; }

    static int BinIndex(double wavelength) {
        return static_cast<int>((wavelength - kLambdaMin) / kLambdaStep);
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
        constexpr double hc = constants::physical::kPlanck * constants::physical::kSpeedOfLight;
        constexpr double hc2 = constants::physical::kPlanck * constants::physical::kSpeedOfLight *
                               constants::physical::kSpeedOfLight;

        SpectralRadiance result;

        for (int i = 0; i < kNumWavelengthBins; ++i) {
            double lambda = Wavelength(i) * 1e-9;  // nm to m.
            double x = hc / (lambda * constants::physical::kBoltzmann * temperature);

            if (x > 700) {
                result.L[i] = 0;  // Avoid overflow.
            } else {
                double B = (2 * hc2 / std::pow(lambda, 5)) / (std::exp(x) - 1);
                result.L[i] = B * 1e-9;  // Convert to per-nm.
            }
        }

        return result;
    }

    // Gravitational/Doppler redshift with g^4 intensity scaling: g = nu_obs/nu_emit,
    // lambda_obs = lambda_emit/g, I_obs = g^4 I_emit.
    SpectralRadiance ApplyRedshift(double g) const {
        SpectralRadiance result;
        double g4 = g * g * g * g;

        for (int i = 0; i < kNumWavelengthBins; ++i) {
            double lambda_emit = Wavelength(i);
            double lambda_obs = lambda_emit / g;

            int j = BinIndex(lambda_obs);

            if (j >= 0 && j < kNumWavelengthBins) {
                result.L[j] += L[i] * g4;
            }
            // Else: shifted out of the visible range (UV or IR).
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

    // Spectrum to CIE XYZ.
    Xyz ToXyz() const {
        Xyz xyz = {0, 0, 0};

        for (int i = 0; i < kNumWavelengthBins; ++i) {
            xyz.X += L[i] * cie1931::kXBar[i] * kLambdaStep;
            xyz.Y += L[i] * cie1931::kYBar[i] * kLambdaStep;
            xyz.Z += L[i] * cie1931::kZBar[i] * kLambdaStep;
        }

        return xyz;
    }

    struct Aces {
        double r, g, b;
    };

    // XYZ to ACES AP0 (D60 white point).
    Aces ToAces() const {
        Xyz xyz = ToXyz();

        Aces aces;
        aces.r = 1.0498110175 * xyz.X + 0.0000000000 * xyz.Y - 0.0000974845 * xyz.Z;
        aces.g = -0.4959030231 * xyz.X + 1.3733130458 * xyz.Y + 0.0982400361 * xyz.Z;
        aces.b = 0.0000000000 * xyz.X + 0.0000000000 * xyz.Y + 0.9912520182 * xyz.Z;

        return aces;
    }

    struct Rgb {
        double r, g, b;
    };

    // XYZ to gamma-corrected sRGB (D65), for preview.
    Rgb ToSrgb() const {
        Xyz xyz = ToXyz();

        double r_lin = 3.2404542 * xyz.X - 1.5371385 * xyz.Y - 0.4985314 * xyz.Z;
        double g_lin = -0.9692660 * xyz.X + 1.8760108 * xyz.Y + 0.0415560 * xyz.Z;
        double b_lin = 0.0556434 * xyz.X - 0.2040259 * xyz.Y + 1.0572252 * xyz.Z;

        auto gamma = [](double u) {
            if (u <= 0.0031308) return 12.92 * u;
            return 1.055 * std::pow(u, 1.0 / 2.4) - 0.055;
        };

        Rgb rgb;
        rgb.r = std::clamp(gamma(r_lin), 0.0, 1.0);
        rgb.g = std::clamp(gamma(g_lin), 0.0, 1.0);
        rgb.b = std::clamp(gamma(b_lin), 0.0, 1.0);

        return rgb;
    }
};

}  // namespace sirius::core
