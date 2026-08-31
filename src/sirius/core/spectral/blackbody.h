#pragma once

// Physically-based spectral rendering primitives: Planck blackbody radiance,
// Wien and Stefan-Boltzmann laws, CIE 1931 colour matching (Gaussian fits),
// wavelength/XYZ/sRGB conversion and redshift.
// Ported from PHSP001A.h.
// Reference: Planck (1901); CIE 1931 observer; sRGB IEC 61966-2-1:1999.

#include "sirius/core/constants.h"
#include "sirius/core/srgb_transfer.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace sirius::core::spectral {

// Planck radiation constants sourced from the constants authority.
constexpr double kPlanckC1 = constants::physical::kPlanckC1;
constexpr double kPlanckC2 = constants::physical::kPlanckC2;

// Linear RGB colour, 0-1 range.
struct Rgb {
    float r, g, b;

    Rgb() : r(0), g(0), b(0) {}
    Rgb(float r_, float g_, float b_) : r(r_), g(g_), b(b_) {}

    Rgb operator*(float s) const { return Rgb(r * s, g * s, b * s); }
    Rgb operator+(const Rgb& o) const { return Rgb(r + o.r, g + o.g, b + o.b); }
    Rgb& operator+=(const Rgb& o) {
        r += o.r;
        g += o.g;
        b += o.b;
        return *this;
    }
};

// CIE 1931 XYZ colour.
struct Xyz {
    float X, Y, Z;

    Xyz() : X(0), Y(0), Z(0) {}
    Xyz(float x, float y, float z) : X(x), Y(y), Z(z) {}

    Xyz operator*(float s) const { return Xyz(X * s, Y * s, Z * s); }
    Xyz operator+(const Xyz& o) const { return Xyz(X + o.X, Y + o.Y, Z + o.Z); }
    Xyz& operator+=(const Xyz& o) {
        X += o.X;
        Y += o.Y;
        Z += o.Z;
        return *this;
    }
};

// Planck blackbody radiance B(lambda, T) = (2hc^2/lambda^5)/(exp(hc/lambda k T) - 1),
// lambda in metres, T in Kelvin; returns W sr^-1 m^-3.
inline double PlanckRadiance(double lambda, double T) {
    if (T <= 0 || lambda <= 0) return 0.0;

    double x = kPlanckC2 / (lambda * T);
    if (x > 700) return 0.0;  // Prevent overflow.

    return kPlanckC1 / (std::pow(lambda, 5) * std::expm1(x));
}

// Wien's displacement law: peak wavelength (m) for temperature T (K).
inline double WienPeakWavelength(double T) {
    constexpr double kWienConstant = 2.897771955e-3;  // m K.
    return kWienConstant / T;
}

// Stefan-Boltzmann total radiance sigma T^4 (W/m^2).
inline double StefanBoltzmannRadiance(double T) {
    constexpr double kSigma = constants::physical::kStefanBoltzmann;
    return kSigma * std::pow(T, 4);
}

// CIE 1931 X colour matching function (Gaussian fit), lambda in nm.
inline double CieX(double lambda_nm) {
    double t1 = (lambda_nm - 442.0) * ((lambda_nm < 442.0) ? 0.0624 : 0.0374);
    double t2 = (lambda_nm - 599.8) * ((lambda_nm < 599.8) ? 0.0264 : 0.0323);
    double t3 = (lambda_nm - 501.1) * ((lambda_nm < 501.1) ? 0.0490 : 0.0382);
    return 0.362 * std::exp(-0.5 * t1 * t1) + 1.056 * std::exp(-0.5 * t2 * t2) -
           0.065 * std::exp(-0.5 * t3 * t3);
}

// CIE 1931 Y colour matching function (Gaussian fit), lambda in nm.
inline double CieY(double lambda_nm) {
    double t1 = (lambda_nm - 568.8) * ((lambda_nm < 568.8) ? 0.0213 : 0.0247);
    double t2 = (lambda_nm - 530.9) * ((lambda_nm < 530.9) ? 0.0613 : 0.0322);
    return 0.821 * std::exp(-0.5 * t1 * t1) + 0.286 * std::exp(-0.5 * t2 * t2);
}

// CIE 1931 Z colour matching function (Gaussian fit), lambda in nm.
inline double CieZ(double lambda_nm) {
    double t1 = (lambda_nm - 437.0) * ((lambda_nm < 437.0) ? 0.0845 : 0.0278);
    double t2 = (lambda_nm - 459.0) * ((lambda_nm < 459.0) ? 0.0385 : 0.0725);
    return 1.217 * std::exp(-0.5 * t1 * t1) + 0.681 * std::exp(-0.5 * t2 * t2);
}

// Wavelength (nm) to CIE XYZ; zero outside the 380-780 nm visible range.
inline Xyz WavelengthToXyz(double lambda_nm) {
    if (lambda_nm < 380.0 || lambda_nm > 780.0) {
        return Xyz(0, 0, 0);
    }
    return Xyz(static_cast<float>(CieX(lambda_nm)), static_cast<float>(CieY(lambda_nm)),
               static_cast<float>(CieZ(lambda_nm)));
}

// XYZ to linear RGB via the sRGB D65 primaries.
inline Rgb XyzToLinearRgb(const Xyz& xyz) {
    float r = 3.2404542f * xyz.X - 1.5371385f * xyz.Y - 0.4985314f * xyz.Z;
    float g = -0.9692660f * xyz.X + 1.8760108f * xyz.Y + 0.0415560f * xyz.Z;
    float b = 0.0556434f * xyz.X - 0.2040259f * xyz.Y + 1.0572252f * xyz.Z;
    return Rgb(r, g, b);
}

// sRGB gamma correction of a single linear channel.
inline float SrgbGamma(float linear) { return colour::EncodeSrgbChannel(linear); }

// Linear RGB to gamma-corrected sRGB.
inline Rgb LinearToSrgb(const Rgb& linear) {
    return Rgb(colour::EncodeClippedSrgbChannel(linear.r),
               colour::EncodeClippedSrgbChannel(linear.g),
               colour::EncodeClippedSrgbChannel(linear.b));
}

// Blackbody temperature T (K) to normalised, display-gamut linear RGB
// (brightest channel = 1). Negative out-of-gamut primary values are clipped;
// they are not negative radiance.
inline Rgb BlackbodyToRgb(double T) {
    if (!std::isfinite(T) || T <= 0) return Rgb(0, 0, 0);

    // Integrate over the visible spectrum.
    Xyz xyz;
    constexpr int kNSamples = 32;
    constexpr double kDLambda =
        (constants::spectral::kLambdaMax - constants::spectral::kLambdaMin) / kNSamples;

    for (int i = 0; i < kNSamples; i++) {
        double lambda = constants::spectral::kLambdaMin + (i + 0.5) * kDLambda;
        double lambda_nm = lambda * 1e9;

        double radiance = PlanckRadiance(lambda, T);
        Xyz sample = WavelengthToXyz(lambda_nm);

        xyz += sample * static_cast<float>(radiance * kDLambda);
    }

    // Convert to RGB.
    Rgb rgb = XyzToLinearRgb(xyz);
    rgb.r = std::max(rgb.r, 0.0f);
    rgb.g = std::max(rgb.g, 0.0f);
    rgb.b = std::max(rgb.b, 0.0f);

    // Normalise to max = 1.
    float maxVal = std::max({rgb.r, rgb.g, rgb.b, 0.001f});
    return Rgb(rgb.r / maxVal, rgb.g / maxVal, rgb.b / maxVal);
}

// Relativistic Doppler factor nu_obs/nu_emit for radial velocity (positive =
// receding); c defaults to geometric units.
inline double DopplerFactor(double velocity, double c = 1.0) {
    if (!std::isfinite(velocity) || !std::isfinite(c) || c <= 0.0) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    const double numerator = c - velocity;
    const double denominator = c + velocity;
    if (!std::isfinite(numerator) || !std::isfinite(denominator) || numerator <= 0.0 ||
        denominator <= 0.0) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    return std::sqrt(numerator / denominator);
}

// Combined gravitational and Doppler redshift z, given g_tt at emission and
// observation and the radial velocity at emission.
inline double TotalRedshift(double g_tt_emit, double g_tt_obs, double velocity) {
    // For static observers in a stationary metric, the measured frequency is
    // proportional to 1/sqrt(-g_tt).  Therefore
    //
    //   g_grav = nu_obs / nu_emit = sqrt(g_tt_emit / g_tt_obs),
    //
    // with both metric components negative.  The previous quotient was the
    // reciprocal and turned a pure gravitational redshift into a blueshift.
    // This utility represents only the static-observer domain; ergoregions and
    // moving observers require the invariant k.u transfer authority used by the
    // live renderer.
    if (!std::isfinite(g_tt_emit) || !std::isfinite(g_tt_obs) || g_tt_emit >= 0.0 ||
        g_tt_obs >= 0.0) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    double grav_factor = std::sqrt(g_tt_emit / g_tt_obs);

    double dopp_factor = DopplerFactor(velocity);

    double total_factor = grav_factor * dopp_factor;
    if (!std::isfinite(total_factor) || total_factor <= 0.0) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    // z = (lambda_obs - lambda_emit) / lambda_emit = 1/factor - 1.
    return (1.0 / total_factor) - 1.0;
}

}  // namespace sirius::core::spectral
