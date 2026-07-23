#pragma once

// Physically-based spectral rendering primitives: Planck blackbody radiance,
// Wien and Stefan-Boltzmann laws, CIE 1931 colour matching (Gaussian fits),
// wavelength/XYZ/sRGB conversion, redshift, and limb darkening.
// Ported from PHSP001A.h.
// Reference: Planck (1901); CIE 1931 observer; sRGB IEC 61966-2-1:1999.

#include "sirius/core/constants.h"

#include <algorithm>
#include <cmath>

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

    return kPlanckC1 / (std::pow(lambda, 5) * (std::exp(x) - 1.0));
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
inline float SrgbGamma(float linear) {
    if (linear <= 0.0031308f) {
        return 12.92f * linear;
    }
    return 1.055f * std::pow(linear, 1.0f / 2.4f) - 0.055f;
}

// Linear RGB to gamma-corrected sRGB.
inline Rgb LinearToSrgb(const Rgb& linear) {
    return Rgb(SrgbGamma(std::clamp(linear.r, 0.0f, 1.0f)),
               SrgbGamma(std::clamp(linear.g, 0.0f, 1.0f)),
               SrgbGamma(std::clamp(linear.b, 0.0f, 1.0f)));
}

// Blackbody temperature T (K) to normalised RGB (brightest channel = 1).
inline Rgb BlackbodyToRgb(double T) {
    if (T <= 0) return Rgb(0, 0, 0);

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

    // Normalise to max = 1.
    float maxVal = std::max({rgb.r, rgb.g, rgb.b, 0.001f});
    return Rgb(rgb.r / maxVal, rgb.g / maxVal, rgb.b / maxVal);
}

// Approximate hue shift for a redshift z (z > 0 redshift, z < 0 blueshift).
// A simplified colour model; the full spectral shift is more involved.
inline Rgb ApplyRedshift(const Rgb& color, float z) {
    if (std::abs(z) < 0.001f) return color;

    float factor = 1.0f / (1.0f + z);

    if (z > 0) {
        // Redshift: move towards red.
        return Rgb(color.r, color.g * factor, color.b * factor * factor);
    } else {
        // Blueshift: move towards blue.
        float inv = 1.0f + z;  // < 1 for blueshift.
        return Rgb(color.r * inv * inv, color.g * inv, color.b);
    }
}

// Relativistic Doppler factor nu_obs/nu_emit for radial velocity (positive =
// receding); c defaults to geometric units.
inline double DopplerFactor(double velocity, double c = 1.0) {
    double beta = velocity / c;
    if (std::abs(beta) >= 1.0) return 1.0;  // Invalid velocity.

    return std::sqrt((1.0 - beta) / (1.0 + beta));
}

// Combined gravitational and Doppler redshift z, given g_tt at emission and
// observation and the radial velocity at emission.
inline double TotalRedshift(double g_tt_emit, double g_tt_obs, double velocity) {
    // Gravitational redshift: sqrt(-g_tt_obs / -g_tt_emit).
    double grav_factor = std::sqrt(std::abs(g_tt_obs / g_tt_emit));

    double dopp_factor = DopplerFactor(velocity);

    double total_factor = grav_factor * dopp_factor;

    // z = (lambda_obs - lambda_emit) / lambda_emit = 1/factor - 1.
    return (1.0 / total_factor) - 1.0;
}

// Wavelength-dependent limb darkening coefficient u(lambda), lambda in nm.
// Law: I(theta)/I(0) = (1 + u cos(theta)) / (1 + u), with u decreasing towards
// red. Reference: Wade & Rucinski (1985), Claret (2000).
inline double LimbDarkeningCoeff(double lambda_nm) {
    // Linear interpolation between empirical values,
    // u(lambda) ~ 1.2 - 0.00114 * lambda(nm) over 400-800 nm.
    constexpr double kUBlue = 0.9;  // u at 400 nm.
    constexpr double kURed = 0.4;   // u at 700 nm.
    constexpr double kLambdaBlue = 400.0;
    constexpr double kLambdaRed = 700.0;

    if (lambda_nm <= kLambdaBlue) return kUBlue;
    if (lambda_nm >= kLambdaRed) return kURed;

    double t = (lambda_nm - kLambdaBlue) / (kLambdaRed - kLambdaBlue);
    return kUBlue + t * (kURed - kUBlue);
}

// Apply wavelength-dependent limb darkening to a colour; cos_theta = 0 at the
// edge, 1 face-on. Effective channel wavelengths R=650, G=550, B=450 nm.
inline Rgb ApplyLimbDarkening(const Rgb& color, double cos_theta) {
    if (cos_theta <= 0) return Rgb(0, 0, 0);

    double u_r = LimbDarkeningCoeff(650.0);
    double u_g = LimbDarkeningCoeff(550.0);
    double u_b = LimbDarkeningCoeff(450.0);

    // Limb darkening law I(theta)/I(0) = (1 + u cos(theta)) / (1 + u).
    float limb_r = static_cast<float>((1.0 + u_r * cos_theta) / (1.0 + u_r));
    float limb_g = static_cast<float>((1.0 + u_g * cos_theta) / (1.0 + u_g));
    float limb_b = static_cast<float>((1.0 + u_b * cos_theta) / (1.0 + u_b));

    return Rgb(color.r * limb_r, color.g * limb_g, color.b * limb_b);
}

}  // namespace sirius::core::spectral
