// =============================================================================
// PHSP001A.h - Spectral Rendering Utilities
// Component ID: PHSP001A (Physics/Spectral)
// =============================================================================
//
// MATHEMATICAL FOUNDATION
// =======================
// This module implements physically-based spectral rendering for relativistic
// ray tracing. Key features:
//
// 1. Blackbody radiation (Planck's law):
//    B(λ,T) = (2hc²/λ⁵) × 1/(exp(hc/λkT) - 1)
//
// 2. Relativistic frequency shift:
//    ν_obs/ν_emit = (k·u)_emit / (k·u)_obs
//    where k = photon 4-momentum, u = observer 4-velocity
//
// 3. Color conversion pipeline:
//    λ → XYZ (CIE 1931) → sRGB (with gamma)
//
// REFERENCES:
//   - Planck, M. (1901). On the Law of Distribution of Energy
//   - CIE 1931 Color Matching Functions
//   - sRGB IEC 61966-2-1:1999
// =============================================================================
#pragma once

#include <cmath>
#include <algorithm>

namespace Sirius {
namespace Spectral {

// =============================================================================
// Physical Constants
// =============================================================================
constexpr double h_PLANCK = 6.62607015e-34;    // Planck constant (J·s)
constexpr double c_LIGHT = 2.99792458e8;       // Speed of light (m/s)
constexpr double k_BOLTZMANN = 1.380649e-23;   // Boltzmann constant (J/K)

// Derived constants for Planck function
constexpr double PLANCK_C1 = 2.0 * h_PLANCK * c_LIGHT * c_LIGHT;  // 2hc²
constexpr double PLANCK_C2 = h_PLANCK * c_LIGHT / k_BOLTZMANN;    // hc/k

// =============================================================================
// Color Structures
// =============================================================================

/// RGB color (linear, 0-1 range)
struct RGB {
    float r, g, b;
    
    RGB() : r(0), g(0), b(0) {}
    RGB(float r_, float g_, float b_) : r(r_), g(g_), b(b_) {}
    
    RGB operator*(float s) const { return RGB(r*s, g*s, b*s); }
    RGB operator+(const RGB& o) const { return RGB(r+o.r, g+o.g, b+o.b); }
    RGB& operator+=(const RGB& o) { r+=o.r; g+=o.g; b+=o.b; return *this; }
};

/// XYZ color (CIE 1931)
struct XYZ {
    float X, Y, Z;
    
    XYZ() : X(0), Y(0), Z(0) {}
    XYZ(float x, float y, float z) : X(x), Y(y), Z(z) {}
    
    XYZ operator*(float s) const { return XYZ(X*s, Y*s, Z*s); }
    XYZ operator+(const XYZ& o) const { return XYZ(X+o.X, Y+o.Y, Z+o.Z); }
    XYZ& operator+=(const XYZ& o) { X+=o.X; Y+=o.Y; Z+=o.Z; return *this; }
};

// =============================================================================
// Blackbody Radiation
// =============================================================================

/// @brief Planck blackbody radiance at given wavelength and temperature
/// @param lambda Wavelength in meters (e.g., 550e-9 for green)
/// @param T Temperature in Kelvin
/// @return Spectral radiance B(λ,T) in W·sr⁻¹·m⁻³
inline double planckRadiance(double lambda, double T) {
    if (T <= 0 || lambda <= 0) return 0.0;
    
    double x = PLANCK_C2 / (lambda * T);
    if (x > 700) return 0.0;  // Prevent overflow
    
    return PLANCK_C1 / (std::pow(lambda, 5) * (std::exp(x) - 1.0));
}

/// @brief Wien's displacement law - peak wavelength for temperature
/// @param T Temperature in Kelvin
/// @return Peak wavelength in meters
inline double wienPeakWavelength(double T) {
    constexpr double WIEN_CONSTANT = 2.897771955e-3;  // m·K
    return WIEN_CONSTANT / T;
}

/// @brief Stefan-Boltzmann total radiance
/// @param T Temperature in Kelvin
/// @return Total radiance in W/m²
inline double stefanBoltzmannRadiance(double T) {
    constexpr double SIGMA = 5.670374419e-8;  // Stefan-Boltzmann constant
    return SIGMA * std::pow(T, 4);
}

// =============================================================================
// CIE Color Matching Functions (Approximate)
// =============================================================================

/// @brief CIE 1931 X color matching function (Gaussian fit)
inline double cie_x(double lambda_nm) {
    double t1 = (lambda_nm - 442.0) * ((lambda_nm < 442.0) ? 0.0624 : 0.0374);
    double t2 = (lambda_nm - 599.8) * ((lambda_nm < 599.8) ? 0.0264 : 0.0323);
    double t3 = (lambda_nm - 501.1) * ((lambda_nm < 501.1) ? 0.0490 : 0.0382);
    return 0.362 * std::exp(-0.5 * t1 * t1) 
         + 1.056 * std::exp(-0.5 * t2 * t2) 
         - 0.065 * std::exp(-0.5 * t3 * t3);
}

/// @brief CIE 1931 Y color matching function (Gaussian fit)
inline double cie_y(double lambda_nm) {
    double t1 = (lambda_nm - 568.8) * ((lambda_nm < 568.8) ? 0.0213 : 0.0247);
    double t2 = (lambda_nm - 530.9) * ((lambda_nm < 530.9) ? 0.0613 : 0.0322);
    return 0.821 * std::exp(-0.5 * t1 * t1) 
         + 0.286 * std::exp(-0.5 * t2 * t2);
}

/// @brief CIE 1931 Z color matching function (Gaussian fit)
inline double cie_z(double lambda_nm) {
    double t1 = (lambda_nm - 437.0) * ((lambda_nm < 437.0) ? 0.0845 : 0.0278);
    double t2 = (lambda_nm - 459.0) * ((lambda_nm < 459.0) ? 0.0385 : 0.0725);
    return 1.217 * std::exp(-0.5 * t1 * t1) 
         + 0.681 * std::exp(-0.5 * t2 * t2);
}

// =============================================================================
// Color Conversion Functions
// =============================================================================

/// @brief Convert wavelength to XYZ color
/// @param lambda_nm Wavelength in nanometers (380-780 visible range)
/// @return XYZ color coordinates
inline XYZ wavelengthToXYZ(double lambda_nm) {
    if (lambda_nm < 380.0 || lambda_nm > 780.0) {
        return XYZ(0, 0, 0);
    }
    return XYZ(
        static_cast<float>(cie_x(lambda_nm)),
        static_cast<float>(cie_y(lambda_nm)),
        static_cast<float>(cie_z(lambda_nm))
    );
}

/// @brief Convert XYZ to linear RGB (sRGB primaries)
inline RGB xyzToLinearRGB(const XYZ& xyz) {
    // sRGB D65 matrix
    float r =  3.2404542f * xyz.X - 1.5371385f * xyz.Y - 0.4985314f * xyz.Z;
    float g = -0.9692660f * xyz.X + 1.8760108f * xyz.Y + 0.0415560f * xyz.Z;
    float b =  0.0556434f * xyz.X - 0.2040259f * xyz.Y + 1.0572252f * xyz.Z;
    return RGB(r, g, b);
}

/// @brief Apply sRGB gamma correction
inline float srgbGamma(float linear) {
    if (linear <= 0.0031308f) {
        return 12.92f * linear;
    }
    return 1.055f * std::pow(linear, 1.0f/2.4f) - 0.055f;
}

/// @brief Convert linear RGB to sRGB (gamma corrected)
inline RGB linearToSRGB(const RGB& linear) {
    return RGB(
        srgbGamma(std::clamp(linear.r, 0.0f, 1.0f)),
        srgbGamma(std::clamp(linear.g, 0.0f, 1.0f)),
        srgbGamma(std::clamp(linear.b, 0.0f, 1.0f))
    );
}

// =============================================================================
// Blackbody to RGB Conversion
// =============================================================================

/// @brief Convert blackbody temperature to RGB color
/// @param T Temperature in Kelvin
/// @return Normalized RGB color (brightest channel = 1)
inline RGB blackbodyToRGB(double T) {
    if (T <= 0) return RGB(0, 0, 0);
    
    // Integrate over visible spectrum
    XYZ xyz;
    constexpr int N_SAMPLES = 32;
    constexpr double LAMBDA_MIN = 380e-9;  // 380 nm
    constexpr double LAMBDA_MAX = 780e-9;  // 780 nm
    constexpr double D_LAMBDA = (LAMBDA_MAX - LAMBDA_MIN) / N_SAMPLES;
    
    for (int i = 0; i < N_SAMPLES; i++) {
        double lambda = LAMBDA_MIN + (i + 0.5) * D_LAMBDA;
        double lambda_nm = lambda * 1e9;
        
        double radiance = planckRadiance(lambda, T);
        XYZ sample = wavelengthToXYZ(lambda_nm);
        
        xyz += sample * static_cast<float>(radiance * D_LAMBDA);
    }
    
    // Convert to RGB
    RGB rgb = xyzToLinearRGB(xyz);
    
    // Normalize to max = 1
    float maxVal = std::max({rgb.r, rgb.g, rgb.b, 0.001f});
    return RGB(rgb.r / maxVal, rgb.g / maxVal, rgb.b / maxVal);
}

// =============================================================================
// Redshift Functions
// =============================================================================

/// @brief Apply redshift to a color
/// @param color Original color
/// @param z Redshift (z > 0 = redshift, z < 0 = blueshift)
/// @return Shifted color (approximate)
inline RGB applyRedshift(const RGB& color, float z) {
    if (std::abs(z) < 0.001f) return color;
    
    // Simple approximation: shift hue based on z
    // This is a simplified model - full spectral shift is more complex
    float factor = 1.0f / (1.0f + z);
    
    if (z > 0) {
        // Redshift: move towards red
        return RGB(
            color.r,
            color.g * factor,
            color.b * factor * factor
        );
    } else {
        // Blueshift: move towards blue
        float inv = 1.0f + z;  // < 1 for blueshift
        return RGB(
            color.r * inv * inv,
            color.g * inv,
            color.b
        );
    }
}

/// @brief Calculate relativistic Doppler factor
/// @param velocity Radial velocity (positive = away)
/// @param c Speed of light (use 1.0 for geometric units)
/// @return Doppler factor (ν_obs/ν_emit)
inline double dopplerFactor(double velocity, double c = 1.0) {
    double beta = velocity / c;
    if (std::abs(beta) >= 1.0) return 1.0;  // Invalid velocity
    
    return std::sqrt((1.0 - beta) / (1.0 + beta));
}

/// @brief Combined gravitational + Doppler redshift
/// @param g_tt Metric component g_tt at emission point
/// @param velocity Radial velocity at emission point
/// @return Total redshift z
inline double totalRedshift(double g_tt_emit, double g_tt_obs, double velocity) {
    // Gravitational redshift: sqrt(-g_tt_obs / -g_tt_emit)
    double grav_factor = std::sqrt(std::abs(g_tt_obs / g_tt_emit));
    
    // Doppler factor
    double dopp_factor = dopplerFactor(velocity);
    
    // Combined
    double total_factor = grav_factor * dopp_factor;
    
    // z = (λ_obs - λ_emit) / λ_emit = 1/factor - 1
    return (1.0 / total_factor) - 1.0;
}

} // namespace Spectral
} // namespace Sirius
