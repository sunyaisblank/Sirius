// MTSB001A.h - Spectral Radiance Representation
// Component ID: MTSB001A
// Purpose: 32-bin spectral radiance for accurate colour rendering
//
// MATHEMATICAL BASIS:
// Spectral radiance L(λ) [W/(m²·sr·nm)] is discretized into 32 bins
// covering the visible spectrum 380-780nm (12.5nm per bin).
//
// Under gravitational redshift with g-factor = ν_obs/ν_emit:
// - Observed wavelength: λ_obs = λ_emit / g
// - Intensity scaling: I_obs = g⁴ × I_emit (relativistic invariant)
//
// REFERENCE: James et al. (2015) "DNGR" Section 3.4

#pragma once

#include <cmath>
#include <algorithm>

namespace sirius::spectral {

//==============================================================================
// Spectral Constants
//==============================================================================

constexpr int NUM_WAVELENGTH_BINS = 32;
constexpr double LAMBDA_MIN = 380.0;  // nm
constexpr double LAMBDA_MAX = 780.0;  // nm
constexpr double LAMBDA_STEP = (LAMBDA_MAX - LAMBDA_MIN) / NUM_WAVELENGTH_BINS;  // 12.5nm

//==============================================================================
// CIE 1931 2° Observer Colour Matching Functions (sampled at bin centres)
//==============================================================================

namespace cie1931 {

// x̄(λ) colour matching function
constexpr double X_BAR[32] = {
    0.014310, 0.043510, 0.134380, 0.283900, 0.348280, 0.336200, 0.290800, 0.195360,
    0.095640, 0.032010, 0.004900, 0.009300, 0.063270, 0.165500, 0.290400, 0.433450,
    0.594500, 0.762100, 0.916300, 1.026300, 1.062200, 1.002600, 0.854450, 0.642400,
    0.447900, 0.283500, 0.164900, 0.087400, 0.046770, 0.022700, 0.011359, 0.005790
};

// ȳ(λ) colour matching function (luminance)
constexpr double Y_BAR[32] = {
    0.000396, 0.001210, 0.004000, 0.011600, 0.023000, 0.038000, 0.060000, 0.091000,
    0.139020, 0.208020, 0.323000, 0.503000, 0.710000, 0.862000, 0.954000, 0.995000,
    0.995000, 0.952000, 0.870000, 0.757000, 0.631000, 0.503000, 0.381000, 0.265000,
    0.175000, 0.107000, 0.061000, 0.032000, 0.017000, 0.008210, 0.004102, 0.002091
};

// z̄(λ) colour matching function
constexpr double Z_BAR[32] = {
    0.067850, 0.207400, 0.645600, 1.385600, 1.747060, 1.772110, 1.669200, 1.287640,
    0.813000, 0.465180, 0.272000, 0.158200, 0.078200, 0.042200, 0.020300, 0.008700,
    0.003900, 0.002100, 0.001650, 0.001100, 0.000800, 0.000340, 0.000190, 0.000050,
    0.000020, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000
};

} // namespace cie1931

//==============================================================================
// SpectralRadiance: 32-bin visible spectrum representation
//==============================================================================

struct SpectralRadiance {
    double L[NUM_WAVELENGTH_BINS];  // Radiance per bin [W/(m²·sr·nm)]
    
    //--------------------------------------------------------------------------
    // Constructors
    //--------------------------------------------------------------------------
    
    SpectralRadiance() {
        for (int i = 0; i < NUM_WAVELENGTH_BINS; ++i) L[i] = 0;
    }
    
    static SpectralRadiance zero() {
        return SpectralRadiance();
    }
    
    //--------------------------------------------------------------------------
    // Wavelength utilities
    //--------------------------------------------------------------------------
    
    static double wavelength(int bin) {
        return LAMBDA_MIN + (bin + 0.5) * LAMBDA_STEP;
    }
    
    static int binIndex(double wavelength) {
        return static_cast<int>((wavelength - LAMBDA_MIN) / LAMBDA_STEP);
    }
    
    //--------------------------------------------------------------------------
    // Arithmetic operators
    //--------------------------------------------------------------------------
    
    SpectralRadiance operator+(const SpectralRadiance& o) const {
        SpectralRadiance result;
        for (int i = 0; i < NUM_WAVELENGTH_BINS; ++i) {
            result.L[i] = L[i] + o.L[i];
        }
        return result;
    }
    
    SpectralRadiance operator*(double s) const {
        SpectralRadiance result;
        for (int i = 0; i < NUM_WAVELENGTH_BINS; ++i) {
            result.L[i] = L[i] * s;
        }
        return result;
    }
    
    SpectralRadiance& operator+=(const SpectralRadiance& o) {
        for (int i = 0; i < NUM_WAVELENGTH_BINS; ++i) {
            L[i] += o.L[i];
        }
        return *this;
    }
    
    SpectralRadiance& operator*=(double s) {
        for (int i = 0; i < NUM_WAVELENGTH_BINS; ++i) {
            L[i] *= s;
        }
        return *this;
    }
    
    //--------------------------------------------------------------------------
    // Blackbody spectrum (Planck function)
    //--------------------------------------------------------------------------
    
    static SpectralRadiance blackbody(double temperature) {
        // B(λ,T) = (2hc²/λ⁵) × 1/(exp(hc/λkT) - 1)
        // Constants in SI: h = 6.626e-34, c = 3e8, k = 1.381e-23
        // Working in nm, we use scaled constants
        
        constexpr double h = 6.62607015e-34;  // J·s
        constexpr double c = 2.99792458e8;    // m/s
        constexpr double k = 1.380649e-23;    // J/K
        constexpr double hc = h * c;
        constexpr double hc2 = h * c * c;
        
        SpectralRadiance result;
        
        for (int i = 0; i < NUM_WAVELENGTH_BINS; ++i) {
            double lambda = wavelength(i) * 1e-9;  // Convert nm to m
            double x = hc / (lambda * k * temperature);
            
            if (x > 700) {
                result.L[i] = 0;  // Avoid overflow
            } else {
                double B = (2 * hc2 / std::pow(lambda, 5)) / (std::exp(x) - 1);
                result.L[i] = B * 1e-9;  // Convert to per-nm
            }
        }
        
        return result;
    }
    
    //--------------------------------------------------------------------------
    // Gravitational/Doppler redshift with g⁴ intensity scaling
    //--------------------------------------------------------------------------
    
    SpectralRadiance applyRedshift(double g) const {
        // g = ν_obs / ν_emit = λ_emit / λ_obs
        // Observed wavelength: λ_obs = λ_emit / g
        // Intensity scaling: I_obs = g⁴ × I_emit (relativistic beaming)
        
        SpectralRadiance result;
        double g4 = g * g * g * g;
        
        for (int i = 0; i < NUM_WAVELENGTH_BINS; ++i) {
            double lambda_emit = wavelength(i);
            double lambda_obs = lambda_emit / g;
            
            int j = binIndex(lambda_obs);
            
            if (j >= 0 && j < NUM_WAVELENGTH_BINS) {
                result.L[j] += L[i] * g4;
            }
            // Else: shifted out of visible range (UV or IR)
        }
        
        return result;
    }
    
    //--------------------------------------------------------------------------
    // Total energy (integral over spectrum)
    //--------------------------------------------------------------------------
    
    double totalEnergy() const {
        double total = 0;
        for (int i = 0; i < NUM_WAVELENGTH_BINS; ++i) {
            total += L[i] * LAMBDA_STEP;
        }
        return total;
    }
    
    //--------------------------------------------------------------------------
    // Colour conversion: Spectrum → CIE XYZ
    //--------------------------------------------------------------------------
    
    struct XYZ {
        double X, Y, Z;
    };
    
    XYZ toXYZ() const {
        XYZ xyz = {0, 0, 0};
        
        for (int i = 0; i < NUM_WAVELENGTH_BINS; ++i) {
            xyz.X += L[i] * cie1931::X_BAR[i] * LAMBDA_STEP;
            xyz.Y += L[i] * cie1931::Y_BAR[i] * LAMBDA_STEP;
            xyz.Z += L[i] * cie1931::Z_BAR[i] * LAMBDA_STEP;
        }
        
        return xyz;
    }
    
    //--------------------------------------------------------------------------
    // Colour conversion: XYZ → ACES AP0
    //--------------------------------------------------------------------------
    
    struct ACES {
        double r, g, b;
    };
    
    ACES toACES() const {
        XYZ xyz = toXYZ();
        
        // ACES AP0 primaries, D60 white point
        // Matrix from XYZ to ACES (Academy Colour Encoding System)
        ACES aces;
        aces.r =  1.0498110175 * xyz.X + 0.0000000000 * xyz.Y - 0.0000974845 * xyz.Z;
        aces.g = -0.4959030231 * xyz.X + 1.3733130458 * xyz.Y + 0.0982400361 * xyz.Z;
        aces.b =  0.0000000000 * xyz.X + 0.0000000000 * xyz.Y + 0.9912520182 * xyz.Z;
        
        return aces;
    }
    
    //--------------------------------------------------------------------------
    // sRGB output (for preview)
    //--------------------------------------------------------------------------
    
    struct RGB {
        double r, g, b;
    };
    
    RGB toSRGB() const {
        XYZ xyz = toXYZ();
        
        // XYZ to linear sRGB (D65)
        double r_lin =  3.2404542 * xyz.X - 1.5371385 * xyz.Y - 0.4985314 * xyz.Z;
        double g_lin = -0.9692660 * xyz.X + 1.8760108 * xyz.Y + 0.0415560 * xyz.Z;
        double b_lin =  0.0556434 * xyz.X - 0.2040259 * xyz.Y + 1.0572252 * xyz.Z;
        
        // Gamma correction (sRGB transfer function)
        auto gamma = [](double u) {
            if (u <= 0.0031308) return 12.92 * u;
            return 1.055 * std::pow(u, 1.0/2.4) - 0.055;
        };
        
        RGB rgb;
        rgb.r = std::clamp(gamma(r_lin), 0.0, 1.0);
        rgb.g = std::clamp(gamma(g_lin), 0.0, 1.0);
        rgb.b = std::clamp(gamma(b_lin), 0.0, 1.0);
        
        return rgb;
    }
};

} // namespace sirius::spectral
