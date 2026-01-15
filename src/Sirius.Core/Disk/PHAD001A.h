// PHAD001A.h - Accretion Disk Physics
// Component ID: PHAD001A
// Purpose: Novikov-Thorne thin disk model for Kerr black holes
//
// MATHEMATICAL BASIS:
// The Novikov-Thorne model describes a geometrically thin, optically thick
// accretion disk around a Kerr black hole with radiative efficiency:
//
//   η = 1 - E(r_ISCO) where E is specific energy at ISCO
//
// Temperature profile (Stefan-Boltzmann):
//   T(r) = [3GMṀ/(8πσr³) × f(r)]^(1/4)
//
// where f(r) is the Page & Thorne (1974) relativistic correction factor.
//
// REFERENCES:
// - Novikov & Thorne (1973) "Astrophysics of Black Holes"
// - Page & Thorne (1974) "Disk-Accretion onto a Black Hole"
// - James et al. (2015) "DNGR" Section 3.2

#pragma once

#include "../Geodesic/PHMT000B.h"
#include "../Metric/PHMT100B.h"
#include "../Symplectic/MTSB001A.h"
#include <cmath>

namespace sirius::physics {

using namespace sirius::spectral;

//==============================================================================
// Physical Constants (SI units)
//==============================================================================

namespace constants {

constexpr double G = 6.67430e-11;       // Gravitational constant [m³/(kg·s²)]
constexpr double c = 2.99792458e8;      // Speed of light [m/s]
constexpr double sigma_SB = 5.670374e-8; // Stefan-Boltzmann [W/(m²·K⁴)]
constexpr double M_sun = 1.98892e30;    // Solar mass [kg]
constexpr double pc = 3.08567758e16;    // Parsec [m]

// Derived
constexpr double GM_sun = G * M_sun;    // GM for 1 solar mass
constexpr double rs_sun = 2 * GM_sun / (c * c);  // Schwarzschild radius for 1 M_sun [m]

} // namespace constants

//==============================================================================
// AccretionDiskD: Novikov-Thorne thin disk model
//==============================================================================

class AccretionDiskD {
public:
    struct Config {
        double M = 1.0;                     // Black hole mass [M_sun]
        double a_star = 0.0;                // Dimensionless spin a/M ∈ [-1, 1]
        double Mdot = 1e-8;                 // Accretion rate [M_sun/year]
        double r_inner = 0;                 // Inner edge (0 = ISCO, auto-computed)
        double r_outer = 500;               // Outer edge [GM/c²]
        double inclination = M_PI/4;        // Disk inclination to observer [rad]
    };
    
    AccretionDiskD() : m_config() { init(); }

    explicit AccretionDiskD(const Config& config)
        : m_config(config) { init(); }

private:
    void init() {
        // Compute derived quantities
        m_M_kg = m_config.M * constants::M_sun;
        m_GM = constants::G * m_M_kg;
        m_rs = 2 * m_GM / (constants::c * constants::c);  // Schwarzschild radius
        m_a = m_config.a_star;  // a/M in geometric units

        // Compute ISCO radius
        m_r_isco = computeISCO(m_a);

        // Set inner edge
        m_r_inner = (m_config.r_inner > 0) ? m_config.r_inner : m_r_isco;
        m_r_outer = m_config.r_outer;

        // Convert accretion rate to SI
        m_Mdot_SI = m_config.Mdot * constants::M_sun / (365.25 * 24 * 3600);  // M_sun/yr → kg/s
    }

public:
    
    //--------------------------------------------------------------------------
    // ISCO Radius Computation
    // r_ISCO = M {3 + Z₂ ∓ √[(3 - Z₁)(3 + Z₁ + 2Z₂)]}
    // where ± is for prograde/retrograde orbits
    //--------------------------------------------------------------------------
    
    static double computeISCO(double a_star) {
        double a = std::abs(a_star);
        if (a > 1) a = 1;  // Clamp to extremal
        
        double Z1 = 1 + std::cbrt(1 - a*a) * (std::cbrt(1 + a) + std::cbrt(1 - a));
        double Z2 = std::sqrt(3*a*a + Z1*Z1);
        
        double r_isco;
        if (a_star >= 0) {
            // Prograde orbit
            r_isco = 3 + Z2 - std::sqrt((3 - Z1)*(3 + Z1 + 2*Z2));
        } else {
            // Retrograde orbit
            r_isco = 3 + Z2 + std::sqrt((3 - Z1)*(3 + Z1 + 2*Z2));
        }
        
        return r_isco;  // In units of M
    }
    
    //--------------------------------------------------------------------------
    // Radiative Flux (Page & Thorne 1974)
    // F(r) = (3GMṀ)/(8πr³) × f(r) where f(r) is the relativistic factor
    //--------------------------------------------------------------------------
    
    double flux(double r) const {
        if (r < m_r_inner || r > m_r_outer) return 0;
        
        double r_M = r;  // r is in units of M
        
        // Classical (Newtonian) flux prefactor
        // F_Newt = 3GMṀ / (8π r³) × [1 - √(r_in/r)]
        
        // For the relativistic correction, we use the simplified form:
        // f(r) ≈ [1 - √(r_isco/r)] × relativistic_correction
        //
        // Full Page & Thorne requires integrals, using simplified form for now
        
        double sqrt_ratio = std::sqrt(m_r_isco / r_M);
        double inner_torque = 1 - sqrt_ratio;
        
        if (inner_torque <= 0) return 0;
        
        // Convert r to physical units
        double r_physical = r_M * m_rs / 2;  // r in meters (rs = 2GM/c²)
        
        // Newtonian flux formula
        double F = (3 * m_GM * m_Mdot_SI) / (8 * M_PI * r_physical * r_physical * r_physical);
        F *= inner_torque;
        
        // Relativistic corrections (approximate)
        // Full treatment requires E(r), L(r), Ω(r) from geodesic equations
        double y = std::sqrt(r_M);
        double y_isco = std::sqrt(m_r_isco);
        
        // Correction factor (simplified approximation)
        double correction = 1.0 / (y * y * y * (y - 3.0/y + 2*m_a/(y*y)));
        if (!std::isfinite(correction) || correction <= 0) correction = 1;
        
        return F * std::min(correction, 10.0);  // Cap correction factor
    }
    
    //--------------------------------------------------------------------------
    // Temperature Profile
    // T(r) = [F(r) / σ_SB]^(1/4)
    //--------------------------------------------------------------------------
    
    double temperature(double r) const {
        double F = flux(r);
        if (F <= 0) return 0;
        
        return std::pow(F / constants::sigma_SB, 0.25);
    }
    
    //--------------------------------------------------------------------------
    // Peak Temperature (for normalisation)
    //--------------------------------------------------------------------------
    
    double peakTemperature() const {
        // Peak is typically around 1.5-2× ISCO radius
        double r_peak = 1.5 * m_r_isco;
        return temperature(r_peak);
    }
    
    //--------------------------------------------------------------------------
    // Spectral Radiance at Radius
    // Blackbody spectrum at local temperature
    //--------------------------------------------------------------------------
    
    SpectralRadiance emissionSpectrum(double r) const {
        double T = temperature(r);
        if (T <= 0) return SpectralRadiance::zero();
        
        return SpectralRadiance::blackbody(T);
    }
    
    //--------------------------------------------------------------------------
    // Limb Darkening
    // I(θ) = I₀ × (1 + u·cos(θ)) / (1 + u)
    // where u is the limb darkening coefficient
    //--------------------------------------------------------------------------
    
    static double limbDarkening(double cos_theta, double u = 0.6) {
        if (cos_theta <= 0) return 0;
        return (1 + u * cos_theta) / (1 + u);
    }
    
    //--------------------------------------------------------------------------
    // Check if Point is Within Disk
    //--------------------------------------------------------------------------
    
    bool isInDisk(double r, double theta) const {
        // Disk is in equatorial plane: |θ - π/2| < ε
        double equatorial_distance = std::abs(theta - M_PI/2);
        if (equatorial_distance > 0.01) return false;  // ~0.5° from equator
        
        return (r >= m_r_inner && r <= m_r_outer);
    }
    
    //--------------------------------------------------------------------------
    // Disk Intersection (for ray/beam)
    // Returns parameter λ at intersection, or -1 if no intersection
    //--------------------------------------------------------------------------
    
    double intersectRay(const Vec4d& x, const Vec4d& k, double lambda_max = 1000) const {
        // Simple: check when θ crosses π/2
        // More sophisticated would solve θ(λ) = π/2
        
        double theta0 = x.theta;
        double dtheta = k.theta;  // Approximate: k_θ ∝ dθ/dλ
        
        if (std::abs(dtheta) < 1e-10) return -1;  // Parallel to equator
        
        double lambda_cross = (M_PI/2 - theta0) / dtheta;
        
        if (lambda_cross < 0 || lambda_cross > lambda_max) return -1;
        
        // Estimate r at crossing
        double r_cross = x.r + k.r * lambda_cross;
        
        if (r_cross >= m_r_inner && r_cross <= m_r_outer) {
            return lambda_cross;
        }
        
        return -1;
    }
    
    //--------------------------------------------------------------------------
    // Accessors
    //--------------------------------------------------------------------------
    
    double iscoRadius() const { return m_r_isco; }
    double innerRadius() const { return m_r_inner; }
    double outerRadius() const { return m_r_outer; }
    double spinParameter() const { return m_a; }
    double massKg() const { return m_M_kg; }
    double schwarzschildRadius() const { return m_rs; }
    
private:
    Config m_config;
    
    // Derived quantities
    double m_M_kg;          // Mass in kg
    double m_GM;            // GM in SI
    double m_rs;            // Schwarzschild radius in meters
    double m_a;             // Dimensionless spin
    double m_r_isco;        // ISCO radius in units of M
    double m_r_inner;       // Inner edge in units of M
    double m_r_outer;       // Outer edge in units of M
    double m_Mdot_SI;       // Accretion rate in kg/s
};

} // namespace sirius::physics
