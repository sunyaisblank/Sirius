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

#include "../Metric/PHMT000B.h"
#include "../Metric/PHMT100B.h"
#include "../Spectral/MTSB001A.h"
#include "PHAD000A.h"
#include <PHCN001A.h>
#include <cmath>

namespace Sirius {

//==============================================================================
// AccretionDiskD: Novikov-Thorne thin disk model
// Implements: Sirius::IDiskModel interface for polymorphic disk rendering
//==============================================================================

class AccretionDiskD : public Sirius::IDiskModel {
public:
    struct Config {
        double M = 1.0;                     // Black hole mass [M_sun]
        double a_star = 0.0;                // Dimensionless spin a/M ∈ [-1, 1]
        double Mdot = 1e-8;                 // Accretion rate [M_sun/year]
        double r_inner = 0;                 // Inner edge (0 = ISCO, auto-computed)
        double r_outer = 500;               // Outer edge [GM/c²]
        double inclination = M_PI/4;        // Disk inclination to observer [rad]

        void validate() {
            if (M <= 0) M = 1.0;
            a_star = std::clamp(a_star, -1.0, 1.0);
            if (Mdot <= 0) Mdot = 1e-8;
            if (r_outer <= r_inner) r_outer = r_inner + 100.0;
        }
    };
    
    AccretionDiskD() : m_config() { init(); }

    explicit AccretionDiskD(const Config& config)
        : m_config(config) { init(); }

private:
    void init() {
        m_config.validate();
        // Compute derived quantities
        m_M_kg = m_config.M * Constants::Physical::M_SUN;
        m_GM = Constants::Physical::G_NEWTON * m_M_kg;
        m_rs = 2 * m_GM / (Constants::Physical::c_LIGHT * Constants::Physical::c_LIGHT);  // Schwarzschild radius
        m_a = m_config.a_star;  // a/M in geometric units

        // Compute ISCO radius
        m_r_isco = computeISCO(m_a);

        // Set inner edge
        m_r_inner = (m_config.r_inner > 0) ? m_config.r_inner : m_r_isco;
        m_r_outer = m_config.r_outer;

        // Convert accretion rate to SI
        m_Mdot_SI = m_config.Mdot * Constants::Physical::M_SUN / (365.25 * 24 * 3600);  // M_sun/yr → kg/s
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
    // Specific Energy for Circular Orbit (Page & Thorne 1974, Eq. 15.11)
    // E(r) = (r² - 2Mr + a√(Mr)) / (r √(r² - 3Mr + 2a√(Mr)))
    //--------------------------------------------------------------------------

    static double specificEnergy(double r, double M, double a) {
        double sqrtMr = std::sqrt(M * r);
        double numerator = r*r - 2*M*r + a * sqrtMr;
        double denominator = r * std::sqrt(r*r - 3*M*r + 2*a*sqrtMr);
        if (std::abs(denominator) < Constants::Tolerances::DIVISION_SAFE_EPS) return 1.0;
        return numerator / denominator;
    }

    //--------------------------------------------------------------------------
    // Specific Angular Momentum for Circular Orbit (Page & Thorne 1974, Eq. 15.12)
    // L(r) = √M (r² - 2a√(Mr) + a²) / (r √(r² - 3Mr + 2a√(Mr)))
    //--------------------------------------------------------------------------

    static double specificAngularMomentum(double r, double M, double a) {
        double sqrtM = std::sqrt(M);
        double sqrtMr = std::sqrt(M * r);
        double numerator = sqrtM * (r*r - 2*a*sqrtMr + a*a);
        double denominator = r * std::sqrt(r*r - 3*M*r + 2*a*sqrtMr);
        if (std::abs(denominator) < Constants::Tolerances::DIVISION_SAFE_EPS) return 0.0;
        return numerator / denominator;
    }

    //--------------------------------------------------------------------------
    // Angular Velocity for Circular Orbit (Bardeen 1970)
    // Ω(r) = √M / (r^(3/2) + a√M)
    //--------------------------------------------------------------------------

    static double angularVelocity(double r, double M, double a) {
        double sqrtM = std::sqrt(M);
        return sqrtM / (std::pow(r, 1.5) + a * sqrtM);
    }

    //--------------------------------------------------------------------------
    // dL/dr - ANALYTIC Derivative of Specific Angular Momentum
    // Required for Page-Thorne flux integral
    //
    // MATHEMATICAL DERIVATION:
    // L(r) = √M × N / (r × √D)  where:
    //   N = r² - 2a√(Mr) + a²
    //   D = r² - 3Mr + 2a√(Mr)
    //
    // Using quotient rule:
    // dL/dr = √M × [dN/dr × (r√D) - N × d(r√D)/dr] / (r√D)²
    //
    // where:
    //   dN/dr = 2r - a√(M/r)
    //   dD/dr = 2r - 3M + a√(M/r)
    //   d(r√D)/dr = √D + r × dD/dr / (2√D) = (2D + r×dD/dr) / (2√D)
    //
    // Reference: Derived analytically from Page & Thorne (1974) Eq. 15.12
    //--------------------------------------------------------------------------

    static double dLdr(double r, double M, double a) {
        if (r <= 0 || M <= 0) return 0;

        double sqrtM = std::sqrt(M);
        double sqrtMr = std::sqrt(M * r);
        double sqrtR = std::sqrt(r);

        // Numerator and denominator inner expressions
        double N = r*r - 2*a*sqrtMr + a*a;           // r² - 2a√(Mr) + a²
        double D = r*r - 3*M*r + 2*a*sqrtMr;         // r² - 3Mr + 2a√(Mr)

        // Protect against D ≤ 0 (inside ISCO)
        if (D <= Constants::Tolerances::DIVISION_SAFE_EPS) return 0;
        double sqrtD = std::sqrt(D);

        // Derivatives of N and D
        double aSqrtMoverR = a * sqrtM / sqrtR;      // a√(M/r)
        double dN_dr = 2*r - aSqrtMoverR;            // 2r - a√(M/r)
        double dD_dr = 2*r - 3*M + aSqrtMoverR;      // 2r - 3M + a√(M/r)

        // d(r√D)/dr = (2D + r·dD/dr) / (2√D)
        double d_rsqrtD_dr = (2*D + r*dD_dr) / (2*sqrtD);

        // Denominator: (r√D)²
        double rSqrtD = r * sqrtD;
        double denom = rSqrtD * rSqrtD;
        if (std::abs(denom) < Constants::Tolerances::FLUX_INTEGRAL_EPS) return 0;

        // Quotient rule: dL/dr = √M × [dN/dr × (r√D) - N × d(r√D)/dr] / (r√D)²
        double result = sqrtM * (dN_dr * rSqrtD - N * d_rsqrtD_dr) / denom;

        return result;
    }

    //--------------------------------------------------------------------------
    // dΩ/dr - ANALYTIC Derivative of Angular Velocity
    //
    // MATHEMATICAL DERIVATION:
    // Ω(r) = √M / (r^(3/2) + a√M)
    //
    // Using chain rule:
    // dΩ/dr = √M × (-(3/2)r^(1/2)) / (r^(3/2) + a√M)²
    //       = -3√(Mr) / [2(r^(3/2) + a√M)²]
    //
    // Reference: Derived analytically from Bardeen (1970)
    //--------------------------------------------------------------------------

    static double dOmegadr(double r, double M, double a) {
        if (r <= 0 || M <= 0) return 0;

        double sqrtM = std::sqrt(M);
        double sqrtMr = std::sqrt(M * r);
        double r32 = std::pow(r, 1.5);              // r^(3/2)
        double denom = r32 + a * sqrtM;             // r^(3/2) + a√M

        if (std::abs(denom) < Constants::Tolerances::FLUX_INTEGRAL_EPS) return 0;

        // dΩ/dr = -3√(Mr) / [2(r^(3/2) + a√M)²]
        return -3.0 * sqrtMr / (2.0 * denom * denom);
    }

    //--------------------------------------------------------------------------
    // Full Page & Thorne (1974) Relativistic Flux
    // F(r) = (Ṁ / 4π√(-g)) × (-Ω,r) / (E - ΩL)² × ∫[r_ISCO to r] (E - ΩL) L,r dr
    //
    // Simplified form using Q-factor:
    // F(r) = (3GMṀ / 8πr³) × Q(r)
    // where Q(r) encodes all relativistic corrections
    //--------------------------------------------------------------------------

    double fullPageThorneFlux(double r) const {
        if (r <= m_r_isco || r > m_r_outer) return 0;

        double M = 1.0;  // Working in units of M=1
        double a = m_a;

        // Get orbital quantities at radius r
        double E_r = specificEnergy(r, M, a);
        double L_r = specificAngularMomentum(r, M, a);
        double Omega_r = angularVelocity(r, M, a);

        // Compute the Page-Thorne flux integral from r_ISCO to r
        // ∫[r_ISCO to r] (E - ΩL) dL/dr' dr'
        //
        // 16-point Gauss-Legendre quadrature on [r_ISCO, r].
        // Nodes are interior to the interval, avoiding the ISCO boundary
        // where the integrand is singular.
        //
        // Nodes and weights for [-1, 1] (Abramowitz & Stegun, Table 25.4):
        static constexpr double gl_nodes[16] = {
            -0.9894009349916499, -0.9445750230732326,
            -0.8656312023878318, -0.7554044083550030,
            -0.6178762444026438, -0.4580167776572274,
            -0.2816035507792589, -0.0950125098376374,
             0.0950125098376374,  0.2816035507792589,
             0.4580167776572274,  0.6178762444026438,
             0.7554044083550030,  0.8656312023878318,
             0.9445750230732326,  0.9894009349916499
        };
        static constexpr double gl_weights[16] = {
            0.0271524594117541, 0.0622535239386479,
            0.0951585116824928, 0.1246289712555339,
            0.1495959888165767, 0.1691565193950025,
            0.1826034150449236, 0.1894506104550685,
            0.1894506104550685, 0.1826034150449236,
            0.1691565193950025, 0.1495959888165767,
            0.1246289712555339, 0.0951585116824928,
            0.0622535239386479, 0.0271524594117541
        };

        // Map [-1, 1] → [r_ISCO, r]
        double half_width = 0.5 * (r - m_r_isco);
        double midpoint = 0.5 * (r + m_r_isco);

        double integral = 0;
        for (int i = 0; i < 16; ++i) {
            double r_i = midpoint + half_width * gl_nodes[i];

            double E_i = specificEnergy(r_i, M, a);
            double L_i = specificAngularMomentum(r_i, M, a);
            double Omega_i = angularVelocity(r_i, M, a);
            double dL_i = dLdr(r_i, M, a);

            double integrand = (E_i - Omega_i * L_i) * dL_i;
            integral += gl_weights[i] * integrand;
        }
        integral *= half_width;

        // Compute (E - ΩL)² at r
        double E_minus_OmegaL = E_r - Omega_r * L_r;
        if (std::abs(E_minus_OmegaL) < Constants::Tolerances::DIVISION_SAFE_EPS) return 0;

        // Compute dΩ/dr at r using ANALYTIC derivative
        double dOmega_dr = dOmegadr(r, M, a);

        // Page-Thorne Q-factor (relativistic correction)
        // Q = (-Ω,r / (E - ΩL)²) × ∫(E - ΩL) L,r dr
        double Q = (-dOmega_dr / (E_minus_OmegaL * E_minus_OmegaL)) * integral;

        // Classical Newtonian flux prefactor in SI
        double r_physical = r * m_rs / 2;  // r in meters
        double F_newt = (3 * m_GM * m_Mdot_SI) / (8 * M_PI * r_physical * r_physical * r_physical);

        // Apply Page-Thorne correction
        double F = F_newt * Q;

        // Validate output
        if (!std::isfinite(F) || F < 0) return 0;

        return F;
    }

    //--------------------------------------------------------------------------
    // Radiative Flux (Page & Thorne 1974) - Implements IDiskModel::flux
    // F(r) = (3GMṀ)/(8πr³) × Q(r) where Q(r) is the full relativistic factor
    //--------------------------------------------------------------------------

    double flux(double r) const override {
        if (r < m_r_inner || r > m_r_outer) return 0;

        double r_M = r;  // r is in units of M

        // Must be outside ISCO for emission
        if (r_M <= m_r_isco) return 0;

        // Inner torque factor: ensures F → 0 at ISCO
        double sqrt_ratio = std::sqrt(m_r_isco / r_M);
        double inner_torque = 1 - sqrt_ratio;
        if (inner_torque <= 0) return 0;

        // Convert r to physical units
        double r_physical = r_M * m_rs / 2;  // r in meters (rs = 2GM/c²)

        // Newtonian flux prefactor
        double F_newt = (3 * m_GM * m_Mdot_SI) / (8 * M_PI * r_physical * r_physical * r_physical);

        // Try full Page-Thorne calculation for r > 1.5 * r_isco
        if (r_M > m_r_isco * 1.5) {
            double F_PT = fullPageThorneFlux(r);
            if (std::isfinite(F_PT) && F_PT > 0) {
                return F_PT;
            }
        }

        // Fallback: simplified relativistic correction
        double y = std::sqrt(r_M);
        double y3 = y * y * y;

        // Q ≈ 1 / (y³(y - 3/y + 2a/y²)) for simplified Page & Thorne
        double denom = y3 * (y - 3.0/y + 2*m_a/(y*y));
        double Q = (std::abs(denom) > 1e-10) ? 1.0 / denom : 1.0;

        // Clamp Q to reasonable range
        if (!std::isfinite(Q) || Q <= 0) Q = 1.0;
        Q = std::clamp(Q, 0.1, 10.0);

        return F_newt * inner_torque * Q;
    }
    
    //--------------------------------------------------------------------------
    // Temperature Profile (implements IDiskModel::temperature)
    // T(r) = [F(r) / σ_SB]^(1/4)
    //--------------------------------------------------------------------------

    double temperature(double r, [[maybe_unused]] double z = 0) const override {
        double F = flux(r);
        if (F <= 0) return 0;

        return std::pow(F / Constants::Physical::SIGMA_SB, 0.25);
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
    // Check if Point is Within Disk - Implements IDiskModel::isInDisk
    //--------------------------------------------------------------------------

    bool isInDisk(double r, double theta) const override {
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

    const Config& config() const { return m_config; }
    double iscoRadius() const { return m_r_isco; }
    double spinParameter() const { return m_a; }
    double massKg() const { return m_M_kg; }
    double schwarzschildRadius() const { return m_rs; }

    // Alias for temperature() for API compatibility
    double effectiveTemperature(double r) const { return temperature(r); }

    //==========================================================================
    // IDiskModel Interface Implementation
    //==========================================================================

    /// @brief Get model name
    const char* modelName() const override { return "Novikov-Thorne"; }

    /// @brief Inner radius in units of M (implements IDiskModel)
    double innerRadius() const override { return m_r_inner; }

    /// @brief Outer radius in units of M (implements IDiskModel)
    double outerRadius() const override { return m_r_outer; }

    /// @brief Half-thickness at radius (thin disk: H = 0)
    double halfThickness([[maybe_unused]] double r) const override {
        return 0.0;  // Novikov-Thorne is infinitely thin
    }

    /// @brief Density at location (thin disk approximation)
    double density([[maybe_unused]] double r, [[maybe_unused]] double z = 0) const override {
        // Thin disk: surface density varies, volume density is delta-function
        // Return normalized value for optical depth calculations
        return 1.0;
    }

    /// @brief Angular velocity for IDiskModel interface
    double angularVelocity(double r) const override {
        return angularVelocity(r, 1.0, m_a);  // Use static version with M=1
    }

    /// @brief Set black hole parameters (implements IDiskModel)
    void setBlackHoleParameters(double mass, double spin) override {
        m_config.M = mass;
        m_config.a_star = spin;
        init();  // Recompute derived quantities
    }

    /// @brief Get black hole mass (implements IDiskModel)
    double blackHoleMass() const override { return m_config.M; }

    /// @brief Get black hole spin (implements IDiskModel)
    double blackHoleSpin() const override { return m_config.a_star; }

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

} // namespace Sirius
