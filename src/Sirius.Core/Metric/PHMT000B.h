// PHMT000B.h - Double-Precision Metric Interface
// Component ID: PHMT000B
// Purpose: Abstract base class for all double-precision spacetime metrics
//
// MATHEMATICAL BASIS:
// The metric tensor g_μν defines spacetime geometry and must satisfy:
// - Symmetry: g_μν = g_νμ
// - Non-degeneracy: det(g) ≠ 0
// - Lorentzian signature: (-,+,+,+)
//
// STRICT REQUIREMENTS:
// - det(g) must be negative (timelike signature)
// - Inverse must satisfy g^μα g_αν = δ^μ_ν to <10^-14
// - Christoffel symmetry: Γ^μ_νρ = Γ^μ_ρν
//
// REFERENCE: James et al. (2015) "Gravitational Lensing by Spinning Black Holes"

#pragma once

#include "../Transport/MTTP001A.h"

namespace Sirius {

//==============================================================================
// IMetricD: Abstract interface for double-precision metric evaluation
//==============================================================================

class IMetricD {
public:
    virtual ~IMetricD() = default;
    
    //--------------------------------------------------------------------------
    // Core Metric Evaluation
    //--------------------------------------------------------------------------
    
    /// Evaluate metric tensor g_μν and its inverse g^μν at position x
    /// @param x Position in Boyer-Lindquist coordinates (t, r, θ, φ)
    /// @param g Output: covariant metric tensor g[μ][ν]
    /// @param g_inv Output: contravariant metric tensor g^μ^ν
    /// @pre x must be in valid coordinate range (r > r_horizon, 0 < θ < π)
    /// @post g and g_inv are symmetric
    /// @post g^μα g_αν = δ^μ_ν within tolerance 10^-14
    virtual void evaluate(
        const Vec4d& x,
        double g[4][4],
        double g_inv[4][4]
    ) const = 0;
    
    //--------------------------------------------------------------------------
    // Connection Coefficients
    //--------------------------------------------------------------------------
    
    /// Christoffel symbols Γ^μ_νρ = (1/2) g^μσ (∂_ν g_σρ + ∂_ρ g_νσ - ∂_σ g_νρ)
    /// @param x Position in Boyer-Lindquist coordinates
    /// @param Gamma Output: Gamma[μ][ν][ρ] with index order (upper, lower, lower)
    /// @post Gamma[μ][ν][ρ] = Gamma[μ][ρ][ν] (symmetric in lower indices)
    virtual void christoffel(
        const Vec4d& x,
        double Gamma[4][4][4]
    ) const = 0;
    
    //--------------------------------------------------------------------------
    // Curvature (Required for Geodesic Deviation / Beam Propagation)
    //--------------------------------------------------------------------------
    
    /// Riemann curvature tensor R^μ_νρσ
    /// Required for beam propagation via geodesic deviation equation:
    ///   D²ξ^μ/dλ² + R^μ_νρσ k^ν ξ^ρ k^σ = 0
    /// @param x Position in Boyer-Lindquist coordinates
    /// @param R Output: R[μ][ν][ρ][σ] (first index up, rest down)
    /// @post Standard Riemann symmetries hold:
    ///       R_μνρσ = -R_νμρσ, R_μνρσ = -R_μνσρ, R_μνρσ = R_ρσμν
    virtual void riemann(
        const Vec4d& x,
        double R[4][4][4][4]
    ) const = 0;
    
    //--------------------------------------------------------------------------
    // Hamiltonian for Symplectic Integration
    //--------------------------------------------------------------------------
    
    /// Compute Hamiltonian H = (1/2) g^μν p_μ p_ν
    /// For null geodesics (photons): H = 0
    /// For timelike geodesics: H = -1/2
    /// @param q Generalised coordinates (position)
    /// @param p Conjugate momenta (covariant wave vector)
    /// @return Hamiltonian value (should be ~0 for null rays, < 10^-10)
    virtual double hamiltonian(
        const Vec4d& q,
        const Vec4d& p
    ) const = 0;
    
    /// Compute ∂H/∂p (velocity from momentum)
    /// dq^μ/dλ = ∂H/∂p_μ = g^μν p_ν
    /// @param q Position
    /// @param p Momentum
    /// @return Contravariant velocity dx/dλ
    virtual Vec4d dHdp(const Vec4d& q, const Vec4d& p) const = 0;
    
    /// Compute ∂H/∂q (momentum change from position)
    /// dp_μ/dλ = -∂H/∂q^μ = -(1/2) (∂g^νρ/∂q^μ) p_ν p_ρ
    /// @param q Position
    /// @param p Momentum
    /// @return Rate of change of covariant momentum
    virtual Vec4d dHdq(const Vec4d& q, const Vec4d& p) const = 0;
    
    //--------------------------------------------------------------------------
    // Coordinate Validity and Boundaries
    //--------------------------------------------------------------------------
    
    /// Check if coordinates are in valid range
    /// @return false for: inside horizon, at coordinate singularity, θ = 0 or π
    virtual bool isValid(const Vec4d& x) const = 0;
    
    /// Return event horizon radius (outer horizon r+ for Kerr)
    virtual double horizonRadius() const = 0;
    
    /// Return inner horizon radius (r- for Kerr, 0 for Schwarzschild)
    virtual double innerHorizonRadius() const = 0;
    
    /// Return ergosphere radius at given theta
    virtual double ergosphereRadius(double theta) const = 0;
    
    /// Return ISCO (innermost stable circular orbit) radius
    /// For prograde orbits in Kerr: r_ISCO → M as a → M
    virtual double iscoRadius() const = 0;
    
    /// Return photon sphere radius (unstable circular photon orbit)
    virtual double photonSphereRadius() const = 0;
    
    //--------------------------------------------------------------------------
    // Metric Parameters
    //--------------------------------------------------------------------------
    
    /// Black hole mass parameter M (geometric units: G = c = 1)
    virtual double mass() const = 0;
    
    /// Spin parameter a = J/M (|a| ≤ M, returns 0 for non-rotating)
    virtual double spin() const = 0;
    
    /// Charge parameter Q (returns 0 for uncharged)
    virtual double charge() const = 0;
    
    //--------------------------------------------------------------------------
    // Utility Functions
    //--------------------------------------------------------------------------
    
    /// Time-transformation function g(q) for TTESI regularisation
    /// Default: g(q) = Σ = r² + a² cos²θ (Kerr)
    /// Used to regularise integration near horizon
    virtual double timeTransformationFunction(const Vec4d& q) const = 0;
    
    /// Compute g-factor (frequency ratio) for a photon with momentum k
    /// observed by a static observer at position x
    /// g = (k · u_obs) / (k · u_emit)
    /// For disk emission: g encodes gravitational + Doppler shift
    virtual double gFactor(
        const Vec4d& x,
        const Vec4d& k,
        const Vec4d& u_emitter
    ) const = 0;
};

//==============================================================================
// MetricParamsD: Double-precision metric parameters
//==============================================================================

struct MetricParamsD {
    double M;       // Mass
    double a;       // Spin parameter (|a| ≤ M)
    double Q;       // Charge parameter
    double Lambda;  // Cosmological constant
    
    // Derived quantities (computed from parameters)
    double rs;      // Schwarzschild radius 2M
    double rplus;   // Outer horizon radius
    double rminus;  // Inner horizon radius
    
    MetricParamsD() : M(1.0), a(0), Q(0), Lambda(0), rs(2.0), rplus(2.0), rminus(0) {}
    
    MetricParamsD(double mass, double spin = 0, double charge = 0)
        : M(mass), a(spin), Q(charge), Lambda(0) {
        rs = 2.0 * M;
        // r± = M ± √(M² - a² - Q²)
        double discriminant = M*M - a*a - Q*Q;
        if (discriminant >= 0) {
            rplus = M + std::sqrt(discriminant);
            rminus = M - std::sqrt(discriminant);
        } else {
            // Naked singularity (invalid for physical black holes)
            rplus = rminus = 0;
        }
    }
};

} // namespace Sirius
