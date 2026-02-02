// =============================================================================
// PHCN001A.h - Unified Numerical Constants and Tolerances
// Component ID: PHCN001A (Physics/Constants)
// =============================================================================
//
// PURPOSE:
// Centralizes all numerical tolerances, physical constants, and precision
// requirements for the Sirius ray tracing engine. All tolerances align with
// the requirements specified in docs/specification.md.
//
// DESIGN PHILOSOPHY:
// "A result is trustworthy only if the mathematical invariants of the
// underlying physics are preserved." - docs/philosophy.md
//
// PRECISION TIERS:
// - GPU (single precision): 23 significand bits, ε ≈ 1.2×10⁻⁷
// - CPU (double precision): 52 significand bits, ε ≈ 2.2×10⁻¹⁶
//
// TESTS: All tolerance values validated by TSDG tests
// =============================================================================

#ifndef PHCN001A_H
#define PHCN001A_H

#include <cmath>
#include <limits>

namespace Sirius {
namespace Constants {

// =============================================================================
// Machine Precision Constants
// =============================================================================

/// @brief Machine epsilon for single precision
constexpr float EPSILON_F = std::numeric_limits<float>::epsilon();  // ~1.19e-7

/// @brief Machine epsilon for double precision
constexpr double EPSILON_D = std::numeric_limits<double>::epsilon();  // ~2.22e-16

/// @brief Safe minimum for avoiding division by zero (float)
constexpr float SAFE_MIN_F = 1e-10f;

/// @brief Safe minimum for avoiding division by zero (double)
constexpr double SAFE_MIN_D = 1e-20;

// =============================================================================
// Geodesic Integration Tolerances
// Reference: docs/specification.md Section "Error Bounds"
// =============================================================================

namespace Geodesic {

/// @brief Null condition tolerance: g_μν k^μ k^ν = 0
/// GPU paths use relaxed tolerance due to single precision accumulation
constexpr float NULL_CONDITION_TOL_GPU = 1e-5f;

/// @brief Null condition tolerance for CPU/double paths
constexpr double NULL_CONDITION_TOL_CPU = 1e-10;

/// @brief Conservation law relative tolerance (E, L_z drift)
/// Per specification: < 10^-4 relative drift over integration
constexpr double CONSERVATION_TOL = 1e-4;

/// @brief Hamiltonian constraint tolerance for symplectic integrators
/// H = (1/2) g^μν p_μ p_ν should equal 0 for null geodesics
constexpr double HAMILTONIAN_TOL = 1e-10;

/// @brief Maximum allowed integration steps per ray
constexpr int MAX_INTEGRATION_STEPS = 10000;

/// @brief Default initial step size (affine parameter units)
constexpr float DEFAULT_STEP_SIZE = 0.01f;

/// @brief Minimum step size before declaring numerical failure
constexpr float MIN_STEP_SIZE = 1e-8f;

/// @brief Maximum step size for stability
constexpr float MAX_STEP_SIZE = 1.0f;

/// @brief Adaptive step safety factor (conservative)
constexpr float STEP_SAFETY_FACTOR = 0.9f;

} // namespace Geodesic

// =============================================================================
// Coordinate Domain Tolerances
// Reference: docs/specification.md Section "Domain Constraints"
// =============================================================================

namespace Coordinates {

/// @brief Horizon buffer: minimum r/r_+ ratio for valid integration
/// Rays crossing inside this buffer are terminated as "horizon crossing"
constexpr double HORIZON_BUFFER = 1.001;

/// @brief Pole avoidance: minimum |sin(θ)| before clamping
/// Prevents coordinate singularity at θ = 0, π
constexpr double POLE_EPSILON = 1e-6;

/// @brief Maximum radial coordinate before declaring "escaped"
/// In units of M (gravitational radius)
constexpr double ESCAPE_RADIUS = 1e6;

/// @brief Coordinate comparison tolerance
constexpr double COORDINATE_TOL = 1e-12;

/// @brief Azimuthal angle wrap tolerance
constexpr double PHI_WRAP_TOL = 1e-10;

} // namespace Coordinates

// =============================================================================
// Metric Tensor Tolerances
// Reference: docs/specification.md Section "Invariants"
// =============================================================================

namespace Metric {

/// @brief Metric symmetry tolerance: max |g_μν - g_νμ|
constexpr double SYMMETRY_TOL = 1e-15;

/// @brief Inverse accuracy: max |g^μα g_αν - δ^μ_ν|
constexpr double INVERSE_TOL = 1e-14;

/// @brief Christoffel symmetry tolerance: max |Γ^λ_μν - Γ^λ_νμ|
constexpr double CHRISTOFFEL_SYMMETRY_TOL = 1e-15;

/// @brief Determinant non-degeneracy threshold
constexpr double DETERMINANT_TOL = 1e-30;

/// @brief Lorentzian signature check threshold
constexpr double SIGNATURE_TOL = 1e-10;

} // namespace Metric

// =============================================================================
// Numerical Differentiation Tolerances
// Used only when analytic derivatives unavailable
// =============================================================================

namespace Differentiation {

/// @brief Relative step size for central differences
/// h = RELATIVE_H * |x| + ABSOLUTE_H
constexpr double RELATIVE_H = 1e-6;

/// @brief Absolute step size for near-zero x values
constexpr double ABSOLUTE_H = 1e-10;

/// @brief Richardson extrapolation convergence threshold
constexpr double RICHARDSON_TOL = 1e-10;

} // namespace Differentiation

// =============================================================================
// Spectral and Radiative Transfer
// =============================================================================

namespace Spectral {

/// @brief Visible spectrum range (meters)
constexpr double LAMBDA_MIN = 380e-9;
constexpr double LAMBDA_MAX = 780e-9;

/// @brief Number of spectral bins
constexpr int N_SPECTRAL_BINS = 32;

/// @brief Minimum intensity for non-zero emission
constexpr double INTENSITY_FLOOR = 1e-30;

/// @brief Maximum allowed redshift factor before clamping
constexpr double MAX_G_FACTOR = 100.0;
constexpr double MIN_G_FACTOR = 0.01;

/// @brief Limb darkening coefficient bounds
constexpr double LIMB_DARK_MIN = 0.0;
constexpr double LIMB_DARK_MAX = 1.0;

} // namespace Spectral

// =============================================================================
// Disk Model Tolerances
// =============================================================================

namespace Disk {

/// @brief Equatorial plane thickness for thin disk intersection
/// |θ - π/2| < EQUATORIAL_TOL defines the disk plane
constexpr double EQUATORIAL_TOL = 0.01;  // ~0.5 degrees

/// @brief Minimum radius for disk emission (relative to ISCO)
constexpr double INNER_EDGE_BUFFER = 1.001;

/// @brief Simpson's rule integration points for Page-Thorne Q-factor
constexpr int INTEGRATION_POINTS = 100;

/// @brief Q-factor validity range (clamp outliers)
constexpr double Q_FACTOR_MIN = 0.01;
constexpr double Q_FACTOR_MAX = 100.0;

} // namespace Disk

// =============================================================================
// Physical Constants (SI Units)
// Reference: CODATA 2018
// =============================================================================

namespace Physical {

/// @brief Planck constant [J·s]
constexpr double h_PLANCK = 6.62607015e-34;

/// @brief Speed of light [m/s]
constexpr double c_LIGHT = 2.99792458e8;

/// @brief Boltzmann constant [J/K]
constexpr double k_BOLTZMANN = 1.380649e-23;

/// @brief Stefan-Boltzmann constant [W/(m²·K⁴)]
constexpr double SIGMA_SB = 5.670374419e-8;

/// @brief Gravitational constant [m³/(kg·s²)]
constexpr double G_NEWTON = 6.67430e-11;

/// @brief Solar mass [kg]
constexpr double M_SUN = 1.98892e30;

/// @brief Parsec [m]
constexpr double PARSEC = 3.08567758e16;

/// @brief Wien displacement constant [m·K]
constexpr double WIEN_B = 2.897771955e-3;

// Derived constants
constexpr double PLANCK_C1 = 2.0 * h_PLANCK * c_LIGHT * c_LIGHT;  // 2hc²
constexpr double PLANCK_C2 = h_PLANCK * c_LIGHT / k_BOLTZMANN;    // hc/k
constexpr double GM_SUN = G_NEWTON * M_SUN;
constexpr double RS_SUN = 2.0 * GM_SUN / (c_LIGHT * c_LIGHT);  // ~2953 m

} // namespace Physical

// =============================================================================
// Mathematical Constants
// =============================================================================

namespace Math {

constexpr double PI = 3.14159265358979323846;
constexpr double TWO_PI = 2.0 * PI;
constexpr double HALF_PI = PI / 2.0;
constexpr double SQRT_2 = 1.41421356237309504880;
constexpr double SQRT_3 = 1.73205080756887729353;
constexpr double INV_SQRT_2 = 0.70710678118654752440;

/// @brief Cube root of 2 (used in Yoshida coefficients)
constexpr double CBRT_2 = 1.25992104989487316477;

} // namespace Math

// =============================================================================
// Utility Functions
// =============================================================================

/// @brief Select appropriate tolerance based on precision type
template<typename T>
constexpr T nullConditionTolerance() {
    if constexpr (std::is_same_v<T, float>) {
        return static_cast<T>(Geodesic::NULL_CONDITION_TOL_GPU);
    } else {
        return static_cast<T>(Geodesic::NULL_CONDITION_TOL_CPU);
    }
}

/// @brief Safe division with floor on denominator
template<typename T>
inline T safeDivide(T numerator, T denominator, T floor = SAFE_MIN_D) {
    if (std::abs(denominator) < floor) {
        denominator = std::copysign(floor, denominator);
    }
    return numerator / denominator;
}

/// @brief Clamp value to valid coordinate range
inline double clampCoordinate(double theta, double epsilon = Coordinates::POLE_EPSILON) {
    return std::clamp(theta, epsilon, Math::PI - epsilon);
}

/// @brief Wrap azimuthal angle to [0, 2π)
inline double wrapPhi(double phi) {
    while (phi < 0) phi += Math::TWO_PI;
    while (phi >= Math::TWO_PI) phi -= Math::TWO_PI;
    return phi;
}

} // namespace Constants
} // namespace Sirius

#endif // PHCN001A_H
