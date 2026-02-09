// =============================================================================
// TSFT001A.h - Unified Test Constants and Tolerances
// Component ID: TSFT001A (Test/Fixtures/Constants)
// =============================================================================
//
// PURPOSE:
// Centralizes all test tolerance constants by importing from PHCN001A.h.
// Ensures all tests use consistent, specification-compliant tolerances.
//
// DESIGN PRINCIPLE:
// "A result is trustworthy only if the mathematical invariants of the
// underlying physics are preserved." - docs/philosophy.md
//
// USAGE:
// All test files should include this header instead of defining local tolerances.
// This ensures consistency with docs/specification.md requirements.
//
// TESTS: Used by all diagnostic and unit tests
// =============================================================================

#ifndef TSFT001A_H
#define TSFT001A_H

#include <PHCN001A.h>
#include <cmath>

namespace sirius::test {
using namespace Sirius;

// =============================================================================
// Geodesic Integration Test Tolerances
// Reference: docs/specification.md Section "Error Bounds"
// =============================================================================

namespace Tolerances {

/// @brief Null condition tolerance for GPU paths
/// Specification: |g_μν k^μ k^ν| < 10^-5 for single precision
constexpr float NULL_CONDITION_GPU = Sirius::Constants::Geodesic::NULL_CONDITION_TOL_GPU;  // 1e-5

/// @brief Null condition tolerance for CPU paths
/// Specification: |g_μν k^μ k^ν| < 10^-10 for double precision
constexpr double NULL_CONDITION_CPU = Sirius::Constants::Geodesic::NULL_CONDITION_TOL_CPU;  // 1e-10

/// @brief Conservation law tolerance (energy, angular momentum)
/// Specification: |ΔE/E| < 10^-4, |ΔL/L| < 10^-4 relative drift
constexpr double CONSERVATION = Sirius::Constants::Geodesic::CONSERVATION_TOL;  // 1e-4

/// @brief Hamiltonian constraint tolerance (symplectic integrators)
/// Specification: |H| < 10^-10 for null geodesics
constexpr double HAMILTONIAN = Sirius::Constants::Geodesic::HAMILTONIAN_TOL;  // 1e-10

/// @brief Metric symmetry tolerance
/// Specification: max|g_μν - g_νμ| < 10^-15 (exact for analytic metrics)
constexpr double METRIC_SYMMETRY = Sirius::Constants::Metric::SYMMETRY_TOL;  // 1e-15

/// @brief Inverse metric accuracy
/// Specification: max|g^μα g_αν - δ^μ_ν| < 10^-14
constexpr double INVERSE_METRIC = Sirius::Constants::Metric::INVERSE_TOL;  // 1e-14

/// @brief Christoffel symmetry tolerance
/// Specification: max|Γ^λ_μν - Γ^λ_νμ| < 10^-15 (torsion-free connection)
constexpr double CHRISTOFFEL_SYMMETRY = Sirius::Constants::Metric::CHRISTOFFEL_SYMMETRY_TOL;  // 1e-15

/// @brief Metric determinant non-degeneracy
/// Specification: |det(g)| > 10^-30
constexpr double DETERMINANT = Sirius::Constants::Metric::DETERMINANT_TOL;  // 1e-30

/// @brief Lorentzian signature tolerance
/// Specification: eigenvalue signs must be (-,+,+,+) to 10^-10
constexpr double SIGNATURE = Sirius::Constants::Metric::SIGNATURE_TOL;  // 1e-10

/// @brief Coordinate comparison tolerance
constexpr double COORDINATE = Sirius::Constants::Coordinates::COORDINATE_TOL;  // 1e-12

/// @brief General numerical equality tolerance (double precision)
/// For tests where exact match is expected but floating-point errors exist
constexpr double GENERAL_DOUBLE = 1e-12;

/// @brief General numerical equality tolerance (single precision)
constexpr float GENERAL_FLOAT = 1e-6f;

/// @brief Weak field approximation tolerance (O(M/r) corrections)
/// For far-field asymptotic tests where 1/r² corrections are acceptable
constexpr double WEAK_FIELD = 1e-4;

/// @brief Analytic reference value tolerance
/// For comparison with known analytic solutions (textbook values)
constexpr double ANALYTIC_REFERENCE = 1e-10;

} // namespace Tolerances

// =============================================================================
// Physical Constants for Test Validation
// =============================================================================

namespace PhysicalRef {

/// @brief Schwarzschild photon sphere radius
constexpr double PHOTON_SPHERE_SCHWARZSCHILD = 3.0;  // r = 3M

/// @brief Schwarzschild ISCO radius
constexpr double ISCO_SCHWARZSCHILD = 6.0;  // r = 6M

/// @brief Extremal Kerr prograde ISCO radius
constexpr double ISCO_KERR_EXTREMAL = 1.0;  // r = M for a = M

/// @brief Einstein light deflection angle (weak field, grazing)
/// Δφ = 4M/b ≈ 1.75 arcsec for Sun
constexpr double EINSTEIN_DEFLECTION_COEFF = 4.0;  // Δφ = 4M/b

/// @brief Schwarzschild horizon radius
constexpr double HORIZON_SCHWARZSCHILD = 2.0;  // r = 2M

} // namespace PhysicalRef

// =============================================================================
// Test Integration Parameters
// =============================================================================

namespace IntegrationParams {

/// @brief Default number of integration steps for conservation tests
constexpr int CONSERVATION_STEPS = 1000;

/// @brief Minimum steps for meaningful statistical tests
constexpr int MINIMUM_STEPS = 100;

/// @brief Steps for long-term stability tests
constexpr int LONG_TERM_STEPS = 10000;

/// @brief Default initial step size
constexpr float DEFAULT_STEP = 0.01f;

/// @brief Tight tolerance step size for precision tests
constexpr float TIGHT_STEP = 0.001f;

/// @brief Maximum radial coordinate before declaring escape
constexpr double ESCAPE_RADIUS = 1000.0;

/// @brief Minimum radial coordinate buffer above horizon
constexpr double HORIZON_BUFFER = 1.001;

} // namespace IntegrationParams

// =============================================================================
// Mathematical Constants
// =============================================================================

namespace MathConst {

using namespace Sirius::Constants::Math;

/// @brief Pi
constexpr double PI = Sirius::Constants::Math::PI;

/// @brief 2π
constexpr double TWO_PI = Sirius::Constants::Math::TWO_PI;

/// @brief π/2
constexpr double HALF_PI = Sirius::Constants::Math::HALF_PI;

/// @brief Standard test angles for parametric tests
constexpr double TEST_ANGLES[] = {0.0, PI/6, PI/4, PI/3, HALF_PI, 2*PI/3, 3*PI/4, PI};
constexpr int NUM_TEST_ANGLES = 8;

/// @brief Standard test radii (in units of M)
constexpr double TEST_RADII[] = {3.0, 5.0, 6.0, 10.0, 20.0, 50.0, 100.0};
constexpr int NUM_TEST_RADII = 7;

/// @brief Standard spin parameters for Kerr tests
constexpr double TEST_SPINS[] = {0.0, 0.3, 0.5, 0.7, 0.9, 0.99};
constexpr int NUM_TEST_SPINS = 6;

} // namespace MathConst

// =============================================================================
// Utility Functions
// =============================================================================

/// @brief Select appropriate null condition tolerance based on type
template<typename T>
constexpr T nullConditionTolerance() {
    return Sirius::Constants::nullConditionTolerance<T>();
}

/// @brief Check if value is within tolerance of expected
template<typename T>
inline bool isNear(T actual, T expected, T tolerance) {
    return std::abs(actual - expected) < tolerance;
}

/// @brief Check relative error is within tolerance
template<typename T>
inline bool isRelativelyNear(T actual, T expected, T tolerance) {
    if (std::abs(expected) < static_cast<T>(1e-30)) {
        return std::abs(actual) < tolerance;
    }
    return std::abs((actual - expected) / expected) < tolerance;
}

/// @brief Compute relative error
template<typename T>
inline T relativeError(T actual, T expected) {
    if (std::abs(expected) < static_cast<T>(1e-30)) {
        return std::abs(actual);
    }
    return std::abs((actual - expected) / expected);
}

} // namespace sirius::test

#endif // TSFT001A_H
