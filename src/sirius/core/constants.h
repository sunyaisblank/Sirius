#pragma once

// Single authority for tolerances, physical constants, and precision
// policy (ported from PHCN001A.h; values bit-identical, names per
// docs/STYLE.md). Tests and validators cite these directly; a numeric
// literal standing in for one of these at a use site is a defect.
//
// Precision tiers: kernels accumulate in single precision (epsilon ~1.2e-7),
// the CPU path and oracle in double (epsilon ~2.2e-16); tolerances that
// differ by tier say so in their names.

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <type_traits>

namespace sirius::core::constants {

// --- Machine precision -------------------------------------------------------

inline constexpr float kEpsilonF = std::numeric_limits<float>::epsilon();
inline constexpr double kEpsilonD = std::numeric_limits<double>::epsilon();

// Safe minima for division floors.
inline constexpr float kSafeMinF = 1e-10f;
inline constexpr double kSafeMinD = 1e-20;

// --- Geodesic integration ----------------------------------------------------

namespace geodesic {

// Null condition g_uv k^u k^v = 0; the kernel tier is relaxed because
// single-precision accumulation absorbs ~3 decimal digits over a ray.
inline constexpr float kNullConditionTolKernel = 1e-5f;
inline constexpr double kNullConditionTolCpu = 1e-10;

// Conservation drift (E, L_z), relative, over a full integration; the
// Mandatory suite fails the build beyond this.
inline constexpr double kConservationTol = 1e-4;

// H = (1/2) g^uv p_u p_v for null geodesics.
inline constexpr double kHamiltonianTol = 1e-10;

inline constexpr int kMaxIntegrationSteps = 10000;
inline constexpr float kDefaultStepSize = 0.01f;
inline constexpr float kMinStepSize = 1e-8f;
inline constexpr float kMaxStepSize = 1.0f;
inline constexpr float kStepSafetyFactor = 0.9f;

}  // namespace geodesic

// --- Coordinate domains ------------------------------------------------------

namespace coordinates {

// Minimum r/r_+ ratio; inside this buffer a ray is a horizon crossing.
inline constexpr double kHorizonBuffer = 1.001;

// Minimum |sin(theta)| before clamping away from the axis singularity.
inline constexpr double kPoleEpsilon = 1e-6;

// Radius (units of M) beyond which a ray has escaped.
inline constexpr double kEscapeRadius = 1e6;

inline constexpr double kCoordinateTol = 1e-12;
inline constexpr double kPhiWrapTol = 1e-10;

}  // namespace coordinates

// --- Ray termination (Kerr-Schild Cartesian, units of M) ---------------------

namespace termination {

// Beyond this and moving outward the geodesic is asymptotically straight.
inline constexpr double kEscapeRadius = 50.0;

// Beyond this the ray has escaped regardless of direction.
inline constexpr double kBackgroundRadius = 100.0;

// Terminate at r <= r_horizon * (1 + kCaptureMargin).
inline constexpr double kCaptureMargin = 0.05;

// Affine-parameter budget per ray; prevents runaway integration.
inline constexpr double kMaxAffineParameter = 100.0;

// Velocity magnitude below which a ray is considered stuck.
inline constexpr double kStalledVelocity = 1e-8;

}  // namespace termination

// --- Metric tensor invariants ------------------------------------------------

namespace metric {

inline constexpr double kSymmetryTol = 1e-15;
inline constexpr double kInverseTol = 1e-14;
inline constexpr double kChristoffelSymmetryTol = 1e-15;
inline constexpr double kDeterminantTol = 1e-30;
inline constexpr double kSignatureTol = 1e-10;

}  // namespace metric

// --- Numerical differentiation (only where analytic forms are absent) --------

namespace differentiation {

// Central-difference step h = kRelativeH * |x| + kAbsoluteH.
inline constexpr double kRelativeH = 1e-6;
inline constexpr double kAbsoluteH = 1e-10;
inline constexpr double kRichardsonTol = 1e-10;

}  // namespace differentiation

// --- Spectral and radiative transfer -----------------------------------------

namespace spectral {

// Visible range [m].
inline constexpr double kLambdaMin = 380e-9;
inline constexpr double kLambdaMax = 780e-9;

inline constexpr int kSpectralBins = 32;
inline constexpr double kIntensityFloor = 1e-30;

// Redshift factor clamp.
inline constexpr double kMaxGFactor = 100.0;
inline constexpr double kMinGFactor = 0.01;

inline constexpr double kLimbDarkMin = 0.0;
inline constexpr double kLimbDarkMax = 1.0;

}  // namespace spectral

// --- Disk models --------------------------------------------------------------

namespace disk {

// |theta - pi/2| below this defines the thin-disk plane (~0.5 degrees).
inline constexpr double kEquatorialTol = 0.01;

// Inner emission edge relative to the ISCO.
inline constexpr double kInnerEdgeBuffer = 1.001;

// Operator disk-temperature scale is the effective temperature at this
// multiple of the inner edge for every selectable radial profile. A
// zero-torque profile is identically cold at the edge itself, so calling the
// scale an "inner temperature" would misstate its physical meaning.
inline constexpr double kTemperatureReferenceRadiusRatio = 1.5;

// Simpson points for the Page-Thorne flux integral.
inline constexpr int kIntegrationPoints = 100;

// Q-factor validity clamp.
inline constexpr double kQFactorMin = 0.01;
inline constexpr double kQFactorMax = 100.0;

}  // namespace disk

// --- Numerical safety ----------------------------------------------------------

namespace tolerances {

inline constexpr double kMetricInversionEps = 1e-12;
inline constexpr double kAngularClampEps = 1e-6;
inline constexpr double kDivisionSafeEps = 1e-15;
inline constexpr double kFluxIntegralEps = 1e-30;

}  // namespace tolerances

// --- Physical constants (SI; CODATA 2018) ------------------------------------

namespace physical {

inline constexpr double kPlanck = 6.62607015e-34;           // [J s]
inline constexpr double kSpeedOfLight = 2.99792458e8;       // [m/s]
inline constexpr double kBoltzmann = 1.380649e-23;          // [J/K]
inline constexpr double kStefanBoltzmann = 5.670374419e-8;  // [W/(m^2 K^4)]
inline constexpr double kGravitation = 6.67430e-11;         // [m^3/(kg s^2)]
inline constexpr double kSolarMass = 1.98892e30;            // [kg]
inline constexpr double kParsec = 3.08567758e16;            // [m]
inline constexpr double kWienB = 2.897771955e-3;            // [m K]

inline constexpr double kPlanckC1 = 2.0 * kPlanck * kSpeedOfLight * kSpeedOfLight;
inline constexpr double kPlanckC2 = kPlanck * kSpeedOfLight / kBoltzmann;
inline constexpr double kGmSun = kGravitation * kSolarMass;
inline constexpr double kSchwarzschildRadiusSun = 2.0 * kGmSun / (kSpeedOfLight * kSpeedOfLight);

}  // namespace physical

// --- Mathematical constants ----------------------------------------------------

namespace math {

inline constexpr double kPi = std::numbers::pi;
inline constexpr double kTwoPi = 2.0 * kPi;
inline constexpr double kHalfPi = kPi / 2.0;
inline constexpr double kSqrt2 = std::numbers::sqrt2;
inline constexpr double kSqrt3 = std::numbers::sqrt3;
inline constexpr double kInvSqrt2 = 0.70710678118654752440;

// Cube root of 2, used by the Yoshida composition coefficients.
inline constexpr double kCbrt2 = 1.25992104989487316477;

}  // namespace math

// --- Utility functions -----------------------------------------------------------

// Tolerance for the null condition selected by accumulation precision.
template <typename T>
[[nodiscard]] constexpr T NullConditionTolerance() {
    if constexpr (std::is_same_v<T, float>) {
        return geodesic::kNullConditionTolKernel;
    } else {
        return static_cast<T>(geodesic::kNullConditionTolCpu);
    }
}

// Division with a floored denominator; the floor preserves sign.
template <typename T>
[[nodiscard]] inline T SafeDivide(T numerator, T denominator, T floor = kSafeMinD) {
    if (std::abs(denominator) < floor) {
        denominator = std::copysign(floor, denominator);
    }
    return numerator / denominator;
}

// Clamps theta away from the axis singularity.
[[nodiscard]] inline double ClampTheta(double theta, double epsilon = coordinates::kPoleEpsilon) {
    return std::clamp(theta, epsilon, math::kPi - epsilon);
}

// Wraps phi to [0, 2 pi).
[[nodiscard]] inline double WrapPhi(double phi) {
    phi = std::fmod(phi, math::kTwoPi);
    if (phi < 0) {
        phi += math::kTwoPi;
    }
    return phi;
}

}  // namespace sirius::core::constants
