// PHPL001A.h - Stokes Parameters for Polarised Light Transport
// Component ID: PHPL001A (Physics/Polarisation/StokesVector)
//
// Implements the Stokes formalism for describing partially polarised light.
// Stokes vector S = (I, Q, U, V) fully characterizes polarisation state.
//
// MATHEMATICAL FOUNDATION:
// ========================
// I = total intensity = |E_x|² + |E_y|²
// Q = linear polarisation (horizontal - vertical) = |E_x|² - |E_y|²
// U = linear polarisation (+45° - -45°) = 2 Re(E_x* E_y)
// V = circular polarisation (RCP - LCP) = 2 Im(E_x* E_y)
//
// Degree of polarisation: p = sqrt(Q² + U² + V²) / I
// Linear polarisation: p_L = sqrt(Q² + U²) / I
// Circular polarisation: p_C = |V| / I
// EVPA (polarisation angle): χ = 0.5 × arctan(U/Q)
//
// PHYSICAL CONSTRAINTS:
// =====================
// 1. I ≥ 0 (non-negative intensity)
// 2. Q² + U² + V² ≤ I² (physical Stokes vectors)
// 3. For fully polarised light: Q² + U² + V² = I²
//
// Reference: Chandrasekhar (1960) "Radiative Transfer"
//            Rybicki & Lightman (1979) Chapter 2

#pragma once

#include <cmath>
#include <algorithm>
#include <array>

namespace Sirius {

//==============================================================================
// Stokes Vector
//==============================================================================
struct StokesVector {
    float I = 0.0f;  ///< Total intensity
    float Q = 0.0f;  ///< Linear polarisation (horizontal - vertical)
    float U = 0.0f;  ///< Linear polarisation (+45° - -45°)
    float V = 0.0f;  ///< Circular polarisation (RCP - LCP)

    //--------------------------------------------------------------------------
    // Constructors
    //--------------------------------------------------------------------------
    StokesVector() = default;
    StokesVector(float i, float q, float u, float v) : I(i), Q(q), U(u), V(v) {}

    /// @brief Create unpolarised light of given intensity
    static StokesVector unpolarised(float intensity) {
        return StokesVector(intensity, 0.0f, 0.0f, 0.0f);
    }

    /// @brief Create horizontally polarised light
    static StokesVector horizontal(float intensity) {
        return StokesVector(intensity, intensity, 0.0f, 0.0f);
    }

    /// @brief Create vertically polarised light
    static StokesVector vertical(float intensity) {
        return StokesVector(intensity, -intensity, 0.0f, 0.0f);
    }

    /// @brief Create +45° linearly polarised light
    static StokesVector diagonal45(float intensity) {
        return StokesVector(intensity, 0.0f, intensity, 0.0f);
    }

    /// @brief Create right circularly polarised light
    static StokesVector rightCircular(float intensity) {
        return StokesVector(intensity, 0.0f, 0.0f, intensity);
    }

    /// @brief Create left circularly polarised light
    static StokesVector leftCircular(float intensity) {
        return StokesVector(intensity, 0.0f, 0.0f, -intensity);
    }

    //--------------------------------------------------------------------------
    // Polarisation Properties
    //--------------------------------------------------------------------------

    /// @brief Total degree of polarisation p ∈ [0, 1]
    float polarisationDegree() const {
        if (I <= 0.0f) return 0.0f;
        float p_sq = (Q*Q + U*U + V*V) / (I*I);
        return std::sqrt(std::min(p_sq, 1.0f));  // Clamp for numerical safety
    }

    /// @brief Linear polarisation degree
    float linearPolarisationDegree() const {
        if (I <= 0.0f) return 0.0f;
        float p_L_sq = (Q*Q + U*U) / (I*I);
        return std::sqrt(std::min(p_L_sq, 1.0f));
    }

    /// @brief Circular polarisation degree (signed)
    float circularPolarisationDegree() const {
        if (I <= 0.0f) return 0.0f;
        return V / I;
    }

    /// @brief Electric Vector Position Angle (EVPA) in radians
    /// @return EVPA χ ∈ [-π/2, π/2]
    float EVPA() const {
        return 0.5f * std::atan2(U, Q);
    }

    /// @brief Check if this is a physical Stokes vector
    bool isPhysical() const {
        if (I < 0.0f) return false;
        return (Q*Q + U*U + V*V) <= I*I * 1.001f;  // Small tolerance
    }

    /// @brief Normalise to ensure physicality (project onto Poincaré sphere)
    void normalise() {
        if (I <= 0.0f) {
            I = Q = U = V = 0.0f;
            return;
        }
        float p_sq = Q*Q + U*U + V*V;
        if (p_sq > I*I) {
            float scale = I / std::sqrt(p_sq);
            Q *= scale;
            U *= scale;
            V *= scale;
        }
    }

    //--------------------------------------------------------------------------
    // Arithmetic Operators
    //--------------------------------------------------------------------------
    StokesVector operator+(const StokesVector& other) const {
        return StokesVector(I + other.I, Q + other.Q, U + other.U, V + other.V);
    }

    StokesVector& operator+=(const StokesVector& other) {
        I += other.I; Q += other.Q; U += other.U; V += other.V;
        return *this;
    }

    StokesVector operator*(float scalar) const {
        return StokesVector(I * scalar, Q * scalar, U * scalar, V * scalar);
    }

    StokesVector& operator*=(float scalar) {
        I *= scalar; Q *= scalar; U *= scalar; V *= scalar;
        return *this;
    }
};

//==============================================================================
// Mueller Matrix (4x4 transformation of Stokes vectors)
//==============================================================================
// S_out = M × S_in
// Used for optical elements: polarisers, retarders, rotators, etc.
//==============================================================================
struct MuellerMatrix {
    std::array<std::array<float, 4>, 4> m;

    MuellerMatrix() {
        // Default to identity
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                m[i][j] = (i == j) ? 1.0f : 0.0f;
    }

    /// @brief Apply Mueller matrix to Stokes vector
    StokesVector apply(const StokesVector& s) const {
        StokesVector out;
        out.I = m[0][0]*s.I + m[0][1]*s.Q + m[0][2]*s.U + m[0][3]*s.V;
        out.Q = m[1][0]*s.I + m[1][1]*s.Q + m[1][2]*s.U + m[1][3]*s.V;
        out.U = m[2][0]*s.I + m[2][1]*s.Q + m[2][2]*s.U + m[2][3]*s.V;
        out.V = m[3][0]*s.I + m[3][1]*s.Q + m[3][2]*s.U + m[3][3]*s.V;
        return out;
    }

    /// @brief Matrix multiplication: M_result = M1 × M2
    MuellerMatrix operator*(const MuellerMatrix& other) const {
        MuellerMatrix result;
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                result.m[i][j] = 0.0f;
                for (int k = 0; k < 4; k++) {
                    result.m[i][j] += m[i][k] * other.m[k][j];
                }
            }
        }
        return result;
    }

    //--------------------------------------------------------------------------
    // Standard Optical Elements
    //--------------------------------------------------------------------------

    /// @brief Horizontal linear polariser
    static MuellerMatrix horizontalPolariser() {
        MuellerMatrix M;
        M.m = {{
            {0.5f, 0.5f, 0.0f, 0.0f},
            {0.5f, 0.5f, 0.0f, 0.0f},
            {0.0f, 0.0f, 0.0f, 0.0f},
            {0.0f, 0.0f, 0.0f, 0.0f}
        }};
        return M;
    }

    /// @brief Vertical linear polariser
    static MuellerMatrix verticalPolariser() {
        MuellerMatrix M;
        M.m = {{
            {0.5f, -0.5f, 0.0f, 0.0f},
            {-0.5f, 0.5f, 0.0f, 0.0f},
            {0.0f, 0.0f, 0.0f, 0.0f},
            {0.0f, 0.0f, 0.0f, 0.0f}
        }};
        return M;
    }

    /// @brief Linear polariser at angle θ
    static MuellerMatrix linearPolariser(float theta) {
        float c2 = std::cos(2.0f * theta);
        float s2 = std::sin(2.0f * theta);
        MuellerMatrix M;
        M.m = {{
            {0.5f, 0.5f * c2, 0.5f * s2, 0.0f},
            {0.5f * c2, 0.5f * c2 * c2, 0.5f * s2 * c2, 0.0f},
            {0.5f * s2, 0.5f * s2 * c2, 0.5f * s2 * s2, 0.0f},
            {0.0f, 0.0f, 0.0f, 0.0f}
        }};
        return M;
    }

    /// @brief Quarter-wave plate (fast axis horizontal)
    static MuellerMatrix quarterWavePlate() {
        MuellerMatrix M;
        M.m = {{
            {1.0f, 0.0f, 0.0f, 0.0f},
            {0.0f, 1.0f, 0.0f, 0.0f},
            {0.0f, 0.0f, 0.0f, 1.0f},
            {0.0f, 0.0f, -1.0f, 0.0f}
        }};
        return M;
    }

    /// @brief Half-wave plate (fast axis horizontal)
    static MuellerMatrix halfWavePlate() {
        MuellerMatrix M;
        M.m = {{
            {1.0f, 0.0f, 0.0f, 0.0f},
            {0.0f, 1.0f, 0.0f, 0.0f},
            {0.0f, 0.0f, -1.0f, 0.0f},
            {0.0f, 0.0f, 0.0f, -1.0f}
        }};
        return M;
    }

    /// @brief Rotation matrix (rotates reference frame by angle θ)
    static MuellerMatrix rotation(float theta) {
        float c2 = std::cos(2.0f * theta);
        float s2 = std::sin(2.0f * theta);
        MuellerMatrix M;
        M.m = {{
            {1.0f, 0.0f, 0.0f, 0.0f},
            {0.0f, c2, s2, 0.0f},
            {0.0f, -s2, c2, 0.0f},
            {0.0f, 0.0f, 0.0f, 1.0f}
        }};
        return M;
    }

    /// @brief Faraday rotation matrix (rotation of polarisation plane)
    /// @param rotation_angle RM × λ² (rotation measure times wavelength squared)
    static MuellerMatrix faradayRotation(float rotation_angle) {
        return rotation(rotation_angle);
    }

    /// @brief Depolariser (partial depolarisation with factor p)
    /// @param p Fraction of polarisation retained (0 = fully depolarised, 1 = unchanged)
    static MuellerMatrix depolariser(float p) {
        MuellerMatrix M;
        M.m = {{
            {1.0f, 0.0f, 0.0f, 0.0f},
            {0.0f, p, 0.0f, 0.0f},
            {0.0f, 0.0f, p, 0.0f},
            {0.0f, 0.0f, 0.0f, p}
        }};
        return M;
    }
};

//==============================================================================
// Polarised Emission Models
//==============================================================================
namespace PolarisedEmission {

/// @brief Synchrotron emission polarisation degree
/// For a power-law electron distribution n(E) ∝ E^(-p):
///   π_L = (p + 1) / (p + 7/3)
/// Typical values: p = 2-3 gives π_L ≈ 0.69-0.75
/// @param spectral_index Electron spectral index p
/// @return Maximum linear polarisation degree
inline float synchrotronPolarisationDegree(float spectral_index) {
    float p = spectral_index;
    return (p + 1.0f) / (p + 7.0f/3.0f);
}

/// @brief Synchrotron Stokes vector for emission from ordered magnetic field
/// @param intensity Total intensity
/// @param pol_degree Polarisation degree (from synchrotronPolarisationDegree)
/// @param evpa Electric vector position angle (perpendicular to B-field)
/// @return Stokes vector for synchrotron emission
inline StokesVector synchrotronEmission(float intensity, float pol_degree, float evpa) {
    StokesVector s;
    s.I = intensity;
    s.Q = intensity * pol_degree * std::cos(2.0f * evpa);
    s.U = intensity * pol_degree * std::sin(2.0f * evpa);
    s.V = 0.0f;  // No circular polarisation from synchrotron (in vacuum)
    return s;
}

/// @brief Thomson scattering polarisation
/// Single scattering from electrons produces polarisation depending on
/// scattering angle θ:
///   π_L = sin²θ / (1 + cos²θ)
/// @param cos_theta Cosine of scattering angle
/// @return Linear polarisation degree from Thomson scattering
inline float thomsonPolarisationDegree(float cos_theta) {
    float sin2 = 1.0f - cos_theta * cos_theta;
    float cos2 = cos_theta * cos_theta;
    if (1.0f + cos2 < 1e-10f) return 0.0f;
    return sin2 / (1.0f + cos2);
}

/// @brief Thermal emission polarisation (from magnetic atmosphere)
/// Thermal emission from a magnetised atmosphere can be partially polarised
/// @param temperature Temperature in K
/// @param B_field_strength Magnetic field strength in Gauss
/// @param angle_to_B Angle between line of sight and magnetic field
/// @return Estimated polarisation degree (simplified model)
inline float thermalMagneticPolarisation(float temperature, float B_field_strength, float angle_to_B) {
    // Simplified model: polarisation scales with B-field and sin²(angle)
    // More accurate model requires solving Stokes-I radiative transfer
    constexpr float T_REF = 10000.0f;  // Reference temperature
    constexpr float B_REF = 1000.0f;   // Reference B-field (Gauss)

    float B_factor = std::min(B_field_strength / B_REF, 1.0f);
    float T_factor = std::min(T_REF / temperature, 1.0f);
    float sin2_angle = 1.0f - std::cos(angle_to_B) * std::cos(angle_to_B);

    return B_factor * T_factor * sin2_angle * 0.3f;  // Max ~30% for strong B
}

} // namespace PolarisedEmission

//==============================================================================
// Parallel Transport of Polarisation
//==============================================================================
// In curved spacetime, the polarisation basis must be parallel transported
// along the geodesic. This causes:
// 1. Gravitational Faraday rotation (frame dragging in Kerr spacetime)
// 2. Geometric rotation from non-flat connection
//
// Reference: Connors, Piran & Stark (1980), ApJ 235, 224
//==============================================================================
namespace ParallelTransport {

/// @brief Gravitational Faraday rotation angle
/// In Kerr spacetime, photon polarisation rotates due to frame dragging
/// @param a Kerr spin parameter (a/M)
/// @param r_impact Impact radius (closest approach)
/// @param inclination Observer inclination angle
/// @return Rotation angle in radians
inline float gravitationalFaradayRotation(float a, float r_impact, float inclination) {
    // Simplified formula for small a
    // Full calculation requires integrating along the geodesic
    // Δχ ≈ 2a / r_impact² × sin(inclination)
    float factor = 2.0f * a / (r_impact * r_impact);
    return factor * std::sin(inclination);
}

/// @brief Apply parallel transport rotation to Stokes vector
/// @param s Input Stokes vector
/// @param rotation_angle Accumulated rotation from parallel transport
/// @return Rotated Stokes vector
inline StokesVector applyParallelTransport(const StokesVector& s, float rotation_angle) {
    auto M = MuellerMatrix::rotation(rotation_angle);
    return M.apply(s);
}

} // namespace ParallelTransport

} // namespace Sirius
