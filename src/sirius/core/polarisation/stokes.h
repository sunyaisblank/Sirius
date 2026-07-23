#pragma once

// Stokes formalism (I, Q, U, V) for partially polarised light: polarisation
// degree and EVPA, Mueller matrices for standard optical elements, synchrotron
// and Thomson emission models, and parallel-transport Faraday rotation.
// Ported from PHPL001A.h.
// Reference: Chandrasekhar (1960) "Radiative Transfer"; Rybicki & Lightman
// (1979) Chapter 2.

#include <algorithm>
#include <array>
#include <cmath>

namespace sirius::core {

// Stokes vector S = (I, Q, U, V) fully characterising a polarisation state.
//   I = |E_x|^2 + |E_y|^2                  total intensity
//   Q = |E_x|^2 - |E_y|^2                  linear (horizontal - vertical)
//   U = 2 Re(E_x* E_y)                     linear (+45 deg - -45 deg)
//   V = 2 Im(E_x* E_y)                     circular (RCP - LCP)
// Physical vectors satisfy Q^2 + U^2 + V^2 <= I^2.
struct StokesVector {
    float I = 0.0f;  // Total intensity.
    float Q = 0.0f;  // Linear polarisation (horizontal - vertical).
    float U = 0.0f;  // Linear polarisation (+45 deg - -45 deg).
    float V = 0.0f;  // Circular polarisation (RCP - LCP).

    StokesVector() = default;
    StokesVector(float i, float q, float u, float v) : I(i), Q(q), U(u), V(v) {}

    // Unpolarised light of the given intensity.
    static StokesVector Unpolarised(float intensity) {
        return StokesVector(intensity, 0.0f, 0.0f, 0.0f);
    }

    // Horizontally polarised light.
    static StokesVector Horizontal(float intensity) {
        return StokesVector(intensity, intensity, 0.0f, 0.0f);
    }

    // Vertically polarised light.
    static StokesVector Vertical(float intensity) {
        return StokesVector(intensity, -intensity, 0.0f, 0.0f);
    }

    // +45 deg linearly polarised light.
    static StokesVector Diagonal45(float intensity) {
        return StokesVector(intensity, 0.0f, intensity, 0.0f);
    }

    // Right circularly polarised light.
    static StokesVector RightCircular(float intensity) {
        return StokesVector(intensity, 0.0f, 0.0f, intensity);
    }

    // Left circularly polarised light.
    static StokesVector LeftCircular(float intensity) {
        return StokesVector(intensity, 0.0f, 0.0f, -intensity);
    }

    // Total degree of polarisation p = sqrt(Q^2 + U^2 + V^2) / I, in [0, 1].
    float PolarisationDegree() const {
        if (I <= 0.0f) return 0.0f;
        float p_sq = (Q * Q + U * U + V * V) / (I * I);
        return std::sqrt(std::min(p_sq, 1.0f));  // Clamp for numerical safety.
    }

    // Linear polarisation degree p_L = sqrt(Q^2 + U^2) / I.
    float LinearPolarisationDegree() const {
        if (I <= 0.0f) return 0.0f;
        float p_L_sq = (Q * Q + U * U) / (I * I);
        return std::sqrt(std::min(p_L_sq, 1.0f));
    }

    // Circular polarisation degree (signed) p_C = V / I.
    float CircularPolarisationDegree() const {
        if (I <= 0.0f) return 0.0f;
        return V / I;
    }

    // Electric Vector Position Angle chi = 0.5 * arctan(U/Q), in [-pi/2, pi/2].
    float Evpa() const { return 0.5f * std::atan2(U, Q); }

    // True when Q^2 + U^2 + V^2 <= I^2 (small tolerance).
    bool IsPhysical() const {
        if (I < 0.0f) return false;
        return (Q * Q + U * U + V * V) <= I * I * 1.001f;  // Small tolerance.
    }

    // Project onto the Poincare sphere to restore physicality.
    void Normalise() {
        if (I <= 0.0f) {
            I = Q = U = V = 0.0f;
            return;
        }
        float p_sq = Q * Q + U * U + V * V;
        if (p_sq > I * I) {
            float scale = I / std::sqrt(p_sq);
            Q *= scale;
            U *= scale;
            V *= scale;
        }
    }

    StokesVector operator+(const StokesVector& other) const {
        return StokesVector(I + other.I, Q + other.Q, U + other.U, V + other.V);
    }

    StokesVector& operator+=(const StokesVector& other) {
        I += other.I;
        Q += other.Q;
        U += other.U;
        V += other.V;
        return *this;
    }

    StokesVector operator*(float scalar) const {
        return StokesVector(I * scalar, Q * scalar, U * scalar, V * scalar);
    }

    StokesVector& operator*=(float scalar) {
        I *= scalar;
        Q *= scalar;
        U *= scalar;
        V *= scalar;
        return *this;
    }
};

// 4x4 Mueller matrix transforming Stokes vectors: S_out = M S_in.
// Used for optical elements (polarisers, retarders, rotators).
struct MuellerMatrix {
    std::array<std::array<float, 4>, 4> m;

    MuellerMatrix() {
        // Default to identity.
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++) m[i][j] = (i == j) ? 1.0f : 0.0f;
    }

    // Apply this Mueller matrix to a Stokes vector.
    StokesVector Apply(const StokesVector& s) const {
        StokesVector out;
        out.I = m[0][0] * s.I + m[0][1] * s.Q + m[0][2] * s.U + m[0][3] * s.V;
        out.Q = m[1][0] * s.I + m[1][1] * s.Q + m[1][2] * s.U + m[1][3] * s.V;
        out.U = m[2][0] * s.I + m[2][1] * s.Q + m[2][2] * s.U + m[2][3] * s.V;
        out.V = m[3][0] * s.I + m[3][1] * s.Q + m[3][2] * s.U + m[3][3] * s.V;
        return out;
    }

    // Matrix product M_result = M1 M2.
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

    // Horizontal linear polariser.
    static MuellerMatrix HorizontalPolariser() {
        MuellerMatrix M;
        M.m = {{{0.5f, 0.5f, 0.0f, 0.0f},
                {0.5f, 0.5f, 0.0f, 0.0f},
                {0.0f, 0.0f, 0.0f, 0.0f},
                {0.0f, 0.0f, 0.0f, 0.0f}}};
        return M;
    }

    // Vertical linear polariser.
    static MuellerMatrix VerticalPolariser() {
        MuellerMatrix M;
        M.m = {{{0.5f, -0.5f, 0.0f, 0.0f},
                {-0.5f, 0.5f, 0.0f, 0.0f},
                {0.0f, 0.0f, 0.0f, 0.0f},
                {0.0f, 0.0f, 0.0f, 0.0f}}};
        return M;
    }

    // Linear polariser at angle theta.
    static MuellerMatrix LinearPolariser(float theta) {
        float c2 = std::cos(2.0f * theta);
        float s2 = std::sin(2.0f * theta);
        MuellerMatrix M;
        M.m = {{{0.5f, 0.5f * c2, 0.5f * s2, 0.0f},
                {0.5f * c2, 0.5f * c2 * c2, 0.5f * s2 * c2, 0.0f},
                {0.5f * s2, 0.5f * s2 * c2, 0.5f * s2 * s2, 0.0f},
                {0.0f, 0.0f, 0.0f, 0.0f}}};
        return M;
    }

    // Quarter-wave plate (fast axis horizontal).
    static MuellerMatrix QuarterWavePlate() {
        MuellerMatrix M;
        M.m = {{{1.0f, 0.0f, 0.0f, 0.0f},
                {0.0f, 1.0f, 0.0f, 0.0f},
                {0.0f, 0.0f, 0.0f, 1.0f},
                {0.0f, 0.0f, -1.0f, 0.0f}}};
        return M;
    }

    // Half-wave plate (fast axis horizontal).
    static MuellerMatrix HalfWavePlate() {
        MuellerMatrix M;
        M.m = {{{1.0f, 0.0f, 0.0f, 0.0f},
                {0.0f, 1.0f, 0.0f, 0.0f},
                {0.0f, 0.0f, -1.0f, 0.0f},
                {0.0f, 0.0f, 0.0f, -1.0f}}};
        return M;
    }

    // Rotation of the reference frame by angle theta.
    static MuellerMatrix Rotation(float theta) {
        float c2 = std::cos(2.0f * theta);
        float s2 = std::sin(2.0f * theta);
        MuellerMatrix M;
        M.m = {{{1.0f, 0.0f, 0.0f, 0.0f},
                {0.0f, c2, s2, 0.0f},
                {0.0f, -s2, c2, 0.0f},
                {0.0f, 0.0f, 0.0f, 1.0f}}};
        return M;
    }

    // Faraday rotation of the polarisation plane.
    // rotation_angle = RM * lambda^2 (rotation measure times wavelength squared).
    static MuellerMatrix FaradayRotation(float rotation_angle) { return Rotation(rotation_angle); }

    // Partial depolariser retaining fraction p of the polarisation
    // (0 = fully depolarised, 1 = unchanged).
    static MuellerMatrix Depolariser(float p) {
        MuellerMatrix M;
        M.m = {{{1.0f, 0.0f, 0.0f, 0.0f},
                {0.0f, p, 0.0f, 0.0f},
                {0.0f, 0.0f, p, 0.0f},
                {0.0f, 0.0f, 0.0f, p}}};
        return M;
    }
};

// Polarised emission models.
namespace polarised_emission {

// Synchrotron polarisation degree for a power-law electron distribution
// n(E) ~ E^(-p): pi_L = (p + 1) / (p + 7/3). Typical p = 2-3 gives ~0.69-0.75.
inline float SynchrotronPolarisationDegree(float spectral_index) {
    float p = spectral_index;
    return (p + 1.0f) / (p + 7.0f / 3.0f);
}

// Synchrotron Stokes vector for emission from an ordered magnetic field;
// evpa is the electric vector position angle (perpendicular to B).
inline StokesVector SynchrotronEmission(float intensity, float pol_degree, float evpa) {
    StokesVector s;
    s.I = intensity;
    s.Q = intensity * pol_degree * std::cos(2.0f * evpa);
    s.U = intensity * pol_degree * std::sin(2.0f * evpa);
    s.V = 0.0f;  // No circular polarisation from synchrotron (in vacuum).
    return s;
}

// Thomson single-scattering polarisation degree as a function of scattering
// angle theta: pi_L = sin^2(theta) / (1 + cos^2(theta)).
inline float ThomsonPolarisationDegree(float cos_theta) {
    float sin2 = 1.0f - cos_theta * cos_theta;
    float cos2 = cos_theta * cos_theta;
    if (1.0f + cos2 < 1e-10f) return 0.0f;
    return sin2 / (1.0f + cos2);
}

// Thermal emission polarisation from a magnetised atmosphere (simplified:
// scales with B-field and sin^2 of the angle to B). An accurate model requires
// solving the Stokes-I radiative transfer.
inline float ThermalMagneticPolarisation(float temperature, float B_field_strength,
                                         float angle_to_B) {
    constexpr float kTRef = 10000.0f;  // Reference temperature.
    constexpr float kBRef = 1000.0f;   // Reference B-field (Gauss).

    float B_factor = std::min(B_field_strength / kBRef, 1.0f);
    float T_factor = std::min(kTRef / temperature, 1.0f);
    float sin2_angle = 1.0f - std::cos(angle_to_B) * std::cos(angle_to_B);

    return B_factor * T_factor * sin2_angle * 0.3f;  // Max ~30% for strong B.
}

}  // namespace polarised_emission

// Parallel transport of the polarisation basis along a geodesic. In curved
// spacetime this produces gravitational Faraday rotation (frame dragging in
// Kerr) and geometric rotation from the non-flat connection.
// Reference: Connors, Piran & Stark (1980), ApJ 235, 224.
namespace parallel_transport {

// Gravitational Faraday rotation angle. Simplified small-a formula
// delta_chi ~ 2a / r_impact^2 * sin(inclination); the full result integrates
// along the geodesic.
inline float GravitationalFaradayRotation(float a, float r_impact, float inclination) {
    float factor = 2.0f * a / (r_impact * r_impact);
    return factor * std::sin(inclination);
}

// Apply an accumulated parallel-transport rotation to a Stokes vector.
inline StokesVector ApplyParallelTransport(const StokesVector& s, float rotation_angle) {
    auto M = MuellerMatrix::Rotation(rotation_angle);
    return M.Apply(s);
}

}  // namespace parallel_transport

}  // namespace sirius::core
