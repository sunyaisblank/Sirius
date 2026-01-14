// RDMetricKerr.cuh - Kerr Metric Utilities
// Component ID: RDTR002A (Transport/MetricKerr)

#pragma once

#include "../Integration/RDMath.cuh"
#include "RDGeodesic.cuh"

namespace Sirius {

//==============================================================================
// Kerr Metric Helper Functions for Hamiltonian Formulation
//==============================================================================

// Compute Kerr metric quantities at (r, θ)
struct KerrHamiltonianParams {
    float Sigma;    // Σ = r² + a²cos²θ
    float Delta;    // Δ = r² - 2Mr + a²
    float A_func;   // A = (r² + a²)² - a²Δsin²θ
    float sin_th;   // sin(θ)
    float cos_th;   // cos(θ)
    float sin2_th;  // sin²(θ)
    float cos2_th;  // cos²(θ)
};

__device__ __forceinline__ KerrHamiltonianParams computeKerrParams(
    float r, float theta, float M, float a)
{
    KerrHamiltonianParams kp;
    
    kp.sin_th = sinf(theta);
    kp.cos_th = cosf(theta);
    kp.sin2_th = kp.sin_th * kp.sin_th;
    kp.cos2_th = kp.cos_th * kp.cos_th;
    
    float r2 = r * r;
    float a2 = a * a;
    
    kp.Sigma = r2 + a2 * kp.cos2_th;
    kp.Delta = r2 - 2.0f * M * r + a2;
    
    float r2_a2 = r2 + a2;
    kp.A_func = r2_a2 * r2_a2 - a2 * kp.Delta * kp.sin2_th;
    
    return kp;
}

// Compute the time transformation function g = Σ
__device__ __forceinline__ float computeTimeTransformG(
    float r, float theta, float a)
{
    float a2 = a * a;
    float cos_th = cosf(theta);
    return r * r + a2 * cos_th * cos_th;  // Σ
}

//==============================================================================
// HAMILTONIAN ERROR ESTIMATOR
//==============================================================================
__device__ float computeHamiltonianError(
    const HamiltonianState& state, float M, float a)
{
    float r = state.q[1];
    float theta = state.q[2];
    float p_r = state.p[1];
    float p_th = state.p[2];
    float E = -state.p[0];
    float Phi = state.p[3];
    
    KerrHamiltonianParams kp = computeKerrParams(r, theta, M, a);
    
    float a2 = a * a;
    float r2 = r * r;
    float r2_a2 = r2 + a2;
    float Delta_safe = fmaxf(kp.Delta, 0.01f);
    float sin2_safe = fmaxf(kp.sin2_th, 0.0001f);
    
    // Kerr Hamiltonian: H = (1/2Σ)[...] = 0 for null geodesics
    // Compute 2Σ·H to avoid small denominator
    float term1 = -(r2_a2 * E - a * Phi) * (r2_a2 * E - a * Phi) / Delta_safe;
    float term2 = kp.Delta * p_r * p_r;
    float term3 = p_th * p_th;
    float term4 = (a * E - Phi / kp.sin_th) * (a * E - Phi / kp.sin_th) / sin2_safe;
    
    float H_2Sigma = term1 + term2 + term3 + term4;
    
    return fabsf(H_2Sigma) / (1.0f + E * E);  // Normalized error
}

//==============================================================================
// CONVERT BETWEEN GEODESIC STATE AND HAMILTONIAN STATE
//==============================================================================

// GeodesicState (q, u = dq/dλ) → HamiltonianState (q, p = g·u)
__device__ HamiltonianState geodesicToHamiltonian(
    const GeodesicState& gs, float M, float a)
{
    HamiltonianState hs;
    
    // Copy positions
    hs.q[0] = gs.x.t;
    hs.q[1] = gs.x.r;
    hs.q[2] = gs.x.theta;
    hs.q[3] = gs.x.phi;
    
    // Compute metric at position
    float r = gs.x.r;
    float theta = gs.x.theta;
    float sin_th = sinf(theta);
    float cos_th = cosf(theta);
    float r2 = r * r;
    float a2 = a * a;
    
    float Sigma = r2 + a2 * cos_th * cos_th;
    float Delta = r2 - 2.0f * M * r + a2;
    float sin2_th = sin_th * sin_th;
    
    // Covariant metric components for Boyer-Lindquist Kerr
    float g_tt = -(1.0f - 2.0f * M * r / Sigma);
    float g_rr = Sigma / fmaxf(Delta, 0.01f);
    float g_thth = Sigma;
    float g_phph = (r2 + a2 + 2.0f * M * r * a2 * sin2_th / Sigma) * sin2_th;
    float g_tph = -2.0f * M * r * a * sin2_th / Sigma;
    
    // Lower indices: p_μ = g_μν u^ν
    float u[4] = {gs.u.t, gs.u.r, gs.u.theta, gs.u.phi};
    hs.p[0] = g_tt * u[0] + g_tph * u[3];
    hs.p[1] = g_rr * u[1];
    hs.p[2] = g_thth * u[2];
    hs.p[3] = g_tph * u[0] + g_phph * u[3];
    
    return hs;
}

// HamiltonianState (q, p) → GeodesicState (q, u = g^(-1)·p)
__device__ GeodesicState hamiltonianToGeodesic(
    const HamiltonianState& hs, float M, float a)
{
    GeodesicState gs;
    
    // Copy positions
    gs.x.t = hs.q[0];
    gs.x.r = hs.q[1];
    gs.x.theta = hs.q[2];
    gs.x.phi = hs.q[3];
    
    float r = hs.q[1];
    float theta = hs.q[2];
    float sin_th = sinf(theta);
    float cos_th = cosf(theta);
    float r2 = r * r;
    float a2 = a * a;
    
    float Sigma = r2 + a2 * cos_th * cos_th;
    float Delta = fmaxf(r2 - 2.0f * M * r + a2, 0.01f);
    float sin2_th = fmaxf(sin_th * sin_th, 0.0001f);
    
    // Contravariant metric components for Boyer-Lindquist Kerr
    float r2_a2 = r2 + a2;
    float A = r2_a2 * r2_a2 - a2 * Delta * sin2_th;
    
    float g_inv_tt = -A / (Sigma * Delta);
    float g_inv_rr = Delta / Sigma;
    float g_inv_thth = 1.0f / Sigma;
    float g_inv_phph = (Delta - a2 * sin2_th) / (Sigma * Delta * sin2_th);
    float g_inv_tph = -2.0f * M * r * a / (Sigma * Delta);
    
    // Raise indices: u^μ = g^μν p_ν
    gs.u.t = g_inv_tt * hs.p[0] + g_inv_tph * hs.p[3];
    gs.u.r = g_inv_rr * hs.p[1];
    gs.u.theta = g_inv_thth * hs.p[2];
    gs.u.phi = g_inv_tph * hs.p[0] + g_inv_phph * hs.p[3];
    
    return gs;
}

} // namespace Sirius
