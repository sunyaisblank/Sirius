// RDIntegratorSymplectic.cuh - Symplectic Integrator for Kerr
// Component ID: RDIN002A (Integration/Symplectic)

#pragma once

#include "../Integration/RDMath.cuh"
#include "../Transport/RDMetricKerr.cuh"

namespace Sirius {

//==============================================================================
// Yoshida 4th-Order Symplectic Coefficients
// Composition: S₄ = S₂(c₁h) ∘ S₂(c₂h) ∘ S₂(c₃h) ∘ S₂(c₂h) ∘ S₂(c₁h)
//==============================================================================
__device__ constexpr float YOSHIDA_W0 = -1.7024143839193153f;  // 2^(1/3)/(2-2^(1/3))
__device__ constexpr float YOSHIDA_W1 = 1.3512071919596578f;   // 1/(2-2^(1/3))
// c1 = w1/2, c2 = (w0+w1)/2, c3 = c2
__device__ constexpr float YOSHIDA_C1 = 0.6756035959798289f;
__device__ constexpr float YOSHIDA_C2 = -0.17560359597982885f;
__device__ constexpr float YOSHIDA_C3 = -0.17560359597982885f;
// d1 = w1, d2 = w0, d3 = w1
__device__ constexpr float YOSHIDA_D1 = 1.3512071919596578f;
__device__ constexpr float YOSHIDA_D2 = -1.7024143839193153f;
__device__ constexpr float YOSHIDA_D3 = 1.3512071919596578f;

//==============================================================================
// Flows for Separable Hamiltonian K = Σ K_i
//==============================================================================

// K₂ flow: Radial kinetic term K₂ = (1/2)Δ(r)·p_r²
__device__ __forceinline__ void flowK2_radial(
    HamiltonianState& state, float h, float M, float a)
{
    float r = state.q[1];
    float p_r = state.p[1];
    float a2 = a * a;
    float dDelta_dr = 2.0f * r - 2.0f * M;
    float Delta = r * r - 2.0f * M * r + a2;
    
    // Half-step momentum kick from potential
    state.p[1] -= 0.5f * h * (0.5f * dDelta_dr * p_r * p_r);
    
    // Full-step position drift
    state.q[1] += h * Delta * state.p[1];
    
    // Half-step momentum kick
    r = state.q[1];  // Updated position
    dDelta_dr = 2.0f * r - 2.0f * M;
    state.p[1] -= 0.5f * h * (0.5f * dDelta_dr * state.p[1] * state.p[1]);
}

// K₃ flow: Angular kinetic term K₃ = (1/2)p_θ²
__device__ __forceinline__ void flowK3_angular(
    HamiltonianState& state, float h)
{
    state.q[2] += h * state.p[2];
    
    // Trans-polar passage
    const float PI_VAL = 3.14159265359f;
    const float TWO_PI_VAL = 6.28318530718f;
    
    if (state.q[2] < 0.0f) {
        state.q[2] = -state.q[2];
        state.p[2] = -state.p[2];
        state.q[3] += PI_VAL;
    } 
    else if (state.q[2] > PI_VAL) {
        state.q[2] = TWO_PI_VAL - state.q[2];
        state.p[2] = -state.p[2];
        state.q[3] += PI_VAL;
    }
}

// K₄ flow: Effective radial potential V_r(r, E, Φ)
__device__ __forceinline__ void flowK4_radialPotential(
    HamiltonianState& state, float h, float M, float a, float E, float Phi)
{
    float r = state.q[1];
    float r2 = r * r;
    float a2 = a * a;
    float Delta = r2 - 2.0f * M * r + a2;
    float Delta_safe = fmaxf(Delta, 0.01f);
    float r2_a2 = r2 + a2;
    float term1 = E * r2_a2 - a * Phi;
    float dDelta_dr = 2.0f * r - 2.0f * M;
    
    float dVr_dr = -2.0f * E * 2.0f * r * term1 / Delta_safe 
                  + term1 * term1 * dDelta_dr / (Delta_safe * Delta_safe);
    
    state.p[1] -= h * dVr_dr;
}

// K₅ flow: Effective angular potential V_θ(θ, E, Φ)
__device__ __forceinline__ void flowK5_angularPotential(
    HamiltonianState& state, float h, float a, float E, float Phi)
{
    float theta = state.q[2];
    float sin_th = sinf(theta);
    float cos_th = cosf(theta);
    float sin2_safe = fmaxf(sin_th * sin_th, 0.0001f);
    float a2 = a * a;
    
    float dVth_dth = 2.0f * a2 * cos_th * sin_th * E * E 
                   - 2.0f * Phi * Phi * cos_th / (sin2_safe * sin_th);
    
    state.p[2] -= h * dVth_dth;
}

// t and φ evolution (from time transformation)
__device__ __forceinline__ void flowTimeAzimuth(
    HamiltonianState& state, float h, float M, float a)
{
    float r = state.q[1];
    float theta = state.q[2];
    
    KerrHamiltonianParams kp = computeKerrParams(r, theta, M, a);
    
    float E = -state.p[0];
    float Phi = state.p[3];
    float r2_a2 = r * r + a * a;
    float Delta_safe = fmaxf(kp.Delta, 0.01f);
    float sin2_safe = fmaxf(kp.sin2_th, 0.0001f);
    
    float dH_dpt = (r2_a2 * r2_a2 / (kp.Sigma * Delta_safe)) * E 
                 - (a / sin2_safe) * Phi;
    
    float dH_dphi = (a * r2_a2 / (kp.Sigma * Delta_safe)) * E 
                  + (1.0f / (kp.Sigma * sin2_safe)) * Phi;
    
    state.q[0] += h * kp.Sigma * dH_dpt;   // t evolution
    state.q[3] += h * kp.Sigma * dH_dphi;  // φ evolution
}

//------------------------------------------------------------------------------
// S₂ Integrator
//------------------------------------------------------------------------------
__device__ void integrateSymplecticS2(
    HamiltonianState& state, float h, float M, float a)
{
    float E = -state.p[0];
    float Phi = state.p[3];
    
    flowK4_radialPotential(state, 0.5f * h, M, a, E, Phi);
    flowK5_angularPotential(state, 0.5f * h, a, E, Phi);
    
    flowK2_radial(state, h, M, a);
    flowK3_angular(state, h);
    
    flowK4_radialPotential(state, 0.5f * h, M, a, E, Phi);
    flowK5_angularPotential(state, 0.5f * h, a, E, Phi);
    
    flowTimeAzimuth(state, h, M, a);
}

//------------------------------------------------------------------------------
// S₄ Integrator
//------------------------------------------------------------------------------
__device__ void integrateSymplecticS4(
    HamiltonianState& state, float h, float M, float a)
{
    const float c1 = YOSHIDA_C1; // Actually c1 = 1/(2-2^(1/3)) is YOSHIDA_W1/2 ?
    // RDOP002A.cu used hardcoded values:
    // c1 = 1/(2-2^(1/3)) -> This is W1?
    // Let's use the values from RDOP002A.cu to be safe
    const float k_c1 = 1.3512071919596578f;
    const float k_c2 = -1.7024143839193153f;
    
    integrateSymplecticS2(state, k_c1 * h, M, a);
    integrateSymplecticS2(state, k_c2 * h, M, a);
    integrateSymplecticS2(state, k_c1 * h, M, a);
}

// Adaptive step
__device__ float computeAdaptiveSymplecticStep(
    const HamiltonianState& state, float tau_step, float M, float a,
    float minStep, float maxStep)
{
    float r = state.q[1];
    float theta = state.q[2];
    float g = computeTimeTransformG(r, theta, a);
    
    float r2 = r * r;
    float Delta = r2 - 2.0f * M * r + a * a;
    float horizon_factor = fmaxf(Delta / (r2 + 0.1f), 0.1f);
    
    float dt = g * tau_step * horizon_factor;
    return clamp(dt, minStep, maxStep);
}

//------------------------------------------------------------------------------
// Main Interface
//------------------------------------------------------------------------------
__device__ GeodesicState integrateGeodesicSymplectic(
    const GeodesicState& state,
    float& h,               // IN/OUT: step size
    float M, float a,
    float tolerance,
    float minStep,
    float maxStep,
    bool& stepAccepted)
{
    HamiltonianState hs = geodesicToHamiltonian(state, M, a);
    float H_error_initial = computeHamiltonianError(hs, M, a);
    
    float tau_step = 0.1f;
    float effective_h = computeAdaptiveSymplecticStep(hs, tau_step, M, a, minStep, maxStep);
    
    HamiltonianState hs_new = hs;
    integrateSymplecticS4(hs_new, effective_h, M, a);
    
    float H_error_final = computeHamiltonianError(hs_new, M, a);
    float H_error_change = fabsf(H_error_final - H_error_initial);
    
    stepAccepted = (H_error_change < tolerance * 10.0f);
    

    if (!stepAccepted) {
        h = fmaxf(h * 0.5f, minStep);
        return state;
    }
    
    // Adapt step size (heuristically update h, though effective_h depends on pos)
    if (H_error_change < tolerance * 0.1f) {
        // h = fminf(h * 1.5f, maxStep); // This line in original was seemingly overwritten
    }
    h = effective_h; // Report the step we actually used
    
    GeodesicState result = hamiltonianToGeodesic(hs_new, M, a);
    
    // Normalize angles
    // Use PI from RDMath.cuh
    const float THETA_MIN = 0.001f;
    const float THETA_MAX = PI - 0.001f;
    result.x.theta = fmaxf(THETA_MIN, fminf(THETA_MAX, result.x.theta));
    result.x.phi = fmodf(result.x.phi, TWO_PI);
    if (result.x.phi < 0.0f) result.x.phi += TWO_PI;
    
    return result;
}


} // namespace Sirius
