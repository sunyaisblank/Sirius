// RDOP005A.h - GPU Debug Buffer Interface
// Component ID: RDOP005A
// Purpose: Debug buffer interface for GPU conservation validation
//
// MATHEMATICAL BASIS:
// This module provides GPU-accessible buffers for validating conservation laws:
// - Hamiltonian constraint: H = ½g^μν p_μ p_ν = 0 (null condition)
// - Killing energy: E = -g_tμ k^μ (conserved for stationary spacetimes)
// - Killing angular momentum: L = g_φμ k^μ (conserved for axisymmetric spacetimes)
//
// SPECIFICATION TARGETS (docs/specification.md):
// - Null condition GPU: |H| < 10^-5
// - Energy conservation GPU: |ΔE/E| < 10^-4
// - Angular momentum GPU: |ΔL/L| < 10^-4
//
// TESTS: TSDG005A.cpp

#pragma once

#include <cstdint>

#ifdef __CUDACC__
#include <cuda_runtime.h>
#else
struct float4 { float x, y, z, w; };
using CUdeviceptr = unsigned long long;
#endif

namespace Sirius {

//==============================================================================
// Debug Ray State - Intermediate values during GPU geodesic integration
//==============================================================================
struct DebugRayState {
    // Position and momentum at sample point
    float4 position;      // (t, x, y, z) in Kerr-Schild Cartesian
    float4 momentum;      // Covariant 4-momentum p_μ
    float4 velocity;      // Contravariant 4-velocity k^μ = dx^μ/dλ

    // Conservation quantities
    float hamiltonian;    // H = ½g^μν p_μ p_ν (should be ~0)
    float killing_energy; // E = -g_tμ k^μ
    float killing_L;      // L = g_φμ k^μ (z-component of angular momentum)
    float affine_param;   // λ - current affine parameter

    // Metric components at this point
    float g_tt;           // g_{00}
    float g_tphi;         // g_{03} (frame-dragging)
    float g_phiphi;       // g_{33}
    float radius;         // Coordinate radius r

    // Step information
    int step_number;
    int terminated;       // 0 = active, 1 = horizon, 2 = escaped, -1 = error
    float step_size;
    float error_estimate; // Local truncation error
};

//==============================================================================
// Debug Buffer Configuration
//==============================================================================
struct DebugBufferConfig {
    // Enable debug output (impacts performance)
    bool enabled;

    // Sampling configuration
    int sample_interval;      // Record every N steps (1 = all steps)
    int max_samples_per_ray;  // Maximum samples to store per ray
    int num_debug_rays;       // Number of rays to track (0 = all)

    // Ray selection (for targeted debugging)
    int debug_pixel_x;        // -1 = all pixels
    int debug_pixel_y;        // -1 = all pixels

    // Output buffer
    DebugRayState* buffer;    // Device pointer to output buffer
    int buffer_capacity;      // Total capacity in DebugRayState entries
    int* buffer_count;        // Device pointer to atomic counter
};

//==============================================================================
// Conservation Statistics (Host-side reduction)
//==============================================================================
struct ConservationStats {
    // Null condition statistics
    double max_hamiltonian;       // Maximum |H| observed
    double mean_hamiltonian;      // Mean |H|
    double rms_hamiltonian;       // RMS |H|

    // Energy conservation
    double max_energy_drift;      // Maximum |ΔE/E|
    double mean_energy_drift;
    double rms_energy_drift;

    // Angular momentum conservation
    double max_L_drift;           // Maximum |ΔL/L|
    double mean_L_drift;
    double rms_L_drift;

    // Integration statistics
    int total_rays;
    int terminated_horizon;
    int terminated_escaped;
    int terminated_error;
    int total_steps;
    double mean_steps_per_ray;
};

//==============================================================================
// GPU Debug Functions (Device-side)
//==============================================================================
#ifdef __CUDACC__

/// @brief Compute Hamiltonian (null condition check)
/// @param pos Position (t, x, y, z)
/// @param mom Covariant momentum p_μ
/// @param g_inv Inverse metric tensor
/// @return H = ½g^μν p_μ p_ν (should be ~0 for null rays)
__device__ __forceinline__ float computeHamiltonian(
    float4 mom,
    const float g_inv[4][4])
{
    float H = 0.0f;
    float p[4] = {mom.x, mom.y, mom.z, mom.w};
    for (int mu = 0; mu < 4; ++mu) {
        for (int nu = 0; nu < 4; ++nu) {
            H += g_inv[mu][nu] * p[mu] * p[nu];
        }
    }
    return 0.5f * H;
}

/// @brief Compute Killing energy E = -g_tμ k^μ
/// @param vel Contravariant 4-velocity k^μ
/// @param g Covariant metric tensor
/// @return Killing energy (conserved for stationary spacetimes)
__device__ __forceinline__ float computeKillingEnergy(
    float4 vel,
    const float g[4][4])
{
    // E = -g_{t μ} k^μ = -(g_{tt} k^t + g_{tx} k^x + g_{ty} k^y + g_{tz} k^z)
    // For Kerr-Schild: g_{ti} = H l_t l_i where l is null
    float v[4] = {vel.x, vel.y, vel.z, vel.w};
    float E = 0.0f;
    for (int mu = 0; mu < 4; ++mu) {
        E -= g[0][mu] * v[mu];
    }
    return E;
}

/// @brief Compute z-component of Killing angular momentum
/// For axisymmetric spacetimes: L_z = x p_y - y p_x (in Cartesian)
/// @param pos Position (t, x, y, z)
/// @param vel Contravariant 4-velocity k^μ
/// @param g Covariant metric tensor
/// @return Angular momentum L_z
__device__ __forceinline__ float computeKillingAngularMomentum(
    float4 pos,
    float4 vel,
    const float g[4][4])
{
    // Killing vector for axisymmetry: ξ^μ = (0, -y, x, 0) in Cartesian
    // L = g_{μν} ξ^μ k^ν
    float xi[4] = {0.0f, -pos.z, pos.y, 0.0f};  // (t, x, y, z) -> ξ = (0, -y, x, 0)
    float v[4] = {vel.x, vel.y, vel.z, vel.w};

    float L = 0.0f;
    for (int mu = 0; mu < 4; ++mu) {
        for (int nu = 0; nu < 4; ++nu) {
            L += g[mu][nu] * xi[mu] * v[nu];
        }
    }
    return L;
}

/// @brief Record debug state for a ray
/// Thread-safe atomic insertion into debug buffer
__device__ __forceinline__ void recordDebugState(
    const DebugBufferConfig& config,
    const DebugRayState& state)
{
    if (!config.enabled || config.buffer == nullptr) return;

    // Atomic increment of buffer count
    int idx = atomicAdd(config.buffer_count, 1);

    // Check capacity
    if (idx < config.buffer_capacity) {
        config.buffer[idx] = state;
    }
}

#endif // __CUDACC__

//==============================================================================
// Host-side Analysis Functions
//==============================================================================
#ifndef __CUDACC__

/// @brief Compute conservation statistics from debug buffer
/// @param buffer Host copy of debug buffer
/// @param count Number of valid entries
/// @return Statistics structure
inline ConservationStats analyzeConservation(
    const DebugRayState* buffer,
    int count)
{
    ConservationStats stats = {};
    if (count == 0) return stats;

    // First pass: find initial values for each ray
    // (Assumes rays are stored contiguously)
    double sum_H = 0, sum_H2 = 0;
    double sum_dE = 0, sum_dE2 = 0;
    double sum_dL = 0, sum_dL2 = 0;

    int current_ray = -1;
    double E_initial = 0, L_initial = 0;

    for (int i = 0; i < count; ++i) {
        const DebugRayState& s = buffer[i];

        // Track Hamiltonian
        double H = std::abs(s.hamiltonian);
        sum_H += H;
        sum_H2 += H * H;
        stats.max_hamiltonian = std::max(stats.max_hamiltonian, H);

        // Detect new ray (step 0 or different pixel)
        if (s.step_number == 0) {
            E_initial = s.killing_energy;
            L_initial = s.killing_L;
            current_ray++;
            stats.total_rays++;
        }

        // Track energy drift
        if (std::abs(E_initial) > 1e-10) {
            double dE = std::abs((s.killing_energy - E_initial) / E_initial);
            sum_dE += dE;
            sum_dE2 += dE * dE;
            stats.max_energy_drift = std::max(stats.max_energy_drift, dE);
        }

        // Track angular momentum drift
        if (std::abs(L_initial) > 1e-10) {
            double dL = std::abs((s.killing_L - L_initial) / L_initial);
            sum_dL += dL;
            sum_dL2 += dL * dL;
            stats.max_L_drift = std::max(stats.max_L_drift, dL);
        }

        // Termination statistics
        stats.total_steps++;
        if (s.terminated == 1) stats.terminated_horizon++;
        else if (s.terminated == 2) stats.terminated_escaped++;
        else if (s.terminated < 0) stats.terminated_error++;
    }

    // Compute means and RMS
    stats.mean_hamiltonian = sum_H / count;
    stats.rms_hamiltonian = std::sqrt(sum_H2 / count);

    stats.mean_energy_drift = sum_dE / count;
    stats.rms_energy_drift = std::sqrt(sum_dE2 / count);

    stats.mean_L_drift = sum_dL / count;
    stats.rms_L_drift = std::sqrt(sum_dL2 / count);

    if (stats.total_rays > 0) {
        stats.mean_steps_per_ray = static_cast<double>(stats.total_steps) / stats.total_rays;
    }

    return stats;
}

#endif // !__CUDACC__

} // namespace Sirius
