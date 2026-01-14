// =============================================================================
// RDOP004A.cu - Ray Bundle Module (Phase 6.3 - DNGR Anti-Aliasing)
// Component ID: RDOP004A
// =============================================================================
//
// MATHEMATICAL FOUNDATION (DNGR Paper Appendix A.2, Eq A.18-A.27)
// ================================================================
//
// Ray bundles track how a pencil of initially parallel light rays evolves
// as it propagates through curved spacetime. This is governed by the
// geodesic deviation equation:
//
//   D²ξ^μ/dλ² = R^μ_νρσ k^ν ξ^ρ k^σ
//
// where:
//   ξ^μ = deviation vector (separation between neighboring rays)
//   k^μ = null geodesic tangent (photon 4-momentum)
//   R^μ_νρσ = Riemann curvature tensor
//   D/dλ = covariant derivative along geodesic
//
// The bundle projects to an ellipse on the celestial sphere, characterized
// by semi-major axis a, semi-minor axis b, and orientation ψ.
//
// REFERENCES:
//   James et al. (2015) "Gravitational Lensing by Spinning Black Holes"
//   Sachs (1961) "Gravitational waves in general relativity VI"
//
// =============================================================================

#pragma once

#include <cuda_runtime.h>
#include <optix.h>
#include "RDOP003A.h"

namespace Sirius {

//==============================================================================
// Riemann Curvature Tensor Computation
// R^μ_νρσ = ∂_ρ Γ^μ_νσ - ∂_σ Γ^μ_νρ + Γ^μ_λρ Γ^λ_νσ - Γ^μ_λσ Γ^λ_νρ
//==============================================================================

// Forward declaration - getChristoffelSymbols is defined in RDOP002A.cu
template<MetricType type>
__device__ void getChristoffelSymbols(const float4& x,
                                       const MetricParams& mp,
                                       float Gamma[4][4][4]);

// Compute Riemann tensor at position x using finite differences of Christoffels
// This is the key curvature quantity needed for geodesic deviation
template<MetricType type>
__device__ void getRiemannTensor(
    const float4& x,
    const MetricParams& mp,
    float R[4][4][4][4])
{
    // Initialize to zero
    for (int mu = 0; mu < 4; mu++)
        for (int nu = 0; nu < 4; nu++)
            for (int rho = 0; rho < 4; rho++)
                for (int sigma = 0; sigma < 4; sigma++)
                    R[mu][nu][rho][sigma] = 0.0f;
    
    // Get Christoffels at current position
    float Gamma[4][4][4];
    getChristoffelSymbols<type>(x, mp, Gamma);
    
    // Finite difference step size (adaptive based on position)
    const float eps = 0.001f * fmaxf(1.0f, x.y);  // x.y = r in spherical
    
    // Compute ∂Γ/∂x^rho by finite differences
    float Gamma_plus[4][4][4], Gamma_minus[4][4][4];
    float dGamma[4][4][4][4]; // dGamma[rho][mu][nu][sigma] = ∂_rho Γ^mu_nu,sigma
    
    for (int rho = 0; rho < 4; rho++) {
        // Perturb in direction rho
        float4 xp = x, xm = x;
        
        switch (rho) {
            case 0: xp.x += eps; xm.x -= eps; break;  // t
            case 1: xp.y += eps; xm.y -= eps; break;  // r
            case 2: xp.z += eps; xm.z -= eps; break;  // theta
            case 3: xp.w += eps; xm.w -= eps; break;  // phi
        }
        
        getChristoffelSymbols<type>(xp, mp, Gamma_plus);
        getChristoffelSymbols<type>(xm, mp, Gamma_minus);
        
        // Central difference
        for (int mu = 0; mu < 4; mu++) {
            for (int nu = 0; nu < 4; nu++) {
                for (int sigma = 0; sigma < 4; sigma++) {
                    dGamma[rho][mu][nu][sigma] = 
                        (Gamma_plus[mu][nu][sigma] - Gamma_minus[mu][nu][sigma]) / (2.0f * eps);
                }
            }
        }
    }
    
    // Compute Riemann tensor
    // R^μ_νρσ = ∂_ρ Γ^μ_νσ - ∂_σ Γ^μ_νρ + Γ^μ_λρ Γ^λ_νσ - Γ^μ_λσ Γ^λ_νρ
    for (int mu = 0; mu < 4; mu++) {
        for (int nu = 0; nu < 4; nu++) {
            for (int rho = 0; rho < 4; rho++) {
                for (int sigma = 0; sigma < 4; sigma++) {
                    // ∂_ρ Γ^μ_νσ - ∂_σ Γ^μ_νρ
                    R[mu][nu][rho][sigma] = dGamma[rho][mu][nu][sigma] 
                                          - dGamma[sigma][mu][nu][rho];
                    
                    // + Γ^μ_λρ Γ^λ_νσ - Γ^μ_λσ Γ^λ_νρ
                    for (int lambda = 0; lambda < 4; lambda++) {
                        R[mu][nu][rho][sigma] += Gamma[mu][lambda][rho] * Gamma[lambda][nu][sigma]
                                               - Gamma[mu][lambda][sigma] * Gamma[lambda][nu][rho];
                    }
                }
            }
        }
    }
}

//==============================================================================
// Geodesic Deviation Acceleration (DNGR Eq A.19)
// D²ξ^μ/dλ² = R^μ_νρσ k^ν ξ^ρ k^σ
//
// This computes the "tidal acceleration" that causes neighboring rays to
// diverge or converge due to spacetime curvature.
//==============================================================================
template<MetricType type>
__device__ float4 geodesicDeviationAccel(
    const float4& x,          // Position
    const float4& k,          // Photon 4-momentum (null geodesic tangent)
    const float4& xi,         // Deviation vector
    const MetricParams& mp)
{
    // Get Riemann tensor
    float R[4][4][4][4];
    getRiemannTensor<type>(x, mp, R);
    
    // Contract: R^μ_νρσ k^ν ξ^ρ k^σ
    float k_arr[4] = {k.x, k.y, k.z, k.w};
    float xi_arr[4] = {xi.x, xi.y, xi.z, xi.w};
    float accel[4] = {0.0f, 0.0f, 0.0f, 0.0f};
    
    for (int mu = 0; mu < 4; mu++) {
        for (int nu = 0; nu < 4; nu++) {
            for (int rho = 0; rho < 4; rho++) {
                for (int sigma = 0; sigma < 4; sigma++) {
                    accel[mu] += R[mu][nu][rho][sigma] * k_arr[nu] * xi_arr[rho] * k_arr[sigma];
                }
            }
        }
    }
    
    return make_float4(accel[0], accel[1], accel[2], accel[3]);
}

//==============================================================================
// Initialize Ray Bundle at Camera
// Sets up initial deviation vectors based on pixel size
//==============================================================================
__device__ RayBundleState initRayBundle(
    const float4& x0,           // Initial position
    const float4& k0,           // Initial 4-momentum
    const float3& screenX,      // Screen X direction (local)
    const float3& screenY,      // Screen Y direction (local)
    float pixelSize)            // Angular size of pixel
{
    RayBundleState bundle;
    
    // Central ray
    bundle.x = x0;
    bundle.u = k0;
    
    // Initial deviation vectors aligned with screen
    // These represent rays at the edges of the pixel
    // In flat space at camera, spatial deviations are simply along screen axes
    bundle.xi1 = make_float4(0.0f, screenX.x * pixelSize, screenX.y * pixelSize, screenX.z * pixelSize);
    bundle.dxi1 = make_float4(0.0f, 0.0f, 0.0f, 0.0f);  // Initially parallel rays
    
    bundle.xi2 = make_float4(0.0f, screenY.x * pixelSize, screenY.y * pixelSize, screenY.z * pixelSize);
    bundle.dxi2 = make_float4(0.0f, 0.0f, 0.0f, 0.0f);
    
    // Initial ellipse is circular (before any curvature effects)
    bundle.semiMajor = pixelSize;
    bundle.semiMinor = pixelSize;
    bundle.orientation = 0.0f;
    
    bundle.affineParam = 0.0f;
    bundle.terminated = false;
    
    return bundle;
}

//==============================================================================
// Compute Ellipse Parameters from Deviation Vectors
// Extracts semi-axes and orientation from ξ₁, ξ₂ (DNGR Eq A.24-A.25)
//==============================================================================
__device__ void computeEllipseParameters(RayBundleState& bundle)
{
    // Project deviation vectors to screen plane (ignore time component)
    float3 xi1_spatial = make_float3(bundle.xi1.y, bundle.xi1.z, bundle.xi1.w);
    float3 xi2_spatial = make_float3(bundle.xi2.y, bundle.xi2.z, bundle.xi2.w);
    
    // Compute 2x2 matrix M = [xi1 | xi2] in screen coords
    // The ellipse axes are eigenvalues of M^T M
    float a11 = dot(xi1_spatial, xi1_spatial);
    float a12 = dot(xi1_spatial, xi2_spatial);
    float a22 = dot(xi2_spatial, xi2_spatial);
    
    // Eigenvalues of symmetric 2x2 matrix [a11 a12; a12 a22]
    float trace = a11 + a22;
    float det = a11 * a22 - a12 * a12;
    float disc = sqrtf(fmaxf(0.0f, trace * trace - 4.0f * det));
    
    float lambda1 = 0.5f * (trace + disc);
    float lambda2 = 0.5f * (trace - disc);
    
    // Semi-axes are sqrt of eigenvalues
    bundle.semiMajor = sqrtf(fmaxf(lambda1, 1e-10f));
    bundle.semiMinor = sqrtf(fmaxf(lambda2, 1e-10f));
    
    // Orientation is direction of major axis
    if (fabsf(a12) > 1e-10f) {
        bundle.orientation = 0.5f * atan2f(2.0f * a12, a11 - a22);
    } else {
        bundle.orientation = (a11 > a22) ? 0.0f : 1.5707963f; // 0 or π/2
    }
}

//==============================================================================
// Geodesic Acceleration from Christoffel Symbols
// du/dλ = -Γ^μ_νρ u^ν u^ρ
//==============================================================================
template<MetricType type>
__device__ __forceinline__ float4 computeGeodesicAccel4(
    const float4& x,
    const float4& u,
    const MetricParams& mp)
{
    float Gamma[4][4][4];
    getChristoffelSymbols<type>(x, mp, Gamma);
    
    float u_arr[4] = {u.x, u.y, u.z, u.w};
    float du_arr[4] = {0.0f, 0.0f, 0.0f, 0.0f};
    for (int mu = 0; mu < 4; mu++) {
        for (int nu = 0; nu < 4; nu++) {
            for (int rho = 0; rho < 4; rho++) {
                du_arr[mu] -= Gamma[mu][nu][rho] * u_arr[nu] * u_arr[rho];
            }
        }
    }
    return make_float4(du_arr[0], du_arr[1], du_arr[2], du_arr[3]);
}

//==============================================================================
// RK4 Integration Step for Ray Bundle
// Integrates both central ray AND deviation vectors
//==============================================================================
template<MetricType type>
__device__ void integrateRayBundleRK4(
    RayBundleState& bundle,
    float dt,
    const MetricParams& mp)
{
    // k1 - evaluate at current state
    float4 dx1 = bundle.u;
    float4 du1 = computeGeodesicAccel4<type>(bundle.x, bundle.u, mp);
    float4 dxi1_1 = bundle.dxi1, dxi2_1 = bundle.dxi2;
    float4 ddxi1_1 = geodesicDeviationAccel<type>(bundle.x, bundle.u, bundle.xi1, mp);
    float4 ddxi2_1 = geodesicDeviationAccel<type>(bundle.x, bundle.u, bundle.xi2, mp);
    
    // k2 - midpoint
    float4 x2 = bundle.x + dx1 * (dt * 0.5f);
    float4 u2 = bundle.u + du1 * (dt * 0.5f);
    float4 xi1_2 = bundle.xi1 + dxi1_1 * (dt * 0.5f);
    float4 xi2_2 = bundle.xi2 + dxi2_1 * (dt * 0.5f);
    float4 dx2 = u2;
    float4 du2 = computeGeodesicAccel4<type>(x2, u2, mp);
    float4 dxi1_2 = bundle.dxi1 + ddxi1_1 * (dt * 0.5f);
    float4 dxi2_2 = bundle.dxi2 + ddxi2_1 * (dt * 0.5f);
    float4 ddxi1_2 = geodesicDeviationAccel<type>(x2, u2, xi1_2, mp);
    float4 ddxi2_2 = geodesicDeviationAccel<type>(x2, u2, xi2_2, mp);
    
    // k3 - second midpoint
    float4 x3 = bundle.x + dx2 * (dt * 0.5f);
    float4 u3 = bundle.u + du2 * (dt * 0.5f);
    float4 xi1_3 = bundle.xi1 + dxi1_2 * (dt * 0.5f);
    float4 xi2_3 = bundle.xi2 + dxi2_2 * (dt * 0.5f);
    float4 dx3 = u3;
    float4 du3 = computeGeodesicAccel4<type>(x3, u3, mp);
    float4 dxi1_3 = bundle.dxi1 + ddxi1_2 * (dt * 0.5f);
    float4 dxi2_3 = bundle.dxi2 + ddxi2_2 * (dt * 0.5f);
    float4 ddxi1_3 = geodesicDeviationAccel<type>(x3, u3, xi1_3, mp);
    float4 ddxi2_3 = geodesicDeviationAccel<type>(x3, u3, xi2_3, mp);
    
    // k4 - endpoint
    float4 x4 = bundle.x + dx3 * dt;
    float4 u4 = bundle.u + du3 * dt;
    float4 xi1_4 = bundle.xi1 + dxi1_3 * dt;
    float4 xi2_4 = bundle.xi2 + dxi2_3 * dt;
    float4 dx4 = u4;
    float4 du4 = computeGeodesicAccel4<type>(x4, u4, mp);
    float4 ddxi1_4 = geodesicDeviationAccel<type>(x4, u4, xi1_4, mp);
    float4 ddxi2_4 = geodesicDeviationAccel<type>(x4, u4, xi2_4, mp);
    
    // Combine RK4 weights
    float inv6 = 1.0f / 6.0f;
    bundle.x = bundle.x + (dx1 + dx2 * 2.0f + dx3 * 2.0f + dx4) * (dt * inv6);
    bundle.u = bundle.u + (du1 + du2 * 2.0f + du3 * 2.0f + du4) * (dt * inv6);
    bundle.xi1 = bundle.xi1 + (dxi1_1 + dxi1_2 * 2.0f + dxi1_3 * 2.0f + bundle.dxi1) * (dt * inv6);
    bundle.dxi1 = bundle.dxi1 + (ddxi1_1 + ddxi1_2 * 2.0f + ddxi1_3 * 2.0f + ddxi1_4) * (dt * inv6);
    bundle.xi2 = bundle.xi2 + (dxi2_1 + dxi2_2 * 2.0f + dxi2_3 * 2.0f + bundle.dxi2) * (dt * inv6);
    bundle.dxi2 = bundle.dxi2 + (ddxi2_1 + ddxi2_2 * 2.0f + ddxi2_3 * 2.0f + ddxi2_4) * (dt * inv6);
    bundle.affineParam += dt;
    
    computeEllipseParameters(bundle);
}

//==============================================================================
// Spatial Filter for Ray Bundle (DNGR Eq A.27)
// Gaussian-weighted average over the ellipse footprint
//==============================================================================
__device__ float3 spatialFilterBundle(
    const RayBundleState& bundle,
    const float3& centralColor,
    cudaTextureObject_t backgroundTexture,
    float filterWidth)
{
    // For now, return central color
    // Full implementation would sample neighboring rays within ellipse
    // and apply Gaussian weighting
    
    // Simple approximation: blend with gaussian-blurred neighbors
    // This is a placeholder - full DNGR samples the ellipse footprint
    return centralColor;
}

} // namespace Sirius
