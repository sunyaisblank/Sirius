// RDGeodesic.cuh - Geodesic State Structures
// Component ID: RDTR001A (Transport/Geodesic)

#pragma once

#include "../Integration/RDMath.cuh"

namespace Sirius {

//==============================================================================
// 4-Vector Structure for Geodesic Integration (Spherical: t, r, θ, φ)
//==============================================================================
struct Vec4 {
    float t, r, theta, phi;
    
    __device__ Vec4() : t(0), r(0), theta(0), phi(0) {}
    __device__ Vec4(float t_, float r_, float theta_, float phi_) 
        : t(t_), r(r_), theta(theta_), phi(phi_) {}
    
    __device__ Vec4 operator+(const Vec4& v) const {
        return Vec4(t + v.t, r + v.r, theta + v.theta, phi + v.phi);
    }
    
    __device__ Vec4 operator*(float s) const {
        return Vec4(t * s, r * s, theta * s, phi * s);
    }
};

//==============================================================================
// 4-Vector Structure for Cartesian Coordinates (t, x, y, z)
// Used by Kerr-Schild unified metric family
//==============================================================================
struct Vec4Cart {
    float t, x, y, z;
    
    __device__ Vec4Cart() : t(0), x(0), y(0), z(0) {}
    __device__ Vec4Cart(float t_, float x_, float y_, float z_) 
        : t(t_), x(x_), y(y_), z(z_) {}
    
    __device__ Vec4Cart operator+(const Vec4Cart& v) const {
        return Vec4Cart(t + v.t, x + v.x, y + v.y, z + v.z);
    }
    
    __device__ Vec4Cart operator*(float s) const {
        return Vec4Cart(t * s, x * s, y * s, z * s);
    }
};

// Additional Vec4Cart operators
__device__ __forceinline__ Vec4Cart operator*(float s, const Vec4Cart& v) {
    return v * s;
}
__device__ __forceinline__ Vec4Cart operator-(const Vec4Cart& a, const Vec4Cart& b) {
    return Vec4Cart(a.t - b.t, a.x - b.x, a.y - b.y, a.z - b.z);
}

//==============================================================================
// Coordinate Conversions
//==============================================================================
__device__ __forceinline__ Vec4Cart vec4SphToCart(const Vec4& sph) {
    float sin_theta, cos_theta, sin_phi, cos_phi;
    sincosf(sph.theta, &sin_theta, &cos_theta);
    sincosf(sph.phi, &sin_phi, &cos_phi);
    return Vec4Cart(
        sph.t,
        sph.r * sin_theta * cos_phi,
        sph.r * sin_theta * sin_phi,
        sph.r * cos_theta
    );
}

__device__ __forceinline__ Vec4 vec4CartToSph(const Vec4Cart& cart) {
    float r = sqrtf(cart.x * cart.x + cart.y * cart.y + cart.z * cart.z);
    float r_safe = fmaxf(r, 1e-10f);
    float theta = acosf(cart.z / r_safe);
    float phi = atan2f(cart.y, cart.x);
    return Vec4(cart.t, r, theta, phi);
}

//==============================================================================
// GeodesicState: Position and Velocity in Boyer-Lindquist coordinates
//==============================================================================
struct GeodesicState {
    Vec4 x;   // Position (t, r, θ, φ)
    Vec4 u;   // 4-velocity (dt/dλ, dr/dλ, dθ/dλ, dφ/dλ)
};

//==============================================================================
// GeodesicStateCart: Position and Velocity in Cartesian coordinates
//==============================================================================
struct GeodesicStateCart {
    Vec4Cart x;   // Position (t, x, y, z)
    Vec4Cart u;   // 4-velocity (dt/dλ, dx/dλ, dy/dλ, dz/dλ)
};

//==============================================================================
// Hamiltonian State: Canonical coordinates (q, p)
//==============================================================================
struct HamiltonianState {
    float q[4];   // Positions: t, r, θ, φ
    float p[4];   // Momenta: p_t, p_r, p_θ, p_φ
    
    __device__ HamiltonianState() {
        for (int i = 0; i < 4; i++) { q[i] = 0.0f; p[i] = 0.0f; }
    }
};

} // namespace Sirius
