// =============================================================================
// PHCT001A.h - Coordinate Transformation Utilities
// Component ID: PHCT001A (Physics/Coordinate/Transformations)
// =============================================================================
//
// PURPOSE
// =======
// Provides conversion utilities between Cartesian and spherical coordinates
// for use with the unified metric families. The Kerr-Schild family operates
// in Cartesian coordinates (t, x, y, z), but some operations (background
// texture sampling, accretion disk geometry, UI display) require spherical.
//
// COORDINATE CONVENTIONS
// ======================
// Cartesian: (t, x, y, z) with standard orientation
// Spherical: (t, r, θ, φ) where:
//   - r = √(x² + y² + z²)
//   - θ = arccos(z/r)     (polar angle from +z axis)
//   - φ = atan2(y, x)     (azimuthal angle in xy-plane)
//
// Both CUDA device and host versions are provided.
// =============================================================================

#ifndef PHCT001A_H
#define PHCT001A_H

#include <cmath>

#ifdef __CUDACC__
#define PHCT_DEVICE __device__ __host__
#define PHCT_DEVICE_ONLY __device__
#else
#define PHCT_DEVICE inline
#define PHCT_DEVICE_ONLY inline
#endif

namespace Sirius {
namespace Coordinates {

// =============================================================================
// Cartesian to Spherical Conversion
// =============================================================================

template<typename T>
PHCT_DEVICE void cartesianToSpherical(T x, T y, T z, T& r, T& theta, T& phi) {
    r = std::sqrt(x*x + y*y + z*z);
    
    // Avoid division by zero at origin
    const T r_safe = (r > T(1e-10)) ? r : T(1e-10);
    
    // θ ∈ [0, π] - polar angle from +z axis
    theta = std::acos(z / r_safe);
    
    // φ ∈ [-π, π] - azimuthal angle in xy-plane
    phi = std::atan2(y, x);
}

// 4-vector version (time component unchanged)
template<typename T>
PHCT_DEVICE void cartesianToSpherical4(const T cart[4], T sph[4]) {
    sph[0] = cart[0];  // t unchanged
    cartesianToSpherical(cart[1], cart[2], cart[3], sph[1], sph[2], sph[3]);
}

// =============================================================================
// Spherical to Cartesian Conversion
// =============================================================================

template<typename T>
PHCT_DEVICE void sphericalToCartesian(T r, T theta, T phi, T& x, T& y, T& z) {
    T sin_theta = std::sin(theta);
    T cos_theta = std::cos(theta);
    T sin_phi = std::sin(phi);
    T cos_phi = std::cos(phi);
    
    x = r * sin_theta * cos_phi;
    y = r * sin_theta * sin_phi;
    z = r * cos_theta;
}

// 4-vector version (time component unchanged)
template<typename T>
PHCT_DEVICE void sphericalToCartesian4(const T sph[4], T cart[4]) {
    cart[0] = sph[0];  // t unchanged
    sphericalToCartesian(sph[1], sph[2], sph[3], cart[1], cart[2], cart[3]);
}

// =============================================================================
// Velocity Transformation: Cartesian to Spherical
// =============================================================================
// Transform 4-velocity from Cartesian to spherical coordinates
// Uses Jacobian: u^i_sph = (∂x_sph^i / ∂x_cart^j) u^j_cart

template<typename T>
PHCT_DEVICE void velocityCartToSph(
    T x, T y, T z,           // Position (Cartesian)
    T ut, T ux, T uy, T uz,  // Velocity (Cartesian)
    T& ut_s, T& ur, T& utheta, T& uphi)  // Velocity (Spherical)
{
    T r_cyl = std::sqrt(x*x + y*y);        // Cylindrical radius
    T r = std::sqrt(x*x + y*y + z*z);      // Spherical radius
    
    const T EPS = T(1e-10);
    r = (r > EPS) ? r : EPS;
    r_cyl = (r_cyl > EPS) ? r_cyl : EPS;
    
    T r2 = r * r;
    
    // Time component unchanged
    ut_s = ut;
    
    // u^r = (x·u_x + y·u_y + z·u_z) / r
    ur = (x*ux + y*uy + z*uz) / r;
    
    // u^θ = (z·(x·u_x + y·u_y) - (x²+y²)·u_z) / (r²·√(x²+y²))
    // Simplified: u^θ = (z·r_cyl·u_r_cyl - r_cyl²·u_z) / (r²·r_cyl)
    //           where u_r_cyl = (x·u_x + y·u_y) / r_cyl
    T u_r_cyl = (x*ux + y*uy) / r_cyl;
    utheta = (z * u_r_cyl - r_cyl * uz) / (r * r);
    
    // u^φ = (x·u_y - y·u_x) / (x² + y²)
    uphi = (x*uy - y*ux) / (r_cyl * r_cyl);
}

// =============================================================================
// Velocity Transformation: Spherical to Cartesian
// =============================================================================

template<typename T>
PHCT_DEVICE void velocitySphToCart(
    T r, T theta, T phi,          // Position (Spherical)
    T ut, T ur, T utheta, T uphi, // Velocity (Spherical)
    T& ut_c, T& ux, T& uy, T& uz) // Velocity (Cartesian)
{
    T sin_theta = std::sin(theta);
    T cos_theta = std::cos(theta);
    T sin_phi = std::sin(phi);
    T cos_phi = std::cos(phi);
    
    // Time component unchanged
    ut_c = ut;
    
    // Jacobian: ∂(x,y,z)/∂(r,θ,φ)
    // ∂x/∂r = sinθ cosφ,  ∂x/∂θ = r cosθ cosφ,  ∂x/∂φ = -r sinθ sinφ
    // ∂y/∂r = sinθ sinφ,  ∂y/∂θ = r cosθ sinφ,  ∂y/∂φ = r sinθ cosφ
    // ∂z/∂r = cosθ,       ∂z/∂θ = -r sinθ,      ∂z/∂φ = 0
    
    ux = sin_theta * cos_phi * ur 
       + r * cos_theta * cos_phi * utheta 
       - r * sin_theta * sin_phi * uphi;
    
    uy = sin_theta * sin_phi * ur 
       + r * cos_theta * sin_phi * utheta 
       + r * sin_theta * cos_phi * uphi;
    
    uz = cos_theta * ur 
       - r * sin_theta * utheta;
}

// =============================================================================
// Kerr Radius from Cartesian Coordinates
// =============================================================================
// For Kerr-Schild coordinates, r is defined implicitly by:
//   r⁴ - (x² + y² + z² - a²)r² - a²z² = 0
// This function solves for r given (x, y, z) and spin parameter a.

template<typename T>
PHCT_DEVICE T computeKerrRadius(T x, T y, T z, T a) {
    T a2 = a * a;
    T R2 = x*x + y*y + z*z;
    
    // For a = 0, r = R directly (Schwarzschild)
    if (std::abs(a) < T(1e-12)) {
        return std::sqrt(R2 > T(1e-20) ? R2 : T(1e-20));
    }
    
    // Solve: r² = [(R² - a²) + √((R² - a²)² + 4a²z²)] / 2
    T Rm2 = R2 - a2;
    T disc = Rm2 * Rm2 + T(4.0) * a2 * z * z;
    T r2 = (Rm2 + std::sqrt(disc > T(0) ? disc : T(0))) / T(2.0);
    
    return std::sqrt(r2 > T(1e-20) ? r2 : T(1e-20));
}

// =============================================================================
// Equirectangular Texture Mapping (for background sphere)
// =============================================================================
// Maps a Cartesian ray direction to (u, v) texture coordinates

template<typename T>
PHCT_DEVICE void directionToEquirectangular(T dx, T dy, T dz, T& u, T& v) {
    // Normalize direction
    T len = std::sqrt(dx*dx + dy*dy + dz*dz);
    const T EPS = T(1e-10);
    len = (len > EPS) ? len : EPS;
    dx /= len; dy /= len; dz /= len;
    
    // Convert to spherical
    T phi = std::atan2(dy, dx);                    // [-π, π]
    T theta = std::acos(dz > T(1) ? T(1) : (dz < T(-1) ? T(-1) : dz));  // [0, π]
    
    // Map to [0, 1] texture coordinates
    const T PI = T(3.14159265358979323846);
    u = (phi + PI) / (T(2) * PI);  // [0, 1]
    v = theta / PI;                 // [0, 1]
}

// =============================================================================
// Accretion Disk Geometry Helpers (Cartesian)
// =============================================================================
// Equatorial plane is z = 0

template<typename T>
PHCT_DEVICE bool inAccretionDisk(T x, T y, T z, T r_inner, T r_outer, T height) {
    // Cylindrical radius
    T r_cyl = std::sqrt(x*x + y*y);
    
    // Height above equatorial plane
    T h = std::abs(z);
    
    return (h < height) && (r_cyl > r_inner) && (r_cyl < r_outer);
}

template<typename T>
PHCT_DEVICE T cylindricalRadius(T x, T y) {
    return std::sqrt(x*x + y*y);
}

template<typename T>
PHCT_DEVICE T azimuthalAngle(T x, T y) {
    return std::atan2(y, x);
}

} // namespace Coordinates
} // namespace Sirius

#endif // PHCT001A_H
