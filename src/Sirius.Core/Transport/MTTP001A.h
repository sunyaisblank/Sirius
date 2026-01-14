// MTTP001A.h - Double-Precision Types for Geodesic Integration
// Component ID: MTTP001A
// Purpose: Provides Vec4d, GeodesicStateD, HamiltonianStateD for offline rendering
//
// MATHEMATICAL BASIS:
// These types represent coordinates and momenta in 4D spacetime with double precision
// Required for O(10^-12) conservation over 10^6 integration steps
//
// STRICT REQUIREMENTS:
// - All components must be IEEE 754 double precision (64-bit)
// - Operators must be __host__ __device__ for CPU/GPU compatibility
// - No implicit conversions from float types

#pragma once

#ifdef __CUDACC__
#include <cuda_runtime.h>
#define SIRIUS_HD __host__ __device__
#else
#define SIRIUS_HD
#include <cmath>
#endif

namespace sirius::math {

//==============================================================================
// Vec4d: Double-precision 4-vector for spacetime coordinates
// Indices: 0=t, 1=r, 2=theta, 3=phi (Boyer-Lindquist)
//==============================================================================

struct Vec4d {
    double t, r, theta, phi;
    
    SIRIUS_HD Vec4d() : t(0), r(0), theta(0), phi(0) {}
    
    SIRIUS_HD Vec4d(double t_, double r_, double th_, double ph_)
        : t(t_), r(r_), theta(th_), phi(ph_) {}
    
    // Indexed access (0=t, 1=r, 2=theta, 3=phi)
    SIRIUS_HD double& operator[](int i) {
        return (&t)[i];
    }
    
    SIRIUS_HD const double& operator[](int i) const {
        return (&t)[i];
    }
    
    // Arithmetic operators
    SIRIUS_HD Vec4d operator+(const Vec4d& o) const {
        return Vec4d(t + o.t, r + o.r, theta + o.theta, phi + o.phi);
    }
    
    SIRIUS_HD Vec4d operator-(const Vec4d& o) const {
        return Vec4d(t - o.t, r - o.r, theta - o.theta, phi - o.phi);
    }
    
    SIRIUS_HD Vec4d operator*(double s) const {
        return Vec4d(t * s, r * s, theta * s, phi * s);
    }
    
    SIRIUS_HD Vec4d operator/(double s) const {
        double inv = 1.0 / s;
        return Vec4d(t * inv, r * inv, theta * inv, phi * inv);
    }
    
    SIRIUS_HD Vec4d& operator+=(const Vec4d& o) {
        t += o.t; r += o.r; theta += o.theta; phi += o.phi;
        return *this;
    }
    
    SIRIUS_HD Vec4d& operator-=(const Vec4d& o) {
        t -= o.t; r -= o.r; theta -= o.theta; phi -= o.phi;
        return *this;
    }
    
    SIRIUS_HD Vec4d& operator*=(double s) {
        t *= s; r *= s; theta *= s; phi *= s;
        return *this;
    }
    
    // Negation
    SIRIUS_HD Vec4d operator-() const {
        return Vec4d(-t, -r, -theta, -phi);
    }
    
    // Euclidean norm squared (for debugging, not covariant)
    SIRIUS_HD double norm2() const {
        return t*t + r*r + theta*theta + phi*phi;
    }
    
    // Zero check
    SIRIUS_HD bool isZero() const {
        return t == 0 && r == 0 && theta == 0 && phi == 0;
    }
};

SIRIUS_HD inline Vec4d operator*(double s, const Vec4d& v) {
    return v * s;
}

//==============================================================================
// GeodesicStateD: Complete state for double-precision geodesic
// Contains position, covariant momentum, and conserved quantities
//==============================================================================

struct GeodesicStateD {
    Vec4d x;        // Position: (t, r, θ, φ) in Boyer-Lindquist
    Vec4d k;        // Wave 4-vector: (k_t, k_r, k_θ, k_φ) covariant components
    double lambda;  // Affine parameter
    
    // Killing conserved quantities (computed once at initialisation)
    double E;       // Energy: E = -k_t (for stationary spacetime)
    double Lz;      // z-angular momentum: L_z = k_φ (for axisymmetric spacetime)
    double Q;       // Carter constant (Kerr spacetime only)
    
    SIRIUS_HD GeodesicStateD() 
        : x(), k(), lambda(0), E(0), Lz(0), Q(0) {}
    
    SIRIUS_HD GeodesicStateD(const Vec4d& pos, const Vec4d& mom) 
        : x(pos), k(mom), lambda(0), E(-mom.t), Lz(mom.phi), Q(0) {}
    
    // Compute conserved quantities from current state (for Kerr)
    SIRIUS_HD void computeConservedQuantities(double a) {
        E = -k.t;
        Lz = k.phi;
        // Carter constant: Q = k_θ² + cos²θ(a²(μ² - E²) + L_z²/sin²θ)
        // For null geodesics (μ² = 0):
        double costh = cos(x.theta);
        double sinth = sin(x.theta);
        double sin2th = sinth * sinth;
        Q = k.theta * k.theta + costh * costh * (-a*a*E*E + Lz*Lz/sin2th);
    }
};

//==============================================================================
// HamiltonianStateD: Phase space state for symplectic integration
// Uses generalised coordinates and conjugate momenta
//==============================================================================

struct HamiltonianStateD {
    Vec4d q;    // Generalised coordinates (= x)
    Vec4d p;    // Conjugate momenta (= covariant k)
    double H;   // Hamiltonian value (should be ~0 for null geodesics)
    
    SIRIUS_HD HamiltonianStateD() : q(), p(), H(0) {}
    
    SIRIUS_HD HamiltonianStateD(const Vec4d& q_, const Vec4d& p_)
        : q(q_), p(p_), H(0) {}
    
    // Convert from GeodesicStateD
    SIRIUS_HD explicit HamiltonianStateD(const GeodesicStateD& geo)
        : q(geo.x), p(geo.k), H(0) {}
    
    // Convert to GeodesicStateD
    SIRIUS_HD GeodesicStateD toGeodesicState(double lambda = 0) const {
        GeodesicStateD geo;
        geo.x = q;
        geo.k = p;
        geo.lambda = lambda;
        geo.E = -p.t;
        geo.Lz = p.phi;
        geo.Q = 0;  // Caller must compute if needed
        return geo;
    }
};

//==============================================================================
// Conversion utilities between float and double types
//==============================================================================

// Convert from float4 position (OptiX format) to Vec4d
SIRIUS_HD inline Vec4d toVec4d(float t, float r, float theta, float phi) {
    return Vec4d(static_cast<double>(t), 
                 static_cast<double>(r), 
                 static_cast<double>(theta), 
                 static_cast<double>(phi));
}

// Convert Vec4d back to floats for display
SIRIUS_HD inline void toFloat4(const Vec4d& v, float& t, float& r, float& theta, float& phi) {
    t = static_cast<float>(v.t);
    r = static_cast<float>(v.r);
    theta = static_cast<float>(v.theta);
    phi = static_cast<float>(v.phi);
}

//==============================================================================
// 4x4 Matrix for Jacobian and metric tensor storage
//==============================================================================

struct Mat4d {
    double m[4][4];
    
    SIRIUS_HD Mat4d() {
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
                m[i][j] = 0.0;
    }
    
    // Identity matrix
    SIRIUS_HD static Mat4d identity() {
        Mat4d I;
        for (int i = 0; i < 4; ++i)
            I.m[i][i] = 1.0;
        return I;
    }
    
    // Zero matrix
    SIRIUS_HD static Mat4d zero() {
        return Mat4d();
    }
    
    SIRIUS_HD double& operator()(int i, int j) { return m[i][j]; }
    SIRIUS_HD const double& operator()(int i, int j) const { return m[i][j]; }
    
    // Matrix-vector multiplication
    SIRIUS_HD Vec4d operator*(const Vec4d& v) const {
        Vec4d result;
        for (int i = 0; i < 4; ++i) {
            result[i] = m[i][0]*v.t + m[i][1]*v.r + m[i][2]*v.theta + m[i][3]*v.phi;
        }
        return result;
    }
    
    // Matrix addition
    SIRIUS_HD Mat4d operator+(const Mat4d& o) const {
        Mat4d result;
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
                result.m[i][j] = m[i][j] + o.m[i][j];
        return result;
    }
    
    // Scalar multiplication
    SIRIUS_HD Mat4d operator*(double s) const {
        Mat4d result;
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
                result.m[i][j] = m[i][j] * s;
        return result;
    }
    
    // Determinant
    SIRIUS_HD double determinant() const;
    
    // Trace
    SIRIUS_HD double trace() const {
        return m[0][0] + m[1][1] + m[2][2] + m[3][3];
    }
};

// 4x4 determinant implementation (Laplace expansion)
SIRIUS_HD inline double Mat4d::determinant() const {
    // Expand along first row
    double det = 0;
    for (int j = 0; j < 4; ++j) {
        // 3x3 minor
        double minor[3][3];
        for (int ii = 1; ii < 4; ++ii) {
            int col = 0;
            for (int jj = 0; jj < 4; ++jj) {
                if (jj == j) continue;
                minor[ii-1][col++] = m[ii][jj];
            }
        }
        double det3 = minor[0][0]*(minor[1][1]*minor[2][2] - minor[1][2]*minor[2][1])
                    - minor[0][1]*(minor[1][0]*minor[2][2] - minor[1][2]*minor[2][0])
                    + minor[0][2]*(minor[1][0]*minor[2][1] - minor[1][1]*minor[2][0]);
        det += ((j % 2 == 0) ? 1 : -1) * m[0][j] * det3;
    }
    return det;
}

} // namespace sirius::math
