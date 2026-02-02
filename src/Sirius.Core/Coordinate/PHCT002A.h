// =============================================================================
// PHCT002A.h - Coordinate Transformation Module
// Component ID: PHCT002A (Physics/Coordinate Transform)
// =============================================================================
//
// PURPOSE:
// Explicit, validated coordinate transformations between:
// - Boyer-Lindquist (BL) coordinates: (t, r, θ, φ) - natural for Kerr physics
// - Kerr-Schild Cartesian (KS): (t, x, y, z) - natural for GPU integration
// - Observer frame (spherical): used for camera and background sampling
//
// MATHEMATICAL FOUNDATION:
// The transformation between BL and KS Cartesian is:
//   x = √(r² + a²) sin(θ) cos(φ')
//   y = √(r² + a²) sin(θ) sin(φ')
//   z = r cos(θ)
// where φ' = φ + ∫ a/(r² - 2Mr + a²) dr for the full transformation
//
// For practical purposes at large r, we use the simplified form:
//   x = r sin(θ) cos(φ)
//   y = r sin(θ) sin(φ)
//   z = r cos(θ)
//
// The Jacobian matrices are provided for proper tensor transformation.
//
// TESTS: TSPH006A.cpp (coordinate transformation tests)
// =============================================================================

#ifndef PHCT002A_H
#define PHCT002A_H

#include "MTTN001A.h"
#include <cmath>
#include <array>

namespace Sirius {
namespace Coordinates {

// =============================================================================
// Coordinate Systems Enumeration
// =============================================================================

enum class CoordinateSystem {
    BoyerLindquist,    // (t, r, θ, φ) - spherical-like
    KerrSchildCart,    // (t, x, y, z) - Cartesian
    Spherical,         // (t, r, θ, φ) - standard spherical (a=0)
    Cartesian          // (t, x, y, z) - flat space Cartesian
};

// =============================================================================
// 4-Vector Types for Different Coordinate Systems
// =============================================================================

/// @brief Boyer-Lindquist coordinates (t, r, θ, φ)
struct Vec4BL {
    double t = 0;      ///< Coordinate time
    double r = 0;      ///< Boyer-Lindquist radius
    double theta = 0;  ///< Polar angle [0, π]
    double phi = 0;    ///< Azimuthal angle [0, 2π)

    Vec4BL() = default;
    Vec4BL(double t_, double r_, double th_, double ph_)
        : t(t_), r(r_), theta(th_), phi(ph_) {}

    double& operator[](int i) {
        switch(i) {
            case 0: return t;
            case 1: return r;
            case 2: return theta;
            default: return phi;
        }
    }
    double operator[](int i) const {
        switch(i) {
            case 0: return t;
            case 1: return r;
            case 2: return theta;
            default: return phi;
        }
    }
};

/// @brief Kerr-Schild Cartesian coordinates (t, x, y, z)
struct Vec4Cart {
    double t = 0;  ///< Coordinate time
    double x = 0;  ///< Cartesian x
    double y = 0;  ///< Cartesian y
    double z = 0;  ///< Cartesian z

    Vec4Cart() = default;
    Vec4Cart(double t_, double x_, double y_, double z_)
        : t(t_), x(x_), y(y_), z(z_) {}

    double& operator[](int i) {
        switch(i) {
            case 0: return t;
            case 1: return x;
            case 2: return y;
            default: return z;
        }
    }
    double operator[](int i) const {
        switch(i) {
            case 0: return t;
            case 1: return x;
            case 2: return y;
            default: return z;
        }
    }

    double radius() const { return std::sqrt(x*x + y*y + z*z); }
};

// =============================================================================
// Jacobian Matrix Type
// =============================================================================

/// @brief 4x4 Jacobian matrix for coordinate transformation
/// J[μ][ν] = ∂x'^μ / ∂x^ν
using Jacobian4x4 = std::array<std::array<double, 4>, 4>;

// =============================================================================
// Coordinate Transformation Functions
// =============================================================================

/// @brief Transform Boyer-Lindquist to Cartesian (simplified for a=0 or large r)
/// @param bl Boyer-Lindquist coordinates
/// @return Cartesian coordinates
inline Vec4Cart BLToCartesian(const Vec4BL& bl) {
    Vec4Cart cart;
    cart.t = bl.t;

    double sin_theta = std::sin(bl.theta);
    double cos_theta = std::cos(bl.theta);
    double sin_phi = std::sin(bl.phi);
    double cos_phi = std::cos(bl.phi);

    cart.x = bl.r * sin_theta * cos_phi;
    cart.y = bl.r * sin_theta * sin_phi;
    cart.z = bl.r * cos_theta;

    return cart;
}

/// @brief Transform Boyer-Lindquist to Kerr-Schild Cartesian (full Kerr)
/// @param bl Boyer-Lindquist coordinates
/// @param a Kerr spin parameter
/// @return Kerr-Schild Cartesian coordinates
inline Vec4Cart BLToKerrSchildCart(const Vec4BL& bl, double a) {
    Vec4Cart cart;
    cart.t = bl.t;

    double sin_theta = std::sin(bl.theta);
    double cos_theta = std::cos(bl.theta);
    double sin_phi = std::sin(bl.phi);
    double cos_phi = std::cos(bl.phi);

    // For Kerr: the oblate spheroidal coordinates have
    // x = √(r² + a²) sin(θ) cos(φ)
    // y = √(r² + a²) sin(θ) sin(φ)
    // z = r cos(θ)
    double rho = std::sqrt(bl.r * bl.r + a * a);

    cart.x = rho * sin_theta * cos_phi;
    cart.y = rho * sin_theta * sin_phi;
    cart.z = bl.r * cos_theta;

    return cart;
}

/// @brief Transform Cartesian to Boyer-Lindquist (simplified for a=0)
/// @param cart Cartesian coordinates
/// @return Boyer-Lindquist coordinates
inline Vec4BL CartesianToBL(const Vec4Cart& cart) {
    Vec4BL bl;
    bl.t = cart.t;

    bl.r = std::sqrt(cart.x*cart.x + cart.y*cart.y + cart.z*cart.z);
    if (bl.r < 1e-15) {
        bl.theta = 0;
        bl.phi = 0;
        return bl;
    }

    bl.theta = std::acos(std::clamp(cart.z / bl.r, -1.0, 1.0));
    bl.phi = std::atan2(cart.y, cart.x);
    if (bl.phi < 0) bl.phi += 2.0 * M_PI;

    return bl;
}

/// @brief Transform Kerr-Schild Cartesian to Boyer-Lindquist (full Kerr)
/// Solves: r⁴ - (R² - a²)r² - a²z² = 0
/// @param cart Kerr-Schild Cartesian coordinates
/// @param a Kerr spin parameter
/// @return Boyer-Lindquist coordinates
inline Vec4BL KerrSchildCartToBL(const Vec4Cart& cart, double a) {
    Vec4BL bl;
    bl.t = cart.t;

    double x = cart.x;
    double y = cart.y;
    double z = cart.z;
    double a2 = a * a;

    // Solve for r: r⁴ - (R² - a²)r² - a²z² = 0
    double R2 = x*x + y*y + z*z;

    if (std::abs(a) < 1e-12) {
        // Schwarzschild case
        bl.r = std::sqrt(R2);
    } else {
        double Rm2 = R2 - a2;
        double disc = Rm2*Rm2 + 4*a2*z*z;
        double r2 = (Rm2 + std::sqrt(disc)) / 2.0;
        bl.r = std::sqrt(std::max(r2, 1e-20));
    }

    if (bl.r < 1e-15) {
        bl.theta = M_PI / 2.0;
        bl.phi = 0;
        return bl;
    }

    // θ = arccos(z/r)
    bl.theta = std::acos(std::clamp(z / bl.r, -1.0, 1.0));

    // φ = atan2(y, x) for simplified transform
    // Full Kerr would need: φ_BL = atan2(y, x) - atan2(a, r)
    bl.phi = std::atan2(y, x);
    if (std::abs(a) > 1e-12) {
        bl.phi -= std::atan2(a, bl.r);
    }
    while (bl.phi < 0) bl.phi += 2.0 * M_PI;
    while (bl.phi >= 2.0 * M_PI) bl.phi -= 2.0 * M_PI;

    return bl;
}

// =============================================================================
// Jacobian Matrices
// =============================================================================

/// @brief Compute Jacobian J[μ][ν] = ∂x^μ_Cart / ∂x^ν_BL
/// For transforming contravariant vectors: V'^μ = J[μ][ν] V^ν
inline Jacobian4x4 JacobianBLToCartesian(const Vec4BL& bl) {
    Jacobian4x4 J = {};

    double sin_th = std::sin(bl.theta);
    double cos_th = std::cos(bl.theta);
    double sin_ph = std::sin(bl.phi);
    double cos_ph = std::cos(bl.phi);
    double r = bl.r;

    // ∂t'/∂t = 1
    J[0][0] = 1.0;

    // ∂x/∂r = sin(θ)cos(φ)
    // ∂x/∂θ = r cos(θ)cos(φ)
    // ∂x/∂φ = -r sin(θ)sin(φ)
    J[1][1] = sin_th * cos_ph;
    J[1][2] = r * cos_th * cos_ph;
    J[1][3] = -r * sin_th * sin_ph;

    // ∂y/∂r = sin(θ)sin(φ)
    // ∂y/∂θ = r cos(θ)sin(φ)
    // ∂y/∂φ = r sin(θ)cos(φ)
    J[2][1] = sin_th * sin_ph;
    J[2][2] = r * cos_th * sin_ph;
    J[2][3] = r * sin_th * cos_ph;

    // ∂z/∂r = cos(θ)
    // ∂z/∂θ = -r sin(θ)
    // ∂z/∂φ = 0
    J[3][1] = cos_th;
    J[3][2] = -r * sin_th;
    J[3][3] = 0.0;

    return J;
}

/// @brief Compute inverse Jacobian J^{-1}[μ][ν] = ∂x^μ_BL / ∂x^ν_Cart
inline Jacobian4x4 JacobianCartesianToBL(const Vec4Cart& cart) {
    Jacobian4x4 Jinv = {};

    double x = cart.x;
    double y = cart.y;
    double z = cart.z;
    double r2 = x*x + y*y + z*z;
    double r = std::sqrt(r2);
    double rho2 = x*x + y*y;
    double rho = std::sqrt(rho2);

    if (r < 1e-15 || rho < 1e-15) {
        // At origin or on z-axis: return identity-ish
        Jinv[0][0] = 1.0;
        Jinv[1][1] = 1.0;
        Jinv[2][2] = 1.0;
        Jinv[3][3] = 1.0;
        return Jinv;
    }

    // ∂t/∂t' = 1
    Jinv[0][0] = 1.0;

    // ∂r/∂x = x/r, ∂r/∂y = y/r, ∂r/∂z = z/r
    Jinv[1][1] = x / r;
    Jinv[1][2] = y / r;
    Jinv[1][3] = z / r;

    // ∂θ/∂x = xz/(r²ρ), ∂θ/∂y = yz/(r²ρ), ∂θ/∂z = -ρ/r²
    Jinv[2][1] = x * z / (r2 * rho);
    Jinv[2][2] = y * z / (r2 * rho);
    Jinv[2][3] = -rho / r2;

    // ∂φ/∂x = -y/ρ², ∂φ/∂y = x/ρ², ∂φ/∂z = 0
    Jinv[3][1] = -y / rho2;
    Jinv[3][2] = x / rho2;
    Jinv[3][3] = 0.0;

    return Jinv;
}

// =============================================================================
// Vector Transformation
// =============================================================================

/// @brief Transform contravariant vector from BL to Cartesian
/// @param v_BL Vector in Boyer-Lindquist coordinates
/// @param bl Position in Boyer-Lindquist coordinates
/// @return Vector in Cartesian coordinates
inline Vec4Cart transformVectorBLToCart(const Vec4BL& v_BL, const Vec4BL& pos_BL) {
    Jacobian4x4 J = JacobianBLToCartesian(pos_BL);
    Vec4Cart v_Cart;

    for (int mu = 0; mu < 4; ++mu) {
        double sum = 0;
        for (int nu = 0; nu < 4; ++nu) {
            sum += J[mu][nu] * v_BL[nu];
        }
        v_Cart[mu] = sum;
    }

    return v_Cart;
}

/// @brief Transform contravariant vector from Cartesian to BL
inline Vec4BL transformVectorCartToBL(const Vec4Cart& v_Cart, const Vec4Cart& pos_Cart) {
    Jacobian4x4 Jinv = JacobianCartesianToBL(pos_Cart);
    Vec4BL v_BL;

    for (int mu = 0; mu < 4; ++mu) {
        double sum = 0;
        for (int nu = 0; nu < 4; ++nu) {
            sum += Jinv[mu][nu] * v_Cart[nu];
        }
        v_BL[mu] = sum;
    }

    return v_BL;
}

// =============================================================================
// Validation
// =============================================================================

/// @brief Verify Jacobian determinant (should be r² sin(θ) for BL→Cart)
inline double JacobianDeterminant(const Jacobian4x4& J) {
    // For 4x4, compute det of spatial 3x3 block (J00 = 1)
    double det3x3 =
        J[1][1] * (J[2][2]*J[3][3] - J[2][3]*J[3][2]) -
        J[1][2] * (J[2][1]*J[3][3] - J[2][3]*J[3][1]) +
        J[1][3] * (J[2][1]*J[3][2] - J[2][2]*J[3][1]);
    return det3x3;  // Full det = J[0][0] * det3x3 = det3x3
}

/// @brief Round-trip test: BL → Cart → BL
/// @return Maximum coordinate deviation
inline double validateRoundTrip(const Vec4BL& original, double a = 0) {
    Vec4Cart cart = (std::abs(a) < 1e-12)
        ? BLToCartesian(original)
        : BLToKerrSchildCart(original, a);

    Vec4BL recovered = (std::abs(a) < 1e-12)
        ? CartesianToBL(cart)
        : KerrSchildCartToBL(cart, a);

    double max_dev = 0;
    max_dev = std::max(max_dev, std::abs(original.t - recovered.t));
    max_dev = std::max(max_dev, std::abs(original.r - recovered.r));
    max_dev = std::max(max_dev, std::abs(original.theta - recovered.theta));

    // Handle phi wrapping
    double dphi = original.phi - recovered.phi;
    while (dphi > M_PI) dphi -= 2*M_PI;
    while (dphi < -M_PI) dphi += 2*M_PI;
    max_dev = std::max(max_dev, std::abs(dphi));

    return max_dev;
}

} // namespace Coordinates
} // namespace Sirius

#endif // PHCT002A_H
