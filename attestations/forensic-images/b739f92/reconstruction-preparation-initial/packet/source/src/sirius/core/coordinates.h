#pragma once

// Explicit, validated transforms between Boyer-Lindquist (t, r, theta, phi),
// Kerr-Schild Cartesian (t, x, y, z), and the observer frame, with the Jacobian
// matrices for tensor transformation. Ported from PHCT002A.h.
//
// The full Kerr embedding is
//   x = sqrt(r^2 + a^2) sin(theta) cos(phi)
//   y = sqrt(r^2 + a^2) sin(theta) sin(phi)
//   z = r cos(theta),
// which degenerates to the ordinary spherical map at a = 0 or large r.

#include "sirius/core/tensor.h"

#include <algorithm>
#include <array>
#include <cmath>

namespace sirius::core::coordinates {

enum class CoordinateSystem {
    BoyerLindquist,  // (t, r, theta, phi), spherical-like.
    KerrSchildCart,  // (t, x, y, z), Cartesian.
    Spherical,       // (t, r, theta, phi), standard spherical (a = 0).
    Cartesian        // (t, x, y, z), flat-space Cartesian.
};

// Boyer-Lindquist coordinates (t, r, theta, phi).
struct Vec4Bl {
    double t = 0;
    double r = 0;
    double theta = 0;  // Polar angle [0, pi].
    double phi = 0;    // Azimuthal angle [0, 2 pi).

    Vec4Bl() = default;
    Vec4Bl(double t_, double r_, double th_, double ph_) : t(t_), r(r_), theta(th_), phi(ph_) {}

    double& operator[](int i) {
        switch (i) {
            case 0: return t;
            case 1: return r;
            case 2: return theta;
            default: return phi;
        }
    }
    double operator[](int i) const {
        switch (i) {
            case 0: return t;
            case 1: return r;
            case 2: return theta;
            default: return phi;
        }
    }
};

// Kerr-Schild Cartesian coordinates (t, x, y, z).
struct Vec4Cart {
    double t = 0;
    double x = 0;
    double y = 0;
    double z = 0;

    Vec4Cart() = default;
    Vec4Cart(double t_, double x_, double y_, double z_) : t(t_), x(x_), y(y_), z(z_) {}

    double& operator[](int i) {
        switch (i) {
            case 0: return t;
            case 1: return x;
            case 2: return y;
            default: return z;
        }
    }
    double operator[](int i) const {
        switch (i) {
            case 0: return t;
            case 1: return x;
            case 2: return y;
            default: return z;
        }
    }

    double Radius() const { return std::sqrt(x * x + y * y + z * z); }
};

// 4x4 Jacobian J[mu][nu] = d x'^mu / d x^nu.
using Jacobian4x4 = std::array<std::array<double, 4>, 4>;

// Boyer-Lindquist to Cartesian, simplified for a = 0 or large r.
inline Vec4Cart BlToCartesian(const Vec4Bl& bl) {
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

// Boyer-Lindquist to Kerr-Schild Cartesian (full Kerr oblate spheroid).
inline Vec4Cart BlToKerrSchildCart(const Vec4Bl& bl, double a) {
    Vec4Cart cart;
    cart.t = bl.t;

    double sin_theta = std::sin(bl.theta);
    double cos_theta = std::cos(bl.theta);
    double sin_phi = std::sin(bl.phi);
    double cos_phi = std::cos(bl.phi);

    double rho = std::sqrt(bl.r * bl.r + a * a);

    cart.x = rho * sin_theta * cos_phi;
    cart.y = rho * sin_theta * sin_phi;
    cart.z = bl.r * cos_theta;

    return cart;
}

// Cartesian to Boyer-Lindquist, simplified for a = 0.
inline Vec4Bl CartesianToBl(const Vec4Cart& cart) {
    Vec4Bl bl;
    bl.t = cart.t;

    bl.r = std::sqrt(cart.x * cart.x + cart.y * cart.y + cart.z * cart.z);
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

// Kerr-Schild Cartesian to Boyer-Lindquist (full Kerr), solving
// r^4 - (R^2 - a^2) r^2 - a^2 z^2 = 0.
inline Vec4Bl KerrSchildCartToBl(const Vec4Cart& cart, double a) {
    Vec4Bl bl;
    bl.t = cart.t;

    double x = cart.x;
    double y = cart.y;
    double z = cart.z;
    double a2 = a * a;

    double R2 = x * x + y * y + z * z;

    if (std::abs(a) < 1e-12) {
        bl.r = std::sqrt(R2);
    } else {
        double Rm2 = R2 - a2;
        double disc = Rm2 * Rm2 + 4 * a2 * z * z;
        double r2 = (Rm2 + std::sqrt(disc)) / 2.0;
        bl.r = std::sqrt(std::max(r2, 1e-20));
    }

    if (bl.r < 1e-15) {
        bl.theta = M_PI / 2.0;
        bl.phi = 0;
        return bl;
    }

    bl.theta = std::acos(std::clamp(z / bl.r, -1.0, 1.0));

    // phi_BL = atan2(y, x) - atan2(a, r) for the full Kerr transform.
    bl.phi = std::atan2(y, x);
    if (std::abs(a) > 1e-12) {
        bl.phi -= std::atan2(a, bl.r);
    }
    while (bl.phi < 0) bl.phi += 2.0 * M_PI;
    while (bl.phi >= 2.0 * M_PI) bl.phi -= 2.0 * M_PI;

    return bl;
}

// Jacobian J[mu][nu] = d x^mu_Cart / d x^nu_BL; transforms contravariant
// vectors via V'^mu = J[mu][nu] V^nu.
inline Jacobian4x4 JacobianBlToCartesian(const Vec4Bl& bl) {
    Jacobian4x4 J = {};

    double sin_th = std::sin(bl.theta);
    double cos_th = std::cos(bl.theta);
    double sin_ph = std::sin(bl.phi);
    double cos_ph = std::cos(bl.phi);
    double r = bl.r;

    J[0][0] = 1.0;

    J[1][1] = sin_th * cos_ph;
    J[1][2] = r * cos_th * cos_ph;
    J[1][3] = -r * sin_th * sin_ph;

    J[2][1] = sin_th * sin_ph;
    J[2][2] = r * cos_th * sin_ph;
    J[2][3] = r * sin_th * cos_ph;

    J[3][1] = cos_th;
    J[3][2] = -r * sin_th;
    J[3][3] = 0.0;

    return J;
}

// Inverse Jacobian J^{-1}[mu][nu] = d x^mu_BL / d x^nu_Cart.
inline Jacobian4x4 JacobianCartesianToBl(const Vec4Cart& cart) {
    Jacobian4x4 Jinv = {};

    double x = cart.x;
    double y = cart.y;
    double z = cart.z;
    double r2 = x * x + y * y + z * z;
    double r = std::sqrt(r2);
    double rho2 = x * x + y * y;
    double rho = std::sqrt(rho2);

    if (r < 1e-15 || rho < 1e-15) {
        // At the origin or on the z-axis: return an identity-like Jacobian.
        Jinv[0][0] = 1.0;
        Jinv[1][1] = 1.0;
        Jinv[2][2] = 1.0;
        Jinv[3][3] = 1.0;
        return Jinv;
    }

    Jinv[0][0] = 1.0;

    Jinv[1][1] = x / r;
    Jinv[1][2] = y / r;
    Jinv[1][3] = z / r;

    Jinv[2][1] = x * z / (r2 * rho);
    Jinv[2][2] = y * z / (r2 * rho);
    Jinv[2][3] = -rho / r2;

    Jinv[3][1] = -y / rho2;
    Jinv[3][2] = x / rho2;
    Jinv[3][3] = 0.0;

    return Jinv;
}

// Transform a contravariant vector from BL to Cartesian at a BL position.
inline Vec4Cart TransformVectorBlToCart(const Vec4Bl& v_bl, const Vec4Bl& pos_bl) {
    Jacobian4x4 J = JacobianBlToCartesian(pos_bl);
    Vec4Cart v_cart;

    for (int mu = 0; mu < 4; ++mu) {
        double sum = 0;
        for (int nu = 0; nu < 4; ++nu) {
            sum += J[mu][nu] * v_bl[nu];
        }
        v_cart[mu] = sum;
    }

    return v_cart;
}

// Transform a contravariant vector from Cartesian to BL at a Cartesian position.
inline Vec4Bl TransformVectorCartToBl(const Vec4Cart& v_cart, const Vec4Cart& pos_cart) {
    Jacobian4x4 Jinv = JacobianCartesianToBl(pos_cart);
    Vec4Bl v_bl;

    for (int mu = 0; mu < 4; ++mu) {
        double sum = 0;
        for (int nu = 0; nu < 4; ++nu) {
            sum += Jinv[mu][nu] * v_cart[nu];
        }
        v_bl[mu] = sum;
    }

    return v_bl;
}

// Jacobian determinant; equals r^2 sin(theta) for the BL->Cart map.
inline double JacobianDeterminant(const Jacobian4x4& J) {
    // J[0][0] = 1, so the full determinant reduces to the spatial 3x3 block.
    double det3x3 =
        J[1][1] * (J[2][2] * J[3][3] - J[2][3] * J[3][2]) -
        J[1][2] * (J[2][1] * J[3][3] - J[2][3] * J[3][1]) +
        J[1][3] * (J[2][1] * J[3][2] - J[2][2] * J[3][1]);
    return det3x3;
}

// Round-trip BL -> Cart -> BL, returning the maximum coordinate deviation.
inline double ValidateRoundTrip(const Vec4Bl& original, double a = 0) {
    Vec4Cart cart = (std::abs(a) < 1e-12)
        ? BlToCartesian(original)
        : BlToKerrSchildCart(original, a);

    Vec4Bl recovered = (std::abs(a) < 1e-12)
        ? CartesianToBl(cart)
        : KerrSchildCartToBl(cart, a);

    double max_dev = 0;
    max_dev = std::max(max_dev, std::abs(original.t - recovered.t));
    max_dev = std::max(max_dev, std::abs(original.r - recovered.r));
    max_dev = std::max(max_dev, std::abs(original.theta - recovered.theta));

    double dphi = original.phi - recovered.phi;
    while (dphi > M_PI) dphi -= 2 * M_PI;
    while (dphi < -M_PI) dphi += 2 * M_PI;
    max_dev = std::max(max_dev, std::abs(dphi));

    return max_dev;
}

}  // namespace sirius::core::coordinates
