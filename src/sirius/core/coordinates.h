#pragma once

// Explicit, validated transforms between Boyer-Lindquist (t, r, theta, phi),
// Kerr-Schild Cartesian (t, x, y, z), and the observer frame, with the Jacobian
// matrices for tensor transformation. Ported from PHCT002A.h.
//
// The full Kerr-Schild embedding is
//   x = (r cos(phi) - a sin(phi)) sin(theta)
//   y = (r sin(phi) + a cos(phi)) sin(theta)
//   z = r cos(theta),
// which degenerates to the ordinary spherical map at a = 0 or large r.

#include "sirius/base/contracts.h"
#include "sirius/core/tensor.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numbers>
#include <optional>

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
            case 0:
                return t;
            case 1:
                return r;
            case 2:
                return theta;
            case 3:
                return phi;
            default:
                SIRIUS_PRE(i >= 0 && i < 4);
                return phi;
        }
    }
    double operator[](int i) const {
        switch (i) {
            case 0:
                return t;
            case 1:
                return r;
            case 2:
                return theta;
            case 3:
                return phi;
            default:
                SIRIUS_PRE(i >= 0 && i < 4);
                return phi;
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
            case 0:
                return t;
            case 1:
                return x;
            case 2:
                return y;
            case 3:
                return z;
            default:
                SIRIUS_PRE(i >= 0 && i < 4);
                return z;
        }
    }
    double operator[](int i) const {
        switch (i) {
            case 0:
                return t;
            case 1:
                return x;
            case 2:
                return y;
            case 3:
                return z;
            default:
                SIRIUS_PRE(i >= 0 && i < 4);
                return z;
        }
    }

    double Radius() const { return std::sqrt(x * x + y * y + z * z); }
};

// The non-negative oblate radial root and its Cartesian gradient.  The root is
// defined by
//   r^4 - (x^2 + y^2 + z^2 - a^2) r^2 - a^2 z^2 = 0.
// Its differential is absent at r=0, where the positive-root sheet terminates
// on the Kerr disk (or at the spherical origin when a=0).  Callers that need
// only the coordinate radius may still observe the exact zero root.
struct KerrSchildRadiusDifferential {
    double radius = 0.0;
    double dx = 0.0;
    double dy = 0.0;
    double dz = 0.0;
};

namespace detail {

struct KerrSchildRadiusSolution {
    double radius = 0.0;
    double scaled_radius = 0.0;
    double scaled_discriminant_root = 0.0;
    double scaled_x = 0.0;
    double scaled_y = 0.0;
    double scaled_z = 0.0;
    double scaled_a = 0.0;
};

[[nodiscard]] inline std::optional<KerrSchildRadiusSolution> TrySolveKerrSchildRadius(
    const Vec4Cart& cart, double a) {
    if (!std::isfinite(cart.t) || !std::isfinite(cart.x) || !std::isfinite(cart.y) ||
        !std::isfinite(cart.z) || !std::isfinite(a)) {
        return std::nullopt;
    }

    const double scale =
        std::max({std::abs(cart.x), std::abs(cart.y), std::abs(cart.z), std::abs(a)});
    if (scale == 0.0) return KerrSchildRadiusSolution{};

    const double x = cart.x / scale;
    const double y = cart.y / scale;
    const double z = cart.z / scale;
    const double spin = a / scale;
    const double spin2 = spin * spin;
    const double reduced = x * x + y * y + z * z - spin2;
    const double discriminant_root = std::hypot(reduced, 2.0 * spin * z);

    // The ordinary quadratic expression loses the entire positive root when
    // reduced < 0 and |z| is small.  Compute sqrt(r^2) directly from the
    // conjugate expression in that branch, avoiding both cancellation and the
    // underflow-prone square of a*z.
    double scaled_radius = 0.0;
    if (reduced >= 0.0) {
        scaled_radius = std::sqrt(0.5 * (reduced + discriminant_root));
    } else if (spin != 0.0 && z != 0.0) {
        scaled_radius =
            std::numbers::sqrt2 * std::abs(spin * z) / std::sqrt(discriminant_root - reduced);
    }

    const double radius = scale * scaled_radius;
    if (!std::isfinite(radius)) return std::nullopt;
    return KerrSchildRadiusSolution{radius, scaled_radius, discriminant_root, x, y, z, spin};
}

}  // namespace detail

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

    // x + i y = (r + i a) exp(i phi) sin(theta).
    cart.x = (bl.r * cos_phi - a * sin_phi) * sin_theta;
    cart.y = (bl.r * sin_phi + a * cos_phi) * sin_theta;
    cart.z = bl.r * cos_theta;

    return cart;
}

// Boyer-Lindquist/Kerr spheroidal radius represented by a Kerr-Schild
// Cartesian point. This is not the Euclidean cylindrical or spherical radius
// when a != 0. Disk orbit, ISCO, Page-Thorne, and horizon calculations must use
// this authority.
inline double KerrSchildRadius(const Vec4Cart& cart, double a) {
    const auto solution = detail::TrySolveKerrSchildRadius(cart, a);
    SIRIUS_PRE(solution.has_value());
    return solution ? solution->radius : std::numeric_limits<double>::quiet_NaN();
}

// Fallible radius-and-gradient authority for metric derivatives and inverse
// chart Jacobians.  It accepts each numerically representable finite
// differentiable point, including non-zero spin below the former epsilon
// branch, and declines the exact r=0 sheet boundary instead of manufacturing
// an epsilon radius.
[[nodiscard]] inline std::optional<KerrSchildRadiusDifferential> TryKerrSchildRadiusDifferential(
    const Vec4Cart& cart, double a) {
    const auto solution = detail::TrySolveKerrSchildRadius(cart, a);
    if (!solution || !(solution->radius > 0.0) || !(solution->scaled_radius > 0.0) ||
        !(solution->scaled_discriminant_root > 0.0)) {
        return std::nullopt;
    }

    const double scaled_r = solution->scaled_radius;
    const double scaled_d = solution->scaled_discriminant_root;
    const double scaled_r2 = scaled_r * scaled_r;
    const double scaled_a2 = solution->scaled_a * solution->scaled_a;
    KerrSchildRadiusDifferential differential;
    differential.radius = solution->radius;
    differential.dx = solution->scaled_x * scaled_r / scaled_d;
    differential.dy = solution->scaled_y * scaled_r / scaled_d;
    differential.dz = solution->scaled_z * (scaled_r2 + scaled_a2) / (scaled_r * scaled_d);
    if (!std::isfinite(differential.dx) || !std::isfinite(differential.dy) ||
        !std::isfinite(differential.dz)) {
        return std::nullopt;
    }
    return differential;
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
    if (bl.phi < 0) bl.phi += 2.0 * std::numbers::pi;

    return bl;
}

// Kerr-Schild Cartesian to Boyer-Lindquist (full Kerr), solving
// r^4 - (R^2 - a^2) r^2 - a^2 z^2 = 0.
inline Vec4Bl KerrSchildCartToBl(const Vec4Cart& cart, double a) {
    Vec4Bl bl;
    bl.t = cart.t;

    const double x = cart.x;
    const double y = cart.y;
    const double z = cart.z;
    bl.r = KerrSchildRadius(cart, a);

    const bool represented = std::isfinite(bl.r) && bl.r > 0.0;
    SIRIUS_PRE(represented);
    if (!represented) {
        const double nan = std::numeric_limits<double>::quiet_NaN();
        return Vec4Bl{nan, nan, nan, nan};
    }

    const double cylindrical_radius = std::hypot(x, y);
    SIRIUS_PRE(cylindrical_radius > 0.0);
    if (!(cylindrical_radius > 0.0)) {
        const double nan = std::numeric_limits<double>::quiet_NaN();
        return Vec4Bl{nan, nan, nan, nan};
    }
    bl.theta = std::acos(std::clamp(z / bl.r, -1.0, 1.0));

    // phi_BL = atan2(y, x) - atan2(a, r) for the full Kerr transform.
    bl.phi = std::atan2(y, x);
    if (a != 0.0) {
        bl.phi -= std::atan2(a, bl.r);
    }
    while (bl.phi < 0) bl.phi += 2.0 * std::numbers::pi;
    while (bl.phi >= 2.0 * std::numbers::pi) bl.phi -= 2.0 * std::numbers::pi;

    return bl;
}

// Transform a contravariant vector from Cartesian Kerr-Schild to
// Boyer-Lindquist components. Besides the inverse spatial oblate Jacobian, the
// time and azimuth components include the Kerr-Schild twist
//   dt_BL = dt_KS - 2Mr/Delta dr,
//   dphi_BL = dphi_KS - a/Delta dr.
// This is the chart authority used by live-path Carter and Walker-Penrose
// monitors; copying a spherical inverse Jacobian here loses both invariants.
inline Vec4Bl TransformVectorKerrSchildCartToBl(const Vec4Cart& vector, const Vec4Cart& position,
                                                double mass, double spin) {
    const double x = position.x;
    const double y = position.y;
    const double z = position.z;
    const bool finite = std::isfinite(vector.t) && std::isfinite(vector.x) &&
                        std::isfinite(vector.y) && std::isfinite(vector.z) && std::isfinite(mass) &&
                        std::isfinite(spin);
    const auto radial_geometry =
        finite ? TryKerrSchildRadiusDifferential(position, spin) : std::nullopt;
    SIRIUS_PRE(radial_geometry.has_value());
    if (!radial_geometry) {
        const double nan = std::numeric_limits<double>::quiet_NaN();
        return Vec4Bl{nan, nan, nan, nan};
    }
    const double r = radial_geometry->radius;
    const double dr_dx = radial_geometry->dx;
    const double dr_dy = radial_geometry->dy;
    const double dr_dz = radial_geometry->dz;

    const double scale =
        std::max({std::abs(x), std::abs(y), std::abs(z), std::abs(r), std::abs(spin)});
    SIRIUS_ASSERT(scale > 0.0);
    const double scaled_x = x / scale;
    const double scaled_y = y / scale;
    const double scaled_z = z / scale;
    const double scaled_r = r / scale;
    const double scaled_a = spin / scale;
    const double scaled_r2 = scaled_r * scaled_r;
    const double sin_theta = std::hypot(scaled_x, scaled_y) / std::hypot(scaled_r, scaled_a);
    SIRIUS_PRE(std::isfinite(sin_theta) && sin_theta > 0.0);
    if (!(std::isfinite(sin_theta) && sin_theta > 0.0)) {
        const double nan = std::numeric_limits<double>::quiet_NaN();
        return Vec4Bl{nan, nan, nan, nan};
    }
    const double inverse_theta_denominator = 1.0 / (scale * scaled_r2 * sin_theta);
    const double dtheta_dx = scaled_z * dr_dx * inverse_theta_denominator;
    const double dtheta_dy = scaled_z * dr_dy * inverse_theta_denominator;
    const double dtheta_dz = -(scaled_r - scaled_z * dr_dz) * inverse_theta_denominator;

    const double numerator = scaled_r * scaled_y - scaled_a * scaled_x;
    const double denominator = scaled_r * scaled_x + scaled_a * scaled_y;
    const double phi_denominator = numerator * numerator + denominator * denominator;
    SIRIUS_PRE(std::isfinite(phi_denominator) && phi_denominator > 0.0);
    if (!(std::isfinite(phi_denominator) && phi_denominator > 0.0)) {
        const double nan = std::numeric_limits<double>::quiet_NaN();
        return Vec4Bl{nan, nan, nan, nan};
    }
    const double numerator_dx = scaled_y * dr_dx - scaled_a;
    const double numerator_dy = scaled_y * dr_dy + scaled_r;
    const double numerator_dz = scaled_y * dr_dz;
    const double denominator_dx = scaled_x * dr_dx + scaled_r;
    const double denominator_dy = scaled_x * dr_dy + scaled_a;
    const double denominator_dz = scaled_x * dr_dz;
    const double dphi_dx =
        (denominator * numerator_dx - numerator * denominator_dx) / (scale * phi_denominator);
    const double dphi_dy =
        (denominator * numerator_dy - numerator * denominator_dy) / (scale * phi_denominator);
    const double dphi_dz =
        (denominator * numerator_dz - numerator * denominator_dz) / (scale * phi_denominator);

    const double radial = dr_dx * vector.x + dr_dy * vector.y + dr_dz * vector.z;
    const double polar = dtheta_dx * vector.x + dtheta_dy * vector.y + dtheta_dz * vector.z;
    const double azimuth_ks = dphi_dx * vector.x + dphi_dy * vector.y + dphi_dz * vector.z;
    const double delta_scale = std::max({std::abs(r), std::abs(mass), std::abs(spin)});
    const double scaled_mass = mass / delta_scale;
    const double delta_radius = r / delta_scale;
    const double delta_spin = spin / delta_scale;
    const double scaled_delta =
        delta_radius * delta_radius - 2.0 * scaled_mass * delta_radius + delta_spin * delta_spin;
    SIRIUS_PRE(std::isfinite(scaled_delta) && scaled_delta != 0.0);
    if (!std::isfinite(scaled_delta) || scaled_delta == 0.0) {
        const double nan = std::numeric_limits<double>::quiet_NaN();
        return Vec4Bl{nan, nan, nan, nan};
    }

    const double time_twist = 2.0 * scaled_mass * delta_radius / scaled_delta;
    const double azimuth_twist = delta_spin / (delta_scale * scaled_delta);

    return Vec4Bl{vector.t - time_twist * radial, radial, polar,
                  azimuth_ks - azimuth_twist * radial};
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

// Exact spatial Jacobian of BlToKerrSchildCart. This includes the spin-coupled
// azimuth terms; using the spherical a=0 Jacobian changes a Kerr camera ray's
// conserved angular momentum before integration begins.
inline Jacobian4x4 JacobianBlToKerrSchildCart(const Vec4Bl& bl, double a) {
    Jacobian4x4 J = {};
    const double sin_th = std::sin(bl.theta);
    const double cos_th = std::cos(bl.theta);
    const double sin_ph = std::sin(bl.phi);
    const double cos_ph = std::cos(bl.phi);
    const double r = bl.r;

    J[0][0] = 1.0;
    J[1][1] = cos_ph * sin_th;
    J[1][2] = (r * cos_ph - a * sin_ph) * cos_th;
    J[1][3] = (-r * sin_ph - a * cos_ph) * sin_th;
    J[2][1] = sin_ph * sin_th;
    J[2][2] = (r * sin_ph + a * cos_ph) * cos_th;
    J[2][3] = (r * cos_ph - a * sin_ph) * sin_th;
    J[3][1] = cos_th;
    J[3][2] = -r * sin_th;
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
    double det3x3 = J[1][1] * (J[2][2] * J[3][3] - J[2][3] * J[3][2]) -
                    J[1][2] * (J[2][1] * J[3][3] - J[2][3] * J[3][1]) +
                    J[1][3] * (J[2][1] * J[3][2] - J[2][2] * J[3][1]);
    return det3x3;
}

// Round-trip BL -> Cart -> BL, returning the maximum coordinate deviation.
inline double ValidateRoundTrip(const Vec4Bl& original, double a = 0) {
    Vec4Cart cart = (a == 0.0) ? BlToCartesian(original) : BlToKerrSchildCart(original, a);

    Vec4Bl recovered = (a == 0.0) ? CartesianToBl(cart) : KerrSchildCartToBl(cart, a);

    double max_dev = 0;
    max_dev = std::max(max_dev, std::abs(original.t - recovered.t));
    max_dev = std::max(max_dev, std::abs(original.r - recovered.r));
    max_dev = std::max(max_dev, std::abs(original.theta - recovered.theta));

    double dphi = original.phi - recovered.phi;
    while (dphi > std::numbers::pi) dphi -= 2 * std::numbers::pi;
    while (dphi < -std::numbers::pi) dphi += 2 * std::numbers::pi;
    max_dev = std::max(max_dev, std::abs(dphi));

    return max_dev;
}

}  // namespace sirius::core::coordinates
