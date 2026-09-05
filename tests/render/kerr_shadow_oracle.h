#pragma once

#include "sirius/core/coordinates.h"
#include "sirius/core/metrics/kerr_schild_family.h"
#include "sirius/core/observer_frame.h"

#include <array>
#include <cmath>
#include <optional>

namespace sirius::render::test {

struct KerrShadowScreenPoint {
    double alpha;
    double beta;
};

struct KerrStationaryObserver {
    core::coordinates::Vec4Cart position;
    core::Metric4d metric;
    core::relativity::ObserverFrame frame;
    // Velocity of the stationary observer relative to the Eulerian frame,
    // expressed in the renderer's screen-forward/up/right convention.
    std::array<double, 3> screen_beta;
};

// Bardeen celestial coordinates for an unstable spherical Kerr photon orbit
// with M=1. The result is the asymptotic upper shadow curve.
[[nodiscard]] inline std::optional<KerrShadowScreenPoint> BardeenShadowPoint(double photon_radius,
                                                                             double spin,
                                                                             double inclination) {
    const double radius = photon_radius;
    const double spin_squared = spin * spin;
    const double xi = (radius * radius * (radius - 3.0) + spin_squared * (radius + 1.0)) /
                      (spin * (1.0 - radius));
    const double eta = radius * radius * radius *
                       (4.0 * spin_squared - radius * (radius - 3.0) * (radius - 3.0)) /
                       (spin_squared * (1.0 - radius) * (1.0 - radius));
    const double sin_inclination = std::sin(inclination);
    const double cos_inclination = std::cos(inclination);
    const double alpha = -xi / sin_inclination;
    const double beta_squared =
        eta + spin_squared * cos_inclination * cos_inclination -
        xi * xi * cos_inclination * cos_inclination / (sin_inclination * sin_inclination);
    if (beta_squared < 0.0) return std::nullopt;
    return KerrShadowScreenPoint{alpha, std::sqrt(beta_squared)};
}

[[nodiscard]] inline std::optional<KerrStationaryObserver> MakeStationaryKerrObserver(
    core::KerrSchildFamily& family, double spin, double radius, double inclination) {
    const core::coordinates::Vec4Bl observer_bl{0.0, radius, inclination, 0.0};
    const auto observer_cart = core::coordinates::BlToKerrSchildCart(observer_bl, spin);
    core::Vec4 event;
    event(0) = observer_cart.t;
    event(1) = observer_cart.x;
    event(2) = observer_cart.y;
    event(3) = observer_cart.z;

    core::Metric4d metric;
    core::Tensor<core::Dual<double>, 4, 4, 4> derivatives;
    family.Evaluate(event, metric, derivatives);
    core::Metric4d inverse;
    if (!family.InverseMetric(event, inverse)) return std::nullopt;

    const double sin_theta = std::sin(inclination);
    const double cos_theta = std::cos(inclination);
    const double cartesian_phi = std::atan2(observer_cart.y, observer_cart.x);
    const double sin_phi = std::sin(cartesian_phi);
    const double cos_phi = std::cos(cartesian_phi);
    std::array<core::Vec4, 3> seeds;
    seeds[0](1) = sin_theta * cos_phi;
    seeds[0](2) = sin_theta * sin_phi;
    seeds[0](3) = cos_theta;
    seeds[1](1) = cos_theta * cos_phi;
    seeds[1](2) = cos_theta * sin_phi;
    seeds[1](3) = -sin_theta;
    seeds[2](1) = -sin_phi;
    seeds[2](2) = cos_phi;
    const auto eulerian = core::relativity::EulerianObserverFrame(metric, inverse, seeds);
    if (!eulerian.has_value()) return std::nullopt;

    core::Vec4 stationary;
    stationary(0) = 1.0 / std::sqrt(-metric(0, 0).real);
    const double gamma = -core::TensorOps::InnerProduct(stationary, eulerian->time, metric);
    if (!(gamma > 0.0)) return std::nullopt;
    std::array<double, 3> radial_theta_phi_beta{};
    for (std::size_t axis = 0; axis < radial_theta_phi_beta.size(); ++axis) {
        radial_theta_phi_beta[axis] =
            core::TensorOps::InnerProduct(stationary, eulerian->spatial[axis], metric) / gamma;
    }
    const auto stationary_frame =
        core::relativity::BoostObserverFrame(*eulerian, radial_theta_phi_beta);
    if (!stationary_frame.has_value()) return std::nullopt;

    return KerrStationaryObserver{
        observer_cart,
        metric,
        *stationary_frame,
        {-radial_theta_phi_beta[0], -radial_theta_phi_beta[1], radial_theta_phi_beta[2]}};
}

// Map an asymptotic Bardeen point to the exact sky of a finite stationary
// observer. This constructs the critical photon from its independent Kerr
// constants (E=1, L_z=xi, Q=eta), transforms it through the repository's full
// BL/Kerr-Schild chart authority, and resolves it in the observer tetrad.
[[nodiscard]] inline std::optional<KerrShadowScreenPoint> ProjectBardeenAtFiniteObserver(
    const KerrShadowScreenPoint& asymptotic, const KerrStationaryObserver& observer, double mass,
    double spin, double radius, double inclination) {
    const double sin_theta = std::sin(inclination);
    const double cos_theta = std::cos(inclination);
    const double sin_theta_squared = sin_theta * sin_theta;
    const double xi = -asymptotic.alpha * sin_theta;
    const double eta = asymptotic.beta * asymptotic.beta - spin * spin * cos_theta * cos_theta +
                       xi * xi * cos_theta * cos_theta / sin_theta_squared;
    const double sigma = radius * radius + spin * spin * cos_theta * cos_theta;
    const double delta = radius * radius - 2.0 * mass * radius + spin * spin;
    const double p = radius * radius + spin * spin - spin * xi;
    const double radial_potential = p * p - delta * ((xi - spin) * (xi - spin) + eta);
    const double polar_potential = eta + spin * spin * cos_theta * cos_theta -
                                   xi * xi * cos_theta * cos_theta / sin_theta_squared;
    if (!(radial_potential >= 0.0) || !(polar_potential >= 0.0)) return std::nullopt;

    core::coordinates::Vec4Bl photon_bl;
    photon_bl.t =
        (-spin * (spin * sin_theta_squared - xi) + (radius * radius + spin * spin) * p / delta) /
        sigma;
    photon_bl.r = std::sqrt(radial_potential) / sigma;
    photon_bl.theta = std::sqrt(polar_potential) / sigma;
    photon_bl.phi = (-(spin - xi / sin_theta_squared) + spin * p / delta) / sigma;

    // Invert the full Kerr-Schild-Cartesian -> BL tangent map at the observer.
    // The spatial oblate Jacobian by itself omits the ingoing Kerr-Schild time
    // and azimuth twists and is not a finite-radius conserved-quantity oracle.
    core::Metric4d cart_to_bl;
    for (int column = 0; column < 4; ++column) {
        core::coordinates::Vec4Cart basis;
        basis[column] = 1.0;
        const auto transformed = core::coordinates::TransformVectorKerrSchildCartToBl(
            basis, observer.position, mass, spin);
        for (int row = 0; row < 4; ++row) cart_to_bl(row, column) = transformed[row];
    }
    const core::Metric4d bl_to_cart = core::TensorOps::Inverse(cart_to_bl);
    core::Vec4 photon;
    for (int row = 0; row < 4; ++row) {
        for (int column = 0; column < 4; ++column) {
            photon(row) += bl_to_cart(row, column).real * photon_bl[column];
        }
    }

    const double radial =
        core::TensorOps::InnerProduct(photon, observer.frame.spatial[0], observer.metric);
    const double polar =
        core::TensorOps::InnerProduct(photon, observer.frame.spatial[1], observer.metric);
    const double azimuthal =
        core::TensorOps::InnerProduct(photon, observer.frame.spatial[2], observer.metric);
    if (!(radial > 0.0) || !std::isfinite(polar) || !std::isfinite(azimuthal)) {
        return std::nullopt;
    }
    return KerrShadowScreenPoint{-radius * azimuthal / radial, radius * polar / radial};
}

}  // namespace sirius::render::test
