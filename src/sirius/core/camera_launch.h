#pragma once

#include "sirius/core/camera.h"
#include "sirius/core/coordinates.h"
#include "sirius/core/geodesic_integrator.h"
#include "sirius/core/observer_frame.h"

#include <array>
#include <cmath>
#include <limits>
#include <optional>

namespace sirius::core {

struct CameraLaunch {
    Vec4 position;
    Vec4 tangent;
    relativity::ObserverFrame observer;
    GeodesicVariations variations;
};

// Convert transported film columns back to the existing launch-angular basis
// only at a legacy beam/source-map consumer. This changes units, not the
// physical trajectory or the four-column acceptance decision.
[[nodiscard]] inline std::optional<std::array<GeodesicVariation, 2>> CameraAngularVariations(
    const GeodesicVariations& columns, const std::optional<CameraPhaseSpaceDifferential>& camera) {
    if (!camera) return std::array<GeodesicVariation, 2>{columns[0], columns[1]};
    const auto& p = camera->film_to_angle;
    const double positive = p[0][0] * p[1][1], negative = p[0][1] * p[1][0];
    const double determinant = std::fma(p[0][0], p[1][1], -negative);
    if (!std::isfinite(determinant) ||
        !(std::abs(determinant) >
          64 * std::numeric_limits<double>::epsilon() * (std::abs(positive) + std::abs(negative))))
        return std::nullopt;
    const std::array<std::array<double, 2>, 2> inverse{
        {{p[1][1] / determinant, -p[0][1] / determinant},
         {-p[1][0] / determinant, p[0][0] / determinant}}};
    std::array<GeodesicVariation, 2> result;
    for (int column = 0; column < 2; ++column) {
        for (int film = 0; film < 2; ++film) {
            result[column].displacement += columns[film].displacement * inverse[film][column];
            result[column].derivative += columns[film].derivative * inverse[film][column];
        }
        if (!relativity::IsFinite(result[column].displacement) ||
            !relativity::IsFinite(result[column].derivative))
            return std::nullopt;
    }
    return result;
}

// One launch authority for the smooth camera and its four film/pupil columns.
// Seeds describe the fixed camera orientation at the central observer. A pupil
// offset displaces that event in the boosted central frame, then rebuilds the
// metric-orthonormal frame at the displaced event. Differentiating only the
// direction, or freezing the second frame, misses physical pupil derivatives.
[[nodiscard]] inline std::optional<CameraLaunch> LaunchCameraRay(IMetric& metric,
                                                                 double absolute_spin,
                                                                 const CameraRay& camera) {
    if (!IsRepresentedCameraRay(camera) || !camera.active || !std::isfinite(absolute_spin))
        return std::nullopt;
    const coordinates::Vec4Bl origin(camera.origin(0), camera.origin(1), camera.origin(2),
                                     camera.origin(3));
    const auto cart = coordinates::BlToKerrSchildCart(origin, absolute_spin);
    CameraLaunch result;
    result.position(0) = cart.t;
    result.position(1) = cart.x;
    result.position(2) = cart.y;
    result.position(3) = cart.z;
    if (!relativity::IsFinite(result.position)) return std::nullopt;

    const double theta = camera.origin(2), phi = std::atan2(cart.y, cart.x);
    const double st = std::sin(theta), ct = std::cos(theta);
    const double sp = std::sin(phi), cp = std::cos(phi);
    std::array<Vec4, 3> seeds;
    seeds[0](1) = st * cp;
    seeds[0](2) = st * sp;
    seeds[0](3) = ct;
    seeds[1](1) = ct * cp;
    seeds[1](2) = ct * sp;
    seeds[1](3) = -st;
    seeds[2](1) = -sp;
    seeds[2](2) = cp;
    const std::array<double, 3> beta{-camera.beta_forward, -camera.beta_up, camera.beta_right};
    Metric4d values, inverse;
    Tensor<Dual<double>, 4, 4, 4> derivatives;
    const auto frame_at = [&](const Vec4& position) -> std::optional<relativity::ObserverFrame> {
        metric.Evaluate(position, values, derivatives);
        if (!metric.InverseMetric(position, inverse)) inverse = TensorOps::Inverse(values);
        const auto eulerian = relativity::EulerianObserverFrame(values, inverse, seeds);
        if (!eulerian) return std::nullopt;
        return relativity::BoostObserverFrame(*eulerian, beta);
    };
    const auto central = frame_at(result.position);
    if (!central) return std::nullopt;
    result.position +=
        central->spatial[2] * camera.aperture_right - central->spatial[1] * camera.aperture_up;
    const auto displaced = frame_at(result.position);
    if (!displaced) return std::nullopt;
    result.observer = *displaced;
    const std::array<double, 3> direction{camera.direction(1), camera.direction(2),
                                          camera.direction(3)};
    const auto tangent = relativity::PastDirectedCameraRay(result.observer, direction);
    if (!tangent) return std::nullopt;
    result.tangent = *tangent;
    if (!camera.phase_space) return result;

    const auto connection = TensorOps::Christoffel(values, derivatives);
    for (std::size_t column = 0; column < result.variations.size(); ++column) {
        auto& variation = result.variations[column];
        const double dr = camera.phase_space->pupil_right[column];
        const double du = camera.phase_space->pupil_up[column];
        if (!std::isfinite(dr) || !std::isfinite(du)) return std::nullopt;
        variation.displacement = central->spatial[2] * dr - central->spatial[1] * du;
        Metric4d varied, varied_inverse;
        Tensor<double, 4, 4> delta;
        for (int mu = 0; mu < 4; ++mu)
            for (int nu = 0; nu < 4; ++nu) {
                for (int axis = 0; axis < 4; ++axis)
                    delta(mu, nu) += derivatives(axis, mu, nu).real * variation.displacement(axis);
                varied(mu, nu) = Dual<double>(values(mu, nu).real, delta(mu, nu));
            }
        for (int mu = 0; mu < 4; ++mu)
            for (int nu = 0; nu < 4; ++nu) {
                double d_inverse = 0;
                for (int a = 0; a < 4; ++a)
                    for (int b = 0; b < 4; ++b)
                        d_inverse -= inverse(mu, a).real * delta(a, b) * inverse(b, nu).real;
                varied_inverse(mu, nu) = Dual<double>(inverse(mu, nu).real, d_inverse);
            }
        const auto frame = relativity::detail::BuildEulerianFrame(varied, varied_inverse, seeds);
        if (!frame) return std::nullopt;
        relativity::ObserverFrame frame_derivative;
        for (int mu = 0; mu < 4; ++mu) {
            frame_derivative.time(mu) = frame->time(mu).dual;
            for (int axis = 0; axis < 3; ++axis)
                frame_derivative.spatial[axis](mu) = frame->spatial[axis](mu).dual;
        }
        // At fixed observer velocity a Lorentz boost is linear in every frame
        // vector, so the same boost authority applies to the derivative frame.
        const auto boosted_derivative = relativity::BoostObserverFrame(frame_derivative, beta);
        if (!boosted_derivative) return std::nullopt;
        std::array<Dual<double>, 3> n;
        Dual<double> norm_squared(0);
        for (int axis = 0; axis < 3; ++axis) {
            n[axis] = Dual<double>(direction[axis], camera.phase_space->direction[axis][column]);
            norm_squared += n[axis] * n[axis];
        }
        if (!isfinite(norm_squared) || !(norm_squared.real > 0)) return std::nullopt;
        const auto norm = sqrt(norm_squared);
        variation.derivative = boosted_derivative->time * -1.0;
        for (int axis = 0; axis < 3; ++axis) {
            const auto component = n[axis] / norm;
            if (!isfinite(component)) return std::nullopt;
            variation.derivative += boosted_derivative->spatial[axis] * component.real +
                                    result.observer.spatial[axis] * component.dual;
        }
        // K is a coordinate derivative; V=K+Gamma(k,X) is the transported
        // Jacobi derivative and transforms as a vector at a chart boundary.
        for (int mu = 0; mu < 4; ++mu)
            for (int nu = 0; nu < 4; ++nu)
                for (int axis = 0; axis < 4; ++axis)
                    variation.derivative(mu) += connection.gamma(mu, nu, axis).real *
                                                result.tangent(nu) * variation.displacement(axis);
        if (!relativity::IsFinite(variation.displacement) ||
            !relativity::IsFinite(variation.derivative))
            return std::nullopt;
    }
    return result;
}

}  // namespace sirius::core
