#pragma once

// Metric-orthonormal camera frames and past-directed ray launch.
//
// The screen direction is specified in the instantaneous rest frame of the
// camera. An Eulerian reference frame (the unit future normal to the coordinate
// time slice) is first constructed from the requested
// spatial axes by Lorentzian Gram-Schmidt, then boosted by the camera's local
// three-velocity.  A backward ray is
//
//   k_past = -u_camera + n^i e_i_camera,
//
// so g(k,k)=0 and the corresponding physical future photon -k has unit
// frequency -(-k).u_camera=1 at launch.  This avoids treating coordinate
// components as an orthonormal frame and avoids replacing backward ray tracing
// with a distinct future-directed ingoing null family.

#include "sirius/core/tensor.h"

#include <array>
#include <cmath>
#include <optional>

namespace sirius::core::relativity {

struct ObserverFrame {
    Vec4 time;
    std::array<Vec4, 3> spatial;
};

[[nodiscard]] inline bool IsFinite(const Vec4& vector) {
    for (int component = 0; component < 4; ++component) {
        if (!std::isfinite(vector(component))) return false;
    }
    return true;
}

// Project a trial vector into the physical two-plane seen by a timelike
// observer transverse to a ray.  For unit future observer u, first form the
// observer-spatial propagation direction from
//
//   k_perp = k + (k.u) u,
//
// then remove both u and the unit direction n=k_perp/|k_perp| from the trial.
// This construction is independent of whether k is future- or past-directed,
// remains exactly transverse when a numerical ray has small null drift, and
// never treats a chart basis vector as an observer worldline.
[[nodiscard]] inline std::optional<Vec4> ProjectToObserverScreen(
    const Metric4d& metric, const Vec4& observer, const Vec4& ray, const Vec4& trial) {
    const auto dot = [&metric](const Vec4& lhs, const Vec4& rhs) {
        return TensorOps::InnerProduct(lhs, rhs, metric);
    };
    if (!IsFinite(observer) || !IsFinite(ray) || !IsFinite(trial)) return std::nullopt;

    const double observer_norm = dot(observer, observer);
    if (!std::isfinite(observer_norm) || std::abs(observer_norm + 1.0) > 1.0e-10) {
        return std::nullopt;
    }

    Vec4 propagation = ray + observer * dot(ray, observer);
    const double propagation_norm = dot(propagation, propagation);
    if (!std::isfinite(propagation_norm) || !(propagation_norm > 1.0e-20)) {
        return std::nullopt;
    }
    propagation = propagation / std::sqrt(propagation_norm);

    Vec4 screen = trial + observer * dot(trial, observer);
    screen -= propagation * dot(screen, propagation);
    const double screen_norm = dot(screen, screen);
    if (!std::isfinite(screen_norm) || !(screen_norm > 1.0e-20)) {
        return std::nullopt;
    }
    screen = screen / std::sqrt(screen_norm);
    if (!IsFinite(screen)) return std::nullopt;
    return screen;
}

[[nodiscard]] inline std::optional<ObserverFrame> EulerianObserverFrame(
    const Metric4d& metric, const Metric4d& inverse_metric,
    const std::array<Vec4, 3>& spatial_seeds) {
    const double inverse_g_tt = inverse_metric(0, 0).real;
    if (!std::isfinite(inverse_g_tt) || !(inverse_g_tt < 0.0)) return std::nullopt;

    ObserverFrame frame;
    const double lapse = 1.0 / std::sqrt(-inverse_g_tt);
    for (int component = 0; component < 4; ++component) {
        frame.time(component) = -lapse * inverse_metric(component, 0).real;
    }
    if (!IsFinite(frame.time) || !(frame.time(0) > 0.0)) return std::nullopt;
    const auto dot = [&metric](const Vec4& lhs, const Vec4& rhs) {
        return TensorOps::InnerProduct(lhs, rhs, metric);
    };

    for (std::size_t index = 0; index < frame.spatial.size(); ++index) {
        Vec4 basis = spatial_seeds[index];
        basis += frame.time * dot(basis, frame.time);
        for (std::size_t previous = 0; previous < index; ++previous) {
            basis -= frame.spatial[previous] * dot(basis, frame.spatial[previous]);
        }
        const double norm_squared = dot(basis, basis);
        if (!std::isfinite(norm_squared) || !(norm_squared > 1.0e-20)) {
            return std::nullopt;
        }
        frame.spatial[index] = basis / std::sqrt(norm_squared);
    }
    return frame;
}

[[nodiscard]] inline std::optional<ObserverFrame> BoostObserverFrame(
    const ObserverFrame& reference, const std::array<double, 3>& beta) {
    const double beta_squared =
        beta[0] * beta[0] + beta[1] * beta[1] + beta[2] * beta[2];
    if (!std::isfinite(beta_squared) || beta_squared >= 1.0) return std::nullopt;
    if (beta_squared == 0.0) return reference;

    const double gamma = 1.0 / std::sqrt(1.0 - beta_squared);
    Vec4 beta_vector;
    for (std::size_t index = 0; index < reference.spatial.size(); ++index) {
        beta_vector += reference.spatial[index] * beta[index];
    }

    ObserverFrame boosted;
    boosted.time = (reference.time + beta_vector) * gamma;
    const double spatial_coefficient = (gamma - 1.0) / beta_squared;
    for (std::size_t index = 0; index < boosted.spatial.size(); ++index) {
        boosted.spatial[index] = reference.spatial[index] +
                                 beta_vector * (spatial_coefficient * beta[index]) +
                                 reference.time * (gamma * beta[index]);
    }
    return boosted;
}

[[nodiscard]] inline std::optional<Vec4> PastDirectedCameraRay(
    const ObserverFrame& frame, const std::array<double, 3>& rest_direction) {
    const double direction_norm = rest_direction[0] * rest_direction[0] +
                                  rest_direction[1] * rest_direction[1] +
                                  rest_direction[2] * rest_direction[2];
    if (!std::isfinite(direction_norm) || !(direction_norm > 0.0)) return std::nullopt;

    Vec4 ray = frame.time * -1.0;
    const double inverse_norm = 1.0 / std::sqrt(direction_norm);
    for (std::size_t index = 0; index < frame.spatial.size(); ++index) {
        ray += frame.spatial[index] * (rest_direction[index] * inverse_norm);
    }
    if (!IsFinite(ray)) return std::nullopt;
    return ray;
}

// Orthonormal Sachs screen at an observer. The input direction is resolved in
// the observer's orthonormal rest frame. Both returned vectors are spacelike,
// unit normal, mutually orthogonal, and orthogonal to the observer and the null
// ray -u+n. Constructing the screen in rest-frame components avoids treating
// chart components as Euclidean vectors in a curved or shifted metric.
[[nodiscard]] inline std::optional<std::array<Vec4, 2>> ObserverScreenBasis(
    const ObserverFrame& frame, const std::array<double, 3>& rest_direction) {
    double direction_norm_squared = 0.0;
    for (const double component : rest_direction) {
        direction_norm_squared += component * component;
    }
    if (!std::isfinite(direction_norm_squared) || !(direction_norm_squared > 0.0)) {
        return std::nullopt;
    }

    std::array<double, 3> direction{};
    const double inverse_norm = 1.0 / std::sqrt(direction_norm_squared);
    for (std::size_t component = 0; component < direction.size(); ++component) {
        direction[component] = rest_direction[component] * inverse_norm;
    }

    std::size_t reference_index = 0;
    for (std::size_t component = 1; component < direction.size(); ++component) {
        if (std::abs(direction[component]) < std::abs(direction[reference_index])) {
            reference_index = component;
        }
    }

    std::array<double, 3> first_components{};
    first_components[reference_index] = 1.0;
    const double projection = direction[reference_index];
    double first_norm_squared = 0.0;
    for (std::size_t component = 0; component < direction.size(); ++component) {
        first_components[component] -= projection * direction[component];
        first_norm_squared += first_components[component] * first_components[component];
    }
    if (!std::isfinite(first_norm_squared) || !(first_norm_squared > 0.0)) {
        return std::nullopt;
    }
    const double first_inverse_norm = 1.0 / std::sqrt(first_norm_squared);
    for (double& component : first_components) component *= first_inverse_norm;

    const std::array<double, 3> second_components{
        direction[1] * first_components[2] - direction[2] * first_components[1],
        direction[2] * first_components[0] - direction[0] * first_components[2],
        direction[0] * first_components[1] - direction[1] * first_components[0]};

    std::array<Vec4, 2> screen{};
    for (std::size_t component = 0; component < direction.size(); ++component) {
        screen[0] += frame.spatial[component] * first_components[component];
        screen[1] += frame.spatial[component] * second_components[component];
    }
    if (!IsFinite(screen[0]) || !IsFinite(screen[1])) return std::nullopt;
    return screen;
}

}  // namespace sirius::core::relativity
