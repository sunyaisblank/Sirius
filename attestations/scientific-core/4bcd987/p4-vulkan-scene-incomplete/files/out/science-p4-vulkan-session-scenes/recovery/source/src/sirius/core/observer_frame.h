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

#include "sirius/core/celestial_tangent_basis.h"
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
[[nodiscard]] inline std::optional<Vec4> ProjectToObserverScreen(const Metric4d& metric,
                                                                 const Vec4& observer,
                                                                 const Vec4& ray,
                                                                 const Vec4& trial) {
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

namespace detail {

template <typename Scalar>
struct EulerianFrameValues {
    Tensor<Scalar, 4> time;
    std::array<Tensor<Scalar, 4>, 3> spatial;
};

inline double FrameValue(double value) { return value; }
inline double FrameValue(const Dual<double>& value) { return value.real; }

template <typename Scalar>
[[nodiscard]] bool FiniteFrameValue(const Scalar& value) {
    using std::isfinite;
    return isfinite(value);
}

template <typename Scalar>
[[nodiscard]] Scalar FrameInnerProduct(const Tensor<Scalar, 4>& lhs, const Tensor<Scalar, 4>& rhs,
                                       const Tensor<Scalar, 4, 4>& metric) {
    Scalar result(0.0);
    for (int mu = 0; mu < 4; ++mu) {
        for (int nu = 0; nu < 4; ++nu) {
            result += metric(mu, nu) * lhs(mu) * rhs(nu);
        }
    }
    return result;
}

// The nominal frame and its directional derivative use exactly the same
// lapse and Gram-Schmidt operations. Seeds are fixed chart vectors; their
// derivatives are zero. A dual scalar carries the metric/frame variation.
template <typename Scalar>
[[nodiscard]] std::optional<EulerianFrameValues<Scalar>> BuildEulerianFrame(
    const Tensor<Scalar, 4, 4>& metric, const Tensor<Scalar, 4, 4>& inverse_metric,
    const std::array<Vec4, 3>& spatial_seeds) {
    const Scalar inverse_g_tt = inverse_metric(0, 0);
    if (!FiniteFrameValue(inverse_g_tt) || !(FrameValue(inverse_g_tt) < 0.0)) {
        return std::nullopt;
    }
    using std::sqrt;
    EulerianFrameValues<Scalar> frame;
    const Scalar lapse = Scalar(1.0) / sqrt(-inverse_g_tt);
    for (int component = 0; component < 4; ++component) {
        frame.time(component) = -lapse * inverse_metric(component, 0);
        if (!FiniteFrameValue(frame.time(component))) return std::nullopt;
    }
    if (!(FrameValue(frame.time(0)) > 0.0)) return std::nullopt;
    const auto dot = [&metric](const auto& lhs, const auto& rhs) {
        return FrameInnerProduct(lhs, rhs, metric);
    };
    for (std::size_t index = 0; index < frame.spatial.size(); ++index) {
        Tensor<Scalar, 4> basis;
        for (int component = 0; component < 4; ++component) {
            basis(component) = Scalar(spatial_seeds[index](component));
        }
        const Scalar time_projection = dot(basis, frame.time);
        for (int component = 0; component < 4; ++component) {
            basis(component) += frame.time(component) * time_projection;
        }
        for (std::size_t previous = 0; previous < index; ++previous) {
            const Scalar projection = dot(basis, frame.spatial[previous]);
            for (int component = 0; component < 4; ++component) {
                basis(component) -= frame.spatial[previous](component) * projection;
            }
        }
        const Scalar norm_squared = dot(basis, basis);
        if (!FiniteFrameValue(norm_squared) || !(FrameValue(norm_squared) > 1.0e-20)) {
            return std::nullopt;
        }
        const Scalar norm = sqrt(norm_squared);
        for (int component = 0; component < 4; ++component) {
            frame.spatial[index](component) = basis(component) / norm;
            if (!FiniteFrameValue(frame.spatial[index](component))) return std::nullopt;
        }
    }
    return frame;
}

}  // namespace detail

[[nodiscard]] inline std::optional<ObserverFrame> EulerianObserverFrame(
    const Metric4d& metric, const Metric4d& inverse_metric,
    const std::array<Vec4, 3>& spatial_seeds) {
    Tensor<double, 4, 4> values;
    Tensor<double, 4, 4> inverse_values;
    for (int mu = 0; mu < 4; ++mu) {
        for (int nu = 0; nu < 4; ++nu) {
            values(mu, nu) = metric(mu, nu).real;
            inverse_values(mu, nu) = inverse_metric(mu, nu).real;
        }
    }
    const auto frame = detail::BuildEulerianFrame(values, inverse_values, spatial_seeds);
    if (!frame) return std::nullopt;
    return ObserverFrame{frame->time, frame->spatial};
}

struct SourceSkyDifferential {
    std::array<double, 3> direction;
    std::array<double, 3> derivative;
    double frequency;
    double frequency_derivative;
};

// Differentiate the normalized direction seen in the source Eulerian frame.
// metric_derivative.real is dg[X]; tangent_derivative is the ordinary
// coordinate variation K, not a covariant Jacobi derivative V. The caller
// applies endpoint/chart changes before supplying K = V - Gamma(k,X).
// Existing Dual metadata on metric inputs is deliberately reseeded here.
[[nodiscard]] inline std::optional<SourceSkyDifferential> EulerianSourceSkyDifferential(
    const Metric4d& metric, const Metric4d& inverse_metric, const Metric4d& metric_derivative,
    const Vec4& past_tangent, const Vec4& tangent_derivative,
    const std::array<Vec4, 3>& spatial_seeds) {
    if (!IsFinite(past_tangent) || !IsFinite(tangent_derivative)) return std::nullopt;
    Metric4d varied_metric;
    Metric4d varied_inverse;
    for (int mu = 0; mu < 4; ++mu) {
        for (int nu = 0; nu < 4; ++nu) {
            const double value = metric(mu, nu).real;
            const double variation = metric_derivative(mu, nu).real;
            if (!std::isfinite(value) || !std::isfinite(variation) ||
                !std::isfinite(inverse_metric(mu, nu).real)) {
                return std::nullopt;
            }
            varied_metric(mu, nu) = Dual<double>(value, variation);
            double inverse_variation = 0.0;
            for (int alpha = 0; alpha < 4; ++alpha) {
                for (int beta = 0; beta < 4; ++beta) {
                    inverse_variation -= inverse_metric(mu, alpha).real *
                                         metric_derivative(alpha, beta).real *
                                         inverse_metric(beta, nu).real;
                }
            }
            varied_inverse(mu, nu) = Dual<double>(inverse_metric(mu, nu).real, inverse_variation);
        }
    }
    const auto frame = detail::BuildEulerianFrame(varied_metric, varied_inverse, spatial_seeds);
    if (!frame) return std::nullopt;
    Tensor<Dual<double>, 4> tangent;
    for (int component = 0; component < 4; ++component) {
        tangent(component) = Dual<double>(past_tangent(component), tangent_derivative(component));
    }
    const Dual<double> frequency = detail::FrameInnerProduct(tangent, frame->time, varied_metric);
    if (!isfinite(frequency) || !(frequency.real > 0.0)) return std::nullopt;
    std::array<Dual<double>, 3> local;
    Dual<double> squared_norm(0.0);
    for (std::size_t component = 0; component < local.size(); ++component) {
        local[component] =
            detail::FrameInnerProduct(tangent, frame->spatial[component], varied_metric) /
            frequency;
        squared_norm += local[component] * local[component];
    }
    if (!isfinite(squared_norm) || !(squared_norm.real > 0.0)) return std::nullopt;
    const Dual<double> norm = sqrt(squared_norm);
    SourceSkyDifferential result{};
    for (std::size_t component = 0; component < local.size(); ++component) {
        const Dual<double> direction = local[component] / norm;
        if (!isfinite(direction)) return std::nullopt;
        result.direction[component] = direction.real;
        result.derivative[component] = direction.dual;
    }
    result.frequency = frequency.real;
    result.frequency_derivative = frequency.dual;
    return result;
}

[[nodiscard]] inline std::optional<ObserverFrame> BoostObserverFrame(
    const ObserverFrame& reference, const std::array<double, 3>& beta) {
    const double beta_squared = beta[0] * beta[0] + beta[1] * beta[1] + beta[2] * beta[2];
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
    const auto components = MakeCelestialTangentBasis(rest_direction);
    if (!components.has_value()) return std::nullopt;

    std::array<Vec4, 2> screen{};
    for (std::size_t component = 0; component < rest_direction.size(); ++component) {
        screen[0] += frame.spatial[component] * components->first[component];
        screen[1] += frame.spatial[component] * components->second[component];
    }
    if (!IsFinite(screen[0]) || !IsFinite(screen[1])) return std::nullopt;
    return screen;
}

}  // namespace sirius::core::relativity
