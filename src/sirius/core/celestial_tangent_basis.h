#pragma once

// One deterministic orientation authority for a tangent plane on the unit
// celestial sphere. Ray-bundle ellipse angles and point-catalogue filtering
// must use the same basis or an anisotropic footprint is silently rotated.

#include <array>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <optional>

namespace sirius::core::relativity {

template <std::floating_point Scalar>
struct CelestialSeparationValue {
    Scalar angle;
    Scalar sine;
    std::array<Scalar, 3> normal;
};

// Unlike acos(dot), this retains first-order angular information when the
// float dot product rounds to one. Subtracting the nearly parallel vectors
// before the cross product also avoids cancellation between order-one terms.
// The atan2 ratio is invariant to the catalogue's permitted norm error.
template <std::floating_point Scalar>
[[nodiscard]] inline CelestialSeparationValue<Scalar> MeasureCelestialSeparation(
    std::array<Scalar, 3> direction, const std::array<Scalar, 3>& star) {
    const Scalar norm_squared =
        direction[0] * direction[0] + direction[1] * direction[1] + direction[2] * direction[2];
    // Public callers may supply any represented nonzero direction magnitude.
    // Keep subtraction on comparable scales without renormalizing the already
    // near-unit renderer inputs (which would perturb their subpixel direction).
    if (norm_squared < 0.25f || norm_squared > 4.0f) {
        const Scalar norm = std::sqrt(norm_squared);
        for (Scalar& component : direction) component /= norm;
    }
    const std::array<Scalar, 3> delta{star[0] - direction[0], star[1] - direction[1],
                                      star[2] - direction[2]};
    const std::array<Scalar, 3> normal{direction[1] * delta[2] - direction[2] * delta[1],
                                       direction[2] * delta[0] - direction[0] * delta[2],
                                       direction[0] * delta[1] - direction[1] * delta[0]};
    const Scalar sine =
        std::sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);
    const Scalar cosine = direction[0] * star[0] + direction[1] * star[1] + direction[2] * star[2];
    return {std::atan2(sine, cosine), sine, normal};
}

// Preserve the established float API, including calls with braced vectors.
using CelestialSeparation = CelestialSeparationValue<float>;
[[nodiscard]] inline CelestialSeparation MeasureCelestialSeparation(
    std::array<float, 3> direction, const std::array<float, 3>& star) {
    return MeasureCelestialSeparation<float>(direction, star);
}

template <std::floating_point Scalar>
struct CelestialTangentBasis {
    std::array<Scalar, 3> first;
    std::array<Scalar, 3> second;
};

template <std::floating_point Scalar>
[[nodiscard]] inline std::optional<CelestialTangentBasis<Scalar>> MakeCelestialTangentBasis(
    const std::array<Scalar, 3>& input_direction) {
    Scalar norm_squared = Scalar(0);
    for (const Scalar component : input_direction) {
        if (!std::isfinite(component)) return std::nullopt;
        norm_squared += component * component;
    }
    if (!std::isfinite(norm_squared) || !(norm_squared > Scalar(0))) return std::nullopt;

    const Scalar inverse_norm = Scalar(1) / std::sqrt(norm_squared);
    std::array<Scalar, 3> direction{};
    for (std::size_t component = 0; component < direction.size(); ++component) {
        direction[component] = input_direction[component] * inverse_norm;
    }

    // Project the coordinate axis least aligned with the direction. This keeps
    // the projection well conditioned and fixes ties to x, then y, then z.
    std::size_t reference_index = 0;
    for (std::size_t component = 1; component < direction.size(); ++component) {
        if (std::abs(input_direction[component]) < std::abs(input_direction[reference_index])) {
            reference_index = component;
        }
    }

    CelestialTangentBasis<Scalar> basis{};
    basis.first[reference_index] = Scalar(1);
    const Scalar projection = direction[reference_index];
    Scalar first_norm_squared = Scalar(0);
    for (std::size_t component = 0; component < direction.size(); ++component) {
        basis.first[component] -= projection * direction[component];
        first_norm_squared += basis.first[component] * basis.first[component];
    }
    if (!std::isfinite(first_norm_squared) || !(first_norm_squared > Scalar(0))) {
        return std::nullopt;
    }
    const Scalar first_inverse_norm = Scalar(1) / std::sqrt(first_norm_squared);
    for (Scalar& component : basis.first) component *= first_inverse_norm;

    basis.second = {direction[1] * basis.first[2] - direction[2] * basis.first[1],
                    direction[2] * basis.first[0] - direction[0] * basis.first[2],
                    direction[0] * basis.first[1] - direction[1] * basis.first[0]};
    for (const Scalar component : basis.second) {
        if (!std::isfinite(component)) return std::nullopt;
    }
    return basis;
}

}  // namespace sirius::core::relativity
