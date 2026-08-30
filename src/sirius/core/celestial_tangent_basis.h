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
        if (std::abs(direction[component]) < std::abs(direction[reference_index])) {
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
