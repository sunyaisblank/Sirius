#pragma once

#include "sirius/core/observer_frame.h"

#include <array>
#include <cmath>
#include <limits>
#include <optional>

namespace sirius::core::relativity {

struct SphericalEndpointVariation {
    Vec4 displacement;
    // |x.k| / sum |x_i k_i| measures cancellation in the boundary derivative.
    // It is a numerical conditioning diagnostic, not a lens magnification.
    double radial_condition;
};

// Move a fixed-affine deviation to the same spherical boundary as the central
// ray: X = xi - k (x.xi)/(x.k). For a geodesic the covariant tangent variation
// remains V. The event itself must already have been localized by the tracer.
// A radial derivative indistinguishable from zero at roundoff cannot define
// this implicit derivative. No denominator floor fabricates a finite map.
[[nodiscard]] inline std::optional<SphericalEndpointVariation> VarySphericalEndpoint(
    const Vec4& position, const Vec4& tangent, const Vec4& deviation) {
    if (!IsFinite(position) || !IsFinite(tangent) || !IsFinite(deviation)) return std::nullopt;
    double radial_derivative = 0.0;
    double radial_scale = 0.0;
    double displacement_derivative = 0.0;
    for (int axis = 1; axis < 4; ++axis) {
        const double product = position(axis) * tangent(axis);
        radial_derivative += product;
        radial_scale += std::abs(product);
        displacement_derivative += position(axis) * deviation(axis);
    }
    constexpr double roundoff_factor = 8.0 * std::numeric_limits<double>::epsilon();
    if (!std::isfinite(radial_derivative) || !std::isfinite(radial_scale) ||
        !(std::abs(radial_derivative) > roundoff_factor * radial_scale)) {
        return std::nullopt;
    }
    const double shift = -displacement_derivative / radial_derivative;
    const Vec4 displacement = deviation + tangent * shift;
    if (!IsFinite(displacement)) return std::nullopt;
    return SphericalEndpointVariation{displacement, std::abs(radial_derivative) / radial_scale};
}

struct SourceSkyAngularMap {
    std::array<double, 3> direction;
    // Rows: fixed central source celestial tangent basis. Columns: the two
    // launch angular directions whose transported amplitudes contain seed.
    std::array<std::array<double, 2>, 2> jacobian;
    double determinant;
};

// Local angular derivative of the actual normalized Eulerian source direction.
// All inputs must be in the same source chart. Endpoint displacements include
// the event shift; covariant tangent variations transform as geometric vectors.
// Singular matrices remain represented here: filter validity/refinement is a
// separate decision and cannot be repaired by clipping this determinant.
[[nodiscard]] inline std::optional<SourceSkyAngularMap> MeasureSourceSkyAngularMap(
    const Metric4d& metric, const Metric4d& inverse_metric,
    const Tensor<Dual<double>, 4, 4, 4>& metric_derivatives, const Vec4& tangent,
    const std::array<Vec4, 2>& endpoint_displacements,
    const std::array<Vec4, 2>& covariant_tangent_variations, double seed) {
    if (!std::isfinite(seed) || !(seed > 0.0)) return std::nullopt;
    const auto connection = TensorOps::Christoffel(metric, metric_derivatives);
    std::array<Vec4, 3> spatial_seeds;
    spatial_seeds[0](1) = 1.0;
    spatial_seeds[1](2) = 1.0;
    spatial_seeds[2](3) = 1.0;
    SourceSkyAngularMap result{};
    std::optional<CelestialTangentBasis<double>> basis;
    for (std::size_t column = 0; column < endpoint_displacements.size(); ++column) {
        const Vec4& displacement = endpoint_displacements[column];
        Vec4 coordinate_variation = covariant_tangent_variations[column];
        if (!IsFinite(displacement) || !IsFinite(coordinate_variation)) return std::nullopt;
        Metric4d metric_variation;
        for (int mu = 0; mu < 4; ++mu) {
            for (int nu = 0; nu < 4; ++nu) {
                for (int axis = 0; axis < 4; ++axis) {
                    metric_variation(mu, nu).real +=
                        metric_derivatives(axis, mu, nu).real * displacement(axis);
                    coordinate_variation(mu) -=
                        connection.gamma(mu, nu, axis).real * tangent(nu) * displacement(axis);
                }
            }
        }
        const auto differential = EulerianSourceSkyDifferential(
            metric, inverse_metric, metric_variation, tangent, coordinate_variation, spatial_seeds);
        if (!differential) return std::nullopt;
        if (column == 0) {
            result.direction = differential->direction;
            basis = MakeCelestialTangentBasis(result.direction);
            if (!basis) return std::nullopt;
        }
        for (std::size_t axis = 0; axis < result.direction.size(); ++axis) {
            result.jacobian[0][column] +=
                basis->first[axis] * differential->derivative[axis] / seed;
            result.jacobian[1][column] +=
                basis->second[axis] * differential->derivative[axis] / seed;
        }
        if (!std::isfinite(result.jacobian[0][column]) ||
            !std::isfinite(result.jacobian[1][column])) {
            return std::nullopt;
        }
    }
    result.determinant = result.jacobian[0][0] * result.jacobian[1][1] -
                         result.jacobian[0][1] * result.jacobian[1][0];
    if (!std::isfinite(result.determinant)) return std::nullopt;
    return result;
}

}  // namespace sirius::core::relativity
