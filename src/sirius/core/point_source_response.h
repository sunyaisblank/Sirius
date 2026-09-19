#pragma once

#include <array>
#include <cmath>
#include <expected>
#include <limits>
#include <numbers>

namespace sirius::core {

using AngularMatrix2 = std::array<std::array<double, 2>, 2>;

enum class PointResponseFailure { InvalidInput, NeedsRefinement, Arithmetic };

// Local affine response only. Its caller must separately establish that the
// source map, camera angular-area density and visibility are sufficiently
// constant over the detector support. The area diagnostics below do not bound
// pointwise relative density changes across the discontinuous support boundary.
// The film stores radiance: integrating this density over image solid angle
// yields source flux / abs(det(source angular map)), without another area factor.
struct AffinePointResponse {
    AngularMatrix2 source_to_standard_normal;
    double density_scale;          // Per source tangent-plane angular area.
    double major_sigma;            // Radians; conservative index query radius / 4.
    double support_radius;         // Radians; exact support is elliptical q <= 16.
    double spherical_area_bound;   // Bounds 1 - sin(radius) / radius on support.
    double arithmetic_area_bound;  // Relative area uncertainty from represented matrix arithmetic.
    double signed_lens_determinant;

    [[nodiscard]] std::expected<double, PointResponseFailure> Density(
        const std::array<double, 2>& source_offset) const {
        if (!std::isfinite(source_offset[0]) || !std::isfinite(source_offset[1])) {
            return std::unexpected(PointResponseFailure::InvalidInput);
        }
        const double a = source_to_standard_normal[0][0] * source_offset[0] +
                         source_to_standard_normal[0][1] * source_offset[1];
        const double b = source_to_standard_normal[1][0] * source_offset[0] +
                         source_to_standard_normal[1][1] * source_offset[1];
        // hypot avoids overflowing a squared offset outside the compact support.
        const double radius = std::hypot(a, b);
        if (std::isnan(radius)) return std::unexpected(PointResponseFailure::Arithmetic);
        if (radius > 4.0) return 0.0;
        const double result = density_scale * std::exp(-0.5 * radius * radius);
        if (!std::isfinite(result)) return std::unexpected(PointResponseFailure::Arithmetic);
        return result;
    }
};

// Density at an already traced image root. Discovery has its own finite-cell
// geometry checks; a local density does not approximate a finite sky footprint.
struct PointImageResponse {
    double density_scale;
    double arithmetic_area_bound;

    [[nodiscard]] std::expected<double, PointResponseFailure> DensityAtOriginalRoot(
        const std::array<double, 2>& original) const {
        if (!std::isfinite(original[0]) || !std::isfinite(original[1]))
            return std::unexpected(PointResponseFailure::InvalidInput);
        const double radius = std::hypot(original[0], original[1]);
        if (radius > 4.0) return 0.0;
        const double result = density_scale * std::exp(-0.5 * radius * radius);
        if (!(result > 0.0) || !std::isfinite(result))
            return std::unexpected(PointResponseFailure::Arithmetic);
        return result;
    }
};

namespace point_response_detail {

// P maps a film-pixel displacement at fixed pupil to the launch rest-frame
// angular basis. J maps that SAME basis to source angles and already includes
// observer aberration. A circular film Gaussian of sigma pixels has source
// covariance B B^T, B = sigma J P. No source-plane axis/determinant floor is used.
[[nodiscard]] inline std::expected<AffinePointResponse, PointResponseFailure> BuildPointResponse(
    const AngularMatrix2& source_map, const AngularMatrix2& film_map, double film_sigma,
    double geometry_area_budget, bool finite_footprint) {
    if (!std::isfinite(film_sigma) || !(film_sigma > 0.0) || !std::isfinite(geometry_area_budget) ||
        !(geometry_area_budget > 0.0) || !(geometry_area_budget < 1.0)) {
        return std::unexpected(PointResponseFailure::InvalidInput);
    }
    for (const auto* matrix : {&source_map, &film_map}) {
        for (const auto& row : *matrix) {
            for (double value : row) {
                if (!std::isfinite(value)) {
                    return std::unexpected(PointResponseFailure::InvalidInput);
                }
            }
        }
    }
    AngularMatrix2 beam{};
    AngularMatrix2 composition_error{};
    constexpr double epsilon = std::numeric_limits<double>::epsilon();
    constexpr double tiny = std::numeric_limits<double>::denorm_min();
    double scale = 0.0;
    for (int row = 0; row < 2; ++row) {
        for (int col = 0; col < 2; ++col) {
            const double first = source_map[row][0] * film_map[0][col];
            const double second = source_map[row][1] * film_map[1][col];
            beam[row][col] = film_sigma * (first + second);
            // Bound both products, their sum and sigma multiplication. Checking
            // only the final determinant misses cancellation while composing J P.
            composition_error[row][col] =
                film_sigma * (8 * epsilon * (std::abs(first) + std::abs(second)) + 8 * tiny) +
                4 * epsilon * std::abs(beam[row][col]) + 8 * tiny;
            if (!std::isfinite(beam[row][col]) || !std::isfinite(composition_error[row][col])) {
                return std::unexpected(PointResponseFailure::Arithmetic);
            }
            scale = std::fmax(scale, std::abs(beam[row][col]));
        }
    }
    if (!(scale > 0.0)) return std::unexpected(PointResponseFailure::NeedsRefinement);
    double maximum_error = 0.0;
    for (int row = 0; row < 2; ++row) {
        for (int col = 0; col < 2; ++col) {
            beam[row][col] /= scale;
            maximum_error = std::fmax(maximum_error, composition_error[row][col] / scale +
                                                         2 * epsilon * std::abs(beam[row][col]));
        }
    }
    const double positive = beam[0][0] * beam[1][1];
    const double negative = beam[0][1] * beam[1][0];
    const double determinant = std::fma(beam[0][0], beam[1][1], -negative);
    const double determinant_uncertainty =
        8.0 * std::numeric_limits<double>::epsilon() * (std::abs(positive) + std::abs(negative));
    if (!(std::abs(determinant) > determinant_uncertainty)) {
        return std::unexpected(PointResponseFailure::NeedsRefinement);
    }
    // Largest singular value from the eigenvalues of B B^T, without forming
    // the cancellation-prone small eigenvalue. Inversion retains the full map.
    const double aa = beam[0][0] * beam[0][0] + beam[0][1] * beam[0][1];
    const double bb = beam[1][0] * beam[1][0] + beam[1][1] * beam[1][1];
    const double ab = beam[0][0] * beam[1][0] + beam[0][1] * beam[1][1];
    // ||B^-1 E||_2 <= ||B^-1||_F ||E||_F. For this scaled 2x2 matrix
    // ||B^-1||_F=sqrt(aa+bb)/|det B| and ||E||_F<=2 max|E_ij|.
    // If rho<1, determinant perturbation is bounded by (1 +/- rho)^2.
    const double rho =
        (2 * maximum_error * std::sqrt(aa + bb) + determinant_uncertainty) / std::abs(determinant) +
        16 * epsilon;
    if (!std::isfinite(rho) || !(rho < 1)) {
        return std::unexpected(PointResponseFailure::NeedsRefinement);
    }
    const double arithmetic_area_bound = std::expm1(-2 * std::log1p(-rho));
    const double major = scale * std::sqrt(0.5 * (aa + bb + std::hypot(aa - bb, 2 * ab)));
    const double support = 4.0 * major;
    const double area_bound = support * support / 6.0;
    if (arithmetic_area_bound > geometry_area_budget ||
        (finite_footprint &&
         (!std::isfinite(area_bound) ||
          area_bound + arithmetic_area_bound + area_bound * arithmetic_area_bound >
              geometry_area_budget))) {
        return std::unexpected(PointResponseFailure::NeedsRefinement);
    }
    // The normalization is exact for the elliptical truncation, not for the
    // major-axis circle used only as a conservative spatial-index bound.
    const double area = std::abs(determinant) * scale * scale;
    const double normalization = 2.0 * std::numbers::pi * area * (-std::expm1(-8.0));
    const double lens_determinant =
        std::fma(source_map[0][0], source_map[1][1], -source_map[0][1] * source_map[1][0]);
    AffinePointResponse response{
        {{{beam[1][1] / determinant / scale, -beam[0][1] / determinant / scale},
          {-beam[1][0] / determinant / scale, beam[0][0] / determinant / scale}}},
        1.0 / normalization,
        major,
        support,
        area_bound,
        arithmetic_area_bound,
        lens_determinant};
    if (!(normalization > 0.0) || !std::isfinite(normalization) ||
        !std::isfinite(response.density_scale) || !std::isfinite(lens_determinant) ||
        lens_determinant == 0) {
        return std::unexpected(PointResponseFailure::Arithmetic);
    }
    for (const auto& row : response.source_to_standard_normal) {
        for (double value : row) {
            if (!std::isfinite(value)) return std::unexpected(PointResponseFailure::Arithmetic);
        }
    }
    return response;
}

}  // namespace point_response_detail

[[nodiscard]] inline std::expected<AffinePointResponse, PointResponseFailure>
MakeAffinePointResponse(const AngularMatrix2& source_map, const AngularMatrix2& film_map,
                        double film_sigma, double geometry_area_budget) {
    return point_response_detail::BuildPointResponse(source_map, film_map, film_sigma,
                                                     geometry_area_budget, true);
}

// source_map is d(source angle)/d(chart coordinate) at the actual image.
// chart_from_standard maps this original sample's unit Gaussian into that
// shared chart. The determinant includes both maps exactly once.
[[nodiscard]] inline std::expected<PointImageResponse, PointResponseFailure> MakePointImageResponse(
    const AngularMatrix2& source_map, const AngularMatrix2& chart_from_standard,
    double arithmetic_area_budget) {
    const auto response = point_response_detail::BuildPointResponse(
        source_map, chart_from_standard, 1.0, arithmetic_area_budget, false);
    if (!response) return std::unexpected(response.error());
    return PointImageResponse{response->density_scale, response->arithmetic_area_bound};
}

// One affine leaf of the ORIGINAL packet's standardized Gaussian support.
// Its source map is d(source angle)/dz, where film = original_film + L*z.
// The query ellipse encloses this rectangular cell; its Gaussian is never
// accumulated. Density uses the unchanged parent Gaussian and normalization.
// The caller owns nonlinear/visibility estimates and uncertain root ownership.
struct RestrictedAffinePointResponse {
    AffinePointResponse query;
    std::array<double, 2> lower;
    std::array<double, 2> upper;
    std::array<double, 2> centre;
    double query_scale;

    [[nodiscard]] std::expected<double, PointResponseFailure> Density(
        const std::array<double, 2>& source_offset) const {
        if (!std::isfinite(source_offset[0]) || !std::isfinite(source_offset[1])) {
            return std::unexpected(PointResponseFailure::InvalidInput);
        }
        std::array<double, 2> original{}, local{}, term_size{}, error{};
        for (int axis = 0; axis < 2; ++axis) {
            const double first = query.source_to_standard_normal[axis][0] * source_offset[0];
            const double second = query.source_to_standard_normal[axis][1] * source_offset[1];
            local[axis] = first + second;
            term_size[axis] = std::abs(first) + std::abs(second);
            original[axis] = std::fma(query_scale, local[axis], centre[axis]);
            if (!std::isfinite(original[axis]) || !std::isfinite(term_size[axis])) {
                return std::unexpected(PointResponseFailure::Arithmetic);
            }
        }
        const double inverse_error =
            query.arithmetic_area_bound * query_scale * std::hypot(local[0], local[1]);
        constexpr double epsilon = std::numeric_limits<double>::epsilon();
        for (int axis = 0; axis < 2; ++axis) {
            // The query's determinant perturbation bound also dominates its
            // inverse relative perturbation. Add represented dot/FMA rounding.
            error[axis] = inverse_error +
                          32 * epsilon * (std::abs(centre[axis]) + query_scale * term_size[axis]) +
                          64 * std::numeric_limits<double>::denorm_min();
            if (!std::isfinite(error[axis])) {
                return std::unexpected(PointResponseFailure::Arithmetic);
            }
            if (original[axis] < lower[axis] - error[axis] ||
                original[axis] > upper[axis] + error[axis])
                return 0.0;
        }
        for (int axis = 0; axis < 2; ++axis) {
            if (std::abs(original[axis] - lower[axis]) <= error[axis] ||
                std::abs(original[axis] - upper[axis]) <= error[axis]) {
                return std::unexpected(PointResponseFailure::NeedsRefinement);
            }
        }
        const double radius = std::hypot(original[0], original[1]);
        const double radius_error = std::hypot(error[0], error[1]) + 8 * epsilon * radius;
        if (std::abs(radius - 4.0) <= radius_error) {
            return std::unexpected(PointResponseFailure::NeedsRefinement);
        }
        return DensityAtOriginalRoot(original);
    }

    // Use one shared root expressed in original z coordinates for adjacent
    // cells. The caller must resolve physical/root uncertainty before supplying
    // this represented value; independently rounded local inverses cannot decide
    // ownership at a shared edge. No epsilon overlap or filter renormalization.
    [[nodiscard]] std::expected<double, PointResponseFailure> DensityAtOriginalRoot(
        const std::array<double, 2>& original) const {
        if (!std::isfinite(original[0]) || !std::isfinite(original[1])) {
            return std::unexpected(PointResponseFailure::InvalidInput);
        }
        for (int axis = 0; axis < 2; ++axis) {
            // Internal upper edges belong to the adjacent cell. The parent's
            // outer edge remains included, matching its compact |z| <= 4 PSF.
            if (original[axis] < lower[axis] || original[axis] > upper[axis] ||
                (original[axis] == upper[axis] && upper[axis] != 4.0)) {
                return 0.0;
            }
        }
        const double radius = std::hypot(original[0], original[1]);
        if (radius > 4.0) return 0.0;
        // query.density_scale includes 1/query_scale^2. Undo only that local
        // query scale: every leaf retains the original unit-Gaussian weight.
        const double parent_scale = query.density_scale * query_scale * query_scale;
        const double density = parent_scale * std::exp(-0.5 * radius * radius);
        if (!(density > 0.0) || !std::isfinite(density)) {
            return std::unexpected(PointResponseFailure::Arithmetic);
        }
        return density;
    }
};

[[nodiscard]] inline std::expected<RestrictedAffinePointResponse, PointResponseFailure>
MakeRestrictedAffinePointResponse(const AngularMatrix2& source_map,
                                  const std::array<double, 2>& lower,
                                  const std::array<double, 2>& upper, double geometry_area_budget) {
    std::array<double, 2> centre{}, half_width{};
    for (int axis = 0; axis < 2; ++axis) {
        if (!std::isfinite(lower[axis]) || !std::isfinite(upper[axis]) || lower[axis] < -4.0 ||
            upper[axis] > 4.0 || !(lower[axis] < upper[axis])) {
            return std::unexpected(PointResponseFailure::InvalidInput);
        }
        half_width[axis] = (upper[axis] - lower[axis]) * 0.5;
        centre[axis] = lower[axis] + half_width[axis];
        if (!(centre[axis] > lower[axis] && centre[axis] < upper[axis])) {
            return std::unexpected(PointResponseFailure::NeedsRefinement);
        }
        // A rounded midpoint need not bisect the represented endpoints (for
        // example, a three-ULP interval). Enclose both actual displacements.
        half_width[axis] =
            std::nextafter(std::fmax(centre[axis] - lower[axis], upper[axis] - centre[axis]),
                           std::numeric_limits<double>::infinity());
    }
    // 4*scale is the cell's half-diagonal. The largest singular value of
    // source_map therefore bounds its affine image about the source centre.
    const double scale = std::hypot(half_width[0], half_width[1]) / 4.0;
    constexpr AngularMatrix2 identity{{{1.0, 0.0}, {0.0, 1.0}}};
    auto query = MakeAffinePointResponse(source_map, identity, scale, geometry_area_budget);
    if (!query) return std::unexpected(query.error());
    // Outward query rounding covers corners whose norm rounds above the SVD
    // radius. This padding changes candidate lookup only, never parent density.
    constexpr double epsilon = std::numeric_limits<double>::epsilon();
    constexpr double tiny = std::numeric_limits<double>::denorm_min();
    double corner_bound = 0.0;
    for (double sign_x : {-1.0, 1.0}) {
        for (double sign_y : {-1.0, 1.0}) {
            std::array<double, 2> corner{}, roundoff{};
            for (int row = 0; row < 2; ++row) {
                const double first = source_map[row][0] * sign_x * half_width[0];
                const double second = source_map[row][1] * sign_y * half_width[1];
                corner[row] = first + second;
                roundoff[row] = 8 * epsilon * (std::abs(first) + std::abs(second)) + 8 * tiny;
            }
            corner_bound = std::fmax(corner_bound, std::hypot(corner[0], corner[1]) +
                                                       std::hypot(roundoff[0], roundoff[1]));
        }
    }
    query->support_radius = std::nextafter(
        std::fmax(query->support_radius, corner_bound) * (1 + 8 * epsilon) + 8 * tiny,
        std::numeric_limits<double>::infinity());
    query->major_sigma =
        std::nextafter(query->support_radius / 4, std::numeric_limits<double>::infinity());
    query->spherical_area_bound = query->support_radius * query->support_radius / 6;
    const double area_error = query->spherical_area_bound + query->arithmetic_area_bound +
                              query->spherical_area_bound * query->arithmetic_area_bound;
    if (!std::isfinite(area_error) || area_error > geometry_area_budget) {
        return std::unexpected(PointResponseFailure::NeedsRefinement);
    }
    return RestrictedAffinePointResponse{*query, lower, upper, centre, scale};
}

}  // namespace sirius::core
