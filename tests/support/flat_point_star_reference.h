#pragma once

#include <array>
#include <cmath>
#include <numbers>
#include <stdexcept>

namespace sirius::test::flat_point_star_reference {

using Coordinate = std::array<long double, 2>;

// Independent square-film gnomonic geometry for an unboosted pinhole at a
// fixed flat-space event. Angles and the tangent coefficient are the measured
// public lens coefficients; no camera projection, tetrad, beam or density
// routine supplies an expected answer. The local source components are in the
// analytic radial/polar/azimuthal orthonormal frame. Roll is zero in this corpus.
struct Geometry {
    int edge;
    long double tangent;
    long double yaw;
    long double pitch;
    std::array<long double, 3> local_source;

    [[nodiscard]] long double PixelScale() const { return 2 * tangent / edge; }

    [[nodiscard]] Coordinate Image() const {
        // Undo yaw and then pitch on (right, -up, forward), independently of
        // the product's forward film projection and least-aligned-axis basis.
        const auto [radial, polar, azimuthal] = local_source;
        const long double x = azimuthal * std::cos(yaw) - radial * std::sin(yaw);
        const long double after_yaw_z = azimuthal * std::sin(yaw) + radial * std::cos(yaw);
        const long double y = -polar * std::cos(pitch) + after_yaw_z * std::sin(pitch);
        const long double z = polar * std::sin(pitch) + after_yaw_z * std::cos(pitch);
        if (!(z < 0)) throw std::runtime_error("reference star is outside the forward film");
        return {edge / 2.L - x / z / PixelScale(), edge / 2.L + y / z / PixelScale()};
    }

    [[nodiscard]] Coordinate Plane(Coordinate film) const {
        return {PixelScale() * (film[0] - edge / 2.L), -PixelScale() * (film[1] - edge / 2.L)};
    }

    [[nodiscard]] long double AngularAreaDensity(Coordinate film) const {
        const auto plane = Plane(film);
        return PixelScale() * PixelScale() /
               std::pow(1 + plane[0] * plane[0] + plane[1] * plane[1], 1.5L);
    }

    [[nodiscard]] long double PixelSolidAngle(Coordinate centre) const {
        const auto plane = Plane(centre);
        const long double half_pixel = PixelScale() / 2;
        const auto primitive = [](long double u, long double v) {
            return std::atan2(u * v, std::sqrt(1 + u * u + v * v));
        };
        return primitive(plane[0] + half_pixel, plane[1] + half_pixel) -
               primitive(plane[0] - half_pixel, plane[1] + half_pixel) -
               primitive(plane[0] + half_pixel, plane[1] - half_pixel) +
               primitive(plane[0] - half_pixel, plane[1] - half_pixel);
    }

    struct Response {
        long double density;
        long double standard_radius;
    };

    [[nodiscard]] Response OriginalResponse(Coordinate image, Coordinate centre,
                                            long double angular_sigma) const {
        const auto plane = Plane(centre);
        const Coordinate displacement{PixelScale() * (image[0] - centre[0]),
                                      -PixelScale() * (image[1] - centre[1])};
        const long double squared_norm = 1 + plane[0] * plane[0] + plane[1] * plane[1];
        const long double longitudinal = plane[0] * displacement[0] + plane[1] * displacement[1];
        // D^T D is the angular Gram matrix of normalize(u,v,-1).
        // z = P(centre)*(image-centre)/sigma, so its squared norm needs
        // neither the product's tangent basis nor its matrix inverse.
        const long double radius_squared =
            (displacement[0] * displacement[0] + displacement[1] * displacement[1] -
             longitudinal * longitudinal / squared_norm) /
            squared_norm / (angular_sigma * angular_sigma);
        const long double radius = std::sqrt(radius_squared);
        if (radius > 4) return {0, radius};
        // The source derivative is evaluated at the actual image; each
        // original L=sigma*P(centre)^-1 is preserved. In flat space the sky
        // map is a rotation, giving det(G(image)*L)=sigma^2*A(image)/A(centre).
        const long double determinant =
            angular_sigma * angular_sigma * AngularAreaDensity(image) / AngularAreaDensity(centre);
        const long double density =
            std::exp(-radius_squared / 2) /
            (2 * std::numbers::pi_v<long double> * -std::expm1(-8.L) * determinant);
        return {density, radius};
    }
};

}  // namespace sirius::test::flat_point_star_reference
