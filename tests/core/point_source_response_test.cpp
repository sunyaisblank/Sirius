#include "sirius/core/point_source_response.h"

#include "sirius/core/camera.h"
#include "sirius/core/camera_sampling.h"
#include "sirius/core/source_sky_map.h"
#include "sirius/core/starfield.h"

#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <limits>
#include <numbers>
#include <vector>

namespace sirius::core {
namespace {

constexpr AngularMatrix2 kFilm{{{0.001, 0.0}, {0.0, -0.002}}};
struct LensCase {
    AngularMatrix2 map;
    double absolute_determinant;
};
constexpr std::array<LensCase, 6> kLenses{{
    {{{{1, 0}, {0, 1}}}, 1},
    {{{{-1, 0}, {0, 1}}}, 1},
    {{{{2, 0}, {0, 3}}}, 6},
    {{{{0.25, 0}, {0, 0.5}}}, 0.125},
    {{{{1, 2}, {0, 1}}}, 1},
    {{{{0, -2}, {3, 0}}}, 6},
}};

// Forward source offsets for these independent affine lens fixtures. The
// product response builds its inverse internally; the oracle never calls it.
std::array<double, 2> Forward(const AngularMatrix2& lens, double x, double y) {
    const double ax = 0.001 * x;
    const double ay = -0.002 * y;
    return {lens[0][0] * ax + lens[0][1] * ay, lens[1][0] * ax + lens[1][1] * ay};
}

TEST(PointSourceResponseTest, IntegratedDetectorFluxFollowsSignedAffineLensMeasure) {
    constexpr int radial_cells = 4096;
    constexpr int angles = 16;
    constexpr double dr = 4.0 / radial_cells;
    constexpr double dphi = 2 * std::numbers::pi / angles;
    constexpr double pixel_solid_angle = 0.000002;
    for (const auto& lens : kLenses) {
        for (double sigma : {0.3, 1.0}) {
            const auto response = MakeAffinePointResponse(lens.map, kFilm, sigma, 0.0025);
            ASSERT_TRUE(response.has_value());
            double integral = 0.0;
            for (int radial = 0; radial < radial_cells; ++radial) {
                const double radius = (radial + 0.5) * dr;
                for (int angle = 0; angle < angles; ++angle) {
                    const double phi = (angle + 0.5) * dphi;
                    const auto source = Forward(lens.map, sigma * radius * std::cos(phi),
                                                sigma * radius * std::sin(phi));
                    const auto density = response->Density(source);
                    ASSERT_TRUE(density.has_value());
                    integral += *density * pixel_solid_angle * sigma * sigma * radius * dr * dphi;
                }
            }
            // Catalogue source flux=3.7; the received flux includes parity only
            // through its absolute area ratio. No additional pixel-area factor.
            EXPECT_NEAR(3.7 * integral, 3.7 / lens.absolute_determinant,
                        3.7e-7 / lens.absolute_determinant);
        }
    }
}

TEST(PointSourceResponseTest, UnresolvedImageShapeStaysInDetectorCoordinates) {
    for (const auto& lens : kLenses) {
        const auto response = MakeAffinePointResponse(lens.map, kFilm, 1.0, 0.0025);
        ASSERT_TRUE(response.has_value());
        const auto center = response->Density({0, 0});
        ASSERT_TRUE(center.has_value());
        for (const std::array<double, 2> offset :
             {std::array{0.25, -0.75}, std::array{1.5, 0.2}, std::array{-0.6, 2.1}}) {
            const auto density = response->Density(Forward(lens.map, offset[0], offset[1]));
            ASSERT_TRUE(density.has_value());
            const double radial_profile =
                std::exp(-0.5 * (offset[0] * offset[0] + offset[1] * offset[1]));
            EXPECT_NEAR(*density / *center, radial_profile, 2e-15);
        }
        const auto outside = response->Density(Forward(lens.map, 0, 4.01));
        ASSERT_TRUE(outside.has_value());
        EXPECT_EQ(*outside, 0.0);
    }
}

TEST(PointSourceResponseTest, EllipticalSupportRejectsInsideMajorCircleOutsideEllipse) {
    const AngularMatrix2 lens{{{4, 0}, {0, 0.25}}};
    const auto response = MakeAffinePointResponse(lens, kFilm, 1, 0.0025);
    ASSERT_TRUE(response.has_value());
    ASSERT_LT(0.003, response->support_radius);
    const auto outside_minor = response->Density({0, 0.003});
    ASSERT_TRUE(outside_minor.has_value());
    EXPECT_EQ(*outside_minor, 0);
    const auto inside_major = response->Density({0.003, 0});
    ASSERT_TRUE(inside_major.has_value());
    EXPECT_GT(*inside_major, 0);
}

TEST(PointSourceResponseTest, SingularAndBroadMapsRequestRefinementWithoutFloors) {
    for (const AngularMatrix2 lens :
         {AngularMatrix2{{{0, 0}, {0, 1}}}, AngularMatrix2{{{1, 2}, {2, 4}}},
          AngularMatrix2{{{500, 0}, {0, 1}}}}) {
        const auto response = MakeAffinePointResponse(lens, kFilm, 1, 0.0025);
        ASSERT_FALSE(response.has_value());
        EXPECT_EQ(response.error(), PointResponseFailure::NeedsRefinement);
    }
    const AngularMatrix2 identity{{{1, 0}, {0, 1}}};
    const auto response = MakeAffinePointResponse(identity, kFilm, 1, 0.0025);
    ASSERT_TRUE(response.has_value());
    const double exact_spherical_area_error =
        1 - std::sin(response->support_radius) / response->support_radius;
    EXPECT_GE(response->spherical_area_bound, exact_spherical_area_error);
    EXPECT_LT(response->spherical_area_bound, 0.0025);
}

TEST(PointSourceResponseTest, InvalidArithmeticIsDistinctFromEmptyContribution) {
    const AngularMatrix2 identity{{{1, 0}, {0, 1}}};
    for (double sigma : {0.0, -1.0, std::numeric_limits<double>::infinity(),
                         std::numeric_limits<double>::quiet_NaN()}) {
        const auto response = MakeAffinePointResponse(identity, kFilm, sigma, 0.0025);
        ASSERT_FALSE(response.has_value());
        EXPECT_EQ(response.error(), PointResponseFailure::InvalidInput);
    }
    const auto response = MakeAffinePointResponse(identity, kFilm, 1, 0.0025);
    ASSERT_TRUE(response.has_value());
    const auto invalid = response->Density({std::numeric_limits<double>::quiet_NaN(), 0});
    ASSERT_FALSE(invalid.has_value());
    EXPECT_EQ(invalid.error(), PointResponseFailure::InvalidInput);
}

TEST(PointSourceResponseTest, CatalogueAccumulationPreservesFluxColorAndExactSupport) {
    StarfieldConfig config;
    config.brightness_scale = 7.0f;
    StarfieldGenerator generator(config);
    const StarEntry center{0, 0, 1, 10, 0, 0.65f, 5800, 0};
    // The second star is inside the major circle but outside the narrow axis
    // of the exact ellipse. It must contribute nothing despite its large flux.
    const StarEntry excluded{0,
                             static_cast<float>(std::sin(0.003)),
                             static_cast<float>(std::cos(0.003)),
                             10,
                             -5,
                             0.65f,
                             5800,
                             0};
    const StarfieldSpatialIndex center_index({center});
    const StarfieldSpatialIndex both_index({center, excluded});
    const AngularMatrix2 lens{{{4, 0}, {0, 0.25}}};
    const auto response = MakeAffinePointResponse(lens, kFilm, 1, 0.0025);
    ASSERT_TRUE(response.has_value());
    const auto actual = generator.AccumulateThroughResponse({0, 0, 1}, *response, both_index, 1.0);
    const auto only_center =
        generator.AccumulateThroughResponse({0, 0, 1}, *response, center_index, 1.0);
    ASSERT_TRUE(actual.has_value());
    ASSERT_TRUE(only_center.has_value());
    EXPECT_EQ(*actual, *only_center);
    // The separately verified continuous Planck transfer defines reference
    // colour; this test isolates actual density, index support and flux units.
    const auto color = spectral::TransferPointSourceBand(center.temperature_K, 1, 1, 1);
    ASSERT_TRUE(color);
    // Known |det(JP)|=2e-6 and unit film sigma; independent on-axis integral.
    const double intensity = 7.0 / (2 * std::numbers::pi * 2e-6 * (1 - std::exp(-8.0)));
    for (int channel = 0; channel < 3; ++channel) {
        EXPECT_NEAR((*actual)[channel], (*color)[channel] * intensity, intensity * 2e-15);
    }
    const auto invalid = generator.AccumulateThroughResponse({0, 0, 0}, *response, both_index, 1.0);
    ASSERT_FALSE(invalid.has_value());
    EXPECT_EQ(invalid.error(), PointResponseFailure::InvalidInput);
}

// Integration oracle uses the independently verified band helper only for
// colour. Geometry, apparent flux and candidate support are calculated here
// from stored catalogue directions and the exact affine Gaussian definition.
TEST(PointSourceResponseTest, CatalogueBandTransferComposesDensityAndMovingFrequency) {
    StarfieldConfig config;
    config.brightness_scale = 7;
    StarfieldGenerator generator(config);
    const std::vector<StarEntry> stars{
        {0, 0, 1, 10, 0, .2f, 3000, 0},
        {static_cast<float>(std::sin(.0003)), 0, static_cast<float>(std::cos(.0003)), 100, 2.5f,
         .4f, 6500, 0},
        {0, static_cast<float>(std::sin(.0007)), static_cast<float>(std::cos(.0007)), 10000, -2.5f,
         -.1f, 20000, 0},
        {static_cast<float>(std::sin(.012)), 0, static_cast<float>(std::cos(.012)), 25, -5, .5f,
         5800, 0}};
    const StarfieldSpatialIndex index(stars);
    const std::array<double, 4> reference_flux{1, static_cast<double>(.1f), 10, 100};
    constexpr double beta = .6, gamma = 1.25;
    for (double cosine : {-1., 0., 1.}) {
        const double g = gamma * (1 + beta * cosine);  // .5,1.25,2.
        for (double parity : {-1., 1.}) {
            const AngularMatrix2 lens{{{parity * g, 0}, {0, g}}};
            const auto response = MakeAffinePointResponse(lens, kFilm, 1, .0025);
            ASSERT_TRUE(response);
            const auto actual = generator.AccumulateThroughResponse({0, 0, 1}, *response, index, g);
            ASSERT_TRUE(actual);
            std::array<double, 3> expected{};
            for (std::size_t i = 0; i < stars.size(); ++i) {
                const auto& star = stars[i];
                const double transverse = std::hypot(star.direction_x, star.direction_y);
                const double angle = std::atan2(transverse, static_cast<double>(star.direction_z));
                const double x = transverse > 0 ? angle * star.direction_x / transverse : 0;
                const double y = transverse > 0 ? angle * star.direction_y / transverse : 0;
                const double squared = std::pow(x / (g * .001), 2) + std::pow(y / (g * .002), 2);
                const double density =
                    squared > 16 ? 0
                                 : std::exp(-.5 * squared) /
                                       (2 * std::numbers::pi * g * g * 2e-6 * (1 - std::exp(-8.)));
                const auto band = spectral::TransferPointSourceBand(star.temperature_K, g, 1, 1);
                ASSERT_TRUE(band);
                for (int c = 0; c < 3; ++c)
                    expected[c] += (*band)[c] * reference_flux[i] * 7 * density;
            }
            for (int c = 0; c < 3; ++c)
                EXPECT_NEAR((*actual)[c], expected[c], std::abs(expected[c]) * 5e-14);
        }
    }
}

TEST(PointSourceResponseTest, CatalogueBandTransferKeepsFluxAndDensityLinear) {
    StarfieldConfig config;
    config.brightness_scale = 7;
    StarfieldGenerator generator(config);
    const std::vector<StarEntry> stars{{0, 0, 1, 10, 0, .2f, 3000, 0},
                                       {0, 0, 1, 10000, 2.5f, .4f, 6500, 0}};
    const StarfieldSpatialIndex index(stars);
    const AngularMatrix2 lens{{{2, 0}, {0, .5}}};
    const AngularMatrix2 wider_film{{{.002, 0}, {0, -.004}}};
    const auto response = MakeAffinePointResponse(lens, kFilm, 1, .0025);
    const auto wider = MakeAffinePointResponse(lens, wider_film, 1, .0025);
    ASSERT_TRUE(response);
    ASSERT_TRUE(wider);
    for (double g : {.5, 1., 1.25, 2.}) {
        const auto original = generator.AccumulateThroughResponse({0, 0, 1}, *response, index, g);
        const auto spread = generator.AccumulateThroughResponse({0, 0, 1}, *wider, index, g);
        ASSERT_TRUE(original);
        ASSERT_TRUE(spread);
        auto brighter_config = config;
        brighter_config.brightness_scale = 14;
        const StarfieldGenerator brighter(brighter_config);
        const auto doubled = brighter.AccumulateThroughResponse({0, 0, 1}, *response, index, g);
        ASSERT_TRUE(doubled);
        auto other_distance = stars;
        for (auto& star : other_distance) star.distance_pc *= 4;
        const StarfieldSpatialIndex relocated(other_distance);
        const auto same = generator.AccumulateThroughResponse({0, 0, 1}, *response, relocated, g);
        ASSERT_TRUE(same);
        EXPECT_EQ(*same, *original);  // Apparent magnitude already contains distance.
        for (int c = 0; c < 3; ++c) {
            EXPECT_DOUBLE_EQ((*doubled)[c], 2 * (*original)[c]);
            EXPECT_DOUBLE_EQ((*spread)[c], .25 * (*original)[c]);
        }
    }
}

TEST(PointSourceResponseTest, CatalogueBandTransferDeclinesInvalidFrequencyAndPartialFailure) {
    StarfieldConfig config;
    config.brightness_scale = 1;
    const StarfieldGenerator generator(config);
    const StarEntry cold{0, 0, 1, 10, 0, .2f, 100, 0};
    const StarEntry hot{0, 0, 1, 10, 0, .2f, 1e6f, 0};
    const StarfieldSpatialIndex empty(std::vector<StarEntry>{});
    const StarfieldSpatialIndex cold_only({cold});
    const StarfieldSpatialIndex both({cold, hot});
    const AngularMatrix2 identity{{{1, 0}, {0, 1}}};
    const auto response = MakeAffinePointResponse(identity, kFilm, 1, .0025);
    ASSERT_TRUE(response);
    auto zero_config = config;
    zero_config.brightness_scale = 0;
    const StarfieldGenerator zero_flux(zero_config);
    for (double invalid : {0., -1., std::numeric_limits<double>::infinity(),
                           std::numeric_limits<double>::quiet_NaN()}) {
        for (const auto* catalogue : {&empty, &both}) {
            const auto invalid_result =
                generator.AccumulateThroughResponse({0, 0, 1}, *response, *catalogue, invalid);
            ASSERT_FALSE(invalid_result);
            EXPECT_EQ(invalid_result.error(), PointResponseFailure::InvalidInput);
            const auto also_invalid =
                zero_flux.AccumulateThroughResponse({0, 0, 1}, *response, *catalogue, invalid);
            ASSERT_FALSE(also_invalid);
            EXPECT_EQ(also_invalid.error(), PointResponseFailure::InvalidInput);
        }
    }
    const auto no_stars = generator.AccumulateThroughResponse({0, 0, 1}, *response, empty, 2);
    ASSERT_TRUE(no_stars);
    EXPECT_EQ(*no_stars, (std::array<double, 3>{0, 0, 0}));
    // The first star transfers finitely; the second overflows its shifted
    // Planck spectrum. Neither its omission nor publication of the prefix is valid.
    const auto prefix = generator.AccumulateThroughResponse({0, 0, 1}, *response, cold_only, 1e292);
    ASSERT_TRUE(prefix);
    EXPECT_GT(std::max({(*prefix)[0], (*prefix)[1], (*prefix)[2]}), 0);
    auto excluded_hot = hot;
    excluded_hot.direction_x = static_cast<float>(std::sin(.006));
    excluded_hot.direction_z = static_cast<float>(std::cos(.006));
    const StarfieldSpatialIndex outside_ellipse({cold, excluded_hot});
    const auto irrelevant =
        generator.AccumulateThroughResponse({0, 0, 1}, *response, outside_ellipse, 1e292);
    ASSERT_TRUE(irrelevant);
    EXPECT_EQ(*irrelevant, *prefix);
    const auto failed = generator.AccumulateThroughResponse({0, 0, 1}, *response, both, 1e292);
    ASSERT_FALSE(failed);
    EXPECT_EQ(failed.error(), PointResponseFailure::Arithmetic);
    auto invalid_density = *response;
    invalid_density.density_scale = std::numeric_limits<double>::infinity();
    const auto density_error =
        generator.AccumulateThroughResponse({0, 0, 1}, invalid_density, both, 1);
    ASSERT_FALSE(density_error);
    EXPECT_EQ(density_error.error(), PointResponseFailure::Arithmetic);
    const auto recovered = generator.AccumulateThroughResponse({0, 0, 1}, *response, both, 1);
    ASSERT_TRUE(recovered);
}

TEST(PointSourceResponseTest, CompositionCancellationAndUnrepresentedParityDecline) {
    // The exact determinant of the product of these represented inputs is
    // 2e16 * 2^-52. Rounded ordinary dot products give determinant 4 instead,
    // a 9.9% error even though the final B is small and apparently well conditioned.
    const AngularMatrix2 lens{{{1e16, 1e16}, {1e16, 1e16 + 2}}};
    const AngularMatrix2 film{{{1, -1}, {-1, std::nextafter(1.0, 2.0)}}};
    const auto cancelled = MakeAffinePointResponse(lens, film, 0.001, 0.0025);
    ASSERT_FALSE(cancelled.has_value());
    EXPECT_EQ(cancelled.error(), PointResponseFailure::NeedsRefinement);
    const AngularMatrix2 tiny_lens{{{1e-200, 0}, {0, 1e-200}}};
    const AngularMatrix2 huge_film{{{1e200, 0}, {0, 1e200}}};
    const auto unrepresented_parity = MakeAffinePointResponse(tiny_lens, huge_film, 0.001, 0.0025);
    ASSERT_FALSE(unrepresented_parity.has_value());
    EXPECT_EQ(unrepresented_parity.error(), PointResponseFailure::Arithmetic);
}

TEST(PointSourceResponseTest, MovingThinLensAndSourceMapComposeIntoReceivedFlux) {
    CameraConfig config;
    config.width = 1920;
    config.height = 1080;
    config.r = 30;
    config.theta = 80 * std::numbers::pi / 180;
    config.fov = 60;
    config.focus_distance = 30;
    config.beta_x = 0.1;
    config.beta_y = 0.02;
    config.beta_z = -0.01;
    ThinLensCamera camera(config);
    relativity::ObserverFrame rest;
    rest.time(0) = 1;
    for (int axis = 0; axis < 3; ++axis) rest.spatial[axis](axis + 1) = 1;
    const std::array<double, 3> velocity{-config.beta_x, -config.beta_y, config.beta_z};
    const auto moving = relativity::BoostObserverFrame(rest, velocity);
    ASSERT_TRUE(moving.has_value());
    Metric4d metric;
    metric(0, 0) = -1;
    for (int axis = 1; axis < 4; ++axis) metric(axis, axis) = 1;
    Tensor<Dual<double>, 4, 4, 4> derivatives;
    unsigned samples = 0;
    ForEachCameraSample(4, [&](const CameraSample& sample) {
        ++samples;
        const auto ray = camera.GenerateRayForObserver(169, 229, sample.image_u, sample.image_v,
                                                       sample.pupil_u, sample.pupil_v);
        const auto film = camera.FilmDifferentialForObserver(
            169, 229, sample.image_u, sample.image_v, sample.pupil_u, sample.pupil_v);
        ASSERT_TRUE(film.has_value());
        std::array<double, 3> direction{ray.direction(1), ray.direction(2), ray.direction(3)};
        const double norm = std::hypot(direction[0], direction[1], direction[2]);
        for (double& value : direction) value /= norm;
        const auto tangent = relativity::PastDirectedCameraRay(*moving, direction);
        const auto screen = relativity::ObserverScreenBasis(*moving, direction);
        ASSERT_TRUE(tangent.has_value());
        ASSERT_TRUE(screen.has_value());
        const auto source = relativity::MeasureSourceSkyAngularMap(
            metric, metric, derivatives, *tangent, {Vec4{}, Vec4{}}, *screen, 1);
        ASSERT_TRUE(source.has_value());
        const auto response =
            MakeAffinePointResponse(source->jacobian, film->angular_jacobian, 1, 0.0025);
        ASSERT_TRUE(response.has_value());
        // Independent Lorentz solid-angle law for a past ray: the source
        // frequency is gamma(1-beta.n), and dOmega_source=dOmega_camera/omega^2.
        double beta_squared = 0, beta_dot_n = 0;
        for (int axis = 0; axis < 3; ++axis) {
            beta_squared += velocity[axis] * velocity[axis];
            beta_dot_n += velocity[axis] * direction[axis];
        }
        const double omega = (1 - beta_dot_n) / std::sqrt(1 - beta_squared);
        EXPECT_NEAR(std::abs(source->determinant), 1 / (omega * omega), 2e-15);
        double integral = 0;
        constexpr int radial_cells = 1024, angles = 8;
        constexpr double dr = 4.0 / radial_cells, dphi = 2 * std::numbers::pi / angles;
        for (int radial = 0; radial < radial_cells; ++radial) {
            const double radius = (radial + 0.5) * dr;
            for (int angle = 0; angle < angles; ++angle) {
                const double phi = (angle + 0.5) * dphi;
                const std::array<double, 2> pixel{radius * std::cos(phi), radius * std::sin(phi)};
                std::array<double, 2> launch{}, offset{};
                for (int row = 0; row < 2; ++row)
                    for (int col = 0; col < 2; ++col)
                        launch[row] += film->angular_jacobian[row][col] * pixel[col];
                for (int row = 0; row < 2; ++row)
                    for (int col = 0; col < 2; ++col)
                        offset[row] += source->jacobian[row][col] * launch[col];
                const auto density = response->Density(offset);
                ASSERT_TRUE(density.has_value());
                integral += *density * film->solid_angle_density * radius * dr * dphi;
            }
        }
        EXPECT_NEAR(integral, omega * omega, 1e-6);
    });
    EXPECT_EQ(samples, 4u);
}

std::array<double, 2> RestrictedSourceOffset(const AngularMatrix2& map, double x, double y) {
    return {map[0][0] * x + map[0][1] * y, map[1][0] * x + map[1][1] * y};
}

TEST(PointSourceResponseTest, RestrictedCellsPreserveOriginalGaussianAndNormalization) {
    constexpr AngularMatrix2 map{{{1.0 / 64, 1.0 / 256}, {1.0 / 512, -1.0 / 128}}};
    constexpr double determinant = std::abs(map[0][0] * map[1][1] - map[0][1] * map[1][0]);
    std::vector<RestrictedAffinePointResponse> cells;
    for (int y = 0; y < 2; ++y) {
        for (int x = 0; x < 2; ++x) {
            const auto cell = MakeRestrictedAffinePointResponse(map, {-4.0 + 4 * x, -4.0 + 4 * y},
                                                                {4.0 * x, 4.0 * y}, 0.01);
            ASSERT_TRUE(cell);
            cells.push_back(*cell);
        }
    }
    const auto sum_at = [&](const std::array<double, 2>& z) {
        double sum = 0.0;
        for (const auto& cell : cells) {
            const auto offset =
                RestrictedSourceOffset(map, z[0] - cell.centre[0], z[1] - cell.centre[1]);
            const auto density = cell.Density(offset);
            EXPECT_TRUE(density);
            if (density) sum += *density;
        }
        return sum;
    };
    for (const auto& z : std::array<std::array<double, 2>, 5>{
             {{-3.5, -0.5}, {-1.5, 1.5}, {0.25, -0.75}, {2.5, 2.0}, {3.5, 3.5}}}) {
        const double radius = std::hypot(z[0], z[1]);
        const double expected =
            radius > 4.0 ? 0.0
                         : std::exp(-0.5 * radius * radius) /
                               (2 * std::numbers::pi * determinant * (-std::expm1(-8.0)));
        EXPECT_NEAR(sum_at(z), expected, 8e-14 * std::max(expected, 1.0));
    }
    // Independent polar quadrature in original coordinates. The source-area
    // Jacobian is applied once; changing child width cannot change this flux.
    constexpr int radii = 1024, angles = 8;
    constexpr double dr = 4.0 / radii, dphi = 2 * std::numbers::pi / angles;
    double integral = 0.0;
    for (int i = 0; i < radii; ++i) {
        const double radius = (i + 0.5) * dr;
        for (int j = 0; j < angles; ++j) {
            const double phi = (j + 0.5) * dphi;
            integral += sum_at({radius * std::cos(phi), radius * std::sin(phi)}) * determinant *
                        radius * dr * dphi;
        }
    }
    EXPECT_NEAR(integral, 1.0, 1e-6);
}

TEST(PointSourceResponseTest, RestrictedCellsOwnSharedEdgesAndKeepParentCutoff) {
    constexpr AngularMatrix2 map{{{1.0 / 64, 0.0}, {0.0, 1.0 / 64}}};
    const auto left = MakeRestrictedAffinePointResponse(map, {-4, -4}, {0, 4}, 0.01);
    const auto right = MakeRestrictedAffinePointResponse(map, {0, -4}, {4, 4}, 0.01);
    ASSERT_TRUE(left && right);
    const auto at = [](const RestrictedAffinePointResponse& cell, double x, double y) {
        return cell.DensityAtOriginalRoot({x, y});
    };
    const auto left_edge = at(*left, 0, 0), right_edge = at(*right, 0, 0);
    ASSERT_TRUE(left_edge && right_edge);
    EXPECT_EQ(*left_edge, 0.0);
    EXPECT_GT(*right_edge, 0.0);
    const auto cutoff = at(*right, 4, 0), corner = at(*right, 3.5, 3.5);
    ASSERT_TRUE(cutoff && corner);
    EXPECT_GT(*cutoff, 0.0);
    EXPECT_EQ(*corner, 0.0);
    const auto beyond = at(*right, 4.001, 0);
    ASSERT_TRUE(beyond);
    EXPECT_EQ(*beyond, 0.0);

    // Independent affine inverses round to opposite sides of this true shared
    // edge. They must report uncertainty, not both silently discard the image.
    constexpr AngularMatrix2 skew{{{1.0 / 64, 1.0 / 256}, {1.0 / 512, -1.0 / 128}}};
    const auto skew_left = MakeRestrictedAffinePointResponse(skew, {-4, -4}, {0, 4}, 0.01);
    const auto skew_right = MakeRestrictedAffinePointResponse(skew, {0, -4}, {4, 4}, 0.01);
    ASSERT_TRUE(skew_left && skew_right);
    for (const auto* cell : {&*skew_left, &*skew_right}) {
        const auto offset =
            RestrictedSourceOffset(skew, -cell->centre[0], -3.89922 - cell->centre[1]);
        const auto uncertain = cell->Density(offset);
        ASSERT_FALSE(uncertain);
        EXPECT_EQ(uncertain.error(), PointResponseFailure::NeedsRefinement);
    }
    const auto owned_left = skew_left->DensityAtOriginalRoot({0, -3.89922});
    const auto owned_right = skew_right->DensityAtOriginalRoot({0, -3.89922});
    ASSERT_TRUE(owned_left && owned_right);
    EXPECT_EQ(*owned_left, 0.0);
    EXPECT_GT(*owned_right, 0.0);

    // A rounded singular-value radius formerly excluded this owned corner by
    // one ULP. Query-only outward padding must include it without changing W.
    constexpr AngularMatrix2 rotated{{{0.0022468547126300972, 0.008799056649877072},
                                      {-0.008799056649877072, 0.0022468547126300972}}};
    const auto rotated_cell = MakeRestrictedAffinePointResponse(rotated, {1, 1}, {1.5, 1.5}, 0.01);
    ASSERT_TRUE(rotated_cell);
    const auto corner_offset = RestrictedSourceOffset(rotated, -0.25, -0.25);
    EXPECT_GE(rotated_cell->query.support_radius, std::hypot(corner_offset[0], corner_offset[1]));
    const auto corner_density = rotated_cell->DensityAtOriginalRoot({1, 1});
    ASSERT_TRUE(corner_density);
    EXPECT_GT(*corner_density, 0.0);

    double three_ulps = 1.0;
    for (int i = 0; i < 3; ++i) three_ulps = std::nextafter(three_ulps, 2.0);
    for (const auto& bounds :
         std::array<std::array<double, 2>, 2>{{{1.0, three_ulps}, {0.1, 0.23456789}}}) {
        const auto narrow = MakeRestrictedAffinePointResponse(map, {bounds[0], bounds[0]},
                                                              {bounds[1], bounds[1]}, 0.01);
        ASSERT_TRUE(narrow);
        for (double x : bounds) {
            for (double y : bounds) {
                const auto offset =
                    RestrictedSourceOffset(map, x - narrow->centre[0], y - narrow->centre[1]);
                EXPECT_GE(narrow->query.support_radius, std::hypot(offset[0], offset[1]));
            }
        }
    }
}

TEST(PointSourceResponseTest, RestrictedCellsDeclineInvalidOrUnrepresentedGeometry) {
    constexpr AngularMatrix2 map{{{1.0 / 64, 0.0}, {0.0, 1.0 / 64}}};
    EXPECT_FALSE(MakeRestrictedAffinePointResponse(map, {-5, -4}, {4, 4}, 0.01));
    EXPECT_FALSE(MakeRestrictedAffinePointResponse(map, {0, -4}, {0, 4}, 0.01));
    const double next = std::nextafter(1.0, 2.0);
    const auto lost = MakeRestrictedAffinePointResponse(map, {1, -1}, {next, 1}, 0.01);
    ASSERT_FALSE(lost);
    EXPECT_EQ(lost.error(), PointResponseFailure::NeedsRefinement);
    constexpr AngularMatrix2 singular{{{1.0 / 64, 0.0}, {0.0, 0.0}}};
    EXPECT_FALSE(MakeRestrictedAffinePointResponse(singular, {-4, -4}, {4, 4}, 0.01));
    const auto cell = MakeRestrictedAffinePointResponse(map, {-4, -4}, {4, 4}, 0.01);
    ASSERT_TRUE(cell);
    const auto invalid = cell->Density({std::numeric_limits<double>::quiet_NaN(), 0});
    ASSERT_FALSE(invalid);
    EXPECT_EQ(invalid.error(), PointResponseFailure::InvalidInput);
    for (double bad :
         {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity()}) {
        for (const auto& root : std::array<std::array<double, 2>, 2>{{{5, bad}, {bad, 5}}}) {
            const auto bad_root = cell->DensityAtOriginalRoot(root);
            ASSERT_FALSE(bad_root);
            EXPECT_EQ(bad_root.error(), PointResponseFailure::InvalidInput);
        }
    }
}

}  // namespace
}  // namespace sirius::core
