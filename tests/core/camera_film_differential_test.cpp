#include "sirius/core/camera.h"
#include "sirius/core/camera_launch.h"
#include "sirius/core/camera_sampling.h"
#include "sirius/core/metrics/kerr_schild_family.h"
#include "sirius/core/observer_frame.h"

#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <numbers>

namespace sirius::core::test {
namespace {
using Direction = std::array<double, 3>;

Direction Unit(Direction v) {
    const double length = std::hypot(v[0], v[1], v[2]);
    for (double& value : v) value /= length;
    return v;
}

Direction DirectionOf(const CameraRay& ray) {
    return Unit({ray.direction(1), ray.direction(2), ray.direction(3)});
}

// Independent continuous geometry oracle. The public camera stores float lens
// coefficients, but its arithmetic quantization is not a differentiable lens.
// Here film coordinates, rotations and normalization are evaluated in double.
// A 2e-9 + 2e-6 relative derivative budget covers the nominal float anchor and
// stored rotation coefficients; it is not an integration/photometry tolerance.
Direction GeometricDirection(const ICamera& camera, double film_x, double film_y,
                             const CameraRay& fixed_pupil) {
    const auto& c = camera.GetConfig();
    const double aspect = static_cast<float>(c.width) / c.height;
    const double x = (2.0 * film_x / c.width - 1.0) * aspect;
    const double y = 1.0 - 2.0 * film_y / c.height;
    const double t = std::tan(c.fov * static_cast<float>(std::numbers::pi) / 360.0f);
    if (camera.GetLensType() == LensType::ThinLens) {
        return Unit({-c.focus_distance, fixed_pupil.aperture_up - y * t * c.focus_distance,
                     x * t * c.focus_distance - fixed_pupil.aperture_right});
    }
    if (camera.GetLensType() == LensType::Fisheye) {
        const double radius = std::hypot(x, y);
        const double angle = radius * c.fov * static_cast<float>(std::numbers::pi) / 360.0;
        if (radius == 0) return {-1, 0, 0};
        return {-std::cos(angle), -std::sin(angle) * y / radius, std::sin(angle) * x / radius};
    }
    // Rotate a camera Cartesian vector independently as three plane rotations.
    Direction v{x * t, y * t, -1};
    const auto plane_rotate = [&](unsigned a, unsigned b, float angle) {
        const double cosine = std::cos(angle), sine = std::sin(angle);
        const double old_a = v[a];
        v[a] = cosine * old_a - sine * v[b];
        v[b] = sine * old_a + cosine * v[b];
    };
    plane_rotate(0, 1, c.roll);
    plane_rotate(1, 2, c.pitch);
    plane_rotate(2, 0, c.yaw);
    return Unit({v[2], -v[1], v[0]});
}

void CheckDifferential(const ICamera& camera, int x, int y, const CameraSample& sample) {
    SCOPED_TRACE(::testing::Message() << camera.GetName() << ' ' << x << ',' << y);
    const auto ray = camera.GenerateRayForObserver(x, y, sample.image_u, sample.image_v,
                                                   sample.pupil_u, sample.pupil_v);
    const auto map = camera.FilmDifferentialForObserver(x, y, sample.image_u, sample.image_v,
                                                        sample.pupil_u, sample.pupil_v);
    ASSERT_TRUE(map);
    const auto n = DirectionOf(ray);
    const auto basis = relativity::MakeCelestialTangentBasis(n);
    ASSERT_TRUE(basis);
    constexpr double h = 0.01;
    for (unsigned column = 0; column < 2; ++column) {
        const double fx = x + static_cast<double>(sample.image_u);
        const double fy = y + static_cast<double>(sample.image_v);
        const auto plus =
            GeometricDirection(camera, fx + (column == 0 ? h : 0), fy + (column == 1 ? h : 0), ray);
        const auto minus =
            GeometricDirection(camera, fx - (column == 0 ? h : 0), fy - (column == 1 ? h : 0), ray);
        double longitudinal = 0;
        for (unsigned i = 0; i < 3; ++i) {
            const double derivative = (plus[i] - minus[i]) / (2 * h);
            EXPECT_NEAR(map->direction_derivative[i][column], derivative,
                        2e-9 + 2e-6 * std::abs(derivative));
            longitudinal += n[i] * map->direction_derivative[i][column];
            const double reconstructed = basis->first[i] * map->angular_jacobian[0][column] +
                                         basis->second[i] * map->angular_jacobian[1][column];
            EXPECT_NEAR(reconstructed, map->direction_derivative[i][column], 2e-17);
        }
        EXPECT_NEAR(longitudinal, 0, 2e-17);
        // Independently differentiate actual public float rays with a fourth
        // order stencil. Two-pixel spacing keeps float rounding below 8e-8
        // rad/pixel; truncation is negligible at these >=800-pixel resolutions.
        std::array<Direction, 4> values;
        for (unsigned j = 0; j < 4; ++j) {
            constexpr std::array<int, 4> offsets{-4, -2, 2, 4};
            const auto neighbour = camera.GenerateRayForObserver(
                x + (column == 0 ? offsets[j] : 0), y + (column == 1 ? offsets[j] : 0),
                sample.image_u, sample.image_v, sample.pupil_u, sample.pupil_v);
            EXPECT_EQ(neighbour.aperture_up, ray.aperture_up);
            EXPECT_EQ(neighbour.aperture_right, ray.aperture_right);
            values[j] = DirectionOf(neighbour);
        }
        for (unsigned i = 0; i < 3; ++i) {
            const double derivative =
                (values[0][i] - 8 * values[1][i] + 8 * values[2][i] - values[3][i]) / 24;
            EXPECT_NEAR(map->direction_derivative[i][column], derivative, 8e-8);
        }
    }
    EXPECT_GT(map->solid_angle_density, 0);
    EXPECT_EQ(map->solid_angle_density, std::abs(map->signed_solid_angle_density));
}

void ExpectSameRay(const CameraRay& a, const CameraRay& b) {
    for (int i = 0; i < 4; ++i) {
        EXPECT_EQ(a.origin(i), b.origin(i));
        EXPECT_EQ(a.direction(i), b.direction(i));
    }
    EXPECT_EQ(a.beta_forward, b.beta_forward);
    EXPECT_EQ(a.beta_up, b.beta_up);
    EXPECT_EQ(a.beta_right, b.beta_right);
    EXPECT_EQ(a.aperture_up, b.aperture_up);
    EXPECT_EQ(a.aperture_right, b.aperture_right);
    EXPECT_EQ(a.active, b.active);
}
}  // namespace

TEST(CameraFilmDifferential, IndependentGeometryAndActualRayStencils) {
    CameraConfig config;
    config.width = 1200;
    config.height = 800;
    config.yaw = 0.31f;
    config.pitch = -0.27f;
    config.roll = 0.19f;
    for (const auto lens : {LensType::Pinhole, LensType::ThinLens, LensType::Fisheye}) {
        const auto camera = CreateCamera(lens, config);
        for (const auto pixel :
             {std::array{121, 231}, std::array{731, 437}, std::array{990, 681}}) {
            ASSERT_NO_FATAL_FAILURE(CheckDifferential(*camera, pixel[0], pixel[1],
                                                      CameraSample{0.37f, 0.63f, 0.21f, 0.82f}));
        }
    }
}

TEST(CameraFilmDifferential, PerspectiveAndFixedPupilSolidAngleLaws) {
    CameraConfig config;
    const double aspect = static_cast<float>(config.width) / config.height;
    const double t = std::tan(config.fov * static_cast<float>(std::numbers::pi) / 360.0f);
    for (const auto lens : {LensType::Pinhole, LensType::ThinLens}) {
        const auto camera = CreateCamera(lens, config);
        for (const auto pupil : {std::array{0.5f, 0.0f}, std::array{0.19f, 0.93f}}) {
            constexpr int x = 283, y = 729;
            constexpr float u = 0.29f, v = 0.71f;
            const auto ray = camera->GenerateRayForObserver(x, y, u, v, pupil[0], pupil[1]);
            const auto map = camera->FilmDifferentialForObserver(x, y, u, v, pupil[0], pupil[1]);
            ASSERT_TRUE(map);
            const double focus = lens == LensType::ThinLens ? config.focus_distance : 1.0;
            const double px = (2 * (x + double(u)) / config.width - 1) * aspect;
            const double py = 1 - 2 * (y + double(v)) / config.height;
            const double length = std::hypot(focus, py * t * focus - ray.aperture_up,
                                             px * t * focus - ray.aperture_right);
            // Solid angle subtended by an infinitesimal area on a plane:
            // dOmega = (normal distance / distance^3) dA.
            const double expected = 4 * aspect / (config.width * double(config.height)) * t * t *
                                    focus * focus * focus / (length * length * length);
            EXPECT_NEAR(map->solid_angle_density, expected, expected * 3e-6);
            EXPECT_GT(map->signed_solid_angle_density, 0);
        }
    }
    // A centred pupil reduces exactly to the same continuous perspective law.
    PinholeCamera pinhole(config);
    ThinLensCamera lens(config);
    const auto a = pinhole.FilmDifferentialForObserver(127, 349);
    const auto b = lens.FilmDifferentialForObserver(127, 349);
    ASSERT_TRUE(a);
    ASSERT_TRUE(b);
    EXPECT_NEAR(a->solid_angle_density, b->solid_angle_density, a->solid_angle_density * 3e-6);
}

TEST(CameraFilmDifferential, EquidistantCentreLimitAndAntipodalMask) {
    CameraConfig config;
    config.width = config.height = 1000;
    config.fov = 360;
    FisheyeCamera camera(config);
    const double a = config.fov * double(static_cast<float>(std::numbers::pi)) / 360;
    for (const auto pixel : {std::array{500, 500}, std::array{501, 500}, std::array{672, 793}}) {
        const auto map = camera.FilmDifferentialForObserver(pixel[0], pixel[1], 0, 0);
        ASSERT_TRUE(map);
        const double x = 2.0 * pixel[0] / config.width - 1;
        const double y = 1 - 2.0 * pixel[1] / config.height;
        const double radius = std::hypot(x, y);
        const double expected = 4.0 / (config.width * config.height) *
                                (radius == 0 ? a * a : a * std::sin(a * radius) / radius);
        EXPECT_NEAR(map->solid_angle_density, expected, expected * 3e-6);
        EXPECT_GT(map->signed_solid_angle_density, 0);
    }
    const auto centre = camera.FilmDifferentialForObserver(500, 500, 0, 0);
    ASSERT_TRUE(centre);
    EXPECT_NEAR(centre->direction_derivative[2][0], 2 * a / config.width, 1e-15);
    EXPECT_NEAR(centre->direction_derivative[1][1], 2 * a / config.height, 1e-15);
    EXPECT_EQ(centre->direction_derivative[0][0], 0);
    EXPECT_EQ(centre->direction_derivative[0][1], 0);
    EXPECT_TRUE(camera.GenerateRay(0, 500, 0, 0).active);
    EXPECT_FALSE(camera.FilmDifferentialForObserver(0, 500, 0, 0));
    EXPECT_FALSE(camera.GenerateRay(0, 0, 0, 0).active);
    EXPECT_FALSE(camera.FilmDifferentialForObserver(0, 0, 0, 0));
}

TEST(CameraFilmDifferential, PixelUnitsAndIntegratedRectangleSolidAngle) {
    CameraConfig config;
    config.width = 320;
    config.height = 240;
    PinholeCamera camera(config);
    // Composite midpoint integration of the production local density. The
    // O(h^2) quadrature error is bounded separately from float lens arithmetic.
    double integral = 0;
    for (int y = 0; y < config.height; ++y) {
        for (int x = 0; x < config.width; ++x) {
            const auto value = camera.FilmDifferentialForObserver(x, y);
            ASSERT_TRUE(value);
            integral += value->solid_angle_density;
        }
    }
    const double t = std::tan(config.fov * static_cast<float>(std::numbers::pi) / 360.0f);
    const double a = t * (static_cast<float>(config.width) / config.height), b = t;
    const double exact = 4 * std::atan(a * b / std::sqrt(1 + a * a + b * b));
    EXPECT_NEAR(integral, exact, 2e-5);
    // Use binary-exact film divisions for the exact scaling identity. With
    // other dimensions the unchanged nominal float ray can differ by an ULP
    // after compiler constant folding, as the geometric tests allow explicitly.
    config.width = 512;
    config.height = 256;
    PinholeCamera scale_camera(config);
    config.width *= 2;
    config.height *= 2;
    PinholeCamera double_resolution(config);
    const auto coarse = scale_camera.FilmDifferentialForObserver(123, 77, 0.25f, 0.25f);
    const auto fine = double_resolution.FilmDifferentialForObserver(246, 154, 0.5f, 0.5f);
    ASSERT_TRUE(coarse);
    ASSERT_TRUE(fine);
    EXPECT_DOUBLE_EQ(fine->solid_angle_density * 4, coarse->solid_angle_density);
    for (unsigned i = 0; i < 3; ++i) {
        for (unsigned j = 0; j < 2; ++j) {
            EXPECT_DOUBLE_EQ(fine->direction_derivative[i][j] * 2,
                             coarse->direction_derivative[i][j]);
        }
    }
}

TEST(CameraFilmDifferential, PreservesNominalRaysAndExistingOrientationSemantics) {
    CameraConfig config;
    CameraConfig rotated = config;
    rotated.yaw = 0.31f;
    rotated.pitch = -0.27f;
    rotated.roll = 0.19f;
    for (const auto lens : {LensType::Pinhole, LensType::ThinLens, LensType::Fisheye}) {
        const auto camera = CreateCamera(lens, config);
        const auto other = CreateCamera(lens, rotated);
        int samples = 0;
        ForEachCameraSample(4, [&](const CameraSample& sample) {
            ++samples;
            const auto before = camera->GenerateRayForObserver(
                297, 381, sample.image_u, sample.image_v, sample.pupil_u, sample.pupil_v);
            const auto differential = camera->FilmDifferentialForObserver(
                297, 381, sample.image_u, sample.image_v, sample.pupil_u, sample.pupil_v);
            ASSERT_TRUE(differential);
            ExpectSameRay(before,
                          camera->GenerateRayForObserver(297, 381, sample.image_u, sample.image_v,
                                                         sample.pupil_u, sample.pupil_v));
            const auto oriented = other->GenerateRayForObserver(
                297, 381, sample.image_u, sample.image_v, sample.pupil_u, sample.pupil_v);
            const auto oriented_map = other->FilmDifferentialForObserver(
                297, 381, sample.image_u, sample.image_v, sample.pupil_u, sample.pupil_v);
            ASSERT_TRUE(oriented_map);
            if (lens == LensType::Pinhole) {
                EXPECT_NE(oriented.direction(1), before.direction(1));
                EXPECT_NEAR(oriented_map->solid_angle_density, differential->solid_angle_density,
                            differential->solid_angle_density * 3e-6);
            } else {
                ExpectSameRay(before, oriented);
                EXPECT_EQ(oriented_map->direction_derivative, differential->direction_derivative);
                EXPECT_EQ(oriented_map->angular_jacobian, differential->angular_jacobian);
            }
        });
        EXPECT_EQ(samples, 4);
    }
}

TEST(CameraFilmDifferential, GovernedFourSampleLaunchAndSingleObserverBoost) {
    CameraConfig config;
    config.r = 30;
    config.theta = 80 * std::numbers::pi / 180;
    config.focus_distance = 30;
    config.beta_x = 0.1;
    config.beta_y = 0.02;
    config.beta_z = -0.01;
    ThinLensCamera camera(config);
    CameraConfig stationary_config = config;
    stationary_config.beta_x = stationary_config.beta_y = stationary_config.beta_z = 0;
    ThinLensCamera stationary(stationary_config);
    relativity::ObserverFrame rest;
    rest.time(0) = 1;
    for (int i = 0; i < 3; ++i) rest.spatial[i](i + 1) = 1;
    const auto moving =
        relativity::BoostObserverFrame(rest, {-config.beta_x, -config.beta_y, config.beta_z});
    ASSERT_TRUE(moving);
    int launches = 0;
    for (const auto pixel : {std::array{169, 229}, std::array{1167, 586}, std::array{874, 555}}) {
        ForEachCameraSample(4, [&](const CameraSample& sample) {
            ++launches;
            ASSERT_NO_FATAL_FAILURE(CheckDifferential(camera, pixel[0], pixel[1], sample));
            const auto map = camera.FilmDifferentialForObserver(
                pixel[0], pixel[1], sample.image_u, sample.image_v, sample.pupil_u, sample.pupil_v);
            const auto unmoving = stationary.FilmDifferentialForObserver(
                pixel[0], pixel[1], sample.image_u, sample.image_v, sample.pupil_u, sample.pupil_v);
            ASSERT_TRUE(map);
            ASSERT_TRUE(unmoving);
            EXPECT_EQ(map->direction_derivative, unmoving->direction_derivative);
            EXPECT_EQ(map->angular_jacobian, unmoving->angular_jacobian);
            EXPECT_EQ(map->solid_angle_density, unmoving->solid_angle_density);
            const auto ray = camera.GenerateRayForObserver(
                pixel[0], pixel[1], sample.image_u, sample.image_v, sample.pupil_u, sample.pupil_v);
            const auto screen = relativity::ObserverScreenBasis(*moving, DirectionOf(ray));
            ASSERT_TRUE(screen);
            for (unsigned column = 0; column < 2; ++column) {
                Vec4 from_rest;
                for (unsigned i = 0; i < 3; ++i)
                    from_rest += moving->spatial[i] * map->direction_derivative[i][column];
                const Vec4 from_screen = (*screen)[0] * map->angular_jacobian[0][column] +
                                         (*screen)[1] * map->angular_jacobian[1][column];
                for (int component = 0; component < 4; ++component)
                    EXPECT_NEAR(from_rest(component), from_screen(component), 2e-17);
                // Actual boosted null launch FD is independent of screen-J
                // reconstruction and catches a second aberration application.
                constexpr double h = 0.01;
                const double fx = pixel[0] + double(sample.image_u);
                const double fy = pixel[1] + double(sample.image_v);
                const auto plus = relativity::PastDirectedCameraRay(
                    *moving, GeometricDirection(camera, fx + (column == 0 ? h : 0),
                                                fy + (column == 1 ? h : 0), ray));
                const auto minus = relativity::PastDirectedCameraRay(
                    *moving, GeometricDirection(camera, fx - (column == 0 ? h : 0),
                                                fy - (column == 1 ? h : 0), ray));
                ASSERT_TRUE(plus);
                ASSERT_TRUE(minus);
                for (int component = 0; component < 4; ++component) {
                    const double expected = ((*plus)(component) - (*minus)(component)) / (2 * h);
                    EXPECT_NEAR(from_screen(component), expected, 2e-9 + 2e-6 * std::abs(expected));
                }
            }
        });
    }
    EXPECT_EQ(launches, 12);
}

TEST(CameraFilmDifferential, CustomCameraDeclinesUnavailableDifferential) {
    // Existing custom projections need not add a made-up Jacobian to remain
    // concrete implementations of ICamera.
    class CustomCamera final : public ICamera {
      public:
        CameraRay GenerateRay(int x, int y, float u, float v, float pu, float pv) const override {
            return delegate.GenerateRay(x, y, u, v, pu, pv);
        }
        LensType GetLensType() const override { return LensType::Pinhole; }
        const char* GetName() const override { return "Custom"; }
        const CameraConfig& GetConfig() const override { return delegate.GetConfig(); }
        void SetConfig(const CameraConfig& c) override { delegate.SetConfig(c); }

      private:
        PinholeCamera delegate;
    } camera;
    EXPECT_TRUE(camera.GenerateRayForObserver(50, 60).active);
    EXPECT_FALSE(camera.FilmDifferentialForObserver(50, 60));
}
TEST(CameraFilmDifferential, FullLaunchDifferentiatesDisplacedMovingKerrObserver) {
    for (const double spin : {-0.9, 0.0, 0.9}) {
        KerrSchildFamily metric(KerrSchildParams::Kerr(1.0, spin));
        CameraConfig c;
        c.r = 12;
        c.theta = 1.1;
        c.phi = 0.7;
        c.beta_x = 0.31;
        c.beta_y = -0.17;
        c.beta_z = 0.21;
        c.width = 1920;
        c.height = 1080;
        c.focus_distance = 20;
        ThinLensCamera camera(c);
        constexpr double film_x = 704.25, film_y = 367.75;
        const auto projection = camera.ProjectFilmForObserver(film_x, film_y, 0.73f, 0.29f);
        ASSERT_TRUE(projection);
        ASSERT_TRUE(projection->ray.phase_space);
        const auto launch = LaunchCameraRay(metric, spin, projection->ray);
        ASSERT_TRUE(launch);
        Metric4d g;
        Tensor<Dual<double>, 4, 4, 4> dg;
        metric.Evaluate(launch->position, g, dg);
        const auto connection = TensorOps::Christoffel(g, dg);
        EXPECT_NEAR(TensorOps::InnerProduct(launch->tangent, launch->tangent, g), 0, 2e-14);
        EXPECT_NEAR(TensorOps::InnerProduct(launch->tangent, launch->observer.time, g), 1, 2e-14);
        double frame_contribution = 0;
        for (int column = 0; column < 4; ++column) {
            // Five-point differences of the complete launch: film coordinates
            // or physical Cartesian pupil coordinates are varied independently.
            // This reference does not consume the analytic derivative fields.
            constexpr double h = 0.002;
            std::array<CameraLaunch, 4> samples;
            int index = 0;
            for (const double offset : {-2 * h, -h, h, 2 * h}) {
                auto varied = camera.ProjectFilmForObserver(film_x + (column == 0 ? offset : 0),
                                                            film_y + (column == 1 ? offset : 0),
                                                            0.73f, 0.29f);
                ASSERT_TRUE(varied);
                CameraRay ray = varied->ray;
                if (column >= 2) {
                    if (column == 2)
                        ray.aperture_right += offset;
                    else
                        ray.aperture_up += offset;
                    // Recover the fixed focus-plane point from the central
                    // ray. Re-evaluating tanf in a separately optimized context
                    // can change a stored lens coefficient by one float ULP.
                    const double distance =
                        -static_cast<double>(c.focus_distance) / projection->ray.direction(1);
                    std::array<double, 3> q{};
                    for (int axis = 0; axis < 3; ++axis)
                        q[axis] = projection->ray.direction(axis + 1) * distance;
                    if (column == 2)
                        q[2] -= offset;
                    else
                        q[1] += offset;
                    const double norm = std::hypot(q[0], q[1], q[2]);
                    for (int axis = 0; axis < 3; ++axis) ray.direction(axis + 1) = q[axis] / norm;
                }
                ray.phase_space.reset();
                const auto sample = LaunchCameraRay(metric, spin, ray);
                ASSERT_TRUE(sample);
                samples[index++] = *sample;
            }
            Vec4 coordinate = launch->variations[column].derivative;
            for (int mu = 0; mu < 4; ++mu)
                for (int nu = 0; nu < 4; ++nu)
                    for (int axis = 0; axis < 4; ++axis)
                        coordinate(mu) -= connection.gamma(mu, nu, axis).real *
                                          launch->tangent(nu) *
                                          launch->variations[column].displacement(axis);
            double differentiated_null = 0;
            for (int mu = 0; mu < 4; ++mu) {
                const double dx = (samples[0].position(mu) - 8 * samples[1].position(mu) +
                                   8 * samples[2].position(mu) - samples[3].position(mu)) /
                                  (12 * h);
                const double dk = (samples[0].tangent(mu) - 8 * samples[1].tangent(mu) +
                                   8 * samples[2].tangent(mu) - samples[3].tangent(mu)) /
                                  (12 * h);
                EXPECT_NEAR(launch->variations[column].displacement(mu), dx, 2e-10);
                EXPECT_NEAR(coordinate(mu), dk, 2e-10);
                double frozen = 0;
                for (int axis = 0; axis < 3; ++axis)
                    frozen += launch->observer.spatial[axis](mu) *
                              projection->ray.phase_space->direction[axis][column];
                if (column >= 2) frame_contribution += std::abs(coordinate(mu) - frozen);
                for (int nu = 0; nu < 4; ++nu) {
                    differentiated_null +=
                        2 * g(mu, nu).real * launch->tangent(mu) * coordinate(nu);
                    for (int axis = 0; axis < 4; ++axis)
                        differentiated_null += dg(axis, mu, nu).real *
                                               launch->variations[column].displacement(axis) *
                                               launch->tangent(mu) * launch->tangent(nu);
                }
            }
            EXPECT_NEAR(differentiated_null, 0, 2e-13);
        }
        EXPECT_GT(frame_contribution, 1e-4);
    }
}

TEST(CameraFilmDifferential, FilmColumnsRecoverAngularScreenAndPinholeHasNoPupilColumns) {
    KerrSchildFamily metric(KerrSchildParams::Kerr(1.0, 0.6));
    CameraConfig c;
    c.beta_x = 0.2;
    c.beta_y = -0.1;
    c.beta_z = 0.3;
    c.yaw = 0.3f;
    c.pitch = 0.2f;
    c.roll = -0.4f;
    PinholeCamera camera(c);
    const auto projection = camera.ProjectFilmForObserver(437.25, 720.125);
    ASSERT_TRUE(projection);
    const auto launch = LaunchCameraRay(metric, 0.6, projection->ray);
    ASSERT_TRUE(launch);
    const auto angular = CameraAngularVariations(launch->variations, projection->ray.phase_space);
    ASSERT_TRUE(angular);
    const auto screen = relativity::ObserverScreenBasis(
        launch->observer,
        {projection->ray.direction(1), projection->ray.direction(2), projection->ray.direction(3)});
    ASSERT_TRUE(screen);
    for (int column = 0; column < 2; ++column)
        for (int mu = 0; mu < 4; ++mu) {
            EXPECT_DOUBLE_EQ((*angular)[column].displacement(mu), 0);
            EXPECT_NEAR((*angular)[column].derivative(mu), (*screen)[column](mu), 2e-14);
            EXPECT_DOUBLE_EQ(launch->variations[column + 2].displacement(mu), 0);
            EXPECT_DOUBLE_EQ(launch->variations[column + 2].derivative(mu), 0);
        }
    auto invalid = projection->ray;
    invalid.phase_space->direction[0][2] = std::numeric_limits<double>::quiet_NaN();
    EXPECT_FALSE(LaunchCameraRay(metric, 0.6, invalid));
    invalid = projection->ray;
    invalid.phase_space->film_to_angle = {};
    EXPECT_FALSE(CameraAngularVariations(launch->variations, invalid.phase_space));
}
}  // namespace sirius::core::test
