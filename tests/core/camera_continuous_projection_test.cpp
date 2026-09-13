#include "sirius/core/camera.h"
#include "sirius/core/camera_sampling.h"

#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <limits>

namespace sirius::core::test {
namespace {
using Vector = std::array<double, 3>;

Vector Unit(Vector v) {
    const double norm = std::hypot(v[0], v[1], v[2]);
    for (double& value : v) value /= norm;
    return v;
}

// Independent geometric oracle: spherical angles for fisheye, intersection of
// the fixed pupil/focus plane for ThinLens, successive plane rotations for
// Pinhole. Stored float lens coefficients remain the represented intrinsics.
Vector Reference(const ICamera& camera, double x, double y, const CameraRay& pupil) {
    const auto& c = camera.GetConfig();
    const double aspect = static_cast<float>(c.width) / c.height;
    const double right = (2 * x / c.width - 1) * aspect, up = 1 - 2 * y / c.height;
    const double tangent = std::tan(c.fov * static_cast<float>(std::numbers::pi) / 360.0f);
    if (camera.GetLensType() == LensType::Fisheye) {
        const double theta = std::hypot(right, up) * c.fov *
                             static_cast<double>(static_cast<float>(std::numbers::pi)) / 360;
        const double phi = std::atan2(up, right);
        return {-std::cos(theta), -std::sin(theta) * std::sin(phi),
                std::sin(theta) * std::cos(phi)};
    }
    if (camera.GetLensType() == LensType::ThinLens)
        return Unit({-c.focus_distance, pupil.aperture_up - up * tangent * c.focus_distance,
                     right * tangent * c.focus_distance - pupil.aperture_right});
    Vector v{right * tangent, up * tangent, -1};
    const auto rotate = [&](unsigned i, unsigned j, float angle) {
        const double old = v[i], cosine = std::cos(angle), sine = std::sin(angle);
        v[i] = cosine * old - sine * v[j];
        v[j] = sine * old + cosine * v[j];
    };
    rotate(0, 1, c.roll);
    rotate(1, 2, c.pitch);
    rotate(2, 0, c.yaw);
    return Unit({v[2], -v[1], v[0]});
}

void SameRay(const CameraRay& a, const CameraRay& b) {
    for (int i = 0; i < 4; ++i) {
        EXPECT_EQ(a.origin(i), b.origin(i));
        EXPECT_EQ(a.direction(i), b.direction(i));
    }
    EXPECT_EQ(a.active, b.active);
    EXPECT_EQ(a.aperture_up, b.aperture_up);
    EXPECT_EQ(a.aperture_right, b.aperture_right);
    EXPECT_EQ(a.beta_forward, b.beta_forward);
    EXPECT_EQ(a.beta_up, b.beta_up);
    EXPECT_EQ(a.beta_right, b.beta_right);
}
}  // namespace

TEST(CameraContinuousProjection, GeometryAndDifferentialAcrossOutputCrop) {
    CameraConfig c;
    c.width = 1280;
    c.height = 640;
    c.yaw = .13f;
    c.pitch = -.21f;
    c.roll = .09f;
    for (auto lens : {LensType::Pinhole, LensType::ThinLens, LensType::Fisheye}) {
        auto lens_config = c;
        if (lens == LensType::ThinLens) lens_config.focus_distance = 30;
        const auto camera = CreateCamera(lens, lens_config);
        for (const auto film : {std::array<double, 2>{-12.25, -8.75},
                                {1291.3, 649.1},
                                {437.375, 255.125},
                                {640, 320}}) {
            SCOPED_TRACE(::testing::Message()
                         << camera->GetName() << ' ' << film[0] << ',' << film[1]);
            const auto projected = camera->ProjectFilmForObserver(film[0], film[1], .2f, .7f);
            ASSERT_TRUE(projected);
            ASSERT_TRUE(projected->ray.active);
            ASSERT_TRUE(projected->differential);
            const auto expected = Reference(*camera, film[0], film[1], projected->ray);
            for (unsigned i = 0; i < 3; ++i)
                EXPECT_NEAR(projected->ray.direction(i + 1), expected[i], 5e-15);
            constexpr double h = .002;
            for (unsigned col = 0; col < 2; ++col) {
                auto plus = film, minus = film;
                plus[col] += h;
                minus[col] -= h;
                const auto p = camera->ProjectFilmForObserver(plus[0], plus[1], .2f, .7f);
                const auto m = camera->ProjectFilmForObserver(minus[0], minus[1], .2f, .7f);
                ASSERT_TRUE(p);
                ASSERT_TRUE(m);
                for (unsigned i = 0; i < 3; ++i) {
                    const double fd = (p->ray.direction(i + 1) - m->ray.direction(i + 1)) / (2 * h);
                    EXPECT_NEAR(projected->differential->direction_derivative[i][col], fd, 3e-12);
                }
            }
            const auto& d = projected->differential->direction_derivative;
            const Vector cross{d[1][0] * d[2][1] - d[2][0] * d[1][1],
                               d[2][0] * d[0][1] - d[0][0] * d[2][1],
                               d[0][0] * d[1][1] - d[1][0] * d[0][1]};
            EXPECT_NEAR(projected->differential->solid_angle_density,
                        std::hypot(cross[0], cross[1], cross[2]), 2e-20);
        }
    }
}

TEST(CameraContinuousProjection, FixedPupilAndObserverArePreserved) {
    CameraConfig c;
    c.width = 1280;
    c.height = 640;
    c.focus_distance = 30;
    c.beta_x = .2;
    c.beta_y = -.1;
    c.beta_z = .3;
    ThinLensCamera camera(c);
    ForEachCameraSample(4, [&](const CameraSample& sample) {
        const auto legacy = camera.GenerateRayForObserver(0, 0, sample.image_u, sample.image_v,
                                                          sample.pupil_u, sample.pupil_v);
        const auto projection =
            camera.ProjectFilmForObserver(-3.25, 650.75, sample.pupil_u, sample.pupil_v);
        ASSERT_TRUE(projection);
        const auto& ray = projection->ray;
        EXPECT_EQ(ray.aperture_up, legacy.aperture_up);
        EXPECT_EQ(ray.aperture_right, legacy.aperture_right);
        EXPECT_EQ(ray.beta_forward, c.beta_x);
        EXPECT_EQ(ray.beta_up, c.beta_y);
        EXPECT_EQ(ray.beta_right, c.beta_z);
        for (int i = 0; i < 4; ++i) EXPECT_EQ(ray.origin(i), legacy.origin(i));
        const double distance = -c.focus_distance / ray.direction(1);
        const double t = std::tan(c.fov * static_cast<float>(std::numbers::pi) / 360.0f);
        EXPECT_NEAR(ray.aperture_right + distance * ray.direction(3),
                    (2 * -3.25 / c.width - 1) * 2 * t * c.focus_distance, 2e-13);
        EXPECT_NEAR(ray.aperture_up - distance * ray.direction(2),
                    (1 - 2 * 650.75 / c.height) * t * c.focus_distance, 2e-13);
    });
}

TEST(CameraContinuousProjection, FisheyeCropMaskCentreAndAntipodalRim) {
    CameraConfig c;
    c.width = 100;
    c.height = 100;
    c.fov = 60;
    FisheyeCamera camera(c);
    const auto outside_crop = camera.ProjectFilmForObserver(150, 50);
    ASSERT_TRUE(outside_crop);
    EXPECT_TRUE(outside_crop->ray.active);
    EXPECT_TRUE(outside_crop->differential);
    // theta=pi occurs at normalized radius six for FOV60, rather than one.
    const auto rim = camera.ProjectFilmForObserver(350, 50);
    ASSERT_TRUE(rim);
    EXPECT_TRUE(rim->ray.active);
    EXPECT_FALSE(rim->differential);
    const auto outside_optics = camera.ProjectFilmForObserver(350.1, 50);
    ASSERT_TRUE(outside_optics);
    EXPECT_FALSE(outside_optics->ray.active);
    EXPECT_FALSE(outside_optics->differential);
    EXPECT_TRUE(IsRepresentedCameraRay(outside_optics->ray));
    c.fov = 360;
    camera.SetConfig(c);
    const auto full_rim = camera.ProjectFilmForObserver(100, 50);
    ASSERT_TRUE(full_rim);
    EXPECT_TRUE(full_rim->ray.active);
    EXPECT_FALSE(full_rim->differential);
    const auto centre = camera.ProjectFilmForObserver(50, 50);
    ASSERT_TRUE(centre);
    ASSERT_TRUE(centre->differential);
    EXPECT_EQ(centre->ray.direction(1), -1);
    const double a = static_cast<float>(std::numbers::pi);
    EXPECT_NEAR(centre->differential->solid_angle_density, a * a * .02 * .02, 2e-18);
}

TEST(CameraContinuousProjection, RefinementDeclinesLostCoordinatesAndDirections) {
    PinholeCamera camera;
    const auto small = camera.ProjectFilmOffsetForObserver(960, 540, 1e-8, -1e-8);
    ASSERT_TRUE(small);
    EXPECT_TRUE(small->differential);
    const auto zero = camera.ProjectFilmOffsetForObserver(960, 540, 0, 0);
    ASSERT_TRUE(zero);
    EXPECT_NE(small->ray.direction(3), zero->ray.direction(3));
    const auto lost = camera.ProjectFilmOffsetForObserver(960, 540, 1e-20, 0);
    ASSERT_FALSE(lost);
    EXPECT_EQ(lost.error(), CameraProjectionFailure::Unrepresentable);
    const auto partially_lost = camera.ProjectFilmOffsetForObserver(960, 540, 1e-20, .25);
    ASSERT_FALSE(partially_lost);
    EXPECT_EQ(partially_lost.error(), CameraProjectionFailure::Unrepresentable);
    const auto lost_direction =
        camera.ProjectFilmOffsetForObserver(0, 540, std::numeric_limits<double>::denorm_min(), 0);
    ASSERT_FALSE(lost_direction);
    EXPECT_EQ(lost_direction.error(), CameraProjectionFailure::Unrepresentable);
}

TEST(CameraContinuousProjection, CropTranslationDoesNotRenormalizeOrMaskProjection) {
    CameraConfig c;
    c.width = 1280;
    c.height = 640;
    for (auto lens : {LensType::Pinhole, LensType::ThinLens, LensType::Fisheye}) {
        const auto original = CreateCamera(lens, c);
        auto larger_config = c;
        larger_config.width = 1600;
        const auto larger = CreateCamera(lens, larger_config);
        for (const double x : {-16.0, .25, 1279.75, 1296.0}) {
            const auto a = original->ProjectFilmForObserver(x, 170.25, .2f, .7f);
            const auto b = larger->ProjectFilmForObserver(x + 160, 170.25, .2f, .7f);
            ASSERT_TRUE(a);
            ASSERT_TRUE(b);
            ASSERT_TRUE(a->differential);
            ASSERT_TRUE(b->differential);
            for (int i = 1; i < 4; ++i)
                EXPECT_NEAR(a->ray.direction(i), b->ray.direction(i), 5e-16);
            EXPECT_NEAR(a->differential->solid_angle_density, b->differential->solid_angle_density,
                        2e-20);
        }
    }
}

TEST(CameraContinuousProjection, OriginalPacketsKeepLegacyValuesWithExplicitSmoothRoundoff) {
    CameraConfig c;
    c.width = 1920;
    c.height = 1080;
    c.yaw = .17f;
    c.pitch = -.11f;
    c.roll = .07f;
    c.beta_x = .2;
    double maximum_difference = 0;
    for (auto lens : {LensType::Pinhole, LensType::ThinLens, LensType::Fisheye}) {
        auto lens_config = c;
        if (lens == LensType::ThinLens) lens_config.focus_distance = 30;
        const auto camera = CreateCamera(lens, lens_config);
        ForEachCameraSample(4, [&](const CameraSample& sample) {
            for (const auto pixel : {std::array<int, 2>{0, 0}, {1919, 1079}, {713, 431}}) {
                const auto before =
                    camera->GenerateRayForObserver(pixel[0], pixel[1], sample.image_u,
                                                   sample.image_v, sample.pupil_u, sample.pupil_v);
                const auto smooth = camera->ProjectFilmForObserver(
                    pixel[0] + static_cast<double>(sample.image_u),
                    pixel[1] + static_cast<double>(sample.image_v), sample.pupil_u, sample.pupil_v);
                ASSERT_TRUE(smooth);
                ASSERT_TRUE(smooth->ray.active);
                const auto after =
                    camera->GenerateRayForObserver(pixel[0], pixel[1], sample.image_u,
                                                   sample.image_v, sample.pupil_u, sample.pupil_v);
                SameRay(before, after);
                for (int i = 1; i < 4; ++i)
                    maximum_difference =
                        std::max(maximum_difference,
                                 std::abs(before.direction(i) - smooth->ray.direction(i)));
                EXPECT_EQ(before.aperture_up, smooth->ray.aperture_up);
                EXPECT_EQ(before.aperture_right, smooth->ray.aperture_right);
            }
        });
    }
    // An explicit compatibility budget for legacy float projection arithmetic,
    // not a detector/transport accuracy allowance. No coordinate snapping.
    EXPECT_GT(maximum_difference, 1e-10);
    EXPECT_LT(maximum_difference, 1e-6);
    RecordProperty("maximum_legacy_smooth_direction_difference", maximum_difference);
}

TEST(CameraContinuousProjection, InvalidAndUnsupportedProjectionPublishesNoRay) {
    PinholeCamera camera;
    for (double bad :
         {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity()}) {
        const auto result = camera.ProjectFilmForObserver(bad, 0);
        ASSERT_FALSE(result);
        EXPECT_EQ(result.error(), CameraProjectionFailure::InvalidInput);
    }
    for (float bad : {-1.0f, 1.0f, std::numeric_limits<float>::quiet_NaN()}) {
        const auto result = camera.ProjectFilmForObserver(0, 0, bad, .5f);
        ASSERT_FALSE(result);
        EXPECT_EQ(result.error(), CameraProjectionFailure::InvalidInput);
    }
    const auto overflow = camera.ProjectFilmForObserver(std::numeric_limits<double>::max(), 0);
    ASSERT_FALSE(overflow);
    EXPECT_EQ(overflow.error(), CameraProjectionFailure::Arithmetic);
    class Custom final : public ICamera {
      public:
        CameraRay GenerateRay(int x, int y, float u, float v, float pu, float pv) const override {
            return delegate.GenerateRay(x, y, u, v, pu, pv);
        }
        LensType GetLensType() const override { return LensType::Pinhole; }
        const char* GetName() const override { return "custom"; }
        const CameraConfig& GetConfig() const override { return delegate.GetConfig(); }
        void SetConfig(const CameraConfig& c) override { delegate.SetConfig(c); }

      private:
        PinholeCamera delegate;
    } custom;
    const auto unsupported = custom.ProjectFilmForObserver(0, 0);
    ASSERT_FALSE(unsupported);
    EXPECT_EQ(unsupported.error(), CameraProjectionFailure::Unsupported);
}
}  // namespace sirius::core::test
