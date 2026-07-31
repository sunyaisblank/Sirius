#include "sirius/backend/cpu/geodesic_tracer.h"
#include "sirius/core/camera.h"
#include "sirius/core/metrics/kerr_schild_family.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <numbers>
#include <optional>
#include <utility>
#include <vector>

namespace sirius::render::test {

namespace {

struct ScreenPoint {
    double alpha;
    double beta;
};

// Bardeen celestial coordinates for an unstable spherical Kerr photon orbit
// (M=1). Bardeen, Press & Teukolsky (1972), equations for xi(r), eta(r), and
// the observer-screen projection.
std::optional<ScreenPoint> BardeenPoint(double photon_radius, double spin, double inclination) {
    const double r = photon_radius;
    const double a2 = spin * spin;
    const double xi = (r * r * (r - 3.0) + a2 * (r + 1.0)) / (spin * (1.0 - r));
    const double eta =
        r * r * r * (4.0 * a2 - r * (r - 3.0) * (r - 3.0)) / (a2 * (1.0 - r) * (1.0 - r));
    const double sin_i = std::sin(inclination);
    const double cos_i = std::cos(inclination);
    const double alpha = -xi / sin_i;
    const double beta_squared =
        eta + a2 * cos_i * cos_i - xi * xi * cos_i * cos_i / (sin_i * sin_i);
    if (beta_squared < 0.0) return std::nullopt;
    return ScreenPoint{alpha, std::sqrt(beta_squared)};
}

class ShadowClassifier {
  public:
    static constexpr int kWidth = 1920;
    static constexpr int kHeight = 1080;
    static constexpr double kObserverRadius = 50.0;
    static constexpr double kInclination = 60.0 * std::numbers::pi / 180.0;
    static constexpr double kFovDegrees = 100.0;

    explicit ShadowClassifier(double spin)
        : metric_(MetricParameters(spin)),
          tracer_(&metric_, TracerParameters()),
          camera_(CameraParameters()) {}

    bool IsCaptured(const ScreenPoint& point) {
        const double tan_half_fov = std::tan(kFovDegrees * std::numbers::pi / 360.0);
        const double aspect = static_cast<double>(kWidth) / kHeight;
        const double image_x =
            0.5 * kWidth * (1.0 + point.alpha / (kObserverRadius * tan_half_fov * aspect));
        const double image_y =
            0.5 * kHeight * (1.0 - point.beta / (kObserverRadius * tan_half_fov));
        const int pixel_x = static_cast<int>(std::floor(image_x));
        const int pixel_y = static_cast<int>(std::floor(image_y));
        const float u = static_cast<float>(image_x - pixel_x);
        const float v = static_cast<float>(image_y - pixel_y);
        const core::CameraRay ray = camera_.GenerateRay(pixel_x, pixel_y, u, v);
        const auto result = tracer_.Trace(ray);
        EXPECT_FALSE(result.numerical_failure);
        EXPECT_NE(result.outcome, backend::TraceResult::Outcome::MaxSteps)
            << "P1 classifier reached its work bound instead of a physical outcome";
        const auto outcome = result.outcome;
        return outcome == backend::TraceResult::Outcome::Horizon ||
               outcome == backend::TraceResult::Outcome::Spiraling;
    }

    double PixelDistance(const ScreenPoint& lhs, const ScreenPoint& rhs) const {
        const double tan_half_fov = std::tan(kFovDegrees * std::numbers::pi / 360.0);
        const double pixels_per_screen_unit = 0.5 * kHeight / (kObserverRadius * tan_half_fov);
        return std::hypot(lhs.alpha - rhs.alpha, lhs.beta - rhs.beta) * pixels_per_screen_unit;
    }

  private:
    static core::KerrSchildParams MetricParameters(double spin) {
        core::KerrSchildParams parameters;
        parameters.M = 1.0;
        parameters.a = spin;
        return parameters;
    }

    static backend::TracerConfig TracerParameters() {
        backend::TracerConfig config;
        config.enable_disk = false;
        config.enable_spiral_termination = false;
        config.escape_radius = 200.0f;
        config.horizon_factor = 1.0f;
        config.max_steps = 30000;
        config.integrator.initial_step = 0.02f;
        config.integrator.max_step = 0.25f;
        config.integrator.min_step = 1e-6f;
        config.integrator.abs_tolerance = 1e-9f;
        config.integrator.rel_tolerance = 1e-9f;
        config.strong_field_radius = 5.0f;
        config.strong_field_max_step = 0.002f;
        return config;
    }

    static core::CameraConfig CameraParameters() {
        core::CameraConfig config;
        config.r = kObserverRadius;
        config.theta = kInclination;
        config.phi = 0.0;
        config.fov = kFovDegrees;
        config.width = kWidth;
        config.height = kHeight;
        return config;
    }

    core::KerrSchildFamily metric_;
    backend::GeodesicTracer tracer_;
    core::PinholeCamera camera_;
};

}  // namespace

TEST(ShadowBoundary, KerrNearExtremalMatchesBardeenWithinOnePixelAt1080p) {
    constexpr double kSpin = 0.998;
    constexpr ScreenPoint kAnalyticCentre{-2.1573218480479185, 0.0};
    ShadowClassifier classifier(kSpin);

    // Sample the full visible upper curve, including both near-equatorial
    // endpoints, rather than a few central witnesses.
    for (const double photon_radius : {1.12, 1.2, 1.35, 1.5, 1.8, 2.1, 2.5, 3.0, 3.5, 3.7}) {
        const auto analytic = BardeenPoint(photon_radius, kSpin, ShadowClassifier::kInclination);
        ASSERT_TRUE(analytic.has_value());

        const ScreenPoint camera_convention{-analytic->alpha, analytic->beta};
        const ScreenPoint delta{camera_convention.alpha - kAnalyticCentre.alpha,
                                camera_convention.beta - kAnalyticCentre.beta};
        auto scaled = [&](double scale) {
            return ScreenPoint{kAnalyticCentre.alpha + scale * delta.alpha,
                               kAnalyticCentre.beta + scale * delta.beta};
        };

        double inside = 0.70;
        double outside = 1.30;
        ASSERT_TRUE(classifier.IsCaptured(scaled(inside)));
        ASSERT_FALSE(classifier.IsCaptured(scaled(outside)));
        for (int iteration = 0; iteration < 14; ++iteration) {
            const double middle = 0.5 * (inside + outside);
            if (classifier.IsCaptured(scaled(middle))) {
                inside = middle;
            } else {
                outside = middle;
            }
        }
        const ScreenPoint measured = scaled(0.5 * (inside + outside));
        EXPECT_LT(classifier.PixelDistance(measured, camera_convention), 1.0)
            << "photon orbit r/M=" << photon_radius << ", scale bracket=[" << inside << ", "
            << outside << "], analytic=(" << camera_convention.alpha << ", "
            << camera_convention.beta << "), measured=(" << measured.alpha << ", " << measured.beta
            << ")";
    }
}

TEST(ShadowBoundary, SchwarzschildCriticalImpactParameterMatchesAnalyticAt1080p) {
    ShadowClassifier classifier(0.0);
    const double critical = 3.0 * std::sqrt(3.0);
    double inside = 0.7;
    double outside = 1.3;
    for (int iteration = 0; iteration < 14; ++iteration) {
        const double middle = 0.5 * (inside + outside);
        if (classifier.IsCaptured({middle * critical, 0.0})) {
            inside = middle;
        } else {
            outside = middle;
        }
    }
    const ScreenPoint measured{0.5 * (inside + outside) * critical, 0.0};
    EXPECT_LT(classifier.PixelDistance(measured, {critical, 0.0}), 1.0);
}

}  // namespace sirius::render::test
