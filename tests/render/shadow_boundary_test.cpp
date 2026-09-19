#include "sirius/backend/cpu/geodesic_tracer.h"
#include "sirius/core/camera.h"
#include "sirius/core/metrics/kerr_schild_family.h"
#include "sirius/core/observer_frame.h"

#include <gtest/gtest.h>

#include "kerr_shadow_oracle.h"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <optional>
#include <utility>
#include <vector>

namespace sirius::render::test {

namespace {

using ScreenPoint = KerrShadowScreenPoint;

class ShadowClassifier {
  public:
    static constexpr int kWidth = 1920;
    static constexpr int kHeight = 1080;
    static constexpr double kObserverRadius = 50.0;
    static constexpr double kInclination = 60.0 * std::numbers::pi / 180.0;
    static constexpr double kFovDegrees = 100.0;

    explicit ShadowClassifier(double spin)
        : metric_(MetricParameters(spin)),
          observer_(MakeStationaryKerrObserver(metric_, spin, kObserverRadius, kInclination)),
          tracer_(&metric_, TracerParameters()),
          camera_(CameraParameters()) {}

    std::optional<ScreenPoint> FiniteObserverPoint(const ScreenPoint& asymptotic) {
        if (!observer_.has_value()) return std::nullopt;
        return ProjectBardeenAtFiniteObserver(asymptotic, *observer_, 1.0, metric_.GetParams().a,
                                              kObserverRadius, kInclination);
    }

    std::optional<bool> Classify(const ScreenPoint& point) {
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
        core::CameraRay ray = camera_.GenerateRay(pixel_x, pixel_y, u, v);

        // Match the stationary observer used by Bardeen rather than the
        // renderer's default Kerr-Schild Eulerian slicing observer.
        if (!observer_.has_value()) {
            ADD_FAILURE() << "stationary Kerr observer is not represented";
            return std::nullopt;
        }
        ray.beta_forward = observer_->screen_beta[0];
        ray.beta_up = observer_->screen_beta[1];
        ray.beta_right = observer_->screen_beta[2];
        const auto result = tracer_.Trace(ray);
        // A marginal photon can exhaust the finite work budget without choosing
        // a physical side of the shadow. It must not move either bracket edge.
        if (!result.cancelled && result.numerical_failure &&
            result.outcome == backend::TraceResult::Outcome::MaxSteps &&
            result.coupled_failure == core::CoupledStepFailure::WorkLimit &&
            result.steps_taken == TracerParameters().max_steps &&
            result.integrator_termination == 0)
            return std::nullopt;
        EXPECT_FALSE(result.cancelled);
        EXPECT_FALSE(result.numerical_failure);
        EXPECT_NE(result.outcome, backend::TraceResult::Outcome::MaxSteps)
            << "P1 classifier did not reach a physical outcome; attempts=" << result.steps_taken
            << ", min_r=" << result.min_radius << ", final=(" << result.final_position(1) << ", "
            << result.final_position(2) << ", " << result.final_position(3)
            << "), integrator_termination=" << result.integrator_termination
            << ", coupled_failure=" << core::CoupledStepFailureName(result.coupled_failure)
            << ", accepted_affine_distance=" << result.affine_length;
        if (result.cancelled || result.numerical_failure) return std::nullopt;
        if (result.outcome != backend::TraceResult::Outcome::Horizon &&
            result.outcome != backend::TraceResult::Outcome::Escaped) {
            ADD_FAILURE() << "shadow classifier reached an unexpected nonphysical outcome";
            return std::nullopt;
        }
        if (result.outcome == backend::TraceResult::Outcome::Horizon) {
            const double terminal_radius = metric_.ComputeKerrRadius(
                result.final_position(1), result.final_position(2), result.final_position(3));
            EXPECT_NEAR(terminal_radius, metric_.OuterHorizonRadius(), 2.0e-12)
                << "captured ray was not localised on the exact Kerr horizon";
        }
        return result.outcome == backend::TraceResult::Outcome::Horizon;
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
    std::optional<KerrStationaryObserver> observer_;
    backend::GeodesicTracer tracer_;
    core::PinholeCamera camera_;
};

}  // namespace

TEST(ShadowBoundary, KerrNearExtremalMatchesBardeenWithinOnePixelAt1080p) {
    constexpr double kSpin = 0.998;
    constexpr ScreenPoint kAnalyticCentre{2.1573218480479185, 0.0};
    ShadowClassifier classifier(kSpin);

    // Sample the full visible upper curve, including both near-equatorial
    // endpoints, rather than a few central witnesses.
    for (const double photon_radius : {1.12, 1.2, 1.35, 1.5, 1.8, 2.1, 2.5, 3.0, 3.5, 3.7}) {
        const auto analytic =
            BardeenShadowPoint(photon_radius, kSpin, ShadowClassifier::kInclination);
        ASSERT_TRUE(analytic.has_value());

        // Camera +right launches the past-directed ray toward +phi.  The
        // corresponding future photon has negative Lz, so Bardeen alpha is
        // positive: camera-right and alpha have the same sign.
        const auto finite_observer = classifier.FiniteObserverPoint(*analytic);
        ASSERT_TRUE(finite_observer.has_value());
        const ScreenPoint camera_convention = *finite_observer;
        const ScreenPoint delta{camera_convention.alpha - kAnalyticCentre.alpha,
                                camera_convention.beta - kAnalyticCentre.beta};
        auto scaled = [&](double scale) {
            return ScreenPoint{kAnalyticCentre.alpha + scale * delta.alpha,
                               kAnalyticCentre.beta + scale * delta.beta};
        };

        double inside = 0.70;
        double outside = 1.30;
        SCOPED_TRACE(photon_radius);
        const auto initial_inside = classifier.Classify(scaled(inside));
        const auto initial_outside = classifier.Classify(scaled(outside));
        ASSERT_TRUE(initial_inside.has_value());
        ASSERT_TRUE(initial_outside.has_value());
        ASSERT_TRUE(*initial_inside);
        ASSERT_FALSE(*initial_outside);
        // Certify a physical bracket at subpixel resolution. Alternating
        // off-centre probes avoids requiring the analytic separatrix itself to
        // terminate, while retaining the same fourteen-ray refinement budget.
        constexpr double kBracketPixels = 0.25;
        for (int iteration = 0; iteration < 14; ++iteration) {
            if (classifier.PixelDistance(scaled(inside), scaled(outside)) <= kBracketPixels) break;
            const double fraction = iteration % 2 == 0 ? 1.0 / 3.0 : 2.0 / 3.0;
            const double trial = inside + fraction * (outside - inside);
            const auto captured = classifier.Classify(scaled(trial));
            if (!captured.has_value()) continue;
            if (*captured)
                inside = trial;
            else
                outside = trial;
        }
        const double width = classifier.PixelDistance(scaled(inside), scaled(outside));
        const double error = std::max(classifier.PixelDistance(scaled(inside), camera_convention),
                                      classifier.PixelDistance(scaled(outside), camera_convention));
        ASSERT_LE(width, kBracketPixels);
        EXPECT_LT(error, 1.0) << "confirmed capture/escape bracket=[" << inside << ", " << outside
                              << "]";
    }
}

TEST(ShadowBoundary, SchwarzschildCriticalImpactParameterMatchesAnalyticAt1080p) {
    ShadowClassifier classifier(0.0);
    const double critical = 3.0 * std::sqrt(3.0);
    double inside = 0.7;
    double outside = 1.3;
    for (int iteration = 0; iteration < 14; ++iteration) {
        const double middle = 0.5 * (inside + outside);
        const auto captured = classifier.Classify({middle * critical, 0.0});
        ASSERT_TRUE(captured.has_value());
        if (*captured) {
            inside = middle;
        } else {
            outside = middle;
        }
    }
    const ScreenPoint measured{0.5 * (inside + outside) * critical, 0.0};
    EXPECT_LT(classifier.PixelDistance(measured, {critical, 0.0}), 1.0);
}

}  // namespace sirius::render::test
