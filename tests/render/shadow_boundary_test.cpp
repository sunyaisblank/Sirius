#include "sirius/backend/cpu/geodesic_tracer.h"
#include "sirius/core/camera.h"
#include "sirius/core/metrics/kerr_schild_family.h"
#include "sirius/core/observer_frame.h"

#ifdef SIRIUS_HAS_RETAINED_COMPUTE
#include "sirius/backend/retained_trace_executor.h"
#endif

#include <gtest/gtest.h>

#include "kerr_shadow_oracle.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <format>
#include <future>
#include <numbers>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace sirius::render::test {

namespace {

using ScreenPoint = KerrShadowScreenPoint;

constexpr std::array kPhotonRadii{1.12, 1.2, 1.35, 1.5, 1.8, 2.1, 2.5, 3.0, 3.5, 3.7};
constexpr int kBracketTrials = 14;
constexpr double kBracketPixels = .25;

struct ClassificationWork {
    std::uint64_t rays = 0;
    std::uint64_t attempts = 0;
    std::uint64_t maximum_attempts = 0;
    std::uint64_t captures = 0;
    std::uint64_t escapes = 0;
    std::uint64_t work_limits = 0;
};

class ShadowClassifier {
  public:
    static constexpr int kWidth = 1920;
    static constexpr int kHeight = 1080;
    static constexpr double kObserverRadius = 50.0;
    static constexpr double kInclination = 60.0 * std::numbers::pi / 180.0;
    static constexpr double kFovDegrees = 100.0;
    static constexpr int kMaximumTraceAttempts = 30000;

    explicit ShadowClassifier(double spin, backend::TraceStepExecutor* executor = nullptr)
        : metric_(MetricParameters(spin)),
          observer_(MakeStationaryKerrObserver(metric_, spin, kObserverRadius, kInclination)),
          tracer_(&metric_, TracerParameters()),
          camera_(CameraParameters()) {
        tracer_.SetStepExecutor(executor);
    }

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
        ++work_.rays;
        const auto attempts = static_cast<std::uint64_t>(std::max(result.steps_taken, 0));
        work_.attempts += attempts;
        work_.maximum_attempts = std::max(work_.maximum_attempts, attempts);
        EXPECT_GE(result.steps_taken, 0);
        EXPECT_LE(result.steps_taken, kMaximumTraceAttempts);
        EXPECT_TRUE(std::isfinite(result.affine_length));
        // A marginal photon can exhaust the finite work budget without choosing
        // a physical side of the shadow. It must not move either bracket edge.
        if (!result.cancelled && result.numerical_failure &&
            result.outcome == backend::TraceResult::Outcome::MaxSteps &&
            result.coupled_failure == core::CoupledStepFailure::WorkLimit &&
            result.steps_taken == TracerParameters().max_steps &&
            result.integrator_termination == 0) {
            ++work_.work_limits;
            return std::nullopt;
        }
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
            ++work_.captures;
            const double terminal_radius = metric_.ComputeKerrRadius(
                result.final_position(1), result.final_position(2), result.final_position(3));
            EXPECT_NEAR(terminal_radius, metric_.OuterHorizonRadius(), 2.0e-12)
                << "captured ray was not localised on the exact Kerr horizon";
        } else {
            ++work_.escapes;
        }
        return result.outcome == backend::TraceResult::Outcome::Horizon;
    }

    ClassificationWork Work() const { return work_; }

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
        config.max_steps = kMaximumTraceAttempts;
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
    ClassificationWork work_;
};

struct BoundaryMeasurement {
    bool complete = false;
    double width_pixels = 0;
    double displacement_pixels = 0;
    ClassificationWork work;
};

void CheckBardeenBracket(ShadowClassifier& classifier, double photon_radius,
                         BoundaryMeasurement& measurement) {
    constexpr ScreenPoint kAnalyticCentre{2.1573218480479185, 0.0};
    SCOPED_TRACE(photon_radius);
    const auto analytic = BardeenShadowPoint(photon_radius, .998, ShadowClassifier::kInclination);
    ASSERT_TRUE(analytic.has_value());

    // Camera +right launches the past-directed ray toward +phi. The
    // corresponding future photon has negative Lz, so camera-right and
    // Bardeen alpha have the same sign.
    const auto finite_observer = classifier.FiniteObserverPoint(*analytic);
    ASSERT_TRUE(finite_observer.has_value());
    const ScreenPoint camera_convention = *finite_observer;
    const ScreenPoint delta{camera_convention.alpha - kAnalyticCentre.alpha,
                            camera_convention.beta - kAnalyticCentre.beta};
    const auto scaled = [&](double scale) {
        return ScreenPoint{kAnalyticCentre.alpha + scale * delta.alpha,
                           kAnalyticCentre.beta + scale * delta.beta};
    };
    double inside = .70, outside = 1.30;
    const auto initial_inside = classifier.Classify(scaled(inside));
    const auto initial_outside = classifier.Classify(scaled(outside));
    ASSERT_TRUE(initial_inside.has_value());
    ASSERT_TRUE(initial_outside.has_value());
    ASSERT_TRUE(*initial_inside);
    ASSERT_FALSE(*initial_outside);
    // Unresolved work-limit rays never move a bracket edge. The retained and
    // CPU gates use the same physical sides, sample sites and finite budget.
    for (int iteration = 0; iteration < kBracketTrials; ++iteration) {
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
    measurement.width_pixels = classifier.PixelDistance(scaled(inside), scaled(outside));
    measurement.displacement_pixels =
        std::max(classifier.PixelDistance(scaled(inside), camera_convention),
                 classifier.PixelDistance(scaled(outside), camera_convention));
    ASSERT_LE(measurement.width_pixels, kBracketPixels);
    EXPECT_LT(measurement.displacement_pixels, 1.0)
        << "confirmed capture/escape bracket=[" << inside << ", " << outside << "]";
    measurement.complete = true;
}

#ifdef SIRIUS_HAS_RETAINED_COMPUTE
void CheckRetainedBardeenBoundary(bool fp64_products) {
    const std::string mode = fp64_products ? "fp64_products" : "fp32_products";
    ::testing::Test::RecordProperty("product_mode", mode);
    const auto inventory = backend::EnumerateVulkanDevices();
    ASSERT_TRUE(inventory) << inventory.error().Description();
    if (inventory->empty()) {
#ifdef SIRIUS_TEST_REQUIRE_VULKAN_RUNTIME
        FAIL() << "required retained Bardeen gate has no Vulkan device";
#else
        GTEST_SKIP() << "no Vulkan device; retained Bardeen boundary unqualified";
#endif
    }
    const auto index = backend::ResolveVulkanDeviceIndex(*inventory);
    ASSERT_TRUE(index) << index.error().Description();
    auto opened = backend::CreateVulkanDevice(*index);
    ASSERT_TRUE(opened) << opened.error().Description();
    auto& device = **opened;
    ::testing::Test::RecordProperty("device", device.Info().name);
    auto compute = backend::RetainedCompute::Create(device, kPhotonRadii.size(), fp64_products);
    ::testing::Test::RecordProperty("device_admission",
                                    compute ? "admitted" : "rejected_comparison_not_run");
    if (!compute)
        ::testing::Test::RecordProperty("device_admission_error", compute.error().Description());
    if (fp64_products && (!device.Info().supports_fp64 || !device.Info().rounds_fp64_to_nearest)) {
        ASSERT_FALSE(compute);
        EXPECT_EQ(compute.error().domain(), base::ErrorDomain::kDevice);
        EXPECT_NE(compute.error().detail().find("binary64"), std::string::npos);
        ::testing::Test::RecordProperty(
            "comparison_disposition",
            "unsupported binary64 products; rejection proved; comparison_not_run");
        return;
    }
    ASSERT_TRUE(compute) << compute.error().Description();
    backend::RetainedTraceExecutor executor(**compute);
    const auto started = std::chrono::steady_clock::now();
    std::vector<std::future<BoundaryMeasurement>> workers;
    workers.reserve(kPhotonRadii.size());
    // Private classifier state per photon radius, shared production executor:
    // independent wavefronts can coalesce into ten-row camera/interval batches.
    for (const double photon_radius : kPhotonRadii) {
        workers.push_back(std::async(std::launch::async, [&executor, &mode, photon_radius] {
            SCOPED_TRACE(mode);
            ShadowClassifier classifier(.998, &executor);
            BoundaryMeasurement measurement;
            CheckBardeenBracket(classifier, photon_radius, measurement);
            measurement.work = classifier.Work();
            return measurement;
        }));
    }
    ClassificationWork total;
    double maximum_width = 0, maximum_displacement = 0;
    for (std::size_t row = 0; row < workers.size(); ++row) {
        SCOPED_TRACE(mode);
        SCOPED_TRACE(kPhotonRadii[row]);
        const auto measurement = workers[row].get();
        const auto& work = measurement.work;
        const auto prefix = std::format("photon_radius_{:.3f}", kPhotonRadii[row]);
        ::testing::Test::RecordProperty(prefix + "_bracket_pixels",
                                        std::format("{:.17g}", measurement.width_pixels));
        ::testing::Test::RecordProperty(prefix + "_maximum_displacement_pixels",
                                        std::format("{:.17g}", measurement.displacement_pixels));
        ::testing::Test::RecordProperty(prefix + "_rays", std::to_string(work.rays));
        ::testing::Test::RecordProperty(prefix + "_attempts", std::to_string(work.attempts));
        ::testing::Test::RecordProperty(prefix + "_captures", std::to_string(work.captures));
        ::testing::Test::RecordProperty(prefix + "_escapes", std::to_string(work.escapes));
        ::testing::Test::RecordProperty(prefix + "_unresolved_work_limits",
                                        std::to_string(work.work_limits));
        ASSERT_TRUE(measurement.complete);
        EXPECT_LE(work.rays, kBracketTrials + 2u);
        EXPECT_LE(work.maximum_attempts, ShadowClassifier::kMaximumTraceAttempts);
        EXPECT_GT(work.captures, 0u);
        EXPECT_GT(work.escapes, 0u);
        EXPECT_EQ(work.rays, work.captures + work.escapes + work.work_limits);
        total.rays += work.rays;
        total.attempts += work.attempts;
        total.maximum_attempts = std::max(total.maximum_attempts, work.maximum_attempts);
        total.captures += work.captures;
        total.escapes += work.escapes;
        total.work_limits += work.work_limits;
        maximum_width = std::max(maximum_width, measurement.width_pixels);
        maximum_displacement = std::max(maximum_displacement, measurement.displacement_pixels);
    }
    ASSERT_FALSE(executor.Error());
    const auto statistics = executor.Statistics();
    EXPECT_GT(statistics.camera_batches, 0u);
    EXPECT_GT(statistics.interval_batches, 0u);
    EXPECT_GT(statistics.reused_phases, 0u);
    EXPECT_LE(statistics.maximum_batch_rows, kPhotonRadii.size());
    ::testing::Test::RecordProperty("confirmed_photon_radii", int(kPhotonRadii.size()));
    ::testing::Test::RecordProperty("classification_rays", std::to_string(total.rays));
    ::testing::Test::RecordProperty("classification_attempts", std::to_string(total.attempts));
    ::testing::Test::RecordProperty("maximum_ray_attempts", std::to_string(total.maximum_attempts));
    ::testing::Test::RecordProperty("unresolved_work_limits", std::to_string(total.work_limits));
    ::testing::Test::RecordProperty("maximum_bracket_pixels",
                                    std::format("{:.17g}", maximum_width));
    ::testing::Test::RecordProperty("maximum_displacement_pixels",
                                    std::format("{:.17g}", maximum_displacement));
    ::testing::Test::RecordProperty("camera_batches", std::to_string(statistics.camera_batches));
    ::testing::Test::RecordProperty("interval_batches",
                                    std::to_string(statistics.interval_batches));
    ::testing::Test::RecordProperty("maximum_batch_rows",
                                    std::to_string(statistics.maximum_batch_rows));
    ::testing::Test::RecordProperty("reused_phases", std::to_string(statistics.reused_phases));
    ::testing::Test::RecordProperty("accepted_intervals",
                                    std::to_string(statistics.accepted_intervals));
    ::testing::Test::RecordProperty("rejected_intervals",
                                    std::to_string(statistics.rejected_intervals));
    ::testing::Test::RecordProperty(
        "elapsed_seconds",
        std::format(
            "{:.17g}",
            std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count()));
}
#endif

}  // namespace

TEST(ShadowBoundary, KerrNearExtremalMatchesBardeenWithinOnePixelAt1080p) {
    ShadowClassifier classifier(.998);

    // Sample the full visible upper curve, including both near-equatorial
    // endpoints, rather than a few central witnesses.
    for (const double photon_radius : kPhotonRadii) {
        BoundaryMeasurement measurement;
        ASSERT_NO_FATAL_FAILURE(CheckBardeenBracket(classifier, photon_radius, measurement));
    }
}

TEST(ShadowBoundary, RetainedFp32ProductsKerrNearExtremalBardeenBoundaryAt1080p) {
#ifdef SIRIUS_HAS_RETAINED_COMPUTE
    ASSERT_NO_FATAL_FAILURE(CheckRetainedBardeenBoundary(false));
#elif defined(SIRIUS_TEST_REQUIRE_VULKAN_RUNTIME)
    FAIL() << "required retained Bardeen gate has no compiled retained kernels";
#else
    GTEST_SKIP() << "retained kernels unavailable; Bardeen boundary unqualified";
#endif
}

TEST(ShadowBoundary, RetainedFp64ProductsKerrNearExtremalBardeenBoundaryAt1080p) {
#ifdef SIRIUS_HAS_RETAINED_COMPUTE
    ASSERT_NO_FATAL_FAILURE(CheckRetainedBardeenBoundary(true));
#elif defined(SIRIUS_TEST_REQUIRE_VULKAN_RUNTIME)
    FAIL() << "required retained Bardeen gate has no compiled retained kernels";
#else
    GTEST_SKIP() << "retained kernels unavailable; Bardeen boundary unqualified";
#endif
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
