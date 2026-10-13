// Finite public-camera CPU ray witnesses, not image or all-domain acceptance.
#include "sirius/backend/cpu/geodesic_tracer.h"
#include "sirius/core/camera.h"
#include "sirius/core/camera_launch.h"
#include "sirius/core/metrics/cpu_metric_factory.h"
#include "sirius/render/session/render_session.h"
#include "sirius/render/trace_domain.h"

#include <gtest/gtest.h>

#include "support/alcubierre_axial_reference.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <numbers>
#include <sstream>
#include <string>

namespace {
namespace reference = sirius::test::alcubierre_axial_reference;
using sirius::backend::GeodesicTracer;
using sirius::backend::TracerConfig;
using sirius::backend::TraceResult;

struct Case {
    const char* name;
    double velocity, radius, sigma;
    unsigned turns;
};
constexpr std::array<Case, 6> kCases{{
    {"flat", 0, 1, 4, 0},
    {"forward", 0.6, 1, 4, 0},
    {"reversed", -0.6, 1, 4, 0},
    {"superluminal", 1.5, 1, 4, 2},
    {"small", 0.6, 0.1, 40, 0},
    {"large", 0.6, 100, 0.04, 0},
}};
constexpr long double kEnvelope = 1.0e-4L;

sirius::render::SessionConfig PublicConfig(const Case& input) {
    sirius::render::SessionConfig config;
    config.backend = sirius::render::RenderBackend::Cpu;
    config.metric_id = sirius::core::MetricId::Alcubierre;
    config.black_hole_mass = 0;
    config.black_hole_spin = config.black_hole_charge = config.cosmological_constant = 0;
    config.warp_velocity = input.velocity;
    config.bubble_radius = input.radius;
    config.bubble_sigma = input.sigma;
    config.observer_distance = 5 * std::max(input.radius, 1 / input.sigma);
    config.observer_inclination = std::numbers::pi / 2;
    config.observer_azimuth = std::numbers::pi;
    config.width = config.height = 3;
    config.tile_size = config.samples_per_pixel = 1;
    config.lens_type = sirius::core::LensType::Pinhole;
    config.camera_beta_forward = config.camera_beta_up = config.camera_beta_right = 0;
    config.enable_disk = config.ray_bundles = config.enable_polarisation = false;
    config.enable_parallel_rendering = config.write_output = config.enable_bloom = false;
    return config;
}

void Record(const std::string& name, long double value) {
    std::ostringstream text;
    text << std::setprecision(17) << value;
    ::testing::Test::RecordProperty(name, text.str());
}
}  // namespace

TEST(CpuAlcubierreAcceptance, PublicCentreRaysMatchIndependentFiniteAxialOrbits) {
    long double maximum_gap = 0, maximum_position = 0, maximum_tangent = 0;
    long double maximum_affine = 0, maximum_sky = 0, maximum_killing = 0, maximum_null = 0;
    int trace_count = 0, maximum_steps = 0;
    for (const auto& input : kCases) {
        SCOPED_TRACE(input.name);
        const auto session = PublicConfig(input);
        ASSERT_FALSE(sirius::render::SessionConfigIssue(session).has_value());
        sirius::core::MetricConstructionParameters parameters;
        parameters.mass = 0;
        parameters.warp_velocity = input.velocity;
        parameters.bubble_radius = input.radius;
        parameters.bubble_sigma = input.sigma;
        auto metric = sirius::core::CreateCpuMetric(session.metric_id, parameters);
        ASSERT_NE(metric, nullptr);
        const auto domain = sirius::render::BuildTraceDomainParameters(
            {session.metric_id, parameters.mass, parameters.cosmological_constant,
             session.observer_distance, parameters.throat_radius, input.radius, input.sigma});
        ASSERT_FALSE(domain.finite_causal_boundary);
        const double scale = std::max(input.radius, 1 / input.sigma);
        sirius::core::CameraConfig camera_config;
        camera_config.r = session.observer_distance;
        camera_config.theta = session.observer_inclination;
        camera_config.phi = session.observer_azimuth;
        camera_config.width = camera_config.height = 3;
        camera_config.fov = session.camera_fov;
        sirius::core::PinholeCamera camera(camera_config);
        // Six admitted centre samples, not 54 pixels. Preserve the real ray:
        // double sin/cos(pi) introduces ~1e-15 scaled plane offsets. The
        // reference itself is axial; no live launch output supplies its answer.
        const auto ray = camera.GenerateRay(1, 1, 0.5f, 0.5f);
        ASSERT_EQ(ray.direction(1), -1);
        ASSERT_EQ(ray.direction(2), 0);
        ASSERT_EQ(ray.direction(3), 0);
        ASSERT_EQ(ray.origin(0), 0);
        const reference::Input ref_input{input.velocity,
                                         input.radius,
                                         input.sigma,
                                         0,
                                         -session.observer_distance,
                                         domain.escape_radius,
                                         -1};
        const auto refined = reference::Refine(ref_input);
        const auto& expected = refined.fine;
        ASSERT_LT(refined.maximum_normalized_gap, 1.0e-6L);
        ASSERT_EQ(refined.coarse.coordinate_turning_points, input.turns);
        ASSERT_EQ(expected.coordinate_turning_points, input.turns);
        maximum_gap = std::max(maximum_gap, refined.maximum_normalized_gap);
        Record(std::string(input.name) + "_reference_gap", refined.maximum_normalized_gap);
        const auto launch = sirius::core::LaunchCameraRay(*metric, 0, ray);
        ASSERT_TRUE(launch.has_value());
        for (int axis = 0; axis < 4; ++axis) {
            ASSERT_TRUE(std::isfinite(launch->position(axis)));
            ASSERT_TRUE(std::isfinite(launch->tangent(axis)));
            EXPECT_LE(
                std::abs(launch->position(axis) / scale - expected.initial_position[axis] / scale),
                1.0e-12L);
            EXPECT_LE(std::abs(launch->tangent(axis) - expected.initial_tangent[axis]), 1.0e-12L);
        }
        for (int refinement = 0; refinement < 3; ++refinement) {
            SCOPED_TRACE(refinement);
            TracerConfig config;
            config.enable_disk = config.enable_ray_bundles = config.enable_polarisation = false;
            config.escape_radius = domain.escape_radius;
            config.finite_causal_boundary = domain.finite_causal_boundary;
            config.horizon_factor = sirius::render::kRenderTraceCaptureFactor;
            config.max_steps = sirius::render::kRenderTraceMaximumAttempts;
            config.integrator.initial_step = domain.cpu_initial_step;
            config.integrator.min_step = domain.cpu_min_step;
            config.integrator.max_step = std::ldexp(domain.max_step, -refinement);
            config.integrator.abs_tolerance = std::ldexp(5.0e-6f, -refinement);
            config.integrator.rel_tolerance = config.integrator.abs_tolerance;
            GeodesicTracer tracer(metric.get(), config);
            const auto actual = tracer.Trace(ray);
            ++trace_count;
            ASSERT_EQ(actual.outcome, TraceResult::Outcome::Escaped);
            ASSERT_FALSE(actual.numerical_failure);
            ASSERT_FALSE(actual.cancelled);
            ASSERT_EQ(actual.integrator_termination, 0);
            ASSERT_EQ(actual.terminal_chart, TraceResult::TerminalChart::MetricNative);
            ASSERT_EQ(actual.asymptotic_sheet, TraceResult::AsymptoticSheet::Observer);
            ASSERT_TRUE(actual.final_tangent.has_value());
            ASSERT_TRUE(std::isfinite(actual.affine_length));
            ASSERT_GT(actual.affine_length, 0);
            ASSERT_GT(actual.steps_taken, 0);
            ASSERT_LE(actual.steps_taken, config.max_steps);
            ASSERT_GT(actual.central_stages, 0);
            ASSERT_GT(actual.variation_stages, 0);
            ASSERT_GT(actual.variation_metric_evaluations, 0);
            maximum_steps = std::max(maximum_steps, actual.steps_taken);
            const long double affine_error =
                std::abs(actual.affine_length - expected.affine) / scale;
            long double position_error = 0, tangent_error = 0, sky_error = 0;
            for (int axis = 0; axis < 4; ++axis) {
                ASSERT_TRUE(std::isfinite(actual.final_position(axis)));
                ASSERT_TRUE(std::isfinite((*actual.final_tangent)(axis)));
                ASSERT_TRUE(std::isfinite(actual.final_direction(axis)));
                position_error = std::max(
                    position_error,
                    std::abs(actual.final_position(axis) - expected.position[axis]) / scale);
                tangent_error = std::max(
                    tangent_error, std::abs((*actual.final_tangent)(axis)-expected.tangent[axis]));
                sky_error = std::max(sky_error,
                                     std::abs(actual.final_direction(axis) - expected.sky[axis]));
            }
            // Include observed reference uncertainty in every acceptance bound.
            EXPECT_LT(position_error + refined.maximum_normalized_gap, kEnvelope);
            EXPECT_LT(tangent_error + refined.maximum_normalized_gap, kEnvelope);
            EXPECT_LT(affine_error + refined.maximum_normalized_gap, kEnvelope);
            EXPECT_LT(sky_error + refined.maximum_normalized_gap, kEnvelope);
            const auto& k = *actual.final_tangent;
            const long double q =
                actual.final_position(1) - input.velocity * actual.final_position(0);
            const long double f = reference::Shape(ref_input, q);
            const long double spatial = k(1) - input.velocity * f * k(0);
            const long double defect = -k(0) * k(0) + spatial * spatial + k(2) * k(2) + k(3) * k(3);
            const long double null_error =
                std::abs(defect) /
                std::max(1.0L, k(0) * k(0) + spatial * spatial + k(2) * k(2) + k(3) * k(3));
            const long double killing_error =
                std::abs(reference::Denominator(ref_input, q) * k(0) - expected.killing_constant) /
                std::abs(expected.killing_constant);
            ASSERT_TRUE(std::isfinite(null_error));
            ASSERT_TRUE(std::isfinite(killing_error));
            EXPECT_LT(null_error, 1.0e-6L);
            EXPECT_LT(killing_error + refined.maximum_normalized_gap, kEnvelope);
            // The radial derivative must be outward at the published sphere.
            EXPECT_GT(actual.final_position(1) * k(1) + actual.final_position(2) * k(2) +
                          actual.final_position(3) * k(3),
                      0);
            maximum_position = std::max(maximum_position, position_error);
            maximum_tangent = std::max(maximum_tangent, tangent_error);
            maximum_affine = std::max(maximum_affine, affine_error);
            maximum_sky = std::max(maximum_sky, sky_error);
            maximum_killing = std::max(maximum_killing, killing_error);
            maximum_null = std::max(maximum_null, null_error);
            const std::string prefix =
                std::string(input.name) + "_r" + std::to_string(refinement) + "_";
            Record(prefix + "position_error_per_L", position_error);
            Record(prefix + "tangent_error", tangent_error);
            Record(prefix + "affine_error_per_L", affine_error);
            Record(prefix + "sky_error", sky_error);
            Record(prefix + "killing_error", killing_error);
            Record(prefix + "null_error", null_error);
            Record(prefix + "steps", actual.steps_taken);
        }
    }
    EXPECT_EQ(trace_count, 18);
    RecordProperty("physical_case_count", int(kCases.size()));
    RecordProperty("actual_trace_count", trace_count);
    RecordProperty("maximum_steps", maximum_steps);
    RecordProperty("reference_coarse_panels", 512);
    RecordProperty("reference_fine_panels", 1024);
    Record("reference_max_normalized_refinement_gap", maximum_gap);
    Record("max_normalized_position_error", maximum_position);
    Record("max_tangent_error", maximum_tangent);
    Record("max_normalized_affine_error", maximum_affine);
    Record("max_finite_sky_error", maximum_sky);
    Record("max_relative_killing_error", maximum_killing);
    Record("max_normalized_null_error", maximum_null);
    RecordProperty("reference_scope", "finite_CPU_axial_escape_no_beam_or_transfer_or_horizon");
}
