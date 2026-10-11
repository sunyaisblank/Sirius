#include "sirius/backend/cpu/geodesic_tracer.h"
#include "sirius/core/camera.h"
#include "sirius/core/metrics/kerr_schild_family.h"
#include "sirius/render/session/point_source_detector.h"
#include "sirius/render/session/ray_work_queue.h"

#ifdef SIRIUS_HAS_RETAINED_COMPUTE
#include "sirius/backend/retained_compute.h"
#include "sirius/backend/retained_trace_executor.h"
#endif

#include <gtest/gtest.h>

#include "../support/flat_point_star_reference.h"
#include "../support/point_source_band_reference.h"
#include "../support/separated_geodesic_reference.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <format>
#include <limits>
#include <memory>
#include <numbers>
#include <stdexcept>
#include <vector>

namespace sirius::test {
namespace {
using namespace sirius::render;

void CheckIndependentMovingKerrFlux(backend::TraceStepExecutor* executor = nullptr) {
    namespace reference = separated_reference;
    constexpr double mass = 1, spin = .7, sigma_pixels = .3;
    constexpr double film_x = 64.5, film_y = 48.5;
    core::CameraConfig camera_config;
    camera_config.width = 128;
    camera_config.height = 96;
    camera_config.r = 8;
    camera_config.theta = 1.1;
    camera_config.phi = .31;
    camera_config.fov = 12;
    camera_config.yaw = static_cast<float>(std::numbers::pi);
    camera_config.beta_x = .08;
    camera_config.beta_y = .025;
    camera_config.beta_z = -.015;
    core::PinholeCamera camera(camera_config);
    core::KerrSchildFamily metric(core::KerrSchildParams::Kerr(mass, spin));
    const auto independent = [&](const DetectorCoordinate& z, long double refinement) {
        const auto film = camera.ProjectFilmForObserver(film_x + sigma_pixels * z[0],
                                                        film_y + sigma_pixels * z[1]);
        if (!film) throw std::runtime_error("independent detector camera unavailable");
        const auto launch = core::LaunchCameraRay(metric, spin, film->ray);
        if (!launch) throw std::runtime_error("independent detector launch unavailable");
        const auto result = reference::Trace(*launch, mass, spin, 40, refinement);
        if (result.fate != reference::Fate::Escape)
            throw std::runtime_error("independent detector patch is not outward escaping");
        return result;
    };
    const DetectorCoordinate specified_image{.35, -.25};
    const auto source = independent(specified_image, .001L);
    core::StarEntry star{};
    star.direction_x = static_cast<float>(source.infinity_direction[0]);
    star.direction_y = static_cast<float>(source.infinity_direction[1]);
    star.direction_z = static_cast<float>(source.infinity_direction[2]);
    star.distance_pc = 10;
    star.magnitude = 4;
    star.temperature_K = 6500;
    const auto target =
        reference::Unit(reference::Three{star.direction_x, star.direction_y, star.direction_z});
    const auto basis = reference::Basis(target);
    // Quantization moves the represented catalogue direction slightly. Solve
    // its image independently; a production centre never positions this star.
    const auto independent_flux = [&](long double refinement, double difference_step) {
        DetectorCoordinate root = specified_image;
        core::AngularMatrix2 derivative{};
        reference::Result centre{};
        for (int iteration = 0; iteration < 3; ++iteration) {
            centre = independent(root, refinement);
            derivative = {};
            for (int column = 0; column < 2; ++column) {
                auto lower = root, upper = root;
                lower[column] -= difference_step;
                upper[column] += difference_step;
                const auto a = independent(lower, refinement);
                const auto b = independent(upper, refinement);
                for (int row = 0; row < 2; ++row)
                    for (int axis = 0; axis < 3; ++axis)
                        derivative[row][column] += static_cast<double>(
                            basis[row][axis] *
                            (b.infinity_direction[axis] - a.infinity_direction[axis]) /
                            (2 * difference_step));
            }
            const double det =
                derivative[0][0] * derivative[1][1] - derivative[0][1] * derivative[1][0];
            if (!(std::abs(det) > 1e-10))
                throw std::runtime_error("independent detector map is ill conditioned");
            if (iteration < 2) {
                std::array<double, 2> residual{};
                for (int row = 0; row < 2; ++row)
                    for (int axis = 0; axis < 3; ++axis)
                        residual[row] += static_cast<double>(
                            basis[row][axis] * (target[axis] - centre.infinity_direction[axis]));
                root[0] += (derivative[1][1] * residual[0] - derivative[0][1] * residual[1]) / det;
                root[1] += (derivative[0][0] * residual[1] - derivative[1][0] * residual[0]) / det;
            }
        }
        double direction_error = 0;
        for (int axis = 0; axis < 3; ++axis)
            direction_error = std::max(
                direction_error,
                static_cast<double>(std::abs(centre.infinity_direction[axis] - target[axis])));
        if (direction_error > 2e-10)
            throw std::runtime_error("independent detector image root did not converge");
        const double determinant =
            std::abs(derivative[0][0] * derivative[1][1] - derivative[0][1] * derivative[1][0]);
        const double power =
            derivative[0][0] * derivative[0][0] + derivative[0][1] * derivative[0][1] +
            derivative[1][0] * derivative[1][0] + derivative[1][1] * derivative[1][1];
        // power/|det|=kappa+1/kappa bounds this specific image's conditioning.
        if (!(power / determinant < 3))
            throw std::runtime_error("independent detector patch exceeds conditioning bound");
        const double frequency_ratio = static_cast<double>(-1 / centre.energy);
        const auto band =
            point_source_band_reference::Reference(star.temperature_K, frequency_ratio);
        const long double density =
            std::exp(-.5L * (root[0] * root[0] + root[1] * root[1])) /
            (2 * std::numbers::pi_v<long double> * -std::expm1(-8.L) * determinant);
        std::array<double, 3> rgb{};
        for (int channel = 0; channel < 3; ++channel)
            rgb[channel] = static_cast<double>(band[channel] *
                                               std::pow(10.L, -.4L * star.magnitude) * density);
        return rgb;
    };
    const auto coarse_reference = independent_flux(.002L, 2e-3);
    const auto expected = independent_flux(.001L, 1e-3);
    double reference_gap = 0;
    for (int channel = 0; channel < 3; ++channel) {
        ASSERT_GT(expected[channel], 0);
        reference_gap =
            std::max(reference_gap,
                     std::abs(coarse_reference[channel] - expected[channel]) / expected[channel]);
    }
    ASSERT_LT(reference_gap, 2e-6);
    core::StarfieldSpatialIndex catalogue({star});
    std::array<double, 3> errors{};
    constexpr std::array maximum_steps{2.f, .5f, .125f};
    constexpr std::array tolerances{1e-4f, 1e-7f, 1e-10f};
    constexpr std::array root_tolerances{2e-5, 2e-6, 2e-7};
    constexpr std::array flux_tolerances{1e-3, 1e-4, 1e-5};
    for (int refinement = 0; refinement < 3; ++refinement) {
        backend::TracerConfig config;
        config.enable_disk = false;
        config.escape_radius = 40;
        config.max_steps = 10000;
        config.horizon_factor = 1;
        config.integrator.initial_step = config.integrator.max_step = maximum_steps[refinement];
        config.integrator.min_step = 1e-6f;
        config.integrator.abs_tolerance = config.integrator.rel_tolerance = tolerances[refinement];
        backend::GeodesicTracer tracer(&metric, config);
        tracer.SetStepExecutor(executor);
        std::uint64_t central_stages = 0, variation_stages = 0;
        const PointDetectorSampler sample = [&](const DetectorCoordinate& z)
            -> std::expected<PointDetectorProbe, PointDetectorFailure> {
            const auto film = camera.ProjectFilmForObserver(film_x + sigma_pixels * z[0],
                                                            film_y + sigma_pixels * z[1]);
            if (!film || !film->differential || !film->ray.active)
                return std::unexpected(PointDetectorFailure::ProjectionUnavailable);
            const auto traced = tracer.TracePointSource(film->ray);
            central_stages += traced.central_stages;
            variation_stages += traced.variation_stages;
            if (traced.numerical_failure || !traced.beam.infinity_source_map)
                return std::unexpected(PointDetectorFailure::TraceFailed);
            const auto& sky = *traced.beam.infinity_source_map;
            PointDetectorProbe point;
            point.visible = true;
            point.direction = sky.map.direction;
            point.camera_over_source_frequency = 1 / sky.frequency;
            point.inner_attempts = static_cast<std::size_t>(traced.steps_taken);
            point.tail_attempts = sky.attempted_steps;
            for (int row = 0; row < 2; ++row)
                for (int column = 0; column < 2; ++column)
                    for (int angular = 0; angular < 2; ++angular)
                        point.source_derivative[row][column] +=
                            sigma_pixels * sky.map.jacobian[row][angular] *
                            film->differential->angular_jacobian[angular][column];
            return point;
        };
        PointDetectorPolicy policy;
        policy.root_error = root_tolerances[refinement];
        policy.relative_rgb_error = flux_tolerances[refinement];
        // The production detector validates the minimum-depth grid again at
        // its grandchildren: its regular base packet already takes >512 rays.
        // Bound this one-star witness at 1024 without changing accuracy policy.
        policy.maximum_probes = 1024;
        const auto result = EvaluatePointDetector(catalogue, 1, sample, {}, policy);
        const auto& statistics = result ? result->statistics : result.error().statistics;
        ::testing::Test::RecordProperty(std::format("detector_probes_{}", refinement),
                                        static_cast<int>(statistics.probes));
        ::testing::Test::RecordProperty(std::format("detector_cells_{}", refinement),
                                        static_cast<int>(statistics.cells));
        ::testing::Test::RecordProperty(std::format("detector_newton_steps_{}", refinement),
                                        static_cast<int>(statistics.newton_steps));
        ::testing::Test::RecordProperty(std::format("central_stages_{}", refinement),
                                        std::to_string(central_stages));
        ::testing::Test::RecordProperty(std::format("variation_stages_{}", refinement),
                                        std::to_string(variation_stages));
        ::testing::Test::RecordProperty(std::format("tail_attempts_{}", refinement),
                                        std::to_string(statistics.tail_attempts));
        ASSERT_TRUE(result) << "failure=" << static_cast<int>(result.error().reason)
                            << " refinement=" << refinement << " probes=" << statistics.probes
                            << " cells=" << statistics.cells << " roots=" << statistics.roots
                            << " Newton=" << statistics.newton_steps
                            << " central_stages=" << central_stages
                            << " variation_stages=" << variation_stages;
        ASSERT_EQ(result->statistics.roots, 1u);
        for (int channel = 0; channel < 3; ++channel)
            errors[refinement] =
                std::max(errors[refinement],
                         std::abs(result->rgb[channel] - expected[channel]) / expected[channel]);
        ::testing::Test::RecordProperty(std::format("rgb_relative_error_{}", refinement),
                                        std::format("{:.17g}", errors[refinement]));
    }
    // The inherited 1e-4 angular-map target allows 1e-4 relative flux here:
    // bounded kappa<2.62, 2 ppm reference/colour uncertainty, and a Gaussian
    // root error below 2e-7 at a root well inside its compact support.
    EXPECT_LT(errors[0], 5e-4);
    EXPECT_LT(errors[1], 1e-4);
    EXPECT_LT(errors[2], 1e-4);
    EXPECT_LE(errors[2], errors[0] + 2e-6);
    EXPECT_LE(errors[2], errors[1] + 2e-6);
    ::testing::Test::RecordProperty("independent_reference_relative_gap",
                                    std::format("{:.17g}", reference_gap));
    ::testing::Test::RecordProperty(
        "reference_rgb",
        std::format("{:.17g},{:.17g},{:.17g}", expected[0], expected[1], expected[2]));
}

TEST(PointSourceDetector, MovingKerrFluxMatchesIndependentSeparatedImageAndBand) {
    RecordProperty("numerical_backend", "cpu_binary64");
    ASSERT_NO_FATAL_FAILURE(CheckIndependentMovingKerrFlux());
}

TEST(PointSourceDetector, RetainedMovingKerrFluxMatchesIndependentSeparatedImageAndBand) {
#ifdef SIRIUS_HAS_RETAINED_COMPUTE
    const auto inventory = backend::EnumerateVulkanDevices();
    ASSERT_TRUE(inventory) << inventory.error().Description();
    if (inventory->empty()) GTEST_SKIP() << "no Vulkan device; numerical backend unqualified";
    const auto index = backend::ResolveVulkanDeviceIndex(*inventory);
    ASSERT_TRUE(index) << index.error().Description();
    auto device = backend::CreateVulkanDevice(*index);
    ASSERT_TRUE(device) << device.error().Description();
    auto compute = backend::RetainedCompute::Create(**device, 1);
    ASSERT_TRUE(compute) << compute.error().Description();
    backend::RetainedTraceExecutor executor(**compute);
    RecordProperty("numerical_backend", "vulkan_retained_fp32");
    RecordProperty("device", (*device)->Info().name);
    ASSERT_NO_FATAL_FAILURE(CheckIndependentMovingKerrFlux(&executor));
    EXPECT_FALSE(executor.Error());
    const auto statistics = executor.Statistics();
    EXPECT_GT(statistics.camera_batches, 0u);
    EXPECT_GT(statistics.interval_batches, 0u);
    RecordProperty("camera_batches", std::to_string(statistics.camera_batches));
    RecordProperty("interval_batches", std::to_string(statistics.interval_batches));
#else
    GTEST_SKIP() << "retained compute kernels unavailable; numerical backend unqualified";
#endif
}

void CheckRotatingIsolatedStarBrightness(backend::TraceStepExecutor* executor = nullptr) {
    namespace reference = flat_point_star_reference;
    constexpr int edge = 12;
    constexpr std::size_t worker_count = 4;
    constexpr std::array phases{.125, .375, .625, .875};
    constexpr std::array widths{1.0, .3};
    core::CameraConfig camera_config;
    camera_config.width = camera_config.height = edge;
    camera_config.r = 8;
    camera_config.fov = 1;
    core::PinholeCamera camera(camera_config);
    core::KerrSchildFamily metric(core::KerrSchildParams::Kerr(0, 0));
    backend::TracerConfig trace_config;
    trace_config.enable_disk = false;
    trace_config.escape_radius = 40;
    trace_config.max_steps = 64;
    trace_config.integrator.initial_step = trace_config.integrator.max_step = 1;
    trace_config.integrator.min_step = 1e-6f;
    trace_config.integrator.abs_tolerance = trace_config.integrator.rel_tolerance = 1e-7f;
    std::vector<std::unique_ptr<backend::GeodesicTracer>> tracers;
    for (std::size_t worker = 0; worker < worker_count; ++worker) {
        auto tracer = std::make_unique<backend::GeodesicTracer>(&metric, trace_config);
        tracer->SetStepExecutor(executor);
        tracers.push_back(std::move(tracer));
    }
    std::array<std::uint64_t, worker_count> central_stages{}, variation_stages{};
    // Queue ownership makes each worker's tracer and counters exclusive. Every
    // Execute joins its entire batch before geometry or counters are inspected.
    RayWorkQueue workers(worker_count, [&](std::size_t worker, const core::CameraRay& ray) {
        auto result = tracers[worker]->TracePointSource(ray);
        central_stages[worker] += result.central_stages;
        variation_stages[worker] += result.variation_stages;
        return result;
    });
    core::StarEntry star{};
    star.direction_x = 1;
    star.distance_pc = 10;
    star.magnitude = 4;
    star.temperature_K = 6500;
    core::StarfieldSpatialIndex catalogue({star});
    const auto band = point_source_band_reference::Reference(star.temperature_K, 1);
    const long double source_flux = std::pow(10.L, -.4L * star.magnitude);
    // Match the declared stored tangent coefficient, not a quantized ray or a
    // traced centre. CameraFilmDifferential tests qualify this lens convention.
    const float stored_tangent =
        std::tan(camera_config.fov * static_cast<float>(std::numbers::pi) / 360.0f);
    const double tangent = stored_tangent;
    const double pixel_scale = 2 * tangent / edge;
    const double angular_pixel = camera_config.fov * std::numbers::pi / 180 / edge;
    PointDetectorPolicy policy;  // Preserve the production root and RGB allowances.
    std::array<std::array<double, phases.size()>, widths.size()> brightness{},
        reference_brightness{};
    const auto sum_stages = [](const auto& values) {
        std::uint64_t sum = 0;
        for (const auto value : values) sum += value;
        return sum;
    };
    for (std::size_t frame = 0; frame < phases.size(); ++frame) {
        // A quarter-pixel diagonal rotation at a fixed event and zero velocity
        // changes neither physical magnification nor g=1. The eighth-pixel
        // offset keeps this smooth-footprint corpus off the hard |z|=4 edge;
        // the detector's separate ownership tests exercise that discontinuity.
        camera_config.pitch = static_cast<float>(std::atan(pixel_scale * (phases[frame] - .5)));
        camera_config.yaw =
            static_cast<float>(std::numbers::pi + std::atan(pixel_scale * (phases[frame] - .5) *
                                                            std::cos(camera_config.pitch)));
        camera.SetConfig(camera_config);
        const reference::Geometry geometry{
            edge,
            tangent,
            camera_config.yaw,
            camera_config.pitch,
            {std::sin(static_cast<long double>(camera_config.theta)),
             std::cos(static_cast<long double>(camera_config.theta)), 0}};
        const auto image = geometry.Image();
        const double phase_rounding =
            4 * std::numeric_limits<float>::epsilon() * std::abs(camera_config.yaw) / pixel_scale;
        for (int axis = 0; axis < 2; ++axis) {
            ASSERT_NEAR(static_cast<double>(image[axis]), 5.5 + phases[frame], phase_rounding);
            ASSERT_GT(image[axis], 4.5L);
            ASSERT_LT(image[axis], edge - 4.5L);
        }
        ::testing::Test::RecordProperty(
            std::format("frame_{}_image", frame),
            std::format("{:.17g},{:.17g}", double(image[0]), double(image[1])));
        for (std::size_t mode = 0; mode < widths.size(); ++mode) {
            SCOPED_TRACE(::testing::Message() << "frame=" << frame << " width=" << widths[mode]);
            const auto prefix = std::format("{}_frame_{}", mode == 0 ? "beam" : "pinhole", frame);
            const double sigma = widths[mode] * angular_pixel;
            trace_config.enable_ray_bundles = mode == 0;
            trace_config.bundle_point_source = mode == 0;
            trace_config.bundle_angular_size = mode == 0
                                                   ? static_cast<float>(angular_pixel)
                                                   : backend::TracerConfig{}.bundle_angular_size;
            for (auto& tracer : tracers) tracer->SetConfig(trace_config);
            std::vector<PointDetectorFootprint> footprints;
            std::vector<std::array<double, 3>> expected;
            std::vector<long double> solid_angles;
            long double minimum_support_margin = 4;
            for (int y = 0; y < edge; ++y)
                for (int x = 0; x < edge; ++x) {
                    const auto film = camera.ProjectFilmForObserver(x + .5, y + .5);
                    ASSERT_TRUE(film && film->ray.active && film->differential);
                    const auto& p = film->differential->angular_jacobian;
                    const double determinant = std::fma(p[0][0], p[1][1], -p[0][1] * p[1][0]);
                    ASSERT_NE(determinant, 0);
                    // These are the production original camera footprints,
                    // not an envelope or an expected answer from the oracle.
                    footprints.push_back(
                        {{double(x), double(y)},
                         {{{sigma * p[1][1] / determinant, -sigma * p[0][1] / determinant},
                           {-sigma * p[1][0] / determinant, sigma * p[0][0] / determinant}}}});
                    const reference::Coordinate centre{x + .5L, y + .5L};
                    const auto response = geometry.OriginalResponse(image, centre, sigma);
                    minimum_support_margin =
                        std::min(minimum_support_margin, std::abs(response.standard_radius - 4));
                    std::array<double, 3> rgb{};
                    for (int channel = 0; channel < 3; ++channel)
                        rgb[channel] =
                            static_cast<double>(band[channel] * source_flux * response.density);
                    expected.push_back(rgb);
                    solid_angles.push_back(geometry.PixelSolidAngle(centre));
                }
            ASSERT_GT(minimum_support_margin, .01L);
            // Every original in the immediate outside rim has zero response:
            // the complete compact star image lies inside the measured patch.
            for (int y = -1; y <= edge; ++y) {
                for (int x = -1; x <= edge; ++x) {
                    if (x == -1 || x == edge || y == -1 || y == edge) {
                        ASSERT_EQ(
                            geometry.OriginalResponse(image, {x + .5L, y + .5L}, sigma).density, 0);
                    }
                }
            }
            double maximum_frequency_error = 0;
            const PointDetectorProbeBatchSampler batch =
                [&](std::span<const DetectorCoordinate> coordinates) {
                    std::vector<core::CameraFilmProjection> films;
                    std::vector<core::CameraRay> rays;
                    for (const auto& q : coordinates) {
                        const auto film = camera.ProjectFilmOffsetForObserver(.5, .5, q[0], q[1]);
                        if (!film || !film->ray.active || !film->differential)
                            return PointDetectorProbeBatch(
                                coordinates.size(),
                                std::unexpected(PointDetectorFailure::ProjectionUnavailable));
                        films.push_back(*film);
                        rays.push_back(film->ray);
                    }
                    const auto traced = workers.Execute(rays);
                    if (!traced)
                        return PointDetectorProbeBatch(
                            coordinates.size(), std::unexpected(PointDetectorFailure::TraceFailed));
                    PointDetectorProbeBatch values;
                    for (std::size_t i = 0; i < films.size(); ++i) {
                        const auto& ray = (*traced)[i];
                        if (ray.cancelled || ray.numerical_failure ||
                            ray.outcome != backend::TraceResult::Outcome::Escaped ||
                            !ray.beam.infinity_source_map || (mode == 0 && !ray.beam.valid)) {
                            values.push_back(std::unexpected(PointDetectorFailure::TraceFailed));
                            continue;
                        }
                        const auto& sky = *ray.beam.infinity_source_map;
                        PointDetectorProbe point;
                        point.visible = true;
                        point.direction = sky.map.direction;
                        point.camera_over_source_frequency = 1 / sky.frequency;
                        point.inner_attempts = static_cast<std::size_t>(ray.steps_taken);
                        point.tail_attempts = sky.attempted_steps;
                        maximum_frequency_error =
                            std::max(maximum_frequency_error,
                                     std::abs(point.camera_over_source_frequency - 1));
                        for (int row = 0; row < 2; ++row)
                            for (int column = 0; column < 2; ++column)
                                for (int angular = 0; angular < 2; ++angular)
                                    point.source_derivative[row][column] +=
                                        sky.map.jacobian[row][angular] *
                                        films[i].differential->angular_jacobian[angular][column];
                        values.push_back(point);
                    }
                    return values;
                };
            const PointDetectorSampler sample = [&](const DetectorCoordinate& q) {
                const std::array<DetectorCoordinate, 1> coordinates{q};
                return batch(coordinates).front();
            };
            const auto central_before = sum_stages(central_stages);
            const auto variation_before = sum_stages(variation_stages);
            const auto result =
                EvaluatePointDetectorGroup(catalogue, 1, footprints, sample, {}, policy, batch);
            const auto& statistics = result ? result->statistics : result.error().statistics;
            ::testing::Test::RecordProperty(prefix + "_probes", std::to_string(statistics.probes));
            ::testing::Test::RecordProperty(prefix + "_probe_batches",
                                            std::to_string(statistics.probe_batches));
            ::testing::Test::RecordProperty(
                prefix + "_central_stages",
                std::to_string(sum_stages(central_stages) - central_before));
            ::testing::Test::RecordProperty(
                prefix + "_variation_stages",
                std::to_string(sum_stages(variation_stages) - variation_before));
            ::testing::Test::RecordProperty(prefix + "_tail_attempts",
                                            std::to_string(statistics.tail_attempts));
            ::testing::Test::RecordProperty(prefix + "_inner_attempts",
                                            std::to_string(statistics.inner_attempts));
            ASSERT_TRUE(result) << "failure=" << static_cast<int>(result.error().reason)
                                << " probes=" << statistics.probes << " cells=" << statistics.cells
                                << " roots=" << statistics.roots;
            ASSERT_EQ(result->samples.size(), edge * edge);
            EXPECT_GT(statistics.roots, 0u);
            EXPECT_GT(statistics.inner_attempts, 0u);
            EXPECT_GT(sum_stages(central_stages) - central_before, 0u);
            EXPECT_GT(sum_stages(variation_stages) - variation_before, 0u);
            EXPECT_GT(statistics.probe_batches, 0u);
            EXPECT_EQ(statistics.maximum_probe_batch, kPointDetectorProbeBatchSize);
            long double observed_flux = 0, expected_flux = 0, total_solid_angle = 0;
            double maximum_relative_error = 0;
            for (std::size_t i = 0; i < expected.size(); ++i) {
                total_solid_angle += solid_angles[i];
                for (int channel = 0; channel < 3; ++channel) {
                    const double actual = result->samples[i].rgb[channel];
                    ASSERT_TRUE(std::isfinite(actual));
                    EXPECT_NEAR(actual, expected[i][channel],
                                policy.absolute_rgb_error +
                                    policy.relative_rgb_error * expected[i][channel])
                        << "original=" << i << " channel=" << channel;
                    if (expected[i][channel] > 0)
                        maximum_relative_error = std::max(
                            maximum_relative_error,
                            std::abs(actual - expected[i][channel]) / expected[i][channel]);
                    observed_flux += actual * solid_angles[i];
                    expected_flux += expected[i][channel] * solid_angles[i];
                }
            }
            ASSERT_GT(expected_flux, 0);
            ASSERT_GT(observed_flux, 0);
            brightness[mode][frame] = static_cast<double>(observed_flux);
            reference_brightness[mode][frame] = static_cast<double>(expected_flux);
            EXPECT_NEAR(
                brightness[mode][frame], reference_brightness[mode][frame],
                static_cast<double>(policy.relative_rgb_error * reference_brightness[mode][frame] +
                                    policy.absolute_rgb_error * 3 * total_solid_angle));
            ::testing::Test::RecordProperty(prefix + "_brightness",
                                            std::format("{:.17g}", brightness[mode][frame]));
            ::testing::Test::RecordProperty(
                prefix + "_reference_brightness",
                std::format("{:.17g}", reference_brightness[mode][frame]));
            ::testing::Test::RecordProperty(prefix + "_maximum_relative_rgb_error",
                                            std::format("{:.17g}", maximum_relative_error));
            ::testing::Test::RecordProperty(prefix + "_maximum_frequency_error",
                                            std::format("{:.17g}", maximum_frequency_error));
            ::testing::Test::RecordProperty(prefix + "_minimum_support_margin",
                                            std::format("{:.17g}", double(minimum_support_margin)));
        }
    }
    const auto cv = [](const auto& sequence) {
        double mean = 0, variance = 0;
        for (const auto value : sequence) mean += value / sequence.size();
        for (const auto value : sequence)
            variance += (value - mean) * (value - mean) / sequence.size();
        return std::sqrt(variance) / mean;
    };
    const auto [minimum, maximum] = std::minmax_element(brightness[0].begin(), brightness[0].end());
    const auto [reference_minimum, reference_maximum] =
        std::minmax_element(reference_brightness[0].begin(), reference_brightness[0].end());
    // SPECIFICATION's DNGR fidelity bar: James et al. 2015, Appendix A.3.1,
    // arxiv.org/pdf/1502.03808, allows at most 2% brightest/dimmest variation
    // for small smoothly changing footprints. This finite flat-space check
    // does not demand constant flux when actual lensing or Doppler shift varies.
    ASSERT_LE(*reference_maximum / *reference_minimum, 1.02);
    EXPECT_LE(*maximum / *minimum, 1.02);
    const double beam_cv = cv(brightness[0]), pinhole_cv = cv(brightness[1]);
    const double reference_beam_cv = cv(reference_brightness[0]);
    const double reference_pinhole_cv = cv(reference_brightness[1]);
    ASSERT_GT(reference_pinhole_cv / std::max(reference_beam_cv, 1e-12), 1.5);
    EXPECT_LT(beam_cv, pinhole_cv);
    EXPECT_GT(pinhole_cv / std::max(beam_cv, 1e-12), 1.5);
    ::testing::Test::RecordProperty("tracked_star_count", 1);
    ::testing::Test::RecordProperty("complete_original_footprints_per_frame", edge * edge);
    ::testing::Test::RecordProperty("frame_count", static_cast<int>(phases.size()));
    ::testing::Test::RecordProperty("trace_workers", static_cast<int>(worker_count));
    ::testing::Test::RecordProperty("detector_probe_batch_limit",
                                    static_cast<int>(kPointDetectorProbeBatchSize));
    ::testing::Test::RecordProperty("beam_brightest_over_dimmest",
                                    std::format("{:.17g}", *maximum / *minimum));
    ::testing::Test::RecordProperty(
        "reference_beam_brightest_over_dimmest",
        std::format("{:.17g}", *reference_maximum / *reference_minimum));
    ::testing::Test::RecordProperty("beam_brightness_cv", std::format("{:.17g}", beam_cv));
    ::testing::Test::RecordProperty("pinhole_brightness_cv", std::format("{:.17g}", pinhole_cv));
    ::testing::Test::RecordProperty("reference_beam_brightness_cv",
                                    std::format("{:.17g}", reference_beam_cv));
    ::testing::Test::RecordProperty("reference_pinhole_brightness_cv",
                                    std::format("{:.17g}", reference_pinhole_cv));
    ::testing::Test::RecordProperty(
        "evidence_scope",
        "one complete isolated star, four optical rotations at a fixed unboosted flat-space event, "
        "original physical Gaussian responses and mapped escape with vacuum infinity handoff; "
        "conditional smooth-footprint photometry, no "
        "full-frame Kerr or interactive qualification");
}

TEST(PointSourceDetector, RotatingIsolatedStarMatchesIndependentSharedRegionBrightness) {
    RecordProperty("numerical_backend", "cpu_binary64");
    ASSERT_NO_FATAL_FAILURE(CheckRotatingIsolatedStarBrightness());
}

TEST(PointSourceDetector, RetainedRotatingIsolatedStarMatchesIndependentSharedRegionBrightness) {
#ifdef SIRIUS_HAS_RETAINED_COMPUTE
    const auto inventory = backend::EnumerateVulkanDevices();
    ASSERT_TRUE(inventory) << inventory.error().Description();
    if (inventory->empty()) GTEST_SKIP() << "no Vulkan device; numerical backend unqualified";
    const auto index = backend::ResolveVulkanDeviceIndex(*inventory);
    ASSERT_TRUE(index) << index.error().Description();
    auto device = backend::CreateVulkanDevice(*index);
    ASSERT_TRUE(device) << device.error().Description();
    const auto buffer_bytes =
        backend::RetainedCompute::RequiredBufferBytes(kPointDetectorProbeBatchSize);
    ASSERT_LE(buffer_bytes, 8u * 1024 * 1024);
    auto compute = backend::RetainedCompute::Create(**device, kPointDetectorProbeBatchSize, false);
    ASSERT_TRUE(compute) << compute.error().Description();
    backend::RetainedTraceExecutor executor(**compute);
    RecordProperty("numerical_backend", "vulkan_retained_fp32");
    RecordProperty("precision_products", "binary32");
    RecordProperty("retained_storage", "three_component_binary32");
    RecordProperty("device", (*device)->Info().name);
    RecordProperty("retained_buffer_bytes", std::to_string(buffer_bytes));
    ASSERT_NO_FATAL_FAILURE(CheckRotatingIsolatedStarBrightness(&executor));
    EXPECT_FALSE(executor.Error());
    const auto statistics = executor.Statistics();
    EXPECT_GT(statistics.camera_batches, 0u);
    EXPECT_GT(statistics.interval_batches, 0u);
    EXPECT_LE(statistics.maximum_batch_rows, kPointDetectorProbeBatchSize);
    RecordProperty("camera_batches", std::to_string(statistics.camera_batches));
    RecordProperty("interval_batches", std::to_string(statistics.interval_batches));
    RecordProperty("maximum_retained_batch_rows", std::to_string(statistics.maximum_batch_rows));
#else
    GTEST_SKIP() << "retained compute kernels unavailable; numerical backend unqualified";
#endif
}

TEST(PointSourceDetector, LargeSharedRegionMatchesOriginalMovingKerrPackets) {
    core::CameraConfig camera_config;
    camera_config.width = 192;
    camera_config.height = 128;
    camera_config.r = 50;
    camera_config.theta = 1.5708;
    camera_config.fov = 2;
    camera_config.beta_x = .1;
    camera_config.beta_y = .8;
    camera_config.focus_distance = 50;
    core::ThinLensCamera camera(camera_config);
    core::KerrSchildFamily metric(core::KerrSchildParams::Kerr(1, .7));
    backend::TracerConfig config;
    config.enable_disk = false;
    config.escape_radius = 200;
    config.max_steps = 20000;
    config.horizon_factor = 1;
    config.integrator.initial_step = .1f;
    config.integrator.max_step = 2;
    config.integrator.min_step = 1e-5f;
    config.integrator.abs_tolerance = 5e-6f;
    config.integrator.rel_tolerance = 5e-6f;
    backend::GeodesicTracer tracer(&metric, config);
    constexpr float pupil_u = .2f, pupil_v = 1.0f / 7.0f;
    const PointDetectorSampler sample = [&](const DetectorCoordinate& q)
        -> std::expected<PointDetectorProbe, PointDetectorFailure> {
        const auto film =
            camera.ProjectFilmOffsetForObserver(80.5, 48.5, q[0], q[1], pupil_u, pupil_v);
        if (!film || !film->differential || !film->ray.active)
            return std::unexpected(PointDetectorFailure::ProjectionUnavailable);
        const auto traced = tracer.TracePointSource(film->ray);
        if (traced.numerical_failure || !traced.beam.infinity_source_map)
            return std::unexpected(PointDetectorFailure::TraceFailed);
        const auto& sky = *traced.beam.infinity_source_map;
        PointDetectorProbe point;
        point.visible = true;
        point.direction = sky.map.direction;
        point.camera_over_source_frequency = 1 / sky.frequency;
        point.inner_attempts = static_cast<std::size_t>(traced.steps_taken);
        point.tail_attempts = sky.attempted_steps;
        for (int row = 0; row < 2; ++row)
            for (int column = 0; column < 2; ++column)
                for (int angular = 0; angular < 2; ++angular)
                    point.source_derivative[row][column] +=
                        sky.map.jacobian[row][angular] *
                        film->differential->angular_jacobian[angular][column];
        return point;
    };
    // Two known physical images give nonzero flux without relying on a random
    // catalogue landing inside this small patch. Both detectors solve the same
    // represented float directions, including their small rounding displacement.
    std::vector<core::StarEntry> stars;
    for (const DetectorCoordinate q : {DetectorCoordinate{14.17, 14.31}, {19.23, 17.19}}) {
        const auto point = sample(q);
        ASSERT_TRUE(point);
        auto star = core::StarEntry{};
        star.direction_x = static_cast<float>(point->direction[0]);
        star.direction_y = static_cast<float>(point->direction[1]);
        star.direction_z = static_cast<float>(point->direction[2]);
        star.distance_pc = 10;
        star.magnitude = 4;
        star.temperature_K = 5800;
        stars.push_back(star);
    }
    core::StarfieldSpatialIndex catalogue(std::move(stars));
    const double sigma = .3 * 2 * std::numbers::pi / 180 / 128;
    std::vector<PointDetectorFootprint> footprints;
    for (int y = 0; y < 32; ++y)
        for (int x = 0; x < 32; ++x) {
            const auto film = camera.ProjectFilmForObserver(80.5 + x, 48.5 + y, pupil_u, pupil_v);
            ASSERT_TRUE(film && film->differential);
            const auto& p = film->differential->angular_jacobian;
            const double determinant = p[0][0] * p[1][1] - p[0][1] * p[1][0];
            footprints.push_back(
                {{double(x), double(y)},
                 {{{sigma * p[1][1] / determinant, -sigma * p[0][1] / determinant},
                   {-sigma * p[1][0] / determinant, sigma * p[0][0] / determinant}}}});
        }
    const PointDetectorProbeBatchSampler batch =
        [&](std::span<const DetectorCoordinate> coordinates) {
            PointDetectorProbeBatch values;
            values.reserve(coordinates.size());
            for (const auto& q : coordinates) values.push_back(sample(q));
            return values;
        };
    const auto group = EvaluatePointDetectorGroup(catalogue, 1, footprints, sample, {}, {}, batch);
    ASSERT_TRUE(group) << static_cast<int>(group.error().reason) << ' '
                       << group.error().statistics.probes;
    ASSERT_EQ(group->samples.size(), 1024u);
    EXPECT_LT(group->statistics.probes, 2048u);
    EXPECT_LT(group->statistics.probe_batches * 4, group->statistics.probes);
    EXPECT_EQ(group->statistics.maximum_probe_batch, kPointDetectorProbeBatchSize);
    double maximum_relative_difference = 0;
    for (const auto index : {14 * 32 + 14, 17 * 32 + 19}) {
        const auto& footprint = footprints[index];
        const auto& m = footprint.chart_from_standard;
        const auto original = EvaluatePointDetector(
            catalogue, 1,
            [&](const DetectorCoordinate& z)
                -> std::expected<PointDetectorProbe, PointDetectorFailure> {
                auto point = sample({footprint.centre[0] + m[0][0] * z[0] + m[0][1] * z[1],
                                     footprint.centre[1] + m[1][0] * z[0] + m[1][1] * z[1]});
                if (!point) return point;
                const auto derivative = point->source_derivative;
                for (int row = 0; row < 2; ++row)
                    for (int column = 0; column < 2; ++column)
                        point->source_derivative[row][column] =
                            derivative[row][0] * m[0][column] + derivative[row][1] * m[1][column];
                return point;
            },
            {});
        ASSERT_TRUE(original);
        for (int channel = 0; channel < 3; ++channel) {
            ASSERT_GT(original->rgb[channel], 0);
            const double difference =
                std::abs(group->samples[index].rgb[channel] - original->rgb[channel]) /
                original->rgb[channel];
            maximum_relative_difference = std::max(maximum_relative_difference, difference);
            EXPECT_LT(difference, 2e-6);
        }
    }
    RecordProperty("shared_probes", static_cast<int>(group->statistics.probes));
    RecordProperty("probe_batches", static_cast<int>(group->statistics.probe_batches));
    RecordProperty("maximum_probe_batch", static_cast<int>(group->statistics.maximum_probe_batch));
    RecordProperty("shared_inner_attempts", static_cast<int>(group->statistics.inner_attempts));
    RecordProperty("shared_tail_attempts", static_cast<int>(group->statistics.tail_attempts));
    RecordProperty("maximum_relative_difference",
                   std::format("{:.17g}", maximum_relative_difference));
}

}  // namespace
}  // namespace sirius::test
