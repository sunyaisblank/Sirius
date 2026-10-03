// Complete production camera-to-event rays against an independently evolved
// Carter flow. This is numerical acceptance, with no image or throughput claim.
#include "sirius/backend/cpu/geodesic_tracer.h"

#include <gtest/gtest.h>

#include "support/separated_geodesic_reference.h"

#include <array>
#include <cmath>
#include <format>
#include <future>
#include <numbers>
#include <string>
#include <vector>

#if defined(SIRIUS_HAS_RETAINED_COMPUTE) && defined(SIRIUS_RETAINED_CAMERA_TEST_DIR)
#include "sirius/backend/device.h"
#include "sirius/backend/retained_compute.h"
#include "sirius/backend/retained_trace_executor.h"
#endif

namespace {
using namespace sirius::core;
using namespace sirius::backend;
namespace reference = sirius::test::separated_reference;
using reference::Four;
using reference::Matrix;
using reference::Scalar;
using reference::Three;

// The existing P1 contract specifies 1e-4 relative position and 1e-4 radians.
// The same 1e-4 fractional envelope judges the represented canonical maps and
// beam covariance. Independently refined expectations must stabilize below
// 1% of that envelope before any product result can pass. These are bounded
// numerical observations, not a global interval or uniform separatrix bound.
constexpr double kAcceptance = 1e-4;
constexpr double kReferenceFraction = .01;

struct Case {
    std::string name;
    double mass, spin, boundary;
    CameraRay camera;
    reference::Fate fate;
    bool radial_turn;
    double disk_inner = 0, disk_outer = 0;
};

CameraRay Ray(double radius, double theta, Three n) {
    CameraRay ray;
    ray.origin(1) = radius;
    ray.origin(2) = theta;
    ray.origin(3) = .31;
    n = reference::Unit(n);
    for (int axis = 0; axis < 3; ++axis) ray.direction(axis + 1) = static_cast<double>(n[axis]);
    return ray;
}
CameraRay StaticSchwarzschildRay(double mass, double impact) {
    constexpr double radius = 12;
    const double f = 1 - 2 / radius;
    const double transverse = impact * std::sqrt(f) / radius;
    auto ray = Ray(radius * mass, std::numbers::pi / 2,
                   {-std::sqrt(1 - transverse * transverse), 0, transverse});
    // A static worldline has beta_r=2M/r relative to the ingoing KS normal.
    // The public camera's forward axis is minus its radial tetrad axis.
    ray.beta_forward = -2 / radius;
    return ray;
}

std::vector<Case> Cases() {
    std::vector<Case> cases;
    auto flat = Ray(8, 1.1, {.6, .3, -.74});
    flat.beta_forward = .2;
    flat.beta_up = -.1;
    flat.beta_right = .05;
    flat.aperture_right = -.017;
    flat.aperture_up = .021;
    cases.push_back({"flat_moving_pupil", 0, 0, 32, flat, reference::Fate::Escape, false});
    cases.push_back({"schwarzschild_capture", 1, 0, 40, StaticSchwarzschildRay(1, 2),
                     reference::Fate::Capture, false});
    cases.push_back({"schwarzschild_turn_mass_01", .1, 0, 4, StaticSchwarzschildRay(.1, 8),
                     reference::Fate::Escape, true});
    cases.push_back({"schwarzschild_critical_plus_005", 1, 0, 40,
                     StaticSchwarzschildRay(1, 3 * std::sqrt(3.) + .05), reference::Fate::Escape,
                     true});
    auto moving = Ray(8, 1.1, {.5, .65, .57});
    moving.beta_forward = .08;
    moving.beta_up = .025;
    moving.beta_right = -.015;
    moving.aperture_right = -.017;
    moving.aperture_up = .021;
    cases.push_back({"kerr_moving_off_plane", 1, .7, 40, moving, reference::Fate::Escape, false});
    auto strong = Ray(12, 1.1, {-.6, .25, std::sqrt(.5775L)});
    cases.push_back({"kerr_0998_inner_turn", 1, .998, 40, strong, reference::Fate::Escape, true});
    auto retrograde = Ray(1200, 1.1, {-.6, -.2, std::sqrt(.60L)});
    cases.push_back(
        {"kerr_minus_0998_mass_100", 100, -99.8, 4000, retrograde, reference::Fate::Escape, true});
    cases.push_back({"kerr_extremal_capture", 1, 1, 40, Ray(8, 1.1, {-1, 0, 0}),
                     reference::Fate::Capture, false});
    cases.push_back({"kerr_minus_extremal_outward", 1, -1, 40, Ray(8, 1.1, {.5, .65, .57}),
                     reference::Fate::Escape, false});
    auto disk = Ray(50, std::numbers::pi / 3, {-std::sqrt(1 - .2L * .2L - .08L * .08L), .2, .08});
    cases.push_back(
        {"schwarzschild_first_disk", 1, 0, 200, disk, reference::Fate::Disk, false, 6, 20});
    disk.beta_forward = .03;
    disk.beta_up = -.01;
    disk.beta_right = .02;
    cases.push_back(
        {"kerr_0998_moving_first_disk", 1, .998, 200, disk, reference::Fate::Disk, false, 6, 20});
    return cases;
}

reference::Result IndependentTrace(const Case& c, const CameraRay& camera, Scalar step) {
    KerrSchildFamily metric(c.mass == 0 ? KerrSchildParams::Minkowski()
                                        : KerrSchildParams::Kerr(c.mass, c.spin));
    const auto launch = LaunchCameraRay(metric, c.spin, camera);
    if (!launch) throw std::runtime_error("unrepresented measured reference launch");
    return reference::Trace(*launch, c.mass, c.spin, c.boundary, step, c.disk_inner, c.disk_outer);
}

double MatrixError(const Matrix& actual, const Matrix& expected) {
    Scalar difference = 0, scale = 0;
    for (unsigned row = 0; row < 2; ++row)
        for (unsigned col = 0; col < 2; ++col) {
            difference = std::max(difference, std::abs(actual[row][col] - expected[row][col]));
            scale = std::max(scale, std::abs(expected[row][col]));
        }
    return static_cast<double>(scale > 0 ? difference / scale : difference);
}
double Angle(const Three& first, const Three& second) {
    const Three cross{first[1] * second[2] - first[2] * second[1],
                      first[2] * second[0] - first[0] * second[2],
                      first[0] * second[1] - first[1] * second[0]};
    return static_cast<double>(
        std::atan2(std::sqrt(reference::Dot(cross, cross)), reference::Dot(first, second)));
}
double FourError(const Four& first, const Four& second, double scale) {
    double worst = 0;
    for (unsigned i = 0; i < 4; ++i)
        worst = std::max(worst, static_cast<double>(std::abs(first[i] - second[i])) / scale);
    return worst;
}

struct Maps {
    Matrix finite{}, infinity{}, covariance{};
};
Maps IndependentMaps(const Case& c, const reference::Result& central, Scalar step, double spacing) {
    const Three n =
        reference::Unit({c.camera.direction(1), c.camera.direction(2), c.camera.direction(3)});
    const auto input_basis = reference::Basis(n);
    const auto finite_basis = reference::Basis(central.finite_direction);
    const auto infinity_basis =
        reference::Basis(central.infinity_direction == Three{} ? central.finite_direction
                                                               : central.infinity_direction);
    const auto geometry = reference::At(central.x, c.mass, c.spin, central.outgoing);
    const auto physical_screen = geometry.Screen(central.k);
    Matrix displacement{};
    Maps result;
    for (unsigned column = 0; column < 2; ++column) {
        std::array<reference::Result, 4> neighbours;
        const std::array<int, 4> offsets{-2, -1, 1, 2};
        for (unsigned index = 0; index < 4; ++index) {
            auto ray = c.camera;
            const Scalar angle = offsets[index] * spacing;
            for (unsigned axis = 0; axis < 3; ++axis)
                ray.direction(static_cast<int>(axis) + 1) = static_cast<double>(
                    std::cos(angle) * n[axis] + std::sin(angle) * input_basis[column][axis]);
            neighbours[index] = IndependentTrace(c, ray, step);
            if (neighbours[index].fate != central.fate)
                throw std::runtime_error("reference neighbour crosses the capture separatrix");
        }
        Three finite{}, infinity{};
        Four position{};
        const std::array<Scalar, 4> weights{1, -8, 8, -1};
        for (unsigned index = 0; index < 4; ++index) {
            for (unsigned axis = 0; axis < 3; ++axis) {
                finite[axis] +=
                    weights[index] * neighbours[index].finite_direction[axis] / (12 * spacing);
                infinity[axis] +=
                    weights[index] * neighbours[index].infinity_direction[axis] / (12 * spacing);
            }
            for (unsigned mu = 0; mu < 4; ++mu)
                position[mu] += weights[index] * neighbours[index].x[mu] / (12 * spacing);
        }
        for (unsigned row = 0; row < 2; ++row) {
            displacement[row][column] = geometry.Inner(physical_screen[row], position);
            if (central.fate == reference::Fate::Escape) {
                result.finite[row][column] = reference::Dot(finite_basis[row], finite);
                result.infinity[row][column] = reference::Dot(infinity_basis[row], infinity);
            }
        }
    }
    for (unsigned row = 0; row < 2; ++row)
        for (unsigned col = 0; col < 2; ++col)
            for (unsigned input = 0; input < 2; ++input)
                result.covariance[row][col] += displacement[row][input] * displacement[col][input];
    return result;
}

struct Witness {
    Case input;
    reference::Result fine;
    Maps maps;
    double endpoint_uncertainty = 0, map_uncertainty = 0;
    Scalar reference_step = 0;
    double reference_spacing = 0;
};

const std::vector<Witness>& Witnesses() {
    static const auto values = [] {
        std::vector<Witness> rows;
        for (const auto& c : Cases()) {
            // Critical scattering amplifies the reference's fourth-order
            // inner-flow error. Refine the expectation rather than reducing
            // its witness set or admitting a larger uncertainty budget.
            const Scalar coarse_step = c.name == "schwarzschild_critical_plus_005" ? .001L : .002L;
            const Scalar fine_step = coarse_step / 2;
            const double narrow_spacing = c.name == "schwarzschild_critical_plus_005" ? 5e-5 : 1e-4;
            const double wide_spacing = 2 * narrow_spacing;
            const auto coarse = IndependentTrace(c, c.camera, coarse_step);
            const auto fine = IndependentTrace(c, c.camera, fine_step);
            const auto coarse_maps = IndependentMaps(c, coarse, coarse_step, wide_spacing);
            const auto maps = IndependentMaps(c, fine, fine_step, narrow_spacing);
            const auto spacing = IndependentMaps(c, fine, fine_step, wide_spacing);
            const double scale = c.mass > 0 ? c.mass : 1;
            double uncertainty =
                std::max({FourError(coarse.x, fine.x, scale), FourError(coarse.k, fine.k, 1),
                          std::abs(static_cast<double>(coarse.affine - fine.affine)) / scale});
            if (fine.fate == reference::Fate::Escape)
                uncertainty = std::max(uncertainty,
                                       Angle(coarse.infinity_direction, fine.infinity_direction));
            const double map_uncertainty =
                std::max({MatrixError(coarse_maps.finite, maps.finite),
                          MatrixError(coarse_maps.infinity, maps.infinity),
                          MatrixError(coarse_maps.covariance, maps.covariance),
                          MatrixError(spacing.finite, maps.finite),
                          MatrixError(spacing.infinity, maps.infinity),
                          MatrixError(spacing.covariance, maps.covariance)});
            rows.push_back(
                {c, fine, maps, uncertainty, map_uncertainty, fine_step, narrow_spacing});
        }
        return rows;
    }();
    return values;
}

TracerConfig Control(const Case& c, unsigned refinement) {
    TracerConfig config;
    config.enable_disk = c.disk_inner > 0;
    if (config.enable_disk) {
        config.disk_inner = c.disk_inner;
        config.disk_outer = c.disk_outer;
    }
    config.enable_ray_bundles = true;
    config.bundle_point_source = true;
    config.escape_radius = static_cast<float>(c.boundary);
    config.max_steps = 20000;
    const float scale = static_cast<float>(c.mass > 0 ? c.mass : 1);
    const std::array<float, 3> cap{2, .5f, .125f};
    const std::array<float, 3> tolerance{5e-6f, 5e-8f, 5e-10f};
    config.integrator.initial_step = .1f * scale;
    config.integrator.max_step = cap[refinement] * scale;
    config.integrator.min_step = 1e-6f * scale;
    config.integrator.abs_tolerance = config.integrator.rel_tolerance = tolerance[refinement];
    return config;
}

struct Errors {
    double endpoint = 0, angle = 0, finite_map = 0, infinity_map = 0, covariance = 0;
};
Errors Check(const Witness& witness, const TraceResult& actual) {
    const auto& c = witness.input;
    const auto& expected = witness.fine;
    SCOPED_TRACE(c.name);
    EXPECT_EQ(expected.fate, c.fate);
    if (c.radial_turn) {
        EXPECT_GE(expected.radial_turns, 1U);
    }
    EXPECT_LE(witness.endpoint_uncertainty, kAcceptance * kReferenceFraction);
    EXPECT_LE(witness.map_uncertainty, kAcceptance * kReferenceFraction);
    EXPECT_FALSE(actual.cancelled);
    EXPECT_FALSE(actual.numerical_failure)
        << CoupledStepFailureName(actual.coupled_failure)
        << " termination=" << actual.integrator_termination << " attempts=" << actual.steps_taken;
    const auto outcome = c.fate == reference::Fate::Capture ? TraceResult::Outcome::Horizon
                         : c.fate == reference::Fate::Disk  ? TraceResult::Outcome::DiskHit
                                                            : TraceResult::Outcome::Escaped;
    EXPECT_EQ(actual.outcome, outcome);
    EXPECT_EQ(actual.terminal_chart, expected.outgoing
                                         ? TraceResult::TerminalChart::OutgoingKerrSchild
                                         : TraceResult::TerminalChart::MetricNative);
    EXPECT_TRUE(actual.final_tangent);
    EXPECT_TRUE(actual.beam.valid);
    if (actual.numerical_failure || !actual.final_tangent || !actual.beam.valid) return {};
    Four position{}, tangent{};
    for (unsigned mu = 0; mu < 4; ++mu) {
        position[mu] = actual.final_position(static_cast<int>(mu));
        tangent[mu] = (*actual.final_tangent)(static_cast<int>(mu));
        EXPECT_TRUE(std::isfinite(position[mu]));
        EXPECT_TRUE(std::isfinite(tangent[mu]));
        EXPECT_TRUE(std::isfinite(expected.x[mu]));
        EXPECT_TRUE(std::isfinite(expected.k[mu]));
    }
    const double scale = c.mass > 0 ? c.mass : 1;
    Errors error;
    error.endpoint =
        std::max({FourError(position, expected.x, scale), FourError(tangent, expected.k, 1),
                  std::abs(actual.affine_length - static_cast<double>(expected.affine)) / scale});
    EXPECT_LE(error.endpoint + witness.endpoint_uncertainty, kAcceptance);
    const double major =
        actual.beam.semi_major / static_cast<double>(Control(c, 0).bundle_angular_size);
    const double minor =
        actual.beam.semi_minor / static_cast<double>(Control(c, 0).bundle_angular_size);
    const double cosine = std::cos(actual.beam.orientation),
                 sine = std::sin(actual.beam.orientation);
    EXPECT_TRUE(std::isfinite(major));
    EXPECT_TRUE(std::isfinite(minor));
    EXPECT_TRUE(std::isfinite(actual.beam.orientation));
    const Matrix covariance{{{major * major * cosine * cosine + minor * minor * sine * sine,
                              (major * major - minor * minor) * cosine * sine},
                             {(major * major - minor * minor) * cosine * sine,
                              major * major * sine * sine + minor * minor * cosine * cosine}}};
    error.covariance = MatrixError(covariance, witness.maps.covariance);
    EXPECT_LE(error.covariance + witness.map_uncertainty, kAcceptance);
    if (expected.fate == reference::Fate::Disk) {
        EXPECT_EQ(actual.num_disk_crossings, 1);
        EXPECT_LE(std::abs(actual.final_position(3)) / scale, 1e-10);
        EXPECT_NEAR(actual.disk_radius,
                    static_cast<double>(reference::At(expected.x, c.mass, c.spin).radius),
                    scale * kAcceptance);
    }
    if (expected.fate == reference::Fate::Escape) {
        EXPECT_TRUE(actual.beam.finite_source_map);
        EXPECT_TRUE(actual.beam.infinity_source_map);
        if (!actual.beam.finite_source_map || !actual.beam.infinity_source_map) return error;
        Matrix finite{}, infinity{};
        for (unsigned row = 0; row < 2; ++row)
            for (unsigned col = 0; col < 2; ++col) {
                finite[row][col] = actual.beam.finite_source_map->jacobian[row][col];
                infinity[row][col] = actual.beam.infinity_source_map->map.jacobian[row][col];
                EXPECT_TRUE(std::isfinite(finite[row][col]));
                EXPECT_TRUE(std::isfinite(infinity[row][col]));
            }
        error.finite_map = MatrixError(finite, witness.maps.finite);
        error.infinity_map = MatrixError(infinity, witness.maps.infinity);
        Three sky{};
        for (unsigned axis = 0; axis < 3; ++axis)
            sky[axis] = actual.beam.infinity_source_map->map.direction[axis];
        error.angle = Angle(sky, expected.infinity_direction);
        Three finite_sky{}, finite_map_direction{};
        for (unsigned axis = 0; axis < 3; ++axis) {
            finite_sky[axis] = actual.final_direction(static_cast<int>(axis) + 1);
            finite_map_direction[axis] = actual.beam.finite_source_map->direction[axis];
        }
        error.angle = std::max({error.angle, Angle(finite_sky, expected.finite_direction),
                                Angle(finite_map_direction, expected.finite_direction)});
        EXPECT_LE(error.angle + witness.endpoint_uncertainty, kAcceptance);
        EXPECT_LE(error.finite_map + witness.map_uncertainty, kAcceptance);
        EXPECT_LE(error.infinity_map + witness.map_uncertainty, kAcceptance);
        EXPECT_NEAR(actual.beam.infinity_source_map->frequency,
                    -static_cast<double>(expected.energy), kAcceptance);
    } else {
        EXPECT_FALSE(actual.beam.infinity_source_map);
    }
    return error;
}

void Record(const Witness& witness, unsigned refinement, const Errors& error,
            const TraceResult& actual) {
    ::testing::Test::RecordProperty(
        witness.input.name + "_r" + std::to_string(refinement),
        std::format("endpoint_per_mass={:.17g};direction_angle_rad={:.17g};finite_map_relative={:."
                    "17g};infinity_map_relative={:.17g};beam_covariance_relative={:.17g};reference_"
                    "endpoint_gap={:.17g};reference_map_gap={:.17g};attempts={};reference_steps={};"
                    "radial_turns={};polar_turns={};reference_step_fraction={:.17g};reference_"
                    "angular_spacing={:.17g}",
                    error.endpoint, error.angle, error.finite_map, error.infinity_map,
                    error.covariance, witness.endpoint_uncertainty, witness.map_uncertainty,
                    actual.steps_taken, witness.fine.steps, witness.fine.radial_turns,
                    witness.fine.polar_turns, static_cast<double>(witness.reference_step),
                    witness.reference_spacing));
}

TEST(FullPathAcceptance, CpuIndependentCarterEventsMapsAndRefinement) {
    for (const auto& witness : Witnesses()) {
        KerrSchildFamily metric(
            witness.input.mass == 0
                ? KerrSchildParams::Minkowski()
                : KerrSchildParams::Kerr(witness.input.mass, witness.input.spin));
        std::array<Errors, 3> errors{};
        for (unsigned refinement = 0; refinement < 3; ++refinement) {
            SCOPED_TRACE(refinement);
            GeodesicTracer tracer(&metric, Control(witness.input, refinement));
            const auto result = tracer.Trace(witness.input.camera);
            errors[refinement] = Check(witness, result);
            Record(witness, refinement, errors[refinement], result);
        }
        // Coupled transport can already dominate the looser controller. A
        // roundoff/reference plateau is accepted and reported explicitly;
        // tightening the controller may not worsen a resolved observable.
        EXPECT_LE(errors[2].endpoint, errors[0].endpoint + 2 * witness.endpoint_uncertainty + 2e-8);
        EXPECT_LE(errors[2].infinity_map,
                  errors[0].infinity_map + 2 * witness.map_uncertainty + 2e-8);
    }
}

TEST(FullPathAcceptance, VulkanRetainedIndependentCarterEventsMapsAndRefinement) {
#if defined(SIRIUS_HAS_RETAINED_COMPUTE) && defined(SIRIUS_RETAINED_CAMERA_TEST_DIR)
    const auto inventory = EnumerateVulkanDevices();
    ASSERT_TRUE(inventory) << inventory.error().Description();
    if (inventory->empty()) GTEST_SKIP() << "no Vulkan device; numerical backend unqualified";
    const auto index = ResolveVulkanDeviceIndex(*inventory);
    ASSERT_TRUE(index) << index.error().Description();
    auto opened = CreateVulkanDevice(*index);
    ASSERT_TRUE(opened) << opened.error().Description();
    auto& device = **opened;
    RecordProperty("device", device.Info().name);
    RecordProperty("evidence_scope",
                   "actual retained camera/interval backend; shared host events and observables; "
                   "no throughput qualification");
    for (const bool wide : {false, true}) {
        if (wide && (!device.Info().supports_fp64 || !device.Info().rounds_fp64_to_nearest)) {
            RecordProperty("fp64_evidence", "unsupported; unqualified");
            continue;
        }
        auto compute = RetainedCompute::Create(device, Witnesses().size(), wide);
        ASSERT_TRUE(compute) << compute.error().Description();
        RetainedTraceExecutor executor(**compute);
        std::vector<std::array<Errors, 3>> convergence(Witnesses().size());
        for (unsigned refinement = 0; refinement < 3; ++refinement) {
            std::vector<std::future<TraceResult>> tasks;
            for (const auto& witness : Witnesses()) {
                tasks.push_back(std::async(std::launch::async, [&executor, &witness, refinement] {
                    KerrSchildFamily metric(
                        witness.input.mass == 0
                            ? KerrSchildParams::Minkowski()
                            : KerrSchildParams::Kerr(witness.input.mass, witness.input.spin));
                    GeodesicTracer tracer(&metric, Control(witness.input, refinement));
                    tracer.SetStepExecutor(&executor);
                    return tracer.Trace(witness.input.camera);
                }));
            }
            for (std::size_t row = 0; row < tasks.size(); ++row) {
                SCOPED_TRACE(wide ? "retained_fp64" : "retained_fp32");
                SCOPED_TRACE(refinement);
                const auto result = tasks[row].get();
                const auto errors = Check(Witnesses()[row], result);
                convergence[row][refinement] = errors;
                auto named = Witnesses()[row];
                named.input.name += wide ? "_fp64" : "_fp32";
                Record(named, refinement, errors, result);
            }
        }
        for (std::size_t row = 0; row < Witnesses().size(); ++row) {
            SCOPED_TRACE(Witnesses()[row].input.name);
            EXPECT_LE(
                convergence[row][2].endpoint,
                convergence[row][0].endpoint + 2 * Witnesses()[row].endpoint_uncertainty + 2e-8);
            EXPECT_LE(
                convergence[row][2].infinity_map,
                convergence[row][0].infinity_map + 2 * Witnesses()[row].map_uncertainty + 2e-8);
        }
        EXPECT_FALSE(executor.Error());
        EXPECT_GT(executor.Statistics().camera_batches, 0U);
        EXPECT_GT(executor.Statistics().interval_batches, 0U);
        RecordProperty(wide ? "fp64_evidence" : "fp32_evidence",
                       "full independent numerical corpus executed");
    }
#else
    GTEST_SKIP() << "retained compute kernels unavailable; Vulkan numerical backend unqualified";
#endif
}
}  // namespace
