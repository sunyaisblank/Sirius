// Complete production camera-to-event rays against an independently evolved
// Carter flow. This is numerical acceptance, with no image or throughput claim.
#include "sirius/backend/cpu/geodesic_tracer.h"
#include "sirius/core/metrics/cpu_metric_factory.h"

#include <gtest/gtest.h>

#include "support/separated_geodesic_reference.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <format>
#include <future>
#include <limits>
#include <numbers>
#include <optional>
#include <string>
#include <vector>

#if defined(SIRIUS_HAS_RETAINED_COMPUTE) && defined(SIRIUS_TEST_HAS_RETAINED_CAMERA)
#include "sirius/backend/device.h"
#include "sirius/backend/retained_compute.h"
#include "sirius/backend/retained_trace_executor.h"
#include "sirius/core/metrics/outgoing_kerr_schild.h"
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
    double charge_ratio = 0;
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

CameraRay StaticChargedRay(double mass, double charge_ratio, double impact) {
    constexpr double radius = 12;
    const double h = 2 / radius - charge_ratio * charge_ratio / (radius * radius);
    const double transverse = impact * std::sqrt(1 - h) / radius;
    auto ray = Ray(radius * mass, std::numbers::pi / 2,
                   {-std::sqrt(1 - transverse * transverse), 0, transverse});
    ray.beta_forward = -h;
    return ray;
}

std::vector<Case> ChargedCases() {
    std::vector<Case> cases;
    cases.push_back({"rn_charge_06_capture", 1, 0, 40, StaticChargedRay(1, .6, 2),
                     reference::Fate::Capture, false, 0, 0, .6});
    cases.push_back({"rn_charge_09_turn_mass_01", .1, 0, 4, StaticChargedRay(.1, .9, 8),
                     reference::Fate::Escape, true, 0, 0, .9});
    constexpr double charge = .999;
    const double photon_radius = (3 + std::sqrt(9 - 8 * charge * charge)) / 2;
    const double photon_f =
        1 - 2 / photon_radius + charge * charge / (photon_radius * photon_radius);
    const double critical_impact = photon_radius / std::sqrt(photon_f);
    cases.push_back({"rn_charge_0999_critical_plus_005_mass_100", 100, 0, 4000,
                     StaticChargedRay(100, charge, critical_impact + .05), reference::Fate::Escape,
                     true, 0, 0, charge});
    auto moving = Ray(8, 1.1, {.5, .65, .57});
    moving.beta_forward = .08;
    moving.beta_up = .025;
    moving.beta_right = -.015;
    moving.aperture_right = -.017;
    moving.aperture_up = .021;
    cases.push_back({"kn_spin_07_charge_06_moving_pupil", 1, .7, 40, moving,
                     reference::Fate::Escape, false, 0, 0, .6});
    cases.push_back({"kn_spin_07_charge_07_inner_turn", 1, .7, 40,
                     Ray(12, 1.1, {-.6, .25, std::sqrt(.5775L)}), reference::Fate::Escape, true, 0,
                     0, .7});
    cases.push_back({"kn_spin_02_charge_097_capture", 1, .2, 40, Ray(8, 1.1, {-1, 0, 0}),
                     reference::Fate::Capture, false, 0, 0, .97});
    return cases;
}

reference::Result IndependentTrace(const Case& c, const CameraRay& camera, Scalar step) {
    const double charge = c.mass * c.charge_ratio;
    KerrSchildFamily metric(c.mass == 0 ? KerrSchildParams::Minkowski()
                                        : KerrSchildParams::KerrNewman(c.mass, c.spin, charge));
    const auto launch = LaunchCameraRay(metric, c.spin, camera);
    if (!launch) throw std::runtime_error("unrepresented measured reference launch");
    return reference::Trace(*launch, c.mass, c.spin, c.boundary, step, c.disk_inner, c.disk_outer,
                            charge);
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
    const auto geometry =
        reference::At(central.x, c.mass, c.spin, central.outgoing, c.mass * c.charge_ratio);
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

std::vector<Witness> BuildWitnesses(const std::vector<Case>& cases) {
    std::vector<Witness> rows;
    for (const auto& c : cases) {
        // Critical scattering amplifies the reference's fourth-order
        // inner-flow error. Refine the expectation rather than reducing
        // its witness set or admitting a larger uncertainty budget.
        const bool critical = c.name.find("critical") != std::string::npos;
        const Scalar coarse_step = critical ? (c.charge_ratio != 0 ? .0005L : .001L) : .002L;
        const Scalar fine_step = coarse_step / 2;
        const double narrow_spacing = critical ? 5e-5 : 1e-4;
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
            uncertainty =
                std::max(uncertainty, Angle(coarse.infinity_direction, fine.infinity_direction));
        const double map_uncertainty = std::max(
            {MatrixError(coarse_maps.finite, maps.finite),
             MatrixError(coarse_maps.infinity, maps.infinity),
             MatrixError(coarse_maps.covariance, maps.covariance),
             MatrixError(spacing.finite, maps.finite), MatrixError(spacing.infinity, maps.infinity),
             MatrixError(spacing.covariance, maps.covariance)});
        rows.push_back({c, fine, maps, uncertainty, map_uncertainty, fine_step, narrow_spacing});
    }
    return rows;
}

const std::vector<Witness>& Witnesses() {
    static const auto values = BuildWitnesses(Cases());
    return values;
}

const std::vector<Witness>& ChargedWitnesses() {
    static const auto values = BuildWitnesses(ChargedCases());
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
        EXPECT_EQ(actual.beam.infinity_source_map.has_value(), c.charge_ratio == 0);
        if (!actual.beam.finite_source_map ||
            (c.charge_ratio == 0 && !actual.beam.infinity_source_map))
            return error;
        Matrix finite{}, infinity{};
        for (unsigned row = 0; row < 2; ++row)
            for (unsigned col = 0; col < 2; ++col) {
                finite[row][col] = actual.beam.finite_source_map->jacobian[row][col];
                if (c.charge_ratio == 0)
                    infinity[row][col] = actual.beam.infinity_source_map->map.jacobian[row][col];
                EXPECT_TRUE(std::isfinite(finite[row][col]));
                EXPECT_TRUE(std::isfinite(infinity[row][col]));
            }
        error.finite_map = MatrixError(finite, witness.maps.finite);
        if (c.charge_ratio == 0) {
            error.infinity_map = MatrixError(infinity, witness.maps.infinity);
            Three sky{};
            for (unsigned axis = 0; axis < 3; ++axis)
                sky[axis] = actual.beam.infinity_source_map->map.direction[axis];
            error.angle = Angle(sky, expected.infinity_direction);
            EXPECT_NEAR(actual.beam.infinity_source_map->frequency,
                        -static_cast<double>(expected.energy), kAcceptance);
        }
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

#if defined(SIRIUS_HAS_RETAINED_COMPUTE) && defined(SIRIUS_TEST_HAS_RETAINED_CAMERA)
struct RayInvariants {
    Scalar energy, angular_momentum, carter, normalized_null;
    Scalar radius, radial_velocity, absolute_sine;
};

// Contract the represented physical view with an independent analytic metric.
// E=-p_t and Lz=x*p_y-y*p_x retain their signs for the past-directed tangent;
// p_theta=Sigma*dtheta/dlambda gives Q=p_theta^2+cos(theta)^2
// *(Lz^2/sin(theta)^2-a^2*E^2). These Killing quantities are unchanged by the
// radial time/azimuth chart shift. No production metric or BL chart map is used.
std::optional<RayInvariants> IndependentInvariants(const Vec4& position, const Vec4& tangent,
                                                   double mass, double spin, bool outgoing) {
    Four x{}, k{}, p{};
    for (unsigned mu = 0; mu < 4; ++mu) {
        x[mu] = position(static_cast<int>(mu));
        k[mu] = tangent(static_cast<int>(mu));
        if (!std::isfinite(x[mu]) || !std::isfinite(k[mu])) return std::nullopt;
    }
    const auto geometry = reference::At(x, mass, spin, outgoing);
    const Scalar r = geometry.radius, sine = std::sin(geometry.theta),
                 cosine = std::cos(geometry.theta), a = spin;
    if (!(std::abs(sine) > 1e-5L)) return std::nullopt;
    Scalar null_defect = 0, null_scale = 0;
    for (unsigned mu = 0; mu < 4; ++mu)
        for (unsigned nu = 0; nu < 4; ++nu) {
            p[mu] += geometry.metric[mu][nu] * k[nu];
            const Scalar term = geometry.metric[mu][nu] * k[mu] * k[nu];
            null_defect += term;
            null_scale += std::abs(term);
        }
    if (!(null_scale > 0)) return std::nullopt;
    const Scalar sigma = r * r + a * a * cosine * cosine;
    const Scalar radial =
        (r * (x[1] * k[1] + x[2] * k[2]) + (r * r + a * a) * x[3] * k[3] / r) / sigma;
    const Scalar polar_momentum = sigma * (cosine * radial - k[3]) / (r * sine);
    const Scalar energy = -p[0], angular = x[1] * p[2] - x[2] * p[1];
    const Scalar carter =
        polar_momentum * polar_momentum +
        cosine * cosine * (angular * angular / (sine * sine) - a * a * energy * energy);
    RayInvariants result{energy, angular, carter,        std::abs(null_defect) / null_scale,
                         r,      radial,  std::abs(sine)};
    for (const Scalar value :
         {result.energy, result.angular_momentum, result.carter, result.normalized_null,
          result.radius, result.radial_velocity, result.absolute_sine})
        if (!std::isfinite(value)) return std::nullopt;
    return result;
}

class ConservationObserver final : public TraceStepExecutor {
  public:
    struct Measurements {
        std::optional<RayInvariants> launch;
        std::array<Scalar, 3> maximum_relative_drift{};
        Scalar maximum_normalized_null = 0;
        Scalar minimum_radius = std::numeric_limits<Scalar>::infinity();
        Scalar minimum_absolute_sine = 1;
        unsigned accepted_intervals = 0, successful_candidates = 0, rejected_candidates = 0;
        unsigned native_samples = 0, outgoing_samples = 0, radial_turns = 0;
        bool finite = true;
    } measured;

    ConservationObserver(RetainedTraceExecutor& executor, double mass, double spin)
        : executor_(executor), mass_(mass), spin_(spin) {}

    void BeginTrace() override { executor_.BeginTrace(); }
    void EndTrace() override { executor_.EndTrace(); }
    std::optional<CameraLaunch> Launch(IMetric& metric, double spin,
                                       const CameraRay& camera) override {
        auto launch = executor_.Launch(metric, spin, camera);
        if (launch) {
            measured.launch = IndependentInvariants(launch->position, launch->tangent, mass_, spin_,
                                                    IsOutgoing(metric));
            Observe(launch->position, launch->tangent, IsOutgoing(metric));
        }
        return launch;
    }
    bool Step(Lightray& ray, IMetric& metric, const IntegratorConfig& config,
              Rk45CoupledState& coupled, Rk45CoupledComparison& comparison) override {
        // Step returns a private device candidate. Only the next input proves
        // that the host committed it, after event/source checks and any clipping.
        // The first input additionally checks the native-to-outgoing launch join.
        if (pending_ || !saw_initial_input_) {
            Observe(ray.position, ray.velocity, IsOutgoing(metric));
            if (pending_) ++measured.accepted_intervals;
            pending_ = false;
            saw_initial_input_ = true;
        }
        const bool success = executor_.Step(ray, metric, config, coupled, comparison);
        if (success) {
            pending_ = true;
            ++measured.successful_candidates;
        }
        return success;
    }
    void RejectLastInterval() override {
        // A rejected success contributes no sample or maximum. Failed device
        // attempts can also call this hook, with no pending candidate to remove.
        if (pending_) ++measured.rejected_candidates;
        pending_ = false;
        executor_.RejectLastInterval();
    }
    bool ObservePhysicalTerminal(const TraceResult& result) {
        if (!pending_ || result.cancelled || result.numerical_failure || !result.final_tangent ||
            result.outcome != TraceResult::Outcome::Escaped)
            return false;
        // The final candidate has no following Step. Sample the actual localized
        // published event in its declared chart, rather than the device overshoot.
        Observe(result.final_position, *result.final_tangent,
                result.terminal_chart == TraceResult::TerminalChart::OutgoingKerrSchild);
        ++measured.accepted_intervals;
        pending_ = false;
        return true;
    }

  private:
    bool IsOutgoing(IMetric& metric) {
        if (dynamic_cast<OutgoingKerrSchild*>(&metric)) return true;
        if (!dynamic_cast<KerrSchildFamily*>(&metric)) measured.finite = false;
        return false;
    }
    void Observe(const Vec4& position, const Vec4& tangent, bool outgoing) {
        const auto sample = IndependentInvariants(position, tangent, mass_, spin_, outgoing);
        if (!sample || !measured.launch || !(std::abs(measured.launch->energy) > 1e-8L) ||
            !(std::abs(measured.launch->angular_momentum) > 1e-8L) ||
            !(std::abs(measured.launch->carter) > 1e-8L)) {
            measured.finite = false;
            return;
        }
        const std::array<Scalar, 3> initial{
            measured.launch->energy, measured.launch->angular_momentum, measured.launch->carter};
        const std::array<Scalar, 3> actual{sample->energy, sample->angular_momentum,
                                           sample->carter};
        for (unsigned i = 0; i < 3; ++i)
            measured.maximum_relative_drift[i] =
                std::max(measured.maximum_relative_drift[i],
                         std::abs((actual[i] - initial[i]) / initial[i]));
        measured.maximum_normalized_null =
            std::max(measured.maximum_normalized_null, sample->normalized_null);
        measured.minimum_radius = std::min(measured.minimum_radius, sample->radius);
        measured.minimum_absolute_sine =
            std::min(measured.minimum_absolute_sine, sample->absolute_sine);
        if (previous_radial_ && *previous_radial_ < 0 && sample->radial_velocity > 0)
            ++measured.radial_turns;
        previous_radial_ = sample->radial_velocity;
        if (outgoing)
            ++measured.outgoing_samples;
        else
            ++measured.native_samples;
    }
    RetainedTraceExecutor& executor_;
    double mass_, spin_;
    bool pending_ = false, saw_initial_input_ = false;
    std::optional<Scalar> previous_radial_;
};
#endif

std::array<Scalar, 2> QuadratureShift(Scalar radius, Scalar mass, Scalar spin, Scalar charge,
                                      unsigned panels) {
    const Scalar horizon = mass + std::sqrt(mass * mass - spin * spin - charge * charge);
    const Scalar anchor = 2 * horizon;
    const Scalar step = (radius - anchor) / panels;
    std::array<Scalar, 2> integral{};
    for (unsigned index = 0; index <= panels; ++index) {
        const Scalar r = anchor + step * index;
        const Scalar delta = r * r - 2 * mass * r + spin * spin + charge * charge;
        const Scalar weight = index == 0 || index == panels ? 1 : index % 2 == 0 ? 2 : 4;
        integral[0] += weight * -2 * (2 * mass * r - charge * charge) / delta;
        integral[1] += weight * (-2 * spin / delta + 2 * spin / (r * r + spin * spin));
    }
    for (auto& value : integral) value *= step / 3;
    return integral;
}

TEST(ChargedReference, ExteriorChartShiftMatchesIndependentQuadrature) {
    for (const auto parameters :
         {std::array<Scalar, 3>{1, 0, .6L}, {1, .7L, .6L}, {1, .2L, .97L}}) {
        const auto [mass, spin, charge] = parameters;
        const Scalar radius = 12 * mass;
        const reference::detail::Constants constants{mass, spin, 0, 0, 0, charge};
        const auto analytic = reference::detail::Shift(radius, constants);
        const auto coarse = QuadratureShift(radius, mass, spin, charge, 2048);
        const auto fine = QuadratureShift(radius, mass, spin, charge, 4096);
        for (unsigned component = 0; component < 2; ++component) {
            EXPECT_LE(std::abs(coarse[component] - fine[component]), 1e-10L);
            EXPECT_LE(std::abs(analytic[component] - fine[component]), 1e-10L);
        }
    }
}

TEST(ChargedReference, SphericalRadialCaptureMatchesExactAffineAndTangent) {
    for (const auto parameters : {std::array<double, 2>{.1, .9}, {1, .6}, {100, .999}}) {
        const auto [mass, ratio] = parameters;
        const double charge = mass * ratio;
        KerrSchildFamily metric(KerrSchildParams::ReissnerNordstrom(mass, charge));
        const auto launch = LaunchCameraRay(metric, 0, StaticChargedRay(mass, ratio, 0));
        ASSERT_TRUE(launch);
        const Scalar radius = 12 * mass;
        const Scalar horizon = mass + std::sqrt(Scalar(mass) * mass - Scalar(charge) * charge);
        const Scalar h = 2 * mass / radius - Scalar(charge) * charge / (radius * radius);
        Scalar radial_velocity = 0;
        Three radial{};
        for (unsigned axis = 0; axis < 3; ++axis) {
            radial[axis] = launch->position(static_cast<int>(axis) + 1) / radius;
            radial_velocity += radial[axis] * launch->tangent(static_cast<int>(axis) + 1);
        }
        const Scalar speed = (-1 + h) * launch->tangent(0) + h * radial_velocity;
        ASSERT_GT(speed, 0);
        const auto expected_shift = QuadratureShift(radius, mass, 0, charge, 4096);
        const auto result = reference::Trace(*launch, mass, 0, 40 * mass, .0005L, 0, 0, charge);
        ASSERT_EQ(result.fate, reference::Fate::Capture);
        EXPECT_TRUE(result.outgoing);
        EXPECT_EQ(result.infinity_direction, Three{});
        EXPECT_LE(std::abs(result.affine - (radius - horizon) / speed) / mass, 1e-8L);
        EXPECT_LE(
            std::abs(result.x[0] - launch->position(0) - expected_shift[0] - (horizon - radius)) /
                mass,
            1e-8L);
        EXPECT_LE(std::abs(result.k[0] + speed), 1e-8L);
        for (unsigned axis = 0; axis < 3; ++axis) {
            EXPECT_LE(std::abs(result.x[axis + 1] - horizon * radial[axis]) / mass, 1e-8L);
            EXPECT_LE(std::abs(result.k[axis + 1] + speed * radial[axis]), 1e-8L);
        }
    }
}

TEST(FullPathAcceptance, CpuChargedFiniteEventsMapsAndRefinement) {
    // The direct coupled tracer monitors Jacobi columns for charged rays.
    // Public charged beam rendering and vacuum infinity transfer still decline.
    for (const auto& witness : ChargedWitnesses()) {
        const auto& c = witness.input;
        MetricConstructionParameters parameters;
        parameters.mass = c.mass;
        parameters.dimensionless_spin = c.spin / c.mass;
        parameters.dimensionless_charge = c.charge_ratio;
        auto metric = CreateCpuMetric(
            c.spin == 0 ? MetricId::ReissnerNordstrom : MetricId::KerrNewman, parameters);
        ASSERT_TRUE(metric) << c.name;
        std::array<Errors, 3> errors{};
        for (unsigned refinement = 0; refinement < 3; ++refinement) {
            SCOPED_TRACE(refinement);
            GeodesicTracer tracer(metric.get(), Control(c, refinement));
            const auto result = tracer.Trace(c.camera);
            errors[refinement] = Check(witness, result);
            Record(witness, refinement, errors[refinement], result);
        }
        EXPECT_LE(errors[2].endpoint, errors[0].endpoint + 2 * witness.endpoint_uncertainty + 2e-8);
        EXPECT_LE(errors[2].finite_map, errors[0].finite_map + 2 * witness.map_uncertainty + 2e-8);
    }
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
#if defined(SIRIUS_HAS_RETAINED_COMPUTE) && defined(SIRIUS_TEST_HAS_RETAINED_CAMERA)
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

        if (device.Info().kind != DeviceKind::kIntegratedGpu &&
            device.Info().kind != DeviceKind::kDiscreteGpu)
            continue;

        // Compare the same mixed physical rays at both renderer capacities.
        // The original independent witnesses judge every row, including capture,
        // disk, turning and moving-pupil events. A start gate supplies real
        // concurrent demand; the histogram records the batches actually executed.
        constexpr std::size_t cohort_rows = 128;
        for (const std::size_t capacity : {64U, 128U}) {
            SCOPED_TRACE(capacity);
            const auto resident_before = device.BufferAllocationBytes();
            const auto allocation = RetainedCompute::RequiredAllocationBytes(device, capacity);
            ASSERT_TRUE(allocation) << allocation.error().Description();
            const auto creation_started = std::chrono::steady_clock::now();
            auto cohort_compute = RetainedCompute::Create(device, capacity, wide);
            ASSERT_TRUE(cohort_compute) << cohort_compute.error().Description();
            const auto creation_ms = std::chrono::duration<double, std::milli>(
                                         std::chrono::steady_clock::now() - creation_started)
                                         .count();
            const auto resident = device.BufferAllocationBytes() - resident_before;
            EXPECT_EQ(resident, *allocation);
            RetainedTraceExecutor cohort_executor(**cohort_compute, {}, 1000);
            std::vector<std::future<TraceResult>> tasks;
            tasks.reserve(cohort_rows);
            // Destroy the promise before the futures if thread creation fails:
            // waiting workers then wake instead of blocking future destruction.
            std::promise<void> launch_signal;
            const auto launch_gate = launch_signal.get_future().share();
            const auto started = std::chrono::steady_clock::now();
            for (std::size_t row = 0; row < cohort_rows; ++row) {
                tasks.push_back(
                    std::async(std::launch::async, [&cohort_executor, launch_gate, row] {
                        const auto& witness = Witnesses()[row % Witnesses().size()];
                        KerrSchildFamily metric(
                            witness.input.mass == 0
                                ? KerrSchildParams::Minkowski()
                                : KerrSchildParams::Kerr(witness.input.mass, witness.input.spin));
                        GeodesicTracer tracer(&metric, Control(witness.input, 0));
                        tracer.SetStepExecutor(&cohort_executor);
                        launch_gate.wait();
                        return tracer.Trace(witness.input.camera);
                    }));
            }
            launch_signal.set_value();
            for (std::size_t row = 0; row < tasks.size(); ++row) {
                SCOPED_TRACE(row);
                Check(Witnesses()[row % Witnesses().size()], tasks[row].get());
            }
            EXPECT_FALSE(cohort_executor.Error());
            const auto timing = cohort_executor.Statistics();
            EXPECT_LE(timing.maximum_batch_rows, capacity);
            const auto prefix = std::format("fanout_{}_{}", wide ? "fp64" : "fp32", capacity);
            RecordProperty(prefix + "_host_wall_ms",
                           std::format("{:.17g}", std::chrono::duration<double, std::milli>(
                                                      std::chrono::steady_clock::now() - started)
                                                      .count()));
            RecordProperty(prefix + "_creation_ms", std::format("{:.17g}", creation_ms));
            RecordProperty(prefix + "_resident_bytes", std::to_string(resident));
            RecordProperty(prefix + "_submission_guard_ms", "1000");
            RecordProperty(prefix + "_larger_batches_executed",
                           timing.maximum_batch_rows > 64 ? "yes" : "no");
            RecordProperty(prefix + "_full_batches", std::to_string(timing.full_batches));
            RecordProperty(prefix + "_maximum_rows", std::to_string(timing.maximum_batch_rows));
            RecordProperty(prefix + "_accepted_intervals",
                           std::to_string(timing.accepted_intervals));
            RecordProperty(prefix + "_rejected_intervals",
                           std::to_string(timing.rejected_intervals));
            RecordProperty(prefix + "_batch_subdivisions",
                           std::to_string(timing.batch_subdivisions));
            RecordProperty(prefix + "_safety_reductions", std::to_string(timing.safety_fallbacks));
            RecordProperty(prefix + "_paired_projection_retries",
                           std::to_string(timing.paired_projection_retries));
            RecordProperty(prefix + "_execute_ms", std::format("{:.17g}", timing.execute_ms));
            std::string histogram;
            for (std::size_t rows = 1; rows < timing.batch_row_counts.size(); ++rows) {
                if (timing.batch_row_counts[rows] != 0)
                    histogram += std::format("{}:{};", rows, timing.batch_row_counts[rows]);
            }
            RecordProperty(prefix + "_batch_histogram", histogram);
            const auto stages = (*cohort_compute)->Statistics();
            for (std::size_t stage = 0; stage < stages.size(); ++stage) {
                const auto& stat = stages[stage];
                RecordProperty(
                    prefix + "_" +
                        RetainedCompute::StageName(
                            static_cast<RetainedCompute::KernelStage>(stage)),
                    std::format(
                        "submissions={};max_submit_wait_ms={:.17g};submit_wait_ms={:.17g};"
                        "pipeline_setup_ms={:.17g};command_setup_ms={:.17g};cleanup_ms={:.17g};"
                        "dispatch_total_ms={:.17g};write_ms={:.17g};read_ms={:.17g};"
                        "pipeline_creations={};target_overshoots={}",
                        stat.submissions, stat.maximum_submit_wait_ms, stat.submit_wait_ms,
                        stat.pipeline_setup_ms, stat.command_setup_ms, stat.cleanup_ms,
                        stat.dispatch_total_ms, stat.write_buffer_ms, stat.read_buffer_ms,
                        stat.pipeline_creations, stat.target_overshoots));
            }
        }
    }
#else
    GTEST_SKIP() << "retained compute kernels unavailable; Vulkan numerical backend unqualified";
#endif
}

TEST(FullPathAcceptance, VulkanRetainedNearExtremalRayConservesIndependentInvariants) {
#if defined(SIRIUS_HAS_RETAINED_COMPUTE) && defined(SIRIUS_TEST_HAS_RETAINED_CAMERA)
    const auto& witnesses = Witnesses();
    const auto selected = std::find_if(witnesses.begin(), witnesses.end(), [](const Witness& item) {
        return item.input.name == "kerr_0998_inner_turn";
    });
    ASSERT_NE(selected, witnesses.end());
    const auto& witness = *selected;
    const auto& c = witness.input;
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
                   "actual retained camera and intervals; committed physical host views and "
                   "localized terminal; independently contracted E/Lz/Q and normalized null; "
                   "one off-plane inward-turn-escape ray per supported product mode; no "
                   "unprojected-stage or uniform-domain claim");
    RecordProperty("reference_scalar_mantissa_bits", std::numeric_limits<Scalar>::digits);
    unsigned executed_modes = 0;
    for (const bool wide : {false, true}) {
        SCOPED_TRACE(wide ? "retained_binary64_products" : "retained_binary32_products");
        if (wide && (!device.Info().supports_fp64 || !device.Info().rounds_fp64_to_nearest)) {
            RecordProperty("binary64_products_evidence", "unsupported; unqualified");
            continue;
        }
        auto compute = RetainedCompute::Create(device, 1, wide);
        ASSERT_TRUE(compute) << compute.error().Description();
        RetainedTraceExecutor executor(**compute);
        ConservationObserver observer(executor, c.mass, c.spin);
        KerrSchildFamily metric(KerrSchildParams::Kerr(c.mass, c.spin));
        // Preserve the coarsest existing coupled controller and its full ray,
        // beam and source-map contract. This gate adds interval maxima only.
        GeodesicTracer tracer(&metric, Control(c, 0));
        tracer.SetStepExecutor(&observer);
        const auto result = tracer.Trace(c.camera);
        ASSERT_FALSE(result.cancelled);
        ASSERT_FALSE(result.numerical_failure) << CoupledStepFailureName(result.coupled_failure)
                                               << " termination=" << result.integrator_termination;
        ASSERT_EQ(result.outcome, TraceResult::Outcome::Escaped);
        ASSERT_EQ(result.terminal_chart, TraceResult::TerminalChart::MetricNative);
        ASSERT_TRUE(observer.ObservePhysicalTerminal(result));
        const auto& measured = observer.measured;
        ASSERT_TRUE(measured.finite);
        ASSERT_TRUE(measured.launch);
        EXPECT_GT(measured.accepted_intervals, 2U);
        EXPECT_EQ(measured.accepted_intervals,
                  measured.successful_candidates - measured.rejected_candidates);
        EXPECT_EQ(measured.native_samples, 2U);  // Measured launch and published terminal.
        EXPECT_EQ(measured.outgoing_samples, measured.accepted_intervals);
        EXPECT_GT(measured.radial_turns, 0U);
        EXPECT_LT(measured.minimum_radius, measured.launch->radius);
        EXPECT_GT(measured.minimum_radius, c.mass + std::sqrt(c.mass * c.mass - c.spin * c.spin));
        EXPECT_GT(measured.minimum_absolute_sine, 1e-5L);
        for (const Scalar drift : measured.maximum_relative_drift) EXPECT_LE(drift, kAcceptance);
        EXPECT_LE(measured.maximum_normalized_null, 1e-6L);
        const auto error = Check(witness, result);
        auto named = witness;
        named.input.name +=
            wide ? "_conservation_binary64_products" : "_conservation_binary32_products";
        Record(named, 0, error, result);
        RecordProperty(
            wide ? "binary64_products_conservation" : "binary32_products_conservation",
            std::format(
                "state=expanded_binary32;products={};E_initial={:.17g};Lz_initial={:.17g};"
                "Q_initial={:.17g};max_E_relative={:.17g};max_Lz_relative={:.17g};"
                "max_Q_relative={:.17g};max_abs_gkk_over_sum_abs_terms={:.17g};"
                "accepted_intervals={};successful_candidates={};rolled_back_candidates={};"
                "native_samples={};outgoing_samples={};radial_turns={};min_r_per_mass={:.17g};"
                "min_abs_sin_theta={:.17g};terminal=physical_escape;affine_per_mass={:.17g}",
                wide ? "binary64" : "binary32", static_cast<double>(measured.launch->energy),
                static_cast<double>(measured.launch->angular_momentum),
                static_cast<double>(measured.launch->carter),
                static_cast<double>(measured.maximum_relative_drift[0]),
                static_cast<double>(measured.maximum_relative_drift[1]),
                static_cast<double>(measured.maximum_relative_drift[2]),
                static_cast<double>(measured.maximum_normalized_null), measured.accepted_intervals,
                measured.successful_candidates, measured.rejected_candidates,
                measured.native_samples, measured.outgoing_samples, measured.radial_turns,
                static_cast<double>(measured.minimum_radius / c.mass),
                static_cast<double>(measured.minimum_absolute_sine),
                result.affine_length / c.mass));
        EXPECT_FALSE(executor.Error());
        EXPECT_GT(executor.Statistics().camera_batches, 0U);
        EXPECT_GT(executor.Statistics().interval_batches, 0U);
        EXPECT_GT(executor.Statistics().reused_phases, 0U);
        ++executed_modes;
    }
    RecordProperty("executed_product_modes", static_cast<int>(executed_modes));
#else
    GTEST_SKIP() << "retained compute kernels unavailable; Vulkan numerical backend unqualified";
#endif
}
}  // namespace
