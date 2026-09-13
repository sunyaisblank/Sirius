// Direct coupled interval controls and actual CPU tracing. No render/device path.
#include "sirius/backend/cpu/geodesic_tracer.h"
#include "sirius/core/disk/novikov_thorne_disk.h"
#include "sirius/core/metrics/kerr_schild_family.h"
#include "sirius/core/metrics/outgoing_kerr_schild.h"
#include "sirius/core/observer_frame.h"
#include "sirius/core/trace_boundary.h"

#include <gtest/gtest.h>

#include <iomanip>
#include <limits>
#include <numbers>
#include <sstream>

namespace {
using namespace sirius::core;
using namespace sirius::backend;

Lightray FlatRay() {
    Lightray ray{};
    ray.position(1) = 5.0;
    ray.velocity(0) = -1.0;
    ray.velocity(1) = 1.0;
    ray.ku_uobsu = 1.0;
    ray.step_size = 2.0f;
    return ray;
}
Rk45CoupledState Columns() {
    Rk45CoupledState state;
    state.length_scale = 5.0;
    state.frequency_scale = 1.0;
    state.tolerance = 1.0e-9;
    state.stationary = true;
    state.variations[0].derivative(2) = 1.0;
    state.variations[1].derivative(3) = 1.0;
    state.variations[2].displacement(2) = 5.0;
    state.variations[3].displacement(3) = 5.0;
    return state;
}
IntegratorConfig Control() {
    IntegratorConfig config;
    config.max_step = config.initial_step = 2.0f;
    return config;
}

TEST(CoupledTransport, FlatFourColumnsUseActualProjectedAndInteriorTrials) {
    KerrSchildFamily metric(KerrSchildParams::Minkowski());
    auto ray = FlatRay();
    auto columns = Columns();
    Rk45CoupledComparison comparison;
    auto config = Control();
    config.min_step = ray.step_size;
    ASSERT_TRUE(Geodesic::IntegrateStepRk45(ray, &metric, config, &columns, &comparison));
    EXPECT_NEAR(ray.position(1), 7.0, 1.0e-13);
    EXPECT_NEAR(columns.variations[0].displacement(2), 2.0, 1.0e-13);
    EXPECT_NEAR(columns.variations[1].displacement(3), 2.0, 1.0e-13);
    EXPECT_NEAR(columns.variations[2].displacement(2), 5.0, 1.0e-13);
    EXPECT_NEAR(columns.variations[3].displacement(3), 5.0, 1.0e-13);
    EXPECT_EQ(columns.central_stages, 21u);
    EXPECT_EQ(columns.variation_stages, 21u);
    EXPECT_GT(columns.variation_metric_evaluations, 14u);
    EXPECT_LE(comparison.error_ratio, 1.0);
}

TEST(CoupledTransport, CameraColumnUnitsPreserveTheJointErrorDecision) {
    const auto ray = FlatRay();
    auto control = Columns();
    const auto original = control.variations;
    auto perturbed = original;
    perturbed[0].displacement(2) += 2e-8;
    perturbed[3].derivative(3) -= 4e-9;
    const double reference =
        Geodesic::CoupledStateError(ray, original, ray, perturbed, Control(), control);
    ASSERT_GT(reference, 1);
    const std::array<double, 4> units{0.0007, 0.002, 0.2, 0.2};
    auto scaled_original = original;
    auto scaled_perturbed = perturbed;
    for (int column = 0; column < 4; ++column) {
        scaled_original[column].displacement *= units[column];
        scaled_original[column].derivative *= units[column];
        scaled_perturbed[column].displacement *= units[column];
        scaled_perturbed[column].derivative *= units[column];
    }
    control.column_scale = units;
    const double scaled = Geodesic::CoupledStateError(ray, scaled_original, ray, scaled_perturbed,
                                                      Control(), control);
    EXPECT_NEAR(scaled, reference, 1e-12);
    control.column_scale = {1, 1, 1, 1};
    EXPECT_LT(Geodesic::CoupledStateError(ray, scaled_original, ray, scaled_perturbed, Control(),
                                          control),
              1);
}

TEST(CoupledTransport, PhysicalCameraColumnsReachTheLiveCpuSourceMap) {
    KerrSchildFamily metric(KerrSchildParams::Minkowski());
    CameraConfig config;
    config.width = 191;
    config.height = 107;
    config.beta_x = 0.2;
    config.beta_y = -0.1;
    config.beta_z = 0.3;
    ThinLensCamera camera(config);
    const auto projection = camera.ProjectFilmForObserver(75.25, 62.75, 0.71f, 0.23f);
    ASSERT_TRUE(projection);
    ASSERT_TRUE(projection->ray.phase_space);
    TracerConfig controls;
    controls.enable_disk = false;
    controls.escape_radius = 100;
    controls.max_steps = 1000;
    controls.integrator = Control();
    GeodesicTracer tracer(&metric, controls);
    const auto physical = tracer.Trace(projection->ray);
    ASSERT_FALSE(physical.numerical_failure);
    ASSERT_EQ(physical.outcome, TraceResult::Outcome::Escaped);
    ASSERT_TRUE(physical.beam.finite_source_map);
    auto angular_ray = projection->ray;
    angular_ray.phase_space.reset();
    const auto angular = tracer.Trace(angular_ray);
    ASSERT_FALSE(angular.numerical_failure);
    ASSERT_TRUE(angular.beam.finite_source_map);
    for (int row = 0; row < 2; ++row)
        for (int column = 0; column < 2; ++column)
            EXPECT_NEAR(physical.beam.finite_source_map->jacobian[row][column],
                        angular.beam.finite_source_map->jacobian[row][column], 1e-11);
    for (int axis = 0; axis < 3; ++axis)
        EXPECT_NEAR(physical.beam.finite_source_map->direction[axis],
                    angular.beam.finite_source_map->direction[axis], 1e-12);
}

TEST(CoupledTransport, NonfiniteColumnsAndUnrepresentedInteriorDeclineWithoutCommit) {
    KerrSchildFamily metric(KerrSchildParams::Minkowski());
    for (const double invalid :
         {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity()}) {
        auto ray = FlatRay();
        auto columns = Columns();
        columns.variations[3].derivative(2) = invalid;
        Rk45CoupledComparison comparison;
        EXPECT_FALSE(Geodesic::IntegrateStepRk45(ray, &metric, Control(), &columns, &comparison));
        EXPECT_EQ(ray.position(1), 5.0);
        EXPECT_NE(ray.terminated, 0);
    }
    auto ray = FlatRay();
    auto columns = Columns();
    auto config = Control();
    ray.step_size = std::numeric_limits<float>::denorm_min();
    config.min_step = ray.step_size;
    Rk45CoupledComparison comparison;
    EXPECT_FALSE(Geodesic::IntegrateStepRk45(ray, &metric, config, &columns, &comparison));
    EXPECT_EQ(ray.position(1), 5.0);
    EXPECT_EQ(columns.variations[0].displacement(2), 0.0);
    EXPECT_NE(ray.terminated, 0);
    EXPECT_EQ(columns.failure, CoupledStepFailure::Interpolation);
}

TEST(CoupledTransport, EndpointAndZeroFractionIncludeEventTimeVariation) {
    KerrSchildFamily metric(KerrSchildParams::Minkowski());
    auto start = FlatRay();
    auto end = start;
    end.position(0) -= 2.0;
    end.position(1) += 2.0;
    auto initial = Columns().variations;
    initial[2].displacement(1) = 3.0;
    auto final = initial;
    for (auto& column : final) column.displacement += column.derivative * 2.0;
    Vec4 normal;
    normal(1) = 1.0;
    for (const double fraction : {0.0, 0.5, 1.0}) {
        const auto event = Geodesic::SampleCoupledSegment(&metric, start, initial, end, final, 2.0,
                                                          fraction, &normal);
        ASSERT_TRUE(event);
        EXPECT_NEAR(event->variations[2].displacement(1), 0.0, 1.0e-13);
        EXPECT_NEAR(event->variations[2].displacement(0), 3.0, 1.0e-13);
        EXPECT_NEAR(event->variations[0].displacement(2), 2.0 * fraction, 1.0e-13);
    }
    normal(1) = 0.0;
    normal(2) = 1.0;
    EXPECT_FALSE(
        Geodesic::SampleCoupledSegment(&metric, start, initial, end, final, 2.0, 1.0, &normal));
}

TEST(CoupledTransport, EveryCanonicalColumnCanRejectProjectedAndInteriorCurvedError) {
    KerrSchildFamily metric(KerrSchildParams::Schwarzschild(1.0));
    auto initial = FlatRay();
    initial.velocity(1) = 0.8;
    initial.velocity(2) = 0.6;
    Metric4d values, inverse;
    Tensor<Dual<double>, 4, 4, 4> derivatives;
    metric.Evaluate(initial.position, values, derivatives);
    const auto projected = Geodesic::ProjectNullTangentPreservingBranch(initial.velocity, values);
    ASSERT_TRUE(projected);
    initial.velocity = *projected;
    ASSERT_TRUE(metric.InverseMetric(initial.position, inverse));
    std::array<Vec4, 3> seeds;
    for (int axis = 0; axis < 3; ++axis) seeds[axis](axis + 1) = 1.0;
    const auto observer = relativity::EulerianObserverFrame(values, inverse, seeds);
    ASSERT_TRUE(observer);
    const double frequency = TensorOps::InnerProduct(initial.velocity, observer->time, values);
    std::array<double, 3> direction{};
    for (int axis = 0; axis < 3; ++axis)
        direction[axis] =
            TensorOps::InnerProduct(initial.velocity, observer->spatial[axis], values) / frequency;
    const auto screen = relativity::ObserverScreenBasis(*observer, direction);
    ASSERT_TRUE(screen);
    for (int active = 0; active < 4; ++active) {
        for (const float interval : {0.125f, 0.5f}) {
            SCOPED_TRACE(active);
            SCOPED_TRACE(interval);
            auto ray = initial;
            ray.step_size = interval;
            auto unmonitored = ray;
            ASSERT_TRUE(Geodesic::IntegrateStepRk45(unmonitored, &metric, Control()));
            auto columns = Columns();
            columns.length_scale = 1.0;
            columns.variations = {};
            if (active < 2)
                columns.variations[active].derivative = (*screen)[active];
            else
                columns.variations[active].displacement = (*screen)[active - 2];
            const auto previous = columns.variations;
            Rk45CoupledComparison comparison;
            EXPECT_FALSE(
                Geodesic::IntegrateStepRk45(ray, &metric, Control(), &columns, &comparison));
            EXPECT_EQ(columns.failure, interval == 0.5f ? CoupledStepFailure::Projection
                                                        : CoupledStepFailure::Interpolation);
            EXPECT_EQ(ray.terminated, 0);
            EXPECT_LT(ray.step_size, interval);
            for (int component = 0; component < 4; ++component) {
                EXPECT_EQ(ray.position(component), initial.position(component));
                EXPECT_EQ(ray.velocity(component), initial.velocity(component));
                for (int column = 0; column < 4; ++column) {
                    EXPECT_EQ(columns.variations[column].displacement(component),
                              previous[column].displacement(component));
                    EXPECT_EQ(columns.variations[column].derivative(component),
                              previous[column].derivative(component));
                }
            }
            ray.step_size = 0.03125f;
            EXPECT_TRUE(
                Geodesic::IntegrateStepRk45(ray, &metric, Control(), &columns, &comparison));
            EXPECT_GT(ray.position(1), initial.position(1));
        }
    }
}

TEST(CoupledTransport, AcceptedIncrementsKeepRefinedEventDerivativesTranslationInvariant) {
    KerrSchildFamily metric(KerrSchildParams::Minkowski());
    for (const double origin : {5.0, 1.0e6}) {
        const double interval = 1.0e-5;
        auto start = FlatRay();
        start.position(1) = origin;
        auto end = start;
        CoupledSegmentIncrement increment;
        increment.position = start.velocity * interval;
        end.position += increment.position;
        auto initial = Columns().variations;
        initial[2].displacement(1) = 3.0;
        auto final = initial;
        for (std::size_t column = 0; column < initial.size(); ++column) {
            increment.displacement[column] = initial[column].derivative * interval;
            final[column].displacement += increment.displacement[column];
        }
        Vec4 normal;
        normal(1) = 1.0;
        const double radius = origin + 0.5 * interval;
        const auto root = FindSphericalBoundaryEvent(
            start.position, start.velocity, end.position, end.velocity, interval, radius,
            SphericalBoundarySense::IncreasingRadius, &increment.position);
        ASSERT_TRUE(root);
        const auto rooted =
            Geodesic::SampleCoupledSegment(&metric, start, initial, end, final, interval,
                                           root->fraction, &normal, nullptr, &increment);
        ASSERT_TRUE(rooted);
        // Isolation and sampling share retained polynomial coefficients;
        // the sphere rounding is bounded separately at represented-coordinate
        // precision. This is not a global large-offset accuracy claim.
        EXPECT_LE(std::abs(rooted->ray.position(1) - root->position(1)),
                  8.0 * std::numeric_limits<double>::epsilon() * origin);
        EXPECT_NEAR(rooted->ray.velocity(1), 1.0, 1.0e-12);
        for (int component = 0; component < 4; ++component)
            EXPECT_EQ(rooted->ray.velocity(component), root->tangent(component));
        for (const double fraction : {0.0, 0.5, 1.0}) {
            const auto event =
                Geodesic::SampleCoupledSegment(&metric, start, initial, end, final, interval,
                                               fraction, &normal, nullptr, &increment);
            ASSERT_TRUE(event);
            for (int component = 0; component < 4; ++component) {
                EXPECT_NEAR(event->ray.acceleration(component), 0.0, 1.0e-12);
                EXPECT_NEAR(event->variations[2].derivative(component), 0.0, 1.0e-12);
            }
            EXPECT_NEAR(event->variations[2].displacement(1), 0.0, 1.0e-12);
            EXPECT_NEAR(event->variations[2].displacement(0), 3.0, 1.0e-12);
            EXPECT_NEAR(event->variations[0].displacement(2), interval * fraction, 1.0e-12);
        }
    }
}

TEST(CoupledTransport, SkyMapIsIndependentOfPhysicalBundleOutputAndWorkLimitDeclines) {
    KerrSchildFamily metric(KerrSchildParams::Minkowski());
    CameraRay ray;
    ray.origin(1) = 5.0;
    ray.origin(2) = std::numbers::pi / 2.0;
    ray.direction(1) = -1.0;
    TracerConfig config;
    config.escape_radius = 12.0f;
    config.enable_disk = false;
    config.integrator.max_step = 7.0f;
    config.integrator.initial_step = 7.0f;
    config.max_steps = 100;
    GeodesicTracer plain(&metric, config);
    const auto reference = plain.Trace(ray);
    ASSERT_FALSE(reference.numerical_failure);
    ASSERT_TRUE(reference.beam.infinity_source_map);
    EXPECT_FALSE(reference.beam.valid);
    for (const bool point : {false, true}) {
        config.enable_ray_bundles = true;
        config.bundle_point_source = point;
        GeodesicTracer with_bundle(&metric, config);
        const auto result = with_bundle.Trace(ray);
        ASSERT_FALSE(result.numerical_failure);
        ASSERT_TRUE(result.beam.valid);
        ASSERT_TRUE(result.beam.infinity_source_map);
        EXPECT_EQ(result.steps_taken, reference.steps_taken);
        EXPECT_EQ(result.central_stages, reference.central_stages);
        EXPECT_EQ(result.affine_length, reference.affine_length);
        for (int component = 0; component < 4; ++component)
            EXPECT_EQ(result.final_position(component), reference.final_position(component));
        EXPECT_EQ(result.beam.infinity_source_map->map.jacobian,
                  reference.beam.infinity_source_map->map.jacobian);
    }
    config.max_steps = 1;
    GeodesicTracer incomplete(&metric, config);
    const auto declined = incomplete.Trace(ray);
    EXPECT_TRUE(declined.numerical_failure);
    EXPECT_EQ(declined.outcome, TraceResult::Outcome::MaxSteps);
    EXPECT_EQ(declined.coupled_failure, CoupledStepFailure::WorkLimit);
    EXPECT_FALSE(declined.beam.valid);
    EXPECT_FALSE(declined.beam.infinity_source_map);
}

// Capture's public beam is a physical screen ellipse at the central event's
// affine parameter. It is not a derivative of arrival time on the horizon.
TEST(CoupledTransport, CapturedKerrBeamMatchesNeighboursAndStepRefinement) {
    for (const double spin : {0.9, -0.9}) {
        KerrSchildParams parameters;
        parameters.M = 1.0;
        parameters.a = spin;
        KerrSchildFamily metric(parameters);
        OutgoingKerrSchild outgoing(metric);
        CameraConfig camera_config;
        camera_config.r = 50.0;
        camera_config.theta = std::numbers::pi / 2.0;
        camera_config.phi = 0.0;
        camera_config.fov = 60.0f;
        camera_config.width = camera_config.height = 64;
        const auto camera_ray =
            PinholeCamera(camera_config).GenerateRay(spin > 0.0 ? 38 : 25, 32, 0.5f, 0.5f);
        TracerConfig config;
        config.escape_radius = 100.0f;
        config.horizon_factor = 1.05f;
        config.max_steps = 20000;
        config.enable_disk = false;
        config.enable_ray_bundles = true;
        config.bundle_point_source = true;
        config.bundle_angular_size = 1.0e-3f;
        config.integrator.initial_step = 0.05f;
        config.integrator.max_step = 1.0f;
        config.integrator.abs_tolerance = config.integrator.rel_tolerance = 1.0e-7f;
        const auto trace = [&](const CameraRay& launch, const TracerConfig& settings) {
            GeodesicTracer tracer(&metric, settings);
            return tracer.Trace(launch);
        };
        const auto central = trace(camera_ray, config);
        ASSERT_FALSE(central.numerical_failure) << spin;
        ASSERT_EQ(central.outcome, TraceResult::Outcome::Horizon);
        ASSERT_EQ(central.terminal_chart, TraceResult::TerminalChart::OutgoingKerrSchild);
        ASSERT_TRUE(central.final_tangent);
        ASSERT_TRUE(central.beam.valid);
        EXPECT_FALSE(central.beam.finite_source_map);
        EXPECT_FALSE(central.beam.infinity_source_map);
        const double capture_radius = metric.OuterHorizonRadius() * config.horizon_factor;
        EXPECT_NEAR(metric.ComputeKerrRadius(central.final_position(1), central.final_position(2),
                                             central.final_position(3)),
                    capture_radius, 2.0e-12);
        const auto public_event = outgoing.ToIngoing(central.final_position);
        ASSERT_TRUE(public_event);
        const Vec4 tangent = public_event->Apply(*central.final_tangent);
        Metric4d values, inverse;
        Tensor<Dual<double>, 4, 4, 4> derivatives;
        metric.Evaluate(public_event->position, values, derivatives);
        ASSERT_TRUE(metric.InverseMetric(public_event->position, inverse));
        std::array<Vec4, 3> seeds{};
        for (int axis = 0; axis < 3; ++axis) seeds[axis](axis + 1) = 1.0;
        const auto observer = relativity::EulerianObserverFrame(values, inverse, seeds);
        ASSERT_TRUE(observer);
        const double frequency = TensorOps::InnerProduct(tangent, observer->time, values);
        ASSERT_GT(frequency, 0.0);
        std::array<double, 3> direction{};
        for (int axis = 0; axis < 3; ++axis)
            direction[axis] =
                TensorOps::InnerProduct(tangent, observer->spatial[axis], values) / frequency;
        const auto screen = relativity::ObserverScreenBasis(*observer, direction);
        ASSERT_TRUE(screen);
        EXPECT_LT(
            std::abs(TensorOps::InnerProduct(tangent, tangent, values)) / (frequency * frequency),
            1.0e-6);
        for (const auto& axis : *screen)
            EXPECT_LT(std::abs(TensorOps::InnerProduct(axis, tangent, values)) / frequency,
                      1.0e-10);
        const auto covariance = [&](const TraceResult::Beam& beam) {
            const double major = beam.semi_major / static_cast<double>(config.bundle_angular_size);
            const double minor = beam.semi_minor / static_cast<double>(config.bundle_angular_size);
            const double c = std::cos(beam.orientation), s = std::sin(beam.orientation);
            return std::array<double, 3>{major * major * c * c + minor * minor * s * s,
                                         (major * major - minor * minor) * c * s,
                                         major * major * s * s + minor * minor * c * c};
        };
        const auto expected = covariance(central.beam);
        const double scale =
            std::max({1.0, std::abs(expected[0]), std::abs(expected[1]), std::abs(expected[2])});
        // Independent orthonormal input basis; a rotation of that basis cannot
        // change the output covariance, so the production basis is not copied.
        Vec4 direction_vector = camera_ray.direction;
        double norm = std::hypot(direction_vector(1), direction_vector(2), direction_vector(3));
        direction_vector = direction_vector / norm;
        Vec4 first;
        first(2) = 1.0;
        double projection = direction_vector(2);
        first -= direction_vector * projection;
        first = first / std::hypot(first(1), first(2), first(3));
        Vec4 second;
        second(1) = direction_vector(2) * first(3) - direction_vector(3) * first(2);
        second(2) = direction_vector(3) * first(1) - direction_vector(1) * first(3);
        second(3) = direction_vector(1) * first(2) - direction_vector(2) * first(1);
        for (const double delta : {2.0e-5, 1.0e-5, 5.0e-6}) {
            double matrix[2][2]{};
            for (int column = 0; column < 2; ++column) {
                Vec4 endpoints[2];
                for (int side = 0; side < 2; ++side) {
                    auto neighbour = camera_ray;
                    neighbour.direction = direction_vector * std::cos(delta) +
                                          (column == 0 ? first : second) *
                                              ((side == 0 ? -1.0 : 1.0) * std::sin(delta));
                    auto settings = config;
                    settings.enable_ray_bundles = false;
                    settings.bundle_point_source = false;
                    settings.bundle_angular_size = TracerConfig{}.bundle_angular_size;
                    const auto result = trace(neighbour, settings);
                    ASSERT_FALSE(result.numerical_failure) << spin << ", " << delta;
                    ASSERT_EQ(result.outcome, TraceResult::Outcome::Horizon);
                    const auto mapped = outgoing.ToIngoing(result.final_position);
                    ASSERT_TRUE(mapped);
                    endpoints[side] = mapped->position;
                }
                const Vec4 deviation = (endpoints[1] - endpoints[0]) / (2.0 * delta);
                for (int row = 0; row < 2; ++row) {
                    matrix[row][column] =
                        TensorOps::InnerProduct((*screen)[row], deviation, values);
                    const double shifted = TensorOps::InnerProduct(
                        (*screen)[row], deviation + tangent * 250.0, values);
                    EXPECT_NEAR(shifted, matrix[row][column], 1.0e-8);
                }
                double component_scale = 1.0;
                for (int mu = 0; mu < 4; ++mu)
                    for (int nu = 0; nu < 4; ++nu)
                        component_scale +=
                            std::abs(values(mu, nu).real * tangent(mu) * deviation(nu));
                EXPECT_LT(
                    std::abs(TensorOps::InnerProduct(tangent, deviation, values)) / component_scale,
                    1.0e-4);
            }
            const std::array<double, 3> measured{
                matrix[0][0] * matrix[0][0] + matrix[0][1] * matrix[0][1],
                matrix[0][0] * matrix[1][0] + matrix[0][1] * matrix[1][1],
                matrix[1][0] * matrix[1][0] + matrix[1][1] * matrix[1][1]};
            for (int component = 0; component < 3; ++component)
                EXPECT_LT(std::abs(measured[component] - expected[component]) / scale, 1.0e-4)
                    << spin << ", " << delta << ", " << component;
        }
        for (const float maximum_step : {1.0f, 0.5f, 0.25f}) {
            auto settings = config;
            settings.integrator.max_step = maximum_step;
            const auto refined = trace(camera_ray, settings);
            ASSERT_FALSE(refined.numerical_failure);
            ASSERT_EQ(refined.outcome, TraceResult::Outcome::Horizon);
            ASSERT_TRUE(refined.beam.valid);
            const auto refined_covariance = covariance(refined.beam);
            for (int component = 0; component < 3; ++component)
                EXPECT_LT(std::abs(refined_covariance[component] - expected[component]) / scale,
                          1.0e-4);
            for (const bool point : {false, true}) {
                settings.enable_ray_bundles = true;
                settings.bundle_point_source = point;
                const auto with_bundle = trace(camera_ray, settings);
                auto plain_settings = settings;
                plain_settings.enable_ray_bundles = false;
                plain_settings.bundle_point_source = false;
                plain_settings.bundle_angular_size = TracerConfig{}.bundle_angular_size;
                const auto plain = trace(camera_ray, plain_settings);
                ASSERT_FALSE(with_bundle.numerical_failure);
                ASSERT_FALSE(plain.numerical_failure);
                ASSERT_TRUE(with_bundle.beam.valid);
                EXPECT_FALSE(plain.beam.valid);
                EXPECT_FALSE(with_bundle.beam.infinity_source_map);
                EXPECT_EQ(with_bundle.affine_length, plain.affine_length);
                EXPECT_EQ(with_bundle.steps_taken, plain.steps_taken);
                EXPECT_EQ(with_bundle.central_stages, plain.central_stages);
                ASSERT_TRUE(with_bundle.final_tangent);
                ASSERT_TRUE(plain.final_tangent);
                for (int component = 0; component < 4; ++component) {
                    EXPECT_EQ(with_bundle.final_position(component),
                              plain.final_position(component));
                    EXPECT_EQ((*with_bundle.final_tangent)(component),
                              (*plain.final_tangent)(component));
                }
            }
        }
    }
}

// A moving event samples the physical geodesic flow. A few ULPs in the
// endpoint tangent must not turn its cubic interpolant's second derivative
// into a fictitious covariant acceleration, amplified by arrival time.
TEST(CoupledTransport, PhysicalEventVariationRetainsArrivalPositionWithoutDenseAccelerationNoise) {
    KerrSchildFamily metric(KerrSchildParams::Minkowski());
    auto start = FlatRay();
    constexpr double interval = 1.0e-6;
    auto end = start;
    end.position = start.position + start.velocity * interval;
    end.velocity(1) = std::nextafter(start.velocity(1), 2.0);
    GeodesicVariations initial{}, final{};
    // Translating the starting event 200 units along x changes the arrival
    // time at a fixed x-plane by -200, without changing the photon tangent.
    initial[0].displacement(1) = 200.0;
    initial[1].derivative(2) = 1.0;
    final = initial;
    final[1].displacement(2) = interval;
    CoupledSegmentIncrement increment;
    increment.position = start.velocity * interval;
    increment.displacement[1](2) = interval;
    Vec4 normal;
    normal(1) = 1.0;
    for (const double fraction : {0.0, 0.25, 0.75, 1.0}) {
        SCOPED_TRACE(fraction);
        const auto event = Geodesic::SampleCoupledSegment(
            &metric, start, initial, end, final, interval, fraction, &normal, nullptr, &increment);
        ASSERT_TRUE(event);
        // Analytic neighboring flat lines hit x=start.x+fraction*interval.
        // Their time derivative is +200 and their x derivative is zero.
        EXPECT_NEAR(event->variations[0].displacement(0), 200.0, 1.0e-12);
        EXPECT_NEAR(event->variations[0].displacement(1), 0.0, 1.0e-12);
        EXPECT_NEAR(event->variations[1].displacement(2), fraction * interval, 1.0e-18);
        for (int component = 0; component < 4; ++component) {
            EXPECT_NEAR(event->variations[0].derivative(component), 0.0, 1.0e-12);
            EXPECT_NEAR(event->variations[1].derivative(component), component == 2 ? 1.0 : 0.0,
                        1.0e-12);
        }
    }
}

TEST(CoupledTransport, OriginalKerrDiskEventMatchesIndependentNeighboursAndRefinement) {
    KerrSchildFamily metric(KerrSchildParams::Kerr(1.0, 0.9));
    // Exact original 64x64, four-sample RenderSession ray (30,24), sample3.
    // Store the actual packet's ray: separately inlining the legacy float
    // pinhole arithmetic can differ by one ULP and miss this regression.
    CameraRay launch;
    launch.origin(1) = 50.0;
    launch.origin(2) = 1.5708;
    launch.direction(1) = -0.99130529165267944;
    launch.direction(2) = -0.12966860830783844;
    launch.direction(3) = -0.022356659173965454;
    TracerConfig config;
    config.escape_radius = 200.0f;
    config.horizon_factor = 1.0f;
    config.max_steps = 20000;
    config.enable_disk = true;
    config.disk_inner = AccretionDiskD::ComputeIsco(0.9);
    config.disk_outer = 20.0;
    config.integrator.initial_step = 0.1f;
    config.integrator.max_step = 2.0f;
    config.integrator.min_step = 1.0e-5f;
    config.integrator.abs_tolerance = config.integrator.rel_tolerance = 5.0e-6f;
    const auto trace = [&](const CameraRay& ray, const TracerConfig& settings) {
        GeodesicTracer tracer(&metric, settings);
        return tracer.Trace(ray);
    };
    const auto central = trace(launch, config);
    ASSERT_FALSE(central.numerical_failure);
    ASSERT_EQ(central.outcome, TraceResult::Outcome::DiskHit);
    ASSERT_TRUE(central.final_tangent);
    ASSERT_EQ(central.terminal_chart, TraceResult::TerminalChart::MetricNative);
    const auto event_roundoff = [](const Vec4& position) {
        // The polynomial locator returns a represented root, not a snapped
        // coordinate. Bound its double arithmetic at the local chart scale.
        return 256.0 * std::numeric_limits<double>::epsilon() *
               (1.0 + std::hypot(position(1), position(2), position(3)));
    };
    EXPECT_LE(std::abs(central.final_position(3)), event_roundoff(central.final_position));
    EXPECT_EQ(central.num_disk_crossings, 1);
    EXPECT_GT(central.disk_radius, config.disk_inner);
    EXPECT_LT(central.disk_radius, config.disk_outer);

    auto beam_config = config;
    beam_config.enable_ray_bundles = true;
    beam_config.bundle_point_source = true;
    beam_config.bundle_angular_size = 1.0e-3f;
    const auto beam = trace(launch, beam_config);
    ASSERT_FALSE(beam.numerical_failure);
    ASSERT_TRUE(beam.beam.valid);
    ASSERT_TRUE(beam.final_tangent);
    EXPECT_EQ(beam.steps_taken, central.steps_taken);
    EXPECT_EQ(beam.affine_length, central.affine_length);
    for (int component = 0; component < 4; ++component) {
        EXPECT_EQ(beam.final_position(component), central.final_position(component));
        EXPECT_EQ((*beam.final_tangent)(component), (*central.final_tangent)(component));
    }
    Metric4d values, inverse;
    Tensor<Dual<double>, 4, 4, 4> derivatives;
    metric.Evaluate(central.final_position, values, derivatives);
    ASSERT_TRUE(metric.InverseMetric(central.final_position, inverse));
    std::array<Vec4, 3> seeds{};
    for (int axis = 0; axis < 3; ++axis) seeds[axis](axis + 1) = 1.0;
    const auto observer = relativity::EulerianObserverFrame(values, inverse, seeds);
    ASSERT_TRUE(observer);
    const auto& tangent = *central.final_tangent;
    const double frequency = TensorOps::InnerProduct(tangent, observer->time, values);
    ASSERT_GT(frequency, 0.0);
    std::array<double, 3> direction{};
    for (int axis = 0; axis < 3; ++axis)
        direction[axis] =
            TensorOps::InnerProduct(tangent, observer->spatial[axis], values) / frequency;
    const auto screen = relativity::ObserverScreenBasis(*observer, direction);
    ASSERT_TRUE(screen);
    const double major =
        beam.beam.semi_major / static_cast<double>(beam_config.bundle_angular_size);
    const double minor =
        beam.beam.semi_minor / static_cast<double>(beam_config.bundle_angular_size);
    const double c = std::cos(beam.beam.orientation), s = std::sin(beam.beam.orientation);
    const std::array<double, 3> covariance{major * major * c * c + minor * minor * s * s,
                                           (major * major - minor * minor) * c * s,
                                           major * major * s * s + minor * minor * c * c};
    const double scale =
        std::max({1.0, std::abs(covariance[0]), std::abs(covariance[1]), std::abs(covariance[2])});
    // Choose an independent input basis. Its rotation does not affect the
    // output covariance; no production Jacobi values seed the neighboring rays.
    Vec4 n = launch.direction;
    n = n / std::hypot(n(1), n(2), n(3));
    Vec4 first;
    first(2) = 1.0;
    first -= n * n(2);
    first = first / std::hypot(first(1), first(2), first(3));
    Vec4 second;
    second(1) = n(2) * first(3) - n(3) * first(2);
    second(2) = n(3) * first(1) - n(1) * first(3);
    second(3) = n(1) * first(2) - n(2) * first(1);
    double previous_covariance_error = -1.0;
    for (const double delta : {1.0e-4, 5.0e-5}) {
        SCOPED_TRACE(delta);
        double matrix[2][2]{};
        for (int column = 0; column < 2; ++column) {
            Vec4 endpoints[2];
            for (int side = 0; side < 2; ++side) {
                auto neighbour = launch;
                neighbour.direction =
                    n * std::cos(delta) +
                    (column == 0 ? first : second) * ((side == 0 ? -1.0 : 1.0) * std::sin(delta));
                const auto result = trace(neighbour, config);
                ASSERT_FALSE(result.numerical_failure);
                ASSERT_EQ(result.outcome, TraceResult::Outcome::DiskHit);
                ASSERT_TRUE(result.final_tangent);
                EXPECT_LE(std::abs(result.final_position(3)),
                          event_roundoff(result.final_position));
                endpoints[side] = result.final_position;
            }
            const Vec4 derivative = (endpoints[1] - endpoints[0]) / (2.0 * delta);
            for (int row = 0; row < 2; ++row)
                matrix[row][column] = TensorOps::InnerProduct((*screen)[row], derivative, values);
        }
        const std::array<double, 3> measured{
            matrix[0][0] * matrix[0][0] + matrix[0][1] * matrix[0][1],
            matrix[0][0] * matrix[1][0] + matrix[0][1] * matrix[1][1],
            matrix[1][0] * matrix[1][0] + matrix[1][1] * matrix[1][1]};
        double covariance_error = 0.0;
        for (int component = 0; component < 3; ++component)
            covariance_error = std::max(
                covariance_error, std::abs(measured[component] - covariance[component]) / scale);
        // Same declared physical-beam accuracy as the captured-neighbor test.
        EXPECT_LT(covariance_error, 1.0e-4);
        std::ostringstream measured_error;
        measured_error << std::setprecision(17) << covariance_error;
        RecordProperty(
            delta == 1.0e-4 ? "coarse_covariance_relative_error" : "fine_covariance_relative_error",
            measured_error.str());
        // Axes/orientation are published as float. Require actual refinement
        // improvement when the previous error exceeds their roundoff floor.
        if (previous_covariance_error > 8.0 * std::numeric_limits<float>::epsilon()) {
            EXPECT_LT(covariance_error, previous_covariance_error);
        }
        previous_covariance_error = covariance_error;
    }
    for (const float cap : {0.1f, 0.025f}) {
        auto settings = config;
        settings.integrator.initial_step = settings.integrator.max_step = cap;
        const auto refined = trace(launch, settings);
        ASSERT_FALSE(refined.numerical_failure);
        ASSERT_EQ(refined.outcome, TraceResult::Outcome::DiskHit);
        ASSERT_TRUE(refined.final_tangent);
        for (int component = 0; component < 4; ++component) {
            EXPECT_NEAR(refined.final_position(component), central.final_position(component),
                        1.0e-10);
            EXPECT_NEAR((*refined.final_tangent)(component), tangent(component), 1.0e-10);
        }
    }
}

}  // namespace
