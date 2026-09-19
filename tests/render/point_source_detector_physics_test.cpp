#include "sirius/backend/cpu/geodesic_tracer.h"
#include "sirius/core/camera.h"
#include "sirius/core/metrics/kerr_schild_family.h"
#include "sirius/render/session/point_source_detector.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <format>
#include <numbers>
#include <vector>

namespace sirius::test {
namespace {
using namespace sirius::render;

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
