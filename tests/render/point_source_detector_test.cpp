#include "sirius/render/session/point_source_detector.h"

#include "sirius/core/observer_frame.h"
#include "sirius/core/spectral/point_source_transfer.h"

#include <gtest/gtest.h>

#include <cmath>
#include <numbers>
#include <set>
#include <vector>

namespace sirius::test {
namespace {
using namespace sirius::render;
constexpr double kScale = .001;

PointDetectorProbe Gnomonic(double x, double y, double dx, double dy) {
    const double norm = std::hypot(1.0, x, y);
    PointDetectorProbe result;
    result.visible = true;
    result.direction = {1 / norm, x / norm, y / norm};
    const auto basis = core::relativity::MakeCelestialTangentBasis(result.direction);
    for (int column = 0; column < 2; ++column) {
        const double derivative = column == 0 ? dx : dy;
        const double projection = result.direction[column + 1] * derivative;
        for (int axis = 0; axis < 3; ++axis) {
            const double component =
                ((axis == column + 1 ? derivative : 0) - result.direction[axis] * projection) /
                norm;
            result.source_derivative[0][column] += basis->first[axis] * component;
            result.source_derivative[1][column] += basis->second[axis] * component;
        }
    }
    return result;
}

core::StarEntry Star(double x, double y) {
    const auto p = Gnomonic(x, y, 1, 1);
    return {static_cast<float>(p.direction[0]),
            static_cast<float>(p.direction[1]),
            static_cast<float>(p.direction[2]),
            10,
            4,
            .5f,
            5800,
            0};
}

std::array<double, 3> Expected(const core::StarEntry& star, DetectorCoordinate z,
                               double determinant, double g = 1, double transmission = 1) {
    const double radius = std::hypot(z[0], z[1]);
    if (radius > 4) return {};
    const double density = std::exp(-.5 * radius * radius) /
                           (2 * std::numbers::pi * -std::expm1(-8.0) * std::abs(determinant));
    auto result =
        *core::spectral::TransferPointSourceBand(star.temperature_K, g, star.Intensity(), density);
    for (auto& value : result) value *= transmission;
    return result;
}

TEST(PointSourceDetector, SharedDiscoveryPreservesDistinctGaussianShapesAndTransfer) {
    core::StarfieldSpatialIndex catalogue(
        {Star(.1 * kScale, .2 * kScale), Star(-1.7 * kScale, -.4 * kScale),
         Star(.6 * kScale, 2.2 * kScale), Star(4.6 * kScale, 0), Star(-5 * kScale, 0)});
    std::vector<PointDetectorFootprint> footprints;
    for (int y = 0; y < 4; ++y)
        for (int x = 0; x < 4; ++x)
            footprints.push_back(
                {{x - 1.5, y - 1.5}, {{{.7 + .05 * x, .3}, {-.2, 1.1 - .03 * y}}}});
    const PointDetectorSampler sample = [](DetectorCoordinate q) {
        auto point = Gnomonic(kScale * q[0], kScale * q[1], kScale, kScale);
        point.camera_over_source_frequency = std::exp(.1 * q[0]);
        point.transmission = .5 + .02 * q[1];
        return point;
    };
    const auto result =
        EvaluatePointDetectorBatch(catalogue, 1, footprints, sample, [] { return false; });
    ASSERT_TRUE(result) << static_cast<int>(result.error().reason) << ' '
                        << result.error().statistics.probes;
    ASSERT_EQ(result->samples.size(), footprints.size());
    std::size_t separate_probes = 0;
    for (std::size_t i = 0; i < footprints.size(); ++i) {
        const auto& footprint = footprints[i];
        const auto& m = footprint.chart_from_standard;
        const double determinant = m[0][0] * m[1][1] - m[0][1] * m[1][0];
        std::array<double, 3> expected{};
        for (const auto& star : catalogue.Stars()) {
            const double x = double(star.direction_y) / star.direction_x;
            const double y = double(star.direction_z) / star.direction_x;
            const double u = x / kScale - footprint.centre[0];
            const double v = y / kScale - footprint.centre[1];
            const DetectorCoordinate original{(m[1][1] * u - m[0][1] * v) / determinant,
                                              (-m[1][0] * u + m[0][0] * v) / determinant};
            const auto contribution = Expected(
                star, original, determinant * kScale * kScale / std::pow(1 + x * x + y * y, 1.5),
                std::exp(.1 * x / kScale), .5 + .02 * y / kScale);
            for (int channel = 0; channel < 3; ++channel)
                expected[channel] += contribution[channel];
        }
        const auto separate = EvaluatePointDetector(
            catalogue, 1,
            [&](DetectorCoordinate z) -> std::expected<PointDetectorProbe, PointDetectorFailure> {
                auto point = sample({footprint.centre[0] + m[0][0] * z[0] + m[0][1] * z[1],
                                     footprint.centre[1] + m[1][0] * z[0] + m[1][1] * z[1]});
                if (!point) return point;
                const auto original = point->source_derivative;
                for (int row = 0; row < 2; ++row)
                    for (int column = 0; column < 2; ++column)
                        point->source_derivative[row][column] =
                            original[row][0] * m[0][column] + original[row][1] * m[1][column];
                return point;
            },
            [] { return false; });
        ASSERT_TRUE(separate);
        separate_probes += separate->statistics.probes;
        for (int channel = 0; channel < 3; ++channel) {
            EXPECT_NEAR(result->samples[i].rgb[channel], expected[channel],
                        2e-6 * expected[channel]);
            EXPECT_NEAR(result->samples[i].rgb[channel], separate->rgb[channel],
                        2e-6 * expected[channel]);
            EXPECT_LE(result->samples[i].estimated_error[channel],
                      1e-7 + 1e-3 * result->samples[i].rgb[channel]);
        }
    }
    EXPECT_LT(result->statistics.probes * 4, separate_probes);
    RecordProperty("shared_probes", static_cast<int>(result->statistics.probes));
    RecordProperty("separate_probes", static_cast<int>(separate_probes));
}

TEST(PointSourceDetector, SharedDiscoveryDeclinesMalformedExhaustedAndCancelledBatches) {
    core::StarfieldSpatialIndex catalogue({Star(0, 0)});
    std::array<PointDetectorFootprint, 2> footprints{};
    footprints[1].centre = {.5, -.25};
    const PointDetectorSampler sample = [](DetectorCoordinate z) {
        return Gnomonic(kScale * z[0], kScale * z[1], kScale, kScale);
    };
    auto policy = PointDetectorPolicy{};
    policy.maximum_probes = 1;
    const auto exhausted =
        EvaluatePointDetectorBatch(catalogue, 1, footprints, sample, [] { return false; }, policy);
    ASSERT_FALSE(exhausted);
    EXPECT_EQ(exhausted.error().reason, PointDetectorFailure::WorkLimit);
    const auto cancelled =
        EvaluatePointDetectorBatch(catalogue, 1, footprints, sample, [] { return true; });
    ASSERT_FALSE(cancelled);
    EXPECT_EQ(cancelled.error().reason, PointDetectorFailure::Cancelled);
    const auto failed = EvaluatePointDetectorBatch(
        catalogue, 1, footprints,
        [](DetectorCoordinate) -> std::expected<PointDetectorProbe, PointDetectorFailure> {
            return std::unexpected(PointDetectorFailure::TraceFailed);
        },
        [] { return false; });
    ASSERT_FALSE(failed);
    EXPECT_EQ(failed.error().reason, PointDetectorFailure::TraceFailed);
    footprints[1].chart_from_standard = {{{1, 1}, {1, 1}}};
    const auto singular =
        EvaluatePointDetectorBatch(catalogue, 1, footprints, sample, [] { return false; });
    ASSERT_FALSE(singular);
    EXPECT_EQ(singular.error().reason, PointDetectorFailure::InvalidInput);
    EXPECT_EQ(singular.error().statistics.probes, 0);
    EXPECT_FALSE(EvaluatePointDetectorBatch(catalogue, 1, {}, sample, [] { return false; }));
}

TEST(PointSourceDetector, SharedDiscoveryDeclinesUncertainOriginalSupportOwnership) {
    auto star = Star(0, 0);
    // This represented catalogue direction has the exact gnomonic ratio
    // (1/128)/1, hence an image at z=(4,0) when the map scale is 1/512.
    star.direction_y = 1.0f / 128;
    core::StarfieldSpatialIndex catalogue({star});
    std::array<PointDetectorFootprint, 2> footprints{};
    footprints[1].centre = {1, 0};
    const auto result = EvaluatePointDetectorBatch(
        catalogue, 1, footprints,
        [](DetectorCoordinate z) { return Gnomonic(z[0] / 512, z[1] / 512, 1.0 / 512, 1.0 / 512); },
        [] { return false; });
    ASSERT_FALSE(result);
    EXPECT_EQ(result.error().reason, PointDetectorFailure::Unresolved);
}

TEST(PointSourceDetector, OriginalGaussianOwnsSharedEdgesAndRejectsOutsideSupport) {
    core::StarfieldSpatialIndex catalogue({Star(0, 0), Star(.3 * kScale, .4 * kScale),
                                           Star(2 * kScale, -kScale), Star(4.2 * kScale, 0)});
    const auto result = EvaluatePointDetector(
        catalogue, 1,
        [](DetectorCoordinate z) { return Gnomonic(kScale * z[0], kScale * z[1], kScale, kScale); },
        [] { return false; });
    ASSERT_TRUE(result) << static_cast<int>(result.error().reason) << ' '
                        << result.error().statistics.probes;
    std::array<double, 3> expected{};
    for (const auto& star : catalogue.Stars()) {
        const double x = double(star.direction_y) / star.direction_x;
        const double y = double(star.direction_z) / star.direction_x;
        const auto value = Expected(star, {x / kScale, y / kScale},
                                    kScale * kScale / std::pow(1 + x * x + y * y, 1.5));
        for (int channel = 0; channel < 3; ++channel) expected[channel] += value[channel];
    }
    for (int channel = 0; channel < 3; ++channel)
        EXPECT_NEAR(result->rgb[channel], expected[channel], 2e-7 * expected[channel]);
}

TEST(PointSourceDetector, CacheRetainsClosePhysicalCoordinatesWithBoundedLookupWork) {
    const std::vector<core::StarEntry> stars{Star(0, 0), Star(1e-10 * kScale, -1e-10 * kScale)};
    core::StarfieldSpatialIndex catalogue(stars);
    std::set<DetectorCoordinate> sampled;
    const auto result = EvaluatePointDetector(
        catalogue, 1,
        [&](DetectorCoordinate z) {
            EXPECT_TRUE(sampled.insert(z).second) << "a shared ray was traced twice";
            return Gnomonic(kScale * z[0], kScale * z[1], kScale, kScale);
        },
        [] { return false; });
    ASSERT_TRUE(result);
    EXPECT_EQ(result->statistics.probes, sampled.size());
    EXPECT_GT(result->statistics.probe_requests, result->statistics.probes);
    EXPECT_GE(result->statistics.probe_cache_comparisons,
              result->statistics.probe_requests - result->statistics.probes);
    // Count work instead of timing this machine: a scan of every prior sample
    // grows quadratically even though the physical rays and result are equal.
    EXPECT_LT(result->statistics.probe_cache_comparisons, 8 * result->statistics.probe_requests);
    EXPECT_TRUE(sampled.contains({0, 0}));
    EXPECT_TRUE(std::any_of(sampled.begin(), sampled.end(), [](DetectorCoordinate z) {
        const double separation = std::hypot(z[0], z[1]);
        return separation > 0 && separation < 1e-8;
    })) << "a distinct corrected image ray was merged into the centre";
    std::array<double, 3> expected{};
    for (const auto& star : stars) {
        const double x = double(star.direction_y) / star.direction_x;
        const double y = double(star.direction_z) / star.direction_x;
        const auto image = Expected(star, {x / kScale, y / kScale},
                                    kScale * kScale / std::pow(1 + x * x + y * y, 1.5));
        for (int channel = 0; channel < 3; ++channel) expected[channel] += image[channel];
    }
    for (int channel = 0; channel < 3; ++channel)
        EXPECT_NEAR(result->rgb[channel], expected[channel], 2e-7 * expected[channel]);
}

TEST(PointSourceDetector, ImageFrequencyAndTransmissionAreAppliedOnce) {
    const auto star = Star(.7 * kScale, -.3 * kScale);
    core::StarfieldSpatialIndex catalogue({star});
    const auto result = EvaluatePointDetector(
        catalogue, 1,
        [](DetectorCoordinate z) {
            auto p = Gnomonic(kScale * z[0], kScale * z[1], kScale, kScale);
            p.camera_over_source_frequency = std::exp(.2 * z[0]);
            p.transmission = std::exp(-2 - .1 * z[1]);
            p.visible = z[0] > .2;  // Original centre is captured.
            return p;
        },
        [] { return false; });
    ASSERT_TRUE(result) << static_cast<int>(result.error().reason) << ' '
                        << result.error().statistics.probes;
    const double x = double(star.direction_y) / star.direction_x;
    const double y = double(star.direction_z) / star.direction_x;
    const auto expected =
        Expected(star, {x / kScale, y / kScale}, kScale * kScale / std::pow(1 + x * x + y * y, 1.5),
                 std::exp(.2 * x / kScale), std::exp(-2 - .1 * y / kScale));
    for (int channel = 0; channel < 3; ++channel)
        EXPECT_NEAR(result->rgb[channel], expected[channel], 2e-7 * expected[channel]);
}

TEST(PointSourceDetector, FoldKeepsBothImagesOfTheSameCatalogueEntry) {
    const auto star = Star(.2 * kScale, .3 * kScale);
    core::StarfieldSpatialIndex catalogue({star});
    auto policy = PointDetectorPolicy{};
    const auto result = EvaluatePointDetector(
        catalogue, 1,
        [](DetectorCoordinate z) {
            return Gnomonic(kScale * (z[0] * z[0] - .4), kScale * z[1], 2 * kScale * z[0], kScale);
        },
        [] { return false; }, policy);
    ASSERT_TRUE(result) << static_cast<int>(result.error().reason) << ' '
                        << result.error().statistics.probes;
    const double x = double(star.direction_y) / star.direction_x;
    const double y = double(star.direction_z) / star.direction_x;
    const double root = std::sqrt(x / kScale + .4);
    const auto image = Expected(star, {root, y / kScale},
                                2 * kScale * kScale * root / std::pow(1 + x * x + y * y, 1.5));
    for (int channel = 0; channel < 3; ++channel)
        EXPECT_NEAR(result->rgb[channel], 2 * image[channel], 2e-6 * image[channel]);
}

TEST(PointSourceDetector, FoldOwnershipSurvivesRotationTranslationAndTighterRefinement) {
    const auto star = Star(.2 * kScale, .3 * kScale);
    core::StarfieldSpatialIndex catalogue({star});
    for (const double angle : {.23, .71}) {
        SCOPED_TRACE(angle);
        const double c = std::cos(angle), s = std::sin(angle);
        const DetectorCoordinate shift{.137, -.219};
        const PointDetectorSampler sample = [&](DetectorCoordinate z) {
            const double u = c * (z[0] - shift[0]) + s * (z[1] - shift[1]);
            const double v = -s * (z[0] - shift[0]) + c * (z[1] - shift[1]);
            auto point = Gnomonic(kScale * (u * u - .4), kScale * v, 2 * kScale * u, kScale);
            const auto b = point.source_derivative;
            for (int row = 0; row < 2; ++row) {
                point.source_derivative[row][0] = b[row][0] * c - b[row][1] * s;
                point.source_derivative[row][1] = b[row][0] * s + b[row][1] * c;
            }
            return point;
        };
        const double x = double(star.direction_y) / star.direction_x;
        const double y = double(star.direction_z) / star.direction_x;
        const double root = std::sqrt(x / kScale + .4);
        std::array<double, 3> expected{};
        for (double u : {-root, root}) {
            const DetectorCoordinate z{shift[0] + c * u - s * y / kScale,
                                       shift[1] + s * u + c * y / kScale};
            const auto image =
                Expected(star, z, 2 * kScale * kScale * u / std::pow(1 + x * x + y * y, 1.5));
            for (int channel = 0; channel < 3; ++channel) expected[channel] += image[channel];
        }
        for (bool tighter : {false, true}) {
            SCOPED_TRACE(tighter);
            PointDetectorPolicy policy;
            if (tighter) {
                policy.geometry_error *= .5;
                policy.maximum_linearization_residual *= .5;
                policy.relative_rgb_error *= .5;
                policy.root_error *= .5;
                policy.maximum_probes *= 2;
                policy.maximum_cells *= 2;
                policy.maximum_depth += 1;
            }
            const auto result =
                EvaluatePointDetector(catalogue, 1, sample, [] { return false; }, policy);
            ASSERT_TRUE(result) << static_cast<int>(result.error().reason) << ' '
                                << result.error().statistics.probes;
            for (int channel = 0; channel < 3; ++channel)
                EXPECT_NEAR(result->rgb[channel], expected[channel], 2e-6 * expected[channel]);
        }
    }
}

TEST(PointSourceDetector, ExhaustionCancellationAndFailedTracesHaveNoPartialRadiance) {
    core::StarfieldSpatialIndex catalogue({Star(0, 0)});
    auto policy = PointDetectorPolicy{};
    policy.maximum_probes = 1;
    const PointDetectorSampler sample = [](DetectorCoordinate z) {
        return Gnomonic(kScale * z[0], kScale * z[1], kScale, kScale);
    };
    const auto exhausted =
        EvaluatePointDetector(catalogue, 1, sample, [] { return false; }, policy);
    ASSERT_FALSE(exhausted);
    EXPECT_EQ(exhausted.error().reason, PointDetectorFailure::WorkLimit);
    const auto cancelled = EvaluatePointDetector(catalogue, 1, sample, [] { return true; });
    ASSERT_FALSE(cancelled);
    EXPECT_EQ(cancelled.error().reason, PointDetectorFailure::Cancelled);
    const auto failed = EvaluatePointDetector(
        catalogue, 1,
        [](DetectorCoordinate) -> std::expected<PointDetectorProbe, PointDetectorFailure> {
            return std::unexpected(PointDetectorFailure::TraceFailed);
        },
        [] { return false; });
    ASSERT_FALSE(failed);
    EXPECT_EQ(failed.error().reason, PointDetectorFailure::TraceFailed);
}

TEST(PointSourceDetector, PolynomialFoldResolvesCloseImagesBetweenTheCoarseNodes) {
    const auto star = Star(0, .2 * kScale);
    core::StarfieldSpatialIndex catalogue({star});
    // The close positive pair is hidden between the parent's regular nodes.
    // The independent oracle uses the three exact polynomial roots.
    for (double separation : {.3, .08}) {
        SCOPED_TRACE(separation);
        const double a = -.83, b = .413, c = b + separation;
        const double angle = .43, cosine = std::cos(angle), sine = std::sin(angle);
        const PointDetectorSampler sample = [&](DetectorCoordinate z) {
            const double u = cosine * z[0] + sine * z[1];
            const double v = -sine * z[0] + cosine * z[1];
            const double f = (u - a) * (u - b) * (u - c);
            const double df = (u - b) * (u - c) + (u - a) * (u - c) + (u - a) * (u - b);
            auto point = Gnomonic(kScale * f, kScale * v, kScale * df, kScale);
            const auto derivative = point.source_derivative;
            for (int row = 0; row < 2; ++row) {
                point.source_derivative[row][0] =
                    derivative[row][0] * cosine - derivative[row][1] * sine;
                point.source_derivative[row][1] =
                    derivative[row][0] * sine + derivative[row][1] * cosine;
            }
            return point;
        };
        const auto result = EvaluatePointDetector(catalogue, 1, sample, [] { return false; });
        ASSERT_TRUE(result) << static_cast<int>(result.error().reason) << ' '
                            << result.error().statistics.probes;
        const double y = double(star.direction_z) / star.direction_x;
        std::array<double, 3> expected{};
        for (const double u : {a, b, c}) {
            const double df = (u - b) * (u - c) + (u - a) * (u - c) + (u - a) * (u - b);
            const DetectorCoordinate z{cosine * u - sine * y / kScale,
                                       sine * u + cosine * y / kScale};
            const auto image = Expected(star, z, kScale * kScale * df / std::pow(1 + y * y, 1.5));
            for (int channel = 0; channel < 3; ++channel) expected[channel] += image[channel];
        }
        for (int channel = 0; channel < 3; ++channel) {
            EXPECT_NEAR(result->rgb[channel], expected[channel], 2e-6 * expected[channel]);
            EXPECT_LE(result->estimated_error[channel], 1e-7 + 1e-3 * result->rgb[channel]);
        }
        std::array<PointDetectorFootprint, 3> footprints{};
        footprints[1].centre = {-.3, .2};
        footprints[2].centre = {.7, -.4};
        const auto batch =
            EvaluatePointDetectorBatch(catalogue, 1, footprints, sample, [] { return false; });
        ASSERT_TRUE(batch) << static_cast<int>(batch.error().reason) << ' '
                           << batch.error().statistics.probes;
        for (std::size_t packet = 0; packet < footprints.size(); ++packet) {
            expected = {};
            for (const double u : {a, b, c}) {
                const double df = (u - b) * (u - c) + (u - a) * (u - c) + (u - a) * (u - b);
                const DetectorCoordinate z{
                    cosine * u - sine * y / kScale - footprints[packet].centre[0],
                    sine * u + cosine * y / kScale - footprints[packet].centre[1]};
                const auto image =
                    Expected(star, z, kScale * kScale * df / std::pow(1 + y * y, 1.5));
                for (int channel = 0; channel < 3; ++channel) expected[channel] += image[channel];
            }
            for (int channel = 0; channel < 3; ++channel)
                EXPECT_NEAR(batch->samples[packet].rgb[channel], expected[channel],
                            2e-6 * expected[channel]);
        }
    }
}

TEST(PointSourceDetector, DisconnectedVisibilityIslandAndStripeUseActualImageVisibility) {
    const std::vector<core::StarEntry> stars{
        Star(.273 * kScale, -.331 * kScale), Star(-.73 * kScale, .41 * kScale),
        Star(.273 * kScale, .31 * kScale), Star(.83 * kScale, -.91 * kScale)};
    core::StarfieldSpatialIndex catalogue(stars);
    const auto visible = [](DetectorCoordinate z) {
        return std::hypot(z[0] - .273, z[1] + .331) < .22 ||
               std::abs(.8 * z[0] + .6 * z[1] + .338) < .09;
    };
    const PointDetectorSampler sample = [&](DetectorCoordinate z) {
        auto point = Gnomonic(kScale * z[0], kScale * z[1], kScale, kScale);
        point.visible = visible(z);
        point.inner_attempts = 7;
        point.tail_attempts = 3;
        return point;
    };
    const auto result = EvaluatePointDetector(catalogue, 1, sample, [] { return false; });
    ASSERT_TRUE(result) << static_cast<int>(result.error().reason) << ' '
                        << result.error().statistics.probes;
    std::array<double, 3> expected{};
    for (const auto& star : stars) {
        const double x = double(star.direction_y) / star.direction_x;
        const double y = double(star.direction_z) / star.direction_x;
        const DetectorCoordinate z{x / kScale, y / kScale};
        if (!visible(z)) continue;
        const auto image = Expected(star, z, kScale * kScale / std::pow(1 + x * x + y * y, 1.5));
        for (int channel = 0; channel < 3; ++channel) expected[channel] += image[channel];
    }
    for (int channel = 0; channel < 3; ++channel)
        EXPECT_NEAR(result->rgb[channel], expected[channel], 2e-7 * expected[channel]);
    EXPECT_EQ(result->statistics.inner_attempts, 7 * result->statistics.probes);
    EXPECT_EQ(result->statistics.tail_attempts, 3 * result->statistics.probes);
    std::array<PointDetectorFootprint, 2> footprints{};
    footprints[1].centre = {.2, .1};
    const auto batch =
        EvaluatePointDetectorBatch(catalogue, 1, footprints, sample, [] { return false; });
    ASSERT_TRUE(batch) << static_cast<int>(batch.error().reason) << ' '
                       << batch.error().statistics.probes;
    for (std::size_t packet = 0; packet < footprints.size(); ++packet) {
        expected = {};
        for (const auto& star : stars) {
            const double x = double(star.direction_y) / star.direction_x;
            const double y = double(star.direction_z) / star.direction_x;
            const DetectorCoordinate q{x / kScale, y / kScale};
            if (!visible(q)) continue;
            const auto image = Expected(
                star, {q[0] - footprints[packet].centre[0], q[1] - footprints[packet].centre[1]},
                kScale * kScale / std::pow(1 + x * x + y * y, 1.5));
            for (int channel = 0; channel < 3; ++channel) expected[channel] += image[channel];
        }
        for (int channel = 0; channel < 3; ++channel)
            EXPECT_NEAR(batch->samples[packet].rgb[channel], expected[channel],
                        2e-7 * expected[channel]);
    }
    EXPECT_EQ(batch->statistics.inner_attempts, 7 * batch->statistics.probes);
    EXPECT_EQ(batch->statistics.tail_attempts, 3 * batch->statistics.probes);
}
}  // namespace
}  // namespace sirius::test
