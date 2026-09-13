#include "sirius/core/metrics/kerr_schild_family.h"
#include "sirius/core/observer_frame.h"

#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <limits>
#include <optional>

namespace sirius::test {
namespace {

using namespace sirius::core;
using namespace sirius::core::relativity;
using Direction = std::array<double, 3>;

Vec4 Vector(double t, double x, double y, double z) {
    Vec4 result;
    result(0) = t;
    result(1) = x;
    result(2) = y;
    result(3) = z;
    return result;
}

double Dot(const Direction& a, const Direction& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

Direction Unit(Direction direction) {
    const double norm = std::sqrt(Dot(direction, direction));
    for (double& component : direction) component /= norm;
    return direction;
}

std::array<Vec4, 3> CartesianSeeds() {
    return {Vector(0, 1, 0, 0), Vector(0, 0, 1, 0), Vector(0, 0, 0, 1)};
}

struct MetricSample {
    Metric4d metric;
    Metric4d inverse;
    Metric4d variation;
};

MetricSample Evaluate(IMetric& authority, const Vec4& position, const Vec4& displacement) {
    MetricSample sample;
    Tensor<Dual<double>, 4, 4, 4> derivative;
    authority.Evaluate(position, sample.metric, derivative);
    if (!authority.InverseMetric(position, sample.inverse)) {
        sample.inverse = TensorOps::Inverse(sample.metric);
    }
    for (int mu = 0; mu < 4; ++mu) {
        for (int nu = 0; nu < 4; ++nu) {
            double value = 0.0;
            for (int axis = 0; axis < 4; ++axis) {
                value += derivative(axis, mu, nu).real * displacement(axis);
            }
            sample.variation(mu, nu) = Dual<double>(value);
        }
    }
    return sample;
}

// Nominal evaluations only: this oracle never calls the differentiated frame
// or the new source-direction/JVP implementation.
std::optional<std::array<double, 4>> NominalSource(IMetric& authority, const Vec4& position,
                                                   const Vec4& tangent) {
    const auto sample = Evaluate(authority, position, Vec4{});
    const auto frame = EulerianObserverFrame(sample.metric, sample.inverse, CartesianSeeds());
    if (!frame) return std::nullopt;
    const double frequency = TensorOps::InnerProduct(tangent, frame->time, sample.metric);
    if (!(frequency > 0.0)) return std::nullopt;
    Direction local{};
    for (int axis = 0; axis < 3; ++axis) {
        local[axis] =
            TensorOps::InnerProduct(tangent, frame->spatial[axis], sample.metric) / frequency;
    }
    local = Unit(local);
    return std::array{local[0], local[1], local[2], frequency};
}

// A four-evaluation, fourth-order finite difference of the actual nominal
// source function. Both x and k vary; the neighboring k need not remain exactly
// null, because shading explicitly normalizes its locally measured direction.
std::optional<std::array<double, 4>> DifferenceOracle(IMetric& authority, const Vec4& position,
                                                      const Vec4& tangent, const Vec4& displacement,
                                                      const Vec4& tangent_variation, double h) {
    std::array<std::array<double, 4>, 4> values{};
    const std::array<double, 4> offsets{-2.0, -1.0, 1.0, 2.0};
    for (std::size_t index = 0; index < offsets.size(); ++index) {
        const double offset = offsets[index] * h;
        const auto value = NominalSource(authority, position + displacement * offset,
                                         tangent + tangent_variation * offset);
        if (!value) return std::nullopt;
        values[index] = *value;
    }
    std::array<double, 4> derivative{};
    for (int component = 0; component < 4; ++component) {
        derivative[component] = (values[0][component] - 8.0 * values[1][component] +
                                 8.0 * values[2][component] - values[3][component]) /
                                (12.0 * h);
    }
    return derivative;
}

TEST(SourceSkyDifferential, KerrPositionAndTangentVariationsMatchIndependentDifferences) {
    for (double spin : {-0.85, 0.85}) {
        KerrSchildFamily authority(KerrSchildParams::KerrNewman(1.0, spin, 0.2));
        for (const Vec4 position : {Vector(2.0, 4.0, -3.0, 2.0), Vector(-7.0, 12.0, 5.0, -7.0)}) {
            const auto initial = Evaluate(authority, position, Vec4{});
            const auto reference =
                EulerianObserverFrame(initial.metric, initial.inverse, CartesianSeeds());
            ASSERT_TRUE(reference);
            const auto camera = BoostObserverFrame(*reference, {0.27, -0.11, 0.08});
            ASSERT_TRUE(camera);
            const auto tangent = PastDirectedCameraRay(*camera, Unit({-0.4, 0.7, 0.2}));
            ASSERT_TRUE(tangent);
            const Vec4 displacement = Vector(0.4, 0.31, -0.27, 0.19);
            const Vec4 tangent_variation = Vector(0.13, -0.07, 0.11, 0.05);
            for (int variation = 0; variation < 3; ++variation) {
                SCOPED_TRACE(::testing::Message() << "spin=" << spin << " x=" << position(1)
                                                  << " variation=" << variation);
                const Vec4 dx = variation == 1 ? Vec4{} : displacement;
                const Vec4 dk = variation == 0 ? Vec4{} : tangent_variation;
                const auto sample = Evaluate(authority, position, dx);
                const auto actual =
                    EulerianSourceSkyDifferential(sample.metric, sample.inverse, sample.variation,
                                                  *tangent, dk, CartesianSeeds());
                ASSERT_TRUE(actual);
                const auto nominal = NominalSource(authority, position, *tangent);
                ASSERT_TRUE(nominal);
                const auto coarse = DifferenceOracle(authority, position, *tangent, dx, dk, 1e-3);
                const auto fine = DifferenceOracle(authority, position, *tangent, dx, dk, 5e-4);
                ASSERT_TRUE(coarse);
                ASSERT_TRUE(fine);
                for (int axis = 0; axis < 3; ++axis) {
                    EXPECT_NEAR(actual->direction[axis], (*nominal)[axis], 2e-13);
                    EXPECT_NEAR((*coarse)[axis], (*fine)[axis], 2e-9);
                    EXPECT_NEAR(actual->derivative[axis], (*fine)[axis], 2e-9);
                }
                EXPECT_NEAR(actual->frequency, (*nominal)[3], 2e-13);
                EXPECT_NEAR((*coarse)[3], (*fine)[3], 2e-9);
                EXPECT_NEAR(actual->frequency_derivative, (*fine)[3], 2e-9);
                EXPECT_NEAR(Dot(actual->direction, actual->derivative), 0.0, 2e-13);
            }
        }
    }
}

TEST(SourceSkyDifferential, MovingFlatCameraMatchesAnalyticLorentzDirectionDerivative) {
    KerrSchildFamily authority(KerrSchildParams::Minkowski());
    const Direction beta{0.3, -0.2, 0.1};
    const Direction n = Unit({-0.4, 0.7, 0.2});
    Direction e{0.2, 0.1, -0.6};
    const double along = Dot(e, n);
    for (int axis = 0; axis < 3; ++axis) e[axis] -= along * n[axis];
    e = Unit(e);
    const double beta_squared = Dot(beta, beta);
    const double gamma = 1.0 / std::sqrt(1.0 - beta_squared);
    const double omega = gamma * (1.0 - Dot(beta, n));
    const double domega = -gamma * Dot(beta, e);
    Direction spatial{}, variation{}, expected{}, derivative{};
    for (int axis = 0; axis < 3; ++axis) {
        spatial[axis] =
            n[axis] + ((gamma - 1.0) * Dot(beta, n) / beta_squared - gamma) * beta[axis];
        variation[axis] = e[axis] + (gamma - 1.0) * Dot(beta, e) / beta_squared * beta[axis];
        expected[axis] = spatial[axis] / omega;
        derivative[axis] = variation[axis] / omega - spatial[axis] * domega / (omega * omega);
    }
    for (const Vec4 position : {Vector(0, 1, 2, 3), Vector(37, -120, 45, 71)}) {
        const auto sample = Evaluate(authority, position, Vector(0.7, -2.0, 3.0, 1.0));
        const auto frame = EulerianObserverFrame(sample.metric, sample.inverse, CartesianSeeds());
        ASSERT_TRUE(frame);
        const auto camera = BoostObserverFrame(*frame, beta);
        ASSERT_TRUE(camera);
        const auto k = PastDirectedCameraRay(*camera, n);
        ASSERT_TRUE(k);
        Vec4 dk;
        for (int axis = 0; axis < 3; ++axis) dk += camera->spatial[axis] * e[axis];
        const auto actual = EulerianSourceSkyDifferential(
            sample.metric, sample.inverse, sample.variation, *k, dk, CartesianSeeds());
        ASSERT_TRUE(actual);
        for (int axis = 0; axis < 3; ++axis) {
            EXPECT_NEAR(actual->direction[axis], expected[axis], 3e-13);
            EXPECT_NEAR(actual->derivative[axis], derivative[axis], 3e-13);
        }
        EXPECT_NEAR(actual->frequency, omega, 3e-13);
        EXPECT_NEAR(actual->frequency_derivative, domega, 3e-13);
    }
}

TEST(SourceSkyDifferential, NormalizesNonNullSpatialDirectionAndItsDerivative) {
    KerrSchildFamily authority(KerrSchildParams::Minkowski());
    const auto sample = Evaluate(authority, Vector(0, 2, 3, 4), Vec4{});
    const auto result = EulerianSourceSkyDifferential(
        sample.metric, sample.inverse, sample.variation, Vector(-2, 3, 4, 0),
        Vector(0.7, 0.3, -0.2, 0.4), CartesianSeeds());
    ASSERT_TRUE(result);
    const Direction expected{0.6, 0.8, 0.0};
    const Direction derivative{0.0576, -0.0432, 0.08};
    for (int axis = 0; axis < 3; ++axis) {
        EXPECT_NEAR(result->direction[axis], expected[axis], 2e-14);
        EXPECT_NEAR(result->derivative[axis], derivative[axis], 2e-14);
    }
    EXPECT_DOUBLE_EQ(result->frequency, 2.0);
    EXPECT_DOUBLE_EQ(result->frequency_derivative, -0.7);
}

TEST(SourceSkyDifferential, DirectionIsInvariantUnderFrequencyAndLongitudinalVariation) {
    KerrSchildFamily authority(KerrSchildParams::Kerr(1.0, 0.7));
    const Vec4 position = Vector(0, 6, -3, 2);
    const auto sample = Evaluate(authority, position, Vector(0.2, 0.3, -0.1, 0.4));
    const auto frame = EulerianObserverFrame(sample.metric, sample.inverse, CartesianSeeds());
    ASSERT_TRUE(frame);
    const auto k = PastDirectedCameraRay(*frame, Unit({-0.3, 0.4, 0.8}));
    ASSERT_TRUE(k);
    const Vec4 dk = Vector(0.1, -0.07, 0.11, 0.05);
    const auto reference = EulerianSourceSkyDifferential(
        sample.metric, sample.inverse, sample.variation, *k, dk, CartesianSeeds());
    ASSERT_TRUE(reference);
    for (double scale : {1e-6, 1.0, 1e6}) {
        for (double longitudinal : {-2.0, 0.0, 3.0}) {
            SCOPED_TRACE(::testing::Message() << scale << ", " << longitudinal);
            const auto actual = EulerianSourceSkyDifferential(
                sample.metric, sample.inverse, sample.variation, *k * scale,
                (dk + *k * longitudinal) * scale, CartesianSeeds());
            ASSERT_TRUE(actual);
            for (int axis = 0; axis < 3; ++axis) {
                EXPECT_NEAR(actual->direction[axis], reference->direction[axis], 3e-13);
                EXPECT_NEAR(actual->derivative[axis], reference->derivative[axis], 3e-13);
            }
            EXPECT_NEAR(actual->frequency / scale, reference->frequency, 3e-13);
            EXPECT_NEAR(actual->frequency_derivative / scale,
                        reference->frequency_derivative + longitudinal * reference->frequency,
                        3e-13);
        }
    }
}

TEST(SourceSkyDifferential, ReseedsMetricVariationWithoutInheritedDualMetadata) {
    KerrSchildFamily authority(KerrSchildParams::Kerr(1.0, -0.7));
    auto sample = Evaluate(authority, Vector(0, 4, -3, 2), Vector(0.1, 0.2, -0.3, 0.4));
    const auto frame = EulerianObserverFrame(sample.metric, sample.inverse, CartesianSeeds());
    ASSERT_TRUE(frame);
    const auto k = PastDirectedCameraRay(*frame, Unit({-0.3, 0.4, 0.8}));
    ASSERT_TRUE(k);
    const Vec4 dk = Vector(0.1, -0.07, 0.11, 0.05);
    const auto reference = EulerianSourceSkyDifferential(
        sample.metric, sample.inverse, sample.variation, *k, dk, CartesianSeeds());
    ASSERT_TRUE(reference);
    // Metric4d is reused by other automatic-differentiation consumers. Its
    // existing dual lane is not the directional derivative supplied to this API.
    for (int mu = 0; mu < 4; ++mu) {
        for (int nu = 0; nu < 4; ++nu) {
            sample.metric(mu, nu).dual = std::numeric_limits<double>::quiet_NaN();
            sample.inverse(mu, nu).dual = 1234.0;
            sample.variation(mu, nu).dual = -9876.0;
        }
    }
    const auto actual = EulerianSourceSkyDifferential(sample.metric, sample.inverse,
                                                      sample.variation, *k, dk, CartesianSeeds());
    ASSERT_TRUE(actual);
    EXPECT_EQ(actual->direction, reference->direction);
    EXPECT_EQ(actual->derivative, reference->derivative);
    EXPECT_DOUBLE_EQ(actual->frequency, reference->frequency);
    EXPECT_DOUBLE_EQ(actual->frequency_derivative, reference->frequency_derivative);
}

TEST(SourceSkyDifferential, DeclinesUnrepresentedFrequencyGeometryAndVariation) {
    KerrSchildFamily authority(KerrSchildParams::Minkowski());
    const auto sample = Evaluate(authority, Vector(0, 2, 3, 4), Vec4{});
    const auto evaluate = [&](const Metric4d& metric, const Metric4d& inverse,
                              const Metric4d& derivative, const Vec4& k, const Vec4& dk,
                              const std::array<Vec4, 3>& seeds) {
        return EulerianSourceSkyDifferential(metric, inverse, derivative, k, dk, seeds);
    };
    for (const Vec4 invalid : {Vector(0, 1, 0, 0), Vector(1, 1, 0, 0), Vector(-1, 0, 0, 0)}) {
        EXPECT_FALSE(evaluate(sample.metric, sample.inverse, sample.variation, invalid, Vec4{},
                              CartesianSeeds()));
    }
    auto degenerate_seeds = CartesianSeeds();
    degenerate_seeds[1] = degenerate_seeds[0];
    EXPECT_FALSE(evaluate(sample.metric, sample.inverse, sample.variation, Vector(-1, 1, 0, 0),
                          Vec4{}, degenerate_seeds));
    const double nan = std::numeric_limits<double>::quiet_NaN();
    EXPECT_FALSE(evaluate(sample.metric, sample.inverse, sample.variation, Vector(-1, 1, 0, 0),
                          Vector(0, nan, 0, 0), CartesianSeeds()));
    auto invalid_derivative = sample.variation;
    invalid_derivative(2, 1) = Dual<double>(nan);
    EXPECT_FALSE(evaluate(sample.metric, sample.inverse, invalid_derivative, Vector(-1, 1, 0, 0),
                          Vec4{}, CartesianSeeds()));
    auto invalid_inverse = sample.inverse;
    invalid_inverse(0, 0) = Dual<double>(1.0);
    EXPECT_FALSE(evaluate(sample.metric, invalid_inverse, sample.variation, Vector(-1, 1, 0, 0),
                          Vec4{}, CartesianSeeds()));
}

}  // namespace
}  // namespace sirius::test
