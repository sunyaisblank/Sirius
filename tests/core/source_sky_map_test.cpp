#include "sirius/core/source_sky_map.h"

#include "sirius/core/metrics/kerr_schild_family.h"
#include "sirius/core/metrics/outgoing_kerr_schild.h"

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
    Vec4 value;
    value(0) = t;
    value(1) = x;
    value(2) = y;
    value(3) = z;
    return value;
}

double Dot(const Direction& a, const Direction& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

Direction Unit(Direction n) {
    const double length = std::sqrt(Dot(n, n));
    for (double& component : n) component /= length;
    return n;
}

Vec4 Spatial(const Direction& n) { return Vector(0, n[0], n[1], n[2]); }

struct EvaluatedMetric {
    Metric4d g;
    Metric4d inverse;
    Tensor<Dual<double>, 4, 4, 4> derivatives;
};

EvaluatedMetric Evaluate(IMetric& authority, const Vec4& x) {
    EvaluatedMetric value;
    authority.Evaluate(x, value.g, value.derivatives);
    if (!authority.InverseMetric(x, value.inverse)) {
        value.inverse = TensorOps::Inverse(value.g);
    }
    return value;
}

std::array<Vec4, 3> Seeds() { return {Vector(0, 1, 0, 0), Vector(0, 0, 1, 0), Vector(0, 0, 0, 1)}; }

std::optional<Direction> NominalSky(IMetric& authority, const Vec4& x, const Vec4& k) {
    const auto metric = Evaluate(authority, x);
    const auto frame = EulerianObserverFrame(metric.g, metric.inverse, Seeds());
    if (!frame) return std::nullopt;
    const double omega = TensorOps::InnerProduct(k, frame->time, metric.g);
    if (!(omega > 0.0)) return std::nullopt;
    Direction direction{};
    for (int axis = 0; axis < 3; ++axis) {
        direction[axis] = TensorOps::InnerProduct(k, frame->spatial[axis], metric.g) / omega;
    }
    return Unit(direction);
}

// The quadratic intersection of an exact straight null ray with the sphere is
// independent of the endpoint-variation formula under test.
Vec4 FlatSphereEvent(const Vec4& origin, const Direction& direction, double radius) {
    const Direction a{origin(1), origin(2), origin(3)};
    const double along = Dot(a, direction);
    const double distance = -along + std::sqrt(along * along + radius * radius - Dot(a, a));
    return origin + Vector(-1, direction[0], direction[1], direction[2]) * distance;
}

TEST(SourceSkyMap, LocalizedFlatSphereVariationMatchesAnalyticRayFamily) {
    const Vec4 origin = Vector(3.0, 2.0, -1.0, 0.5);
    const Direction n = Unit({0.7, 0.2, -0.4});
    const auto basis = MakeCelestialTangentBasis(n);
    ASSERT_TRUE(basis);
    const Vec4 k = Vector(-1, n[0], n[1], n[2]);
    constexpr double seed = 1e-3;
    constexpr double h = 1e-5;
    for (double radius : {8.0, 40.0, 200.0}) {
        const Vec4 event = FlatSphereEvent(origin, n, radius);
        const double affine = origin(0) - event(0);
        for (const Direction e : {basis->first, basis->second}) {
            const Vec4 xi = Spatial(e) * (affine * seed);
            const auto actual = VarySphericalEndpoint(event, k, xi);
            ASSERT_TRUE(actual);
            Direction plus{}, minus{};
            for (int axis = 0; axis < 3; ++axis) {
                plus[axis] = std::cos(h) * n[axis] + std::sin(h) * e[axis];
                minus[axis] = std::cos(h) * n[axis] - std::sin(h) * e[axis];
            }
            const Vec4 expected =
                (FlatSphereEvent(origin, plus, radius) - FlatSphereEvent(origin, minus, radius)) /
                (2.0 * h);
            for (int component = 0; component < 4; ++component) {
                EXPECT_NEAR(actual->displacement(component) / seed, expected(component), 2e-8);
            }
            double radial = 0.0;
            for (int axis = 1; axis < 4; ++axis) radial += event(axis) * actual->displacement(axis);
            EXPECT_NEAR(radial, 0.0, 2e-12);
        }
    }
}

TEST(SourceSkyMap, EndpointVariationPreservesGaugeAndDeclinesGrazingContact) {
    const Vec4 x = Vector(0, 3, 4, 0);
    const Vec4 k = Vector(-1, 0.6, 0.8, 0);
    const Vec4 xi = Vector(0.2, 0.4, -0.1, 0.3);
    const auto reference = VarySphericalEndpoint(x, k, xi);
    ASSERT_TRUE(reference);
    for (double scale : {1e-6, 1.0, 1e6}) {
        for (double gauge : {-3.0, 0.0, 7.0}) {
            const auto actual = VarySphericalEndpoint(x, k * scale, xi + k * gauge);
            ASSERT_TRUE(actual);
            for (int component = 0; component < 4; ++component) {
                EXPECT_NEAR(actual->displacement(component), reference->displacement(component),
                            3e-14);
            }
        }
    }
    EXPECT_FALSE(VarySphericalEndpoint(x, Vector(-1, 0.8, -0.6, 0), xi));
    EXPECT_FALSE(VarySphericalEndpoint(
        Vector(0, 1, 1, 0), Vector(-1, 1, -1 + std::numeric_limits<double>::epsilon(), 0), xi));
    EXPECT_FALSE(VarySphericalEndpoint(Vec4{}, k, xi));
    EXPECT_FALSE(
        VarySphericalEndpoint(x, k, Vector(0, std::numeric_limits<double>::quiet_NaN(), 0, 0)));
}

TEST(SourceSkyMap, FlatAngularMapIsIndependentOfSphereRadiusSeedAndFrequency) {
    KerrSchildFamily authority(KerrSchildParams::Minkowski());
    const Vec4 origin = Vector(3, 2, -1, 0.5);
    const Direction n = Unit({0.7, 0.2, -0.4});
    const auto basis = MakeCelestialTangentBasis(n);
    ASSERT_TRUE(basis);
    const Vec4 k = Vector(-1, n[0], n[1], n[2]);
    for (double radius : {8.0, 200.0, 3200.0}) {
        const Vec4 event = FlatSphereEvent(origin, n, radius);
        const double affine = origin(0) - event(0);
        const auto metric = Evaluate(authority, event);
        for (double seed : {1e-6, 1e-3, 0.1}) {
            for (double frequency : {0.01, 1.0, 100.0}) {
                std::array<Vec4, 2> displacements{}, variations{};
                const std::array<Direction, 2> axes{basis->first, basis->second};
                for (int column = 0; column < 2; ++column) {
                    const auto boundary =
                        VarySphericalEndpoint(event, k, Spatial(axes[column]) * (affine * seed));
                    ASSERT_TRUE(boundary);
                    displacements[column] = boundary->displacement;
                    variations[column] = Spatial(axes[column]) * (seed * frequency);
                }
                const auto actual =
                    MeasureSourceSkyAngularMap(metric.g, metric.inverse, metric.derivatives,
                                               k * frequency, displacements, variations, seed);
                ASSERT_TRUE(actual);
                for (int row = 0; row < 2; ++row) {
                    for (int column = 0; column < 2; ++column) {
                        EXPECT_NEAR(actual->jacobian[row][column], row == column ? 1.0 : 0.0,
                                    3e-13);
                    }
                }
                EXPECT_NEAR(actual->determinant, 1.0, 4e-13);
            }
        }
    }
}

// Independent connection: finite differences of nominal metric values rather
// than the metric derivative tensor or TensorOps::Christoffel used by the API.
std::array<Vec4, 2> CovariantVariations(IMetric& authority, const Vec4& x, const Vec4& k,
                                        const std::array<Vec4, 2>& displacements,
                                        const std::array<Vec4, 2>& coordinate_variations) {
    const auto central = Evaluate(authority, x);
    double derivative[4][4][4]{};
    constexpr double h = 2e-4;
    for (int axis = 0; axis < 4; ++axis) {
        std::array<Metric4d, 4> values;
        const std::array<double, 4> offsets{-2.0, -1.0, 1.0, 2.0};
        for (int index = 0; index < 4; ++index) {
            Vec4 q = x;
            q(axis) += offsets[index] * h;
            values[index] = Evaluate(authority, q).g;
        }
        for (int mu = 0; mu < 4; ++mu) {
            for (int nu = 0; nu < 4; ++nu) {
                derivative[axis][mu][nu] = (values[0](mu, nu).real - 8.0 * values[1](mu, nu).real +
                                            8.0 * values[2](mu, nu).real - values[3](mu, nu).real) /
                                           (12.0 * h);
            }
        }
    }
    auto result = coordinate_variations;
    for (int column = 0; column < 2; ++column) {
        for (int mu = 0; mu < 4; ++mu) {
            for (int a = 0; a < 4; ++a) {
                for (int b = 0; b < 4; ++b) {
                    double gamma = 0.0;
                    for (int c = 0; c < 4; ++c) {
                        gamma += 0.5 * central.inverse(mu, c).real *
                                 (derivative[a][c][b] + derivative[b][c][a] - derivative[c][a][b]);
                    }
                    result[column](mu) += gamma * k(a) * displacements[column](b);
                }
            }
        }
    }
    return result;
}

TEST(SourceSkyMap, KerrCovariantChartVariationsMatchMappedNominalFamily) {
    for (double spin : {-0.85, 0.85}) {
        KerrSchildFamily incoming(KerrSchildParams::KerrNewman(1, spin, 0.2));
        OutgoingKerrSchild outgoing(incoming);
        const Vec4 x = Vector(0.3, 5.0, -3.0, 2.0);
        const auto local_metric = Evaluate(outgoing, x);
        const auto frame = EulerianObserverFrame(local_metric.g, local_metric.inverse, Seeds());
        ASSERT_TRUE(frame);
        const auto k = PastDirectedCameraRay(*frame, Unit({0.3, -0.7, 0.2}));
        ASSERT_TRUE(k);
        const std::array<Vec4, 2> dx{Vector(0.3, 0.2, -0.1, 0.05), Vector(-0.2, 0.1, 0.3, -0.15)};
        const std::array<Vec4, 2> dk{Vector(0.07, -0.1, 0.03, 0.06),
                                     Vector(-0.04, 0.08, -0.02, 0.09)};
        const auto covariant = CovariantVariations(outgoing, x, *k, dx, dk);
        const auto mapping = outgoing.ToIngoing(x);
        ASSERT_TRUE(mapping);
        const auto metric = Evaluate(incoming, mapping->position);
        constexpr double seed = 1e-3;
        std::array<Vec4, 2> mapped_dx{}, mapped_covariant{};
        for (int column = 0; column < 2; ++column) {
            mapped_dx[column] = mapping->Apply(dx[column]) * seed;
            mapped_covariant[column] = mapping->Apply(covariant[column]) * seed;
        }
        const auto actual =
            MeasureSourceSkyAngularMap(metric.g, metric.inverse, metric.derivatives,
                                       mapping->Apply(*k), mapped_dx, mapped_covariant, seed);
        ASSERT_TRUE(actual);
        const auto nominal = NominalSky(incoming, mapping->position, mapping->Apply(*k));
        ASSERT_TRUE(nominal);
        const auto basis = MakeCelestialTangentBasis(*nominal);
        ASSERT_TRUE(basis);
        for (int axis = 0; axis < 3; ++axis) {
            EXPECT_NEAR(actual->direction[axis], (*nominal)[axis], 2e-13);
        }
        for (int column = 0; column < 2; ++column) {
            for (double h : {1e-3, 5e-4}) {
                std::array<Direction, 4> values{};
                const std::array<double, 4> offsets{-2.0, -1.0, 1.0, 2.0};
                for (int index = 0; index < 4; ++index) {
                    const double offset = offsets[index] * h;
                    const auto varied_map = outgoing.ToIngoing(x + dx[column] * offset);
                    ASSERT_TRUE(varied_map);
                    const auto value = NominalSky(incoming, varied_map->position,
                                                  varied_map->Apply(*k + dk[column] * offset));
                    ASSERT_TRUE(value);
                    values[index] = *value;
                }
                Direction expected{};
                for (int axis = 0; axis < 3; ++axis) {
                    expected[axis] = (values[0][axis] - 8.0 * values[1][axis] +
                                      8.0 * values[2][axis] - values[3][axis]) /
                                     (12.0 * h);
                }
                EXPECT_NEAR(actual->jacobian[0][column], Dot(basis->first, expected), 2e-8);
                EXPECT_NEAR(actual->jacobian[1][column], Dot(basis->second, expected), 2e-8);
            }
        }
    }
}

TEST(SourceSkyMap, RetainsZeroRankDeficiencyAndSignedParity) {
    KerrSchildFamily authority(KerrSchildParams::Minkowski());
    const auto metric = Evaluate(authority, Vector(0, 2, 3, 4));
    const Vec4 k = Vector(-1, 0, 0, 1);
    const std::array<Vec4, 2> zero{};
    const auto constant = MeasureSourceSkyAngularMap(metric.g, metric.inverse, metric.derivatives,
                                                     k, zero, zero, 1.0);
    ASSERT_TRUE(constant);
    EXPECT_DOUBLE_EQ(constant->determinant, 0.0);
    for (const auto& row : constant->jacobian) {
        for (double value : row) EXPECT_DOUBLE_EQ(value, 0.0);
    }
    const auto rank_one =
        MeasureSourceSkyAngularMap(metric.g, metric.inverse, metric.derivatives, k, zero,
                                   {Vector(0, 1, 0, 0), Vector(0, 2, 0, 0)}, 1.0);
    ASSERT_TRUE(rank_one);
    EXPECT_DOUBLE_EQ(rank_one->determinant, 0.0);
    EXPECT_DOUBLE_EQ(rank_one->jacobian[0][1], 2.0);
    const auto reversed =
        MeasureSourceSkyAngularMap(metric.g, metric.inverse, metric.derivatives, k, zero,
                                   {Vector(0, 1, 0, 0), Vector(0, 0, -1, 0)}, 1.0);
    ASSERT_TRUE(reversed);
    EXPECT_DOUBLE_EQ(reversed->determinant, -1.0);
    for (double invalid : {0.0, -1.0, std::numeric_limits<double>::infinity()}) {
        EXPECT_FALSE(MeasureSourceSkyAngularMap(metric.g, metric.inverse, metric.derivatives, k,
                                                zero, zero, invalid));
    }
}

}  // namespace
}  // namespace sirius::test
