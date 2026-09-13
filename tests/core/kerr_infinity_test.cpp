// Independent Cartesian flat rays and scalar Carter equations are the oracles.
// No production tail RHS, initial-state extraction or radial validator is used
// to construct reference directions or their derivatives.
#include "sirius/core/kerr_infinity.h"

#include "sirius/core/metrics/kerr_schild_family.h"

#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <limits>
#include <numbers>

namespace {
using namespace sirius::core;
using namespace sirius::core::relativity;
using Vector = std::array<double, 3>;

double Dot(const Vector& a, const Vector& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
Vector Unit(Vector a) {
    const double length = std::hypot(a[0], a[1], a[2]);
    for (double& x : a) x /= length;
    return a;
}
Vector Cross(const Vector& a, const Vector& b) {
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}
std::array<Vector, 2> Basis(const Vector& n) {
    unsigned least = 0;
    for (unsigned i = 1; i < 3; ++i)
        if (std::abs(n[i]) < std::abs(n[least])) least = i;
    Vector e{};
    for (unsigned i = 0; i < 3; ++i) e[i] = (i == least ? 1.0 : 0.0) - n[least] * n[i];
    e = Unit(e);
    return {e, Cross(n, e)};
}
Vec4 Four(double t, const Vector& x) {
    Vec4 result;
    result(0) = t;
    for (unsigned i = 0; i < 3; ++i) result(i + 1) = x[i];
    return result;
}
Metric4d FlatMetric() {
    Metric4d metric;
    metric(0, 0).real = -1;
    for (int i = 1; i < 4; ++i) metric(i, i).real = 1;
    return metric;
}
KerrInfinityConfig Accurate() {
    KerrInfinityConfig config;
    config.absolute_tolerance = 1e-12;
    config.relative_tolerance = 1e-11;
    return config;
}

TEST(KerrInfinity, FlatOblateRaysAndAxesReachExactCartesianSky) {
    const auto metric = FlatMetric();
    const Tensor<Dual<double>, 4, 4, 4> dg;
    // Nonzero spin at zero mass is an oblate-coordinate control, not a claim
    // that the live KerrSchildFamily accepts a massless spinning black hole.
    for (double a : {-0.9, 0.0, 0.9}) {
        for (Vector x :
             {Vector{20, -7, 5}, Vector{-4, 18, 8}, Vector{0, 0, 20}, Vector{0, 0, -20}}) {
            Vector n = Unit({0.8, 0.3, std::copysign(0.6, x[2])});
            const auto result = TraceKerrInfinity(0, a, metric, dg, Four(0.2, x), Four(-1, n), {},
                                                  {}, 1e-3, Accurate());
            ASSERT_TRUE(result.has_value()) << a;
            for (unsigned i = 0; i < 3; ++i) EXPECT_NEAR(result->map.direction[i], n[i], 2e-10);
            EXPECT_DOUBLE_EQ(result->map.determinant, 0);
            EXPECT_DOUBLE_EQ(result->frequency, 1);
            EXPECT_GT(result->accepted_steps, 0u);
            EXPECT_GE(result->attempted_steps, result->accepted_steps);
            EXPECT_LE(result->maximum_local_error_ratio, 1);
        }
    }
}

TEST(KerrInfinity, VaryingHandoffAndTangentMatchExactFlatDerivative) {
    const Vector x{20, -7, 5}, k{0.8, 0.3, 0.6};
    const double frequency = std::sqrt(Dot(k, k));
    const Vector n = Unit(k);
    const auto screen = Basis(n);
    const std::array<Vector, 2> dx{Vector{1, 0.2, -0.3}, Vector{-0.7, 0.4, 0.1}};
    const std::array<Vector, 2> dk{Vector{0.03, -0.02, 0.015}, Vector{-0.02, 0.04, 0.03}};
    for (double seed : {1e-6, 1e-3, 0.2}) {
        std::array<Vec4, 2> X, K;
        std::array<Vector, 2> expected;
        for (unsigned column = 0; column < 2; ++column) {
            X[column] = Four(0.37, dx[column]) * seed;
            K[column] = Four(-Dot(n, dk[column]), dk[column]) * seed;
            for (unsigned i = 0; i < 3; ++i)
                expected[column][i] = (dk[column][i] - n[i] * Dot(n, dk[column])) / frequency;
        }
        const auto result = TraceKerrInfinity(0, 0.9, FlatMetric(), {}, Four(0, x),
                                              Four(-frequency, k), X, K, seed, Accurate());
        ASSERT_TRUE(result.has_value());
        for (unsigned column = 0; column < 2; ++column) {
            for (unsigned row = 0; row < 2; ++row)
                EXPECT_NEAR(result->map.jacobian[row][column], Dot(screen[row], expected[column]),
                            2e-10);
            EXPECT_NEAR(result->frequency_derivative[column], Dot(n, dk[column]), 2e-14);
        }
    }
}

TEST(KerrInfinity, FrequencyScaleAndSignedParityArePreserved) {
    const Vector x{20, -7, 5}, n = Unit({0.8, 0.3, 0.6});
    const auto screen = Basis(n);
    constexpr double seed = 1e-3;
    for (double frequency : {1e-6, 1.0, 1e6}) {
        const Vec4 k = Four(-1, n) * frequency;
        std::array<Vec4, 2> K;
        for (unsigned column = 0; column < 2; ++column) {
            const double parity = column == 0 ? 1.0 : -1.0;
            K[column] = Four(0, screen[column]) * (parity * frequency * seed) +
                        k * ((column == 0 ? 0.7 : -0.2) * seed);
        }
        const auto result =
            TraceKerrInfinity(0, -0.9, FlatMetric(), {}, Four(0, x), k, {}, K, seed, Accurate());
        ASSERT_TRUE(result.has_value());
        EXPECT_NEAR(result->map.jacobian[0][0], 1, 2e-9);
        EXPECT_NEAR(result->map.jacobian[1][1], -1, 2e-9);
        EXPECT_NEAR(result->map.jacobian[1][0], 0, 2e-9);
        EXPECT_NEAR(result->map.jacobian[0][1], 0, 2e-9);
        EXPECT_NEAR(result->map.determinant, -1, 4e-9);
        EXPECT_NEAR(result->frequency, frequency, frequency * 1e-14);
        EXPECT_NEAR(result->frequency_derivative[0], 0.7 * frequency, frequency * 1e-14);
        EXPECT_NEAR(result->frequency_derivative[1], -0.2 * frequency, frequency * 1e-14);
    }
}

TEST(KerrInfinity, AxisPrincipalRaysHaveAnalyticAngularVariation) {
    // On either principal axis, E=-1,L=0,Q=-a^2. The linear transverse
    // Carter system integrates to dn_inf=(r+i*a)/(r^2+a^2) dv_initial.
    // The covector conversion gives dv_initial=(r-i*a) dk_transverse,
    // so the exact limiting transverse derivative is dk_transverse.
    for (double a : {-1.0, -0.9, 0.0, 0.9, 1.0}) {
        KerrSchildFamily authority(KerrSchildParams::Kerr(1, a));
        for (double sign : {-1.0, 1.0}) {
            for (double radius : {3.0, 20.0}) {
                const Vector n{0, 0, sign};
                const Vec4 x = Four(0, {0, 0, sign * radius}), k = Four(-1, n);
                Metric4d metric;
                Tensor<Dual<double>, 4, 4, 4> dg;
                authority.Evaluate(x, metric, dg);
                const auto screen = Basis(n);
                const std::array<Vec4, 2> K{Four(0, screen[0]) * 1e-3, Four(0, screen[1]) * 1e-3};
                const auto result =
                    TraceKerrInfinity(1, a, metric, dg, x, k, {}, K, 1e-3, Accurate());
                ASSERT_TRUE(result.has_value()) << a << " " << sign << " " << radius;
                for (unsigned i = 0; i < 3; ++i) EXPECT_NEAR(result->map.direction[i], n[i], 1e-14);
                for (unsigned row = 0; row < 2; ++row)
                    for (unsigned col = 0; col < 2; ++col)
                        EXPECT_NEAR(result->map.jacobian[row][col], row == col ? 1 : 0, 2e-10);
            }
        }
    }
}

struct CarterInput {
    double mass, spin, radius, theta, phi, energy, angular_momentum, polar_momentum;
};
struct Event {
    Vec4 x, k;
};
CarterInput Vary(CarterInput c, const std::array<double, 6>& d, double amount) {
    c.radius += amount * d[0];
    c.theta += amount * d[1];
    c.phi += amount * d[2];
    c.energy += amount * d[3];
    c.angular_momentum += amount * d[4];
    c.polar_momentum += amount * d[5];
    return c;
}
Event CarterEvent(const CarterInput& c) {
    const double r = c.radius, a = c.spin, E = c.energy, L = c.angular_momentum;
    const double st = std::sin(c.theta), ct = std::cos(c.theta);
    const Vector n{st * std::cos(c.phi), st * std::sin(c.phi), ct};
    const Vector et{ct * std::cos(c.phi), ct * std::sin(c.phi), -st};
    const Vector ep{-std::sin(c.phi), std::cos(c.phi), 0};
    const auto transform = [r, a](const Vector& vector) {
        return Vector{r * vector[0] - a * vector[1], r * vector[1] + a * vector[0], r * vector[2]};
    };
    const double delta = r * r - 2 * c.mass * r + a * a, sigma = r * r + a * a * ct * ct;
    const double Q =
        c.polar_momentum * c.polar_momentum + ct * ct * (L * L / (st * st) - a * a * E * E);
    const double P = E * (r * r + a * a) - a * L;
    const double radial = std::sqrt(P * P - delta * (Q + (L - a * E) * (L - a * E)));
    const double phi_rate = L / (st * st) - a * E + a * P / delta + a * radial / delta;
    Vector angular;
    for (unsigned i = 0; i < 3; ++i) angular[i] = et[i] * c.polar_momentum + ep[i] * st * phi_rate;
    Vector spatial = transform(angular);
    for (unsigned i = 0; i < 3; ++i) spatial[i] = (spatial[i] + n[i] * radial) / sigma;
    const double time = ((r * r + a * a) * P / delta + a * (L - a * E * st * st) +
                         2 * c.mass * r * radial / delta) /
                        sigma;
    return {Four(0, transform(n)), Four(time, spatial)};
}

// Independent fixed-step, long-double scalar Carter integration. In particular
// this does not share the product's sphere-vector RHS or adaptive integrator.
Vector Reference(const CarterInput& c, unsigned steps = 4096) {
    using State = std::array<long double, 3>;
    const long double a = c.spin, E = c.energy, L = c.angular_momentum, M = c.mass;
    const long double ct = std::cos(static_cast<long double>(c.theta));
    const long double st = std::sin(static_cast<long double>(c.theta));
    const long double Q =
        c.polar_momentum * c.polar_momentum + ct * ct * (L * L / (st * st) - a * a * E * E);
    const auto f = [&](long double u, const State& y) {
        const long double sine = std::sin(y[0]), cosine = std::cos(y[0]);
        const long double D = 1 - 2 * M * u + a * a * u * u, P = E + (a * a * E - a * L) * u * u;
        const long double root = std::sqrt(P * P - D * u * u * (Q + (L - a * E) * (L - a * E)));
        return State{-y[1] / root,
                     (a * a * E * E * sine * cosine - L * L * cosine / (sine * sine * sine)) / root,
                     -(L / (sine * sine) - a * E + a * P / D) / root - a / D};
    };
    const auto add = [](State y, const State& d, long double h) {
        for (unsigned i = 0; i < 3; ++i) y[i] += h * d[i];
        return y;
    };
    State y{c.theta, c.polar_momentum, c.phi};
    const long double u0 = 1 / static_cast<long double>(c.radius), h = -u0 / steps;
    for (unsigned step = 0; step < steps; ++step) {
        const long double u = u0 + step * h;
        const auto k1 = f(u, y), k2 = f(u + h / 2, add(y, k1, h / 2));
        const auto k3 = f(u + h / 2, add(y, k2, h / 2)), k4 = f(u + h, add(y, k3, h));
        for (unsigned i = 0; i < 3; ++i) y[i] += h * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) / 6;
    }
    return {static_cast<double>(std::sin(y[0]) * std::cos(y[2])),
            static_cast<double>(std::sin(y[0]) * std::sin(y[2])),
            static_cast<double>(std::cos(y[0]))};
}

TEST(KerrInfinity, KerrTailMatchesIndependentSeparatedReference) {
    const std::array<std::array<double, 6>, 2> variation{
        std::array<double, 6>{0.7, 0.03, -0.02, 0.04, 0.11, -0.08},
        std::array<double, 6>{-0.4, -0.02, 0.04, -0.03, -0.05, 0.09}};
    for (double spin : {-1.0, -0.9, 0.0, 0.9, 1.0}) {
        const CarterInput input{1, spin, 8, 0.9, 1.1, -1, 2, 0.3};
        const Event event = CarterEvent(input);
        KerrSchildFamily authority(KerrSchildParams::Kerr(1, spin));
        Metric4d metric;
        Tensor<Dual<double>, 4, 4, 4> dg;
        authority.Evaluate(event.x, metric, dg);
        std::array<Vec4, 2> X, K;
        constexpr double h = 1e-4, seed = 1e-3;
        for (unsigned column = 0; column < 2; ++column) {
            const Event mm = CarterEvent(Vary(input, variation[column], -2 * h));
            const Event m = CarterEvent(Vary(input, variation[column], -h));
            const Event p = CarterEvent(Vary(input, variation[column], h));
            const Event pp = CarterEvent(Vary(input, variation[column], 2 * h));
            X[column] = (mm.x - m.x * 8 + p.x * 8 - pp.x) * (seed / (12 * h));
            K[column] = (mm.k - m.k * 8 + p.k * 8 - pp.k) * (seed / (12 * h));
        }
        const auto result =
            TraceKerrInfinity(1, spin, metric, dg, event.x, event.k, X, K, seed, Accurate());
        ASSERT_TRUE(result.has_value()) << spin;
        const auto reference = Reference(input), finer = Reference(input, 8192);
        const auto screen = Basis(reference);
        for (unsigned i = 0; i < 3; ++i) {
            EXPECT_NEAR(reference[i], finer[i], 2e-12);
            EXPECT_NEAR(result->map.direction[i], reference[i], 2e-10);
        }
        for (unsigned column = 0; column < 2; ++column) {
            std::array<Vector, 2> derivatives;
            for (unsigned refinement = 0; refinement < 2; ++refinement) {
                const double dh = refinement == 0 ? 2e-4 : 1e-4;
                const auto m = Reference(Vary(input, variation[column], -dh));
                const auto p = Reference(Vary(input, variation[column], dh));
                for (unsigned i = 0; i < 3; ++i)
                    derivatives[refinement][i] = (p[i] - m[i]) / (2 * dh);
            }
            for (unsigned i = 0; i < 3; ++i)
                EXPECT_NEAR(derivatives[0][i], derivatives[1][i], 2e-9);
            for (unsigned row = 0; row < 2; ++row)
                EXPECT_NEAR(result->map.jacobian[row][column], Dot(screen[row], derivatives[1]),
                            2e-9);
            EXPECT_NEAR(result->frequency_derivative[column], -variation[column][3], 2e-10);
        }
    }
}

TEST(KerrInfinity, RejectsHiddenRadialTurningPoint) {
    // r=2.5 is inside the photon sphere. This ray initially moves outward,
    // but b=5.3 exceeds sqrt(27), producing a forbidden radial interval before
    // infinity. H is positive at both endpoints and negative at r=3.
    const CarterInput input{1, 0, 2.5, std::numbers::pi / 2, 0, -1, 5.3, 0};
    const Event event = CarterEvent(input);
    EXPECT_GT(event.k(1), 0);
    EXPECT_GT(1 - 5.3 * 5.3 / (2.5 * 2.5) * (1 - 2 / 2.5), 0);
    EXPECT_LT(1 - 5.3 * 5.3 / 27, 0);
    KerrSchildFamily authority(KerrSchildParams::Schwarzschild(1));
    Metric4d metric;
    Tensor<Dual<double>, 4, 4, 4> dg;
    authority.Evaluate(event.x, metric, dg);
    const auto result = TraceKerrInfinity(1, 0, metric, dg, event.x, event.k, {}, {}, 1e-3);
    ASSERT_FALSE(result.has_value());
    EXPECT_EQ(result.error(), KerrInfinityFailure::RadialUnresolved);
}

TEST(KerrInfinity, MassSpinAndRadiusScalingPreserveDirectionAndJacobian) {
    std::optional<KerrInfinityResult> reference;
    for (double length : {1.0, 1e-3, 1e3}) {
        const CarterInput input{length, -0.9 * length, 8 * length,  0.9, 1.1,
                                -1,     2 * length,    0.3 * length};
        const Event event = CarterEvent(input);
        KerrSchildFamily authority(KerrSchildParams::Kerr(input.mass, input.spin));
        Metric4d metric;
        Tensor<Dual<double>, 4, 4, 4> dg;
        authority.Evaluate(event.x, metric, dg);
        std::array<Vec4, 2> X, K;
        const std::array<std::array<double, 6>, 2> variation{
            std::array<double, 6>{0.7 * length, 0.03, -0.02, 0.04, 0.11 * length, -0.08 * length},
            std::array<double, 6>{-0.4 * length, -0.02, 0.04, -0.03, -0.05 * length,
                                  0.09 * length}};
        constexpr double h = 1e-4, seed = 1e-3;
        for (unsigned column = 0; column < 2; ++column) {
            const auto mm = CarterEvent(Vary(input, variation[column], -2 * h));
            const auto m = CarterEvent(Vary(input, variation[column], -h));
            const auto p = CarterEvent(Vary(input, variation[column], h));
            const auto pp = CarterEvent(Vary(input, variation[column], 2 * h));
            X[column] = (mm.x - m.x * 8 + p.x * 8 - pp.x) * (seed / (12 * h));
            K[column] = (mm.k - m.k * 8 + p.k * 8 - pp.k) * (seed / (12 * h));
        }
        const auto result = TraceKerrInfinity(input.mass, input.spin, metric, dg, event.x, event.k,
                                              X, K, seed, Accurate());
        ASSERT_TRUE(result.has_value()) << length;
        if (!reference) reference = *result;
        for (unsigned i = 0; i < 3; ++i)
            EXPECT_NEAR(result->map.direction[i], reference->map.direction[i], 2e-10);
        for (unsigned row = 0; row < 2; ++row)
            for (unsigned column = 0; column < 2; ++column)
                EXPECT_NEAR(result->map.jacobian[row][column], reference->map.jacobian[row][column],
                            2e-9);
        EXPECT_NEAR(result->map.determinant, reference->map.determinant, 2e-10);
    }
}

TEST(KerrInfinity, NearCriticalSchwarzschildDerivativeMatchesIndependentRadialIntegral) {
    // Equatorial Schwarzschild gives phi_inf=phi0+integral b/sqrt(H) du,
    // dphi_inf/db=integral H^(-3/2) du. This is an independent scalar
    // quadrature with no angular ODE or transported derivative state.
    constexpr double radius = 2.5, seed = 1e-3;
    const double b = std::sqrt(27.0) - 1e-3;
    const double f = 1 - 2 / radius, kr = std::sqrt(1 - f * b * b / (radius * radius));
    const Vec4 x = Four(0, {radius, 0, 0}),
               k = Four((-1 + 2 / radius * kr) / f, {kr, b / radius, 0});
    const double dkr = -f * b / (radius * radius * kr);
    std::array<Vec4, 2> K{};
    K[0] = Four((2 / radius) * dkr / f, {dkr, 1 / radius, 0}) * seed;
    KerrSchildFamily authority(KerrSchildParams::Schwarzschild(1));
    Metric4d metric;
    Tensor<Dual<double>, 4, 4, 4> dg;
    authority.Evaluate(x, metric, dg);
    const auto integral = [b](unsigned intervals) {
        std::array<long double, 2> sum{};
        const long double width = 1 / (static_cast<long double>(radius) * intervals);
        for (unsigned i = 0; i <= intervals; ++i) {
            const long double u = i * width, bb = b;
            const long double H = 1 - bb * bb * u * u * (1 - 2 * u), root = std::sqrt(H);
            const long double weight = (i == 0 || i == intervals) ? 1 : (i % 2 ? 4 : 2);
            sum[0] += weight * bb / root;
            sum[1] += weight / (H * root);
        }
        for (auto& value : sum) value *= width / 3;
        return sum;
    };
    const auto expected = integral(16384), finer = integral(32768);
    EXPECT_NEAR(static_cast<double>(expected[0]), static_cast<double>(finer[0]), 1e-10);
    EXPECT_NEAR(static_cast<double>(expected[1]), static_cast<double>(finer[1]), 1e-7);
    const double phi = static_cast<double>(finer[0]), slope = static_cast<double>(finer[1]);
    const Vector n{std::cos(phi), std::sin(phi), 0},
        derivative{-std::sin(phi) * slope, std::cos(phi) * slope, 0};
    const auto screen = Basis(n);
    const auto result = TraceKerrInfinity(1, 0, metric, dg, x, k, {}, K, seed, Accurate());
    ASSERT_TRUE(result.has_value());
    for (unsigned i = 0; i < 3; ++i) EXPECT_NEAR(result->map.direction[i], n[i], 2e-9);
    for (unsigned row = 0; row < 2; ++row)
        EXPECT_NEAR(result->map.jacobian[row][0], Dot(screen[row], derivative), 2e-6);
    EXPECT_DOUBLE_EQ(result->map.determinant, 0);
}

TEST(KerrInfinity, RejectsUnrepresentedInputsAndWorkExhaustion) {
    const Vec4 x = Four(0, {20, -7, 5}), k = Four(-1, Unit({0.8, 0.3, 0.6}));
    auto config = Accurate();
    config.maximum_attempts = 1;
    auto result = TraceKerrInfinity(0, 0.9, FlatMetric(), {}, x, k, {}, {}, 1e-3, config);
    ASSERT_FALSE(result.has_value());
    EXPECT_EQ(result.error(), KerrInfinityFailure::WorkLimit);
    result = TraceKerrInfinity(0, 0.9, FlatMetric(), {}, x, k * -1, {}, {}, 1e-3);
    ASSERT_FALSE(result.has_value());
    EXPECT_EQ(result.error(), KerrInfinityFailure::InvalidInput);
    auto nonnull = k;
    nonnull(0) = -2;
    result = TraceKerrInfinity(0, 0.9, FlatMetric(), {}, x, nonnull, {}, {}, 1e-3);
    ASSERT_FALSE(result.has_value());
    EXPECT_EQ(result.error(), KerrInfinityFailure::InvalidInput);
    std::array<Vec4, 2> bad_variation{};
    bad_variation[0](0) = 1;
    result = TraceKerrInfinity(0, 0.9, FlatMetric(), {}, x, k, {}, bad_variation, 1e-3);
    ASSERT_FALSE(result.has_value());
    EXPECT_EQ(result.error(), KerrInfinityFailure::InvalidInput);
    bad_variation[0](0) = std::numeric_limits<double>::quiet_NaN();
    result = TraceKerrInfinity(0, 0.9, FlatMetric(), {}, x, k, {}, bad_variation, 1e-3);
    ASSERT_FALSE(result.has_value());
    EXPECT_EQ(result.error(), KerrInfinityFailure::InvalidInput);
}
}  // namespace
