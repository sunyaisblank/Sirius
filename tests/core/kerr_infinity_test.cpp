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
Event CarterEvent(const CarterInput& c, double radial_sign = 1.0) {
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
    const double radial = radial_sign * std::sqrt(P * P - delta * (Q + (L - a * E) * (L - a * E)));
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

// Continue the scalar Carter system in Mino time with an independent radial
// velocity variable. It crosses the turning point without the product's root
// isolation, factorization, sphere-vector representation, or RK45 controller.
std::optional<Vector> ScatteringReference(const CarterInput& c, long double step) {
    using State = std::array<long double, 5>;
    const long double a = c.spin, E = c.energy, L = c.angular_momentum, M = c.mass;
    const long double ct = std::cos(static_cast<long double>(c.theta));
    const long double st = std::sin(static_cast<long double>(c.theta));
    const long double Q =
        c.polar_momentum * c.polar_momentum + ct * ct * (L * L / (st * st) - a * a * E * E);
    const long double C = Q + (L - a * E) * (L - a * E);
    const long double Z = a * a * E - a * L;
    const long double A = 2 * E * Z - C, B = 2 * M * C, D = Z * Z - a * a * C;
    const auto field = [&](const State& y) {
        const auto [theta, polar, phi, u, velocity] = y;
        const long double sine = std::sin(theta), cosine = std::cos(theta);
        const long double delta = 1 - 2 * M * u + a * a * u * u;
        const long double phi_rate =
            L / (sine * sine) - a * E + a * (E + Z * u * u - velocity) / delta;
        return State{polar, (L * L / (sine * sine * sine) - a * a * E * E * sine) * cosine,
                     phi_rate, velocity, u * (A + u * (1.5L * B + 2 * D * u))};
    };
    const auto add = [](State y, const State& d, long double dt) {
        for (unsigned i = 0; i < y.size(); ++i) y[i] += dt * d[i];
        return y;
    };
    const auto advance = [&](const State& y, long double dt) {
        const auto k1 = field(y), k2 = field(add(y, k1, dt / 2));
        const auto k3 = field(add(y, k2, dt / 2)), k4 = field(add(y, k3, dt));
        auto next = y;
        for (unsigned i = 0; i < y.size(); ++i)
            next[i] += dt * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) / 6;
        return next;
    };
    const long double u = 1 / static_cast<long double>(c.radius);
    State y{c.theta, c.polar_momentum, c.phi, u, std::sqrt(E * E + u * u * (A + u * (B + u * D)))};
    for (unsigned i = 0; i < 100000; ++i) {
        const auto next = advance(y, step);
        if (next[3] <= 0) {
            long double lo = 0, hi = step;
            for (unsigned refinement = 0; refinement < 48; ++refinement) {
                const auto mid = std::midpoint(lo, hi);
                if (advance(y, mid)[3] > 0)
                    lo = mid;
                else
                    hi = mid;
            }
            y = advance(y, std::midpoint(lo, hi));
            return Vector{static_cast<double>(std::sin(y[0]) * std::cos(y[2])),
                          static_cast<double>(std::sin(y[0]) * std::sin(y[2])),
                          static_cast<double>(std::cos(y[0]))};
        }
        y = next;
    }
    return {};
}

TEST(KerrInfinity, RadialTurnMatchesIndependentMinoTimeScatteringAndColumns) {
    const std::array<std::array<double, 6>, 2> variations{
        std::array<double, 6>{0.3, 0.02, -0.01, 0.03, 0.08, -0.05},
        std::array<double, 6>{-0.4, -0.02, 0.03, -0.02, -0.07, 0.09}};
    constexpr double seed = 1e-3, delta = 1e-4;
    struct ScatteringCase {
        CarterInput input;
        double minimum_radius, maximum_radius;
    };
    const double equator = std::numbers::pi / 2;
    const std::array cases{
        ScatteringCase{{1, -0.9, 20, 1.2, 1.1, -1, 8, 0.8}, 4, 20},
        ScatteringCase{{1, 0, 20, 1.2, 1.1, -1, 8, 0.8}, 4, 20},
        ScatteringCase{{1, 0.9, 20, 1.2, 1.1, -1, 8, 0.8}, 4, 20},
        ScatteringCase{{1, 0.9, 20, equator, 1.1, -1, -4, std::sqrt(0.45)}, 2.5, 3.1},
        ScatteringCase{{1, -0.9, 20, equator, 1.1, -1, 4, std::sqrt(0.45)}, 2.5, 3.1},
        // The potential is positive at both 2M and 4M and negative between
        // two roots, so endpoint-only probing will miss this scattering ray.
        ScatteringCase{{1, 0.7, 20, equator, 1.1, -1, -1, std::sqrt(28.6)}, 3.7, 4},
        ScatteringCase{{1, 0, 20, equator, 1.1, -1, 5.3, 0}, 3.4, 3.5},
        ScatteringCase{{1, 0.9, 20, equator, 1.1, -1, -2.9, 0}, 1.70, 1.72},
    };
    for (const auto& [input, min_radius, max_radius] : cases) {
        const double spin = input.spin;
        SCOPED_TRACE(spin);
        SCOPED_TRACE(input.angular_momentum);
        const auto event = CarterEvent(input, -1);
        KerrSchildFamily authority(KerrSchildParams::Kerr(1, spin));
        Metric4d metric;
        Tensor<Dual<double>, 4, 4, 4> dg;
        authority.Evaluate(event.x, metric, dg);
        std::array<Vec4, 2> X, K;
        for (unsigned column = 0; column < 2; ++column) {
            const auto mm = CarterEvent(Vary(input, variations[column], -2 * delta), -1);
            const auto m = CarterEvent(Vary(input, variations[column], -delta), -1);
            const auto p = CarterEvent(Vary(input, variations[column], delta), -1);
            const auto pp = CarterEvent(Vary(input, variations[column], 2 * delta), -1);
            X[column] = (mm.x - m.x * 8 + p.x * 8 - pp.x) * (seed / (12 * delta));
            K[column] = (mm.k - m.k * 8 + p.k * 8 - pp.k) * (seed / (12 * delta));
        }
        const auto ordinary = TraceKerrInfinity(1, spin, metric, dg, event.x, event.k, X, K, seed);
        ASSERT_FALSE(ordinary);
        EXPECT_EQ(ordinary.error(), KerrInfinityFailure::InvalidInput);
        auto config = Accurate();
        config.allow_radial_turn = true;
        double turn_radius = 0.0;
        const auto result = TraceKerrInfinity(1, spin, metric, dg, event.x, event.k, X, K, seed,
                                              config, &turn_radius);
        ASSERT_TRUE(result) << spin << ": " << static_cast<int>(result.error());
        EXPECT_GT(turn_radius, min_radius);
        EXPECT_LT(turn_radius, max_radius);
        const auto reference = ScatteringReference(input, 1e-4L);
        const auto finer = ScatteringReference(input, 5e-5L);
        ASSERT_TRUE(reference && finer);
        for (unsigned i = 0; i < 3; ++i) {
            EXPECT_NEAR((*reference)[i], (*finer)[i], 2e-12);
            EXPECT_NEAR(result->map.direction[i], (*finer)[i], 2e-10);
        }
        const auto screen = Basis(*finer);
        for (unsigned column = 0; column < 2; ++column) {
            const auto mm = ScatteringReference(Vary(input, variations[column], -2 * delta), 5e-5L);
            const auto m = ScatteringReference(Vary(input, variations[column], -delta), 5e-5L);
            const auto p = ScatteringReference(Vary(input, variations[column], delta), 5e-5L);
            const auto pp = ScatteringReference(Vary(input, variations[column], 2 * delta), 5e-5L);
            ASSERT_TRUE(mm && m && p && pp);
            Vector derivative;
            for (unsigned i = 0; i < 3; ++i)
                derivative[i] = ((*mm)[i] - 8 * (*m)[i] + 8 * (*p)[i] - (*pp)[i]) / (12 * delta);
            for (unsigned row = 0; row < 2; ++row)
                EXPECT_NEAR(result->map.jacobian[row][column], Dot(screen[row], derivative), 2e-9);
            EXPECT_NEAR(result->frequency_derivative[column], -variations[column][3], 2e-10);
        }
    }
}

TEST(KerrInfinity, RadialTurnDeclinesCaptureAndUnresolvedTurns) {
    const double equator = std::numbers::pi / 2;
    const std::array inputs{
        CarterInput{1, 0, 20, equator, 0, -1, 2.0, 0},
        CarterInput{1, 0, 20, equator, 0, -1, std::sqrt(27.0), 0},
        CarterInput{1, 0.9, 20, equator, 0, -1, -2.0, 0},
    };
    for (const auto& input : inputs) {
        const auto event = CarterEvent(input, -1);
        KerrSchildFamily authority(KerrSchildParams::Kerr(1, input.spin));
        Metric4d metric;
        Tensor<Dual<double>, 4, 4, 4> dg;
        authority.Evaluate(event.x, metric, dg);
        auto config = Accurate();
        config.allow_radial_turn = true;
        double turn_radius = -1.0;
        const auto result = TraceKerrInfinity(1, input.spin, metric, dg, event.x, event.k, {}, {},
                                              1e-3, config, &turn_radius);
        ASSERT_FALSE(result) << input.angular_momentum;
        EXPECT_EQ(result.error(), KerrInfinityFailure::RadialUnresolved);
        EXPECT_EQ(turn_radius, -1.0);

        const double impact = std::abs(input.angular_momentum);
        const auto parameters = authority.GetParams();
        EXPECT_EQ(
            sirius::core::relativity::CertifyKerrDetectorCapture(parameters, event.x, event.k, 64),
            impact == 2.0);
        const auto outward = CarterEvent(input, 1);
        EXPECT_FALSE(sirius::core::relativity::CertifyKerrDetectorCapture(parameters, outward.x,
                                                                          outward.k, 64));
        if (impact == 2.0) {
            using sirius::core::relativity::finite_tail_detail::Range;
            using sirius::core::relativity::finite_tail_detail::Widen;
            std::array<Range, 4> position_bounds, tangent_bounds;
            for (int i = 0; i < 4; ++i) {
                position_bounds[i] = Widen(event.x(i) - 1e-3, event.x(i) + 1e-3);
                tangent_bounds[i] = Widen(event.k(i) - 1e-3, event.k(i) + 1e-3);
            }
            EXPECT_TRUE(sirius::core::relativity::CertifyKerrDetectorCapture(
                parameters, position_bounds, tangent_bounds, 64));
            auto turning_input = input;
            turning_input.angular_momentum = input.spin == 0 ? 5.4 : -8.0;
            const auto turning = CarterEvent(turning_input, -1);
            for (int i = 0; i < 4; ++i) {
                auto& bounds = tangent_bounds[i];
                bounds.lo = std::min(bounds.lo, turning.k(i));
                bounds.hi = std::max(bounds.hi, turning.k(i));
            }
            EXPECT_FALSE(sirius::core::relativity::CertifyKerrDetectorCapture(
                parameters, position_bounds, tangent_bounds, 64));
        }
    }
    // The first rejected Radeon launches in the governed 1080p scene.
    // Hexadecimal literals retain the actual public ingoing binary64 state.
    const auto position = Four(0x1.24a6d18db4p-10,
                               {0x1.d884a21ad7b96p+4, 0x1.d04c7d7c1169p-1, 0x1.51733439618ap+2});
    const std::array tangents{
        Four(-0x1.deb0a674e7f3p-1,
             {-0x1.935b61994b898p-1, -0x1.858031d03bdb1p-4, -0x1.9f7b1fe6305p-3}),
        Four(-0x1.deb9d8491266p-1,
             {-0x1.93366da6c8fd1p-1, -0x1.7a9ddae0f2a43p-4, -0x1.a4992d0defep-3}),
    };
    const auto kerr = KerrSchildParams::Kerr(1, 0.9);
    KerrSchildFamily kerr_metric(kerr);
    Metric4d device_metric;
    Tensor<Dual<double>, 4, 4, 4> device_derivatives;
    kerr_metric.Evaluate(position, device_metric, device_derivatives);
    const std::array scattering{
        Four(-0x1.deb5b5c6e8e38p-1,
             {-0x1.92f7f407dc8c7p-1, -0x1.93ad4880817b2p-4, -0x1.a2f4ea78b027p-3}),
        Four(-0x1.de917a1a791bbp-1,
             {-0x1.942f58adb117fp-1, -0x1.941f5365512dfp-4, -0x1.8cbe679530dcp-3}),
        Four(-0x1.de914d60819c1p-1,
             {-0x1.94319243b7847p-1, -0x1.93f086614a59p-4, -0x1.8ca0f3bf482dp-3}),
        Four(-0x1.de9114dc37fafp-1,
             {-0x1.94334b86ad86p-1, -0x1.93fc4dac70a0bp-4, -0x1.8c7dbf5d817cp-3}),
        Four(-0x1.de90e141a4c6fp-1,
             {-0x1.9434cc8dfbae1p-1, -0x1.940b9b8b24eb6p-4, -0x1.8c5db9300432p-3}),
        Four(-0x1.de90c471408dp-1,
             {-0x1.9435b18a1a608p-1, -0x1.9410965fa77dp-4, -0x1.8c4bbc809922p-3}),
    };
    constexpr std::array decimal_radii{2.075103142542746, 1.712572390557077, 1.694091803789609,
                                      1.689357681104494, 1.686593027729126, 1.682462964901301};
    auto through_turn = Accurate();
    through_turn.allow_radial_turn = true;
    for (std::size_t i = 0; i < scattering.size(); ++i) {
        const auto& tangent = scattering[i];
        EXPECT_FALSE(CertifyKerrDetectorCapture(kerr, position, tangent, 64));
        double radius = 0.0;
        const auto resumed = TraceKerrInfinity(1, 0.9, device_metric, device_derivatives, position,
                                                tangent, {}, {}, 1e-3, through_turn, &radius);
        ASSERT_TRUE(resumed);
        EXPECT_NEAR(radius, decimal_radii[i], 2e-12);
    }
    through_turn.minimum_turn_radius = 2.0;
    const auto governed = TraceKerrInfinity(1, 0.9, device_metric, device_derivatives, position,
                                             scattering.back(), {}, {}, 1e-3, through_turn);
    ASSERT_FALSE(governed);
    EXPECT_EQ(governed.error(), KerrInfinityFailure::RadialUnresolved);
    for (const auto& tangent : tangents) {
        EXPECT_TRUE(CertifyKerrDetectorCapture(kerr, position, tangent, 64));
        EXPECT_FALSE(CertifyKerrDetectorCapture(kerr, position, tangent, 0));
    }

    // The exact Schwarzschild potential has its minimum at u=1/3 and
    // its sign is fixed by b^2 versus 27. Test both sides of the boundary.
    for (double shift : {-1e-9, 1e-9}) {
        const double impact = std::sqrt(27.0) * (1.0 + shift);
        const auto close = CarterEvent({1, 0, 20, equator, 0, -1, impact, 0}, -1);
        EXPECT_EQ(CertifyKerrDetectorCapture(KerrSchildParams::Schwarzschild(1), close.x, close.k,
                                             64),
                  1.0 - impact * impact / 27.0 > 0.0);
    }
}

TEST(KerrInfinity, RejectsHiddenRadialTurningPoint) {
    using kerr_infinity_detail::Constants;
    using kerr_infinity_detail::ExactRadialSign;
    Constants critical;
    critical.energy = kerr_infinity_detail::Jet(-1.0);
    // At u=1/4 in Schwarzschild, H=1-C/32. Adjacent binary64 values
    // of C must yield opposite signs, while C=32 gives exactly zero.
    for (const auto& [carter, sign] : std::array{std::pair{std::nextafter(32.0, 0.0), 1},
                                               std::pair{32.0, 0},
                                               std::pair{std::nextafter(32.0, 64.0), -1}}) {
        critical.carter_combination = kerr_infinity_detail::Jet(carter);
        EXPECT_EQ(ExactRadialSign(1, 0, critical, .25), sign);
    }
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
