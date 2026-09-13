#pragma once
// Independent analytic and scalar-Carter continuation oracle. No production
// infinity helper, angular-vector RHS or radial certificate supplies expected
// direction, Jacobian, frequency or status. The existing metric/connection is
// used only to adapt independent coordinate variations to covariant inputs.
//
// Carter equations: Gralla & Lupsasca (2020), arXiv:1910.12881, equations 3--9.
// Flat Cartesian and principal-axis references are analytic. Non-axis values
// use fixed-step scalar Carter integration and independently refined angular
// differences. These local angular tests do not establish detector accuracy.
#include "sirius/core/metrics/kerr_schild_family.h"
#include "sirius/core/tensor.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <numbers>
#include <stdexcept>
#include <string>
#include <vector>

namespace sirius::test::kerr_infinity_oracle {
struct Case {
    std::string name;
    double mass, spin, seed;
    sirius::core::Vec4 x, k;
    std::array<sirius::core::Vec4, 2> X, V;
    std::array<double, 3> direction;
    std::array<std::array<double, 2>, 2> jacobian;
    double frequency;
    std::array<double, 2> frequency_derivative;
    unsigned expected_status = 1;
    unsigned maximum_attempts = 512;
};
namespace detail {
using namespace sirius::core;
using Vector = std::array<double, 3>;
inline double Dot(const Vector& a, const Vector& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
inline Vector Unit(Vector a) {
    const double length = std::hypot(a[0], a[1], a[2]);
    for (double& x : a) x /= length;
    return a;
}
inline Vector Cross(const Vector& a, const Vector& b) {
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}
inline std::array<Vector, 2> Basis(const Vector& n) {
    unsigned least = 0;
    for (unsigned i = 1; i < 3; ++i)
        if (std::abs(n[i]) < std::abs(n[least])) least = i;
    Vector e{};
    for (unsigned i = 0; i < 3; ++i) e[i] = (i == least ? 1.0 : 0.0) - n[least] * n[i];
    e = Unit(e);
    return {e, Cross(n, e)};
}
inline Vec4 Four(double t, const Vector& x) {
    Vec4 result;
    result(0) = t;
    for (unsigned i = 0; i < 3; ++i) result(i + 1) = x[i];
    return result;
}
struct CarterInput {
    double mass, spin, radius, theta, phi, energy, angular_momentum, polar_momentum;
};
struct Event {
    Vec4 x, k;
};
inline CarterInput Vary(CarterInput c, const std::array<double, 6>& d, double amount) {
    c.radius += amount * d[0];
    c.theta += amount * d[1];
    c.phi += amount * d[2];
    c.energy += amount * d[3];
    c.angular_momentum += amount * d[4];
    c.polar_momentum += amount * d[5];
    return c;
}
inline Event CarterEvent(const CarterInput& c) {
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
inline Vector Reference(const CarterInput& c, unsigned steps = 4096) {
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

struct ReferenceAudit {
    std::string name;
    double nominal_step_refinement_error = 0;
    double angular_difference_refinement_error = 0;
};
inline double DeviceParameter(double value) {
    return static_cast<double>(static_cast<float>(value));
}

// The sole source-authority use is adapting an independently specified
// coordinate variation into the covariant variation accepted by the shader.
// The massless oblate fixture has exactly Cartesian Minkowski g and Gamma=0.
inline void CovariantInput(Case& c, const std::array<Vec4, 2>& K) {
    c.V = K;
    if (c.mass == 0) return;
    KerrSchildFamily authority(KerrSchildParams::Kerr(c.mass, c.spin));
    Metric4d metric;
    Tensor<Dual<double>, 4, 4, 4> dg;
    authority.Evaluate(c.x, metric, dg);
    const auto connection = TensorOps::Christoffel(metric, dg);
    for (unsigned column = 0; column < 2; ++column)
        for (unsigned mu = 0; mu < 4; ++mu)
            for (unsigned alpha = 0; alpha < 4; ++alpha)
                for (unsigned beta = 0; beta < 4; ++beta)
                    c.V[column](mu) +=
                        connection.gamma(mu, alpha, beta).real * c.k(alpha) * c.X[column](beta);
}

inline std::vector<Case> BuildCases(std::vector<ReferenceAudit>* audits = nullptr) {
    std::vector<Case> cases;
    const std::array<Vector, 4> flat_positions{Vector{20, -7, 5}, Vector{-4, 18, 8},
                                               Vector{0, 0, 20}, Vector{0, 0, -20}};
    const std::array<double, 4> frequencies{1e-6, 1, 1e6, 2};
    const std::array<double, 4> seeds{1e-6, 1e-3, 0.1, 1e-3};
    for (int sign : {-1, 1}) {
        for (unsigned index = 0; index < flat_positions.size(); ++index) {
            Case c{};
            c.name = "flat_" + std::to_string(sign) + "_" + std::to_string(index);
            c.mass = 0;
            c.spin = DeviceParameter(sign * 0.9);
            c.seed = seeds[index];
            c.x = Four(0.2, flat_positions[index]);
            c.direction = Unit({0.8, 0.3, std::copysign(0.6, flat_positions[index][2])});
            c.frequency = frequencies[index];
            c.k = Four(-1, c.direction) * c.frequency;
            const auto screen = Basis(c.direction);
            const std::array<Vector, 2> dx{Vector{1, 0.2, -0.3}, Vector{-0.7, 0.4, 0.1}};
            std::array<Vec4, 2> K;
            for (unsigned column = 0; column < 2; ++column) {
                const double parity = (column == 1 && index % 2) ? -1.0 : 1.0;
                const double rate = column == 0 ? 0.7 : -0.2;
                c.X[column] = Four(0.37, dx[column]) * c.seed;
                K[column] = Four(0, screen[column]) * (parity * c.frequency * c.seed) +
                            c.k * (rate * c.seed);
                c.jacobian[column][column] = parity;
                c.frequency_derivative[column] = rate * c.frequency;
            }
            CovariantInput(c, K);
            cases.push_back(c);
        }
    }
    // Exact principal-axis first variation: dn_inf=(r+ia)/(r^2+a^2)dv0,
    // while dv0=(r-ia)dk_transverse. Thus the angular matrix is identity.
    for (double requested_spin : {-1.0, -0.9, 0.0, 0.9, 1.0}) {
        for (int pole : {-1, 1}) {
            Case c{};
            c.name = "axis_" + std::to_string(requested_spin) + "_" + std::to_string(pole);
            c.mass = DeviceParameter(1);
            c.spin = DeviceParameter(requested_spin);
            c.seed = 1e-3;
            c.direction = {0, 0, static_cast<double>(pole)};
            c.x = Four(0, {0, 0, pole * (pole < 0 ? 3.0 : 20.0)});
            c.k = Four(-1, c.direction);
            c.frequency = 1;
            c.jacobian = {{{1, 0}, {0, 1}}};
            const auto screen = Basis(c.direction);
            const std::array<Vec4, 2> K{Four(0, screen[0]) * c.seed, Four(0, screen[1]) * c.seed};
            CovariantInput(c, K);
            cases.push_back(c);
        }
    }
    const std::array<std::array<double, 3>, 8> definitions{
        std::array<double, 3>{1, -1, 8},    std::array<double, 3>{1, -0.9, 8},
        std::array<double, 3>{1, 0, 8},     std::array<double, 3>{1, 0.9, 8},
        std::array<double, 3>{1, 1, 8},     std::array<double, 3>{1, -0.9, 40},
        std::array<double, 3>{1, 0.9, 200}, std::array<double, 3>{1e-3, -0.9e-3, 8e-3}};
    for (unsigned index = 0; index < definitions.size(); ++index) {
        Case c{};
        c.name = "carter_" + std::to_string(index);
        c.mass = DeviceParameter(definitions[index][0]);
        c.spin = DeviceParameter(definitions[index][1]);
        c.seed = 1e-3;
        const double length = c.mass;
        const CarterInput input{c.mass, c.spin,     definitions[index][2], 0.9, 1.1,
                                -1,     2 * length, 0.3 * length};
        const Event event = CarterEvent(input);
        c.x = event.x;
        c.k = event.k;
        const std::array<std::array<double, 6>, 2> variation{
            std::array<double, 6>{0.7 * length, 0.03, -0.02, 0.04, 0.11 * length, -0.08 * length},
            std::array<double, 6>{-0.4 * length, -0.02, 0.04, -0.03, -0.05 * length,
                                  0.09 * length}};
        std::array<Vec4, 2> K;
        constexpr double h = 1e-4;
        for (unsigned column = 0; column < 2; ++column) {
            const auto mm = CarterEvent(Vary(input, variation[column], -2 * h));
            const auto m = CarterEvent(Vary(input, variation[column], -h));
            const auto p = CarterEvent(Vary(input, variation[column], h));
            const auto pp = CarterEvent(Vary(input, variation[column], 2 * h));
            c.X[column] = (mm.x - m.x * 8 + p.x * 8 - pp.x) * (c.seed / (12 * h));
            K[column] = (mm.k - m.k * 8 + p.k * 8 - pp.k) * (c.seed / (12 * h));
        }
        CovariantInput(c, K);
        c.direction = Reference(input, 8192);
        const auto coarse = Reference(input, 4096);
        const auto screen = Basis(c.direction);
        ReferenceAudit audit{};
        audit.name = c.name;
        for (unsigned i = 0; i < 3; ++i)
            audit.nominal_step_refinement_error =
                std::max(audit.nominal_step_refinement_error, std::abs(c.direction[i] - coarse[i]));
        for (unsigned column = 0; column < 2; ++column) {
            std::array<Vector, 2> derivatives;
            for (unsigned refinement = 0; refinement < 2; ++refinement) {
                const double dh = refinement == 0 ? 2e-4 : 1e-4;
                const auto m = Reference(Vary(input, variation[column], -dh), 8192);
                const auto p = Reference(Vary(input, variation[column], dh), 8192);
                for (unsigned i = 0; i < 3; ++i)
                    derivatives[refinement][i] = (p[i] - m[i]) / (2 * dh);
            }
            for (unsigned i = 0; i < 3; ++i)
                audit.angular_difference_refinement_error =
                    std::max(audit.angular_difference_refinement_error,
                             std::abs(derivatives[0][i] - derivatives[1][i]));
            for (unsigned row = 0; row < 2; ++row)
                c.jacobian[row][column] = Dot(screen[row], derivatives[1]);
            c.frequency_derivative[column] = -variation[column][3];
        }
        c.frequency = -input.energy;
        if (audit.nominal_step_refinement_error > 2e-12 ||
            audit.angular_difference_refinement_error > 2e-9)
            throw std::runtime_error("Independent Carter reference did not converge: " + c.name);
        if (audits) audits->push_back(audit);
        cases.push_back(c);
    }
    {
        Case c{};
        c.name = "hidden_radial_turn";
        c.mass = 1;
        c.spin = 0;
        c.seed = 1e-3;
        c.expected_status = 3;
        const CarterInput input{1, 0, 2.5, std::numbers::pi / 2, 0, -1, 5.3, 0};
        const auto event = CarterEvent(input);
        c.x = event.x;
        c.k = event.k;
        CovariantInput(c, {});
        cases.push_back(c);
    }
    {
        Case c = cases.front();
        c.name = "invalid_seed";
        c.seed = 0;
        c.expected_status = 2;
        cases.push_back(c);
    }
    {
        Case c = cases.front();
        c.name = "zero_work_budget";
        c.maximum_attempts = 0;
        c.expected_status = 4;
        cases.push_back(c);
    }
    return cases;
}
}  // namespace detail
inline std::vector<Case> Cases() { return detail::BuildCases(); }
}  // namespace sirius::test::kerr_infinity_oracle
