// The actual tracer transports screens on accepted segments. The independent
// algebraic reference reconstructs their disk angle from conserved
// Walker-Penrose constants, using the test-only Boyer-Lindquist formulation.
// References: Walker & Penrose (1970); James et al. (2015), arXiv:1502.03808.
// Committed-segment invariant: Gelles et al. (2021), PRD104,044060, Eq11-12,
// https://link.aps.org/accepted/10.1103/PhysRevD.104.044060.
#include "sirius/backend/cpu/geodesic_tracer.h"
#include "sirius/core/camera_launch.h"
#include "sirius/core/constants.h"
#include "sirius/core/metrics/kerr_schild_family.h"
#include "sirius/oracle/polarisation_transport.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <format>
#include <limits>
#include <numbers>
#include <string>

namespace {
using namespace sirius::core;
using namespace sirius::backend;
using Bl = sirius::oracle::Vec4d;

struct BlEvent {
    Bl position;
    Bl tangent;
};

// Differentiate the defining oblate-coordinate quartic and the azimuthal
// twist explicitly. No production coordinate adapter supplies the reference.
BlEvent ToBl(const Vec4& x, const Vec4& v, double a,
             TraceResult::TerminalChart chart = TraceResult::TerminalChart::MetricNative) {
    const double X = x(1), Y = x(2), Z = x(3);
    const double s = X * X + Y * Y + Z * Z - a * a;
    const double u = (s + std::sqrt(s * s + 4 * a * a * Z * Z)) / 2;
    const double r = std::sqrt(u);
    const double theta = std::acos(Z / r);
    const double ds = 2 * (X * v(1) + Y * v(2) + Z * v(3));
    const double dr = (u * ds + 2 * a * a * Z * v(3)) / ((2 * u - s) * 2 * r);
    const double dtheta = (std::cos(theta) * dr - v(3)) / (r * std::sin(theta));
    // Outgoing KS reverses the radial time/azimuth terms. Differentiate the
    // oblate definitions and dt_KS=dt_BL +/- 2r dr/Delta directly; no live map.
    const double sign = chart == TraceResult::TerminalChart::OutgoingKerrSchild ? -1 : 1;
    const double dphi = (X * v(2) - Y * v(1)) / (X * X + Y * Y) + sign * a * dr / (r * r + a * a) -
                        sign * a * dr / (r * r - 2 * r + a * a);
    return {{0, r, theta, std::atan2(Y, X) - sign * std::atan2(a, r)},
            {v(0) - sign * 2 * r * dr / (r * r - 2 * r + a * a), dr, dtheta, dphi}};
}

// Gelles et al., PRD104,044060 (2021), Eq11-12. Antisymmetric k/f products
// also make additions proportional to k gauge-invariant. Independent of the
// production transport and of oracle::WalkerPenroseConstant.
std::complex<double> BlConstant(const BlEvent& event, const Bl& f, double a) {
    const auto& k = event.tangent;
    const double r = event.position.r, theta = event.position.theta;
    const double sine = std::sin(theta);
    const double first = k.t * f.r - k.r * f.t + a * sine * sine * (k.r * f.phi - k.phi * f.r);
    const double second = ((r * r + a * a) * (k.phi * f.theta - k.theta * f.phi) -
                           a * (k.t * f.theta - k.theta * f.t)) *
                          sine;
    return std::complex<double>(first, -second) * std::complex<double>(r, -a * std::cos(theta));
}

// Integrate the independently derived chart derivatives from the launch radius.
// Relative t/phi origins are anchored at the measured launch, so neither the
// production chart's additive constants nor its coordinate helpers are reused.
Bl BlPosition(const Vec4& x, double a, double launch_radius, TraceResult::TerminalChart chart) {
    auto result = ToBl(x, Vec4{}, a, chart).position;
    const double outer = 1 + std::sqrt(1 - a * a), inner = 1 - std::sqrt(1 - a * a);
    const double outer_log = std::log((result.r - outer) / (launch_radius - outer));
    const double inner_log = std::log((result.r - inner) / (launch_radius - inner));
    const double inverse_delta = (outer_log - inner_log) / (outer - inner);
    const double time_shift = 2 * (outer * outer_log - inner * inner_log) / (outer - inner);
    const double sign = chart == TraceResult::TerminalChart::OutgoingKerrSchild ? -1 : 1;
    result.t = x(0) - sign * time_shift;
    result.phi -= sign * a * inverse_delta;
    return result;
}

class CommittedPolarisationMonitor {
  public:
    CommittedPolarisationMonitor(double spin, const CameraLaunch& launch, const CameraRay& camera,
                                 double maximum_step)
        : spin_(spin), launch_(launch), camera_(camera), maximum_step_(maximum_step) {}

    void Observe(const std::array<PolarisedRay, 2>& basis, TraceResult::TerminalChart chart) {
        if (::testing::Test::HasFatalFailure()) return;
        EXPECT_EQ(chart, TraceResult::TerminalChart::OutgoingKerrSchild);
        for (const auto& state : basis) {
            ASSERT_TRUE(std::isfinite(state.affine));
            for (int axis = 0; axis < 4; ++axis) {
                ASSERT_TRUE(std::isfinite(state.position(axis)));
                ASSERT_TRUE(std::isfinite(state.velocity(axis)));
                ASSERT_TRUE(std::isfinite(state.polarisation(axis)));
            }
        }
        EXPECT_EQ(basis[0].affine, basis[1].affine);
        if (samples_ == 0)
            EXPECT_EQ(basis[0].affine, 0);
        else
            EXPECT_GT(basis[0].affine, last_[0].affine);  // No rollback or duplicate observations.
        const auto event = ToBl(basis[0].position, basis[0].velocity, spin_, chart);
        ASSERT_GT(event.position.r, 1 + std::sqrt(1 - spin_ * spin_));
        ASSERT_GT(event.position.theta, sirius::oracle::kBoyerLindquistPoleMargin);
        ASSERT_LT(event.position.theta,
                  std::numbers::pi - sirius::oracle::kBoyerLindquistPoleMargin);
        std::array<Bl, 2> vectors;
        for (int column = 0; column < 2; ++column) {
            for (int axis = 0; axis < 4; ++axis) {
                EXPECT_EQ(basis[column].position(axis), basis[0].position(axis));
                EXPECT_EQ(basis[column].velocity(axis), basis[0].velocity(axis));
            }
            vectors[column] =
                ToBl(basis[column].position, basis[column].polarisation, spin_, chart).tangent;
            for (int axis = 0; axis < 4; ++axis) {
                ASSERT_TRUE(std::isfinite(event.tangent[axis]));
                ASSERT_TRUE(std::isfinite(vectors[column][axis]));
            }
            const auto constant = BlConstant(event, vectors[column], spin_);
            ASSERT_TRUE(std::isfinite(constant.real()));
            ASSERT_TRUE(std::isfinite(constant.imag()));
            if (samples_ == 0) {
                initial_[column] = constant;
                ASSERT_TRUE(std::isfinite(std::abs(constant)));
                ASSERT_GT(std::abs(constant), 0);
                const std::array<double, 3> rest{camera_.direction(1), camera_.direction(2),
                                                 camera_.direction(3)};
                const auto screens = relativity::ObserverScreenBasis(launch_.observer, rest);
                ASSERT_TRUE(screens);
                const auto expected =
                    BlConstant(ToBl(launch_.position, launch_.tangent, spin_),
                               ToBl(launch_.position, (*screens)[column], spin_).tangent, spin_);
                ASSERT_TRUE(std::isfinite(expected.real()) && std::isfinite(expected.imag()));
                ASSERT_TRUE(std::isfinite(std::abs(expected)));
                ASSERT_GT(std::abs(expected), 0);
                const double launch_defect = std::abs(constant - expected) / std::abs(expected);
                ASSERT_TRUE(std::isfinite(launch_defect));
                EXPECT_LT(launch_defect, constants::geodesic::kConservationTol);
            }
            const double relative_drift =
                std::abs(constant - initial_[column]) / std::abs(initial_[column]);
            ASSERT_TRUE(std::isfinite(relative_drift));
            maximum_[column] = std::max(maximum_[column], relative_drift);
        }
        sirius::oracle::KerrMetricD metric(1, spin_);
        double g[4][4], inverse[4][4];
        metric.Evaluate(event.position, g, inverse);
        const auto dot = [&](const Bl& first, const Bl& second) {
            return sirius::oracle::InnerProductD(g, first, second);
        };
        const double energy = std::abs(-(g[0][0] * event.tangent.t + g[0][3] * event.tangent.phi));
        ASSERT_TRUE(std::isfinite(energy));
        ASSERT_GT(energy, 0);
        for (const auto& f : vectors) {
            const double norm_defect = std::abs(dot(f, f) - 1);
            const double tangent_defect = std::abs(dot(f, event.tangent)) / energy;
            ASSERT_TRUE(std::isfinite(norm_defect) && std::isfinite(tangent_defect));
            maximum_screen_ = std::max({maximum_screen_, norm_defect, tangent_defect});
        }
        const double cross_defect = std::abs(dot(vectors[0], vectors[1]));
        ASSERT_TRUE(std::isfinite(cross_defect));
        maximum_screen_ = std::max(maximum_screen_, cross_defect);
        if (samples_ == 0) {
            launch_radius_ = event.position.r;
            const auto native = BlPosition(launch_.position, spin_, launch_radius_,
                                           TraceResult::TerminalChart::MetricNative);
            const auto current = BlPosition(basis[0].position, spin_, launch_radius_, chart);
            time_origin_ = native.t - current.t;
            angle_origin_ = native.phi - current.phi;
        }
        last_ = basis;
        last_chart_ = chart;
        ++samples_;
    }

    void Check(const TraceResult& result, const std::string& prefix) const {
        ASSERT_GT(samples_, 1U) << "outgoing worker or committed observations are absent";
        EXPECT_LE(samples_ - 1, static_cast<unsigned>(result.steps_taken));
        // Every accepted interval is bounded by this unchanged public cap.
        // An initial/final-only hook cannot qualify a long transport path.
        EXPECT_GE(static_cast<double>(samples_ - 1), result.affine_length / maximum_step_);
        EXPECT_EQ(last_[0].affine, result.affine_length) << "final clipped frame was not observed";
        ASSERT_TRUE(result.final_tangent);
        const auto observed = ToBl(last_[0].position, last_[0].velocity, spin_, last_chart_);
        const auto expected =
            ToBl(result.final_position, *result.final_tangent, spin_, result.terminal_chart);
        auto position = BlPosition(last_[0].position, spin_, launch_radius_, last_chart_);
        position.t += time_origin_;
        position.phi += angle_origin_;
        const auto terminal =
            BlPosition(result.final_position, spin_, launch_radius_, result.terminal_chart);
        for (int axis = 0; axis < 4; ++axis) {
            ASSERT_TRUE(std::isfinite(position[axis]) && std::isfinite(terminal[axis]));
            ASSERT_TRUE(std::isfinite(observed.tangent[axis]) &&
                        std::isfinite(expected.tangent[axis]));
            const double difference =
                axis == 3 ? std::remainder(position[axis] - terminal[axis], 2 * std::numbers::pi)
                          : position[axis] - terminal[axis];
            EXPECT_LT(std::abs(difference) / std::max(1., std::abs(terminal[axis])),
                      constants::geodesic::kConservationTol);
            EXPECT_LT(std::abs(observed.tangent[axis] - expected.tangent[axis]) /
                          std::max(1., std::abs(expected.tangent[axis])),
                      constants::geodesic::kConservationTol);
        }
        for (int column = 0; column < 2; ++column) {
            EXPECT_LT(maximum_[column], constants::geodesic::kConservationTol);
            ::testing::Test::RecordProperty(
                prefix + "_wp_max_relative_drift_" + std::to_string(column),
                std::format("{:.17g}", maximum_[column]));
        }
        EXPECT_LT(maximum_screen_, constants::geodesic::kConservationTol);
        ::testing::Test::RecordProperty(prefix + "_committed_samples", std::to_string(samples_));
        ::testing::Test::RecordProperty(prefix + "_screen_max_defect",
                                        std::format("{:.17g}", maximum_screen_));
        ::testing::Test::RecordProperty(prefix + "_terminal_affine",
                                        std::format("{:.17g}", last_[0].affine));
    }

  private:
    double spin_;
    const CameraLaunch& launch_;
    const CameraRay& camera_;
    double maximum_step_;
    unsigned samples_ = 0;
    std::array<std::complex<double>, 2> initial_{};
    std::array<double, 2> maximum_{};
    double maximum_screen_ = 0, launch_radius_ = 0, time_origin_ = 0, angle_origin_ = 0;
    std::array<PolarisedRay, 2> last_{};
    TraceResult::TerminalChart last_chart_ = TraceResult::TerminalChart::MetricNative;
};

struct DiskReference {
    double angle, degree, intensity, redshift;
};

DiskReference DiskAngle(const CameraLaunch& launch, const CameraRay& camera,
                        const TraceResult& result, double a) {
    const std::array<double, 3> rest{camera.direction(1), camera.direction(2), camera.direction(3)};
    const auto screens = relativity::ObserverScreenBasis(launch.observer, rest);
    const auto initial = ToBl(launch.position, launch.tangent, a);
    std::array<std::complex<double>, 2> constants;
    for (int column = 0; column < 2; ++column) {
        const auto f = ToBl(launch.position, (*screens)[column], a).tangent;
        constants[column] = sirius::oracle::WalkerPenroseConstant(
            sirius::oracle::PolarisedStateD(initial.position, initial.tangent, f), a);
    }
    const auto terminal = ToBl(result.final_position, *result.final_tangent, a);
    sirius::oracle::KerrMetricD metric(1, a);
    double g[4][4], inv[4][4];
    metric.Evaluate(terminal.position, g, inv);
    const auto dot = [&](const Bl& lhs, const Bl& rhs) {
        return sirius::oracle::InnerProductD(g, lhs, rhs);
    };
    // Circular emitter, directly normalized with the independent BL metric.
    const double omega = 1 / (std::pow(terminal.position.r, 1.5) + a);
    Bl emitter(1, 0, 0, omega);
    emitter = emitter * (1 / std::sqrt(-dot(emitter, emitter)));
    const double frequency = dot(terminal.tangent, emitter);
    Bl n = terminal.tangent * (-1 / frequency) - emitter;
    Bl normal(0, 0, -1 / std::sqrt(g[2][2]), 0);
    const double mu = std::abs(dot(n, normal));
    Bl meridian = normal - n * dot(normal, n);
    meridian = meridian * (1 / std::sqrt(dot(meridian, meridian)));
    const auto endpoint = sirius::oracle::WalkerPenroseConstant(
        sirius::oracle::PolarisedStateD(terminal.position, terminal.tangent, meridian), a);
    const double det =
        constants[0].real() * constants[1].imag() - constants[1].real() * constants[0].imag();
    const double first =
        (endpoint.real() * constants[1].imag() - constants[1].real() * endpoint.imag()) / det;
    const double second =
        (constants[0].real() * endpoint.imag() - endpoint.real() * constants[0].imag()) / det;
    return {std::remainder(std::atan2(second, first) + std::numbers::pi / 2, std::numbers::pi),
            .1171 * (1 - mu) / (1 + 3.582 * mu), (1 + 2.06 * mu) / (1 + (2.0 / 3) * 2.06),
            1 / frequency};
}

TEST(PolarisationAcceptance, NormalEmissionHasFiniteUnpolarisedStokes) {
    KerrSchildFamily metric(KerrSchildParams::Schwarzschild(1));
    TracerConfig config;
    config.enable_disk = true;
    config.enable_polarisation = true;
    config.disk_inner = 6;
    config.disk_outer = 20;
    for (double tilt : {0., 5e-4, 1e-3}) {
        CameraRay ray;
        ray.origin(1) = 10;
        ray.origin(2) = std::numbers::pi / 2;
        ray.direction(1) = tilt;
        ray.direction(2) = std::sqrt(1 - tilt * tilt);
        ray.beta_forward = -.2;
        // Match the represented circular emitter, including its fp32 omega.
        ray.beta_right =
            static_cast<double>(static_cast<float>(std::pow(10., -1.5))) * 10 * std::sqrt(1.2);
        GeodesicTracer tracer(&metric, config);
        const auto result = tracer.Trace(ray);
        ASSERT_EQ(result.outcome, TraceResult::Outcome::DiskHit);
        ASSERT_FALSE(result.numerical_failure);
        ASSERT_EQ(result.num_disk_crossings, 1);
        const auto& crossing = result.disk_crossings[0];
        ASSERT_TRUE(crossing.polarisation_valid);
        EXPECT_TRUE(std::isfinite(crossing.polarisation_evpa));
        // The atmosphere law gives p=O(tilt^2); EVPA is unobservable at p=0.
        EXPECT_LE(crossing.polarisation_degree, tilt * tilt * .1);
        if (tilt == 0) {
            EXPECT_FLOAT_EQ(crossing.polarisation_degree, 0);
            EXPECT_FLOAT_EQ(crossing.polarisation_evpa, 0);
        }
    }
}

TEST(PolarisationAcceptance, ActualDiskStokesMatchWalkerPenroseReconstruction) {
    // Existing launch tests separately judge the measured input frame. This
    // comparison covers real accepted-segment transport and disk extraction;
    // independent full-path trajectory tests judge the central event itself.
    for (double a : {0., .7, .998}) {
        KerrSchildFamily metric(KerrSchildParams::Kerr(1, a));
        CameraRay ray;
        ray.origin(1) = 50;
        ray.origin(2) = std::numbers::pi / 3;
        ray.direction(1) = -std::sqrt(1 - .2 * .2 - .08 * .08);
        ray.direction(2) = .2;
        ray.direction(3) = .08;
        const auto launch = LaunchCameraRay(metric, a, ray);
        ASSERT_TRUE(launch);
        double previous_error = 1;
        double first_error = 0;
        for (int refinement = 0; refinement < 3; ++refinement) {
            SCOPED_TRACE(std::format("a={} refinement={}", a, refinement));
            TracerConfig config;
            config.enable_disk = true;
            config.enable_polarisation = true;
            config.disk_inner = 6;
            config.disk_outer = 20;
            config.integrator.initial_step = .02f;
            config.integrator.max_step = .4f / static_cast<float>(1 << refinement);
            config.integrator.abs_tolerance = 1e-7f / static_cast<float>(1 << (3 * refinement));
            config.integrator.rel_tolerance = config.integrator.abs_tolerance;
            GeodesicTracer tracer(&metric, config);
            CommittedPolarisationMonitor monitor(a, *launch, ray, config.integrator.max_step);
            tracer.SetPolarisationObserver(
                [&](const auto& basis, auto chart) { monitor.Observe(basis, chart); });
            const auto result = tracer.Trace(ray);
            ASSERT_EQ(result.outcome, TraceResult::Outcome::DiskHit);
            ASSERT_FALSE(result.numerical_failure);
            ASSERT_TRUE(result.final_tangent);
            ASSERT_NO_FATAL_FAILURE(
                monitor.Check(result, std::format("spin_{:.3f}_r{}", a, refinement)));
            const auto reference = DiskAngle(*launch, ray, result, a);
            const auto& crossing = result.disk_crossings[0];
            ASSERT_TRUE(crossing.polarisation_valid);
            const double error = std::abs(
                std::remainder(crossing.polarisation_evpa - reference.angle, std::numbers::pi));
            // Same 1e-4 live transport budget as the existing WP gate. The
            // thermal angle/degree public fields have binary32 precision.
            EXPECT_LT(error, 1e-4);
            if (refinement == 0) first_error = error;
            EXPECT_LE(error, previous_error + 2e-7);
            EXPECT_NEAR(crossing.polarisation_degree, reference.degree, 2e-7);
            EXPECT_NEAR(crossing.polarisation_intensity_scale, reference.intensity, 2e-6);
            EXPECT_NEAR(crossing.redshift, reference.redshift, 2e-6);
            for (bool sine : {false, true}) {
                const double observed = crossing.polarisation_degree *
                                        (sine ? std::sin(2 * crossing.polarisation_evpa)
                                              : std::cos(2 * crossing.polarisation_evpa));
                const double expected = reference.degree * (sine ? std::sin(2 * reference.angle)
                                                                 : std::cos(2 * reference.angle));
                EXPECT_NEAR(observed, expected, 2e-5);
            }
            const auto prefix = std::format("spin_{:.3f}_r{}", a, refinement);
            RecordProperty(prefix + "_evpa_error_rad", std::format("{:.17g}", error));
            RecordProperty(prefix + "_disk_radius", std::format("{:.17g}", crossing.r));
            previous_error = error;
        }
        // An error already below the public binary32 angle's rounding scale
        // cannot exhibit a reliable convergence slope. Above that scale,
        // refinement must reduce the independently measured error.
        const double angle_rounding = std::numeric_limits<float>::epsilon() * std::numbers::pi / 2;
        EXPECT_LE(previous_error, std::max(first_error / 2, angle_rounding));
        RecordProperty(std::format("spin_{:.3f}_refinement_behavior", a),
                       first_error <= angle_rounding ? "binary32_angle_plateau" : "error_reduced");
    }
}

TEST(PolarisationAcceptance, CommittedOutwardTransportReachesPhysicalEscapeAndObserverIsPassive) {
    constexpr double spin = .998;
    KerrSchildFamily metric(KerrSchildParams::Kerr(1, spin));
    CameraRay ray;
    ray.origin(1) = 50;
    ray.origin(2) = std::numbers::pi / 3;
    ray.direction(1) = std::sqrt(1 - .2 * .2 - .08 * .08);
    ray.direction(2) = .2;
    ray.direction(3) = .08;
    const auto launch = LaunchCameraRay(metric, spin, ray);
    ASSERT_TRUE(launch);
    TracerConfig config;
    config.enable_polarisation = true;
    config.disk_inner = 6;
    config.disk_outer = 20;
    config.escape_radius = 80;  // Declared finite sphere enclosing this exterior observer.
    config.integrator.initial_step = .02f;
    config.integrator.max_step = .4f;
    config.integrator.abs_tolerance = config.integrator.rel_tolerance = 1e-7f;
    GeodesicTracer tracer(&metric, config);
    CommittedPolarisationMonitor monitor(spin, *launch, ray, config.integrator.max_step);
    tracer.SetPolarisationObserver(
        [&](const auto& basis, auto chart) { monitor.Observe(basis, chart); });
    const auto observed = tracer.Trace(ray);
    ASSERT_FALSE(observed.numerical_failure);
    ASSERT_EQ(observed.outcome, TraceResult::Outcome::Escaped);
    EXPECT_EQ(observed.num_disk_crossings, 0);
    EXPECT_LT(std::abs(std::hypot(observed.final_position(1), observed.final_position(2),
                                  observed.final_position(3)) -
                       config.escape_radius) /
                  config.escape_radius,
              constants::geodesic::kConservationTol);
    ASSERT_NO_FATAL_FAILURE(monitor.Check(observed, "outward"));
    tracer.SetPolarisationObserver({});
    const auto passive = tracer.Trace(ray);
    ASSERT_FALSE(passive.numerical_failure);
    EXPECT_EQ(passive.outcome, observed.outcome);
    EXPECT_EQ(passive.steps_taken, observed.steps_taken);
    EXPECT_EQ(passive.affine_length, observed.affine_length);
    ASSERT_TRUE(passive.final_tangent && observed.final_tangent);
    for (int axis = 0; axis < 4; ++axis) {
        EXPECT_EQ(passive.final_position(axis), observed.final_position(axis));
        EXPECT_EQ((*passive.final_tangent)(axis), (*observed.final_tangent)(axis));
    }
    config.enable_polarisation = false;
    GeodesicTracer unpolarised(&metric, config);
    unsigned unexpected_samples = 0;
    unpolarised.SetPolarisationObserver([&](const auto&, auto) { ++unexpected_samples; });
    const auto without_transport = unpolarised.Trace(ray);
    ASSERT_FALSE(without_transport.numerical_failure);
    EXPECT_EQ(without_transport.outcome, observed.outcome);
    EXPECT_EQ(unexpected_samples, 0U);
}
}  // namespace
