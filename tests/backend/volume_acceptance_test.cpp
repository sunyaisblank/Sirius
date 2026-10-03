// Numerical acceptance of the represented relative grey disk, without a
// renderer. A Schwarzschild radial null ray reduces the reference transport to
// a one-dimensional ODE in radius; no production metric, source, Planck,
// opacity, interpolation or transfer helper supplies its expected values.
#include "sirius/backend/cpu/geodesic_tracer.h"
#include "sirius/core/camera_launch.h"
#include "sirius/core/metrics/kerr_schild_family.h"

#ifdef SIRIUS_HAS_RETAINED_COMPUTE
#include "sirius/backend/retained_trace_executor.h"
#endif

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <format>
#include <numbers>
#include <string>
#include <vector>

namespace {
using sirius::backend::DiskTemperatureModel;
using sirius::backend::GeodesicTracer;
using sirius::backend::TracerConfig;
using sirius::backend::TraceResult;
using sirius::core::CameraRay;
using sirius::core::KerrSchildFamily;
using sirius::core::KerrSchildParams;
using sirius::core::color_modes::Mode;

constexpr double kObserverRadius = 30;
constexpr double kInner = 6, kOuter = 20, kHOverR = .1;
constexpr double kTemperatureKelvin = 30000;
// These are separate controls: the integrator limits accepted affine segment
// length, while volume samples set midpoint quadrature inside each segment.
// Geodesic error tolerances do not control a disk-support crossing or source
// quadrature. Qualify this explicit finest profile to 0.05%, not every profile.
constexpr std::array<float, 3> kMaximumSteps{.5f, .25f, .125f};
constexpr std::array<int, 3> kVolumeSamples{8, 32, 128};
using Channels = std::array<double, 3>;
using State = std::array<double, 4>;  // tau, relative observed RGB.

CameraRay RadialRay() {
    CameraRay ray;
    ray.origin(1) = kObserverRadius;
    ray.origin(2) = std::numbers::pi / 2;
    ray.direction(1) = -1;
    return ray;
}

double PrincipalEnergy(double spin) {
    const double r = kObserverRadius;
    return (r * r - 2 * r + spin * spin) * std::sqrt(1 + 2 / r) / (r * r + 2 * r + spin * spin);
}

CameraRay PrincipalRay(double spin) {
    auto ray = RadialRay();
    // At launch phi_KS=0, (x,y)=(r,a), the Kerr-Schild null one-form has
    // spatial direction (1,0,0). Its metric has g_tx=2/r, g_xx=1+2/r.
    // An outgoing principal photon has L=aE, Q=0 and the following tangent.
    // Project its past tangent onto the independently normalized radial and
    // azimuthal Eulerian axes to specify the public camera rest direction.
    const double r = kObserverRadius, a = spin, q = 1 + 2 / r;
    const double delta = r * r - 2 * r + a * a, d = r * r + 2 * r + a * a;
    const double e = PrincipalEnergy(a);
    const double t = -e * d / delta, x = -e * (1 - 2 * a * a / delta);
    const double y = -2 * a * r * e / delta;
    const double px = 2 / r * t + q * x;
    ray.direction(1) = (r * px + a * y) / std::sqrt(d);
    ray.direction(3) = (-a * px + (r + 2) * y) / std::sqrt(d * q);
    return ray;
}

void CheckPrincipalLaunch(KerrSchildFamily& metric, const CameraRay& ray, double spin) {
    const auto launch = sirius::core::LaunchCameraRay(metric, spin, ray);
    ASSERT_TRUE(launch);
    const double r = kObserverRadius, q = 1 + 2 / r;
    const auto& k = launch->tangent;
    // Independent equatorial metric contractions verify the actual camera
    // family; measured production constants never set the reference inputs.
    const double energy = (q - 2) * k(0) + 2 / r * k(1);
    const double px = -(2 / r * k(0) + q * k(1)), py = -k(2);
    const double angular_momentum = r * py - spin * px;
    EXPECT_NEAR(energy, PrincipalEnergy(spin), 3e-13);
    EXPECT_NEAR(angular_momentum, spin * energy, 3e-13);
    EXPECT_NEAR(k(1) + spin / r * k(2), -energy, 3e-13);
    EXPECT_NEAR(k(3), 0, 3e-13);
    EXPECT_NEAR(ray.direction(1) * ray.direction(1) + ray.direction(3) * ray.direction(3), 1,
                3e-13);
}

TracerConfig VolumeConfig(
    float maximum_step, int samples,
    DiskTemperatureModel temperature_model = DiskTemperatureModel::ShakuraSunyaev) {
    TracerConfig config;
    config.escape_radius = 80;
    config.horizon_factor = 1;
    config.max_steps = 10000;
    config.enable_disk = config.enable_volumetric = true;
    config.disk_inner = kInner;
    config.disk_outer = kOuter;
    config.disk_temperature_model = temperature_model;
    config.disk_temperature_scale_kelvin = kTemperatureKelvin;
    config.volumetric_scale_height_ratio = kHOverR;
    config.volumetric_flare_power = 0;
    config.volumetric_tau_midplane = .2f;
    config.volumetric_tau_max = 10;
    config.volumetric_samples = samples;
    config.integrator.initial_step = maximum_step;
    config.integrator.max_step = maximum_step;
    config.integrator.min_step = 1e-6f;
    config.integrator.abs_tolerance = config.integrator.rel_tolerance = 1e-9f;
    return config;
}

// The ingoing Kerr-Schild Eulerian camera is not a static Schwarzschild
// observer. For the physical outward radial photon, E=-p_t and unit launch
// frequency imply E=(1-2/r0)/sqrt(1+2/r0). Along the past ray dr/dlambda=-E;
// circular matter has u^t=1/sqrt(1-3/r), so dl_comoving=u^t |dr| and
// g=nu_camera/nu_emitter=1/(E u^t). This fixes signs and affine normalisation
// independently of the production camera and outgoing chart.
double Frequency(double radius) {
    const double energy = (1 - 2 / kObserverRadius) / std::sqrt(1 + 2 / kObserverRadius);
    return std::sqrt(1 - 3 / radius) / energy;
}

double FluxShape(double radius, DiskTemperatureModel model) {
    if (radius <= kInner) return 0;
    if (model == DiskTemperatureModel::ShakuraSunyaev)
        return std::pow(kInner / radius, 3) * (1 - std::sqrt(kInner / radius));
    // Page & Thorne (1974), Eqs. 11b and 15n, specialised analytically to a=0:
    // (E-Omega L)L'=(r-6)/(2 sqrt(r)(r-3)). Substitution t=sqrt(r)
    // integrates this as t-sqrt(3)/2 log[(t-sqrt(3))/(t+sqrt(3))].
    // The explicit 0.1% zero-torque edge buffer belongs to the current model.
    if (radius <= kInner * 1.001) return 0;
    const auto primitive = [](double r) {
        const double t = std::sqrt(r), a = std::sqrt(3.);
        return t - a / 2 * std::log((t - a) / (t + a));
    };
    return 1.5 * (primitive(radius) - primitive(kInner)) /
           (std::pow(radius, 3.5) * (1 - 3 / radius));
}

double Temperature(double radius, const TracerConfig& config) {
    return config.disk_temperature_inner *
           std::pow(FluxShape(radius, config.disk_temperature_model) /
                        FluxShape(1.5 * kInner, config.disk_temperature_model),
                    .25);
}

Channels IndependentColour(double temperature) {
    // Same declared 32-midpoint approximate visible response, evaluated from
    // the defining SI Planck constants in long double. These outputs are
    // brightest-channel-normalised display colour, not SI visible radiance.
    std::array<long double, 3> xyz{};
    constexpr long double h = 6.62607015e-34L, c = 299792458.L, k = 1.380649e-23L;
    constexpr long double step = 400e-9L / 32;
    for (int i = 0; i < 32; ++i) {
        const long double wavelength = 380e-9L + (i + .5L) * step;
        const long double nm = wavelength * 1e9L;
        const auto gaussian = [nm](long double centre, long double left, long double right) {
            const long double t = (nm - centre) * (nm < centre ? left : right);
            return std::exp(-t * t / 2);
        };
        const std::array<long double, 3> matching{
            .362L * gaussian(442, .0624L, .0374L) + 1.056L * gaussian(599.8L, .0264L, .0323L) -
                .065L * gaussian(501.1L, .0490L, .0382L),
            .821L * gaussian(568.8L, .0213L, .0247L) + .286L * gaussian(530.9L, .0613L, .0322L),
            1.217L * gaussian(437, .0845L, .0278L) + .681L * gaussian(459, .0385L, .0725L)};
        const long double radiance =
            2 * h * c * c /
            (std::pow(wavelength, 5) * std::expm1(h * c / (wavelength * k * temperature)));
        for (int channel = 0; channel < 3; ++channel)
            xyz[channel] += matching[channel] * radiance * step;
    }
    const std::array<long double, 3> rgb{
        std::max(0.L, 12831.L / 3959 * xyz[0] - 329.L / 214 * xyz[1] - 1974.L / 3959 * xyz[2]),
        std::max(0.L, -851781.L / 878810 * xyz[0] + 1648619.L / 878810 * xyz[1] +
                          36519.L / 878810 * xyz[2]),
        std::max(0.L, 705.L / 12673 * xyz[0] - 2585.L / 12673 * xyz[1] + 705.L / 667 * xyz[2])};
    const long double normalizer = std::max({rgb[0], rgb[1], rgb[2], .001L});
    return {double(rgb[0] / normalizer), double(rgb[1] / normalizer), double(rgb[2] / normalizer)};
}

Channels RelativeSource(double t, double g, const TracerConfig& config) {
    const double emitted = std::pow(t, 4);
    if (t == 0) return {};
    if (config.color_mode == Mode::RedshiftMap) {
        const double diagnostic_g = std::clamp(g, .1, 3.);
        if (diagnostic_g < 1)
            return {emitted, emitted * (.1 + .9 * diagnostic_g),
                    emitted * (.05 + .95 * diagnostic_g)};
        return {emitted * (1 - .45 * (diagnostic_g - 1)), emitted * (1 - .25 * (diagnostic_g - 1)),
                emitted};
    }
    auto colour = IndependentColour(t * g * config.disk_temperature_scale_kelvin);
    for (auto& channel : colour) channel *= emitted * std::pow(g, 4);
    return colour;
}

Channels Source(double radius, const TracerConfig& config) {
    return RelativeSource(Temperature(radius, config), Frequency(radius), config);
}

long double CircularUt(long double r, long double a) {
    const long double r32 = r * std::sqrt(r);
    return (1 + a / r32) / std::sqrt(1 - 3 / r + 2 * a / r32);
}

long double PageThorneIntegrand(long double r, long double a, long double derivative_step) {
    // Differentiate the defining circular-orbit L with a five-point stencil;
    // production uses a quotient-rule derivative and Gauss-Legendre integral.
    const auto angular_momentum = [a](long double radius) {
        const long double root = std::sqrt(radius), r32 = radius * root;
        return (radius * radius - 2 * a * root + a * a) /
               (r32 * std::sqrt(1 - 3 / radius + 2 * a / r32));
    };
    const long double h = r * derivative_step;
    const long double derivative = (angular_momentum(r - 2 * h) - 8 * angular_momentum(r - h) +
                                    8 * angular_momentum(r + h) - angular_momentum(r + 2 * h)) /
                                   (12 * h);
    return derivative / CircularUt(r, a);  // E-Omega L=1/u^t.
}

template <class Function>
long double SimpsonIntegral(const Function& function, long double lower, long double upper,
                            int panels) {
    const long double step = (upper - lower) / panels;
    long double value = function(lower) + function(upper);
    for (int i = 1; i < panels; ++i) value += (i % 2 == 0 ? 2 : 4) * function(lower + i * step);
    return value * step / 3;
}

long double KerrFlux(long double r, long double a, long double integral) {
    const long double denominator = r * std::sqrt(r) + a;
    // Page-Thorne F~[-Omega,r/(r*(E-Omega L)^2)] integral.
    return 1.5L * std::sqrt(r) * CircularUt(r, a) * CircularUt(r, a) * integral /
           (r * denominator * denominator);
}

State KerrReference(const TracerConfig& config, int panels, double spin) {
    // The equatorial principal branch has L=aE,Q=0 and dr/dlambda=-E.
    // Circular matter measures nu=E*u^t*(1-a*Omega). Thus its comoving
    // length per |dr| is u^t*(1-a*Omega), and observed g is its inverse /E.
    // The model uses spheroidal BL r and Cartesian z; at z=0 its finite
    // Gaussian normalization and H=.1r are the same declared radial law.
    const double column = std::sqrt(2 * std::numbers::pi) * std::erf(3 / std::sqrt(2.));
    const auto path_per_radius = [spin](double radius) {
        const double omega = 1 / (radius * std::sqrt(radius) + spin);
        return double(CircularUt(radius, spin)) * (1 - spin * omega);
    };
    const auto rate = [&](double radius) {
        return config.volumetric_tau_midplane * std::pow(radius / kInner, -1.5) /
               (column * config.volumetric_scale_height_ratio * radius) * path_per_radius(radius);
    };
    const bool page_thorne = config.disk_temperature_model == DiskTemperatureModel::NovikovThorne;
    const double lower = page_thorne ? kInner * 1.001 : kInner;
    const long double derivative_step = .1L / panels;
    const auto integrand = [=](long double radius) {
        return PageThorneIntegrand(radius, spin, derivative_step);
    };
    const long double normalizer =
        page_thorne ? KerrFlux(1.5L * kInner, spin,
                               SimpsonIntegral(integrand, kInner, 1.5L * kInner, panels / 8))
                    : FluxShape(1.5 * kInner, DiskTemperatureModel::ShakuraSunyaev);
    const double half_step = (kOuter - lower) / (2 * panels);
    long double integral =
        page_thorne ? SimpsonIntegral(integrand, kInner, lower, panels / 128) : 0;
    std::vector<State> coefficients(2 * panels + 1);
    for (int j = 0; j <= 2 * panels; ++j) {
        const double radius = lower + j * half_step;
        if (page_thorne && j > 0) {
            const double previous = radius - half_step;
            integral += half_step *
                        (integrand(previous) + 4 * integrand(previous + half_step / 2) +
                         integrand(radius)) /
                        6;
        }
        const long double flux = page_thorne
                                     ? KerrFlux(radius, spin, integral)
                                     : FluxShape(radius, DiskTemperatureModel::ShakuraSunyaev);
        const double temperature =
            config.disk_temperature_inner * std::pow(double(flux / normalizer), .25);
        const double g = 1 / (PrincipalEnergy(spin) * path_per_radius(radius));
        const auto source = RelativeSource(temperature, g, config);
        const double opacity_rate = rate(radius);
        coefficients[j] = {opacity_rate, source[0] * opacity_rate, source[1] * opacity_rate,
                           source[2] * opacity_rate};
    }
    const auto rhs = [&](int index, const State& state) {
        const auto& c = coefficients[index];
        const double attenuation = std::exp(-state[0]);
        return State{c[0], attenuation * c[1], attenuation * c[2], attenuation * c[3]};
    };
    const auto add = [](const State& a, const State& b, double weight) {
        State result{};
        for (int j = 0; j < 4; ++j) result[j] = a[j] + weight * b[j];
        return result;
    };
    State state{};
    const double step = 2 * half_step;
    for (int i = 0; i < panels; ++i) {
        const int outer_index = 2 * (panels - i);
        const auto first = rhs(outer_index, state);
        const auto second = rhs(outer_index - 1, add(state, first, step / 2));
        const auto third = rhs(outer_index - 1, add(state, second, step / 2));
        const auto fourth = rhs(outer_index - 2, add(state, third, step));
        for (int j = 0; j < 4; ++j)
            state[j] += step * (first[j] + 2 * second[j] + 2 * third[j] + fourth[j]) / 6;
    }
    // Locate the declared source-buffer discontinuity exactly in this oracle.
    // Its outer one-sided source was integrated above; extinction continues
    // through the dark inner material while accumulated observer emission stays.
    if (page_thorne) state[0] += double(SimpsonIntegral(rate, kInner, lower, panels / 128));
    return state;
}

State Reference(const TracerConfig& config, int panels) {
    // Integrate dI/dr=-exp(-tau) S kappa u^t from the observer-nearest
    // outer edge to the inner edge, by RK4 on the reduced transfer ODE. This
    // uses neither accepted geodesic segments nor the production homogeneous
    // layer recurrence. The unmodulated, equatorial reference stays below its
    // tau budget; cap behaviour is tested separately.
    const double column_integral = std::sqrt(2 * std::numbers::pi) * std::erf(3 / std::sqrt(2.));
    const auto rhs = [&](double radius, const State& state) {
        const double opacity = config.volumetric_tau_midplane * std::pow(radius / kInner, -1.5) /
                               (column_integral * config.volumetric_scale_height_ratio * radius);
        const double rate = opacity / std::sqrt(1 - 3 / radius);
        const auto source = Source(radius, config);
        return State{rate, std::exp(-state[0]) * source[0] * rate,
                     std::exp(-state[0]) * source[1] * rate,
                     std::exp(-state[0]) * source[2] * rate};
    };
    const auto add = [](const State& a, const State& b, double weight) {
        State result{};
        for (int i = 0; i < 4; ++i) result[i] = a[i] + weight * b[i];
        return result;
    };
    const double step = (kOuter - kInner) / panels;
    State state{};
    for (int i = 0; i < panels; ++i) {
        const double r = kOuter - i * step;
        const auto first = rhs(r, state);
        const auto second = rhs(r - step / 2, add(state, first, step / 2));
        const auto third = rhs(r - step / 2, add(state, second, step / 2));
        const auto fourth = rhs(r - step, add(state, third, step));
        for (int j = 0; j < 4; ++j)
            state[j] += step * (first[j] + 2 * second[j] + 2 * third[j] + fourth[j]) / 6;
    }
    return state;
}

double RelativeError(double actual, double expected) { return std::abs(actual / expected - 1); }

void CheckReferenceRefinement(const State& fine, const State& coarse) {
    for (int i = 0; i < 4; ++i) {
        ASSERT_GT(fine[i], 0);
        EXPECT_LT(RelativeError(coarse[i], fine[i]), 1e-6);
    }
}

double CheckTrace(const TraceResult& trace, const State& reference) {
    EXPECT_FALSE(trace.numerical_failure);
    EXPECT_EQ(trace.outcome, TraceResult::Outcome::Horizon);
    EXPECT_TRUE(trace.volumetric_hit);
    EXPECT_TRUE(std::isfinite(trace.optical_depth));
    double error = RelativeError(trace.optical_depth, reference[0]);
    for (int c = 0; c < 3; ++c) {
        EXPECT_TRUE(std::isfinite(trace.volumetric_emission[c]));
        const double channel_error = RelativeError(trace.volumetric_emission[c], reference[c + 1]);
        EXPECT_TRUE(std::isfinite(channel_error));
        error = std::max(error, channel_error);
    }
    return error;
}

void RecordTraceMetrics(const TraceResult& trace, const State& reference,
                        const std::string& prefix) {
    const State actual{trace.optical_depth, trace.volumetric_emission[0],
                       trace.volumetric_emission[1], trace.volumetric_emission[2]};
    constexpr std::array<const char*, 4> names{"tau", "rgb_0", "rgb_1", "rgb_2"};
    for (int channel = 0; channel < 4; ++channel) {
        ::testing::Test::RecordProperty(std::format("{}_{}", prefix, names[channel]),
                                        std::format("{:.17g}", actual[channel]));
        ::testing::Test::RecordProperty(
            std::format("{}_{}_relative_error", prefix, names[channel]),
            std::format("{:.17g}", RelativeError(actual[channel], reference[channel])));
    }
}

template <class Trace>
void CheckConvergence(const Trace& trace, DiskTemperatureModel model, const char* prefix,
                      double spin = 0) {
    const auto config = VolumeConfig(kMaximumSteps.back(), kVolumeSamples.back(), model);
    const auto reference = spin == 0 ? Reference(config, 8192) : KerrReference(config, 8192, spin);
    const auto coarse_reference =
        spin == 0 ? Reference(config, 4096) : KerrReference(config, 4096, spin);
    CheckReferenceRefinement(reference, coarse_reference);
    double reference_error = 0;
    for (int channel = 0; channel < 4; ++channel)
        reference_error =
            std::max(reference_error, RelativeError(coarse_reference[channel], reference[channel]));
    ::testing::Test::RecordProperty(std::format("{}_reference_refinement_relative_error", prefix),
                                    std::format("{:.17g}", reference_error));
    std::array<double, 3> errors{};
    for (int level = 0; level < 3; ++level) {
        const auto settings = VolumeConfig(kMaximumSteps[level], kVolumeSamples[level], model);
        const auto result = trace(settings);
        errors[level] = CheckTrace(result, reference);
        const auto level_prefix = std::format("{}_level_{}", prefix, level);
        RecordTraceMetrics(result, reference, level_prefix);
        ::testing::Test::RecordProperty(std::format("{}_level_{}_relative_error", prefix, level),
                                        std::format("{:.17g}", errors[level]));
        ::testing::Test::RecordProperty(std::format("{}_level_{}_attempts", prefix, level),
                                        result.steps_taken);
        ::testing::Test::RecordProperty(std::format("{}_max_affine_step", level_prefix),
                                        std::format("{:.17g}", settings.integrator.max_step));
        ::testing::Test::RecordProperty(std::format("{}_midpoint_samples", level_prefix),
                                        settings.volumetric_samples);
    }
    // Smooth midpoint quadrature is second order, while finite-support edge
    // sampling has a first-order envelope and can change sign as its midpoint
    // phase changes. Each level reduces the maximum affine sample spacing by
    // eight. Require aggregate improvement and the unchanged 0.05% bound.
    EXPECT_LE(errors[1], errors[0] + 1e-5);
    EXPECT_LE(errors[2], errors[1] + 1e-5);
    EXPECT_LE(errors[2], .5 * errors[0] + 1e-5);
    EXPECT_LT(errors[2], 5e-4);
    ::testing::Test::RecordProperty(std::format("{}_geodesic_abs_tolerance", prefix),
                                    std::format("{:.17g}", config.integrator.abs_tolerance));
    ::testing::Test::RecordProperty(std::format("{}_geodesic_rel_tolerance", prefix),
                                    std::format("{:.17g}", config.integrator.rel_tolerance));
    ::testing::Test::RecordProperty(std::format("{}_reference_tau", prefix),
                                    std::format("{:.17g}", reference[0]));
    for (int c = 0; c < 3; ++c)
        ::testing::Test::RecordProperty(std::format("{}_reference_rgb_{}", prefix, c),
                                        std::format("{:.17g}", reference[c + 1]));
}

TEST(VolumeAcceptance, TemperatureAmplitudeRetainsRelativeStefanBoltzmannScale) {
    KerrSchildFamily metric(KerrSchildParams::Schwarzschild(1));
    auto config = VolumeConfig(.25f, 8);
    config.color_mode = Mode::RedshiftMap;  // Independent of blackbody chromaticity.
    GeodesicTracer original(&metric, config);
    const auto baseline = original.Trace(RadialRay());
    ASSERT_FALSE(baseline.numerical_failure);
    ASSERT_TRUE(baseline.volumetric_hit);
    for (float amplitude : {.5f, 2.f}) {
        config.disk_temperature_inner = amplitude;
        GeodesicTracer changed(&metric, config);
        const auto result = changed.Trace(RadialRay());
        ASSERT_FALSE(result.numerical_failure);
        ASSERT_TRUE(result.volumetric_hit);
        EXPECT_FLOAT_EQ(result.optical_depth, baseline.optical_depth);
        for (int c = 0; c < 3; ++c) {
            ASSERT_GT(baseline.volumetric_emission[c], 0);
            EXPECT_NEAR(result.volumetric_emission[c] / baseline.volumetric_emission[c],
                        std::pow(amplitude, 4), 3e-6 * std::pow(amplitude, 4));
        }
    }
}

TEST(VolumeAcceptance, SchwarzschildRadialObservablesConvergeToIndependentTransfer) {
    KerrSchildFamily metric(KerrSchildParams::Schwarzschild(1));
    const auto trace = [&](const TracerConfig& settings) {
        GeodesicTracer tracer(&metric, settings);
        return tracer.Trace(RadialRay());
    };
    ASSERT_NO_FATAL_FAILURE(
        CheckConvergence(trace, DiskTemperatureModel::ShakuraSunyaev, "newtonian"));
    ASSERT_NO_FATAL_FAILURE(
        CheckConvergence(trace, DiskTemperatureModel::NovikovThorne, "page_thorne"));
}

TEST(VolumeAcceptance, KerrPrincipalObservablesConvergeToIndependentTransfer) {
    for (const double spin : {.7, .998}) {
        KerrSchildFamily metric(KerrSchildParams::Kerr(1, spin));
        const auto ray = PrincipalRay(spin);
        ASSERT_NO_FATAL_FAILURE(CheckPrincipalLaunch(metric, ray, spin));
        const auto trace = [&](const TracerConfig& settings) {
            GeodesicTracer tracer(&metric, settings);
            return tracer.Trace(ray);
        };
        for (const auto model :
             {DiskTemperatureModel::ShakuraSunyaev, DiskTemperatureModel::NovikovThorne}) {
            const auto prefix = std::format(
                "spin_{:.3f}_{}", spin,
                model == DiskTemperatureModel::ShakuraSunyaev ? "newtonian" : "page_thorne");
            ASSERT_NO_FATAL_FAILURE(CheckConvergence(trace, model, prefix.c_str(), spin));
        }
    }
}

TEST(VolumeAcceptance, DarkPageThorneInnerBufferRetainsIndependentExtinction) {
    KerrSchildFamily metric(KerrSchildParams::Schwarzschild(1));
    auto config = VolumeConfig(.001f, 8, DiskTemperatureModel::NovikovThorne);
    // An entirely dark but material interval, wholly inside [6,6.006]. The
    // finite absorbing boundary makes this a few-step check rather than
    // evolving a small step all the way to the horizon.
    config.escape_radius = 6.005f;
    config.finite_causal_boundary = true;
    auto ray = RadialRay();
    constexpr double lower = 6.002;
    ray.origin(1) = lower;
    ray.direction(1) = 1;
    GeodesicTracer tracer(&metric, config);
    const auto result = tracer.Trace(ray);
    ASSERT_FALSE(result.numerical_failure);
    EXPECT_EQ(result.outcome, TraceResult::Outcome::Escaped);
    EXPECT_TRUE(result.volumetric_hit);

    const double column_integral = std::sqrt(2 * std::numbers::pi) * std::erf(3 / std::sqrt(2.));
    const auto rate = [&](double radius) {
        return config.volumetric_tau_midplane * std::pow(radius / kInner, -1.5) /
               (column_integral * config.volumetric_scale_height_ratio * radius *
                std::sqrt(1 - 3 / radius));
    };
    const double upper = config.escape_radius;
    const double expected_tau =
        (upper - lower) * (rate(lower) + 4 * rate((lower + upper) / 2) + rate(upper)) / 6;
    ASSERT_GT(expected_tau, 0);
    EXPECT_LT(RelativeError(result.optical_depth, expected_tau), 5e-5);
    for (const auto channel : result.volumetric_emission) EXPECT_FLOAT_EQ(channel, 0);
    RecordProperty("independent_tau", std::format("{:.17g}", expected_tau));
    RecordProperty("relative_tau_error",
                   std::format("{:.17g}", RelativeError(result.optical_depth, expected_tau)));
}

TEST(VolumeAcceptance, OpticallyThinAndCappedLayersPreservePhysicalRayFate) {
    KerrSchildFamily metric(KerrSchildParams::Schwarzschild(1));
    auto config = VolumeConfig(kMaximumSteps.back(), kVolumeSamples.back());
    config.volumetric_tau_midplane = 1e-5f;
    const auto reference = Reference(config, 8192);
    CheckReferenceRefinement(reference, Reference(config, 4096));
    GeodesicTracer thin(&metric, config);
    const auto small = thin.Trace(RadialRay());
    ASSERT_LT(small.optical_depth, .01f);
    EXPECT_LT(CheckTrace(small, reference), 5e-4);
    RecordTraceMetrics(small, reference, "optically_thin");

    config.volumetric_tau_midplane = 100;
    config.volumetric_tau_max = .25f;
    GeodesicTracer capped(&metric, config);
    const auto thick = capped.Trace(RadialRay());
    ASSERT_FALSE(thick.numerical_failure);
    EXPECT_EQ(thick.outcome, TraceResult::Outcome::Horizon);
    EXPECT_FLOAT_EQ(thick.optical_depth, .25f);
    EXPECT_TRUE(thick.volumetric_hit);
    for (float value : thick.volumetric_emission) {
        EXPECT_TRUE(std::isfinite(value));
        EXPECT_GT(value, 0);
    }
}

#ifdef SIRIUS_HAS_RETAINED_COMPUTE
TEST(VolumeAcceptance, RetainedVulkanRadialObservablesMatchIndependentTransfer) {
    using namespace sirius::backend;
    const auto inventory = EnumerateVulkanDevices();
    ASSERT_TRUE(inventory) << inventory.error().Description();
    if (inventory->empty()) GTEST_SKIP() << "No Vulkan device available";
    const auto index = ResolveVulkanDeviceIndex(*inventory);
    ASSERT_TRUE(index) << index.error().Description();
    auto device = CreateVulkanDevice(*index);
    ASSERT_TRUE(device) << device.error().Description();
    auto compute = RetainedCompute::Create(**device, 1);
    RecordProperty("device_admission", compute ? "admitted" : "rejected_comparison_not_run");
    if (!compute) RecordProperty("device_admission_error", compute.error().Description());
    ASSERT_TRUE(compute) << compute.error().Description();
    RetainedTraceExecutor executor(**compute);
    KerrSchildFamily metric(KerrSchildParams::Schwarzschild(1));
    const auto trace = [&](const TracerConfig& settings) {
        GeodesicTracer tracer(&metric, settings);
        tracer.SetStepExecutor(&executor);
        return tracer.Trace(RadialRay());
    };
    CheckConvergence(trace, DiskTemperatureModel::ShakuraSunyaev, "newtonian");
    CheckConvergence(trace, DiskTemperatureModel::NovikovThorne, "page_thorne");
    for (const double spin : {.7, .998}) {
        KerrSchildFamily spinning_metric(KerrSchildParams::Kerr(1, spin));
        const auto ray = PrincipalRay(spin);
        ASSERT_NO_FATAL_FAILURE(CheckPrincipalLaunch(spinning_metric, ray, spin));
        const auto spinning_trace = [&](const TracerConfig& settings) {
            GeodesicTracer tracer(&spinning_metric, settings);
            tracer.SetStepExecutor(&executor);
            return tracer.Trace(ray);
        };
        for (const auto model :
             {DiskTemperatureModel::ShakuraSunyaev, DiskTemperatureModel::NovikovThorne}) {
            const auto prefix = std::format(
                "spin_{:.3f}_{}", spin,
                model == DiskTemperatureModel::ShakuraSunyaev ? "newtonian" : "page_thorne");
            ASSERT_NO_FATAL_FAILURE(CheckConvergence(spinning_trace, model, prefix.c_str(), spin));
        }
    }
    EXPECT_FALSE(executor.Error());
    const auto statistics = executor.Statistics();
    EXPECT_GT(statistics.camera_batches, 0u);
    EXPECT_GT(statistics.interval_batches, 0u);
    EXPECT_GT(statistics.reused_phases, 0u);
    ::testing::Test::RecordProperty("device_name", (*device)->Info().name);
    ::testing::Test::RecordProperty("device_interval_batches",
                                    std::to_string(statistics.interval_batches));
}
#endif
}  // namespace
