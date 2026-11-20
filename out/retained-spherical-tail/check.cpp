#include "sirius/core/source_sky_map.h"

#include "coupled_handoff.h"
#include "support/separated_geodesic_reference.h"
#include <boost/math/quadrature/gauss_kronrod.hpp>

#include <algorithm>
#include <array>
#include <cfenv>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <limits>

namespace {
using sirius::core::CameraLaunch;
using sirius::core::Dual;
using sirius::core::GeodesicVariations;
using sirius::core::IMetric;
using sirius::core::KerrSchildFamily;
using sirius::core::KerrSchildParams;
using sirius::core::Metric4d;
using sirius::core::OutgoingKerrSchild;
using sirius::core::Tensor;
using sirius::core::TensorOps;
using sirius::core::Vec4;
using namespace sirius::core::relativity;
namespace reference = sirius::test::separated_reference;
using reference::Four;
using reference::Scalar;
using reference::Three;

struct Event {
    Four position;
    Four tangent;
    Scalar affine = 0;
};

struct Case {
    Scalar mass;
    Scalar radius;
    Scalar phi;
    Scalar transverse;
    Scalar radial;
    Scalar boundary;
};

double max_event_error = 0;
double max_variation_error = 0;
double max_covariant_roundtrip_error = 0;
double max_source_sky_error = 0;
unsigned panels = 0;
unsigned most_panels = 0;
unsigned enclosed_integrals = 0;

void Check(bool condition, const char* name) {
    if (!condition) {
        std::fprintf(stderr, "Failed: %s\n", name);
        std::abort();
    }
}

void Agree(double actual, Scalar expected, double tolerance, double& maximum, const char* name) {
    const double error =
        static_cast<double>(std::abs(actual - expected) / std::max(1.L, std::abs(expected)));
    maximum = std::max(maximum, error);
    if (!(error <= tolerance)) {
        std::fprintf(stderr, "%s: %.17g expected %.21Lg relative error %.9g limit %.9g\n", name,
                     actual, expected, error, tolerance);
        std::abort();
    }
}

Vec4 Convert(const Four& values) {
    Vec4 result;
    for (int i = 0; i < 4; ++i) result(i) = static_cast<double>(values[i]);
    return result;
}

Case Offset(Case c, unsigned column, Scalar delta) {
    switch (column) {
        case 0:
            c.radius += delta;
            break;
        case 1:
            c.phi += delta;
            break;
        case 2:
            c.transverse += delta;
            break;
        case 3:
            c.radial += delta;
            break;
    }
    return c;
}

Event Launch(const Case& c, bool outgoing) {
    // Non-equatorial launch. The past-directed time component comes from
    // the independent Cartesian metric quadratic, using arbitrary spatial k.
    const Scalar sine = std::sin(c.phi), cosine = std::cos(c.phi);
    const Three direction{.8L * cosine, .8L * sine, .6L};
    Event event;
    event.position = {1.2L, c.radius * direction[0], c.radius * direction[1],
                      c.radius * direction[2]};
    event.tangent = {0, c.radial * direction[0] - c.transverse * sine,
                     c.radial * direction[1] + c.transverse * cosine, c.radial * direction[2]};
    const auto geometry = reference::At(event.position, c.mass, 0, outgoing);
    Scalar mixed = 0;
    for (unsigned i = 1; i < 4; ++i) mixed += geometry.metric[0][i] * event.tangent[i];
    const Scalar spatial = geometry.Inner(event.tangent, event.tangent);
    const Scalar temporal = geometry.metric[0][0];
    event.tangent[0] = (-mixed + std::sqrt(mixed * mixed - temporal * spatial)) / temporal;
    return event;
}

// The production solver uses planar rotation and a quadrature in r. The
// established independent long-double oracle integrates scalar Carter
// equations in Mino time and bisects the first radial event.
Event Oracle(const Case& c, bool outgoing, Scalar refinement) {
    const Event launch = Launch(c, outgoing);
    auto x = launch.position;
    auto k = launch.tangent;
    const auto to_outgoing = [&](Scalar radius) {
        // Spherical specialization of the existing chart's gauge, with t
        // origins equal at r=4M. Only the test converts the independent
        // oracle's public ingoing event to the outgoing chart.
        return -4 * c.mass * std::log((radius - 2 * c.mass) / (2 * c.mass));
    };
    const auto shift_tangent = [&](Scalar radius, const Four& position, const Four& tangent) {
        const Scalar radial =
            (position[1] * tangent[1] + position[2] * tangent[2] + position[3] * tangent[3]) /
            radius;
        return 4 * c.mass * radial / (radius - 2 * c.mass);
    };
    if (outgoing && c.mass > 0) {
        x[0] -= to_outgoing(c.radius);
        k[0] += shift_tangent(c.radius, x, k);
    }
    CameraLaunch packet{};
    packet.position = Convert(x);
    packet.tangent = Convert(k);
    const auto result = reference::Trace(packet, c.mass, 0, c.boundary, refinement);
    Check(result.fate == reference::Fate::Escape && result.radial_turns == 0,
          "independent first outward event");
    Event event{result.x, result.k, result.affine};
    if (outgoing && c.mass > 0) {
        event.position[0] += to_outgoing(c.boundary);
        event.tangent[0] -= shift_tangent(c.boundary, event.position, event.tangent);
    }
    return event;
}

void CheckCovariantResult(const Case& test, bool outgoing, const Vec4& position,
                          const Vec4& tangent, const std::array<SphericalTailVariation, 4>& columns,
                          const SphericalTailResult& direct, const Event& oracle_terminal,
                          Scalar epsilon) {
    KerrSchildFamily family(KerrSchildParams::Schwarzschild(test.mass));
    OutgoingKerrSchild outgoing_metric(family);
    IMetric& metric = outgoing ? static_cast<IMetric&>(outgoing_metric) : family;
    Metric4d g;
    Tensor<Dual<double>, 4, 4, 4> dg;
    metric.Evaluate(position, g, dg);
    const auto start_connection = TensorOps::Christoffel(g, dg);
    GeodesicVariations covariant;
    for (unsigned i = 0; i < columns.size(); ++i) {
        const auto& seed = columns[i];
        covariant[i] = {
            seed.position,
            seed.tangent - NegConnectionContraction(start_connection, tangent, seed.position)};
    }
    const auto trial = PrepareSphericalCoupledTail(metric, position, tangent, covariant,
                                                   test.radius * .9L, test.boundary);
    Check(trial.has_value(), "coupled terminal handoff");
    metric.Evaluate(trial->position, g, dg);
    const auto end_connection = TensorOps::Christoffel(g, dg);
    std::array<Vec4, 2> displacements, derivatives;
    for (unsigned i = 0; i < covariant.size(); ++i) {
        const auto& variation = trial->variations[i];
        const auto event =
            VarySphericalEndpoint(trial->position, trial->tangent, variation.displacement);
        Check(event.has_value(), "ordinary finite event reconstruction");
        if (i < displacements.size()) {
            displacements[i] = event->displacement;
            derivatives[i] = variation.derivative;
        }
        const Vec4 derivative =
            variation.derivative +
            NegConnectionContraction(end_connection, trial->tangent, event->displacement);
        for (int j = 0; j < 4; ++j) {
            Agree(event->displacement(j), direct.event_variations[i].position(j), 5e-12,
                  max_covariant_roundtrip_error, "fixed-affine deviation roundtrip");
            Agree(derivative(j), direct.event_variations[i].tangent(j), 5e-12,
                  max_covariant_roundtrip_error, "covariant momentum roundtrip");
        }
    }
    Metric4d inverse;
    Check(metric.InverseMetric(trial->position, inverse), "source-sky inverse metric");
    const auto map =
        MeasureSourceSkyAngularMap(g, inverse, dg, trial->tangent, displacements, derivatives, 1.0);
    Check(map.has_value(), "source-sky angular map");
    const auto measured_sky = [&](const Event& event) {
        return reference::At(event.position, test.mass, 0, outgoing).Sky(event.tangent);
    };
    const auto expected_direction = measured_sky(oracle_terminal);
    const auto basis = reference::Basis(expected_direction);
    for (unsigned axis = 0; axis < expected_direction.size(); ++axis)
        Agree(map->direction[axis], expected_direction[axis], 3e-10, max_source_sky_error,
              "finite source-sky direction");
    for (unsigned column = 0; column < displacements.size(); ++column) {
        const auto plus = measured_sky(Oracle(Offset(test, column, epsilon), outgoing, .001L));
        const auto minus = measured_sky(Oracle(Offset(test, column, -epsilon), outgoing, .001L));
        Three difference;
        for (unsigned axis = 0; axis < difference.size(); ++axis)
            difference[axis] = (plus[axis] - minus[axis]) / (2 * epsilon);
        for (unsigned row = 0; row < basis.size(); ++row)
            Agree(map->jacobian[row][column], reference::Dot(basis[row], difference), 2e-7,
                  max_source_sky_error, "finite source-sky Jacobian");
    }
}

void Exercise(const Case& test, bool outgoing) {
    const auto input = Launch(test, outgoing);
    const auto chart =
        outgoing ? SphericalKerrSchildChart::Outgoing : SphericalKerrSchildChart::Ingoing;
    const auto position = Convert(input.position);
    const auto tangent = Convert(input.tangent);
    std::array<SphericalTailVariation, 4> variations;
    const Scalar epsilon = 0.00004L;
    const auto delta = [&](const Event& a, const Event& b) {
        Event result;
        for (unsigned i = 0; i < 4; ++i) {
            result.position[i] = (a.position[i] - b.position[i]) / (2 * epsilon);
            result.tangent[i] = (a.tangent[i] - b.tangent[i]) / (2 * epsilon);
        }
        result.affine = (a.affine - b.affine) / (2 * epsilon);
        return result;
    };
    for (unsigned column = 0; column < variations.size(); ++column) {
        const auto varied = delta(Launch(Offset(test, column, epsilon), outgoing),
                                  Launch(Offset(test, column, -epsilon), outgoing));
        variations[column] = {Convert(varied.position), Convert(varied.tangent)};
    }
    const auto result = ContinueSphericalTail(position, tangent, variations, test.mass, chart,
                                              test.radius * .9L, test.boundary);
    if (!result) {
        std::fprintf(stderr, "tail failure %u at radius %.9Lg, boundary %.9Lg\n",
                     static_cast<unsigned>(result.error()), test.radius, test.boundary);
        std::abort();
    }
    panels += result->panels;
    most_panels = std::max(most_panels, result->panels);

    const auto observed = Oracle(test, outgoing, .001L);
    const auto refined = Oracle(test, outgoing, .0005L);
    CheckCovariantResult(test, outgoing, position, tangent, variations, *result, refined, epsilon);
    for (int i = 0; i < 4; ++i) {
        Agree(observed.position[i], refined.position[i], 2e-10, max_event_error,
              "reference position convergence");
        Agree(observed.tangent[i], refined.tangent[i], 2e-10, max_event_error,
              "reference tangent convergence");
        Agree(result->position(i), refined.position[i], 3e-10, max_event_error,
              "tail event position");
        Agree(result->tangent(i), refined.tangent[i], 3e-10, max_event_error, "tail event tangent");
    }
    Agree(result->affine_length, refined.affine, 3e-10, max_event_error, "affine length");
    for (unsigned column = 0; column < variations.size(); ++column) {
        const auto independent = delta(Oracle(Offset(test, column, epsilon), outgoing, .001L),
                                       Oracle(Offset(test, column, -epsilon), outgoing, .001L));
        for (int i = 0; i < 4; ++i) {
            Agree(result->event_variations[column].position(i), independent.position[i], 2e-7,
                  max_variation_error, "fixed-radius position derivative");
            Agree(result->event_variations[column].tangent(i), independent.tangent[i], 2e-7,
                  max_variation_error, "fixed-radius tangent derivative");
        }
        Agree(result->affine_variations[column], independent.affine, 2e-7, max_variation_error,
              "affine derivative");
    }
}

void DeclinesAndRadial() {
    const auto chart = SphericalKerrSchildChart::Outgoing;
    const auto x = Convert({0, 0, 0, 25});
    const auto k = Convert({-1, 0, 0, 1});
    std::array<SphericalTailVariation, 4> columns{};
    // This seed opens angular momentum from zero along the polar axis.
    columns[0].tangent(1) = 1;
    const auto radial = ContinueSphericalTail(x, k, columns, 0, chart, 0, 80);
    Check(radial.has_value(), "axis has no polar singularity");
    Agree(radial->position(0), -55, 2e-12, max_event_error, "flat axial time");
    Agree(radial->position(3), 80, 2e-12, max_event_error, "flat axial boundary");
    Agree(radial->event_variations[0].position(1), 55, 2e-10, max_variation_error, "radial seed");
    Agree(radial->event_variations[0].tangent(1), 1, 2e-11, max_variation_error,
          "radial tangent seed");

    auto guarded = ContinueSphericalTail(x, k, columns, 0, chart, 25, 80);
    Check(!guarded && guarded.error() == FiniteTailFailure::Domain, "source support boundary");
    guarded = ContinueSphericalTail(x, k, columns, 0, chart, 0, 25);
    Check(!guarded && guarded.error() == FiniteTailFailure::Domain,
          "target sphere already reached");
    guarded = ContinueSphericalTail(x, Convert({-1, 0, 0, -1}), columns, 0, chart, 0, 80);
    Check(!guarded && guarded.error() == FiniteTailFailure::Domain, "inward tangent");
    columns[3].tangent(0) = 1;
    guarded = ContinueSphericalTail(x, k, columns, 0, chart, 0, 80);
    Check(!guarded && guarded.error() == FiniteTailFailure::Constraint,
          "differentiated null constraint");
    columns[3].tangent(0) = 0;
    guarded = ContinueSphericalTail(x, Convert({-1, 0, 0, .9}), columns, 0, chart, 0, 80);
    Check(!guarded && guarded.error() == FiniteTailFailure::Constraint, "central null constraint");
    guarded = ContinueSphericalTail(x, Convert({-1, 1, 0, 1e-5}), columns, 0, chart, 0, 80);
    Check(!guarded && guarded.error() == FiniteTailFailure::Constraint,
          "null residual cannot conceal a poorly matched radial root");
    guarded = ContinueSphericalTail(x, k, columns, 0, chart, 0, 80, {.maximum_panels = 0});
    Check(!guarded && guarded.error() == FiniteTailFailure::WorkLimit,
          "zero quadrature work allowance");
    guarded =
        ContinueSphericalTail(x, k, columns, 0, chart, 0, std::numeric_limits<double>::infinity());
    Check(!guarded && guarded.error() == FiniteTailFailure::Input, "finite target required");

    KerrSchildFamily rotating(KerrSchildParams::Kerr(1, .1));
    const auto ineligible = PrepareSphericalCoupledTail(rotating, x, k, {}, 0, 80);
    Check(!ineligible && ineligible.error() == FiniteTailFailure::Domain,
          "full rotating metric declines");
    Check(std::fesetround(FE_UPWARD) == 0, "set alternate rounding for control");
    guarded = ContinueSphericalTail(x, k, columns, 0, chart, 0, 80);
    Check(std::fesetround(FE_TONEAREST) == 0, "restore original rounding");
    Check(!guarded && guarded.error() == FiniteTailFailure::Input, "rounding environment decline");
}

void CheckIntegralEnclosures() {
    using Complex = std::complex<Scalar>;
    using OracleQuadrature = boost::math::quadrature::gauss_kronrod<Scalar, 31>;
    using namespace finite_tail_detail;
    const SphericalTailControl control;
    // Exact input words for the interval calculation. Complex-step
    // differentiation supplies sensitivities independently of its partials.
    for (const auto& fixture :
         {std::array{22., 200., 1., 144., 1.}, std::array{7., 120., 1., 60., -1.},
          std::array{21., 200., 1., 0., -1.}, std::array{40., 200., 0., 670., 1.}}) {
        const auto [start, end, mass, impact, orientation] = fixture;
        Quadrature trial{mass, orientation, impact, start, end, control};
        const auto bounds = trial.Run();
        Check(bounds.has_value(), "radial enclosure quadrature");
        constexpr Scalar perturbation = 1e-30L;
        for (unsigned i = 0; i < 3; ++i) {
            const auto kernel = [=](Scalar radius) -> Complex {
                const Scalar beta = 2 * mass / radius;
                const Scalar lapse = 1 - beta;
                const Complex q =
                    Complex(1) - Complex(impact, perturbation) * lapse / (radius * radius);
                const Complex speed = std::sqrt(q);
                if (i == 0) return Scalar(1) / speed;
                if (i == 1) return Scalar(1) / (radius * radius * speed);
                return (Complex(orientation * beta) - Scalar(1) / speed) / lapse;
            };
            const auto direct = OracleQuadrature::integrate(
                [&](Scalar r) { return kernel(r).real(); }, start, end, 12, 1e-16L);
            const auto derivative = OracleQuadrature::integrate(
                [&](Scalar r) { return kernel(r).imag() / perturbation; }, start, end, 12, 1e-16L);
            Check((*bounds)[i].lo <= direct && direct <= (*bounds)[i].hi,
                  "independent integral inside certified interval");
            Check((*bounds)[i + 3].lo <= derivative && derivative <= (*bounds)[i + 3].hi,
                  "complex-step sensitivity inside certified interval");
            enclosed_integrals += 2;
        }
    }
}

void MeasureProvisionalCost() {
    constexpr unsigned repetitions = 32;
    using Clock = std::chrono::steady_clock;
    const auto start = Clock::now();
    unsigned attempted_panels = 0;
    for (unsigned i = 0; i < repetitions; ++i) {
        const Case test{1, 25 + (i % 4) * 10.L, .35L, .04L * i, .85L, 200};
        const auto input = Launch(test, false);
        const auto result =
            ContinueSphericalTail(Convert(input.position), Convert(input.tangent), {}, test.mass,
                                  SphericalKerrSchildChart::Ingoing, 20, test.boundary);
        Check(result.has_value(), "representative exterior cost fixture");
        attempted_panels += result->panels;
    }
    const double ms = std::chrono::duration<double, std::milli>(Clock::now() - start).count();
    std::printf("diagnostic_tails=%u diagnostic_panels=%u diagnostic_ms=%.6f\n", repetitions,
                attempted_panels, ms);
}
}  // namespace

int main() {
    DeclinesAndRadial();
    CheckIntegralEnclosures();
    constexpr std::array cases{
        Case{0, 22, .35L, .2L, .8L, 200},    Case{1, 21, .35L, 0, .8L, 200},
        Case{1, 22, .35L, .2L, .8L, 200},    Case{1, 50, 1.6L, -.54L, .7L, 500},
        Case{1, 60, -2.1L, .8L, .21L, 140},  Case{1, 7, -.8L, .54L, .16L, 120},
        Case{8, 168, -.8L, .33L, .6L, 1600},
    };
    for (const auto& test : cases)
        for (const bool outgoing : {false, true}) Exercise(test, outgoing);
    std::printf(
        "fixture cases=%zu charts=2 max_scaled_event_error=%.9g max_scaled_variation_error=%.9g "
        "covariant_roundtrip_error=%.9g source_sky_error=%.9g enclosed_integrals=%u "
        "total_panels=%u "
        "max_panels=%u\n",
        cases.size(), max_event_error, max_variation_error, max_covariant_roundtrip_error,
        max_source_sky_error, enclosed_integrals, panels, most_panels);
    MeasureProvisionalCost();
}
