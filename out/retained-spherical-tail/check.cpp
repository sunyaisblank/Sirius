#include "sirius/core/polarisation/walker_penrose.h"
#include "sirius/core/source_sky_map.h"
#include "sirius/core/spherical_continuation/coupled_handoff.h"

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
using Complex = std::complex<Scalar>;
using ComplexFour = std::array<Complex, 4>;
using ComplexThree = std::array<Complex, 3>;
using ComplexMetric = std::array<ComplexFour, 4>;
using OracleQuadrature = boost::math::quadrature::gauss_kronrod<Scalar, 31>;
constexpr Scalar kComplexStep = 1e-28L;

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
    Scalar source_support = 0;
};

double SourceSupport(const Case& test) {
    return static_cast<double>(test.source_support > 0 ? test.source_support : test.radius * .9L);
}

SphericalTailPath Path(const Case& test) {
    return test.radial < 0 ? SphericalTailPath::ExteriorTurn : SphericalTailPath::Outward;
}

double max_event_error = 0;
double max_variation_error = 0;
double max_covariant_roundtrip_error = 0;
double max_source_sky_error = 0;
unsigned panels = 0;
unsigned most_panels = 0;
unsigned enclosed_integrals = 0;
unsigned enclosed_events = 0;
double max_scaled_half_width = 0;

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

Four Promote(const Vec4& values) { return {values(0), values(1), values(2), values(3)}; }

void Contains(SphericalTailBound bound, Scalar expected, const char* name) {
    if (!(bound.lower <= expected && expected <= bound.upper)) {
        std::fprintf(stderr, "%s: %.21Lg outside [%.17g, %.17g]\n", name, expected, bound.lower,
                     bound.upper);
        std::abort();
    }
    const double width = (bound.upper - bound.lower) * 0.5;
    max_scaled_half_width = std::max(
        max_scaled_half_width, width / std::max(1.0, std::abs(static_cast<double>(expected))));
    ++enclosed_events;
}

Complex Radius(const ComplexFour& x) { return std::sqrt(x[1] * x[1] + x[2] * x[2] + x[3] * x[3]); }

ComplexMetric ReferenceMetric(const ComplexFour& x, Scalar mass, bool outgoing) {
    const Complex r = Radius(x);
    const ComplexFour ell{outgoing ? -1.L : 1.L, x[1] / r, x[2] / r, x[3] / r};
    ComplexMetric metric;
    for (unsigned i = 0; i < ell.size(); ++i)
        for (unsigned j = 0; j < ell.size(); ++j) {
            metric[i][j] = 2 * mass / r * ell[i] * ell[j];
            if (i == j) metric[i][j] += i == 0 ? -1.L : 1.L;
        }
    return metric;
}

Four ReferenceConnection(const Four& x, const Four& a, const Four& b, Scalar mass, bool outgoing) {
    const ComplexFour event{x[0], x[1], x[2], x[3]};
    const auto metric = ReferenceMetric(event, mass, outgoing);
    std::array<ComplexMetric, 4> differentiated;
    for (unsigned i = 0; i < event.size(); ++i) {
        auto varied = event;
        varied[i] += Complex(0, kComplexStep);
        differentiated[i] = ReferenceMetric(varied, mass, outgoing);
    }
    // Solve g C = Γ_first-kind(a,b), using the complex-step metric derivatives.
    std::array<std::array<Scalar, 5>, 4> system{};
    for (unsigned i = 0; i < event.size(); ++i)
        for (unsigned j = 0; j < event.size(); ++j) {
            system[i][j] = metric[i][j].real();
            for (unsigned k = 0; k < event.size(); ++k) {
                const auto derivatives =
                    differentiated[j][i][k] + differentiated[k][i][j] - differentiated[i][j][k];
                system[i][4] += derivatives.imag() * a[j] * b[k] / (2 * kComplexStep);
            }
        }
    for (unsigned i = 0; i < event.size(); ++i) {
        const Scalar pivot = system[i][i];
        for (auto& entry : system[i]) entry /= pivot;
        for (unsigned j = 0; j < event.size(); ++j) {
            if (i == j) continue;
            const Scalar factor = system[j][i];
            for (unsigned k = 0; k < system[j].size(); ++k) system[j][k] -= factor * system[i][k];
        }
    }
    return {system[0][4], system[1][4], system[2][4], system[3][4]};
}

ComplexThree Cross(const ComplexThree& a, const ComplexThree& b) {
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}

struct ComplexEvent {
    ComplexFour position;
    ComplexFour tangent;
    Complex affine;
};

ComplexEvent ExactInputReference(const Four& x0, const Four& k0, const Four& dx, const Four& dk,
                                 Scalar mass, Scalar target, bool outgoing, bool through_turn) {
    ComplexFour x, k;
    for (unsigned i = 0; i < x.size(); ++i) {
        x[i] = Complex(x0[i], kComplexStep * dx[i]);
        k[i] = Complex(k0[i], kComplexStep * dk[i]);
    }
    const Complex r0 = Radius(x);
    const auto metric = ReferenceMetric(x, mass, outgoing);
    Complex energy = 0;
    for (unsigned i = 0; i < x.size(); ++i) energy += metric[0][i] * k[i];
    const auto angular_momentum = Cross({x[1], x[2], x[3]}, {k[1], k[2], k[3]});
    const Complex length = std::sqrt(angular_momentum[0] * angular_momentum[0] +
                                     angular_momentum[1] * angular_momentum[1] +
                                     angular_momentum[2] * angular_momentum[2]);
    const Complex impact = length / energy;
    const Scalar orientation = outgoing ? -1.L : 1.L;
    const auto speed = [&](Complex r) {
        return std::sqrt(1.L - (1.L - 2 * mass / r) * impact * impact / (r * r));
    };
    const auto integrate_complex = [&](const auto& integrand) -> Complex {
        return {OracleQuadrature::integrate([&](Scalar u) { return integrand(u).real(); }, 0, 1, 12,
                                            1e-16L),
                kComplexStep * OracleQuadrature::integrate(
                                   [&](Scalar u) { return integrand(u).imag() / kComplexStep; }, 0,
                                   1, 12, 1e-16L)};
    };
    const auto radial_integral = [&](const auto& kernel) {
        return integrate_complex(
            [&](Scalar u) { return (target - r0) * kernel(r0 + (target - r0) * u); });
    };
    std::array<Complex, 3> integrals;
    if (through_turn) {
        // The reference finds the root by high-precision complex Newton and
        // regularizes it in inverse radius. The candidate uses an interval
        // bracket and the Cartesian radius difference instead.
        Complex turn = impact;
        const Complex b2 = impact * impact;
        for (unsigned i = 0; i < 20; ++i)
            turn -= (turn * turn * turn / b2 - turn + 2 * mass) / (3.L * turn * turn / b2 - 1.L);
        const Complex q = 2 * mass / turn;
        for (unsigned i = 0; i < integrals.size(); ++i) {
            const auto leg = [&](Complex end, Scalar sign) {
                const Complex scale = std::sqrt(1.L - turn / end);
                return integrate_complex([&](Scalar t) {
                    const Complex z = t * scale;
                    const Complex u = 1.L - z * z;
                    const Complex p = ((1.L + u) - q * (1.L + u + u * u)) / (1.L - q);
                    const Complex jacobian = 2.L * turn * scale / (u * u);
                    const Complex measure = jacobian / std::sqrt(p);
                    if (i == 0) return measure;
                    if (i == 1) return measure * u * u / (turn * turn);
                    return (sign * q * u * jacobian * z - measure) / (1.L - q * u);
                });
            };
            integrals[i] = leg(r0, -orientation) + leg(target, orientation);
        }
    } else {
        integrals[0] = radial_integral([&](Complex r) { return 1.L / speed(r); });
        integrals[1] = radial_integral([&](Complex r) { return 1.L / (r * r * speed(r)); });
        integrals[2] = radial_integral([&](Complex r) {
            const Complex h = 2 * mass / r;
            return (orientation * h - 1.L / speed(r)) / (1.L - h);
        });
    }
    ComplexEvent result;
    result.affine = integrals[0] / energy;
    result.position[0] = x[0] + integrals[2];
    const Complex angle = impact * integrals[1];
    const ComplexThree initial{x[1] / r0, x[2] / r0, x[3] / r0};
    const auto angular = Cross(angular_momentum, initial);
    ComplexThree normal;
    const Complex rate = energy * speed(target);
    const Scalar h = 2 * mass / target;
    result.tangent[0] = (orientation * h * rate - energy) / (1 - h);
    for (unsigned i = 0; i < normal.size(); ++i) {
        // sin(angle)/length has a finite limit at a radial launch.
        normal[i] = initial[i] * std::cos(angle);
        if (length != Complex{}) normal[i] += angular[i] * std::sin(angle) / length;
        result.position[i + 1] = target * normal[i];
    }
    const auto terminal_angular = Cross(angular_momentum, normal);
    for (unsigned i = 0; i < normal.size(); ++i)
        result.tangent[i + 1] = rate * normal[i] + terminal_angular[i] / target;
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
    const unsigned expected_turns = c.mass > 0 && c.radial < 0 ? 1U : 0U;
    if (result.fate != reference::Fate::Escape || result.radial_turns != expected_turns) {
        std::fprintf(stderr,
                     "Carter event mismatch at M=%Lg r=%Lg kr=%Lg outgoing=%d fate=%d turns=%u\n",
                     c.mass, c.radius, c.radial, outgoing, static_cast<int>(result.fate),
                     result.radial_turns);
        std::abort();
    }
    Event event{result.x, result.k, result.affine};
    if (outgoing && c.mass > 0) {
        event.position[0] += to_outgoing(c.boundary);
        event.tangent[0] -= shift_tangent(c.boundary, event.position, event.tangent);
    }
    return event;
}

void CheckEventBounds(const SphericalTailResult& candidate, const Case& test, bool outgoing,
                      const Vec4& position, const Vec4& tangent,
                      const std::array<SphericalTailVariation, 4>& columns) {
    const auto x = Promote(position), k = Promote(tangent);
    const auto center =
        ExactInputReference(x, k, {}, {}, test.mass, test.boundary, outgoing, test.radial < 0);
    Contains(candidate.event_bounds.affine_length, center.affine.real(), "affine enclosure");
    for (unsigned j = 0; j < x.size(); ++j) {
        Contains(candidate.event_bounds.position[j], center.position[j].real(),
                 "position enclosure");
        Contains(candidate.event_bounds.tangent[j], center.tangent[j].real(), "tangent enclosure");
    }
    for (unsigned i = 0; i < columns.size(); ++i) {
        const auto varied =
            ExactInputReference(x, k, Promote(columns[i].position), Promote(columns[i].tangent),
                                test.mass, test.boundary, outgoing, test.radial < 0);
        const auto& bound = candidate.variation_bounds[i];
        Contains(bound.affine_length, varied.affine.imag() / kComplexStep,
                 "affine sensitivity enclosure");
        for (unsigned j = 0; j < x.size(); ++j) {
            Contains(bound.position[j], varied.position[j].imag() / kComplexStep,
                     "position sensitivity enclosure");
            Contains(bound.tangent[j], varied.tangent[j].imag() / kComplexStep,
                     "tangent sensitivity enclosure");
        }
    }
}

void CheckCoupledBounds(const SphericalCoupledTail& candidate, const Case& test, bool outgoing,
                        const Vec4& position, const Vec4& tangent,
                        const GeodesicVariations& covariant) {
    const auto x = Promote(position), k = Promote(tangent);
    const auto center =
        ExactInputReference(x, k, {}, {}, test.mass, test.boundary, outgoing, test.radial < 0);
    Four end_x, end_k;
    for (unsigned j = 0; j < x.size(); ++j) {
        end_x[j] = center.position[j].real();
        end_k[j] = center.tangent[j].real();
    }
    for (unsigned i = 0; i < covariant.size(); ++i) {
        const Four dx = Promote(covariant[i].displacement);
        Four dk = Promote(covariant[i].derivative);
        const auto start_connection = ReferenceConnection(x, k, dx, test.mass, outgoing);
        for (unsigned j = 0; j < x.size(); ++j) dk[j] -= start_connection[j];
        const auto varied =
            ExactInputReference(x, k, dx, dk, test.mass, test.boundary, outgoing, test.radial < 0);
        Four event_dx;
        for (unsigned j = 0; j < x.size(); ++j)
            event_dx[j] = varied.position[j].imag() / kComplexStep;
        const auto end_connection =
            ReferenceConnection(end_x, end_k, event_dx, test.mass, outgoing);
        for (unsigned j = 0; j < x.size(); ++j) {
            const Scalar fixed_affine =
                event_dx[j] - end_k[j] * varied.affine.imag() / kComplexStep;
            const Scalar derivative = varied.tangent[j].imag() / kComplexStep + end_connection[j];
            Contains(candidate.variation_bounds[i].displacement[j], fixed_affine,
                     "fixed-affine displacement enclosure");
            Contains(candidate.variation_bounds[i].derivative[j], derivative,
                     "covariant derivative enclosure");
        }
    }
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
    const auto trial = PrepareSphericalCoupledTail(
        metric, position, tangent, covariant, SourceSupport(test), test.boundary, {}, Path(test));
    if (!trial) {
        std::fprintf(stderr, "coupled tail failure %u at radius %.9Lg, boundary %.9Lg, chart %s\n",
                     static_cast<unsigned>(trial.error()), test.radius, test.boundary,
                     outgoing ? "outgoing" : "ingoing");
        std::abort();
    }
    CheckCoupledBounds(*trial, test, outgoing, position, tangent, covariant);
    if (test.radial < 0) {
        // Exercise the input enclosure used when the exact launch may lie
        // between two binary64 words. The complete event and four covariant
        // columns must still contain the independent result at this launch.
        using namespace finite_tail_detail;
        const auto surround = [](const Vec4& value) {
            return finite_tail_detail::Four{Widen(value(0), value(0)), Widen(value(1), value(1)),
                                            Widen(value(2), value(2)), Widen(value(3), value(3))};
        };
        CoupledLaunch enclosed;
        enclosed.position = surround(position);
        enclosed.tangent = surround(tangent);
        for (unsigned i = 0; i < covariant.size(); ++i) {
            enclosed.displacement[i] = surround(covariant[i].displacement);
            enclosed.derivative[i] = surround(covariant[i].derivative);
        }
        const auto bounded =
            ContinueCoupled(metric, enclosed, SourceSupport(test), test.boundary, {}, Path(test));
        Check(bounded.has_value(), "enclosed launch state is accepted");
        CheckCoupledBounds(*bounded, test, outgoing, position, tangent, covariant);
    }
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
                                              SourceSupport(test), test.boundary, {}, Path(test));
    if (!result) {
        std::fprintf(stderr, "tail failure %u at radius %.9Lg, boundary %.9Lg\n",
                     static_cast<unsigned>(result.error()), test.radius, test.boundary);
        std::abort();
    }
    panels += result->panels;
    most_panels = std::max(most_panels, result->panels);
    Check(result->periapsis.has_value() == (test.radial < 0), "expected radial branch");
    if (result->periapsis)
        Check(result->periapsis->lower > SourceSupport(test), "turn excludes source support");
    CheckEventBounds(*result, test, outgoing, position, tangent, variations);

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
        ContinueSphericalTail(x, k, columns, 0, chart, 0, 80,
                              {.event_absolute_tolerance = 1e-18, .event_relative_tolerance = 0});
    Check(!guarded && guarded.error() == FiniteTailFailure::EventAccuracy,
          "event enclosure must satisfy its budget");
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
    using namespace finite_tail_detail;
    const SphericalTailControl control;
    // Exact input words for the interval calculation. Complex-step
    // differentiation supplies sensitivities independently of its partials.
    for (const auto& fixture :
         {std::array{22., 200., 1., 144., 1.}, std::array{7., 120., 1., 60., -1.},
          std::array{21., 200., 1., 0., -1.}, std::array{40., 200., 0., 670., 1.}}) {
        const auto [start, end, mass, impact, orientation] = fixture;
        Quadrature trial{mass, orientation, Range(impact), Range(start), end, control};
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

void CheckCoupledCancellationAndWorkLimit() {
    KerrSchildFamily metric(KerrSchildParams::Schwarzschild(1));
    for (const auto& fixture :
         {Case{1, 25, .35L, .55L, .8L, 200, 20}, Case{1, 50, .35L, .55L, -.8L, 200, 20}}) {
        const auto launch = Launch(fixture, false);
        const auto x = Convert(launch.position), k = Convert(launch.tangent);
        const auto solve = [&](const SphericalTailControl& control) {
            return PrepareSphericalCoupledTail(metric, x, k, {}, SourceSupport(fixture),
                                               fixture.boundary, control, Path(fixture));
        };
        unsigned checks = 0;
        const std::function<bool()> cancel = [&] { return ++checks >= 8; };
        const auto interrupted = solve({.should_cancel = &cancel});
        Check(!interrupted && interrupted.error() == FiniteTailFailure::Cancelled && checks == 8,
              "cancellation exits in-progress continuation");
        const auto retry = solve({});
        Check(retry.has_value() && retry->panels > 8, "fresh continuation after cancellation");
        const auto limited = solve({.maximum_panels = retry->panels - 1});
        Check(!limited && limited.error() == FiniteTailFailure::WorkLimit,
              "work budget includes rejected quadrature panels");
        const auto complete = solve({.maximum_panels = retry->panels});
        Check(complete.has_value() && complete->panels == retry->panels,
              "continuation accepts the exact budget");
    }
    const auto launch = Launch({1, 50, .35L, .55L, -.8L, 200, 20}, false);
    const auto would_enter_sources =
        PrepareSphericalCoupledTail(metric, Convert(launch.position), Convert(launch.tangent), {},
                                    40, 200, {}, SphericalTailPath::ExteriorTurn);
    Check(!would_enter_sources && would_enter_sources.error() == FiniteTailFailure::Domain,
          "nearer periapsis cannot pass through source support");
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
    CheckCoupledCancellationAndWorkLimit();
    constexpr std::array cases{
        Case{0, 22, .35L, .2L, .8L, 200},        Case{1, 21, .35L, 0, .8L, 200},
        Case{1, 22, .35L, .2L, .8L, 200},        Case{1, 50, 1.6L, -.54L, .7L, 500},
        Case{1, 60, -2.1L, .8L, .21L, 140},      Case{1, 7, -.8L, .54L, .16L, 120},
        Case{8, 168, -.8L, .33L, .6L, 1600},     Case{0, 50, .65L, .65L, -.8L, 200, 20},
        Case{1, 50, -.35L, .75L, -.6L, 200, 20}, Case{1, 30, -.8L, .46L, -.9L, 160, 8},
        Case{8, 168, .49L, .6L, -.7L, 1200, 75},
    };
    for (const auto& test : cases)
        for (const bool outgoing : {false, true}) Exercise(test, outgoing);
    std::printf(
        "fixture cases=%zu charts=2 max_scaled_event_error=%.9g max_scaled_variation_error=%.9g "
        "covariant_roundtrip_error=%.9g source_sky_error=%.9g enclosed_integrals=%u "
        "enclosed_event_components=%u max_scaled_half_width=%.9g total_panels=%u "
        "max_panels=%u\n",
        cases.size(), max_event_error, max_variation_error, max_covariant_roundtrip_error,
        max_source_sky_error, enclosed_integrals, enclosed_events, max_scaled_half_width, panels,
        most_panels);
    MeasureProvisionalCost();
}
