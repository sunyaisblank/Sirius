// Independent field-equation gates for the production Cartesian
// Kerr-Schild family. The curvature side differentiates the live analytic
// connection; the matter side differentiates the standard Kerr-Newman
// electromagnetic potential and never reuses the production H function.

#include "sirius/core/metrics/kerr_schild_family.h"

#include <gtest/gtest.h>

#include "metric_curvature_oracle.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <string_view>

namespace sirius::test {
using namespace sirius::core;

namespace {

using Matrix4 = std::array<std::array<double, 4>, 4>;

Vec4 KerrSchildPoint(double radius, double theta, double phi, double spin) {
    Vec4 point;
    const double cylindrical = std::sqrt(radius * radius + spin * spin) * std::sin(theta);
    point(1) = cylindrical * std::cos(phi);
    point(2) = cylindrical * std::sin(phi);
    point(3) = radius * std::cos(theta);
    return point;
}

double IndependentKerrRadius(const Vec4& point, double spin) {
    const double x = point(1);
    const double y = point(2);
    const double z = point(3);
    const double spin_squared = spin * spin;
    const double reduced = x * x + y * y + z * z - spin_squared;
    const double radius_squared =
        0.5 * (reduced + std::sqrt(reduced * reduced + 4.0 * spin_squared * z * z));
    return std::sqrt(radius_squared);
}

// A_mu = -Q r l_mu / Sigma in Cartesian Kerr-Schild coordinates, where
// Sigma=r^2+a^2 cos^2(theta). A global sign convention is immaterial to the
// stress tensor and to the source-free Maxwell equation.
std::array<double, 4> KerrNewmanPotential(const Vec4& point, double spin, double charge) {
    if (charge == 0.0) return {};
    const double radius = IndependentKerrRadius(point, spin);
    const double radius_squared = radius * radius;
    const double spin_squared = spin * spin;
    const double denominator = radius_squared + spin_squared;
    const double cosine = point(3) / radius;
    const double sigma = radius_squared + spin_squared * cosine * cosine;
    const double l[4] = {
        1.0,
        (radius * point(1) + spin * point(2)) / denominator,
        (radius * point(2) - spin * point(1)) / denominator,
        point(3) / radius,
    };
    std::array<double, 4> potential{};
    for (int component = 0; component < 4; ++component) {
        potential[component] = -charge * radius * l[component] / sigma;
    }
    return potential;
}

double PotentialDerivative(const Vec4& point, double spin, double charge, int derivative,
                           int component, double h) {
    if (derivative == 0) return 0.0;
    Vec4 plus_one = point;
    Vec4 minus_one = point;
    Vec4 plus_two = point;
    Vec4 minus_two = point;
    plus_one(derivative) += h;
    minus_one(derivative) -= h;
    plus_two(derivative) += 2.0 * h;
    minus_two(derivative) -= 2.0 * h;
    return (-KerrNewmanPotential(plus_two, spin, charge)[component] +
            8.0 * KerrNewmanPotential(plus_one, spin, charge)[component] -
            8.0 * KerrNewmanPotential(minus_one, spin, charge)[component] +
            KerrNewmanPotential(minus_two, spin, charge)[component]) /
           (12.0 * h);
}

Matrix4 MaxwellField(const Vec4& point, double spin, double charge, double h) {
    Matrix4 field{};
    for (int mu = 0; mu < 4; ++mu) {
        for (int nu = 0; nu < 4; ++nu) {
            field[mu][nu] = PotentialDerivative(point, spin, charge, mu, nu, h) -
                            PotentialDerivative(point, spin, charge, nu, mu, h);
        }
    }
    return field;
}

Metric4d MetricAt(KerrSchildFamily& metric, const Vec4& point) {
    Metric4d value;
    Tensor<Dual<double>, 4, 4, 4> derivatives;
    metric.Evaluate(point, value, derivatives);
    return value;
}

Matrix4 RaisedMaxwellField(KerrSchildFamily& metric, const KerrSchildParams& parameters,
                           const Vec4& point, double h) {
    const Matrix4 covariant = MaxwellField(point, parameters.a, parameters.Q, h);
    const Metric4d inverse = TensorOps::Inverse(MetricAt(metric, point));
    Matrix4 raised{};
    for (int mu = 0; mu < 4; ++mu) {
        for (int nu = 0; nu < 4; ++nu) {
            for (int alpha = 0; alpha < 4; ++alpha) {
                for (int beta = 0; beta < 4; ++beta) {
                    raised[mu][nu] +=
                        inverse(mu, alpha).real * inverse(nu, beta).real * covariant[alpha][beta];
                }
            }
        }
    }
    return raised;
}

Matrix4 EightPiMaxwellStress(KerrSchildFamily& metric, const KerrSchildParams& parameters,
                             const Vec4& point, double h) {
    const Metric4d covariant_metric = MetricAt(metric, point);
    const Metric4d inverse_metric = TensorOps::Inverse(covariant_metric);
    const Matrix4 field = MaxwellField(point, parameters.a, parameters.Q, h);
    const Matrix4 raised = RaisedMaxwellField(metric, parameters, point, h);

    double invariant = 0.0;
    for (int alpha = 0; alpha < 4; ++alpha) {
        for (int beta = 0; beta < 4; ++beta) {
            invariant += field[alpha][beta] * raised[alpha][beta];
        }
    }

    Matrix4 stress{};
    for (int mu = 0; mu < 4; ++mu) {
        for (int nu = 0; nu < 4; ++nu) {
            double contraction = 0.0;
            for (int alpha = 0; alpha < 4; ++alpha) {
                for (int beta = 0; beta < 4; ++beta) {
                    contraction +=
                        field[mu][alpha] * inverse_metric(alpha, beta).real * field[nu][beta];
                }
            }
            stress[mu][nu] = 2.0 * contraction - 0.5 * covariant_metric(mu, nu).real * invariant;
        }
    }
    return stress;
}

struct FieldEquationCase {
    std::string_view name;
    KerrSchildParams parameters;
    double radius;
    double theta;
    double phi;
    bool require_extremal = false;
    bool require_horizonless = false;
};

const std::array<FieldEquationCase, 13>& FieldEquationCases() {
    static const std::array<FieldEquationCase, 13> cases = {{
        {"Minkowski", KerrSchildParams::Minkowski(), 3.0, 0.8, 0.4},
        {"Schwarzschild", KerrSchildParams::Schwarzschild(1.0), 5.0, 0.9, -0.3},
        {"Kerr prograde", KerrSchildParams::Kerr(1.0, 0.7), 5.0, 0.7, 0.6},
        {"Kerr retrograde", KerrSchildParams::Kerr(1.0, -0.7), 5.0, 1.1, -0.4},
        {"Reissner-Nordstrom +Q", KerrSchildParams::ReissnerNordstrom(1.0, 0.6), 5.0, 0.8, 0.5},
        {"Reissner-Nordstrom -Q", KerrSchildParams::ReissnerNordstrom(1.0, -0.6), 4.0, 1.2, -0.7},
        {"Kerr-Newman +a,+Q", KerrSchildParams::KerrNewman(1.0, 0.5, 0.4), 5.0, 0.9, 0.3},
        {"Kerr-Newman -a,-Q", KerrSchildParams::KerrNewman(1.0, -0.5, -0.4), 4.0, 1.0, -0.5},
        {"extremal Kerr-Newman", KerrSchildParams::KerrNewman(1.0, 0.8, 0.6), 4.0, 0.7, -0.6, true,
         false},
        {"horizonless Kerr-Newman", KerrSchildParams::KerrNewman(1.0, 0.9, -0.7), 4.0, 1.2, 0.2,
         false, true},
        {"scaled Kerr-Newman", KerrSchildParams::KerrNewman(0.25, 0.15, 0.1), 1.5, 0.6, 0.8},
        {"de Sitter", KerrSchildParams::DeSitter(0.02), 2.0, 0.8, -0.2},
        {"Kottler", KerrSchildParams{1.0, 0.0, 0.0, 0.01}, 5.0, 1.0, 0.4},
    }};
    return cases;
}

}  // namespace

TEST(KerrSchildFieldEquations, LiveCartesianFamilySatisfiesEinsteinMaxwell) {
    for (const FieldEquationCase& sample : FieldEquationCases()) {
        KerrSchildFamily metric(sample.parameters);
        if (sample.require_extremal) {
            EXPECT_TRUE(metric.HasHorizon()) << sample.name;
            EXPECT_NEAR(metric.ExtremalityParameter(), 1.0, 2.0e-15) << sample.name;
        }
        if (sample.require_horizonless) {
            EXPECT_FALSE(metric.HasHorizon()) << sample.name;
        }
        const Vec4 point =
            KerrSchildPoint(sample.radius, sample.theta, sample.phi, sample.parameters.a);
        const double curvature_step = 2.0e-4 * sample.radius;
        const double field_step = 5.0e-5 * sample.radius;
        const sirius::test_support::RicciSample ricci =
            sirius::test_support::RicciFromConnectionFiniteDifference(metric, point,
                                                                      curvature_step);
        const Metric4d covariant_metric = MetricAt(metric, point);
        const Matrix4 matter = EightPiMaxwellStress(metric, sample.parameters, point, field_step);

        double equation_scale = std::abs(sample.parameters.Lambda);
        double matter_scale = 0.0;
        double cosmological_scale = 0.0;
        for (int mu = 0; mu < 4; ++mu) {
            for (int nu = 0; nu < 4; ++nu) {
                const double left = ricci.covariant[mu][nu] -
                                    0.5 * covariant_metric(mu, nu).real * ricci.scalar +
                                    sample.parameters.Lambda * covariant_metric(mu, nu).real;
                equation_scale =
                    std::max({equation_scale, std::abs(left), std::abs(matter[mu][nu])});
                matter_scale = std::max(matter_scale, std::abs(matter[mu][nu]));
                cosmological_scale =
                    std::max(cosmological_scale,
                             std::abs(sample.parameters.Lambda * covariant_metric(mu, nu).real));
            }
        }
        const double tolerance =
            std::max(2.0e-8 / (sample.radius * sample.radius), 3.0e-4 * equation_scale);
        EXPECT_NEAR(ricci.scalar, 4.0 * sample.parameters.Lambda, tolerance) << sample.name;
        for (int mu = 0; mu < 4; ++mu) {
            for (int nu = 0; nu < 4; ++nu) {
                const double left = ricci.covariant[mu][nu] -
                                    0.5 * covariant_metric(mu, nu).real * ricci.scalar +
                                    sample.parameters.Lambda * covariant_metric(mu, nu).real;
                EXPECT_NEAR(left, matter[mu][nu], tolerance)
                    << sample.name << " component (" << mu << "," << nu << ")";
            }
        }
        if (sample.parameters.Q != 0.0) {
            EXPECT_GT(matter_scale, 10.0 * tolerance)
                << sample.name << " would not reject omission of the Maxwell source";
        }
        if (sample.parameters.Lambda != 0.0) {
            EXPECT_GT(cosmological_scale, 10.0 * tolerance)
                << sample.name << " would not reject omission of the cosmological source";
        }
    }
}

TEST(KerrSchildFieldEquations, KerrNewmanPotentialIsSourceFreeOutsideTheRing) {
    for (const FieldEquationCase& sample : FieldEquationCases()) {
        if (sample.parameters.Q == 0.0) continue;
        KerrSchildFamily metric(sample.parameters);
        const Vec4 point =
            KerrSchildPoint(sample.radius, sample.theta, sample.phi, sample.parameters.a);
        const double determinant = TensorOps::Determinant(MetricAt(metric, point));
        ASSERT_NEAR(determinant, -1.0, 2.0e-12) << sample.name;
        const double volume_density = std::sqrt(-determinant);
        const double derivative_step = 2.0e-4 * sample.radius;
        const double field_step = 4.0e-5 * sample.radius;
        for (int nu = 0; nu < 4; ++nu) {
            double divergence = 0.0;
            for (int mu = 1; mu < 4; ++mu) {
                Vec4 plus_one = point;
                Vec4 minus_one = point;
                Vec4 plus_two = point;
                Vec4 minus_two = point;
                plus_one(mu) += derivative_step;
                minus_one(mu) -= derivative_step;
                plus_two(mu) += 2.0 * derivative_step;
                minus_two(mu) -= 2.0 * derivative_step;
                const Matrix4 fp1 =
                    RaisedMaxwellField(metric, sample.parameters, plus_one, field_step);
                const Matrix4 fm1 =
                    RaisedMaxwellField(metric, sample.parameters, minus_one, field_step);
                const Matrix4 fp2 =
                    RaisedMaxwellField(metric, sample.parameters, plus_two, field_step);
                const Matrix4 fm2 =
                    RaisedMaxwellField(metric, sample.parameters, minus_two, field_step);
                const double dp2 = std::sqrt(-TensorOps::Determinant(MetricAt(metric, plus_two)));
                const double dp1 = std::sqrt(-TensorOps::Determinant(MetricAt(metric, plus_one)));
                const double dm1 = std::sqrt(-TensorOps::Determinant(MetricAt(metric, minus_one)));
                const double dm2 = std::sqrt(-TensorOps::Determinant(MetricAt(metric, minus_two)));
                divergence += (-dp2 * fp2[mu][nu] + 8.0 * dp1 * fp1[mu][nu] -
                               8.0 * dm1 * fm1[mu][nu] + dm2 * fm2[mu][nu]) /
                              (12.0 * derivative_step * volume_density);
            }
            const double scale =
                std::abs(sample.parameters.Q) / (sample.radius * sample.radius * sample.radius);
            EXPECT_NEAR(divergence, 0.0, std::max(2.0e-8, 2.0e-4 * scale))
                << sample.name << " divergence component " << nu;
        }
    }
}

}  // namespace sirius::test
