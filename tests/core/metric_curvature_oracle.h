#pragma once

// Test-only curvature oracle shared by live Cartesian metric gates.  The
// production metric supplies g and its analytic first derivatives; this oracle
// independently differentiates the resulting Levi-Civita connection with a
// fourth-order centred stencil and assembles the Ricci tensor.  Closed-form
// curvature or stress-energy expressions under test never enter this path.

#include "sirius/core/metrics/metric.h"

namespace sirius::test_support {

struct RicciSample {
    double covariant[4][4] = {};
    double scalar = 0.0;
};

inline sirius::core::ChristoffelSymbols ConnectionAt(sirius::core::IMetric& metric,
                                                     const sirius::core::Vec4& position) {
    sirius::core::Metric4d g;
    sirius::core::Tensor<sirius::core::Dual<double>, 4, 4, 4> dg;
    metric.Evaluate(position, g, dg);
    return sirius::core::TensorOps::Christoffel(g, dg);
}

inline RicciSample RicciFromConnectionFiniteDifference(sirius::core::IMetric& metric,
                                                       const sirius::core::Vec4& position,
                                                       double h) {
    const sirius::core::ChristoffelSymbols gamma = ConnectionAt(metric, position);
    double d_gamma[4][4][4][4] = {};
    for (int derivative = 0; derivative < 4; ++derivative) {
        sirius::core::Vec4 plus_one = position;
        sirius::core::Vec4 minus_one = position;
        sirius::core::Vec4 plus_two = position;
        sirius::core::Vec4 minus_two = position;
        plus_one(derivative) += h;
        minus_one(derivative) -= h;
        plus_two(derivative) += 2.0 * h;
        minus_two(derivative) -= 2.0 * h;
        const sirius::core::ChristoffelSymbols gp1 = ConnectionAt(metric, plus_one);
        const sirius::core::ChristoffelSymbols gm1 = ConnectionAt(metric, minus_one);
        const sirius::core::ChristoffelSymbols gp2 = ConnectionAt(metric, plus_two);
        const sirius::core::ChristoffelSymbols gm2 = ConnectionAt(metric, minus_two);
        for (int upper = 0; upper < 4; ++upper) {
            for (int lower_one = 0; lower_one < 4; ++lower_one) {
                for (int lower_two = 0; lower_two < 4; ++lower_two) {
                    d_gamma[derivative][upper][lower_one][lower_two] =
                        (-gp2.gamma(upper, lower_one, lower_two).real +
                         8.0 * gp1.gamma(upper, lower_one, lower_two).real -
                         8.0 * gm1.gamma(upper, lower_one, lower_two).real +
                         gm2.gamma(upper, lower_one, lower_two).real) /
                        (12.0 * h);
                }
            }
        }
    }

    RicciSample result;
    for (int sigma = 0; sigma < 4; ++sigma) {
        for (int nu = 0; nu < 4; ++nu) {
            for (int rho = 0; rho < 4; ++rho) {
                double component = d_gamma[rho][rho][nu][sigma] - d_gamma[nu][rho][rho][sigma];
                for (int lambda = 0; lambda < 4; ++lambda) {
                    component +=
                        gamma.gamma(rho, rho, lambda).real * gamma.gamma(lambda, nu, sigma).real -
                        gamma.gamma(rho, nu, lambda).real * gamma.gamma(lambda, rho, sigma).real;
                }
                result.covariant[sigma][nu] += component;
            }
        }
    }

    sirius::core::Metric4d g;
    sirius::core::Tensor<sirius::core::Dual<double>, 4, 4, 4> dg;
    metric.Evaluate(position, g, dg);
    sirius::core::Metric4d inverse;
    if (!metric.InverseMetric(position, inverse)) {
        inverse = sirius::core::TensorOps::Inverse(g);
    }
    for (int mu = 0; mu < 4; ++mu) {
        for (int nu = 0; nu < 4; ++nu) {
            result.scalar += inverse(mu, nu).real * result.covariant[mu][nu];
        }
    }
    return result;
}

}  // namespace sirius::test_support
