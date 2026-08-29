// Implementation of the metric operations declared in tensor.h. Ported from
// MTTN001A.cpp; arithmetic order and literals bit-identical.

#include "sirius/core/tensor.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace sirius::core {

// 2x2 minors shared by the determinant and the inverse (Cramer's rule).
// Naming: A<r1><r2><c1><c2> is the determinant of rows {r1,r2}, columns {c1,c2}
// read from the lower-right 3x3 blocks of m.
namespace {

struct Minors2x2 {
    double A2323, A1323, A1223, A0323, A0223, A0123;
    double A2313, A1313, A1213, A2312, A1312, A1212;
    double A0313, A0213, A0312, A0212, A0113, A0112;
};

inline void ExtractReal(const Metric4d& g, double m[4][4]) {
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++) m[i][j] = g(i, j).real;
}

inline Minors2x2 ComputeMinors(const double m[4][4]) {
    Minors2x2 A;
    A.A2323 = m[2][2] * m[3][3] - m[2][3] * m[3][2];
    A.A1323 = m[2][1] * m[3][3] - m[2][3] * m[3][1];
    A.A1223 = m[2][1] * m[3][2] - m[2][2] * m[3][1];
    A.A0323 = m[2][0] * m[3][3] - m[2][3] * m[3][0];
    A.A0223 = m[2][0] * m[3][2] - m[2][2] * m[3][0];
    A.A0123 = m[2][0] * m[3][1] - m[2][1] * m[3][0];
    A.A2313 = m[1][2] * m[3][3] - m[1][3] * m[3][2];
    A.A1313 = m[1][1] * m[3][3] - m[1][3] * m[3][1];
    A.A1213 = m[1][1] * m[3][2] - m[1][2] * m[3][1];
    A.A2312 = m[1][2] * m[2][3] - m[1][3] * m[2][2];
    A.A1312 = m[1][1] * m[2][3] - m[1][3] * m[2][1];
    A.A1212 = m[1][1] * m[2][2] - m[1][2] * m[2][1];
    A.A0313 = m[1][0] * m[3][3] - m[1][3] * m[3][0];
    A.A0213 = m[1][0] * m[3][2] - m[1][2] * m[3][0];
    A.A0312 = m[1][0] * m[2][3] - m[1][3] * m[2][0];
    A.A0212 = m[1][0] * m[2][2] - m[1][2] * m[2][0];
    A.A0113 = m[1][0] * m[3][1] - m[1][1] * m[3][0];
    A.A0112 = m[1][0] * m[2][1] - m[1][1] * m[2][0];
    return A;
}

inline double DeterminantFromMinors(const double m[4][4], const Minors2x2& A) {
    return m[0][0] * (m[1][1] * A.A2323 - m[1][2] * A.A1323 + m[1][3] * A.A1223) -
           m[0][1] * (m[1][0] * A.A2323 - m[1][2] * A.A0323 + m[1][3] * A.A0223) +
           m[0][2] * (m[1][0] * A.A1323 - m[1][1] * A.A0323 + m[1][3] * A.A0123) -
           m[0][3] * (m[1][0] * A.A1223 - m[1][1] * A.A0223 + m[1][2] * A.A0123);
}

}  // namespace

Metric4d TensorOps::Inverse(const Metric4d& g) {
    // Full Cramer-rule inverse. Only the real parts are inverted; the dual parts
    // of the result are zero by contract, because no consumer differentiates
    // through the inverse (the geodesic formulation contracts g^mu_nu with
    // d g_mu_nu, never with d g^mu_nu).
    //
    // Precondition: g is non-degenerate. A degenerate metric produces non-finite
    // entries here, which the integrator's invalid-state guard turns into ray
    // termination. No flat-space stand-in is fabricated.
    double m[4][4];
    ExtractReal(g, m);
    const Minors2x2 A = ComputeMinors(m);
    const double det = DeterminantFromMinors(m, A);
    const double inv_det = 1.0 / det;

    Metric4d g_inv;
    g_inv(0, 0) =
        Dual<double>(inv_det * (m[1][1] * A.A2323 - m[1][2] * A.A1323 + m[1][3] * A.A1223));
    g_inv(0, 1) =
        Dual<double>(-inv_det * (m[0][1] * A.A2323 - m[0][2] * A.A1323 + m[0][3] * A.A1223));
    g_inv(0, 2) =
        Dual<double>(inv_det * (m[0][1] * A.A2313 - m[0][2] * A.A1313 + m[0][3] * A.A1213));
    g_inv(0, 3) =
        Dual<double>(-inv_det * (m[0][1] * A.A2312 - m[0][2] * A.A1312 + m[0][3] * A.A1212));
    g_inv(1, 0) =
        Dual<double>(-inv_det * (m[1][0] * A.A2323 - m[1][2] * A.A0323 + m[1][3] * A.A0223));
    g_inv(1, 1) =
        Dual<double>(inv_det * (m[0][0] * A.A2323 - m[0][2] * A.A0323 + m[0][3] * A.A0223));
    g_inv(1, 2) =
        Dual<double>(-inv_det * (m[0][0] * A.A2313 - m[0][2] * A.A0313 + m[0][3] * A.A0213));
    g_inv(1, 3) =
        Dual<double>(inv_det * (m[0][0] * A.A2312 - m[0][2] * A.A0312 + m[0][3] * A.A0212));
    g_inv(2, 0) =
        Dual<double>(inv_det * (m[1][0] * A.A1323 - m[1][1] * A.A0323 + m[1][3] * A.A0123));
    g_inv(2, 1) =
        Dual<double>(-inv_det * (m[0][0] * A.A1323 - m[0][1] * A.A0323 + m[0][3] * A.A0123));
    g_inv(2, 2) =
        Dual<double>(inv_det * (m[0][0] * A.A1313 - m[0][1] * A.A0313 + m[0][3] * A.A0113));
    g_inv(2, 3) =
        Dual<double>(-inv_det * (m[0][0] * A.A1312 - m[0][1] * A.A0312 + m[0][3] * A.A0112));
    g_inv(3, 0) =
        Dual<double>(-inv_det * (m[1][0] * A.A1223 - m[1][1] * A.A0223 + m[1][2] * A.A0123));
    g_inv(3, 1) =
        Dual<double>(inv_det * (m[0][0] * A.A1223 - m[0][1] * A.A0223 + m[0][2] * A.A0123));
    g_inv(3, 2) =
        Dual<double>(-inv_det * (m[0][0] * A.A1213 - m[0][1] * A.A0213 + m[0][2] * A.A0113));
    g_inv(3, 3) =
        Dual<double>(inv_det * (m[0][0] * A.A1212 - m[0][1] * A.A0212 + m[0][2] * A.A0112));
    return g_inv;
}

double TensorOps::Determinant(const Metric4d& g) {
    double m[4][4];
    ExtractReal(g, m);
    const Minors2x2 A = ComputeMinors(m);
    return DeterminantFromMinors(m, A);
}

ChristoffelSymbols TensorOps::Christoffel(const Metric4d& g,
                                          const Tensor<Dual<double>, 4, 4, 4>& dg) {
    ChristoffelSymbols gamma;
    Metric4d g_inv = Inverse(g);

    // Gamma^mu_nu_rho = 1/2 g^mu_sigma (d_rho g_sigma_nu + d_nu g_sigma_rho
    //                                   - d_sigma g_nu_rho)
    for (int mu = 0; mu < 4; mu++) {
        for (int nu = 0; nu < 4; nu++) {
            for (int rho = 0; rho < 4; rho++) {
                gamma.gamma(mu, nu, rho) = Dual<double>(0.0, 0.0);

                for (int sigma = 0; sigma < 4; sigma++) {
                    Dual<double> term =
                        dg(rho, sigma, nu) + dg(nu, sigma, rho) - dg(sigma, nu, rho);
                    gamma.gamma(mu, nu, rho) =
                        gamma.gamma(mu, nu, rho) + g_inv(mu, sigma) * term * 0.5;
                }
            }
        }
    }

    return gamma;
}

Vec4 TensorOps::GeodesicAcceleration(const Vec4& velocity, const ChristoffelSymbols& gamma) {
    Vec4 acceleration;

    // a^mu = -Gamma^mu_nu_rho v^nu v^rho (geodesic equation).
    for (int mu = 0; mu < 4; mu++) {
        double accel = 0.0;
        for (int nu = 0; nu < 4; nu++) {
            for (int rho = 0; rho < 4; rho++) {
                accel -= gamma.gamma(mu, nu, rho).real * velocity(nu) * velocity(rho);
            }
        }
        acceleration(mu) = accel;
    }

    return acceleration;
}

Vec4 TensorOps::GeodesicAccelerationDirect(const Vec4& velocity, const Metric4d& g_inv,
                                           const Tensor<Dual<double>, 4, 4, 4>& dg) {
    // Acceleration without building the Christoffel tensor:
    //   a^mu = -(1/2) g^mu_sigma (d_nu g_sigma_rho + d_rho g_sigma_nu
    //                             - d_sigma g_nu_rho) v^nu v^rho.
    // The caller supplies g_inv so families with a closed-form inverse
    // (Kerr-Schild) can pass the exact one; Inverse is the generic fallback.

    // Velocity outer product vv[nu][rho] = v^nu v^rho.
    double vv[4][4];
    for (int nu = 0; nu < 4; nu++) {
        for (int rho = 0; rho < 4; rho++) {
            vv[nu][rho] = velocity(nu) * velocity(rho);
        }
    }

    Vec4 acceleration;

    for (int mu = 0; mu < 4; mu++) {
        double accel = 0.0;

        for (int sigma = 0; sigma < 4; sigma++) {
            double g_inv_mu_sigma = g_inv(mu, sigma).real;
            if (std::abs(g_inv_mu_sigma) < 1e-15) continue;  // Skip zero entries.

            double sum = 0.0;
            for (int nu = 0; nu < 4; nu++) {
                for (int rho = 0; rho < 4; rho++) {
                    double term =
                        dg(nu, sigma, rho).real + dg(rho, sigma, nu).real - dg(sigma, nu, rho).real;
                    sum += term * vv[nu][rho];
                }
            }
            accel += g_inv_mu_sigma * sum;
        }

        acceleration(mu) = -0.5 * accel;
    }

    return acceleration;
}

Vec4 TensorOps::LowerIndex(const Vec4& vector, const Metric4d& g) {
    Vec4 lowered;
    // v_mu = g_mu_nu v^nu.
    for (int mu = 0; mu < 4; mu++) {
        for (int nu = 0; nu < 4; nu++) {
            lowered(mu) += g(mu, nu).real * vector(nu);
        }
    }
    return lowered;
}

Vec4 TensorOps::RaiseIndex(const Vec4& vector, const Metric4d& g_inv) {
    Vec4 raised;
    // v^mu = g^mu_nu v_nu.
    for (int mu = 0; mu < 4; mu++) {
        for (int nu = 0; nu < 4; nu++) {
            raised(mu) += g_inv(mu, nu).real * vector(nu);
        }
    }
    return raised;
}

double TensorOps::InnerProduct(const Vec4& u, const Vec4& v, const Metric4d& g) {
    double result = 0.0;
    // g_mu_nu u^mu v^nu.
    for (int mu = 0; mu < 4; mu++) {
        for (int nu = 0; nu < 4; nu++) {
            result += g(mu, nu).real * u(mu) * v(nu);
        }
    }
    return result;
}

Vec4 TensorOps::NormalizeNull(const Vec4& velocity, const Metric4d& g) {
    Vec4 normalized = velocity;

    // Solve for k^0 from spatial components k^i via the null condition
    //   g_00 (k^0)^2 + 2 g_0i k^0 k^i + g_ij k^i k^j = 0,
    // a quadratic A (k^0)^2 + B k^0 + C = 0 with
    //   A = g_00, B = 2 g_0i k^i, C = g_ij k^i k^j.

    double g00 = g(0, 0).real;

    double B = 0.0;
    for (int i = 1; i < 4; i++) {
        B += 2.0 * g(0, i).real * velocity(i);
    }

    double C = 0.0;
    for (int i = 1; i < 4; i++) {
        for (int j = 1; j < 4; j++) {
            C += g(i, j).real * velocity(i) * velocity(j);
        }
    }

    if (std::abs(g00) < 1e-15) {
        // Degenerate g_00 ~ 0: use the exact linear branch when it exists.
        if (std::abs(B) > 1e-15) {
            normalized(0) = -C / B;
        } else {
            normalized(0) = std::numeric_limits<double>::quiet_NaN();
        }
    } else {
        double discriminant = B * B - 4.0 * g00 * C;

        if (discriminant < 0.0) {
            normalized(0) = std::numeric_limits<double>::quiet_NaN();
        } else {
            double sqrt_disc = std::sqrt(discriminant);
            double k0_plus = (-B + sqrt_disc) / (2.0 * g00);
            double k0_minus = (-B - sqrt_disc) / (2.0 * g00);
            // A nonzero input k^0 identifies the null family and must survive a
            // projection (in particular, past-directed render rays must not be
            // replaced by the other physical light cone). A zero input is the
            // legacy explicit-initialisation request and selects the larger
            // coordinate-time root outside an ergoregion.
            normalized(0) =
                std::abs(velocity(0)) > 1.0e-15
                    ? (std::abs(k0_plus - velocity(0)) <= std::abs(k0_minus - velocity(0))
                           ? k0_plus
                           : k0_minus)
                    : std::max(k0_plus, k0_minus);
        }
    }

    return normalized;
}

}  // namespace sirius::core
