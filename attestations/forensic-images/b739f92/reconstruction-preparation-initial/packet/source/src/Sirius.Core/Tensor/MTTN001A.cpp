// MTTN001A.cpp - Tensor Operations Implementation

#include "MTTN001A.h"
#include <cmath>

namespace Sirius {

// -----------------------------------------------------------------------------
// 2x2 minors shared by determinant and inverse (Cramer's rule).
// Naming: A<r1><r2><c1><c2> is the determinant of rows {r1,r2}, columns {c1,c2}
// read from the lower-right 3x3 blocks of m.
// -----------------------------------------------------------------------------
namespace {

struct Minors2x2 {
    double A2323, A1323, A1223, A0323, A0223, A0123;
    double A2313, A1313, A1213, A2312, A1312, A1212;
    double A0313, A0213, A0312, A0212, A0113, A0112;
};

inline void extractReal(const Metric4D& g, double m[4][4]) {
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            m[i][j] = g(i, j).real;
}

inline Minors2x2 computeMinors(const double m[4][4]) {
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

inline double determinantFromMinors(const double m[4][4], const Minors2x2& A) {
    return m[0][0] * (m[1][1] * A.A2323 - m[1][2] * A.A1323 + m[1][3] * A.A1223)
         - m[0][1] * (m[1][0] * A.A2323 - m[1][2] * A.A0323 + m[1][3] * A.A0223)
         + m[0][2] * (m[1][0] * A.A1323 - m[1][1] * A.A0323 + m[1][3] * A.A0123)
         - m[0][3] * (m[1][0] * A.A1223 - m[1][1] * A.A0223 + m[1][2] * A.A0123);
}

} // namespace

Metric4D TensorOps::inverse(const Metric4D& g) {
    // Full Cramer-rule inverse. Only the real parts are inverted; the dual
    // parts of the result are zero by contract, because no consumer
    // differentiates through the inverse (the geodesic formulation contracts
    // g^μν with ∂g_μν, never with ∂g^μν).
    //
    // Precondition: g is non-degenerate. A degenerate metric produces
    // non-finite entries here, which the integrator's invalid-state guard
    // turns into ray termination. No flat-space stand-in is fabricated.
    double m[4][4];
    extractReal(g, m);
    const Minors2x2 A = computeMinors(m);
    const double det = determinantFromMinors(m, A);
    const double invDet = 1.0 / det;

    Metric4D g_inv;
    g_inv(0, 0) = Dual<double>( invDet * (m[1][1] * A.A2323 - m[1][2] * A.A1323 + m[1][3] * A.A1223));
    g_inv(0, 1) = Dual<double>(-invDet * (m[0][1] * A.A2323 - m[0][2] * A.A1323 + m[0][3] * A.A1223));
    g_inv(0, 2) = Dual<double>( invDet * (m[0][1] * A.A2313 - m[0][2] * A.A1313 + m[0][3] * A.A1213));
    g_inv(0, 3) = Dual<double>(-invDet * (m[0][1] * A.A2312 - m[0][2] * A.A1312 + m[0][3] * A.A1212));
    g_inv(1, 0) = Dual<double>(-invDet * (m[1][0] * A.A2323 - m[1][2] * A.A0323 + m[1][3] * A.A0223));
    g_inv(1, 1) = Dual<double>( invDet * (m[0][0] * A.A2323 - m[0][2] * A.A0323 + m[0][3] * A.A0223));
    g_inv(1, 2) = Dual<double>(-invDet * (m[0][0] * A.A2313 - m[0][2] * A.A0313 + m[0][3] * A.A0213));
    g_inv(1, 3) = Dual<double>( invDet * (m[0][0] * A.A2312 - m[0][2] * A.A0312 + m[0][3] * A.A0212));
    g_inv(2, 0) = Dual<double>( invDet * (m[1][0] * A.A1323 - m[1][1] * A.A0323 + m[1][3] * A.A0123));
    g_inv(2, 1) = Dual<double>(-invDet * (m[0][0] * A.A1323 - m[0][1] * A.A0323 + m[0][3] * A.A0123));
    g_inv(2, 2) = Dual<double>( invDet * (m[0][0] * A.A1313 - m[0][1] * A.A0313 + m[0][3] * A.A0113));
    g_inv(2, 3) = Dual<double>(-invDet * (m[0][0] * A.A1312 - m[0][1] * A.A0312 + m[0][3] * A.A0112));
    g_inv(3, 0) = Dual<double>(-invDet * (m[1][0] * A.A1223 - m[1][1] * A.A0223 + m[1][2] * A.A0123));
    g_inv(3, 1) = Dual<double>( invDet * (m[0][0] * A.A1223 - m[0][1] * A.A0223 + m[0][2] * A.A0123));
    g_inv(3, 2) = Dual<double>(-invDet * (m[0][0] * A.A1213 - m[0][1] * A.A0213 + m[0][2] * A.A0113));
    g_inv(3, 3) = Dual<double>( invDet * (m[0][0] * A.A1212 - m[0][1] * A.A0212 + m[0][2] * A.A0112));
    return g_inv;
}

double TensorOps::determinant(const Metric4D& g) {
    double m[4][4];
    extractReal(g, m);
    const Minors2x2 A = computeMinors(m);
    return determinantFromMinors(m, A);
}

ChristoffelSymbols TensorOps::christoffel(const Metric4D& g, const Tensor<Dual<double>, 4, 4, 4>& dg) {
    ChristoffelSymbols gamma;
    Metric4D g_inv = inverse(g);

    // Γ^μ_νρ = ½ g^μσ (∂_ρ g_σν + ∂_ν g_σρ - ∂_σ g_νρ)
    for (int mu = 0; mu < 4; mu++) {
        for (int nu = 0; nu < 4; nu++) {
            for (int rho = 0; rho < 4; rho++) {
                gamma.gamma(mu, nu, rho) = Dual<double>(0.0, 0.0);

                for (int sigma = 0; sigma < 4; sigma++) {
                    Dual<double> term = dg(rho, sigma, nu) + dg(nu, sigma, rho) - dg(sigma, nu, rho);
                    gamma.gamma(mu, nu, rho) = gamma.gamma(mu, nu, rho) + g_inv(mu, sigma) * term * 0.5;
                }
            }
        }
    }

    return gamma;
}

Vec4 TensorOps::geodesicAcceleration(const Vec4& velocity, const ChristoffelSymbols& gamma) {
    Vec4 acceleration;

    // a^μ = -Γ^μ_νρ v^ν v^ρ (geodesic equation)
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

Vec4 TensorOps::geodesicAccelerationDirect(const Vec4& velocity, const Metric4D& g_inv,
                                            const Tensor<Dual<double>, 4, 4, 4>& dg) {
    // Compute acceleration directly without building the Christoffel tensor:
    // a^μ = -(1/2) g^μσ (∂_ν g_σρ + ∂_ρ g_σν - ∂_σ g_νρ) v^ν v^ρ
    //
    // The caller supplies g_inv so that metric families with a closed-form
    // inverse (Kerr-Schild) can pass the exact one; TensorOps::inverse is the
    // generic fallback.

    // Pre-compute velocity outer product vv[ν][ρ] = v^ν v^ρ
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
            if (std::abs(g_inv_mu_sigma) < 1e-15) continue;  // Skip zero entries

            double sum = 0.0;
            for (int nu = 0; nu < 4; nu++) {
                for (int rho = 0; rho < 4; rho++) {
                    // Γ^μ_νρ contribution via σ:
                    // (1/2) g^μσ (dg(ν,σ,ρ) + dg(ρ,σ,ν) - dg(σ,ν,ρ))
                    double term = dg(nu, sigma, rho).real
                                + dg(rho, sigma, nu).real
                                - dg(sigma, nu, rho).real;
                    sum += term * vv[nu][rho];
                }
            }
            accel += g_inv_mu_sigma * sum;
        }

        acceleration(mu) = -0.5 * accel;
    }

    return acceleration;
}

Vec4 TensorOps::lowerIndex(const Vec4& vector, const Metric4D& g) {
    Vec4 lowered;
    // v_μ = g_μν v^ν
    for (int mu = 0; mu < 4; mu++) {
        for (int nu = 0; nu < 4; nu++) {
            lowered(mu) += g(mu, nu).real * vector(nu);
        }
    }
    return lowered;
}

Vec4 TensorOps::raiseIndex(const Vec4& vector, const Metric4D& g_inv) {
    Vec4 raised;
    // v^μ = g^μν v_ν
    for (int mu = 0; mu < 4; mu++) {
        for (int nu = 0; nu < 4; nu++) {
            raised(mu) += g_inv(mu, nu).real * vector(nu);
        }
    }
    return raised;
}

double TensorOps::innerProduct(const Vec4& u, const Vec4& v, const Metric4D& g) {
    double result = 0.0;
    // g_μν u^μ v^ν
    for (int mu = 0; mu < 4; mu++) {
        for (int nu = 0; nu < 4; nu++) {
            result += g(mu, nu).real * u(mu) * v(nu);
        }
    }
    return result;
}

Vec4 TensorOps::normalizeNull(const Vec4& velocity, const Metric4D& g) {
    Vec4 normalized = velocity;

    // Solve for k^0 given spatial components k^i using null condition:
    // g_00 (k^0)² + 2 g_0i k^0 k^i + g_ij k^i k^j = 0
    //
    // This is a quadratic equation: A (k^0)² + B k^0 + C = 0
    // where:
    //   A = g_00
    //   B = 2 * (g_01 k^1 + g_02 k^2 + g_03 k^3)
    //   C = g_ij k^i k^j (spatial part of inner product)

    double g00 = g(0, 0).real;

    // Compute B = 2 * g_0i k^i
    double B = 0.0;
    for (int i = 1; i < 4; i++) {
        B += 2.0 * g(0, i).real * velocity(i);
    }

    // Compute C = g_ij k^i k^j (only spatial components i,j = 1,2,3)
    double C = 0.0;
    for (int i = 1; i < 4; i++) {
        for (int j = 1; j < 4; j++) {
            C += g(i, j).real * velocity(i) * velocity(j);
        }
    }

    // Solve A*(k^0)^2 + B*k^0 + C = 0
    // k^0 = (-B ± sqrt(B^2 - 4AC)) / (2A)

    if (std::abs(g00) < 1e-15) {
        // Degenerate: g_00 ≈ 0 (near horizon), use linear solution
        if (std::abs(B) > 1e-15) {
            normalized(0) = -C / B;
        } else {
            normalized(0) = 1.0; // Fallback
        }
    } else {
        double discriminant = B * B - 4.0 * g00 * C;

        if (discriminant < 0.0) {
            // No real solution (shouldn't happen for valid spatial directions)
            // Just return 1.0 as fallback
            normalized(0) = 1.0;
        } else {
            double sqrt_disc = std::sqrt(discriminant);

            // Two solutions: (-B + sqrt_disc)/(2*g00) and (-B - sqrt_disc)/(2*g00)
            // For Lorentzian signature with g00 < 0 (timelike), we need k^0 > 0 (future-directed)
            double k0_plus = (-B + sqrt_disc) / (2.0 * g00);
            double k0_minus = (-B - sqrt_disc) / (2.0 * g00);

            // Choose the positive, future-directed root
            if (k0_plus > 0.0) {
                normalized(0) = k0_plus;
            } else if (k0_minus > 0.0) {
                normalized(0) = k0_minus;
            } else {
                // Both negative - take the larger magnitude (less negative)
                normalized(0) = std::max(k0_plus, k0_minus);
            }
        }
    }

    return normalized;
}

} // namespace Sirius
