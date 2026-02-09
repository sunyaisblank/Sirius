// MTTN001A.cpp - Tensor Operations Implementation

#include "MTTN001A.h"
#include <cmath>

namespace Sirius {

Metric4D TensorOps::inverse(const Metric4D& g) {
    Metric4D g_inv{};
    
    double g00 = g(0, 0).real, g11 = g(1, 1).real, g22 = g(2, 2).real, g33 = g(3, 3).real;
    double g01 = g(0, 1).real, g02 = g(0, 2).real, g03 = g(0, 3).real;
    double g12 = g(1, 2).real, g13 = g(1, 3).real, g23 = g(2, 3).real;
    
    // Check if metric is diagonal (common case for most spacetimes)
    if (std::abs(g01) < 1e-10 && std::abs(g02) < 1e-10 && std::abs(g03) < 1e-10 &&
        std::abs(g12) < 1e-10 && std::abs(g13) < 1e-10 && std::abs(g23) < 1e-10) {
        g_inv(0, 0) = Dual<double>(1.0 / g00, 0.0);
        g_inv(1, 1) = Dual<double>(1.0 / g11, 0.0);
        g_inv(2, 2) = Dual<double>(1.0 / g22, 0.0);
        g_inv(3, 3) = Dual<double>(1.0 / g33, 0.0);
    } else {
        // General case with off-diagonal elements
        double det = determinant(g);
        if (std::abs(det) < 1e-15) return g_inv;
        
        // Simplified inverse for metrics with off-diagonal terms
        g_inv(0, 0) = Dual<double>((g(1, 1).real * g(2, 2).real * g(3, 3).real) / det, 0.0);
        g_inv(1, 1) = Dual<double>((g(0, 0).real * g(2, 2).real * g(3, 3).real) / det, 0.0);
        g_inv(2, 2) = Dual<double>((g(0, 0).real * g(1, 1).real * g(3, 3).real) / det, 0.0);
        g_inv(3, 3) = Dual<double>((g(0, 0).real * g(1, 1).real * g(2, 2).real) / det, 0.0);
    }
    
    return g_inv;
}

double TensorOps::determinant(const Metric4D& g) {
    // Fast path for diagonal metrics
    if (std::abs(g(0, 1).real) < 1e-10 && std::abs(g(0, 2).real) < 1e-10 && std::abs(g(0, 3).real) < 1e-10 &&
        std::abs(g(1, 2).real) < 1e-10 && std::abs(g(1, 3).real) < 1e-10 && std::abs(g(2, 3).real) < 1e-10) {
        return g(0, 0).real * g(1, 1).real * g(2, 2).real * g(3, 3).real;
    }
    
    // General 4x4 determinant via cofactor expansion
    double det = g(0, 0).real * (g(1, 1).real * (g(2, 2).real * g(3, 3).real - g(2, 3).real * g(3, 2).real) -
                                  g(1, 2).real * (g(2, 1).real * g(3, 3).real - g(2, 3).real * g(3, 1).real) +
                                  g(1, 3).real * (g(2, 1).real * g(3, 2).real - g(2, 2).real * g(3, 1).real));
    return det;
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

Vec4 TensorOps::geodesicAccelerationDirect(const Vec4& velocity, const Metric4D& g,
                                            const Tensor<Dual<double>, 4, 4, 4>& dg) {
    // Compute acceleration directly without building Christoffel tensor
    // a^μ = -Γ^μ_νρ v^ν v^ρ
    //     = -(1/2) g^μσ (∂_ν g_σρ + ∂_ρ g_σν - ∂_σ g_νρ) v^ν v^ρ
    //
    // Optimization: pre-compute v^ν v^ρ products (symmetric, so 10 unique values)
    // Then contract with metric derivatives

    // Get inverse metric
    Metric4D g_inv = inverse(g);

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
            normalized(0) = static_cast<float>(-C / B);
        } else {
            normalized(0) = 1.0f; // Fallback
        }
    } else {
        double discriminant = B * B - 4.0 * g00 * C;
        
        if (discriminant < 0.0) {
            // No real solution (shouldn't happen for valid spatial directions)
            // Just return 1.0 as fallback
            normalized(0) = 1.0f;
        } else {
            double sqrt_disc = std::sqrt(discriminant);
            
            // Two solutions: (-B + sqrt_disc)/(2*g00) and (-B - sqrt_disc)/(2*g00)
            // For Lorentzian signature with g00 < 0 (timelike), we need k^0 > 0 (future-directed)
            double k0_plus = (-B + sqrt_disc) / (2.0 * g00);
            double k0_minus = (-B - sqrt_disc) / (2.0 * g00);
            
            // Choose the positive, future-directed root
            if (k0_plus > 0.0) {
                normalized(0) = static_cast<float>(k0_plus);
            } else if (k0_minus > 0.0) {
                normalized(0) = static_cast<float>(k0_minus);
            } else {
                // Both negative - take the larger magnitude (less negative)
                normalized(0) = static_cast<float>(std::max(k0_plus, k0_minus));
            }
        }
    }
    
    return normalized;
}

} // namespace Sirius
