#pragma once

// Tensor algebra for general relativity: a dimension-parameterised Tensor with
// specialisations for vectors, matrices, and rank-3 Christoffel stacks, plus
// TensorOps for the metric operations the integrator needs. Ported from
// MTTN001A.h.
//
// Governing relations:
//   Christoffel:   Gamma^mu_nu_rho = 1/2 g^mu_sigma (d_rho g_sigma_nu
//                                     + d_nu g_sigma_rho - d_sigma g_nu_rho)
//   Null vector:   solve g_mu_nu k^mu k^nu = 0 for k^0
//   Inner product: <u, v>_g = g_mu_nu u^mu v^nu

#include "sirius/core/dual_number.h"

#include <array>
#include <cmath>
#include <map>

namespace sirius::core {

class TensorOps;

// Generic rank-3 tensor; Dim2 = Dim3 = 1 recovers a vector, Dim3 = 1 a matrix.
template <typename T, int Dim1, int Dim2 = 1, int Dim3 = 1>
class Tensor {
  public:
    std::array<std::array<std::array<T, Dim3>, Dim2>, Dim1> data;

    Tensor() { Zero(); }

    T& operator()(int d1, int d2 = 0, int d3 = 0) { return data[d1][d2][d3]; }
    const T& operator()(int d1, int d2 = 0, int d3 = 0) const { return data[d1][d2][d3]; }

    void Zero() {
        for (int i = 0; i < Dim1; ++i) {
            for (int j = 0; j < Dim2; ++j) {
                for (int k = 0; k < Dim3; ++k) {
                    data[i][j][k] = T(0.0);
                }
            }
        }
    }
};

// Matrix specialisation.
template <typename T, int Dim1, int Dim2>
class Tensor<T, Dim1, Dim2, 1> {
  public:
    std::array<std::array<T, Dim2>, Dim1> data;

    Tensor() { Zero(); }

    T& operator()(int d1, int d2 = 0) { return data[d1][d2]; }
    const T& operator()(int d1, int d2 = 0) const { return data[d1][d2]; }

    void Zero() {
        for (int i = 0; i < Dim1; ++i) {
            for (int j = 0; j < Dim2; ++j) {
                data[i][j] = T(0.0);
            }
        }
    }
};

// Vector specialisation, with the arithmetic a 4-vector needs.
template <typename T, int Dim1>
class Tensor<T, Dim1, 1, 1> {
  public:
    std::array<T, Dim1> data;

    Tensor() { Zero(); }

    T& operator()(int d1) { return data[d1]; }
    const T& operator()(int d1) const { return data[d1]; }

    void Zero() {
        for (int i = 0; i < Dim1; ++i) {
            data[i] = T(0.0);
        }
    }

    Tensor& operator+=(const Tensor& other) {
        for (int i = 0; i < Dim1; ++i) data[i] += other.data[i];
        return *this;
    }

    Tensor& operator-=(const Tensor& other) {
        for (int i = 0; i < Dim1; ++i) data[i] -= other.data[i];
        return *this;
    }

    Tensor& operator*=(double scalar) {
        for (int i = 0; i < Dim1; ++i) data[i] *= scalar;
        return *this;
    }

    Tensor& operator/=(double scalar) {
        for (int i = 0; i < Dim1; ++i) data[i] /= scalar;
        return *this;
    }

    Tensor operator-() const {
        Tensor result;
        for (int i = 0; i < Dim1; ++i) result.data[i] = -data[i];
        return result;
    }

    double Length() const {
        double sum_sq = 0.0;
        for (int i = 0; i < Dim1; ++i) sum_sq += data[i] * data[i];
        return std::sqrt(sum_sq);
    }
};

template <typename T, int Dim1>
Tensor<T, Dim1, 1, 1> operator+(Tensor<T, Dim1, 1, 1> lhs, const Tensor<T, Dim1, 1, 1>& rhs) {
    return lhs += rhs;
}

template <typename T, int Dim1>
Tensor<T, Dim1, 1, 1> operator-(Tensor<T, Dim1, 1, 1> lhs, const Tensor<T, Dim1, 1, 1>& rhs) {
    return lhs -= rhs;
}

template <typename T, int Dim1>
Tensor<T, Dim1, 1, 1> operator*(Tensor<T, Dim1, 1, 1> lhs, double rhs) {
    return lhs *= rhs;
}

template <typename T, int Dim1>
Tensor<T, Dim1, 1, 1> operator*(double lhs, Tensor<T, Dim1, 1, 1> rhs) {
    return rhs *= lhs;
}

template <typename T, int Dim1>
Tensor<T, Dim1, 1, 1> operator/(Tensor<T, Dim1, 1, 1> lhs, double rhs) {
    return lhs /= rhs;
}

// --- General-relativity type aliases -----------------------------------------

// Contravariant spacetime 4-vector.
using Vec4 = Tensor<double, 4>;

// 4x4 metric tensor carrying dual numbers for autodiff compatibility.
using Metric4d = Tensor<Dual<double>, 4, 4>;

// Christoffel symbols Gamma^mu_nu_rho.
struct ChristoffelSymbols {
    Tensor<Dual<double>, 4, 4, 4> gamma;

    ChristoffelSymbols() { gamma.Zero(); }
};

// --- Metric and geodesic operations ------------------------------------------

// Static tensor operations in curved spacetime.
class TensorOps {
  public:
    // Inverse metric by the full Cramer rule.
    // Precondition (documented, deliberately not enforced): g is non-degenerate.
    // A degenerate g yields non-finite entries here, which the integrator's
    // invalid-state guard converts into ray termination; a contract that halted
    // would defeat that designed NaN propagation, and no flat-space stand-in is
    // fabricated.
    static Metric4d Inverse(const Metric4d& g);

    // Determinant by the full 4x4 expansion.
    static double Determinant(const Metric4d& g);

    // Christoffel symbols from the metric and its derivatives dg[sigma,mu,nu] =
    // d_sigma g_mu_nu.
    static ChristoffelSymbols Christoffel(const Metric4d& g,
                                          const Tensor<Dual<double>, 4, 4, 4>& dg);

    // Geodesic acceleration a^mu = -Gamma^mu_nu_rho v^nu v^rho.
    static Vec4 GeodesicAcceleration(const Vec4& velocity, const ChristoffelSymbols& gamma);

    // Geodesic acceleration straight from metric derivatives, bypassing explicit
    // Christoffel construction:
    //   a^mu = -(1/2) g^mu_sigma (d_nu g_sigma_rho + d_rho g_sigma_nu
    //                             - d_sigma g_nu_rho) v^nu v^rho.
    // The caller supplies g_inv so a family with a closed-form inverse can pass
    // the exact one; Inverse is the generic fallback.
    static Vec4 GeodesicAccelerationDirect(const Vec4& velocity, const Metric4d& g_inv,
                                           const Tensor<Dual<double>, 4, 4, 4>& dg);

    // Lower an index: v_mu = g_mu_nu v^nu.
    static Vec4 LowerIndex(const Vec4& vector, const Metric4d& g);

    // Raise an index: v^mu = g^mu_nu v_nu.
    static Vec4 RaiseIndex(const Vec4& vector, const Metric4d& g_inv);

    // Inner product g_mu_nu u^mu v^nu.
    static double InnerProduct(const Vec4& u, const Vec4& v, const Metric4d& g);

    // Normalise a null vector to satisfy g_mu_nu k^mu k^nu = 0.
    static Vec4 NormalizeNull(const Vec4& velocity, const Metric4d& g);
};

}  // namespace sirius::core
