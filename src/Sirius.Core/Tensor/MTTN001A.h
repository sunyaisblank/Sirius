// MTTN001A.h - Tensor Algebra for General Relativity
//
// Tensor class with specializations for vectors (1D), matrices (2D), and
// 3-tensors (Christoffel symbols). TensorOps provides metric operations.
//
// Key formulas:
//   Christoffel: Γᵘₘₙ = ½ gᵘᵅ (∂ₘ gₐₙ + ∂ₙ gₐₘ - ∂ₐ gₘₙ)
//   Null normalization: solve gₘₙ kᵘ kᵛ = 0 for k⁰
//   Inner product: ⟨u, v⟩_g = gₘₙ uᵘ vᵛ
//
// Tests: TSMT002A.cpp (tensors), TSPH002A.cpp (Christoffel)

#pragma once

#include "MTDL001A.h"
#include <array>
#include <map>
#include <cmath>

namespace Sirius {

// Forward declaration for TensorOps
class TensorOps;

/// @brief Generic N-dimensional tensor
/// @tparam T Element type
/// @tparam Dim1 First dimension size
/// @tparam Dim2 Second dimension size (1 for vectors)
/// @tparam Dim3 Third dimension size (1 for vectors/matrices)
template<typename T, int Dim1, int Dim2 = 1, int Dim3 = 1>
class Tensor {
public:
    std::array<std::array<std::array<T, Dim3>, Dim2>, Dim1> data;

    Tensor() {
        zero();
    }

    T& operator()(int d1, int d2 = 0, int d3 = 0) { return data[d1][d2][d3]; }
    const T& operator()(int d1, int d2 = 0, int d3 = 0) const { return data[d1][d2][d3]; }

    void zero() {
        for (int i = 0; i < Dim1; ++i) {
            for (int j = 0; j < Dim2; ++j) {
                for (int k = 0; k < Dim3; ++k) {
                    data[i][j][k] = T(0.0);
                }
            }
        }
    }
};

// =============================================================================
// Specialization for 2D Tensors (Matrices)
// =============================================================================

template<typename T, int Dim1, int Dim2>
class Tensor<T, Dim1, Dim2, 1> {
public:
    std::array<std::array<T, Dim2>, Dim1> data;

    Tensor() {
        zero();
    }

    T& operator()(int d1, int d2 = 0) { return data[d1][d2]; }
    const T& operator()(int d1, int d2 = 0) const { return data[d1][d2]; }

    void zero() {
        for (int i = 0; i < Dim1; ++i) {
            for (int j = 0; j < Dim2; ++j) {
                data[i][j] = T(0.0);
            }
        }
    }
};

// =============================================================================
// Specialization for 1D Tensors (Vectors)
// =============================================================================

template<typename T, int Dim1>
class Tensor<T, Dim1, 1, 1> {
public:
    std::array<T, Dim1> data;

    Tensor() {
        zero();
    }

    T& operator()(int d1) { return data[d1]; }
    const T& operator()(int d1) const { return data[d1]; }

    void zero() {
        for (int i = 0; i < Dim1; ++i) {
            data[i] = T(0.0);
        }
    }

    // Vector operations
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

    // Unary minus
    Tensor operator-() const {
        Tensor result;
        for (int i = 0; i < Dim1; ++i) result.data[i] = -data[i];
        return result;
    }

    // Length (magnitude)
    double length() const {
        double sum_sq = 0.0;
        for (int i = 0; i < Dim1; ++i) sum_sq += data[i] * data[i];
        return std::sqrt(sum_sq);
    }
};

// =============================================================================
// Global operators for Tensor (vector)
// =============================================================================

template<typename T, int Dim1>
Tensor<T, Dim1, 1, 1> operator+(Tensor<T, Dim1, 1, 1> lhs, const Tensor<T, Dim1, 1, 1>& rhs) { return lhs += rhs; }

template<typename T, int Dim1>
Tensor<T, Dim1, 1, 1> operator-(Tensor<T, Dim1, 1, 1> lhs, const Tensor<T, Dim1, 1, 1>& rhs) { return lhs -= rhs; }

template<typename T, int Dim1>
Tensor<T, Dim1, 1, 1> operator*(Tensor<T, Dim1, 1, 1> lhs, double rhs) { return lhs *= rhs; }

template<typename T, int Dim1>
Tensor<T, Dim1, 1, 1> operator*(double lhs, Tensor<T, Dim1, 1, 1> rhs) { return rhs *= lhs; }

template<typename T, int Dim1>
Tensor<T, Dim1, 1, 1> operator/(Tensor<T, Dim1, 1, 1> lhs, double rhs) { return lhs /= rhs; }

// =============================================================================
// Type Aliases for General Relativity
// =============================================================================

/// @brief 4-vector in spacetime (contravariant)
using Vec4 = Tensor<double, 4>;

/// @brief 4x4 metric tensor with dual number support for autodiff
using Metric4D = Tensor<Dual<double>, 4, 4>;

/// @brief Christoffel symbols Γ^μ_νρ
struct ChristoffelSymbols {
    Tensor<Dual<double>, 4, 4, 4> gamma;
    
    ChristoffelSymbols() {
        gamma.zero();
    }
};

// =============================================================================
// TensorOps - Metric and Geodesic Operations
// =============================================================================

/// @brief Static utility class for tensor operations in curved spacetime
class TensorOps {
public:
    /// @brief Compute inverse of metric tensor
    static Metric4D inverse(const Metric4D& g);
    
    /// @brief Compute determinant of metric tensor
    static double determinant(const Metric4D& g);
    
    /// @brief Compute Christoffel symbols from metric and its derivatives
    /// @param g Metric tensor
    /// @param dg Metric derivatives: dg[σ,μ,ν] = ∂_σ g_μν
    static ChristoffelSymbols christoffel(const Metric4D& g, const Tensor<Dual<double>, 4, 4, 4>& dg);
    
    /// @brief Compute geodesic acceleration: a^μ = -Γ^μ_νρ v^ν v^ρ
    static Vec4 geodesicAcceleration(const Vec4& velocity, const ChristoffelSymbols& gamma);

    /// @brief Compute geodesic acceleration directly from metric derivatives
    /// This bypasses explicit Christoffel symbol construction for better performance.
    /// a^μ = -(1/2) g^μσ (∂_ν g_σρ + ∂_ρ g_σν - ∂_σ g_νρ) v^ν v^ρ
    /// @param velocity 4-velocity v^μ
    /// @param g Metric tensor
    /// @param dg Metric derivatives: dg[σ,μ,ν] = ∂_σ g_μν
    /// @return Geodesic acceleration a^μ
    static Vec4 geodesicAccelerationDirect(const Vec4& velocity, const Metric4D& g,
                                            const Tensor<Dual<double>, 4, 4, 4>& dg);

    /// @brief Lower index: v_μ = g_μν v^ν
    static Vec4 lowerIndex(const Vec4& vector, const Metric4D& g);
    
    /// @brief Raise index: v^μ = g^μν v_ν
    static Vec4 raiseIndex(const Vec4& vector, const Metric4D& g_inv);
    
    /// @brief Compute inner product: g_μν u^μ v^ν
    static double innerProduct(const Vec4& u, const Vec4& v, const Metric4D& g);
    
    /// @brief Normalize null vector to satisfy g_μν k^μ k^ν = 0
    static Vec4 normalizeNull(const Vec4& velocity, const Metric4D& g);
};

} // namespace Sirius
