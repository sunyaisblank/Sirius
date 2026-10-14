// MTDL001A.h - Dual Numbers for Automatic Differentiation
//
// Dual numbers: D = {a + bε : ε² = 0}
// Implements forward-mode autodiff: f(a + bε) = f(a) + b·f'(a)·ε
//
// Key derivatives: sin→cos, sqrt→1/(2√), exp→exp, log→1/x
// Tests: TSMT001A.cpp

#ifndef MTDL001A_H
#define MTDL001A_H

#include <cmath>
#include <limits>

namespace Sirius {

/// @brief Dual number for automatic differentiation
/// @tparam T Underlying numeric type (typically double)
template<typename T>
struct Dual {
    T real;  ///< Real part (function value)
    T dual;  ///< Dual part (derivative)

    Dual(T real = 0, T dual = 0) : real(real), dual(dual) {}

    Dual& operator+=(const Dual& other) {
        real += other.real;
        dual += other.dual;
        return *this;
    }

    Dual& operator-=(const Dual& other) {
        real -= other.real;
        dual -= other.dual;
        return *this;
    }

    Dual& operator*=(const Dual& other) {
        // Product rule: (a + bε)(c + dε) = ac + (ad + bc)ε
        dual = real * other.dual + dual * other.real;
        real *= other.real;
        return *this;
    }

    Dual& operator/=(const Dual& other) {
        // Quotient rule: (a + bε)/(c + dε) = a/c + (bc - ad)/c² ε
        dual = (dual * other.real - real * other.dual) / (other.real * other.real);
        real /= other.real;
        return *this;
    }

    // Scalar multiplication
    Dual& operator*=(T scalar) {
        dual *= scalar;
        real *= scalar;
        return *this;
    }

    // Scalar division
    Dual& operator/=(T scalar) {
        dual /= scalar;
        real /= scalar;
        return *this;
    }

    // Unary minus
    Dual operator-() const {
        return Dual(-real, -dual);
    }
};

// Binary operators
template<typename T>
Dual<T> operator+(Dual<T> lhs, const Dual<T>& rhs) { return lhs += rhs; }

template<typename T>
Dual<T> operator-(Dual<T> lhs, const Dual<T>& rhs) { return lhs -= rhs; }

template<typename T>
Dual<T> operator*(Dual<T> lhs, const Dual<T>& rhs) { return lhs *= rhs; }

// Scalar multiplication (Dual * scalar)
template<typename T>
Dual<T> operator*(Dual<T> lhs, T rhs) { return lhs *= rhs; }

// Scalar multiplication (scalar * Dual)
template<typename T>
Dual<T> operator*(T lhs, Dual<T> rhs) { return rhs *= lhs; }

template<typename T>
Dual<T> operator/(Dual<T> lhs, const Dual<T>& rhs) { return lhs /= rhs; }

// Scalar division (Dual / scalar)
template<typename T>
Dual<T> operator/(Dual<T> lhs, T rhs) { return lhs /= rhs; }

// =============================================================================
// Transcendental functions with automatic differentiation
// =============================================================================

/// @brief sin(a + bε) = sin(a) + b·cos(a)·ε
template<typename T>
Dual<T> sin(const Dual<T>& d) {
    return Dual<T>(std::sin(d.real), d.dual * std::cos(d.real));
}

/// @brief cos(a + bε) = cos(a) - b·sin(a)·ε
template<typename T>
Dual<T> cos(const Dual<T>& d) {
    return Dual<T>(std::cos(d.real), -d.dual * std::sin(d.real));
}

/// @brief sqrt(a + bε) = √a + b/(2√a)·ε
template<typename T>
Dual<T> sqrt(const Dual<T>& d) {
    T sqrt_real = std::sqrt(d.real);
    return Dual<T>(sqrt_real, d.dual / (T(2.0) * sqrt_real));
}

/// @brief exp(a + bε) = exp(a) + b·exp(a)·ε
template<typename T>
Dual<T> exp(const Dual<T>& d) {
    T exp_real = std::exp(d.real);
    return Dual<T>(exp_real, d.dual * exp_real);
}

/// @brief log(a + bε) = log(a) + b/a·ε
template<typename T>
Dual<T> log(const Dual<T>& d) {
    return Dual<T>(std::log(d.real), d.dual / d.real);
}

/// @brief pow(a + bε, n) = a^n + b·n·a^(n-1)·ε
template<typename T>
Dual<T> pow(const Dual<T>& base, T exponent) {
    T pow_real = std::pow(base.real, exponent);
    return Dual<T>(pow_real, base.dual * exponent * std::pow(base.real, exponent - T(1)));
}

/// @brief atan2(y, x) with automatic differentiation
template<typename T>
Dual<T> atan2(const Dual<T>& y, const Dual<T>& x) {
    T denom = x.real * x.real + y.real * y.real;
    return Dual<T>(
        std::atan2(y.real, x.real),
        (x.real * y.dual - y.real * x.dual) / denom
    );
}

/// @brief fabs(a + bε) = |a| + sign(a)·b·ε
template<typename T>
Dual<T> fabs(const Dual<T>& d) {
    T sign = (d.real >= T(0)) ? T(1) : T(-1);
    return Dual<T>(std::fabs(d.real), sign * d.dual);
}

// =============================================================================
// Extended Transcendental Functions (Phase 2)
// =============================================================================

/// @brief tan(a + bε) = tan(a) + b·sec²(a)·ε
/// Note: sec²(a) = 1 + tan²(a) = 1/cos²(a)
template<typename T>
Dual<T> tan(const Dual<T>& d) {
    T tan_real = std::tan(d.real);
    T cos_real = std::cos(d.real);
    T sec2 = T(1) / (cos_real * cos_real);  // sec²(x) = 1/cos²(x)
    return Dual<T>(tan_real, d.dual * sec2);
}

/// @brief sinh(a + bε) = sinh(a) + b·cosh(a)·ε
template<typename T>
Dual<T> sinh(const Dual<T>& d) {
    return Dual<T>(std::sinh(d.real), d.dual * std::cosh(d.real));
}

/// @brief cosh(a + bε) = cosh(a) + b·sinh(a)·ε
template<typename T>
Dual<T> cosh(const Dual<T>& d) {
    return Dual<T>(std::cosh(d.real), d.dual * std::sinh(d.real));
}

/// @brief tanh(a + bε) = tanh(a) + b·sech²(a)·ε
/// Note: sech²(a) = 1 - tanh²(a) = 1/cosh²(a)
template<typename T>
Dual<T> tanh(const Dual<T>& d) {
    T tanh_real = std::tanh(d.real);
    T sech2 = T(1) - tanh_real * tanh_real;  // sech²(x) = 1 - tanh²(x)
    return Dual<T>(tanh_real, d.dual * sech2);
}

/// @brief asin(a + bε) = asin(a) + b/√(1-a²)·ε
/// Domain: |a| < 1
template<typename T>
Dual<T> asin(const Dual<T>& d) {
    T one_minus_x2 = T(1) - d.real * d.real;
    // Protect against domain boundary
    if (one_minus_x2 <= T(0)) {
        one_minus_x2 = std::numeric_limits<T>::epsilon();
    }
    T deriv = T(1) / std::sqrt(one_minus_x2);
    return Dual<T>(std::asin(d.real), d.dual * deriv);
}

/// @brief acos(a + bε) = acos(a) - b/√(1-a²)·ε
/// Domain: |a| < 1
template<typename T>
Dual<T> acos(const Dual<T>& d) {
    T one_minus_x2 = T(1) - d.real * d.real;
    // Protect against domain boundary
    if (one_minus_x2 <= T(0)) {
        one_minus_x2 = std::numeric_limits<T>::epsilon();
    }
    T deriv = T(-1) / std::sqrt(one_minus_x2);
    return Dual<T>(std::acos(d.real), d.dual * deriv);
}

/// @brief atan(a + bε) = atan(a) + b/(1+a²)·ε
template<typename T>
Dual<T> atan(const Dual<T>& d) {
    T deriv = T(1) / (T(1) + d.real * d.real);
    return Dual<T>(std::atan(d.real), d.dual * deriv);
}

/// @brief hypot(x, y) = √(x² + y²) with derivatives
/// Safe computation avoiding overflow: d(hypot)/dx = x/hypot, d(hypot)/dy = y/hypot
template<typename T>
Dual<T> hypot(const Dual<T>& x, const Dual<T>& y) {
    T hypot_real = std::hypot(x.real, y.real);
    if (hypot_real < std::numeric_limits<T>::epsilon()) {
        return Dual<T>(T(0), T(0));
    }
    // d(hypot)/dx = x/hypot, d(hypot)/dy = y/hypot
    T dual = (x.real * x.dual + y.real * y.dual) / hypot_real;
    return Dual<T>(hypot_real, dual);
}

/// @brief copysign(x, y) = |x| × sign(y)
/// Derivative follows x if signs match, negates otherwise
template<typename T>
Dual<T> copysign(const Dual<T>& x, const Dual<T>& y) {
    T mag_sign = (y.real >= T(0)) ? T(1) : T(-1);
    T x_sign = (x.real >= T(0)) ? T(1) : T(-1);
    T result_sign = mag_sign;

    // If x was negative and y positive (or vice versa), negate the derivative
    T deriv_factor = x_sign * mag_sign;

    return Dual<T>(std::copysign(x.real, y.real), deriv_factor * x.dual);
}

/// @brief abs(a + bε) = |a| + sign(a)·b·ε
/// Alias for fabs for consistency with std::abs
template<typename T>
Dual<T> abs(const Dual<T>& d) {
    return fabs(d);
}

/// @brief max(x, y) with derivative tracking
template<typename T>
Dual<T> max(const Dual<T>& x, const Dual<T>& y) {
    if (x.real >= y.real) {
        return x;
    }
    return y;
}

/// @brief min(x, y) with derivative tracking
template<typename T>
Dual<T> min(const Dual<T>& x, const Dual<T>& y) {
    if (x.real <= y.real) {
        return x;
    }
    return y;
}

/// @brief clamp(x, lo, hi) with derivative tracking
template<typename T>
Dual<T> clamp(const Dual<T>& x, const Dual<T>& lo, const Dual<T>& hi) {
    if (x.real < lo.real) {
        return Dual<T>(lo.real, T(0));  // At boundary, derivative is 0
    }
    if (x.real > hi.real) {
        return Dual<T>(hi.real, T(0));  // At boundary, derivative is 0
    }
    return x;  // In range, derivative passes through
}

/// @brief clamp with scalar bounds
template<typename T>
Dual<T> clamp(const Dual<T>& x, T lo, T hi) {
    if (x.real < lo) {
        return Dual<T>(lo, T(0));
    }
    if (x.real > hi) {
        return Dual<T>(hi, T(0));
    }
    return x;
}

// =============================================================================
// Comparison Operators
// =============================================================================

template<typename T>
bool operator<(const Dual<T>& lhs, const Dual<T>& rhs) {
    return lhs.real < rhs.real;
}

template<typename T>
bool operator>(const Dual<T>& lhs, const Dual<T>& rhs) {
    return lhs.real > rhs.real;
}

template<typename T>
bool operator<=(const Dual<T>& lhs, const Dual<T>& rhs) {
    return lhs.real <= rhs.real;
}

template<typename T>
bool operator>=(const Dual<T>& lhs, const Dual<T>& rhs) {
    return lhs.real >= rhs.real;
}

template<typename T>
bool operator==(const Dual<T>& lhs, const Dual<T>& rhs) {
    return lhs.real == rhs.real;
}

template<typename T>
bool operator!=(const Dual<T>& lhs, const Dual<T>& rhs) {
    return lhs.real != rhs.real;
}

// Scalar comparisons
template<typename T>
bool operator<(const Dual<T>& lhs, T rhs) { return lhs.real < rhs; }

template<typename T>
bool operator>(const Dual<T>& lhs, T rhs) { return lhs.real > rhs; }

template<typename T>
bool operator<(T lhs, const Dual<T>& rhs) { return lhs < rhs.real; }

template<typename T>
bool operator>(T lhs, const Dual<T>& rhs) { return lhs > rhs.real; }

// =============================================================================
// Utility Functions
// =============================================================================

/// @brief Check if dual number has finite values
template<typename T>
bool isfinite(const Dual<T>& d) {
    return std::isfinite(d.real) && std::isfinite(d.dual);
}

/// @brief Check if dual number is NaN
template<typename T>
bool isnan(const Dual<T>& d) {
    return std::isnan(d.real) || std::isnan(d.dual);
}

/// @brief Check if dual number is infinite
template<typename T>
bool isinf(const Dual<T>& d) {
    return std::isinf(d.real) || std::isinf(d.dual);
}

} // namespace Sirius

#endif // MTDL001A_H
