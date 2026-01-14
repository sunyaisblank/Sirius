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
// Backward compatibility alias
// =============================================================================
// Legacy code may include <Dual.h> - provide compatibility
#ifndef SIRIUS_NO_LEGACY_ALIASES
// Dual is already named Dual, no alias needed
#endif

#endif // MTDL001A_H
