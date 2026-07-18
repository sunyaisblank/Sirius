#pragma once

// Forward-mode automatic differentiation over a dual number D = a + b epsilon
// with epsilon^2 = 0, so f(a + b eps) = f(a) + b f'(a) eps. The transcendental
// overloads carry the analytic derivative of each function in the dual part.
// Ported from MTDL001A.h; algorithms and literals bit-identical.
//
// The free math functions mirror the std spelling deliberately: they are ADL
// customisation points so generic code that writes sin(x) resolves to this
// overload for Dual and to std::sin for a plain scalar.

#include <cmath>
#include <limits>

namespace sirius::core {

// Dual number for automatic differentiation; T is the underlying scalar.
template <typename T>
struct Dual {
    T real;  // Function value.
    T dual;  // Derivative.

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
        // Product rule: (a + b eps)(c + d eps) = ac + (ad + bc) eps.
        dual = real * other.dual + dual * other.real;
        real *= other.real;
        return *this;
    }

    Dual& operator/=(const Dual& other) {
        // Quotient rule: (a + b eps)/(c + d eps) = a/c + (bc - ad)/c^2 eps.
        dual = (dual * other.real - real * other.dual) / (other.real * other.real);
        real /= other.real;
        return *this;
    }

    Dual& operator*=(T scalar) {
        dual *= scalar;
        real *= scalar;
        return *this;
    }

    Dual& operator/=(T scalar) {
        dual /= scalar;
        real /= scalar;
        return *this;
    }

    Dual operator-() const { return Dual(-real, -dual); }
};

template <typename T>
Dual<T> operator+(Dual<T> lhs, const Dual<T>& rhs) {
    return lhs += rhs;
}

template <typename T>
Dual<T> operator-(Dual<T> lhs, const Dual<T>& rhs) {
    return lhs -= rhs;
}

template <typename T>
Dual<T> operator*(Dual<T> lhs, const Dual<T>& rhs) {
    return lhs *= rhs;
}

template <typename T>
Dual<T> operator*(Dual<T> lhs, T rhs) {
    return lhs *= rhs;
}

template <typename T>
Dual<T> operator*(T lhs, Dual<T> rhs) {
    return rhs *= lhs;
}

template <typename T>
Dual<T> operator/(Dual<T> lhs, const Dual<T>& rhs) {
    return lhs /= rhs;
}

template <typename T>
Dual<T> operator/(Dual<T> lhs, T rhs) {
    return lhs /= rhs;
}

// --- Transcendental functions with automatic differentiation ------------------

// sin(a + b eps) = sin(a) + b cos(a) eps.
template <typename T>
Dual<T> sin(const Dual<T>& d) {
    return Dual<T>(std::sin(d.real), d.dual * std::cos(d.real));
}

// cos(a + b eps) = cos(a) - b sin(a) eps.
template <typename T>
Dual<T> cos(const Dual<T>& d) {
    return Dual<T>(std::cos(d.real), -d.dual * std::sin(d.real));
}

// sqrt(a + b eps) = sqrt(a) + b/(2 sqrt(a)) eps.
template <typename T>
Dual<T> sqrt(const Dual<T>& d) {
    T sqrt_real = std::sqrt(d.real);
    return Dual<T>(sqrt_real, d.dual / (T(2.0) * sqrt_real));
}

// exp(a + b eps) = exp(a) + b exp(a) eps.
template <typename T>
Dual<T> exp(const Dual<T>& d) {
    T exp_real = std::exp(d.real);
    return Dual<T>(exp_real, d.dual * exp_real);
}

// log(a + b eps) = log(a) + b/a eps.
template <typename T>
Dual<T> log(const Dual<T>& d) {
    return Dual<T>(std::log(d.real), d.dual / d.real);
}

// pow(a + b eps, n) = a^n + b n a^(n-1) eps.
template <typename T>
Dual<T> pow(const Dual<T>& base, T exponent) {
    T pow_real = std::pow(base.real, exponent);
    return Dual<T>(pow_real, base.dual * exponent * std::pow(base.real, exponent - T(1)));
}

// atan2(y, x) with automatic differentiation.
template <typename T>
Dual<T> atan2(const Dual<T>& y, const Dual<T>& x) {
    T denom = x.real * x.real + y.real * y.real;
    return Dual<T>(std::atan2(y.real, x.real), (x.real * y.dual - y.real * x.dual) / denom);
}

// fabs(a + b eps) = |a| + sign(a) b eps.
template <typename T>
Dual<T> fabs(const Dual<T>& d) {
    T sign = (d.real >= T(0)) ? T(1) : T(-1);
    return Dual<T>(std::fabs(d.real), sign * d.dual);
}

// tan(a + b eps) = tan(a) + b sec^2(a) eps, with sec^2 = 1/cos^2.
template <typename T>
Dual<T> tan(const Dual<T>& d) {
    T tan_real = std::tan(d.real);
    T cos_real = std::cos(d.real);
    T sec2 = T(1) / (cos_real * cos_real);
    return Dual<T>(tan_real, d.dual * sec2);
}

// sinh(a + b eps) = sinh(a) + b cosh(a) eps.
template <typename T>
Dual<T> sinh(const Dual<T>& d) {
    return Dual<T>(std::sinh(d.real), d.dual * std::cosh(d.real));
}

// cosh(a + b eps) = cosh(a) + b sinh(a) eps.
template <typename T>
Dual<T> cosh(const Dual<T>& d) {
    return Dual<T>(std::cosh(d.real), d.dual * std::sinh(d.real));
}

// tanh(a + b eps) = tanh(a) + b sech^2(a) eps, with sech^2 = 1 - tanh^2.
template <typename T>
Dual<T> tanh(const Dual<T>& d) {
    T tanh_real = std::tanh(d.real);
    T sech2 = T(1) - tanh_real * tanh_real;
    return Dual<T>(tanh_real, d.dual * sech2);
}

// asin(a + b eps) = asin(a) + b/sqrt(1 - a^2) eps; domain |a| < 1.
template <typename T>
Dual<T> asin(const Dual<T>& d) {
    T one_minus_x2 = T(1) - d.real * d.real;
    if (one_minus_x2 <= T(0)) {
        one_minus_x2 = std::numeric_limits<T>::epsilon();
    }
    T deriv = T(1) / std::sqrt(one_minus_x2);
    return Dual<T>(std::asin(d.real), d.dual * deriv);
}

// acos(a + b eps) = acos(a) - b/sqrt(1 - a^2) eps; domain |a| < 1.
template <typename T>
Dual<T> acos(const Dual<T>& d) {
    T one_minus_x2 = T(1) - d.real * d.real;
    if (one_minus_x2 <= T(0)) {
        one_minus_x2 = std::numeric_limits<T>::epsilon();
    }
    T deriv = T(-1) / std::sqrt(one_minus_x2);
    return Dual<T>(std::acos(d.real), d.dual * deriv);
}

// atan(a + b eps) = atan(a) + b/(1 + a^2) eps.
template <typename T>
Dual<T> atan(const Dual<T>& d) {
    T deriv = T(1) / (T(1) + d.real * d.real);
    return Dual<T>(std::atan(d.real), d.dual * deriv);
}

// hypot(x, y) = sqrt(x^2 + y^2); derivatives x/hypot, y/hypot avoid overflow.
template <typename T>
Dual<T> hypot(const Dual<T>& x, const Dual<T>& y) {
    T hypot_real = std::hypot(x.real, y.real);
    if (hypot_real < std::numeric_limits<T>::epsilon()) {
        return Dual<T>(T(0), T(0));
    }
    T dual = (x.real * x.dual + y.real * y.dual) / hypot_real;
    return Dual<T>(hypot_real, dual);
}

// copysign(x, y) = |x| sign(y); the derivative follows x when the signs match.
template <typename T>
Dual<T> copysign(const Dual<T>& x, const Dual<T>& y) {
    T mag_sign = (y.real >= T(0)) ? T(1) : T(-1);
    T x_sign = (x.real >= T(0)) ? T(1) : T(-1);
    T deriv_factor = x_sign * mag_sign;
    return Dual<T>(std::copysign(x.real, y.real), deriv_factor * x.dual);
}

// abs(a + b eps): alias for fabs, matching std::abs spelling.
template <typename T>
Dual<T> abs(const Dual<T>& d) {
    return fabs(d);
}

// max(x, y) with derivative tracking.
template <typename T>
Dual<T> max(const Dual<T>& x, const Dual<T>& y) {
    if (x.real >= y.real) {
        return x;
    }
    return y;
}

// min(x, y) with derivative tracking.
template <typename T>
Dual<T> min(const Dual<T>& x, const Dual<T>& y) {
    if (x.real <= y.real) {
        return x;
    }
    return y;
}

// clamp(x, lo, hi): the derivative is zero at a clamped boundary, else passes
// through.
template <typename T>
Dual<T> clamp(const Dual<T>& x, const Dual<T>& lo, const Dual<T>& hi) {
    if (x.real < lo.real) {
        return Dual<T>(lo.real, T(0));
    }
    if (x.real > hi.real) {
        return Dual<T>(hi.real, T(0));
    }
    return x;
}

// clamp with scalar bounds.
template <typename T>
Dual<T> clamp(const Dual<T>& x, T lo, T hi) {
    if (x.real < lo) {
        return Dual<T>(lo, T(0));
    }
    if (x.real > hi) {
        return Dual<T>(hi, T(0));
    }
    return x;
}

// --- Comparison operators (real part only) -----------------------------------

template <typename T>
bool operator<(const Dual<T>& lhs, const Dual<T>& rhs) {
    return lhs.real < rhs.real;
}

template <typename T>
bool operator>(const Dual<T>& lhs, const Dual<T>& rhs) {
    return lhs.real > rhs.real;
}

template <typename T>
bool operator<=(const Dual<T>& lhs, const Dual<T>& rhs) {
    return lhs.real <= rhs.real;
}

template <typename T>
bool operator>=(const Dual<T>& lhs, const Dual<T>& rhs) {
    return lhs.real >= rhs.real;
}

template <typename T>
bool operator==(const Dual<T>& lhs, const Dual<T>& rhs) {
    return lhs.real == rhs.real;
}

template <typename T>
bool operator!=(const Dual<T>& lhs, const Dual<T>& rhs) {
    return lhs.real != rhs.real;
}

template <typename T>
bool operator<(const Dual<T>& lhs, T rhs) {
    return lhs.real < rhs;
}

template <typename T>
bool operator>(const Dual<T>& lhs, T rhs) {
    return lhs.real > rhs;
}

template <typename T>
bool operator<(T lhs, const Dual<T>& rhs) {
    return lhs < rhs.real;
}

template <typename T>
bool operator>(T lhs, const Dual<T>& rhs) {
    return lhs > rhs.real;
}

// --- Finiteness predicates ---------------------------------------------------

template <typename T>
bool isfinite(const Dual<T>& d) {
    return std::isfinite(d.real) && std::isfinite(d.dual);
}

template <typename T>
bool isnan(const Dual<T>& d) {
    return std::isnan(d.real) || std::isnan(d.dual);
}

template <typename T>
bool isinf(const Dual<T>& d) {
    return std::isinf(d.real) || std::isinf(d.dual);
}

}  // namespace sirius::core
