#pragma once

#include <array>
#include <cmath>

namespace sirius::core {

// Differentiate three Cartesian coordinates while retaining the arithmetic
// precision of Scalar in both the value and its gradient.
template <typename Scalar>
struct FirstOrder3 {
    Scalar value;
    std::array<Scalar, 3> gradient{};

    FirstOrder3(Scalar scalar = 0) : value(scalar) {}
    static FirstOrder3 Variable(Scalar scalar, int axis) {
        FirstOrder3 result(scalar);
        result.gradient[axis] = 1;
        return result;
    }
    friend FirstOrder3 operator+(const FirstOrder3& a, const FirstOrder3& b) {
        FirstOrder3 result(a.value + b.value);
        for (int i = 0; i < 3; ++i) result.gradient[i] = a.gradient[i] + b.gradient[i];
        return result;
    }
    friend FirstOrder3 operator-(const FirstOrder3& a, const FirstOrder3& b) {
        FirstOrder3 result(a.value - b.value);
        for (int i = 0; i < 3; ++i) result.gradient[i] = a.gradient[i] - b.gradient[i];
        return result;
    }
    friend FirstOrder3 operator*(const FirstOrder3& a, const FirstOrder3& b) {
        FirstOrder3 result(a.value * b.value);
        for (int i = 0; i < 3; ++i)
            result.gradient[i] = a.gradient[i] * b.value + a.value * b.gradient[i];
        return result;
    }
    friend FirstOrder3 operator/(const FirstOrder3& a, const FirstOrder3& b) {
        FirstOrder3 result(a.value / b.value);
        for (int i = 0; i < 3; ++i)
            result.gradient[i] = (a.gradient[i] - result.value * b.gradient[i]) / b.value;
        return result;
    }
    friend FirstOrder3 sqrt(const FirstOrder3& a) {
        using std::sqrt;
        FirstOrder3 result(sqrt(a.value));
        for (int i = 0; i < 3; ++i) result.gradient[i] = a.gradient[i] / (result.value * Scalar(2));
        return result;
    }
};

}  // namespace sirius::core
