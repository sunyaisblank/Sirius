#pragma once

#include <array>

namespace sirius::core {

// Value, gradient and Hessian of a scalar function of three spatial
// coordinates. Used to differentiate the same metric scalar/null-vector
// authorities as the nominal path without neighboring metric evaluations.
struct SecondOrder3 {
    double value = 0;
    std::array<double, 3> gradient{};
    std::array<std::array<double, 3>, 3> hessian{};

    SecondOrder3(double scalar = 0) : value(scalar) {}
    static SecondOrder3 Variable(double scalar, int axis) {
        SecondOrder3 result(scalar);
        result.gradient[axis] = 1;
        return result;
    }
};

inline SecondOrder3 operator+(const SecondOrder3& a, const SecondOrder3& b) {
    SecondOrder3 result(a.value + b.value);
    for (int i = 0; i < 3; ++i) {
        result.gradient[i] = a.gradient[i] + b.gradient[i];
        for (int j = 0; j < 3; ++j) result.hessian[i][j] = a.hessian[i][j] + b.hessian[i][j];
    }
    return result;
}

inline SecondOrder3 operator-(const SecondOrder3& a, const SecondOrder3& b) {
    SecondOrder3 result(a.value - b.value);
    for (int i = 0; i < 3; ++i) {
        result.gradient[i] = a.gradient[i] - b.gradient[i];
        for (int j = 0; j < 3; ++j) result.hessian[i][j] = a.hessian[i][j] - b.hessian[i][j];
    }
    return result;
}

inline SecondOrder3 operator*(const SecondOrder3& a, const SecondOrder3& b) {
    SecondOrder3 result(a.value * b.value);
    for (int i = 0; i < 3; ++i) {
        result.gradient[i] = a.gradient[i] * b.value + a.value * b.gradient[i];
        for (int j = 0; j < 3; ++j)
            result.hessian[i][j] = a.hessian[i][j] * b.value + a.value * b.hessian[i][j] +
                                   a.gradient[i] * b.gradient[j] + a.gradient[j] * b.gradient[i];
    }
    return result;
}

inline SecondOrder3 operator/(const SecondOrder3& a, const SecondOrder3& b) {
    SecondOrder3 inverse(1 / b.value);
    for (int i = 0; i < 3; ++i) {
        inverse.gradient[i] = -(b.gradient[i] / b.value) * inverse.value;
        for (int j = 0; j < 3; ++j)
            inverse.hessian[i][j] = (2 * (b.gradient[i] / b.value) * (b.gradient[j] / b.value) -
                                     b.hessian[i][j] / b.value) *
                                    inverse.value;
    }
    SecondOrder3 result = a * inverse;
    result.value = a.value / b.value;
    return result;
}

}  // namespace sirius::core
