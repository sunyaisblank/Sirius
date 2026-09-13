#pragma once
#include <cmath>
namespace sirius::core {
// A normalized sum of two binary64 values. Error-free sums/products retain
// arithmetic lost when large coordinate terms cancel into small covariant
// variations. This is extra working precision, not an interval enclosure.
// Operations require ordinary IEEE arithmetic (no reassociation/fast-math);
// callers reject nonfinite results before publishing a physical state.
struct Twofold {
    double hi = 0, lo = 0;
    Twofold(double value = 0) : hi(value) {}
    Twofold(double high, double low) : hi(high), lo(low) {}
    double Rounded() const { return hi + lo; }
    static Twofold Sum(double a, double b) {
        const double s = a + b, z = s - a;
        return {s, (a - (s - z)) + (b - z)};
    }
    static Twofold Product(double a, double b) {
        const double p = a * b;
        return {p, std::fma(a, b, -p)};
    }
    friend Twofold operator+(Twofold a, Twofold b) {
        const auto s = Sum(a.hi, b.hi), t = Sum(a.lo, b.lo);
        const auto u = Sum(s.lo, t.hi), v = Sum(s.hi, u.hi);
        return Sum(v.hi, ((v.lo + u.lo) + t.lo));
    }
    friend Twofold operator-(Twofold a) { return {-a.hi, -a.lo}; }
    friend Twofold operator-(Twofold a, Twofold b) { return a + -b; }
    friend Twofold operator*(Twofold a, Twofold b) {
        return (Product(a.hi, b.hi) + Product(a.hi, b.lo)) +
               (Product(a.lo, b.hi) + Product(a.lo, b.lo));
    }
    friend Twofold operator*(Twofold a, double b) { return Product(a.hi, b) + Product(a.lo, b); }
    friend Twofold operator*(double a, Twofold b) { return b * a; }
    friend Twofold operator/(Twofold a, Twofold b) {
        Twofold q(a.hi / b.hi);
        const auto residual = a - b * q;
        return q + Twofold(residual.Rounded() / b.hi);
    }
    friend Twofold operator/(Twofold a, double b) {
        const Twofold q(a.hi / b);
        const auto residual = a - q * b;
        return q + Twofold(residual.Rounded() / b);
    }
    Twofold& operator+=(Twofold b) { return *this = *this + b; }
    Twofold& operator-=(Twofold b) { return *this = *this - b; }
};

inline Twofold sqrt(Twofold value) {
    const Twofold initial(std::sqrt(value.Rounded()));
    if (initial.hi == 0.0 && value.hi == 0.0 && value.lo == 0.0) return initial;
    return initial + (value - initial * initial) / (initial * 2.0);
}
}  // namespace sirius::core
