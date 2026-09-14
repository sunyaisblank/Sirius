#pragma once
#include <cmath>

#if defined(_MSC_VER)
#define SIRIUS_TWOFOLD_INLINE __forceinline
#elif defined(__GNUC__) || defined(__clang__)
#define SIRIUS_TWOFOLD_INLINE inline __attribute__((always_inline))
#else
#define SIRIUS_TWOFOLD_INLINE inline
#endif

namespace sirius::core {
// A normalized sum of two binary64 values. Error-free sums/products retain
// arithmetic lost when large coordinate terms cancel into small covariant
// variations. This is extra working precision, not an interval enclosure.
// Operations require ordinary IEEE arithmetic (no reassociation/fast-math);
// callers reject nonfinite results before publishing a physical state.
struct Twofold {
    double hi = 0, lo = 0;
    Twofold(double value = 0) : hi(value) {}
    // Establish the normalized representation once, including a leading zero.
    Twofold(double high, double low) { *this = Sum(high, low); }

  private:
    struct RawTag {};
    Twofold(RawTag, double high, double low) : hi(high), lo(low) {}
    static SIRIUS_TWOFOLD_INLINE Twofold Raw(double high, double low) {
        return {RawTag{}, high, low};
    }

  public:
    double Rounded() const { return hi + lo; }
    static SIRIUS_TWOFOLD_INLINE Twofold Sum(double a, double b) {
        const double s = a + b, z = s - a;
        return Raw(s, (a - (s - z)) + (b - z));
    }
    static SIRIUS_TWOFOLD_INLINE Twofold Product(double a, double b) {
        const double p = a * b;
        // Exact zero factors have no residual. Keep the signed high zero and
        // avoid an FMA call for the many zero metric/variation coefficients.
        // Nonzero factors that underflow, and nonfinite products, still use
        // the error-free transform below.
        if (p == 0 && (a == 0 || b == 0)) return Raw(p, 0);
        return Raw(p, std::fma(a, b, -p));
    }

  private:
    static SIRIUS_TWOFOLD_INLINE Twofold FastSum(double a, double b) {
        const double s = a + b;
        const double e = b - (s - a);
        if (s == 0 && e == 0) return {};
        return Raw(s, e);
    }
    // The short product keeps its leading error-free transform away from
    // overflow and underflow. Extreme products retain the general fallback.
    static bool ProductRange(double p) { return std::abs(p) >= 0x1p-900 && std::abs(p) <= 0x1p900; }

  public:
    // Accurate double-word addition: Joldes, Muller and Popescu (2017),
    // Algorithm 5. Renormalize both the high and low error-free sums.
    friend SIRIUS_TWOFOLD_INLINE Twofold operator+(Twofold a, Twofold b) {
        if (a.hi == 0 && std::isfinite(b.hi)) return b.hi == 0 ? Twofold{} : b;
        if (b.hi == 0 && std::isfinite(a.hi)) return a;
        const auto s = Sum(a.hi, b.hi), t = Sum(a.lo, b.lo);
        const auto v = FastSum(s.hi, s.lo + t.hi);
        return FastSum(v.hi, t.lo + v.lo);
    }
    friend SIRIUS_TWOFOLD_INLINE Twofold operator-(Twofold a) { return Raw(-a.hi, -a.lo); }
    friend SIRIUS_TWOFOLD_INLINE Twofold operator-(Twofold a, Twofold b) { return a + -b; }
    friend SIRIUS_TWOFOLD_INLINE Twofold operator*(Twofold a, Twofold b) {
        if ((a.hi == 0 || b.hi == 0) && std::isfinite(a.hi) && std::isfinite(b.hi)) return {};
        const auto leading = Product(a.hi, b.hi);
        if (ProductRange(leading.hi)) {
            if (a.lo == 0 && b.lo == 0) return leading;
            if (a.lo == 0) return FastSum(leading.hi, std::fma(a.hi, b.lo, leading.lo));
            if (b.lo == 0) return FastSum(leading.hi, std::fma(a.lo, b.hi, leading.lo));
            // Keep the low*low term even when the two cross terms cancel:
            // sparse limbs can carry a represented tail far below 106 bits.
            const double cross = a.hi * b.lo + a.lo * b.hi;
            return FastSum(leading.hi, std::fma(a.lo, b.lo, leading.lo + cross));
        }
        return (Product(a.hi, b.hi) + Product(a.hi, b.lo)) +
               (Product(a.lo, b.hi) + Product(a.lo, b.lo));
    }
    friend SIRIUS_TWOFOLD_INLINE Twofold operator*(Twofold a, double b) {
        if ((a.hi == 0 || b == 0) && std::isfinite(a.hi) && std::isfinite(b)) return {};
        if (b == 1 && std::isfinite(a.hi)) return a;
        if (b == -1 && std::isfinite(a.hi)) return -a;
        const double scaled = a.hi * b;
        // These binary scalings are exact in the guarded normal range.
        if ((b == .5 || b == -.5 || b == 2 || b == -2) && ProductRange(scaled)) {
            const double low = a.lo * b;
            return Raw(scaled, low == 0 ? 0 : low);
        }
        const auto leading = Product(a.hi, b);
        if (ProductRange(leading.hi)) return FastSum(leading.hi, std::fma(a.lo, b, leading.lo));
        return leading + Product(a.lo, b);
    }
    friend SIRIUS_TWOFOLD_INLINE Twofold operator*(double a, Twofold b) { return b * a; }
    friend SIRIUS_TWOFOLD_INLINE Twofold operator/(Twofold a, Twofold b) {
        Twofold q(a.hi / b.hi);
        const auto residual = a - b * q;
        return q + Twofold(residual.Rounded() / b.hi);
    }
    friend SIRIUS_TWOFOLD_INLINE Twofold operator/(Twofold a, double b) {
        const Twofold q(a.hi / b);
        const auto residual = a - q * b;
        return q + Twofold(residual.Rounded() / b);
    }
    SIRIUS_TWOFOLD_INLINE Twofold& operator+=(Twofold b) { return *this = *this + b; }
    SIRIUS_TWOFOLD_INLINE Twofold& operator-=(Twofold b) { return *this = *this - b; }
};

inline Twofold sqrt(Twofold value) {
    const Twofold initial(std::sqrt(value.Rounded()));
    if (initial.hi == 0.0 && value.hi == 0.0 && value.lo == 0.0) return initial;
    return initial + (value - initial * initial) / (initial * 2.0);
}
}  // namespace sirius::core

#undef SIRIUS_TWOFOLD_INLINE
