#pragma once

#include "sirius/core/spherical_continuation/radial_quadrature.h"

namespace sirius::core::relativity::finite_tail_detail {

// Enclose a fixed-pupil calculation as c + a*s + b*t + r, s,t in [-1,1].
// The same s and t are used in every component. Nonlinear terms enter r,
// while cancellation of shared linear terms is preserved.
struct PupilAffine {
    Range centre, first, second, remainder;

    PupilAffine() = default;
    explicit PupilAffine(double value) : centre(value) {}
    explicit PupilAffine(Range constant) : centre(constant) {}
    PupilAffine(Range c, Range a, Range b, Range r)
        : centre(c), first(a), second(b), remainder(r) {}

    static PupilAffine Coordinate(Range bound, bool first_axis) {
        if (!bound.Finite()) return PupilAffine(InvalidRange());
        const Range middle(std::midpoint(bound.lo, bound.hi));
        const Range radius(MidpointError(bound));
        return {middle, first_axis ? radius : Range{}, first_axis ? Range{} : radius, {}};
    }

    Range Variation() const {
        const Range magnitude = Abs(first) + Abs(second);
        return Range(-magnitude.hi, magnitude.hi) + remainder;
    }

    Range Bounds() const { return centre + Variation(); }
};

inline Range Enclosure(Range x) { return x; }
inline Range Enclosure(const PupilAffine& x) { return x.Bounds(); }

inline PupilAffine operator+(const PupilAffine& a, const PupilAffine& b) {
    return {a.centre + b.centre, a.first + b.first, a.second + b.second, a.remainder + b.remainder};
}

inline PupilAffine operator-(const PupilAffine& a, const PupilAffine& b) {
    return {a.centre - b.centre, a.first - b.first, a.second - b.second, a.remainder - b.remainder};
}

inline PupilAffine operator*(const PupilAffine& a, const PupilAffine& b) {
    return {a.centre * b.centre, a.centre * b.first + b.centre * a.first,
            a.centre * b.second + b.centre * a.second,
            a.centre * b.remainder + b.centre * a.remainder + a.Variation() * b.Variation()};
}

inline PupilAffine Reciprocal(const PupilAffine& x) {
    const Range bound = x.Bounds();
    if (!(bound.lo > 0.0)) return PupilAffine(InvalidRange());
    // 1/(c+d) = 1/c - d/c² + d²/[c²(c+d)].
    const Range c2 = Square(x.centre);
    const Range gradient = Range(-1.0) / c2;
    return {Range(1.0) / x.centre, gradient * x.first, gradient * x.second,
            gradient * x.remainder + Square(x.Variation()) / (c2 * bound)};
}

inline PupilAffine operator/(const PupilAffine& a, const PupilAffine& b) {
    return a * Reciprocal(b);
}

inline PupilAffine sqrt(const PupilAffine& x) {
    const Range central_root = sqrt(x.centre);
    const Range full_root = sqrt(x.Bounds());
    // sqrt(c+d) = sqrt(c) + d/(2 sqrt(c))
    //             - d²/[2 sqrt(c) (sqrt(c+d) + sqrt(c))²].
    const Range twice = Range(2.0) * central_root;
    const Range gradient = Range(1.0) / twice;
    return {central_root, gradient * x.first, gradient * x.second,
            gradient * x.remainder -
                Square(x.Variation()) / (twice * Square(full_root + central_root))};
}

inline PupilAffine Square(const PupilAffine& x) {
    const Range twice = Range(2.0) * x.centre;
    return {Square(x.centre), twice * x.first, twice * x.second,
            twice * x.remainder + Square(x.Variation())};
}

}  // namespace sirius::core::relativity::finite_tail_detail
