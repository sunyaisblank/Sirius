#pragma once

#include "coupled_handoff.h"

#include <optional>

namespace sirius::core::relativity {

// Outside 3M, 1 - b²(1 - 2M/r)/r² is strictly increasing for b² > 0.
// A negative value at the exterior support boundary and a positive value at
// launch therefore certify a unique simple periapsis outside all sources.
[[nodiscard]] inline std::optional<SphericalTailBound> BracketVacuumTurn(
    const Vec4& position, const Vec4& tangent, double mass, SphericalKerrSchildChart chart,
    double source_support_radius, double boundary, const SphericalTailControl& control = {}) {
    using namespace finite_tail_detail;
    if (!RoundingEnvironment() || !std::isfinite(mass) || mass < 0.0 || !std::isfinite(boundary) ||
        !std::isfinite(control.constraint_tolerance) || control.constraint_tolerance < 0.0 ||
        !std::isfinite(source_support_radius))
        return {};
    const Range two_mass = Range(2.0) * Range(mass);
    double lower = source_support_radius;
    if (!(lower > (Range(3.0) * Range(mass)).hi)) return {};
    std::array<Jet, 4> x, k;
    for (unsigned i = 0; i < x.size(); ++i) {
        x[i] = Jet(position(i));
        k[i] = Jet(tangent(i));
    }
    const double orientation = chart == SphericalKerrSchildChart::Outgoing ? -1.0 : 1.0;
    const auto ray = Extract(x, k, mass, orientation, lower, boundary, control.constraint_tolerance,
                             RadialBranch::Inward);
    if (!ray) return {};
    const auto potential = [&](double radius) {
        const Range r(radius);
        return Range(1.0) - (Range(1.0) - two_mass / r) * ray->impact_squared.real / Square(r);
    };
    double upper = ray->radius.real.lo;
    if (!(potential(lower).hi < 0.0 && potential(upper).lo > 0.0)) return {};
    // Keep only endpoints with certified, opposite signs. Stop as soon as
    // rounding in the launch invariants makes the sign at a midpoint unclear.
    for (unsigned i = 0; i < std::numeric_limits<double>::digits; ++i) {
        const double middle = std::midpoint(lower, upper);
        if (!(lower < middle && middle < upper)) break;
        const Range value = potential(middle);
        if (value.hi < 0.0)
            lower = middle;
        else if (value.lo > 0.0)
            upper = middle;
        else
            break;
    }
    return SphericalTailBound{lower, upper};
}
}  // namespace sirius::core::relativity
