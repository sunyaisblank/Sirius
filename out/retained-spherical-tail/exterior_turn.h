#pragma once

#include "finite_tail.h"

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
    std::array<Jet, 4> x, k;
    for (unsigned i = 0; i < x.size(); ++i) {
        x[i] = Jet(position(i));
        k[i] = Jet(tangent(i));
    }
    const double orientation = chart == SphericalKerrSchildChart::Outgoing ? -1.0 : 1.0;
    const auto ray = Extract(x, k, mass, orientation, source_support_radius, boundary,
                             control.constraint_tolerance, SphericalTailPath::ExteriorTurn);
    if (!ray) return {};
    const auto turn =
        VacuumTurnRange(ray->radius.real, ray->impact_squared.real, mass, source_support_radius);
    if (!turn) return {};
    return Bound(*turn);
}
}  // namespace sirius::core::relativity
