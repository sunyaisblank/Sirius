#pragma once

#include "sirius/core/geodesic_integrator.h"
#include "sirius/core/metrics/outgoing_kerr_schild.h"

#include "finite_tail.h"

namespace sirius::core::relativity {

struct SphericalCoupledVariationBounds {
    std::array<SphericalTailBound, 4> displacement;
    std::array<SphericalTailBound, 4> derivative;
};

struct SphericalCoupledTail {
    Vec4 position;
    Vec4 tangent;
    GeodesicVariations variations;
    double affine_length;
    unsigned panels;
    SphericalTailEventBounds event_bounds;
    std::array<SphericalCoupledVariationBounds, 4> variation_bounds;
    std::optional<SphericalTailBound> periapsis;
};

namespace finite_tail_detail {

inline Range Enclosure(SphericalTailBound bound) { return {bound.lower, bound.upper}; }

inline Four Components(const Vec4& vector) {
    return {Range(vector(0)), Range(vector(1)), Range(vector(2)), Range(vector(3))};
}

inline Four Components(const std::array<SphericalTailBound, 4>& bounds) {
    return {Enclosure(bounds[0]), Enclosure(bounds[1]), Enclosure(bounds[2]), Enclosure(bounds[3])};
}

struct CoupledLaunch {
    Four position, tangent;
    std::array<Four, 4> displacement, derivative;
};

inline Range Dot(const Four& a, const Four& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2] + a[3] * b[3];
}

// Contract the first-kind Christoffel symbols against two vectors using
// g=eta+h l*l, then raise the remaining index with eta-h l^sharp*l^sharp.
inline Four SphericalConnection(const Four& position, const Four& a, const Four& b, double mass,
                                double orientation, Range radius) {
    if (mass == 0.0) return {};
    const Range h = Range(2.0) * Range(mass) / radius;
    const Four n{Range(0.0), position[1] / radius, position[2] / radius, position[3] / radius};
    const Range sign(orientation);
    const Four l{Range(1.0), sign * n[1], sign * n[2], sign * n[3]};
    const Four raised{Range(-1.0), l[1], l[2], l[3]};
    const Range la = Dot(l, a), lb = Dot(l, b);
    Four first_kind;
    for (unsigned axis = 1; axis < n.size(); ++axis) {
        const Range dh = (Range(0.0) - h) * n[axis] / radius;
        Four dl;
        for (unsigned i = 1; i < n.size(); ++i)
            dl[i] = sign * (Range(i == axis ? 1.0 : 0.0) - n[i] * n[axis]) / radius;
        const Range dla = Dot(dl, a), dlb = Dot(dl, b);
        for (unsigned i = 0; i < n.size(); ++i) {
            const Range dgb = dh * l[i] * lb + h * (dl[i] * lb + l[i] * dlb);
            const Range dga = dh * l[i] * la + h * (dl[i] * la + l[i] * dla);
            first_kind[i] = first_kind[i] + a[axis] * dgb + b[axis] * dga;
        }
        first_kind[axis] = first_kind[axis] - (dh * la * lb + h * (dla * lb + la * dlb));
    }
    const Range inverse_correction = h * Dot(raised, first_kind);
    Four result;
    for (unsigned i = 0; i < n.size(); ++i) {
        const Range flat = Range(i == 0 ? -1.0 : 1.0) * first_kind[i];
        result[i] = Range(0.5) * (flat - raised[i] * inverse_correction);
    }
    return result;
}

inline std::expected<SphericalCoupledTail, FiniteTailFailure> ContinueCoupled(
    IMetric& metric, const CoupledLaunch& launch, double source_support_radius, double boundary,
    const SphericalTailControl& control, SphericalTailPath path) {
    if (Cancelled(control)) return std::unexpected(FiniteTailFailure::Cancelled);
    const auto* outgoing = dynamic_cast<const OutgoingKerrSchild*>(&metric);
    const auto* family =
        outgoing ? &outgoing->Source() : dynamic_cast<const KerrSchildFamily*>(&metric);
    if (!family) return std::unexpected(FiniteTailFailure::Domain);
    const auto parameters = family->GetParams();
    if (parameters.a != 0 || parameters.Q != 0 || parameters.Lambda != 0)
        return std::unexpected(FiniteTailFailure::Domain);
    Vec4 nominal_position;
    for (unsigned i = 0; i < launch.position.size(); ++i) {
        if (!launch.position[i].Finite() || !launch.tangent[i].Finite())
            return std::unexpected(FiniteTailFailure::Input);
        nominal_position(i) = std::midpoint(launch.position[i].lo, launch.position[i].hi);
    }
    if (!RoundingEnvironment() || !metric.IsValidEvent(nominal_position))
        return std::unexpected(FiniteTailFailure::Input);

    const double orientation = outgoing ? -1.0 : 1.0;
    const auto& start_x = launch.position;
    const auto& start_k = launch.tangent;
    const Range start_radius = sqrt(Square(start_x[1]) + Square(start_x[2]) + Square(start_x[3]));
    Columns x, k;
    for (unsigned i = 0; i < launch.displacement.size(); ++i) {
        const auto& seed_x = launch.displacement[i];
        const auto& seed_k = launch.derivative[i];
        const auto connection =
            SphericalConnection(start_x, start_k, seed_x, parameters.M, orientation, start_radius);
        for (unsigned j = 0; j < seed_x.size(); ++j) {
            x[i][j] = Jet(start_x[j], seed_x[j]);
            k[i][j] = Jet(start_k[j], seed_k[j] - connection[j]);
        }
    }
    const auto tail =
        Continue(x, k, parameters.M,
                 outgoing ? SphericalKerrSchildChart::Outgoing : SphericalKerrSchildChart::Ingoing,
                 source_support_radius, boundary, control, path);
    if (!tail) return std::unexpected(tail.error());
    if (!metric.IsValidEvent(tail->position)) return std::unexpected(FiniteTailFailure::Domain);

    const auto end_x = Components(tail->event_bounds.position);
    const auto end_k = Components(tail->event_bounds.tangent);
    SphericalCoupledTail result{tail->position, tail->tangent,      {}, tail->affine_length,
                                tail->panels,   tail->event_bounds, {}, tail->periapsis};
    for (unsigned i = 0; i < result.variations.size(); ++i) {
        const auto& event = tail->variation_bounds[i];
        const Four event_x = Components(event.position), event_k = Components(event.tangent);
        const auto connection =
            SphericalConnection(end_x, end_k, event_x, parameters.M, orientation, Range(boundary));
        const Range shift = Enclosure(event.affine_length);
        // X_event = X_fixed + k d(lambda_event). The covariant derivative
        // needs no affine shift since the covariant acceleration of k vanishes.
        auto& variation = result.variations[i];
        auto& bounds = result.variation_bounds[i];
        for (unsigned j = 0; j < event_x.size(); ++j) {
            const Range displacement = event_x[j] - end_k[j] * shift;
            const Range derivative = event_k[j] + connection[j];
            if (!Within(displacement, control.event_absolute_tolerance,
                        control.event_relative_tolerance) ||
                !Within(derivative, control.event_absolute_tolerance,
                        control.event_relative_tolerance))
                return std::unexpected(FiniteTailFailure::EventAccuracy);
            variation.displacement(j) = std::midpoint(displacement.lo, displacement.hi);
            variation.derivative(j) = std::midpoint(derivative.lo, derivative.hi);
            bounds.displacement[j] = Bound(displacement);
            bounds.derivative[j] = Bound(derivative);
        }
    }
    if (Cancelled(control)) return std::unexpected(FiniteTailFailure::Cancelled);
    return result;
}
}  // namespace finite_tail_detail

// The receiver can still run the original spherical event variation and
// source-sky code on the result. The fixed-affine state is supplied in the
// input metric's active chart and all four covariant derivatives travel with
// the central ray. No part of the trial is written to caller-owned state.
[[nodiscard]] inline std::expected<SphericalCoupledTail, FiniteTailFailure>
PrepareSphericalCoupledTail(IMetric& metric, const Vec4& position, const Vec4& tangent,
                            const GeodesicVariations& variations, double source_support_radius,
                            double boundary, const SphericalTailControl& control = {},
                            SphericalTailPath path = SphericalTailPath::Outward) {
    using namespace finite_tail_detail;
    CoupledLaunch launch{Components(position), Components(tangent), {}, {}};
    for (unsigned i = 0; i < variations.size(); ++i) {
        launch.displacement[i] = Components(variations[i].displacement);
        launch.derivative[i] = Components(variations[i].derivative);
    }
    return ContinueCoupled(metric, launch, source_support_radius, boundary, control, path);
}
}  // namespace sirius::core::relativity
