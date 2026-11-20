#pragma once

#include "sirius/core/geodesic_integrator.h"
#include "sirius/core/metrics/outgoing_kerr_schild.h"
#include "sirius/core/observer_frame.h"
#include "sirius/core/polarisation/walker_penrose.h"

#include "finite_tail.h"

namespace sirius::core::relativity {

struct SphericalCoupledTail {
    Vec4 position;
    Vec4 tangent;
    GeodesicVariations variations;
    double affine_length;
    unsigned panels;
};

// The receiver can still run the original spherical event variation and
// source-sky code on the result. The fixed-affine state is supplied in the
// input metric's active chart and all four covariant derivatives travel with
// the central ray. No part of the trial is written to caller-owned state.
[[nodiscard]] inline std::expected<SphericalCoupledTail, FiniteTailFailure>
PrepareSphericalCoupledTail(IMetric& metric, const Vec4& position, const Vec4& tangent,
                            const GeodesicVariations& variations, double source_support_radius,
                            double boundary, const SphericalTailControl& control = {}) {
    const auto* outgoing = dynamic_cast<const OutgoingKerrSchild*>(&metric);
    const auto* family =
        outgoing ? &outgoing->Source() : dynamic_cast<const KerrSchildFamily*>(&metric);
    if (!family) return std::unexpected(FiniteTailFailure::Domain);
    const auto parameters = family->GetParams();
    if (parameters.a != 0 || parameters.Q != 0 || parameters.Lambda != 0)
        return std::unexpected(FiniteTailFailure::Domain);
    if (!IsFinite(position) || !IsFinite(tangent) || !metric.IsValidEvent(position))
        return std::unexpected(FiniteTailFailure::Input);

    Metric4d values;
    Tensor<Dual<double>, 4, 4, 4> derivatives;
    metric.Evaluate(position, values, derivatives);
    const auto start_connection = TensorOps::Christoffel(values, derivatives);
    std::array<SphericalTailVariation, 4> columns;
    for (unsigned i = 0; i < columns.size(); ++i) {
        const auto& variation = variations[i];
        columns[i] = {variation.displacement,
                      variation.derivative + NegConnectionContraction(start_connection, tangent,
                                                                      variation.displacement)};
    }
    const auto tail = ContinueSphericalTail(
        position, tangent, columns, parameters.M,
        outgoing ? SphericalKerrSchildChart::Outgoing : SphericalKerrSchildChart::Ingoing,
        source_support_radius, boundary, control);
    if (!tail) return std::unexpected(tail.error());
    if (!metric.IsValidEvent(tail->position)) return std::unexpected(FiniteTailFailure::Domain);

    metric.Evaluate(tail->position, values, derivatives);
    const auto end_connection = TensorOps::Christoffel(values, derivatives);
    SphericalCoupledTail result{
        tail->position, tail->tangent, {}, tail->affine_length, tail->panels};
    for (unsigned i = 0; i < columns.size(); ++i) {
        const auto& event = tail->event_variations[i];
        // X_event = X_fixed + k d(lambda_event). The covariant derivative
        // needs no affine shift since the covariant acceleration of k vanishes.
        auto& variation = result.variations[i];
        variation.displacement = event.position - tail->tangent * tail->affine_variations[i];
        variation.derivative =
            event.tangent - NegConnectionContraction(end_connection, tail->tangent, event.position);
        if (!IsFinite(variation.displacement) || !IsFinite(variation.derivative))
            return std::unexpected(FiniteTailFailure::Arithmetic);
    }
    return result;
}
}  // namespace sirius::core::relativity
