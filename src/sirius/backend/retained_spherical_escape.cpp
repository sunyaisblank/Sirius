#include "sirius/backend/retained_trace_executor.h"
#include "sirius/core/spherical_continuation/coupled_admission.h"

namespace sirius::backend {

void RetainedTraceExecutor::RecordSphericalEscape(const SphericalEscapeTrial& trial) {
    std::lock_guard lock(mutex_);
    ++stats_.certified_spherical_escapes;
    stats_.certified_spherical_turns += trial.periapsis.has_value();
    stats_.certified_spherical_panels += trial.panels;
}

std::optional<SphericalEscapeTrial> RetainedTraceExecutor::TrySphericalEscape(
    core::IMetric& metric, const core::Lightray& ray, const core::IntegratorConfig& config,
    const core::Rk45CoupledState& coupled, const SphericalEscapeRequest& request,
    const std::function<bool()>& should_cancel) {
    using namespace core::relativity;
    using namespace finite_tail_detail;
    if (!SupportsSphericalEscape() || request.maximum_panels == 0) return {};
    const std::function<bool()> cancelled = [&] {
        return (should_cancel && should_cancel()) || (should_cancel_ && should_cancel_());
    };
    if (cancelled()) return {};
    const auto* outgoing = dynamic_cast<const core::OutgoingKerrSchild*>(&metric);
    const auto* family =
        outgoing ? &outgoing->Source() : dynamic_cast<const core::KerrSchildFamily*>(&metric);
    if (!family) return {};
    const auto parameters = family->GetParams();
    if (parameters.a != 0 || parameters.Q != 0 || parameters.Lambda != 0) return {};

    const auto snapshot = PhysicalSnapshot(ray, coupled.variations, parameters, outgoing ? -1 : 1);
    std::array<RetainedValue, 40> physical;
    {
        std::lock_guard lock(mutex_);
        const auto thread = std::this_thread::get_id();
        const auto found = continuations_.find(thread);
        if (error_ || stopping_ || !active_traces_.contains(thread) ||
            found == continuations_.end() || !found->second.finish.valid || !found->second.dopri ||
            !(snapshot == found->second.after))
            return {};
        physical = found->second.finish.physical;
    }

    // The committed expansion defines the next local initial value, as in
    // Step. Preserve all three limbs of each projected physical component.
    const auto expansion = [&](std::size_t i) {
        const auto& value = physical[i];
        return value.IsRepresented() ? Range(value.high) + Range(value.low) + Range(value.tail)
                                     : InvalidRange();
    };
    CoupledLaunch launch;
    for (unsigned i = 0; i < launch.position.size(); ++i) {
        launch.position[i] = expansion(i);
        launch.tangent[i] = expansion(4 + i);
        for (unsigned column = 0; column < launch.displacement.size(); ++column) {
            launch.displacement[column][i] = expansion(8 + 8 * column + i);
            launch.derivative[column][i] = expansion(12 + 8 * column + i);
        }
    }
    SphericalTailControl control;
    control.maximum_panels = std::min(control.maximum_panels, request.maximum_panels);
    control.should_cancel = &cancelled;
    const auto path =
        request.through_turn ? SphericalTailPath::ExteriorTurn : SphericalTailPath::Outward;
    const auto tail = ContinueCoupled(metric, launch, request.source_support_radius,
                                      request.boundary, control, path);
    if (!tail || !(SphericalCoupledStateErrorBound(*tail, config, coupled) <= 1.0) || cancelled())
        return {};
    const auto periapsis =
        tail->periapsis
            ? std::optional(std::midpoint(tail->periapsis->lower, tail->periapsis->upper))
            : std::nullopt;
    return SphericalEscapeTrial{tail->position,      tail->tangent, tail->variations,
                                tail->affine_length, periapsis,     tail->panels};
}
}  // namespace sirius::backend
