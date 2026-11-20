#pragma once

#include "sirius/core/camera_launch.h"
#include "sirius/core/geodesic_integrator.h"

#include <functional>

namespace sirius::backend {

struct SphericalEscapeRequest {
    double source_support_radius;
    double boundary;
    unsigned maximum_panels;
    bool through_turn;
};

struct SphericalEscapeTrial {
    core::Vec4 position;
    core::Vec4 tangent;
    core::GeodesicVariations variations;
    double affine_length;
    std::optional<double> periapsis;
    unsigned panels;
};

// Optional execution seam for the shared trace/event/source coordinator.
// Its default remains the CPU integrator. An installed executor owns every
// attempted joint interval and must fail explicitly instead of falling back.
// Successful full/lower/midpoint/refined endpoints are projected physical views
// of that executor's working state. A fixed-affine endpoint must not be rebuilt
// from rounded host increments before its next interval.
class TraceStepExecutor {
  public:
    virtual ~TraceStepExecutor() = default;
    virtual void BeginTrace() {}
    virtual void EndTrace() {}
    virtual void RejectLastInterval() {}
    virtual std::optional<core::CameraLaunch> Launch(core::IMetric& metric, double absolute_spin,
                                                     const core::CameraRay& camera) = 0;
    virtual bool Step(core::Lightray& ray, core::IMetric& metric,
                      const core::IntegratorConfig& config, core::Rk45CoupledState& coupled,
                      core::Rk45CoupledComparison& comparison) = 0;
    // Sample the executor-owned trial at its original fraction. A declined
    // sample never authorizes a host reconstruction of the retained state.
    virtual std::optional<core::CoupledSegmentSample> Sample(core::CoupledTrial, double,
                                                             const core::Vec4* = nullptr) {
        return std::nullopt;
    }
    // An executor can offer a certified terminal event from its committed
    // state. A declined proposal leaves that state available for the next Step.
    virtual bool SupportsSphericalEscape() const { return false; }
    virtual std::optional<SphericalEscapeTrial> TrySphericalEscape(core::IMetric&,
                                                                   const core::Lightray&,
                                                                   const core::IntegratorConfig&,
                                                                   const core::Rk45CoupledState&,
                                                                   const SphericalEscapeRequest&,
                                                                   const std::function<bool()>&) {
        return std::nullopt;
    }
    virtual void RecordSphericalEscape(const SphericalEscapeTrial&) {}
};

}  // namespace sirius::backend
