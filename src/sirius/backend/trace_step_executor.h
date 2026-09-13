#pragma once

#include "sirius/core/camera_launch.h"
#include "sirius/core/geodesic_integrator.h"

namespace sirius::backend {

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
};

}  // namespace sirius::backend
