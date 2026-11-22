#pragma once

#include "sirius/core/camera_launch.h"
#include "sirius/core/geodesic_integrator.h"
#include "sirius/core/spherical_continuation/radial_quadrature.h"

#include <functional>

namespace sirius::backend {

struct DetectorLaunchFrame {
    using Range = core::relativity::finite_tail_detail::Range;
    using Column = std::array<Range, 4>;
    using Variations = std::array<Column, 4>;
    Column position;
    std::array<Column, 4> tetrad;
    Variations displacement;
    // Covariant derivative of -e_0 and each e_i, indexed by basis, seed, component.
    std::array<Variations, 4> covariant;
    std::array<double, 4> metric;
    core::CameraRay pupil;
};

struct DetectorCameraLaunch {
    core::CameraLaunch central;
    DetectorLaunchFrame::Column position;
    DetectorLaunchFrame::Column tangent;
};

[[nodiscard]] std::optional<DetectorCameraLaunch> ProjectDetectorLaunch(const DetectorLaunchFrame&,
                                                                        const core::CameraRay&);

struct SphericalEscapeRequest {
    struct ThinDiskClearance {
        double inner_radius;
        double outer_radius;
        double minimum_radius;
    };

    double source_support_radius;
    double boundary;
    unsigned maximum_panels;
    bool through_turn;
    std::optional<ThinDiskClearance> thin_disk = std::nullopt;
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
    // The launch event and tetrad stay fixed as the film coordinates vary at
    // a fixed pupil. Only an executor that encloses both may certify them.
    virtual std::optional<DetectorLaunchFrame> EncloseDetectorLaunchFrame(core::IMetric&, double,
                                                                          const core::CameraRay&) {
        return std::nullopt;
    }
};

}  // namespace sirius::backend
