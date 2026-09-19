#pragma once

#include "sirius/backend/retained_compute.h"
#include "sirius/core/geodesic_integrator.h"

namespace sirius::backend {

struct RetainedIntervalControl {
    core::IntegratorConfig integrator;
    double length_scale = 0;
    double frequency_scale = 0;
    double tolerance = 0;
    std::array<double, 4> column_scale{1, 1, 1, 1};
};

struct RetainedIntervalInput {
    std::array<RetainedValue, 4> metric;
    RetainedEndpointOutput start;
    double chart = 0;
    double interval = 0;
    RetainedIntervalControl control;
};

struct RetainedIntervalOutput {
    // Private until all embedded, midpoint and refined comparisons admit.
    // The tracer must still compare localized events before committing it.
    RetainedEndpointOutput full, lower, midpoint, refined;
    std::array<RetainedValue, 20> full_increment{}, lower_increment{}, midpoint_increment{},
        refined_increment{};
    std::uint32_t attempted_stages = 0;
    double error_ratio = 0;
    core::CoupledStepFailure failure = core::CoupledStepFailure::None;
    bool admissible = false;
};

// One bounded batch of private coupled attempts. RK, projection and dense
// evaluation run on the device. The host forms retained lower increments and
// compares centers with the existing component budgets; no CPU ray is substituted.
// Arithmetic radii remain attached to each result and private substage. The
// embedded/refinement center comparisons are local numerical error controls,
// not interval certification. An admitted expansion may define the exact
// numerical initial value of a new local problem without carrying the prior
// problem's arithmetic box as an independent global uncertainty.
[[nodiscard]] base::Expected<std::vector<RetainedIntervalOutput>> AttemptRetainedIntervals(
    RetainedCompute& compute, std::span<const RetainedIntervalInput> inputs);

[[nodiscard]] double RetainedPhysicalError(const std::array<RetainedValue, 40>& first,
                                           const std::array<RetainedValue, 40>& second,
                                           const RetainedIntervalControl& control);

}  // namespace sirius::backend
