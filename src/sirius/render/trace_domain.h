#pragma once

// Shared geometric trace-domain authority for the CPU and Vulkan render paths.
// The escape sphere and largest step follow the full scene scale; the initial
// and minimum steps follow its smallest feature. The escape sphere must enclose
// the launch point: the Vulkan kernel tests it before taking its first step, and
// an absolute 200-unit sphere used to make every more-distant observer escape
// without ever reaching the lens.

#include "sirius/base/contracts.h"
#include "sirius/core/metrics/registry.h"

#include <algorithm>
#include <cmath>

namespace sirius::render {

struct TraceDomainParameters {
    float escape_radius;
    float cpu_initial_step;
    float cpu_min_step;
    float vulkan_min_step;
    float max_step;
};

struct TraceDomainRequest {
    core::MetricId metric_id;
    double metric_mass;
    double observer_radius;
    double throat_radius;
    double bubble_radius;
    double bubble_sigma;
};

[[nodiscard]] inline TraceDomainParameters BuildTraceDomainParameters(
    const TraceDomainRequest& request) {
    SIRIUS_PRE(std::isfinite(request.metric_mass) && request.metric_mass >= 0.0);
    SIRIUS_PRE(std::isfinite(request.observer_radius) && request.observer_radius > 0.0);
    SIRIUS_PRE(std::isfinite(request.throat_radius) && request.throat_radius > 0.0);
    SIRIUS_PRE(std::isfinite(request.bubble_radius) && request.bubble_radius > 0.0);
    SIRIUS_PRE(std::isfinite(request.bubble_sigma) && request.bubble_sigma > 0.0);
    SIRIUS_PRE(core::MetricUsesMass(request.metric_id) ? request.metric_mass > 0.0
                                                       : request.metric_mass == 0.0);
    SIRIUS_PRE(
        request.metric_id != core::MetricId::Alcubierre ||
        !core::AlcubierreScaleIssue(request.bubble_radius, request.bubble_sigma).has_value());

    const double scene_scale =
        core::MetricSceneLengthScale(request.metric_id, request.metric_mass, request.throat_radius,
                                     request.bubble_radius, request.bubble_sigma);
    const double feature_scale = core::MetricFeatureLengthScale(
        request.metric_id, request.metric_mass, request.throat_radius, request.bubble_radius,
        request.bubble_sigma);
    SIRIUS_PRE(std::isfinite(scene_scale) && scene_scale > 0.0);
    SIRIUS_PRE(std::isfinite(feature_scale) && feature_scale > 0.0);
    const double escape_radius = std::max(200.0 * scene_scale, 1.25 * request.observer_radius);

    return {
        static_cast<float>(escape_radius),          static_cast<float>(0.1 * feature_scale),
        static_cast<float>(1.0e-5 * feature_scale), static_cast<float>(0.02 * feature_scale),
        static_cast<float>(2.0 * scene_scale),
    };
}

}  // namespace sirius::render
