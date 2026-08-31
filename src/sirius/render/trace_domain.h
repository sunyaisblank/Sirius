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

// One bounded attempt envelope for both production backends. Rejected
// adaptive attempts count on each path, so this is a work and termination
// semantic rather than an integrator-specific tuning constant.
inline constexpr int kRenderTraceMaximumAttempts = 20000;
inline constexpr float kRenderTraceCaptureFactor = 1.0f;

struct TraceDomainParameters {
    float escape_radius;
    float cpu_initial_step;
    float cpu_min_step;
    float vulkan_min_step;
    float max_step;
    bool finite_causal_boundary;
};

struct TraceDomainRequest {
    core::MetricId metric_id;
    double metric_mass;
    double cosmological_constant;
    double observer_radius;
    double throat_radius;
    double bubble_radius;
    double bubble_sigma;
};

[[nodiscard]] inline TraceDomainParameters BuildTraceDomainParameters(
    const TraceDomainRequest& request) {
    SIRIUS_PRE(std::isfinite(request.metric_mass) && request.metric_mass >= 0.0);
    SIRIUS_PRE(std::isfinite(request.cosmological_constant) &&
               request.cosmological_constant >= 0.0);
    SIRIUS_PRE(std::isfinite(request.observer_radius) && request.observer_radius > 0.0);
    SIRIUS_PRE(std::isfinite(request.throat_radius) &&
               request.throat_radius >= core::kMinMorrisThorneThroatRadius &&
               request.throat_radius <= core::kMaxMorrisThorneThroatRadius);
    SIRIUS_PRE(std::isfinite(request.bubble_radius) && request.bubble_radius > 0.0);
    SIRIUS_PRE(std::isfinite(request.bubble_sigma) && request.bubble_sigma > 0.0);
    SIRIUS_PRE(core::MetricUsesMass(request.metric_id) ? request.metric_mass > 0.0
                                                       : request.metric_mass == 0.0);
    SIRIUS_PRE(
        !core::MetricParameterIssue(request.metric_id, 0.0, 0.0, request.cosmological_constant)
             .has_value());
    SIRIUS_PRE(!core::MetricHorizonIssue(request.metric_id, request.metric_mass,
                                         request.cosmological_constant)
                    .has_value());
    SIRIUS_PRE(
        request.metric_id != core::MetricId::Alcubierre ||
        !core::AlcubierreScaleIssue(request.bubble_radius, request.bubble_sigma).has_value());
    SIRIUS_PRE(core::MetricObserverRadiusIssueFor(
                   request.metric_id, request.metric_mass, request.cosmological_constant,
                   request.observer_radius, request.throat_radius, request.bubble_radius,
                   request.bubble_sigma) == core::MetricObserverRadiusIssue::None);

    const double scene_scale =
        core::MetricSceneLengthScale(request.metric_id, request.metric_mass, request.throat_radius,
                                     request.bubble_radius, request.bubble_sigma);
    const double feature_scale = core::MetricFeatureLengthScale(
        request.metric_id, request.metric_mass, request.throat_radius, request.bubble_radius,
        request.bubble_sigma);
    SIRIUS_PRE(std::isfinite(scene_scale) && scene_scale > 0.0);
    SIRIUS_PRE(std::isfinite(feature_scale) && feature_scale > 0.0);
    double escape_radius = std::max(200.0 * scene_scale, 1.25 * request.observer_radius);
    bool finite_causal_boundary = false;
    if (const auto cosmological_horizon = core::MetricCosmologicalHorizonRadius(
            request.metric_id, request.metric_mass, request.cosmological_constant);
        cosmological_horizon.has_value()) {
        // de Sitter has no asymptotic spatial infinity.  The directional sky is
        // imposed no later than the causal-patch boundary, never on an arbitrary
        // sphere beyond the cosmological horizon.
        escape_radius = std::min(escape_radius, *cosmological_horizon);
        finite_causal_boundary = true;
    }
    SIRIUS_POST(std::isfinite(escape_radius) && escape_radius > request.observer_radius);

    float stored_escape_radius = static_cast<float>(escape_radius);
    if (finite_causal_boundary && static_cast<double>(stored_escape_radius) > escape_radius) {
        stored_escape_radius = std::nextafter(stored_escape_radius, 0.0f);
    }
    SIRIUS_POST(std::isfinite(stored_escape_radius) &&
                stored_escape_radius > static_cast<float>(request.observer_radius));
    SIRIUS_POST(!finite_causal_boundary ||
                static_cast<double>(stored_escape_radius) <= escape_radius);

    return {
        stored_escape_radius,
        static_cast<float>(0.1 * feature_scale),
        static_cast<float>(1.0e-5 * feature_scale),
        static_cast<float>(0.02 * feature_scale),
        static_cast<float>(2.0 * scene_scale),
        finite_causal_boundary,
    };
}

}  // namespace sirius::render
