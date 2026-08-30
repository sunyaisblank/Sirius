#pragma once

// Shared geometric trace-domain authority for the CPU and Vulkan render paths.
// Black-hole coordinates scale with M, while massless scenes use one geometric
// unit. The escape sphere must enclose the launch point: the Vulkan kernel tests
// it before taking its first step, and an absolute 200-unit sphere used to make
// every observer beyond r=200 escape without ever reaching the lens.

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

[[nodiscard]] inline double MetricTraceLengthScale(core::MetricId metric_id, double metric_mass) {
    SIRIUS_PRE(std::isfinite(metric_mass) && metric_mass >= 0.0);

    switch (metric_id) {
        case core::MetricId::Schwarzschild:
        case core::MetricId::Kerr:
        case core::MetricId::ReissnerNordstrom:
        case core::MetricId::KerrNewman:
        case core::MetricId::SchwarzschildDeSitter:
            SIRIUS_PRE(metric_mass > 0.0);
            return metric_mass;
        case core::MetricId::Minkowski:
        case core::MetricId::DeSitter:
        case core::MetricId::MorrisThorne:
        case core::MetricId::Alcubierre:
            return 1.0;
    }
    SIRIUS_ASSERT(false);
    return 1.0;
}

[[nodiscard]] inline TraceDomainParameters BuildTraceDomainParameters(core::MetricId metric_id,
                                                                      double metric_mass,
                                                                      double observer_radius) {
    SIRIUS_PRE(std::isfinite(observer_radius) && observer_radius > 0.0);

    // Metrics without a mass parameter retain the established unit numerical
    // scale; their own geometric parameters remain inside their metric authority.
    const double scale = MetricTraceLengthScale(metric_id, metric_mass);
    const double escape_radius = std::max(200.0 * scale, 1.25 * observer_radius);

    return {
        static_cast<float>(escape_radius),  static_cast<float>(0.1 * scale),
        static_cast<float>(1.0e-5 * scale), static_cast<float>(0.02 * scale),
        static_cast<float>(2.0 * scale),
    };
}

}  // namespace sirius::render
