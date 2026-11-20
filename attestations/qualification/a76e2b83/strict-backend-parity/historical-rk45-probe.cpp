#include "sirius/backend/cpu/geodesic_tracer.h"
#include "sirius/core/camera.h"
#include "sirius/core/metrics/kerr_schild_family.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <numbers>

using sirius::backend::GeodesicTracer;
using sirius::backend::TraceResult;
using sirius::backend::TracerConfig;
using namespace sirius::core;

namespace {

struct RejectionObservation {
    int retries = 0;
    int changed_phase = 0;
};

RejectionObservation observation;

// At the equatorial observer's origin (radius, 0, 0), the conserved
// angular momentum is radius times transverse velocity. The positive
// conserved energy is the discriminant root of the null condition.
double impact_parameter(const CameraRay& ray) {
    const double radius = ray.origin(1);
    const double radial = ray.direction(1);
    const double transverse = std::hypot(ray.direction(2), ray.direction(3));
    const double h = 2.0 / radius;
    const double a = h - 1.0;
    const double b = h * radial;
    const double c = (1.0 + h) * radial * radial + transverse * transverse;
    return radius * transverse / std::sqrt(b * b - a * c);
}

const char* outcome(TraceResult::Outcome result) {
    switch (result) {
        case TraceResult::Outcome::Escaped: return "escaped";
        case TraceResult::Outcome::Horizon: return "captured";
        case TraceResult::Outcome::DiskHit: return "disk";
        case TraceResult::Outcome::MaxSteps: return "max_steps";
        case TraceResult::Outcome::Spiraling: return "spiraling";
    }
    return "unknown";
}

}  // namespace

extern "C" bool real_step(Lightray&, IMetric*, const IntegratorConfig&)
    asm("__real_" SIRIUS_RK45_LINK_SYMBOL);
extern "C" bool observed_step(Lightray&, IMetric*, const IntegratorConfig&)
    asm("__wrap_" SIRIUS_RK45_LINK_SYMBOL);

bool observed_step(Lightray& ray, IMetric* metric, const IntegratorConfig& config) {
    const auto position = ray.position;
    const auto velocity = ray.velocity;
    const float step = ray.step_size;
    const bool accepted = real_step(ray, metric, config);
    if (!accepted && !ray.terminated) {
        ++observation.retries;
        bool unchanged = ray.step_size < step;
        for (int axis = 0; axis < 4; ++axis) {
            unchanged &= position(axis) == ray.position(axis) && velocity(axis) == ray.velocity(axis);
        }
        if (!unchanged) ++observation.changed_phase;
    }
    return accepted;
}

int main() {
    KerrSchildParams parameters;
    parameters.M = 1.0;
    parameters.a = 0.0;
    KerrSchildFamily metric(parameters);
    CameraConfig camera_config;
    camera_config.r = 50.0;
    camera_config.theta = std::numbers::pi / 2.0;
    camera_config.phi = 0.0;
    camera_config.fov = 60.0f;
    camera_config.width = 64;
    camera_config.height = 64;
    PinholeCamera camera(camera_config);

    std::cout << std::setprecision(17);
    std::cout << "profile,x,y,impact,analytic,outcome,attempts,numerical_failure,rejections,changed_phase\n";
    for (const bool stress : {false, true}) {
        TracerConfig config;
        config.escape_radius = 100.0f;
        config.horizon_factor = 1.05f;
        config.max_steps = 5000;
        config.enable_disk = false;
        config.integrator.abs_tolerance = 1.0e-6f;
        config.integrator.rel_tolerance = 1.0e-6f;
        if (stress) {
            config.integrator.initial_step = 8.0f;
            config.integrator.max_step = 8.0f;
        }
        GeodesicTracer tracer(&metric, config);
        const double critical = std::sqrt(27.0);
        for (int x = 20; x < 44; ++x) {
            constexpr int y = 32;
            const CameraRay ray = camera.GenerateRay(x, y);
            const double impact = impact_parameter(ray);
            if (impact < 2.5 || impact > 8.0 || std::abs(impact - critical) < 0.2) continue;
            observation = {};
            const auto trace = tracer.Trace(ray);
            std::cout << (stress ? "rejection_stress" : "default_step") << ',' << x << ',' << y << ','
                      << impact << ',' << (impact < critical ? "captured" : "escaped") << ','
                      << outcome(trace.outcome) << ',' << trace.steps_taken << ','
                      << trace.numerical_failure << ',' << observation.retries << ','
                      << observation.changed_phase << '\n';
        }
    }
}
