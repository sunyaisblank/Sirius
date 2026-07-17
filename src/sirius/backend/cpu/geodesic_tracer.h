#pragma once

// Geodesic ray tracer for batch rendering: wraps the core RK45 geodesic
// integrator, converts camera rays from Boyer-Lindquist into Kerr-Schild
// Cartesian, and classifies each ray's termination (escape, horizon, disk hit).
// Ported from GTRC001A.h.
//
// Traces null geodesics d^2 x^mu/dlambda^2 + Gamma^mu_ab (dx^a/dlambda)
// (dx^b/dlambda) = 0 under the null constraint g_mu_nu k^mu k^nu = 0. The camera
// generates rays in BL coordinates; integration runs in Kerr-Schild Cartesian,
// and this class owns the transformation. Reference for the photon-orbit and
// image-order machinery: Bardeen, Press & Teukolsky (1972); Gralla, Lupsasca &
// Marrone (2020).

#include "sirius/core/camera.h"
#include "sirius/core/coordinates.h"
#include "sirius/core/geodesic_integrator.h"
#include "sirius/core/metrics/metric.h"
#include "sirius/core/tensor.h"

#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace sirius::backend {

// Outcome and per-ray data produced by a trace.
struct TraceResult {
    // Ray termination classification.
    enum class Outcome {
        Escaped,    // Escaped to infinity (sample background).
        Horizon,    // Captured by the black hole.
        DiskHit,    // Intersected the accretion disk.
        MaxSteps,   // Terminated at the step limit.
        Spiraling   // Detected as spiralling near the photon sphere.
    };

    Outcome outcome = Outcome::MaxSteps;

    // Final state (valid for Escaped).
    sirius::core::Vec4 final_position;
    sirius::core::Vec4 final_direction;

    // Disk-hit data (valid for DiskHit).
    float disk_radius = 0.0f;
    float disk_temperature = 0.0f;
    float disk_phi = 0.0f;

    // Motion-blur decomposition of the g-factor (valid for DiskHit); enables
    // analytic g(phi + dphi) = grav / (gamma (1 - v_orb (A cos dphi + B sin dphi))).
    float gfactor_grav = 1.0f;
    float gfactor_gamma = 1.0f;
    float gfactor_v_orb = 0.0f;
    float gfactor_A = 0.0f;
    float gfactor_B = 0.0f;

    // Diagnostics.
    int steps_taken = 0;
    float redshift = 1.0f;
    bool numerical_failure = false;

    // Magnification and photon ring.
    float magnification = 1.0f;
    float min_radius = 1e10f;
    bool photon_ring = false;

    // Higher-order imaging (winding number and image order).
    int equatorial_crossings = 0;
    float total_phi_change = 0.0f;
    int image_order = 0;
    float impact_parameter = 0.0f;

    // Track up to n = 3 images.
    static constexpr int kMaxDiskCrossings = 4;

    // Per-intersection disk crossing data.
    struct DiskCrossing {
        float r = 0.0f;
        float phi = 0.0f;
        float temperature = 0.0f;
        float redshift = 1.0f;
        int crossing_index = 0;
        bool valid = false;
    };

    DiskCrossing disk_crossings[kMaxDiskCrossings] = {};
    int num_disk_crossings = 0;

    // Volumetric disk accumulation (ray marching through a 3D disk).
    float volumetric_emission[3] = {0.0f, 0.0f, 0.0f};
    float optical_depth = 0.0f;
    float volumetric_path_length = 0.0f;
    bool volumetric_hit = false;
};

// Termination, disk, and integration parameters for the tracer.
struct TracerConfig {
    // Termination.
    float escape_radius = 100.0f;
    float horizon_factor = 1.05f;
    int max_steps = 5000;

    // Thin accretion disk.
    bool enable_disk = true;
    float disk_inner = 6.0f;   // ISCO for Schwarzschild = 6M.
    float disk_outer = 20.0f;
    float disk_thickness = 0.01f;
    float disk_temperature_inner = 1.0f;

    // Volumetric disk.
    bool enable_volumetric = false;
    float volumetric_H_over_r = 0.1f;
    float volumetric_H_power = 0.25f;
    float volumetric_tau_midplane = 10.0f;
    int volumetric_samples = 32;
    float volumetric_tau_max = 10.0f;

    // Integration (core RK45 configuration).
    sirius::core::IntegratorConfig integrator;

    TracerConfig() {
        integrator = sirius::core::Geodesic::GetDefaultConfig();
        integrator.abs_tolerance = 1e-6f;
        integrator.rel_tolerance = 1e-6f;
        integrator.min_step = 1e-6f;
        integrator.max_step = 0.1f;
        integrator.initial_step = 0.01f;
        integrator.use_rk45 = true;
    }
};

// Traces camera rays through a spacetime and reports each ray's fate.
class GeodesicTracer {
  public:
    // Construct with a metric (must outlive the tracer) and configuration.
    GeodesicTracer(sirius::core::IMetric* metric, const TracerConfig& config);

    ~GeodesicTracer() = default;

    // Trace one camera ray (BL origin, BL-basis direction). Postcondition:
    // result.outcome describes the termination condition.
    TraceResult Trace(const sirius::core::CameraRay& camera_ray);

    // Trace with per-call metric parameters (mass M, spin a, with spin passed as
    // an absolute a = spin so that a/M = spin/mass); updates the disk inner edge
    // to the ISCO for the new spin.
    TraceResult Trace(const sirius::core::CameraRay& camera_ray, double mass, double spin);

    void SetConfig(const TracerConfig& config) { config_ = config; }
    const TracerConfig& GetConfig() const { return config_; }
    void SetMetric(sirius::core::IMetric* metric) { metric_ = metric; }

  private:
    sirius::core::IMetric* metric_;
    TracerConfig config_;

    // Metric parameters cached once per trace.
    double cached_m_ = 1.0;
    double cached_a_ = 0.0;  // a/M.

    void CacheMetricParameters();

    // Convert a camera ray (BL) into a Cartesian Lightray, normalised to null.
    sirius::core::Lightray InitializeLightray(const sirius::core::CameraRay& camera_ray);

    // Whether the segment (pos_old, pos_new) crosses the equatorial disk; writes
    // the radius and azimuth of the crossing when it does.
    bool CheckDiskIntersection(const sirius::core::Vec4& pos_old,
                               const sirius::core::Vec4& pos_new, float& intersection_r,
                               float& intersection_phi);

    // Novikov-Thorne thin-disk temperature T(r) ~ r^(-3/4) normalised at the edge.
    float ComputeDiskTemperature(float r);

    // Whether any ray component is NaN/Inf.
    bool HasInvalidState(const sirius::core::Lightray& ray);

    // g-factor (gravitational + Doppler) for disk emission; intensity ~ g^4.
    float ComputeGFactor(float r, float phi, const sirius::core::Vec4& ray_vel);

    // g-factor with the decomposed components stored for analytic motion blur.
    void ComputeGFactorWithComponents(float r, float phi, const sirius::core::Vec4& ray_vel,
                                      TraceResult& result);

    // Keplerian angular velocity Omega(r) for the Kerr prograde orbit.
    float ComputeOrbitalVelocity(float r);

    // Volumetric disk helpers.
    bool IsInVolumetricDisk(float r, float z);
    float ComputeScaleHeight(float r);
    float ComputeVolumetricOpacityDensity(float r, float z);
    float ComputeVolumetricTemperature(float r, float z);
    void AccumulateVolumetricEmission(const sirius::core::Lightray& ray,
                                      const sirius::core::Vec4& entry_pos,
                                      const sirius::core::Vec4& exit_pos, TraceResult& result);
};

// ---- Inline implementations -------------------------------------------------

inline GeodesicTracer::GeodesicTracer(sirius::core::IMetric* metric, const TracerConfig& config)
    : metric_(metric), config_(config) {}

inline void GeodesicTracer::CacheMetricParameters() {
    const auto& params = metric_->GetParameters();
    cached_m_ = params.count("mass") ? params.at("mass").value : 1.0;
    cached_a_ = params.count("spin") ? params.at("spin").value : 0.0;
}

inline float GeodesicTracer::ComputeDiskTemperature(float r) {
    // Novikov-Thorne T(r) ~ r^(-3/4), normalised so T(disk_inner) = the inner
    // temperature.
    float r_in = config_.disk_inner;
    float T_in = config_.disk_temperature_inner;
    return T_in * std::pow(r_in / r, 0.75f);
}

inline bool GeodesicTracer::HasInvalidState(const sirius::core::Lightray& ray) {
    for (int i = 0; i < 4; ++i) {
        if (std::isnan(ray.position(i)) || std::isinf(ray.position(i))) return true;
        if (std::isnan(ray.velocity(i)) || std::isinf(ray.velocity(i))) return true;
    }
    return false;
}

inline bool GeodesicTracer::CheckDiskIntersection(const sirius::core::Vec4& pos_old,
                                                  const sirius::core::Vec4& pos_new,
                                                  float& intersection_r, float& intersection_phi) {
    // Equatorial-plane crossing (z = 0 in Cartesian).
    double z_old = pos_old(3);
    double z_new = pos_new(3);

    if (z_old * z_new >= 0) return false;

    // Linear interpolation to the crossing point.
    double alpha = -z_old / (z_new - z_old);

    double x_cross = pos_old(1) + alpha * (pos_new(1) - pos_old(1));
    double y_cross = pos_old(2) + alpha * (pos_new(2) - pos_old(2));

    double r_cross = std::sqrt(x_cross * x_cross + y_cross * y_cross);

    if (r_cross < config_.disk_inner || r_cross > config_.disk_outer) {
        return false;
    }

    intersection_r = static_cast<float>(r_cross);
    intersection_phi = static_cast<float>(std::atan2(y_cross, x_cross));

    return true;
}

}  // namespace sirius::backend
