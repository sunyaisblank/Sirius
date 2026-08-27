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
#include "sirius/core/disk/corona.h"
#include "sirius/core/disk/novikov_thorne_disk.h"
#include "sirius/core/disk/turbulence.h"
#include "sirius/core/geodesic_integrator.h"
#include "sirius/core/metrics/metric.h"
#include "sirius/core/polarisation/walker_penrose.h"
#include "sirius/core/tensor.h"

#include <array>
#include <cmath>
#include <memory>
#include <numbers>

namespace sirius::backend {

enum class DiskTemperatureModel {
    NovikovThorne,
    ShakuraSunyaev,
};

// Outcome and per-ray data produced by a trace.
struct TraceResult {
    // Ray termination classification.
    enum class Outcome {
        Escaped,   // Escaped to infinity (sample background).
        Horizon,   // Captured by the black hole.
        DiskHit,   // Intersected the accretion disk.
        MaxSteps,  // Terminated at the step limit.
        Spiraling  // Detected as spiralling near the photon sphere.
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
    float gfactor_cosine_coefficient = 0.0f;
    float gfactor_sine_coefficient = 0.0f;

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
        // Thermal single-scattering disk polarisation, expressed in the
        // observer's transported screen basis. Valid only when the tracer's
        // polarisation path was enabled.
        float polarisation_evpa = 0.0f;
        float polarisation_degree = 0.0f;
        bool polarisation_valid = false;
    };

    std::array<DiskCrossing, kMaxDiskCrossings> disk_crossings{};
    int num_disk_crossings = 0;

    // Volumetric disk accumulation (ray marching through a 3D disk).
    std::array<float, 3> volumetric_emission{};
    float optical_depth = 0.0f;
    float volumetric_path_length = 0.0f;
    bool volumetric_hit = false;

    // Ray-bundle footprint (P2): the beam ellipse the geodesic-deviation solution
    // produces at escape or disk-hit, DNGR's distinguishing technique (James et
    // al. 2015, CQG 32 065001, section 3 and appendix B). Populated only when the
    // tracer runs with ray bundles enabled; otherwise valid stays false and the
    // point-sampled behaviour is unchanged. Angles are in the same units as the
    // initial bundle angular size.
    struct Beam {
        bool valid = false;
        float semi_major = 0.0f;       // Ellipse semi-major axis.
        float semi_minor = 0.0f;       // Ellipse semi-minor axis.
        float orientation = 0.0f;      // Position angle of the major axis [rad].
        float area_ratio = 1.0f;       // Final transverse area / initial area.
        float magnification = 1.0f;    // Initial / final area (oracle 1/|det J|).
        float transverse_area = 0.0f;  // Final transverse area (angular units^2).
        // Angular footprint on the celestial sphere (radians), the DNGR beam that
        // filters the star field (P3): the transverse extent divided by the
        // affine length converges to the ray-direction spread as the ray escapes.
        // Meaningful when the pupil (point-source) bundle mode is used.
        float footprint_major = 0.0f;
        float footprint_minor = 0.0f;
    };

    Beam beam;
};

// Termination, disk, and integration parameters for the tracer.
struct TracerConfig {
    // Termination.
    float escape_radius = 100.0f;
    // Kerr-Schild is horizon penetrating: the physical event horizon is the
    // capture boundary. Tests may opt into a larger diagnostic surface.
    float horizon_factor = 1.0f;
    int max_steps = 5000;

    // Thin accretion disk.
    bool enable_disk = true;
    float disk_inner = 6.0f;  // ISCO for Schwarzschild = 6M.
    float disk_outer = 20.0f;
    float disk_thickness = 0.01f;
    float disk_temperature_inner = 1.0f;
    DiskTemperatureModel disk_temperature_model = DiskTemperatureModel::NovikovThorne;

    // Cinematic work bound for rays that orbit near the photon separatrix.
    // Production may terminate these exponentially dim higher-order paths as
    // Spiraling; analytic shadow classifiers disable the heuristic because an
    // escaping near-critical ray is not a captured ray.
    bool enable_spiral_termination = true;

    // Optional strong-field integration cap for accuracy-critical probes. When
    // both values are positive, only steps inside the named radius are capped;
    // asymptotic travel retains the ordinary adaptive maximum so escape remains
    // reachable within max_steps.
    float strong_field_radius = 0.0f;
    float strong_field_max_step = 0.0f;

    // Doppler beaming toggle (P4). When true (default, full physics) the disk
    // g-factor carries the orbital Doppler asymmetry; when false the orbital
    // velocity is treated as zero for the brightness/colour factor while the
    // gravitational redshift is retained, mirroring the Interstellar artistic
    // choice (James et al. 2015, CQG 32 065001, section 5). True reproduces the
    // pinned render exactly. Reference: MTW section 22.3 (relativistic beaming).
    bool doppler_beaming = true;

    // Ray bundles (P2). When true the tracer propagates two geodesic-deviation
    // vectors alongside the central ray and reports the beam ellipse at
    // termination. Default false leaves the point-sampled path untouched (the
    // pinned render). bundle_angular_size is the initial half-extent of the
    // parallel bundle in the transverse plane (radians for a pixel footprint);
    // it cancels in the magnification and only sets the ellipse's absolute scale.
    bool enable_ray_bundles = false;
    float bundle_angular_size = 1.0e-3f;

    // Bundle initial condition. False (default) is a parallel bundle (identity
    // Jacobian, xi = orthonormal transverse pair, D xi / d lambda = 0), matching
    // the oracle's BeamStateD and giving the lensing magnification. True is a
    // pupil (point-source) bundle (xi = 0, D xi / d lambda = orthonormal), whose
    // transverse extent divided by affine length is the celestial-sphere
    // footprint the star filter needs (P3).
    bool bundle_point_source = false;

    // E2: transport an observer screen basis along the ray and project the
    // disk's single-scattering polarisation into it at every thin-disk crossing.
    // Off by default so ordinary renders pay no transport cost.
    bool enable_polarisation = false;

    // Volumetric disk.
    bool enable_volumetric = false;
    float volumetric_scale_height_ratio = 0.1f;
    float volumetric_flare_power = 0.25f;
    float volumetric_tau_midplane = 10.0f;
    int volumetric_samples = 32;
    float volumetric_tau_max = 10.0f;
    bool enable_turbulence = false;
    sirius::core::TurbulenceConfig turbulence;
    bool enable_corona = false;
    sirius::core::CoronaConfig corona;

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

    void SetConfig(const TracerConfig& config) {
        config_ = config;
        config_.turbulence.Validate();
        config_.corona.Validate();
    }
    const TracerConfig& GetConfig() const { return config_; }
    void SetMetric(sirius::core::IMetric* metric) { metric_ = metric; }

    // Geodesic-deviation (tidal) acceleration a^mu = -R^mu_nu_rho_sigma k^nu
    // xi^rho k^sigma in the tracer's Kerr-Schild Cartesian chart (P2). Exposed so
    // the beam physics can be validated against the oracle's analytic
    // Boyer-Lindquist Riemann tensor at a matched event: the contraction is a
    // genuine vector, so its metric-invariant magnitude is chart-independent.
    sirius::core::Vec4 TidalAcceleration(const sirius::core::Vec4& pos, const sirius::core::Vec4& k,
                                         const sirius::core::Vec4& xi);

    // Kretschmann scalar K = R_mu_nu_rho_sig R^mu_nu_rho_sig from the same
    // Cartesian Riemann the bundle rides on. A coordinate invariant, so it equals
    // the analytic value (48 M^2 / r^6 for Schwarzschild), which pins the
    // finite-differenced curvature against a chart-independent ground truth.
    double KretschmannScalar(const sirius::core::Vec4& pos);

  private:
    sirius::core::IMetric* metric_;
    TracerConfig config_;

    // Metric parameters cached once per trace.
    double cached_m_ = 1.0;
    double cached_a_ = 0.0;  // a/M.
    std::unique_ptr<sirius::core::AccretionDiskD> page_thorne_disk_;
    double page_thorne_reference_temperature_ = 0.0;
    double page_thorne_cached_m_ = -1.0;
    double page_thorne_cached_a_ = -2.0;

    void CacheMetricParameters();

    // Convert a camera ray (BL) into a Cartesian Lightray, normalised to null.
    sirius::core::Lightray InitializeLightray(const sirius::core::CameraRay& camera_ray);

    // Whether the segment (pos_old, pos_new) crosses the equatorial disk; writes
    // the radius and azimuth of the crossing when it does.
    bool CheckDiskIntersection(const sirius::core::Vec4& pos_old, const sirius::core::Vec4& pos_new,
                               float& intersection_r, float& intersection_phi);

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
    float ComputeVolumetricOpacityDensity(float r, float z, float phi);
    float ComputeVolumetricTemperature(float r, float z);
    void AccumulateVolumetricEmission(const sirius::core::Lightray& ray,
                                      const sirius::core::Vec4& entry_pos,
                                      const sirius::core::Vec4& exit_pos, TraceResult& result);

    // --- Ray-bundle (geodesic deviation) machinery (P2) ----------------------
    // A parallel bundle carried alongside the central ray: two deviation vectors
    // xi and their covariant derivatives V = D xi / d lambda, integrated by the
    // Jacobi equation D^2 xi^mu / d lambda^2 = -R^mu_nu_rho_sigma k^nu xi^rho
    // k^sigma (MTW eq 11.10; James et al. 2015 appendix B).
    struct RayBundle {
        sirius::core::Vec4 xi[2];  // Deviation vectors.
        sirius::core::Vec4 V[2];   // Covariant derivatives D xi / d lambda.
    };

    // Christoffel Gamma^mu_nu_rho in the tracer's Kerr-Schild Cartesian chart,
    // from the metric's analytic derivatives (single authority, no re-derivation).
    void ComputeChristoffelCart(const sirius::core::Vec4& pos, double Gamma[4][4][4]);

    // Riemann R^mu_nu_rho_sigma by central differences of the Christoffels in the
    // same chart the render integrates in, so k, xi, and R never mix charts
    // (mirrors kernels/gr_deviation.slang GetRiemannTensorCart).
    void ComputeRiemannCart(const sirius::core::Vec4& pos, double R[4][4][4][4]);

    // Initialise a parallel bundle at the ray start: xi_0, xi_1 an orthonormal
    // transverse pair scaled by the angular size, V = 0 (identity Jacobian).
    void InitBundle(const sirius::core::Vec4& k, RayBundle& bundle) const;

    // Advance the bundle by one affine step d_lambda, freezing Gamma and R at the
    // step's start point (they vary slowly over one adaptive step); RK4 on the
    // linear deviation system.
    void StepBundle(const sirius::core::Vec4& pos, const sirius::core::Vec4& k, double d_lambda,
                    RayBundle& bundle);

    // Extract the beam ellipse (semi-axes, orientation, magnification) from the
    // final bundle, projected onto the plane transverse to the ray direction; the
    // affine length lambda converts the transverse extent to the angular sky
    // footprint for the pupil bundle.
    void FinaliseBundle(const RayBundle& bundle, const sirius::core::Vec4& k, double lambda,
                        TraceResult::Beam& out) const;

    // --- Polarisation transport (E2) -----------------------------------------
    struct PolarisationFrame {
        sirius::core::PolarisedRay reference;
        sirius::core::PolarisedRay perpendicular;
    };

    void InitPolarisationFrame(const sirius::core::Lightray& ray, PolarisationFrame& frame);
    void AdvancePolarisationFrame(PolarisationFrame& frame, double d_lambda);
    void ReconditionPolarisationFrame(PolarisationFrame& frame, const sirius::core::Vec4& position,
                                      const sirius::core::Vec4& velocity);
    void SetDiskPolarisation(const PolarisationFrame& frame, TraceResult::DiskCrossing& crossing);
};

// ---- Inline implementations -------------------------------------------------

inline GeodesicTracer::GeodesicTracer(sirius::core::IMetric* metric, const TracerConfig& config)
    : metric_(metric), config_(config) {
    config_.turbulence.Validate();
    config_.corona.Validate();
}

inline void GeodesicTracer::CacheMetricParameters() {
    const auto& params = metric_->GetParameters();
    cached_m_ = params.count("mass") ? params.at("mass").value : 1.0;
    cached_a_ = params.count("spin") ? params.at("spin").value : 0.0;
    if (page_thorne_disk_ == nullptr || cached_m_ != page_thorne_cached_m_ ||
        cached_a_ != page_thorne_cached_a_) {
        sirius::core::AccretionDiskD::Config disk_config;
        disk_config.M = std::max(std::abs(cached_m_), 0.1);
        disk_config.a_star = cached_a_;
        disk_config.r_outer =
            std::max(500.0, static_cast<double>(config_.disk_outer) / std::max(cached_m_, 0.1));
        page_thorne_disk_ = std::make_unique<sirius::core::AccretionDiskD>(disk_config);
        const double reference_radius = 1.5 * page_thorne_disk_->IscoRadius();
        page_thorne_reference_temperature_ = page_thorne_disk_->Temperature(reference_radius);
        page_thorne_cached_m_ = cached_m_;
        page_thorne_cached_a_ = cached_a_;
    }
}

inline float GeodesicTracer::ComputeDiskTemperature(float r) {
    float r_in = config_.disk_inner;
    float T_in = config_.disk_temperature_inner;
    if (config_.disk_temperature_model == DiskTemperatureModel::ShakuraSunyaev) {
        return T_in * std::pow(r_in / r, 0.75f);
    }

    // Full Page-Thorne flux authority, normalised at 1.5 r_ISCO so
    // disk_temperature_inner remains the operator's Kelvin scale rather than
    // silently becoming an accretion-rate input.
    if (page_thorne_disk_ == nullptr || page_thorne_reference_temperature_ <= 0.0 ||
        cached_m_ <= 0.0) {
        return 0.0f;
    }
    const double model_temperature =
        page_thorne_disk_->Temperature(static_cast<double>(r) / cached_m_);
    return T_in * static_cast<float>(std::max(model_temperature, 0.0) /
                                     page_thorne_reference_temperature_);
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

    sirius::core::coordinates::Vec4Cart crossing_cart{0.0, x_cross, y_cross, 0.0};
    const double absolute_spin = cached_a_ * cached_m_;
    double r_cross = sirius::core::coordinates::KerrSchildRadius(crossing_cart, absolute_spin);

    if (r_cross < config_.disk_inner || r_cross > config_.disk_outer) {
        return false;
    }

    intersection_r = static_cast<float>(r_cross);
    intersection_phi = static_cast<float>(std::atan2(y_cross, x_cross));

    return true;
}

}  // namespace sirius::backend
