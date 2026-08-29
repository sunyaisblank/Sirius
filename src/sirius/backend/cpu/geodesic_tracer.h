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
#include "sirius/core/disk/novikov_thorne_disk.h"
#include "sirius/core/disk/turbulence.h"
#include "sirius/core/geodesic_integrator.h"
#include "sirius/core/metrics/metric.h"
#include "sirius/core/observer_frame.h"
#include "sirius/core/polarisation/walker_penrose.h"
#include "sirius/core/relativistic_transfer.h"
#include "sirius/core/spectral/colour_modes.h"
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
        Escaped,  // Escaped to infinity (sample background).
        Horizon,  // Captured by the black hole.
        DiskHit,  // Intersected the accretion disk.
        MaxSteps  // Terminated at the step limit.
    };

    Outcome outcome = Outcome::MaxSteps;

    // Terminal event/state position for every outcome. final_direction is the
    // Eulerian sky direction and is valid only for Escaped.
    sirius::core::Vec4 final_position;
    sirius::core::Vec4 final_direction;

    // Disk-hit data (valid for DiskHit).
    float disk_radius = 0.0f;
    float disk_temperature = 0.0f;
    float disk_phi = 0.0f;

    // Exact first-crossing frequency transfers. full_disk_redshift uses the
    // circular geodesic emitter; zamo_disk_redshift is the locally non-rotating
    // diagnostic branch selected when orbital Doppler transfer is suppressed.
    float full_disk_redshift = 1.0f;
    float zamo_disk_redshift = 1.0f;

    // Diagnostics.
    int steps_taken = 0;
    float redshift = 1.0f;
    bool numerical_failure = false;

    // Exact minimum Kerr-Schild chart radius reached by the accepted samples.
    float min_radius = 1e10f;

    // An optically thick thin disk terminates the ray at its nearest crossing.
    static constexpr int kMaxDiskCrossings = 1;

    // Per-intersection disk crossing data.
    struct DiskCrossing {
        float r = 0.0f;
        float phi = 0.0f;
        float temperature = 0.0f;
        float redshift = 1.0f;
        float full_redshift = 1.0f;
        float zamo_redshift = 1.0f;
        float emission_cosine = 1.0f;
        int crossing_index = 0;
        bool valid = false;
        // Semi-infinite pure electron-scattering atmosphere, expressed in the
        // observer's transported screen basis. Valid only when the tracer's
        // polarisation path was enabled.
        float polarisation_evpa = 0.0f;
        float polarisation_degree = 0.0f;
        float polarisation_intensity_scale = 0.0f;
        bool polarisation_valid = false;
    };

    std::array<DiskCrossing, kMaxDiskCrossings> disk_crossings{};
    int num_disk_crossings = 0;

    // Volumetric disk accumulation (ray marching through a 3D disk).
    std::array<float, 3> volumetric_emission{};
    float optical_depth = 0.0f;
    // Observer-normalised affine length over samples containing extinction.
    float volumetric_affine_length = 0.0f;
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
    float disk_temperature_inner = 1.0f;
    float disk_temperature_scale_kelvin = 30000.0f;
    sirius::core::color_modes::Mode color_mode = sirius::core::color_modes::Mode::TrueColor;
    DiskTemperatureModel disk_temperature_model = DiskTemperatureModel::NovikovThorne;

    // Optional strong-field integration cap for accuracy-critical probes. When
    // both values are positive, only steps inside the named radius are capped;
    // asymptotic travel retains the ordinary adaptive maximum so escape remains
    // reachable within max_steps.
    float strong_field_radius = 0.0f;
    float strong_field_max_step = 0.0f;

    // Doppler transfer toggle (P4). True (the physical default) uses the
    // circular-emitter invariant frequency. False substitutes the exact ZAMO
    // branch as an explicitly nonphysical frame-dragging/gravitational
    // diagnostic inspired by the film treatment (James et al. 2015, section 5).
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
    // disk's semi-infinite pure electron-scattering atmosphere into it at the
    // observer-nearest thin-disk crossing. Off by default so ordinary renders
    // pay no transport cost.
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

    // Convert a camera ray (BL) into a Cartesian Lightray, normalised to null,
    // and return the same metric-orthonormal camera frame used for launch.
    sirius::core::Lightray InitializeLightray(
        const sirius::core::CameraRay& camera_ray,
        sirius::core::relativity::ObserverFrame* launch_frame);

    // Whether an accepted central-ray segment crosses the equatorial disk.
    // Solves z(lambda)=0 on the cubic Hermite segment and returns the crossing
    // event, tangent, affine fraction, radius, and azimuth.
    bool FindDiskIntersection(const sirius::core::Vec4& start_position,
                              const sirius::core::Vec4& start_tangent,
                              const sirius::core::Vec4& end_position,
                              const sirius::core::Vec4& end_tangent, double d_lambda,
                              float& intersection_r, float& intersection_phi,
                              double& intersection_fraction,
                              sirius::core::Vec4& intersection_position,
                              sirius::core::Vec4& intersection_tangent);

    // Novikov-Thorne thin-disk temperature T(r) ~ r^(-3/4) normalised at the edge.
    float ComputeDiskTemperature(float r);

    // Whether any ray component is NaN/Inf.
    bool HasInvalidState(const sirius::core::Lightray& ray);

    struct DiskFrequencySample {
        sirius::core::relativity::KerrDiskFrequencyTransfer transfer;
        float emission_cosine = 1.0f;
    };

    // Invariant g-factor (gravitational, frame-dragging, and orbital Doppler)
    // for disk emission; intensity ~ g^4. The observer frequency and the
    // photon's Killing quantities retain one affine scaling end to end.
    DiskFrequencySample ComputeDiskFrequencySample(float r, float cartesian_phi,
                                                   const sirius::core::Vec4& ray_vel,
                                                   float observer_frequency);
    void SetDiskFrequencyTransfer(float r, float phi, const sirius::core::Vec4& ray_vel,
                                  float observer_frequency, TraceResult& result,
                                  TraceResult::DiskCrossing& crossing);

    // Keplerian angular velocity Omega(r) for the Kerr prograde orbit.
    float ComputeOrbitalVelocity(float r);

    // Volumetric disk helpers.
    bool IsInVolumetricDisk(float r, float z);
    float ComputeScaleHeight(float r);
    float ComputeVolumetricOpacityDensity(float r, float z, float phi);
    float ComputeVolumetricTemperature(float r, float z);
    void AccumulateVolumetricEmission(const sirius::core::Vec4& entry_velocity,
                                      const sirius::core::Vec4& exit_velocity, double affine_length,
                                      float observer_frequency, const sirius::core::Vec4& entry_pos,
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

    // Initialise in the observer's Sachs screen, not a Euclidean chart plane.
    void InitBundle(const std::array<sirius::core::Vec4, 2>& launch_screen,
                    RayBundle& bundle) const;

    // Advance the bundle by one accepted affine step. RK4 evaluates connection
    // and curvature at its own start/midpoint/end stages; a cubic Hermite
    // interpolant of the accepted central ray supplies the midpoint event and
    // tangent without introducing a second geodesic authority.
    void StepBundle(const sirius::core::Vec4& start_position,
                    const sirius::core::Vec4& start_tangent, const sirius::core::Vec4& end_position,
                    const sirius::core::Vec4& end_tangent, double d_lambda, RayBundle& bundle);

    // Extract the beam ellipse (semi-axes, orientation, magnification) from the
    // final bundle, projected onto the plane transverse to the ray direction; the
    // affine length lambda converts the transverse extent to the angular sky
    // footprint for the pupil bundle.
    void FinaliseBundle(const RayBundle& bundle, const sirius::core::Vec4& position,
                        const sirius::core::Vec4& k, TraceResult::Beam& out) const;

    // --- Polarisation transport (E2) -----------------------------------------
    struct PolarisationFrame {
        sirius::core::PolarisedRay reference;
        sirius::core::PolarisedRay perpendicular;
    };

    void InitPolarisationFrame(const sirius::core::Lightray& ray,
                               const std::array<sirius::core::Vec4, 2>& launch_screen,
                               PolarisationFrame& frame);
    void AdvancePolarisationFrame(PolarisationFrame& frame, const sirius::core::Vec4& end_position,
                                  const sirius::core::Vec4& end_tangent, double d_lambda);
    void ReconditionPolarisationFrame(PolarisationFrame& frame, const sirius::core::Vec4& position,
                                      const sirius::core::Vec4& velocity);
    void SetDiskPolarisation(const PolarisationFrame& frame, TraceResult::DiskCrossing& crossing);
};

// ---- Inline implementations -------------------------------------------------

inline GeodesicTracer::GeodesicTracer(sirius::core::IMetric* metric, const TracerConfig& config)
    : metric_(metric), config_(config) {
    config_.turbulence.Validate();
}

inline void GeodesicTracer::CacheMetricParameters() {
    const auto& params = metric_->GetParameters();
    cached_m_ = params.count("mass") ? params.at("mass").value : 1.0;
    cached_a_ = params.count("spin") ? params.at("spin").value : 0.0;
    if (page_thorne_disk_ == nullptr || cached_m_ != page_thorne_cached_m_ ||
        cached_a_ != page_thorne_cached_a_) {
        sirius::core::AccretionDiskD::Config disk_config;
        // Only the dimensionless Page-Thorne radial shape is consumed below;
        // its physical flux amplitude is divided out at the reference radius.
        // Fixing M=1 makes that normalisation explicit and prevents the chart's
        // geometric mass scale from masquerading as a solar-mass input.
        disk_config.M = 1.0;
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

}  // namespace sirius::backend
