// Render session implementation. Ported from SNRS001A.cpp.
//
// The CPU render path (metric construction, spiral tile scheduling, the
// thread-pool, per-pixel shading, and output writing) is preserved exactly. The
// legacy OptiX/GPU branches are removed; each removal site carries a one-line
// decision comment. The Vulkan compute backend will re-enter through
// sirius::backend::device, not from here.

#include "sirius/render/session/render_session.h"

#include "sirius/base/resource_locator.h"
#include "sirius/render/exr_writer.h"
#include "sirius/render/film_pipeline.h"
#include "sirius/render/image_buffer.h"
#include "sirius/render/pixel_sampling.h"
#include "sirius/render/png_writer.h"
#include "sirius/render/trace_domain.h"

#ifdef SIRIUS_HAS_VULKAN_BACKEND
#include "sirius/backend/device.h"  // Device enumeration for backend auto-resolution.
#include "sirius/render/vulkan_renderer.h"
#endif

#include "sirius/core/constants.h"
#include "sirius/core/metrics/cpu_metric_factory.h"
#include "sirius/core/spectral/blackbody.h"  // Spectral blackbody colour.

// stb_image decoder for the starfield texture; the implementation lives once in
// stb_impl.cpp within this same library.
#include "stb_image.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <chrono>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <locale>
#include <memory>
#include <sstream>
#include <stdexcept>

namespace sirius::render {

using backend::GeodesicTracer;
using backend::TracerConfig;
using backend::TraceResult;
using core::AccretionDiskD;
using core::CameraConfig;
using core::CameraRay;
using core::MetricId;
using core::Vec4;

namespace math = core::constants::math;

namespace {
constexpr int kBloomRadius = 12;
constexpr float kShadowLift = 0.02f;

core::MetricConstructionParameters MetricConstructionParametersFor(const SessionConfig& config) {
    core::MetricConstructionParameters parameters;
    parameters.mass = config.black_hole_mass;
    parameters.dimensionless_spin = config.black_hole_spin;
    parameters.dimensionless_charge = config.black_hole_charge;
    parameters.cosmological_constant = config.cosmological_constant;
    parameters.throat_radius = config.throat_radius;
    parameters.wormhole_topology = config.wormhole_topology;
    parameters.warp_velocity = config.warp_velocity;
    parameters.bubble_radius = config.bubble_radius;
    parameters.bubble_sigma = config.bubble_sigma;
    return parameters;
}
}  // namespace

RenderSession::~RenderSession() {
    (void)Cancel();
    WaitForCompletion();
}

base::Expected<void> RenderSession::Configure(const SessionConfig& config) {
    std::lock_guard<std::mutex> lock(lifecycle_mutex_);
    if (fsm_.GetState() != SessionState::Idle || render_thread_.joinable() || join_in_progress_) {
        return base::Fail(base::ErrorDomain::kConfiguration, "configure render session",
                          "session is not idle");
    }
    if (const auto issue = SessionConfigIssue(config); issue.has_value()) {
        return base::Fail(base::ErrorDomain::kConfiguration, "configure render session", *issue);
    }
    config_ = config;
    return {};
}

bool RenderSession::Start() {
    std::lock_guard<std::mutex> lock(lifecycle_mutex_);
    if (fsm_.GetState() != SessionState::Idle || render_thread_.joinable() || join_in_progress_) {
        return false;
    }
    progress_.GetCancellationToken().Reset();
    stop_workers_ = false;
    try {
        render_thread_ = std::thread([this] {
            try {
                if (!fsm_.Process(SessionEvent::Start)) {
                    error_message_ = "render session could not enter Initialising";
                    fsm_.ForceState(SessionState::Failed);
                    OnSessionEnd(SessionState::Failed);
                }
            } catch (const std::exception& error) {
                error_message_ = std::string("unhandled render exception: ") + error.what();
                if (!IsTerminal(fsm_.GetState())) {
                    fsm_.ForceState(SessionState::Failed);
                    OnSessionEnd(SessionState::Failed);
                }
            } catch (...) {
                error_message_ = "unhandled non-standard render exception";
                if (!IsTerminal(fsm_.GetState())) {
                    fsm_.ForceState(SessionState::Failed);
                    OnSessionEnd(SessionState::Failed);
                }
            }
        });
    } catch (...) {
        return false;
    }
    return true;
}

bool RenderSession::Cancel() {
    const SessionState state = fsm_.GetState();
    if (IsTerminal(state) || state == SessionState::Completing) {
        return false;
    }
    if (state == SessionState::Idle) {
        std::lock_guard<std::mutex> lock(lifecycle_mutex_);
        if (!render_thread_.joinable() && !join_in_progress_) {
            return false;
        }
    }
    progress_.GetCancellationToken().Cancel();
    stop_workers_ = true;
    return true;
}

void RenderSession::WaitForCompletion() {
    std::unique_lock<std::mutex> lock(lifecycle_mutex_);
    const std::thread::id caller = std::this_thread::get_id();
    if (join_in_progress_) {
        if (caller == render_thread_id_) {
            return;
        }
        lifecycle_cv_.wait(lock, [this] { return !join_in_progress_; });
        return;
    }
    if (!render_thread_.joinable() || render_thread_.get_id() == caller) {
        return;
    }

    render_thread_id_ = render_thread_.get_id();
    std::thread completed_thread = std::move(render_thread_);
    join_in_progress_ = true;
    lock.unlock();
    completed_thread.join();
    lock.lock();
    join_in_progress_ = false;
    render_thread_id_ = {};
    lock.unlock();
    lifecycle_cv_.notify_all();
}

std::string SessionSceneEvidenceJson(const SessionConfig& config, std::size_t point_star_count) {
    const char* backend = nullptr;
    switch (config.backend) {
        case RenderBackend::Cpu:
            backend = "Cpu";
            break;
        case RenderBackend::Vulkan:
            backend = "Vulkan";
            break;
        default:
            SIRIUS_ASSERT(false);
            backend = "Invalid";
            break;
    }

    const char* lens = nullptr;
    switch (config.lens_type) {
        case core::LensType::Pinhole:
            lens = "Pinhole";
            break;
        case core::LensType::ThinLens:
            lens = "ThinLens";
            break;
        case core::LensType::Fisheye:
            lens = "Fisheye";
            break;
        default:
            SIRIUS_ASSERT(false);
            lens = "Invalid";
            break;
    }

    std::ostringstream evidence;
    evidence.imbue(std::locale::classic());
    evidence << std::setprecision(std::numeric_limits<double>::max_digits10)
             << "{\"schema\":\"sirius-render-scene-v1\"" << ",\"backend\":\"" << backend << "\""
             << ",\"metric\":\"" << core::MetricInfoFor(config.metric_id).canonical_name << "\""
             << ",\"spin\":" << config.black_hole_spin << ",\"width\":" << config.width
             << ",\"height\":" << config.height
             << ",\"samples_per_pixel\":" << config.samples_per_pixel
             << ",\"field_of_view\":" << static_cast<double>(config.camera_fov)
             << ",\"disk_enabled\":" << (config.enable_disk ? "true" : "false")
             << ",\"ray_bundles\":" << (config.ray_bundles ? "true" : "false")
             << ",\"point_starfield\":" << (config.point_starfield ? "true" : "false")
             << ",\"point_star_count\":" << point_star_count << ",\"point_brightness_scale\":"
             << static_cast<double>(config.point_starfield_config.brightness_scale)
             << ",\"camera_beta\":[" << config.camera_beta_forward << ',' << config.camera_beta_up
             << ',' << config.camera_beta_right << "]" << ",\"lens\":\"" << lens << "\""
             << ",\"focal_length\":" << static_cast<double>(config.camera_focal_length)
             << ",\"aperture\":" << static_cast<double>(config.camera_aperture)
             << ",\"focus_distance\":" << static_cast<double>(config.camera_focus_distance) << '}';
    return evidence.str();
}

std::optional<std::string> SessionConfigIssue(const SessionConfig& config) {
    const SessionConfig defaults;
    const auto finite = [](double value) { return std::isfinite(value); };
    const auto in_range = [&finite](double value, double minimum, double maximum) {
        return finite(value) && value >= minimum && value <= maximum;
    };

    if (config.width < 1 || config.width > 8192 || config.height < 1 || config.height > 8192) {
        return "width and height must each be between 1 and 8192";
    }
    if (config.tile_size < 1 || config.tile_size > 256 ||
        (config.tile_size & (config.tile_size - 1)) != 0) {
        return "tile size must be a power of two between 1 and 256";
    }
    if (config.samples_per_pixel < 1 || config.samples_per_pixel > 4096) {
        return "samples per pixel must be between 1 and 4096";
    }
    if (config.ray_bundles && config.metric_id != MetricId::Minkowski &&
        config.metric_id != MetricId::Schwarzschild && config.metric_id != MetricId::Kerr) {
        return "ray-bundle curvature transport is represented only for Minkowski, Schwarzschild, "
               "and Kerr";
    }
    if (config.thread_count < 0 || config.thread_count > 1024) {
        return "thread count must be between 0 and 1024";
    }
    if (!config.enable_parallel_rendering && config.thread_count != 0) {
        return "thread count requires parallel CPU rendering";
    }
    if (config.point_starfield) {
        if (!core::IsRepresentedPointStarfieldConfig(config.point_starfield_config)) {
            return "point-starfield parameters are outside the represented domain";
        }
    } else if (config.point_starfield_config != defaults.point_starfield_config) {
        return "point-starfield parameters require point-starfield mode";
    }
    if (config.write_output) {
        if (config.output_path.empty() || config.output_path.find('\0') != std::string::npos) {
            return "output path must not be empty or contain a null byte";
        }
        std::string extension;
        if (const auto dot = config.output_path.rfind('.'); dot != std::string::npos) {
            extension = config.output_path.substr(dot);
            std::transform(extension.begin(), extension.end(), extension.begin(),
                           [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
        }
        if (extension != ".ppm" && extension != ".png" && extension != ".exr") {
            return "output path extension must be .ppm, .png, or .exr";
        }
    } else if (config.output_path != defaults.output_path) {
        return "output path requires file output";
    }

    switch (config.backend) {
        case RenderBackend::Cpu:
        case RenderBackend::Vulkan:
            break;
        default:
            return "invalid render backend";
    }
    if (config.backend == RenderBackend::Vulkan &&
        (config.tile_size != defaults.tile_size || config.thread_count != defaults.thread_count ||
         config.enable_parallel_rendering != defaults.enable_parallel_rendering)) {
        return "CPU tile and thread controls are not represented by the device-governed Vulkan "
               "backend";
    }
    switch (config.wormhole_topology) {
        case WormholeTopology::OneSheetCapture:
        case WormholeTopology::TwoSheet:
            break;
        default:
            return "invalid wormhole topology";
    }
    const core::MetricConstructionParameters metric_parameters =
        MetricConstructionParametersFor(config);
    if (const auto issue =
            core::MetricConstructionParameterIssue(config.metric_id, metric_parameters);
        issue.has_value()) {
        return std::string(*issue);
    }

    const core::MetricObserverRadiusIssue observer_radius_issue =
        core::MetricObserverRadiusIssueFor(config.metric_id, config.black_hole_mass,
                                           config.cosmological_constant, config.observer_distance,
                                           config.throat_radius, config.bubble_radius,
                                           config.bubble_sigma);
    if (observer_radius_issue == core::MetricObserverRadiusIssue::CosmologicalHorizon) {
        return "positive-lambda observer must remain inside the cosmological trace boundary";
    }
    if (observer_radius_issue != core::MetricObserverRadiusIssue::None ||
        !finite(config.observer_inclination) ||
        config.observer_inclination <= 0.1 * math::kPi / 180.0 ||
        config.observer_inclination >= 179.9 * math::kPi / 180.0 ||
        !finite(config.observer_azimuth) || std::abs(config.observer_azimuth) > 2.0 * math::kPi ||
        !finite(config.camera_fov) || config.camera_fov < 1.0f || config.camera_fov > 170.0f) {
        return "observer distance, angles, or field of view is outside the represented domain";
    }
    const double bundle_angular_size =
        (static_cast<double>(config.camera_fov) * math::kPi / 180.0) /
        static_cast<double>(config.height);
    if (config.ray_bundles && bundle_angular_size > 0.1) {
        return "ray-bundle angular extent exceeds the represented linear-deviation domain";
    }
    const double beta_squared = config.camera_beta_forward * config.camera_beta_forward +
                                config.camera_beta_up * config.camera_beta_up +
                                config.camera_beta_right * config.camera_beta_right;
    if (!finite(beta_squared) || beta_squared >= 1.0) {
        return "camera beta magnitude must be finite and below one";
    }
    switch (config.lens_type) {
        case core::LensType::Pinhole:
        case core::LensType::ThinLens:
        case core::LensType::Fisheye:
            break;
        default:
            return "invalid lens model";
    }
    if (!in_range(config.camera_focal_length, std::numeric_limits<float>::min(), 10000.0) ||
        !in_range(config.camera_aperture, std::numeric_limits<float>::min(), 128.0) ||
        !in_range(config.camera_focus_distance, std::numeric_limits<float>::min(), 1.0e6)) {
        return "camera focal length, aperture, or focus distance is outside the represented domain";
    }
    if (const auto issue =
            core::LensSpecificParameterIssue(config.lens_type, config.camera_focal_length,
                                             config.camera_aperture, config.camera_focus_distance);
        issue.has_value()) {
        return std::string(*issue);
    }
    if (const auto issue = core::ThinLensGeometryIssue(
            config.lens_type, config.observer_distance, config.camera_focal_length,
            config.camera_aperture, config.camera_focus_distance);
        issue.has_value()) {
        return std::string(*issue);
    }

    switch (config.temperature_model) {
        case DiskTemperatureModel::NovikovThorne:
        case DiskTemperatureModel::ShakuraSunyaev:
            break;
        default:
            return "invalid disk temperature model";
    }
    if (!finite(config.disk_temperature_scale) || config.disk_temperature_scale < 100.0f ||
        config.disk_temperature_scale > 1.0e8f) {
        return "disk temperature scale must be between 100 and 100000000 Kelvin";
    }
    if (config.enable_disk &&
        core::DiskSupportFor(config.metric_id) != core::DiskSupport::PageThorne) {
        return "the selected metric has no represented accretion-disk emission model";
    }
    if (!config.enable_disk &&
        (config.temperature_model != DiskTemperatureModel::NovikovThorne ||
         config.disk_temperature_scale != core::kDefaultDiskTemperatureKelvin)) {
        return "disk temperature model and scale require disk emission";
    }
    if (!config.enable_disk && !config.doppler_beaming) {
        return "Doppler-beaming control requires disk emission";
    }
    if (config.enable_corona) {
        return "corona emission requires frequency-dependent covariant Compton transfer, which "
               "is not represented";
    }
    if (!config.enable_disk && (config.enable_volumetric_disk || config.enable_turbulence)) {
        return "volumetric disk and turbulence require the disk";
    }
    if (!config.enable_volumetric_disk && config.enable_turbulence) {
        return "turbulence requires volumetric transfer";
    }
    if (config.volumetric_samples < 1 || config.volumetric_samples > 4096 ||
        !finite(config.volumetric_h_over_r) || config.volumetric_h_over_r < 0.01f ||
        config.volumetric_h_over_r > 0.5f || !finite(config.volumetric_h_power) ||
        config.volumetric_h_power < -2.0f || config.volumetric_h_power > 4.0f ||
        !finite(config.volumetric_tau_midplane) || config.volumetric_tau_midplane < 0.0f ||
        config.volumetric_tau_midplane > 1.0e6f) {
        return "volumetric transfer parameters are outside the represented domain";
    }
    if (!config.enable_volumetric_disk &&
        (config.volumetric_h_over_r != defaults.volumetric_h_over_r ||
         config.volumetric_h_power != defaults.volumetric_h_power ||
         config.volumetric_tau_midplane != defaults.volumetric_tau_midplane ||
         config.volumetric_samples != defaults.volumetric_samples)) {
        return "volumetric parameters require volumetric transfer";
    }

    switch (config.color_mode) {
        case core::color_modes::Mode::TrueColor:
        case core::color_modes::Mode::TemperatureMap:
        case core::color_modes::Mode::RedshiftMap:
        case core::color_modes::Mode::Polarisation:
            break;
        default:
            return "invalid colour mode";
    }
    const bool polarisation_mode = config.color_mode == core::color_modes::Mode::Polarisation;
    if (!config.enable_disk && config.color_mode != core::color_modes::Mode::TrueColor) {
        return "diagnostic colour modes require disk emission";
    }
    if (config.enable_polarisation != polarisation_mode) {
        return "polarisation transport and colour mode must be enabled together";
    }
    if (polarisation_mode) {
        if (!config.enable_disk) {
            return "polarisation requires the thin accretion disk";
        }
        if (config.enable_volumetric_disk) {
            return "polarisation is not represented for volumetric transfer";
        }
        if (config.metric_id != MetricId::Schwarzschild && config.metric_id != MetricId::Kerr) {
            return "polarisation is represented only for Schwarzschild and Kerr";
        }
        if (config.enable_motion_blur) {
            return "polarisation is not represented with temporal disk motion blur";
        }
    }
    switch (config.tonemapper) {
        case core::TonemapType::None:
        case core::TonemapType::Reinhard:
        case core::TonemapType::AcesFit:
        case core::TonemapType::Filmic:
        case core::TonemapType::Exposure:
            break;
        default:
            return "invalid tonemapper";
    }
    if (!in_range(config.exposure, std::numeric_limits<float>::min(), 100.0) ||
        !in_range(config.bloom_intensity, 0.0, 5.0) ||
        !in_range(config.bloom_threshold, 0.0, 100.0) || !in_range(config.contrast, 0.0, 4.0) ||
        !in_range(config.saturation, 0.0, 4.0)) {
        return "display-pipeline parameters are outside the represented domain";
    }
    if (!config.enable_bloom && (config.bloom_intensity != defaults.bloom_intensity ||
                                 config.bloom_threshold != defaults.bloom_threshold)) {
        return "bloom intensity and threshold require bloom";
    }
    if (!config.enable_motion_blur &&
        (config.shutter_time != defaults.shutter_time ||
         config.motion_blur_samples != defaults.motion_blur_samples)) {
        return "motion-blur parameters require temporal integration";
    }
    if (config.enable_motion_blur) {
        if (!in_range(config.shutter_time, 0.0, 1000.0) || config.motion_blur_samples < 1 ||
            config.motion_blur_samples > 4096) {
            return "motion-blur parameters are outside the represented domain";
        }
        return "temporal disk motion blur requires a represented time-dependent emissivity model";
    }
    if (config.enable_film_finish) {
        if (!IsRepresentedFilmConfig(config.film_config)) {
            return "film-finish parameters are outside the represented domain";
        }
    } else if (config.film_config != defaults.film_config) {
        return "film-finish parameters require the film finish";
    }
    if (config.enable_jets) {
        return "relativistic jets require covariant geodesic radiative transfer, which is not "
               "represented";
    }
    return std::nullopt;
}

// =============================================================================
// FSM wiring.
// =============================================================================
void RenderSession::SetupActions() {
    fsm_.SetEntryAction([this](SessionState state) { OnEnterState(state); });

    fsm_.SetTransitionAction([](SessionState from, SessionEvent event, SessionState to) {
        std::cout << "[Session] " << StateName(from) << " --" << EventName(event) << "--> "
                  << StateName(to) << std::endl;
    });
}

void RenderSession::OnEnterState(SessionState state) {
    switch (state) {
        case SessionState::Initialising: {
            auto initialised = Initialise();
            if (!initialised) {
                error_message_ = initialised.error().Description();
                fsm_.Process(progress_.GetCancellationToken().IsCancelled() ? SessionEvent::Cancel
                                                                            : SessionEvent::Error);
            }
            break;
        }
        case SessionState::Scheduling:
            ScheduleNextTile();
            break;
        case SessionState::Rendering:
            // Tile is being processed.
            break;
        case SessionState::Completing: {
            auto written = WriteOutput();
            if (!written) {
                error_message_ = written.error().Description();
                fsm_.Process(SessionEvent::Error);
            }
            break;
        }
        case SessionState::Complete:
        case SessionState::Failed:
        case SessionState::Cancelled:
            OnSessionEnd(state);
            break;
        case SessionState::Idle:
            break;
        default:
            SIRIUS_ASSERT(false);
            break;
    }
}

// =============================================================================
// Initialisation.
// =============================================================================
base::Expected<void> RenderSession::Initialise() {
    if (progress_.GetCancellationToken().IsCancelled()) {
        fsm_.Process(SessionEvent::Cancel);
        return {};
    }
    std::cout << "[Session] Initialising render..." << std::endl;
    std::cout << "  Resolution: " << config_.width << " x " << config_.height << std::endl;
    if (config_.backend == RenderBackend::Cpu) {
        std::cout << "  Tile size:  " << config_.tile_size << std::endl;
    }
    std::cout << "  Samples:    " << config_.samples_per_pixel << std::endl;

    // CPU tiles are operator-sized and spiral ordered. Vulkan derives a separate
    // device-budgeted tile plan and reports its actual total at dispatch time.
    if (config_.backend == RenderBackend::Cpu) {
        tiles_.Initialise(config_.width, config_.height, config_.tile_size);
        std::cout << "  Tiles:      " << tiles_.GetTileCount() << " (spiral order)" << std::endl;
    } else {
        std::cout << "  Tiles:      device-budget governed (Vulkan)" << std::endl;
    }

    // Display buffer.
    display_.Initialise(config_.width, config_.height);

    // Progress tracker (Start() before SetTotals to avoid a reset).
    progress_.Start();
    if (config_.backend == RenderBackend::Cpu) {
        progress_.SetTotals(tiles_.GetTileCount(), config_.samples_per_pixel);
    } else {
        progress_.SetTotals(1, 1);
    }

    // The Vulkan renderer owns device scene construction, resource loading,
    // catalogue upload, and dispatch. Do not construct an unused CPU metric,
    // camera, tracer pool, texture, or duplicate point catalogue first.
    if (config_.backend == RenderBackend::Vulkan) {
        const std::size_t point_star_count =
            config_.point_starfield ? config_.point_starfield_config.star_count : 0;
        std::cout << "[Session] Scene evidence: "
                  << SessionSceneEvidenceJson(config_, point_star_count) << std::endl;
        if (progress_.GetCancellationToken().IsCancelled()) {
            fsm_.Process(SessionEvent::Cancel);
            return {};
        }
        fsm_.Process(SessionEvent::Ready);
        return {};
    }

    // Construct through the same authority probed by the non-render backend
    // estate. Charge and spin remain dimensionless at this boundary; the
    // factory performs their one owned conversion to geometric lengths.
    const core::MetricConstructionParameters metric_parameters =
        MetricConstructionParametersFor(config_);
    metric_ = core::CreateCpuMetric(config_.metric_id, metric_parameters);
    if (metric_) {
        std::cout << "  Metric:     " << metric_->GetName() << std::endl;
    } else {
        std::cout << "  Metric:     " << core::MetricInfoFor(config_.metric_id).canonical_name
                  << " (GPU backend only)" << std::endl;
    }

    // Camera.
    CameraConfig camera_config;
    camera_config.r = config_.observer_distance;
    camera_config.theta = config_.observer_inclination;
    camera_config.phi = config_.observer_azimuth;
    camera_config.fov = config_.camera_fov;
    camera_config.width = config_.width;
    camera_config.height = config_.height;
    // Camera worldline aberration (P5): beta components in the local
    // ray-component frame. Zero keeps every ray unaberrated (byte-pin).
    camera_config.beta_x = config_.camera_beta_forward;
    camera_config.beta_y = config_.camera_beta_up;
    camera_config.beta_z = config_.camera_beta_right;
    camera_config.focal_length = config_.camera_focal_length;
    camera_config.aperture = config_.camera_aperture;
    camera_config.focus_distance = config_.camera_focus_distance;
    camera_ = core::CreateCamera(config_.lens_type, camera_config);
    std::cout << "  Observer:   r=" << camera_config.r << " coordinate units";
    if (core::MetricUsesMass(config_.metric_id)) {
        std::cout << " (r/M=" << camera_config.r / config_.black_hole_mass << ')';
    } else if (config_.metric_id == MetricId::MorrisThorne) {
        std::cout << " (rho/b0=" << camera_config.r / config_.throat_radius << ')';
    } else if (config_.metric_id == MetricId::Alcubierre) {
        const double scale = core::MetricSceneLengthScale(
            config_.metric_id, config_.black_hole_mass, config_.throat_radius,
            config_.bubble_radius, config_.bubble_sigma);
        std::cout << " (r/L=" << camera_config.r / scale << ", L=max(R,1/sigma))";
    }
    std::cout << ", theta=" << (camera_config.theta * 180.0 / math::kPi) << " deg" << std::endl;

    // Geodesic tracer.
    TracerConfig tracer_config;
    const TraceDomainParameters trace_domain = BuildTraceDomainParameters({
        .metric_id = config_.metric_id,
        .metric_mass = config_.black_hole_mass,
        .cosmological_constant = config_.cosmological_constant,
        .observer_radius = config_.observer_distance,
        .throat_radius = config_.throat_radius,
        .bubble_radius = config_.bubble_radius,
        .bubble_sigma = config_.bubble_sigma,
    });
    tracer_config.escape_radius = trace_domain.escape_radius;
    tracer_config.finite_causal_boundary = trace_domain.finite_causal_boundary;
    // Kerr-Schild coordinates are horizon-penetrating, so the exact capture
    // surface is numerically safe. Enlarging it inflates the near-extremal
    // shadow.
    tracer_config.horizon_factor = 1.0f;
    tracer_config.max_steps = kRenderTraceMaximumAttempts;
    tracer_config.wormhole_topology = config_.wormhole_topology;

    // Large steps far from the hole, small near the horizon.
    tracer_config.integrator.initial_step = trace_domain.cpu_initial_step;
    tracer_config.integrator.max_step = trace_domain.max_step;
    tracer_config.integrator.min_step = trace_domain.cpu_min_step;
    tracer_config.integrator.abs_tolerance = 5e-6f;
    tracer_config.integrator.rel_tolerance = 5e-6f;

    // Disk inner edge from the ISCO. The thin disk is a black-hole construct;
    // horizonless spacetimes render lensing and background only.
    const auto disk_support = core::DiskSupportFor(config_.metric_id);
    const bool disk_capable = disk_support == core::DiskSupport::PageThorne;
    SIRIUS_ASSERT(!config_.enable_disk || disk_support == core::DiskSupport::PageThorne);
    SIRIUS_ASSERT(config_.enable_disk ||
                  (!config_.enable_volumetric_disk && !config_.enable_turbulence));
    auto isco_radius = AccretionDiskD::ComputeIsco(config_.black_hole_spin);
    tracer_config.enable_disk = config_.enable_disk && disk_capable;
    if (tracer_config.enable_disk) {
        tracer_config.disk_inner = isco_radius * config_.black_hole_mass;
        tracer_config.disk_outer = 20.0 * config_.black_hole_mass;
        tracer_config.enable_polarisation = config_.enable_polarisation;
    }

    // Doppler beaming toggle (P4); true keeps the full physics and the
    // pinned render byte-for-byte.
    tracer_config.doppler_beaming = config_.doppler_beaming;
    switch (config_.temperature_model) {
        case DiskTemperatureModel::NovikovThorne:
            tracer_config.disk_temperature_model = backend::DiskTemperatureModel::NovikovThorne;
            break;
        case DiskTemperatureModel::ShakuraSunyaev:
            tracer_config.disk_temperature_model = backend::DiskTemperatureModel::ShakuraSunyaev;
            break;
        default:
            return base::Fail(base::ErrorDomain::kConfiguration, "initialise render session",
                              "invalid disk temperature model reached initialisation");
    }

    // Ray bundles (P2/P3). Enabled for the beam footprint the filtered star
    // field consumes; the pupil (point-source) mode gives the celestial-sphere
    // footprint. Default off leaves the point-sampled path and the pinned
    // render untouched.
    pixel_angular_size_ = (config_.camera_fov * math::kPi / 180.0) / std::max(1, config_.height);
    tracer_config.enable_ray_bundles = config_.ray_bundles;
    if (tracer_config.enable_ray_bundles) {
        tracer_config.bundle_point_source = true;
        tracer_config.bundle_angular_size = static_cast<float>(pixel_angular_size_);
    }

    // Volumetric disk configuration.
    tracer_config.enable_volumetric = config_.enable_volumetric_disk;
    if (tracer_config.enable_volumetric) {
        tracer_config.volumetric_scale_height_ratio = config_.volumetric_h_over_r;
        tracer_config.volumetric_flare_power = config_.volumetric_h_power;
        tracer_config.volumetric_tau_midplane = config_.volumetric_tau_midplane;
        tracer_config.volumetric_samples = config_.volumetric_samples;
        tracer_config.enable_turbulence = config_.enable_turbulence;
        tracer_config.disk_temperature_scale_kelvin = config_.disk_temperature_scale;
        tracer_config.color_mode = config_.color_mode;
    }
    if (metric_) {
        tracer_ = std::make_unique<GeodesicTracer>(metric_.get(), tracer_config);
        if (tracer_config.enable_disk) {
            std::cout << "  Disk:       r_in=" << tracer_config.disk_inner
                      << "M, r_out=" << tracer_config.disk_outer << "M" << std::endl;
        } else if (disk_capable) {
            std::cout << "  Disk:       disabled" << std::endl;
        }
        std::cout << "[Session] Physics engine initialized (geodesic integration enabled)"
                  << std::endl;
    }

    // Relativistic jets.
    // Colour mode.
    const char* mode_name = "TrueColor";
    switch (config_.color_mode) {
        case core::color_modes::Mode::TrueColor:
            mode_name = "TrueColor (Physical)";
            break;
        case core::color_modes::Mode::TemperatureMap:
            mode_name = "TemperatureMap (False Color)";
            break;
        case core::color_modes::Mode::RedshiftMap:
            mode_name = "RedshiftMap (g-factor)";
            break;
        case core::color_modes::Mode::Polarisation:
            mode_name = "Polarisation";
            break;
    }
    std::cout << "[Session] Color mode: " << mode_name << std::endl;

    // Multi-threaded rendering setup.
    if (config_.enable_parallel_rendering) {
        num_threads_ = config_.thread_count;
        if (num_threads_ <= 0) {
            // Auto-detect: hardware concurrency, leaving one core for the system.
            num_threads_ = static_cast<int>(std::thread::hardware_concurrency());
            if (num_threads_ > 1) num_threads_--;
            if (num_threads_ < 1) num_threads_ = 1;
        }

        // Per-thread tracers share the metric but keep independent state.
        thread_tracers_.clear();
        thread_tracers_.reserve(num_threads_);
        if (metric_) {
            for (int i = 0; i < num_threads_; ++i) {
                thread_tracers_.push_back(
                    std::make_unique<GeodesicTracer>(metric_.get(), tracer_config));
            }
        }

        std::cout << "[Session] Multi-threaded rendering enabled: " << num_threads_ << " threads"
                  << std::endl;
    } else {
        num_threads_ = 1;
        std::cout << "[Session] Single-threaded rendering" << std::endl;
    }

    // The background texture is a physics input, not decorative packaging:
    // substituting grey changes every escaped ray. Missing resources
    // therefore fail initialisation instead of quietly changing the scene.
    const auto starfield = base::ResolveResource("assets/Starfield.png");
    if (!starfield || !LoadStarfieldTexture(starfield->string())) {
        return base::Fail(
            base::ErrorDomain::kIo, "load render resource",
#if SIRIUS_RELEASE_RESOURCE_LOCKED
            "assets/Starfield.png is missing or unreadable in the packaged Sirius volume");
#else
            "assets/Starfield.png is missing or unreadable (set SIRIUS_RESOURCE_DIR to an "
            "installed share/sirius directory)");
#endif
    }

    // Filtered point-source star field (P3): build the deterministic catalogue
    // once. The beam footprint (ray bundles on) or a pinhole sigma (bundles
    // off) filters it per escaping ray.
    if (config_.point_starfield) {
        const core::StarfieldConfig scfg =
            core::ExpandPointStarfieldConfig(config_.point_starfield_config);
        star_generator_ = std::make_unique<core::StarfieldGenerator>(scfg);
        star_index_ =
            std::make_unique<core::StarfieldSpatialIndex>(star_generator_->GenerateCatalogue());
        std::cout << "[Session] Point-source star field: " << star_index_->Size() << " stars, "
                  << (star_index_->MemoryBytes() / 1024) << " KiB index, beams "
                  << (config_.ray_bundles ? "on" : "off") << std::endl;
    }
    std::cout << "[Session] Scene evidence: "
              << SessionSceneEvidenceJson(config_, star_index_ ? star_index_->Size() : 0)
              << std::endl;

    // GPU acceleration removed: the legacy OptiX backend init and starfield
    // upload lived here. OptiX is retired; the Vulkan compute path arrives
    // through sirius::backend::device later. The CPU path renders directly.

    if (progress_.GetCancellationToken().IsCancelled()) {
        fsm_.Process(SessionEvent::Cancel);
        return {};
    }

    // Transition to Scheduling.
    fsm_.Process(SessionEvent::Ready);
    return {};
}

// =============================================================================
// Tile scheduling.
// =============================================================================
void RenderSession::ScheduleNextTile() {
    if (progress_.GetCancellationToken().IsCancelled()) {
        fsm_.Process(SessionEvent::Cancel);
        return;
    }

    // Vulkan render path: dispatches the Slang trace kernel through the device
    // seam (replacing the retired OptiX single-launch that once lived here). It
    // fills the display buffer directly and fires AllTilesComplete so WriteOutput
    // applies the host display pipeline, exactly as the CPU path does.
    if (config_.backend == RenderBackend::Vulkan) {
        RenderVulkanPath();
        return;
    }

    // The CPU path needs a CPU-representable metric; refuse rather than render a
    // different spacetime.
    if (!metric_) {
        error_message_ = "Metric '" +
                         std::string(core::MetricInfoFor(config_.metric_id).canonical_name) +
                         "' is not supported on the CPU backend (GPU/Vulkan required)";
        std::cerr << "[Session] " << error_message_ << std::endl;
        fsm_.Process(SessionEvent::Error);
        return;
    }

    // Parallel rendering when enabled and multiple threads are available.
    if (config_.enable_parallel_rendering && num_threads_ > 1) {
        RenderTilesParallel();
        return;
    }

    // Single-threaded: sequential tile processing.
    Tile* tile = tiles_.GetNextTile();

    if (tile == nullptr) {
        fsm_.Process(SessionEvent::AllTilesComplete);
        return;
    }

    fsm_.Process(SessionEvent::TileAvailable);

    RenderTile(tile);
}

// =============================================================================
// Pixel shading helpers (unified for single- and multi-threaded paths).
// =============================================================================
RenderSession::PixelResult RenderSession::ShadeDiskHit(const TraceResult& result) const {
    PixelResult px;

    // Volumetric disk: use the pre-integrated emission from ray marching.
    if (result.volumetric_hit) {
        // Samples were coloured at their own temperature and g-factor before
        // invariant transfer. Reconstructing one effective blackbody here would
        // apply the spectral mapping twice.
        px.r = result.volumetric_emission[0];
        px.g = result.volumetric_emission[1];
        px.b = result.volumetric_emission[2];
        return px;
    }

    // Thin disk: accumulate emission from all disk crossings. Relativistic
    // beaming is applied exactly once: emitted T^4 becomes observed g^4 T^4.
    float total_r = 0.0f, total_g = 0.0f, total_b = 0.0f;
    core::StokesVector total_stokes;
    const bool polarisation_mode = config_.color_mode == core::color_modes::Mode::Polarisation;

    for (int crossing_idx = 0; crossing_idx < result.num_disk_crossings; crossing_idx++) {
        const auto& crossing = result.disk_crossings[crossing_idx];
        if (!crossing.valid) continue;

        float T_emit = crossing.temperature;
        float g = crossing.redshift;

        // The stationary axisymmetric disk has one covariant transfer sample at
        // a crossing. Temporal requests fail at validation until an evolving
        // emissivity (rather than an azimuth-shifted steady flow) is represented.
        const std::array<float, 1> temporal_redshifts{g};

        const float emitted_intensity = std::pow(T_emit, 4.0f);
        const float observed_intensity =
            core::color_modes::ObservedBolometricIntensity(emitted_intensity, g);

        if (polarisation_mode) {
            SIRIUS_ASSERT(crossing.polarisation_valid);
            if (!crossing.polarisation_valid) continue;

            const float chi = crossing.polarisation_evpa;
            const float degree = crossing.polarisation_degree;
            const float atmosphere_intensity =
                observed_intensity * crossing.polarisation_intensity_scale;
            core::StokesVector crossing_stokes{
                atmosphere_intensity, atmosphere_intensity * degree * std::cos(2.0f * chi),
                atmosphere_intensity * degree * std::sin(2.0f * chi), 0.0f};
            total_stokes += crossing_stokes;
            continue;
        }

        const core::spectral::Rgb disk_color = core::color_modes::AverageTemporalColorMode(
            config_.color_mode, T_emit, temporal_redshifts, emitted_intensity,
            config_.disk_temperature_scale);

        total_r += disk_color.r;
        total_g += disk_color.g;
        total_b += disk_color.b;
    }

    // Lensing changes the ray-to-solid-angle map, not radiance along one ray.
    constexpr float output_scale = 1.0f;

    if (polarisation_mode) {
        total_stokes *= output_scale;
        total_stokes.Normalise();
        const core::spectral::Rgb visualised =
            core::color_modes::polarisation_vis::StokesToRgbHsv(total_stokes);
        px.r = visualised.r;
        px.g = visualised.g;
        px.b = visualised.b;
    } else {
        px.r = total_r * output_scale;
        px.g = total_g * output_scale;
        px.b = total_b * output_scale;
    }

    return px;
}

RenderSession::PixelResult RenderSession::ShadeEscaped(const TraceResult& result) const {
    PixelResult px;
    if (config_.point_starfield && star_index_ && star_index_->Size() > 0) {
        SampleStarfieldPoints(result, px.r, px.g, px.b);
    } else {
        SampleStarfield(result.final_direction, px.r, px.g, px.b);
    }

    return px;
}

RenderSession::PixelResult RenderSession::ShadePixel(int px_coord, int py_coord,
                                                     GeodesicTracer* tracer) const {
    PixelResult result;
    float r_acc = 0.0f, g_acc = 0.0f, b_acc = 0.0f;

    const int samples_taken =
        ForEachCameraSample(config_.samples_per_pixel, [&](const CameraSample& sample) {
            CameraRay camera_ray = camera_->GenerateRayForObserver(
                px_coord, py_coord, sample.image_u, sample.image_v, sample.pupil_u, sample.pupil_v);
            SIRIUS_ASSERT(core::IsRepresentedCameraRay(camera_ray));
            if (!camera_ray.active) return;
            TraceResult trace_result = tracer->Trace(camera_ray);

            float sr = 0.0f, sg = 0.0f, sb = 0.0f;

            switch (trace_result.outcome) {
                case TraceResult::Outcome::Horizon:
                case TraceResult::Outcome::Throat:
                    break;

                case TraceResult::Outcome::DiskHit: {
                    PixelResult disk = ShadeDiskHit(trace_result);
                    sr = disk.r;
                    sg = disk.g;
                    sb = disk.b;
                    break;
                }

                case TraceResult::Outcome::Escaped: {
                    PixelResult esc = ShadeEscaped(trace_result);
                    sr = esc.r;
                    sg = esc.g;
                    sb = esc.b;
                    break;
                }

                case TraceResult::Outcome::MaxSteps:
                    sr = sg = sb = 0.0f;
                    break;
                default:
                    SIRIUS_ASSERT(false);
                    sr = 1.0f;
                    sg = 0.0f;
                    sb = 1.0f;
                    break;
            }

            // Volumetric transfer composes with the terminal surface/background;
            // it is not a terminal ray outcome. Apply I = I_bg exp(-tau) + I_vol
            // after shading the actual fate of the central ray.
            if (trace_result.volumetric_hit) {
                PixelResult volume = ShadeDiskHit(trace_result);
                const float transmission = std::exp(-std::max(trace_result.optical_depth, 0.0f));
                sr = sr * transmission + volume.r;
                sg = sg * transmission + volume.g;
                sb = sb * transmission + volume.b;
            }

            r_acc += sr;
            g_acc += sg;
            b_acc += sb;
        });

    float inv_samples = 1.0f / static_cast<float>(samples_taken);
    result.r = r_acc * inv_samples;
    result.g = g_acc * inv_samples;
    result.b = b_acc * inv_samples;

    return result;
}

// =============================================================================
// Tile rendering.
// =============================================================================
void RenderSession::RenderTile(Tile* tile) {
    if (!tile) return;

    std::vector<float> tileBuffer(tile->width * tile->height * 4, 0.0f);

    for (int ty = 0; ty < tile->height; ++ty) {
        if (progress_.GetCancellationToken().IsCancelled()) {
            fsm_.Process(SessionEvent::Cancel);
            return;
        }
        for (int tx = 0; tx < tile->width; ++tx) {
            int px = tile->x + tx;
            int py = tile->y + ty;

            PixelResult pixel = ShadePixel(px, py, tracer_.get());

            int idx = (ty * tile->width + tx) * 4;
            tileBuffer[idx + 0] = pixel.r;
            tileBuffer[idx + 1] = pixel.g;
            tileBuffer[idx + 2] = pixel.b;
            tileBuffer[idx + 3] = 1.0f;
        }
    }

    display_.UpdateTile(tile->x, tile->y, tile->width, tile->height, tileBuffer.data());

    tiles_.CompleteTile(tile->id);
    // The ProgressTracker callback is the single progress surface (the CLI
    // renders it); a second raw carriage-return writer here would fight it.
    progress_.CompleteTile(tile->PixelCount());

    fsm_.Process(SessionEvent::TileComplete);
}

// =============================================================================
// Vulkan render path.
// =============================================================================
void RenderSession::RenderVulkanPath() {
#ifdef SIRIUS_HAS_VULKAN_BACKEND
    std::cout << "[Session] Dispatching Vulkan render path..." << std::endl;

    auto stats = RenderVulkanToDisplay(
        config_, display_,
        [this](int /*done*/, int total) {
            // The scheduler's tile count is the CPU spiral grid; the Vulkan path
            // governs its own tile count, so retotal the tracker to it on the
            // first report and then count governed tiles through CompleteTile.
            if (progress_.GetTilesTotal() != total) {
                progress_.SetTotals(total, 1);
            }
            progress_.CompleteTile(1);
        },
        [this] { return progress_.GetCancellationToken().IsCancelled(); });

    if (!stats) {
        if (progress_.GetCancellationToken().IsCancelled()) {
            fsm_.Process(SessionEvent::Cancel);
            return;
        }
        error_message_ = stats.error().Description();
        std::cerr << "[Session] Vulkan render declined: " << error_message_ << std::endl;
        fsm_.Process(SessionEvent::Error);
        return;
    }

    std::cout << "[Session] Vulkan render complete: " << stats->metric_name << " on "
              << stats->device_name << ", " << stats->tiles_rendered << " tile(s) of "
              << stats->tile_plan.tile_edge << "px in " << stats->band_dispatches
              << " governed dispatch(es), " << stats->seconds << "s" << std::endl;
    fsm_.Process(SessionEvent::AllTilesComplete);
#else
    error_message_ = "Vulkan backend not compiled in (build without Vulkan development files)";
    std::cerr << "[Session] " << error_message_ << std::endl;
    fsm_.Process(SessionEvent::Error);
#endif
}

// =============================================================================
// Starfield texture loading.
// =============================================================================
bool RenderSession::LoadStarfieldTexture(const std::string& path) {
    int channels;
    auto deleter = [](unsigned char* p) {
        if (p) stbi_image_free(p);
    };
    std::unique_ptr<unsigned char[], decltype(deleter)> data(
        stbi_load(path.c_str(), &starfield_width_, &starfield_height_, &channels, 4), deleter);

    if (!data) {
        std::cerr << "[Session] Failed to load starfield texture: " << path << std::endl;
        std::cerr << "[Session] STB error: " << stbi_failure_reason() << std::endl;
        return false;
    }

    const std::size_t size = static_cast<std::size_t>(starfield_width_) *
                             static_cast<std::size_t>(starfield_height_) * 4;
    starfield_data_.assign(data.get(), data.get() + size);

    starfield_loaded_ = true;
    std::cout << "[Session] Loaded starfield texture: " << starfield_width_ << "x"
              << starfield_height_ << " RGBA (" << (size / 1024 / 1024) << " MB)" << std::endl;

    return true;
}

// =============================================================================
// Starfield background sampling (equirectangular projection).
// =============================================================================
void RenderSession::SampleStarfield(const Vec4& direction, float& r, float& g, float& b) const {
    if (!starfield_loaded_ || starfield_data_.empty()) {
        SIRIUS_ASSERT(starfield_loaded_ && !starfield_data_.empty());
        r = g = b = 0.0f;
        return;
    }

    // Direction is Cartesian (x, y, z) from the geodesic tracer.
    double dx = direction(1);
    double dy = direction(2);
    double dz = direction(3);

    double len = std::sqrt(dx * dx + dy * dy + dz * dz);
    if (len < 1e-10) {
        SIRIUS_ASSERT(len >= 1e-10);
        r = g = b = 0.0f;
        return;
    }

    dx /= len;
    dy /= len;
    dz /= len;

    // Spherical coordinates: theta 0 at +Z to pi at -Z; phi 0 at +X toward +Y.
    double theta = std::acos(std::clamp(dz, -1.0, 1.0));
    double phi = std::atan2(dy, dx);
    if (phi < 0) phi += math::kTwoPi;

    double u = phi / math::kTwoPi;
    double v = theta / math::kPi;

    double px = u * (starfield_width_ - 1);
    double py = v * (starfield_height_ - 1);

    int x0 = static_cast<int>(std::floor(px));
    int y0 = static_cast<int>(std::floor(py));
    int x1 = std::min(x0 + 1, starfield_width_ - 1);
    int y1 = std::min(y0 + 1, starfield_height_ - 1);

    double fx = px - x0;
    double fy = py - y0;

    auto sample = [this](int x, int y) -> std::array<float, 3> {
        int idx = (y * starfield_width_ + x) * 4;  // 4 bytes per pixel (RGBA).
        return {starfield_data_[idx + 0] / 255.0f, starfield_data_[idx + 1] / 255.0f,
                starfield_data_[idx + 2] / 255.0f};
    };

    auto c00 = sample(x0, y0);
    auto c10 = sample(x1, y0);
    auto c01 = sample(x0, y1);
    auto c11 = sample(x1, y1);

    float w00 = static_cast<float>((1.0 - fx) * (1.0 - fy));
    float w10 = static_cast<float>(fx * (1.0 - fy));
    float w01 = static_cast<float>((1.0 - fx) * fy);
    float w11 = static_cast<float>(fx * fy);

    r = c00[0] * w00 + c10[0] * w10 + c01[0] * w01 + c11[0] * w11;
    g = c00[1] * w00 + c10[1] * w10 + c01[1] * w01 + c11[1] * w11;
    b = c00[2] * w00 + c10[2] * w10 + c01[2] * w01 + c11[2] * w11;
}

// =============================================================================
// Filtered point-source star field sampling (P3).
// =============================================================================
void RenderSession::SampleStarfieldPoints(const TraceResult& result, float& r, float& g,
                                          float& b) const {
    // Beam footprint on the sky. With ray bundles the tracer supplies the full
    // lensed ellipse; a pinhole (bundles off) samples at a fraction of the pixel
    // angular size, so a star pops in and out as the camera rotates (the flicker
    // the beam filter removes). Both ellipse axes are floored at the pixel size so
    // a pixel always integrates at least its own solid angle.
    constexpr float kPinholeFraction = 0.3f;  // Pinhole sigma as a fraction of a pixel.
    float pixel = static_cast<float>(pixel_angular_size_);
    float sigma_major;
    float sigma_minor;
    float orientation = 0.0f;
    if (config_.ray_bundles && result.beam.valid) {
        sigma_major = std::max(result.beam.footprint_major, pixel);
        sigma_minor = std::max(result.beam.footprint_minor, pixel);
        orientation = result.beam.orientation;
    } else {
        sigma_major = kPinholeFraction * pixel;
        sigma_minor = sigma_major;
    }

    const auto& d = result.final_direction;
    star_generator_->AccumulateThroughBeam(static_cast<float>(d(1)), static_cast<float>(d(2)),
                                           static_cast<float>(d(3)), sigma_major, sigma_minor,
                                           orientation, *star_index_, r, g, b);
}

// =============================================================================
// Output writing.
// =============================================================================
base::Expected<void> RenderSession::WriteOutput() {
    std::cout << "\n[Session] Writing output..." << std::endl;

    const int width = config_.width;
    const int height = config_.height;

    if (const auto bad = display_.FirstNonFiniteIndex(); bad.has_value()) {
        return base::Fail(
            base::ErrorDomain::kPhysics, "write render output",
            "linear radiance contains a non-finite sample at channel " + std::to_string(*bad));
    }

    // Format follows the output extension: .exr | .png | .ppm (default).
    std::string ext;
    {
        std::size_t dot = config_.output_path.rfind('.');
        if (dot != std::string::npos) ext = config_.output_path.substr(dot + 1);
        std::transform(ext.begin(), ext.end(), ext.begin(),
                       [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    }

    if (config_.write_output && ext == "exr") {
        // EXR is the HDR interchange path: the linear radiance buffer is
        // written untouched. Tonemapping, grading, and transfer encoding are
        // display concerns and deliberately do not apply here.
        ImageBufferRGBA hdr(width, height);
        hdr.pixels = display_.SnapshotFloatData();

        EXRMetadata meta;
        meta.metric_type = core::MetricInfoFor(config_.metric_id).canonical_name;
        meta.black_hole_mass = config_.black_hole_mass;
        meta.black_hole_spin = config_.black_hole_spin;
        meta.observer_distance = config_.observer_distance;
        meta.observer_inclination = config_.observer_inclination * 180.0 / math::kPi;
        meta.samples_per_pixel = config_.samples_per_pixel;

        if (!EXRWriter::WriteExr(config_.output_path, hdr, meta)) {
            return base::Fail(base::ErrorDomain::kIo, "write EXR", config_.output_path);
        }
    } else {
        // Display pipeline: tonemap and grade to display-linear values. The
        // transfer encode belongs to the writer, so it happens exactly once
        // per format (EXRWriter's PPM boundary and PNGWriter each apply sRGB).
        core::PostProcessConfig postprocess_config;
        postprocess_config.tonemapper = config_.tonemapper;
        postprocess_config.exposure = config_.exposure;
        postprocess_config.gamma = 1.0f;  // Writers encode.

        postprocess_config.enable_bloom = config_.enable_bloom;
        if (config_.enable_bloom) {
            postprocess_config.bloom_intensity = config_.bloom_intensity;
            postprocess_config.bloom_threshold = config_.bloom_threshold;
            postprocess_config.bloom_radius = kBloomRadius;
        }

        postprocess_config.saturation = config_.saturation;
        postprocess_config.contrast = config_.contrast;
        postprocess_config.lift = kShadowLift;
        postprocess_config.gain = 1.0f;

        display_.MutateFloatData([&](std::vector<float>& pixels) {
            core::PostProcessor::Process(pixels, width, height, postprocess_config);

            // Film is a display-referred finishing pipeline. It is intentionally
            // absent from the EXR branch above, which must retain untouched
            // linear HDR radiance.
            if (config_.enable_film_finish) {
                FilmPipeline film(config_.film_config);
                film.Apply(pixels.data(), width, height, 0);
            }
        });

        if (const auto bad = display_.FirstNonFiniteIndex(); bad.has_value()) {
            return base::Fail(base::ErrorDomain::kInternal, "apply display pipeline",
                              "produced a non-finite sample at channel " + std::to_string(*bad));
        }

        const std::vector<float> output = display_.SnapshotFloatData();
        const float* data = output.data();

        if (!config_.write_output) {
            // In-memory progressive previews consume the display-linear
            // buffer directly and must not overwrite the configured output.
        } else if (ext == "png") {
            if (!PNGWriter::WriteRgba(config_.output_path, width, height, data)) {
                return base::Fail(base::ErrorDomain::kIo, "write PNG", config_.output_path);
            }
        } else {
            if (ext != "ppm") {
                std::cout << "[Session] Unrecognised output extension '" << ext
                          << "'; writing PPM data" << std::endl;
            }
            if (!EXRWriter::WritePpmRgba(config_.output_path, width, height, data)) {
                return base::Fail(base::ErrorDomain::kIo, "write PPM", config_.output_path);
            }
        }
    }

    if (config_.write_output) {
        std::cout << "[Session] Wrote: " << config_.output_path << std::endl;
    }

    fsm_.Process(SessionEvent::OutputWritten);
    return {};
}

// =============================================================================
// Session end.
// =============================================================================
void RenderSession::OnSessionEnd(SessionState state) {
    std::cout << "\n[Session] Finished with state: " << StateName(state) << std::endl;

    std::string message;
    switch (state) {
        case SessionState::Complete:
            message = "Render completed successfully";
            break;
        case SessionState::Failed:
            message = error_message_;
            break;
        case SessionState::Cancelled:
            message = "Render cancelled by user";
            break;
        default:
            SIRIUS_ASSERT(false);
            message = "Unknown end state";
    }

    CompletionCallback callback;
    {
        std::lock_guard<std::mutex> lock(callback_mutex_);
        callback = completion_callback_;
    }
    if (callback) {
        try {
            callback(state, message);
        } catch (...) {
            std::cerr << "[Session] completion callback threw; terminal state retained"
                      << std::endl;
        }
    }

    stop_workers_ = true;
    for (auto& thread : worker_threads_) {
        if (thread.joinable()) {
            thread.join();
        }
    }
    worker_threads_.clear();
}

// =============================================================================
// Multi-threaded tile rendering.
// =============================================================================
void RenderSession::RenderTilesParallel() {
    stop_workers_ = false;
    active_workers_ = 0;

    worker_threads_.clear();
    worker_threads_.reserve(num_threads_);

    std::cout << "[Session] Starting " << num_threads_ << " worker threads..." << std::endl;

    for (int i = 0; i < num_threads_; ++i) {
        worker_threads_.emplace_back(&RenderSession::WorkerThread, this, i);
    }

    for (auto& thread : worker_threads_) {
        if (thread.joinable()) {
            thread.join();
        }
    }
    worker_threads_.clear();

    std::cout << "\n[Session] All worker threads completed" << std::endl;

    if (progress_.GetCancellationToken().IsCancelled()) {
        fsm_.Process(SessionEvent::Cancel);
    } else {
        fsm_.Process(SessionEvent::AllTilesComplete);
    }
}

void RenderSession::WorkerThread(int thread_id) {
    active_workers_++;

    while (!stop_workers_) {
        Tile* tile = nullptr;
        {
            std::lock_guard<std::mutex> lock(tile_mutex_);
            tile = tiles_.GetNextTile();
        }

        if (tile == nullptr) {
            break;
        }

        if (progress_.GetCancellationToken().IsCancelled()) {
            stop_workers_ = true;
            break;
        }

        if (!RenderTileThreaded(tile, thread_id)) {
            stop_workers_ = true;
            break;
        }

        {
            std::lock_guard<std::mutex> lock(tile_mutex_);
            tiles_.CompleteTile(tile->id);
            // Single progress surface: the ProgressTracker callback.
            progress_.CompleteTile(tile->PixelCount());
        }
    }

    active_workers_--;
}

bool RenderSession::RenderTileThreaded(Tile* tile, int thread_id) {
    if (!tile) return false;

    GeodesicTracer* tracer = nullptr;
    if (thread_id >= 0 && thread_id < static_cast<int>(thread_tracers_.size())) {
        tracer = thread_tracers_[thread_id].get();
    } else {
        tracer = tracer_.get();
    }

    std::vector<float> tileBuffer(tile->width * tile->height * 4, 0.0f);

    for (int ty = 0; ty < tile->height; ++ty) {
        if (ty % 8 == 0 && (stop_workers_ || progress_.GetCancellationToken().IsCancelled())) {
            return false;
        }

        for (int tx = 0; tx < tile->width; ++tx) {
            int px = tile->x + tx;
            int py = tile->y + ty;

            PixelResult pixel = ShadePixel(px, py, tracer);

            int idx = (ty * tile->width + tx) * 4;
            tileBuffer[idx + 0] = pixel.r;
            tileBuffer[idx + 1] = pixel.g;
            tileBuffer[idx + 2] = pixel.b;
            tileBuffer[idx + 3] = 1.0f;
        }
    }

    {
        std::lock_guard<std::mutex> lock(display_mutex_);
        display_.UpdateTile(tile->x, tile->y, tile->width, tile->height, tileBuffer.data());
    }
    return true;
}

}  // namespace sirius::render
