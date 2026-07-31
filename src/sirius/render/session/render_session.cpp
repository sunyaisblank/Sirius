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

#ifdef SIRIUS_HAS_VULKAN_BACKEND
#include "sirius/backend/device.h"  // Device enumeration for backend auto-resolution.
#include "sirius/render/vulkan_renderer.h"
#endif

#include "sirius/core/constants.h"
#include "sirius/core/metrics/morris_thorne_family.h"
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
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>

namespace sirius::render {

using backend::GeodesicTracer;
using backend::TracerConfig;
using backend::TraceResult;
using core::AccretionDiskD;
using core::CameraConfig;
using core::CameraRay;
using core::JetConfig;
using core::KerrSchildFamily;
using core::KerrSchildParams;
using core::MetricId;
using core::RelativisticJet;
using core::Vec4;
using core::WarpDriveFamily;
using core::WarpDriveParams;

namespace math = core::constants::math;

namespace {
// Shading and grading constants (kept from the legacy reference path).
constexpr float kDiskIntensityBoost = 5.0f;
constexpr float kPhotonRingBoostDisk = 2.0f;
constexpr float kPhotonRingBoostVolumetric = 1.5f;
constexpr float kPhotonRingBoostEscaped = 1.2f;
constexpr float kSpiralingBrightness = 0.02f;
constexpr float kMaxStepsBrightness = 0.01f;
constexpr float kDopplerClampMin = 0.1f;
constexpr float kDopplerClampMax = 10.0f;
constexpr float kGFactorClampMin = 0.1f;
constexpr float kGFactorClampMax = 5.0f;
constexpr float kJetNormalisation = 0.1f;
constexpr int kBloomRadius = 12;
constexpr float kShadowLift = 0.02f;
}  // namespace

RenderSession::~RenderSession() {
    (void)Cancel();
    WaitForCompletion();
}

bool RenderSession::Configure(const SessionConfig& config) {
    std::lock_guard<std::mutex> lock(lifecycle_mutex_);
    if (fsm_.GetState() != SessionState::Idle || render_thread_.joinable() || join_in_progress_) {
        return false;
    }
    config_ = config;
    return true;
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

std::optional<std::string> SessionConfigIssue(const SessionConfig& config) {
    const auto finite = [](double value) { return std::isfinite(value); };
    const auto in_range = [&finite](double value, double minimum, double maximum) {
        return finite(value) && value >= minimum && value <= maximum;
    };

    if (config.width < 1 || config.width > 8192 || config.height < 1 || config.height > 8192) {
        return "width and height must each be between 1 and 8192";
    }
    if (config.tileSize < 1 || config.tileSize > 256 ||
        (config.tileSize & (config.tileSize - 1)) != 0) {
        return "tile size must be a power of two between 1 and 256";
    }
    if (config.samplesPerPixel < 1 || config.samplesPerPixel > 4096) {
        return "samples per pixel must be between 1 and 4096";
    }
    if (config.threadCount < 0 || config.threadCount > 1024) {
        return "thread count must be between 0 and 1024";
    }
    if (config.pointStarfield) {
        const core::StarfieldConfig& stars = config.starfieldConfig;
        if (!stars.enabled || stars.star_count < 1 || stars.star_count > 10000000u ||
            !finite(stars.min_distance_pc) || stars.min_distance_pc < 0.1f ||
            !finite(stars.max_distance_pc) ||
            stars.max_distance_pc < stars.min_distance_pc + 1.0f ||
            !finite(stars.magnitude_limit) || stars.magnitude_limit < 0.0f ||
            stars.magnitude_limit > 20.0f || !finite(stars.aperture_mm) ||
            stars.aperture_mm < 0.0f || stars.aperture_mm > 1000.0f ||
            !finite(stars.focus_distance_pc) || stars.focus_distance_pc < stars.min_distance_pc ||
            stars.focus_distance_pc > stars.max_distance_pc || !finite(stars.brightness_scale) ||
            stars.brightness_scale < 0.0f) {
            return "point-starfield parameters are outside the represented domain";
        }
    }
    if (config.writeOutput) {
        if (config.outputPath.empty() || config.outputPath.find('\0') != std::string::npos) {
            return "output path must not be empty or contain a null byte";
        }
        std::string extension;
        if (const auto dot = config.outputPath.rfind('.'); dot != std::string::npos) {
            extension = config.outputPath.substr(dot);
            std::transform(extension.begin(), extension.end(), extension.begin(),
                           [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
        }
        if (extension != ".ppm" && extension != ".png" && extension != ".exr") {
            return "output path extension must be .ppm, .png, or .exr";
        }
    }

    switch (config.backend) {
        case RenderBackend::Cpu:
        case RenderBackend::Vulkan:
            break;
        default:
            return "invalid render backend";
    }
    if (const auto issue =
            core::MetricParameterIssue(config.metricId, config.blackHoleSpin,
                                       config.blackHoleCharge, config.cosmologicalConstant);
        issue.has_value()) {
        return std::string(*issue);
    }
    const bool massless_metric =
        config.metricId == MetricId::Minkowski || config.metricId == MetricId::DeSitter;
    if (!finite(config.blackHoleMass) || config.blackHoleMass < 0.0 ||
        config.blackHoleMass > 100.0 || !finite(config.blackHoleSpin) ||
        config.blackHoleSpin < 0.0 || config.blackHoleSpin > 0.998 ||
        !finite(config.blackHoleCharge) || config.blackHoleCharge < 0.0 ||
        config.blackHoleCharge > 0.999 || !finite(config.cosmologicalConstant) ||
        std::abs(config.cosmologicalConstant) > 0.1) {
        return "metric mass, spin, charge, or lambda is outside the represented domain";
    }
    if (config.blackHoleSpin * config.blackHoleSpin +
            config.blackHoleCharge * config.blackHoleCharge >=
        0.999) {
        return "combined spin squared plus charge squared must be below 0.999";
    }
    if ((massless_metric && config.blackHoleMass != 0.0) ||
        (!massless_metric && config.blackHoleMass < 0.1)) {
        return massless_metric ? "massless metrics require mass to be zero"
                               : "mass must be between 0.1 and 100 for this metric";
    }

    const double distance_scale = massless_metric ? 1.0 : config.blackHoleMass;
    if (!finite(config.observerDistance) || config.observerDistance < 5.0 * distance_scale ||
        config.observerDistance > 1000.0 * distance_scale || !finite(config.observerInclination) ||
        config.observerInclination <= 0.1 * math::kPi / 180.0 ||
        config.observerInclination >= 179.9 * math::kPi / 180.0 ||
        !finite(config.observerAzimuth) || std::abs(config.observerAzimuth) > 2.0 * math::kPi ||
        !finite(config.cameraFOV) || config.cameraFOV < 1.0f || config.cameraFOV > 170.0f) {
        return "observer distance, angles, or field of view is outside the represented domain";
    }
    const double beta_squared = config.cameraBetaForward * config.cameraBetaForward +
                                config.cameraBetaUp * config.cameraBetaUp +
                                config.cameraBetaRight * config.cameraBetaRight;
    if (!finite(beta_squared) || beta_squared >= 1.0) {
        return "camera beta magnitude must be finite and below one";
    }
    switch (config.lensType) {
        case core::LensType::Pinhole:
        case core::LensType::ThinLens:
        case core::LensType::Fisheye:
            break;
        default:
            return "invalid lens model";
    }
    if (!in_range(config.cameraFocalLength, std::numeric_limits<float>::min(), 10000.0) ||
        !in_range(config.cameraAperture, std::numeric_limits<float>::min(), 128.0) ||
        !in_range(config.cameraFocusDistance, std::numeric_limits<float>::min(), 1.0e6)) {
        return "camera focal length, aperture, or focus distance is outside the represented domain";
    }

    switch (config.temperatureModel) {
        case DiskTemperatureModel::NovikovThorne:
        case DiskTemperatureModel::ShakuraSunyaev:
            break;
        default:
            return "invalid disk temperature model";
    }
    if (!finite(config.diskTemperatureScale) || config.diskTemperatureScale < 100.0f ||
        config.diskTemperatureScale > 1.0e8f) {
        return "disk temperature scale must be between 100 and 100000000 Kelvin";
    }
    if (config.enableDisk &&
        core::DiskSupportFor(config.metricId) != core::DiskSupport::PageThorne) {
        return "the selected metric has no represented accretion-disk emission model";
    }
    if (!config.enableDisk &&
        (config.enableVolumetricDisk || config.enableTurbulence || config.enableCorona)) {
        return "volumetric disk, turbulence, and corona require the disk";
    }
    if (!config.enableVolumetricDisk && (config.enableTurbulence || config.enableCorona)) {
        return "turbulence and corona require volumetric transfer";
    }
    if (config.volumetricSamples < 1 || config.volumetricSamples > 4096 ||
        !finite(config.volumetricHOverR) || config.volumetricHOverR <= 0.0f ||
        config.volumetricHOverR > 2.0f || !finite(config.volumetricHPower) ||
        config.volumetricHPower < -2.0f || config.volumetricHPower > 4.0f ||
        !finite(config.volumetricTauMidplane) || config.volumetricTauMidplane < 0.0f ||
        config.volumetricTauMidplane > 1.0e6f) {
        return "volumetric transfer parameters are outside the represented domain";
    }

    switch (config.colorMode) {
        case core::color_modes::Mode::TrueColor:
        case core::color_modes::Mode::TemperatureMap:
        case core::color_modes::Mode::RedshiftMap:
        case core::color_modes::Mode::Narrowband:
        case core::color_modes::Mode::Polarisation:
            break;
        default:
            return "invalid colour mode";
    }
    const bool polarisation_mode = config.colorMode == core::color_modes::Mode::Polarisation;
    if (config.enablePolarisation != polarisation_mode) {
        return "polarisation transport and colour mode must be enabled together";
    }
    if (polarisation_mode) {
        if (!config.enableDisk) {
            return "polarisation requires the thin accretion disk";
        }
        if (config.enableVolumetricDisk) {
            return "polarisation is not represented for volumetric transfer";
        }
        if (config.metricId != MetricId::Schwarzschild && config.metricId != MetricId::Kerr) {
            return "polarisation is represented only for Schwarzschild and Kerr";
        }
        if (config.enableMotionBlur) {
            return "polarisation is not represented with temporal disk motion blur";
        }
    }
    switch (config.tonemapper) {
        case core::TonemapType::None:
        case core::TonemapType::Reinhard:
        case core::TonemapType::Aces:
        case core::TonemapType::Filmic:
        case core::TonemapType::Exposure:
            break;
        default:
            return "invalid tonemapper";
    }
    if (!in_range(config.exposure, std::numeric_limits<float>::min(), 100.0) ||
        !in_range(config.bloomIntensity, 0.0, 5.0) ||
        !in_range(config.bloomThreshold, 0.0, 100.0) || !in_range(config.contrast, 0.0, 4.0) ||
        !in_range(config.saturation, 0.0, 4.0)) {
        return "display-pipeline parameters are outside the represented domain";
    }
    if (config.enableMotionBlur) {
        if (!config.enableDisk) {
            return "temporal disk motion blur requires the disk";
        }
        if (!in_range(config.shutterTime, 0.0, 1000.0) || config.motionBlurSamples < 1 ||
            config.motionBlurSamples > 4096) {
            return "motion-blur parameters are outside the represented domain";
        }
    }
    if (config.enableFilmSimulation) {
        const FilmConfig& film = config.filmConfig;
        switch (film.format) {
            case FilmFormat::IMAX70mm_15perf:
            case FilmFormat::IMAX70mm_5perf:
            case FilmFormat::Panavision70mm:
            case FilmFormat::VistaVision:
            case FilmFormat::Academy35mm:
            case FilmFormat::Anamorphic235:
            case FilmFormat::Digital:
                break;
            default:
                return "invalid film format";
        }
        switch (film.stock) {
            case FilmStock::KodakVision3_500T:
            case FilmStock::KodakVision3_250D:
            case FilmStock::KodakVision3_50D:
            case FilmStock::KodakEktachrome:
            case FilmStock::FujiEterna:
            case FilmStock::Custom:
                break;
            default:
                return "invalid film stock";
        }
        if (!in_range(film.aspect_ratio, 0.1, 10.0) || film.width < 1 || film.width > 8192 ||
            film.height < 1 || film.height > 8192 || !in_range(film.iso, 1.0, 100000.0) ||
            !in_range(film.grain_intensity, 0.0, 1.0) || !in_range(film.grain_size, 0.0, 256.0) ||
            !in_range(film.grain_uniformity, 0.0, 1.0) ||
            !in_range(film.halation_radius, 0.0, 256.0) ||
            !in_range(film.halation_strength, 0.0, 5.0) ||
            !in_range(film.halation_threshold, 0.0, 100.0) ||
            !in_range(film.halation_color_r, 0.0, 10.0) ||
            !in_range(film.halation_color_g, 0.0, 10.0) ||
            !in_range(film.halation_color_b, 0.0, 10.0) || !in_range(film.saturation, 0.0, 4.0) ||
            !in_range(film.contrast, 0.0, 4.0) || !in_range(film.exposure, -100.0, 100.0) ||
            !in_range(film.color_temperature_K, 100.0, 100000.0) ||
            !in_range(film.tint, -10.0, 10.0) || !in_range(film.toe_strength, 0.0, 10.0) ||
            !in_range(film.shoulder_strength, 0.0, 10.0) ||
            !in_range(film.midtone_point, 0.01, 0.99) ||
            !in_range(film.weave_amplitude_x, 0.0, 256.0) ||
            !in_range(film.weave_amplitude_y, 0.0, 256.0) ||
            !in_range(film.weave_frequency, 0.0, 1000.0) ||
            !in_range(film.vignette_strength, 0.0, 2.0) ||
            !in_range(film.vignette_radius, 0.0, 10.0) ||
            !in_range(film.vignette_softness, std::numeric_limits<float>::min(), 10.0) ||
            !in_range(film.bloom_intensity, 0.0, 5.0) ||
            !in_range(film.bloom_threshold, 0.0, 100.0) ||
            !in_range(film.bloom_radius, 0.0, 256.0)) {
            return "film-simulation parameters are outside the represented domain";
        }
    }
    if (!finite(config.throatRadius) || config.throatRadius <= 0.0 ||
        config.throatRadius > 1000.0 || !finite(config.warpVelocity) ||
        std::abs(config.warpVelocity) > 10.0 || !finite(config.bubbleRadius) ||
        config.bubbleRadius <= 0.0 || config.bubbleRadius > 1000.0 || !finite(config.bubbleSigma) ||
        config.bubbleSigma <= 0.0 || config.bubbleSigma > 1000.0) {
        return "wormhole or warp-drive parameters are outside the represented domain";
    }
    switch (config.wormholeTopology) {
        case WormholeTopology::OneSheetCapture:
            break;
        case WormholeTopology::TwoSheet:
            return "two-sheet wormhole continuation and a second environment are not represented";
        default:
            return "invalid wormhole topology";
    }
    if (config.enableJets &&
        (!finite(config.jetLorentzFactor) || config.jetLorentzFactor < 1.0f ||
         !finite(config.jetOpeningAngle) || config.jetOpeningAngle <= 0.0f ||
         !finite(config.jetLaunchRadius) || config.jetLaunchRadius <= 0.0f ||
         !finite(config.jetMaxExtent) || config.jetMaxExtent <= config.jetLaunchRadius ||
         !finite(config.jetCollimation) || !finite(config.jetSpectralIndex) ||
         !finite(config.jetIntensity) || config.jetIntensity < 0.0f)) {
        return "relativistic-jet parameters are outside the represented domain";
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
        case SessionState::Initialising:
            Initialise();
            break;
        case SessionState::Scheduling:
            ScheduleNextTile();
            break;
        case SessionState::Rendering:
            // Tile is being processed.
            break;
        case SessionState::Completing:
            WriteOutput();
            break;
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
void RenderSession::Initialise() {
    try {
        if (progress_.GetCancellationToken().IsCancelled()) {
            fsm_.Process(SessionEvent::Cancel);
            return;
        }
        if (const auto issue = SessionConfigIssue(config_); issue.has_value()) {
            throw std::invalid_argument("SessionConfig: " + *issue);
        }
        std::cout << "[Session] Initialising render..." << std::endl;
        std::cout << "  Resolution: " << config_.width << " x " << config_.height << std::endl;
        std::cout << "  Tile size:  " << config_.tileSize << std::endl;
        std::cout << "  Samples:    " << config_.samplesPerPixel << std::endl;

        // Tiles.
        tiles_.Initialise(config_.width, config_.height, config_.tileSize);
        std::cout << "  Tiles:      " << tiles_.GetTileCount() << " (spiral order)" << std::endl;

        // Display buffer.
        display_.Initialise(config_.width, config_.height);

        // Progress tracker (Start() before SetTotals to avoid a reset).
        progress_.Start();
        progress_.SetTotals(tiles_.GetTileCount(), config_.samplesPerPixel);

        // Construct the spacetime from its registry identity. Charge and the
        // cosmological constant flow through; a metric the CPU tracer cannot
        // represent leaves metric_ null, and the scheduler then refuses the CPU
        // path instead of substituting a different spacetime.
        switch (config_.metricId) {
            case MetricId::Minkowski:
                metric_ = std::make_unique<KerrSchildFamily>(KerrSchildParams::Minkowski());
                break;
            case MetricId::Schwarzschild:
                if (config_.blackHoleSpin != 0.0 || config_.blackHoleCharge != 0.0 ||
                    config_.cosmologicalConstant != 0.0) {
                    throw std::invalid_argument(
                        "SessionConfig: Schwarzschild received spin, charge, or lambda");
                }
                metric_ = std::make_unique<KerrSchildFamily>(
                    KerrSchildParams::Schwarzschild(config_.blackHoleMass));
                break;
            case MetricId::Kerr:
                if (config_.blackHoleCharge != 0.0 || config_.cosmologicalConstant != 0.0) {
                    throw std::invalid_argument("SessionConfig: Kerr received charge or lambda");
                }
                metric_ = std::make_unique<KerrSchildFamily>(KerrSchildParams::Kerr(
                    config_.blackHoleMass, config_.blackHoleSpin * config_.blackHoleMass));
                break;
            case MetricId::ReissnerNordstrom:
                if (config_.blackHoleSpin != 0.0 || config_.cosmologicalConstant != 0.0) {
                    throw std::invalid_argument(
                        "SessionConfig: Reissner-Nordstrom received spin or lambda");
                }
                metric_ = std::make_unique<KerrSchildFamily>(KerrSchildParams::ReissnerNordstrom(
                    config_.blackHoleMass, config_.blackHoleCharge * config_.blackHoleMass));
                break;
            case MetricId::KerrNewman:
                if (config_.cosmologicalConstant != 0.0) {
                    throw std::invalid_argument("SessionConfig: Kerr-Newman received lambda");
                }
                metric_ = std::make_unique<KerrSchildFamily>(KerrSchildParams::KerrNewman(
                    config_.blackHoleMass, config_.blackHoleSpin * config_.blackHoleMass,
                    config_.blackHoleCharge * config_.blackHoleMass));
                break;
            case MetricId::DeSitter: {
                if (config_.blackHoleSpin != 0.0 || config_.blackHoleCharge != 0.0) {
                    throw std::invalid_argument("SessionConfig: de-Sitter received spin or charge");
                }
                metric_ = std::make_unique<KerrSchildFamily>(
                    KerrSchildParams::DeSitter(config_.cosmologicalConstant));
                break;
            }
            case MetricId::SchwarzschildDeSitter: {
                if (config_.blackHoleSpin != 0.0 || config_.blackHoleCharge != 0.0) {
                    throw std::invalid_argument(
                        "SessionConfig: Schwarzschild-de-Sitter received spin or charge");
                }
                KerrSchildParams params = KerrSchildParams::Schwarzschild(config_.blackHoleMass);
                params.Lambda = config_.cosmologicalConstant;
                metric_ = std::make_unique<KerrSchildFamily>(params);
                break;
            }
            case MetricId::Alcubierre: {
                WarpDriveParams wp;
                wp.vs = config_.warpVelocity;
                wp.R = config_.bubbleRadius;
                wp.sigma = config_.bubbleSigma;
                metric_ = std::make_unique<WarpDriveFamily>(wp);
                break;
            }
            case MetricId::MorrisThorne: {
                // The Cartesian embedding of the spherical family; one sheet,
                // throat as the capture surface (see morris_thorne_family.h).
                metric_ = std::make_unique<core::MorrisThorneCartesian>(
                    core::MorrisThorneParams::Ellis(config_.throatRadius));
                break;
            }
        }
        if (metric_) {
            std::cout << "  Metric:     " << metric_->GetName() << std::endl;
        } else {
            std::cout << "  Metric:     " << core::MetricInfoFor(config_.metricId).canonical_name
                      << " (GPU backend only)" << std::endl;
        }

        // Camera.
        CameraConfig camConfig;
        camConfig.r = config_.observerDistance;
        camConfig.theta = config_.observerInclination;
        camConfig.phi = config_.observerAzimuth;
        camConfig.fov = config_.cameraFOV;
        camConfig.width = config_.width;
        camConfig.height = config_.height;
        // Camera worldline aberration (P5): beta components in the local
        // ray-component frame. Zero keeps every ray unaberrated (byte-pin).
        camConfig.beta_x = config_.cameraBetaForward;
        camConfig.beta_y = config_.cameraBetaUp;
        camConfig.beta_z = config_.cameraBetaRight;
        camConfig.focal_length = config_.cameraFocalLength;
        camConfig.aperture = config_.cameraAperture;
        camConfig.focus_distance = config_.cameraFocusDistance;
        camera_ = core::CreateCamera(config_.lensType, camConfig);
        std::cout << "  Observer:   r=" << camConfig.r
                  << "M, theta=" << (camConfig.theta * 180.0 / math::kPi) << " deg" << std::endl;

        // Geodesic tracer.
        TracerConfig tracerConfig;
        tracerConfig.escape_radius = 200.0f;
        // Kerr-Schild coordinates are horizon-penetrating, so the exact capture
        // surface is numerically safe. Enlarging it inflates the near-extremal
        // shadow.
        tracerConfig.horizon_factor = 1.0f;
        tracerConfig.max_steps = 20000;

        // Large steps far from the hole, small near the horizon.
        tracerConfig.integrator.initial_step = 0.1f;
        tracerConfig.integrator.max_step = 2.0f;
        tracerConfig.integrator.min_step = 1e-5f;
        tracerConfig.integrator.abs_tolerance = 5e-6f;
        tracerConfig.integrator.rel_tolerance = 5e-6f;

        // Disk inner edge from the ISCO. The thin disk is a black-hole construct;
        // horizonless spacetimes render lensing and background only.
        const auto diskSupport = core::DiskSupportFor(config_.metricId);
        const bool diskCapable = diskSupport == core::DiskSupport::PageThorne;
        if (config_.enableDisk && diskSupport != core::DiskSupport::PageThorne) {
            throw std::invalid_argument(
                "SessionConfig: the selected metric has no represented accretion-disk emission "
                "model; disable the disk");
        }
        if (!config_.enableDisk &&
            (config_.enableVolumetricDisk || config_.enableTurbulence || config_.enableCorona)) {
            throw std::invalid_argument(
                "SessionConfig: volumetric disk, turbulence, and corona require the disk");
        }
        auto rISCO = AccretionDiskD::ComputeIsco(config_.blackHoleSpin);
        tracerConfig.disk_inner = static_cast<float>(rISCO * config_.blackHoleMass);
        tracerConfig.disk_outer = static_cast<float>(20.0 * config_.blackHoleMass);
        tracerConfig.enable_disk = config_.enableDisk && diskCapable;
        tracerConfig.enable_polarisation = config_.enablePolarisation;

        // Doppler beaming toggle (P4); true keeps the full physics and the
        // pinned render byte-for-byte.
        tracerConfig.doppler_beaming = config_.dopplerBeaming;
        switch (config_.temperatureModel) {
            case DiskTemperatureModel::NovikovThorne:
                tracerConfig.disk_temperature_model = backend::DiskTemperatureModel::NovikovThorne;
                break;
            case DiskTemperatureModel::ShakuraSunyaev:
                tracerConfig.disk_temperature_model = backend::DiskTemperatureModel::ShakuraSunyaev;
                break;
            default:
                throw std::invalid_argument(
                    "SessionConfig: invalid disk temperature model reached CPU initialisation");
        }

        // Ray bundles (P2/P3). Enabled for the beam footprint the filtered star
        // field consumes; the pupil (point-source) mode gives the celestial-sphere
        // footprint. Default off leaves the point-sampled path and the pinned
        // render untouched.
        pixel_angular_size_ = (config_.cameraFOV * math::kPi / 180.0) / std::max(1, config_.height);
        tracerConfig.enable_ray_bundles = config_.rayBundles;
        tracerConfig.bundle_point_source = true;
        tracerConfig.bundle_angular_size = static_cast<float>(pixel_angular_size_);

        // Volumetric disk configuration.
        tracerConfig.enable_volumetric = config_.enableVolumetricDisk;
        tracerConfig.volumetric_H_over_r = config_.volumetricHOverR;
        tracerConfig.volumetric_H_power = config_.volumetricHPower;
        tracerConfig.volumetric_tau_midplane = config_.volumetricTauMidplane;
        tracerConfig.volumetric_samples = config_.volumetricSamples;
        tracerConfig.enable_turbulence = config_.enableTurbulence;
        tracerConfig.turbulence.enabled = config_.enableTurbulence;
        tracerConfig.enable_corona = config_.enableCorona;
        tracerConfig.corona.enabled = config_.enableCorona;
        tracerConfig.corona.inner_radius_M = tracerConfig.disk_inner;
        tracerConfig.corona.outer_radius_M = tracerConfig.disk_outer;

        if (metric_) {
            tracer_ = std::make_unique<GeodesicTracer>(metric_.get(), tracerConfig);
            if (diskCapable) {
                std::cout << "  Disk:       r_in=" << tracerConfig.disk_inner
                          << "M, r_out=" << tracerConfig.disk_outer << "M" << std::endl;
            }
            std::cout << "[Session] Physics engine initialized (geodesic integration enabled)"
                      << std::endl;
        }

        // Relativistic jets.
        if (config_.enableJets) {
            JetConfig jetConfig;
            jetConfig.lorentz_factor = config_.jetLorentzFactor;
            jetConfig.opening_angle = config_.jetOpeningAngle;
            jetConfig.r_launch = config_.jetLaunchRadius;
            jetConfig.r_max = config_.jetMaxExtent;
            jetConfig.collimation = config_.jetCollimation;
            jetConfig.spectral_index = config_.jetSpectralIndex;
            jet_ = std::make_unique<RelativisticJet>(jetConfig);
            std::cout << "[Session] Relativistic jets enabled: Gamma=" << jetConfig.lorentz_factor
                      << ", theta_open=" << (jetConfig.opening_angle * 180.0 / math::kPi) << " deg"
                      << std::endl;
        }

        // Colour mode.
        const char* modeName = "TrueColor";
        switch (config_.colorMode) {
            case core::color_modes::Mode::TrueColor:
                modeName = "TrueColor (Physical)";
                break;
            case core::color_modes::Mode::TemperatureMap:
                modeName = "TemperatureMap (False Color)";
                break;
            case core::color_modes::Mode::RedshiftMap:
                modeName = "RedshiftMap (g-factor)";
                break;
            case core::color_modes::Mode::Narrowband:
                modeName = "Narrowband (Hubble Palette)";
                break;
            case core::color_modes::Mode::Polarisation:
                modeName = "Polarisation";
                break;
        }
        std::cout << "[Session] Color mode: " << modeName << std::endl;

        // Multi-threaded rendering setup.
        if (config_.enableParallelRendering) {
            num_threads_ = config_.threadCount;
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
                        std::make_unique<GeodesicTracer>(metric_.get(), tracerConfig));
                }
            }

            std::cout << "[Session] Multi-threaded rendering enabled: " << num_threads_
                      << " threads" << std::endl;
        } else {
            num_threads_ = 1;
            std::cout << "[Session] Single-threaded rendering" << std::endl;
        }

        // The background texture is a physics input, not decorative packaging:
        // substituting grey changes every escaped ray. Missing resources
        // therefore fail initialisation instead of quietly changing the scene.
        const auto starfield = base::ResolveResource("assets/Starfield.png");
        if (!starfield || !LoadStarfieldTexture(starfield->string())) {
            throw std::runtime_error(
                "required runtime resource assets/Starfield.png is missing or unreadable "
                "(set SIRIUS_RESOURCE_DIR to an installed share/sirius directory)");
        }

        // Filtered point-source star field (P3): build the deterministic catalogue
        // once. The beam footprint (ray bundles on) or a pinhole sigma (bundles
        // off) filters it per escaping ray.
        if (config_.pointStarfield) {
            core::StarfieldConfig scfg = config_.starfieldConfig;
            star_generator_ = std::make_unique<core::StarfieldGenerator>(scfg);
            star_catalogue_ = star_generator_->GenerateCatalogue();
            star_index_ = std::make_unique<core::StarfieldSpatialIndex>(star_catalogue_);
            std::cout << "[Session] Point-source star field: " << star_catalogue_.size()
                      << " stars, " << (star_index_->MemoryBytes() / 1024) << " KiB index, beams "
                      << (config_.rayBundles ? "on" : "off") << std::endl;
        }

        // GPU acceleration removed: the legacy OptiX backend init and starfield
        // upload lived here. OptiX is retired; the Vulkan compute path arrives
        // through sirius::backend::device later. The CPU path renders directly.

        if (progress_.GetCancellationToken().IsCancelled()) {
            fsm_.Process(SessionEvent::Cancel);
            return;
        }

        // Transition to Scheduling.
        fsm_.Process(SessionEvent::Ready);
    } catch (const std::exception& e) {
        error_message_ = e.what();
        fsm_.Process(progress_.GetCancellationToken().IsCancelled() ? SessionEvent::Cancel
                                                                    : SessionEvent::Error);
    }
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
                         std::string(core::MetricInfoFor(config_.metricId).canonical_name) +
                         "' is not supported on the CPU backend (GPU/Vulkan required)";
        std::cerr << "[Session] " << error_message_ << std::endl;
        fsm_.Process(SessionEvent::Error);
        return;
    }

    // Parallel rendering when enabled and multiple threads are available.
    if (config_.enableParallelRendering && num_threads_ > 1) {
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
        float vol_intensity = result.volumetric_emission[0];
        float T_effective = std::pow(vol_intensity, 0.25f);
        float T_kelvin = std::clamp(T_effective * config_.diskTemperatureScale, 1000.0f, 100000.0f);

        core::spectral::Rgb bbColor = core::spectral::BlackbodyToRgb(static_cast<double>(T_kelvin));
        float mag = result.magnification;
        px.r = bbColor.r * result.volumetric_emission[0] * mag;
        px.g = bbColor.g * result.volumetric_emission[1] * mag;
        px.b = bbColor.b * result.volumetric_emission[2] * mag;

        if (result.photon_ring) {
            px.r *= kPhotonRingBoostVolumetric;
            px.g *= kPhotonRingBoostVolumetric;
            px.b *= kPhotonRingBoostVolumetric;
        }
        return px;
    }

    // Thin disk: accumulate emission from all disk crossings. Relativistic
    // beaming is applied exactly once: emitted T^4 becomes observed g^4 T^4.
    float total_r = 0.0f, total_g = 0.0f, total_b = 0.0f;
    core::StokesVector total_stokes;
    const bool polarisation_mode = config_.colorMode == core::color_modes::Mode::Polarisation;

    for (int crossing_idx = 0; crossing_idx < result.num_disk_crossings; crossing_idx++) {
        const auto& crossing = result.disk_crossings[crossing_idx];
        if (!crossing.valid) continue;

        float T_emit = crossing.temperature;
        float g = crossing.redshift;
        float r_cross = crossing.r;

        // Higher-order demagnification ~exp(-n pi).
        float order_demag = std::exp(-static_cast<float>(math::kPi) * crossing_idx);

        // Motion blur (primary crossing only). Retain the temporal samples:
        // radiance and colour are nonlinear in g, so averaging g itself is not
        // a represented temporal integral.
        std::vector<float> temporal_redshifts{g};
        if (crossing_idx == 0 && config_.enableMotionBlur && config_.motionBlurSamples > 1) {
            float grav = result.gfactor_grav;
            float gamma = result.gfactor_gamma;
            float v_orb = result.gfactor_v_orb;
            float A = result.gfactor_A;
            float B = result.gfactor_B;

            float M = static_cast<float>(config_.blackHoleMass);
            float a = static_cast<float>(config_.blackHoleSpin * M);
            float sqrtM = std::sqrt(M);
            float Omega = sqrtM / (std::pow(r_cross, 1.5f) + a * sqrtM);

            float delta_phi_max = Omega * config_.shutterTime;
            int N = config_.motionBlurSamples;
            temporal_redshifts.clear();
            temporal_redshifts.reserve(static_cast<std::size_t>(N));

            for (int i = 0; i < N; i++) {
                float t = (N > 1) ? static_cast<float>(i) / (N - 1) : 0.5f;
                float delta_phi = delta_phi_max * (t - 0.5f);

                float cos_dphi = std::cos(delta_phi);
                float sin_dphi = std::sin(delta_phi);
                float v_dot_n_offset = v_orb * (A * cos_dphi + B * sin_dphi);
                float doppler_denom = gamma * (1.0f - v_dot_n_offset);
                doppler_denom = std::clamp(doppler_denom, kDopplerClampMin, kDopplerClampMax);

                float g_offset = grav / doppler_denom;
                g_offset = std::clamp(g_offset, kGFactorClampMin, kGFactorClampMax);
                temporal_redshifts.push_back(g_offset);
            }
        }

        const float emitted_intensity = std::pow(T_emit, 4.0f);
        const float g2 = g * g;
        const float observed_intensity = emitted_intensity * g2 * g2;

        // Limb darkening (primary crossing only). It is scalar for the
        // Chandrasekhar atmosphere, so it applies coherently to I, Q, and U.
        float limb_scale = 1.0f;
        if (crossing_idx == 0) {
            float A = result.gfactor_A;
            float B = result.gfactor_B;
            float n_xy_squared = A * A + B * B;
            float cos_theta = std::sqrt(std::max(0.0f, 1.0f - n_xy_squared));
            const core::spectral::Rgb darkened = core::spectral::ApplyLimbDarkening(
                core::spectral::Rgb(1.0f, 1.0f, 1.0f), cos_theta);
            limb_scale = darkened.r;
        }

        if (polarisation_mode) {
            SIRIUS_ASSERT(crossing.polarisation_valid);
            if (!crossing.polarisation_valid) continue;

            const float chi = crossing.polarisation_evpa;
            const float degree = crossing.polarisation_degree;
            core::StokesVector crossing_stokes{
                observed_intensity, observed_intensity * degree * std::cos(2.0f * chi),
                observed_intensity * degree * std::sin(2.0f * chi), 0.0f};
            crossing_stokes *= limb_scale * order_demag;
            total_stokes += crossing_stokes;
            continue;
        }

        const core::spectral::Rgb diskColor = core::color_modes::AverageTemporalColorMode(
            config_.colorMode, T_emit, temporal_redshifts, emitted_intensity,
            config_.diskTemperatureScale);

        total_r += diskColor.r * limb_scale * order_demag;
        total_g += diskColor.g * limb_scale * order_demag;
        total_b += diskColor.b * limb_scale * order_demag;
    }

    // Gravitational lensing magnification and cinematic boost.
    float output_scale = result.magnification * kDiskIntensityBoost;
    if (result.photon_ring) {
        output_scale *= kPhotonRingBoostDisk;
    }

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
    if (config_.pointStarfield && !star_catalogue_.empty()) {
        SampleStarfieldPoints(result, px.r, px.g, px.b);
    } else {
        SampleStarfield(result.final_direction, px.r, px.g, px.b);
    }

    if (result.magnification > 1.0f) {
        px.r *= result.magnification;
        px.g *= result.magnification;
        px.b *= result.magnification;
    }

    if (result.photon_ring) {
        px.r *= kPhotonRingBoostEscaped;
        px.g *= kPhotonRingBoostEscaped;
        px.b *= kPhotonRingBoostEscaped;
    }

    // Relativistic jet emission.
    if (jet_ && config_.enableJets) {
        double obs_r = config_.observerDistance;
        double obs_theta = config_.observerInclination;
        float obs_x = static_cast<float>(obs_r * std::sin(obs_theta));
        float obs_y = 0.0f;
        float obs_z = static_cast<float>(obs_r * std::cos(obs_theta));

        float end_x = static_cast<float>(result.final_position(1));
        float end_y = static_cast<float>(result.final_position(2));
        float end_z = static_cast<float>(result.final_position(3));

        float jet_emission = core::jet_ray_marching::IntegrateJetEmission(
            *jet_, obs_x, obs_y, obs_z, end_x, end_y, end_z, obs_x, obs_y, obs_z, 64);

        if (jet_emission > 0.0f) {
            float jet_scale = config_.jetIntensity * kJetNormalisation;
            px.r += jet_emission * jet_scale * 0.8f;
            px.g += jet_emission * jet_scale * 0.9f;
            px.b += jet_emission * jet_scale * 1.0f;
        }
    }

    return px;
}

RenderSession::PixelResult RenderSession::ShadePixel(int px_coord, int py_coord,
                                                     GeodesicTracer* tracer) const {
    PixelResult result;
    float r_acc = 0.0f, g_acc = 0.0f, b_acc = 0.0f;

    const int samples_taken = ForEachPixelSample(config_.samplesPerPixel, [&](float u, float v) {
        CameraRay camRay = camera_->GenerateRayAberrated(px_coord, py_coord, u, v);
        TraceResult traceResult = tracer->Trace(camRay);

        float sr = 0.0f, sg = 0.0f, sb = 0.0f;

        switch (traceResult.outcome) {
            case TraceResult::Outcome::Horizon:
                break;

            case TraceResult::Outcome::DiskHit: {
                PixelResult disk = ShadeDiskHit(traceResult);
                sr = disk.r;
                sg = disk.g;
                sb = disk.b;
                break;
            }

            case TraceResult::Outcome::Escaped: {
                PixelResult esc = ShadeEscaped(traceResult);
                sr = esc.r;
                sg = esc.g;
                sb = esc.b;
                break;
            }

            case TraceResult::Outcome::Spiraling:
                sr = sg = sb = kSpiralingBrightness;
                break;

            case TraceResult::Outcome::MaxSteps:
                sr = sg = sb = kMaxStepsBrightness;
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
        if (traceResult.volumetric_hit) {
            PixelResult volume = ShadeDiskHit(traceResult);
            const float transmission = std::exp(-std::max(traceResult.optical_depth, 0.0f));
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

    const size_t size =
        static_cast<size_t>(starfield_width_) * static_cast<size_t>(starfield_height_) * 4;
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
    if (config_.rayBundles && result.beam.valid) {
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
                                           orientation, star_catalogue_, *star_index_, r, g, b);
}

// =============================================================================
// Output writing.
// =============================================================================
void RenderSession::WriteOutput() {
    try {
        std::cout << "\n[Session] Writing output..." << std::endl;

        const int width = config_.width;
        const int height = config_.height;
        const size_t pixelCount = static_cast<size_t>(width) * height;

        if (const auto bad = display_.FirstNonFiniteIndex(); bad.has_value()) {
            throw std::runtime_error("linear radiance contains a non-finite sample at channel " +
                                     std::to_string(*bad));
        }

        // Format follows the output extension: .exr | .png | .ppm (default).
        std::string ext;
        {
            size_t dot = config_.outputPath.rfind('.');
            if (dot != std::string::npos) ext = config_.outputPath.substr(dot + 1);
            std::transform(ext.begin(), ext.end(), ext.begin(),
                           [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
        }

        if (config_.writeOutput && ext == "exr") {
            // EXR is the HDR interchange path: the linear radiance buffer is
            // written untouched. Tonemapping, grading, and transfer encoding are
            // display concerns and deliberately do not apply here.
            ImageBufferRGBA hdr(width, height);
            std::memcpy(hdr.pixels.data(), display_.GetFloatData(), pixelCount * 4 * sizeof(float));

            EXRMetadata meta;
            meta.metricType = core::MetricInfoFor(config_.metricId).canonical_name;
            meta.blackHoleMass = config_.blackHoleMass;
            meta.blackHoleSpin = config_.blackHoleSpin;
            meta.observerDistance = config_.observerDistance;
            meta.observerInclination = config_.observerInclination * 180.0 / math::kPi;
            meta.samplesPerPixel = config_.samplesPerPixel;

            if (!EXRWriter::WriteExr(config_.outputPath, hdr, meta)) {
                throw std::runtime_error("Failed to write EXR: " + config_.outputPath);
            }
        } else {
            // Display pipeline: tonemap and grade to display-linear values. The
            // transfer encode belongs to the writer, so it happens exactly once
            // per format (PPM applies gamma below; PNG applies sRGB inside
            // PNGWriter).
            core::PostProcessConfig ppConfig;
            ppConfig.tonemapper = config_.tonemapper;
            ppConfig.exposure = config_.exposure;
            ppConfig.gamma = 1.0f;  // Writers encode.

            ppConfig.enable_bloom = config_.enableBloom;
            ppConfig.bloom_intensity = config_.bloomIntensity;
            ppConfig.bloom_threshold = config_.bloomThreshold;
            ppConfig.bloom_radius = kBloomRadius;

            ppConfig.saturation = config_.saturation;
            ppConfig.contrast = config_.contrast;
            ppConfig.lift = kShadowLift;
            ppConfig.gain = 1.0f;

            core::PostProcessor::Process(display_.GetFloatBuffer(), width, height, ppConfig);

            // Film is a display-referred finishing pipeline. It is intentionally
            // absent from the EXR branch above, which must retain untouched
            // linear HDR radiance.
            if (config_.enableFilmSimulation) {
                FilmPipeline film(config_.filmConfig);
                film.Apply(display_.GetFloatBuffer().data(), width, height, 0);
            }

            if (const auto bad = display_.FirstNonFiniteIndex(); bad.has_value()) {
                throw std::runtime_error(
                    "display pipeline produced a non-finite sample at channel " +
                    std::to_string(*bad));
            }

            const float* data = display_.GetFloatData();

            if (!config_.writeOutput) {
                // In-memory progressive previews consume the display-linear
                // buffer directly and must not overwrite the configured output.
            } else if (ext == "png") {
                if (!PNGWriter::WriteRgba(config_.outputPath, width, height, data)) {
                    throw std::runtime_error("Failed to write PNG: " + config_.outputPath);
                }
            } else {
                if (ext != "ppm") {
                    std::cout << "[Session] Unrecognised output extension '" << ext
                              << "'; writing PPM data" << std::endl;
                }
                std::ofstream file(config_.outputPath, std::ios::binary);
                if (!file) {
                    throw std::runtime_error("Failed to open output file: " + config_.outputPath);
                }

                file << "P6\n" << width << " " << height << "\n255\n";

                constexpr float kDisplayGamma = 2.2f;
                for (int y = 0; y < height; ++y) {
                    for (int x = 0; x < width; ++x) {
                        int idx = (y * width + x) * 4;
                        unsigned char rgb[3];
                        for (int c = 0; c < 3; ++c) {
                            float val = std::clamp(data[idx + c], 0.0f, 1.0f);
                            val = std::pow(val, 1.0f / kDisplayGamma);
                            rgb[c] =
                                static_cast<unsigned char>(std::clamp(val * 255.0f, 0.0f, 255.0f));
                        }
                        file.write(reinterpret_cast<char*>(rgb), 3);
                    }
                }
            }
        }

        if (config_.writeOutput) {
            std::cout << "[Session] Wrote: " << config_.outputPath << std::endl;
        }

        fsm_.Process(SessionEvent::OutputWritten);
    } catch (const std::exception& e) {
        error_message_ = e.what();
        fsm_.Process(SessionEvent::Error);
    }
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
            message = "Render failed: " + error_message_;
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

// =============================================================================
// Configuration conversion.
// =============================================================================
SessionConfig SessionConfig::FromSiriusConfig(const SiriusConfig& config) {
    SessionConfig sc;
    sc.width = config.render.width;
    sc.height = config.render.height;
    sc.samplesPerPixel = config.render.samplesPerPixel;
    sc.tileSize = config.render.tileSize;
    sc.threadCount = config.render.threadCount;
    sc.outputPath = config.render.outputPath;
    sc.enableDisk = config.diskEnabled;

    // The validator has already accepted the name; a parse failure here is an
    // invariant violation, so it halts rather than substituting a default.
    auto metricId = core::ParseMetricName(config.metric.name);
    if (!metricId.has_value()) {
        throw std::invalid_argument("SessionConfig: unvalidated metric name '" +
                                    config.metric.name + "' reached the session boundary");
    }
    sc.metricId = *metricId;
    sc.blackHoleMass = config.metric.mass;
    // The config validator is the authority on the spin range; no silent clamp
    // is duplicated here.
    sc.blackHoleSpin = config.metric.spin;
    // Tonemapper string was validated at the config boundary; parse here.
    {
        const std::string& tm = config.postprocess.tonemapper;
        if (tm == "ACES")
            sc.tonemapper = core::TonemapType::Aces;
        else if (tm == "Reinhard")
            sc.tonemapper = core::TonemapType::Reinhard;
        else if (tm == "Filmic" || tm == "Uncharted2")
            sc.tonemapper = core::TonemapType::Filmic;
        else if (tm == "None" || tm == "Linear")
            sc.tonemapper = core::TonemapType::None;
        else
            throw std::invalid_argument("SessionConfig: unvalidated tonemapper '" + tm +
                                        "' reached the session boundary");
    }
    sc.blackHoleCharge = config.metric.charge;
    sc.cosmologicalConstant = config.metric.lambda;
    if (config.metric.temperatureModel == "NovikovThorne" ||
        config.metric.temperatureModel == "NT") {
        sc.temperatureModel = DiskTemperatureModel::NovikovThorne;
    } else if (config.metric.temperatureModel == "ShakuraSunyaev" ||
               config.metric.temperatureModel == "SS") {
        sc.temperatureModel = DiskTemperatureModel::ShakuraSunyaev;
    } else {
        throw std::invalid_argument("SessionConfig: unvalidated disk temperature model '" +
                                    config.metric.temperatureModel +
                                    "' reached the session boundary");
    }
    sc.diskTemperatureScale = config.metric.diskTemperature;
    sc.dopplerBeaming = config.dopplerBeaming;
    sc.pointStarfield = config.pointStarfield;
    sc.rayBundles = config.rayBundles;
    if (config.colorMode == "TrueColor") {
        sc.colorMode = core::color_modes::Mode::TrueColor;
    } else if (config.colorMode == "TemperatureMap") {
        sc.colorMode = core::color_modes::Mode::TemperatureMap;
    } else if (config.colorMode == "RedshiftMap") {
        sc.colorMode = core::color_modes::Mode::RedshiftMap;
    } else if (config.colorMode == "Narrowband") {
        sc.colorMode = core::color_modes::Mode::Narrowband;
    } else if (config.colorMode == "Polarisation") {
        sc.colorMode = core::color_modes::Mode::Polarisation;
    } else {
        throw std::invalid_argument("SessionConfig: unvalidated color mode '" + config.colorMode +
                                    "' reached the session boundary");
    }
    sc.enablePolarisation = sc.colorMode == core::color_modes::Mode::Polarisation;
    sc.throatRadius = config.metric.throatRadius;
    if (config.metric.wormholeTopology == "OneSheetCapture") {
        sc.wormholeTopology = WormholeTopology::OneSheetCapture;
    } else if (config.metric.wormholeTopology == "TwoSheet") {
        sc.wormholeTopology = WormholeTopology::TwoSheet;
    } else {
        throw std::invalid_argument("SessionConfig: unvalidated wormhole topology '" +
                                    config.metric.wormholeTopology +
                                    "' reached the session boundary");
    }
    sc.warpVelocity = config.metric.warpVelocity;
    sc.bubbleRadius = config.metric.bubbleRadius;
    sc.bubbleSigma = config.metric.bubbleSigma;
    sc.observerDistance = config.observer.distance;
    sc.observerInclination = config.observer.inclination * math::kPi / 180.0;
    sc.observerAzimuth = config.observer.azimuth * math::kPi / 180.0;
    sc.cameraFOV = static_cast<float>(config.observer.fov);
    sc.cameraBetaForward = config.observer.cameraBetaForward;
    sc.cameraBetaUp = config.observer.cameraBetaUp;
    sc.cameraBetaRight = config.observer.cameraBetaRight;
    if (config.observer.lensModel == "Pinhole") {
        sc.lensType = core::LensType::Pinhole;
    } else if (config.observer.lensModel == "ThinLens") {
        sc.lensType = core::LensType::ThinLens;
    } else if (config.observer.lensModel == "Fisheye") {
        sc.lensType = core::LensType::Fisheye;
    } else {
        throw std::invalid_argument("SessionConfig: unvalidated lens model '" +
                                    config.observer.lensModel + "' reached the session boundary");
    }
    sc.cameraFocalLength = config.observer.focalLength;
    sc.cameraAperture = config.observer.aperture;
    sc.cameraFocusDistance = config.observer.focusDistance;

    // Post-processing.
    sc.enableBloom = config.postprocess.enableBloom;
    sc.bloomIntensity = config.postprocess.bloomIntensity;
    sc.bloomThreshold = config.postprocess.bloomThreshold;
    sc.exposure = config.postprocess.exposure;
    sc.contrast = config.postprocess.contrast;
    sc.saturation = config.postprocess.saturation;

    // Volumetric disk.
    sc.enableVolumetricDisk = config.volumetric.enabled;
    sc.volumetricHOverR = config.volumetric.hOverR;
    sc.volumetricHPower = config.volumetric.hPower;
    sc.volumetricTauMidplane = config.volumetric.tauMidplane;
    sc.volumetricSamples = config.volumetric.samples;
    sc.enableTurbulence = config.volumetric.enableTurbulence;
    sc.enableCorona = config.volumetric.enableCorona;

    // Temporal thin-disk integration.
    sc.enableMotionBlur = config.motionBlur.enabled;
    sc.shutterTime = config.motionBlur.shutterTime;
    sc.motionBlurSamples = config.motionBlur.samples;

    // Film simulation.
    sc.enableFilmSimulation = config.film.enabled;
    if (config.film.preset == "Interstellar") {
        sc.filmConfig = FilmConfig::Interstellar();
    } else if (config.film.preset == "SpaceOdyssey2001") {
        sc.filmConfig = FilmConfig::SpaceOdyssey2001();
    } else if (config.film.preset == "DigitalClean") {
        sc.filmConfig = FilmConfig::DigitalClean();
    } else {
        throw std::invalid_argument("SessionConfig: unvalidated film preset '" +
                                    config.film.preset + "' reached the session boundary");
    }
    if (config.film.enabled) {
        sc.filmConfig.grain_intensity = config.film.grainIntensity;
        sc.filmConfig.halation_strength = config.film.halationStrength;
        sc.filmConfig.vignette_strength = config.film.vignetteStrength;
    }

    // Backend selection from the fully layered SiriusConfig. ConfigLoader has
    // already applied file and environment values; CLI parsing, when present,
    // has applied the highest-priority override. GO-LIVE (owner decision,
    // 2026-07-18, specification
    // section 1.5): 'auto' now resolves to Vulkan when the backend is compiled
    // in, a device is present, and the registry marks the metric
    // gpu-dispatchable; otherwise it falls back to the CPU path with the
    // reason logged — a fallback, never a decline. 'cpu' pins the CPU path;
    // 'vulkan' selects the Vulkan path unconditionally (device absence then
    // surfaces as the render's loud decline, not a silent CPU switch).
    sc.backend = RenderBackend::Cpu;
    if (config.backend.preferred == "vulkan") {
        sc.backend = RenderBackend::Vulkan;
    } else if (config.backend.preferred == "auto") {
#ifdef SIRIUS_HAS_VULKAN_BACKEND
        const bool gpu_metric = core::MetricInfoFor(sc.metricId).gpu_supported;
        if (gpu_metric) {
            if (auto compatible = ValidateVulkanRenderConfig(sc); !compatible) {
                std::cout << "[Session] backend auto: " << compatible.error().detail()
                          << "; using the CPU path" << std::endl;
            } else {
                if (auto devices = backend::EnumerateVulkanDevices();
                    devices.has_value() && !devices->empty()) {
                    sc.backend = RenderBackend::Vulkan;
                } else {
                    std::cout << "[Session] backend auto: no Vulkan device visible; "
                                 "using the CPU path"
                              << std::endl;
                }
            }
        } else {
            std::cout << "[Session] backend auto: metric '"
                      << core::MetricInfoFor(sc.metricId).canonical_name
                      << "' is CPU-only (registry); using the CPU path" << std::endl;
        }
#endif
    }

    return sc;
}

}  // namespace sirius::render
