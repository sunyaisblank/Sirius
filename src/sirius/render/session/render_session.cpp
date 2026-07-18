// Render session implementation. Ported from SNRS001A.cpp.
//
// The CPU render path (metric construction, spiral tile scheduling, the
// thread-pool, per-pixel shading, and output writing) is preserved exactly. The
// legacy OptiX/GPU branches are removed; each removal site carries a one-line
// decision comment. The Vulkan compute backend will re-enter through
// sirius::backend::device, not from here.

#include "sirius/render/session/render_session.h"

#include "sirius/render/exr_writer.h"
#include "sirius/render/image_buffer.h"
#include "sirius/render/png_writer.h"

#ifdef SIRIUS_HAS_VULKAN_BACKEND
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
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>

namespace sirius::render {

using backend::GeodesicTracer;
using backend::TraceResult;
using backend::TracerConfig;
using core::AccretionDiskD;
using core::CameraConfig;
using core::CameraRay;
using core::JetConfig;
using core::KerrSchildFamily;
using core::KerrSchildParams;
using core::MetricId;
using core::PinholeCamera;
using core::RelativisticJet;
using core::StokesVector;
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
constexpr float kSynchrotronPolDegree = 0.7f;
constexpr float kSpiralingBrightness = 0.02f;
constexpr float kMaxStepsBrightness = 0.01f;
constexpr float kBackgroundFallback = 0.001f;
constexpr float kTInnerKelvin = 30000.0f;
constexpr float kDopplerClampMin = 0.1f;
constexpr float kDopplerClampMax = 10.0f;
constexpr float kGFactorClampMin = 0.1f;
constexpr float kGFactorClampMax = 5.0f;
constexpr float kJetNormalisation = 0.1f;
constexpr int kBloomRadius = 12;
constexpr float kShadowLift = 0.02f;
}  // namespace

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
        default:
            break;
    }
}

// =============================================================================
// Initialisation.
// =============================================================================
void RenderSession::Initialise() {
    try {
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
            default: {
                KerrSchildParams params;
                params.M = config_.blackHoleMass;
                params.a = config_.blackHoleSpin * config_.blackHoleMass;    // spin given as a/M.
                params.Q = config_.blackHoleCharge * config_.blackHoleMass;  // charge given as Q/M.
                params.Lambda = config_.cosmologicalConstant;
                metric_ = std::make_unique<KerrSchildFamily>(params);
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
        camConfig.phi = 0.0;
        camConfig.fov = config_.cameraFOV;
        camConfig.width = config_.width;
        camConfig.height = config_.height;
        // Camera worldline aberration (P5): beta components in the local
        // ray-component frame. Zero keeps every ray unaberrated (byte-pin).
        camConfig.beta_x = config_.cameraBetaForward;
        camConfig.beta_y = config_.cameraBetaUp;
        camConfig.beta_z = config_.cameraBetaRight;
        camera_ = std::make_unique<PinholeCamera>(camConfig);
        std::cout << "  Observer:   r=" << camConfig.r
                  << "M, theta=" << (camConfig.theta * 180.0 / math::kPi) << " deg" << std::endl;

        // Geodesic tracer.
        TracerConfig tracerConfig;
        tracerConfig.escape_radius = 200.0f;
        tracerConfig.horizon_factor = 1.05f;
        tracerConfig.max_steps = 20000;

        // Large steps far from the hole, small near the horizon.
        tracerConfig.integrator.initial_step = 0.1f;
        tracerConfig.integrator.max_step = 2.0f;
        tracerConfig.integrator.min_step = 1e-5f;
        tracerConfig.integrator.abs_tolerance = 5e-6f;
        tracerConfig.integrator.rel_tolerance = 5e-6f;

        // Disk inner edge from the ISCO. The thin disk is a black-hole construct;
        // horizonless spacetimes render lensing and background only.
        const bool diskCapable = config_.metricId == MetricId::Schwarzschild ||
                                 config_.metricId == MetricId::Kerr ||
                                 config_.metricId == MetricId::ReissnerNordstrom ||
                                 config_.metricId == MetricId::KerrNewman ||
                                 config_.metricId == MetricId::SchwarzschildDeSitter;
        auto rISCO = AccretionDiskD::ComputeIsco(config_.blackHoleSpin);
        tracerConfig.disk_inner = static_cast<float>(rISCO * config_.blackHoleMass);
        tracerConfig.disk_outer = static_cast<float>(20.0 * config_.blackHoleMass);
        tracerConfig.enable_disk = diskCapable;

        // Doppler beaming toggle (P4); true keeps the full physics and the
        // pinned render byte-for-byte.
        tracerConfig.doppler_beaming = config_.dopplerBeaming;

        // Ray bundles (P2/P3). Enabled for the beam footprint the filtered star
        // field consumes; the pupil (point-source) mode gives the celestial-sphere
        // footprint. Default off leaves the point-sampled path and the pinned
        // render untouched.
        pixel_angular_size_ =
            (config_.cameraFOV * math::kPi / 180.0) / std::max(1, config_.height);
        tracerConfig.enable_ray_bundles = config_.rayBundles;
        tracerConfig.bundle_point_source = true;
        tracerConfig.bundle_angular_size = static_cast<float>(pixel_angular_size_);

        // Volumetric disk configuration.
        tracerConfig.enable_volumetric = config_.enableVolumetricDisk;
        tracerConfig.volumetric_H_over_r = config_.volumetricHOverR;
        tracerConfig.volumetric_H_power = config_.volumetricHPower;
        tracerConfig.volumetric_tau_midplane = config_.volumetricTauMidplane;
        tracerConfig.volumetric_samples = config_.volumetricSamples;

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

        // Polarisation buffer.
        if (config_.enablePolarisation) {
            polarisation_buffer_.resize(config_.width * config_.height);
            std::cout << "[Session] Polarisation tracking enabled" << std::endl;
        }

        // Colour mode.
        const char* modeName = "TrueColor";
        switch (config_.colorMode) {
            case core::color_modes::Mode::TrueColor: modeName = "TrueColor (Physical)"; break;
            case core::color_modes::Mode::TemperatureMap: modeName = "TemperatureMap (False Color)"; break;
            case core::color_modes::Mode::RedshiftMap: modeName = "RedshiftMap (g-factor)"; break;
            case core::color_modes::Mode::Narrowband: modeName = "Narrowband (Hubble Palette)"; break;
            case core::color_modes::Mode::Polarisation: modeName = "Polarisation"; break;
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

        // Load the starfield background texture (CPU background sampling needs it).
        std::vector<std::string> texturePaths = {
            "assets/Starfield.png",
            "../assets/Starfield.png",
            "../../assets/Starfield.png",
            "../../../assets/Starfield.png",
            "../../../../assets/Starfield.png"};

        bool textureLoaded = false;
        for (const auto& path : texturePaths) {
            if (LoadStarfieldTexture(path)) {
                textureLoaded = true;
                break;
            }
        }

        if (!textureLoaded) {
            std::cerr << "[Session] Warning: Could not load starfield texture from any path"
                      << std::endl;
            std::cerr << "[Session] Background will render as faint grey" << std::endl;
        }

        // Filtered point-source star field (P3): build the deterministic catalogue
        // once. The beam footprint (ray bundles on) or a pinhole sigma (bundles
        // off) filters it per escaping ray.
        if (config_.pointStarfield) {
            core::StarfieldConfig scfg = config_.starfieldConfig;
            scfg.star_count = std::max(scfg.star_count, 100000u);  // P3 catalogue floor.
            star_generator_ = std::make_unique<core::StarfieldGenerator>(scfg);
            star_catalogue_ = star_generator_->GenerateCatalogue();
            std::cout << "[Session] Point-source star field: " << star_catalogue_.size()
                      << " stars, beams " << (config_.rayBundles ? "on" : "off") << std::endl;
        }

        // GPU acceleration removed: the legacy OptiX backend init and starfield
        // upload lived here. OptiX is retired; the Vulkan compute path arrives
        // through sirius::backend::device later. The CPU path renders directly.

        // Transition to Scheduling.
        fsm_.Process(SessionEvent::Ready);
    } catch (const std::exception& e) {
        error_message_ = e.what();
        fsm_.Process(SessionEvent::Error);
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
        float T_kelvin = std::clamp(T_effective * kTInnerKelvin, 1000.0f, 100000.0f);

        core::spectral::Rgb bbColor = core::spectral::BlackbodyToRgb(static_cast<double>(T_kelvin));
        float mag = result.magnification;
        px.r = bbColor.r * vol_intensity * mag;
        px.g = bbColor.g * vol_intensity * mag;
        px.b = bbColor.b * vol_intensity * mag;

        if (result.photon_ring) {
            px.r *= kPhotonRingBoostVolumetric;
            px.g *= kPhotonRingBoostVolumetric;
            px.b *= kPhotonRingBoostVolumetric;
        }
        return px;
    }

    // Thin disk: accumulate emission from all disk crossings.
    float total_r = 0.0f, total_g = 0.0f, total_b = 0.0f;

    for (int crossing_idx = 0; crossing_idx < result.num_disk_crossings; crossing_idx++) {
        const auto& crossing = result.disk_crossings[crossing_idx];
        if (!crossing.valid) continue;

        float T_emit = crossing.temperature;
        float g = crossing.redshift;
        float r_cross = crossing.r;

        // Higher-order demagnification ~exp(-n pi).
        float order_demag = std::exp(-static_cast<float>(math::kPi) * crossing_idx);
        float T_obs = T_emit * g;
        float intensity = std::pow(T_obs, 4.0f);

        core::spectral::Rgb diskColor =
            core::color_modes::ApplyColorMode(config_.colorMode, T_emit, g, intensity, nullptr);

        float cr = diskColor.r;
        float cg = diskColor.g;
        float cb = diskColor.b;

        // Limb darkening (primary crossing only).
        if (crossing_idx == 0) {
            float A = result.gfactor_A;
            float B = result.gfactor_B;
            float n_xy_squared = A * A + B * B;
            float cos_theta = std::sqrt(std::max(0.0f, 1.0f - n_xy_squared));

            core::spectral::Rgb limbInput(cr, cg, cb);
            core::spectral::Rgb darkened = core::spectral::ApplyLimbDarkening(limbInput, cos_theta);
            cr = darkened.r;
            cg = darkened.g;
            cb = darkened.b;
        }

        // Motion blur (primary crossing only).
        float g_blur = g;
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
            float g_sum = 0.0f;

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
                g_sum += g_offset;
            }

            g_blur = g_sum / N;
        }

        // Relativistic beaming: I_obs = g^4 I_emit.
        float g4 = g_blur * g_blur * g_blur * g_blur;
        cr *= g4;
        cg *= g4;
        cb *= g4;

        total_r += cr * order_demag;
        total_g += cg * order_demag;
        total_b += cb * order_demag;
    }

    // Gravitational lensing magnification and cinematic boost.
    px.r = total_r * result.magnification * kDiskIntensityBoost;
    px.g = total_g * result.magnification * kDiskIntensityBoost;
    px.b = total_b * result.magnification * kDiskIntensityBoost;

    if (result.photon_ring) {
        px.r *= kPhotonRingBoostDisk;
        px.g *= kPhotonRingBoostDisk;
        px.b *= kPhotonRingBoostDisk;
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
    StokesVector stokes_acc;

    int spp = std::max(1, config_.samplesPerPixel);
    int grid_size = static_cast<int>(std::sqrt(static_cast<float>(spp)));
    if (grid_size < 1) grid_size = 1;

    for (int sy = 0; sy < grid_size; ++sy) {
        for (int sx = 0; sx < grid_size; ++sx) {
            float u = (sx + 0.5f) / grid_size;
            float v = (sy + 0.5f) / grid_size;

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
                default:
                    sr = sg = sb = kMaxStepsBrightness;
                    break;
            }

            r_acc += sr;
            g_acc += sg;
            b_acc += sb;

            // Polarisation accumulation.
            if (config_.enablePolarisation && traceResult.outcome == TraceResult::Outcome::DiskHit) {
                float disk_phi = traceResult.disk_phi;
                float evpa = disk_phi + static_cast<float>(math::kHalfPi);
                float I_sample = std::sqrt(sr * sr + sg * sg + sb * sb);
                StokesVector sample_stokes = core::polarised_emission::SynchrotronEmission(
                    I_sample, kSynchrotronPolDegree, evpa);
                stokes_acc += sample_stokes;
            }
        }
    }

    int total_samples = grid_size * grid_size;
    float inv_samples = 1.0f / static_cast<float>(total_samples);
    result.r = r_acc * inv_samples;
    result.g = g_acc * inv_samples;
    result.b = b_acc * inv_samples;

    if (config_.enablePolarisation) {
        stokes_acc *= inv_samples;
        stokes_acc.Normalise();
        result.stokes = stokes_acc;
    }

    return result;
}

// =============================================================================
// Tile rendering.
// =============================================================================
void RenderSession::RenderTile(Tile* tile) {
    if (!tile) return;

    std::vector<float> tileBuffer(tile->width * tile->height * 4, 0.0f);

    for (int ty = 0; ty < tile->height; ++ty) {
        for (int tx = 0; tx < tile->width; ++tx) {
            int px = tile->x + tx;
            int py = tile->y + ty;

            PixelResult pixel = ShadePixel(px, py, tracer_.get());

            int idx = (ty * tile->width + tx) * 4;
            tileBuffer[idx + 0] = pixel.r;
            tileBuffer[idx + 1] = pixel.g;
            tileBuffer[idx + 2] = pixel.b;
            tileBuffer[idx + 3] = 1.0f;

            if (config_.enablePolarisation && !polarisation_buffer_.empty()) {
                int pixel_idx = py * config_.width + px;
                if (pixel_idx >= 0 && pixel_idx < static_cast<int>(polarisation_buffer_.size())) {
                    polarisation_buffer_[pixel_idx] = pixel.stokes;
                }
            }
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
        config_, display_, [this](int /*done*/, int total) {
            // The scheduler's tile count is the CPU spiral grid; the Vulkan path
            // governs its own tile count, so retotal the tracker to it on the
            // first report and then count governed tiles through CompleteTile.
            if (progress_.GetTilesTotal() != total) {
                progress_.SetTotals(total, 1);
            }
            progress_.CompleteTile(1);
        });

    if (!stats) {
        error_message_ = stats.error().Description();
        std::cerr << "[Session] Vulkan render declined: " << error_message_ << std::endl;
        fsm_.Process(SessionEvent::Error);
        return;
    }

    std::cout << "[Session] Vulkan render complete: " << stats->metric_name << " on "
              << stats->device_name << ", " << stats->tiles_rendered << " tile(s) of "
              << stats->tile_plan.tile_edge << "px in " << stats->seconds << "s" << std::endl;
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

    size_t size = static_cast<size_t>(starfield_width_ * starfield_height_ * 4);
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
        r = g = b = kBackgroundFallback;
        return;
    }

    // Direction is Cartesian (x, y, z) from the geodesic tracer.
    double dx = direction(1);
    double dy = direction(2);
    double dz = direction(3);

    double len = std::sqrt(dx * dx + dy * dy + dz * dz);
    if (len < 1e-10) {
        r = g = b = 0.001f;
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
    // Beam footprint radius on the sky. With ray bundles the tracer supplies the
    // lensed footprint; a pinhole (bundles off) samples at a fraction of the pixel
    // angular size, so a star pops in and out as the camera rotates (the flicker
    // the beam filter removes). The footprint is floored at the pixel size so a
    // pixel always integrates at least its own solid angle.
    constexpr float kPinholeFraction = 0.3f;  // Pinhole sigma as a fraction of a pixel.
    float pixel = static_cast<float>(pixel_angular_size_);
    float sigma;
    if (config_.rayBundles && result.beam.valid) {
        float footprint = std::sqrt(std::max(0.0f, result.beam.footprint_major) *
                                    std::max(0.0f, result.beam.footprint_minor));
        sigma = std::max(footprint, pixel);
    } else {
        sigma = kPinholeFraction * pixel;
    }

    const auto& d = result.final_direction;
    star_generator_->AccumulateThroughBeam(static_cast<float>(d(1)), static_cast<float>(d(2)),
                                           static_cast<float>(d(3)), sigma, star_catalogue_, r, g,
                                           b);
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

        // Format follows the output extension: .exr | .png | .ppm (default).
        std::string ext;
        {
            size_t dot = config_.outputPath.rfind('.');
            if (dot != std::string::npos) ext = config_.outputPath.substr(dot + 1);
            std::transform(ext.begin(), ext.end(), ext.begin(),
                           [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
        }

        if (ext == "exr") {
            // EXR is the HDR interchange path: the linear radiance buffer is
            // written untouched. Tonemapping, grading, and transfer encoding are
            // display concerns and deliberately do not apply here.
            ImageBufferRGBA hdr(width, height);
            std::memcpy(hdr.pixels.data(), display_.GetFloatData(),
                        pixelCount * 4 * sizeof(float));

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

            const float* data = display_.GetFloatData();

            if (ext == "png") {
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
                            rgb[c] = static_cast<unsigned char>(std::clamp(val * 255.0f, 0.0f, 255.0f));
                        }
                        file.write(reinterpret_cast<char*>(rgb), 3);
                    }
                }
            }
        }

        std::cout << "[Session] Wrote: " << config_.outputPath << std::endl;

        // Polarisation map: a separate image showing polarisation degree and
        // angle (Hue = EVPA, Saturation = degree, Value = intensity).
        if (config_.outputPolarisationMap && !polarisation_buffer_.empty()) {
            std::string polPath = config_.outputPath;
            size_t dotPos = polPath.rfind('.');
            if (dotPos != std::string::npos) {
                polPath.insert(dotPos, "_polarisation");
            } else {
                polPath += "_polarisation.ppm";
            }

            std::ofstream polFile(polPath, std::ios::binary);
            if (polFile) {
                polFile << "P6\n" << config_.width << " " << config_.height << "\n255\n";

                for (int y = 0; y < config_.height; ++y) {
                    for (int x = 0; x < config_.width; ++x) {
                        int idx = y * config_.width + x;
                        const StokesVector& stokes = polarisation_buffer_[idx];

                        core::spectral::Rgb color =
                            core::color_modes::polarisation_vis::StokesToRgbHsv(stokes);

                        unsigned char rgb[3] = {
                            static_cast<unsigned char>(std::clamp(color.r, 0.0f, 1.0f) * 255.0f),
                            static_cast<unsigned char>(std::clamp(color.g, 0.0f, 1.0f) * 255.0f),
                            static_cast<unsigned char>(std::clamp(color.b, 0.0f, 1.0f) * 255.0f)};
                        polFile.write(reinterpret_cast<char*>(rgb), 3);
                    }
                }

                std::cout << "[Session] Wrote polarisation map: " << polPath << std::endl;
            }
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
            message = "Unknown end state";
    }

    if (completion_callback_) {
        completion_callback_(state, message);
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

    fsm_.Process(SessionEvent::AllTilesComplete);
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

        RenderTileThreaded(tile, thread_id);

        {
            std::lock_guard<std::mutex> lock(tile_mutex_);
            tiles_.CompleteTile(tile->id);
            // Single progress surface: the ProgressTracker callback.
            progress_.CompleteTile(tile->PixelCount());
        }
    }

    active_workers_--;
}

void RenderSession::RenderTileThreaded(Tile* tile, int thread_id) {
    if (!tile) return;

    GeodesicTracer* tracer = nullptr;
    if (thread_id >= 0 && thread_id < static_cast<int>(thread_tracers_.size())) {
        tracer = thread_tracers_[thread_id].get();
    } else {
        tracer = tracer_.get();
    }

    std::vector<float> tileBuffer(tile->width * tile->height * 4, 0.0f);

    for (int ty = 0; ty < tile->height; ++ty) {
        if (ty % 8 == 0 && stop_workers_) {
            return;
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

            if (config_.enablePolarisation && !polarisation_buffer_.empty()) {
                int pixel_idx = py * config_.width + px;
                if (pixel_idx >= 0 && pixel_idx < static_cast<int>(polarisation_buffer_.size())) {
                    std::lock_guard<std::mutex> lock(display_mutex_);
                    polarisation_buffer_[pixel_idx] = pixel.stokes;
                }
            }
        }
    }

    {
        std::lock_guard<std::mutex> lock(display_mutex_);
        display_.UpdateTile(tile->x, tile->y, tile->width, tile->height, tileBuffer.data());
    }
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
    sc.outputPath = config.render.outputPath;

    // The validator has already accepted the name; a parse failure here is an
    // invariant violation, so it halts rather than substituting a default.
    auto metricId = core::ParseMetricName(config.metric.name);
    if (!metricId.has_value()) {
        throw std::invalid_argument("SessionConfig: unvalidated metric name '" + config.metric.name +
                                    "' reached the session boundary");
    }
    sc.metricId = *metricId;
    sc.blackHoleMass = (sc.metricId == MetricId::Minkowski || sc.metricId == MetricId::DeSitter)
                           ? 0.0
                           : config.metric.mass;
    // The config validator is the authority on the spin range; no silent clamp
    // is duplicated here.
    sc.blackHoleSpin = config.metric.spin;
    // Tonemapper string was validated at the config boundary; parse here.
    {
        const std::string& tm = config.postprocess.tonemapper;
        if (tm == "ACES") sc.tonemapper = core::TonemapType::Aces;
        else if (tm == "Reinhard") sc.tonemapper = core::TonemapType::Reinhard;
        else if (tm == "Filmic" || tm == "Uncharted2") sc.tonemapper = core::TonemapType::Filmic;
        else if (tm == "None" || tm == "Linear") sc.tonemapper = core::TonemapType::None;
        else
            throw std::invalid_argument("SessionConfig: unvalidated tonemapper '" + tm +
                                        "' reached the session boundary");
    }
    sc.blackHoleCharge = config.metric.charge;
    sc.cosmologicalConstant = config.metric.lambda;
    sc.temperatureModel = config.metric.temperatureModel;
    sc.diskTemperatureScale = config.metric.diskTemperature;
    sc.dopplerBeaming = config.dopplerBeaming;
    sc.pointStarfield = config.pointStarfield;
    sc.rayBundles = config.rayBundles;
    sc.throatRadius = config.metric.throatRadius;
    sc.warpVelocity = config.metric.warpVelocity;
    sc.bubbleRadius = config.metric.bubbleRadius;
    sc.bubbleSigma = config.metric.bubbleSigma;
    sc.observerDistance = config.observer.distance;
    sc.observerInclination = config.observer.inclination * math::kPi / 180.0;
    sc.cameraFOV = static_cast<float>(config.observer.fov);
    sc.cameraBetaForward = config.observer.cameraBetaForward;
    sc.cameraBetaUp = config.observer.cameraBetaUp;
    sc.cameraBetaRight = config.observer.cameraBetaRight;

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

    // Film simulation.
    sc.enableFilmSimulation = config.film.enabled;
    if (config.film.enabled) {
        if (config.film.preset == "Interstellar") {
            sc.filmConfig = FilmConfig::Interstellar();
        } else if (config.film.preset == "SpaceOdyssey2001") {
            sc.filmConfig = FilmConfig::SpaceOdyssey2001();
        } else {
            sc.filmConfig = FilmConfig::DigitalClean();
        }
        sc.filmConfig.grain_intensity = config.film.grainIntensity;
        sc.filmConfig.halation_strength = config.film.halationStrength;
        sc.filmConfig.vignette_strength = config.film.vignetteStrength;
    }

    // Backend selection. 'auto' resolves to Cpu today (the go-live default flip
    // is an owner milestone, specification section 1.5); SIRIUS_RENDER_BACKEND=
    // vulkan opts into the Vulkan render path explicitly. An explicit --backend
    // vulkan / --gpu request is applied by the CLI after this conversion, because
    // the config validator's backend.preferred set is {auto, cpu} only.
    sc.useGPU = (config.backend.preferred != "cpu");
    sc.backend = RenderBackend::Cpu;
    if (const char* rb = std::getenv("SIRIUS_RENDER_BACKEND"); rb != nullptr) {
        std::string value(rb);
        std::transform(value.begin(), value.end(), value.begin(),
                       [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
        if (value == "vulkan") {
            sc.backend = RenderBackend::Vulkan;
        }
    }

    return sc;
}

}  // namespace sirius::render
