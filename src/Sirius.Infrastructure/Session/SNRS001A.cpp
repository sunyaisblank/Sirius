// SNRS001A.cpp - Render Session Implementation
// Component ID: SNRS001A (Session/RenderSession)
//
// Modified: Priority 1 Fix - Integrate geodesic physics engine
// Previously used geometric approximations; now uses proper geodesic integration.

#include "SNRS001A.h"
#include <PPOP001A.h>
#include <PHSP001A.h>   // Phase 2 P2: Spectral blackbody functions
#include <PHSC001A.h>   // Phase 7: Astronomical color modes
#include <fstream>
#include <cmath>
#include <algorithm>
#include <array>
#include <chrono>
#include <cstring>

// Unified constants (replaces M_PI macro)
#include <PHCN001A.h>

// GPU acceleration factory
namespace Sirius::Acceleration {
    extern std::unique_ptr<IAccelerator> createAccelerator(BackendType type);
    extern std::vector<BackendType> getAvailableBackends();
    extern BackendType getBestBackend();
}

// stb_image for texture loading (implementation defined in RDRT001A.cpp)
#include <stb_image.h>

// Alias for backward compatibility
using namespace Sirius::Constants;

namespace {
    constexpr float DISK_INTENSITY_BOOST = 5.0f;
    constexpr float PHOTON_RING_BOOST_DISK = 2.0f;
    constexpr float PHOTON_RING_BOOST_VOLUMETRIC = 1.5f;
    constexpr float PHOTON_RING_BOOST_ESCAPED = 1.2f;
    constexpr float SYNCHROTRON_POL_DEGREE = 0.7f;
    constexpr float SPIRALING_BRIGHTNESS = 0.02f;
    constexpr float MAX_STEPS_BRIGHTNESS = 0.01f;
    constexpr float BACKGROUND_FALLBACK = 0.001f;
    constexpr float T_INNER_KELVIN = 30000.0f;
    constexpr float DOPPLER_CLAMP_MIN = 0.1f;
    constexpr float DOPPLER_CLAMP_MAX = 10.0f;
    constexpr float G_FACTOR_CLAMP_MIN = 0.1f;
    constexpr float G_FACTOR_CLAMP_MAX = 5.0f;
    constexpr float JET_NORMALISATION = 0.1f;
    constexpr int   BLOOM_RADIUS = 12;
    constexpr float SHADOW_LIFT = 0.02f;
}

namespace Sirius {

//==============================================================================
// Initialisation
//==============================================================================
void RenderSession::initialise() {
    try {
        std::cout << "[Session] Initialising render..." << std::endl;
        std::cout << "  Resolution: " << m_Config.width << " x " << m_Config.height << std::endl;
        std::cout << "  Tile size:  " << m_Config.tileSize << std::endl;
        std::cout << "  Samples:    " << m_Config.samplesPerPixel << std::endl;

        // Initialise tiles
        m_Tiles.initialise(m_Config.width, m_Config.height, m_Config.tileSize);
        std::cout << "  Tiles:      " << m_Tiles.getTileCount() << " (spiral order)" << std::endl;

        // Initialise display buffer
        m_Display.initialise(m_Config.width, m_Config.height);

        // Initialise progress tracker (start() before setTotals to avoid reset)
        m_Progress.start();
        m_Progress.setTotals(m_Tiles.getTileCount(), m_Config.samplesPerPixel);

        // =================================================================
        // Priority 1 Fix: Initialize physics components
        // =================================================================

        // Create metric from config (Kerr-Schild family)
        KerrSchildParams params;
        params.M = m_Config.blackHoleMass;
        params.a = m_Config.blackHoleSpin * m_Config.blackHoleMass;  // spin as a/M
        params.Q = 0.0;  // No charge for now
        params.Lambda = 0.0;  // No cosmological constant
        m_Metric = std::make_unique<KerrSchildFamily>(params);
        std::cout << "  Metric:     " << m_Metric->getName() << " (M=" << params.M
                  << ", a=" << params.a << ")" << std::endl;

        // Create camera
        CameraConfig camConfig;
        camConfig.r = m_Config.observerDistance;
        camConfig.theta = m_Config.observerInclination;
        camConfig.phi = 0.0;  // Azimuthal angle
        camConfig.fov = m_Config.cameraFOV;
        camConfig.width = m_Config.width;
        camConfig.height = m_Config.height;
        m_Camera = std::make_unique<PinholeCamera>(camConfig);
        std::cout << "  Observer:   r=" << camConfig.r << "M, θ="
                  << (camConfig.theta * 180.0 / Math::PI) << "°" << std::endl;

        // Create geodesic tracer
        TracerConfig tracerConfig;
        tracerConfig.escape_radius = 200.0f;  // Increased for better escape detection
        tracerConfig.horizon_factor = 1.05f;
        tracerConfig.max_steps = 20000;  // Increased for long geodesics near photon sphere

        // Optimize integrator for faster rendering while maintaining visual accuracy
        // Balance: large steps far from BH, small steps near horizon
        tracerConfig.integrator.initial_step = 0.1f;   // 10x faster start (was 0.01)
        tracerConfig.integrator.max_step = 2.0f;       // Allow larger steps far from BH (was 0.1)
        tracerConfig.integrator.min_step = 1e-5f;      // Precise near horizon
        tracerConfig.integrator.abs_tolerance = 5e-6f; // Moderate precision
        tracerConfig.integrator.rel_tolerance = 5e-6f;

        // Compute ISCO for disk inner radius
        auto rISCO = AccretionDiskD::computeISCO(m_Config.blackHoleSpin);
        tracerConfig.disk_inner = static_cast<float>(rISCO * m_Config.blackHoleMass);
        tracerConfig.disk_outer = static_cast<float>(20.0 * m_Config.blackHoleMass);
        tracerConfig.enable_disk = true;

        // Phase 6: Volumetric disk configuration
        tracerConfig.enable_volumetric = m_Config.enableVolumetricDisk;
        tracerConfig.volumetric_H_over_r = m_Config.volumetricHOverR;
        tracerConfig.volumetric_H_power = m_Config.volumetricHPower;
        tracerConfig.volumetric_tau_midplane = m_Config.volumetricTauMidplane;
        tracerConfig.volumetric_samples = m_Config.volumetricSamples;

        m_Tracer = std::make_unique<GeodesicTracer>(m_Metric.get(), tracerConfig);
        std::cout << "  Disk:       r_in=" << tracerConfig.disk_inner << "M, r_out="
                  << tracerConfig.disk_outer << "M" << std::endl;
        std::cout << "[Session] Physics engine initialized (geodesic integration enabled)" << std::endl;

        // =================================================================
        // Phase 7: Initialize Relativistic Jets
        // =================================================================
        if (m_Config.enableJets) {
            JetConfig jetConfig;
            jetConfig.lorentz_factor = m_Config.jetLorentzFactor;
            jetConfig.opening_angle = m_Config.jetOpeningAngle;
            jetConfig.r_launch = m_Config.jetLaunchRadius;
            jetConfig.r_max = m_Config.jetMaxExtent;
            jetConfig.collimation = m_Config.jetCollimation;
            jetConfig.spectral_index = m_Config.jetSpectralIndex;
            m_Jet = std::make_unique<RelativisticJet>(jetConfig);
            std::cout << "[Session] Relativistic jets enabled: Γ=" << jetConfig.lorentz_factor
                      << ", θ_open=" << (jetConfig.opening_angle * 180.0 / Math::PI) << "°" << std::endl;
        }

        // =================================================================
        // Phase 7: Initialize Polarisation Buffer
        // =================================================================
        if (m_Config.enablePolarisation) {
            m_PolarisationBuffer.resize(m_Config.width * m_Config.height);
            std::cout << "[Session] Polarisation tracking enabled" << std::endl;
        }

        // =================================================================
        // Color Mode
        // =================================================================
        const char* modeName = "TrueColor";
        switch (m_Config.colorMode) {
            case ColorModes::Mode::TrueColor: modeName = "TrueColor (Physical)"; break;
            case ColorModes::Mode::TemperatureMap: modeName = "TemperatureMap (False Color)"; break;
            case ColorModes::Mode::RedshiftMap: modeName = "RedshiftMap (g-factor)"; break;
            case ColorModes::Mode::Narrowband: modeName = "Narrowband (Hubble Palette)"; break;
            case ColorModes::Mode::Polarisation: modeName = "Polarisation"; break;
        }
        std::cout << "[Session] Color mode: " << modeName << std::endl;

        // =================================================================
        // Multi-threaded Rendering Setup
        // =================================================================
        if (m_Config.enableParallelRendering) {
            m_NumThreads = m_Config.threadCount;
            if (m_NumThreads <= 0) {
                // Auto-detect: use hardware concurrency, leave 1 core for system
                m_NumThreads = static_cast<int>(std::thread::hardware_concurrency());
                if (m_NumThreads > 1) m_NumThreads--;
                if (m_NumThreads < 1) m_NumThreads = 1;
            }

            // Create per-thread tracers (they share the metric but each has its own state)
            m_ThreadTracers.clear();
            m_ThreadTracers.reserve(m_NumThreads);
            for (int i = 0; i < m_NumThreads; ++i) {
                m_ThreadTracers.push_back(std::make_unique<GeodesicTracer>(m_Metric.get(), tracerConfig));
            }

            std::cout << "[Session] Multi-threaded rendering enabled: " << m_NumThreads << " threads" << std::endl;
        } else {
            m_NumThreads = 1;
            std::cout << "[Session] Single-threaded rendering" << std::endl;
        }

        // =================================================================
        // Load Starfield Background Texture (must happen BEFORE GPU init)
        // =================================================================
        // Try multiple paths to find the starfield texture
        std::vector<std::string> texturePaths = {
            "Sirius.Render/Texture/Starfield.png",           // Relative to binary
            "../Sirius.Render/Texture/Starfield.png",        // One level up
            "src/Sirius.Render/Texture/Starfield.png",       // From project root
            "../src/Sirius.Render/Texture/Starfield.png",    // From build dir
            "../../src/Sirius.Render/Texture/Starfield.png"  // From bin/bin
        };

        bool textureLoaded = false;
        for (const auto& path : texturePaths) {
            if (loadStarfieldTexture(path)) {
                textureLoaded = true;
                break;
            }
        }

        if (!textureLoaded) {
            std::cerr << "[Session] Warning: Could not load starfield texture from any path" << std::endl;
            std::cerr << "[Session] Background will render as faint grey" << std::endl;
        }

        // =================================================================
        // GPU Acceleration (OptiX)
        // =================================================================
        m_UseGPU = false;
        if (m_Config.useGPU) {
            try {
                auto bestBackend = Acceleration::getBestBackend();
                if (bestBackend != Acceleration::BackendType::None) {
                    m_Accelerator = Acceleration::createAccelerator(bestBackend);
                    if (m_Accelerator && m_Accelerator->initialise(m_Config.width, m_Config.height)) {
                        // Upload starfield texture to GPU (now loaded above)
                        if (m_StarfieldLoaded) {
                            if (m_Accelerator->uploadBackground(m_StarfieldData.data(),
                                                                m_StarfieldWidth, m_StarfieldHeight)) {
                                std::cout << "[Session] Uploaded starfield texture to GPU" << std::endl;
                            } else {
                                std::cerr << "[Session] Warning: Failed to upload starfield to GPU" << std::endl;
                            }
                        }
                        m_UseGPU = true;
                        std::cout << "[Session] GPU acceleration enabled: "
                                  << Acceleration::backendName(bestBackend)
                                  << " (" << m_Accelerator->getCapabilities().deviceName << ")" << std::endl;
                    } else {
                        std::cout << "[Session] GPU initialization failed, falling back to CPU" << std::endl;
                        if (m_Accelerator) {
                            std::cout << "[Session] Error: " << m_Accelerator->getLastError() << std::endl;
                        }
                        m_Accelerator.reset();
                    }
                } else {
                    std::cout << "[Session] No GPU backend available, using CPU" << std::endl;
                }
            } catch (const std::exception& e) {
                std::cout << "[Session] GPU initialization exception: " << e.what() << std::endl;
                std::cout << "[Session] Falling back to CPU rendering" << std::endl;
                m_Accelerator.reset();
            }
        } else {
            std::cout << "[Session] GPU disabled by config, using CPU rendering" << std::endl;
        }

        // Transition to Scheduling
        m_FSM.process(SessionEvent::Ready);
    }
    catch (const std::exception& e) {
        m_ErrorMessage = e.what();
        m_FSM.process(SessionEvent::Error);
    }
}

//==============================================================================
// Tile Scheduling
//==============================================================================
void RenderSession::scheduleNextTile() {
    // Check for cancellation
    if (m_Progress.getCancellationToken().isCancelled()) {
        m_FSM.process(SessionEvent::Cancel);
        return;
    }

    // GPU rendering (single launch for entire frame)
    if (m_UseGPU && m_Accelerator && m_Accelerator->isInitialised()) {
        renderGPU();
        return;
    }

    // Use parallel rendering if enabled and multiple threads available
    if (m_Config.enableParallelRendering && m_NumThreads > 1) {
        // Parallel mode: launch all workers and process all tiles
        renderTilesParallel();
        return;
    }

    // Single-threaded mode: sequential tile processing
    Tile* tile = m_Tiles.getNextTile();

    if (tile == nullptr) {
        // All tiles complete
        m_FSM.process(SessionEvent::AllTilesComplete);
        return;
    }

    // Signal tile available and render
    m_FSM.process(SessionEvent::TileAvailable);

    // Render the tile synchronously
    renderTile(tile);
}

//==============================================================================
// Pixel Shading Helpers (unified for single-threaded and multi-threaded)
//==============================================================================

RenderSession::PixelResult RenderSession::shadeDiskHit(const TraceResult& result) const {
    PixelResult px;

    // Volumetric disk: use pre-integrated emission from ray marching
    if (result.volumetric_hit) {
        float vol_intensity = result.volumetric_emission[0];
        float T_effective = std::pow(vol_intensity, 0.25f);
        float T_kelvin = std::clamp(T_effective * T_INNER_KELVIN, 1000.0f, 100000.0f);

        Spectral::RGB bbColor = Spectral::blackbodyToRGB(static_cast<double>(T_kelvin));
        float mag = result.magnification;
        px.r = bbColor.r * vol_intensity * mag;
        px.g = bbColor.g * vol_intensity * mag;
        px.b = bbColor.b * vol_intensity * mag;

        if (result.photon_ring) {
            px.r *= PHOTON_RING_BOOST_VOLUMETRIC;
            px.g *= PHOTON_RING_BOOST_VOLUMETRIC;
            px.b *= PHOTON_RING_BOOST_VOLUMETRIC;
        }
        return px;
    }

    // Thin disk: accumulate emission from all disk crossings
    float total_r = 0.0f, total_g = 0.0f, total_b = 0.0f;

    for (int crossing_idx = 0; crossing_idx < result.num_disk_crossings; crossing_idx++) {
        const auto& crossing = result.disk_crossings[crossing_idx];
        if (!crossing.valid) continue;

        float T_emit = crossing.temperature;
        float g = crossing.redshift;
        float r_cross = crossing.r;

        // Higher-order demagnification: ~exp(-nπ)
        float order_demag = std::exp(-static_cast<float>(Math::PI) * crossing_idx);
        float T_obs = T_emit * g;
        float intensity = std::pow(T_obs, 4.0f);

        Spectral::RGB diskColor = ColorModes::applyColorMode(
            m_Config.colorMode, T_emit, g, intensity, nullptr);

        float cr = diskColor.r;
        float cg = diskColor.g;
        float cb = diskColor.b;

        // Limb darkening (primary crossing only)
        if (crossing_idx == 0) {
            float A = result.gfactor_A;
            float B = result.gfactor_B;
            float n_xy_squared = A * A + B * B;
            float cos_theta = std::sqrt(std::max(0.0f, 1.0f - n_xy_squared));

            Spectral::RGB limbInput(cr, cg, cb);
            Spectral::RGB darkened = Spectral::applyLimbDarkening(limbInput, cos_theta);
            cr = darkened.r;
            cg = darkened.g;
            cb = darkened.b;
        }

        // Motion blur (primary crossing only)
        float g_blur = g;
        if (crossing_idx == 0 && m_Config.enableMotionBlur && m_Config.motionBlurSamples > 1) {
            float grav = result.gfactor_grav;
            float gamma = result.gfactor_gamma;
            float v_orb = result.gfactor_v_orb;
            float A = result.gfactor_A;
            float B = result.gfactor_B;

            float M = static_cast<float>(m_Config.blackHoleMass);
            float a = static_cast<float>(m_Config.blackHoleSpin * M);
            float sqrtM = std::sqrt(M);
            float Omega = sqrtM / (std::pow(r_cross, 1.5f) + a * sqrtM);

            float delta_phi_max = Omega * m_Config.shutterTime;
            int N = m_Config.motionBlurSamples;
            float g_sum = 0.0f;

            for (int i = 0; i < N; i++) {
                float t = (N > 1) ? static_cast<float>(i) / (N - 1) : 0.5f;
                float delta_phi = delta_phi_max * (t - 0.5f);

                float cos_dphi = std::cos(delta_phi);
                float sin_dphi = std::sin(delta_phi);
                float v_dot_n_offset = v_orb * (A * cos_dphi + B * sin_dphi);
                float doppler_denom = gamma * (1.0f - v_dot_n_offset);
                doppler_denom = std::clamp(doppler_denom, DOPPLER_CLAMP_MIN, DOPPLER_CLAMP_MAX);

                float g_offset = grav / doppler_denom;
                g_offset = std::clamp(g_offset, G_FACTOR_CLAMP_MIN, G_FACTOR_CLAMP_MAX);
                g_sum += g_offset;
            }

            g_blur = g_sum / N;
        }

        // Relativistic beaming: I_obs = g^4 × I_emit
        float g4 = g_blur * g_blur * g_blur * g_blur;
        cr *= g4;
        cg *= g4;
        cb *= g4;

        total_r += cr * order_demag;
        total_g += cg * order_demag;
        total_b += cb * order_demag;
    }

    // Apply gravitational lensing magnification and cinematic boost
    px.r = total_r * result.magnification * DISK_INTENSITY_BOOST;
    px.g = total_g * result.magnification * DISK_INTENSITY_BOOST;
    px.b = total_b * result.magnification * DISK_INTENSITY_BOOST;

    if (result.photon_ring) {
        px.r *= PHOTON_RING_BOOST_DISK;
        px.g *= PHOTON_RING_BOOST_DISK;
        px.b *= PHOTON_RING_BOOST_DISK;
    }

    return px;
}

RenderSession::PixelResult RenderSession::shadeEscaped(const TraceResult& result) const {
    PixelResult px;
    sampleStarfield(result.final_direction, px.r, px.g, px.b);

    if (result.magnification > 1.0f) {
        px.r *= result.magnification;
        px.g *= result.magnification;
        px.b *= result.magnification;
    }

    if (result.photon_ring) {
        px.r *= PHOTON_RING_BOOST_ESCAPED;
        px.g *= PHOTON_RING_BOOST_ESCAPED;
        px.b *= PHOTON_RING_BOOST_ESCAPED;
    }

    // Relativistic jet emission
    if (m_Jet && m_Config.enableJets) {
        double obs_r = m_Config.observerDistance;
        double obs_theta = m_Config.observerInclination;
        float obs_x = static_cast<float>(obs_r * std::sin(obs_theta));
        float obs_y = 0.0f;
        float obs_z = static_cast<float>(obs_r * std::cos(obs_theta));

        float end_x = static_cast<float>(result.final_position(1));
        float end_y = static_cast<float>(result.final_position(2));
        float end_z = static_cast<float>(result.final_position(3));

        float jet_emission = JetRayMarching::integrateJetEmission(
            *m_Jet, obs_x, obs_y, obs_z,
            end_x, end_y, end_z,
            obs_x, obs_y, obs_z, 64);

        if (jet_emission > 0.0f) {
            float jet_scale = m_Config.jetIntensity * JET_NORMALISATION;
            px.r += jet_emission * jet_scale * 0.8f;
            px.g += jet_emission * jet_scale * 0.9f;
            px.b += jet_emission * jet_scale * 1.0f;
        }
    }

    return px;
}

RenderSession::PixelResult RenderSession::shadePixel(int px_coord, int py_coord, GeodesicTracer* tracer) const {
    PixelResult result;
    float r_acc = 0.0f, g_acc = 0.0f, b_acc = 0.0f;
    StokesVector stokes_acc;

    int spp = std::max(1, m_Config.samplesPerPixel);
    int grid_size = static_cast<int>(std::sqrt(static_cast<float>(spp)));
    if (grid_size < 1) grid_size = 1;

    for (int sy = 0; sy < grid_size; ++sy) {
        for (int sx = 0; sx < grid_size; ++sx) {
            float u = (sx + 0.5f) / grid_size;
            float v = (sy + 0.5f) / grid_size;

            CameraRay camRay = m_Camera->generateRay(px_coord, py_coord, u, v);
            TraceResult traceResult = tracer->trace(camRay);

            float sr = 0.0f, sg = 0.0f, sb = 0.0f;

            switch (traceResult.outcome) {
                case TraceResult::Outcome::HORIZON:
                    break;

                case TraceResult::Outcome::DISK_HIT: {
                    PixelResult disk = shadeDiskHit(traceResult);
                    sr = disk.r;
                    sg = disk.g;
                    sb = disk.b;
                    break;
                }

                case TraceResult::Outcome::ESCAPED: {
                    PixelResult esc = shadeEscaped(traceResult);
                    sr = esc.r;
                    sg = esc.g;
                    sb = esc.b;
                    break;
                }

                case TraceResult::Outcome::SPIRALING:
                    sr = sg = sb = SPIRALING_BRIGHTNESS;
                    break;

                case TraceResult::Outcome::MAX_STEPS:
                default:
                    sr = sg = sb = MAX_STEPS_BRIGHTNESS;
                    break;
            }

            r_acc += sr;
            g_acc += sg;
            b_acc += sb;

            // Polarisation accumulation
            if (m_Config.enablePolarisation && traceResult.outcome == TraceResult::Outcome::DISK_HIT) {
                float disk_phi = traceResult.disk_phi;
                float evpa = disk_phi + static_cast<float>(Math::HALF_PI);
                float I_sample = std::sqrt(sr*sr + sg*sg + sb*sb);
                StokesVector sample_stokes = PolarisedEmission::synchrotronEmission(
                    I_sample, SYNCHROTRON_POL_DEGREE, evpa);
                stokes_acc += sample_stokes;
            }
        }
    }

    int total_samples = grid_size * grid_size;
    float inv_samples = 1.0f / static_cast<float>(total_samples);
    result.r = r_acc * inv_samples;
    result.g = g_acc * inv_samples;
    result.b = b_acc * inv_samples;

    if (m_Config.enablePolarisation) {
        stokes_acc *= inv_samples;
        stokes_acc.normalise();
        result.stokes = stokes_acc;
    }

    return result;
}

//==============================================================================
// Tile Rendering
//==============================================================================
void RenderSession::renderTile(Tile* tile) {
    if (!tile) return;

    std::vector<float> tileBuffer(tile->width * tile->height * 4, 0.0f);

    for (int ty = 0; ty < tile->height; ++ty) {
        for (int tx = 0; tx < tile->width; ++tx) {
            int px = tile->x + tx;
            int py = tile->y + ty;

            PixelResult pixel = shadePixel(px, py, m_Tracer.get());

            int idx = (ty * tile->width + tx) * 4;
            tileBuffer[idx + 0] = pixel.r;
            tileBuffer[idx + 1] = pixel.g;
            tileBuffer[idx + 2] = pixel.b;
            tileBuffer[idx + 3] = 1.0f;

            if (m_Config.enablePolarisation && !m_PolarisationBuffer.empty()) {
                int pixel_idx = py * m_Config.width + px;
                if (pixel_idx >= 0 && pixel_idx < static_cast<int>(m_PolarisationBuffer.size())) {
                    m_PolarisationBuffer[pixel_idx] = pixel.stokes;
                }
            }
        }
    }

    m_Display.updateTile(tile->x, tile->y, tile->width, tile->height, tileBuffer.data());

    m_Tiles.completeTile(tile->id);
    m_Progress.completeTile(tile->pixelCount());

    int percent = static_cast<int>(m_Progress.getProgress() * 100);
    std::cout << "\r[Session] Progress: " << percent << "% | ETA: "
              << m_Progress.getETAString() << "    " << std::flush;

    m_FSM.process(SessionEvent::TileComplete);
}

//==============================================================================
// Starfield Texture Loading
//==============================================================================
bool RenderSession::loadStarfieldTexture(const std::string& path) {
    int channels;
    auto deleter = [](unsigned char* p) { if (p) stbi_image_free(p); };
    std::unique_ptr<unsigned char[], decltype(deleter)> data(
        stbi_load(path.c_str(), &m_StarfieldWidth, &m_StarfieldHeight, &channels, 4),
        deleter);

    if (!data) {
        std::cerr << "[Session] Failed to load starfield texture: " << path << std::endl;
        std::cerr << "[Session] STB error: " << stbi_failure_reason() << std::endl;
        return false;
    }

    size_t size = static_cast<size_t>(m_StarfieldWidth * m_StarfieldHeight * 4);
    m_StarfieldData.assign(data.get(), data.get() + size);

    m_StarfieldLoaded = true;
    std::cout << "[Session] Loaded starfield texture: " << m_StarfieldWidth << "x"
              << m_StarfieldHeight << " RGBA (" << (size / 1024 / 1024) << " MB)" << std::endl;

    return true;
}

//==============================================================================
// Starfield Background Sampling (Equirectangular Projection)
//==============================================================================
void RenderSession::sampleStarfield(const Vec4& direction, float& r, float& g, float& b) const {
    // Default to black if texture not loaded
    if (!m_StarfieldLoaded || m_StarfieldData.empty()) {
        r = g = b = BACKGROUND_FALLBACK;
        return;
    }

    // Convert direction to spherical angles for equirectangular lookup
    // direction is in Cartesian (x, y, z) from geodesic tracer
    double dx = direction(1);
    double dy = direction(2);
    double dz = direction(3);

    // Normalize
    double len = std::sqrt(dx*dx + dy*dy + dz*dz);
    if (len < 1e-10) {
        r = g = b = 0.001f;
        return;
    }

    dx /= len;
    dy /= len;
    dz /= len;

    // Convert to spherical coordinates (theta = polar, phi = azimuthal)
    // theta: 0 at +Z (north pole) to PI at -Z (south pole)
    // phi: 0 at +X, increasing toward +Y
    double theta = std::acos(std::clamp(dz, -1.0, 1.0));
    double phi = std::atan2(dy, dx);
    if (phi < 0) phi += Math::TWO_PI;

    // Map to UV coordinates for equirectangular projection
    // u = phi / (2*PI), range [0, 1]
    // v = theta / PI, range [0, 1]
    double u = phi / Math::TWO_PI;
    double v = theta / Math::PI;

    // Convert to pixel coordinates with bilinear sampling
    double px = u * (m_StarfieldWidth - 1);
    double py = v * (m_StarfieldHeight - 1);

    int x0 = static_cast<int>(std::floor(px));
    int y0 = static_cast<int>(std::floor(py));
    int x1 = std::min(x0 + 1, m_StarfieldWidth - 1);
    int y1 = std::min(y0 + 1, m_StarfieldHeight - 1);

    double fx = px - x0;
    double fy = py - y0;

    // Sample 4 neighboring pixels (RGBA format, 4 bytes per pixel)
    auto sample = [this](int x, int y) -> std::array<float, 3> {
        int idx = (y * m_StarfieldWidth + x) * 4;  // 4 bytes per pixel (RGBA)
        return {
            m_StarfieldData[idx + 0] / 255.0f,
            m_StarfieldData[idx + 1] / 255.0f,
            m_StarfieldData[idx + 2] / 255.0f
            // idx + 3 is alpha, which we ignore
        };
    };

    auto c00 = sample(x0, y0);
    auto c10 = sample(x1, y0);
    auto c01 = sample(x0, y1);
    auto c11 = sample(x1, y1);

    // Bilinear interpolation
    float w00 = static_cast<float>((1.0 - fx) * (1.0 - fy));
    float w10 = static_cast<float>(fx * (1.0 - fy));
    float w01 = static_cast<float>((1.0 - fx) * fy);
    float w11 = static_cast<float>(fx * fy);

    r = c00[0] * w00 + c10[0] * w10 + c01[0] * w01 + c11[0] * w11;
    g = c00[1] * w00 + c10[1] * w10 + c01[1] * w01 + c11[1] * w11;
    b = c00[2] * w00 + c10[2] * w10 + c01[2] * w01 + c11[2] * w11;
}

//==============================================================================
// Output Writing
//==============================================================================
void RenderSession::writeOutput() {
    try {
        std::cout << "\n[Session] Writing output..." << std::endl;

        // Apply post-processing (cinematic pipeline)
        PostProcessConfig ppConfig;
        ppConfig.tonemapper = TonemapType::ACES;       // Filmic ACES curve
        ppConfig.exposure = m_Config.exposure;          // Default: 3.0
        ppConfig.gamma = 2.2f;

        // Bloom for accretion disk glow
        ppConfig.enableBloom = m_Config.enableBloom;
        ppConfig.bloomIntensity = m_Config.bloomIntensity;  // Default: 0.5
        ppConfig.bloomThreshold = m_Config.bloomThreshold;  // Default: 0.3
        ppConfig.bloomRadius = BLOOM_RADIUS;

        // Color grading for cinematic look
        ppConfig.saturation = m_Config.saturation;      // Default: 1.15
        ppConfig.contrast = m_Config.contrast;          // Default: 1.1
        ppConfig.lift = SHADOW_LIFT;
        ppConfig.gain = 1.0f;

        PostProcessor::process(m_Display.getFloatBuffer(), m_Config.width, m_Config.height, ppConfig);

        // Write main PPM output
        std::ofstream file(m_Config.outputPath, std::ios::binary);
        if (!file) {
            throw std::runtime_error("Failed to open output file: " + m_Config.outputPath);
        }

        file << "P6\n" << m_Config.width << " " << m_Config.height << "\n255\n";

        const float* data = m_Display.getFloatData();
        for (int y = 0; y < m_Config.height; ++y) {
            for (int x = 0; x < m_Config.width; ++x) {
                int idx = (y * m_Config.width + x) * 4;
                unsigned char rgb[3] = {
                    static_cast<unsigned char>(std::clamp(data[idx + 0], 0.0f, 1.0f) * 255.0f),
                    static_cast<unsigned char>(std::clamp(data[idx + 1], 0.0f, 1.0f) * 255.0f),
                    static_cast<unsigned char>(std::clamp(data[idx + 2], 0.0f, 1.0f) * 255.0f)
                };
                file.write(reinterpret_cast<char*>(rgb), 3);
            }
        }

        std::cout << "[Session] Wrote: " << m_Config.outputPath << std::endl;

        // =================================================================
        // Phase 7: Write Polarisation Map Output
        // =================================================================
        // Outputs a separate image showing polarisation degree and angle.
        // Uses HSV colormap: Hue = EVPA, Saturation = pol degree, Value = intensity
        // =================================================================
        if (m_Config.outputPolarisationMap && !m_PolarisationBuffer.empty()) {
            std::string polPath = m_Config.outputPath;
            size_t dotPos = polPath.rfind('.');
            if (dotPos != std::string::npos) {
                polPath.insert(dotPos, "_polarisation");
            } else {
                polPath += "_polarisation.ppm";
            }

            std::ofstream polFile(polPath, std::ios::binary);
            if (polFile) {
                polFile << "P6\n" << m_Config.width << " " << m_Config.height << "\n255\n";

                for (int y = 0; y < m_Config.height; ++y) {
                    for (int x = 0; x < m_Config.width; ++x) {
                        int idx = y * m_Config.width + x;
                        const StokesVector& stokes = m_PolarisationBuffer[idx];

                        // Convert Stokes to RGB using HSV scheme
                        Spectral::RGB color = ColorModes::PolarisationVis::stokesToRGB_HSV(stokes);

                        unsigned char rgb[3] = {
                            static_cast<unsigned char>(std::clamp(color.r, 0.0f, 1.0f) * 255.0f),
                            static_cast<unsigned char>(std::clamp(color.g, 0.0f, 1.0f) * 255.0f),
                            static_cast<unsigned char>(std::clamp(color.b, 0.0f, 1.0f) * 255.0f)
                        };
                        polFile.write(reinterpret_cast<char*>(rgb), 3);
                    }
                }

                std::cout << "[Session] Wrote polarisation map: " << polPath << std::endl;
            }
        }

        m_FSM.process(SessionEvent::OutputWritten);
    }
    catch (const std::exception& e) {
        m_ErrorMessage = e.what();
        m_FSM.process(SessionEvent::Error);
    }
}

//==============================================================================
// Session End
//==============================================================================
void RenderSession::onSessionEnd(SessionState state) {
    std::cout << "\n[Session] Finished with state: " << stateName(state) << std::endl;
    
    std::string message;
    switch (state) {
        case SessionState::Complete:
            message = "Render completed successfully";
            break;
        case SessionState::Failed:
            message = "Render failed: " + m_ErrorMessage;
            break;
        case SessionState::Cancelled:
            message = "Render cancelled by user";
            break;
        default:
            message = "Unknown end state";
    }
    
    if (m_CompletionCallback) {
        m_CompletionCallback(state, message);
    }

    // Clean up worker threads if any
    m_StopWorkers = true;
    for (auto& thread : m_WorkerThreads) {
        if (thread.joinable()) {
            thread.join();
        }
    }
    m_WorkerThreads.clear();
}

//==============================================================================
// Multi-threaded Tile Rendering
//==============================================================================
void RenderSession::renderTilesParallel() {
    m_StopWorkers = false;
    m_ActiveWorkers = 0;

    // Launch worker threads
    m_WorkerThreads.clear();
    m_WorkerThreads.reserve(m_NumThreads);

    std::cout << "[Session] Starting " << m_NumThreads << " worker threads..." << std::endl;

    for (int i = 0; i < m_NumThreads; ++i) {
        m_WorkerThreads.emplace_back(&RenderSession::workerThread, this, i);
    }

    // Wait for all workers to complete
    for (auto& thread : m_WorkerThreads) {
        if (thread.joinable()) {
            thread.join();
        }
    }
    m_WorkerThreads.clear();

    std::cout << "\n[Session] All worker threads completed" << std::endl;

    // Transition to completion
    m_FSM.process(SessionEvent::AllTilesComplete);
}

void RenderSession::workerThread(int threadId) {
    m_ActiveWorkers++;

    while (!m_StopWorkers) {
        // Get next tile (thread-safe)
        Tile* tile = nullptr;
        {
            std::lock_guard<std::mutex> lock(m_TileMutex);
            tile = m_Tiles.getNextTile();
        }

        if (tile == nullptr) {
            // No more tiles - exit worker
            break;
        }

        // Check for cancellation
        if (m_Progress.getCancellationToken().isCancelled()) {
            m_StopWorkers = true;
            break;
        }

        // Render the tile using per-thread tracer
        renderTileThreaded(tile, threadId);

        // Update progress (thread-safe)
        {
            std::lock_guard<std::mutex> lock(m_TileMutex);
            m_Tiles.completeTile(tile->id);
            m_Progress.completeTile(tile->pixelCount());

            // Print progress (from any thread)
            int percent = static_cast<int>(m_Progress.getProgress() * 100);
            std::cout << "\r[Session] Progress: " << percent << "% | ETA: "
                      << m_Progress.getETAString() << " | Threads: " << m_ActiveWorkers.load()
                      << "    " << std::flush;
        }
    }

    m_ActiveWorkers--;
}

//==============================================================================
// GPU Rendering (OptiX Backend)
//==============================================================================
void RenderSession::renderGPU() {
    if (!m_Accelerator || !m_Accelerator->isInitialised()) {
        std::cerr << "[Session] GPU not initialized" << std::endl;
        m_FSM.process(SessionEvent::Error);
        return;
    }

    std::cout << "[Session] Starting GPU render..." << std::endl;

    // Configure launch parameters
    Acceleration::LaunchConfig config;
    config.width = m_Config.width;
    config.height = m_Config.height;
    config.samplesPerPixel = m_Config.samplesPerPixel;
    config.tileSize = m_Config.tileSize;
    config.blackHoleMass = static_cast<float>(m_Config.blackHoleMass);
    config.blackHoleSpin = static_cast<float>(m_Config.blackHoleSpin);
    config.observerDistance = static_cast<float>(m_Config.observerDistance);
    config.observerInclination = static_cast<float>(m_Config.observerInclination);
    config.cameraFOV = m_Config.cameraFOV;

    // Compute ISCO for disk inner radius
    double M = m_Config.blackHoleMass;
    auto rISCO = AccretionDiskD::computeISCO(m_Config.blackHoleSpin);
    config.diskInnerRadius = static_cast<float>(rISCO * M);
    config.diskOuterRadius = static_cast<float>(20.0 * M);
    config.diskTemperature = 30000.0f;

    // Integration parameters
    config.maxSteps = 20000;
    config.maxStepSize = 2.0f;

    // Metric type: 2 = Kerr
    config.metricType = (std::abs(m_Config.blackHoleSpin) < 0.01) ? 1 : 2;  // 1=Schwarzschild, 2=Kerr
    config.metricFamily = 0;  // Kerr-Schild family

    // Post-processing / exposure
    // For realistic dark space scenes, exposure < 1.0 darkens the scene
    config.exposure = m_Config.exposure;
    config.contrast = m_Config.contrast;
    config.saturation = m_Config.saturation;
    config.bloomIntensity = m_Config.bloomIntensity;

    // Reset accumulation for fresh render
    m_Accelerator->resetAccumulation();

    // Launch GPU render
    auto start_time = std::chrono::high_resolution_clock::now();
    m_Accelerator->launch(config);
    m_Accelerator->synchronise();
    auto end_time = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    std::cout << "[Session] GPU render completed in " << duration.count() << "ms" << std::endl;

    // Copy frame buffer to display
    float* gpuBuffer = m_Accelerator->getFrameBuffer();
    if (gpuBuffer) {
        // Copy GPU buffer to display buffer
        std::vector<float>& displayBuffer = m_Display.getFloatBuffer();
        size_t bufferSize = m_Config.width * m_Config.height * 4;
        std::memcpy(displayBuffer.data(), gpuBuffer, bufferSize * sizeof(float));

        std::cout << "[Session] Frame buffer copied to display" << std::endl;
    } else {
        std::cerr << "[Session] GPU frame buffer is null" << std::endl;
    }

    // Mark all tiles complete
    m_Progress.setTotals(1, 1);  // Treat as single "tile"
    m_Progress.completeTile(m_Config.width * m_Config.height);

    // Transition to completion
    m_FSM.process(SessionEvent::AllTilesComplete);
}

void RenderSession::renderTileThreaded(Tile* tile, int threadId) {
    if (!tile) return;

    GeodesicTracer* tracer = nullptr;
    if (threadId >= 0 && threadId < static_cast<int>(m_ThreadTracers.size())) {
        tracer = m_ThreadTracers[threadId].get();
    } else {
        tracer = m_Tracer.get();
    }

    std::vector<float> tileBuffer(tile->width * tile->height * 4, 0.0f);

    for (int ty = 0; ty < tile->height; ++ty) {
        if (ty % 8 == 0 && m_StopWorkers) {
            return;
        }

        for (int tx = 0; tx < tile->width; ++tx) {
            int px = tile->x + tx;
            int py = tile->y + ty;

            PixelResult pixel = shadePixel(px, py, tracer);

            int idx = (ty * tile->width + tx) * 4;
            tileBuffer[idx + 0] = pixel.r;
            tileBuffer[idx + 1] = pixel.g;
            tileBuffer[idx + 2] = pixel.b;
            tileBuffer[idx + 3] = 1.0f;

            if (m_Config.enablePolarisation && !m_PolarisationBuffer.empty()) {
                int pixel_idx = py * m_Config.width + px;
                if (pixel_idx >= 0 && pixel_idx < static_cast<int>(m_PolarisationBuffer.size())) {
                    std::lock_guard<std::mutex> lock(m_DisplayMutex);
                    m_PolarisationBuffer[pixel_idx] = pixel.stokes;
                }
            }
        }
    }

    {
        std::lock_guard<std::mutex> lock(m_DisplayMutex);
        m_Display.updateTile(tile->x, tile->y, tile->width, tile->height, tileBuffer.data());
    }
}

//==============================================================================
// Configuration Conversion
//==============================================================================
SessionConfig SessionConfig::fromSiriusConfig(const Configuration::SiriusConfig& config) {
    SessionConfig sc;
    sc.width = config.render.width;
    sc.height = config.render.height;
    sc.samplesPerPixel = config.render.samplesPerPixel;
    sc.tileSize = config.render.tileSize;
    sc.outputPath = config.render.outputPath;
    sc.metricName = config.metric.name;
    sc.blackHoleMass = config.metric.mass;
    sc.blackHoleSpin = config.metric.spin;
    sc.observerDistance = config.observer.distance;
    sc.observerInclination = config.observer.inclination * Math::PI / 180.0;
    sc.cameraFOV = static_cast<float>(config.observer.fov);

    // Post-processing
    sc.enableBloom = config.postprocess.enableBloom;
    sc.bloomIntensity = config.postprocess.bloomIntensity;
    sc.bloomThreshold = config.postprocess.bloomThreshold;
    sc.exposure = config.postprocess.exposure;
    sc.contrast = config.postprocess.contrast;
    sc.saturation = config.postprocess.saturation;

    // Volumetric disk
    sc.enableVolumetricDisk = config.volumetric.enabled;
    sc.volumetricHOverR = config.volumetric.hOverR;
    sc.volumetricHPower = config.volumetric.hPower;
    sc.volumetricTauMidplane = config.volumetric.tauMidplane;
    sc.volumetricSamples = config.volumetric.samples;

    // Advanced volumetric (turbulence + corona)
    if (config.volumetric.enableTurbulence || config.volumetric.enableCorona) {
        sc.enableVolumetricDiskAdvanced = true;
        sc.volumetricDiskConfig.turbulence.enabled = config.volumetric.enableTurbulence;
        sc.volumetricDiskConfig.corona.enabled = config.volumetric.enableCorona;
        sc.volumetricDiskConfig.H_over_r = config.volumetric.hOverR;
        sc.volumetricDiskConfig.H_power = config.volumetric.hPower;
        sc.volumetricDiskConfig.volumetric_samples = config.volumetric.samples;
    }

    // Film simulation
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

    // Backend
    sc.useGPU = (config.backend.preferred != "cpu");

    return sc;
}

} // namespace Sirius
