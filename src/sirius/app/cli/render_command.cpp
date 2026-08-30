// The render command. Ported from CRCL002A.cpp.

#include "sirius/app/cli/render_command.h"

#include "sirius/app/cli/cli_output.h"
#include "sirius/app/config/config_loader.h"
#include "sirius/app/config/session_config_adapter.h"
#include "sirius/render/session/render_session.h"

#ifdef SIRIUS_HAS_VULKAN_BACKEND
#include "sirius/backend/device.h"
#endif

#include <nlohmann/json.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace sirius::app {

namespace cli = cli_output;

namespace {

bool IsGpuBackendName(std::string name) {
    std::transform(name.begin(), name.end(), name.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return name == "gpu" || name == "vulkan";
}

int ParseInteger(const std::string& text) {
    std::size_t consumed = 0;
    const int value = std::stoi(text, &consumed);
    if (consumed != text.size()) {
        throw std::invalid_argument("trailing characters");
    }
    return value;
}

double ParseDouble(const std::string& text) {
    std::size_t consumed = 0;
    const double value = std::stod(text, &consumed);
    if (consumed != text.size() || !std::isfinite(value)) {
        throw std::invalid_argument("expected one finite number");
    }
    return value;
}

float ParseFloat(const std::string& text) {
    const double value = ParseDouble(text);
    if (value < -std::numeric_limits<float>::max() || value > std::numeric_limits<float>::max()) {
        throw std::out_of_range("outside the float range");
    }
    return static_cast<float>(value);
}

void ApplyCompatibleLensDefaults(ObserverConfig& observer) {
    const auto lens = core::ParseLensType(observer.lens_model);
    if (!lens.has_value() || *lens == core::LensType::ThinLens) return;
    observer.focal_length = core::kDefaultCameraFocalLength;
    observer.aperture = core::kDefaultCameraAperture;
    observer.focus_distance = core::kDefaultCameraFocusDistance;
}

}  // namespace

std::string RenderCommand::Usage() const {
    return R"(Usage: sirius render [options]

Render a black hole visualization using geodesic ray tracing.

Basic Options:
  -o, --output <path>       Output file path (default: render.ppm)
  -w, --width <n>           Image width (default: 1920)
  -h, --height <n>          Image height (default: 1080)
  -s, --samples <n>         Samples per pixel (default: 64)
  -t, --tile-size <n>       Tile size in pixels (default: 64)
  -m, --metric <name>       Metric (see 'sirius info metrics'): Minkowski, Schwarzschild,
                            Kerr, Reissner-Nordstrom, Kerr-Newman, de-Sitter,
                            Schwarzschild-de-Sitter, Morris-Thorne, Alcubierre
  --mass <length>           Metric mass M in geometric coordinate units (0 when absent;
                            otherwise 0.1-100)
  -d, --distance <r>        Observer coordinate radius (5M-1000M for mass metrics;
                            5L-1000L for governed exotic metrics; default: 50)
  -i, --inclination <deg>   Observer inclination (default: 90)
  -a, --spin <a/M>          Dimensionless black-hole spin 0-0.998 (default: 0)
  --fov <deg>               Camera field of view (default: 60)
  --lens <name>             Lens: Pinhole, ThinLens, or Fisheye (fisheye CPU-only;
                            default: Pinhole)
  --focal-length <mm-eq>    Thin-lens focal length; selects ThinLens (50 mm-eq = 1 lens unit)
  --aperture <f-number>     Thin-lens aperture; selects ThinLens (default: 2.8)
  --focus-distance <length> Thin-lens focus distance in geometric coordinate units
                            and selects ThinLens (default: 50)
  --temperature-model <m>   Disk temperature model: NovikovThorne (NT) or
                            ShakuraSunyaev (SS) (default: NovikovThorne)
  --disk-temperature <T>    Disk temperature at 1.5 times the inner edge in Kelvin
                            (default: 50000)
  --color-mode <name>       TrueColor, TemperatureMap, RedshiftMap, or
                            Polarisation (Narrowband is reserved and declines)
  --throat-radius <b0>      Morris-Thorne-only throat radius (default: 1.0)
  --wormhole-topology <t>   OneSheetCapture (supported) or TwoSheet (declines explicitly)
  --warp-velocity <vs>      Alcubierre-only warp velocity (default: 0.5)
  --bubble-radius <R>       Alcubierre-only bubble radius (default: 1.0)
  --bubble-sigma <1/length> Alcubierre inverse wall scale (default: 0.5;
                            0.1 <= sigma*R <= 100)

Post-Processing:
  --exposure <e>            Exposure value (default: 1.0)
  --bloom <intensity>       Bloom intensity 0-1 (default: 0.3)
  --bloom-threshold <t>     Bloom brightness threshold; enables bloom (default: 0.3)
  --contrast <c>            Contrast adjustment (default: 1.0)
  --saturation <s>          Saturation adjustment (default: 1.0)
  --tonemapper <name>       Tonemapper: ACES, Reinhard, Filmic, Uncharted2,
                            None, Linear (default: ACES)
  --no-bloom                Disable bloom post-processing

Volumetric Disk (3D accretion disk with thickness):
  --volumetric              Enable volumetric disk rendering
  --h-over-r <ratio>        Inner-edge H/r in [0.01,0.5]; enables volume (default: 0.1)
  --h-power <exp>           Flaring exponent; enables volume (default: 0.25)
  --tau <depth>             Inner-edge vertical optical depth; enables volume (default: 10.0)
  --vol-samples <n>         Ray marching samples; enables volume (default: 64)
  --turbulence              Enable deterministic procedural density perturbations (not GRMHD)
  --corona                  Decline: spectral covariant Compton transfer is not represented
  --no-disk                 Disable accretion-disk emission (required for charged/Λ metrics)

Temporal Disk:
  --motion-blur             Enable CPU temporal integration of thin-disk rotation
  --shutter-time <t>        Exposure interval in geometric time units (default: 0.1)
  --motion-samples <n>      Temporal samples (default: 3)

Film-Inspired Display Finish:
  --film                    Enable the bounded display finish
  --film-preset <name>      Preset: Interstellar, SpaceOdyssey2001
                            (default: Interstellar)
  --grain <intensity>       Override preset grain intensity; enables finish
  --halation <strength>     Override preset halation strength; enables finish
  --vignette <strength>     Override preset vignette strength; enables finish

DNGR Physics (default off; the pinned render is unchanged):
  --beams                   Propagate ray bundles (geodesic deviation) on the live path (P2)
  --starfield <mode>        Star field: 'point' (beam-filtered point sources, P3) or
                            'texture' (default equirectangular lookup)
  --doppler-beaming <on|off> Disk Doppler asymmetry: 'on' (full physics, default) or
                            'off' (suppressed, the Interstellar look, P4)
  --camera-beta <f[,u,r]>   Camera four-velocity beta as forward[,up,right], |beta|<1 (P5)

Presets:
  --cinematic               Enable cinematic defaults (volumetric, bloom, high exposure)

Backend:
  --cpu                     Force CPU rendering (the reference path)
  --no-gpu                  Alias of --cpu
  --gpu                     Render on the Vulkan backend (declines if no device)
  --backend <name>          Select backend: auto (Vulkan only when the entire scene is
                            representable), cpu, vulkan

Examples:
  sirius render -o test.ppm -s 32
  sirius render -m Kerr -a 0.9 -s 256 -o kerr.ppm
  sirius render -w 3840 -h 2160 --cinematic --volumetric -o cinematic.ppm
  sirius render --volumetric --h-over-r 0.15 -o volumetric.ppm
)";
}

int RenderCommand::Execute(const std::vector<std::string>& args, const GlobalOptions& globals,
                           SiriusConfig& config) {
    // Command objects are reusable through the router/test seam. Backend intent
    // belongs to this invocation and must never leak from an earlier --gpu run.
    gpu_backend_requested_ = IsGpuBackendName(config.backend.preferred);
    if (!ParseArgs(args, globals, config)) {
        return 1;
    }

    // Backend selection. An explicit request (--gpu / --backend vulkan, or a
    // lower-layer SIRIUS_BACKEND=vulkan) is resolved here, before validation,
    // so a backend decline (not a config-name error) is what the user sees when
    // no device is visible or the backend is not compiled in. 'auto' is resolved
    // by MakeSessionConfig (the single authority): Vulkan when a
    // device and a registry-gpu-dispatchable metric align, CPU otherwise
    // (go-live, owner decision 2026-07-18).
    const bool use_vulkan = gpu_backend_requested_;

    if (use_vulkan) {
#ifdef SIRIUS_HAS_VULKAN_BACKEND
        auto devices = backend::EnumerateVulkanDevices();
        if (!devices.has_value()) {
            cli::Error("Vulkan device enumeration failed: " + devices.error().Description());
            return 1;
        }
        if (devices->empty()) {
            cli::Error("Vulkan backend requested but no Vulkan device is present; use --cpu");
            return 1;
        }
        if (auto selected = backend::ResolveVulkanDeviceIndex(*devices); !selected) {
            cli::Error(selected.error().Description());
            return 1;
        }
#else
        cli::Error("Vulkan backend requested but this build has no Vulkan support; use --cpu");
        return 1;
#endif
    }

    auto errors = ConfigLoader::Validate(config);
    if (!errors.empty()) {
        for (const auto& err : errors) {
            cli::Error(err);
        }
        return 1;
    }

    if (!globals.json_output) {
        PrintConfig(config, globals.verbose);
    }

    return ExecuteSession(config, globals, use_vulkan);
}

bool RenderCommand::ParseArgs(const std::vector<std::string>& args,
                              const GlobalOptions& /*globals*/, SiriusConfig& config) {
    bool metric_overridden = false;
    bool mass_overridden = false;
    for (std::size_t i = 0; i < args.size(); ++i) {
        const std::string& arg = args[i];

        try {
            if ((arg == "-o" || arg == "--output") && i + 1 < args.size()) {
                config.render.output_path = args[++i];
            } else if ((arg == "-w" || arg == "--width") && i + 1 < args.size()) {
                config.render.width = ParseInteger(args[++i]);
            } else if ((arg == "-h" || arg == "--height") && i + 1 < args.size()) {
                config.render.height = ParseInteger(args[++i]);
            } else if ((arg == "-s" || arg == "--samples") && i + 1 < args.size()) {
                config.render.samples_per_pixel = ParseInteger(args[++i]);
            } else if ((arg == "-t" || arg == "--tile-size") && i + 1 < args.size()) {
                config.render.tile_size = ParseInteger(args[++i]);
            } else if ((arg == "-m" || arg == "--metric") && i + 1 < args.size()) {
                config.metric.name = args[++i];
                metric_overridden = true;
            } else if (arg == "--mass" && i + 1 < args.size()) {
                config.metric.mass = ParseDouble(args[++i]);
                mass_overridden = true;
            } else if ((arg == "-d" || arg == "--distance") && i + 1 < args.size()) {
                config.observer.distance = ParseDouble(args[++i]);
            } else if ((arg == "-i" || arg == "--inclination") && i + 1 < args.size()) {
                config.observer.inclination = ParseDouble(args[++i]);
            } else if ((arg == "-a" || arg == "--spin") && i + 1 < args.size()) {
                config.metric.spin = ParseDouble(args[++i]);
            } else if (arg == "--fov" && i + 1 < args.size()) {
                config.observer.fov = ParseDouble(args[++i]);
            } else if (arg == "--lens" && i + 1 < args.size()) {
                config.observer.lens_model = args[++i];
                ApplyCompatibleLensDefaults(config.observer);
            } else if (arg == "--focal-length" && i + 1 < args.size()) {
                config.observer.focal_length = ParseFloat(args[++i]);
                config.observer.lens_model = "ThinLens";
            } else if (arg == "--aperture" && i + 1 < args.size()) {
                config.observer.aperture = ParseFloat(args[++i]);
                config.observer.lens_model = "ThinLens";
            } else if (arg == "--focus-distance" && i + 1 < args.size()) {
                config.observer.focus_distance = ParseFloat(args[++i]);
                config.observer.lens_model = "ThinLens";
            } else if (arg == "--temperature-model" && i + 1 < args.size()) {
                config.metric.temperature_model = args[++i];
                config.disk_enabled = true;
            } else if (arg == "--disk-temperature" && i + 1 < args.size()) {
                config.metric.disk_temperature = ParseFloat(args[++i]);
                config.disk_enabled = true;
            } else if (arg == "--color-mode" && i + 1 < args.size()) {
                config.color_mode = args[++i];
                if (config.color_mode != "TrueColor") config.disk_enabled = true;
            } else if (arg == "--throat-radius" && i + 1 < args.size()) {
                config.metric.throat_radius = ParseDouble(args[++i]);
            } else if (arg == "--wormhole-topology" && i + 1 < args.size()) {
                const std::string topology = args[++i];
                if (topology == "one-sheet" || topology == "OneSheetCapture") {
                    config.metric.wormhole_topology = "OneSheetCapture";
                } else if (topology == "two-sheet" || topology == "TwoSheet") {
                    config.metric.wormhole_topology = "TwoSheet";
                } else {
                    cli::Error(
                        "--wormhole-topology expects OneSheetCapture/one-sheet or "
                        "TwoSheet/two-sheet");
                    return false;
                }
            } else if (arg == "--warp-velocity" && i + 1 < args.size()) {
                config.metric.warp_velocity = ParseDouble(args[++i]);
            } else if (arg == "--bubble-radius" && i + 1 < args.size()) {
                config.metric.bubble_radius = ParseDouble(args[++i]);
            } else if (arg == "--bubble-sigma" && i + 1 < args.size()) {
                config.metric.bubble_sigma = ParseDouble(args[++i]);
            } else if (arg == "--exposure" && i + 1 < args.size()) {
                config.postprocess.exposure = ParseFloat(args[++i]);
            } else if (arg == "--bloom" && i + 1 < args.size()) {
                config.postprocess.enable_bloom = true;
                config.postprocess.bloom_intensity = ParseFloat(args[++i]);
            } else if (arg == "--bloom-threshold" && i + 1 < args.size()) {
                config.postprocess.bloom_threshold = ParseFloat(args[++i]);
                config.postprocess.enable_bloom = true;
            } else if (arg == "--contrast" && i + 1 < args.size()) {
                config.postprocess.contrast = ParseFloat(args[++i]);
            } else if (arg == "--saturation" && i + 1 < args.size()) {
                config.postprocess.saturation = ParseFloat(args[++i]);
            } else if (arg == "--tonemapper" && i + 1 < args.size()) {
                config.postprocess.tonemapper = args[++i];
            } else if (arg == "--no-bloom") {
                config.postprocess.enable_bloom = false;
            } else if (arg == "--volumetric") {
                config.volumetric.enabled = true;
                config.disk_enabled = true;
            } else if (arg == "--h-over-r" && i + 1 < args.size()) {
                config.volumetric.h_over_r = ParseFloat(args[++i]);
                config.volumetric.enabled = true;
                config.disk_enabled = true;
            } else if (arg == "--h-power" && i + 1 < args.size()) {
                config.volumetric.h_power = ParseFloat(args[++i]);
                config.volumetric.enabled = true;
                config.disk_enabled = true;
            } else if (arg == "--tau" && i + 1 < args.size()) {
                config.volumetric.tau_midplane = ParseFloat(args[++i]);
                config.volumetric.enabled = true;
                config.disk_enabled = true;
            } else if (arg == "--vol-samples" && i + 1 < args.size()) {
                config.volumetric.samples = ParseInteger(args[++i]);
                config.volumetric.enabled = true;
                config.disk_enabled = true;
            } else if (arg == "--turbulence") {
                config.volumetric.enabled = true;
                config.volumetric.enable_turbulence = true;
                config.disk_enabled = true;
            } else if (arg == "--corona") {
                config.volumetric.enabled = true;
                config.volumetric.enable_corona = true;
                config.disk_enabled = true;
            } else if (arg == "--no-disk") {
                config.disk_enabled = false;
            } else if (arg == "--motion-blur") {
                config.motion_blur.enabled = true;
                config.disk_enabled = true;
            } else if (arg == "--shutter-time" && i + 1 < args.size()) {
                config.motion_blur.shutter_time = ParseFloat(args[++i]);
                config.motion_blur.enabled = true;
                config.disk_enabled = true;
            } else if (arg == "--motion-samples" && i + 1 < args.size()) {
                config.motion_blur.samples = ParseInteger(args[++i]);
                config.motion_blur.enabled = true;
                config.disk_enabled = true;
            } else if (arg == "--film") {
                config.film.enabled = true;
            } else if (arg == "--film-preset" && i + 1 < args.size()) {
                config.film.preset = args[++i];
                config.film.enabled = true;
            } else if (arg == "--grain" && i + 1 < args.size()) {
                config.film.grain_intensity = ParseFloat(args[++i]);
                config.film.enabled = true;
            } else if (arg == "--halation" && i + 1 < args.size()) {
                config.film.halation_strength = ParseFloat(args[++i]);
                config.film.enabled = true;
            } else if (arg == "--vignette" && i + 1 < args.size()) {
                config.film.vignette_strength = ParseFloat(args[++i]);
                config.film.enabled = true;
            } else if (arg == "--cinematic") {
                config.volumetric.enabled = true;
                config.disk_enabled = true;
                config.volumetric.h_over_r = 0.12f;
                config.volumetric.samples = 64;
                config.postprocess.enable_bloom = true;
                config.postprocess.bloom_intensity = 0.35f;
                config.postprocess.bloom_threshold = 0.4f;
                config.postprocess.exposure = 1.2f;
                config.postprocess.contrast = 1.1f;
                config.postprocess.saturation = 1.15f;
            } else if (arg == "--beams" || arg == "--ray-bundles") {
                config.ray_bundles = true;  // P2: propagate geodesic deviation.
            } else if (arg == "--starfield" && i + 1 < args.size()) {
                std::string mode = args[++i];
                std::transform(mode.begin(), mode.end(), mode.begin(),
                               [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
                if (mode == "point") {
                    config.point_starfield = true;  // P3: filtered point sources.
                } else if (mode == "texture") {
                    config.point_starfield = false;
                } else {
                    cli::Error("--starfield expects 'point' or 'texture'");
                    return false;
                }
            } else if (arg == "--doppler-beaming" && i + 1 < args.size()) {
                std::string mode = args[++i];
                std::transform(mode.begin(), mode.end(), mode.begin(),
                               [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
                if (mode == "off" || mode == "false" || mode == "0") {
                    config.doppler_beaming = false;  // P4: suppress the asymmetry.
                } else if (mode == "on" || mode == "true" || mode == "1") {
                    config.doppler_beaming = true;
                } else {
                    cli::Error("--doppler-beaming expects 'on' or 'off'");
                    return false;
                }
                config.disk_enabled = true;
            } else if (arg == "--camera-beta" && i + 1 < args.size()) {
                // P5: camera four-velocity beta as forward[,up,right] in [0, 1).
                const std::string spec = args[++i];
                double comp[3] = {0.0, 0.0, 0.0};
                std::size_t start = 0;
                int components = 0;
                while (start <= spec.size()) {
                    if (components == 3) {
                        cli::Error("--camera-beta accepts at most three components");
                        return false;
                    }
                    std::size_t comma = spec.find(',', start);
                    std::string tok = spec.substr(start, comma - start);
                    if (tok.empty()) {
                        cli::Error("--camera-beta components must not be empty");
                        return false;
                    }
                    comp[components++] = ParseDouble(tok);
                    if (comma == std::string::npos) break;
                    start = comma + 1;
                }
                double mag = std::sqrt(comp[0] * comp[0] + comp[1] * comp[1] + comp[2] * comp[2]);
                if (mag >= 1.0) {
                    cli::Error("--camera-beta magnitude must be < 1 (sub-luminal worldline)");
                    return false;
                }
                config.observer.camera_beta_forward = comp[0];
                config.observer.camera_beta_up = comp[1];
                config.observer.camera_beta_right = comp[2];
            } else if (arg == "--no-gpu" || arg == "--cpu") {
                config.backend.preferred = "cpu";
                gpu_backend_requested_ = false;
            } else if (arg == "--gpu") {
                // Keep the schema and the explicit execution override aligned.
                // Execute performs the early device probe; the typed session
                // must still see "vulkan" so it cannot emit an auto-fallback
                // diagnostic before ExecuteSession pins the Vulkan path.
                config.backend.preferred = "vulkan";
                gpu_backend_requested_ = true;
            } else if (arg == "--backend" && i + 1 < args.size()) {
                const std::string& name = args[++i];
                if (IsGpuBackendName(name)) {
                    config.backend.preferred = "vulkan";
                    gpu_backend_requested_ = true;
                } else {
                    config.backend.preferred = name;  // cpu or auto; validated later.
                    gpu_backend_requested_ = false;
                }
            } else if (arg.substr(0, 1) == "-") {
                cli::Error("Unknown option: " + arg);
                return false;
            } else {
                cli::Error("Unexpected positional argument: " + arg);
                return false;
            }
        } catch (const std::exception& e) {
            cli::Error("Invalid value for " + arg + ": " + e.what());
            return false;
        }
    }

    if (metric_overridden && !mass_overridden) {
        const auto metric = core::ParseMetricName(config.metric.name);
        if (metric.has_value()) {
            if (!core::MetricUsesMass(*metric)) {
                config.metric.mass = 0.0;
            } else if (config.metric.mass == 0.0) {
                config.metric.mass = MetricConfig{}.mass;
            }
        }
    }
    return true;
}

void RenderCommand::PrintConfig(const SiriusConfig& config, bool verbose) {
    cli::Rule("Sirius Render");
    std::cout << std::endl;

    std::cout << "Configuration:" << std::endl;
    std::cout << "  Resolution:  " << config.render.width << " x " << config.render.height
              << std::endl;
    std::cout << "  Samples:     " << config.render.samples_per_pixel << " spp" << std::endl;
    std::cout << "  Tile size:   " << config.render.tile_size << " px" << std::endl;
    std::cout << "  Metric:      " << config.metric.name;
    if (config.metric.spin > 0) {
        std::cout << " (a=" << std::fixed << std::setprecision(3) << config.metric.spin << ")";
    }
    std::cout << std::endl;
    const auto metric_id = core::ParseMetricName(config.metric.name);
    const bool uses_mass = metric_id.has_value() && core::MetricUsesMass(*metric_id);
    if (uses_mass) {
        std::cout << "  Metric mass: " << config.metric.mass << " coordinate units" << std::endl;
    }
    std::cout << "  Observer:    r=" << std::fixed << std::setprecision(1)
              << config.observer.distance << " coordinate units";
    if (uses_mass) {
        std::cout << " (r/M=" << config.observer.distance / config.metric.mass << ")";
    }
    std::cout << ", θ=" << config.observer.inclination << "°" << std::endl;
    std::cout << "  FOV:         " << config.observer.fov << "°" << std::endl;
    std::cout << "  Color mode:  " << config.color_mode << std::endl;
    std::cout << "  Output:      " << config.render.output_path << std::endl;
    std::cout << std::endl;

    if (config.volumetric.enabled) {
        std::cout << "Volumetric Disk:" << std::endl;
        std::cout << "  Scale height: H/r=" << std::fixed << std::setprecision(2)
                  << config.volumetric.h_over_r << std::endl;
        std::cout << "  Flaring:      " << config.volumetric.h_power << std::endl;
        std::cout << "  Optical depth:" << config.volumetric.tau_midplane << std::endl;
        std::cout << "  Ray samples:  " << config.volumetric.samples << std::endl;
        if (config.volumetric.enable_turbulence)
            std::cout << "  Turbulence:   enabled" << std::endl;
        if (config.volumetric.enable_corona) std::cout << "  Corona:       enabled" << std::endl;
        std::cout << std::endl;
    }

    if (config.film.enabled) {
        const auto finish = ProjectFilmFinishConfig(config.film);
        if (!finish) return;
        std::cout << "Film-Inspired Display Finish:" << std::endl;
        std::cout << "  Preset:      " << config.film.preset << std::endl;
        std::cout << "  Grain:       " << std::fixed << std::setprecision(2)
                  << finish->grain_intensity << std::endl;
        std::cout << "  Halation:    " << finish->halation_strength << std::endl;
        std::cout << "  Vignette:    " << finish->vignette_strength << std::endl;
        std::cout << std::endl;
    }

    if (verbose || config.postprocess.exposure != 1.0f) {
        std::cout << "Post-processing:" << std::endl;
        std::cout << "  Bloom:       "
                  << (config.postprocess.enable_bloom ? "enabled" : "disabled");
        if (config.postprocess.enable_bloom) {
            std::cout << " (intensity=" << std::fixed << std::setprecision(2)
                      << config.postprocess.bloom_intensity << ")";
        }
        std::cout << std::endl;
        std::cout << "  Exposure:    " << std::fixed << std::setprecision(2)
                  << config.postprocess.exposure << std::endl;
        if (config.postprocess.contrast != 1.0f || config.postprocess.saturation != 1.0f) {
            std::cout << "  Contrast:    " << config.postprocess.contrast << std::endl;
            std::cout << "  Saturation:  " << config.postprocess.saturation << std::endl;
        }
        std::cout << std::endl;
    }
}

int RenderCommand::ExecuteSession(const SiriusConfig& config, const GlobalOptions& globals,
                                  bool use_vulkan) {
    auto adapted = MakeSessionConfig(config);
    if (!adapted) {
        cli::Error(adapted.error().Description());
        return 1;
    }
    auto session_config = std::move(*adapted);
    if (use_vulkan) {
        session_config.backend = render::RenderBackend::Vulkan;
    }

    render::RenderSession session;
    if (auto configured = session.Configure(session_config); !configured) {
        cli::Error(configured.error().Description());
        return 1;
    }

    auto start_time = std::chrono::steady_clock::now();

    if (!globals.json_output) {
        session.SetProgressCallback([&](float progress, int completed, int total, double eta) {
            auto now = std::chrono::steady_clock::now();
            double elapsed = std::chrono::duration<double>(now - start_time).count();

            cli::ProgressState state;
            state.progress = progress;
            state.tiles_completed = completed;
            state.tiles_total = total;
            state.elapsed_seconds = elapsed;
            state.eta_seconds = eta;

            cli::PrintProgress(state);
        });
    }

    // The final state carries success; only the failure message needs capturing.
    std::string error_message;
    session.SetCompletionCallback([&](render::SessionState state, const std::string& message) {
        if (!globals.json_output) {
            cli::ClearProgress();
        }
        if (state != render::SessionState::Complete) {
            error_message = message;
        }
    });

    render::SessionState result = session.Execute();

    if (globals.json_output) {
        nlohmann::json j;
        j["success"] = (result == render::SessionState::Complete);
        j["output"] = config.render.output_path;
        j["state"] = std::string(render::StateName(result));
        if (!error_message.empty()) {
            j["error"] = error_message;
        }
        cli::PrintJson(j.dump(2));
    } else {
        std::cout << std::endl;
        cli::Rule();
        if (result == render::SessionState::Complete) {
            cli::Success("Render complete: " + config.render.output_path);
        } else if (result == render::SessionState::Cancelled) {
            cli::Warning("Render cancelled");
        } else {
            cli::Error("Render failed: " + error_message);
        }
    }

    return (result == render::SessionState::Complete) ? 0 : 1;
}

}  // namespace sirius::app
