// Vulkan render path implementation. See vulkan_renderer.h for the contract.
//
// Parameter mapping (Task A.3): the camera frame is built in world Cartesian to
// reproduce the CPU pinhole geometry exactly, so both backends trace the same
// scene. The CPU tracer places the observer at BlToCartesian(r, inclination,
// phi=0) and maps image +x -> +phi, image +y (top) -> -theta; the kernel's NDC
// runs v top-to-bottom, so camUp is +theta_hat. The two integrators differ by
// design (RK45 Hamiltonian on the CPU, Cartesian RK4 on the kernel), so parity
// is statistical, not bitwise (docs/ARCHITECTURE.md section 3).

#include "sirius/render/vulkan_renderer.h"

#include "sirius/backend/device.h"
#include "sirius/base/resource_locator.h"
#include "sirius/core/constants.h"
#include "sirius/core/coordinates.h"
#include "sirius/core/disk/novikov_thorne_disk.h"
#include "sirius/core/metrics/registry.h"
#include "sirius/core/starfield.h"
#include "sirius/render/dispatch_governor.h"
#include "sirius/render/pixel_sampling.h"
#include "sirius/render/session/display_buffer.h"
#include "sirius/render/session/render_session.h"
#include "sirius/render/trace_domain.h"

#include "stb_image.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <format>
#include <fstream>
#include <iostream>
#include <optional>
#include <span>
#include <string>
#include <vector>

namespace sirius::render {

namespace {

using base::ErrorDomain;
using base::Expected;
using base::Fail;
namespace math = core::constants::math;

// --- Trace-kernel integration constants --------------------------------------
// Termination and step control for the Vulkan render dispatch. escape_radius and
// capture_factor mirrors the CPU session's exact Kerr-Schild capture surface.
// so the shadow and background geometry coincide; the step schedule is the
// kernel's own (the integrator differs from the CPU RK45 by design).
constexpr float kTraceStepScale = 0.08f;

// Thin-disk extent in units of M. The outer edge matches the CPU session's
// disk_outer = 20 M and the inner edge is the ISCO. Emission and transfer are
// owned by the shared Page-Thorne/blackbody/invariant-transfer kernel path.
constexpr float kDiskOuterFactor = 20.0f;

// Exact dense params-buffer ABI through the last consumed slot (64).
// Host and kernel indices are kept explicit so a new control must extend both.
constexpr std::uint32_t kParamCount = 65;
constexpr int kMaxVulkanVolumeSamples = 128;

// Kernel dispatch ids (must match trace.slang / gr_types.slang).
constexpr float kDispatchKerrSchild = 0.0f;
constexpr float kDispatchWarpDrive = 2.0f;
constexpr float kDispatchMorrisThorne = 3.0f;

// What the trace kernel needs to render one metric, or a decline.
struct KernelScene {
    float metric_id = kDispatchKerrSchild;
    float M = 0.0f;
    float a = 0.0f;  // absolute spin a = (a/M) * M
    float Q = 0.0f;
    float Lambda = 0.0f;
    float warp_vs = 0.0f;
    float warp_sigma = 0.0f;
    float warp_R = 0.0f;
    float throat_b0 = 0.0f;   // Morris-Thorne throat radius
    float worm_shape = 0.0f;  // 0 Ellis, 1 zero-tidal, 2 absurdly benign
    bool disk_enabled = false;
    std::string metric_name;
};

// Maps a session metric onto the trace kernel, or declines. The gpu-renderable
// set is exactly the registry's gpu_supported metrics the Cartesian render path
// carries: Minkowski/Schwarzschild/Kerr (Kerr-Schild family), Alcubierre, and
// Morris-Thorne through its Cartesian embedding (one sheet, throat capture,
// mirroring the CPU authority). Charge/lambda metrics (registry
// gpu_supported=false) decline loudly.
[[nodiscard]] Expected<KernelScene> MapMetric(const SessionConfig& config) {
    KernelScene scene;
    const auto info = core::MetricInfoFor(config.metric_id);
    scene.metric_name = info.canonical_name;

    switch (config.metric_id) {
        case core::MetricId::Minkowski:
            scene.metric_id = kDispatchKerrSchild;
            scene.M = 0.0f;
            scene.disk_enabled = false;
            return scene;
        case core::MetricId::Schwarzschild:
            scene.metric_id = kDispatchKerrSchild;
            scene.M = static_cast<float>(config.black_hole_mass);
            scene.disk_enabled = true;
            return scene;
        case core::MetricId::Kerr:
            scene.metric_id = kDispatchKerrSchild;
            scene.M = static_cast<float>(config.black_hole_mass);
            scene.a = static_cast<float>(config.black_hole_spin * config.black_hole_mass);
            scene.disk_enabled = true;
            return scene;
        case core::MetricId::Alcubierre:
            scene.metric_id = kDispatchWarpDrive;
            scene.warp_vs = static_cast<float>(config.warp_velocity);
            scene.warp_sigma = static_cast<float>(config.bubble_sigma);
            scene.warp_R = static_cast<float>(config.bubble_radius);
            scene.disk_enabled = false;
            return scene;
        case core::MetricId::MorrisThorne:
            scene.metric_id = kDispatchMorrisThorne;
            scene.throat_b0 = static_cast<float>(config.throat_radius);
            scene.worm_shape = 0.0f;  // Ellis; the session constructs Ellis(b0).
            scene.disk_enabled = false;
            return scene;
        case core::MetricId::ReissnerNordstrom:
        case core::MetricId::KerrNewman:
        case core::MetricId::DeSitter:
        case core::MetricId::SchwarzschildDeSitter:
            return Fail(ErrorDomain::kDevice, "select Vulkan render metric",
                        std::format("'{}' carries charge or a cosmological constant that is not "
                                    "plumbed to the Vulkan render dispatch (registry "
                                    "gpu_supported=false); use --cpu",
                                    info.canonical_name));
    }
    return Fail(ErrorDomain::kDevice, "select Vulkan render metric", "unknown metric id");
}

// Selects the precision-ladder rung (docs/ARCHITECTURE.md section 6).
// SIRIUS_PRECISION=fp64 selects the double-precision trace kernel
// (trace_fp64.spv); it requires the device to report shaderFloat64 and
// declines loudly otherwise, never silently running fp32 (docs/STYLE.md
// section 4). SIRIUS_PRECISION=fp32-comp selects the compensated rung
// (trace_fp32comp.spv, Kahan state accumulation), which any device runs.
// Unset keeps the plain fp32 rung, the default: fp64 costs multiples of fp32
// throughput on consumer GPUs and the fp32 path passes the parity gate. Any
// other value is a loud error, not a silent default.
[[nodiscard]] Expected<PrecisionRung> SelectPrecisionRung(bool device_supports_fp64) {
    const char* precision = std::getenv("SIRIUS_PRECISION");
    if (precision == nullptr || *precision == '\0') {
        return PrecisionRung::Fp32;
    }
    const std::string requested(precision);
    if (requested == "fp32") {
        return PrecisionRung::Fp32;
    }
    if (requested == "fp32-comp") {
        return PrecisionRung::Fp32Comp;
    }
    if (requested == "fp64") {
        if (!device_supports_fp64) {
            return Fail(ErrorDomain::kDevice, "select precision rung",
                        "SIRIUS_PRECISION=fp64 requested but the device lacks shaderFloat64");
        }
        return PrecisionRung::Fp64;
    }
    return Fail(ErrorDomain::kDevice, "select precision rung",
                "unknown SIRIUS_PRECISION '" + requested + "' (accepted: fp32, fp32-comp, fp64)");
}

[[nodiscard]] const char* RungName(PrecisionRung rung) {
    switch (rung) {
        case PrecisionRung::Fp32:
            return "fp32";
        case PrecisionRung::Fp32Comp:
            return "fp32-comp";
        case PrecisionRung::Fp64:
            return "fp64";
    }
    return "fp32";
}

[[nodiscard]] const char* RungKernelName(PrecisionRung rung) {
    switch (rung) {
        case PrecisionRung::Fp32:
            return "trace.spv";
        case PrecisionRung::Fp32Comp:
            return "trace_fp32comp.spv";
        case PrecisionRung::Fp64:
            return "trace_fp64.spv";
    }
    return "trace.spv";
}

[[nodiscard]] std::vector<std::uint32_t> LoadSpirv(PrecisionRung rung) {
    const std::string kernel_name = RungKernelName(rung);
    const auto path = base::ResolveResource("kernels/" + kernel_name);
    if (!path) {
        return {};
    }
    std::ifstream file(*path, std::ios::binary | std::ios::ate);
    if (!file) return {};
    const auto size = static_cast<std::size_t>(file.tellg());
    if (size == 0 || size % sizeof(std::uint32_t) != 0) return {};
    std::vector<std::uint32_t> words(size / sizeof(std::uint32_t));
    file.seekg(0);
    file.read(reinterpret_cast<char*>(words.data()), static_cast<std::streamsize>(size));
    return words;
}

// Loads and packs the equirectangular starfield to one uint32 (RGBA8) per texel,
// matching the CPU path's byte layout. Returns an empty vector when the required
// asset is missing or unreadable; the caller declines the request.
[[nodiscard]] std::vector<std::uint32_t> LoadStarfield(int& out_width, int& out_height) {
    const auto path = base::ResolveResource("assets/Starfield.png");
    if (!path) {
        out_width = 0;
        out_height = 0;
        return {};
    }
    int channels = 0;
    unsigned char* pixels =
        stbi_load(path->string().c_str(), &out_width, &out_height, &channels, 4);
    if (pixels == nullptr) return {};
    const auto texels = static_cast<std::size_t>(out_width) * out_height;
    std::vector<std::uint32_t> packed(texels);
    for (std::size_t i = 0; i < texels; ++i) {
        packed[i] = static_cast<std::uint32_t>(pixels[i * 4 + 0]) |
                    (static_cast<std::uint32_t>(pixels[i * 4 + 1]) << 8) |
                    (static_cast<std::uint32_t>(pixels[i * 4 + 2]) << 16) |
                    (static_cast<std::uint32_t>(pixels[i * 4 + 3]) << 24);
    }
    stbi_image_free(pixels);
    return packed;
}

// Fills the scene-invariant slots of the params buffer (everything except the
// per-tile origin/size in [32..35], written per dispatch).
void FillSceneParams(std::vector<float>& params, const SessionConfig& config,
                     const KernelScene& scene, bool starfield_enabled, int starfield_width,
                     int starfield_height) {
    const double theta = config.observer_inclination;
    const double st = std::sin(theta);
    const double ct = std::cos(theta);
    const double phi = config.observer_azimuth;
    const double sp = std::sin(phi);
    const double cp = std::cos(phi);
    const double r = config.observer_distance;

    params[0] = static_cast<float>(config.width);
    params[1] = static_cast<float>(config.height);
    params[2] = scene.metric_id;
    params[3] = scene.M;
    params[4] = scene.a;
    params[5] = scene.Q;
    params[6] = scene.Lambda;

    // Match GeodesicTracer::InitializeLightray exactly. The observer position
    // is the oblate Boyer-Lindquist -> Kerr-Schild Cartesian map, not a
    // spherical substitution. Camera directions are local orthonormal
    // components rotated about the actual Cartesian azimuth of that position.
    const double x = (r * cp - scene.a * sp) * st;
    const double y = (r * sp + scene.a * cp) * st;
    const double cartesian_phi = std::atan2(y, x);
    const double cart_sp = std::sin(cartesian_phi);
    const double cart_cp = std::cos(cartesian_phi);
    params[7] = static_cast<float>(x);
    params[8] = static_cast<float>(y);
    params[9] = static_cast<float>(r * ct);
    params[10] = static_cast<float>(-st * cart_cp);
    params[11] = static_cast<float>(-st * cart_sp);
    params[12] = static_cast<float>(-ct);
    params[13] = static_cast<float>(-cart_sp);
    params[14] = static_cast<float>(cart_cp);
    params[15] = 0.0f;
    params[16] = static_cast<float>(ct * cart_cp);
    params[17] = static_cast<float>(ct * cart_sp);
    params[18] = static_cast<float>(-st);

    params[19] = static_cast<float>(config.camera_fov * math::kPi / 180.0);  // full-angle radians
    params[20] = static_cast<float>(config.width) / static_cast<float>(config.height);

    const TraceDomainParameters trace_domain = BuildTraceDomainParameters({
        .metric_id = config.metric_id,
        .metric_mass = config.black_hole_mass,
        .cosmological_constant = config.cosmological_constant,
        .observer_radius = config.observer_distance,
        .throat_radius = config.throat_radius,
        .bubble_radius = config.bubble_radius,
        .bubble_sigma = config.bubble_sigma,
    });
    params[21] = static_cast<float>(kRenderTraceMaximumAttempts);
    params[22] = kTraceStepScale;
    params[23] = trace_domain.vulkan_min_step;
    params[24] = trace_domain.max_step;
    params[25] = trace_domain.escape_radius;
    params[26] = kRenderTraceCaptureFactor;

    const double spin = (config.black_hole_mass > 0.0) ? config.black_hole_spin : 0.0;
    const double isco = core::AccretionDiskD::ComputeIsco(spin);
    params[27] = scene.disk_enabled ? 1.0f : 0.0f;
    params[28] = static_cast<float>(isco * config.black_hole_mass);
    params[29] = static_cast<float>(kDiskOuterFactor * config.black_hole_mass);
    params[30] = config.disk_temperature_scale;

    // [31..34] tile origin/size: written per tile.
    params[35] = starfield_enabled ? 1.0f : 0.0f;
    params[36] = static_cast<float>(starfield_width);
    params[37] = static_cast<float>(starfield_height);
    params[38] = scene.warp_vs;
    params[39] = scene.warp_sigma;
    params[40] = scene.warp_R;

    // Ray bundles (P2): propagate two deviation vectors in the trace kernel and
    // carry the geometric-mean expansion in radiance alpha. Off by default so the
    // parity and pinned renders are unchanged.
    params[41] = config.ray_bundles ? 1.0f : 0.0f;

    // Morris-Thorne (dispatch id 3): throat radius and shape selector.
    params[42] = scene.throat_b0;
    params[43] = scene.worm_shape;
    params[44] = 0.5f;  // Per-dispatch sub-pixel sample u.
    params[45] = 0.5f;  // Per-dispatch sub-pixel sample v.
    params[46] = 0.0f;  // Running-average sample index.
    params[47] = static_cast<float>(config.camera_beta_forward);
    params[48] = static_cast<float>(config.camera_beta_up);
    params[49] = static_cast<float>(config.camera_beta_right);
    params[50] = config.temperature_model == DiskTemperatureModel::ShakuraSunyaev ? 1.0f : 0.0f;
    params[51] = config.doppler_beaming ? 1.0f : 0.0f;
    params[52] = config.lens_type == core::LensType::ThinLens ? 1.0f : 0.0f;
    params[53] = config.camera_focal_length;
    params[54] = config.camera_aperture;
    params[55] = config.camera_focus_distance;
    params[56] = config.point_starfield ? 1.0f : 0.0f;
    params[57] =
        static_cast<float>((config.camera_fov * math::kPi / 180.0) / std::max(1, config.height));
    params[58] = config.point_starfield_config.brightness_scale;
    params[59] = config.enable_volumetric_disk ? 1.0f : 0.0f;
    params[60] = config.volumetric_h_over_r;
    params[61] = config.volumetric_h_power;
    params[62] = config.volumetric_tau_midplane;
    params[63] = static_cast<float>(config.volumetric_samples);
    params[64] = config.enable_turbulence ? 1.0f : 0.0f;
}

}  // namespace

Expected<void> ValidateVulkanRenderConfig(const SessionConfig& config) {
    if (const auto issue = SessionConfigIssue(config); issue.has_value()) {
        return Fail(ErrorDomain::kConfiguration, "validate Vulkan render configuration", *issue);
    }
    const SessionConfig defaults;
    if (config.tile_size != defaults.tile_size || config.thread_count != defaults.thread_count ||
        config.enable_parallel_rendering != defaults.enable_parallel_rendering) {
        return Fail(ErrorDomain::kConfiguration, "validate Vulkan render configuration",
                    "CPU tile and thread controls are not represented by the device-governed "
                    "Vulkan backend");
    }
    if (!config.enable_disk && (config.enable_volumetric_disk || config.enable_turbulence)) {
        return Fail(ErrorDomain::kConfiguration, "validate Vulkan render configuration",
                    "volumetric disk and turbulence require the disk");
    }
    if (config.enable_polarisation) {
        return Fail(ErrorDomain::kConfiguration, "validate Vulkan render configuration",
                    "rendered polarisation is not represented by the live Vulkan kernel");
    }
    if (config.color_mode != core::color_modes::Mode::TrueColor) {
        return Fail(ErrorDomain::kConfiguration, "validate Vulkan render configuration",
                    "the requested diagnostic colour mode is not represented by the Vulkan "
                    "render kernel");
    }
    if (config.enable_disk &&
        core::DiskSupportFor(config.metric_id) != core::DiskSupport::PageThorne) {
        return Fail(ErrorDomain::kConfiguration, "validate Vulkan render configuration",
                    "the selected metric has no represented accretion-disk emission model");
    }
    switch (config.temperature_model) {
        case DiskTemperatureModel::NovikovThorne:
        case DiskTemperatureModel::ShakuraSunyaev:
            break;
        default:
            return Fail(ErrorDomain::kConfiguration, "validate Vulkan render configuration",
                        "invalid disk temperature model");
    }
    switch (config.lens_type) {
        case core::LensType::Pinhole:
        case core::LensType::ThinLens:
            break;
        case core::LensType::Fisheye:
            return Fail(ErrorDomain::kConfiguration, "validate Vulkan render configuration",
                        "fisheye lens is not represented by the Vulkan render kernel");
        default:
            return Fail(ErrorDomain::kConfiguration, "validate Vulkan render configuration",
                        "invalid lens model");
    }
    if ((config.enable_volumetric_disk || config.enable_turbulence) &&
        config.volumetric_samples > kMaxVulkanVolumeSamples) {
        return Fail(ErrorDomain::kDevice, "validate Vulkan render configuration",
                    "Vulkan volumetric transfer supports at most " +
                        std::to_string(kMaxVulkanVolumeSamples) +
                        " midpoint samples per geodesic segment; reduce volumetric.samples or "
                        "use --cpu");
    }
    return {};
}

Expected<VulkanRenderStats> RenderVulkanToDisplay(const SessionConfig& config,
                                                  DisplayBuffer& display,
                                                  const std::function<void(int, int)>& on_tile,
                                                  const std::function<bool()>& should_cancel) {
    const auto start = std::chrono::steady_clock::now();

    if (auto compatible = ValidateVulkanRenderConfig(config); !compatible) {
        return std::unexpected(compatible.error());
    }

    auto scene = MapMetric(config);
    if (!scene) {
        return std::unexpected(scene.error());
    }
    scene->disk_enabled = scene->disk_enabled && config.enable_disk;

    // Resolve the same explicit device identity used by readiness, direct
    // backend tests, and physical attestations.
    auto devices = backend::EnumerateVulkanDevices();
    if (!devices) {
        return std::unexpected(devices.error());
    }
    if (devices->empty()) {
        return Fail(ErrorDomain::kDevice, "open Vulkan render device",
                    "no Vulkan device present (no ICD/loader)");
    }
    auto device_index = backend::ResolveVulkanDeviceIndex(*devices);
    if (!device_index) {
        return std::unexpected(device_index.error());
    }
    auto device_result = backend::CreateVulkanDevice(*device_index);
    if (!device_result) {
        return std::unexpected(device_result.error());
    }
    backend::ComputeDevice& device = **device_result;
    const backend::DeviceInfo& info = device.Info();

    auto rung = SelectPrecisionRung(info.supports_fp64);
    if (!rung) {
        return std::unexpected(rung.error());
    }

    const auto spirv = LoadSpirv(*rung);
    if (spirv.empty()) {
        return Fail(
            ErrorDomain::kKernel, "load trace kernel",
#if SIRIUS_RELEASE_RESOURCE_LOCKED
            std::string(RungKernelName(*rung)) + " not found in the packaged Sirius volume");
#else
            std::string(RungKernelName(*rung)) +
                " not found below the Sirius runtime resource root "
                "(set SIRIUS_RESOURCE_DIR to an installed share/sirius directory)");
#endif
    }
    auto kernel = device.LoadKernel(spirv);
    if (!kernel) {
        return std::unexpected(kernel.error());
    }

    // Exactly one escaped-ray scene source is active. Texture mode requires the
    // packaged equirectangular asset; point mode owns its generated catalogue and
    // must not load or require an unused texture.
    int sf_w = 0;
    int sf_h = 0;
    std::vector<std::uint32_t> starfield;
    const bool use_starfield = !config.point_starfield;
    if (use_starfield) {
        starfield = LoadStarfield(sf_w, sf_h);
        if (starfield.empty()) {
            return Fail(ErrorDomain::kIo, "load Vulkan starfield",
                        "required runtime resource assets/Starfield.png is missing or unreadable");
        }
    }
    const std::uint64_t params_overhead = kParamCount * sizeof(float);
    const std::uint64_t starfield_bytes = starfield.size() * sizeof(std::uint32_t);
    std::unique_ptr<core::StarfieldSpatialIndex> point_index;
    if (config.point_starfield) {
        const core::StarfieldConfig catalogue_config =
            core::ExpandPointStarfieldConfig(config.point_starfield_config);
        core::StarfieldGenerator generator(catalogue_config);
        point_index = std::make_unique<core::StarfieldSpatialIndex>(generator.GenerateCatalogue());
    }
    const std::uint64_t point_bytes =
        (point_index ? point_index->Size() : 0) * sizeof(core::StarEntry) +
        (point_index ? point_index->MemoryBytes() : 0);
    const auto resolved_budget = ResolveBudgetBytes(info.render_memory_bytes);
    if (!resolved_budget) {
        return std::unexpected(resolved_budget.error());
    }
    const std::uint64_t budget = *resolved_budget;
    if (budget == 0) {
        return Fail(ErrorDomain::kDevice, "resolve memory budget",
                    "device reports no host-visible coherent render memory and no "
                    "SIRIUS_MEMORY_BUDGET_MB override");
    }

    auto plan =
        DeriveTilePlan(budget, config.width, config.height,
                       params_overhead + (use_starfield ? starfield_bytes : 0) + point_bytes);
    if (!plan) {
        return std::unexpected(plan.error());
    }

    std::cout << "[Vulkan] device[" << *device_index << "]: " << info.name << " ("
              << backend::ToString(info.kind) << ")\n";
    std::cout << "[Vulkan] budget: " << (budget / (1024 * 1024))
              << " MiB, tile: " << plan->tile_edge << "x" << plan->tile_edge << " (working set "
              << (plan->tile_working_set_bytes / 1024) << " KiB)\n";
    std::cout << "[Vulkan] precision: " << RungName(*rung) << " rung\n";
    if (config.point_starfield) {
        std::cout << "[Vulkan] point-source star field: " << point_index->Size() << " stars, "
                  << (point_index->MemoryBytes() / 1024) << " KiB index\n";
    }

    // Dispatch-duration governor: tiles are submitted as adaptive row bands so
    // no single dispatch can outlive the OS GPU watchdog (dispatch_governor.h).
    auto target_ms = ResolveDispatchTargetMs();
    if (!target_ms) {
        return std::unexpected(target_ms.error());
    }
    BandController bands(plan->tile_edge, *target_ms);
    if (bands.Enabled()) {
        std::cout << "[Vulkan] dispatch governor: " << *target_ms << " ms/band target\n";
    } else {
        std::cout << "[Vulkan] dispatch governor disabled: one dispatch per tile\n";
    }

    const int edge = plan->tile_edge;
    const std::uint64_t radiance_capacity =
        static_cast<std::uint64_t>(edge) * edge * 4 * sizeof(float);

    auto radiance_buf = device.CreateBuffer(radiance_capacity, backend::BufferUsage::kStorage);
    auto params_buf =
        device.CreateBuffer(kParamCount * sizeof(float), backend::BufferUsage::kStorage);
    const std::vector<std::uint32_t> dummy_star = {0u};
    const std::span<const std::uint32_t> star_span =
        use_starfield ? std::span<const std::uint32_t>(starfield)
                      : std::span<const std::uint32_t>(dummy_star);
    auto star_buf = device.CreateBuffer(star_span.size_bytes(), backend::BufferUsage::kStorage);
    const std::array<float, 8> dummy_point_star{};
    const std::array<std::uint32_t, 2> dummy_offsets{};
    const std::array<std::uint32_t, 1> dummy_indices{};
    const std::span<const core::StarEntry> point_star_span =
        point_index ? point_index->Stars() : std::span<const core::StarEntry>();
    const std::span<const std::uint32_t> point_offset_span =
        point_index ? point_index->Offsets() : std::span<const std::uint32_t>(dummy_offsets);
    const std::span<const std::uint32_t> point_index_span =
        point_index ? point_index->Indices() : std::span<const std::uint32_t>(dummy_indices);
    auto point_star_buf = device.CreateBuffer(
        point_star_span.empty() ? sizeof(dummy_point_star) : point_star_span.size_bytes(),
        backend::BufferUsage::kStorage);
    auto point_offset_buf =
        device.CreateBuffer(point_offset_span.size_bytes(), backend::BufferUsage::kStorage);
    auto point_index_buf =
        device.CreateBuffer(point_index_span.size_bytes(), backend::BufferUsage::kStorage);
    if (!radiance_buf || !params_buf || !star_buf || !point_star_buf || !point_offset_buf ||
        !point_index_buf) {
        return Fail(ErrorDomain::kDevice, "allocate render buffers",
                    "tile radiance, params, or starfield buffer allocation failed");
    }
    if (auto w = device.WriteBuffer(*star_buf, std::as_bytes(star_span)); !w) {
        return std::unexpected(w.error());
    }
    const auto raw_point_stars = point_star_span.empty()
                                     ? std::as_bytes(std::span<const float>(dummy_point_star))
                                     : std::as_bytes(point_star_span);
    if (auto w = device.WriteBuffer(*point_star_buf, raw_point_stars); !w) {
        return std::unexpected(w.error());
    }
    if (auto w = device.WriteBuffer(*point_offset_buf, std::as_bytes(point_offset_span)); !w) {
        return std::unexpected(w.error());
    }
    if (auto w = device.WriteBuffer(*point_index_buf, std::as_bytes(point_index_span)); !w) {
        return std::unexpected(w.error());
    }

    std::vector<float> params(kParamCount, 0.0f);
    FillSceneParams(params, config, *scene, use_starfield, sf_w, sf_h);

    const int tiles_x = (config.width + edge - 1) / edge;
    const int tiles_y = (config.height + edge - 1) / edge;
    const int tiles_total = tiles_x * tiles_y;
    std::vector<float> tile_pixels(static_cast<std::size_t>(edge) * edge * 4, 0.0f);
    // Transactional host staging: cancellation or a late device error must not
    // leave a partial frame in the caller's display buffer.
    std::vector<float> frame_pixels(static_cast<std::size_t>(config.width) * config.height * 4,
                                    0.0f);

    int tiles_done = 0;
    int band_dispatches = 0;
    const backend::BufferHandle bindings[] = {*radiance_buf,   *params_buf,       *star_buf,
                                              *point_star_buf, *point_offset_buf, *point_index_buf};
    for (int tj = 0; tj < tiles_y; ++tj) {
        for (int ti = 0; ti < tiles_x; ++ti) {
            if (should_cancel && should_cancel()) {
                return Fail(ErrorDomain::kInternal, "render Vulkan frame", "cancelled");
            }
            const int ox = ti * edge;
            const int oy = tj * edge;
            const int tw = std::min(edge, config.width - ox);
            const int th = std::min(edge, config.height - oy);

            params[31] = static_cast<float>(ox);
            params[33] = static_cast<float>(tw);

            // The kernel addresses pixels absolutely (tile origin + thread id),
            // so submitting the tile as row bands changes submission
            // granularity and no pixel value.
            for (int by = 0; by < th;) {
                if (should_cancel && should_cancel()) {
                    return Fail(ErrorDomain::kInternal, "render Vulkan frame", "cancelled");
                }
                const int bh = bands.NextRows(th - by, tw);
                params[32] = static_cast<float>(oy + by);
                params[34] = static_cast<float>(bh);
                const auto gx = static_cast<std::uint32_t>((tw + 7) / 8);
                const auto gy = static_cast<std::uint32_t>((bh + 7) / 8);
                int sample_index = 0;
                std::optional<base::Error> sample_error;
                ForEachPixelSample(config.samples_per_pixel, [&](float sample_u, float sample_v) {
                    if (sample_error.has_value()) return;
                    if (should_cancel && should_cancel()) {
                        sample_error =
                            base::Error{ErrorDomain::kInternal, "render Vulkan frame", "cancelled"};
                        return;
                    }
                    params[44] = sample_u;
                    params[45] = sample_v;
                    params[46] = static_cast<float>(sample_index);
                    if (auto w = device.WriteBuffer(*params_buf,
                                                    std::as_bytes(std::span<const float>(params)));
                        !w) {
                        sample_error = w.error();
                        return;
                    }

                    const auto band_start = std::chrono::steady_clock::now();
                    if (auto d = device.Dispatch(*kernel, bindings, gx, gy, 1); !d) {
                        sample_error = d.error();
                        return;
                    }
                    bands.Record(static_cast<std::int64_t>(tw) * bh,
                                 std::chrono::duration<double, std::milli>(
                                     std::chrono::steady_clock::now() - band_start)
                                     .count());
                    ++sample_index;
                    ++band_dispatches;
                });
                if (sample_error.has_value()) {
                    return std::unexpected(std::move(*sample_error));
                }

                const std::size_t band_floats = static_cast<std::size_t>(tw) * bh * 4;
                if (auto rd = device.ReadBuffer(
                        *radiance_buf,
                        std::as_writable_bytes(std::span<float>(tile_pixels.data(), band_floats)));
                    !rd) {
                    return std::unexpected(rd.error());
                }

                for (int row = 0; row < bh; ++row) {
                    const auto source =
                        tile_pixels.begin() + static_cast<std::ptrdiff_t>(row * tw * 4);
                    const auto destination =
                        frame_pixels.begin() +
                        static_cast<std::ptrdiff_t>(((oy + by + row) * config.width + ox) * 4);
                    std::copy_n(source, static_cast<std::size_t>(tw) * 4, destination);
                }
                by += bh;
            }

            ++tiles_done;
            if (on_tile) {
                on_tile(tiles_done, tiles_total);
            }
        }
    }

    if (should_cancel && should_cancel()) {
        return Fail(ErrorDomain::kInternal, "render Vulkan frame", "cancelled");
    }
    if (const auto bad = std::find_if_not(frame_pixels.begin(), frame_pixels.end(),
                                          [](float value) { return std::isfinite(value); });
        bad != frame_pixels.end()) {
        return Fail(ErrorDomain::kKernel, "validate Vulkan radiance",
                    std::format("non-finite sample at channel {}",
                                std::distance(frame_pixels.begin(), bad)));
    }
    display.UpdateTile(0, 0, config.width, config.height, frame_pixels.data());

    VulkanRenderStats stats;
    stats.device_name = info.name;
    stats.metric_name = scene->metric_name;
    stats.tile_plan = *plan;
    stats.precision = *rung;
    stats.starfield_uploaded = use_starfield;
    stats.point_catalogue_uploaded = config.point_starfield;
    stats.tiles_rendered = tiles_total;
    stats.band_dispatches = band_dispatches;
    stats.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    return stats;
}

}  // namespace sirius::render
