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
#include "sirius/core/constants.h"
#include "sirius/core/coordinates.h"
#include "sirius/core/disk/novikov_thorne_disk.h"
#include "sirius/core/metrics/registry.h"
#include "sirius/render/session/display_buffer.h"
#include "sirius/render/session/render_session.h"

#include "stb_image.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <format>
#include <fstream>
#include <iostream>
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
// capture_factor mirror the CPU session's TracerConfig (200 M, horizon x 1.05)
// so the shadow and background geometry coincide; the step schedule is the
// kernel's own (the integrator differs from the CPU RK45 by design).
constexpr float kTraceMaxSteps = 3000.0f;
constexpr float kTraceStepScale = 0.08f;
constexpr float kTraceMinStep = 0.02f;
constexpr float kTraceMaxStep = 2.0f;
constexpr float kTraceEscapeRadius = 200.0f;
constexpr float kTraceCaptureFactor = 1.05f;

// Thin-disk extent (units of M) and a cinematic radiance scale for the emitted
// blackbody colour. The outer edge matches the CPU session's disk_outer = 20 M;
// the inner edge is the ISCO. The emission scale is a display choice (the disk
// is not a parity-gated quantity), sized so the disk reads under the default
// ACES exposure without clipping the starfield.
constexpr float kDiskOuterFactor = 20.0f;
constexpr float kDiskEmissionScale = 2.0f;

// params buffer: 48 floats (indices documented in trace.slang). A power-of-two
// slack over the 43 slots used leaves room for later fields without an ABI move.
constexpr std::uint32_t kParamCount = 48;

// Kernel dispatch ids (must match trace.slang / gr_types.slang).
constexpr float kDispatchKerrSchild = 0.0f;
constexpr float kDispatchWarpDrive = 2.0f;

// What the trace kernel needs to render one metric, or a decline.
struct KernelScene {
    float metric_id = kDispatchKerrSchild;
    float M = 0.0f;
    float a = 0.0f;       // absolute spin a = (a/M) * M
    float Q = 0.0f;
    float Lambda = 0.0f;
    float warp_vs = 0.0f;
    float warp_sigma = 0.0f;
    float warp_R = 0.0f;
    bool disk_enabled = false;
    std::string metric_name;
};

// Maps a session metric onto the trace kernel, or declines. The gpu-renderable
// set is exactly the registry's gpu_supported metrics the Cartesian render path
// carries: Minkowski/Schwarzschild/Kerr (Kerr-Schild family) and Alcubierre.
// Charge/lambda metrics (registry gpu_supported=false) and Morris-Thorne (a
// spherical chart whose render-session integration is deferred; its device
// physics is covered by the parity probe) decline loudly.
[[nodiscard]] Expected<KernelScene> MapMetric(const SessionConfig& config) {
    KernelScene scene;
    const auto info = core::MetricInfoFor(config.metricId);
    scene.metric_name = info.canonical_name;

    switch (config.metricId) {
        case core::MetricId::Minkowski:
            scene.metric_id = kDispatchKerrSchild;
            scene.M = 0.0f;
            scene.disk_enabled = false;
            return scene;
        case core::MetricId::Schwarzschild:
            scene.metric_id = kDispatchKerrSchild;
            scene.M = static_cast<float>(config.blackHoleMass);
            scene.disk_enabled = true;
            return scene;
        case core::MetricId::Kerr:
            scene.metric_id = kDispatchKerrSchild;
            scene.M = static_cast<float>(config.blackHoleMass);
            scene.a = static_cast<float>(config.blackHoleSpin * config.blackHoleMass);
            scene.disk_enabled = true;
            return scene;
        case core::MetricId::Alcubierre:
            scene.metric_id = kDispatchWarpDrive;
            scene.warp_vs = static_cast<float>(config.warpVelocity);
            scene.warp_sigma = static_cast<float>(config.bubbleSigma);
            scene.warp_R = static_cast<float>(config.bubbleRadius);
            scene.disk_enabled = false;
            return scene;
        case core::MetricId::MorrisThorne:
            return Fail(ErrorDomain::kDevice, "select Vulkan render metric",
                        "Morris-Thorne renders in a spherical chart; its Vulkan render-session "
                        "integration is a recorded follow-up (the device kernel physics is covered "
                        "by the kernel parity probe)");
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

// Selects the precision-ladder rung. The fp64 trace kernel is deferred (the fp32
// path passes the parity gate); an explicit fp64 request therefore declines
// loudly rather than silently running fp32 (docs/STYLE.md section 4).
[[nodiscard]] Expected<PrecisionRung> SelectPrecisionRung(bool device_supports_fp64) {
    const char* precision = std::getenv("SIRIUS_PRECISION");
    const bool wants_fp64 = precision != nullptr && std::string(precision) == "fp64";
    if (!wants_fp64) {
        return PrecisionRung::Fp32;
    }
    if (!device_supports_fp64) {
        return Fail(ErrorDomain::kDevice, "select precision rung",
                    "SIRIUS_PRECISION=fp64 requested but the device lacks shaderFloat64");
    }
    return Fail(ErrorDomain::kDevice, "select precision rung",
                "SIRIUS_PRECISION=fp64 requested and supported, but the fp64 trace kernel is "
                "deferred (the fp32 path passes the parity gate); unset SIRIUS_PRECISION for the "
                "fp32 rung");
}

// Candidate locations for a build artefact, searched relative to the working
// directory (the CLI and tests run from various depths of the tree).
[[nodiscard]] std::vector<std::string> ResourceCandidates(const std::string& relative) {
    std::vector<std::string> paths;
#ifdef SIRIUS_KERNEL_DIR
    if (relative == "trace.spv") {
        paths.push_back(std::string(SIRIUS_KERNEL_DIR) + "/trace.spv");
    }
#endif
    for (const char* prefix : {"", "../", "../../", "../../../", "../../../../"}) {
        paths.push_back(std::string(prefix) + relative);
    }
    return paths;
}

[[nodiscard]] std::vector<std::uint32_t> LoadSpirv() {
    for (const auto& path : ResourceCandidates("trace.spv")) {
        std::ifstream file(path, std::ios::binary | std::ios::ate);
        if (!file) continue;
        const auto size = static_cast<std::size_t>(file.tellg());
        if (size == 0 || size % sizeof(std::uint32_t) != 0) continue;
        std::vector<std::uint32_t> words(size / sizeof(std::uint32_t));
        file.seekg(0);
        file.read(reinterpret_cast<char*>(words.data()), static_cast<std::streamsize>(size));
        return words;
    }
    return {};
}

// Loads and packs the equirectangular starfield to one uint32 (RGBA8) per texel,
// matching the CPU path's byte layout. Returns an empty vector when no asset is
// found; the caller then renders the analytic gradient background.
[[nodiscard]] std::vector<std::uint32_t> LoadStarfield(int& out_width, int& out_height) {
    for (const auto& path : ResourceCandidates("assets/Starfield.png")) {
        int channels = 0;
        unsigned char* pixels = stbi_load(path.c_str(), &out_width, &out_height, &channels, 4);
        if (pixels == nullptr) continue;
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
    out_width = 0;
    out_height = 0;
    return {};
}

// Fills the scene-invariant slots of the params buffer (everything except the
// per-tile origin/size in [32..35], written per dispatch).
void FillSceneParams(std::vector<float>& params, const SessionConfig& config,
                     const KernelScene& scene, bool starfield_enabled, int starfield_width,
                     int starfield_height) {
    const double theta = config.observerInclination;
    const double st = std::sin(theta);
    const double ct = std::cos(theta);
    const double r = config.observerDistance;

    params[0] = static_cast<float>(config.width);
    params[1] = static_cast<float>(config.height);
    params[2] = scene.metric_id;
    params[3] = scene.M;
    params[4] = scene.a;
    params[5] = scene.Q;
    params[6] = scene.Lambda;

    // camPos = BlToCartesian(r, theta, phi=0); camDir = -r_hat; camRight = phi_hat;
    // camUp = +theta_hat (kernel NDC v runs top-to-bottom, matching CPU image up).
    params[7] = static_cast<float>(r * st);
    params[8] = 0.0f;
    params[9] = static_cast<float>(r * ct);
    params[10] = static_cast<float>(-st);
    params[11] = 0.0f;
    params[12] = static_cast<float>(-ct);
    params[13] = 0.0f;
    params[14] = 1.0f;
    params[15] = 0.0f;
    params[16] = static_cast<float>(ct);
    params[17] = 0.0f;
    params[18] = static_cast<float>(-st);

    params[19] = static_cast<float>(config.cameraFOV * math::kPi / 180.0);  // full-angle radians
    params[20] = static_cast<float>(config.width) / static_cast<float>(config.height);

    params[21] = kTraceMaxSteps;
    params[22] = kTraceStepScale;
    params[23] = kTraceMinStep;
    params[24] = kTraceMaxStep;
    params[25] = kTraceEscapeRadius;
    params[26] = kTraceCaptureFactor;

    const double spin = (config.blackHoleMass > 0.0) ? config.blackHoleSpin : 0.0;
    const double isco = core::AccretionDiskD::ComputeIsco(spin);
    params[27] = scene.disk_enabled ? 1.0f : 0.0f;
    params[28] = static_cast<float>(isco * config.blackHoleMass);
    params[29] = static_cast<float>(kDiskOuterFactor * config.blackHoleMass);
    params[30] = config.diskTemperatureScale;
    params[31] = kDiskEmissionScale;

    // [32..35] tile origin/size: written per tile.
    params[36] = starfield_enabled ? 1.0f : 0.0f;
    params[37] = static_cast<float>(starfield_width);
    params[38] = static_cast<float>(starfield_height);
    params[39] = 0.0f;  // reserved (background fallback floor)
    params[40] = scene.warp_vs;
    params[41] = scene.warp_sigma;
    params[42] = scene.warp_R;
}

}  // namespace

Expected<VulkanRenderStats> RenderVulkanToDisplay(
    const SessionConfig& config, DisplayBuffer& display,
    const std::function<void(int, int)>& on_tile) {
    const auto start = std::chrono::steady_clock::now();

    auto scene = MapMetric(config);
    if (!scene) {
        return std::unexpected(scene.error());
    }

    // Open the first Vulkan device (Lavapipe here). An empty enumeration is a
    // decline: no fabricated device, the caller falls back or reports.
    auto devices = backend::EnumerateVulkanDevices();
    if (!devices) {
        return std::unexpected(devices.error());
    }
    if (devices->empty()) {
        return Fail(ErrorDomain::kDevice, "open Vulkan render device",
                    "no Vulkan device present (no ICD/loader)");
    }
    auto device_result = backend::CreateVulkanDevice(0);
    if (!device_result) {
        return std::unexpected(device_result.error());
    }
    backend::ComputeDevice& device = **device_result;
    const backend::DeviceInfo& info = device.Info();

    auto rung = SelectPrecisionRung(info.supports_fp64);
    if (!rung) {
        return std::unexpected(rung.error());
    }

    const auto spirv = LoadSpirv();
    if (spirv.empty()) {
        return Fail(ErrorDomain::kKernel, "load trace kernel",
                    "trace.spv not found (build the kernels or set SIRIUS_KERNEL_DIR)");
    }
    auto kernel = device.LoadKernel(spirv);
    if (!kernel) {
        return std::unexpected(kernel.error());
    }

    // Starfield: uploaded when it fits the governed budget as fixed device
    // residency, else dropped for the analytic gradient background so a tight
    // budget still renders (a loud log, not a decline).
    int sf_w = 0;
    int sf_h = 0;
    std::vector<std::uint32_t> starfield = LoadStarfield(sf_w, sf_h);
    const std::uint64_t params_overhead = kParamCount * sizeof(float);
    const std::uint64_t starfield_bytes = starfield.size() * sizeof(std::uint32_t);
    const std::uint64_t budget = ResolveBudgetBytes(info.device_local_bytes);
    if (budget == 0) {
        return Fail(ErrorDomain::kDevice, "resolve memory budget",
                    "device reports no device-local memory and no SIRIUS_MEMORY_BUDGET_MB override");
    }

    bool use_starfield = !starfield.empty();
    auto plan = DeriveTilePlan(budget, config.width, config.height,
                               params_overhead + (use_starfield ? starfield_bytes : 0));
    if (!plan && use_starfield) {
        // Retry without the starfield: it is optional background, not working set.
        use_starfield = false;
        std::cout << "[Vulkan] starfield dropped under the memory budget; rendering the analytic "
                     "gradient background\n";
        plan = DeriveTilePlan(budget, config.width, config.height, params_overhead);
    }
    if (!plan) {
        return std::unexpected(plan.error());
    }

    std::cout << "[Vulkan] device: " << info.name << " (" << backend::ToString(info.kind) << ")\n";
    std::cout << "[Vulkan] budget: " << (budget / (1024 * 1024)) << " MiB, tile: "
              << plan->tile_edge << "x" << plan->tile_edge << " (working set "
              << (plan->tile_working_set_bytes / 1024) << " KiB)\n";
    std::cout << "[Vulkan] precision: fp32 rung\n";

    const int edge = plan->tile_edge;
    const std::uint64_t radiance_capacity =
        static_cast<std::uint64_t>(edge) * edge * 4 * sizeof(float);

    auto radiance_buf = device.CreateBuffer(radiance_capacity, backend::BufferUsage::kStorage);
    auto params_buf = device.CreateBuffer(kParamCount * sizeof(float),
                                          backend::BufferUsage::kStorage);
    const std::vector<std::uint32_t> dummy_star = {0u};
    const std::span<const std::uint32_t> star_span =
        use_starfield ? std::span<const std::uint32_t>(starfield)
                      : std::span<const std::uint32_t>(dummy_star);
    auto star_buf = device.CreateBuffer(star_span.size_bytes(), backend::BufferUsage::kStorage);
    if (!radiance_buf || !params_buf || !star_buf) {
        return Fail(ErrorDomain::kDevice, "allocate render buffers",
                    "tile radiance, params, or starfield buffer allocation failed");
    }
    if (auto w = device.WriteBuffer(*star_buf, std::as_bytes(star_span)); !w) {
        return std::unexpected(w.error());
    }

    std::vector<float> params(kParamCount, 0.0f);
    FillSceneParams(params, config, *scene, use_starfield, sf_w, sf_h);

    const int tiles_x = (config.width + edge - 1) / edge;
    const int tiles_y = (config.height + edge - 1) / edge;
    const int tiles_total = tiles_x * tiles_y;
    std::vector<float> tile_pixels(static_cast<std::size_t>(edge) * edge * 4, 0.0f);

    int tiles_done = 0;
    const backend::BufferHandle bindings[] = {*radiance_buf, *params_buf, *star_buf};
    for (int tj = 0; tj < tiles_y; ++tj) {
        for (int ti = 0; ti < tiles_x; ++ti) {
            const int ox = ti * edge;
            const int oy = tj * edge;
            const int tw = std::min(edge, config.width - ox);
            const int th = std::min(edge, config.height - oy);

            params[32] = static_cast<float>(ox);
            params[33] = static_cast<float>(oy);
            params[34] = static_cast<float>(tw);
            params[35] = static_cast<float>(th);
            if (auto w = device.WriteBuffer(*params_buf,
                                            std::as_bytes(std::span<const float>(params)));
                !w) {
                return std::unexpected(w.error());
            }

            const auto gx = static_cast<std::uint32_t>((tw + 7) / 8);
            const auto gy = static_cast<std::uint32_t>((th + 7) / 8);
            if (auto d = device.Dispatch(*kernel, bindings, gx, gy, 1); !d) {
                return std::unexpected(d.error());
            }

            const std::size_t tile_floats = static_cast<std::size_t>(tw) * th * 4;
            if (auto rd = device.ReadBuffer(
                    *radiance_buf,
                    std::as_writable_bytes(std::span<float>(tile_pixels.data(), tile_floats)));
                !rd) {
                return std::unexpected(rd.error());
            }

            display.UpdateTile(ox, oy, tw, th, tile_pixels.data());
            ++tiles_done;
            if (on_tile) {
                on_tile(tiles_done, tiles_total);
            }
        }
    }

    VulkanRenderStats stats;
    stats.device_name = info.name;
    stats.metric_name = scene->metric_name;
    stats.tile_plan = *plan;
    stats.precision = *rung;
    stats.starfield_uploaded = use_starfield;
    stats.tiles_rendered = tiles_total;
    stats.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    return stats;
}

}  // namespace sirius::render
