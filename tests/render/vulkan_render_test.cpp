// Vulkan render-path gates (specification programmes 3 and 4). Three concerns:
//   1. a Kerr render through the RenderSession Vulkan path completes under a
//      constrained memory budget and yields finite, non-constant radiance with a
//      bounded horizon shadow (64x64 and 160x120, exercising >1 governed tile);
//   2. the same kernel dispatched directly agrees with the CPU reference tracer
//      on the same scene geometry within statistical bounds (the two integrators
//      differ by design, docs/ARCHITECTURE.md section 3, so the gate is on
//      radiance statistics, not bitwise equality);
//   3. everything skips cleanly where no Vulkan device or no compiled kernel is
//      present.
//
// All comparisons are on LINEAR radiance (the kernel and the CPU tracer both
// return linear radiance; tonemapping is host-side and out of scope here). The
// parity scene disables the disk and the starfield so the comparison isolates
// the integrator + camera geometry against the shared analytic gradient
// background; the disk emission models differ deliberately between the cinematic
// CPU shading and the physical kernel radiance and are not parity-gated.

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <numbers>
#include <span>
#include <string>
#include <vector>

#ifdef SIRIUS_HAS_VULKAN_BACKEND
#include "sirius/backend/cpu/geodesic_tracer.h"
#include "sirius/backend/device.h"
#include "sirius/core/camera.h"
#include "sirius/core/disk/novikov_thorne_disk.h"
#include "sirius/core/metrics/kerr_schild_family.h"
#include "sirius/core/metrics/morris_thorne_family.h"
#include "sirius/render/session/display_buffer.h"
#include "sirius/render/session/render_session.h"
#endif

namespace {

#ifdef SIRIUS_HAS_VULKAN_BACKEND

using sirius::backend::BufferHandle;
using sirius::backend::BufferUsage;
using sirius::backend::ComputeDevice;
using sirius::backend::CreateVulkanDevice;
using sirius::backend::EnumerateVulkanDevices;
using sirius::backend::GeodesicTracer;
using sirius::backend::TracerConfig;
using sirius::backend::TraceResult;
using sirius::core::CameraConfig;
using sirius::core::CameraRay;
using sirius::core::KerrSchildFamily;
using sirius::core::KerrSchildParams;
using sirius::core::PinholeCamera;

constexpr double kPi = std::numbers::pi;

std::vector<std::uint32_t> LoadSpirv(const std::string& path) {
    std::ifstream file(path, std::ios::binary | std::ios::ate);
    if (!file) return {};
    const auto size = static_cast<std::size_t>(file.tellg());
    std::vector<std::uint32_t> words(size / sizeof(std::uint32_t));
    file.seekg(0);
    file.read(reinterpret_cast<char*>(words.data()), static_cast<std::streamsize>(size));
    return words;
}

// A scene shared by the parity gate: Kerr M=1, a=0.9 seen from 30M at 80 deg
// inclination, disk and starfield off (analytic gradient background).
struct Scene {
    int width = 96;
    int height = 96;
    float M = 1.0f;
    float spin = 0.9f;  // a/M
    double distance = 30.0;
    double inclination_deg = 80.0;
    float fov_deg = 50.0f;
};

// The analytic gradient background the kernel and the CPU reference both sample
// (matches trace.slang SampleBackgroundGradient): 0.5 + 0.5 * dir_hat.
void GradientBackground(double dx, double dy, double dz, float& r, float& g, float& b) {
    const double len = std::sqrt(dx * dx + dy * dy + dz * dz);
    if (len < 1e-10) {
        r = g = b = 0.001f;
        return;
    }
    r = static_cast<float>(0.5 + 0.5 * dx / len);
    g = static_cast<float>(0.5 + 0.5 * dy / len);
    b = static_cast<float>(0.5 + 0.5 * dz / len);
}

float Luminance(float r, float g, float b) { return 0.2126f * r + 0.7152f * g + 0.0722f * b; }

// Builds the trace params for `scene` as a full single tile, disk and starfield
// disabled, with the same camera frame the Vulkan renderer builds.
std::vector<float> BuildTraceParams(const Scene& scene) {
    const double theta = scene.inclination_deg * kPi / 180.0;
    const double st = std::sin(theta);
    const double ct = std::cos(theta);
    const double r = scene.distance;

    std::vector<float> p(48, 0.0f);
    p[0] = static_cast<float>(scene.width);
    p[1] = static_cast<float>(scene.height);
    p[2] = 0.0f;  // Kerr-Schild dispatch
    p[3] = scene.M;
    p[4] = scene.spin * scene.M;  // absolute a
    p[7] = static_cast<float>(r * st);
    p[8] = 0.0f;
    p[9] = static_cast<float>(r * ct);
    p[10] = static_cast<float>(-st);
    p[11] = 0.0f;
    p[12] = static_cast<float>(-ct);
    p[13] = 0.0f;
    p[14] = 1.0f;
    p[15] = 0.0f;
    p[16] = static_cast<float>(ct);
    p[17] = 0.0f;
    p[18] = static_cast<float>(-st);
    p[19] = static_cast<float>(scene.fov_deg * kPi / 180.0);
    p[20] = static_cast<float>(scene.width) / static_cast<float>(scene.height);
    p[21] = 3000.0f;  // maxSteps
    p[22] = 0.08f;    // stepScale
    p[23] = 0.02f;    // minStep
    p[24] = 2.0f;     // maxStep
    p[25] = 200.0f;   // escapeRadius
    p[26] = 1.05f;    // captureFactor
    p[27] = 0.0f;     // disk disabled
    p[32] = 0.0f;     // tileOriginX
    p[33] = 0.0f;     // tileOriginY
    p[34] = static_cast<float>(scene.width);
    p[35] = static_cast<float>(scene.height);
    p[36] = 0.0f;  // starfield disabled -> gradient background
    p[37] = 1.0f;
    p[38] = 1.0f;
    return p;
}

// Opens the first Vulkan device and the trace kernel, or leaves `ready` false so
// the caller GTEST_SKIPs.
struct KernelFixture {
    std::unique_ptr<ComputeDevice> device;
    sirius::backend::KernelHandle kernel;
    bool ready = false;
};

KernelFixture OpenKernel() {
    KernelFixture f;
#ifdef SIRIUS_KERNEL_DIR
    const auto devices = EnumerateVulkanDevices();
    if (!devices.has_value() || devices->empty()) return f;
    auto device = CreateVulkanDevice(0);
    if (!device.has_value()) return f;
    const auto spirv = LoadSpirv(std::string(SIRIUS_KERNEL_DIR) + "/trace.spv");
    if (spirv.empty()) return f;
    auto kernel = (*device)->LoadKernel(spirv);
    if (!kernel.has_value()) return f;
    f.device = std::move(*device);
    f.kernel = *kernel;
    f.ready = true;
#endif
    return f;
}

// --- Session Vulkan path smoke ----------------------------------------------

// Renders a Kerr scene through the RenderSession Vulkan path to an EXR (which
// leaves the display buffer's linear radiance untouched) under a constrained
// budget, returning the linear RGBA buffer, or an empty vector when skipped.
std::vector<float> RenderSessionVulkan(int width, int height, const std::string& out_path) {
    const auto devices = EnumerateVulkanDevices();
    if (!devices.has_value() || devices->empty()) return {};

    sirius::render::SessionConfig config;
    config.width = width;
    config.height = height;
    config.samplesPerPixel = 1;
    config.tileSize = 64;
    config.enableParallelRendering = false;
    config.backend = sirius::render::RenderBackend::Vulkan;
    config.metricId = sirius::core::MetricId::Kerr;
    config.blackHoleMass = 1.0;
    config.blackHoleSpin = 0.9;
    config.observerDistance = 30.0;
    config.observerInclination = 80.0 * kPi / 180.0;
    config.cameraFOV = 50.0f;
    config.outputPath = out_path;

    sirius::render::RenderSession session;
    session.Configure(config);
    const sirius::render::SessionState state = session.Execute();
    if (state != sirius::render::SessionState::Complete) return {};

    sirius::render::DisplayBuffer& display = session.GetDisplayBuffer();
    const float* data = display.GetFloatData();
    return std::vector<float>(data, data + static_cast<std::size_t>(width) * height * 4);
}

void ExpectFiniteNonConstantWithShadow(const std::vector<float>& rgba, int width, int height,
                                       const char* tag) {
    ASSERT_EQ(rgba.size(), static_cast<std::size_t>(width) * height * 4) << tag;
    float min_lum = std::numeric_limits<float>::max();
    float max_lum = 0.0f;
    std::size_t shadow = 0;
    for (int p = 0; p < width * height; ++p) {
        const float r = rgba[p * 4 + 0];
        const float g = rgba[p * 4 + 1];
        const float b = rgba[p * 4 + 2];
        ASSERT_TRUE(std::isfinite(r) && std::isfinite(g) && std::isfinite(b))
            << tag << " non-finite radiance";
        const float lum = Luminance(r, g, b);
        min_lum = std::min(min_lum, lum);
        max_lum = std::max(max_lum, lum);
        if (r < 1e-4f && g < 1e-4f && b < 1e-4f) ++shadow;
    }
    const double fraction = static_cast<double>(shadow) / (width * height);
    std::cout << "[ " << tag << " ] shadow fraction=" << fraction << " luminance range=[" << min_lum
              << ", " << max_lum << "]\n";
    EXPECT_GT(max_lum - min_lum, 1e-3f) << tag << " radiance field is constant";
    EXPECT_GE(fraction, 0.005) << tag << " shadow fraction too small";
    EXPECT_LE(fraction, 0.60) << tag << " shadow fraction too large";
}

TEST(VulkanRenderSession, Kerr64CompletesUnderConstrainedBudgetWithFiniteRadiance) {
    if (const auto d = EnumerateVulkanDevices(); !d.has_value() || d->empty()) {
        GTEST_SKIP() << "no Vulkan device present";
    }
    setenv("SIRIUS_MEMORY_BUDGET_MB", "1", 1);  // constrained: drops the starfield
    const std::string out = std::string(std::getenv("TMPDIR") ? std::getenv("TMPDIR") : "/tmp") +
                            "/sirius_vk_smoke_64.exr";
    const auto rgba = RenderSessionVulkan(64, 64, out);
    unsetenv("SIRIUS_MEMORY_BUDGET_MB");
    ASSERT_FALSE(rgba.empty()) << "session Vulkan render did not complete";
    ExpectFiniteNonConstantWithShadow(rgba, 64, 64, "vk-session-64");
}

// The fp64 rung end to end: SIRIUS_PRECISION=fp64 must select trace_fp64.spv
// on a shaderFloat64 device and complete with a finite, non-constant frame
// carrying the same shadow structure as the fp32 rung. On devices without
// shaderFloat64 the session must fail (decline), never silently render fp32.
TEST(VulkanRenderSession, Fp64RungRendersOrDeclinesLoudly) {
    const auto devices = EnumerateVulkanDevices();
    if (!devices.has_value() || devices->empty()) {
        GTEST_SKIP() << "no Vulkan device present";
    }
    auto device = CreateVulkanDevice(0);
    ASSERT_TRUE(device.has_value());
    const bool fp64_supported = (*device)->Info().supports_fp64;
    device->reset();

    setenv("SIRIUS_PRECISION", "fp64", 1);
    const std::string out = std::string(std::getenv("TMPDIR") ? std::getenv("TMPDIR") : "/tmp") +
                            "/sirius_vk_fp64_64.exr";
    const auto rgba = RenderSessionVulkan(64, 64, out);
    unsetenv("SIRIUS_PRECISION");

    if (fp64_supported) {
        ASSERT_FALSE(rgba.empty()) << "fp64 rung requested and supported but did not complete";
        ExpectFiniteNonConstantWithShadow(rgba, 64, 64, "vk-session-fp64");
    } else {
        EXPECT_TRUE(rgba.empty()) << "fp64 request must decline on a device without shaderFloat64";
    }
}

// The compensated rung end to end: SIRIUS_PRECISION=fp32-comp selects
// trace_fp32comp.spv on any device and completes with a finite, non-constant
// frame carrying the shadow structure.
TEST(VulkanRenderSession, CompensatedRungRendersOnAnyDevice) {
    if (const auto d = EnumerateVulkanDevices(); !d.has_value() || d->empty()) {
        GTEST_SKIP() << "no Vulkan device present";
    }
    setenv("SIRIUS_PRECISION", "fp32-comp", 1);
    const std::string out = std::string(std::getenv("TMPDIR") ? std::getenv("TMPDIR") : "/tmp") +
                            "/sirius_vk_comp_64.exr";
    const auto rgba = RenderSessionVulkan(64, 64, out);
    unsetenv("SIRIUS_PRECISION");
    ASSERT_FALSE(rgba.empty()) << "compensated rung did not complete";
    ExpectFiniteNonConstantWithShadow(rgba, 64, 64, "vk-session-comp");
}

TEST(VulkanRenderSession, Kerr160x120CompletesAcrossMultipleGovernedTiles) {
    if (const auto d = EnumerateVulkanDevices(); !d.has_value() || d->empty()) {
        GTEST_SKIP() << "no Vulkan device present";
    }
    setenv("SIRIUS_MEMORY_BUDGET_MB", "1", 1);  // tile edge 120 -> 2 tiles for 160 wide
    const std::string out = std::string(std::getenv("TMPDIR") ? std::getenv("TMPDIR") : "/tmp") +
                            "/sirius_vk_smoke_160.exr";
    const auto rgba = RenderSessionVulkan(160, 120, out);
    unsetenv("SIRIUS_MEMORY_BUDGET_MB");
    ASSERT_FALSE(rgba.empty()) << "session Vulkan render did not complete";
    ExpectFiniteNonConstantWithShadow(rgba, 160, 120, "vk-session-160");
}

// --- CPU-versus-Vulkan parity -----------------------------------------------

TEST(VulkanRenderSession, CpuVulkanAgreeOnKerrGeometryWithinStatisticalBounds) {
    KernelFixture f = OpenKernel();
    if (!f.ready) GTEST_SKIP() << "no Vulkan device or trace kernel absent";

    const Scene scene;
    const int w = scene.width;
    const int h = scene.height;

    // --- Vulkan render (direct kernel dispatch, disk + starfield off) --------
    std::vector<float> params = BuildTraceParams(scene);
    std::vector<float> vk(static_cast<std::size_t>(w) * h * 4, 0.0f);
    const std::vector<std::uint32_t> star_dummy = {0u};

    const auto rbuf = f.device->CreateBuffer(vk.size() * sizeof(float), BufferUsage::kStorage);
    const auto pbuf = f.device->CreateBuffer(params.size() * sizeof(float), BufferUsage::kStorage);
    const auto sbuf =
        f.device->CreateBuffer(star_dummy.size() * sizeof(std::uint32_t), BufferUsage::kStorage);
    ASSERT_TRUE(rbuf && pbuf && sbuf);
    ASSERT_TRUE(f.device->WriteBuffer(*pbuf, std::as_bytes(std::span<const float>(params))));
    ASSERT_TRUE(
        f.device->WriteBuffer(*sbuf, std::as_bytes(std::span<const std::uint32_t>(star_dummy))));
    const BufferHandle bind[] = {*rbuf, *pbuf, *sbuf};
    ASSERT_TRUE(f.device->Dispatch(f.kernel, bind, (w + 7) / 8, (h + 7) / 8, 1).has_value());
    ASSERT_TRUE(f.device->ReadBuffer(*rbuf, std::as_writable_bytes(std::span<float>(vk))));

    // --- CPU reference (the geodesic tracer on the same scene) ---------------
    KerrSchildParams kp;
    kp.M = scene.M;
    kp.a = scene.spin * scene.M;
    KerrSchildFamily metric(kp);

    TracerConfig tc;
    tc.escape_radius = 200.0f;
    tc.horizon_factor = 1.05f;
    tc.max_steps = 20000;
    tc.enable_disk = false;
    tc.integrator.initial_step = 0.1f;
    tc.integrator.max_step = 2.0f;
    tc.integrator.min_step = 1e-5f;
    tc.integrator.abs_tolerance = 5e-6f;
    tc.integrator.rel_tolerance = 5e-6f;
    GeodesicTracer tracer(&metric, tc);

    CameraConfig cc;
    cc.r = scene.distance;
    cc.theta = scene.inclination_deg * kPi / 180.0;
    cc.phi = 0.0;
    cc.fov = scene.fov_deg;
    cc.width = w;
    cc.height = h;
    PinholeCamera camera(cc);

    std::vector<float> cpu(static_cast<std::size_t>(w) * h * 4, 0.0f);
    for (int y = 0; y < h; ++y) {
        for (int x = 0; x < w; ++x) {
            const CameraRay ray = camera.GenerateRay(x, y, 0.5f, 0.5f);
            const TraceResult res = tracer.Trace(ray);
            float r = 0.0f, g = 0.0f, b = 0.0f;
            if (res.outcome == TraceResult::Outcome::Escaped) {
                GradientBackground(res.final_direction(1), res.final_direction(2),
                                   res.final_direction(3), r, g, b);
            } else if (res.outcome == TraceResult::Outcome::Spiraling) {
                r = g = b = 0.02f;
            } else if (res.outcome == TraceResult::Outcome::MaxSteps) {
                r = g = b = 0.01f;
            }  // Horizon -> black
            const int idx = (y * w + x) * 4;
            cpu[idx + 0] = r;
            cpu[idx + 1] = g;
            cpu[idx + 2] = b;
        }
    }

    // --- Statistics ----------------------------------------------------------
    // Shadow fractions per backend, and per-pixel luminance agreement over the
    // background (both non-shadow) region. Disk-off + shared gradient background
    // isolate the integrator + camera geometry.
    std::size_t shadow_vk = 0, shadow_cpu = 0, bg_pixels = 0;
    double abs_sum = 0.0, max_rel = 0.0;
    constexpr float kShadowLum = 0.02f;  // below this luminance a pixel is shadow-like
    constexpr float kRelFloor = 0.05f;   // relative-diff denominator floor
    for (int p = 0; p < w * h; ++p) {
        const float lv = Luminance(vk[p * 4], vk[p * 4 + 1], vk[p * 4 + 2]);
        const float lc = Luminance(cpu[p * 4], cpu[p * 4 + 1], cpu[p * 4 + 2]);
        ASSERT_TRUE(std::isfinite(lv) && std::isfinite(lc));
        if (lv < kShadowLum) ++shadow_vk;
        if (lc < kShadowLum) ++shadow_cpu;
        if (lv >= kShadowLum && lc >= kShadowLum) {
            ++bg_pixels;
            const float d = std::abs(lv - lc);
            abs_sum += d;
            max_rel = std::max<double>(max_rel, d / std::max(lc, kRelFloor));
        }
    }

    const double frac_vk = static_cast<double>(shadow_vk) / (w * h);
    const double frac_cpu = static_cast<double>(shadow_cpu) / (w * h);
    const double mean_abs = bg_pixels > 0 ? abs_sum / static_cast<double>(bg_pixels) : 0.0;
    std::cout << "[ cpu-vulkan-parity ] shadow_vk=" << frac_vk << " shadow_cpu=" << frac_cpu
              << " |dfrac|=" << std::abs(frac_vk - frac_cpu) << " bg_pixels=" << bg_pixels
              << " mean|dlum|=" << mean_abs << " max_rel_dlum=" << max_rel << "\n";

    ASSERT_GT(bg_pixels, static_cast<std::size_t>(w * h) / 4)
        << "too few shared background pixels to compare";

    // Primary parity gate: the shared analytic background agrees on luminance.
    // This is the geometry gate the task specifies (mean abs diff, max rel diff);
    // it is tight because an identical camera frame makes both backends sample
    // the same escape directions where lensing is weak. A camera-mapping
    // regression blows these by orders of magnitude. Achieved on Lavapipe:
    // mean|dlum| ~6e-7, max_rel ~3e-5 (printed above); the bounds carry a wide
    // margin for transcendental-library and integrator drift.
    EXPECT_LT(mean_abs, 5e-3)
        << "mean background luminance difference too large (camera mismatch?)";
    EXPECT_LT(max_rel, 0.2) << "worst-case background luminance disagreement too large";

    // Secondary gate: each backend renders a horizon shadow within the broad
    // band kernel_trace_test uses, and the two fractions agree closely. The
    // historical 0.24-vs-0.11 gap (bounded loosely at 0.25) was not integrator
    // divergence: the CPU trace loop was terminating rays at their first
    // RK45 step rejection (outcome MaxSteps, luminance 0.01, classified as
    // shadow here), which near-critical photon-ring rays always hit. With
    // rejections retried, the measured gap on Lavapipe is ~1e-3; the bound
    // carries a 20x margin for driver and transcendental-library variation
    // while still failing a reintroduced ray-kill decisively.
    EXPECT_GE(frac_vk, 0.005);
    EXPECT_LE(frac_vk, 0.40) << "Vulkan shadow fraction outside the broad band";
    EXPECT_GE(frac_cpu, 0.005);
    EXPECT_LE(frac_cpu, 0.40) << "CPU shadow fraction outside the broad band";
    EXPECT_LT(std::abs(frac_vk - frac_cpu), 0.02)
        << "shadow fractions diverge: a capture-semantics or ray-kill regression";
}

// --- CPU-versus-Vulkan Morris-Thorne parity ----------------------------------
// The wormhole's Cartesian embedding runs on both backends (MorrisThorneCartesian
// on the CPU, gr_metrics GetMorrisThorneCartesian* under dispatch id 3 on the
// device); same one-sheet chart, same throat capture convention, same analytic
// gradient background. The gate mirrors the Kerr parity test: background
// luminance agreement plus throat-shadow fractions within the tightened band.
TEST(VulkanRenderSession, CpuVulkanAgreeOnMorrisThorneGeometryWithinStatisticalBounds) {
    KernelFixture f = OpenKernel();
    if (!f.ready) GTEST_SKIP() << "no Vulkan device or trace kernel absent";

    constexpr int w = 96;
    constexpr int h = 96;
    constexpr float kThroat = 4.0f;
    constexpr double kDistance = 30.0;
    constexpr double kInclinationDeg = 80.0;
    constexpr float kFovDeg = 50.0f;

    // --- Vulkan render (direct kernel dispatch, disk + starfield off) --------
    const double theta = kInclinationDeg * kPi / 180.0;
    const double st = std::sin(theta);
    const double ct = std::cos(theta);
    std::vector<float> params(48, 0.0f);
    params[0] = w;
    params[1] = h;
    params[2] = 3.0f;  // Morris-Thorne dispatch
    params[7] = static_cast<float>(kDistance * st);
    params[9] = static_cast<float>(kDistance * ct);
    params[10] = static_cast<float>(-st);
    params[12] = static_cast<float>(-ct);
    params[14] = 1.0f;
    params[16] = static_cast<float>(ct);
    params[18] = static_cast<float>(-st);
    params[19] = static_cast<float>(kFovDeg * kPi / 180.0);
    params[20] = 1.0f;
    params[21] = 3000.0f;
    params[22] = 0.08f;
    params[23] = 0.02f;
    params[24] = 2.0f;
    params[25] = 200.0f;
    params[26] = 1.05f;
    params[34] = w;
    params[35] = h;
    params[37] = 1.0f;
    params[38] = 1.0f;
    params[44] = kThroat;
    params[45] = 0.0f;  // Ellis

    std::vector<float> vk(static_cast<std::size_t>(w) * h * 4, 0.0f);
    const std::vector<std::uint32_t> star_dummy = {0u};
    const auto rbuf = f.device->CreateBuffer(vk.size() * sizeof(float), BufferUsage::kStorage);
    const auto pbuf = f.device->CreateBuffer(params.size() * sizeof(float), BufferUsage::kStorage);
    const auto sbuf =
        f.device->CreateBuffer(star_dummy.size() * sizeof(std::uint32_t), BufferUsage::kStorage);
    ASSERT_TRUE(rbuf && pbuf && sbuf);
    ASSERT_TRUE(f.device->WriteBuffer(*pbuf, std::as_bytes(std::span<const float>(params))));
    ASSERT_TRUE(
        f.device->WriteBuffer(*sbuf, std::as_bytes(std::span<const std::uint32_t>(star_dummy))));
    const BufferHandle bind[] = {*rbuf, *pbuf, *sbuf};
    ASSERT_TRUE(f.device->Dispatch(f.kernel, bind, (w + 7) / 8, (h + 7) / 8, 1).has_value());
    ASSERT_TRUE(f.device->ReadBuffer(*rbuf, std::as_writable_bytes(std::span<float>(vk))));

    // --- CPU reference (the geodesic tracer on the same scene) ---------------
    sirius::core::MorrisThorneCartesian metric(sirius::core::MorrisThorneParams::Ellis(kThroat));

    TracerConfig tc;
    tc.escape_radius = 200.0f;
    tc.horizon_factor = 1.05f;
    tc.max_steps = 20000;
    tc.enable_disk = false;
    tc.integrator.initial_step = 0.1f;
    tc.integrator.max_step = 2.0f;
    tc.integrator.min_step = 1e-5f;
    tc.integrator.abs_tolerance = 5e-6f;
    tc.integrator.rel_tolerance = 5e-6f;
    GeodesicTracer tracer(&metric, tc);

    CameraConfig cc;
    cc.r = kDistance;
    cc.theta = theta;
    cc.phi = 0.0;
    cc.fov = kFovDeg;
    cc.width = w;
    cc.height = h;
    PinholeCamera camera(cc);

    std::vector<float> cpu(static_cast<std::size_t>(w) * h * 4, 0.0f);
    for (int y = 0; y < h; ++y) {
        for (int x = 0; x < w; ++x) {
            const CameraRay ray = camera.GenerateRay(x, y, 0.5f, 0.5f);
            const TraceResult res = tracer.Trace(ray);
            float r = 0.0f, g = 0.0f, b = 0.0f;
            if (res.outcome == TraceResult::Outcome::Escaped) {
                GradientBackground(res.final_direction(1), res.final_direction(2),
                                   res.final_direction(3), r, g, b);
            } else if (res.outcome == TraceResult::Outcome::Spiraling) {
                r = g = b = 0.02f;
            } else if (res.outcome == TraceResult::Outcome::MaxSteps) {
                r = g = b = 0.01f;
            }  // Horizon (throat) -> black
            const int idx = (y * w + x) * 4;
            cpu[idx + 0] = r;
            cpu[idx + 1] = g;
            cpu[idx + 2] = b;
        }
    }

    // --- Statistics -----------------------------------------------------------
    std::size_t shadow_vk = 0, shadow_cpu = 0, bg_pixels = 0;
    double abs_sum = 0.0, max_rel = 0.0;
    constexpr float kShadowLum = 0.02f;
    constexpr float kRelFloor = 0.05f;
    for (int p = 0; p < w * h; ++p) {
        const float lv = Luminance(vk[p * 4], vk[p * 4 + 1], vk[p * 4 + 2]);
        const float lc = Luminance(cpu[p * 4], cpu[p * 4 + 1], cpu[p * 4 + 2]);
        ASSERT_TRUE(std::isfinite(lv) && std::isfinite(lc));
        if (lv < kShadowLum) ++shadow_vk;
        if (lc < kShadowLum) ++shadow_cpu;
        if (lv >= kShadowLum && lc >= kShadowLum) {
            ++bg_pixels;
            const float d = std::abs(lv - lc);
            abs_sum += d;
            max_rel = std::max<double>(max_rel, d / std::max(lc, kRelFloor));
        }
    }

    const double frac_vk = static_cast<double>(shadow_vk) / (w * h);
    const double frac_cpu = static_cast<double>(shadow_cpu) / (w * h);
    const double mean_abs = bg_pixels > 0 ? abs_sum / static_cast<double>(bg_pixels) : 0.0;
    std::cout << "[ mt-parity ] shadow_vk=" << frac_vk << " shadow_cpu=" << frac_cpu
              << " |dfrac|=" << std::abs(frac_vk - frac_cpu) << " bg_pixels=" << bg_pixels
              << " mean|dlum|=" << mean_abs << " max_rel_dlum=" << max_rel << "\n";

    ASSERT_GT(bg_pixels, static_cast<std::size_t>(w * h) / 4)
        << "too few shared background pixels to compare";
    EXPECT_LT(mean_abs, 5e-3) << "mean background luminance difference too large";
    EXPECT_LT(max_rel, 0.2) << "worst-case background luminance disagreement too large";
    EXPECT_GE(frac_vk, 0.005);
    EXPECT_LE(frac_vk, 0.40);
    EXPECT_GE(frac_cpu, 0.005);
    EXPECT_LE(frac_cpu, 0.40);
    EXPECT_LT(std::abs(frac_vk - frac_cpu), 0.02)
        << "throat-shadow fractions diverge between the backends";
}

#endif  // SIRIUS_HAS_VULKAN_BACKEND

// A always-present placeholder so the suite is never empty when the backend is
// compiled out; it documents why the real gates are absent.
TEST(VulkanRenderSession, BackendCompiledOrSkipped) {
#ifndef SIRIUS_HAS_VULKAN_BACKEND
    GTEST_SKIP() << "Vulkan backend not compiled in";
#else
    SUCCEED();
#endif
}

}  // namespace
