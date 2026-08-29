// Vulkan render-path gates (specification programmes 3 and 4). Three concerns:
//   1. a Kerr render through the RenderSession Vulkan path yields finite,
//      non-constant radiance with a bounded horizon shadow (64x64 and 160x120,
//      exercising >1 governed tile), while a budget below the mandatory scene
//      residency declines instead of changing the scene;
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
// background. Disk physics is independently governed by the non-render P4
// kernel-parity probes, so this scene does not duplicate that programme.

#include <gtest/gtest.h>

#include "support/scoped_environment.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <functional>
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
#include "sirius/render/vulkan_renderer.h"
#endif

namespace {

#ifdef SIRIUS_HAS_VULKAN_BACKEND
using sirius::test::ScopedEnvironmentVariable;

using sirius::backend::BufferHandle;
using sirius::backend::BufferUsage;
using sirius::backend::ComputeDevice;
using sirius::backend::CreateVulkanDevice;
using sirius::backend::EnumerateVulkanDevices;
using sirius::backend::GeodesicTracer;
using sirius::backend::ResolveVulkanDeviceIndex;
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

// A scene shared by the kernel-integration gate: Kerr M=1, a=0.9 seen from 30M
// at 80 degrees, with disk and starfield disabled so the geometric comparison
// isolates camera and trajectory behavior.
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
    const double absolute_spin = scene.spin * scene.M;
    const double x = r * st;
    const double y = absolute_spin * st;
    const double cartesian_phi = std::atan2(y, x);
    const double sp = std::sin(cartesian_phi);
    const double cp = std::cos(cartesian_phi);

    std::vector<float> p(72, 0.0f);
    p[46] = 0.5f;
    p[47] = 0.5f;
    p[53] = 1.0f;
    p[0] = static_cast<float>(scene.width);
    p[1] = static_cast<float>(scene.height);
    p[2] = 0.0f;  // Kerr-Schild dispatch
    p[3] = scene.M;
    p[4] = scene.spin * scene.M;  // absolute a
    p[7] = static_cast<float>(x);
    p[8] = static_cast<float>(y);
    p[9] = static_cast<float>(r * ct);
    p[10] = static_cast<float>(-st * cp);
    p[11] = static_cast<float>(-st * sp);
    p[12] = static_cast<float>(-ct);
    p[13] = static_cast<float>(-sp);
    p[14] = static_cast<float>(cp);
    p[15] = 0.0f;
    p[16] = static_cast<float>(ct * cp);
    p[17] = static_cast<float>(ct * sp);
    p[18] = static_cast<float>(-st);
    p[19] = static_cast<float>(scene.fov_deg * kPi / 180.0);
    p[20] = static_cast<float>(scene.width) / static_cast<float>(scene.height);
    p[21] = 3000.0f;  // max_steps
    p[22] = 0.08f;    // stepScale
    p[23] = 0.02f;    // min_step
    p[24] = 2.0f;     // max_step
    p[25] = 200.0f;   // escape_radius
    p[26] = 1.0f;     // exact Kerr-Schild horizon capture surface
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
    const auto selected = ResolveVulkanDeviceIndex(*devices);
    if (!selected.has_value()) return f;
    auto device = CreateVulkanDevice(*selected);
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

struct VulkanScreenPoint {
    double alpha;
    double beta;
};

std::optional<VulkanScreenPoint> BardeenScreenPoint(double photon_radius, double spin,
                                                    double inclination) {
    const double r = photon_radius;
    const double a2 = spin * spin;
    const double xi = (r * r * (r - 3.0) + a2 * (r + 1.0)) / (spin * (1.0 - r));
    const double eta =
        r * r * r * (4.0 * a2 - r * (r - 3.0) * (r - 3.0)) / (a2 * (1.0 - r) * (1.0 - r));
    const double sin_i = std::sin(inclination);
    const double cos_i = std::cos(inclination);
    const double alpha = -xi / sin_i;
    const double beta_squared =
        eta + a2 * cos_i * cos_i - xi * xi * cos_i * cos_i / (sin_i * sin_i);
    if (beta_squared < 0.0) return std::nullopt;
    return VulkanScreenPoint{alpha, std::sqrt(beta_squared)};
}

// --- Session Vulkan path smoke ----------------------------------------------

// Renders a Kerr scene through the RenderSession Vulkan path to an EXR (which
// leaves the display buffer's linear radiance untouched), returning the linear
// RGBA buffer, or an empty vector when the session declines.
std::vector<float> RenderSessionVulkan(
    int width, int height, const std::string& out_path,
    const std::function<void(sirius::render::SessionConfig&)>& configure = {}) {
    const auto devices = EnumerateVulkanDevices();
    if (!devices.has_value() || devices->empty()) return {};

    sirius::render::SessionConfig config;
    config.width = width;
    config.height = height;
    config.samples_per_pixel = 1;
    config.tile_size = 64;
    config.enable_parallel_rendering = false;
    config.backend = sirius::render::RenderBackend::Vulkan;
    config.metric_id = sirius::core::MetricId::Kerr;
    config.black_hole_mass = 1.0;
    config.black_hole_spin = 0.9;
    config.observer_distance = 30.0;
    config.observer_inclination = 80.0 * kPi / 180.0;
    config.camera_fov = 50.0f;
    config.output_path = out_path;
    if (configure) configure(config);

    sirius::render::RenderSession session;
    const auto configured = session.Configure(config);
    if (!configured) {
        ADD_FAILURE() << configured.error().Description();
        return {};
    }
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

void ExpectFiniteResolvedPointField(const std::vector<float>& rgba, int width, int height,
                                    const char* tag) {
    ASSERT_EQ(rgba.size(), static_cast<std::size_t>(width) * height * 4) << tag;
    float minimum_luminance = std::numeric_limits<float>::max();
    float maximum_luminance = 0.0f;
    std::size_t lit_pixels = 0;
    for (int pixel = 0; pixel < width * height; ++pixel) {
        const float r = rgba[pixel * 4 + 0];
        const float g = rgba[pixel * 4 + 1];
        const float b = rgba[pixel * 4 + 2];
        ASSERT_TRUE(std::isfinite(r) && std::isfinite(g) && std::isfinite(b))
            << tag << " non-finite radiance";
        const float luminance = Luminance(r, g, b);
        minimum_luminance = std::min(minimum_luminance, luminance);
        maximum_luminance = std::max(maximum_luminance, luminance);
        if (luminance > 1.0e-4f) ++lit_pixels;
    }
    const double lit_fraction = static_cast<double>(lit_pixels) / (width * height);
    std::cout << "[ " << tag << " ] lit fraction=" << lit_fraction << " luminance range=["
              << minimum_luminance << ", " << maximum_luminance << "]\n";
    EXPECT_GT(maximum_luminance - minimum_luminance, 1.0e-2f)
        << tag << " point radiance does not survive display quantisation";
    EXPECT_GE(lit_fraction, 0.02) << tag << " has too few resolved point-source pixels";
    EXPECT_LE(lit_fraction, 0.50) << tag << " collapsed into a diffuse background";
}

TEST(VulkanRenderSession, CapabilityBoundaryAcceptsRepresentedSceneSemantics) {
    sirius::render::SessionConfig config;
    config.samples_per_pixel = 7;
    config.temperature_model = sirius::render::DiskTemperatureModel::ShakuraSunyaev;
    config.doppler_beaming = false;
    config.camera_beta_forward = 0.1;
    config.lens_type = sirius::core::LensType::ThinLens;
    config.enable_volumetric_disk = true;
    config.enable_turbulence = true;
    config.volumetric_samples = 64;
    EXPECT_TRUE(sirius::render::ValidateVulkanRenderConfig(config).has_value());
}

TEST(VulkanRenderSession, CapabilityBoundaryRejectsUnrepresentedSceneSemantics) {
    sirius::render::SessionConfig config;
    config.enable_volumetric_disk = true;
    config.volumetric_samples = 129;
    EXPECT_FALSE(sirius::render::ValidateVulkanRenderConfig(config).has_value());

    config.volumetric_samples = 64;
    config.point_starfield = true;
    EXPECT_TRUE(sirius::render::ValidateVulkanRenderConfig(config).has_value());

    config.temperature_model = static_cast<sirius::render::DiskTemperatureModel>(255);
    EXPECT_FALSE(sirius::render::ValidateVulkanRenderConfig(config).has_value());
    config.temperature_model = sirius::render::DiskTemperatureModel::NovikovThorne;
    config.lens_type = sirius::core::LensType::Fisheye;
    EXPECT_FALSE(sirius::render::ValidateVulkanRenderConfig(config).has_value());
    config.lens_type = static_cast<sirius::core::LensType>(255);
    EXPECT_FALSE(sirius::render::ValidateVulkanRenderConfig(config).has_value());
    config.lens_type = sirius::core::LensType::Pinhole;
    config.enable_polarisation = true;
    EXPECT_FALSE(sirius::render::ValidateVulkanRenderConfig(config).has_value());
    config.enable_polarisation = false;
    config.metric_id = sirius::core::MetricId::KerrNewman;
    config.enable_disk = true;
    EXPECT_FALSE(sirius::render::ValidateVulkanRenderConfig(config).has_value());
    config.metric_id = sirius::core::MetricId::Minkowski;
    config.black_hole_mass = 0.0;
    EXPECT_FALSE(sirius::render::ValidateVulkanRenderConfig(config).has_value());
    config.metric_id = sirius::core::MetricId::Schwarzschild;
    config.black_hole_mass = 1.0;
    config.black_hole_spin = 0.4;
    EXPECT_FALSE(sirius::render::ValidateVulkanRenderConfig(config).has_value());
    config.black_hole_spin = 0.0;
    config.enable_jets = true;
    EXPECT_FALSE(sirius::render::ValidateVulkanRenderConfig(config).has_value());
    config.enable_jets = false;
    config.enable_corona = true;
    EXPECT_FALSE(sirius::render::ValidateVulkanRenderConfig(config).has_value());
    config.enable_corona = false;
    config.enable_motion_blur = true;
    EXPECT_FALSE(sirius::render::ValidateVulkanRenderConfig(config).has_value());
    config.enable_motion_blur = false;
    config.color_mode = sirius::core::color_modes::Mode::RedshiftMap;
    EXPECT_FALSE(sirius::render::ValidateVulkanRenderConfig(config).has_value());
}

TEST(VulkanRenderSession, ProceduralVolumetricTurbulenceReachesLiveKernel) {
    if (const auto devices = EnumerateVulkanDevices(); !devices || devices->empty()) {
        GTEST_SKIP() << "no Vulkan device present";
    }
    const std::string root = std::string(std::getenv("TMPDIR") ? std::getenv("TMPDIR") : "/tmp");
    const auto baseline = RenderSessionVulkan(40, 24, root + "/sirius_vk_volume_baseline.exr",
                                              [](sirius::render::SessionConfig& config) {
                                                  config.enable_volumetric_disk = true;
                                                  config.volumetric_samples = 4;
                                              });
    const auto represented = RenderSessionVulkan(40, 24, root + "/sirius_vk_volume_turbulence.exr",
                                                 [](sirius::render::SessionConfig& config) {
                                                     config.enable_volumetric_disk = true;
                                                     config.enable_turbulence = true;
                                                     config.volumetric_samples = 4;
                                                 });
    ASSERT_FALSE(baseline.empty());
    ASSERT_EQ(represented.size(), baseline.size());
    double absoluteDifference = 0.0;
    for (std::size_t i = 0; i < represented.size(); ++i) {
        ASSERT_TRUE(std::isfinite(represented[i]));
        absoluteDifference += std::abs(static_cast<double>(represented[i] - baseline[i]));
    }
    EXPECT_GT(absoluteDifference, 1.0e-5)
        << "Vulkan procedural turbulence did not affect live volumetric transfer";
}

TEST(VulkanRenderSession, ThinAndVolumetricDopplerSuppressionAffectLiveEmission) {
    if (const auto devices = EnumerateVulkanDevices(); !devices || devices->empty()) {
        GTEST_SKIP() << "no Vulkan device present";
    }
    const std::string root = std::string(std::getenv("TMPDIR") ? std::getenv("TMPDIR") : "/tmp");
    auto render = [&](const char* suffix, bool volumetric, bool doppler) {
        return RenderSessionVulkan(40, 24, root + "/sirius_vk_doppler_" + suffix + ".exr",
                                   [=](sirius::render::SessionConfig& config) {
                                       config.enable_volumetric_disk = volumetric;
                                       config.volumetric_samples = 4;
                                       config.doppler_beaming = doppler;
                                   });
    };
    const auto thin_beamed = render("thin_on", false, true);
    const auto thin_gravity = render("thin_off", false, false);
    const auto volume_beamed = render("volume_on", true, true);
    const auto volume_gravity = render("volume_off", true, false);

    auto absolute_difference = [](const std::vector<float>& lhs, const std::vector<float>& rhs) {
        EXPECT_FALSE(lhs.empty());
        EXPECT_EQ(lhs.size(), rhs.size());
        double difference = 0.0;
        for (std::size_t i = 0; i < std::min(lhs.size(), rhs.size()); ++i) {
            EXPECT_TRUE(std::isfinite(lhs[i]));
            EXPECT_TRUE(std::isfinite(rhs[i]));
            difference += std::abs(static_cast<double>(lhs[i] - rhs[i]));
        }
        return difference;
    };

    EXPECT_GT(absolute_difference(thin_beamed, thin_gravity), 1.0e-5)
        << "the Doppler toggle did not reach live Vulkan thin-disk emission";
    EXPECT_GT(absolute_difference(volume_beamed, volume_gravity), 1.0e-5)
        << "the Doppler toggle did not reach live Vulkan volumetric emission";
}

TEST(VulkanRenderSession, NonSquareMultisamplingCameraAndLensReachLiveKernel) {
    if (const auto d = EnumerateVulkanDevices(); !d.has_value() || d->empty()) {
        GTEST_SKIP() << "no Vulkan device present";
    }
    const std::string root = std::string(std::getenv("TMPDIR") ? std::getenv("TMPDIR") : "/tmp");
    const auto baseline = RenderSessionVulkan(48, 32, root + "/sirius_vk_semantics_baseline.exr");
    const auto represented = RenderSessionVulkan(
        48, 32, root + "/sirius_vk_semantics_represented.exr",
        [](sirius::render::SessionConfig& config) {
            config.samples_per_pixel = 3;
            config.camera_beta_forward = 0.1;
            config.lens_type = sirius::core::LensType::ThinLens;
            config.camera_aperture = 2.8f;
            config.camera_focal_length = 50.0f;
            config.camera_focus_distance = 30.0f;
            config.temperature_model = sirius::render::DiskTemperatureModel::ShakuraSunyaev;
            config.doppler_beaming = false;
        });
    ASSERT_FALSE(baseline.empty());
    ASSERT_EQ(represented.size(), baseline.size());
    double absolute_difference = 0.0;
    for (std::size_t i = 0; i < represented.size(); ++i) {
        ASSERT_TRUE(std::isfinite(represented[i]));
        absolute_difference += std::abs(static_cast<double>(represented[i] - baseline[i]));
    }
    EXPECT_GT(absolute_difference, 1.0e-4)
        << "represented Vulkan sampling/camera/lens controls did not affect the live output";
}

TEST(VulkanRenderSession, IndexedPointCatalogueReachesLiveKernel) {
    if (const auto devices = EnumerateVulkanDevices(); !devices || devices->empty()) {
        GTEST_SKIP() << "no Vulkan device present";
    }
    const std::string root = std::string(std::getenv("TMPDIR") ? std::getenv("TMPDIR") : "/tmp");
    const auto texture = RenderSessionVulkan(
        48, 32, root + "/sirius_vk_texture_stars.exr",
        [](sirius::render::SessionConfig& config) { config.enable_disk = false; });
    const auto points = RenderSessionVulkan(48, 32, root + "/sirius_vk_point_stars.exr",
                                            [](sirius::render::SessionConfig& config) {
                                                config.point_starfield = true;
                                                config.ray_bundles = true;
                                                config.enable_disk = false;
                                                config.starfield_config.star_count = 100000;
                                                config.starfield_config.brightness_scale = 100.0f;
                                            });
    ASSERT_FALSE(texture.empty());
    ASSERT_EQ(points.size(), texture.size());
    double point_energy = 0.0;
    double absolute_difference = 0.0;
    double minimum_luminance = std::numeric_limits<double>::infinity();
    double maximum_luminance = 0.0;
    for (std::size_t i = 0; i < points.size(); i += 4) {
        ASSERT_TRUE(std::isfinite(points[i]) && std::isfinite(points[i + 1]) &&
                    std::isfinite(points[i + 2]));
        point_energy += points[i] + points[i + 1] + points[i + 2];
        const double luminance = Luminance(points[i], points[i + 1], points[i + 2]);
        minimum_luminance = std::min(minimum_luminance, luminance);
        maximum_luminance = std::max(maximum_luminance, luminance);
        absolute_difference += std::abs(static_cast<double>(points[i] - texture[i])) +
                               std::abs(static_cast<double>(points[i + 1] - texture[i + 1])) +
                               std::abs(static_cast<double>(points[i + 2] - texture[i + 2]));
    }
    EXPECT_GT(point_energy, 1.0e-2) << "indexed point catalogue produced no visible radiance";
    EXPECT_GT(maximum_luminance - minimum_luminance, 1.0e-4)
        << "indexed point catalogue collapsed to a display-invisible constant field";
    EXPECT_GT(absolute_difference, 1.0e-3)
        << "point-catalogue request appears to have sampled the texture path";
}

TEST(VulkanRenderSession, CombinedParitySceneRetainsResolvedImageStructure) {
    if (const auto devices = EnumerateVulkanDevices(); !devices || devices->empty()) {
        GTEST_SKIP() << "no Vulkan device present";
    }
    const std::string root = std::string(std::getenv("TMPDIR") ? std::getenv("TMPDIR") : "/tmp");
    constexpr int kWidth = 192;
    constexpr int kHeight = 128;
    const auto represented =
        RenderSessionVulkan(kWidth, kHeight, root + "/sirius_vk_combined_p3_p5.exr",
                            [](sirius::render::SessionConfig& config) {
                                config.samples_per_pixel = 3;
                                config.camera_fov = 60.0f;
                                config.camera_beta_forward = 0.1;
                                config.camera_beta_up = 0.02;
                                config.camera_beta_right = -0.01;
                                config.lens_type = sirius::core::LensType::ThinLens;
                                config.camera_focal_length = 50.0f;
                                config.camera_aperture = 2.8f;
                                config.camera_focus_distance = 30.0f;
                                config.point_starfield = true;
                                config.ray_bundles = true;
                                config.enable_disk = false;
                                config.starfield_config.star_count = 100000;
                            });
    ASSERT_FALSE(represented.empty());
    ExpectFiniteResolvedPointField(represented, kWidth, kHeight, "combined P3/P5 scene");
}

TEST(VulkanRenderSession, CpuVulkanPointCatalogueAgreeOnFlatScene) {
    if (const auto devices = EnumerateVulkanDevices(); !devices || devices->empty()) {
        GTEST_SKIP() << "no Vulkan device present";
    }
    const std::string root = std::string(std::getenv("TMPDIR") ? std::getenv("TMPDIR") : "/tmp");
    auto configure_points = [](sirius::render::SessionConfig& config) {
        config.metric_id = sirius::core::MetricId::Minkowski;
        config.black_hole_mass = 0.0;
        config.black_hole_spin = 0.0;
        config.enable_disk = false;
        config.observer_distance = 20.0;
        config.observer_inclination = 1.1;
        config.camera_fov = 30.0f;
        config.point_starfield = true;
        config.ray_bundles = false;
        config.starfield_config.star_count = 100000;
        config.starfield_config.brightness_scale = 1.0f;
    };
    const auto gpu =
        RenderSessionVulkan(32, 20, root + "/sirius_vk_point_flat.exr", configure_points);
    const auto cpu = RenderSessionVulkan(32, 20, root + "/sirius_cpu_point_flat.exr",
                                         [&](sirius::render::SessionConfig& config) {
                                             configure_points(config);
                                             config.backend = sirius::render::RenderBackend::Cpu;
                                         });
    ASSERT_EQ(gpu.size(), cpu.size());
    ASSERT_FALSE(gpu.empty());
    double absolute_sum = 0.0;
    double signal_sum = 0.0;
    double maximum = 0.0;
    std::size_t channels = 0;
    for (std::size_t i = 0; i < gpu.size(); i += 4) {
        for (int channel = 0; channel < 3; ++channel) {
            const double difference = std::abs(static_cast<double>(gpu[i + channel]) -
                                               static_cast<double>(cpu[i + channel]));
            absolute_sum += difference;
            maximum = std::max(maximum, difference);
            signal_sum += std::abs(static_cast<double>(cpu[i + channel]));
            ++channels;
        }
    }
    const double mean = absolute_sum / static_cast<double>(channels);
    std::cout << "[ point-star parity ] mean|d|=" << mean << " max|d|=" << maximum
              << " relative L1=" << absolute_sum / std::max(signal_sum, 1.0e-20) << "\n";
    EXPECT_GT(signal_sum, 1.0e-5);
    EXPECT_LT(mean, 2.0e-5);
    EXPECT_LT(absolute_sum / signal_sum, 0.02);
}

TEST(VulkanRenderSession, ConstrainedBudgetDeclinesRatherThanChangingBackground) {
    if (const auto d = EnumerateVulkanDevices(); !d.has_value() || d->empty()) {
        GTEST_SKIP() << "no Vulkan device present";
    }
    ScopedEnvironmentVariable budget("SIRIUS_MEMORY_BUDGET_MB",
                                     "1");  // Less than the mandatory starfield residency.
    const std::string out = std::string(std::getenv("TMPDIR") ? std::getenv("TMPDIR") : "/tmp") +
                            "/sirius_vk_smoke_64.exr";
    const auto rgba = RenderSessionVulkan(64, 64, out);
    EXPECT_TRUE(rgba.empty())
        << "the Vulkan path must not substitute an analytic background under pressure";
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
    const auto selected = ResolveVulkanDeviceIndex(*devices);
    ASSERT_TRUE(selected.has_value()) << selected.error().Description();
    auto device = CreateVulkanDevice(*selected);
    ASSERT_TRUE(device.has_value());
    const bool fp64_supported = (*device)->Info().supports_fp64;
    device->reset();

    ScopedEnvironmentVariable precision("SIRIUS_PRECISION", "fp64");
    const std::string out = std::string(std::getenv("TMPDIR") ? std::getenv("TMPDIR") : "/tmp") +
                            "/sirius_vk_fp64_64.exr";
    const auto rgba = RenderSessionVulkan(64, 64, out);

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
    ScopedEnvironmentVariable precision("SIRIUS_PRECISION", "fp32-comp");
    const std::string out = std::string(std::getenv("TMPDIR") ? std::getenv("TMPDIR") : "/tmp") +
                            "/sirius_vk_comp_64.exr";
    const auto rgba = RenderSessionVulkan(64, 64, out);
    ASSERT_FALSE(rgba.empty()) << "compensated rung did not complete";
    ExpectFiniteNonConstantWithShadow(rgba, 64, 64, "vk-session-comp");
}

TEST(VulkanRenderSession, Kerr160x120CompletesAcrossMultipleGovernedTiles) {
    if (const auto d = EnumerateVulkanDevices(); !d.has_value() || d->empty()) {
        GTEST_SKIP() << "no Vulkan device present";
    }
    const std::string out = std::string(std::getenv("TMPDIR") ? std::getenv("TMPDIR") : "/tmp") +
                            "/sirius_vk_smoke_160.exr";
    const auto rgba = RenderSessionVulkan(160, 120, out);
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
    const auto psbuf = f.device->CreateBuffer(sizeof(std::uint32_t), BufferUsage::kStorage);
    const auto pobuf = f.device->CreateBuffer(sizeof(std::uint32_t), BufferUsage::kStorage);
    const auto pibuf = f.device->CreateBuffer(sizeof(std::uint32_t), BufferUsage::kStorage);
    ASSERT_TRUE(rbuf && pbuf && sbuf && psbuf && pobuf && pibuf);
    ASSERT_TRUE(f.device->WriteBuffer(*pbuf, std::as_bytes(std::span<const float>(params))));
    ASSERT_TRUE(
        f.device->WriteBuffer(*sbuf, std::as_bytes(std::span<const std::uint32_t>(star_dummy))));
    const BufferHandle bind[] = {*rbuf, *pbuf, *sbuf, *psbuf, *pobuf, *pibuf};
    ASSERT_TRUE(f.device->Dispatch(f.kernel, bind, (w + 7) / 8, (h + 7) / 8, 1).has_value());
    ASSERT_TRUE(f.device->ReadBuffer(*rbuf, std::as_writable_bytes(std::span<float>(vk))));

    // --- CPU reference (the geodesic tracer on the same scene) ---------------
    KerrSchildParams kp;
    kp.M = scene.M;
    kp.a = scene.spin * scene.M;
    KerrSchildFamily metric(kp);

    TracerConfig tc;
    tc.escape_radius = 200.0f;
    tc.horizon_factor = 1.0f;
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
            }  // Horizon or numerical work-bound outcome -> black
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

TEST(VulkanRenderSession, KerrNearExtremalBardeenBoundaryAt1080p) {
    KernelFixture fixture = OpenKernel();
    if (!fixture.ready) GTEST_SKIP() << "no Vulkan device or trace kernel absent";

    Scene scene;
    scene.width = 1920;
    scene.height = 1080;
    scene.spin = 0.998f;
    scene.distance = 50.0;
    scene.inclination_deg = 60.0;
    scene.fov_deg = 100.0f;
    std::vector<float> params = BuildTraceParams(scene);
    params[21] = 30000.0f;
    params[22] = 0.02f;
    params[24] = 0.25f;
    params[34] = 1.0f;
    params[35] = 1.0f;

    std::array<float, 4> radiance{};
    const std::array<std::uint32_t, 1> star_dummy{0u};
    const auto radiance_buffer =
        fixture.device->CreateBuffer(sizeof(radiance), BufferUsage::kStorage);
    const auto params_buffer =
        fixture.device->CreateBuffer(params.size() * sizeof(float), BufferUsage::kStorage);
    const auto star_buffer =
        fixture.device->CreateBuffer(sizeof(star_dummy), BufferUsage::kStorage);
    const auto point_star_buffer =
        fixture.device->CreateBuffer(sizeof(std::uint32_t), BufferUsage::kStorage);
    const auto point_offset_buffer =
        fixture.device->CreateBuffer(sizeof(std::uint32_t), BufferUsage::kStorage);
    const auto point_index_buffer =
        fixture.device->CreateBuffer(sizeof(std::uint32_t), BufferUsage::kStorage);
    ASSERT_TRUE(radiance_buffer && params_buffer && star_buffer && point_star_buffer &&
                point_offset_buffer && point_index_buffer);
    ASSERT_TRUE(fixture.device->WriteBuffer(
        *star_buffer, std::as_bytes(std::span<const std::uint32_t>(star_dummy))));
    const BufferHandle bindings[] = {*radiance_buffer,   *params_buffer,       *star_buffer,
                                     *point_star_buffer, *point_offset_buffer, *point_index_buffer};

    constexpr double kInclination = 60.0 * kPi / 180.0;
    constexpr VulkanScreenPoint kAnalyticCentre{-2.1573218480479185, 0.0};
    const double tan_half_fov = std::tan(scene.fov_deg * kPi / 360.0);
    const double aspect = static_cast<double>(scene.width) / scene.height;
    const double pixels_per_screen_unit = 0.5 * scene.height / (scene.distance * tan_half_fov);

    auto is_captured = [&](const VulkanScreenPoint& point) {
        const double image_x =
            0.5 * scene.width * (1.0 + point.alpha / (scene.distance * tan_half_fov * aspect));
        const double image_y =
            0.5 * scene.height * (1.0 - point.beta / (scene.distance * tan_half_fov));
        const int pixel_x = static_cast<int>(std::floor(image_x));
        const int pixel_y = static_cast<int>(std::floor(image_y));
        params[32] = static_cast<float>(pixel_x);
        params[33] = static_cast<float>(pixel_y);
        params[46] = static_cast<float>(image_x - pixel_x);
        params[47] = static_cast<float>(image_y - pixel_y);
        radiance.fill(0.0f);
        EXPECT_TRUE(fixture.device->WriteBuffer(*radiance_buffer,
                                                std::as_bytes(std::span<const float>(radiance))));
        EXPECT_TRUE(fixture.device->WriteBuffer(*params_buffer,
                                                std::as_bytes(std::span<const float>(params))));
        EXPECT_TRUE(fixture.device->Dispatch(fixture.kernel, bindings, 1, 1, 1));
        EXPECT_TRUE(fixture.device->ReadBuffer(*radiance_buffer,
                                               std::as_writable_bytes(std::span<float>(radiance))));
        return std::max({radiance[0], radiance[1], radiance[2]}) < 1.0e-6f;
    };

    // Same full-curve sample set as the independent CPU classifier.
    for (const double photon_radius : {1.12, 1.2, 1.35, 1.5, 1.8, 2.1, 2.5, 3.0, 3.5, 3.7}) {
        const auto analytic = BardeenScreenPoint(photon_radius, scene.spin, kInclination);
        ASSERT_TRUE(analytic.has_value());
        const VulkanScreenPoint camera_convention{-analytic->alpha, analytic->beta};
        const VulkanScreenPoint delta{camera_convention.alpha - kAnalyticCentre.alpha,
                                      camera_convention.beta - kAnalyticCentre.beta};
        auto scaled = [&](double scale) {
            return VulkanScreenPoint{kAnalyticCentre.alpha + scale * delta.alpha,
                                     kAnalyticCentre.beta + scale * delta.beta};
        };

        double inside = 0.70;
        double outside = 1.30;
        ASSERT_TRUE(is_captured(scaled(inside)));
        ASSERT_FALSE(is_captured(scaled(outside)));
        for (int iteration = 0; iteration < 14; ++iteration) {
            const double middle = 0.5 * (inside + outside);
            if (is_captured(scaled(middle))) {
                inside = middle;
            } else {
                outside = middle;
            }
        }
        const VulkanScreenPoint measured = scaled(0.5 * (inside + outside));
        const double displacement = std::hypot(measured.alpha - camera_convention.alpha,
                                               measured.beta - camera_convention.beta) *
                                    pixels_per_screen_unit;
        EXPECT_LT(displacement, 1.0) << "Vulkan photon orbit r/M=" << photon_radius
                                     << ", scale bracket=[" << inside << ", " << outside << "]";
    }
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
    std::vector<float> params(72, 0.0f);
    params[46] = 0.5f;
    params[47] = 0.5f;
    params[53] = 1.0f;
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
    const auto psbuf = f.device->CreateBuffer(sizeof(std::uint32_t), BufferUsage::kStorage);
    const auto pobuf = f.device->CreateBuffer(sizeof(std::uint32_t), BufferUsage::kStorage);
    const auto pibuf = f.device->CreateBuffer(sizeof(std::uint32_t), BufferUsage::kStorage);
    ASSERT_TRUE(rbuf && pbuf && sbuf && psbuf && pobuf && pibuf);
    ASSERT_TRUE(f.device->WriteBuffer(*pbuf, std::as_bytes(std::span<const float>(params))));
    ASSERT_TRUE(
        f.device->WriteBuffer(*sbuf, std::as_bytes(std::span<const std::uint32_t>(star_dummy))));
    const BufferHandle bind[] = {*rbuf, *pbuf, *sbuf, *psbuf, *pobuf, *pibuf};
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
            }  // Horizon (throat) or numerical work-bound outcome -> black
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

}  // namespace
