#include "support/moving_kerr_detector_scene.h"
#include "support/point_source_band_reference.h"
#include "support/scoped_environment.h"
#include "support/separated_geodesic_reference.h"
// End-to-end CPU render-session probe (session-level gate).
//
// Drives RenderSession CPU-only at 64x64, 4 spp, Kerr a=0.9, writing a PNG and
// an EXR to a temp directory. Asserts: the session reaches Complete, both files
// exist, the PNG decodes (stb) and the EXR loads (tinyexr), and each decoded
// image is finite and non-constant.

#include "sirius/core/camera_sampling.h"
#include "sirius/core/metrics/kerr_schild_family.h"
#include "sirius/render/session/render_session.h"
#include "sirius/render/trace_domain.h"
#ifdef SIRIUS_HAS_VULKAN_BACKEND
#include "sirius/backend/device.h"
#include "sirius/render/vulkan_renderer.h"
#endif

#include <gtest/gtest.h>

#include "support/scoped_temporary_directory.h"
#include <nlohmann/json.hpp>
#include <stb_image.h>
#include <tinyexr.h>

#include <array>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <format>
#include <fstream>
#include <future>
#include <iostream>
#include <limits>
#include <numbers>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using sirius::render::RenderSession;
using sirius::render::SessionConfig;
using sirius::render::SessionState;
using sirius::test::ScopedTemporaryDirectory;

SessionConfig ProbeConfig(const std::string& output_path) {
    SessionConfig cfg;
    cfg.width = 64;
    cfg.height = 64;
    cfg.samples_per_pixel = 4;
    cfg.tile_size = 64;
    cfg.enable_parallel_rendering = false;  // Deterministic single-thread for the probe.
    cfg.metric_id = sirius::core::MetricId::Kerr;
    cfg.black_hole_mass = 1.0;
    cfg.black_hole_spin = 0.9;
    cfg.output_path = output_path;
    return cfg;
}

// Run one CPU render to output_path; returns the terminal session state.
SessionState RenderTo(const std::string& output_path) {
    auto config = ProbeConfig(output_path);
    // Writer coverage uses the complete physical frame. Independent tiles can
    // use the normal worker pool; serial/parallel equality has its own probe.
    config.tile_size = 8;
    config.enable_parallel_rendering = true;
    RenderSession session;
    if (!session.Configure(config)) {
        return SessionState::Failed;
    }
    return session.Execute();
}

}  // namespace

namespace {
namespace disk_reference = sirius::test::separated_reference;
namespace disk_band_reference = sirius::test::point_source_band_reference;

SessionConfig DiskCompositionConfig(double spin) {
    SessionConfig config;
    config.width = config.height = 3;
    config.samples_per_pixel = 1;
    config.tile_size = 4;
    config.enable_parallel_rendering = false;
    config.write_output = false;
    config.metric_id =
        spin == 0 ? sirius::core::MetricId::Schwarzschild : sirius::core::MetricId::Kerr;
    config.black_hole_mass = 1;
    config.black_hole_spin = spin;
    config.observer_distance = 50;
    config.observer_inclination = std::numbers::pi / 3;
    config.observer_azimuth = .31;
    config.camera_fov = 32;
    if (spin != 0) {
        config.camera_beta_forward = .03;
        config.camera_beta_up = -.01;
        config.camera_beta_right = .02;
    }
    config.temperature_model = sirius::render::DiskTemperatureModel::ShakuraSunyaev;
    config.disk_temperature_scale = 12000;
    config.enable_disk = true;
    config.doppler_beaming = true;
    config.enable_bloom = false;
    config.enable_film_finish = false;
    config.tonemapper = sirius::core::TonemapType::None;
    config.exposure = config.contrast = config.saturation = 1;
    return config;
}

// The declared model is a Newtonian zero-torque radial flux profile carried
// by a relativistic circular emitter. This test does not call the production
// ISCO, temperature, frequency-transfer or colour helpers for its expectation.
long double IndependentDiskInner(double spin) {
    const long double a = spin;
    const long double z1 = 1 + std::cbrt(1 - a * a) * (std::cbrt(1 + a) + std::cbrt(1 - a));
    const long double z2 = std::sqrt(3 * a * a + z1 * z1);
    return 3 + z2 - (a >= 0 ? 1 : -1) * std::sqrt((3 - z1) * (3 + z1 + 2 * z2));
}

struct DiskCompositionExpected {
    disk_reference::Result event;
    long double radius, temperature, frequency_ratio;
    disk_band_reference::Channels rgb;
};

DiskCompositionExpected IndependentDiskPixel(const SessionConfig& config,
                                             const sirius::core::CameraLaunch& launch,
                                             long double step_fraction) {
    const long double mass = config.black_hole_mass;
    const long double spin = config.black_hole_spin * mass;
    const long double inner = IndependentDiskInner(config.black_hole_spin) * mass;
    auto event = disk_reference::Trace(launch, static_cast<double>(mass), static_cast<double>(spin),
                                       200 * static_cast<double>(mass), step_fraction,
                                       static_cast<double>(inner), 20 * static_cast<double>(mass));
    if (event.fate != disk_reference::Fate::Disk)
        throw std::runtime_error("declared disk composition witness missed its first disk event");
    const auto geometry = disk_reference::At(event.x, mass, spin);
    const long double r = geometry.radius;
    const long double omega = std::sqrt(mass) / (std::pow(r, 1.5L) + spin * std::sqrt(mass));
    // Fixed-r circular coordinate velocity agrees in BL and Cartesian KS.
    const long double gtt = -1 + 2 * mass / r;
    const long double gtphi = -2 * mass * spin / r;
    const long double gphiphi = r * r + spin * spin + 2 * mass * spin * spin / r;
    const long double norm = -(gtt + 2 * omega * gtphi + omega * omega * gphiphi);
    if (!(r > inner && r < 20 * mass && norm > 0))
        throw std::runtime_error("disk witness outside stable timelike circular-emitter domain");
    // The reference carries the past-directed tangent: E=-p_t, L=p_phi.
    // The physical future photon has constants (-E,-L), camera frequency one.
    const long double emitted_frequency =
        (-event.energy + omega * event.angular_momentum) / std::sqrt(norm);
    if (!(emitted_frequency > 0)) throw std::runtime_error("invalid independent disk frequency");
    const long double g = 1 / emitted_frequency;
    const auto flux = [inner](long double radius) {
        const long double q = inner / radius;
        return q * q * q * (1 - std::sqrt(q));
    };
    const long double temperature = std::pow(flux(r) / flux(1.5L * inner), .25L);
    // The represented thin-disk colour is relative hue (brightest channel one)
    // times observed bolometric intensity, not an absolute visible-band flux.
    auto rgb = disk_band_reference::Rgb(disk_band_reference::Integrate(
        temperature * g * config.disk_temperature_scale, 1, 32, false, false));
    const long double normalizer = std::max({rgb[0], rgb[1], rgb[2], .001L});
    for (auto& channel : rgb) channel = channel / normalizer * std::pow(temperature * g, 4);
    return {std::move(event), r, temperature, g, rgb};
}

struct DiskCompositionWitness {
    int x, y;
    DiskCompositionExpected expected;
    long double reference_gap;
    long double event_gap;
    long double frequency_gap;
};

std::vector<DiskCompositionWitness> DiskCompositionWitnesses(const SessionConfig& config) {
    sirius::core::CameraConfig camera_config;
    camera_config.r = config.observer_distance;
    camera_config.theta = config.observer_inclination;
    camera_config.phi = config.observer_azimuth;
    camera_config.fov = config.camera_fov;
    camera_config.width = config.width;
    camera_config.height = config.height;
    camera_config.beta_x = config.camera_beta_forward;
    camera_config.beta_y = config.camera_beta_up;
    camera_config.beta_z = config.camera_beta_right;
    sirius::core::PinholeCamera camera(camera_config);
    sirius::core::KerrSchildFamily metric(sirius::core::KerrSchildParams::Kerr(
        config.black_hole_mass, config.black_hole_spin * config.black_hole_mass));
    std::vector<DiskCompositionWitness> witnesses;
    for (const int x : {0, 2}) {
        // Public camera projections and LaunchCameraRay are measured boundary
        // inputs. Their independent analytic qualification is P1; all complete
        // trajectories and all emitter/spectral physics below are independent.
        const auto film = camera.ProjectFilmForObserver(x + .5, 2.5, .2f, 1.f / 7);
        if (!film) throw std::runtime_error("disk witness camera projection declined");
        const auto launch = sirius::core::LaunchCameraRay(
            metric, config.black_hole_spin * config.black_hole_mass, film->ray);
        if (!launch) throw std::runtime_error("disk witness measured camera launch declined");
        auto coarse = IndependentDiskPixel(config, *launch, .002L);
        auto fine = IndependentDiskPixel(config, *launch, .001L);
        const long double peak = std::max({fine.rgb[0], fine.rgb[1], fine.rgb[2]});
        long double gap = 0;
        for (int channel = 0; channel < 3; ++channel)
            gap = std::max(gap, std::abs(coarse.rgb[channel] - fine.rgb[channel]) / peak);
        long double event_gap = 0;
        for (int mu = 0; mu < 4; ++mu)
            event_gap = std::max(event_gap, std::abs(coarse.event.x[mu] - fine.event.x[mu]) /
                                                config.black_hole_mass);
        const long double frequency_gap =
            std::abs(coarse.frequency_ratio - fine.frequency_ratio) / fine.frequency_ratio;
        witnesses.push_back({x, 2, std::move(fine), gap, event_gap, frequency_gap});
    }
    return witnesses;
}

void CheckDiskComposition(const SessionConfig& config, const std::vector<float>& pixels,
                          const std::vector<DiskCompositionWitness>& witnesses,
                          const std::string& backend) {
    ASSERT_EQ(pixels.size(), static_cast<std::size_t>(config.width * config.height * 4));
    for (const auto& witness : witnesses) {
        const auto& expected = witness.expected;
        SCOPED_TRACE(backend + " spin=" + std::to_string(config.black_hole_spin) +
                     " pixel=" + std::to_string(witness.x) + "," + std::to_string(witness.y));
        ASSERT_EQ(expected.event.fate, disk_reference::Fate::Disk);
        ASSERT_LT(witness.reference_gap, 1e-6L);
        ASSERT_LT(witness.event_gap, 1e-6L);
        ASSERT_LT(witness.frequency_gap, 1e-6L);
        ASSERT_LE(std::abs(expected.event.x[3]) / config.black_hole_mass, 1e-12L);
        ASSERT_EQ(expected.event.radial_turns, 0u);
        // These are transverse first crossings, away from either disk edge;
        // this finite witness makes no near-edge or caustic conditioning claim.
        ASSERT_GT(
            expected.radius / config.black_hole_mass - IndependentDiskInner(config.black_hole_spin),
            1);
        ASSERT_GT(20 - expected.radius / config.black_hole_mass, 1);
        ASSERT_GT(std::abs(expected.event.k[3]), .1L);
        const long double peak = std::max({expected.rgb[0], expected.rgb[1], expected.rgb[2]});
        ASSERT_GT(peak, 1e-3L);
        const std::size_t index = (witness.y * config.width + witness.x) * 4;
        long double error = 0;
        for (int channel = 0; channel < 3; ++channel) {
            ASSERT_TRUE(std::isfinite(pixels[index + channel]));
            ASSERT_GE(pixels[index + channel], 0);
            error =
                std::max(error, std::abs(pixels[index + channel] - expected.rgb[channel]) / peak);
        }
        // Fixed session controller abs/rel=5e-6, float film inputs and channels;
        // use the existing central-trajectory 1e-4 envelope without fitting it
        // to the observed result. The oracle spends at most 1% of that budget.
        EXPECT_LE(error + witness.reference_gap, 1e-4L);
        EXPECT_EQ(pixels[index + 3], 1);
        const auto key = backend + "_a" + std::to_string(config.black_hole_spin) + "_x" +
                         std::to_string(witness.x);
        ::testing::Test::RecordProperty(
            key, std::format("radius/M={:.12g},temperature/profile={:.12g},g={:.12g},"
                             "rgb=[{:.12g},{:.12g},{:.12g}],relative_channel_error={:.12g},"
                             "reference_relative_gap={:.12g},reference_endpoint_gap/M={:.12g},"
                             "reference_frequency_gap={:.12g},reference_steps={}",
                             expected.radius / config.black_hole_mass, expected.temperature,
                             expected.frequency_ratio, expected.rgb[0], expected.rgb[1],
                             expected.rgb[2], error, witness.reference_gap, witness.event_gap,
                             witness.frequency_gap, expected.event.steps));
    }
}
}  // namespace

TEST(RenderSessionProbe, CpuThinDiskPublishedLinearChannelsMatchIndependentFirstEventPhysics) {
    RecordProperty("scope",
                   "typed CPU RenderSession; published pre-grade tile radiance; "
                   "first opaque disk Carter event; moving Kerr observer; no output");
    for (const double spin : {0., .7}) {
        auto config = DiskCompositionConfig(spin);
        const auto witnesses = DiskCompositionWitnesses(config);
        RenderSession session;
        const auto configured = session.Configure(config);
        ASSERT_TRUE(configured) << configured.error().Description();
        std::vector<float> linear;
        // CPU UpdateTile publishes complete physical radiance before reporting
        // completion. The final in-memory display later receives fixed shadow
        // lift/clipping even when output and tonemapping are disabled.
        session.SetProgressCallback([&](float, int done, int total, double) {
            if (done == total) linear = session.GetDisplayBuffer().SnapshotFloatData();
        });
        ASSERT_EQ(session.Execute(), SessionState::Complete) << session.GetErrorMessage();
        ASSERT_NO_FATAL_FAILURE(CheckDiskComposition(config, linear, witnesses, "cpu"));
        EXPECT_EQ(session.GetTileScheduler().GetCompletedCount(), 1);
    }
}

TEST(RenderSessionProbe, VulkanThinDiskPublishedLinearChannelsMatchIndependentFirstEventPhysics) {
#ifdef SIRIUS_HAS_VULKAN_BACKEND
    const auto devices = sirius::backend::EnumerateVulkanDevices();
    ASSERT_TRUE(devices) << devices.error().Description();
    if (devices->empty()) GTEST_SKIP() << "Vulkan unavailable; disk session unqualified";
    const sirius::test::ScopedEnvironmentVariable precision("SIRIUS_PRECISION", "fp32");
    RecordProperty("scope",
                   "production RenderVulkanToDisplay -> retained executor -> internal "
                   "RenderSession shared shading -> linear publication; no output");
    for (const double spin : {0., .7}) {
        auto config = DiskCompositionConfig(spin);
        config.backend = sirius::render::RenderBackend::Vulkan;
        // Device batching accepts the Vulkan defaults, rather than CPU tile
        // and thread controls from the shared composition fixture.
        const SessionConfig defaults;
        config.tile_size = defaults.tile_size;
        config.enable_parallel_rendering = defaults.enable_parallel_rendering;
        const auto witnesses = DiskCompositionWitnesses(config);
        sirius::render::DisplayBuffer display;
        display.Initialise(config.width, config.height);
        const auto rendered = sirius::render::RenderVulkanToDisplay(config, display);
        ASSERT_TRUE(rendered) << rendered.error().Description();
        ASSERT_TRUE(rendered->retained_intervals);
        EXPECT_GT(rendered->camera_batches, 0u);
        EXPECT_GT(rendered->accepted_intervals, 0u);
        EXPECT_EQ(rendered->precision, sirius::render::PrecisionRung::Fp32);
        RecordProperty("device", rendered->device_name);
        ASSERT_NO_FATAL_FAILURE(
            CheckDiskComposition(config, display.SnapshotFloatData(), witnesses, "vulkan_fp32"));
    }
#else
    GTEST_SKIP() << "Vulkan backend absent; disk session unqualified";
#endif
}

TEST(RenderSessionProbe, CpuKerrRenderProducesValidPngAndExr) {
    namespace fs = std::filesystem;
    const ScopedTemporaryDirectory temporary_directory("sirius-render-probe");
    const fs::path& dir = temporary_directory.path();

    const std::string png_path = (dir / "probe_kerr.png").string();
    const std::string exr_path = (dir / "probe_kerr.exr").string();

    // --- PNG render ---------------------------------------------------------
    ASSERT_EQ(RenderTo(png_path), SessionState::Complete) << "PNG render did not complete";
    ASSERT_TRUE(fs::exists(png_path)) << "PNG file was not written";

    int pw = 0, ph = 0, pc = 0;
    unsigned char* png = stbi_load(png_path.c_str(), &pw, &ph, &pc, 3);
    ASSERT_NE(png, nullptr) << "PNG failed to decode: " << stbi_failure_reason();
    EXPECT_EQ(pw, 64);
    EXPECT_EQ(ph, 64);

    // PNG (8-bit) is always finite; require it to be non-constant.
    bool png_varies = false;
    for (int i = 3; i < pw * ph * 3 && !png_varies; ++i) {
        if (png[i] != png[0]) png_varies = true;
    }
    EXPECT_TRUE(png_varies) << "PNG image is constant (nothing rendered)";
    stbi_image_free(png);

    // --- EXR render ---------------------------------------------------------
    ASSERT_EQ(RenderTo(exr_path), SessionState::Complete) << "EXR render did not complete";
    ASSERT_TRUE(fs::exists(exr_path)) << "EXR file was not written";

    float* exr = nullptr;
    int ew = 0, eh = 0;
    const char* err = nullptr;
    int ret = LoadEXR(&exr, &ew, &eh, exr_path.c_str(), &err);
    ASSERT_EQ(ret, TINYEXR_SUCCESS) << (err ? err : "unknown tinyexr error");
    EXPECT_EQ(ew, 64);
    EXPECT_EQ(eh, 64);

    // EXR carries linear HDR floats; require every sample finite and the frame
    // non-constant.
    bool exr_finite = true;
    bool exr_varies = false;
    const float first = exr[0];
    for (int i = 0; i < ew * eh * 4; ++i) {
        if (!std::isfinite(exr[i])) exr_finite = false;
        if ((i % 4) != 3 && exr[i] != first) exr_varies = true;  // Ignore the alpha channel.
    }
    EXPECT_TRUE(exr_finite) << "EXR contains non-finite samples";
    EXPECT_TRUE(exr_varies) << "EXR image is constant (nothing rendered)";
    std::free(exr);
}

TEST(RenderSessionProbe, CpuKerrRenderProducesValidPpmThroughTheOwnedWriter) {
    namespace fs = std::filesystem;
    const ScopedTemporaryDirectory temporary_directory("sirius-ppm-probe");
    const fs::path output = temporary_directory.path() / "probe_kerr.ppm";

    ASSERT_EQ(RenderTo(output.string()), SessionState::Complete) << "PPM render did not complete";
    std::ifstream file(output, std::ios::binary);
    ASSERT_TRUE(file) << "PPM file was not written";
    std::string magic;
    int width = 0;
    int height = 0;
    int maximum = 0;
    ASSERT_TRUE(file >> magic >> width >> height >> maximum);
    EXPECT_EQ(magic, "P6");
    EXPECT_EQ(width, 64);
    EXPECT_EQ(height, 64);
    EXPECT_EQ(maximum, 255);
    ASSERT_EQ(file.get(), '\n');

    std::vector<unsigned char> pixels(static_cast<std::size_t>(width) * height * 3);
    file.read(reinterpret_cast<char*>(pixels.data()), static_cast<std::streamsize>(pixels.size()));
    ASSERT_EQ(file.gcount(), static_cast<std::streamsize>(pixels.size()));
    EXPECT_EQ(file.peek(), std::char_traits<char>::eof())
        << "PPM has an unexpected trailing payload";
    const bool varies =
        std::any_of(pixels.begin() + 1, pixels.end(),
                    [first = pixels.front()](unsigned char value) { return value != first; });
    EXPECT_TRUE(varies) << "PPM image is constant (nothing rendered)";
}

TEST(RenderSessionProbe, PhysicalPointDetectorCompletesAMovingThinLensKerrFrame) {
    const ScopedTemporaryDirectory temporary_directory("sirius-point-detector-probe");
    SessionConfig config;
    // Keep the original aspect ratio and physical scene while exercising a
    // full 512-pixel shared region rather than the former eight-pixel witness.
    config.width = 32;
    config.height = 16;
    config.samples_per_pixel = 1;
    config.tile_size = 32;
    config.enable_parallel_rendering = false;
    // Inspect linear radiance: display grading can legitimately suppress this
    // deliberately faint catalogue. EXR also exercises the connected writer.
    config.output_path = (temporary_directory.path() / "detector.exr").string();
    sirius::test::ConfigureMovingKerrDetector(config);
    RenderSession session;
    const auto configured = session.Configure(config);
    ASSERT_TRUE(configured) << configured.error().Description();
    ASSERT_EQ(session.Execute(), SessionState::Complete);
    const auto pixels = session.GetDisplayBuffer().SnapshotFloatData();
    ASSERT_EQ(pixels.size(), 32u * 16u * 4u);
    double total = 0;
    bool varies = false;
    for (std::size_t i = 0; i < pixels.size(); ++i) {
        ASSERT_TRUE(std::isfinite(pixels[i]));
        if (i % 4 == 3) continue;
        total += pixels[i];
        varies = varies || pixels[i] != pixels[0];
    }
    EXPECT_GT(total, 0);
    EXPECT_TRUE(varies);

    // The sampler must use each worker's tracer and keep packet caches local.
    config.enable_parallel_rendering = true;
    config.thread_count = 2;
    config.tile_size = 16;
    config.output_path = (temporary_directory.path() / "detector-parallel.exr").string();
    RenderSession parallel;
    ASSERT_TRUE(parallel.Configure(config));
    ASSERT_EQ(parallel.Execute(), SessionState::Complete);
    EXPECT_EQ(parallel.GetDisplayBuffer().SnapshotFloatData(), pixels);
    EXPECT_EQ(parallel.GetTileScheduler().GetTileCount(), (config.width + 15) / 16);
}

TEST(RenderSessionProbe, PhysicalPointBlocksPreservePartialEdgesAndNonSquareSamples) {
    const ScopedTemporaryDirectory directory("sirius-point-block-edges");
    SessionConfig config;
    sirius::test::ConfigureMovingKerrDetector(config);
    // Flat transport keeps this scheduling regression bounded. A full block
    // and a one-pixel edge exercise cache replacement and scalar fallback;
    // three SPP retain distinct film/pupil samples and their accumulation order.
    config.metric_id = sirius::core::MetricId::Minkowski;
    config.black_hole_mass = 0;
    config.black_hole_spin = 0;
    config.width = 33;
    config.height = 1;
    config.camera_fov = 1;
    config.point_starfield_config.star_count = 10000;
    config.samples_per_pixel = 3;
    config.tile_size = 64;
    config.enable_parallel_rendering = false;
    config.output_path = (directory.path() / "serial.exr").string();
    RenderSession serial;
    const auto configured = serial.Configure(config);
    ASSERT_TRUE(configured) << configured.error().Description();
    ASSERT_EQ(serial.Execute(), SessionState::Complete) << serial.GetErrorMessage();
    const auto pixels = serial.GetDisplayBuffer().SnapshotFloatData();
    ASSERT_EQ(pixels.size(), 132u);
    double total = 0;
    for (std::size_t i = 0; i < pixels.size(); ++i) {
        ASSERT_TRUE(std::isfinite(pixels[i]));
        if (i % 4 == 3)
            EXPECT_EQ(pixels[i], 1);
        else
            total += pixels[i];
    }
    EXPECT_GT(total, 0);
    // A tile boundary through the shared block must not change discovery,
    // kernel ownership or the partial block at the image edge.
    config.enable_parallel_rendering = true;
    config.thread_count = 2;
    config.tile_size = 16;
    config.output_path = (directory.path() / "parallel.exr").string();
    RenderSession parallel;
    ASSERT_TRUE(parallel.Configure(config));
    ASSERT_EQ(parallel.Execute(), SessionState::Complete) << parallel.GetErrorMessage();
    EXPECT_EQ(parallel.GetDisplayBuffer().SnapshotFloatData(), pixels);
    EXPECT_EQ(parallel.GetTileScheduler().GetTileCount(), (config.width + 15) / 16);
}

TEST(RenderSessionProbe, FilmAffectsDisplayOutputButNeverLinearExr) {
    namespace fs = std::filesystem;
    const ScopedTemporaryDirectory temporary_directory("sirius-film-probe");
    const fs::path& dir = temporary_directory.path();

    auto render = [&](bool film, const std::string& extension) {
        SessionConfig cfg;
        cfg.width = 16;
        cfg.height = 12;
        cfg.samples_per_pixel = 1;
        cfg.tile_size = 16;
        cfg.enable_parallel_rendering = false;
        cfg.metric_id = sirius::core::MetricId::Minkowski;
        cfg.black_hole_mass = 0.0;
        cfg.enable_disk = false;
        cfg.enable_bloom = false;
        cfg.enable_film_finish = film;
        cfg.film_config = sirius::render::FilmConfig::Interstellar();
        cfg.output_path = (dir / ("film_probe_" + std::to_string(film) + extension)).string();
        RenderSession session;
        EXPECT_TRUE(session.Configure(cfg));
        EXPECT_EQ(session.Execute(), SessionState::Complete);
        std::vector<float> copy = session.GetDisplayBuffer().SnapshotFloatData();
        fs::remove(cfg.output_path);
        return copy;
    };

    const auto display_plain = render(false, ".png");
    const auto display_film = render(true, ".png");
    ASSERT_EQ(display_plain.size(), display_film.size());
    double display_difference = 0.0;
    for (std::size_t i = 0; i < display_plain.size(); ++i) {
        display_difference += std::abs(static_cast<double>(display_plain[i] - display_film[i]));
    }
    EXPECT_GT(display_difference, 1.0e-4) << "enabled film pipeline was inert";

    const auto exr_plain = render(false, ".exr");
    const auto exr_film = render(true, ".exr");
    ASSERT_EQ(exr_plain.size(), exr_film.size());
    EXPECT_EQ(exr_plain, exr_film)
        << "film finish contaminated the untouched linear-HDR EXR branch";
}

TEST(RenderSessionProbe, StartIsAsynchronousAndCancellationIsTerminalWithoutOutput) {
    namespace fs = std::filesystem;
    const ScopedTemporaryDirectory temporary_directory("sirius-cancellation-probe");
    const fs::path output = temporary_directory.path() / "cancelled-render-must-not-exist.ppm";

    SessionConfig config = ProbeConfig(output.string());
    config.width = 512;
    config.height = 512;
    config.samples_per_pixel = 4;
    config.tile_size = 64;
    config.enable_parallel_rendering = false;

    RenderSession session;
    ASSERT_TRUE(session.Configure(config));
    const auto start = std::chrono::steady_clock::now();
    ASSERT_TRUE(session.Start());
    const double launch_seconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    EXPECT_LT(launch_seconds, 0.25) << "Start performed render work synchronously";

    EXPECT_TRUE(session.Cancel());
    session.WaitForCompletion();
    EXPECT_EQ(session.GetState(), SessionState::Cancelled);
    EXPECT_FALSE(fs::exists(output));
    EXPECT_FALSE(session.Start()) << "a terminal session must not be silently restarted";
}

TEST(RenderSessionProbe, CancellationInterruptsAnActivePrivateRayBeforePublication) {
    struct BlockingExecutor final : sirius::backend::TraceStepExecutor {
        std::promise<void> entered, resume;
        std::shared_future<void> resumed = resume.get_future().share();
        int steps = 0, rejected = 0;
        std::optional<sirius::core::CameraLaunch> Launch(
            sirius::core::IMetric& metric, double spin,
            const sirius::core::CameraRay& camera) override {
            return sirius::core::LaunchCameraRay(metric, spin, camera);
        }
        bool Step(sirius::core::Lightray& ray, sirius::core::IMetric& metric,
                  const sirius::core::IntegratorConfig& config,
                  sirius::core::Rk45CoupledState& coupled,
                  sirius::core::Rk45CoupledComparison& comparison) override {
            if (++steps == 1) entered.set_value();
            resumed.wait();
            return sirius::core::Geodesic::IntegrateStepRk45(ray, &metric, config, &coupled,
                                                             &comparison);
        }
        void RejectLastInterval() override { ++rejected; }
    } executor;
    auto entered = executor.entered.get_future();
    SessionConfig config;
    config.width = config.height = 4;
    config.samples_per_pixel = 1;
    config.metric_id = sirius::core::MetricId::Minkowski;
    config.black_hole_mass = config.black_hole_spin = 0;
    config.enable_disk = false;
    config.write_output = false;
    RenderSession session(executor, 1, 1);
    config.backend = sirius::render::RenderBackend::Vulkan;
    const auto configured = session.Configure(config);
    ASSERT_TRUE(configured) << configured.error().Description();
    ASSERT_TRUE(session.Start());
    const bool active = entered.wait_for(std::chrono::seconds(5)) == std::future_status::ready;
    const bool cancelled = session.Cancel();
    executor.resume.set_value();
    session.WaitForCompletion();
    ASSERT_TRUE(active) << "the private ray did not start";
    EXPECT_TRUE(cancelled);
    EXPECT_EQ(session.GetState(), SessionState::Cancelled);
    EXPECT_EQ(executor.steps, 1);
    EXPECT_EQ(executor.rejected, 1);
    EXPECT_EQ(session.GetTileScheduler().GetCompletedCount(), 0);
}

TEST(RenderSessionProbe, OnePointRegionFeedsConcurrentDeviceProbesAndCancelsPrivately) {
    struct BlockingExecutor final : sirius::backend::TraceStepExecutor {
        std::promise<void> entered, resume;
        std::shared_future<void> resumed = resume.get_future().share();
        std::atomic<int> steps{0}, rejected{0};
        std::optional<sirius::core::CameraLaunch> Launch(
            sirius::core::IMetric& metric, double spin,
            const sirius::core::CameraRay& camera) override {
            return sirius::core::LaunchCameraRay(metric, spin, camera);
        }
        bool Step(sirius::core::Lightray& ray, sirius::core::IMetric& metric,
                  const sirius::core::IntegratorConfig& config,
                  sirius::core::Rk45CoupledState& coupled,
                  sirius::core::Rk45CoupledComparison& comparison) override {
            if (++steps == 2) entered.set_value();
            resumed.wait();
            return sirius::core::Geodesic::IntegrateStepRk45(ray, &metric, config, &coupled,
                                                             &comparison);
        }
        void RejectLastInterval() override { ++rejected; }
    } executor;
    struct Capture {
        std::ostringstream output;
        std::streambuf* previous = std::cout.rdbuf(output.rdbuf());
        ~Capture() { std::cout.rdbuf(previous); }
    } capture;
    auto entered = executor.entered.get_future();
    SessionConfig config;
    sirius::test::ConfigureMovingKerrDetector(config);
    config.width = config.height = 4;
    config.samples_per_pixel = 1;
    config.write_output = false;
    config.backend = sirius::render::RenderBackend::Vulkan;
    RenderSession session(executor, 2, sirius::render::kPointDetectorBlockEdge);
    const auto configured = session.Configure(config);
    ASSERT_TRUE(configured) << configured.error().Description();
    ASSERT_TRUE(session.Start());
    const bool active = entered.wait_for(std::chrono::seconds(5)) == std::future_status::ready;
    const bool cancelled = session.Cancel();
    executor.resume.set_value();
    session.WaitForCompletion();
    EXPECT_TRUE(active) << "one detector region did not submit concurrent independent rays";
    EXPECT_TRUE(cancelled);
    EXPECT_EQ(session.GetState(), SessionState::Cancelled);
    EXPECT_EQ(executor.steps.load(), 2);
    EXPECT_EQ(executor.rejected.load(), 2);
    EXPECT_EQ(session.GetTileScheduler().GetCompletedCount(), 0);
    const auto output = capture.output.str();
    EXPECT_EQ(output.find(sirius::render::kSceneEvidencePrefix), std::string::npos);
    EXPECT_EQ(output.find(sirius::render::kVulkanEvidencePrefix), std::string::npos);
    const auto prefix = sirius::render::kSourceSceneEvidencePrefix;
    const auto source = output.find(prefix);
    ASSERT_NE(source, std::string::npos);
    EXPECT_EQ(output.find(prefix, source + prefix.size()), std::string::npos);
    const auto record = nlohmann::json::parse(
        output.substr(source + prefix.size(), output.find('\n', source) - source - prefix.size()));
    EXPECT_EQ(record["point_star_count"], 100000);
    EXPECT_EQ(record["point_seed"], 42);
    EXPECT_EQ(record["backend"], "Vulkan");
    EXPECT_EQ(record["width"], 4);
}

TEST(RenderSessionProbe, CompletionCallbackCanReenterLifecycleWithoutDeadlock) {
    SessionConfig config = ProbeConfig("unused.ppm");
    config.width = 8;
    config.height = 8;
    config.tile_size = 8;
    config.samples_per_pixel = 1;
    config.metric_id = sirius::core::MetricId::Minkowski;
    config.black_hole_mass = 0.0;
    config.black_hole_spin = 0.0;
    config.enable_disk = false;
    config.write_output = false;
    config.output_path = SessionConfig{}.output_path;

    RenderSession session;
    const auto issue = sirius::render::SessionConfigIssue(config);
    ASSERT_FALSE(issue.has_value()) << *issue;
    ASSERT_TRUE(session.Configure(config));
    bool callback_completed = false;
    bool callback_configured = true;
    session.SetCompletionCallback([&](SessionState state, const std::string&) {
        EXPECT_EQ(state, SessionState::Complete);
        callback_configured = session.Configure(config).has_value();
        session.WaitForCompletion();
        callback_completed = true;
    });

    EXPECT_EQ(session.Execute(), SessionState::Complete);
    EXPECT_TRUE(callback_completed);
    EXPECT_FALSE(callback_configured);
}

TEST(RenderSessionProbe, PointStarfieldRejectsValuesItsGeneratorWouldClamp) {
    SessionConfig config;
    config.point_starfield = true;
    config.point_starfield_config.star_count = std::numeric_limits<std::uint32_t>::max();
    const auto issue = sirius::render::SessionConfigIssue(config);
    ASSERT_TRUE(issue.has_value());
    EXPECT_NE(issue->find("point-starfield"), std::string::npos);

    config.point_starfield_config = sirius::core::PointStarfieldConfig{};
    config.point_starfield_config.min_distance_pc = std::numeric_limits<float>::quiet_NaN();
    EXPECT_TRUE(sirius::render::SessionConfigIssue(config).has_value());
}

TEST(RenderSessionProbe, SceneEvidenceBindsCanonicalTypedConfiguration) {
    SessionConfig config;
    config.backend = sirius::render::RenderBackend::Vulkan;
    config.metric_id = sirius::core::MetricId::Kerr;
    config.black_hole_spin = 0.9;
    config.width = 5616;
    config.height = 4096;
    config.samples_per_pixel = 4;
    config.camera_fov = 60.0f;
    config.enable_disk = false;
    config.ray_bundles = true;
    config.point_starfield = true;
    config.point_starfield_config.seed = 4294967295u;
    config.point_starfield_config.min_distance_pc = 2.5f;
    config.point_starfield_config.max_distance_pc = 9000.0f;
    config.camera_beta_forward = 0.1;
    config.camera_beta_up = 0.02;
    config.camera_beta_right = -0.01;
    config.lens_type = sirius::core::LensType::ThinLens;
    config.camera_focal_length = 50.0f;
    config.camera_aperture = 2.8f;
    config.camera_focus_distance = 30.0f;

    const std::string evidence = sirius::render::SessionSceneEvidenceJson(config, 100000);
    EXPECT_TRUE(
        evidence.starts_with("{\"schema\":\"sirius-render-scene-v1\",\"backend\":\"Vulkan\","));
    for (const char* field : {
             "\"metric\":\"Kerr\"",
             "\"width\":5616",
             "\"height\":4096",
             "\"samples_per_pixel\":4",
             "\"field_of_view\":60",
             "\"disk_enabled\":false",
             "\"ray_bundles\":true",
             "\"point_starfield\":true",
             "\"point_star_count\":100000",
             "\"point_brightness_scale\":100",
             "\"point_seed\":4294967295",
             "\"point_min_distance_pc\":2.5",
             "\"point_max_distance_pc\":9000",
             "\"camera_beta\":[",
             "\"lens\":\"ThinLens\"",
             "\"focal_length\":50",
             "\"aperture\":",
             "\"focus_distance\":30",
         }) {
        EXPECT_NE(evidence.find(field), std::string::npos) << field;
    }
    EXPECT_EQ(evidence.back(), '}');
}

TEST(RenderSessionProbe, TypedNumericBoundariesMatchTheExternalConfigurationBoundary) {
    SessionConfig config;

    config.camera_focal_length = 10000.1f;
    EXPECT_TRUE(sirius::render::SessionConfigIssue(config).has_value());
    config.camera_focal_length = 50.0f;

    config.volumetric_tau_midplane = 1.0e6f + 1.0f;
    EXPECT_TRUE(sirius::render::SessionConfigIssue(config).has_value());
    config.volumetric_tau_midplane = 10.0f;

    config.enable_turbulence = true;
    EXPECT_EQ(sirius::render::SessionConfigIssue(config),
              "turbulence requires volumetric transfer");
    config.enable_turbulence = false;

    config.exposure = 100.1f;
    EXPECT_TRUE(sirius::render::SessionConfigIssue(config).has_value());
    config.exposure = 3.0f;

    config.enable_motion_blur = true;
    config.shutter_time = 1000.1f;
    EXPECT_TRUE(sirius::render::SessionConfigIssue(config).has_value());
    config.enable_motion_blur = false;
    config.shutter_time = sirius::core::kDefaultMotionBlurShutterTime;

    config.enable_film_finish = true;
    config.film_config.halation_radius = std::numeric_limits<float>::infinity();
    EXPECT_EQ(sirius::render::SessionConfigIssue(config),
              "film-finish parameters are outside the represented domain");
    config.film_config.halation_radius = 257.0f;
    EXPECT_EQ(sirius::render::SessionConfigIssue(config),
              "film-finish parameters are outside the represented domain");
    config.film_config.halation_radius = 8.0f;
    EXPECT_FALSE(sirius::render::SessionConfigIssue(config).has_value());

    config.enable_disk = false;
    config.metric_id = sirius::core::MetricId::DeSitter;
    config.black_hole_mass = 0.0;
    EXPECT_EQ(sirius::render::SessionConfigIssue(config), "de-Sitter requires positive lambda");
    config.cosmological_constant = 0.001;
    EXPECT_FALSE(sirius::render::SessionConfigIssue(config).has_value());
    config.observer_distance = 55.0;
    EXPECT_EQ(sirius::render::SessionConfigIssue(config),
              "positive-lambda observer must remain inside the cosmological trace boundary");
    config.observer_distance = 50.0;

    config.metric_id = sirius::core::MetricId::SchwarzschildDeSitter;
    config.black_hole_mass = 2.0;
    EXPECT_FALSE(sirius::render::SessionConfigIssue(config).has_value());
    config.cosmological_constant = 0.02;
    EXPECT_EQ(sirius::render::SessionConfigIssue(config),
              "positive-lambda observer must remain inside the cosmological trace boundary");
    config.cosmological_constant = 0.03;
    EXPECT_EQ(sirius::render::SessionConfigIssue(config),
              "Schwarzschild-de-Sitter requires 9*lambda*mass^2 < 1 (sub-Nariai sector)");
}

TEST(RenderSessionProbe, PolarisedRequestsDeclineAndTwoSheetIsRepresented) {
    SessionConfig config;
    config.metric_id = sirius::core::MetricId::Kerr;
    config.black_hole_spin = 0.7;
    config.color_mode = sirius::core::color_modes::Mode::Polarisation;
    config.enable_polarisation = true;

    config.enable_volumetric_disk = true;
    EXPECT_EQ(sirius::render::SessionConfigIssue(config),
              "polarisation is not represented for volumetric transfer");

    config.enable_volumetric_disk = false;
    config.enable_motion_blur = true;
    EXPECT_EQ(sirius::render::SessionConfigIssue(config),
              "polarisation is not represented with temporal disk motion blur");

    config.enable_motion_blur = false;
    config.color_mode = sirius::core::color_modes::Mode::TrueColor;
    config.enable_polarisation = false;
    config.enable_disk = false;
    config.metric_id = sirius::core::MetricId::MorrisThorne;
    config.black_hole_spin = 0.0;
    config.black_hole_mass = 1.0;
    config.wormhole_topology = sirius::render::WormholeTopology::OneSheetCapture;
    EXPECT_EQ(sirius::render::SessionConfigIssue(config),
              "metrics without a mass parameter require mass to be zero");

    config.black_hole_mass = 0.0;
    config.wormhole_topology = sirius::render::WormholeTopology::TwoSheet;
    EXPECT_FALSE(sirius::render::SessionConfigIssue(config).has_value());
}

TEST(DisplayBuffer, NonFiniteRadianceIsIdentifiedBeforeEncoding) {
    sirius::render::DisplayBuffer display;
    display.Initialise(1, 1);
    const float invalid[] = {0.0f, std::numeric_limits<float>::quiet_NaN(), 0.0f, 1.0f};
    display.UpdateTile(0, 0, 1, 1, invalid);
    const auto bad = display.FirstNonFiniteIndex();
    ASSERT_TRUE(bad.has_value());
    EXPECT_EQ(*bad, 1u);
}

TEST(DisplayBuffer, MalformedDimensionsAndTilesFailClosed) {
    sirius::render::DisplayBuffer display;
    EXPECT_DEATH(display.Initialise(-1, 1), "precondition.*enforced, terminating");

    display.Initialise(1, 1);
    EXPECT_DEATH(display.UpdateTile(0, 0, 1, 1, nullptr), "precondition.*enforced, terminating");
    const float pixel[] = {0.0f, 0.0f, 0.0f, 1.0f};
    EXPECT_DEATH(display.UpdateTile(1, 0, 1, 1, pixel), "precondition.*enforced, terminating");
    EXPECT_EQ(display.GetUpdateCounter(), 0u);
}

TEST(RenderSessionProbe, CpuPolarisationModeConsumesTransportedDiskStokes) {
    auto render = [](sirius::core::color_modes::Mode mode) {
        SessionConfig cfg;
        cfg.width = 24;
        cfg.height = 24;
        cfg.tile_size = 32;
        cfg.samples_per_pixel = 1;
        cfg.enable_parallel_rendering = false;
        cfg.write_output = false;
        cfg.metric_id = sirius::core::MetricId::Kerr;
        cfg.black_hole_mass = 1.0;
        cfg.black_hole_spin = 0.7;
        cfg.observer_inclination = 75.0 * std::numbers::pi / 180.0;
        cfg.color_mode = mode;
        cfg.enable_polarisation = mode == sirius::core::color_modes::Mode::Polarisation;
        cfg.enable_bloom = false;
        cfg.tonemapper = sirius::core::TonemapType::None;
        cfg.exposure = 1.0f;
        cfg.contrast = 1.0f;
        cfg.saturation = 1.0f;

        RenderSession session;
        EXPECT_TRUE(session.Configure(cfg));
        EXPECT_EQ(session.Execute(), SessionState::Complete);
        return session.GetDisplayBuffer().SnapshotFloatData();
    };

    const auto true_color = render(sirius::core::color_modes::Mode::TrueColor);
    const auto polarisation = render(sirius::core::color_modes::Mode::Polarisation);
    ASSERT_EQ(true_color.size(), polarisation.size());

    double difference = 0.0;
    bool finite = true;
    for (std::size_t i = 0; i < polarisation.size(); ++i) {
        finite = finite && std::isfinite(polarisation[i]);
        if ((i % 4) != 3) {
            difference += std::abs(static_cast<double>(polarisation[i] - true_color[i]));
        }
    }
    EXPECT_TRUE(finite);
    EXPECT_GT(difference, 1.0e-3)
        << "transported disk Stokes data did not reach the rendered colour branch";
}

// The Morris-Thorne wormhole renders on the CPU path through the exact isotropic
// Cartesian Ellis chart: the session must complete (not decline), and the frame must show
// the one-sheet wormhole structure - some rays captured at the throat (the
// dark centre) and some escaping past it (the lensed background), so the
// image is non-constant with genuinely dark pixels present.
TEST(RenderSessionProbe, CpuMorrisThorneRenderCompletes) {
    namespace fs = std::filesystem;
    const ScopedTemporaryDirectory temporary_directory("sirius-wormhole-probe");
    const fs::path& dir = temporary_directory.path();
    const std::string png_path = (dir / "probe_wormhole.png").string();

    SessionConfig cfg;
    cfg.width = 64;
    cfg.height = 64;
    cfg.samples_per_pixel = 4;
    cfg.tile_size = 64;
    cfg.enable_parallel_rendering = false;
    cfg.metric_id = sirius::core::MetricId::MorrisThorne;
    cfg.black_hole_mass = 0.0;
    cfg.enable_disk = false;
    // Throat large enough that its shadow spans several pixels at 64x64 with
    // the default observer distance; a b0 = 1 throat subtends ~1 pixel and
    // vanishes under sample jitter and tonemapping.
    cfg.throat_radius = 5.0;
    // The probe asserts trace physics (captured versus escaped rays), so the
    // film bloom stays off: at these scales it floods the throat shadow with
    // light from the surrounding Einstein ring.
    cfg.enable_bloom = false;
    cfg.output_path = png_path;

    RenderSession session;
    ASSERT_TRUE(session.Configure(cfg)) << "Session must accept the CPU wormhole config";
    ASSERT_EQ(session.Execute(), SessionState::Complete)
        << "Morris-Thorne must render on the CPU path, not decline";
    ASSERT_TRUE(fs::exists(png_path));

    int pw = 0, ph = 0, pc = 0;
    unsigned char* png = stbi_load(png_path.c_str(), &pw, &ph, &pc, 3);
    ASSERT_NE(png, nullptr) << stbi_failure_reason();
    EXPECT_EQ(pw, 64);
    EXPECT_EQ(ph, 64);

    bool varies = false;
    int dark_pixels = 0;
    for (int p = 0; p < pw * ph; ++p) {
        unsigned char r = png[3 * p], g = png[3 * p + 1], b = png[3 * p + 2];
        if (r != png[0] || g != png[1] || b != png[2]) varies = true;
        if (r < 8 && g < 8 && b < 8) ++dark_pixels;
    }
    EXPECT_TRUE(varies) << "Wormhole frame is constant (nothing rendered)";
    EXPECT_GT(dark_pixels, 0) << "No throat shadow: no rays were captured";
    EXPECT_LT(dark_pixels, pw * ph) << "Frame entirely dark: no rays escaped";
    stbi_image_free(png);
}

TEST(RenderSessionProbe, EveryRegisteredCpuMetricCompletesAFrame) {
    std::size_t advertised_cpu_metrics = 0;
    for (const auto& info : sirius::core::MetricRegistry()) {
        if (!info.cpu_supported) continue;
        ++advertised_cpu_metrics;
        SCOPED_TRACE(info.canonical_name);

        SessionConfig cfg;
        cfg.width = 4;
        cfg.height = 4;
        cfg.samples_per_pixel = 1;
        cfg.tile_size = 8;
        cfg.enable_parallel_rendering = false;
        cfg.write_output = false;
        cfg.enable_disk = false;
        cfg.enable_bloom = false;
        cfg.metric_id = info.id;

        switch (info.id) {
            case sirius::core::MetricId::Minkowski:
                cfg.black_hole_mass = 0.0;
                break;
            case sirius::core::MetricId::Kerr:
                cfg.black_hole_spin = 0.5;
                break;
            case sirius::core::MetricId::ReissnerNordstrom:
                cfg.black_hole_charge = 0.3;
                break;
            case sirius::core::MetricId::KerrNewman:
                cfg.black_hole_spin = 0.3;
                cfg.black_hole_charge = 0.3;
                break;
            case sirius::core::MetricId::DeSitter:
                cfg.black_hole_mass = 0.0;
                cfg.cosmological_constant = 0.001;
                break;
            case sirius::core::MetricId::SchwarzschildDeSitter:
                cfg.cosmological_constant = 0.001;
                break;
            case sirius::core::MetricId::Schwarzschild:
                break;
            case sirius::core::MetricId::MorrisThorne:
            case sirius::core::MetricId::Alcubierre:
                cfg.black_hole_mass = 0.0;
                break;
        }

        RenderSession session;
        ASSERT_TRUE(session.Configure(cfg));
        EXPECT_EQ(session.Execute(), SessionState::Complete);
    }
    EXPECT_EQ(advertised_cpu_metrics, 9u);
}

TEST(RenderSessionProbe, NumericalRayFailureKeepsCpuTilesPrivateAndPreventsOutput) {
    const ScopedTemporaryDirectory temporary_directory("sirius-cpu-ray-failure");
    for (bool parallel : {false, true}) {
        for (bool with_callback : {false, true}) {
            SCOPED_TRACE(with_callback);
            SCOPED_TRACE(parallel);
            auto config = ProbeConfig(
                (temporary_directory.path() / (parallel ? "parallel.ppm" : "sequential.ppm"))
                    .string());
            config.width = 1;
            config.height = 1;
            config.tile_size = 1;
            config.enable_parallel_rendering = parallel;
            config.thread_count = parallel ? 2 : 0;
            config.observer_distance = 30.0;
            config.observer_inclination = 80.0 * std::numbers::pi / 180.0;
            config.enable_disk = false;
            // A represented timelike observer whose highly boosted rays reach the
            // integrator's minimum step. This exercises the real tracer failure,
            // rather than substituting a failed shading callback.
            config.camera_beta_forward = -0.9999999999;
            // Obtain the diagnostic values from the real tracer independently of
            // the session's error formatting and callback routing.
            sirius::core::KerrSchildFamily metric(sirius::core::KerrSchildParams::Kerr(
                config.black_hole_mass, config.black_hole_spin));
            const auto domain = sirius::render::BuildTraceDomainParameters(
                {.metric_id = config.metric_id,
                 .metric_mass = config.black_hole_mass,
                 .cosmological_constant = 0.0,
                 .observer_radius = config.observer_distance,
                 .throat_radius = 1.0,
                 .bubble_radius = 1.0,
                 .bubble_sigma = 1.0});
            sirius::backend::TracerConfig trace;
            trace.enable_disk = false;
            trace.escape_radius = domain.escape_radius;
            trace.horizon_factor = 1.0f;
            trace.max_steps = sirius::render::kRenderTraceMaximumAttempts;
            trace.integrator.initial_step = domain.cpu_initial_step;
            trace.integrator.max_step = domain.max_step;
            trace.integrator.min_step = domain.cpu_min_step;
            trace.integrator.abs_tolerance = 5e-6f;
            trace.integrator.rel_tolerance = 5e-6f;
            sirius::backend::GeodesicTracer tracer(&metric, trace);
            sirius::core::CameraConfig camera_config;
            camera_config.r = config.observer_distance;
            camera_config.theta = config.observer_inclination;
            camera_config.fov = config.camera_fov;
            camera_config.width = config.width;
            camera_config.height = config.height;
            camera_config.beta_x = config.camera_beta_forward;
            sirius::core::PinholeCamera camera(camera_config);
            std::optional<sirius::backend::TraceResult> first;
            sirius::core::ForEachCameraSample(config.samples_per_pixel, [&](const auto& sample) {
                if (first) return;
                const auto projection = camera.ProjectFilmForObserver(
                    sample.image_u, sample.image_v, sample.pupil_u, sample.pupil_v);
                ASSERT_TRUE(projection);
                first = tracer.Trace(projection->ray);
            });
            ASSERT_TRUE(first.has_value());
            ASSERT_TRUE(first->numerical_failure);
            ASSERT_NE(first->coupled_failure, sirius::core::CoupledStepFailure::WorkLimit);
            const auto diagnostic = std::format(
                "numerical ray failure (coupled failure: {}, integrator termination: {}, attempts: "
                "{}, "
                "accepted affine distance: {})",
                sirius::core::CoupledStepFailureName(first->coupled_failure),
                first->integrator_termination, first->steps_taken, first->affine_length);
            RenderSession session;
            SessionState callback_state = SessionState::Idle;
            std::string message;
            if (with_callback) {
                session.SetCompletionCallback([&](SessionState state, const std::string& detail) {
                    callback_state = state;
                    message = detail;
                });
            }
            ASSERT_TRUE(session.Configure(config));
            testing::internal::CaptureStderr();
            const auto terminal = session.Execute();
            const auto stderr_output = testing::internal::GetCapturedStderr();
            EXPECT_EQ(terminal, SessionState::Failed);
            EXPECT_NE(stderr_output.find(diagnostic), std::string::npos) << stderr_output;
            EXPECT_NE(stderr_output.find("sample 0"), std::string::npos) << stderr_output;
            if (with_callback) {
                EXPECT_EQ(callback_state, SessionState::Failed);
                EXPECT_NE(message.find(diagnostic), std::string::npos) << message;
                EXPECT_NE(stderr_output.find(message), std::string::npos) << stderr_output;
            } else {
                EXPECT_EQ(callback_state, SessionState::Idle);
                EXPECT_TRUE(message.empty());
            }
            EXPECT_EQ(session.GetTileScheduler().GetCompletedCount(), 0);
            EXPECT_FALSE(std::filesystem::exists(config.output_path));
            const auto frame = session.GetDisplayBuffer().SnapshotFloatData();
            ASSERT_EQ(frame.size(), 4u);
            EXPECT_EQ(frame, std::vector<float>(4, 0.0f));
        }
    }
}

TEST(RenderSessionProbe, ExhaustedNonKerrRayCannotPublishACompletedBlackFrame) {
    const ScopedTemporaryDirectory temporary_directory("sirius-cpu-work-limit");
    for (const bool parallel : {false, true}) {
        SCOPED_TRACE(parallel);
        SessionConfig config;
        config.width = config.height = config.tile_size = 1;
        config.samples_per_pixel = 1;
        config.enable_parallel_rendering = parallel;
        config.thread_count = parallel ? 2 : 0;
        config.metric_id = sirius::core::MetricId::MorrisThorne;
        config.black_hole_mass = 0.0;
        config.enable_disk = false;
        config.enable_bloom = false;
        config.observer_distance = 1000.0;
        config.camera_beta_forward = 0.9999;
        config.output_path =
            (temporary_directory.path() / (parallel ? "parallel.ppm" : "sequential.ppm")).string();
        // This receding observer gives the central past ray coordinate speed
        // sqrt((1-beta)/(1+beta)), about 0.0071. It cannot reach either the
        // unit throat or enclosing escape sphere in 20000 two-unit attempts.
        // This exercises the production work envelope without reducing it or
        // replacing the tracer with a failed callback.
        RenderSession session;
        SessionState callback_state = SessionState::Idle;
        std::string message;
        session.SetCompletionCallback([&](SessionState state, const std::string& detail) {
            callback_state = state;
            message = detail;
        });
        ASSERT_TRUE(session.Configure(config));
        EXPECT_EQ(session.Execute(), SessionState::Failed);
        EXPECT_EQ(callback_state, SessionState::Failed);
        EXPECT_NE(message.find("ray work limit exhausted"), std::string::npos) << message;
        EXPECT_EQ(session.GetTileScheduler().GetCompletedCount(), 0);
        EXPECT_FALSE(std::filesystem::exists(config.output_path));
        EXPECT_EQ(session.GetDisplayBuffer().SnapshotFloatData(), std::vector<float>(4, 0.0f));

        // A completed capture remains a valid dark image on the same route.
        config.observer_distance = 30.0;
        config.camera_beta_forward = 0.0;
        config.write_output = false;
        config.output_path = SessionConfig{}.output_path;
        RenderSession captured;
        ASSERT_TRUE(captured.Configure(config));
        EXPECT_EQ(captured.Execute(), SessionState::Complete);
        EXPECT_EQ(captured.GetTileScheduler().GetCompletedCount(), 1);
    }
}

TEST(RenderSessionProbe, LaterGoodCameraSampleCannotEraseCpuNumericalFailure) {
    const ScopedTemporaryDirectory temporary_directory("sirius-cpu-sample-failure");
    auto config = ProbeConfig((temporary_directory.path() / "mixed.ppm").string());
    config.width = 1;
    config.height = 1;
    config.tile_size = 1;
    config.samples_per_pixel = 3;
    config.observer_distance = 30.0;
    config.observer_inclination = 80.0 * std::numbers::pi / 180.0;
    config.camera_fov = 170.0f;
    config.camera_beta_forward = 0.999999;
    config.enable_disk = false;

    // Independently trace the unchanged deterministic sample packets. The
    // receding observer makes the central past ray too slow to reach capture
    // within the production attempt envelope; both off-axis rays escape.
    // Thus a later success must not erase the middle ray's numerical failure.
    // The older extreme approaching boost now makes all three rays decline
    // under coupled accuracy checks and cannot witness this ordering contract.
    sirius::core::KerrSchildFamily metric(sirius::core::KerrSchildParams::Kerr(1.0, 0.9));
    const auto domain =
        sirius::render::BuildTraceDomainParameters({.metric_id = config.metric_id,
                                                    .metric_mass = config.black_hole_mass,
                                                    .cosmological_constant = 0.0,
                                                    .observer_radius = config.observer_distance,
                                                    .throat_radius = 1.0,
                                                    .bubble_radius = 1.0,
                                                    .bubble_sigma = 1.0});
    sirius::backend::TracerConfig trace;
    trace.enable_disk = false;
    trace.escape_radius = domain.escape_radius;
    trace.horizon_factor = 1.0f;
    trace.max_steps = sirius::render::kRenderTraceMaximumAttempts;
    trace.integrator.initial_step = domain.cpu_initial_step;
    trace.integrator.max_step = domain.max_step;
    trace.integrator.min_step = domain.cpu_min_step;
    trace.integrator.abs_tolerance = 5e-6f;
    trace.integrator.rel_tolerance = 5e-6f;
    sirius::backend::GeodesicTracer tracer(&metric, trace);
    sirius::core::CameraConfig camera_config;
    camera_config.r = config.observer_distance;
    camera_config.theta = config.observer_inclination;
    camera_config.fov = config.camera_fov;
    camera_config.width = config.width;
    camera_config.height = config.height;
    camera_config.beta_x = config.camera_beta_forward;
    sirius::core::PinholeCamera camera(camera_config);
    std::vector<bool> failures;
    std::vector<sirius::backend::TraceResult::Outcome> outcomes;
    std::vector<sirius::core::CoupledStepFailure> reasons;
    sirius::core::ForEachCameraSample(config.samples_per_pixel, [&](const auto& sample) {
        const auto ray = camera.GenerateRayForObserver(0, 0, sample.image_u, sample.image_v,
                                                       sample.pupil_u, sample.pupil_v);
        const auto result = tracer.Trace(ray);
        failures.push_back(result.numerical_failure);
        outcomes.push_back(result.outcome);
        reasons.push_back(result.coupled_failure);
    });
    ASSERT_EQ(failures, (std::vector<bool>{false, true, false}));
    using Outcome = sirius::backend::TraceResult::Outcome;
    ASSERT_EQ(outcomes,
              (std::vector<Outcome>{Outcome::Escaped, Outcome::MaxSteps, Outcome::Escaped}));
    ASSERT_EQ(reasons[1], sirius::core::CoupledStepFailure::WorkLimit);

    for (bool parallel : {false, true}) {
        SCOPED_TRACE(parallel);
        config.enable_parallel_rendering = parallel;
        config.thread_count = parallel ? 2 : 0;
        RenderSession session;
        std::string message;
        session.SetCompletionCallback(
            [&](SessionState, const std::string& detail) { message = detail; });
        ASSERT_TRUE(session.Configure(config));
        EXPECT_EQ(session.Execute(), SessionState::Failed);
        EXPECT_NE(message.find("sample 1"), std::string::npos) << message;
        EXPECT_NE(message.find("ray work limit exhausted"), std::string::npos) << message;
        EXPECT_EQ(session.GetTileScheduler().GetCompletedCount(), 0);
        EXPECT_FALSE(std::filesystem::exists(config.output_path));
        EXPECT_EQ(session.GetDisplayBuffer().SnapshotFloatData(), std::vector<float>(4, 0.0f));
    }
}
