// End-to-end CPU render-session probe (session-level gate).
//
// Drives RenderSession CPU-only at 64x64, 4 spp, Kerr a=0.9, writing a PNG and
// an EXR to a temp directory. Asserts: the session reaches Complete, both files
// exist, the PNG decodes (stb) and the EXR loads (tinyexr), and each decoded
// image is finite and non-constant.

#include "sirius/render/session/render_session.h"

#include <gtest/gtest.h>

#include "support/scoped_temporary_directory.h"
#include <stb_image.h>
#include <tinyexr.h>

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <limits>
#include <numbers>
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
    RenderSession session;
    if (!session.Configure(ProbeConfig(output_path))) {
        return SessionState::Failed;
    }
    return session.Execute();
}

}  // namespace

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

    RenderSession session;
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
              "turbulence and corona require volumetric transfer");
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

TEST(TileScheduler, ReinitialiseResetsCompletionLedger) {
    sirius::render::TileScheduler scheduler;
    scheduler.Initialise(128, 64, 64);
    ASSERT_EQ(scheduler.GetTileCount(), 2);
    const sirius::render::Tile* tile = scheduler.GetNextTile();
    ASSERT_NE(tile, nullptr);
    scheduler.CompleteTile(tile->id);
    ASSERT_EQ(scheduler.GetCompletedCount(), 1);

    scheduler.Initialise(64, 64, 64);
    EXPECT_EQ(scheduler.GetTileCount(), 1);
    EXPECT_EQ(scheduler.GetCompletedCount(), 0);
    EXPECT_FALSE(scheduler.AllComplete());
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
