// End-to-end CPU render-session probe (session-level gate).
//
// Drives RenderSession CPU-only at 64x64, 4 spp, Kerr a=0.9, writing a PNG and
// an EXR to a temp directory. Asserts: the session reaches Complete, both files
// exist, the PNG decodes (stb) and the EXR loads (tinyexr), and each decoded
// image is finite and non-constant.

#include "sirius/render/render_config.h"
#include "sirius/render/session/render_session.h"

#ifdef SIRIUS_HAS_VULKAN_BACKEND
#include "sirius/backend/device.h"
#endif

#include <gtest/gtest.h>

#include <stb_image.h>
#include <tinyexr.h>

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <limits>
#include <numbers>
#include <string>
#include <vector>

namespace {

using sirius::render::RenderSession;
using sirius::render::SessionConfig;
using sirius::render::SessionState;

SessionConfig ProbeConfig(const std::string& output_path) {
    SessionConfig cfg;
    cfg.width = 64;
    cfg.height = 64;
    cfg.samplesPerPixel = 4;
    cfg.tileSize = 64;
    cfg.enableParallelRendering = false;  // Deterministic single-thread for the probe.
    cfg.metricId = sirius::core::MetricId::Kerr;
    cfg.blackHoleMass = 1.0;
    cfg.blackHoleSpin = 0.9;
    cfg.outputPath = output_path;
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
    const fs::path dir = fs::temp_directory_path() / "sirius_render_probe";
    fs::create_directories(dir);

    const std::string pngPath = (dir / "probe_kerr.png").string();
    const std::string exrPath = (dir / "probe_kerr.exr").string();
    fs::remove(pngPath);
    fs::remove(exrPath);

    // --- PNG render ---------------------------------------------------------
    ASSERT_EQ(RenderTo(pngPath), SessionState::Complete) << "PNG render did not complete";
    ASSERT_TRUE(fs::exists(pngPath)) << "PNG file was not written";

    int pw = 0, ph = 0, pc = 0;
    unsigned char* png = stbi_load(pngPath.c_str(), &pw, &ph, &pc, 3);
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
    ASSERT_EQ(RenderTo(exrPath), SessionState::Complete) << "EXR render did not complete";
    ASSERT_TRUE(fs::exists(exrPath)) << "EXR file was not written";

    float* exr = nullptr;
    int ew = 0, eh = 0;
    const char* err = nullptr;
    int ret = LoadEXR(&exr, &ew, &eh, exrPath.c_str(), &err);
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

    fs::remove(pngPath);
    fs::remove(exrPath);
}

TEST(RenderSessionProbe, FilmAffectsDisplayOutputButNeverLinearExr) {
    namespace fs = std::filesystem;
    const fs::path dir = fs::temp_directory_path() / "sirius_render_probe";
    fs::create_directories(dir);

    auto render = [&](bool film, const std::string& extension) {
        SessionConfig cfg;
        cfg.width = 16;
        cfg.height = 12;
        cfg.samplesPerPixel = 1;
        cfg.tileSize = 16;
        cfg.enableParallelRendering = false;
        cfg.metricId = sirius::core::MetricId::Minkowski;
        cfg.blackHoleMass = 0.0;
        cfg.enableDisk = false;
        cfg.enableBloom = false;
        cfg.enableFilmSimulation = film;
        cfg.filmConfig = sirius::render::FilmConfig::Interstellar();
        cfg.outputPath = (dir / ("film_probe_" + std::to_string(film) + extension)).string();
        RenderSession session;
        EXPECT_TRUE(session.Configure(cfg));
        EXPECT_EQ(session.Execute(), SessionState::Complete);
        const float* pixels = session.GetDisplayBuffer().GetFloatData();
        std::vector<float> copy(pixels, pixels + cfg.width * cfg.height * 4);
        fs::remove(cfg.outputPath);
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
        << "film simulation contaminated the untouched linear-HDR EXR branch";
}

// Backend auto-resolution after the go-live flip (owner decision, 2026-07-18):
// 'auto' selects Vulkan only when the device, metric, and requested scene
// semantics are all represented, including exact non-square multisampling.
TEST(RenderSessionProbe, BackendAutoResolvesByDeviceRegistryAndCapabilities) {
    sirius::render::SiriusConfig cfg = sirius::render::SiriusConfig::defaults();
    cfg.metric.name = "Kerr";
    cfg.backend.preferred = "cpu";
    EXPECT_EQ(SessionConfig::FromSiriusConfig(cfg).backend, sirius::render::RenderBackend::Cpu)
        << "'cpu' must pin the CPU path";

    cfg.backend.preferred = "auto";
#ifdef SIRIUS_HAS_VULKAN_BACKEND
    const auto devices = sirius::backend::EnumerateVulkanDevices();
    const bool device_present = devices.has_value() && !devices->empty();
    EXPECT_EQ(
        SessionConfig::FromSiriusConfig(cfg).backend,
        device_present ? sirius::render::RenderBackend::Vulkan : sirius::render::RenderBackend::Cpu)
        << "auto must follow device presence when the full sampled scene is represented";

    cfg.volumetric.enabled = true;
    cfg.volumetric.samples = 129;
    EXPECT_EQ(SessionConfig::FromSiriusConfig(cfg).backend, sirius::render::RenderBackend::Cpu)
        << "auto must select CPU above the Vulkan volume-sample boundary";
    cfg.volumetric.enabled = false;
    cfg.volumetric.samples = 64;

    cfg.metric.name = "Reissner-Nordstrom";  // registry gpu_supported = false
    EXPECT_EQ(SessionConfig::FromSiriusConfig(cfg).backend, sirius::render::RenderBackend::Cpu)
        << "auto must resolve CPU for a metric the registry marks CPU-only";
#else
    EXPECT_EQ(SessionConfig::FromSiriusConfig(cfg).backend, sirius::render::RenderBackend::Cpu)
        << "auto must resolve CPU when the Vulkan backend is not compiled in";
#endif
}

TEST(RenderSessionProbe, ConfigurationConversionPreservesObserverAndDiskControls) {
    sirius::render::SiriusConfig cfg = sirius::render::SiriusConfig::defaults();
    cfg.backend.preferred = "cpu";
    cfg.observer.azimuth = 37.5;
    cfg.metric.temperatureModel = "ShakuraSunyaev";
    cfg.metric.diskTemperature = 42000.0f;
    cfg.render.threadCount = 7;
    cfg.observer.cameraBetaForward = 0.1;
    cfg.observer.lensModel = "ThinLens";
    cfg.observer.focalLength = 85.0f;
    cfg.observer.aperture = 1.4f;
    cfg.observer.focusDistance = 40.0f;
    cfg.colorMode = "Polarisation";
    cfg.motionBlur.enabled = false;
    cfg.motionBlur.shutterTime = 0.25f;
    cfg.motionBlur.samples = 7;

    const SessionConfig session = SessionConfig::FromSiriusConfig(cfg);
    EXPECT_NEAR(session.observerAzimuth, 37.5 * std::numbers::pi / 180.0, 1e-12);
    EXPECT_EQ(session.temperatureModel, sirius::render::DiskTemperatureModel::ShakuraSunyaev);
    EXPECT_FLOAT_EQ(session.diskTemperatureScale, 42000.0f);
    EXPECT_EQ(session.threadCount, 7);
    EXPECT_DOUBLE_EQ(session.cameraBetaForward, 0.1);
    EXPECT_EQ(session.lensType, sirius::core::LensType::ThinLens);
    EXPECT_FLOAT_EQ(session.cameraFocalLength, 85.0f);
    EXPECT_FLOAT_EQ(session.cameraAperture, 1.4f);
    EXPECT_FLOAT_EQ(session.shutterTime, 0.25f);
    EXPECT_EQ(session.motionBlurSamples, 7);
    EXPECT_FLOAT_EQ(session.cameraFocusDistance, 40.0f);
    EXPECT_EQ(session.colorMode, sirius::core::color_modes::Mode::Polarisation);
    EXPECT_TRUE(session.enablePolarisation);

    cfg.metric.temperatureModel = "UnknownTemperature";
    EXPECT_THROW(SessionConfig::FromSiriusConfig(cfg), std::invalid_argument);
    cfg.metric.temperatureModel = "NovikovThorne";
    cfg.observer.lensModel = "UnknownLens";
    EXPECT_THROW(SessionConfig::FromSiriusConfig(cfg), std::invalid_argument);
    cfg.observer.lensModel = "Pinhole";
    cfg.film.preset = "UnknownFilm";
    EXPECT_THROW(SessionConfig::FromSiriusConfig(cfg), std::invalid_argument);
    cfg.film.preset = "Interstellar";
    cfg.observer.lensModel = "Fisheye";
    EXPECT_EQ(SessionConfig::FromSiriusConfig(cfg).lensType, sirius::core::LensType::Fisheye);

    SessionConfig malformed = session;
    malformed.width = 1;
    malformed.height = 1;
    malformed.tileSize = 1;
    malformed.samplesPerPixel = 1;
    malformed.enableParallelRendering = false;
    malformed.temperatureModel = static_cast<sirius::render::DiskTemperatureModel>(255);
    RenderSession cpu_session;
    ASSERT_TRUE(cpu_session.Configure(malformed));
    EXPECT_EQ(cpu_session.Execute(), SessionState::Failed);

    malformed.temperatureModel = sirius::render::DiskTemperatureModel::NovikovThorne;
    malformed.colorMode = sirius::core::color_modes::Mode::TrueColor;
    malformed.enablePolarisation = true;
    RenderSession polarisation_session;
    ASSERT_TRUE(polarisation_session.Configure(malformed));
    EXPECT_EQ(polarisation_session.Execute(), SessionState::Failed);

    malformed.colorMode = sirius::core::color_modes::Mode::Polarisation;
    malformed.enablePolarisation = false;
    RenderSession polarisation_colour_session;
    ASSERT_TRUE(polarisation_colour_session.Configure(malformed));
    EXPECT_EQ(polarisation_colour_session.Execute(), SessionState::Failed);

    malformed.colorMode = static_cast<sirius::core::color_modes::Mode>(255);
    RenderSession invalid_colour_session;
    ASSERT_TRUE(invalid_colour_session.Configure(malformed));
    EXPECT_EQ(invalid_colour_session.Execute(), SessionState::Failed);

    malformed.colorMode = sirius::core::color_modes::Mode::TrueColor;
    malformed.tonemapper = static_cast<sirius::core::TonemapType>(255);
    RenderSession invalid_tonemapper_session;
    ASSERT_TRUE(invalid_tonemapper_session.Configure(malformed));
    EXPECT_EQ(invalid_tonemapper_session.Execute(), SessionState::Failed);

    malformed.tonemapper = sirius::core::TonemapType::Aces;
    malformed.metricId = sirius::core::MetricId::Schwarzschild;
    malformed.blackHoleSpin = 0.4;
    RenderSession mismatched_metric_session;
    ASSERT_TRUE(mismatched_metric_session.Configure(malformed));
    EXPECT_EQ(mismatched_metric_session.Execute(), SessionState::Failed);

    malformed.metricId = sirius::core::MetricId::KerrNewman;
    malformed.blackHoleSpin = 0.4;
    malformed.blackHoleCharge = 0.2;
    malformed.enableDisk = true;
    RenderSession charged_disk_session;
    ASSERT_TRUE(charged_disk_session.Configure(malformed));
    EXPECT_EQ(charged_disk_session.Execute(), SessionState::Failed);

    malformed.metricId = sirius::core::MetricId::Minkowski;
    malformed.blackHoleMass = 0.0;
    malformed.blackHoleSpin = 0.0;
    malformed.blackHoleCharge = 0.0;
    malformed.enableDisk = true;
    RenderSession inapplicable_disk_session;
    ASSERT_TRUE(inapplicable_disk_session.Configure(malformed));
    EXPECT_EQ(inapplicable_disk_session.Execute(), SessionState::Failed);

    malformed.metricId = sirius::core::MetricId::Schwarzschild;
    malformed.blackHoleMass = 1.0;
    malformed.blackHoleSpin = 0.0;
    malformed.blackHoleCharge = 0.0;
    malformed.enableDisk = true;
    malformed.width = 0;
    RenderSession invalid_dimensions_session;
    ASSERT_TRUE(invalid_dimensions_session.Configure(malformed));
    EXPECT_EQ(invalid_dimensions_session.Execute(), SessionState::Failed);

    const auto preview_path =
        std::filesystem::temp_directory_path() / "sirius_in_memory_preview_must_not_exist.ppm";
    std::filesystem::remove(preview_path);
    SessionConfig preview;
    preview.width = 8;
    preview.height = 8;
    preview.tileSize = 8;
    preview.samplesPerPixel = 1;
    preview.enableParallelRendering = false;
    preview.metricId = sirius::core::MetricId::Minkowski;
    preview.blackHoleMass = 0.0;
    preview.enableDisk = false;
    preview.writeOutput = false;
    preview.outputPath = preview_path.string();
    RenderSession preview_session;
    ASSERT_TRUE(preview_session.Configure(preview));
    EXPECT_EQ(preview_session.Execute(), SessionState::Complete);
    EXPECT_FALSE(std::filesystem::exists(preview_path));
}

TEST(RenderSessionProbe, StartIsAsynchronousAndCancellationIsTerminalWithoutOutput) {
    namespace fs = std::filesystem;
    const fs::path output =
        fs::temp_directory_path() / "sirius_cancelled_render_must_not_exist.ppm";
    fs::remove(output);

    SessionConfig config = ProbeConfig(output.string());
    config.width = 512;
    config.height = 512;
    config.samplesPerPixel = 4;
    config.tileSize = 64;
    config.enableParallelRendering = false;

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
    config.tileSize = 8;
    config.samplesPerPixel = 1;
    config.metricId = sirius::core::MetricId::Minkowski;
    config.blackHoleMass = 0.0;
    config.blackHoleSpin = 0.0;
    config.enableDisk = false;
    config.writeOutput = false;

    RenderSession session;
    ASSERT_TRUE(session.Configure(config));
    bool callback_completed = false;
    bool callback_configured = true;
    session.SetCompletionCallback([&](SessionState state, const std::string&) {
        EXPECT_EQ(state, SessionState::Complete);
        callback_configured = session.Configure(config);
        session.WaitForCompletion();
        callback_completed = true;
    });

    EXPECT_EQ(session.Execute(), SessionState::Complete);
    EXPECT_TRUE(callback_completed);
    EXPECT_FALSE(callback_configured);
}

TEST(RenderSessionProbe, PointStarfieldRejectsValuesItsGeneratorWouldClamp) {
    SessionConfig config;
    config.pointStarfield = true;
    config.starfieldConfig.star_count = std::numeric_limits<std::uint32_t>::max();
    const auto issue = sirius::render::SessionConfigIssue(config);
    ASSERT_TRUE(issue.has_value());
    EXPECT_NE(issue->find("point-starfield"), std::string::npos);

    config.starfieldConfig = sirius::core::StarfieldConfig{};
    config.starfieldConfig.min_distance_pc = std::numeric_limits<float>::quiet_NaN();
    EXPECT_TRUE(sirius::render::SessionConfigIssue(config).has_value());
}

TEST(RenderSessionProbe, TypedNumericBoundariesMatchTheExternalConfigurationBoundary) {
    SessionConfig config;

    config.cameraFocalLength = 10000.1f;
    EXPECT_TRUE(sirius::render::SessionConfigIssue(config).has_value());
    config.cameraFocalLength = 50.0f;

    config.volumetricTauMidplane = 1.0e6f + 1.0f;
    EXPECT_TRUE(sirius::render::SessionConfigIssue(config).has_value());
    config.volumetricTauMidplane = 10.0f;

    config.enableTurbulence = true;
    EXPECT_EQ(sirius::render::SessionConfigIssue(config),
              "turbulence and corona require volumetric transfer");
    config.enableTurbulence = false;

    config.exposure = 100.1f;
    EXPECT_TRUE(sirius::render::SessionConfigIssue(config).has_value());
    config.exposure = 3.0f;

    config.enableMotionBlur = true;
    config.shutterTime = 1000.1f;
    EXPECT_TRUE(sirius::render::SessionConfigIssue(config).has_value());
    config.enableMotionBlur = false;

    config.enableFilmSimulation = true;
    config.filmConfig.halation_radius = std::numeric_limits<float>::infinity();
    EXPECT_EQ(sirius::render::SessionConfigIssue(config),
              "film-simulation parameters are outside the represented domain");
    config.filmConfig.halation_radius = 257.0f;
    EXPECT_EQ(sirius::render::SessionConfigIssue(config),
              "film-simulation parameters are outside the represented domain");
    config.filmConfig.halation_radius = 8.0f;
    EXPECT_FALSE(sirius::render::SessionConfigIssue(config).has_value());
}

TEST(RenderSessionProbe, PolarisedAndTwoSheetRequestsDeclineAtTheTypedBoundary) {
    SessionConfig config;
    config.metricId = sirius::core::MetricId::Kerr;
    config.blackHoleSpin = 0.7;
    config.colorMode = sirius::core::color_modes::Mode::Polarisation;
    config.enablePolarisation = true;

    config.enableVolumetricDisk = true;
    EXPECT_EQ(sirius::render::SessionConfigIssue(config),
              "polarisation is not represented for volumetric transfer");

    config.enableVolumetricDisk = false;
    config.enableMotionBlur = true;
    EXPECT_EQ(sirius::render::SessionConfigIssue(config),
              "polarisation is not represented with temporal disk motion blur");

    config.enableMotionBlur = false;
    config.colorMode = sirius::core::color_modes::Mode::TrueColor;
    config.enablePolarisation = false;
    config.enableDisk = false;
    config.metricId = sirius::core::MetricId::MorrisThorne;
    config.blackHoleSpin = 0.0;
    config.wormholeTopology = sirius::render::WormholeTopology::TwoSheet;
    EXPECT_EQ(sirius::render::SessionConfigIssue(config),
              "two-sheet wormhole continuation and a second environment are not represented");
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

TEST(DisplayBuffer, MalformedDimensionsTilesAndGammaFailClosed) {
    sirius::render::DisplayBuffer display;
    display.Initialise(-1, 1);
    EXPECT_TRUE(display.SnapshotFloatData().empty());
    EXPECT_EQ(display.GetWidth(), 0);
    EXPECT_EQ(display.GetHeight(), 0);

    display.Initialise(1, 1);
    display.UpdateTile(0, 0, 1, 1, nullptr);
    EXPECT_EQ(display.GetUpdateCounter(), 0u);
    EXPECT_EQ(display.GetByteData(0.0f), nullptr);
    EXPECT_EQ(display.GetByteData(std::numeric_limits<float>::quiet_NaN()), nullptr);
}

TEST(RenderSessionProbe, CpuPolarisationModeConsumesTransportedDiskStokes) {
    auto render = [](sirius::core::color_modes::Mode mode) {
        SessionConfig cfg;
        cfg.width = 24;
        cfg.height = 24;
        cfg.tileSize = 32;
        cfg.samplesPerPixel = 1;
        cfg.enableParallelRendering = false;
        cfg.writeOutput = false;
        cfg.metricId = sirius::core::MetricId::Kerr;
        cfg.blackHoleMass = 1.0;
        cfg.blackHoleSpin = 0.7;
        cfg.observerInclination = 75.0 * std::numbers::pi / 180.0;
        cfg.colorMode = mode;
        cfg.enablePolarisation = mode == sirius::core::color_modes::Mode::Polarisation;
        cfg.enableBloom = false;
        cfg.tonemapper = sirius::core::TonemapType::None;
        cfg.exposure = 1.0f;
        cfg.contrast = 1.0f;
        cfg.saturation = 1.0f;

        RenderSession session;
        EXPECT_TRUE(session.Configure(cfg));
        EXPECT_EQ(session.Execute(), SessionState::Complete);
        const float* pixels = session.GetDisplayBuffer().GetFloatData();
        return std::vector<float>(pixels, pixels + cfg.width * cfg.height * 4);
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

// The Morris-Thorne wormhole renders on the CPU path through the Cartesian
// embedding: the session must complete (not decline), and the frame must show
// the one-sheet wormhole structure - some rays captured at the throat (the
// dark centre) and some escaping past it (the lensed background), so the
// image is non-constant with genuinely dark pixels present.
TEST(RenderSessionProbe, CpuMorrisThorneRenderCompletes) {
    namespace fs = std::filesystem;
    const fs::path dir = fs::temp_directory_path() / "sirius_render_probe";
    fs::create_directories(dir);
    const std::string pngPath = (dir / "probe_wormhole.png").string();
    fs::remove(pngPath);

    SessionConfig cfg;
    cfg.width = 64;
    cfg.height = 64;
    cfg.samplesPerPixel = 4;
    cfg.tileSize = 64;
    cfg.enableParallelRendering = false;
    cfg.metricId = sirius::core::MetricId::MorrisThorne;
    cfg.enableDisk = false;
    // Throat large enough that its shadow spans several pixels at 64x64 with
    // the default observer distance; a b0 = 1 throat subtends ~1 pixel and
    // vanishes under sample jitter and tonemapping.
    cfg.throatRadius = 5.0;
    // The probe asserts trace physics (captured versus escaped rays), so the
    // film bloom stays off: at these scales it floods the throat shadow with
    // light from the surrounding Einstein ring.
    cfg.enableBloom = false;
    cfg.outputPath = pngPath;

    RenderSession session;
    ASSERT_TRUE(session.Configure(cfg)) << "Session must accept the CPU wormhole config";
    ASSERT_EQ(session.Execute(), SessionState::Complete)
        << "Morris-Thorne must render on the CPU path, not decline";
    ASSERT_TRUE(fs::exists(pngPath));

    int pw = 0, ph = 0, pc = 0;
    unsigned char* png = stbi_load(pngPath.c_str(), &pw, &ph, &pc, 3);
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

    fs::remove(pngPath);
}
