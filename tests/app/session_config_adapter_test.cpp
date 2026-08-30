#include "sirius/app/config/session_config_adapter.h"

#include "sirius/render/exr_writer.h"
#include "sirius/render/session/render_session.h"
#include "sirius/render/trace_domain.h"

#ifdef SIRIUS_HAS_VULKAN_BACKEND
#include "sirius/backend/device.h"
#endif

#include <gtest/gtest.h>

#include <filesystem>
#include <numbers>

namespace sirius::app::test {

using render::RenderSession;
using render::SessionConfig;
using render::SessionState;

TEST(RenderSessionProbe, TraceDomainScalesWithMassAndEnclosesTheObserver) {
    const render::TraceDomainParameters baseline =
        render::BuildTraceDomainParameters(core::MetricId::Kerr, 1.0, 50.0);
    EXPECT_FLOAT_EQ(baseline.escape_radius, 200.0f);
    EXPECT_FLOAT_EQ(baseline.cpu_initial_step, 0.1f);
    EXPECT_FLOAT_EQ(baseline.cpu_min_step, 1.0e-5f);
    EXPECT_FLOAT_EQ(baseline.vulkan_min_step, 0.02f);
    EXPECT_FLOAT_EQ(baseline.max_step, 2.0f);

    const render::TraceDomainParameters small =
        render::BuildTraceDomainParameters(core::MetricId::Kerr, 0.1, 100.0);
    EXPECT_GT(small.escape_radius, 100.0f);
    EXPECT_FLOAT_EQ(small.escape_radius, 125.0f);
    EXPECT_FLOAT_EQ(small.max_step, 0.2f);

    const render::TraceDomainParameters large =
        render::BuildTraceDomainParameters(core::MetricId::Kerr, 100.0, 100000.0);
    EXPECT_GT(large.escape_radius, 100000.0f);
    EXPECT_FLOAT_EQ(large.escape_radius, 125000.0f);
    EXPECT_FLOAT_EQ(large.cpu_initial_step, 10.0f);
    EXPECT_FLOAT_EQ(large.cpu_min_step, 1.0e-3f);
    EXPECT_FLOAT_EQ(large.vulkan_min_step, 2.0f);
    EXPECT_FLOAT_EQ(large.max_step, 200.0f);

    const render::TraceDomainParameters massless =
        render::BuildTraceDomainParameters(core::MetricId::DeSitter, 0.0, 1000.0);
    EXPECT_GT(massless.escape_radius, 1000.0f);
    EXPECT_FLOAT_EQ(massless.max_step, 2.0f);

    const render::TraceDomainParameters non_mass_geometry =
        render::BuildTraceDomainParameters(core::MetricId::MorrisThorne, 0.0, 50.0);
    EXPECT_FLOAT_EQ(non_mass_geometry.escape_radius, 200.0f);
    EXPECT_FLOAT_EQ(non_mass_geometry.max_step, 2.0f);

    SessionConfig non_mass_session;
    non_mass_session.metric_id = core::MetricId::MorrisThorne;
    non_mass_session.enable_disk = false;
    EXPECT_EQ(render::SessionConfigIssue(non_mass_session),
              "metrics without a mass parameter require mass to be zero");
    non_mass_session.black_hole_mass = 0.0;
    EXPECT_FALSE(render::SessionConfigIssue(non_mass_session).has_value());
    non_mass_session.throat_radius = 2.0;
    EXPECT_FALSE(render::SessionConfigIssue(non_mass_session).has_value());
    non_mass_session.warp_velocity = 1.0;
    EXPECT_EQ(render::SessionConfigIssue(non_mass_session),
              "warp velocity, bubble radius, and bubble sigma apply only to Alcubierre");

    SessionConfig warp_session;
    warp_session.metric_id = core::MetricId::Alcubierre;
    warp_session.black_hole_mass = 0.0;
    warp_session.enable_disk = false;
    warp_session.warp_velocity = 1.0;
    warp_session.bubble_radius = 2.0;
    EXPECT_FALSE(render::SessionConfigIssue(warp_session).has_value());
    warp_session.throat_radius = 2.0;
    EXPECT_EQ(render::SessionConfigIssue(warp_session),
              "throat radius and wormhole topology apply only to Morris-Thorne");
}

TEST(RenderSessionProbe, GeometricMetadataNeverInventsPhysicalLengthUnits) {
    render::EXRMetadata metadata;
    metadata.black_hole_mass = 2.0;
    metadata.black_hole_spin = 0.5;
    metadata.observer_distance = 100.0;
    const std::string header = render::EXRWriter::GenerateMetadataHeader(metadata);

    EXPECT_NE(header.find("Metric Mass Parameter: 2 coordinate units"), std::string::npos);
    EXPECT_NE(header.find("Dimensionless Spin a/M: 0.5"), std::string::npos);
    EXPECT_NE(header.find("Observer Coordinate Radius: 100 coordinate units"), std::string::npos);
    EXPECT_EQ(header.find("M_sun"), std::string::npos);
}

TEST(RenderSessionProbe, FeatureSpecificControlsRequireOwnersAtTypedBoundary) {
    SiriusConfig app_defaults = SiriusConfig::Defaults();
    app_defaults.backend.preferred = "cpu";
    const auto projected_defaults = MakeSessionConfig(app_defaults);
    ASSERT_TRUE(projected_defaults.has_value()) << projected_defaults.error().Description();
    const SessionConfig typed_defaults;
    EXPECT_FLOAT_EQ(projected_defaults->bloom_intensity, typed_defaults.bloom_intensity);
    EXPECT_FLOAT_EQ(projected_defaults->bloom_threshold, typed_defaults.bloom_threshold);
    EXPECT_FLOAT_EQ(projected_defaults->volumetric_h_over_r, typed_defaults.volumetric_h_over_r);
    EXPECT_FLOAT_EQ(projected_defaults->volumetric_h_power, typed_defaults.volumetric_h_power);
    EXPECT_FLOAT_EQ(projected_defaults->volumetric_tau_midplane,
                    typed_defaults.volumetric_tau_midplane);
    EXPECT_EQ(projected_defaults->volumetric_samples, typed_defaults.volumetric_samples);
    EXPECT_EQ(projected_defaults->film_config, typed_defaults.film_config);

    SessionConfig config;
    config.lens_type = core::LensType::Pinhole;
    config.camera_focal_length = 85.0f;
    EXPECT_EQ(render::SessionConfigIssue(config),
              "focal length, aperture, and focus distance apply only to ThinLens");

    config.lens_type = core::LensType::ThinLens;
    EXPECT_FALSE(render::SessionConfigIssue(config).has_value());

    config = SessionConfig{};
    config.enable_disk = false;
    EXPECT_FALSE(render::SessionConfigIssue(config).has_value());
    config.temperature_model = render::DiskTemperatureModel::ShakuraSunyaev;
    EXPECT_EQ(render::SessionConfigIssue(config),
              "disk temperature model and scale require disk emission");

    config.temperature_model = render::DiskTemperatureModel::NovikovThorne;
    config.disk_temperature_scale = 42000.0f;
    EXPECT_EQ(render::SessionConfigIssue(config),
              "disk temperature model and scale require disk emission");

    config.disk_temperature_scale = core::kDefaultDiskTemperatureKelvin;
    config.doppler_beaming = false;
    EXPECT_EQ(render::SessionConfigIssue(config), "Doppler-beaming control requires disk emission");

    config.doppler_beaming = true;
    config.color_mode = core::color_modes::Mode::RedshiftMap;
    EXPECT_EQ(render::SessionConfigIssue(config), "diagnostic colour modes require disk emission");

    config = SessionConfig{};
    config.enable_bloom = false;
    config.bloom_threshold = 0.8f;
    EXPECT_EQ(render::SessionConfigIssue(config), "bloom intensity and threshold require bloom");

    config = SessionConfig{};
    config.volumetric_h_power = 0.5f;
    EXPECT_EQ(render::SessionConfigIssue(config),
              "volumetric parameters require volumetric transfer");

    config = SessionConfig{};
    config.shutter_time = 0.25f;
    EXPECT_EQ(render::SessionConfigIssue(config),
              "motion-blur parameters require temporal integration");

    config = SessionConfig{};
    config.film_config.grain_intensity = 0.4f;
    EXPECT_EQ(render::SessionConfigIssue(config),
              "film-simulation parameters require film simulation");

    config = SessionConfig{};
    config.starfield_config.brightness_scale = 200.0f;
    EXPECT_EQ(render::SessionConfigIssue(config),
              "point-starfield parameters require point-starfield mode");
}

TEST(RenderSessionProbe, BackendAutoResolvesByDeviceRegistryAndCapabilities) {
    SiriusConfig config = SiriusConfig::Defaults();
    config.metric.name = "Kerr";
    config.backend.preferred = "cpu";
    EXPECT_EQ(MakeSessionConfig(config)->backend, render::RenderBackend::Cpu)
        << "'cpu' must pin the CPU path";

    config.backend.preferred = "auto";
#ifdef SIRIUS_HAS_VULKAN_BACKEND
    const auto devices = backend::EnumerateVulkanDevices();
    const bool device_present = devices.has_value() && !devices->empty();
    EXPECT_EQ(MakeSessionConfig(config)->backend,
              device_present ? render::RenderBackend::Vulkan : render::RenderBackend::Cpu)
        << "auto must follow device presence when the full sampled scene is represented";

    config.volumetric.enabled = true;
    config.volumetric.samples = 129;
    EXPECT_EQ(MakeSessionConfig(config)->backend, render::RenderBackend::Cpu)
        << "auto must select CPU above the Vulkan volume-sample boundary";
    config.volumetric.enabled = false;
    config.volumetric.samples = 64;

    config.metric.name = "Reissner-Nordstrom";
    config.disk_enabled = false;
    EXPECT_EQ(MakeSessionConfig(config)->backend, render::RenderBackend::Cpu)
        << "auto must resolve CPU for a metric the registry marks CPU-only";
#else
    EXPECT_EQ(MakeSessionConfig(config)->backend, render::RenderBackend::Cpu)
        << "auto must resolve CPU when the Vulkan backend is not compiled in";
#endif
}

TEST(RenderSessionProbe, ConfigurationConversionPreservesObserverAndDiskControls) {
    SiriusConfig config = SiriusConfig::Defaults();
    config.backend.preferred = "cpu";
    config.observer.azimuth = 37.5;
    config.metric.temperature_model = "ShakuraSunyaev";
    config.metric.disk_temperature = 42000.0f;
    config.render.thread_count = 7;
    config.observer.camera_beta_forward = 0.1;
    config.observer.lens_model = "ThinLens";
    config.observer.focal_length = 85.0f;
    config.observer.aperture = 1.4f;
    config.observer.focus_distance = 40.0f;
    config.color_mode = "Polarisation";

    const auto adapted = MakeSessionConfig(config);
    ASSERT_TRUE(adapted.has_value()) << adapted.error().Description();
    const SessionConfig session = *adapted;
    EXPECT_NEAR(session.observer_azimuth, 37.5 * std::numbers::pi / 180.0, 1e-12);
    EXPECT_EQ(session.temperature_model, render::DiskTemperatureModel::ShakuraSunyaev);
    EXPECT_FLOAT_EQ(session.disk_temperature_scale, 42000.0f);
    EXPECT_EQ(session.thread_count, 7);
    EXPECT_DOUBLE_EQ(session.camera_beta_forward, 0.1);
    EXPECT_EQ(session.lens_type, core::LensType::ThinLens);
    EXPECT_FLOAT_EQ(session.camera_focal_length, 85.0f);
    EXPECT_FLOAT_EQ(session.camera_aperture, 1.4f);
    EXPECT_FLOAT_EQ(session.camera_focus_distance, 40.0f);
    EXPECT_EQ(session.color_mode, core::color_modes::Mode::Polarisation);
    EXPECT_TRUE(session.enable_polarisation);

    config.metric.temperature_model = "UnknownTemperature";
    EXPECT_FALSE(MakeSessionConfig(config).has_value());
    config.metric.temperature_model = "NovikovThorne";
    config.observer.lens_model = "UnknownLens";
    EXPECT_FALSE(MakeSessionConfig(config).has_value());
    config.observer.lens_model = "Pinhole";
    config.film.preset = "UnknownFilm";
    EXPECT_FALSE(MakeSessionConfig(config).has_value());
    config.film.preset = "Interstellar";
    config.observer.lens_model = "Fisheye";
    config.observer.focal_length = core::kDefaultCameraFocalLength;
    config.observer.aperture = core::kDefaultCameraAperture;
    config.observer.focus_distance = core::kDefaultCameraFocusDistance;
    EXPECT_EQ(MakeSessionConfig(config)->lens_type, core::LensType::Fisheye);

    SessionConfig malformed = session;
    malformed.width = 1;
    malformed.height = 1;
    malformed.tile_size = 1;
    malformed.samples_per_pixel = 1;
    malformed.enable_parallel_rendering = false;
    malformed.temperature_model = static_cast<render::DiskTemperatureModel>(255);
    RenderSession cpu_session;
    const auto invalid_temperature = cpu_session.Configure(malformed);
    ASSERT_FALSE(invalid_temperature);
    EXPECT_EQ(invalid_temperature.error().domain(), base::ErrorDomain::kConfiguration);
    EXPECT_NE(invalid_temperature.error().detail().find("temperature model"), std::string::npos);

    malformed.temperature_model = render::DiskTemperatureModel::NovikovThorne;
    malformed.color_mode = core::color_modes::Mode::TrueColor;
    malformed.enable_polarisation = true;
    RenderSession polarisation_session;
    EXPECT_FALSE(polarisation_session.Configure(malformed));

    malformed.color_mode = core::color_modes::Mode::Polarisation;
    malformed.enable_polarisation = false;
    RenderSession polarisation_colour_session;
    EXPECT_FALSE(polarisation_colour_session.Configure(malformed));

    malformed.color_mode = static_cast<core::color_modes::Mode>(255);
    RenderSession invalid_colour_session;
    EXPECT_FALSE(invalid_colour_session.Configure(malformed));

    malformed.color_mode = core::color_modes::Mode::TrueColor;
    malformed.tonemapper = static_cast<core::TonemapType>(255);
    RenderSession invalid_tonemapper_session;
    EXPECT_FALSE(invalid_tonemapper_session.Configure(malformed));

    malformed.tonemapper = core::TonemapType::Aces;
    malformed.metric_id = core::MetricId::Schwarzschild;
    malformed.black_hole_spin = 0.4;
    RenderSession mismatched_metric_session;
    EXPECT_FALSE(mismatched_metric_session.Configure(malformed));

    malformed.metric_id = core::MetricId::KerrNewman;
    malformed.black_hole_spin = 0.4;
    malformed.black_hole_charge = 0.2;
    malformed.enable_disk = true;
    RenderSession charged_disk_session;
    EXPECT_FALSE(charged_disk_session.Configure(malformed));

    malformed.metric_id = core::MetricId::Minkowski;
    malformed.black_hole_mass = 0.0;
    malformed.black_hole_spin = 0.0;
    malformed.black_hole_charge = 0.0;
    malformed.enable_disk = true;
    RenderSession inapplicable_disk_session;
    EXPECT_FALSE(inapplicable_disk_session.Configure(malformed));

    malformed.metric_id = core::MetricId::Schwarzschild;
    malformed.black_hole_mass = 1.0;
    malformed.black_hole_spin = 0.0;
    malformed.black_hole_charge = 0.0;
    malformed.enable_disk = true;
    malformed.width = 0;
    RenderSession invalid_dimensions_session;
    EXPECT_FALSE(invalid_dimensions_session.Configure(malformed));
}

TEST(RenderSessionProbe, InMemoryPreviewCompletesWithoutWritingOutput) {
    const auto preview_path =
        std::filesystem::temp_directory_path() / "sirius_in_memory_preview_must_not_exist.ppm";
    std::filesystem::remove(preview_path);
    SessionConfig preview;
    preview.width = 8;
    preview.height = 8;
    preview.tile_size = 8;
    preview.samples_per_pixel = 1;
    preview.enable_parallel_rendering = false;
    preview.metric_id = core::MetricId::Minkowski;
    preview.black_hole_mass = 0.0;
    preview.enable_disk = false;
    preview.write_output = false;
    preview.output_path = preview_path.string();
    RenderSession preview_session;
    ASSERT_TRUE(preview_session.Configure(preview));
    EXPECT_EQ(preview_session.Execute(), SessionState::Complete);
    EXPECT_FALSE(std::filesystem::exists(preview_path));
}

}  // namespace sirius::app::test
