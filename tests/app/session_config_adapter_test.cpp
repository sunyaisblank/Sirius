#include "sirius/app/config/session_config_adapter.h"

#include "sirius/render/session/render_session.h"

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
    config.motion_blur.enabled = false;
    config.motion_blur.shutter_time = 0.25f;
    config.motion_blur.samples = 7;

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
    EXPECT_FLOAT_EQ(session.shutter_time, 0.25f);
    EXPECT_EQ(session.motion_blur_samples, 7);
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
