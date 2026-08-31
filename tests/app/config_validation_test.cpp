// Configuration validation tests: the ConfigLoader::Validate contract accepts a
// default configuration and rejects each out-of-range or unknown field with a
// message. These exercise the CLI/config objects linked into this target (July
// remediation ledger item 8).

#include "sirius/app/config/config_loader.h"
#include "sirius/app/config/config_schema.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <limits>
#include <string_view>

namespace sirius::app::test {

TEST(ConfigValidation, DefaultConfigurationIsValid) {
    SiriusConfig config = SiriusConfig::Defaults();
    EXPECT_TRUE(ConfigLoader::Validate(config).empty());
}

TEST(ConfigValidation, WidthBelowMinimumRejected) {
    SiriusConfig config = SiriusConfig::Defaults();
    config.render.width = 64;  // Below the 128 floor.
    EXPECT_FALSE(ConfigLoader::Validate(config).empty());
}

TEST(ConfigValidation, NonPowerOfTwoTileSizeRejected) {
    SiriusConfig config = SiriusConfig::Defaults();
    config.render.tile_size = 48;  // In range but not a power of two.
    EXPECT_FALSE(ConfigLoader::Validate(config).empty());
}

TEST(ConfigValidation, UnknownMetricNameRejected) {
    SiriusConfig config = SiriusConfig::Defaults();
    config.metric.name = "Gargantua";
    EXPECT_FALSE(ConfigLoader::Validate(config).empty());
}

TEST(ConfigValidation, KnownMetricAliasAccepted) {
    SiriusConfig config = SiriusConfig::Defaults();
    config.metric.name = "Wormhole";  // Morris-Thorne alias in the registry.
    // The name is valid; validation does not gate on backend support (that is
    // the session's decline), so no metric-name error is produced.
    auto errors = ConfigLoader::Validate(config);
    for (const auto& e : errors) {
        EXPECT_EQ(e.find("is not a known metric"), std::string::npos) << e;
    }
}

TEST(ConfigValidation, SpinAboveNearExtremalRejected) {
    SiriusConfig config = SiriusConfig::Defaults();
    config.metric.name = "Kerr";
    config.metric.spin = 1.5;
    EXPECT_FALSE(ConfigLoader::Validate(config).empty());
}

TEST(ConfigValidation, RotatingLambdaRejected) {
    SiriusConfig config = SiriusConfig::Defaults();
    config.metric.spin = 0.5;
    config.metric.cosmological_constant = 0.01;  // Lambda requires spin = 0.
    EXPECT_FALSE(ConfigLoader::Validate(config).empty());
}

TEST(ConfigValidation, ObserverInsidePoleBufferRejected) {
    SiriusConfig config = SiriusConfig::Defaults();
    config.observer.inclination = 0.0;  // At the pole.
    EXPECT_FALSE(ConfigLoader::Validate(config).empty());
}

TEST(ConfigValidation, UnknownTonemapperRejected) {
    SiriusConfig config = SiriusConfig::Defaults();
    config.postprocess.tonemapper = "Cineon";
    EXPECT_FALSE(ConfigLoader::Validate(config).empty());
}

TEST(ConfigValidation, AcesFitIsExplicitAndBareAcesIsRejected) {
    SiriusConfig config = SiriusConfig::Defaults();
    EXPECT_EQ(config.postprocess.tonemapper, "ACESFit");
    EXPECT_TRUE(ConfigLoader::Validate(config).empty());

    config.postprocess.tonemapper = "ACES";
    const auto errors = ConfigLoader::Validate(config);
    ASSERT_FALSE(errors.empty());
    EXPECT_NE(errors.front().find("ACESFit"), std::string::npos);
}

TEST(ConfigValidation, VulkanBackendNameAccepted) {
    // Go-live (2026-07-18): vulkan joined the accepted set; device absence is
    // the render's decline, not a config-name error.
    SiriusConfig config = SiriusConfig::Defaults();
    config.backend.preferred = "vulkan";
    for (const auto& e : ConfigLoader::Validate(config)) {
        EXPECT_EQ(e.find("backend.preferred"), std::string::npos) << e;
    }
}

TEST(ConfigValidation, UnknownBackendNameRejected) {
    SiriusConfig config = SiriusConfig::Defaults();
    config.backend.preferred = "optix";  // Retired; never silently remapped.
    EXPECT_FALSE(ConfigLoader::Validate(config).empty());
}

TEST(ConfigValidation, CpuBackendAccepted) {
    SiriusConfig config = SiriusConfig::Defaults();
    config.backend.preferred = "cpu";
    EXPECT_TRUE(ConfigLoader::Validate(config).empty());
}

TEST(ConfigValidation, NonFiniteValuesAreRejected) {
    SiriusConfig config = SiriusConfig::Defaults();
    config.metric.spin = std::numeric_limits<double>::quiet_NaN();
    config.postprocess.exposure = std::numeric_limits<float>::infinity();
    EXPECT_FALSE(ConfigLoader::Validate(config).empty());
}

TEST(ConfigValidation, UnknownOutputExtensionIsRejected) {
    SiriusConfig config = SiriusConfig::Defaults();
    config.render.output_path = "render.png.exe";
    EXPECT_FALSE(ConfigLoader::Validate(config).empty());
}

TEST(ConfigValidation, UnrepresentedNarrowbandEmissionDeclines) {
    SiriusConfig config = SiriusConfig::Defaults();
    config.color_mode = "Narrowband";
    const auto errors = ConfigLoader::Validate(config);
    EXPECT_NE(std::find_if(errors.begin(), errors.end(),
                           [](const std::string& error) {
                               return error.find("ionisation") != std::string::npos &&
                                      error.find("frequency-dependent transfer") !=
                                          std::string::npos;
                           }),
              errors.end());
}

TEST(ConfigValidation, TurbulenceRequiresVolumeAndSpectralCoronaDeclines) {
    SiriusConfig config = SiriusConfig::Defaults();
    config.volumetric.enable_turbulence = true;
    EXPECT_FALSE(ConfigLoader::Validate(config).empty());

    config.volumetric.enabled = true;
    EXPECT_TRUE(ConfigLoader::Validate(config).empty());

    config.volumetric.enable_corona = true;
    const auto corona_errors = ConfigLoader::Validate(config);
    EXPECT_NE(std::find_if(corona_errors.begin(), corona_errors.end(),
                           [](const std::string& error) {
                               return error.find("covariant Compton transfer") != std::string::npos;
                           }),
              corona_errors.end());
    config.volumetric.enable_corona = false;

    config.disk_enabled = false;
    EXPECT_FALSE(ConfigLoader::Validate(config).empty());
}

TEST(ConfigValidation, UnimplementedDenoiserRequestIsRejected) {
    SiriusConfig config = SiriusConfig::Defaults();
    config.backend.enable_denoiser = true;
    EXPECT_FALSE(ConfigLoader::Validate(config).empty());
}

TEST(ConfigValidation, ObserverAzimuthRangeIsValidated) {
    SiriusConfig config = SiriusConfig::Defaults();
    config.observer.azimuth = 361.0;
    EXPECT_FALSE(ConfigLoader::Validate(config).empty());
}

TEST(ConfigValidation, CameraWorldlineAndLensAreValidated) {
    SiriusConfig config = SiriusConfig::Defaults();
    config.observer.camera_beta_forward = 1.0;
    config.observer.lens_model = "UnknownLens";
    config.observer.aperture = 0.0f;
    const auto errors = ConfigLoader::Validate(config);
    EXPECT_GE(errors.size(), 3u);

    config = SiriusConfig::Defaults();
    config.observer.camera_beta_forward = std::numeric_limits<double>::max();
    EXPECT_FALSE(ConfigLoader::Validate(config).empty());

    config = SiriusConfig::Defaults();
    config.observer.camera_beta_forward = 0.2;
    config.observer.lens_model = "ThinLens";
    config.observer.aperture = 1.4f;
    EXPECT_TRUE(ConfigLoader::Validate(config).empty());

    config.observer.lens_model = "Fisheye";
    EXPECT_FALSE(ConfigLoader::Validate(config).empty());
    config.observer.aperture = core::kDefaultCameraAperture;
    EXPECT_TRUE(ConfigLoader::Validate(config).empty());

    config = SiriusConfig::Defaults();
    config.metric.name = "Schwarzschild";
    config.metric.spin = 0.4;
    EXPECT_FALSE(ConfigLoader::Validate(config).empty());

    config = SiriusConfig::Defaults();
    config.metric.name = "Kerr-Newman";
    config.metric.spin = 0.4;
    config.metric.charge = 0.2;
    EXPECT_FALSE(ConfigLoader::Validate(config).empty());
    config.disk_enabled = false;
    EXPECT_TRUE(ConfigLoader::Validate(config).empty());
}

TEST(ConfigValidation, FeatureSpecificControlsRequireTheirOwningModels) {
    SiriusConfig config = SiriusConfig::Defaults();
    config.observer.focal_length = 85.0f;
    auto errors = ConfigLoader::Validate(config);
    EXPECT_NE(std::find_if(errors.begin(), errors.end(),
                           [](const std::string& error) {
                               return error.find("only to ThinLens") != std::string::npos;
                           }),
              errors.end());

    config.observer.lens_model = "ThinLens";
    EXPECT_TRUE(ConfigLoader::Validate(config).empty());

    config = SiriusConfig::Defaults();
    config.disk_enabled = false;
    EXPECT_TRUE(ConfigLoader::Validate(config).empty());

    config.metric.temperature_model = "ShakuraSunyaev";
    errors = ConfigLoader::Validate(config);
    EXPECT_NE(std::find_if(errors.begin(), errors.end(),
                           [](const std::string& error) {
                               return error.find("temperature model") != std::string::npos;
                           }),
              errors.end());

    config.metric.temperature_model = "NovikovThorne";
    config.metric.disk_temperature = 42000.0f;
    EXPECT_FALSE(ConfigLoader::Validate(config).empty());

    config.metric.disk_temperature = core::kDefaultDiskTemperatureKelvin;
    config.doppler_beaming = false;
    errors = ConfigLoader::Validate(config);
    EXPECT_NE(std::find_if(errors.begin(), errors.end(),
                           [](const std::string& error) {
                               return error.find("Doppler") != std::string::npos ||
                                      error.find("doppler") != std::string::npos;
                           }),
              errors.end());

    config.doppler_beaming = true;
    config.color_mode = "RedshiftMap";
    errors = ConfigLoader::Validate(config);
    EXPECT_NE(std::find_if(errors.begin(), errors.end(),
                           [](const std::string& error) {
                               return error.find("diagnostic color modes") != std::string::npos;
                           }),
              errors.end());

    config = SiriusConfig::Defaults();
    config.postprocess.enable_bloom = false;
    config.postprocess.bloom_threshold = 0.8f;
    EXPECT_FALSE(ConfigLoader::Validate(config).empty());

    config = SiriusConfig::Defaults();
    config.volumetric.h_power = 0.5f;
    EXPECT_FALSE(ConfigLoader::Validate(config).empty());

    config = SiriusConfig::Defaults();
    config.motion_blur.shutter_time = 0.25f;
    EXPECT_FALSE(ConfigLoader::Validate(config).empty());

    config = SiriusConfig::Defaults();
    config.film.grain_intensity = 0.4f;
    EXPECT_FALSE(ConfigLoader::Validate(config).empty());
}

TEST(ConfigValidation, MetricMassAndObserverCoordinateRadiusAreIdentityAware) {
    for (const std::string_view name : {"Minkowski", "de-Sitter", "Morris-Thorne", "Alcubierre"}) {
        SiriusConfig config = SiriusConfig::Defaults();
        config.metric.name = name;
        config.metric.mass = 0.0;
        config.metric.cosmological_constant = name == "de-Sitter" ? 0.001 : 0.0;
        config.disk_enabled = false;
        EXPECT_TRUE(ConfigLoader::Validate(config).empty()) << name;

        config.metric.mass = 1.0;
        const auto mass_errors = ConfigLoader::Validate(config);
        EXPECT_NE(std::find_if(mass_errors.begin(), mass_errors.end(),
                               [](const std::string& error) {
                                   return error.find("without a mass parameter") !=
                                          std::string::npos;
                               }),
                  mass_errors.end())
            << name;

        config.metric.mass = 0.0;
        config.observer.distance = 4.0;
        const auto distance_errors = ConfigLoader::Validate(config);
        const std::string_view expected_distance_text = name == "Morris-Thorne" ? "5*b0"
                                                        : name == "Alcubierre"
                                                            ? "5*L"
                                                            : "geometric coordinate units";
        EXPECT_NE(std::find_if(distance_errors.begin(), distance_errors.end(),
                               [expected_distance_text](const std::string& error) {
                                   return error.find(expected_distance_text) != std::string::npos;
                               }),
                  distance_errors.end())
            << name;
    }

    SiriusConfig scaled = SiriusConfig::Defaults();
    scaled.metric.mass = 10.0;
    scaled.observer.distance = 49.0;
    const auto below_errors = ConfigLoader::Validate(scaled);
    EXPECT_NE(std::find_if(below_errors.begin(), below_errors.end(),
                           [](const std::string& error) {
                               return error.find("5*M <= r <= 1000*M") != std::string::npos;
                           }),
              below_errors.end());
    scaled.observer.distance = 50.0;
    EXPECT_TRUE(ConfigLoader::Validate(scaled).empty());
    scaled.observer.distance = 10000.0;
    EXPECT_TRUE(ConfigLoader::Validate(scaled).empty());
    scaled.observer.distance = 10001.0;
    EXPECT_FALSE(ConfigLoader::Validate(scaled).empty());

    SiriusConfig metric_specific = SiriusConfig::Defaults();
    metric_specific.metric.throat_radius = 2.0;
    auto specific_errors = ConfigLoader::Validate(metric_specific);
    EXPECT_NE(std::find_if(specific_errors.begin(), specific_errors.end(),
                           [](const std::string& error) {
                               return error.find("only to Morris-Thorne") != std::string::npos;
                           }),
              specific_errors.end());

    metric_specific = SiriusConfig::Defaults();
    metric_specific.metric.name = "Morris-Thorne";
    metric_specific.metric.mass = 0.0;
    metric_specific.metric.throat_radius = 2.0;
    metric_specific.disk_enabled = false;
    EXPECT_TRUE(ConfigLoader::Validate(metric_specific).empty());
    metric_specific.metric.throat_radius = 0.09;
    EXPECT_FALSE(ConfigLoader::Validate(metric_specific).empty());
    metric_specific.metric.throat_radius = 20.0;
    metric_specific.observer.distance = 99.0;
    EXPECT_FALSE(ConfigLoader::Validate(metric_specific).empty());
    metric_specific.observer.distance = 100.0;
    EXPECT_TRUE(ConfigLoader::Validate(metric_specific).empty());
    metric_specific.observer.distance = 20000.0;
    EXPECT_TRUE(ConfigLoader::Validate(metric_specific).empty());
    metric_specific.observer.distance = 20001.0;
    EXPECT_FALSE(ConfigLoader::Validate(metric_specific).empty());
    metric_specific.metric.throat_radius = 2.0;
    metric_specific.observer.distance = 50.0;
    metric_specific.metric.warp_velocity = 1.0;
    specific_errors = ConfigLoader::Validate(metric_specific);
    EXPECT_NE(std::find_if(specific_errors.begin(), specific_errors.end(),
                           [](const std::string& error) {
                               return error.find("only to Alcubierre") != std::string::npos;
                           }),
              specific_errors.end());

    metric_specific = SiriusConfig::Defaults();
    metric_specific.metric.name = "Alcubierre";
    metric_specific.metric.mass = 0.0;
    metric_specific.metric.warp_velocity = 1.0;
    metric_specific.metric.bubble_radius = 2.0;
    metric_specific.disk_enabled = false;
    EXPECT_TRUE(ConfigLoader::Validate(metric_specific).empty());
    metric_specific.metric.bubble_sigma = 0.05;  // sigma*R=0.1, L=20.
    metric_specific.observer.distance = 99.0;
    EXPECT_FALSE(ConfigLoader::Validate(metric_specific).empty());
    metric_specific.observer.distance = 100.0;
    EXPECT_TRUE(ConfigLoader::Validate(metric_specific).empty());
    metric_specific.metric.bubble_sigma = 0.04;
    auto scale_errors = ConfigLoader::Validate(metric_specific);
    EXPECT_NE(std::find_if(scale_errors.begin(), scale_errors.end(),
                           [](const std::string& error) {
                               return error.find("bubble_sigma*bubble_radius") != std::string::npos;
                           }),
              scale_errors.end());
    metric_specific.metric.bubble_sigma = 51.0;
    scale_errors = ConfigLoader::Validate(metric_specific);
    EXPECT_NE(std::find_if(scale_errors.begin(), scale_errors.end(),
                           [](const std::string& error) {
                               return error.find("bubble_sigma*bubble_radius") != std::string::npos;
                           }),
              scale_errors.end());
    metric_specific.metric.bubble_sigma = 0.5;
    metric_specific.observer.distance = 50.0;
    metric_specific.metric.throat_radius = 2.0;
    EXPECT_FALSE(ConfigLoader::Validate(metric_specific).empty());
}

TEST(ConfigValidation, DeSitterRequestsEnforcePositiveLambdaAndSubNariaiBlackHole) {
    SiriusConfig config = SiriusConfig::Defaults();
    config.disk_enabled = false;
    config.metric.name = "de-Sitter";
    config.metric.mass = 0.0;

    config.metric.cosmological_constant = 0.0;
    EXPECT_FALSE(ConfigLoader::Validate(config).empty());
    config.metric.cosmological_constant = -0.001;
    EXPECT_FALSE(ConfigLoader::Validate(config).empty());
    config.metric.cosmological_constant = 0.001;
    EXPECT_TRUE(ConfigLoader::Validate(config).empty());
    config.observer.distance = 55.0;
    auto errors = ConfigLoader::Validate(config);
    EXPECT_NE(std::find_if(errors.begin(), errors.end(),
                           [](const std::string& error) {
                               return error.find("cosmological horizon") != std::string::npos;
                           }),
              errors.end());
    config.observer.distance = 50.0;

    config.metric.name = "Schwarzschild-de-Sitter";
    config.metric.mass = 2.0;
    config.metric.cosmological_constant = 0.001;
    EXPECT_TRUE(ConfigLoader::Validate(config).empty());

    config.metric.cosmological_constant = 0.02;  // Sub-Nariai, but r=50 lies outside r_c.
    errors = ConfigLoader::Validate(config);
    EXPECT_NE(std::find_if(errors.begin(), errors.end(),
                           [](const std::string& error) {
                               return error.find("cosmological horizon") != std::string::npos;
                           }),
              errors.end());

    config.metric.cosmological_constant = 0.03;  // 9 Lambda M^2 = 1.08.
    errors = ConfigLoader::Validate(config);
    EXPECT_NE(std::find_if(errors.begin(), errors.end(),
                           [](const std::string& error) {
                               return error.find("sub-Nariai") != std::string::npos;
                           }),
              errors.end());
}

TEST(ConfigValidation, PolarisationRequiresRepresentedThinBlackHoleDisk) {
    SiriusConfig config = SiriusConfig::Defaults();
    config.color_mode = "Polarisation";
    config.metric.name = "Kerr";
    config.metric.spin = 0.7;
    EXPECT_TRUE(ConfigLoader::Validate(config).empty());

    config.volumetric.enabled = true;
    EXPECT_FALSE(ConfigLoader::Validate(config).empty());

    config.volumetric.enabled = false;
    config.motion_blur.enabled = true;
    EXPECT_FALSE(ConfigLoader::Validate(config).empty());

    config.motion_blur.enabled = false;
    config.disk_enabled = false;
    EXPECT_FALSE(ConfigLoader::Validate(config).empty());

    config.disk_enabled = true;
    config.metric.name = "Morris-Thorne";
    config.metric.spin = 0.0;
    EXPECT_FALSE(ConfigLoader::Validate(config).empty());

    config = SiriusConfig::Defaults();
    config.color_mode = "UnverifiedPolarMap";
    EXPECT_FALSE(ConfigLoader::Validate(config).empty());
}

TEST(ConfigValidation, MotionBlurAndWormholeTopologyHaveExplicitOperatorBoundaries) {
    SiriusConfig config = SiriusConfig::Defaults();
    config.metric.name = "Kerr";
    config.metric.spin = 0.7;
    config.motion_blur.enabled = true;
    config.motion_blur.shutter_time = 0.25f;
    config.motion_blur.samples = 7;
    auto errors = ConfigLoader::Validate(config);
    EXPECT_FALSE(errors.empty());
    EXPECT_NE(std::find_if(errors.begin(), errors.end(),
                           [](const std::string& error) {
                               return error.find("temporal emissivity model") != std::string::npos;
                           }),
              errors.end());

    config.backend.preferred = "vulkan";
    EXPECT_FALSE(ConfigLoader::Validate(config).empty());

    config.motion_blur.samples = 0;
    EXPECT_FALSE(ConfigLoader::Validate(config).empty());
    config.motion_blur.samples = 7;
    config.disk_enabled = false;
    EXPECT_FALSE(ConfigLoader::Validate(config).empty());

    config = SiriusConfig::Defaults();
    config.metric.name = "Morris-Thorne";
    config.metric.mass = 0.0;
    config.disk_enabled = false;
    config.metric.wormhole_topology = "OneSheetCapture";
    EXPECT_TRUE(ConfigLoader::Validate(config).empty());
    config.metric.wormhole_topology = "TwoSheet";
    const auto topology_errors = ConfigLoader::Validate(config);
    ASSERT_FALSE(topology_errors.empty());
    EXPECT_NE(std::find_if(topology_errors.begin(), topology_errors.end(),
                           [](const std::string& error) {
                               return error.find("TwoSheet is not represented") !=
                                      std::string::npos;
                           }),
              topology_errors.end());
}

TEST(ConfigValidation, DiskRequestDeclinesForEveryMetricWithoutAnEmissionModel) {
    SiriusConfig config = SiriusConfig::Defaults();
    config.disk_enabled = true;
    config.metric.spin = 0.0;
    for (const std::string_view metric :
         {"Minkowski", "de-Sitter", "Morris-Thorne", "Alcubierre", "Reissner-Nordstrom",
          "Kerr-Newman", "Schwarzschild-de-Sitter"}) {
        config.metric.name = metric;
        EXPECT_FALSE(ConfigLoader::Validate(config).empty()) << metric;
    }

    config.metric.name = "Schwarzschild";
    EXPECT_TRUE(ConfigLoader::Validate(config).empty());
    config.metric.name = "Kerr";
    EXPECT_TRUE(ConfigLoader::Validate(config).empty());
}

TEST(ConfigValidation, RayBundlesRejectMetricsWithoutStationaryCurvatureTransport) {
    SiriusConfig config = SiriusConfig::Defaults();
    config.disk_enabled = false;
    config.ray_bundles = true;

    for (const std::string_view metric :
         {"de-Sitter", "Morris-Thorne", "Alcubierre", "Reissner-Nordstrom", "Kerr-Newman",
          "Schwarzschild-de-Sitter"}) {
        config.metric.name = metric;
        const auto errors = ConfigLoader::Validate(config);
        EXPECT_NE(std::find_if(errors.begin(), errors.end(),
                               [](const std::string& error) {
                                   return error.find("covariant curvature transport") !=
                                          std::string::npos;
                               }),
                  errors.end())
            << metric;
    }

    for (const std::string_view metric : {"Minkowski", "Schwarzschild", "Kerr"}) {
        config.metric.name = metric;
        config.metric.mass = metric == "Minkowski" ? 0.0 : 1.0;
        config.metric.spin = metric == "Kerr" ? 0.7 : 0.0;
        EXPECT_TRUE(ConfigLoader::Validate(config).empty()) << metric;
    }
}

}  // namespace sirius::app::test
