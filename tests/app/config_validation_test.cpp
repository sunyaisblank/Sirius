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

TEST(ConfigValidation, VolumetricExtensionsRequireTheLiveVolumePath) {
    SiriusConfig config = SiriusConfig::Defaults();
    config.volumetric.enable_corona = true;
    config.volumetric.enable_turbulence = true;
    EXPECT_FALSE(ConfigLoader::Validate(config).empty());

    config.volumetric.enabled = true;
    EXPECT_TRUE(ConfigLoader::Validate(config).empty());

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

TEST(ConfigValidation, MasslessMetricUsesAUnitObserverDistanceScale) {
    SiriusConfig config = SiriusConfig::Defaults();
    config.metric.name = "Minkowski";
    config.metric.mass = 0.0;
    config.metric.spin = 0.0;
    config.disk_enabled = false;
    EXPECT_TRUE(ConfigLoader::Validate(config).empty());

    config.metric.name = "de-Sitter";
    config.metric.cosmological_constant = 0.001;
    EXPECT_TRUE(ConfigLoader::Validate(config).empty());

    config.metric.mass = 1.0;
    EXPECT_FALSE(ConfigLoader::Validate(config).empty());
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
    EXPECT_TRUE(ConfigLoader::Validate(config).empty());

    config.backend.preferred = "vulkan";
    // Backend-specific capability is checked after projection; the external
    // configuration itself remains coherent.
    EXPECT_TRUE(ConfigLoader::Validate(config).empty());

    config.motion_blur.samples = 0;
    EXPECT_FALSE(ConfigLoader::Validate(config).empty());
    config.motion_blur.samples = 7;
    config.disk_enabled = false;
    EXPECT_FALSE(ConfigLoader::Validate(config).empty());

    config = SiriusConfig::Defaults();
    config.metric.name = "Morris-Thorne";
    config.disk_enabled = false;
    config.metric.wormhole_topology = "OneSheetCapture";
    EXPECT_TRUE(ConfigLoader::Validate(config).empty());
    config.metric.wormhole_topology = "TwoSheet";
    const auto errors = ConfigLoader::Validate(config);
    ASSERT_FALSE(errors.empty());
    EXPECT_NE(std::find_if(errors.begin(), errors.end(),
                           [](const std::string& error) {
                               return error.find("TwoSheet is not represented") !=
                                      std::string::npos;
                           }),
              errors.end());
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

}  // namespace sirius::app::test
