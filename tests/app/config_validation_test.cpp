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
    SiriusConfig config = SiriusConfig::defaults();
    EXPECT_TRUE(ConfigLoader::Validate(config).empty());
}

TEST(ConfigValidation, WidthBelowMinimumRejected) {
    SiriusConfig config = SiriusConfig::defaults();
    config.render.width = 64;  // Below the 128 floor.
    EXPECT_FALSE(ConfigLoader::Validate(config).empty());
}

TEST(ConfigValidation, NonPowerOfTwoTileSizeRejected) {
    SiriusConfig config = SiriusConfig::defaults();
    config.render.tileSize = 48;  // In range but not a power of two.
    EXPECT_FALSE(ConfigLoader::Validate(config).empty());
}

TEST(ConfigValidation, UnknownMetricNameRejected) {
    SiriusConfig config = SiriusConfig::defaults();
    config.metric.name = "Gargantua";
    EXPECT_FALSE(ConfigLoader::Validate(config).empty());
}

TEST(ConfigValidation, KnownMetricAliasAccepted) {
    SiriusConfig config = SiriusConfig::defaults();
    config.metric.name = "Wormhole";  // Morris-Thorne alias in the registry.
    // The name is valid; validation does not gate on backend support (that is
    // the session's decline), so no metric-name error is produced.
    auto errors = ConfigLoader::Validate(config);
    for (const auto& e : errors) {
        EXPECT_EQ(e.find("is not a known metric"), std::string::npos) << e;
    }
}

TEST(ConfigValidation, SpinAboveNearExtremalRejected) {
    SiriusConfig config = SiriusConfig::defaults();
    config.metric.name = "Kerr";
    config.metric.spin = 1.5;
    EXPECT_FALSE(ConfigLoader::Validate(config).empty());
}

TEST(ConfigValidation, RotatingLambdaRejected) {
    SiriusConfig config = SiriusConfig::defaults();
    config.metric.spin = 0.5;
    config.metric.lambda = 0.01;  // Lambda requires spin = 0.
    EXPECT_FALSE(ConfigLoader::Validate(config).empty());
}

TEST(ConfigValidation, ObserverInsidePoleBufferRejected) {
    SiriusConfig config = SiriusConfig::defaults();
    config.observer.inclination = 0.0;  // At the pole.
    EXPECT_FALSE(ConfigLoader::Validate(config).empty());
}

TEST(ConfigValidation, UnknownTonemapperRejected) {
    SiriusConfig config = SiriusConfig::defaults();
    config.postprocess.tonemapper = "Cineon";
    EXPECT_FALSE(ConfigLoader::Validate(config).empty());
}

TEST(ConfigValidation, VulkanBackendNameAccepted) {
    // Go-live (2026-07-18): vulkan joined the accepted set; device absence is
    // the render's decline, not a config-name error.
    SiriusConfig config = SiriusConfig::defaults();
    config.backend.preferred = "vulkan";
    for (const auto& e : ConfigLoader::Validate(config)) {
        EXPECT_EQ(e.find("backend.preferred"), std::string::npos) << e;
    }
}

TEST(ConfigValidation, UnknownBackendNameRejected) {
    SiriusConfig config = SiriusConfig::defaults();
    config.backend.preferred = "optix";  // Retired; never silently remapped.
    EXPECT_FALSE(ConfigLoader::Validate(config).empty());
}

TEST(ConfigValidation, CpuBackendAccepted) {
    SiriusConfig config = SiriusConfig::defaults();
    config.backend.preferred = "cpu";
    EXPECT_TRUE(ConfigLoader::Validate(config).empty());
}

TEST(ConfigValidation, NonFiniteValuesAreRejected) {
    SiriusConfig config = SiriusConfig::defaults();
    config.metric.spin = std::numeric_limits<double>::quiet_NaN();
    config.postprocess.exposure = std::numeric_limits<float>::infinity();
    EXPECT_FALSE(ConfigLoader::Validate(config).empty());
}

TEST(ConfigValidation, UnknownOutputExtensionIsRejected) {
    SiriusConfig config = SiriusConfig::defaults();
    config.render.outputPath = "render.png.exe";
    EXPECT_FALSE(ConfigLoader::Validate(config).empty());
}

TEST(ConfigValidation, VolumetricExtensionsRequireTheLiveVolumePath) {
    SiriusConfig config = SiriusConfig::defaults();
    config.volumetric.enableCorona = true;
    config.volumetric.enableTurbulence = true;
    EXPECT_FALSE(ConfigLoader::Validate(config).empty());

    config.volumetric.enabled = true;
    EXPECT_TRUE(ConfigLoader::Validate(config).empty());

    config.diskEnabled = false;
    EXPECT_FALSE(ConfigLoader::Validate(config).empty());
}

TEST(ConfigValidation, UnimplementedDenoiserRequestIsRejected) {
    SiriusConfig config = SiriusConfig::defaults();
    config.backend.enableDenoiser = true;
    EXPECT_FALSE(ConfigLoader::Validate(config).empty());
}

TEST(ConfigValidation, ObserverAzimuthRangeIsValidated) {
    SiriusConfig config = SiriusConfig::defaults();
    config.observer.azimuth = 361.0;
    EXPECT_FALSE(ConfigLoader::Validate(config).empty());
}

TEST(ConfigValidation, CameraWorldlineAndLensAreValidated) {
    SiriusConfig config = SiriusConfig::defaults();
    config.observer.cameraBetaForward = 1.0;
    config.observer.lensModel = "UnknownLens";
    config.observer.aperture = 0.0f;
    const auto errors = ConfigLoader::Validate(config);
    EXPECT_GE(errors.size(), 3u);

    config = SiriusConfig::defaults();
    config.observer.cameraBetaForward = std::numeric_limits<double>::max();
    EXPECT_FALSE(ConfigLoader::Validate(config).empty());

    config = SiriusConfig::defaults();
    config.observer.cameraBetaForward = 0.2;
    config.observer.lensModel = "ThinLens";
    config.observer.aperture = 1.4f;
    EXPECT_TRUE(ConfigLoader::Validate(config).empty());

    config.observer.lensModel = "Fisheye";
    EXPECT_TRUE(ConfigLoader::Validate(config).empty());

    config = SiriusConfig::defaults();
    config.metric.name = "Schwarzschild";
    config.metric.spin = 0.4;
    EXPECT_FALSE(ConfigLoader::Validate(config).empty());

    config = SiriusConfig::defaults();
    config.metric.name = "Kerr-Newman";
    config.metric.spin = 0.4;
    config.metric.charge = 0.2;
    EXPECT_FALSE(ConfigLoader::Validate(config).empty());
    config.diskEnabled = false;
    EXPECT_TRUE(ConfigLoader::Validate(config).empty());
}

TEST(ConfigValidation, MasslessMetricUsesAUnitObserverDistanceScale) {
    SiriusConfig config = SiriusConfig::defaults();
    config.metric.name = "Minkowski";
    config.metric.mass = 0.0;
    config.metric.spin = 0.0;
    config.diskEnabled = false;
    EXPECT_TRUE(ConfigLoader::Validate(config).empty());

    config.metric.name = "de-Sitter";
    config.metric.lambda = 0.001;
    EXPECT_TRUE(ConfigLoader::Validate(config).empty());

    config.metric.mass = 1.0;
    EXPECT_FALSE(ConfigLoader::Validate(config).empty());
}

TEST(ConfigValidation, PolarisationRequiresRepresentedThinBlackHoleDisk) {
    SiriusConfig config = SiriusConfig::defaults();
    config.colorMode = "Polarisation";
    config.metric.name = "Kerr";
    config.metric.spin = 0.7;
    EXPECT_TRUE(ConfigLoader::Validate(config).empty());

    config.volumetric.enabled = true;
    EXPECT_FALSE(ConfigLoader::Validate(config).empty());

    config.volumetric.enabled = false;
    config.motionBlur.enabled = true;
    EXPECT_FALSE(ConfigLoader::Validate(config).empty());

    config.motionBlur.enabled = false;
    config.diskEnabled = false;
    EXPECT_FALSE(ConfigLoader::Validate(config).empty());

    config.diskEnabled = true;
    config.metric.name = "Morris-Thorne";
    config.metric.spin = 0.0;
    EXPECT_FALSE(ConfigLoader::Validate(config).empty());

    config = SiriusConfig::defaults();
    config.colorMode = "UnverifiedPolarMap";
    EXPECT_FALSE(ConfigLoader::Validate(config).empty());
}

TEST(ConfigValidation, MotionBlurAndWormholeTopologyHaveExplicitOperatorBoundaries) {
    SiriusConfig config = SiriusConfig::defaults();
    config.metric.name = "Kerr";
    config.metric.spin = 0.7;
    config.motionBlur.enabled = true;
    config.motionBlur.shutterTime = 0.25f;
    config.motionBlur.samples = 7;
    EXPECT_TRUE(ConfigLoader::Validate(config).empty());

    config.backend.preferred = "vulkan";
    // Backend-specific capability is checked after projection; the external
    // configuration itself remains coherent.
    EXPECT_TRUE(ConfigLoader::Validate(config).empty());

    config.motionBlur.samples = 0;
    EXPECT_FALSE(ConfigLoader::Validate(config).empty());
    config.motionBlur.samples = 7;
    config.diskEnabled = false;
    EXPECT_FALSE(ConfigLoader::Validate(config).empty());

    config = SiriusConfig::defaults();
    config.metric.name = "Morris-Thorne";
    config.diskEnabled = false;
    config.metric.wormholeTopology = "OneSheetCapture";
    EXPECT_TRUE(ConfigLoader::Validate(config).empty());
    config.metric.wormholeTopology = "TwoSheet";
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
    SiriusConfig config = SiriusConfig::defaults();
    config.diskEnabled = true;
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
