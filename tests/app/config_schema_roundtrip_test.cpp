// Configuration JSON tests: the non-intrusive serialisation re-attached in
// config_schema.h round-trips a SiriusConfig and parses the exact legacy field
// spellings, so existing config files load identically. This is the schema
// compatibility check for the app/config port.

#include "sirius/app/config/config_schema.h"

#include <gtest/gtest.h>

#include <nlohmann/json.hpp>

namespace sirius::app::test {

TEST(ConfigSchema, DefaultsRoundTripThroughJson) {
    SiriusConfig original = SiriusConfig::Defaults();
    original.metric.name = "Kerr";
    original.metric.spin = 0.9;
    original.render.samples_per_pixel = 128;
    original.postprocess.tonemapper = "Filmic";
    original.observer.camera_beta_forward = 0.2;
    original.observer.lens_model = "ThinLens";
    original.metric.wormhole_topology = "TwoSheet";
    original.motion_blur.enabled = true;
    original.motion_blur.shutter_time = 0.25f;
    original.motion_blur.samples = 7;
    original.disk_enabled = false;
    original.doppler_beaming = false;
    original.point_starfield = true;
    original.ray_bundles = true;
    original.color_mode = "Polarisation";

    nlohmann::json j = original;
    SiriusConfig restored = j.get<SiriusConfig>();

    EXPECT_EQ(restored.metric.name, "Kerr");
    EXPECT_DOUBLE_EQ(restored.metric.spin, 0.9);
    EXPECT_EQ(restored.render.samples_per_pixel, 128);
    EXPECT_EQ(restored.postprocess.tonemapper, "Filmic");
    EXPECT_DOUBLE_EQ(restored.observer.camera_beta_forward, 0.2);
    EXPECT_EQ(restored.observer.lens_model, "ThinLens");
    EXPECT_EQ(restored.metric.wormhole_topology, "TwoSheet");
    EXPECT_TRUE(restored.motion_blur.enabled);
    EXPECT_FLOAT_EQ(restored.motion_blur.shutter_time, 0.25f);
    EXPECT_EQ(restored.motion_blur.samples, 7);
    EXPECT_FALSE(restored.disk_enabled);
    EXPECT_FALSE(restored.doppler_beaming);
    EXPECT_TRUE(restored.point_starfield);
    EXPECT_TRUE(restored.ray_bundles);
    EXPECT_EQ(restored.color_mode, "Polarisation");
}

TEST(ConfigSchema, LegacyFieldSpellingsParse) {
    // The exact key spellings the legacy CRCF003A schema used; a config file
    // written by the old binary must parse into the same fields.
    const char* legacy = R"({
        "render": { "width": 800, "height": 600, "samplesPerPixel": 16, "tileSize": 32,
                    "threadCount": 4, "outputPath": "out.exr" },
        "metric": { "name": "Kerr", "mass": 2.0, "spin": 0.75, "charge": 0.1,
                    "lambda": 0.0, "temperatureModel": "ShakuraSunyaev",
                    "diskTemperature": 60000.0, "throatRadius": 1.5,
                    "wormholeTopology": "OneSheetCapture",
                    "warpVelocity": 0.8, "bubbleRadius": 2.0, "bubbleSigma": 0.25 },
        "observer": { "distance": 100.0, "inclination": 60.0, "azimuth": 30.0, "fov": 45.0 },
        "postprocess": { "enableBloom": false, "bloomIntensity": 0.5, "bloomThreshold": 0.2,
                         "exposure": 2.0, "contrast": 1.2, "saturation": 0.9,
                         "tonemapper": "Reinhard" },
        "volumetric": { "enabled": true, "hOverR": 0.15, "hPower": 0.3, "tauMidplane": 12.0,
                        "samples": 96, "enableTurbulence": true, "enableCorona": false },
        "motionBlur": { "enabled": true, "shutterTime": 0.5, "samples": 9 },
        "film": { "enabled": true, "preset": "SpaceOdyssey2001", "grainIntensity": 0.2,
                  "halationStrength": 0.1, "vignetteStrength": 0.4 },
        "backend": { "preferred": "cpu", "enableDenoiser": false, "cudaDevice": 0 }
    })";

    SiriusConfig config = nlohmann::json::parse(legacy).get<SiriusConfig>();

    EXPECT_EQ(config.render.width, 800);
    EXPECT_EQ(config.render.output_path, "out.exr");
    EXPECT_EQ(config.metric.name, "Kerr");
    EXPECT_DOUBLE_EQ(config.metric.mass, 2.0);
    EXPECT_EQ(config.metric.temperature_model, "ShakuraSunyaev");
    EXPECT_FLOAT_EQ(config.metric.disk_temperature, 60000.0f);
    EXPECT_DOUBLE_EQ(config.observer.azimuth, 30.0);
    EXPECT_FALSE(config.postprocess.enable_bloom);
    EXPECT_EQ(config.postprocess.tonemapper, "Reinhard");
    EXPECT_TRUE(config.volumetric.enable_turbulence);
    EXPECT_TRUE(config.motion_blur.enabled);
    EXPECT_EQ(config.motion_blur.samples, 9);
    EXPECT_EQ(config.film.preset, "SpaceOdyssey2001");
    EXPECT_EQ(config.backend.preferred, "cpu");
}

TEST(ConfigSchema, PartialJsonKeepsDefaultsForOmittedFields) {
    // WITH_DEFAULT semantics: a document that names only some fields leaves the
    // rest at their struct defaults.
    const char* partial = R"({ "metric": { "name": "Alcubierre", "warpVelocity": 0.6 } })";

    SiriusConfig config = nlohmann::json::parse(partial).get<SiriusConfig>();

    EXPECT_EQ(config.metric.name, "Alcubierre");
    EXPECT_DOUBLE_EQ(config.metric.warp_velocity, 0.6);
    // Omitted fields keep defaults.
    EXPECT_EQ(config.render.width, 1920);
    EXPECT_EQ(config.render.output_path, "render.ppm");
    EXPECT_DOUBLE_EQ(config.metric.charge, 0.0);
}

}  // namespace sirius::app::test
