// Configuration JSON tests: the non-intrusive serialisation re-attached in
// config_schema.h round-trips a SiriusConfig and parses the exact legacy field
// spellings, so existing config files load identically. This is the schema
// compatibility check for the app/config port.

#include <gtest/gtest.h>

#include <nlohmann/json.hpp>

#include "sirius/app/config/config_schema.h"

namespace sirius::app::test {

TEST(ConfigSchema, DefaultsRoundTripThroughJson) {
    SiriusConfig original = SiriusConfig::defaults();
    original.metric.name = "Kerr";
    original.metric.spin = 0.9;
    original.render.samplesPerPixel = 128;
    original.postprocess.tonemapper = "Filmic";

    nlohmann::json j = original;
    SiriusConfig restored = j.get<SiriusConfig>();

    EXPECT_EQ(restored.metric.name, "Kerr");
    EXPECT_DOUBLE_EQ(restored.metric.spin, 0.9);
    EXPECT_EQ(restored.render.samplesPerPixel, 128);
    EXPECT_EQ(restored.postprocess.tonemapper, "Filmic");
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
                    "warpVelocity": 0.8, "bubbleRadius": 2.0, "bubbleSigma": 0.25 },
        "observer": { "distance": 100.0, "inclination": 60.0, "azimuth": 30.0, "fov": 45.0 },
        "postprocess": { "enableBloom": false, "bloomIntensity": 0.5, "bloomThreshold": 0.2,
                         "exposure": 2.0, "contrast": 1.2, "saturation": 0.9,
                         "tonemapper": "Reinhard" },
        "volumetric": { "enabled": true, "hOverR": 0.15, "hPower": 0.3, "tauMidplane": 12.0,
                        "samples": 96, "enableTurbulence": true, "enableCorona": false },
        "film": { "enabled": true, "preset": "SpaceOdyssey2001", "grainIntensity": 0.2,
                  "halationStrength": 0.1, "vignetteStrength": 0.4 },
        "backend": { "preferred": "cpu", "enableDenoiser": false, "cudaDevice": 0 }
    })";

    SiriusConfig config = nlohmann::json::parse(legacy).get<SiriusConfig>();

    EXPECT_EQ(config.render.width, 800);
    EXPECT_EQ(config.render.outputPath, "out.exr");
    EXPECT_EQ(config.metric.name, "Kerr");
    EXPECT_DOUBLE_EQ(config.metric.mass, 2.0);
    EXPECT_EQ(config.metric.temperatureModel, "ShakuraSunyaev");
    EXPECT_FLOAT_EQ(config.metric.diskTemperature, 60000.0f);
    EXPECT_DOUBLE_EQ(config.observer.azimuth, 30.0);
    EXPECT_FALSE(config.postprocess.enableBloom);
    EXPECT_EQ(config.postprocess.tonemapper, "Reinhard");
    EXPECT_TRUE(config.volumetric.enableTurbulence);
    EXPECT_EQ(config.film.preset, "SpaceOdyssey2001");
    EXPECT_EQ(config.backend.preferred, "cpu");
}

TEST(ConfigSchema, PartialJsonKeepsDefaultsForOmittedFields) {
    // WITH_DEFAULT semantics: a document that names only some fields leaves the
    // rest at their struct defaults.
    const char* partial = R"({ "metric": { "name": "Alcubierre", "warpVelocity": 0.6 } })";

    SiriusConfig config = nlohmann::json::parse(partial).get<SiriusConfig>();

    EXPECT_EQ(config.metric.name, "Alcubierre");
    EXPECT_DOUBLE_EQ(config.metric.warpVelocity, 0.6);
    // Omitted fields keep defaults.
    EXPECT_EQ(config.render.width, 1920);
    EXPECT_EQ(config.render.outputPath, "render.ppm");
    EXPECT_DOUBLE_EQ(config.metric.charge, 0.0);
}

}  // namespace sirius::app::test
