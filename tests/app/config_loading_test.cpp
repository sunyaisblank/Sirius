#include "sirius/app/cli/config_command.h"
#include "sirius/app/config/config_loader.h"
#include "sirius/app/config/session_config_adapter.h"

#include <gtest/gtest.h>

#include <chrono>
#include <filesystem>
#include <fstream>
#include <string>

namespace sirius::app::test {

namespace {

class TemporaryConfig {
  public:
    explicit TemporaryConfig(const std::string& contents) {
        const auto stamp = std::chrono::steady_clock::now().time_since_epoch().count();
        path_ = std::filesystem::temp_directory_path() /
                ("sirius-config-adversarial-" + std::to_string(stamp) + ".json");
        std::ofstream output(path_);
        output << contents;
        output.close();
    }

    ~TemporaryConfig() {
        std::error_code ec;
        std::filesystem::remove(path_, ec);
    }

    const std::filesystem::path& Path() const { return path_; }

  private:
    std::filesystem::path path_;
};

}  // namespace

TEST(ConfigLoading, ExplicitMissingFileDeclines) {
    const auto path =
        std::filesystem::temp_directory_path() / "sirius-config-definitely-absent.json";
    std::error_code ec;
    std::filesystem::remove(path, ec);
    EXPECT_FALSE(ConfigLoader::Load(path.string()).has_value());
}

TEST(ConfigLoading, MalformedJsonDeclinesWithoutDefaults) {
    TemporaryConfig config(R"({"render": {"width": 512,})");
    EXPECT_FALSE(ConfigLoader::LoadFromFile(config.Path()).has_value());
}

TEST(ConfigLoading, UnknownNestedFieldDeclines) {
    TemporaryConfig config(R"({"render": {"widht": 512}})");
    EXPECT_FALSE(ConfigLoader::LoadFromFile(config.Path()).has_value());
}

TEST(ConfigLoading, UnknownTopLevelFieldDeclines) {
    TemporaryConfig config(R"({"renderer": {"width": 512}})");
    EXPECT_FALSE(ConfigLoader::LoadFromFile(config.Path()).has_value());
}

TEST(ConfigLoading, DuplicateFieldsDeclineInsteadOfUsingParserOrder) {
    TemporaryConfig top_level(R"({"diskEnabled": true, "diskEnabled": false})");
    EXPECT_FALSE(ConfigLoader::LoadFromFile(top_level.Path()).has_value());

    TemporaryConfig nested(R"({"render": {"width": 512, "width": 1024}})");
    EXPECT_FALSE(ConfigLoader::LoadFromFile(nested.Path()).has_value());
}

TEST(ConfigLoading, ValidateCommandUsesTheSameStrictParserAsStartup) {
    ConfigCommand command;
    GlobalOptions globals;
    SiriusConfig config = SiriusConfig::Defaults();

    TemporaryConfig unknown(R"({"render": {"widht": 512}})");
    EXPECT_EQ(command.Execute({"validate", unknown.Path().string()}, globals, config), 1);

    TemporaryConfig duplicate(R"({"render": {"width": 512, "width": 1024}})");
    EXPECT_EQ(command.Execute({"validate", duplicate.Path().string()}, globals, config), 1);

    TemporaryConfig valid(R"({"render": {"width": 512}})");
    EXPECT_EQ(command.Execute({"validate", valid.Path().string()}, globals, config), 0);
}

TEST(ConfigLoading, InvalidKnownValueDeclines) {
    TemporaryConfig config(R"({"render": {"width": 64}})");
    EXPECT_FALSE(ConfigLoader::LoadFromFile(config.Path()).has_value());
}

TEST(ConfigLoading, ValidPartialFileMergesOverDefaults) {
    {
        TemporaryConfig config(R"({"render": {"width": 512}})");
        const auto loaded = ConfigLoader::LoadFromFile(config.Path());
        ASSERT_TRUE(loaded.has_value()) << loaded.error().Description();
        EXPECT_EQ(loaded->render.width, 512);
        EXPECT_EQ(loaded->render.height, SiriusConfig::Defaults().render.height);
    }
    {
        TemporaryConfig config(R"({"metric": {"name": "Alcubierre"}, "diskEnabled": false})");
        const auto loaded = ConfigLoader::LoadFromFile(config.Path());
        ASSERT_TRUE(loaded.has_value()) << loaded.error().Description();
        EXPECT_DOUBLE_EQ(loaded->metric.mass, 0.0);
    }
    {
        TemporaryConfig config(
            R"({"metric": {"name": "Alcubierre", "mass": 1.0}, "diskEnabled": false})");
        EXPECT_FALSE(ConfigLoader::LoadFromFile(config.Path()).has_value());
    }
}

TEST(ConfigLoading, PointCatalogueSurvivesLoadingSavingAndSessionProjection) {
    TemporaryConfig input(R"({
        "pointStarfield": true,
        "pointStarfieldConfig": {
            "starCount": 1e5, "minDistancePc": 2.5, "maxDistancePc": 9000,
            "brightnessScale": 1e-5, "seed": 4294967295
        }
    })");
    const auto loaded = ConfigLoader::LoadFromFile(input.Path());
    ASSERT_TRUE(loaded) << loaded.error().Description();
    EXPECT_EQ(loaded->point_starfield_config.star_count, 100000u);
    EXPECT_EQ(loaded->point_starfield_config.seed, 4294967295u);
    EXPECT_FLOAT_EQ(loaded->point_starfield_config.min_distance_pc, 2.5f);
    EXPECT_FLOAT_EQ(loaded->point_starfield_config.max_distance_pc, 9000.0f);
    EXPECT_FLOAT_EQ(loaded->point_starfield_config.brightness_scale, 1e-5f);
    const auto session = MakeSessionConfig(*loaded);
    ASSERT_TRUE(session) << session.error().Description();
    EXPECT_EQ(session->point_starfield_config, loaded->point_starfield_config);
    EXPECT_TRUE(session->point_starfield);
    TemporaryConfig saved("{}");
    ASSERT_TRUE(ConfigLoader::SaveToFile(*loaded, saved.Path()));
    const auto reloaded = ConfigLoader::LoadFromFile(saved.Path());
    ASSERT_TRUE(reloaded);
    EXPECT_EQ(reloaded->point_starfield_config, session->point_starfield_config);
}

TEST(ConfigLoading, PointCatalogueRejectsChangedIntegerIdentitiesAndUnknownFields) {
    for (const auto* field : {"starCount", "seed"}) {
        for (const auto* value :
             {"-1", "-4294967295", "4294967296", "4294967297", "1.5", "1e30", "true", "\"42\""}) {
            SCOPED_TRACE(std::string(field) + "=" + value);
            TemporaryConfig config(
                std::string("{\"pointStarfield\":true,\"pointStarfieldConfig\":{\"") + field +
                "\":" + value + "}}");
            const auto loaded = ConfigLoader::LoadFromFile(config.Path());
            ASSERT_FALSE(loaded);
            EXPECT_NE(loaded.error().Description().find(field), std::string::npos);
        }
    }
    TemporaryConfig unknown(R"({"pointStarfield":true,"pointStarfieldConfig":{"starCout":100}})");
    EXPECT_FALSE(ConfigLoader::LoadFromFile(unknown.Path()));
    TemporaryConfig duplicate(
        R"({"pointStarfield":true,"pointStarfieldConfig":{"seed":1,"seed":2}})");
    EXPECT_FALSE(ConfigLoader::LoadFromFile(duplicate.Path()));
    TemporaryConfig malformed(R"({"pointStarfield":true,"pointStarfieldConfig":[]})");
    EXPECT_FALSE(ConfigLoader::LoadFromFile(malformed.Path()));
}

TEST(ConfigLoading, SaveDeclinesWhenParentCannotBeCreated) {
    TemporaryConfig regular_file_parent("{}");
    const auto impossible_child = regular_file_parent.Path() / "child.json";
    EXPECT_FALSE(ConfigLoader::SaveToFile(SiriusConfig::Defaults(), impossible_child).has_value());
    EXPECT_FALSE(std::filesystem::exists(impossible_child));
}

TEST(ConfigLoading, SaveDeclinesInvalidConfiguration) {
    TemporaryConfig output("{}");
    SiriusConfig invalid = SiriusConfig::Defaults();
    invalid.render.width = 0;
    EXPECT_FALSE(ConfigLoader::SaveToFile(invalid, output.Path()).has_value());
}

}  // namespace sirius::app::test
