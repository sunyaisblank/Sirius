#include "sirius/app/cli/config_command.h"
#include "sirius/app/config/config_loader.h"

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
    EXPECT_THROW(ConfigLoader::Load(path.string()), std::runtime_error);
}

TEST(ConfigLoading, MalformedJsonDeclinesWithoutDefaults) {
    TemporaryConfig config(R"({"render": {"width": 512,})");
    EXPECT_THROW(ConfigLoader::LoadFromFile(config.Path()), nlohmann::json::exception);
}

TEST(ConfigLoading, UnknownNestedFieldDeclines) {
    TemporaryConfig config(R"({"render": {"widht": 512}})");
    EXPECT_THROW(ConfigLoader::LoadFromFile(config.Path()), std::invalid_argument);
}

TEST(ConfigLoading, UnknownTopLevelFieldDeclines) {
    TemporaryConfig config(R"({"renderer": {"width": 512}})");
    EXPECT_THROW(ConfigLoader::LoadFromFile(config.Path()), std::invalid_argument);
}

TEST(ConfigLoading, DuplicateFieldsDeclineInsteadOfUsingParserOrder) {
    TemporaryConfig top_level(R"({"diskEnabled": true, "diskEnabled": false})");
    EXPECT_THROW(ConfigLoader::LoadFromFile(top_level.Path()), std::invalid_argument);

    TemporaryConfig nested(R"({"render": {"width": 512, "width": 1024}})");
    EXPECT_THROW(ConfigLoader::LoadFromFile(nested.Path()), std::invalid_argument);
}

TEST(ConfigLoading, ValidateCommandUsesTheSameStrictParserAsStartup) {
    ConfigCommand command;
    GlobalOptions globals;
    SiriusConfig config = SiriusConfig::defaults();

    TemporaryConfig unknown(R"({"render": {"widht": 512}})");
    EXPECT_EQ(command.Execute({"validate", unknown.Path().string()}, globals, config), 1);

    TemporaryConfig duplicate(R"({"render": {"width": 512, "width": 1024}})");
    EXPECT_EQ(command.Execute({"validate", duplicate.Path().string()}, globals, config), 1);

    TemporaryConfig valid(R"({"render": {"width": 512}})");
    EXPECT_EQ(command.Execute({"validate", valid.Path().string()}, globals, config), 0);
}

TEST(ConfigLoading, InvalidKnownValueDeclines) {
    TemporaryConfig config(R"({"render": {"width": 64}})");
    EXPECT_THROW(ConfigLoader::LoadFromFile(config.Path()), std::invalid_argument);
}

TEST(ConfigLoading, ValidPartialFileMergesOverDefaults) {
    TemporaryConfig config(R"({"render": {"width": 512}})");
    const SiriusConfig loaded = ConfigLoader::LoadFromFile(config.Path());
    EXPECT_EQ(loaded.render.width, 512);
    EXPECT_EQ(loaded.render.height, SiriusConfig::defaults().render.height);
}

TEST(ConfigLoading, SaveDeclinesWhenParentCannotBeCreated) {
    TemporaryConfig regular_file_parent("{}");
    const auto impossible_child = regular_file_parent.Path() / "child.json";
    EXPECT_FALSE(ConfigLoader::SaveToFile(SiriusConfig::defaults(), impossible_child));
    EXPECT_FALSE(std::filesystem::exists(impossible_child));
}

TEST(ConfigLoading, SaveDeclinesInvalidConfiguration) {
    TemporaryConfig output("{}");
    SiriusConfig invalid = SiriusConfig::defaults();
    invalid.render.width = 0;
    EXPECT_FALSE(ConfigLoader::SaveToFile(invalid, output.Path()));
}

}  // namespace sirius::app::test
