#include "sirius/app/cli/command_router.h"

#include "sirius/app/cli/config_command.h"
#include "sirius/app/cli/info_command.h"
#include "sirius/app/config/config_schema.h"

#include <gtest/gtest.h>

#include <array>
#include <string>

namespace sirius::app::test {

namespace {

int InvokeRouter(CommandRouter& router, std::initializer_list<std::string> values) {
    std::vector<std::string> storage(values);
    std::vector<char*> argv;
    argv.reserve(storage.size());
    for (std::string& value : storage) {
        argv.push_back(value.data());
    }
    return router.Run(static_cast<int>(argv.size()), argv.data());
}

}  // namespace

TEST(CommandRouter, UnknownCommandReturnsFailure) {
    CommandRouter router;
    EXPECT_EQ(InvokeRouter(router, {"sirius", "definitely-not-a-command"}), 1);
}

TEST(CommandRouter, NoCommandStillPrintsHelpSuccessfully) {
    CommandRouter router;
    EXPECT_EQ(InvokeRouter(router, {"sirius"}), 0);
}

TEST(CommandRouter, ExplicitMissingConfigPropagatesFailure) {
    CommandRouter router;
    EXPECT_THROW(InvokeRouter(router, {"sirius", "--config", "/definitely/absent/sirius.json",
                                       "info", "system"}),
                 std::runtime_error);
}

TEST(CommandRouter, InfoRejectsEmptyAndSurplusArguments) {
    InfoCommand command;
    GlobalOptions globals;
    SiriusConfig config = SiriusConfig::defaults();
    EXPECT_EQ(command.Execute({"system", "ignored"}, globals, config), 1);
    EXPECT_EQ(command.Execute({""}, globals, config), 1);
    EXPECT_EQ(command.Execute({"--unknown"}, globals, config), 1);
}

TEST(CommandRouter, ConfigRejectsMalformedSubcommandArguments) {
    ConfigCommand command;
    GlobalOptions globals;
    SiriusConfig config = SiriusConfig::defaults();
    EXPECT_EQ(command.Execute({"show", "ignored"}, globals, config), 1);
    EXPECT_EQ(command.Execute({"paths", "ignored"}, globals, config), 1);
    EXPECT_EQ(command.Execute({"validate", "--unknown"}, globals, config), 1);
    EXPECT_EQ(command.Execute({"init", "--output"}, globals, config), 1);
    EXPECT_EQ(command.Execute({"init", "--unknown", "path"}, globals, config), 1);
}

}  // namespace sirius::app::test
