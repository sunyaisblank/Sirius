// Platform path tests: the resolver reports a platform name, a non-empty config
// search list, and declines an absent resource without throwing.

#include <gtest/gtest.h>

#include "sirius/app/platform_paths.h"

namespace sirius::app::test {

TEST(PlatformPaths, PlatformNameIsNonEmpty) {
    EXPECT_FALSE(PlatformPaths::PlatformName().empty());
}

TEST(PlatformPaths, ConfigSearchPathsAreOrderedAndNonEmpty) {
    auto paths = PlatformPaths::ConfigSearchPaths();
    ASSERT_FALSE(paths.empty());
    // The current directory is the highest-priority location.
    EXPECT_EQ(paths.front().filename().string(), "sirius.json");
}

TEST(PlatformPaths, AbsentResourceResolvesToNullopt) {
    auto resolved = PlatformPaths::ResolveResource("definitely-not-a-real-resource-xyz.frag");
    EXPECT_FALSE(resolved.has_value());
}

TEST(PlatformPaths, ExecutableDirectoryIsResolved) {
    // On Linux the /proc/self/exe readlink yields a concrete directory.
    EXPECT_FALSE(PlatformPaths::ExecutableDirectory().empty());
}

}  // namespace sirius::app::test
