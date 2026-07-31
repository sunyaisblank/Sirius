// Platform path tests: the resolver reports a platform name, a non-empty config
// search list, and declines an absent resource without throwing.

#include "sirius/app/platform_paths.h"

#include <gtest/gtest.h>

#include "support/scoped_environment.h"

#include <filesystem>
#include <fstream>

namespace sirius::app::test {

TEST(PlatformPaths, PlatformNameIsNonEmpty) { EXPECT_FALSE(PlatformPaths::PlatformName().empty()); }

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

TEST(PlatformPaths, ExplicitResourceRootDisablesFallbacks) {
    const auto empty_root =
        std::filesystem::temp_directory_path() / "sirius-strict-resource-root-does-not-exist";
    std::error_code ec;
    std::filesystem::remove_all(empty_root, ec);
    sirius::test::ScopedEnvironmentVariable strict_root("SIRIUS_RESOURCE_DIR",
                                                        empty_root.string().c_str());
    EXPECT_FALSE(PlatformPaths::ResolveResource("assets/Starfield.png").has_value())
        << "an explicit volume root must not borrow assets from the source tree";
}

TEST(PlatformPaths, ResourceNamesAndSymlinksCannotEscapeTheSelectedVolume) {
    const auto parent = std::filesystem::temp_directory_path() / "sirius-resource-containment-test";
    const auto root = parent / "volume";
    const auto outside = parent / "outside.txt";
    std::error_code ec;
    std::filesystem::remove_all(parent, ec);
    std::filesystem::create_directories(root, ec);
    ASSERT_FALSE(ec);
    {
        std::ofstream file(outside);
        ASSERT_TRUE(file.good());
        file << "not part of the Sirius volume";
    }

    sirius::test::ScopedEnvironmentVariable strict_root("SIRIUS_RESOURCE_DIR",
                                                        root.string().c_str());
    EXPECT_FALSE(PlatformPaths::ResolveResource("../outside.txt").has_value());
    EXPECT_FALSE(PlatformPaths::ResolveResource("./outside.txt").has_value());
    EXPECT_FALSE(PlatformPaths::ResolveResource(outside.string()).has_value());

#if !defined(_WIN32)
    std::filesystem::create_symlink(outside, root / "escaped.txt", ec);
    ASSERT_FALSE(ec);
    EXPECT_FALSE(PlatformPaths::ResolveResource("escaped.txt").has_value());
#endif

    std::filesystem::remove_all(parent, ec);
}

}  // namespace sirius::app::test
