#include "sirius/base/resource_locator.h"

#include <cstdint>
#include <cstdlib>
#include <system_error>
#include <vector>

#if defined(_WIN32)
#include <windows.h>
#elif defined(__linux__)
#include <unistd.h>
#elif defined(__APPLE__)
#include <mach-o/dyld.h>
#endif

namespace sirius::base {

namespace {

[[nodiscard]] bool IsSafeRelativePath(const fs::path& path) {
    if (path.empty() || path.is_absolute()) return false;
    for (const auto& component : path) {
        if (component == "." || component == "..") return false;
    }
    return true;
}

[[nodiscard]] bool IsWithin(const fs::path& path, const fs::path& root) {
    auto path_it = path.begin();
    for (auto root_it = root.begin(); root_it != root.end(); ++root_it, ++path_it) {
        if (path_it == path.end() || *path_it != *root_it) return false;
    }
    return true;
}

[[nodiscard]] fs::path DiscoverExecutableDirectory() {
#if defined(_WIN32)
    std::vector<wchar_t> path(512);
    while (path.size() <= 32768) {
        const DWORD length =
            GetModuleFileNameW(nullptr, path.data(), static_cast<DWORD>(path.size()));
        if (length == 0) return {};
        if (length < path.size()) {
            return fs::path(std::wstring_view(path.data(), length)).parent_path();
        }
        path.resize(path.size() * 2);
    }
#elif defined(__linux__)
    std::vector<char> path(1024);
    while (path.size() <= 1048576) {
        const ssize_t length = readlink("/proc/self/exe", path.data(), path.size());
        if (length < 0) return {};
        if (static_cast<std::size_t>(length) < path.size()) {
            return fs::path(std::string_view(path.data(), static_cast<std::size_t>(length)))
                .parent_path();
        }
        path.resize(path.size() * 2);
    }
#elif defined(__APPLE__)
    std::uint32_t size = 0;
    (void)_NSGetExecutablePath(nullptr, &size);
    std::vector<char> path(size);
    if (size > 0 && _NSGetExecutablePath(path.data(), &size) == 0) {
        std::error_code ec;
        const fs::path canonical = fs::canonical(fs::path(path.data()), ec);
        if (!ec) {
            return canonical.parent_path();
        }
        return fs::path(path.data()).parent_path();
    }
#endif
    return {};
}

}  // namespace

fs::path ExecutableDirectory() {
    static const fs::path directory = DiscoverExecutableDirectory();
    return directory;
}

std::vector<fs::path> ResourceCandidates(std::string_view relative_path) {
    const fs::path relative(relative_path);
    if (!IsSafeRelativePath(relative)) return {};
    if (const char* explicit_root = std::getenv("SIRIUS_RESOURCE_DIR");
        explicit_root != nullptr && *explicit_root != '\0') {
        return {fs::path(explicit_root) / relative};
    }

    const fs::path executable = ExecutableDirectory();
    std::vector<fs::path> candidates;
    candidates.reserve(2);
    if (!executable.empty()) {
        candidates.push_back(executable / "resources" / relative);
        candidates.push_back(executable.parent_path() / "share" / "sirius" / relative);
    }
    return candidates;
}

std::optional<fs::path> ResolveResource(std::string_view relative_path) {
    const fs::path relative(relative_path);
    if (!IsSafeRelativePath(relative)) return std::nullopt;

    for (const fs::path& candidate : ResourceCandidates(relative_path)) {
        std::error_code ec;
        if (!fs::is_regular_file(candidate, ec) || ec) {
            continue;
        }
        const fs::path canonical = fs::canonical(candidate, ec);
        if (ec) continue;

        fs::path root = candidate;
        for (auto it = relative.begin(); it != relative.end(); ++it) {
            root = root.parent_path();
        }
        const fs::path canonical_root = fs::canonical(root, ec);
        if (!ec && IsWithin(canonical, canonical_root)) return canonical;
    }
    return std::nullopt;
}

}  // namespace sirius::base
