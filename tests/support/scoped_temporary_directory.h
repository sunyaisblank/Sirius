#pragma once

#include <chrono>
#include <cstdint>
#include <filesystem>
#include <functional>
#include <random>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>

// Owns one collision-resistant test directory. Fixed filenames directly under
// the system temporary directory are unsafe when compiler estates run in
// parallel: an otherwise independent process can remove or replace them.
namespace sirius::test {

class ScopedTemporaryDirectory {
  public:
    explicit ScopedTemporaryDirectory(std::string_view tag) {
        const auto now = std::chrono::steady_clock::now().time_since_epoch().count();
        const auto working_directory = std::filesystem::current_path().string();
        std::random_device random;
        const std::uint64_t entropy =
            static_cast<std::uint64_t>(now) ^
            static_cast<std::uint64_t>(std::hash<std::string>{}(working_directory)) ^
            (static_cast<std::uint64_t>(random()) << 32U) ^ random();

        for (std::uint64_t attempt = 0; attempt < 128; ++attempt) {
            path_ = std::filesystem::temp_directory_path() /
                    (std::string(tag) + "-" + std::to_string(entropy + attempt));
            std::error_code error;
            if (std::filesystem::create_directory(path_, error)) return;
            if (error && error != std::errc::file_exists) {
                throw std::runtime_error("could not create an isolated test directory: " +
                                         error.message());
            }
        }
        throw std::runtime_error("could not allocate a unique isolated test directory");
    }

    ~ScopedTemporaryDirectory() {
        std::error_code ignored;
        std::filesystem::remove_all(path_, ignored);
    }

    ScopedTemporaryDirectory(const ScopedTemporaryDirectory&) = delete;
    ScopedTemporaryDirectory& operator=(const ScopedTemporaryDirectory&) = delete;
    ScopedTemporaryDirectory(ScopedTemporaryDirectory&&) = delete;
    ScopedTemporaryDirectory& operator=(ScopedTemporaryDirectory&&) = delete;

    [[nodiscard]] const std::filesystem::path& path() const noexcept { return path_; }

  private:
    std::filesystem::path path_;
};

}  // namespace sirius::test
