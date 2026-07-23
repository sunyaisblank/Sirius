#pragma once

// Cross-platform, process-local environment mutation for tests. The original
// value is restored at scope exit so individually discovered tests remain
// order-independent even when the invoking shell supplied the variable.

#include <cstdlib>
#include <optional>
#include <stdexcept>
#include <string>

namespace sirius::test {

inline int MutateEnvironment(const char* name, const char* value) {
#ifdef _WIN32
    return ::_putenv_s(name, value != nullptr ? value : "");
#else
    return value != nullptr ? ::setenv(name, value, 1) : ::unsetenv(name);
#endif
}

class ScopedEnvironmentVariable {
  public:
    ScopedEnvironmentVariable(const char* name, const char* value) : name_(name) {
        if (const char* original = std::getenv(name); original != nullptr) {
            original_ = original;
        }
        if (MutateEnvironment(name_.c_str(), value) != 0) {
            throw std::runtime_error("failed to set test environment variable: " + name_);
        }
    }

    ~ScopedEnvironmentVariable() {
        const char* value = original_.has_value() ? original_->c_str() : nullptr;
        (void)MutateEnvironment(name_.c_str(), value);
    }

    ScopedEnvironmentVariable(const ScopedEnvironmentVariable&) = delete;
    ScopedEnvironmentVariable& operator=(const ScopedEnvironmentVariable&) = delete;
    ScopedEnvironmentVariable(ScopedEnvironmentVariable&&) = delete;
    ScopedEnvironmentVariable& operator=(ScopedEnvironmentVariable&&) = delete;

  private:
    std::string name_;
    std::optional<std::string> original_;
};

}  // namespace sirius::test
