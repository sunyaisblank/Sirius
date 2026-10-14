#pragma once

// The single error channel for fallible operations per docs/STYLE.md
// section 4: origin, operation, offending input, no implementation detail
// across boundaries. Invariant violations are contracts (contracts.h), not
// errors; declined work on live paths is an Error, never a fabricated value.

#include <expected>
#include <format>
#include <string>
#include <utility>

namespace sirius::base {

// Originating subsystem, for message prefixes and routing; callers branch on
// their own outcomes, not on another layer's domain.
enum class ErrorDomain {
    kConfiguration,
    kIo,
    kDevice,
    kKernel,
    kPhysics,
    kInternal,
};

[[nodiscard]] constexpr const char* ToString(ErrorDomain domain) noexcept {
    switch (domain) {
        case ErrorDomain::kConfiguration:
            return "configuration";
        case ErrorDomain::kIo:
            return "io";
        case ErrorDomain::kDevice:
            return "device";
        case ErrorDomain::kKernel:
            return "kernel";
        case ErrorDomain::kPhysics:
            return "physics";
        case ErrorDomain::kInternal:
            return "internal";
    }
    return "unknown";
}

class [[nodiscard]] Error {
  public:
    Error(ErrorDomain domain, std::string operation, std::string detail)
        : domain_(domain), operation_(std::move(operation)), detail_(std::move(detail)) {}

    [[nodiscard]] ErrorDomain domain() const noexcept { return domain_; }
    [[nodiscard]] const std::string& operation() const noexcept { return operation_; }
    [[nodiscard]] const std::string& detail() const noexcept { return detail_; }

    [[nodiscard]] std::string Description() const {
        return std::format("{}: {}: {}", ToString(domain_), operation_, detail_);
    }

  private:
    ErrorDomain domain_;
    std::string operation_;
    std::string detail_;
};

template <typename T>
using Expected = std::expected<T, Error>;

// The one-line failure constructor: Fail(kIo, "open scene file", path).
[[nodiscard]] inline std::unexpected<Error> Fail(ErrorDomain domain, std::string operation,
                                                 std::string detail) {
    return std::unexpected(Error{domain, std::move(operation), std::move(detail)});
}

}  // namespace sirius::base
