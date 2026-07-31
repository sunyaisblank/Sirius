#pragma once

// Single authority for metric identity: every subsystem that names, parses,
// validates, lists, or dispatches on a metric goes through this table. Parsing
// is case-insensitive and failure is explicit; there is no default metric and no
// silent fallthrough. Ported from PHMT200A.h.
//
// Before this registry existed, metric identity was a free-form string compared
// by exact literal at six independent sites whose alias sets had drifted apart;
// the Morris-Thorne routing defect (a validator-accepted name falling through a
// GPU router to the wrong spacetime) was the direct result.

#include "sirius/base/contracts.h"

#include <array>
#include <cctype>
#include <optional>
#include <string>
#include <string_view>

namespace sirius::core {

// Closed set of spacetimes the engine can represent.
enum class MetricId {
    Minkowski,
    Schwarzschild,
    Kerr,
    ReissnerNordstrom,
    KerrNewman,
    DeSitter,
    SchwarzschildDeSitter,
    MorrisThorne,
    Alcubierre,
};

// Registry row: identity, spellings, parameters, backend support.
struct MetricInfo {
    MetricId id;
    const char* canonical_name;          // The one name Sirius prints.
    std::array<const char*, 5> aliases;  // Accepted spellings (nullptr-padded),
                                         // including GetName() outputs so display
                                         // names round-trip.
    const char* parameters;              // Human-readable parameter summary.
    bool cpu_supported;                  // Constructible on the CPU tracer path.
    bool gpu_supported;                  // Dispatchable to the device kernel.
};

enum class DiskSupport {
    PageThorne,
    NotApplicable,
    Unsupported,
};

[[nodiscard]] constexpr DiskSupport DiskSupportFor(MetricId id) noexcept {
    switch (id) {
        case MetricId::Schwarzschild:
        case MetricId::Kerr:
            return DiskSupport::PageThorne;
        case MetricId::ReissnerNordstrom:
        case MetricId::KerrNewman:
        case MetricId::SchwarzschildDeSitter:
            return DiskSupport::Unsupported;
        case MetricId::Minkowski:
        case MetricId::DeSitter:
        case MetricId::MorrisThorne:
        case MetricId::Alcubierre:
            return DiskSupport::NotApplicable;
    }
    SIRIUS_ASSERT(false);
    return DiskSupport::Unsupported;
}

[[nodiscard]] constexpr const char* ToString(DiskSupport support) noexcept {
    switch (support) {
        case DiskSupport::PageThorne:
            return "Page-Thorne";
        case DiskSupport::NotApplicable:
            return "not applicable";
        case DiskSupport::Unsupported:
            return "disabled required";
    }
    SIRIUS_ASSERT(false);
    return "unknown";
}

// Identity/parameter compatibility shared by config, CPU session, and Vulkan
// capability boundaries. Numeric range checks remain the owning boundary's
// responsibility; this authority prevents one metric name from acquiring
// another metric's spin, charge, or cosmological sector.
[[nodiscard]] constexpr std::optional<std::string_view> MetricParameterIssue(
    MetricId id, double spin, double charge, double lambda) noexcept {
    const bool has_spin = spin != 0.0;
    const bool has_charge = charge != 0.0;
    const bool has_lambda = lambda != 0.0;
    switch (id) {
        case MetricId::Schwarzschild:
            if (has_spin || has_charge || has_lambda) {
                return "Schwarzschild requires spin, charge, and lambda to be zero";
            }
            break;
        case MetricId::Kerr:
            if (has_charge || has_lambda) return "Kerr requires charge and lambda to be zero";
            break;
        case MetricId::ReissnerNordstrom:
            if (has_spin || has_lambda) {
                return "Reissner-Nordstrom requires spin and lambda to be zero";
            }
            break;
        case MetricId::KerrNewman:
            if (has_lambda) return "Kerr-Newman requires lambda to be zero";
            break;
        case MetricId::SchwarzschildDeSitter:
            if (has_spin || has_charge) {
                return "Schwarzschild-de-Sitter requires spin and charge to be zero";
            }
            break;
        case MetricId::DeSitter:
            if (has_spin || has_charge) return "de-Sitter requires spin and charge to be zero";
            break;
        case MetricId::Minkowski:
            if (has_spin || has_charge || has_lambda) {
                return "Minkowski requires spin, charge, and lambda to be zero";
            }
            break;
        case MetricId::MorrisThorne:
        case MetricId::Alcubierre:
            if (has_spin || has_charge || has_lambda) {
                return "Morris-Thorne and Alcubierre require spin, charge, and lambda to be zero";
            }
            break;
        default:
            return "unknown metric identity";
    }
    return std::nullopt;
}

// The registry. Backend flags state what each path genuinely renders today;
// requesting an unsupported combination is an error at the session boundary,
// never a substituted spacetime.
inline const std::array<MetricInfo, 9>& MetricRegistry() {
    static const std::array<MetricInfo, 9> registry = {{
        {MetricId::Minkowski,
         "Minkowski",
         {nullptr, nullptr, nullptr, nullptr, nullptr},
         "none (flat spacetime)",
         true,
         true},
        {MetricId::Schwarzschild,
         "Schwarzschild",
         {nullptr, nullptr, nullptr, nullptr, nullptr},
         "mass",
         true,
         true},
        {MetricId::Kerr,
         "Kerr",
         {nullptr, nullptr, nullptr, nullptr, nullptr},
         "mass, spin",
         true,
         true},
        {MetricId::ReissnerNordstrom,
         "Reissner-Nordstrom",
         {"ReissnerNordstrom", "Reissner-Nordström", nullptr, nullptr, nullptr},
         "mass, charge",
         true,
         false},
        {MetricId::KerrNewman,
         "Kerr-Newman",
         {"KerrNewman", nullptr, nullptr, nullptr, nullptr},
         "mass, spin, charge",
         true,
         false},
        {MetricId::DeSitter,
         "de-Sitter",
         {"DeSitter", "de Sitter", nullptr, nullptr, nullptr},
         "lambda",
         true,
         false},
        {MetricId::SchwarzschildDeSitter,
         "Schwarzschild-de-Sitter",
         {"SchwarzschildDeSitter", "Schwarzschild-de Sitter", "SdS", nullptr, nullptr},
         "mass, lambda",
         true,
         false},
        // CPU support arrives through MorrisThorneCartesian, the flat-plus-
        // rank-one embedding of the spherical family (one sheet, throat as the
        // capture surface); the spherical family remains the kernel-parity
        // reference and is never driven by the Cartesian tracer directly.
        {MetricId::MorrisThorne,
         "Morris-Thorne",
         {"MorrisThorne", "Wormhole", "Morris-Thorne Wormhole", "Ellis Drainhole",
          "Zero-Tidal Wormhole"},
         "throat radius",
         true,
         true},
        {MetricId::Alcubierre,
         "Alcubierre",
         {"WarpDrive", "Alcubierre Warp Drive", nullptr, nullptr, nullptr},
         "warp velocity, bubble radius, wall thickness",
         true,
         true},
    }};
    return registry;
}

// The registry row for an id.
inline const MetricInfo& MetricInfoFor(MetricId id) {
    for (const auto& info : MetricRegistry()) {
        if (info.id == id) return info;
    }
    // The registry covers every enumerator; a malformed enum must not acquire
    // Minkowski semantics.
    SIRIUS_ASSERT(false);
    return MetricRegistry()[0];
}

namespace detail {
// Case-insensitive ASCII comparison; multi-byte UTF-8 (the ö in the
// Reissner-Nordström display name) passes through byte-exact.
inline bool NamesEqual(const char* a, const std::string& b) {
    std::size_t i = 0;
    for (; i < b.size() && a[i] != '\0'; ++i) {
        unsigned char ca = static_cast<unsigned char>(a[i]);
        unsigned char cb = static_cast<unsigned char>(b[i]);
        if (std::tolower(ca) != std::tolower(cb)) return false;
    }
    return a[i] == '\0' && i == b.size();
}
}  // namespace detail

// Parse a metric name (canonical or alias, case-insensitive). Returns nullopt
// for an unknown name; callers must treat nullopt as an error and report the
// accepted names, never substitute a default.
inline std::optional<MetricId> ParseMetricName(const std::string& name) {
    for (const auto& info : MetricRegistry()) {
        if (detail::NamesEqual(info.canonical_name, name)) return info.id;
        for (const char* alias : info.aliases) {
            if (alias != nullptr && detail::NamesEqual(alias, name)) return info.id;
        }
    }
    return std::nullopt;
}

// Comma-separated canonical names, for error messages and help text.
inline std::string KnownMetricNames() {
    std::string names;
    for (const auto& info : MetricRegistry()) {
        if (!names.empty()) names += ", ";
        names += info.canonical_name;
    }
    return names;
}

}  // namespace sirius::core
