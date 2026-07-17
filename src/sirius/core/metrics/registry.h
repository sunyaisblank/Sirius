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

#include <array>
#include <cctype>
#include <optional>
#include <string>

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

// The registry. Backend flags state what each path genuinely renders today;
// requesting an unsupported combination is an error at the session boundary,
// never a substituted spacetime.
inline const std::array<MetricInfo, 9>& MetricRegistry() {
    static const std::array<MetricInfo, 9> registry = {{
        {MetricId::Minkowski, "Minkowski",
         {nullptr, nullptr, nullptr, nullptr, nullptr},
         "none (flat spacetime)", true, true},
        {MetricId::Schwarzschild, "Schwarzschild",
         {nullptr, nullptr, nullptr, nullptr, nullptr},
         "mass", true, true},
        {MetricId::Kerr, "Kerr",
         {nullptr, nullptr, nullptr, nullptr, nullptr},
         "mass, spin", true, true},
        {MetricId::ReissnerNordstrom, "Reissner-Nordstrom",
         {"ReissnerNordstrom", "Reissner-Nordström", nullptr, nullptr, nullptr},
         "mass, charge", true, false},
        {MetricId::KerrNewman, "Kerr-Newman",
         {"KerrNewman", nullptr, nullptr, nullptr, nullptr},
         "mass, spin, charge", true, false},
        {MetricId::DeSitter, "de-Sitter",
         {"DeSitter", "de Sitter", nullptr, nullptr, nullptr},
         "lambda", true, false},
        {MetricId::SchwarzschildDeSitter, "Schwarzschild-de-Sitter",
         {"SchwarzschildDeSitter", "Schwarzschild-de Sitter", "SdS", nullptr, nullptr},
         "mass, lambda", true, false},
        // CPU false: the Morris-Thorne family evaluates in spherical coordinates,
        // and the CPU tracer drives Cartesian positions; wiring them together
        // would evaluate the metric with x read as r. A Cartesian embedding is
        // the recorded follow-up; until then the wormhole is GPU-only.
        {MetricId::MorrisThorne, "Morris-Thorne",
         {"MorrisThorne", "Wormhole", "Morris-Thorne Wormhole", "Ellis Drainhole",
          "Zero-Tidal Wormhole"},
         "throat radius", false, true},
        {MetricId::Alcubierre, "Alcubierre",
         {"WarpDrive", "Alcubierre Warp Drive", nullptr, nullptr, nullptr},
         "warp velocity, bubble radius, wall thickness", true, true},
    }};
    return registry;
}

// The registry row for an id.
inline const MetricInfo& MetricInfoFor(MetricId id) {
    for (const auto& info : MetricRegistry()) {
        if (info.id == id) return info;
    }
    // The registry covers every enumerator; reaching here is a logic error.
    return MetricRegistry()[0];
}

namespace detail {
// Case-insensitive ASCII comparison; multi-byte UTF-8 (the ö in the
// Reissner-Nordström display name) passes through byte-exact.
inline bool NamesEqual(const char* a, const std::string& b) {
    std::size_t i = 0;
    for (; a[i] != '\0' && i < b.size(); ++i) {
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
