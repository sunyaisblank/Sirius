// PHMT200A.h - Metric Identity Registry
// Component ID: PHMT200A (Physics/Metric/Identity Registry)
//
// Single authority for metric identification. Every subsystem that needs to
// name, parse, validate, list, or dispatch on a metric goes through this
// table. Before this registry existed, metric identity was a free-form string
// compared by exact literal at six independent sites, whose alias sets had
// drifted apart; the Morris-Thorne routing defect (a validator-accepted name
// falling through a GPU router to render the wrong spacetime) was the direct
// result. Parsing is case-insensitive and failure is explicit: there is no
// default metric and no silent fallthrough anywhere downstream.

#ifndef PHMT200A_H
#define PHMT200A_H

#include <array>
#include <cctype>
#include <optional>
#include <string>

namespace Sirius {

/// @brief Closed set of spacetimes the engine can represent
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

/// @brief Registry row: identity, spellings, parameters, backend support
struct MetricInfo {
    MetricId id;
    const char* canonicalName;              ///< The one name Sirius prints
    std::array<const char*, 5> aliases;     ///< Accepted spellings (nullptr-padded),
                                            ///< including IMetric::getName() outputs
                                            ///< so display names round-trip
    const char* parameters;                 ///< Human-readable parameter summary
    bool cpuSupported;                      ///< Constructible on the CPU tracer path
    bool gpuSupported;                      ///< Dispatchable to the OptiX kernel
};

/// @brief The registry. Backend flags state what each path genuinely renders
/// today; requesting an unsupported combination is an error at the session
/// boundary, never a substituted spacetime.
inline const std::array<MetricInfo, 9>& metricRegistry() {
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
        // CPU false: MorrisThorneFamily evaluates in spherical coordinates,
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

/// @brief Look up the registry row for an id
inline const MetricInfo& metricInfo(MetricId id) {
    for (const auto& info : metricRegistry()) {
        if (info.id == id) return info;
    }
    // The registry covers every enumerator; reaching here is a logic error.
    return metricRegistry()[0];
}

namespace Detail {
/// Case-insensitive ASCII comparison; multi-byte UTF-8 (the ö in the
/// Reissner-Nordström display name) passes through byte-exact.
inline bool namesEqual(const char* a, const std::string& b) {
    std::size_t i = 0;
    for (; a[i] != '\0' && i < b.size(); ++i) {
        unsigned char ca = static_cast<unsigned char>(a[i]);
        unsigned char cb = static_cast<unsigned char>(b[i]);
        if (std::tolower(ca) != std::tolower(cb)) return false;
    }
    return a[i] == '\0' && i == b.size();
}
} // namespace Detail

/// @brief Parse a metric name (canonical or alias, case-insensitive)
/// @return The identity, or std::nullopt for an unknown name. Callers must
///         treat nullopt as an error and report the accepted names; they must
///         not substitute a default.
inline std::optional<MetricId> parseMetricName(const std::string& name) {
    for (const auto& info : metricRegistry()) {
        if (Detail::namesEqual(info.canonicalName, name)) return info.id;
        for (const char* alias : info.aliases) {
            if (alias != nullptr && Detail::namesEqual(alias, name)) return info.id;
        }
    }
    return std::nullopt;
}

/// @brief Comma-separated canonical names, for error messages and help text
inline std::string knownMetricNames() {
    std::string names;
    for (const auto& info : metricRegistry()) {
        if (!names.empty()) names += ", ";
        names += info.canonicalName;
    }
    return names;
}

} // namespace Sirius

#endif // PHMT200A_H
