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

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <limits>
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

inline constexpr double kDefaultMorrisThorneThroatRadius = 1.0;
inline constexpr double kMinMorrisThorneThroatRadius = 0.1;
inline constexpr double kMaxMorrisThorneThroatRadius = 1000.0;
inline constexpr double kDefaultAlcubierreWarpVelocity = 0.5;
inline constexpr double kDefaultAlcubierreBubbleRadius = 1.0;
inline constexpr double kDefaultAlcubierreBubbleSigma = 0.5;
inline constexpr double kMinAlcubierreSigmaRadius = 0.1;
inline constexpr double kMaxAlcubierreSigmaRadius = 100.0;
inline constexpr double kMaxCosmologicalObserverFraction = 0.99;

// Exact spherical positive-Lambda horizon authority.  These functions live
// below the metric class so configuration, tracing, and KerrSchildFamily cannot
// drift onto different roots of the Kottler lapse
//   f(r) = 1 - 2M/r - Lambda r^2/3.
// The scale-safe forms avoid overflowing 3/Lambda or Lambda*r^2 for very small
// positive Lambda values that are otherwise representable in double precision.
[[nodiscard]] inline double KottlerStaticLapse(double mass, double lambda, double radius) noexcept {
    if (!(mass >= 0.0) || !(lambda >= 0.0) || !(radius > 0.0) || !std::isfinite(mass) ||
        !std::isfinite(lambda) || !std::isfinite(radius)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    const double scaled_radius = std::sqrt(lambda) * radius / std::sqrt(3.0);
    return 1.0 - 2.0 * mass / radius - scaled_radius * scaled_radius;
}

[[nodiscard]] inline double KottlerBlackHoleHorizonRadius(double mass, double lambda) noexcept {
    if (!(mass > 0.0) || !(lambda > 0.0) || !std::isfinite(mass) || !std::isfinite(lambda) ||
        !(9.0 * lambda * mass * mass <= 1.0)) {
        return -1.0;
    }
    const double stationary_radius = std::cbrt(3.0 * mass) / std::cbrt(lambda);
    double lower = 2.0 * mass;
    double upper = stationary_radius;
    for (int iteration = 0; iteration < 128; ++iteration) {
        const double midpoint = lower + 0.5 * (upper - lower);
        if (KottlerStaticLapse(mass, lambda, midpoint) > 0.0) {
            upper = midpoint;
        } else {
            lower = midpoint;
        }
    }
    return lower + 0.5 * (upper - lower);
}

[[nodiscard]] inline double KottlerCosmologicalHorizonRadius(double mass, double lambda) noexcept {
    if (!(mass >= 0.0) || !(lambda > 0.0) || !std::isfinite(mass) || !std::isfinite(lambda) ||
        !(9.0 * lambda * mass * mass <= 1.0)) {
        return -1.0;
    }
    const double curvature_radius = std::sqrt(3.0) / std::sqrt(lambda);
    if (mass == 0.0) return curvature_radius;

    const double stationary_radius = std::cbrt(3.0 * mass) / std::cbrt(lambda);
    double lower = stationary_radius;
    double upper = curvature_radius;
    for (int iteration = 0; iteration < 128; ++iteration) {
        const double midpoint = lower + 0.5 * (upper - lower);
        if (KottlerStaticLapse(mass, lambda, midpoint) > 0.0) {
            lower = midpoint;
        } else {
            upper = midpoint;
        }
    }
    return lower + 0.5 * (upper - lower);
}

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

// Whether M is an actual parameter of this spacetime identity. This authority
// is shared by operator validation, typed-session validation, and numerical
// trace scaling so an unrelated value cannot silently alter a metric that has
// no mass parameter.
[[nodiscard]] constexpr bool MetricUsesMass(MetricId id) noexcept {
    switch (id) {
        case MetricId::Schwarzschild:
        case MetricId::Kerr:
        case MetricId::ReissnerNordstrom:
        case MetricId::KerrNewman:
        case MetricId::SchwarzschildDeSitter:
            return true;
        case MetricId::Minkowski:
        case MetricId::DeSitter:
        case MetricId::MorrisThorne:
        case MetricId::Alcubierre:
            return false;
    }
    SIRIUS_ASSERT(false);
    return false;
}

// Characteristic coordinate scales shared by observer validation and numerical
// tracing. For the tanh Alcubierre profile, sigma has inverse-length units: the
// represented scene extends over max(R, 1/sigma), while its narrowest feature is
// min(R, 1/sigma). Callers validate the selected parameters before relying on
// these values.
[[nodiscard]] constexpr double MetricSceneLengthScale(MetricId id, double mass,
                                                      double throat_radius, double bubble_radius,
                                                      double bubble_sigma) noexcept {
    if (MetricUsesMass(id)) return mass;
    switch (id) {
        case MetricId::MorrisThorne:
            return throat_radius;
        case MetricId::Alcubierre:
            return std::max(bubble_radius, 1.0 / bubble_sigma);
        case MetricId::Minkowski:
        case MetricId::DeSitter:
            return 1.0;
        case MetricId::Schwarzschild:
        case MetricId::Kerr:
        case MetricId::ReissnerNordstrom:
        case MetricId::KerrNewman:
        case MetricId::SchwarzschildDeSitter:
            break;
    }
    SIRIUS_ASSERT(false);
    return 1.0;
}

[[nodiscard]] constexpr double MetricFeatureLengthScale(MetricId id, double mass,
                                                        double throat_radius, double bubble_radius,
                                                        double bubble_sigma) noexcept {
    if (id == MetricId::Alcubierre) {
        return std::min(bubble_radius, 1.0 / bubble_sigma);
    }
    return MetricSceneLengthScale(id, mass, throat_radius, bubble_radius, bubble_sigma);
}

// The positive-Lambda renderer represents a causal-patch boundary-value
// problem, not an asymptotically flat sky.  Its observer must remain strictly
// inside the cosmological horizon and the trace boundary may not pass through
// that horizon.  A one-percent separation is the governed fp32 event margin;
// the camera itself remains the chart's timelike ADM/Eulerian observer rather
// than being reinterpreted as a static-coordinate observer.
[[nodiscard]] inline std::optional<double> MetricCosmologicalHorizonRadius(MetricId id, double mass,
                                                                           double lambda) noexcept {
    double horizon = -1.0;
    switch (id) {
        case MetricId::DeSitter:
            horizon = KottlerCosmologicalHorizonRadius(0.0, lambda);
            break;
        case MetricId::SchwarzschildDeSitter:
            horizon = KottlerCosmologicalHorizonRadius(mass, lambda);
            break;
        case MetricId::Minkowski:
        case MetricId::Schwarzschild:
        case MetricId::Kerr:
        case MetricId::ReissnerNordstrom:
        case MetricId::KerrNewman:
        case MetricId::MorrisThorne:
        case MetricId::Alcubierre:
            return std::nullopt;
        default:
            SIRIUS_ASSERT(false);
            return std::nullopt;
    }
    if (std::isfinite(horizon) && horizon > 0.0) return horizon;
    return std::nullopt;
}

enum class MetricObserverRadiusIssue {
    None,
    NaturalScale,
    CosmologicalHorizon,
};

[[nodiscard]] inline MetricObserverRadiusIssue MetricObserverRadiusIssueFor(
    MetricId id, double mass, double lambda, double observer_radius, double throat_radius,
    double bubble_radius, double bubble_sigma) noexcept {
    const double scale =
        MetricSceneLengthScale(id, mass, throat_radius, bubble_radius, bubble_sigma);
    if (!std::isfinite(observer_radius) || !std::isfinite(scale) || !(scale > 0.0) ||
        observer_radius < 5.0 * scale || observer_radius > 1000.0 * scale) {
        return MetricObserverRadiusIssue::NaturalScale;
    }
    if (const auto cosmological_horizon = MetricCosmologicalHorizonRadius(id, mass, lambda);
        cosmological_horizon.has_value() &&
        !(observer_radius <= kMaxCosmologicalObserverFraction * *cosmological_horizon)) {
        return MetricObserverRadiusIssue::CosmologicalHorizon;
    }
    return MetricObserverRadiusIssue::None;
}

// A very diffuse profile makes the two tanh terms nearly coincident in fp32; a
// very sharp profile makes the wall too narrow relative to R. This bound keeps
// the smallest named scale at least one percent of the scene scale used by the
// shared fp32 state. The separate max(R, 1/sigma) authority places the outer
// boundary in the asymptotic region.
[[nodiscard]] constexpr std::optional<std::string_view> AlcubierreScaleIssue(
    double bubble_radius, double bubble_sigma) noexcept {
    const double sigma_radius = bubble_sigma * bubble_radius;
    if (!(sigma_radius >= kMinAlcubierreSigmaRadius && sigma_radius <= kMaxAlcubierreSigmaRadius)) {
        return "Alcubierre requires 0.1 <= bubble_sigma*bubble_radius <= 100 so its radius "
               "and inverse-wall scales remain jointly resolved";
    }
    return std::nullopt;
}

// Identity compatibility for parameters whose schema defaults must remain
// populated even when their owning exotic metric is inactive. A canonical
// inactive default is harmless; changing an irrelevant value is an operator
// error rather than a silently ignored request.
[[nodiscard]] constexpr std::optional<std::string_view> MetricSpecificParameterIssue(
    MetricId id, double throat_radius, bool one_sheet_topology, double warp_velocity,
    double bubble_radius, double bubble_sigma) noexcept {
    switch (id) {
        case MetricId::Minkowski:
        case MetricId::Schwarzschild:
        case MetricId::Kerr:
        case MetricId::ReissnerNordstrom:
        case MetricId::KerrNewman:
        case MetricId::DeSitter:
        case MetricId::SchwarzschildDeSitter:
        case MetricId::MorrisThorne:
        case MetricId::Alcubierre:
            break;
        default:
            return "unknown metric identity";
    }
    if (id != MetricId::MorrisThorne &&
        (throat_radius != kDefaultMorrisThorneThroatRadius || !one_sheet_topology)) {
        return "throat radius and wormhole topology apply only to Morris-Thorne";
    }
    if (id != MetricId::Alcubierre && (warp_velocity != kDefaultAlcubierreWarpVelocity ||
                                       bubble_radius != kDefaultAlcubierreBubbleRadius ||
                                       bubble_sigma != kDefaultAlcubierreBubbleSigma)) {
        return "warp velocity, bubble radius, and bubble sigma apply only to Alcubierre";
    }
    return std::nullopt;
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
            if (!(lambda > 0.0)) {
                return "Schwarzschild-de-Sitter requires positive lambda";
            }
            break;
        case MetricId::DeSitter:
            if (has_spin || has_charge) return "de-Sitter requires spin and charge to be zero";
            if (!(lambda > 0.0)) return "de-Sitter requires positive lambda";
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

// The named Schwarzschild-de Sitter scene is restricted to the black-hole
// Kottler sector with distinct black-hole and cosmological horizons. At the
// Nariai value 9 Lambda M^2 = 1 the static region degenerates; above it there
// is no black-hole exterior to advertise as Schwarzschild-de Sitter rendering.
[[nodiscard]] constexpr std::optional<std::string_view> MetricHorizonIssue(MetricId id, double mass,
                                                                           double lambda) noexcept {
    if (id != MetricId::SchwarzschildDeSitter) return std::nullopt;
    if (!(mass > 0.0) || !(lambda > 0.0)) {
        return "Schwarzschild-de-Sitter requires positive mass and lambda";
    }
    if (!(9.0 * lambda * mass * mass < 1.0)) {
        return "Schwarzschild-de-Sitter requires 9*lambda*mass^2 < 1 (sub-Nariai sector)";
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
         "positive lambda",
         true,
         false},
        {MetricId::SchwarzschildDeSitter,
         "Schwarzschild-de-Sitter",
         {"SchwarzschildDeSitter", "Schwarzschild-de Sitter", "SdS", nullptr, nullptr},
         "mass, positive lambda; 9 lambda mass^2 < 1",
         true,
         false},
        // CPU support arrives through the exact isotropic Cartesian chart of
        // the zero-tidal Ellis member (one output sheet, throat as the capture
        // surface); the areal-radius family remains an independent authority
        // and is never driven by the Cartesian tracer directly.
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
         "warp velocity, bubble radius, inverse wall scale sigma",
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
