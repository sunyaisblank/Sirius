// Memory governor implementation. Pure arithmetic over a reported budget; no
// device handles cross this boundary, so it is unit-testable without Vulkan.

#include "sirius/render/memory_governor.h"

#include "sirius/base/contracts.h"

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <format>
#include <limits>
#include <string>

namespace sirius::render {

using base::ErrorDomain;
using base::Expected;
using base::Fail;

Expected<TilePlan> DeriveTilePlan(std::uint64_t budget_bytes, int image_width, int image_height,
                                  std::uint64_t fixed_overhead_bytes, int tile_edge_cap) {
    SIRIUS_PRE(image_width > 0 && image_height > 0);
    SIRIUS_PRE(tile_edge_cap >= kMinTileEdge);

    const auto usable =
        static_cast<std::uint64_t>(static_cast<double>(budget_bytes) * kResidencyFraction);

    if (usable <= fixed_overhead_bytes) {
        return Fail(ErrorDomain::kDevice, "derive tile plan",
                    std::format("budget {} bytes (usable {} at fraction {:.2f}) cannot seat the "
                                "fixed device residency of {} bytes",
                                budget_bytes, usable, kResidencyFraction, fixed_overhead_bytes));
    }

    const std::uint64_t tile_budget = usable - fixed_overhead_bytes;
    const std::uint64_t max_area = tile_budget / kTileWorkingSetBytesPerPixel;

    // Largest square edge whose area fits, then clamp to the image extent and the
    // sanity cap. floor(sqrt) never over-commits because the square only shrinks.
    auto edge = static_cast<int>(std::floor(std::sqrt(static_cast<double>(max_area))));
    edge = std::min(edge, std::min(image_width, image_height));
    edge = std::min(edge, std::min(kMaxTileEdge, tile_edge_cap));

    if (edge < kMinTileEdge) {
        return Fail(ErrorDomain::kDevice, "derive tile plan",
                    std::format("budget {} bytes seats only a {}x{} tile after {} bytes overhead; "
                                "minimum viable tile edge is {}",
                                budget_bytes, edge, edge, fixed_overhead_bytes, kMinTileEdge));
    }

    TilePlan plan;
    plan.tile_edge = edge;
    plan.budget_bytes = budget_bytes;
    plan.usable_bytes = usable;
    plan.fixed_overhead_bytes = fixed_overhead_bytes;
    plan.tile_working_set_bytes = static_cast<std::uint64_t>(edge) *
                                  static_cast<std::uint64_t>(edge) * kTileWorkingSetBytesPerPixel;
    return plan;
}

Expected<std::uint64_t> ResolveBudgetBytes(std::uint64_t render_memory_bytes) {
    if (const char* override_mb = std::getenv("SIRIUS_MEMORY_BUDGET_MB");
        override_mb != nullptr && override_mb[0] != '\0') {
        char* end = nullptr;
        errno = 0;
        const double mb = std::strtod(override_mb, &end);
        constexpr double kBytesPerMiB = 1024.0 * 1024.0;
        const bool valid =
            end != override_mb && *end == '\0' && errno != ERANGE && std::isfinite(mb) &&
            mb > 0.0 &&
            mb <= static_cast<double>(std::numeric_limits<std::uint64_t>::max()) / kBytesPerMiB;
        if (!valid) {
            return Fail(ErrorDomain::kConfiguration, "resolve memory budget",
                        "SIRIUS_MEMORY_BUDGET_MB='" + std::string(override_mb) +
                            "' is not a finite positive MiB count in range");
        }
        const auto bytes = static_cast<std::uint64_t>(mb * kBytesPerMiB);
        if (bytes == 0) {
            return Fail(ErrorDomain::kConfiguration, "resolve memory budget",
                        "SIRIUS_MEMORY_BUDGET_MB resolves below one byte");
        }
        return bytes;
    }
    return render_memory_bytes;
}

}  // namespace sirius::render
