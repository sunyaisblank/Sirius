// Memory governor implementation. Pure arithmetic over a reported budget; no
// device handles cross this boundary, so it is unit-testable without Vulkan.

#include "sirius/render/memory_governor.h"

#include "sirius/base/contracts.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <format>

namespace sirius::render {

using base::ErrorDomain;
using base::Expected;
using base::Fail;

Expected<TilePlan> DeriveTilePlan(std::uint64_t budget_bytes, int image_width, int image_height,
                                  std::uint64_t fixed_overhead_bytes) {
    SIRIUS_PRE(image_width > 0 && image_height > 0);

    const auto usable = static_cast<std::uint64_t>(static_cast<double>(budget_bytes) *
                                                   kResidencyFraction);

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
    edge = std::min(edge, kMaxTileEdge);

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
    plan.tile_working_set_bytes =
        static_cast<std::uint64_t>(edge) * static_cast<std::uint64_t>(edge) *
        kTileWorkingSetBytesPerPixel;
    return plan;
}

std::uint64_t ResolveBudgetBytes(std::uint64_t device_local_bytes) {
    if (const char* override_mb = std::getenv("SIRIUS_MEMORY_BUDGET_MB");
        override_mb != nullptr && override_mb[0] != '\0') {
        char* end = nullptr;
        const double mb = std::strtod(override_mb, &end);
        if (end != override_mb && mb > 0.0) {
            return static_cast<std::uint64_t>(mb * 1024.0 * 1024.0);
        }
    }
    return device_local_bytes;
}

}  // namespace sirius::render
