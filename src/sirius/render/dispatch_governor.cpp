// Dispatch governor implementation. See dispatch_governor.h for the contract.

#include "sirius/render/dispatch_governor.h"

#include "sirius/base/contracts.h"

#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <string>

namespace sirius::render {

using base::ErrorDomain;
using base::Expected;
using base::Fail;

int BandController::NextRows(int remaining_rows, int band_width) const {
    SIRIUS_PRE(remaining_rows > 0);
    SIRIUS_PRE(band_width > 0);
    if (!Enabled()) {
        return remaining_rows;
    }
    const auto rows =
        static_cast<int>(std::min<std::int64_t>(pixels_ / band_width, remaining_rows));
    return std::max(rows, 1);
}

void BandController::Record(std::int64_t dispatched_pixels, double measured_ms) {
    SIRIUS_PRE(dispatched_pixels > 0);
    if (!Enabled()) {
        return;
    }
    // A measurement below clock resolution carries no rate information; take
    // the capped growth step rather than dividing by (near) zero.
    double ratio = kBandGrowthCap;
    if (measured_ms > 0.0) {
        ratio = std::min(target_ms_ / measured_ms, kBandGrowthCap);
    }
    // The base is the work actually dispatched, never the area the controller
    // wanted: a truncated tail band must only speak for itself.
    const double scaled = std::floor(static_cast<double>(dispatched_pixels) * ratio);
    pixels_ = std::clamp(static_cast<std::int64_t>(scaled), std::int64_t{1}, max_pixels_);
}

Expected<double> ResolveDispatchTargetMs() {
    const char* override_ms = std::getenv("SIRIUS_DISPATCH_TARGET_MS");
    if (override_ms == nullptr || override_ms[0] == '\0') {
        return kDefaultDispatchTargetMs;
    }
    char* end = nullptr;
    errno = 0;
    const double ms = std::strtod(override_ms, &end);
    if (end == override_ms || *end != '\0' || errno == ERANGE || !std::isfinite(ms) || ms < 0.0) {
        return Fail(ErrorDomain::kDevice, "resolve dispatch target",
                    "SIRIUS_DISPATCH_TARGET_MS='" + std::string(override_ms) +
                        "' is not a finite non-negative millisecond count (0 disables banding)");
    }
    return ms;
}

}  // namespace sirius::render
