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
    SIRIUS_PRE(band_width <= max_pixels_);
    const auto allowed = minimum_fallback_ ? std::int64_t{1} : (Enabled() ? pixels_ : max_pixels_);
    const auto rows = std::min<std::int64_t>(std::max<std::int64_t>(1, allowed / band_width),
                                             std::min(remaining_rows, max_rows_));
    return static_cast<int>(rows);
}

bool BandController::Record(std::int64_t dispatched_pixels, double measured_ms) {
    SIRIUS_PRE(dispatched_pixels > 0);
    if (!std::isfinite(measured_ms) || measured_ms < 0.0 || dispatched_pixels > max_pixels_ ||
        (measured_ms > kDispatchStopMs && dispatched_pixels <= minimum_pixels_)) {
        pixels_ = 1;
        return false;
    }
    if (measured_ms == 0.0 || measured_ms > kDispatchStopMs) {
        minimum_fallback_ = true;
        pixels_ = 1;
        return true;
    }
    if (!Enabled()) return true;
    const double ratio = measured_ms < target_ms_   ? kBandGrowthCap
                         : measured_ms > target_ms_ ? 0.5
                                                    : 1.0;
    const double scaled = std::floor(static_cast<double>(dispatched_pixels) * ratio);
    // Clamp before conversion; the next band also honours row and area caps.
    pixels_ = static_cast<std::int64_t>(std::clamp(scaled, 1.0, static_cast<double>(max_pixels_)));
    return true;
}

Expected<double> ResolveDispatchTargetMs(double default_target_ms) {
    const char* override_ms = std::getenv("SIRIUS_DISPATCH_TARGET_MS");
    if (override_ms == nullptr || override_ms[0] == '\0') {
        return default_target_ms;
    }
    char* end = nullptr;
    errno = 0;
    const double ms = std::strtod(override_ms, &end);
    if (end == override_ms || *end != '\0' || errno == ERANGE || !std::isfinite(ms) || ms < 0.0) {
        return Fail(ErrorDomain::kDevice, "resolve dispatch target",
                    "SIRIUS_DISPATCH_TARGET_MS='" + std::string(override_ms) +
                        "' is not a finite non-negative millisecond count (0 disables adaptation)");
    }
    return ms;
}

}  // namespace sirius::render
