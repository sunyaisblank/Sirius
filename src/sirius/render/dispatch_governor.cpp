// Dispatch governor implementation. See dispatch_governor.h for the contract.

#include "sirius/render/dispatch_governor.h"

#include "sirius/base/contracts.h"

#include <cmath>
#include <cstdlib>
#include <string>

namespace sirius::render {

using base::ErrorDomain;
using base::Expected;
using base::Fail;

int BandController::NextRows(int remaining_rows) const {
    SIRIUS_PRE(remaining_rows > 0);
    if (!Enabled()) {
        return remaining_rows;
    }
    return std::clamp(rows_, 1, remaining_rows);
}

void BandController::Record(double measured_ms) {
    if (!Enabled()) {
        return;
    }
    // A measurement below clock resolution carries no rate information; take
    // the capped growth step rather than dividing by (near) zero.
    double ratio = kBandGrowthCap;
    if (measured_ms > 0.0) {
        ratio = std::min(target_ms_ / measured_ms, kBandGrowthCap);
    }
    const double scaled = std::floor(static_cast<double>(rows_) * ratio);
    rows_ = std::clamp(static_cast<int>(scaled), 1, max_rows_);
}

Expected<double> ResolveDispatchTargetMs() {
    const char* override_ms = std::getenv("SIRIUS_DISPATCH_TARGET_MS");
    if (override_ms == nullptr || override_ms[0] == '\0') {
        return kDefaultDispatchTargetMs;
    }
    char* end = nullptr;
    const double ms = std::strtod(override_ms, &end);
    if (end == override_ms || *end != '\0' || std::isnan(ms) || ms < 0.0) {
        return Fail(ErrorDomain::kDevice, "resolve dispatch target",
                    "SIRIUS_DISPATCH_TARGET_MS='" + std::string(override_ms) +
                        "' is not a non-negative millisecond count (0 disables banding)");
    }
    return ms;
}

}  // namespace sirius::render
