// Dispatch governor implementation. See dispatch_governor.h for the contract.

#include "sirius/render/dispatch_governor.h"

#include "sirius/base/contracts.h"

#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <format>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace sirius::render {

using base::ErrorDomain;
using base::Expected;
using base::Fail;

std::array<DispatchRegion, 2> SplitDispatchRegion(const DispatchRegion& region) {
    SIRIUS_PRE(region.width > 0 && region.height > 0 && region.Pixels() > 1);
    DispatchRegion first = region;
    DispatchRegion second = region;
    if (region.height > 1) {
        first.height = region.height / 2;
        second.y += first.height;
        second.height -= first.height;
    } else {
        first.width = region.width / 2;
        second.x += first.width;
        second.width -= first.width;
    }
    return {first, second};
}

int BandController::NextRows(int remaining_rows, int band_width) const {
    SIRIUS_PRE(remaining_rows > 0);
    SIRIUS_PRE(band_width > 0);
    SIRIUS_PRE(band_width <= max_pixels_);
    const auto allowed = safety_fallback_ ? std::int64_t{1} : (Enabled() ? pixels_ : max_pixels_);
    const auto rows = std::min<std::int64_t>(std::max<std::int64_t>(1, allowed / band_width),
                                             std::min(remaining_rows, max_rows_));
    return static_cast<int>(rows);
}

bool BandController::Record(std::int64_t dispatched_pixels, double measured_ms) {
    SIRIUS_PRE(dispatched_pixels > 0);
    if (!std::isfinite(measured_ms) || measured_ms < 0.0 || dispatched_pixels > max_pixels_ ||
        (measured_ms > kDispatchStopMs && dispatched_pixels == 1)) {
        pixels_ = 1;
        return false;
    }
    if (measured_ms == 0.0 || measured_ms > kDispatchStopMs) {
        safety_fallback_ = true;
        const auto reduced_cap =
            measured_ms == 0.0 ? std::int64_t{1} : std::max<std::int64_t>(1, dispatched_pixels / 2);
        safety_pixel_cap_ = std::min(safety_pixel_cap_, reduced_cap);
        pixels_ = safety_pixel_cap_;
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

Expected<void> ExecuteDispatchRegions(
    const DispatchRegion& logical_region, int samples_per_pixel, BandController& bands,
    const std::function<Expected<double>(const DispatchRegion&, const CameraSample&, int)>& submit,
    const std::function<Expected<void>(const DispatchRegion&)>& completed,
    const std::function<bool()>& should_cancel) {
    SIRIUS_PRE(logical_region.width > 0 && logical_region.height > 0);
    SIRIUS_PRE(samples_per_pixel >= 1 && samples_per_pixel <= 4096);
    std::vector<DispatchRegion> pending{logical_region};
    while (!pending.empty()) {
        if (should_cancel && should_cancel()) {
            return Fail(ErrorDomain::kInternal, "render Vulkan frame", "cancelled");
        }
        const DispatchRegion region = pending.back();
        pending.pop_back();
        if (region.Pixels() > bands.SafetyPixelCap()) {
            const auto children = SplitDispatchRegion(region);
            pending.push_back(children[1]);
            pending.push_back(children[0]);
            continue;
        }

        int sample_index = 0;
        bool retry_smaller = false;
        std::optional<base::Error> sample_error;
        ForEachCameraSample(samples_per_pixel, [&](const CameraSample& sample) {
            if (sample_error.has_value() || retry_smaller) return;
            if (should_cancel && should_cancel()) {
                sample_error =
                    base::Error{ErrorDomain::kInternal, "render Vulkan frame", "cancelled"};
                return;
            }
            const auto timing = submit(region, sample, sample_index);
            if (!timing) {
                sample_error = timing.error();
                return;
            }
            if (!bands.Record(region.Pixels(), *timing)) {
                sample_error = base::Error{
                    ErrorDomain::kDevice, "govern Vulkan dispatch",
                    std::format("refusing further submissions after {} ms for {} pixels "
                                "at ({}, {}) in {}x{} sample {} "
                                "(invalid timing or one pixel exceeded {} ms)",
                                *timing, region.Pixels(), region.x, region.y, region.width,
                                region.height, sample_index, kDispatchStopMs)};
                return;
            }
            if (region.Pixels() > bands.SafetyPixelCap()) {
                // Changing width also changes the radiance buffer stride. The
                // partial parent must never reach the completion callback.
                retry_smaller = true;
                return;
            }
            ++sample_index;
        });
        if (sample_error.has_value()) {
            return std::unexpected(std::move(*sample_error));
        }
        if (retry_smaller) {
            pending.push_back(region);
            continue;
        }
        if (should_cancel && should_cancel()) {
            return Fail(ErrorDomain::kInternal, "render Vulkan frame", "cancelled");
        }
        if (auto accepted = completed(region); !accepted) {
            return std::unexpected(accepted.error());
        }
    }
    return {};
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
