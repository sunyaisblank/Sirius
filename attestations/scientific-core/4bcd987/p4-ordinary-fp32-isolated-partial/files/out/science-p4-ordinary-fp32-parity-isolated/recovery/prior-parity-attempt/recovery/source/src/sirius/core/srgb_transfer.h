#pragma once

// Single host authority for the IEC 61966-2-1 sRGB opto-electronic encoding
// curve and its clipped 8-bit display quantisation. Writers own when the
// transfer is applied; they do not own independent copies of its mathematics.

#include <algorithm>
#include <cmath>
#include <concepts>
#include <cstdint>
#include <limits>
#include <optional>

namespace sirius::core::colour {

template <std::floating_point T>
[[nodiscard]] inline T EncodeSrgbChannel(T linear) {
    if (!std::isfinite(linear)) return std::numeric_limits<T>::quiet_NaN();

    constexpr T kLinearBreakpoint = static_cast<T>(0.0031308L);
    constexpr T kLinearSlope = static_cast<T>(12.92L);
    constexpr T kPowerScale = static_cast<T>(1.055L);
    constexpr T kPowerOffset = static_cast<T>(0.055L);
    constexpr T kPowerExponent = static_cast<T>(1.0L / 2.4L);
    if (linear <= kLinearBreakpoint) return kLinearSlope * linear;
    return kPowerScale * std::pow(linear, kPowerExponent) - kPowerOffset;
}

template <std::floating_point T>
[[nodiscard]] inline T EncodeClippedSrgbChannel(T linear) {
    if (!std::isfinite(linear)) return std::numeric_limits<T>::quiet_NaN();
    return EncodeSrgbChannel(std::clamp(linear, T{0}, T{1}));
}

template <std::floating_point T>
[[nodiscard]] inline std::optional<std::uint8_t> TryEncodeSrgb8(T linear) {
    if (!std::isfinite(linear)) return std::nullopt;

    constexpr T kCodeValueMax = static_cast<T>(255);
    constexpr T kRoundToNearest = static_cast<T>(0.5L);
    const T encoded = EncodeClippedSrgbChannel(linear);
    const T quantised = std::clamp(encoded * kCodeValueMax + kRoundToNearest, T{0}, kCodeValueMax);
    return static_cast<std::uint8_t>(quantised);
}

}  // namespace sirius::core::colour
