#pragma once

// Single host authority for converting CIE XYZ, relative to the sRGB D65
// white point, to linear-light sRGB.  The rational coefficients are the exact
// inverse matrix implied by the IEC sRGB primary chromaticities and D65 white;
// retaining them as ratios prevents a rounded decimal copy from becoming a
// second colour-space definition.

#include <concepts>

namespace sirius::core::colour {

template <std::floating_point T>
struct LinearSrgb {
    T r;
    T g;
    T b;
};

template <std::floating_point T>
[[nodiscard]] constexpr LinearSrgb<T> XyzD65ToLinearSrgb(T x, T y, T z) noexcept {
    constexpr T kRedFromX = static_cast<T>(12831.0L / 3959.0L);
    constexpr T kRedFromY = static_cast<T>(-329.0L / 214.0L);
    constexpr T kRedFromZ = static_cast<T>(-1974.0L / 3959.0L);
    constexpr T kGreenFromX = static_cast<T>(-851781.0L / 878810.0L);
    constexpr T kGreenFromY = static_cast<T>(1648619.0L / 878810.0L);
    constexpr T kGreenFromZ = static_cast<T>(36519.0L / 878810.0L);
    constexpr T kBlueFromX = static_cast<T>(705.0L / 12673.0L);
    constexpr T kBlueFromY = static_cast<T>(-2585.0L / 12673.0L);
    constexpr T kBlueFromZ = static_cast<T>(705.0L / 667.0L);

    return {
        kRedFromX * x + kRedFromY * y + kRedFromZ * z,
        kGreenFromX * x + kGreenFromY * y + kGreenFromZ * z,
        kBlueFromX * x + kBlueFromY * y + kBlueFromZ * z,
    };
}

}  // namespace sirius::core::colour
