#pragma once

// Operator-visible accretion-disk defaults shared across the JSON boundary and
// typed render session. Keeping the physical temperature scale here prevents a
// disabled disk from carrying a seemingly meaningful but unconsumed override.

namespace sirius::core {

inline constexpr float kDefaultDiskTemperatureKelvin = 50000.0f;

}  // namespace sirius::core
