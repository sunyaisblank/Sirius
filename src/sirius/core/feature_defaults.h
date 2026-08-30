#pragma once

// Canonical scalar defaults shared by the app schema and typed render session.
// An inactive feature therefore has one unambiguous neutral parameter state at
// every boundary.

namespace sirius::core {

inline constexpr float kDefaultBloomIntensity = 0.3f;
inline constexpr float kDefaultBloomThreshold = 0.3f;

inline constexpr float kDefaultVolumetricHOverR = 0.1f;
inline constexpr float kDefaultVolumetricHPower = 0.25f;
inline constexpr float kDefaultVolumetricTauMidplane = 10.0f;
inline constexpr int kDefaultVolumetricSamples = 64;

inline constexpr float kDefaultMotionBlurShutterTime = 0.1f;
inline constexpr int kDefaultMotionBlurSamples = 3;

}  // namespace sirius::core
