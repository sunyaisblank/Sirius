// MTQT001A.h - Quaternion Mathematics (GLM wrapper)

#pragma once

#include <glm/gtc/quaternion.hpp>

/// @brief Quaternion for 3D rotations (wraps GLM)
using Quat = glm::quat;

#ifndef SIRIUS_NO_LEGACY_ALIASES
// Legacy code may include <Quaternion.h>
#endif
