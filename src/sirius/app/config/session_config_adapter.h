#pragma once

#include "sirius/app/config/config_schema.h"
#include "sirius/base/error.h"
#include "sirius/render/session/render_session.h"

namespace sirius::app {

// Cross the operator/renderer boundary exactly once. Strings and compatibility
// spellings end here; the render layer receives closed enums and typed values.
[[nodiscard]] base::Expected<render::SessionConfig> MakeSessionConfig(const SiriusConfig& config);

}  // namespace sirius::app
