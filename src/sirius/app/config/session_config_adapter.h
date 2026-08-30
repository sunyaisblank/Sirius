#pragma once

#include "sirius/app/config/config_schema.h"
#include "sirius/base/error.h"
#include "sirius/render/session/render_session.h"

namespace sirius::app {

// Resolve a named finish first, then apply only controls the operator actually
// supplied. This keeps preset identity intact across JSON/CLI projection.
[[nodiscard]] base::Expected<render::FilmConfig> ProjectFilmFinishConfig(
    const FilmFinishConfig& config);

// Cross the operator/renderer boundary exactly once. Strings and compatibility
// spellings end here; the render layer receives closed enums and typed values.
[[nodiscard]] base::Expected<render::SessionConfig> MakeSessionConfig(const SiriusConfig& config);

}  // namespace sirius::app
