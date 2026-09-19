#pragma once

#include "sirius/render/session/render_session.h"

namespace sirius::test {
// The existing bounded CPU detector witness, shared with device image parity.
// Critical-ray precision is judged separately by the independent transport fixtures.
inline void ConfigureMovingKerrDetector(render::SessionConfig& config) {
    config.metric_id = core::MetricId::Kerr;
    config.black_hole_mass = 1;
    config.black_hole_spin = .7;
    config.observer_distance = 50;
    config.observer_inclination = render::SessionConfig{}.observer_inclination;
    config.camera_fov = 2;
    config.camera_beta_forward = .1;
    config.camera_beta_up = .8;
    config.camera_beta_right = 0;
    config.lens_type = core::LensType::ThinLens;
    config.camera_focus_distance = 50;
    config.enable_disk = false;
    config.enable_bloom = false;
    config.point_starfield = true;
    config.point_starfield_config.star_count = 100000;
    config.point_starfield_config.brightness_scale = 1e-5f;
}
}  // namespace sirius::test
