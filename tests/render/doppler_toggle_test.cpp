// Doppler beaming toggle tests (specification P4).
//
// The disk transfer records exact circular-emitter and ZAMO frequency ratios. The toggle
// (TracerConfig::doppler_beaming) suppresses the Doppler asymmetry - the film's
// artistic choice (James et al. 2015, CQG 32 065001, section 5) - while keeping
// the gravitational redshift. Viewed near edge-on, the approaching and receding
// halves of the disk carry a large brightness asymmetry with the toggle on that
// collapses with it off. Default true reproduces the pinned render.

#include "sirius/backend/cpu/geodesic_tracer.h"
#include "sirius/core/camera.h"
#include "sirius/core/metrics/kerr_schild_family.h"

#include <gtest/gtest.h>

#include <cmath>
#include <iostream>
#include <numbers>

namespace {

using namespace sirius::core;
using namespace sirius::backend;

struct AsymmetryMeasurement {
    double observed = 0.0;
    double doppler_factor = 0.0;
    int disk_hits = 0;
};

// Measures both the rendered g^4 T^4 asymmetry and the isolated Doppler factor
// (g/g_grav)^4. The latter is averaged per hit on each image half so Kerr
// frame-dragging/lensing cannot masquerade as an unsuppressed emitter velocity
// merely by mapping different radii or hit counts onto the two halves.
AsymmetryMeasurement DiskAsymmetry(bool doppler_beaming) {
    KerrSchildParams p;
    p.M = 1.0;
    p.a = 0.9;
    KerrSchildFamily metric(p);

    TracerConfig tc;
    tc.escape_radius = 100.0f;
    tc.horizon_factor = 1.05f;
    tc.max_steps = 20000;
    tc.enable_disk = true;
    tc.disk_inner = AccretionDiskD::ComputeIsco(0.9);
    tc.disk_outer = 20.0f;
    tc.disk_temperature_inner = 1.0f;
    tc.doppler_beaming = doppler_beaming;
    tc.integrator.initial_step = 0.1f;
    tc.integrator.max_step = 2.0f;
    GeodesicTracer tracer(&metric, tc);

    const int width = 120, height = 68;
    CameraConfig cam;
    cam.r = 30.0;
    cam.theta = 80.0 * std::numbers::pi / 180.0;  // ~80 deg inclination (near edge-on).
    cam.phi = 0.0;
    cam.fov = 55.0f;
    cam.width = width;
    cam.height = height;
    PinholeCamera camera(cam);

    double lum_left = 0.0, lum_right = 0.0;
    double factor_left = 0.0, factor_right = 0.0;
    int hits_left = 0, hits_right = 0;
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            CameraRay ray = camera.GenerateRay(x, y, 0.5f, 0.5f);
            TraceResult r = tracer.Trace(ray);
            if (r.outcome != TraceResult::Outcome::DiskHit) continue;
            const bool left = x < width / 2;
            if (left)
                ++hits_left;
            else
                ++hits_right;
            double g = r.redshift;
            double t_emit = r.disk_temperature;
            double lum = std::pow(g, 4.0) * std::pow(t_emit, 4.0);  // Observed beamed flux.
            // Measure the branch selected by the toggle.  The trace result
            // intentionally retains both physical transfer witnesses even
            // when shading selects the ZAMO diagnostic branch.
            const double doppler_factor = std::pow(r.redshift / r.zamo_disk_redshift, 4.0);
            if (left) {
                lum_left += lum;
                factor_left += doppler_factor;
            } else {
                lum_right += lum;
                factor_right += doppler_factor;
            }
        }
    }
    const double total = lum_left + lum_right;
    const double mean_factor_left = factor_left / std::max(hits_left, 1);
    const double mean_factor_right = factor_right / std::max(hits_right, 1);
    const double factor_total = mean_factor_left + mean_factor_right;
    return AsymmetryMeasurement{
        .observed = (total > 1e-30) ? std::abs(lum_left - lum_right) / total : 0.0,
        .doppler_factor = (factor_total > 1e-30)
                              ? std::abs(mean_factor_left - mean_factor_right) / factor_total
                              : 0.0,
        .disk_hits = hits_left + hits_right,
    };
}

// -----------------------------------------------------------------------------
// With Doppler beaming on the disk halves differ strongly; with it off the
// asymmetry collapses toward zero (only gravitational redshift remains, which is
// symmetric about the spin axis at fixed radius).
// -----------------------------------------------------------------------------
TEST(DopplerToggleTest, SuppressionCollapsesDiskAsymmetry) {
    const AsymmetryMeasurement on = DiskAsymmetry(true);
    const AsymmetryMeasurement off = DiskAsymmetry(false);

    std::cout << "[doppler] observed_on=" << on.observed << " observed_off=" << off.observed
              << " factor_on=" << on.doppler_factor << " factor_off=" << off.doppler_factor
              << " hits=" << on.disk_hits << "\n";

    ASSERT_GT(on.disk_hits, 100) << "too few disk hits to measure asymmetry";
    EXPECT_EQ(on.disk_hits, off.disk_hits)
        << "the toggle must not change the geometry, only the shading";

    EXPECT_GT(on.doppler_factor, 0.1);
    EXPECT_LT(off.doppler_factor, 1.0e-6)
        << "emitter-velocity beaming remains after Doppler suppression";
    EXPECT_LT(off.observed, 0.25 * on.observed)
        << "the rendered image did not substantially reduce its total asymmetry";
}

}  // namespace
