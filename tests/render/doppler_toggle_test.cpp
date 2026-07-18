// Doppler beaming toggle tests (specification P4).
//
// The disk g-factor separates the gravitational redshift from the orbital
// Doppler term (tracer gfactor_grav, gfactor_gamma, gfactor_v_orb). The toggle
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

namespace {

using namespace sirius::core;
using namespace sirius::backend;

// Ratio of the approaching/receding brightness asymmetry: |L_left - L_right| /
// (L_left + L_right), where luminance per disk-hit pixel is the beamed observed
// flux g^4 T_emit^4. Rendered near edge-on at the given Doppler setting.
double DiskAsymmetry(bool doppler_beaming, int& disk_hits) {
    KerrSchildParams p;
    p.M = 1.0;
    p.a = 0.9;
    KerrSchildFamily metric(p);

    TracerConfig tc;
    tc.escape_radius = 100.0f;
    tc.horizon_factor = 1.05f;
    tc.max_steps = 20000;
    tc.enable_disk = true;
    tc.disk_inner = 2.32f;  // Prograde ISCO for a = 0.9.
    tc.disk_outer = 20.0f;
    tc.disk_temperature_inner = 1.0f;
    tc.doppler_beaming = doppler_beaming;
    tc.integrator.initial_step = 0.1f;
    tc.integrator.max_step = 2.0f;
    GeodesicTracer tracer(&metric, tc);

    const int width = 120, height = 68;
    CameraConfig cam;
    cam.r = 30.0;
    cam.theta = 80.0 * M_PI / 180.0;  // ~80 deg inclination (near edge-on).
    cam.phi = 0.0;
    cam.fov = 55.0f;
    cam.width = width;
    cam.height = height;
    PinholeCamera camera(cam);

    double lum_left = 0.0, lum_right = 0.0;
    disk_hits = 0;
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            CameraRay ray = camera.GenerateRay(x, y, 0.5f, 0.5f);
            TraceResult r = tracer.Trace(ray);
            if (r.outcome != TraceResult::Outcome::DiskHit) continue;
            disk_hits++;
            double g = r.redshift;
            double t_emit = r.disk_temperature;
            double lum = std::pow(g, 4.0) * std::pow(t_emit, 4.0);  // Observed beamed flux.
            if (x < width / 2)
                lum_left += lum;
            else
                lum_right += lum;
        }
    }
    double total = lum_left + lum_right;
    return (total > 1e-30) ? std::abs(lum_left - lum_right) / total : 0.0;
}

// -----------------------------------------------------------------------------
// With Doppler beaming on the disk halves differ strongly; with it off the
// asymmetry collapses toward zero (only gravitational redshift remains, which is
// symmetric about the spin axis at fixed radius).
// -----------------------------------------------------------------------------
TEST(DopplerToggleTest, SuppressionCollapsesDiskAsymmetry) {
    int hits_on = 0, hits_off = 0;
    double asym_on = DiskAsymmetry(true, hits_on);
    double asym_off = DiskAsymmetry(false, hits_off);

    std::cout << "[doppler] asym_on=" << asym_on << " (" << hits_on << " hits) asym_off=" << asym_off
              << " (" << hits_off << " hits) collapse=" << (asym_on / std::max(asym_off, 1e-9))
              << "\n";

    ASSERT_GT(hits_on, 100) << "too few disk hits to measure asymmetry";
    EXPECT_EQ(hits_on, hits_off) << "the toggle must not change the geometry, only the shading";

    // The full-physics asymmetry is large.
    EXPECT_GT(asym_on, 0.1);
    // Suppression collapses it by a wide margin.
    EXPECT_LT(asym_off, 0.2 * asym_on);
}

}  // namespace
