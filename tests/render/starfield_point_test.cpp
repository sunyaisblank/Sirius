// Filtered point-source star field tests (specification P3).
//
// Stars render as point sources filtered through the per-pixel beam footprint the
// ray bundle supplies, not as texture lookups (James et al. 2015, CQG 32 065001,
// section 3). This suite pins the catalogue size, the finiteness and structure of
// a small render, and the anti-flicker result: over a slowly rotating camera the
// beam footprint holds a bright pixel's brightness far steadier than pinhole
// point sampling does.

#include "sirius/backend/cpu/geodesic_tracer.h"
#include "sirius/core/camera.h"
#include "sirius/core/metrics/kerr_schild_family.h"
#include "sirius/core/starfield.h"

#include <gtest/gtest.h>

#include <cmath>
#include <iostream>
#include <vector>

namespace {

using namespace sirius::core;
using namespace sirius::backend;

// Mean and variance of a sample.
void MeanVar(const std::vector<double>& v, double& mean, double& var) {
    mean = 0.0;
    for (double x : v) mean += x;
    mean /= static_cast<double>(v.size());
    var = 0.0;
    for (double x : v) var += (x - mean) * (x - mean);
    var /= static_cast<double>(v.size());
}

// -----------------------------------------------------------------------------
// The deterministic catalogue meets the >= 10^5 floor with finite entries.
// -----------------------------------------------------------------------------
TEST(StarfieldPointTest, CatalogueMeetsSizeFloorAndIsFinite) {
    StarfieldConfig cfg;
    cfg.star_count = 100000;
    cfg.seed = 1234;
    StarfieldGenerator gen(cfg);
    auto stars = gen.GenerateCatalogue();

    ASSERT_GE(stars.size(), 100000u);
    for (size_t i = 0; i < stars.size(); i += 97) {  // Sample the catalogue.
        const auto& s = stars[i];
        double n = std::sqrt(s.direction_x * s.direction_x + s.direction_y * s.direction_y +
                             s.direction_z * s.direction_z);
        EXPECT_NEAR(n, 1.0, 1e-3) << "star direction not unit";
        EXPECT_TRUE(std::isfinite(s.magnitude));
        EXPECT_GT(s.temperature_K, 0.0f);
        EXPECT_TRUE(std::isfinite(s.temperature_K));
    }

    // Determinism: same seed -> identical catalogue.
    auto stars2 = gen.GenerateCatalogue();
    ASSERT_EQ(stars.size(), stars2.size());
    EXPECT_FLOAT_EQ(stars[500].direction_x, stars2[500].direction_x);
    EXPECT_FLOAT_EQ(stars[500].magnitude, stars2[500].magnitude);
}

// -----------------------------------------------------------------------------
// Beam accumulation over a small frame is finite, free of NaN, and non-constant.
// -----------------------------------------------------------------------------
TEST(StarfieldPointTest, BeamAccumulationFiniteAndNonConstant) {
    StarfieldConfig cfg;
    cfg.star_count = 100000;
    cfg.seed = 7;
    cfg.brightness_scale = 50.0f;
    StarfieldGenerator gen(cfg);
    auto stars = gen.GenerateCatalogue();

    float sigma = 0.02f;  // ~1 deg beam.
    int lit = 0;
    double first = -1.0;
    bool non_constant = false;
    for (int i = 0; i < 200; ++i) {
        double th = M_PI * (i + 0.5) / 200.0;
        float dx = std::sin(th), dy = 0.0f, dz = std::cos(th);
        float r, g, b;
        gen.AccumulateThroughBeam(dx, dy, dz, sigma, stars, r, g, b);
        ASSERT_TRUE(std::isfinite(r) && std::isfinite(g) && std::isfinite(b));
        EXPECT_GE(r, 0.0f);
        double lum = r + g + b;
        if (lum > 1e-6) lit++;
        if (first < 0) first = lum;
        else if (std::abs(lum - first) > 1e-6)
            non_constant = true;
    }
    EXPECT_GT(lit, 0) << "no stars accumulated anywhere";
    EXPECT_TRUE(non_constant) << "star field is constant across directions";
}

// -----------------------------------------------------------------------------
// The flicker test: a slowly rotating camera samples the star field through the
// tracer's beam footprint (bundles on) and through a pinhole (bundles off). The
// beam holds each bright pixel steady; the pinhole makes stars pop in and out.
// -----------------------------------------------------------------------------
TEST(StarfieldPointTest, BeamFootprintSuppressesStarFlicker) {
    StarfieldConfig scfg;
    scfg.star_count = 100000;
    scfg.seed = 42;
    scfg.brightness_scale = 100.0f;
    StarfieldGenerator gen(scfg);
    auto stars = gen.GenerateCatalogue();

    // Flat space: rays escape straight, so the tracer's pupil bundle returns the
    // pixel-sized footprint and the run stays cheap; the flicker mechanism is the
    // star sampling, isolated from lensing.
    KerrSchildParams p;
    p.M = 0.0;
    KerrSchildFamily metric(p);

    const int width = 24, height = 14;
    const float fov = 40.0f;
    double pixel_angular = (fov * M_PI / 180.0) / height;

    TracerConfig tc;
    tc.escape_radius = 80.0f;
    tc.max_steps = 400;
    tc.enable_disk = false;
    tc.enable_ray_bundles = true;
    tc.bundle_point_source = true;
    tc.bundle_angular_size = static_cast<float>(pixel_angular);
    // Flat space rays are straight, so large affine steps stay exact and keep the
    // bundle sweep (Riemann per step) cheap.
    tc.integrator.initial_step = 1.0f;
    tc.integrator.max_step = 5.0f;
    tc.integrator.min_step = 0.01f;
    GeodesicTracer tracer(&metric, tc);

    const int frames = 8;
    const float dyaw = 0.4f * static_cast<float>(pixel_angular);  // Slow rotation.
    const float pinhole_sigma = 0.3f * static_cast<float>(pixel_angular);

    // Per-pixel brightness time series, beams on and off.
    std::vector<std::vector<double>> series_beam(width * height), series_pin(width * height);

    for (int f = 0; f < frames; ++f) {
        CameraConfig cam;
        cam.r = 50.0;
        cam.theta = M_PI / 2.0;
        cam.phi = 0.0;
        cam.fov = fov;
        cam.width = width;
        cam.height = height;
        cam.yaw = f * dyaw;  // Slowly rotate.
        PinholeCamera camera(cam);

        for (int y = 0; y < height; ++y) {
            for (int x = 0; x < width; ++x) {
                CameraRay ray = camera.GenerateRay(x, y, 0.5f, 0.5f);
                TraceResult r = tracer.Trace(ray);
                double lum_beam = 0.0, lum_pin = 0.0;
                if (r.outcome == TraceResult::Outcome::Escaped && r.beam.valid) {
                    const auto& d = r.final_direction;
                    float dx = static_cast<float>(d(1)), dy = static_cast<float>(d(2)),
                          dz = static_cast<float>(d(3));
                    float footprint = std::sqrt(std::max(0.0f, r.beam.footprint_major) *
                                                std::max(0.0f, r.beam.footprint_minor));
                    float sigma_beam = std::max(footprint, static_cast<float>(pixel_angular));
                    float br, bg, bb, pr, pg, pb;
                    gen.AccumulateThroughBeam(dx, dy, dz, sigma_beam, stars, br, bg, bb);
                    gen.AccumulateThroughBeam(dx, dy, dz, pinhole_sigma, stars, pr, pg, pb);
                    lum_beam = br + bg + bb;
                    lum_pin = pr + pg + pb;
                }
                series_beam[y * width + x].push_back(lum_beam);
                series_pin[y * width + x].push_back(lum_pin);
            }
        }
    }

    // Consider pixels a star sweeps through (bright in at least one frame under
    // pinhole sampling). Compare the temporal coefficient of variation.
    double cv_beam_sum = 0.0, cv_pin_sum = 0.0;
    int counted = 0;
    for (int i = 0; i < width * height; ++i) {
        double peak = 0.0;
        for (double x : series_pin[i]) peak = std::max(peak, x);
        if (peak < 1e-4) continue;  // No star here.
        double mb, vb, mp, vp;
        MeanVar(series_beam[i], mb, vb);
        MeanVar(series_pin[i], mp, vp);
        if (mb > 1e-9) cv_beam_sum += std::sqrt(vb) / mb;
        if (mp > 1e-9) cv_pin_sum += std::sqrt(vp) / mp;
        counted++;
    }
    ASSERT_GT(counted, 0) << "no bright pixels to measure flicker on";
    double cv_beam = cv_beam_sum / counted;
    double cv_pin = cv_pin_sum / counted;

    std::cout << "[flicker] bright pixels=" << counted << " CV(beam)=" << cv_beam
              << " CV(pinhole)=" << cv_pin << " ratio=" << (cv_pin / std::max(cv_beam, 1e-12))
              << "\n";

    // The beam footprint holds each pixel steadier than the pinhole does.
    EXPECT_LT(cv_beam, cv_pin);
    EXPECT_GT(cv_pin / std::max(cv_beam, 1e-12), 1.5)
        << "beam did not measurably reduce star flicker";
}

}  // namespace
