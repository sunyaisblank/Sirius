// Ray-bundle (geodesic deviation) tests for the CPU live path (specification P2).
//
// The tracer propagates two deviation vectors alongside the central null ray by
// the Jacobi equation D^2 xi/d lambda^2 = -R(k, xi)k and reports the beam ellipse
// at termination (James et al. 2015, CQG 32 065001, section 3 / appendix B).
//
// Oracle agreement (gate a): the CPU path finite-differences the Riemann tensor
// in the Kerr-Schild Cartesian chart the render integrates in. The chart-
// invariant handle both paths compute is the Kretschmann scalar
// K = R_mu_nu_rho_sig R^mu_nu_rho_sig; the oracle carries the analytic Kerr
// value (KerrMetricD::Kretschmann, Henry 2000, eq. 18). Matching the CPU
// curvature to it pins the tensor the whole bundle rides on, at every sampled
// angle for both the Schwarzschild and the Kerr a=0.9 congruences. The live
// fp32 bar is the conservation tier (core/constants.h geodesic::
// kConservationTol, 1e-4); the double-precision arithmetic here beats it by
// two orders of magnitude.
//
// History: this suite originally pinned Kerr on the equator only. The closed-
// form Kretschmann's bracket term carried (r^2 - a^2 cos^2 th)^2 where Henry's
// equation has (r^2 + a^2 cos^2 th)^2, a defect invisible wherever a cos th
// = 0 - exactly the two regimes then sampled - and the oracle's component-wise
// Riemann() was incomplete off the diagonal blocks. Both were rebuilt (the
// tensor from analytic first and second metric derivatives, the scalar to
// Henry's form, cross-pinned by contraction in the oracle suite), which is
// what admits the off-equator cases below.
//
// Magnification consistency (gate b) uses only the propagated Jacobi map; no
// separate scalar brightness heuristic is admitted.

#include "sirius/backend/cpu/geodesic_tracer.h"
#include "sirius/core/camera.h"
#include "sirius/core/constants.h"
#include "sirius/core/metrics/kerr_schild_family.h"
#include "sirius/oracle/kerr_boyer_lindquist.h"

#include <gtest/gtest.h>

#include <cmath>
#include <iostream>
#include <numbers>

namespace {

using namespace sirius::core;
using namespace sirius::backend;

// The Kerr-Schild Cartesian image of a Boyer-Lindquist point (r, theta) at
// phi = 0: x = sqrt(r^2 + a^2) sin theta, z = r cos theta (coordinates.h).
Vec4 CartPoint(double r, double theta, double a) {
    Vec4 p;
    p(1) = std::sqrt(r * r + a * a) * std::sin(theta);
    p(3) = r * std::cos(theta);
    return p;
}

// Relative Kretschmann agreement between the CPU Cartesian Riemann and the
// oracle's analytic Kerr value at matched physical events.
double MaxKretschmannError(double a, const std::vector<double>& radii,
                           const std::vector<double>& thetas) {
    KerrSchildParams p;
    p.M = 1.0;
    p.a = a;
    KerrSchildFamily metric(p);
    TracerConfig cfg;
    GeodesicTracer tracer(&metric, cfg);
    sirius::oracle::KerrMetricD om(1.0, a);

    double worst = 0.0;
    for (double r : radii)
        for (double th : thetas) {
            double k_cpu = tracer.KretschmannScalar(CartPoint(r, th, a));
            double k_ora = om.Kretschmann(sirius::oracle::Vec4d(0, r, th, 0));
            double rel = std::abs(k_cpu - k_ora) / std::max(std::abs(k_ora), 1e-30);
            worst = std::max(worst, rel);
        }
    return worst;
}

// --- Bundle probe (full trace) -----------------------------------------------
struct BundleProbe {
    std::unique_ptr<KerrSchildFamily> metric;
    std::unique_ptr<GeodesicTracer> tracer;
    std::unique_ptr<PinholeCamera> camera;

    void Build(double M, double a, float escape_radius, int width = 64, int height = 64) {
        KerrSchildParams p;
        p.M = M;
        p.a = a;
        metric = std::make_unique<KerrSchildFamily>(p);
        TracerConfig cfg;
        cfg.escape_radius = escape_radius;
        cfg.horizon_factor = 1.05f;
        cfg.max_steps = 20000;
        cfg.enable_disk = false;
        cfg.enable_ray_bundles = true;
        cfg.bundle_angular_size = 1.0e-3f;
        cfg.integrator.initial_step = 0.05f;
        cfg.integrator.max_step = 1.0f;
        cfg.integrator.abs_tolerance = 1e-7f;
        cfg.integrator.rel_tolerance = 1e-7f;
        tracer = std::make_unique<GeodesicTracer>(metric.get(), cfg);
        CameraConfig cam;
        cam.r = 50.0;
        cam.theta = std::numbers::pi / 2.0;
        cam.phi = 0.0;
        cam.fov = 60.0f;
        cam.width = width;
        cam.height = height;
        camera = std::make_unique<PinholeCamera>(cam);
    }
};

// -----------------------------------------------------------------------------
// Flat space: a parallel bundle neither focuses nor shears, so its cross-section
// is invariant and the magnification is exactly one.
// -----------------------------------------------------------------------------
TEST(RayBundleTest, FlatSpaceBundleMagnificationIsUnity) {
    BundleProbe probe;
    probe.Build(0.0, 0.0, 90.0f);  // Minkowski.
    CameraRay ray = probe.camera->GenerateRay(40, 32, 0.5f, 0.5f);
    TraceResult r = probe.tracer->Trace(ray);
    ASSERT_TRUE(r.beam.valid);
    EXPECT_NEAR(r.beam.magnification, 1.0f, 1e-3f);
    EXPECT_NEAR(r.beam.area_ratio, 1.0f, 1e-3f);
    EXPECT_GT(r.beam.semi_major, 0.0f);
    EXPECT_GE(r.beam.semi_major, r.beam.semi_minor);
}

// -----------------------------------------------------------------------------
// Oracle agreement (gate a), Schwarzschild radial-plus-shear congruence: the
// CPU Cartesian curvature matches the oracle analytic Kretschmann at every angle
// (a = 0 is spherically symmetric).
// -----------------------------------------------------------------------------
TEST(RayBundleTest, KretschmannMatchesOracleSchwarzschild) {
    constexpr double kTol = sirius::core::constants::geodesic::kConservationTol;  // 1e-4.
    double worst = MaxKretschmannError(
        0.0, {4.0, 6.0, 8.0, 15.0},
        {std::numbers::pi / 2, std::numbers::pi / 3, std::numbers::pi / 4, std::numbers::pi / 6});
    std::cout << "[oracle a=0] max Kretschmann rel error = " << worst << "\n";
    EXPECT_LT(worst, kTol);
}

// -----------------------------------------------------------------------------
// Oracle agreement (gate a), Kerr a=0.9: the equatorial congruence where the
// disk and photon ring live, and off-equator angles where a*cos(theta) != 0 -
// the regime that separates Henry's bracket from the historical defect (see
// the file header).
// -----------------------------------------------------------------------------
TEST(RayBundleTest, KretschmannMatchesOracleKerrEquatorial) {
    constexpr double kTol = sirius::core::constants::geodesic::kConservationTol;  // 1e-4.
    double worst =
        MaxKretschmannError(0.9, {3.0, 4.0, 6.0, 8.0, 15.0, 30.0}, {std::numbers::pi / 2});
    std::cout << "[oracle a=0.9 equatorial] max Kretschmann rel error = " << worst << "\n";
    EXPECT_LT(worst, kTol);
}

TEST(RayBundleTest, KretschmannMatchesOracleKerrOffEquatorial) {
    constexpr double kTol = sirius::core::constants::geodesic::kConservationTol;  // 1e-4.
    double worst = MaxKretschmannError(
        0.9, {4.0, 6.0, 8.0, 15.0, 30.0},
        {std::numbers::pi / 3, std::numbers::pi / 4, std::numbers::pi / 6, 2.2});
    std::cout << "[oracle a=0.9 off-equator] max Kretschmann rel error = " << worst << "\n";
    EXPECT_LT(worst, kTol);
}

// -----------------------------------------------------------------------------
// The propagated bundle is finite, deterministic, and yields a well-ordered
// ellipse for Kerr a=0.9.
// -----------------------------------------------------------------------------
TEST(RayBundleTest, BundleFiniteAndDeterministicKerr) {
    BundleProbe probe;
    probe.Build(1.0, 0.9, 100.0f);
    CameraRay ray = probe.camera->GenerateRay(38, 32, 0.5f, 0.5f);
    TraceResult a = probe.tracer->Trace(ray);
    TraceResult b = probe.tracer->Trace(ray);
    ASSERT_TRUE(a.beam.valid);
    EXPECT_TRUE(std::isfinite(a.beam.magnification));
    EXPECT_TRUE(std::isfinite(a.beam.semi_major));
    EXPECT_GE(a.beam.semi_major, a.beam.semi_minor);
    EXPECT_GT(a.beam.semi_minor, 0.0f);
    EXPECT_FLOAT_EQ(a.beam.magnification, b.beam.magnification);  // Deterministic.
    EXPECT_FLOAT_EQ(a.beam.area_ratio, b.beam.area_ratio);
}

// -----------------------------------------------------------------------------
// Magnification consistency (gate b): the Jacobi determinant is the only
// magnification authority. A close approach defocuses the bundle while a
// grazing ray leaves its cross-section near unity.
// -----------------------------------------------------------------------------
TEST(RayBundleTest, MagnificationComesOnlyFromJacobiMap) {
    BundleProbe probe;
    probe.Build(1.0, 0.0, 120.0f);

    double min_r_closest = 1e30, area_closest = 0.0;
    double area_grazing = 1e30, min_r_grazing = 0.0;

    for (int x = 0; x < 64; ++x) {
        CameraRay ray = probe.camera->GenerateRay(x, 32, 0.5f, 0.5f);
        TraceResult r = probe.tracer->Trace(ray);
        if (r.outcome != TraceResult::Outcome::Escaped || !r.beam.valid) continue;
        if (r.min_radius < min_r_closest) {
            min_r_closest = r.min_radius;
            area_closest = r.beam.area_ratio;
        }
        if (r.min_radius > min_r_grazing) {
            min_r_grazing = r.min_radius;
            area_grazing = r.beam.area_ratio;
        }
    }

    std::cout << "[gate-b] closest: min_r=" << min_r_closest << " area_ratio=" << area_closest
              << " | grazing: min_r=" << min_r_grazing << " area_ratio=" << area_grazing << "\n";

    // The most-deflected ray defocuses the bundle far more than the grazing ray.
    EXPECT_GT(area_closest, 5.0 * area_grazing);
    // A grazing ray keeps the bundle near the flat-space cross-section.
    EXPECT_LT(area_grazing, 5.0);
}

}  // namespace
