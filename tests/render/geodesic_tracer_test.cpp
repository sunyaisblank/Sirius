// GeodesicTracer integration tests. Ported from TSIN005A.cpp.
//
// Verifies that the CPU geodesic tracer traces null geodesics through curved
// spacetime and correctly classifies each ray's termination. Assertions and
// tolerances are identical to the legacy suite.

#include "sirius/backend/cpu/geodesic_tracer.h"

#include "sirius/core/camera.h"
#include "sirius/core/metrics/kerr_schild_family.h"
#include "sirius/core/metrics/morris_thorne_family.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <iostream>
#include <memory>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace {

using namespace sirius::core;
using namespace sirius::backend;

class GeodesicTracerTest : public ::testing::Test {
  protected:
    void SetUp() override {
        // Schwarzschild metric (a = 0).
        KerrSchildParams params;
        params.M = 1.0;
        params.a = 0.0;
        params.Q = 0.0;
        params.Lambda = 0.0;
        m_Metric = std::make_unique<KerrSchildFamily>(params);

        TracerConfig config;
        config.escape_radius = 100.0f;
        config.horizon_factor = 1.05f;
        config.max_steps = 5000;
        config.enable_disk = true;
        config.disk_inner = 6.0f;  // ISCO for Schwarzschild.
        config.disk_outer = 20.0f;
        config.integrator.abs_tolerance = 1e-6f;
        config.integrator.rel_tolerance = 1e-6f;

        m_Tracer = std::make_unique<GeodesicTracer>(m_Metric.get(), config);

        CameraConfig camConfig;
        camConfig.r = 50.0;
        camConfig.theta = M_PI / 2.0;
        camConfig.phi = 0.0;
        camConfig.fov = 60.0f;
        camConfig.width = 64;
        camConfig.height = 64;
        m_Camera = std::make_unique<PinholeCamera>(camConfig);
    }

    std::unique_ptr<KerrSchildFamily> m_Metric;
    std::unique_ptr<GeodesicTracer> m_Tracer;
    std::unique_ptr<PinholeCamera> m_Camera;
};

TEST_F(GeodesicTracerTest, Construction) {
    EXPECT_NE(m_Tracer.get(), nullptr);
    EXPECT_NE(m_Metric.get(), nullptr);
    EXPECT_NE(m_Camera.get(), nullptr);
}

TEST_F(GeodesicTracerTest, CameraRayGeneration) {
    CameraRay ray = m_Camera->GenerateRay(32, 32, 0.5f, 0.5f);

    EXPECT_NEAR(ray.origin(1), 50.0, 0.1);         // r = 50M.
    EXPECT_NEAR(ray.origin(2), M_PI / 2.0, 0.01);  // theta = pi/2.

    double dir_len =
        std::sqrt(ray.direction(1) * ray.direction(1) + ray.direction(2) * ray.direction(2) +
                  ray.direction(3) * ray.direction(3));
    EXPECT_GT(dir_len, 0.0);
}

TEST_F(GeodesicTracerTest, BasicTracing) {
    CameraRay ray = m_Camera->GenerateRay(32, 32, 0.5f, 0.5f);

    TraceResult result = m_Tracer->Trace(ray);

    EXPECT_GT(result.steps_taken, 0);
    EXPECT_FALSE(result.numerical_failure);

    EXPECT_TRUE(result.outcome == TraceResult::Outcome::Escaped ||
                result.outcome == TraceResult::Outcome::Horizon ||
                result.outcome == TraceResult::Outcome::DiskHit ||
                result.outcome == TraceResult::Outcome::MaxSteps);
}

TEST_F(GeodesicTracerTest, HorizonCapture) {
    int horizon_or_disk = 0;
    int total = 0;

    for (int y = 24; y < 40; y += 2) {
        for (int x = 24; x < 40; x += 2) {
            CameraRay ray = m_Camera->GenerateRay(x, y, 0.5f, 0.5f);
            TraceResult result = m_Tracer->Trace(ray);
            total++;

            if (result.outcome == TraceResult::Outcome::Horizon ||
                result.outcome == TraceResult::Outcome::DiskHit) {
                horizon_or_disk++;
            }
        }
    }

    EXPECT_GT(total, 0) << "No rays traced";
    EXPECT_GE(horizon_or_disk + (total - horizon_or_disk), total);
}

TEST_F(GeodesicTracerTest, EscapeToInfinity) {
    CameraRay ray = m_Camera->GenerateRay(0, 0, 0.5f, 0.5f);

    TraceResult result = m_Tracer->Trace(ray);

    EXPECT_TRUE(result.outcome == TraceResult::Outcome::Escaped ||
                result.outcome == TraceResult::Outcome::DiskHit ||
                result.outcome == TraceResult::Outcome::MaxSteps);

    if (result.outcome == TraceResult::Outcome::Escaped) {
        double dir_len = std::sqrt(result.final_direction(1) * result.final_direction(1) +
                                   result.final_direction(2) * result.final_direction(2) +
                                   result.final_direction(3) * result.final_direction(3));
        EXPECT_GT(dir_len, 0.0);
    }
}

TEST_F(GeodesicTracerTest, DiskIntersection) {
    int disk_hits = 0;
    int total_rays = 0;
    bool disk_data_valid = true;

    for (int y = 0; y < 64; y += 4) {
        for (int x = 0; x < 64; x += 4) {
            CameraRay ray = m_Camera->GenerateRay(x, y, 0.5f, 0.5f);
            TraceResult result = m_Tracer->Trace(ray);
            total_rays++;

            if (result.outcome == TraceResult::Outcome::DiskHit) {
                disk_hits++;

                if (result.disk_radius < 6.0f || result.disk_radius > 20.0f) {
                    disk_data_valid = false;
                }
                if (result.disk_temperature <= 0.0f) {
                    disk_data_valid = false;
                }
            }
        }
    }

    if (disk_hits > 0) {
        EXPECT_TRUE(disk_data_valid) << "Invalid disk intersection data";
    }

    EXPECT_GT(total_rays, 0) << "No rays traced";
}

TEST_F(GeodesicTracerTest, DiskTemperatureProfile) {
    TracerConfig config = m_Tracer->GetConfig();

    float T_inner = config.disk_temperature_inner;
    float r_in = config.disk_inner;

    // Novikov-Thorne: T(r) ~ r^(-3/4).
    float T_at_inner = T_inner * std::pow(r_in / r_in, 0.75f);
    EXPECT_NEAR(T_at_inner, T_inner, 0.001f);

    float T_at_2r = T_inner * std::pow(r_in / (2.0f * r_in), 0.75f);
    EXPECT_NEAR(T_at_2r, T_inner * 0.5946f, 0.01f);

    for (float r = r_in; r <= 20.0f; r += 2.0f) {
        float T = T_inner * std::pow(r_in / r, 0.75f);
        EXPECT_GT(T, 0.0f) << "Temperature should be positive at r=" << r;
    }

    EXPECT_GT(config.disk_inner, 0.0f);
    EXPECT_GT(config.disk_outer, config.disk_inner);
}

TEST_F(GeodesicTracerTest, NoNumericalFailures) {
    int failures = 0;
    int total = 0;

    for (int y = 0; y < 64; y += 4) {
        for (int x = 0; x < 64; x += 4) {
            CameraRay ray = m_Camera->GenerateRay(x, y, 0.5f, 0.5f);
            TraceResult result = m_Tracer->Trace(ray);
            total++;

            if (result.numerical_failure) {
                failures++;
            }

            if (result.outcome == TraceResult::Outcome::Escaped) {
                EXPECT_FALSE(std::isnan(result.final_direction(0)));
                EXPECT_FALSE(std::isnan(result.final_direction(1)));
                EXPECT_FALSE(std::isnan(result.final_direction(2)));
                EXPECT_FALSE(std::isnan(result.final_direction(3)));
            }
        }
    }

    double failure_rate = static_cast<double>(failures) / total;
    EXPECT_LT(failure_rate, 0.01) << "Failure rate: " << (failure_rate * 100) << "%";
}

TEST_F(GeodesicTracerTest, KerrMetricTracing) {
    KerrSchildParams kerrParams;
    kerrParams.M = 1.0;
    kerrParams.a = 0.9;
    kerrParams.Q = 0.0;
    kerrParams.Lambda = 0.0;
    auto kerrMetric = std::make_unique<KerrSchildFamily>(kerrParams);

    TracerConfig config;
    config.escape_radius = 100.0f;
    config.horizon_factor = 1.1f;
    config.max_steps = 10000;
    config.enable_disk = true;
    config.disk_inner = 2.32f;  // Approximate prograde ISCO.
    config.disk_outer = 20.0f;

    auto kerrTracer = std::make_unique<GeodesicTracer>(kerrMetric.get(), config);

    for (int x = 24; x < 40; x += 4) {
        CameraRay ray = m_Camera->GenerateRay(x, 32, 0.5f, 0.5f);
        TraceResult result = kerrTracer->Trace(ray);

        EXPECT_FALSE(result.numerical_failure) << "Numerical failure at x=" << x;
        EXPECT_GT(result.steps_taken, 0);
    }
}

TEST_F(GeodesicTracerTest, GFactorDecompositionConsistency) {
    int disk_hits = 0;
    double max_relative_error = 0.0;

    for (int y = 20; y < 44; y += 2) {
        for (int x = 20; x < 44; x += 2) {
            CameraRay ray = m_Camera->GenerateRay(x, y, 0.5f, 0.5f);
            TraceResult result = m_Tracer->Trace(ray);

            if (result.outcome == TraceResult::Outcome::DiskHit) {
                disk_hits++;

                float grav = result.gfactor_grav;
                float gamma = result.gfactor_gamma;
                float v_orb = result.gfactor_v_orb;
                float A = result.gfactor_A;

                float v_dot_n = v_orb * A;
                float doppler_denom = gamma * (1.0f - v_dot_n);
                doppler_denom = std::clamp(doppler_denom, 0.1f, 10.0f);
                float g_reconstructed = grav / doppler_denom;
                g_reconstructed = std::clamp(g_reconstructed, 0.1f, 5.0f);

                float g_stored = result.redshift;
                float rel_error = std::abs(g_reconstructed - g_stored) / g_stored;
                max_relative_error = std::max(max_relative_error, static_cast<double>(rel_error));

                EXPECT_NEAR(g_reconstructed, g_stored, 1e-5f)
                    << "g-factor decomposition mismatch at disk_r=" << result.disk_radius
                    << ", disk_phi=" << result.disk_phi;
            }
        }
    }

    EXPECT_GT(disk_hits, 10) << "Need sufficient disk hits to validate";
    EXPECT_LT(max_relative_error, 1e-4) << "Max relative error in g-factor reconstruction";
}

TEST_F(GeodesicTracerTest, MotionBlurConvergence) {
    TraceResult disk_result;
    bool found = false;
    for (int y = 16; y < 48 && !found; y += 2) {
        for (int x = 16; x < 48 && !found; x += 2) {
            CameraRay ray = m_Camera->GenerateRay(x, y, 0.5f, 0.5f);
            disk_result = m_Tracer->Trace(ray);
            if (disk_result.outcome == TraceResult::Outcome::DiskHit) {
                found = true;
            }
        }
    }

    if (!found) {
        GTEST_SKIP() << "No disk hit found for motion blur test";
        return;
    }

    float grav = disk_result.gfactor_grav;
    float gamma = disk_result.gfactor_gamma;
    float v_orb = disk_result.gfactor_v_orb;
    float A = disk_result.gfactor_A;
    float B = disk_result.gfactor_B;
    float r = disk_result.disk_radius;

    float M = 1.0f;  // Schwarzschild.
    float a = 0.0f;
    float sqrtM = std::sqrt(M);
    float Omega = sqrtM / (std::pow(r, 1.5f) + a * sqrtM);
    float shutter_time = 0.1f;
    float delta_phi_max = Omega * shutter_time;

    auto compute_blur = [&](int N) -> float {
        float g_sum = 0.0f;
        for (int i = 0; i < N; i++) {
            float t = (N > 1) ? static_cast<float>(i) / (N - 1) : 0.5f;
            float delta_phi = delta_phi_max * (t - 0.5f);
            float cos_dphi = std::cos(delta_phi);
            float sin_dphi = std::sin(delta_phi);
            float v_dot_n = v_orb * (A * cos_dphi + B * sin_dphi);
            float denom = gamma * (1.0f - v_dot_n);
            denom = std::clamp(denom, 0.1f, 10.0f);
            float g_offset = std::clamp(grav / denom, 0.1f, 5.0f);
            g_sum += g_offset;
        }
        return g_sum / N;
    };

    float g_N4 = compute_blur(4);
    float g_N8 = compute_blur(8);
    float g_N16 = compute_blur(16);
    float g_N32 = compute_blur(32);
    float g_N64 = compute_blur(64);

    float err_4_to_64 = std::abs(g_N4 - g_N64);
    float err_8_to_64 = std::abs(g_N8 - g_N64);
    float err_16_to_64 = std::abs(g_N16 - g_N64);
    float err_32_to_64 = std::abs(g_N32 - g_N64);
    (void)err_32_to_64;

    EXPECT_LE(err_8_to_64, err_4_to_64 + 1e-6f) << "Error should not increase with more samples";
    EXPECT_LE(err_16_to_64, err_8_to_64 + 1e-6f) << "Error should not increase with more samples";

    float relative_err_8 = (g_N64 > 1e-6f) ? (err_8_to_64 / g_N64) : 0.0f;
    EXPECT_LT(relative_err_8, 0.01f) << "N=8 should be within 1% of converged value";

    EXPECT_GT(g_N4, 0.0f);
    EXPECT_GT(g_N64, 0.0f);
    EXPECT_LT(g_N64, 10.0f);
}

TEST_F(GeodesicTracerTest, GFactorCoefficientNormalization) {
    int disk_hits = 0;

    for (int y = 20; y < 44; y += 3) {
        for (int x = 20; x < 44; x += 3) {
            CameraRay ray = m_Camera->GenerateRay(x, y, 0.5f, 0.5f);
            TraceResult result = m_Tracer->Trace(ray);

            if (result.outcome == TraceResult::Outcome::DiskHit) {
                disk_hits++;

                float A = result.gfactor_A;
                float B = result.gfactor_B;
                float norm_squared = A * A + B * B;

                EXPECT_LE(norm_squared, 1.01f) << "A^2 + B^2 should be <= 1, got " << norm_squared;
                EXPECT_GT(norm_squared, 0.0f)
                    << "A^2 + B^2 should be > 0 for rays hitting equatorial disk";
            }
        }
    }

    EXPECT_GT(disk_hits, 5) << "Need disk hits to validate coefficients";
}

TEST_F(GeodesicTracerTest, TracingPerformance) {
    TracerConfig perfConfig;
    perfConfig.escape_radius = 100.0f;
    perfConfig.horizon_factor = 1.05f;
    perfConfig.max_steps = 500;
    perfConfig.enable_disk = false;
    perfConfig.integrator.abs_tolerance = 1e-3f;
    perfConfig.integrator.rel_tolerance = 1e-3f;
    perfConfig.integrator.initial_step = 0.1f;

    auto perfTracer = std::make_unique<GeodesicTracer>(m_Metric.get(), perfConfig);

    auto start = std::chrono::high_resolution_clock::now();

    int total_rays = 0;
    for (int y = 0; y < 64; y += 4) {
        for (int x = 0; x < 64; x += 4) {
            CameraRay ray = m_Camera->GenerateRay(x, y, 0.5f, 0.5f);
            TraceResult result = perfTracer->Trace(ray);
            (void)result;
            total_rays++;
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    EXPECT_LT(duration.count(), 60000)
        << "Traced " << total_rays << " rays in " << duration.count() << "ms";

    double rays_per_sec = 1000.0 * total_rays / duration.count();
    std::cout << "CPU tracing: " << total_rays << " rays in " << duration.count() << "ms ("
              << rays_per_sec << " rays/sec)" << std::endl;
}

// =============================================================================
// Morris-Thorne (Ellis) wormhole through the same Cartesian tracer, via the
// MorrisThorneCartesian embedding. One sheet; the throat is the capture
// surface, so a central ray terminates as Horizon (the tracer's capture
// outcome) and offset rays escape with gravitational deflection.
// =============================================================================

class MorrisThorneTracerTest : public ::testing::Test {
  protected:
    void SetUp() override {
        m_Metric = std::make_unique<MorrisThorneCartesian>(MorrisThorneParams::Ellis(1.0));

        TracerConfig config;
        config.escape_radius = 200.0f;
        config.horizon_factor = 1.05f;
        config.max_steps = 20000;
        config.enable_disk = false;  // The thin disk is a black-hole construct.
        config.integrator.abs_tolerance = 1e-7f;
        config.integrator.rel_tolerance = 1e-7f;
        config.integrator.initial_step = 0.05f;
        config.integrator.max_step = 1.0f;

        m_Tracer = std::make_unique<GeodesicTracer>(m_Metric.get(), config);

        CameraConfig camConfig;
        camConfig.r = 50.0;
        camConfig.theta = M_PI / 2.0;
        camConfig.phi = 0.0;
        camConfig.fov = 60.0f;
        camConfig.width = 64;
        camConfig.height = 64;
        m_Camera = std::make_unique<PinholeCamera>(camConfig);
    }

    // Deflection angle between the launch direction and the escape direction,
    // both in Cartesian components. At theta = pi/2, phi = 0 the camera's
    // spherical direction slots map to Cartesian as
    //   v_x = dir(1),  v_y = dir(3),  v_z = -dir(2)
    // (the Jacobian in GeodesicTracer::InitializeLightray specialised to the
    // probe's observer placement).
    static double Deflection(const CameraRay& ray, const TraceResult& result) {
        double ix = ray.direction(1);
        double iy = ray.direction(3);
        double iz = -ray.direction(2);
        double il = std::sqrt(ix * ix + iy * iy + iz * iz);

        double fx = result.final_direction(1);
        double fy = result.final_direction(2);
        double fz = result.final_direction(3);
        double fl = std::sqrt(fx * fx + fy * fy + fz * fz);

        double c = (ix * fx + iy * fy + iz * fz) / (il * fl);
        return std::acos(std::clamp(c, -1.0, 1.0));
    }

    std::unique_ptr<MorrisThorneCartesian> m_Metric;
    std::unique_ptr<GeodesicTracer> m_Tracer;
    std::unique_ptr<PinholeCamera> m_Camera;
};

TEST_F(MorrisThorneTracerTest, CentralRayCapturedAtThroat) {
    CameraRay ray = m_Camera->GenerateRay(32, 32, 0.5f, 0.5f);
    TraceResult result = m_Tracer->Trace(ray);
    EXPECT_EQ(result.outcome, TraceResult::Outcome::Horizon)
        << "A radially ingoing ray must terminate at the throat capture surface"
        << " (outcome=" << static_cast<int>(result.outcome) << ", min_radius=" << result.min_radius
        << ", steps=" << result.steps_taken << ", numerical_failure=" << result.numerical_failure
        << ")";
}

TEST_F(MorrisThorneTracerTest, EdgeRayEscapes) {
    CameraRay ray = m_Camera->GenerateRay(2, 32, 0.5f, 0.5f);
    TraceResult result = m_Tracer->Trace(ray);
    EXPECT_EQ(result.outcome, TraceResult::Outcome::Escaped);
}

// Ellis lensing: the leading-order deflection is alpha ~ (pi/4)(b0/rho)^2
// (Chetouani & Clement 1984), so doubling the impact parameter divides the
// deflection by about four. The gate checks sign, monotone decrease, and the
// quadratic ratio within generous bounds that still exclude both zero
// deflection and the 1/rho scaling a mass-like defect would produce.
TEST_F(MorrisThorneTracerTest, DeflectionFallsQuadraticallyWithImpactParameter) {
    // Pixel offsets chosen so the impact parameters are roughly 5, 10, 20 b0.
    struct Probe {
        int px;
        double deflection = 0.0;
        double rho = 0.0;
    };
    std::array<Probe, 3> probes = {{{38, 0, 0}, {44, 0, 0}, {57, 0, 0}}};

    for (auto& p : probes) {
        CameraRay ray = m_Camera->GenerateRay(p.px, 32, 0.5f, 0.5f);
        TraceResult result = m_Tracer->Trace(ray);
        ASSERT_EQ(result.outcome, TraceResult::Outcome::Escaped)
            << "probe pixel " << p.px << " must escape";
        p.deflection = Deflection(ray, result);

        // Impact parameter |x0 x v_hat| from the actual launch geometry.
        double ix = ray.direction(1), iy = ray.direction(3), iz = -ray.direction(2);
        double il = std::sqrt(ix * ix + iy * iy + iz * iz);
        ix /= il;
        iy /= il;
        iz /= il;
        double x0 = 50.0, y0 = 0.0, z0 = 0.0;
        double cx = y0 * iz - z0 * iy;
        double cy = z0 * ix - x0 * iz;
        double cz = x0 * iy - y0 * ix;
        p.rho = std::sqrt(cx * cx + cy * cy + cz * cz);
    }

    EXPECT_GT(probes[0].deflection, probes[1].deflection);
    EXPECT_GT(probes[1].deflection, probes[2].deflection);
    EXPECT_GT(probes[2].deflection, 0.0);

    for (const auto& p : probes) {
        double predicted = (M_PI / 4.0) * (1.0 / (p.rho * p.rho));  // b0 = 1.
        EXPECT_NEAR(p.deflection, predicted, 0.5 * predicted)
            << "Ellis leading-order deflection at rho=" << p.rho << ": traced=" << p.deflection
            << ", predicted=" << predicted;
    }
}

}  // namespace
