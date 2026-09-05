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
#include <numbers>

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
        camConfig.theta = std::numbers::pi / 2.0;
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

    EXPECT_NEAR(ray.origin(1), 50.0, 0.1);                     // r = 50M.
    EXPECT_NEAR(ray.origin(2), std::numbers::pi / 2.0, 0.01);  // theta = pi/2.

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
                result.outcome == TraceResult::Outcome::Throat ||
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
                result.outcome == TraceResult::Outcome::Throat ||
                result.outcome == TraceResult::Outcome::DiskHit) {
                horizon_or_disk++;
            }
        }
    }

    EXPECT_GT(total, 0) << "No rays traced";
    EXPECT_GT(horizon_or_disk, 0) << "the central image region contains no captured or disk rays";
}

TEST_F(GeodesicTracerTest, EscapeToInfinity) {
    CameraRay ray = m_Camera->GenerateRay(0, 0, 0.5f, 0.5f);

    TraceResult result = m_Tracer->Trace(ray);

    ASSERT_EQ(result.outcome, TraceResult::Outcome::Escaped);
    double dir_len = std::sqrt(result.final_direction(1) * result.final_direction(1) +
                               result.final_direction(2) * result.final_direction(2) +
                               result.final_direction(3) * result.final_direction(3));
    EXPECT_GT(dir_len, 0.0);
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
                if (std::abs(result.final_position(3)) > 1.0e-6f) {
                    disk_data_valid = false;
                }
            }
        }
    }

    EXPECT_GT(total_rays, 0) << "No rays traced";
    ASSERT_GT(disk_hits, 0) << "deterministic disk scene produced no disk intersection";
    EXPECT_TRUE(disk_data_valid) << "Invalid disk intersection data";
}

TEST_F(GeodesicTracerTest, LiveDiskTemperatureUsesZeroTorqueShakuraSunyaevProfile) {
    TracerConfig config = m_Tracer->GetConfig();
    config.disk_temperature_model = DiskTemperatureModel::ShakuraSunyaev;
    config.disk_temperature_inner = 10000.0f;
    GeodesicTracer tracer(m_Metric.get(), config);

    int compared = 0;
    for (int y = 8; y < 56 && compared < 12; y += 2) {
        for (int x = 8; x < 56 && compared < 12; x += 2) {
            const TraceResult result = tracer.Trace(m_Camera->GenerateRay(x, y, 0.5f, 0.5f));
            for (int i = 0; i < result.num_disk_crossings; ++i) {
                const auto& crossing = result.disk_crossings[i];
                if (!crossing.valid) continue;
                const double expected = ShakuraSunyaevTemperature(config.disk_temperature_inner,
                                                                  crossing.r, config.disk_inner);
                EXPECT_NEAR(crossing.temperature, expected, 2.0e-3)
                    << "live crossing at r=" << crossing.r
                    << " did not consume the zero-torque Shakura-Sunyaev profile";
                ++compared;
            }
        }
    }
    EXPECT_GE(compared, 6)
        << "insufficient live disk crossings for the Shakura-Sunyaev profile gate";
}

TEST_F(GeodesicTracerTest, LiveDiskTemperatureUsesFullPageThorneProfile) {
    AccretionDiskD::Config disk_config;
    disk_config.M = 1.0;
    disk_config.a_star = 0.0;
    AccretionDiskD oracle(disk_config);
    const double reference = oracle.Temperature(1.5 * oracle.IscoRadius());
    ASSERT_GT(reference, 0.0);

    int compared = 0;
    for (int y = 8; y < 56 && compared < 12; y += 2) {
        for (int x = 8; x < 56 && compared < 12; x += 2) {
            const TraceResult result = m_Tracer->Trace(m_Camera->GenerateRay(x, y, 0.5f, 0.5f));
            for (int i = 0; i < result.num_disk_crossings; ++i) {
                const auto& crossing = result.disk_crossings[i];
                if (!crossing.valid) continue;
                const double expected = oracle.Temperature(crossing.r) / reference;
                EXPECT_NEAR(crossing.temperature, expected, 2.0e-5)
                    << "live crossing at r=" << crossing.r
                    << " did not consume the full Page-Thorne profile";
                ++compared;
            }
        }
    }
    EXPECT_GE(compared, 6) << "insufficient live disk crossings for the Page-Thorne gate";
}

TEST_F(GeodesicTracerTest, LiveDiskCrossingCarriesTransportedPhysicalStokesOrientation) {
    CameraRay disk_ray;
    bool found = false;
    for (int y = 8; y < 56 && !found; y += 2) {
        for (int x = 8; x < 56 && !found; x += 2) {
            const CameraRay candidate = m_Camera->GenerateRay(x, y, 0.5f, 0.5f);
            if (m_Tracer->Trace(candidate).outcome == TraceResult::Outcome::DiskHit) {
                disk_ray = candidate;
                found = true;
            }
        }
    }
    ASSERT_TRUE(found) << "no thin-disk ray available for the live polarisation witness";

    TracerConfig config = m_Tracer->GetConfig();
    config.enable_polarisation = true;
    GeodesicTracer polarised_tracer(m_Metric.get(), config);
    const TraceResult first = polarised_tracer.Trace(disk_ray);
    const TraceResult repeat = polarised_tracer.Trace(disk_ray);

    ASSERT_EQ(first.outcome, TraceResult::Outcome::DiskHit);
    ASSERT_GT(first.num_disk_crossings, 0);
    ASSERT_EQ(first.num_disk_crossings, repeat.num_disk_crossings);
    for (int i = 0; i < first.num_disk_crossings; ++i) {
        const auto& crossing = first.disk_crossings[i];
        const auto& repeated = repeat.disk_crossings[i];
        ASSERT_TRUE(crossing.polarisation_valid);
        EXPECT_TRUE(std::isfinite(crossing.polarisation_evpa));
        EXPECT_GE(crossing.polarisation_degree, 0.0f);
        EXPECT_LE(crossing.polarisation_degree, 0.1171f);
        EXPECT_GT(crossing.polarisation_intensity_scale, 0.0f);
        EXPECT_FLOAT_EQ(crossing.polarisation_evpa, repeated.polarisation_evpa);
        EXPECT_FLOAT_EQ(crossing.polarisation_degree, repeated.polarisation_degree);
        EXPECT_FLOAT_EQ(crossing.polarisation_intensity_scale,
                        repeated.polarisation_intensity_scale);
    }

    const TraceResult unpolarised = m_Tracer->Trace(disk_ray);
    EXPECT_FALSE(unpolarised.disk_crossings[0].polarisation_valid);
}

TEST(GeodesicTracerVolumetric, TransferAccumulatesAcrossEveryTraversedSegment) {
    KerrSchildParams params;
    params.M = 1.0;
    KerrSchildFamily metric(params);

    TracerConfig config;
    config.escape_radius = 80.0f;
    config.max_steps = 20000;
    config.horizon_factor = 1.05f;
    config.enable_disk = true;
    config.enable_volumetric = true;
    config.disk_inner = 6.0f;
    config.disk_outer = 20.0f;
    config.disk_temperature_inner = 1.0f;
    config.disk_temperature_model = DiskTemperatureModel::ShakuraSunyaev;
    config.volumetric_scale_height_ratio = 0.1f;
    config.volumetric_tau_midplane = 2.0f;
    config.volumetric_samples = 4;
    config.integrator.initial_step = 0.5f;
    config.integrator.max_step = 1.0f;
    GeodesicTracer tracer(&metric, config);

    CameraConfig camera_config;
    camera_config.r = 50.0;
    camera_config.theta = std::numbers::pi / 2.0;
    camera_config.width = 3;
    camera_config.height = 3;
    PinholeCamera camera(camera_config);

    const TraceResult result = tracer.Trace(camera.GenerateRay(1, 1, 0.5f, 0.5f));
    EXPECT_EQ(result.outcome, TraceResult::Outcome::Horizon)
        << "volumetric transfer must not replace the ray's terminal outcome";
    EXPECT_TRUE(result.volumetric_hit);
    EXPECT_GT(result.volumetric_affine_length, 10.0f)
        << "only a final RK45 segment appears to have been retained";
    EXPECT_GT(result.optical_depth, 0.1f);
    EXPECT_GT(result.volumetric_emission[0], 0.0f);
    EXPECT_TRUE(std::isfinite(result.volumetric_emission[0]));

    TracerConfig capped_config = config;
    capped_config.volumetric_tau_max = 0.25f;
    GeodesicTracer capped_tracer(&metric, capped_config);
    const TraceResult capped = capped_tracer.Trace(camera.GenerateRay(1, 1, 0.5f, 0.5f));
    EXPECT_TRUE(capped.volumetric_hit);
    EXPECT_NEAR(capped.optical_depth, capped_config.volumetric_tau_max, 1.0e-6f);
    EXPECT_LE(capped.optical_depth, capped_config.volumetric_tau_max + 1.0e-6f)
        << "the final transfer sample overshot the configured optical-depth ceiling";
}

TEST(GeodesicTracerVolumetric, OpticallyThinTransferIsNotDiscardedAtCompositionBoundary) {
    KerrSchildParams params;
    params.M = 1.0;
    KerrSchildFamily metric(params);

    TracerConfig config;
    config.escape_radius = 80.0f;
    config.max_steps = 2000;
    config.enable_disk = true;
    config.enable_volumetric = true;
    config.disk_inner = 6.0f;
    config.disk_outer = 20.0f;
    config.disk_temperature_inner = 1.0f;
    config.disk_temperature_model = DiskTemperatureModel::ShakuraSunyaev;
    config.volumetric_scale_height_ratio = 0.1f;
    config.volumetric_tau_midplane = 1.0e-5f;
    config.volumetric_samples = 4;
    config.integrator.initial_step = 0.5f;
    config.integrator.max_step = 1.0f;
    GeodesicTracer tracer(&metric, config);

    CameraConfig camera_config;
    camera_config.r = 50.0;
    camera_config.theta = std::numbers::pi / 2.0;
    camera_config.width = 3;
    camera_config.height = 3;
    PinholeCamera camera(camera_config);

    const TraceResult result = tracer.Trace(camera.GenerateRay(1, 1, 0.5f, 0.5f));
    ASSERT_GT(result.optical_depth, 0.0f);
    ASSERT_LT(result.optical_depth, 0.01f)
        << "the fixture must remain below the former display-oriented cutoff";
    EXPECT_TRUE(result.volumetric_hit);
    EXPECT_GT(result.volumetric_emission[0], 0.0f);
}

TEST(GeodesicTracerVolumetric, RedshiftAndDopplerReachTheLiveVolumeSource) {
    KerrSchildParams params;
    params.M = 1.0;
    params.a = 0.6;
    KerrSchildFamily metric(params);

    TracerConfig config;
    config.escape_radius = 80.0f;
    config.max_steps = 4000;
    config.enable_disk = true;
    config.enable_volumetric = true;
    config.disk_inner = 4.0f;
    config.disk_outer = 20.0f;
    config.disk_temperature_inner = 1.0f;
    config.disk_temperature_model = DiskTemperatureModel::ShakuraSunyaev;
    config.volumetric_scale_height_ratio = 0.15f;
    config.volumetric_tau_midplane = 1.0f;
    config.volumetric_samples = 2;
    config.integrator.initial_step = 0.25f;
    config.integrator.max_step = 1.0f;

    CameraConfig camera_config;
    camera_config.r = 40.0;
    camera_config.theta = std::numbers::pi / 2.0;
    camera_config.width = 3;
    camera_config.height = 3;
    PinholeCamera camera(camera_config);
    const CameraRay ray = camera.GenerateRay(1, 1, 0.5f, 0.5f);

    config.doppler_beaming = true;
    GeodesicTracer beamed_tracer(&metric, config);
    const TraceResult beamed = beamed_tracer.Trace(ray);
    config.doppler_beaming = false;
    GeodesicTracer gravity_only_tracer(&metric, config);
    const TraceResult gravity_only = gravity_only_tracer.Trace(ray);

    ASSERT_TRUE(beamed.volumetric_hit);
    ASSERT_TRUE(gravity_only.volumetric_hit);
    EXPECT_GT(beamed.volumetric_emission[0], 0.0f);
    EXPECT_GT(gravity_only.volumetric_emission[0], 0.0f);
    EXPECT_LT(beamed.volumetric_emission[0], gravity_only.volumetric_emission[0])
        << "the orbital Doppler term did not reach the live volumetric source";
}

TEST(GeodesicTracerVolumetric, ProceduralTurbulenceAltersLiveTransferDeterministically) {
    KerrSchildParams params;
    params.M = 1.0;
    KerrSchildFamily metric(params);

    TracerConfig baseline_config;
    baseline_config.escape_radius = 80.0f;
    baseline_config.max_steps = 2000;
    baseline_config.enable_disk = true;
    baseline_config.enable_volumetric = true;
    baseline_config.disk_inner = 6.0f;
    baseline_config.disk_outer = 20.0f;
    baseline_config.disk_temperature_inner = 1.0f;
    baseline_config.disk_temperature_model = DiskTemperatureModel::ShakuraSunyaev;
    baseline_config.volumetric_scale_height_ratio = 0.1f;
    baseline_config.volumetric_tau_midplane = 2.0f;
    baseline_config.volumetric_samples = 4;
    baseline_config.integrator.initial_step = 0.5f;
    baseline_config.integrator.max_step = 1.0f;

    CameraConfig camera_config;
    camera_config.r = 50.0;
    camera_config.theta = std::numbers::pi / 2.0;
    camera_config.width = 3;
    camera_config.height = 3;
    PinholeCamera camera(camera_config);
    const CameraRay ray = camera.GenerateRay(1, 1, 0.5f, 0.5f);

    GeodesicTracer baseline_tracer(&metric, baseline_config);
    const TraceResult baseline = baseline_tracer.Trace(ray);

    TracerConfig enhanced_config = baseline_config;
    enhanced_config.enable_turbulence = true;
    GeodesicTracer enhanced_tracer(&metric, enhanced_config);

    const TraceResult first = enhanced_tracer.Trace(ray);
    const TraceResult repeat = enhanced_tracer.Trace(ray);
    ASSERT_TRUE(first.volumetric_hit);
    EXPECT_TRUE(std::isfinite(first.optical_depth));
    EXPECT_TRUE(std::isfinite(first.volumetric_emission[0]));
    EXPECT_NE(first.optical_depth, baseline.optical_depth);
    EXPECT_NE(first.volumetric_emission[0], baseline.volumetric_emission[0]);
    EXPECT_FLOAT_EQ(first.optical_depth, repeat.optical_depth);
    for (int channel = 0; channel < 3; ++channel) {
        EXPECT_FLOAT_EQ(first.volumetric_emission[channel], repeat.volumetric_emission[channel]);
    }
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
    config.disk_inner = AccretionDiskD::ComputeIsco(0.9);
    config.disk_outer = 20.0f;

    auto kerrTracer = std::make_unique<GeodesicTracer>(kerrMetric.get(), config);

    for (int x = 24; x < 40; x += 4) {
        CameraRay ray = m_Camera->GenerateRay(x, 32, 0.5f, 0.5f);
        TraceResult result = kerrTracer->Trace(ray);

        EXPECT_FALSE(result.numerical_failure) << "Numerical failure at x=" << x;
        EXPECT_GT(result.steps_taken, 0);
    }
}

TEST_F(GeodesicTracerTest, InvariantFrequencyTransferSelectsTheCircularEmitterBranch) {
    int disk_hits = 0;

    for (int y = 20; y < 44; y += 2) {
        for (int x = 20; x < 44; x += 2) {
            CameraRay ray = m_Camera->GenerateRay(x, y, 0.5f, 0.5f);
            TraceResult result = m_Tracer->Trace(ray);

            if (result.outcome == TraceResult::Outcome::DiskHit) {
                disk_hits++;
                EXPECT_TRUE(std::isfinite(result.full_disk_redshift));
                EXPECT_TRUE(std::isfinite(result.zamo_disk_redshift));
                EXPECT_GT(result.full_disk_redshift, 0.0f);
                EXPECT_GT(result.zamo_disk_redshift, 0.0f);
                EXPECT_FLOAT_EQ(result.redshift, result.full_disk_redshift);
                ASSERT_GT(result.num_disk_crossings, 0);
                EXPECT_FLOAT_EQ(result.disk_crossings[0].full_redshift, result.full_disk_redshift);
                EXPECT_FLOAT_EQ(result.disk_crossings[0].zamo_redshift, result.zamo_disk_redshift);
            }
        }
    }

    EXPECT_GT(disk_hits, 10) << "Need sufficient disk hits to validate";
}

TEST(GeodesicTracerRedshift, NearExtremalInnerDiskEmissionRemainsFinite) {
    KerrSchildParams params;
    params.M = 1.0;
    params.a = 0.998;
    KerrSchildFamily metric(params);

    TracerConfig config;
    config.escape_radius = 100.0f;
    config.horizon_factor = 1.0f;
    config.max_steps = 20000;
    config.enable_disk = true;
    config.disk_inner = AccretionDiskD::ComputeIsco(0.998);
    config.disk_outer = 20.0f;
    config.integrator.initial_step = 0.1f;
    config.integrator.max_step = 1.0f;
    GeodesicTracer tracer(&metric, config);

    CameraConfig camera_config;
    camera_config.r = 30.0;
    camera_config.theta = 80.0 * std::numbers::pi / 180.0;
    camera_config.fov = 55.0f;
    camera_config.width = 96;
    camera_config.height = 54;
    PinholeCamera camera(camera_config);

    int inner_crossings = 0;
    for (int y = 4; y < camera_config.height - 4 && inner_crossings < 3; y += 2) {
        for (int x = 4; x < camera_config.width - 4 && inner_crossings < 3; x += 2) {
            const TraceResult result = tracer.Trace(camera.GenerateRay(x, y, 0.5f, 0.5f));
            for (int i = 0; i < result.num_disk_crossings; ++i) {
                const auto& crossing = result.disk_crossings[i];
                if (!crossing.valid || crossing.r >= 2.0f) continue;
                EXPECT_TRUE(std::isfinite(crossing.redshift));
                EXPECT_GT(crossing.redshift, 0.0f);
                ++inner_crossings;
            }
        }
    }
    EXPECT_GE(inner_crossings, 1)
        << "the near-extremal probe did not exercise disk emission inside r=2M";
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
// exact MorrisThorneCartesian isotropic chart. This fixture selects the
// explicit dark one-sheet output boundary; a central ray terminates as Throat
// while offset rays escape with gravitational deflection.
// =============================================================================

class MorrisThorneTracerTest : public ::testing::Test {
  protected:
    void SetUp() override {
        m_Metric = std::make_unique<MorrisThorneCartesian>(MorrisThorneParams::Ellis(1.0));

        TracerConfig config;
        config.escape_radius = 200.0f;
        config.horizon_factor = 1.0f;
        config.max_steps = 20000;
        config.enable_disk = false;  // The thin disk is a black-hole construct.
        config.integrator.abs_tolerance = 1e-7f;
        config.integrator.rel_tolerance = 1e-7f;
        config.integrator.initial_step = 0.05f;
        config.integrator.max_step = 1.0f;

        m_Tracer = std::make_unique<GeodesicTracer>(m_Metric.get(), config);

        CameraConfig camConfig;
        camConfig.r = 50.0;
        camConfig.theta = std::numbers::pi / 2.0;
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

TEST_F(MorrisThorneTracerTest, CentralRayTerminatesAtExplicitThroatBoundary) {
    CameraRay ray = m_Camera->GenerateRay(32, 32, 0.5f, 0.5f);
    TraceResult result = m_Tracer->Trace(ray);
    EXPECT_EQ(result.outcome, TraceResult::Outcome::Throat)
        << "A radially ingoing one-sheet ray must terminate at the explicit throat boundary"
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
        double predicted = (std::numbers::pi / 4.0) * (1.0 / (p.rho * p.rho));  // b0 = 1.
        EXPECT_NEAR(p.deflection, predicted, 0.5 * predicted)
            << "Ellis leading-order deflection at rho=" << p.rho << ": traced=" << p.deflection
            << ", predicted=" << predicted;
    }
}

}  // namespace
