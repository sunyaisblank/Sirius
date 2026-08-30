// CPU Geodesic Reference Tests
// Reference-ray construction from Boyer-Lindquist initial conditions, CPU RK45
// integration to horizon capture or escape in Schwarzschild and Kerr, and the
// impact-parameter grid the parity gate specifies. Ported from TSIN006A.cpp;
// the GPU/OptiX comparison arm is retired (OptiX removed this engagement), so
// the surviving content is the CPU baseline plus the tolerance specification
// anchors. Assertions and tolerances unchanged; includes are core-only.

#include "sirius/core/geodesic_integrator.h"
#include "sirius/core/metrics/kerr_schild_family.h"
#include "sirius/core/metrics/metric.h"
#include "sirius/core/tensor.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <numbers>
#include <vector>

namespace sirius::test {
using namespace sirius::core;

// Parity tolerances (GPU is single precision, the CPU reference is double;
// differences accumulate over integration steps).
constexpr double POSITION_RELATIVE_TOLERANCE = 1e-4;  // |Δr/r| < 10^-4.
constexpr double DIRECTION_TOLERANCE_RAD = 1e-4;      // |Δθ| < 10^-4 rad.
constexpr int NUM_REFERENCE_RAYS = 100;               // Grid of impact parameters.

// Reference ray in Boyer-Lindquist spherical initial conditions.
struct ReferenceRay {
    double r_start;
    double theta_start;
    double phi_start;

    double v_r;
    double v_theta;
    double v_phi;

    enum class Outcome { Escaped, Horizon, Unknown };
    Outcome expected;
};

class CpuGeodesicReferenceTests : public ::testing::Test {
  protected:
    IntegratorConfig cpu_config;
    std::vector<ReferenceRay> referenceRays;

    void SetUp() override {
        cpu_config = Geodesic::GetDefaultConfig();
        cpu_config.abs_tolerance = 1e-8f;
        cpu_config.rel_tolerance = 1e-8f;
        cpu_config.min_step = 1e-8f;
        cpu_config.max_step = 0.1f;
        cpu_config.initial_step = 0.01f;

        generateReferenceRays();
    }

    void generateReferenceRays() {
        referenceRays.clear();

        double r_obs = 50.0;
        double theta_obs = std::numbers::pi / 2.0;  // Equatorial plane.
        double phi_obs = 0.0;

        int n_b = 10;    // Impact parameter samples.
        int n_phi = 10;  // Azimuthal samples.

        for (int ib = 0; ib < n_b; ++ib) {
            double b = 1.0 + 14.0 * ib / (n_b - 1);

            for (int ip = 0; ip < n_phi; ++ip) {
                double phi_dir = 2.0 * std::numbers::pi * ip / n_phi;

                ReferenceRay ray;
                ray.r_start = r_obs;
                ray.theta_start = theta_obs;
                ray.phi_start = phi_obs;

                ray.v_r = -1.0;
                ray.v_theta = 0.0;
                ray.v_phi = b / r_obs * std::cos(phi_dir);

                if (b < 5.0) {
                    ray.expected = ReferenceRay::Outcome::Horizon;
                } else {
                    ray.expected = ReferenceRay::Outcome::Escaped;
                }

                referenceRays.push_back(ray);
            }
        }
    }

    // Lightray from a reference configuration, normalised to the null cone.
    Lightray createLightray(const ReferenceRay& ref, IMetric* metric) {
        Lightray ray{};

        double sin_th = std::sin(ref.theta_start);
        double cos_th = std::cos(ref.theta_start);
        double sin_ph = std::sin(ref.phi_start);
        double cos_ph = std::cos(ref.phi_start);

        ray.position(0) = 0.0f;                           // t.
        ray.position(1) = ref.r_start * sin_th * cos_ph;  // x.
        ray.position(2) = ref.r_start * sin_th * sin_ph;  // y.
        ray.position(3) = ref.r_start * cos_th;           // z.

        ray.proper_time = 0.0f;
        ray.coordinate_time = 0.0f;
        ray.terminated = 0;
        ray.step_size = cpu_config.initial_step;

        double r = ref.r_start;
        ray.velocity(1) = ref.v_r * sin_th * cos_ph + r * ref.v_theta * cos_th * cos_ph -
                          r * ref.v_phi * sin_th * sin_ph;
        ray.velocity(2) = ref.v_r * sin_th * sin_ph + r * ref.v_theta * cos_th * sin_ph +
                          r * ref.v_phi * sin_th * cos_ph;
        ray.velocity(3) = ref.v_r * cos_th - r * ref.v_theta * sin_th;

        Metric4d g;
        Tensor<Dual<double>, 4, 4, 4> dg;
        metric->Evaluate(ray.position, g, dg);
        ray.velocity = TensorOps::NormalizeNull(ray.velocity, g);

        ray.acceleration = Geodesic::CalculateAcceleration(ray.velocity, ray.position, metric);

        return ray;
    }

    struct IntegrationResult {
        Vec4 final_position;
        Vec4 final_velocity;
        int steps;
        bool hit_horizon;
        bool escaped;
    };

    IntegrationResult integrateCPU(const ReferenceRay& ref, IMetric* metric, int max_steps = 1000) {
        IntegrationResult result = {};

        Lightray ray = createLightray(ref, metric);

        for (int step = 0; step < max_steps; ++step) {
            bool success = Geodesic::IntegrateStepRk45(ray, metric, cpu_config);

            if (!success || ray.terminated) {
                result.hit_horizon = (ray.terminated == 1);
                break;
            }

            double r =
                std::sqrt(ray.position(1) * ray.position(1) + ray.position(2) * ray.position(2) +
                          ray.position(3) * ray.position(3));
            if (r > 200.0) {
                result.escaped = true;
                break;
            }

            if (r < 2.1) {
                result.hit_horizon = true;
                break;
            }

            result.steps++;
        }

        result.final_position = ray.position;
        result.final_velocity = ray.velocity;

        return result;
    }
};

// --- CPU integration baselines ----------------------------------------------

TEST_F(CpuGeodesicReferenceTests, CPUBaselineSchwarzschildEscaping) {
    KerrSchildFamily metric{KerrSchildParams::Schwarzschild(1.0)};
    metric.SetParameter("mass", 1.0);

    ReferenceRay ref;
    ref.r_start = 50.0;
    ref.theta_start = std::numbers::pi / 2.0;
    ref.phi_start = 0.0;
    ref.v_r = -1.0;
    ref.v_theta = 0.0;
    ref.v_phi = 0.2;  // b ~ 10M.
    ref.expected = ReferenceRay::Outcome::Escaped;

    IntegrationResult result = integrateCPU(ref, &metric);

    EXPECT_TRUE(result.escaped || result.steps > 500)
        << "Ray with large impact parameter should escape or integrate many steps";
}

TEST_F(CpuGeodesicReferenceTests, CPUBaselineSchwarzschildHorizon) {
    KerrSchildFamily metric{KerrSchildParams::Schwarzschild(1.0)};
    metric.SetParameter("mass", 1.0);

    ReferenceRay ref;
    ref.r_start = 10.0;
    ref.theta_start = std::numbers::pi / 2.0;
    ref.phi_start = 0.0;
    ref.v_r = -1.0;
    ref.v_theta = 0.0;
    ref.v_phi = 0.02;  // b ~ 0.2M (well below b_c ~ 5.2M).
    ref.expected = ReferenceRay::Outcome::Horizon;

    IntegrationResult result = integrateCPU(ref, &metric, 5000);

    EXPECT_TRUE(result.hit_horizon) << "Ray with small impact parameter should hit horizon";
}

TEST_F(CpuGeodesicReferenceTests, CPUBaselineKerrPrograde) {
    KerrSchildFamily metric{KerrSchildParams::Kerr(1.0, 0.9)};
    metric.SetParameter("mass", 1.0);
    metric.SetParameter("spin", 0.9);

    ReferenceRay ref;
    ref.r_start = 50.0;
    ref.theta_start = std::numbers::pi / 2.0;
    ref.phi_start = 0.0;
    ref.v_r = -1.0;
    ref.v_theta = 0.0;
    ref.v_phi = 0.15;  // Prograde.
    ref.expected = ReferenceRay::Outcome::Escaped;

    IntegrationResult result = integrateCPU(ref, &metric);

    EXPECT_GT(result.steps, 0) << "Kerr integration should proceed";
    EXPECT_FALSE(result.hit_horizon && result.steps < 10)
        << "Ray should not immediately hit horizon";
}

// --- Specification anchors --------------------------------------------------

TEST_F(CpuGeodesicReferenceTests, ParityToleranceSpecification) {
    EXPECT_LE(POSITION_RELATIVE_TOLERANCE, 1e-4) << "Position parity tolerance should be < 10^-4";
    EXPECT_LE(DIRECTION_TOLERANCE_RAD, 1e-4) << "Direction parity tolerance should be < 10^-4 rad";
}

TEST_F(CpuGeodesicReferenceTests, ReferenceRayGridGeneration) {
    EXPECT_EQ(referenceRays.size(), static_cast<std::size_t>(NUM_REFERENCE_RAYS))
        << "Should generate " << NUM_REFERENCE_RAYS << " reference rays";

    for (const auto& ray : referenceRays) {
        EXPECT_GT(ray.r_start, 0) << "Radius should be positive";
        EXPECT_GE(ray.theta_start, 0) << "Theta should be >= 0";
        EXPECT_LE(ray.theta_start, std::numbers::pi) << "Theta should be <= pi";
    }
}

TEST_F(CpuGeodesicReferenceTests, ReferenceRaysCoverImpactParameters) {
    std::vector<double> impact_params;
    for (const auto& ray : referenceRays) {
        double b = ray.v_phi * ray.r_start;
        impact_params.push_back(std::abs(b));
    }

    double min_b = *std::min_element(impact_params.begin(), impact_params.end());
    double max_b = *std::max_element(impact_params.begin(), impact_params.end());

    EXPECT_LT(min_b, 5.0) << "Should include rays that hit horizon";
    EXPECT_GT(max_b, 6.0) << "Should include rays that escape";
}

}  // namespace sirius::test
