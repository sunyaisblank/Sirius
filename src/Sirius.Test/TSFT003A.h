// =============================================================================
// TSFT003A.h - Geodesic Path Generator and Conservation Test Fixture
// Component ID: TSFT003A (Test/Fixtures/GeodesicGenerator)
// =============================================================================
//
// PURPOSE:
// Provides reusable infrastructure for geodesic integration tests following
// the pattern from Alaris VolatilityTestDataGenerator.
//
// CONSERVATION LAWS TESTED:
// 1. Null condition: g_μν k^μ k^ν = 0
// 2. Killing energy: E = -g_tμ k^μ (stationary spacetimes)
// 3. Killing angular momentum: L_z = g_φμ k^μ (axisymmetric spacetimes)
// 4. Carter constant: Q (Kerr spacetime)
//
// USAGE:
// class MyGeodesicTests : public sirius::test::GeodesicTestFixture { ... };
//
// REFERENCE: docs/specification.md Section "Conservation Laws"
// =============================================================================

#ifndef TSFT003A_H
#define TSFT003A_H

#include <gtest/gtest.h>
#include <PHCN001A.h>
#include <MTTN001A.h>
#include <MTDL001A.h>
#include <PHGD001A.h>
#include "TSFT001A.h"
#include "TSFT002A.h"

#include <cmath>
#include <vector>
#include <random>

namespace sirius::test {
using namespace Sirius;

// =============================================================================
// Conservation Statistics
// =============================================================================

/// @brief Statistics from conservation law monitoring during integration
struct ConservationStats {
    // Integration metadata
    int steps_completed = 0;
    int steps_requested = 0;
    bool terminated_early = false;
    std::string termination_reason;

    // Null condition statistics
    double null_condition_initial = 0.0;
    double null_condition_max = 0.0;
    double null_condition_final = 0.0;
    double null_condition_mean = 0.0;

    // Energy conservation statistics
    double energy_initial = 0.0;
    double energy_final = 0.0;
    double energy_drift_max = 0.0;
    double energy_drift_relative_max = 0.0;

    // Angular momentum conservation statistics
    double L_initial = 0.0;
    double L_final = 0.0;
    double L_drift_max = 0.0;
    double L_drift_relative_max = 0.0;

    // Position tracking
    double r_initial = 0.0;
    double r_final = 0.0;
    double r_min = 0.0;
    double r_max = 0.0;

    /// @brief Check if all conservation laws are satisfied
    bool allConserved(double null_tol = Tolerances::NULL_CONDITION_CPU,
                      double E_tol = Tolerances::CONSERVATION,
                      double L_tol = Tolerances::CONSERVATION) const {
        return (null_condition_max < null_tol) &&
               (energy_drift_relative_max < E_tol) &&
               (L_drift_relative_max < L_tol);
    }

    /// @brief Get human-readable summary
    std::string summary() const {
        std::string s = "ConservationStats (" + std::to_string(steps_completed) + "/" +
                        std::to_string(steps_requested) + " steps):\n";
        if (terminated_early) {
            s += "  Terminated: " + termination_reason + "\n";
        }
        s += "  Null condition: max=" + std::to_string(null_condition_max) + "\n";
        s += "  Energy drift: relative_max=" + std::to_string(energy_drift_relative_max) + "\n";
        s += "  L_z drift: relative_max=" + std::to_string(L_drift_relative_max) + "\n";
        s += "  Radial range: [" + std::to_string(r_min) + ", " + std::to_string(r_max) + "]\n";
        return s;
    }
};

// =============================================================================
// Lightray Configuration
// =============================================================================

/// @brief Configuration for generating test lightrays
struct LightrayConfig {
    // Starting position (spherical coordinates)
    double r = 10.0;
    double theta = MathConst::HALF_PI;
    double phi = 0.0;

    // Initial velocity direction (spherical basis)
    double v_r = -0.5;      // Radial component
    double v_theta = 0.0;   // Polar component
    double v_phi = 0.5;     // Azimuthal component

    // Integration parameters
    int max_steps = IntegrationParams::CONSERVATION_STEPS;
    double r_min_bound = 2.5;    // Terminate if r < this
    double r_max_bound = 100.0;  // Terminate if r > this

    /// @brief Create inward radial ray
    static LightrayConfig radialInward(double r_start = 10.0) {
        LightrayConfig cfg;
        cfg.r = r_start;
        cfg.v_r = -1.0;
        cfg.v_theta = 0.0;
        cfg.v_phi = 0.0;
        return cfg;
    }

    /// @brief Create outward radial ray
    static LightrayConfig radialOutward(double r_start = 5.0) {
        LightrayConfig cfg;
        cfg.r = r_start;
        cfg.v_r = 1.0;
        cfg.v_theta = 0.0;
        cfg.v_phi = 0.0;
        return cfg;
    }

    /// @brief Create equatorial orbit-like ray
    static LightrayConfig equatorialOrbital(double r_start = 10.0) {
        LightrayConfig cfg;
        cfg.r = r_start;
        cfg.theta = MathConst::HALF_PI;
        cfg.v_r = -0.3;
        cfg.v_theta = 0.0;
        cfg.v_phi = 0.8;
        return cfg;
    }

    /// @brief Create off-equatorial ray
    static LightrayConfig offEquatorial(double r_start = 15.0, double theta_start = MathConst::PI/3) {
        LightrayConfig cfg;
        cfg.r = r_start;
        cfg.theta = theta_start;
        cfg.v_r = -0.4;
        cfg.v_theta = 0.1;
        cfg.v_phi = 0.5;
        return cfg;
    }
};

// =============================================================================
// Geodesic Test Fixture
// =============================================================================

/// @brief Base fixture for geodesic integration and conservation tests
class GeodesicTestFixture : public MetricTestFixture {
protected:
    IntegratorConfig integrator_config;

    void SetUp() override {
        MetricTestFixture::SetUp();

        // Configure RK45 integrator with tight tolerances
        integrator_config = Geodesic::getDefaultConfig();
        integrator_config.abs_tolerance = 1e-8f;
        integrator_config.rel_tolerance = 1e-8f;
        integrator_config.min_step = 1e-8f;
        integrator_config.max_step = 0.05f;
        integrator_config.initial_step = 0.005f;
        integrator_config.use_rk45 = true;
    }

    // =========================================================================
    // Lightray Creation
    // =========================================================================

    /// @brief Create a lightray from configuration
    /// @param cfg Lightray configuration
    /// @param metric Metric to use for null normalization
    /// @return Initialized lightray
    template<typename MetricType>
    Lightray createLightray(const LightrayConfig& cfg, MetricType* metric) {
        Lightray ray;

        // Convert position: spherical → Cartesian
        double sin_th = std::sin(cfg.theta);
        double cos_th = std::cos(cfg.theta);
        double sin_ph = std::sin(cfg.phi);
        double cos_ph = std::cos(cfg.phi);

        ray.position(0) = 0.0;  // t = 0
        ray.position(1) = cfg.r * sin_th * cos_ph;  // x
        ray.position(2) = cfg.r * sin_th * sin_ph;  // y
        ray.position(3) = cfg.r * cos_th;           // z

        ray.proper_time = 0.0f;
        ray.coordinate_time = 0.0f;
        ray.terminated = 0;
        ray.step_size = integrator_config.initial_step;
        ray.bounce_count = 0;

        // Convert velocity: spherical basis → Cartesian basis
        // v^x = v^r sin(θ)cos(φ) + r v^θ cos(θ)cos(φ) - r v^φ sin(θ)sin(φ)
        // v^y = v^r sin(θ)sin(φ) + r v^θ cos(θ)sin(φ) + r v^φ sin(θ)cos(φ)
        // v^z = v^r cos(θ)        - r v^θ sin(θ)
        ray.velocity(1) = cfg.v_r * sin_th * cos_ph
                        + cfg.r * cfg.v_theta * cos_th * cos_ph
                        - cfg.r * cfg.v_phi * sin_th * sin_ph;
        ray.velocity(2) = cfg.v_r * sin_th * sin_ph
                        + cfg.r * cfg.v_theta * cos_th * sin_ph
                        + cfg.r * cfg.v_phi * sin_th * cos_ph;
        ray.velocity(3) = cfg.v_r * cos_th
                        - cfg.r * cfg.v_theta * sin_th;

        // Normalize to null geodesic
        Metric4D g;
        Tensor<Dual<double>, 4, 4, 4> dg;
        metric->evaluate(ray.position, g, dg);
        ray.velocity = TensorOps::normalizeNull(ray.velocity, g);

        // Compute initial acceleration
        ray.acceleration = Geodesic::calculateAcceleration(ray.velocity, ray.position, metric);

        return ray;
    }

    /// @brief Generate multiple test rays with varied configurations
    template<typename MetricType>
    std::vector<Lightray> generateTestRays(MetricType* metric, int count = 10) {
        std::vector<Lightray> rays;

        // Fixed test configurations
        rays.push_back(createLightray(LightrayConfig::radialInward(10.0), metric));
        rays.push_back(createLightray(LightrayConfig::radialOutward(5.0), metric));
        rays.push_back(createLightray(LightrayConfig::equatorialOrbital(10.0), metric));
        rays.push_back(createLightray(LightrayConfig::offEquatorial(15.0, MathConst::PI/3), metric));

        // Additional rays if requested
        if (count > 4) {
            std::mt19937 rng(42);  // Fixed seed for reproducibility
            std::uniform_real_distribution<double> r_dist(5.0, 20.0);
            std::uniform_real_distribution<double> theta_dist(MathConst::PI/6, 5*MathConst::PI/6);
            std::uniform_real_distribution<double> v_dist(-1.0, 1.0);

            for (int i = 4; i < count; ++i) {
                LightrayConfig cfg;
                cfg.r = r_dist(rng);
                cfg.theta = theta_dist(rng);
                cfg.phi = 0.0;
                cfg.v_r = v_dist(rng);
                cfg.v_theta = v_dist(rng) * 0.2;  // Smaller polar motion
                cfg.v_phi = v_dist(rng);

                // Ensure inward motion for more interesting geodesics
                if (cfg.v_r > 0) cfg.v_r = -cfg.v_r;

                rays.push_back(createLightray(cfg, metric));
            }
        }

        return rays;
    }

    // =========================================================================
    // Conservation Law Computation
    // =========================================================================

    /// @brief Compute null condition: g_μν k^μ k^ν
    double computeNullCondition(const Vec4& k, const Metric4D& g) const {
        double norm = 0.0;
        for (int mu = 0; mu < 4; ++mu) {
            for (int nu = 0; nu < 4; ++nu) {
                norm += g(mu, nu).real * k(mu) * k(nu);
            }
        }
        return norm;
    }

    /// @brief Compute Killing energy: E = -g_tμ k^μ
    /// Conserved for stationary spacetimes (∂/∂t is Killing)
    double computeKillingEnergy(const Vec4& k, const Metric4D& g) const {
        double E = 0.0;
        for (int mu = 0; mu < 4; ++mu) {
            E -= g(0, mu).real * k(mu);  // ξ^t = (1,0,0,0)
        }
        return E;
    }

    /// @brief Compute Killing angular momentum in Cartesian coords
    /// L_z = ξ_φ · k where ξ^φ = x ∂_y - y ∂_x = (0, -y, x, 0)
    double computeKillingAngularMomentum(const Vec4& k, const Metric4D& g, const Vec4& pos) const {
        double L = 0.0;
        Vec4 xi;  // Rotational Killing vector in Cartesian
        xi(0) = 0.0;
        xi(1) = -pos(2);  // -y
        xi(2) = pos(1);   // x
        xi(3) = 0.0;

        for (int mu = 0; mu < 4; ++mu) {
            for (int nu = 0; nu < 4; ++nu) {
                L += g(mu, nu).real * k(mu) * xi(nu);
            }
        }
        return L;
    }

    /// @brief Compute radial coordinate from Cartesian position
    double computeRadius(const Vec4& pos) const {
        return std::sqrt(pos(1)*pos(1) + pos(2)*pos(2) + pos(3)*pos(3));
    }

    // =========================================================================
    // Integration with Conservation Monitoring
    // =========================================================================

    /// @brief Integrate geodesic while monitoring conservation laws
    template<typename MetricType>
    ConservationStats integrateWithConservation(
        Lightray& ray,
        MetricType* metric,
        const LightrayConfig& cfg = LightrayConfig())
    {
        ConservationStats stats;
        stats.steps_requested = cfg.max_steps;

        // Initial values
        Metric4D g_init;
        Tensor<Dual<double>, 4, 4, 4> dg_init;
        metric->evaluate(ray.position, g_init, dg_init);

        stats.null_condition_initial = computeNullCondition(ray.velocity, g_init);
        stats.energy_initial = computeKillingEnergy(ray.velocity, g_init);
        stats.L_initial = computeKillingAngularMomentum(ray.velocity, g_init, ray.position);
        stats.r_initial = computeRadius(ray.position);

        stats.null_condition_max = std::abs(stats.null_condition_initial);
        stats.r_min = stats.r_initial;
        stats.r_max = stats.r_initial;

        double null_sum = std::abs(stats.null_condition_initial);
        int null_count = 1;

        // Integration loop
        for (int step = 0; step < cfg.max_steps; ++step) {
            bool success = Geodesic::integrateStepRK45(ray, metric, integrator_config);

            if (!success || ray.terminated) {
                stats.terminated_early = true;
                stats.termination_reason = ray.terminated ? "ray terminated" : "integration failed";
                break;
            }

            // Check bounds
            double r = computeRadius(ray.position);
            if (r < cfg.r_min_bound) {
                stats.terminated_early = true;
                stats.termination_reason = "horizon crossing (r < " + std::to_string(cfg.r_min_bound) + ")";
                break;
            }
            if (r > cfg.r_max_bound) {
                stats.terminated_early = true;
                stats.termination_reason = "escaped (r > " + std::to_string(cfg.r_max_bound) + ")";
                break;
            }

            stats.steps_completed++;

            // Update radial tracking
            stats.r_min = std::min(stats.r_min, r);
            stats.r_max = std::max(stats.r_max, r);

            // Evaluate metric at current position
            Metric4D g;
            Tensor<Dual<double>, 4, 4, 4> dg;
            metric->evaluate(ray.position, g, dg);

            // Monitor conservation laws
            double null = computeNullCondition(ray.velocity, g);
            double E = computeKillingEnergy(ray.velocity, g);
            double L = computeKillingAngularMomentum(ray.velocity, g, ray.position);

            stats.null_condition_max = std::max(stats.null_condition_max, std::abs(null));
            null_sum += std::abs(null);
            null_count++;

            // Energy drift
            double E_drift = std::abs(E - stats.energy_initial);
            stats.energy_drift_max = std::max(stats.energy_drift_max, E_drift);
            if (std::abs(stats.energy_initial) > 1e-10) {
                double E_rel = E_drift / std::abs(stats.energy_initial);
                stats.energy_drift_relative_max = std::max(stats.energy_drift_relative_max, E_rel);
            }

            // Angular momentum drift
            double L_drift = std::abs(L - stats.L_initial);
            stats.L_drift_max = std::max(stats.L_drift_max, L_drift);
            if (std::abs(stats.L_initial) > 1e-10) {
                double L_rel = L_drift / std::abs(stats.L_initial);
                stats.L_drift_relative_max = std::max(stats.L_drift_relative_max, L_rel);
            }

            // Update finals
            stats.null_condition_final = null;
            stats.energy_final = E;
            stats.L_final = L;
        }

        stats.r_final = computeRadius(ray.position);
        stats.null_condition_mean = null_sum / null_count;

        return stats;
    }

    // =========================================================================
    // GTest Assertion Helpers
    // =========================================================================

    /// @brief Assert null condition is preserved
    void assertNullCondition(const ConservationStats& stats,
                             double tolerance = Tolerances::NULL_CONDITION_CPU) {
        EXPECT_LT(stats.null_condition_max, tolerance)
            << "Null condition violation: max=" << stats.null_condition_max
            << " (tolerance=" << tolerance << ")";
    }

    /// @brief Assert energy is conserved
    void assertEnergyConservation(const ConservationStats& stats,
                                  double tolerance = Tolerances::CONSERVATION) {
        EXPECT_LT(stats.energy_drift_relative_max, tolerance)
            << "Energy drift: relative_max=" << stats.energy_drift_relative_max
            << " (tolerance=" << tolerance << ")";
    }

    /// @brief Assert angular momentum is conserved
    void assertAngularMomentumConservation(const ConservationStats& stats,
                                           double tolerance = Tolerances::CONSERVATION) {
        if (std::abs(stats.L_initial) > 1e-10) {
            EXPECT_LT(stats.L_drift_relative_max, tolerance)
                << "Angular momentum drift: relative_max=" << stats.L_drift_relative_max
                << " (tolerance=" << tolerance << ")";
        }
    }

    /// @brief Assert all conservation laws
    void assertAllConservationLaws(const ConservationStats& stats,
                                   double null_tol = Tolerances::NULL_CONDITION_CPU,
                                   double E_tol = Tolerances::CONSERVATION,
                                   double L_tol = Tolerances::CONSERVATION) {
        assertNullCondition(stats, null_tol);
        assertEnergyConservation(stats, E_tol);
        assertAngularMomentumConservation(stats, L_tol);
    }
};

} // namespace sirius::test

#endif // TSFT003A_H
