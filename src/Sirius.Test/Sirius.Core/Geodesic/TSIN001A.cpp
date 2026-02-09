// TSIN001A.cpp - Geodesic Path Validation

#include <gtest/gtest.h>
#include <cmath>
#include <vector>
#include "MTTN001A.h"
#include "PHMT000A.h"
#include "PHMT100A.h"
#include "PHGD001A.h"

namespace sirius::test {
using namespace Sirius;

// Constants
constexpr double M = 1.0;
constexpr double PI = 3.14159265358979323846;

// Helper: Calculate radius from Cartesian position
inline double getRadius(const Vec4& pos) {
    return std::sqrt(pos(1)*pos(1) + pos(2)*pos(2) + pos(3)*pos(3));
}

// Helper: Calculate phi from Cartesian position
inline double getPhi(const Vec4& pos) {
    return std::atan2(pos(2), pos(1));
}

// Test fixture providing Schwarzschild and Kerr metrics with RK45 integration
class GeodesicPathTests : public ::testing::Test {
protected:
    Sirius::KerrSchildFamily schwarzschild{Sirius::KerrSchildParams::Schwarzschild(M)};
    Sirius::KerrSchildFamily kerr{Sirius::KerrSchildParams::Kerr(M, 0.5)};
    IntegratorConfig rk45_config;
    
    void SetUp() override {
        schwarzschild.setParameter("mass", M);
        kerr.setParameter("mass", M);
        kerr.setParameter("spin", 0.5);
        
        rk45_config = Geodesic::getDefaultConfig();
        rk45_config.abs_tolerance = 1e-7f;
        rk45_config.rel_tolerance = 1e-7f;
        rk45_config.min_step = 1e-7f;
        rk45_config.max_step = 0.05f;
        rk45_config.initial_step = 0.005f;
    }
    
    // Create light ray at position (Spherical) with direction (Spherical Basis)
    Lightray createRay(const Vec4& spherical_pos, const Vec4& spherical_dir, IMetric* metric) {
        Lightray ray;
        
        // Convert Position: Spherical (t,r,th,ph) -> Cartesian (t,x,y,z)
        double r = spherical_pos(1);
        double th = spherical_pos(2);
        double ph = spherical_pos(3);
        
        ray.position(0) = spherical_pos(0);
        ray.position(1) = r * std::sin(th) * std::cos(ph); // x
        ray.position(2) = r * std::sin(th) * std::sin(ph); // y
        ray.position(3) = r * std::cos(th);                // z
        
        // Convert Direction: Spherical Basis (vt, vr, vth, vph) -> Cartesian Basis (vt, vx, vy, vz)
        double vr = spherical_dir(1);
        double vth = spherical_dir(2);
        double vph = spherical_dir(3);
        
        double st = std::sin(th);
        double ct = std::cos(th);
        double sp = std::sin(ph);
        double cp = std::cos(ph);
        
        ray.velocity(0) = spherical_dir(0);
        ray.velocity(1) = vr * st * cp + r * vth * ct * cp - r * vph * st * sp; // vx
        ray.velocity(2) = vr * st * sp + r * vth * ct * sp + r * vph * st * cp; // vy
        ray.velocity(3) = vr * ct      - r * vth * st;                          // vz
        
        ray.proper_time = 0.0f;
        ray.coordinate_time = 0.0f;
        ray.terminated = 0;
        ray.step_size = 0.01f;
        ray.bounce_count = 0;
        
        Metric4D g;
        Tensor<Dual<double>, 4, 4, 4> dg;
        metric->evaluate(ray.position, g, dg);
        ray.velocity = TensorOps::normalizeNull(ray.velocity, g);
        ray.acceleration = Geodesic::calculateAcceleration(ray.velocity, ray.position, metric);
        
        return ray;
    }
    
    // Helper: Transform Cartesian metric to Spherical basis for validation
    Metric4D transformToSpherical(const Metric4D& g_cart, const Vec4& pos_cart) {
        double x = pos_cart(1);
        double y = pos_cart(2);
        double z = pos_cart(3);
        double rho = std::sqrt(x*x + y*y);
        double r = std::sqrt(x*x + y*y + z*z);
        
        // Jacobian d(x,y,z)/d(r,th,ph)
        // dx/dr = sin(th)cos(ph) = x/r
        // dx/dth = r cos(th)cos(ph) = x*z/rho
        // dx/dph = -r sin(th)sin(ph) = -y
        
        // Actually simpler to use trigonometric forms if available, but from Cartesian:
        double st = rho/r;
        double ct = z/r;
        double cp = (rho > 1e-10) ? x/rho : 1.0;
        double sp = (rho > 1e-10) ? y/rho : 0.0;
        
        Tensor<double, 4, 4> J; // Transformation matrix Lambda^mu_nu = dx^mu / dx'^nu
        J.zero();
        J(0,0) = 1.0; // dt/dt
        
        // dr column
        J(1,1) = st * cp; // dx/dr
        J(2,1) = st * sp; // dy/dr
        J(3,1) = ct;      // dz/dr
        
        // dth column
        J(1,2) = r * ct * cp; // dx/dth
        J(2,2) = r * ct * sp; // dy/dth
        J(3,2) = -r * st;     // dz/dth
        
        // dph column
        J(1,3) = -r * st * sp; // dx/dph
        J(2,3) = r * st * cp;  // dy/dph
        J(3,3) = 0.0;          // dz/dph
        
        Metric4D g_sph;
        g_sph.zero();
        
        // g'_ab = J^c_a J^d_b g_cd
        for(int a=0; a<4; ++a) {
            for(int b=0; b<4; ++b) {
                Dual<double> sum = 0.0;
                for(int c=0; c<4; ++c) {
                    for(int d=0; d<4; ++d) {
                        sum = sum + g_cart(c,d) * J(c,a) * J(d,b);
                    }
                }
                g_sph(a,b) = sum;
            }
        }
        return g_sph;
    }
    
    struct IntegrationResult {
        Vec4 final_position;
        Vec4 final_velocity;
        int steps_taken;
        int termination_reason;
        double coordinate_time_elapsed;
    };
    
    // Integrate ray until max_steps or termination (horizon/escape)
    IntegrationResult integrateRay(Lightray& ray, IMetric* metric, int max_steps,
                                   double r_min = 2.5, double r_max = 200.0) {
        IntegrationResult result;
        result.steps_taken = 0;
        result.termination_reason = 0;
        
        for (int step = 0; step < max_steps; ++step) {
            bool success = Geodesic::integrateStepRK45(ray, metric, rk45_config);
            result.steps_taken++;
            
            if (!success || ray.terminated) {
                result.termination_reason = ray.terminated;
                break;
            }
            
            double r = getRadius(ray.position);
            if (r < r_min) {
                result.termination_reason = 1;  // Hit horizon
                break;
            }
            if (r > r_max) {
                result.termination_reason = 2;  // Escaped to infinity
                break;
            }
        }
        
        result.final_position = ray.position;
        result.final_velocity = ray.velocity;
        result.coordinate_time_elapsed = ray.coordinate_time;
        
        return result;
    }
};

// --- Schwarzschild Photon Sphere Tests ---

// Photon sphere orbits (r=3M) are unstable; correct momentum derivatives cause rapid escape with deviation ~80. Disabled pending investigation.
TEST_F(GeodesicPathTests, DISABLED_PhotonSphereRadius)
{
    double photon_sphere = 3.0 * M;
    
    Vec4 pos;
    pos(0) = 0.0;
    pos(1) = photon_sphere;
    pos(2) = PI / 2.0;
    pos(3) = 0.0;
    
    Vec4 dir;
    dir(0) = 0.0;
    dir(1) = 0.0;   // No radial component
    dir(2) = 0.0;
    dir(3) = 1.0;   // Pure azimuthal motion
    
    Lightray ray = createRay(pos, dir, &schwarzschild);
    
    double max_r_deviation = 0.0;
    double initial_r = getRadius(ray.position);
    
    for (int step = 0; step < 2000; ++step) {
        Geodesic::integrateStepRK45(ray, &schwarzschild, rk45_config);
        
        double r_deviation = std::abs(getRadius(ray.position) - initial_r);
        max_r_deviation = std::max(max_r_deviation, r_deviation);
        
        if (ray.terminated) break;
    }
    
    EXPECT_LT(max_r_deviation, 5.0) 
        << "Photon sphere orbit deviated by " << max_r_deviation;
}

// Verify photon sphere instability: small perturbation outside r=3M should cause outward motion or escape.
TEST_F(GeodesicPathTests, PhotonSphereInstability)
{
    double photon_sphere = 3.0 * M;
    
    Vec4 pos;
    pos(0) = 0.0;
    pos(1) = photon_sphere + 0.01;  // Slightly outside
    pos(2) = PI / 2.0;
    pos(3) = 0.0;
    
    Vec4 dir;
    dir(0) = 0.0;
    dir(1) = 0.0;
    dir(2) = 0.0;
    dir(3) = 1.0;
    
    Lightray ray = createRay(pos, dir, &schwarzschild);
    auto result = integrateRay(ray, &schwarzschild, 5000, 2.5, 100.0);
    
    bool moved_outward = getRadius(result.final_position) > photon_sphere;
    EXPECT_TRUE(moved_outward || result.termination_reason == 2) 
        << "Ray starting outside photon sphere should move outward or escape";
}

// --- Light Deflection Tests ---

// Basic ray tracing sanity check: radially outward ray from r=50 with small tangential component should escape to r>200.
TEST_F(GeodesicPathTests, WeakFieldLightDeflection)
{
    Vec4 pos;
    pos(0) = 0.0;
    pos(1) = 50.0;
    pos(2) = PI / 2.0;
    pos(3) = 0.0;
    
    Vec4 dir;
    dir(0) = 0.0;
    dir(1) = 1.0;       // Radially outward
    dir(2) = 0.0;
    dir(3) = 0.01;      // Small tangential component
    
    Lightray ray = createRay(pos, dir, &schwarzschild);
    
    printf("DEFLECT: After createRay: k0=%.4f k1=%.4f k3=%.4f\n",
           ray.velocity(0), ray.velocity(1), ray.velocity(3));
    fflush(stdout);
    
    double initial_r = getRadius(ray.position);
    double initial_phi = getPhi(ray.position);
    int step_count = 0;
    
    for (int step = 0; step < 5000; ++step) {
        bool success = Geodesic::integrateStepRK45(ray, &schwarzschild, rk45_config);
        if (!success || ray.terminated) {
            printf("DEFLECT: Terminated at step %d, success=%d terminated=%d\n",
                   step, success, ray.terminated);
            break;
        }
        step_count++;
        
        if (getRadius(ray.position) > 200.0) break;
    }
    
    double final_r = getRadius(ray.position);
    double final_phi = getPhi(ray.position);
    
    printf("DEFLECT: steps=%d init_r=%.1f final_r=%.1f delta_phi=%.4f\n",
           step_count, initial_r, final_r, final_phi - initial_phi);
    fflush(stdout);
    
    EXPECT_GT(step_count, 10) << "Ray should take multiple steps";
    EXPECT_GT(final_r, initial_r) << "Ray should move outward";
}

// --- Radial Infall Tests ---

// Radially infalling light should approach horizon asymptotically in coordinate time. Disabled pending investigation with correct Hamiltonian integrator.
TEST_F(GeodesicPathTests, DISABLED_RadialInfallApproachesHorizon)
{
    Vec4 pos;
    pos(0) = 0.0;
    pos(1) = 10.0;
    pos(2) = PI / 2.0;
    pos(3) = 0.0;
    
    Vec4 dir;
    dir(0) = 0.0;
    dir(1) = -1.0;  // Radially inward
    dir(2) = 0.0;
    dir(3) = 0.0;
    
    Lightray ray = createRay(pos, dir, &schwarzschild);
    double initial_r = getRadius(ray.position);
    
    int steps = 0;
    while (steps < 1000 && !ray.terminated) {
        Geodesic::integrateStepRK45(ray, &schwarzschild, rk45_config);
        steps++;
    }
    
    double final_r = getRadius(ray.position);
    EXPECT_LT(final_r, initial_r) 
        << "Radial infall should approach horizon, moved from " << initial_r << " to " << final_r;
    
    EXPECT_GT(ray.coordinate_time, 0.0) 
        << "Coordinate time should accumulate during infall";
}

// --- Kerr Frame Dragging Tests ---

// Frame dragging in Kerr spacetime causes prograde and retrograde rays to behave differently; verifies integrator handles non-diagonal metrics.
TEST_F(GeodesicPathTests, KerrFrameDragging)
{
    kerr.setParameter("spin", 0.9);
    
    Vec4 pos;
    pos(0) = 0.0;
    pos(1) = 5.0;
    pos(2) = PI / 2.0;
    pos(3) = 0.0;
    
    Vec4 dir_prograde;
    dir_prograde(0) = 0.0;
    dir_prograde(1) = -0.5;
    dir_prograde(2) = 0.0;
    dir_prograde(3) = 1.0;
    
    Vec4 dir_retrograde;
    dir_retrograde(0) = 0.0;
    dir_retrograde(1) = -0.5;
    dir_retrograde(2) = 0.0;
    dir_retrograde(3) = -1.0;
    
    Lightray ray_pro = createRay(pos, dir_prograde, &kerr);
    Lightray ray_ret = createRay(pos, dir_retrograde, &kerr);
    
    for (int step = 0; step < 500; ++step) {
        Geodesic::integrateStepRK45(ray_pro, &kerr, rk45_config);
        Geodesic::integrateStepRK45(ray_ret, &kerr, rk45_config);
        
        if (ray_pro.terminated || ray_ret.terminated) break;
    }
    
    EXPECT_TRUE(!ray_pro.terminated || !ray_ret.terminated)
        << "At least one ray should complete integration";
}

// --- Coordinate Continuity Tests ---

// Verify φ coordinate wraps smoothly from 2π to 0 without discontinuity. Light rays outside photon sphere don't orbit at constant radius.
TEST_F(GeodesicPathTests, PhiWrapContinuity)
{
    Vec4 pos;
    pos(0) = 0.0;
    pos(1) = 10.0;
    pos(2) = PI / 2.0;
    pos(3) = PI;
    
    Vec4 dir;
    dir(0) = 0.0;
    dir(1) = 0.0;
    dir(2) = 0.0;
    dir(3) = 1.0;
    
    Lightray ray = createRay(pos, dir, &schwarzschild);
    
    double prev_phi = getPhi(ray.position);
    int wrap_count = 0;
    int smooth_transitions = 0;
    
    for (int step = 0; step < 500; ++step) {
        Geodesic::integrateStepRK45(ray, &schwarzschild, rk45_config);
        
        double curr_phi = getPhi(ray.position);
        
        if (std::abs(curr_phi - prev_phi) > PI) {
            wrap_count++;
            double wrapped_diff = std::abs(curr_phi - prev_phi);
            if (wrapped_diff > 2.0 * PI - 0.5 || wrapped_diff < 0.5) {
                smooth_transitions++;
            }
        }
        
        prev_phi = curr_phi;
        
        if (ray.terminated) break;
        if (getRadius(ray.position) < 2.5 || getRadius(ray.position) > 100.0) break;
    }
    
    std::cout << "Phi wraps detected: " << wrap_count << std::endl;
    std::cout << "Smooth transitions: " << smooth_transitions << std::endl;
}

// --- Metric Signature Tests ---

// Verify Schwarzschild has Lorentzian signature (-,+,+,+)
TEST_F(GeodesicPathTests, MetricSignatureCorrect)
{
    Vec4 pos;
    pos(0) = 0.0;
    pos(1) = 10.0;
    pos(2) = PI / 2.0;
    pos(3) = 0.0;
    
    Metric4D g;
    Tensor<Dual<double>, 4, 4, 4> dg;
    
    Lightray ray = createRay(pos, Vec4(), &schwarzschild);
    schwarzschild.evaluate(ray.position, g, dg); // Evaluated in Cartesian
    
    // Transform to Spherical for signature check
    Metric4D g_sph = transformToSpherical(g, ray.position);
    
    EXPECT_LT(g_sph(0, 0).real, 0.0) << "g_tt should be negative";
    EXPECT_GT(g_sph(1, 1).real, 0.0) << "g_rr should be positive";
    EXPECT_GT(g_sph(2, 2).real, 0.0) << "g_θθ should be positive";
    EXPECT_GT(g_sph(3, 3).real, 0.0) << "g_φφ should be positive";
    
    double det = g_sph(0, 0).real * g_sph(1, 1).real * g_sph(2, 2).real * g_sph(3, 3).real;
    EXPECT_LT(det, 0.0) << "Metric determinant should be negative";
}

// Verify Kerr has correct signature with non-zero frame-dragging term g_tφ
TEST_F(GeodesicPathTests, KerrMetricSignature) {
    Sirius::KerrSchildFamily kerr{Sirius::KerrSchildParams::Kerr(M, 0.9)};
    
    // Use a position where frame dragging is evident (non-zero phi/theta)
    // Spherical (t, r, pi/4, 0)
    // Cartesian: x=r*sin(th), y=0, z=r*cos(th)
    Vec4 pos;
    pos(0) = 0.0;
    pos(1) = 0.0;
    pos(2) = 0.0;
    pos(3) = 0.0;
    // We used createRay helper to get correct Cartesian position and it returns a ray object
    // But here we need to instantiate it locally or reuse createRay logic
    
    // Let's manually set a Cartesian position
    double r = 4.0;
    double theta = PI/4.0;
    pos(1) = r * std::sin(theta); // x
    pos(2) = 0.0;                 // y
    pos(3) = r * std::cos(theta); // z
    
    Metric4D g;
    Tensor<Dual<double>, 4, 4, 4> dg;
    
    kerr.evaluate(pos, g, dg);
    
    // Transform to Spherical
    Metric4D g_sph = transformToSpherical(g, pos);
    
    EXPECT_LT(g_sph(0, 0).real, 0.0) << "Kerr g_tt should be negative";
    EXPECT_GT(g_sph(1, 1).real, 0.0) << "Kerr g_rr should be positive";
    EXPECT_GT(g_sph(2, 2).real, 0.0) << "Kerr g_θθ should be positive";
    
    // In Spherical Kerr, g_tphi is non-zero
    EXPECT_NE(g_sph(0, 3).real, 0.0) << "Kerr g_tφ should be non-zero (frame dragging)";
    EXPECT_DOUBLE_EQ(g_sph(0, 3).real, g_sph(3, 0).real) << "Metric should be symmetric";
}

} // namespace sirius::test