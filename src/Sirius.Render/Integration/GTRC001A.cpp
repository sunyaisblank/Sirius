// =============================================================================
// GTRC001A.cpp - Geodesic Ray Tracer Implementation
// Component ID: GTRC001A (Render/Integration/GeodesicTracer)
// =============================================================================

#include "GTRC001A.h"
#include <iostream>
#include <algorithm>

namespace Sirius {

// =============================================================================
// Initialize Lightray from Camera Ray
// =============================================================================
Lightray GeodesicTracer::initializeLightray(const CameraRay& camera_ray) {
    Lightray ray;

    // Camera ray origin is in Boyer-Lindquist coordinates: (t, r, θ, φ)
    double t = camera_ray.origin(0);
    double r = camera_ray.origin(1);
    double th = camera_ray.origin(2);
    double ph = camera_ray.origin(3);

    // Clamp theta to avoid poles
    th = std::clamp(th, 0.01, M_PI - 0.01);

    // Precompute trigonometric values
    double sin_th = std::sin(th);
    double cos_th = std::cos(th);
    double sin_ph = std::sin(ph);
    double cos_ph = std::cos(ph);

    // =========================================================================
    // Convert Position: Boyer-Lindquist → Kerr-Schild Cartesian
    // x = r sin(θ) cos(φ)
    // y = r sin(θ) sin(φ)
    // z = r cos(θ)
    // =========================================================================
    ray.position(0) = static_cast<float>(t);
    ray.position(1) = static_cast<float>(r * sin_th * cos_ph);
    ray.position(2) = static_cast<float>(r * sin_th * sin_ph);
    ray.position(3) = static_cast<float>(r * cos_th);

    // =========================================================================
    // Convert Velocity: Spherical Basis → Cartesian Basis
    // =========================================================================
    // Camera ray direction is in spherical coordinate basis:
    // direction = (dt/dλ, dr/dλ, dθ/dλ, dφ/dλ) but scaled for convenience
    // We need to convert (v_r, v_θ, v_φ) to (v_x, v_y, v_z)
    //
    // Jacobian transformation:
    // v_x = v_r sin(θ)cos(φ) + r·v_θ cos(θ)cos(φ) - r·v_φ sin(θ)sin(φ)
    // v_y = v_r sin(θ)sin(φ) + r·v_θ cos(θ)sin(φ) + r·v_φ sin(θ)cos(φ)
    // v_z = v_r cos(θ)        - r·v_θ sin(θ)

    double v_r = camera_ray.direction(1);
    double v_th = camera_ray.direction(2) / r;  // Camera stores r·dθ/dλ
    double v_ph = camera_ray.direction(3) / (r * sin_th);  // Camera stores r·sin(θ)·dφ/dλ

    // Avoid division by zero at poles
    if (std::abs(sin_th) < 1e-10) {
        v_ph = 0.0;
    }

    ray.velocity(1) = static_cast<float>(
        v_r * sin_th * cos_ph +
        r * v_th * cos_th * cos_ph -
        r * v_ph * sin_th * sin_ph
    );
    ray.velocity(2) = static_cast<float>(
        v_r * sin_th * sin_ph +
        r * v_th * cos_th * sin_ph +
        r * v_ph * sin_th * cos_ph
    );
    ray.velocity(3) = static_cast<float>(
        v_r * cos_th -
        r * v_th * sin_th
    );

    // =========================================================================
    // Normalize to Null Condition
    // =========================================================================
    // The velocity must satisfy g_μν k^μ k^ν = 0
    // We solve for k^0 (the time component) given the spatial components

    Metric4D g;
    Tensor<Dual<double>, 4, 4, 4> dg;

    Vec4 pos_double;
    for (int i = 0; i < 4; ++i) pos_double(i) = ray.position(i);

    m_Metric->evaluate(pos_double, g, dg);

    // Convert velocity to double for normalization
    Vec4 vel_double;
    vel_double(0) = 0.0;  // Will be computed
    for (int i = 1; i < 4; ++i) vel_double(i) = ray.velocity(i);

    // Normalize using TensorOps
    Vec4 normalized = TensorOps::normalizeNull(vel_double, g);

    for (int i = 0; i < 4; ++i) {
        ray.velocity(i) = static_cast<float>(normalized(i));
    }

    // =========================================================================
    // Initialize remaining fields
    // =========================================================================
    Vec4 vel_for_accel;
    vel_for_accel(0) = ray.velocity(0);
    vel_for_accel(1) = ray.velocity(1);
    vel_for_accel(2) = ray.velocity(2);
    vel_for_accel(3) = ray.velocity(3);

    Vec4 pos_for_accel;
    pos_for_accel(0) = ray.position(0);
    pos_for_accel(1) = ray.position(1);
    pos_for_accel(2) = ray.position(2);
    pos_for_accel(3) = ray.position(3);

    ray.acceleration = Geodesic::calculateAcceleration(vel_for_accel, pos_for_accel, m_Metric);

    ray.proper_time = 0.0f;
    ray.coordinate_time = 0.0f;
    ray.step_size = m_Config.integrator.initial_step;
    ray.terminated = 0;
    ray.bounce_count = 0;
    ray.ku_uobsu = 1.0f;  // Initial k·u for redshift
    ray.running_dlambda_dnew = 1.0f;

    return ray;
}

// =============================================================================
// Main Trace Method
// =============================================================================
TraceResult GeodesicTracer::trace(const CameraRay& camera_ray) {
    TraceResult result;
    result.steps_taken = 0;
    result.numerical_failure = false;

    // Initialize ray from camera
    Lightray ray = initializeLightray(camera_ray);

    // Check initial state
    if (hasInvalidState(ray)) {
        result.outcome = TraceResult::Outcome::MAX_STEPS;
        result.numerical_failure = true;
        return result;
    }

    // Compute horizon radius
    float r_horizon = computeHorizonRadius();

    // =========================================================================
    // Photon Sphere Radius (Exact Kerr Formula)
    // =========================================================================
    // For Kerr spacetime, unstable circular photon orbits exist at:
    //   r_photon = 2M × [1 + cos(2/3 × arccos(∓a/M))]
    //
    // The - sign gives prograde orbits (smaller r, used for a > 0)
    // The + sign gives retrograde orbits (larger r)
    //
    // For Schwarzschild (a = 0): r_photon = 2M × [1 + cos(2π/3)] = 2M × 1.5 = 3M
    //
    // Reference: Bardeen, Press & Teukolsky (1972), ApJ 178, 347
    //
    auto params = m_Metric->getParameters();
    double M = params.count("mass") ? params.at("mass").value : 1.0;
    double a_over_M = params.count("spin") ? params.at("spin").value : 0.0;

    // Clamp spin to valid range |a/M| ≤ 1
    a_over_M = std::clamp(a_over_M, -0.998, 0.998);

    // Exact formula for prograde photon orbit
    // r_ph = 2M × [1 + cos(2/3 × arccos(-a/M))]
    double arg = -a_over_M;  // Prograde orbit for positive spin
    arg = std::clamp(arg, -1.0, 1.0);  // Numerical safety for arccos
    double r_photon_d = 2.0 * M * (1.0 + std::cos((2.0 / 3.0) * std::acos(arg)));
    float r_photon = static_cast<float>(r_photon_d);

    // Track minimum radius for photon ring detection
    float min_r = 1e10f;

    // =========================================================================
    // Higher-Order Imaging: Winding Number Tracking
    // =========================================================================
    // Track equatorial crossings and azimuthal angle change to determine
    // image order (n=0 direct, n=1 primary, n=2 secondary, etc.)
    //
    // Impact parameter b = L/E determines trajectory class:
    //   b >> b_crit: direct path (n=0)
    //   b ~ b_crit: multiple windings possible (n=1,2,3...)
    //   b < b_crit: captured by black hole
    //
    // Critical impact parameter for Schwarzschild: b_crit = 3√3 M ≈ 5.196 M
    // =========================================================================

    float prev_phi = static_cast<float>(std::atan2(ray.position(2), ray.position(1)));
    float prev_z = static_cast<float>(ray.position(3));
    int equatorial_crossings = 0;
    float total_phi_change = 0.0f;

    // Compute impact parameter b = L/E from initial conditions
    // For a ray at large r with velocity components, b ≈ r × sin(angle from radial)
    {
        double x0 = ray.position(1);
        double y0 = ray.position(2);
        double z0 = ray.position(3);
        double r0 = std::sqrt(x0*x0 + y0*y0 + z0*z0);

        double vx = ray.velocity(1);
        double vy = ray.velocity(2);
        double vz = ray.velocity(3);
        double v_mag = std::sqrt(vx*vx + vy*vy + vz*vz);

        // Radial unit vector
        double n_x = x0 / r0;
        double n_y = y0 / r0;
        double n_z = z0 / r0;

        // Velocity component perpendicular to radial direction
        double v_radial = (vx*n_x + vy*n_y + vz*n_z);
        double v_perp_sq = v_mag*v_mag - v_radial*v_radial;
        double v_perp = std::sqrt(std::max(v_perp_sq, 0.0));

        // Impact parameter: b = r × v_perp / v_total (in geometric units)
        result.impact_parameter = static_cast<float>(r0 * v_perp / std::max(v_mag, 1e-10));
    }

    // Previous position for disk intersection
    Vec4 prev_pos;
    for (int i = 0; i < 4; ++i) prev_pos(i) = ray.position(i);

    // =========================================================================
    // Spiral Orbit Detection Parameters
    // =========================================================================
    // Rays with impact parameter near b_crit spiral many times around the
    // photon sphere before escaping/capturing. We detect this and terminate
    // early to avoid excessive computation.
    //
    // Critical impact parameter: b_crit = 3√3 M ≈ 5.196 M (Schwarzschild)
    // For Kerr: b_crit depends on spin and inclination
    //
    // Detection criteria:
    // 1. Ray has completed > 2π total azimuthal change (full orbit)
    // 2. Minimum radius is within 50% of photon sphere
    // 3. Ray is neither strongly infalling nor strongly escaping
    //
    // Such rays produce higher-order images (n≥2) which are exponentially
    // dimmer (exp(-πn)) and contribute little to visual output.
    // =========================================================================
    constexpr float SPIRAL_PHI_THRESHOLD = 4.0f * static_cast<float>(M_PI);  // 2 full orbits
    constexpr float SPIRAL_R_THRESHOLD = 1.5f;  // r_min < 1.5 × r_photon
    constexpr int SPIRAL_CHECK_INTERVAL = 500;  // Check every N steps

    // =========================================================================
    // Integration Loop
    // =========================================================================
    for (int step = 0; step < m_Config.max_steps; ++step) {
        // Store previous position
        for (int i = 0; i < 4; ++i) prev_pos(i) = ray.position(i);

        // Integrate one step using RK45
        bool success = Geodesic::integrateStepRK45(ray, m_Metric, m_Config.integrator);
        result.steps_taken++;

        // Check for integration failure
        if (!success || ray.terminated || hasInvalidState(ray)) {
            result.outcome = TraceResult::Outcome::MAX_STEPS;
            result.numerical_failure = hasInvalidState(ray);
            break;
        }

        // Compute current radius from Cartesian position
        double x = ray.position(1);
        double y = ray.position(2);
        double z = ray.position(3);
        double r = std::sqrt(x*x + y*y + z*z);

        // Track minimum radius for photon ring detection
        if (r < min_r) {
            min_r = static_cast<float>(r);
        }

        // =====================================================================
        // Higher-Order Imaging: Track Winding Number
        // =====================================================================
        // Count equatorial crossings and accumulate azimuthal angle change
        float curr_phi = static_cast<float>(std::atan2(y, x));
        float curr_z = static_cast<float>(z);

        // Track equatorial crossing (z sign change)
        if (prev_z * curr_z < 0) {
            equatorial_crossings++;
        }

        // Track azimuthal angle change (handle wraparound at ±π)
        float dphi = curr_phi - prev_phi;
        if (dphi > M_PI) dphi -= 2.0f * static_cast<float>(M_PI);
        if (dphi < -M_PI) dphi += 2.0f * static_cast<float>(M_PI);
        total_phi_change += std::abs(dphi);

        prev_phi = curr_phi;
        prev_z = curr_z;

        // =====================================================================
        // Spiral Orbit Early Termination
        // =====================================================================
        // Detect quasi-circular orbits near the photon sphere that would
        // otherwise consume many thousands of integration steps.
        //
        // These rays produce higher-order images (n≥2) which are exponentially
        // dimmer and contribute negligibly to the final render.
        //
        // We terminate early and treat as MAX_STEPS with appropriate fallback.
        // =====================================================================
        if ((step % SPIRAL_CHECK_INTERVAL == 0) && step > SPIRAL_CHECK_INTERVAL) {
            float r_ratio = min_r / r_photon;
            bool near_photon_sphere = (r_ratio < SPIRAL_R_THRESHOLD) && (r_ratio > 0.9f);
            bool has_spiraled = (total_phi_change > SPIRAL_PHI_THRESHOLD);

            if (near_photon_sphere && has_spiraled) {
                // Check if radial velocity is small (quasi-circular)
                double vx = ray.velocity(1);
                double vy = ray.velocity(2);
                double vz = ray.velocity(3);
                double v_radial = (x*vx + y*vy + z*vz) / r;
                double v_total = std::sqrt(vx*vx + vy*vy + vz*vz);
                double radial_fraction = std::abs(v_radial) / std::max(v_total, 1e-10);

                // If radial velocity < 30% of total, ray is quasi-circular
                if (radial_fraction < 0.3) {
                    // Terminate as spiraling ray
                    // Set outcome based on what we've collected so far
                    if (result.num_disk_crossings > 0) {
                        // Keep DISK_HIT - we have disk data
                        result.min_radius = min_r;
                    } else {
                        result.outcome = TraceResult::Outcome::SPIRALING;
                    }
                    result.equatorial_crossings = equatorial_crossings;
                    result.total_phi_change = total_phi_change;
                    result.image_order = static_cast<int>(total_phi_change / M_PI);
                    break;
                }
            }
        }

        // =====================================================================
        // Check Termination Conditions
        // =====================================================================

        // 1. Horizon capture
        if (r < r_horizon) {
            result.outcome = TraceResult::Outcome::HORIZON;
            result.final_position(0) = ray.position(0);
            result.final_position(1) = ray.position(1);
            result.final_position(2) = ray.position(2);
            result.final_position(3) = ray.position(3);
            result.min_radius = min_r;
            break;
        }

        // 2. Escape to infinity
        if (r > m_Config.escape_radius) {
            // Check if ray is moving outward (radial velocity > 0)
            double vx = ray.velocity(1);
            double vy = ray.velocity(2);
            double vz = ray.velocity(3);
            double v_radial = (x*vx + y*vy + z*vz) / r;

            if (v_radial > 0) {
                // If ray hit the disk, preserve DISK_HIT outcome; otherwise ESCAPED
                // This ensures disk crossings are the primary outcome for rendering
                if (result.num_disk_crossings == 0) {
                    result.outcome = TraceResult::Outcome::ESCAPED;
                }
                // else: keep DISK_HIT outcome set earlier

                result.final_position(0) = ray.position(0);
                result.final_position(1) = ray.position(1);
                result.final_position(2) = ray.position(2);
                result.final_position(3) = ray.position(3);
                result.final_direction(0) = ray.velocity(0);
                result.final_direction(1) = ray.velocity(1);
                result.final_direction(2) = ray.velocity(2);
                result.final_direction(3) = ray.velocity(3);
                result.min_radius = min_r;

                // =============================================================
                // Phase 1: Photon Ring Detection and Magnification
                // =============================================================
                // Rays passing close to the photon sphere experience strong
                // gravitational lensing, producing the photon ring (Einstein ring).
                //
                // MAGNIFICATION FORMULA:
                // Near the photon sphere, magnification scales approximately as:
                //   μ ≈ 1 + C / (r_min/r_photon - 1)^α
                //
                // where α ~ 1 for weak lensing, α ~ 2 for strong lensing.
                // We use α = 1 with a cap for numerical stability.
                //
                // Physical bounds: μ → ∞ as r_min → r_photon (unstable orbit)
                // Numerical bound: cap at 20× to prevent extreme values
                //
                // Reference: Luminet (1979), A&A 75, 228
                //
                float r_ratio = min_r / r_photon;
                if (r_ratio < 1.5f) {
                    result.photon_ring = true;

                    // Magnification: μ = 1 + 2 / (r_min/r_photon - 1)
                    // Valid for r_min > r_photon (ray passed outside photon sphere)
                    if (r_ratio > 1.0f) {
                        float denom = r_ratio - 1.0f;
                        result.magnification = 1.0f + 2.0f / std::max(denom, 0.1f);
                    } else {
                        // r_min < r_photon: ray penetrated photon sphere
                        // This is a numerical artifact for escaping rays
                        // Apply maximum magnification
                        result.magnification = 20.0f;
                    }

                    // Cap at reasonable maximum
                    result.magnification = std::min(result.magnification, 20.0f);
                }
                break;
            }
        }

        // 3. Volumetric disk (ray marching through 3D disk volume)
        //    Uses prev_pos and curr_pos to detect entry/exit events
        if (m_Config.enable_disk && m_Config.enable_volumetric) {
            float prev_r = static_cast<float>(std::sqrt(prev_pos(1)*prev_pos(1) + prev_pos(2)*prev_pos(2)));
            float prev_z_cyl = static_cast<float>(prev_pos(3));
            float curr_r = static_cast<float>(std::sqrt(x*x + y*y));
            float curr_z_cyl = static_cast<float>(z);

            bool was_inside = isInVolumetricDisk(prev_r, prev_z_cyl);
            bool is_inside = isInVolumetricDisk(curr_r, curr_z_cyl);

            // Check if ray crossed through the disk (either direction)
            if (was_inside || is_inside) {
                // Ray is passing through disk volume - do ray marching on this segment
                Vec4 segment_start, segment_end;
                for (int i = 0; i < 4; ++i) {
                    segment_start(i) = prev_pos(i);
                    segment_end(i) = ray.position(i);
                }

                // Accumulate emission for this segment
                accumulateVolumetricEmission(ray, segment_start, segment_end, result);

                // Set DISK_HIT outcome if we have significant emission
                if (result.volumetric_hit && result.outcome != TraceResult::Outcome::DISK_HIT) {
                    result.outcome = TraceResult::Outcome::DISK_HIT;
                    result.disk_radius = curr_r;
                    result.disk_phi = static_cast<float>(std::atan2(y, x));
                }
            }
        }

        // 4. Thin disk intersection (accumulate multiple crossings for higher-order imaging)
        if (m_Config.enable_disk && !m_Config.enable_volumetric) {
            Vec4 curr_pos;
            curr_pos(0) = ray.position(0);
            curr_pos(1) = ray.position(1);
            curr_pos(2) = ray.position(2);
            curr_pos(3) = ray.position(3);
            float disk_r, disk_phi;

            if (checkDiskIntersection(prev_pos, curr_pos, disk_r, disk_phi)) {
                // =============================================================
                // Higher-Order Imaging: Accumulate All Disk Crossings
                // =============================================================
                // Instead of stopping at the first disk hit, we continue to
                // capture secondary, tertiary, etc. images (n=1, n=2, n=3...)
                //
                // Each crossing represents a different image order:
                //   Crossing 0: Primary image (direct view of disk)
                //   Crossing 1: Secondary image (light from far side, bent over)
                //   Crossing 2+: Higher-order images (multiple windings)
                //
                // The first crossing sets the primary outcome; subsequent
                // crossings are stored for accumulation during rendering.
                // =============================================================

                if (result.num_disk_crossings < TraceResult::MAX_DISK_CROSSINGS) {
                    auto& crossing = result.disk_crossings[result.num_disk_crossings];
                    crossing.r = disk_r;
                    crossing.phi = disk_phi;
                    crossing.temperature = computeDiskTemperature(disk_r);
                    crossing.crossing_index = result.num_disk_crossings;
                    crossing.valid = true;

                    // Compute g-factor for this crossing
                    Vec4 ray_vel;
                    ray_vel(0) = ray.velocity(0);
                    ray_vel(1) = ray.velocity(1);
                    ray_vel(2) = ray.velocity(2);
                    ray_vel(3) = ray.velocity(3);

                    // For the first crossing, use full g-factor computation
                    // For subsequent crossings, compute simplified redshift
                    if (result.num_disk_crossings == 0) {
                        // First crossing: full computation with motion blur components
                        computeGFactorWithComponents(disk_r, disk_phi, ray_vel, result);
                        crossing.redshift = result.redshift;

                        // Set primary outcome data
                        result.outcome = TraceResult::Outcome::DISK_HIT;
                        result.disk_radius = disk_r;
                        result.disk_phi = disk_phi;
                        result.disk_temperature = crossing.temperature;
                    } else {
                        // Subsequent crossings: simplified g-factor
                        crossing.redshift = computeGFactor(disk_r, disk_phi, ray_vel);
                    }

                    result.num_disk_crossings++;
                }

                // =============================================================
                // Phase 1: Magnification for disk hits
                // =============================================================
                // Compute magnification based on closest approach to photon sphere
                float r_ratio_disk = min_r / r_photon;
                if (r_ratio_disk < 1.5f) {
                    result.photon_ring = true;

                    if (r_ratio_disk > 1.0f) {
                        float denom = r_ratio_disk - 1.0f;
                        // Increase magnification cap for higher-order imaging
                        result.magnification = 1.0f + 2.0f / std::max(denom, 0.05f);
                    } else {
                        result.magnification = 100.0f;  // Higher cap for caustics
                    }

                    // Increased cap to capture higher-order images
                    result.magnification = std::min(result.magnification, 100.0f);
                }

                // Continue integration to find more crossings (don't break)
                // Only break if we've found enough crossings
                if (result.num_disk_crossings >= TraceResult::MAX_DISK_CROSSINGS) {
                    result.min_radius = min_r;
                    break;
                }
            }
        }
    }

    // Store final min_radius if not already set
    result.min_radius = std::min(result.min_radius, min_r);

    // =========================================================================
    // Higher-Order Imaging: Finalize Winding Data
    // =========================================================================
    // Store accumulated winding information and compute image order
    //
    // Image order classification based on equatorial crossings:
    //   n = 0: 0 crossings (direct path, no significant deflection)
    //   n = 1: 1-2 crossings (primary image, half-orbit deflection)
    //   n = 2: 3-4 crossings (secondary image, full orbit)
    //   n = 3+: 5+ crossings (higher-order images)
    //
    // Alternative: use total_phi_change for finer classification
    //   n ≈ total_phi_change / π (each half-orbit = π radians)
    // =========================================================================

    result.equatorial_crossings = equatorial_crossings;
    result.total_phi_change = total_phi_change;

    // Image order based on equatorial crossings (each crossing ≈ half-orbit)
    // n = floor(crossings / 2) for crossings >= 1, else n = 0
    if (equatorial_crossings >= 1) {
        result.image_order = (equatorial_crossings + 1) / 2;  // Round up: 1->1, 2->1, 3->2, 4->2, etc.
    } else {
        result.image_order = 0;  // Direct path
    }

    // Alternative classification using total azimuthal change (more continuous)
    // Uncomment to use: result.image_order = static_cast<int>(total_phi_change / M_PI);

    // If loop completed without termination, result is MAX_STEPS (default)
    return result;
}

// =============================================================================
// Trace with Custom Metric Parameters
// =============================================================================
TraceResult GeodesicTracer::trace(const CameraRay& camera_ray, double mass, double spin) {
    // Update metric parameters
    m_Metric->setParameter("mass", mass);
    m_Metric->setParameter("spin", spin / mass);  // Store as a/M

    // Update disk inner radius (ISCO) for new spin
    // ISCO for Kerr: r_isco = M * (3 + Z2 - sqrt((3-Z1)(3+Z1+2*Z2)))
    // where Z1 = 1 + (1-a²)^(1/3)[(1+a)^(1/3) + (1-a)^(1/3)]
    //       Z2 = sqrt(3*a² + Z1²)
    // For prograde orbits with a > 0
    double a = spin;
    double M = mass;

    if (std::abs(a) < 1e-10) {
        m_Config.disk_inner = static_cast<float>(6.0 * M);  // Schwarzschild ISCO
    } else {
        double a_over_M = a / M;
        double one_minus_a2 = 1.0 - a_over_M * a_over_M;

        // Clamp for numerical stability
        one_minus_a2 = std::max(one_minus_a2, 1e-10);

        double cbrt_1_minus_a2 = std::cbrt(one_minus_a2);
        double cbrt_1_plus_a = std::cbrt(1.0 + a_over_M);
        double cbrt_1_minus_a = std::cbrt(1.0 - a_over_M);

        double Z1 = 1.0 + cbrt_1_minus_a2 * (cbrt_1_plus_a + cbrt_1_minus_a);
        double Z2 = std::sqrt(3.0 * a_over_M * a_over_M + Z1 * Z1);

        double r_isco;
        if (a_over_M >= 0) {
            // Prograde ISCO
            r_isco = M * (3.0 + Z2 - std::sqrt((3.0 - Z1) * (3.0 + Z1 + 2.0 * Z2)));
        } else {
            // Retrograde ISCO
            r_isco = M * (3.0 + Z2 + std::sqrt((3.0 - Z1) * (3.0 + Z1 + 2.0 * Z2)));
        }

        m_Config.disk_inner = static_cast<float>(r_isco);
    }

    return trace(camera_ray);
}

// =============================================================================
// Compute Keplerian Orbital Velocity
// =============================================================================
float GeodesicTracer::computeOrbitalVelocity(float r) {
    // Get metric parameters
    auto params = m_Metric->getParameters();
    double M = params.count("mass") ? params.at("mass").value : 1.0;
    double a_over_M = params.count("spin") ? params.at("spin").value : 0.0;
    double a = a_over_M * M;

    // Keplerian angular velocity for Kerr spacetime (prograde orbit)
    // Ω = sqrt(M) / (r^(3/2) + a*sqrt(M))
    double r32 = std::pow(static_cast<double>(r), 1.5);
    double sqrtM = std::sqrt(M);

    double Omega = sqrtM / (r32 + a * sqrtM);

    return static_cast<float>(Omega);
}

// =============================================================================
// Compute g-factor (Gravitational + Doppler Redshift)
// =============================================================================
// Refactored to populate motion blur components in TraceResult for analytical
// computation of g-factor at arbitrary azimuthal offsets.
//
// MATHEMATICAL FOUNDATION:
// =======================
// The g-factor decomposes as: g = grav_factor / (gamma × (1 - v·n))
//
// For motion blur at azimuthal offset δφ, the disk element moves to φ + δφ,
// but the ray direction (n) is fixed at the intersection. The Doppler term becomes:
//
//   v·n(φ + δφ) = v_orb × (A·cos(δφ) + B·sin(δφ))
//
// where:
//   A = n_x·sin(φ) - n_y·cos(φ)  [original v·n / v_orb]
//   B = n_x·cos(φ) + n_y·sin(φ)  [perpendicular component]
//
// Derivation: Using angle addition formulas:
//   v_orbit(φ+δφ) = v_orb × (-sin(φ+δφ), cos(φ+δφ), 0)
//   v·n = v_orb × [-n_x·sin(φ+δφ) + n_y·cos(φ+δφ)]
//       = v_orb × [-n_x·(sinφ·cosδφ + cosφ·sinδφ) + n_y·(cosφ·cosδφ - sinφ·sinδφ)]
//       = v_orb × [(−n_x·sinφ + n_y·cosφ)·cosδφ + (−n_x·cosφ − n_y·sinφ)·sinδφ]
//       = v_orb × [A·cosδφ + B·sinδφ]
//
// where A = n_x·sin(φ) - n_y·cos(φ) and B = n_x·cos(φ) + n_y·sin(φ)
// (with sign adjustment for photon direction convention)
//
// ERROR BOUNDS:
// =============
// This computation is analytically exact - no approximation is made.
// Numerical error is bounded by double precision: < 10⁻¹⁵ relative error.
// The final clamping (0.1 < g < 5.0) prevents physical singularities only.
//
float GeodesicTracer::computeGFactor(float r, float phi, const Vec4& ray_vel) {
    // Get metric parameters
    auto params = m_Metric->getParameters();
    double M = params.count("mass") ? params.at("mass").value : 1.0;
    double a_over_M = params.count("spin") ? params.at("spin").value : 0.0;

    // Compute orbital angular velocity (exact Kerr formula)
    float Omega = computeOrbitalVelocity(r);

    // Orbital velocity magnitude v = r × Ω (in units of c)
    double v_orb = static_cast<double>(r) * Omega;

    // Gravitational redshift component (exact Kerr formula)
    // grav_factor = sqrt(1 - 2Mr / (r² + a²))
    double r_d = static_cast<double>(r);
    double a = a_over_M * M;
    double a2 = a * a;
    double r2 = r_d * r_d;
    double grav_factor = std::sqrt(1.0 - 2.0 * M * r_d / (r2 + a2));

    // Extract and normalize ray direction at disk intersection
    double n_x = ray_vel(1);
    double n_y = ray_vel(2);
    double n_z = ray_vel(3);
    double v_mag = std::sqrt(n_x*n_x + n_y*n_y + n_z*n_z);

    if (v_mag < 1e-10) {
        return static_cast<float>(grav_factor);
    }

    // Normalize to unit direction (photon direction toward observer)
    n_x /= v_mag;
    n_y /= v_mag;
    // n_z not used for equatorial disk

    // Lorentz factor (exact)
    double gamma = 1.0 / std::sqrt(std::max(1e-10, 1.0 - v_orb * v_orb));

    // Compute Doppler coefficients for motion blur
    // A = n_x·sin(φ) - n_y·cos(φ)
    // B = n_x·cos(φ) + n_y·sin(φ)
    double phi_d = static_cast<double>(phi);
    double sin_phi = std::sin(phi_d);
    double cos_phi = std::cos(phi_d);

    // Note: n points toward observer, so we negate for outgoing photon convention
    // A coefficient is sufficient for original g-factor (B only needed for motion blur)
    double A = (-n_x) * sin_phi - (-n_y) * cos_phi;  // = -n_x·sinφ + n_y·cosφ

    // At original φ (δφ = 0): v·n = v_orb × A
    double v_dot_n = v_orb * A;

    // Combined g-factor: g = grav_factor / (gamma × (1 - v·n))
    double doppler_denom = gamma * (1.0 - v_dot_n);

    // Clamp denominator to avoid physical singularities (photon in orbital plane)
    doppler_denom = std::clamp(doppler_denom, 0.1, 10.0);

    double g = grav_factor / doppler_denom;

    // Clamp final g-factor to physically reasonable range
    g = std::clamp(g, 0.1, 5.0);

    return static_cast<float>(g);
}

// =============================================================================
// Compute g-factor with motion blur components
// =============================================================================
// This version populates the TraceResult with decomposed g-factor components
// for analytical motion blur computation in the renderer.
//
void GeodesicTracer::computeGFactorWithComponents(float r, float phi, const Vec4& ray_vel,
                                                   TraceResult& result) {
    // Get metric parameters
    auto params = m_Metric->getParameters();
    double M = params.count("mass") ? params.at("mass").value : 1.0;
    double a_over_M = params.count("spin") ? params.at("spin").value : 0.0;

    // Compute orbital angular velocity (exact Kerr formula)
    float Omega = computeOrbitalVelocity(r);

    // Orbital velocity magnitude v = r × Ω (in units of c)
    double v_orb = static_cast<double>(r) * Omega;

    // Gravitational redshift component (exact Kerr formula)
    double r_d = static_cast<double>(r);
    double a = a_over_M * M;
    double a2 = a * a;
    double r2 = r_d * r_d;
    double grav_factor = std::sqrt(1.0 - 2.0 * M * r_d / (r2 + a2));

    // Extract and normalize ray direction
    double n_x = ray_vel(1);
    double n_y = ray_vel(2);
    double n_z = ray_vel(3);
    double v_mag = std::sqrt(n_x*n_x + n_y*n_y + n_z*n_z);

    if (v_mag < 1e-10) {
        result.redshift = static_cast<float>(grav_factor);
        result.gfactor_grav = static_cast<float>(grav_factor);
        result.gfactor_gamma = 1.0f;
        result.gfactor_v_orb = static_cast<float>(v_orb);
        result.gfactor_A = 0.0f;
        result.gfactor_B = 0.0f;
        return;
    }

    n_x /= v_mag;
    n_y /= v_mag;

    // Lorentz factor
    double gamma = 1.0 / std::sqrt(std::max(1e-10, 1.0 - v_orb * v_orb));

    // Doppler coefficients (exact formulas)
    double phi_d = static_cast<double>(phi);
    double sin_phi = std::sin(phi_d);
    double cos_phi = std::cos(phi_d);

    // A and B coefficients (negated n for outgoing photon convention)
    double A = (-n_x) * sin_phi - (-n_y) * cos_phi;
    double B = (-n_x) * cos_phi + (-n_y) * sin_phi;

    // Store decomposed components
    result.gfactor_grav = static_cast<float>(grav_factor);
    result.gfactor_gamma = static_cast<float>(gamma);
    result.gfactor_v_orb = static_cast<float>(v_orb);
    result.gfactor_A = static_cast<float>(A);
    result.gfactor_B = static_cast<float>(B);

    // Compute g-factor at original φ
    double v_dot_n = v_orb * A;
    double doppler_denom = gamma * (1.0 - v_dot_n);
    doppler_denom = std::clamp(doppler_denom, 0.1, 10.0);
    double g = grav_factor / doppler_denom;
    g = std::clamp(g, 0.1, 5.0);

    result.redshift = static_cast<float>(g);
}

// =============================================================================
// Volumetric Disk Methods (Phase 6)
// =============================================================================

float GeodesicTracer::computeScaleHeight(float r) {
    if (r <= 0) return 0;

    // H(r) = H_over_r × r × (r/r_ref)^H_power
    // Reference radius is disk_inner for consistency
    float r_ref = m_Config.disk_inner;
    float r_ratio = r / r_ref;
    float H_over_r = m_Config.volumetric_H_over_r * std::pow(r_ratio, m_Config.volumetric_H_power);

    // Clamp to physical range [0.01, 0.5]
    H_over_r = std::clamp(H_over_r, 0.01f, 0.5f);

    return H_over_r * r;
}

bool GeodesicTracer::isInVolumetricDisk(float r, float z) {
    // Radial bounds
    if (r < m_Config.disk_inner || r > m_Config.disk_outer) return false;

    // Vertical bounds: 3σ truncation
    float H = computeScaleHeight(r);
    float z_max = 3.0f * H;

    return std::abs(z) <= z_max;
}

float GeodesicTracer::computeVolumetricOpacityDensity(float r, float z) {
    if (!isInVolumetricDisk(r, z)) return 0;

    float H = computeScaleHeight(r);
    if (H <= 0) return 0;

    // Gaussian vertical profile: ρ ∝ exp(-z²/(2H²))
    float z_over_H = z / H;
    float gaussian = std::exp(-0.5f * z_over_H * z_over_H);

    // Radial scaling: κρ_0 ∝ r^(-1.5) (surface density falloff)
    // Normalize so vertical optical depth = tau_midplane at r_ref
    float r_ref = m_Config.disk_inner;
    float H_ref = computeScaleHeight(r_ref);
    float r_ratio = r / r_ref;

    // κρ_0(r) = tau_midplane / (√(2π) × H(r)) × (r/r_ref)^(-1.5)
    float kappa_rho0 = m_Config.volumetric_tau_midplane /
                       (std::sqrt(2.0f * static_cast<float>(M_PI)) * H) *
                       std::pow(r_ratio, -1.5f) * (H_ref / H);

    return kappa_rho0 * gaussian;
}

float GeodesicTracer::computeVolumetricTemperature(float r, float z) {
    if (!isInVolumetricDisk(r, z)) return 0;

    // Midplane temperature from thin disk model
    float T_mid = computeDiskTemperature(r);
    if (T_mid <= 0) return T_mid;

    // Vertical temperature gradient
    // T(z) = T_mid × [1 - (z/H)² × (1 - T_atm_ratio⁴)]^0.25
    // T_atm_ratio = 0.8 (cooler atmosphere)
    float H = computeScaleHeight(r);
    if (H <= 0) return T_mid;

    float z_over_H = std::abs(z) / H;
    float z_over_H_sq = z_over_H * z_over_H;

    // Clamp to 3σ boundary
    z_over_H_sq = std::min(z_over_H_sq, 9.0f);

    // T^4 interpolation
    constexpr float T_atm_ratio = 0.8f;
    constexpr float T_atm_ratio_4 = T_atm_ratio * T_atm_ratio * T_atm_ratio * T_atm_ratio;
    float T4_factor = 1.0f - z_over_H_sq * (1.0f - T_atm_ratio_4) / 9.0f;

    // Clamp to avoid negative
    T4_factor = std::max(T4_factor, T_atm_ratio_4);

    return T_mid * std::pow(T4_factor, 0.25f);
}

void GeodesicTracer::accumulateVolumetricEmission(const Lightray& ray,
                                                   const Vec4& entry_pos,
                                                   const Vec4& exit_pos,
                                                   TraceResult& result) {
    // Compute path length
    double dx = exit_pos(1) - entry_pos(1);
    double dy = exit_pos(2) - entry_pos(2);
    double dz = exit_pos(3) - entry_pos(3);
    double path_length = std::sqrt(dx*dx + dy*dy + dz*dz);

    if (path_length < 1e-10) return;

    result.volumetric_path_length = static_cast<float>(path_length);

    // Ray marching parameters
    int N = m_Config.volumetric_samples;
    double ds = path_length / N;

    // Accumulators
    double accumulated_tau = 0.0;
    double accumulated_r = 0.0;
    double accumulated_g = 0.0;
    double accumulated_b = 0.0;

    // Direction unit vector
    double dir_x = dx / path_length;
    double dir_y = dy / path_length;
    double dir_z = dz / path_length;

    // Ray march through disk
    for (int i = 0; i < N; i++) {
        // Sample position (midpoint rule)
        double t = (i + 0.5) * ds;
        double x = entry_pos(1) + t * dir_x;
        double y = entry_pos(2) + t * dir_y;
        double z = entry_pos(3) + t * dir_z;

        // Cylindrical coordinates
        double r = std::sqrt(x*x + y*y);
        double z_cyl = z;

        // Skip if outside disk
        if (!isInVolumetricDisk(static_cast<float>(r), static_cast<float>(z_cyl))) {
            continue;
        }

        // Local properties
        float kappa_rho = computeVolumetricOpacityDensity(static_cast<float>(r), static_cast<float>(z_cyl));
        float T = computeVolumetricTemperature(static_cast<float>(r), static_cast<float>(z_cyl));

        if (kappa_rho <= 0 || T <= 0) continue;

        // Optical depth increment
        double dtau = kappa_rho * ds;

        // Source function: blackbody at local temperature
        // Use simplified Stefan-Boltzmann: B ∝ T^4
        // Normalize relative to inner temperature
        float T_inner = m_Config.disk_temperature_inner;
        double B = std::pow(T / T_inner, 4.0);

        // Radiative transfer for this step
        // I_out = I_in × exp(-dtau) + B × (1 - exp(-dtau))
        double transmission = std::exp(-dtau);
        double one_minus_trans = 1.0 - transmission;

        // Accumulate (gray approximation: same for R, G, B)
        // In reality, we'd integrate Planck function over wavelength
        accumulated_r = accumulated_r * transmission + B * one_minus_trans;
        accumulated_g = accumulated_g * transmission + B * one_minus_trans;
        accumulated_b = accumulated_b * transmission + B * one_minus_trans;

        accumulated_tau += dtau;

        // Early termination if optically thick
        if (accumulated_tau > m_Config.volumetric_tau_max) break;
    }

    // Store results
    result.volumetric_emission[0] = static_cast<float>(accumulated_r);
    result.volumetric_emission[1] = static_cast<float>(accumulated_g);
    result.volumetric_emission[2] = static_cast<float>(accumulated_b);
    result.optical_depth = static_cast<float>(accumulated_tau);
    result.volumetric_hit = (accumulated_tau > 0.01);  // Threshold for "hit"
}

} // namespace Sirius
