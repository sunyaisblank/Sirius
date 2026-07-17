// Geodesic ray tracer implementation. Ported from GTRC001A.cpp.
//
// The arithmetic is preserved exactly from the legacy reference path; only the
// namespaces, symbol casing, and the core API spellings change. This is the
// reference tracer whose output the image-identity gate pins.

#include "sirius/backend/cpu/geodesic_tracer.h"

#include "sirius/core/disk/novikov_thorne_disk.h"  // AccretionDiskD::ComputeIsco (ISCO authority).

#include <algorithm>
#include <cmath>

namespace sirius::backend {

using namespace sirius::core;

// =============================================================================
// Initialise a Lightray from a camera ray.
// =============================================================================
Lightray GeodesicTracer::InitializeLightray(const CameraRay& camera_ray) {
    Lightray ray;

    // Camera ray origin is Boyer-Lindquist (t, r, theta, phi).
    double t = camera_ray.origin(0);
    double r = camera_ray.origin(1);
    double th = camera_ray.origin(2);
    double ph = camera_ray.origin(3);

    // Clamp theta away from the poles.
    th = std::clamp(th, 0.01, M_PI - 0.01);

    // Position: Boyer-Lindquist -> Kerr-Schild Cartesian (consolidated transforms).
    coordinates::Vec4Bl pos_bl(t, r, th, ph);
    coordinates::Vec4Cart pos_cart = coordinates::BlToCartesian(pos_bl);

    ray.position(0) = static_cast<float>(pos_cart.t);
    ray.position(1) = static_cast<float>(pos_cart.x);
    ray.position(2) = static_cast<float>(pos_cart.y);
    ray.position(3) = static_cast<float>(pos_cart.z);

    double sin_th = std::sin(th);
    double cos_th = std::cos(th);
    double sin_ph = std::sin(ph);
    double cos_ph = std::cos(ph);

    // Velocity: spherical basis -> Cartesian basis via the Jacobian
    //   v_x = v_r sin th cos ph + r v_th cos th cos ph - r v_ph sin th sin ph
    //   v_y = v_r sin th sin ph + r v_th cos th sin ph + r v_ph sin th cos ph
    //   v_z = v_r cos th        - r v_th sin th
    double v_r = camera_ray.direction(1);
    double v_th = camera_ray.direction(2) / r;             // Camera stores r dtheta/dlambda.
    double v_ph = camera_ray.direction(3) / (r * sin_th);  // Camera stores r sin(theta) dphi/dlambda.

    // Avoid division by zero at the poles.
    if (std::abs(sin_th) < 1e-10) {
        v_ph = 0.0;
    }

    ray.velocity(1) = static_cast<float>(v_r * sin_th * cos_ph + r * v_th * cos_th * cos_ph -
                                         r * v_ph * sin_th * sin_ph);
    ray.velocity(2) = static_cast<float>(v_r * sin_th * sin_ph + r * v_th * cos_th * sin_ph +
                                         r * v_ph * sin_th * cos_ph);
    ray.velocity(3) = static_cast<float>(v_r * cos_th - r * v_th * sin_th);

    // Normalise to the null condition g_mu_nu k^mu k^nu = 0; solve for k^0.
    Metric4d g;
    Tensor<Dual<double>, 4, 4, 4> dg;

    Vec4 pos_double;
    for (int i = 0; i < 4; ++i) pos_double(i) = ray.position(i);

    metric_->Evaluate(pos_double, g, dg);

    Vec4 vel_double;
    vel_double(0) = 0.0;  // Computed by the null normalisation.
    for (int i = 1; i < 4; ++i) vel_double(i) = ray.velocity(i);

    Vec4 normalized = TensorOps::NormalizeNull(vel_double, g);

    for (int i = 0; i < 4; ++i) {
        ray.velocity(i) = static_cast<float>(normalized(i));
    }

    // Initialise the remaining fields.
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

    ray.acceleration = Geodesic::CalculateAcceleration(vel_for_accel, pos_for_accel, metric_);

    ray.proper_time = 0.0f;
    ray.coordinate_time = 0.0f;
    ray.step_size = config_.integrator.initial_step;
    ray.terminated = 0;
    ray.bounce_count = 0;
    ray.ku_uobsu = 1.0f;  // Initial k.u for redshift.
    ray.running_dlambda_dnew = 1.0f;

    return ray;
}

// =============================================================================
// Main trace.
// =============================================================================
TraceResult GeodesicTracer::Trace(const CameraRay& camera_ray) {
    TraceResult result;
    result.steps_taken = 0;
    result.numerical_failure = false;

    CacheMetricParameters();

    Lightray ray = InitializeLightray(camera_ray);

    if (HasInvalidState(ray)) {
        result.outcome = TraceResult::Outcome::MaxSteps;
        result.numerical_failure = true;
        return result;
    }

    // Photon-sphere radius, exact Kerr formula:
    //   r_photon = 2M [1 + cos(2/3 arccos(-a/M))]   (prograde, a > 0).
    // Reference: Bardeen, Press & Teukolsky (1972).
    double M = cached_m_;
    double a_over_M = cached_a_;

    a_over_M = std::clamp(a_over_M, -0.998, 0.998);

    double arg = -a_over_M;  // Prograde orbit for positive spin.
    arg = std::clamp(arg, -1.0, 1.0);
    double r_photon_d = 2.0 * M * (1.0 + std::cos((2.0 / 3.0) * std::acos(arg)));
    float r_photon = static_cast<float>(r_photon_d);

    float min_r = 1e10f;

    // Higher-order imaging: winding-number tracking.
    float prev_phi = static_cast<float>(std::atan2(ray.position(2), ray.position(1)));
    float prev_z = static_cast<float>(ray.position(3));
    int equatorial_crossings = 0;
    float total_phi_change = 0.0f;

    // Impact parameter b = L/E from the initial conditions.
    {
        double x0 = ray.position(1);
        double y0 = ray.position(2);
        double z0 = ray.position(3);
        double r0 = std::sqrt(x0 * x0 + y0 * y0 + z0 * z0);

        double vx = ray.velocity(1);
        double vy = ray.velocity(2);
        double vz = ray.velocity(3);
        double v_mag = std::sqrt(vx * vx + vy * vy + vz * vz);

        double n_x = x0 / r0;
        double n_y = y0 / r0;
        double n_z = z0 / r0;

        double v_radial = (vx * n_x + vy * n_y + vz * n_z);
        double v_perp_sq = v_mag * v_mag - v_radial * v_radial;
        double v_perp = std::sqrt(std::max(v_perp_sq, 0.0));

        result.impact_parameter = static_cast<float>(r0 * v_perp / std::max(v_mag, 1e-10));
    }

    Vec4 prev_pos;
    for (int i = 0; i < 4; ++i) prev_pos(i) = ray.position(i);

    // Spiral-orbit early termination: rays near b_crit orbit the photon sphere
    // many times; the higher-order images they produce are exponentially dim
    // (~exp(-pi n)), so they are cut once quasi-circular motion is detected.
    constexpr float SPIRAL_PHI_THRESHOLD = 2.0f * static_cast<float>(M_PI);
    constexpr float SPIRAL_R_THRESHOLD = 1.5f;
    constexpr int SPIRAL_CHECK_INTERVAL = 100;

    for (int step = 0; step < config_.max_steps; ++step) {
        for (int i = 0; i < 4; ++i) prev_pos(i) = ray.position(i);

        bool success = Geodesic::IntegrateStepRk45(ray, metric_, config_.integrator);
        result.steps_taken++;

        if (!success || ray.terminated || HasInvalidState(ray)) {
            result.outcome = TraceResult::Outcome::MaxSteps;
            result.numerical_failure = HasInvalidState(ray);
            break;
        }

        double x = ray.position(1);
        double y = ray.position(2);
        double z = ray.position(3);
        double r = std::sqrt(x * x + y * y + z * z);

        if (r < min_r) {
            min_r = static_cast<float>(r);
        }

        // Higher-order imaging: track winding.
        float curr_phi = static_cast<float>(std::atan2(y, x));
        float curr_z = static_cast<float>(z);

        if (prev_z * curr_z < 0) {
            equatorial_crossings++;
        }

        float dphi = curr_phi - prev_phi;
        if (dphi > M_PI) dphi -= 2.0f * static_cast<float>(M_PI);
        if (dphi < -M_PI) dphi += 2.0f * static_cast<float>(M_PI);
        total_phi_change += std::abs(dphi);

        prev_phi = curr_phi;
        prev_z = curr_z;

        // Spiral-orbit early termination.
        if ((step % SPIRAL_CHECK_INTERVAL == 0) && step > SPIRAL_CHECK_INTERVAL) {
            float r_ratio = min_r / r_photon;
            bool near_photon_sphere = (r_ratio < SPIRAL_R_THRESHOLD) && (r_ratio > 0.9f);
            bool has_spiraled = (total_phi_change > SPIRAL_PHI_THRESHOLD);

            if (near_photon_sphere && has_spiraled) {
                double vx = ray.velocity(1);
                double vy = ray.velocity(2);
                double vz = ray.velocity(3);
                double v_radial = (x * vx + y * vy + z * vz) / r;
                double v_total = std::sqrt(vx * vx + vy * vy + vz * vz);
                double radial_fraction = std::abs(v_radial) / std::max(v_total, 1e-10);

                if (radial_fraction < 0.3) {
                    if (result.num_disk_crossings > 0) {
                        result.min_radius = min_r;
                    } else {
                        result.outcome = TraceResult::Outcome::Spiraling;
                    }
                    result.equatorial_crossings = equatorial_crossings;
                    result.total_phi_change = total_phi_change;
                    result.image_order = static_cast<int>(total_phi_change / M_PI);
                    break;
                }
            }
        }

        // 1. Horizon capture: the metric decides in its own coordinates, so the
        // oblate Kerr horizon is placed exactly even at high spin.
        if (metric_->InsideCaptureSurface(ray.position, config_.horizon_factor - 1.0)) {
            result.outcome = TraceResult::Outcome::Horizon;
            result.final_position(0) = ray.position(0);
            result.final_position(1) = ray.position(1);
            result.final_position(2) = ray.position(2);
            result.final_position(3) = ray.position(3);
            result.min_radius = min_r;
            break;
        }

        // 2. Escape to infinity.
        if (r > config_.escape_radius) {
            double vx = ray.velocity(1);
            double vy = ray.velocity(2);
            double vz = ray.velocity(3);
            double v_radial = (x * vx + y * vy + z * vz) / r;

            if (v_radial > 0) {
                // A disk hit found earlier keeps priority for rendering.
                if (result.num_disk_crossings == 0) {
                    result.outcome = TraceResult::Outcome::Escaped;
                }

                result.final_position(0) = ray.position(0);
                result.final_position(1) = ray.position(1);
                result.final_position(2) = ray.position(2);
                result.final_position(3) = ray.position(3);
                result.final_direction(0) = ray.velocity(0);
                result.final_direction(1) = ray.velocity(1);
                result.final_direction(2) = ray.velocity(2);
                result.final_direction(3) = ray.velocity(3);
                result.min_radius = min_r;

                // Photon-ring detection and magnification (Luminet 1979):
                //   mu = 1 + 2 / (r_min/r_photon - 1), capped at 20.
                float r_ratio = min_r / r_photon;
                if (r_ratio < 1.5f) {
                    result.photon_ring = true;

                    if (r_ratio > 1.0f) {
                        float denom = r_ratio - 1.0f;
                        result.magnification = 1.0f + 2.0f / std::max(denom, 0.1f);
                    } else {
                        // Ray penetrated the photon sphere: numerical artefact
                        // for an escaping ray; apply the maximum magnification.
                        result.magnification = 20.0f;
                    }

                    result.magnification = std::min(result.magnification, 20.0f);
                }
                break;
            }
        }

        // 3. Volumetric disk (ray marching through a 3D disk volume).
        if (config_.enable_disk && config_.enable_volumetric) {
            float prev_r = static_cast<float>(
                std::sqrt(prev_pos(1) * prev_pos(1) + prev_pos(2) * prev_pos(2)));
            float prev_z_cyl = static_cast<float>(prev_pos(3));
            float curr_r = static_cast<float>(std::sqrt(x * x + y * y));
            float curr_z_cyl = static_cast<float>(z);

            bool was_inside = IsInVolumetricDisk(prev_r, prev_z_cyl);
            bool is_inside = IsInVolumetricDisk(curr_r, curr_z_cyl);

            if (was_inside || is_inside) {
                Vec4 segment_start, segment_end;
                for (int i = 0; i < 4; ++i) {
                    segment_start(i) = prev_pos(i);
                    segment_end(i) = ray.position(i);
                }

                AccumulateVolumetricEmission(ray, segment_start, segment_end, result);

                if (result.volumetric_hit && result.outcome != TraceResult::Outcome::DiskHit) {
                    result.outcome = TraceResult::Outcome::DiskHit;
                    result.disk_radius = curr_r;
                    result.disk_phi = static_cast<float>(std::atan2(y, x));
                }
            }
        }

        // 4. Thin disk (accumulate multiple crossings for higher-order imaging).
        if (config_.enable_disk && !config_.enable_volumetric) {
            Vec4 curr_pos;
            curr_pos(0) = ray.position(0);
            curr_pos(1) = ray.position(1);
            curr_pos(2) = ray.position(2);
            curr_pos(3) = ray.position(3);
            float disk_r, disk_phi;

            if (CheckDiskIntersection(prev_pos, curr_pos, disk_r, disk_phi)) {
                // Each crossing is a different image order: crossing 0 primary,
                // 1 secondary, and so on; the first sets the primary outcome.
                if (result.num_disk_crossings < TraceResult::kMaxDiskCrossings) {
                    auto& crossing = result.disk_crossings[result.num_disk_crossings];
                    crossing.r = disk_r;
                    crossing.phi = disk_phi;
                    crossing.temperature = ComputeDiskTemperature(disk_r);
                    crossing.crossing_index = result.num_disk_crossings;
                    crossing.valid = true;

                    Vec4 ray_vel;
                    ray_vel(0) = ray.velocity(0);
                    ray_vel(1) = ray.velocity(1);
                    ray_vel(2) = ray.velocity(2);
                    ray_vel(3) = ray.velocity(3);

                    if (result.num_disk_crossings == 0) {
                        // First crossing: full computation with motion-blur terms.
                        ComputeGFactorWithComponents(disk_r, disk_phi, ray_vel, result);
                        crossing.redshift = result.redshift;

                        result.outcome = TraceResult::Outcome::DiskHit;
                        result.disk_radius = disk_r;
                        result.disk_phi = disk_phi;
                        result.disk_temperature = crossing.temperature;
                    } else {
                        crossing.redshift = ComputeGFactor(disk_r, disk_phi, ray_vel);
                    }

                    result.num_disk_crossings++;
                }

                // Magnification for disk hits (closest approach to photon sphere).
                float r_ratio_disk = min_r / r_photon;
                if (r_ratio_disk < 1.5f) {
                    result.photon_ring = true;

                    if (r_ratio_disk > 1.0f) {
                        float denom = r_ratio_disk - 1.0f;
                        result.magnification = 1.0f + 2.0f / std::max(denom, 0.05f);
                    } else {
                        result.magnification = 100.0f;  // Higher cap for caustics.
                    }

                    result.magnification = std::min(result.magnification, 100.0f);
                }

                if (result.num_disk_crossings >= TraceResult::kMaxDiskCrossings) {
                    result.min_radius = min_r;
                    break;
                }
            }
        }
    }

    result.min_radius = std::min(result.min_radius, min_r);

    // Finalise winding data. Image order from equatorial crossings (each ~half
    // orbit): n = ceil(crossings / 2) for crossings >= 1, else 0.
    result.equatorial_crossings = equatorial_crossings;
    result.total_phi_change = total_phi_change;

    if (equatorial_crossings >= 1) {
        result.image_order = (equatorial_crossings + 1) / 2;
    } else {
        result.image_order = 0;
    }

    // No termination reached: result stays MaxSteps (default).
    return result;
}

// =============================================================================
// Trace with per-call metric parameters.
// =============================================================================
TraceResult GeodesicTracer::Trace(const CameraRay& camera_ray, double mass, double spin) {
    metric_->SetParameter("mass", mass);
    metric_->SetParameter("spin", spin / mass);  // Store as a/M.
    CacheMetricParameters();

    // Update the disk inner edge to the ISCO for the new spin; the Bardeen
    // formula lives in one place (AccretionDiskD::ComputeIsco, units of M).
    config_.disk_inner = static_cast<float>(AccretionDiskD::ComputeIsco(spin / mass) * mass);

    return Trace(camera_ray);
}

// =============================================================================
// Keplerian orbital velocity.
// =============================================================================
float GeodesicTracer::ComputeOrbitalVelocity(float r) {
    double M = cached_m_;
    double a = cached_a_ * M;

    // Prograde Kerr angular velocity Omega = sqrt(M) / (r^(3/2) + a sqrt(M)).
    double r32 = std::pow(static_cast<double>(r), 1.5);
    double sqrtM = std::sqrt(M);

    double Omega = sqrtM / (r32 + a * sqrtM);

    return static_cast<float>(Omega);
}

// =============================================================================
// g-factor (gravitational + Doppler redshift).
// =============================================================================
// g = grav_factor / (gamma (1 - v.n)). The Doppler term is decomposed so the
// renderer can evaluate g at an azimuthal offset dphi analytically:
//   v.n(phi + dphi) = v_orb (A cos dphi + B sin dphi),
//   A = n_x sin phi - n_y cos phi,  B = n_x cos phi + n_y sin phi.
// This is exact; only the final clamp guards physical singularities.
float GeodesicTracer::ComputeGFactor(float r, float phi, const Vec4& ray_vel) {
    double M = cached_m_;
    double a_over_M = cached_a_;

    float Omega = ComputeOrbitalVelocity(r);
    double v_orb = static_cast<double>(r) * Omega;

    // Gravitational redshift sqrt(1 - 2Mr / (r^2 + a^2)).
    double r_d = static_cast<double>(r);
    double a = a_over_M * M;
    double a2 = a * a;
    double r2 = r_d * r_d;
    double grav_factor = std::sqrt(1.0 - 2.0 * M * r_d / (r2 + a2));

    double n_x = ray_vel(1);
    double n_y = ray_vel(2);
    double n_z = ray_vel(3);
    double v_mag = std::sqrt(n_x * n_x + n_y * n_y + n_z * n_z);

    if (v_mag < 1e-10) {
        return static_cast<float>(grav_factor);
    }

    n_x /= v_mag;
    n_y /= v_mag;
    // n_z is unused for the equatorial disk.

    double gamma = 1.0 / std::sqrt(std::max(1e-10, 1.0 - v_orb * v_orb));

    double phi_d = static_cast<double>(phi);
    double sin_phi = std::sin(phi_d);
    double cos_phi = std::cos(phi_d);

    // n points toward the observer, so negate for the outgoing-photon convention.
    double A = (-n_x) * sin_phi - (-n_y) * cos_phi;

    double v_dot_n = v_orb * A;

    double doppler_denom = gamma * (1.0 - v_dot_n);
    doppler_denom = std::clamp(doppler_denom, 0.1, 10.0);

    double g = grav_factor / doppler_denom;
    g = std::clamp(g, 0.1, 5.0);

    return static_cast<float>(g);
}

// =============================================================================
// g-factor with the motion-blur components stored in the result.
// =============================================================================
void GeodesicTracer::ComputeGFactorWithComponents(float r, float phi, const Vec4& ray_vel,
                                                  TraceResult& result) {
    double M = cached_m_;
    double a_over_M = cached_a_;

    float Omega = ComputeOrbitalVelocity(r);
    double v_orb = static_cast<double>(r) * Omega;

    double r_d = static_cast<double>(r);
    double a = a_over_M * M;
    double a2 = a * a;
    double r2 = r_d * r_d;
    double grav_factor = std::sqrt(1.0 - 2.0 * M * r_d / (r2 + a2));

    double n_x = ray_vel(1);
    double n_y = ray_vel(2);
    double n_z = ray_vel(3);
    double v_mag = std::sqrt(n_x * n_x + n_y * n_y + n_z * n_z);

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

    double gamma = 1.0 / std::sqrt(std::max(1e-10, 1.0 - v_orb * v_orb));

    double phi_d = static_cast<double>(phi);
    double sin_phi = std::sin(phi_d);
    double cos_phi = std::cos(phi_d);

    double A = (-n_x) * sin_phi - (-n_y) * cos_phi;
    double B = (-n_x) * cos_phi + (-n_y) * sin_phi;

    result.gfactor_grav = static_cast<float>(grav_factor);
    result.gfactor_gamma = static_cast<float>(gamma);
    result.gfactor_v_orb = static_cast<float>(v_orb);
    result.gfactor_A = static_cast<float>(A);
    result.gfactor_B = static_cast<float>(B);

    double v_dot_n = v_orb * A;
    double doppler_denom = gamma * (1.0 - v_dot_n);
    doppler_denom = std::clamp(doppler_denom, 0.1, 10.0);
    double g = grav_factor / doppler_denom;
    g = std::clamp(g, 0.1, 5.0);

    result.redshift = static_cast<float>(g);
}

// =============================================================================
// Volumetric disk.
// =============================================================================
float GeodesicTracer::ComputeScaleHeight(float r) {
    if (r <= 0) return 0;

    // H(r) = H_over_r r (r/r_ref)^H_power, referenced at the inner edge.
    float r_ref = config_.disk_inner;
    float r_ratio = r / r_ref;
    float H_over_r = config_.volumetric_H_over_r * std::pow(r_ratio, config_.volumetric_H_power);

    H_over_r = std::clamp(H_over_r, 0.01f, 0.5f);

    return H_over_r * r;
}

bool GeodesicTracer::IsInVolumetricDisk(float r, float z) {
    if (r < config_.disk_inner || r > config_.disk_outer) return false;

    // Vertical bounds: 3-sigma truncation.
    float H = ComputeScaleHeight(r);
    float z_max = 3.0f * H;

    return std::abs(z) <= z_max;
}

float GeodesicTracer::ComputeVolumetricOpacityDensity(float r, float z) {
    if (!IsInVolumetricDisk(r, z)) return 0;

    float H = ComputeScaleHeight(r);
    if (H <= 0) return 0;

    // Gaussian vertical profile rho ~ exp(-z^2 / (2 H^2)).
    float z_over_H = z / H;
    float gaussian = std::exp(-0.5f * z_over_H * z_over_H);

    // Radial scaling normalised so the vertical optical depth is tau_midplane
    // at r_ref: kappa rho_0(r) = tau_midplane / (sqrt(2 pi) H) (r/r_ref)^(-1.5).
    float r_ref = config_.disk_inner;
    float H_ref = ComputeScaleHeight(r_ref);
    float r_ratio = r / r_ref;

    float kappa_rho0 = config_.volumetric_tau_midplane /
                       (std::sqrt(2.0f * static_cast<float>(M_PI)) * H) *
                       std::pow(r_ratio, -1.5f) * (H_ref / H);

    return kappa_rho0 * gaussian;
}

float GeodesicTracer::ComputeVolumetricTemperature(float r, float z) {
    if (!IsInVolumetricDisk(r, z)) return 0;

    float T_mid = ComputeDiskTemperature(r);
    if (T_mid <= 0) return T_mid;

    // Vertical temperature gradient with a cooler atmosphere (T_atm_ratio = 0.8):
    //   T(z) = T_mid [1 - (z/H)^2 (1 - T_atm_ratio^4)]^(1/4).
    float H = ComputeScaleHeight(r);
    if (H <= 0) return T_mid;

    float z_over_H = std::abs(z) / H;
    float z_over_H_sq = z_over_H * z_over_H;

    z_over_H_sq = std::min(z_over_H_sq, 9.0f);

    constexpr float T_atm_ratio = 0.8f;
    constexpr float T_atm_ratio_4 = T_atm_ratio * T_atm_ratio * T_atm_ratio * T_atm_ratio;
    float T4_factor = 1.0f - z_over_H_sq * (1.0f - T_atm_ratio_4) / 9.0f;

    T4_factor = std::max(T4_factor, T_atm_ratio_4);

    return T_mid * std::pow(T4_factor, 0.25f);
}

void GeodesicTracer::AccumulateVolumetricEmission([[maybe_unused]] const Lightray& ray,
                                                  const Vec4& entry_pos, const Vec4& exit_pos,
                                                  TraceResult& result) {
    // ray is retained in the signature for parity with the legacy interface; the
    // segment endpoints carry everything the ray march needs.
    double dx = exit_pos(1) - entry_pos(1);
    double dy = exit_pos(2) - entry_pos(2);
    double dz = exit_pos(3) - entry_pos(3);
    double path_length = std::sqrt(dx * dx + dy * dy + dz * dz);

    if (path_length < 1e-10) return;

    result.volumetric_path_length = static_cast<float>(path_length);

    int N = config_.volumetric_samples;
    double ds = path_length / N;

    double accumulated_tau = 0.0;
    double accumulated_r = 0.0;
    double accumulated_g = 0.0;
    double accumulated_b = 0.0;

    double dir_x = dx / path_length;
    double dir_y = dy / path_length;
    double dir_z = dz / path_length;

    for (int i = 0; i < N; i++) {
        // Midpoint sample.
        double t = (i + 0.5) * ds;
        double x = entry_pos(1) + t * dir_x;
        double y = entry_pos(2) + t * dir_y;
        double z = entry_pos(3) + t * dir_z;

        double r = std::sqrt(x * x + y * y);
        double z_cyl = z;

        if (!IsInVolumetricDisk(static_cast<float>(r), static_cast<float>(z_cyl))) {
            continue;
        }

        float kappa_rho =
            ComputeVolumetricOpacityDensity(static_cast<float>(r), static_cast<float>(z_cyl));
        float T = ComputeVolumetricTemperature(static_cast<float>(r), static_cast<float>(z_cyl));

        if (kappa_rho <= 0 || T <= 0) continue;

        double dtau = kappa_rho * ds;

        // Grey source function, Stefan-Boltzmann normalised to the inner edge.
        float T_inner = config_.disk_temperature_inner;
        double B = std::pow(T / T_inner, 4.0);

        // Radiative transfer I_out = I_in exp(-dtau) + B (1 - exp(-dtau)).
        double transmission = std::exp(-dtau);
        double one_minus_trans = 1.0 - transmission;

        accumulated_r = accumulated_r * transmission + B * one_minus_trans;
        accumulated_g = accumulated_g * transmission + B * one_minus_trans;
        accumulated_b = accumulated_b * transmission + B * one_minus_trans;

        accumulated_tau += dtau;

        if (accumulated_tau > config_.volumetric_tau_max) break;
    }

    result.volumetric_emission[0] = static_cast<float>(accumulated_r);
    result.volumetric_emission[1] = static_cast<float>(accumulated_g);
    result.volumetric_emission[2] = static_cast<float>(accumulated_b);
    result.optical_depth = static_cast<float>(accumulated_tau);
    result.volumetric_hit = (accumulated_tau > 0.01);  // "Hit" threshold.
}

}  // namespace sirius::backend
