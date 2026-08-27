// Geodesic ray tracer implementation. Ported from GTRC001A.cpp.
//
// The arithmetic is preserved exactly from the legacy reference path; only the
// namespaces, symbol casing, and the core API spellings change. This is the
// reference tracer whose output the image-identity gate pins.

#include "sirius/backend/cpu/geodesic_tracer.h"

#include "sirius/core/disk/novikov_thorne_disk.h"  // AccretionDiskD::ComputeIsco (ISCO authority).
#include "sirius/core/spectral/colour_modes.h"

#include <algorithm>
#include <cmath>
#include <numbers>

namespace sirius::backend {

using namespace sirius::core;

namespace {

void TransverseBasis(double nx, double ny, double nz, double& ax, double& ay, double& az,
                     double& bx, double& by, double& bz);

}  // namespace

// =============================================================================
// Initialise a Lightray from a camera ray.
// =============================================================================
Lightray GeodesicTracer::InitializeLightray(const CameraRay& camera_ray) {
    // Value-initialize the full device-facing record so screen coordinates and
    // alignment padding never carry indeterminate bytes into copies or hashes.
    Lightray ray{};

    // Camera ray origin is Boyer-Lindquist (t, r, theta, phi).
    double t = camera_ray.origin(0);
    double r = camera_ray.origin(1);
    double th = camera_ray.origin(2);
    double ph = camera_ray.origin(3);

    // Clamp theta away from the poles.
    th = std::clamp(th, 0.01, std::numbers::pi - 0.01);

    // Position: Boyer-Lindquist -> Kerr-Schild Cartesian using the spin-aware
    // oblate transform.
    coordinates::Vec4Bl pos_bl(t, r, th, ph);
    const double absolute_spin = cached_a_ * cached_m_;
    coordinates::Vec4Cart pos_cart = coordinates::BlToKerrSchildCart(pos_bl, absolute_spin);

    ray.position(0) = static_cast<float>(pos_cart.t);
    ray.position(1) = static_cast<float>(pos_cart.x);
    ray.position(2) = static_cast<float>(pos_cart.y);
    ray.position(3) = static_cast<float>(pos_cart.z);

    double sin_th = std::sin(th);
    // Direction components are an orthonormal camera-local triad, not BL
    // coordinate derivatives. Rotate that triad directly into the asymptotic
    // Cartesian frame at the actual oblate position.
    const double cartesian_phi = std::atan2(pos_cart.y, pos_cart.x);
    const double sin_ph = std::sin(cartesian_phi);
    const double cos_ph = std::cos(cartesian_phi);
    const double v_r = camera_ray.direction(1);
    const double v_theta = camera_ray.direction(2);
    const double v_phi = camera_ray.direction(3);
    ray.velocity(1) = static_cast<float>(v_r * sin_th * cos_ph + v_theta * std::cos(th) * cos_ph -
                                         v_phi * sin_ph);
    ray.velocity(2) = static_cast<float>(v_r * sin_th * sin_ph + v_theta * std::cos(th) * sin_ph +
                                         v_phi * cos_ph);
    ray.velocity(3) = static_cast<float>(v_r * std::cos(th) - v_theta * sin_th);

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

void GeodesicTracer::InitPolarisationFrame(const Lightray& ray, PolarisationFrame& frame) {
    Vec4 position = ray.position;
    Vec4 velocity = ray.velocity;

    double nx = velocity(1);
    double ny = velocity(2);
    double nz = velocity(3);
    double length = std::sqrt(nx * nx + ny * ny + nz * nz);
    SIRIUS_PRE(length > 0.0);
    nx /= length;
    ny /= length;
    nz /= length;

    double ax, ay, az, bx, by, bz;
    TransverseBasis(nx, ny, nz, ax, ay, az, bx, by, bz);

    Vec4 reference_trial;
    reference_trial(1) = ax;
    reference_trial(2) = ay;
    reference_trial(3) = az;
    Vec4 perpendicular_trial;
    perpendicular_trial(1) = bx;
    perpendicular_trial(2) = by;
    perpendicular_trial(3) = bz;

    frame.reference.position = position;
    frame.reference.velocity = velocity;
    frame.reference.polarisation =
        MakeOrthonormalPolarisation(*metric_, position, velocity, reference_trial);
    frame.reference.affine = 0.0;

    frame.perpendicular.position = position;
    frame.perpendicular.velocity = velocity;
    frame.perpendicular.polarisation =
        MakeOrthonormalPolarisation(*metric_, position, velocity, perpendicular_trial);
    frame.perpendicular.affine = 0.0;

    ReconditionPolarisationFrame(frame, position, velocity);
}

void GeodesicTracer::AdvancePolarisationFrame(PolarisationFrame& frame, double d_lambda) {
    SIRIUS_PRE(d_lambda > 0.0);
    ParallelTransportStep(*metric_, frame.reference, d_lambda);
    ParallelTransportStep(*metric_, frame.perpendicular, d_lambda);
}

void GeodesicTracer::ReconditionPolarisationFrame(PolarisationFrame& frame, const Vec4& position,
                                                  const Vec4& velocity) {
    // Reset the carrier path to the accepted RK45 endpoint and remove only
    // numerical gauge/norm drift from the transported screen vectors.
    frame.reference.position = position;
    frame.reference.velocity = velocity;
    frame.reference.polarisation =
        MakeOrthonormalPolarisation(*metric_, position, velocity, frame.reference.polarisation);

    Metric4d g;
    Tensor<Dual<double>, 4, 4, 4> dg;
    metric_->Evaluate(position, g, dg);

    frame.perpendicular.position = position;
    frame.perpendicular.velocity = velocity;
    Vec4 perpendicular =
        MakeOrthonormalPolarisation(*metric_, position, velocity, frame.perpendicular.polarisation);
    perpendicular -= frame.reference.polarisation *
                     TensorOps::InnerProduct(perpendicular, frame.reference.polarisation, g);
    const double norm = TensorOps::InnerProduct(perpendicular, perpendicular, g);
    SIRIUS_PRE(norm > 0.0);
    frame.perpendicular.polarisation = perpendicular / std::sqrt(norm);
}

void GeodesicTracer::SetDiskPolarisation(const PolarisationFrame& frame,
                                         TraceResult::DiskCrossing& crossing) {
    const Vec4& position = frame.reference.position;
    const Vec4& wave = frame.reference.velocity;

    Metric4d g;
    Tensor<Dual<double>, 4, 4, 4> dg;
    metric_->Evaluate(position, g, dg);

    // Equatorial Page-Thorne emitter: u is the normalised circular worldline
    // with angular velocity Omega. At z=0 the Cartesian Kerr-Schild azimuth has
    // the same fixed-r derivative as Boyer-Lindquist phi.
    const double x = position(1);
    const double y = position(2);
    const coordinates::Vec4Cart crossing_cart{position(0), x, y, position(3)};
    const double orbital_radius =
        coordinates::KerrSchildRadius(crossing_cart, cached_a_ * cached_m_);
    const double omega = ComputeOrbitalVelocity(static_cast<float>(orbital_radius));
    Vec4 emitter;
    emitter(0) = 1.0;
    emitter(1) = -omega * y;
    emitter(2) = omega * x;
    const double emitter_norm = TensorOps::InnerProduct(emitter, emitter, g);
    SIRIUS_PRE(emitter_norm < 0.0);
    emitter = emitter / std::sqrt(-emitter_norm);

    const double emitted_frequency = -TensorOps::InnerProduct(emitter, wave, g);
    SIRIUS_PRE(emitted_frequency > 0.0);
    const Vec4 photon_direction = wave / emitted_frequency - emitter;

    Vec4 disk_normal;
    disk_normal(3) = 1.0;
    const double normal_norm = TensorOps::InnerProduct(disk_normal, disk_normal, g);
    SIRIUS_PRE(normal_norm > 0.0);
    disk_normal = disk_normal / std::sqrt(normal_norm);

    const float inclination_cosine = static_cast<float>(
        std::clamp(std::abs(TensorOps::InnerProduct(photon_direction, disk_normal, g)), 0.0, 1.0));
    crossing.polarisation_degree =
        polarised_emission::ThomsonPolarisationDegree(inclination_cosine);

    // The electric vector of the single-scattering atmosphere is parallel to
    // the disk plane, hence perpendicular to the projected disk normal. Gauge
    // additions proportional to k vanish against the transported screen basis.
    const Vec4 meridian = MakeOrthonormalPolarisation(*metric_, position, wave, disk_normal);
    const double meridian_reference =
        TensorOps::InnerProduct(meridian, frame.reference.polarisation, g);
    const double meridian_perpendicular =
        TensorOps::InnerProduct(meridian, frame.perpendicular.polarisation, g);
    const double meridian_angle = std::atan2(meridian_perpendicular, meridian_reference);
    crossing.polarisation_evpa = static_cast<float>(
        std::remainder(meridian_angle + std::numbers::pi / 2.0, std::numbers::pi));
    crossing.polarisation_valid = true;
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

    // Ray-bundle state (P2); propagated only when enabled, so the point-sampled
    // path is untouched. prev_vel and prev_pt carry the pre-step velocity and
    // affine parameter the deviation step needs.
    RayBundle bundle;
    Vec4 prev_vel;
    double prev_pt = 0.0;
    if (config_.enable_ray_bundles) {
        InitBundle(ray.velocity, bundle);
    }
    PolarisationFrame polarisation_frame;
    PolarisationFrame previous_polarisation_frame;
    if (config_.enable_polarisation) {
        InitPolarisationFrame(ray, polarisation_frame);
    }

    // Spiral-orbit early termination: rays near b_crit orbit the photon sphere
    // many times; the higher-order images they produce are exponentially dim
    // (~exp(-pi n)), so they are cut once quasi-circular motion is detected.
    constexpr float kSpiralPhiThreshold = 2.0f * static_cast<float>(std::numbers::pi);
    constexpr float kSpiralRadiusThreshold = 1.5f;
    constexpr int kSpiralCheckInterval = 100;

    for (int step = 0; step < config_.max_steps; ++step) {
        for (int i = 0; i < 4; ++i) prev_pos(i) = ray.position(i);
        if (config_.enable_ray_bundles || config_.enable_polarisation) {
            prev_vel = ray.velocity;
            prev_pt = static_cast<double>(ray.proper_time);
        }
        if (config_.enable_polarisation) {
            previous_polarisation_frame = polarisation_frame;
        }

        const double pre_step_radius =
            std::sqrt(static_cast<double>(ray.position(1)) * ray.position(1) +
                      static_cast<double>(ray.position(2)) * ray.position(2) +
                      static_cast<double>(ray.position(3)) * ray.position(3));
        auto step_config = config_.integrator;
        if (config_.strong_field_radius > 0.0f && config_.strong_field_max_step > 0.0f &&
            pre_step_radius < config_.strong_field_radius) {
            step_config.max_step = std::min(step_config.max_step, config_.strong_field_max_step);
        }
        bool success = Geodesic::IntegrateStepRk45(ray, metric_, step_config);
        result.steps_taken++;

        if (ray.terminated || HasInvalidState(ray)) {
            result.outcome = TraceResult::Outcome::MaxSteps;
            result.numerical_failure = HasInvalidState(ray);
            break;
        }
        if (!success) {
            // Error control rejected the step: the state is unchanged and
            // step_size was reduced, so retry. The attempt still counts
            // against max_steps, which bounds the work a stiff region can
            // consume. Treating rejection as termination here used to kill
            // rays at their first rejection with outcome MaxSteps (shaded
            // near-black), which no Kerr-Schild session scene triggered but
            // any metric with a steeper derivative profile did.
            continue;
        }

        const double d_lambda = static_cast<double>(ray.proper_time) - prev_pt;

        // Advance the ray bundle over the affine step just taken, with the
        // connection and curvature frozen at the step's start point.
        if (config_.enable_ray_bundles) {
            if (d_lambda > 0.0) {
                StepBundle(prev_pos, prev_vel, d_lambda, bundle);
            }
        }
        if (config_.enable_polarisation && d_lambda > 0.0) {
            AdvancePolarisationFrame(polarisation_frame, d_lambda);
            ReconditionPolarisationFrame(polarisation_frame, ray.position, ray.velocity);
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
        if (dphi > std::numbers::pi) dphi -= 2.0f * static_cast<float>(std::numbers::pi);
        if (dphi < -std::numbers::pi) dphi += 2.0f * static_cast<float>(std::numbers::pi);
        total_phi_change += std::abs(dphi);

        prev_phi = curr_phi;
        prev_z = curr_z;

        // Spiral-orbit early termination.
        if (config_.enable_spiral_termination && (step % kSpiralCheckInterval == 0) &&
            step > kSpiralCheckInterval) {
            float r_ratio = min_r / r_photon;
            bool near_photon_sphere = (r_ratio < kSpiralRadiusThreshold) && (r_ratio > 0.9f);
            bool has_spiraled = (total_phi_change > kSpiralPhiThreshold);

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
                    result.image_order = static_cast<int>(total_phi_change / std::numbers::pi);
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
        if ((config_.enable_disk && config_.enable_volumetric) || config_.enable_corona) {
            const double absolute_spin = cached_a_ * cached_m_;
            const coordinates::Vec4Cart previous_cart{prev_pos(0), prev_pos(1), prev_pos(2),
                                                      prev_pos(3)};
            const coordinates::Vec4Cart current_cart{ray.position(0), x, y, z};
            float prev_r =
                static_cast<float>(coordinates::KerrSchildRadius(previous_cart, absolute_spin));
            float prev_z_cyl = static_cast<float>(prev_pos(3));
            float curr_r =
                static_cast<float>(coordinates::KerrSchildRadius(current_cart, absolute_spin));
            float curr_z_cyl = static_cast<float>(z);

            bool was_inside = IsInVolumetricDisk(prev_r, prev_z_cyl);
            bool is_inside = IsInVolumetricDisk(curr_r, curr_z_cyl);

            if (was_inside || is_inside || config_.enable_corona) {
                Vec4 segment_start, segment_end;
                for (int i = 0; i < 4; ++i) {
                    segment_start(i) = prev_pos(i);
                    segment_end(i) = ray.position(i);
                }

                AccumulateVolumetricEmission(ray, segment_start, segment_end, result);
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

                    if (config_.enable_polarisation) {
                        const double denominator = ray.position(3) - prev_pos(3);
                        SIRIUS_ASSERT(std::abs(denominator) > 0.0);
                        const double crossing_fraction =
                            std::clamp(-prev_pos(3) / denominator, 0.0, 1.0);
                        PolarisationFrame crossing_frame = previous_polarisation_frame;
                        if (crossing_fraction > 0.0 && d_lambda > 0.0) {
                            AdvancePolarisationFrame(crossing_frame, crossing_fraction * d_lambda);
                        }
                        SetDiskPolarisation(crossing_frame, crossing);
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

    // Beam ellipse from the propagated bundle (P2), evaluated against the ray's
    // terminal direction. Only when ray bundles are enabled.
    if (config_.enable_ray_bundles) {
        FinaliseBundle(bundle, ray.velocity, static_cast<double>(ray.proper_time), result.beam);
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
    double grav_factor = std::sqrt(std::max(0.0, 1.0 - 2.0 * M * r_d / (r2 + a2)));

    double n_x = ray_vel(1);
    double n_y = ray_vel(2);
    double n_z = ray_vel(3);
    double v_mag = std::sqrt(n_x * n_x + n_y * n_y + n_z * n_z);

    if (v_mag < 1e-10) {
        return static_cast<float>(std::clamp(grav_factor, 0.1, 5.0));
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

    // Doppler suppression (P4): drop the orbital-velocity term, keeping only the
    // gravitational redshift, so the disk's approaching/receding asymmetry
    // collapses (the film's artistic choice).
    if (!config_.doppler_beaming) {
        return static_cast<float>(std::clamp(grav_factor, 0.1, 5.0));
    }

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
    double grav_factor = std::sqrt(std::max(0.0, 1.0 - 2.0 * M * r_d / (r2 + a2)));

    double n_x = ray_vel(1);
    double n_y = ray_vel(2);
    double n_z = ray_vel(3);
    double v_mag = std::sqrt(n_x * n_x + n_y * n_y + n_z * n_z);

    if (v_mag < 1e-10) {
        result.redshift = static_cast<float>(std::clamp(grav_factor, 0.1, 5.0));
        result.gfactor_grav = static_cast<float>(grav_factor);
        result.gfactor_gamma = 1.0f;
        result.gfactor_v_orb = static_cast<float>(v_orb);
        result.gfactor_cosine_coefficient = 0.0f;
        result.gfactor_sine_coefficient = 0.0f;
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
    result.gfactor_cosine_coefficient = static_cast<float>(A);
    result.gfactor_sine_coefficient = static_cast<float>(B);

    double v_dot_n = v_orb * A;
    double doppler_denom = gamma * (1.0 - v_dot_n);
    doppler_denom = std::clamp(doppler_denom, 0.1, 10.0);
    double g = grav_factor / doppler_denom;
    g = std::clamp(g, 0.1, 5.0);

    result.redshift = static_cast<float>(g);

    // Doppler suppression (P4): v_orb treated as zero for the brightness factor
    // (gamma -> 1, v.n -> 0, so g -> gravitational redshift), while A and B are
    // retained because they carry the disk-plane viewing geometry that limb
    // darkening needs, not the Doppler asymmetry. The default (true) path above
    // is untouched, so the pinned render does not move.
    if (!config_.doppler_beaming) {
        result.gfactor_gamma = 1.0f;
        result.gfactor_v_orb = 0.0f;
        result.redshift = static_cast<float>(std::clamp(grav_factor, 0.1, 5.0));
    }
}

// =============================================================================
// Volumetric disk.
// =============================================================================
float GeodesicTracer::ComputeScaleHeight(float r) {
    if (r <= 0) return 0;

    // H(r) = H_over_r r (r/r_ref)^H_power, referenced at the inner edge.
    float r_ref = config_.disk_inner;
    float r_ratio = r / r_ref;
    float H_over_r =
        config_.volumetric_scale_height_ratio * std::pow(r_ratio, config_.volumetric_flare_power);

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

float GeodesicTracer::ComputeVolumetricOpacityDensity(float r, float z, float phi) {
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
                       (std::sqrt(2.0f * static_cast<float>(std::numbers::pi)) * H) *
                       std::pow(r_ratio, -1.5f) * (H_ref / H);

    float density = kappa_rho0 * gaussian;
    if (config_.enable_turbulence) {
        const float spherical_r = std::sqrt(r * r + z * z);
        const float theta =
            spherical_r > 0.0f ? std::acos(std::clamp(z / spherical_r, -1.0f, 1.0f)) : 0.0f;
        density *= turbulence_noise::SampleDensityPerturbation(spherical_r, theta, phi,
                                                               config_.turbulence);
    }
    return density;
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

void GeodesicTracer::AccumulateVolumetricEmission(const Lightray& ray, const Vec4& entry_pos,
                                                  const Vec4& exit_pos, TraceResult& result) {
    // ray is retained in the signature for parity with the legacy interface; the
    // segment endpoints carry everything the ray march needs.
    double dx = exit_pos(1) - entry_pos(1);
    double dy = exit_pos(2) - entry_pos(2);
    double dz = exit_pos(3) - entry_pos(3);
    double path_length = std::sqrt(dx * dx + dy * dy + dz * dz);

    if (path_length < 1e-10) return;

    int N = config_.volumetric_samples;
    double ds = path_length / N;

    // Carry the radiative-transfer solution across integration segments. The
    // previous implementation restarted from zero for every RK45 segment and
    // overwrote the result, so the reported optical depth depended on the final
    // step rather than the traversed volume.
    double accumulated_tau = result.optical_depth;
    double accumulated_r = result.volumetric_emission[0];
    double accumulated_g = result.volumetric_emission[1];
    double accumulated_b = result.volumetric_emission[2];
    double active_path_length = 0.0;
    const double max_tau = static_cast<double>(config_.volumetric_tau_max);

    if (accumulated_tau >= max_tau) {
        result.optical_depth = static_cast<float>(max_tau);
        return;
    }

    double dir_x = dx / path_length;
    double dir_y = dy / path_length;
    double dir_z = dz / path_length;
    Vec4 ray_velocity;
    for (int component = 0; component < 4; ++component) {
        ray_velocity(component) = ray.velocity(component);
    }

    for (int i = 0; i < N; i++) {
        // Midpoint sample.
        double t = (i + 0.5) * ds;
        double x = entry_pos(1) + t * dir_x;
        double y = entry_pos(2) + t * dir_y;
        double z = entry_pos(3) + t * dir_z;

        const double cylindrical_r = std::sqrt(x * x + y * y);
        const double spherical_r = std::sqrt(cylindrical_r * cylindrical_r + z * z);
        const coordinates::Vec4Cart sample_cart{0.0, x, y, z};
        const double disk_r = coordinates::KerrSchildRadius(sample_cart, cached_a_ * cached_m_);
        const float phi = static_cast<float>(std::atan2(y, x));
        const float theta =
            spherical_r > 0.0
                ? static_cast<float>(std::acos(std::clamp(z / spherical_r, -1.0, 1.0)))
                : 0.0f;

        float disk_opacity = 0.0f;
        float disk_source = 0.0f;
        if (config_.enable_volumetric &&
            IsInVolumetricDisk(static_cast<float>(disk_r), static_cast<float>(z))) {
            disk_opacity = ComputeVolumetricOpacityDensity(static_cast<float>(disk_r),
                                                           static_cast<float>(z), phi);
            const float temperature =
                ComputeVolumetricTemperature(static_cast<float>(disk_r), static_cast<float>(z));
            if (disk_opacity > 0.0f && temperature > 0.0f) {
                const float emitted_source =
                    std::pow(temperature / config_.disk_temperature_inner, 4.0f);
                const float g = ComputeGFactor(static_cast<float>(disk_r), phi, ray_velocity);
                disk_source = core::color_modes::ObservedBolometricIntensity(emitted_source, g);
            }
        }

        float corona_opacity = 0.0f;
        float corona_source = 0.0f;
        if (config_.enable_corona &&
            corona_physics::IsInsideCorona(static_cast<float>(spherical_r), theta, phi,
                                           config_.corona, config_.disk_inner)) {
            corona_opacity = config_.corona.optical_depth / config_.corona.scale_height_M;
            const float seed = corona_physics::Emissivity(static_cast<float>(spherical_r), theta,
                                                          config_.corona, config_.disk_inner);
            corona_source = corona_physics::ComptonizedSource(seed, config_.corona);
        }

        const float total_opacity = disk_opacity + corona_opacity;
        if (total_opacity <= 0.0f) continue;
        const double remaining_tau = max_tau - accumulated_tau;
        const double dtau = std::min(static_cast<double>(total_opacity) * ds, remaining_tau);
        active_path_length += dtau / total_opacity;
        const double source_r =
            (disk_opacity * disk_source + corona_opacity * corona_source * 0.65f) / total_opacity;
        const double source_g =
            (disk_opacity * disk_source + corona_opacity * corona_source * 0.80f) / total_opacity;
        const double source_b =
            (disk_opacity * disk_source + corona_opacity * corona_source) / total_opacity;

        // Radiative transfer I_out = I_in exp(-dtau) + S (1 - exp(-dtau)).
        const double transmission = std::exp(-dtau);
        double one_minus_trans = 1.0 - transmission;

        accumulated_r = accumulated_r * transmission + source_r * one_minus_trans;
        accumulated_g = accumulated_g * transmission + source_g * one_minus_trans;
        accumulated_b = accumulated_b * transmission + source_b * one_minus_trans;

        accumulated_tau += dtau;

        if (accumulated_tau >= max_tau) break;
    }

    result.volumetric_emission[0] = static_cast<float>(accumulated_r);
    result.volumetric_emission[1] = static_cast<float>(accumulated_g);
    result.volumetric_emission[2] = static_cast<float>(accumulated_b);
    result.optical_depth = static_cast<float>(accumulated_tau);
    result.volumetric_path_length += static_cast<float>(active_path_length);
    result.volumetric_hit = result.volumetric_hit || (accumulated_tau > 0.01);
}

// =============================================================================
// Ray-bundle (geodesic deviation) machinery (P2).
// =============================================================================
namespace {

// Orthonormal transverse basis (e_a, e_b) spanning the plane perpendicular to a
// unit direction n; e_b = n x e_a completes a right-handed triad. The reference
// vector avoids the degenerate cross product when n is near the z axis.
void TransverseBasis(double nx, double ny, double nz, double& ax, double& ay, double& az,
                     double& bx, double& by, double& bz) {
    double rx = 0.0, ry = 0.0, rz = 1.0;
    if (std::abs(nz) > 0.9) {
        rx = 1.0;
        ry = 0.0;
        rz = 0.0;
    }
    double rn = rx * nx + ry * ny + rz * nz;
    ax = rx - rn * nx;
    ay = ry - rn * ny;
    az = rz - rn * nz;
    double al = std::sqrt(ax * ax + ay * ay + az * az);
    if (al < 1e-30) al = 1.0;
    ax /= al;
    ay /= al;
    az /= al;
    bx = ny * az - nz * ay;
    by = nz * ax - nx * az;
    bz = nx * ay - ny * ax;
}

}  // namespace

void GeodesicTracer::ComputeChristoffelCart(const Vec4& pos, double Gamma[4][4][4]) {
    Metric4d g;
    Tensor<Dual<double>, 4, 4, 4> dg;
    metric_->Evaluate(pos, g, dg);
    ChristoffelSymbols cs = TensorOps::Christoffel(g, dg);
    for (int mu = 0; mu < 4; ++mu)
        for (int nu = 0; nu < 4; ++nu)
            for (int rho = 0; rho < 4; ++rho) Gamma[mu][nu][rho] = cs.gamma(mu, nu, rho).real;
}

void GeodesicTracer::ComputeRiemannCart(const Vec4& pos, double R[4][4][4][4]) {
    for (int mu = 0; mu < 4; ++mu)
        for (int nu = 0; nu < 4; ++nu)
            for (int rho = 0; rho < 4; ++rho)
                for (int sig = 0; sig < 4; ++sig) R[mu][nu][rho][sig] = 0.0;

    double Gamma[4][4][4];
    ComputeChristoffelCart(pos, Gamma);

    double rr = std::sqrt(pos(1) * pos(1) + pos(2) * pos(2) + pos(3) * pos(3));
    // Finite-difference step scaled to the local length; the fixed relative size
    // matches kernels/gr_deviation.slang so the two paths difference identically.
    double eps = 1.0e-3 * std::max(1.0, rr);

    // dGamma[rho][mu][nu][sig] = d_rho Gamma^mu_nu_sig; time derivatives vanish
    // (the Kerr-Schild family is stationary).
    double dGamma[4][4][4][4];
    for (int mu = 0; mu < 4; ++mu)
        for (int nu = 0; nu < 4; ++nu)
            for (int sig = 0; sig < 4; ++sig) dGamma[0][mu][nu][sig] = 0.0;

    for (int d = 1; d < 4; ++d) {
        Vec4 xp = pos, xm = pos;
        xp(d) += eps;
        xm(d) -= eps;
        double Gp[4][4][4], Gm[4][4][4];
        ComputeChristoffelCart(xp, Gp);
        ComputeChristoffelCart(xm, Gm);
        for (int mu = 0; mu < 4; ++mu)
            for (int nu = 0; nu < 4; ++nu)
                for (int sig = 0; sig < 4; ++sig)
                    dGamma[d][mu][nu][sig] = (Gp[mu][nu][sig] - Gm[mu][nu][sig]) / (2.0 * eps);
    }

    // R^mu_nu_rho_sig = d_rho Gamma^mu_nu_sig - d_sig Gamma^mu_nu_rho
    //                 + Gamma^mu_lam_rho Gamma^lam_nu_sig
    //                 - Gamma^mu_lam_sig Gamma^lam_nu_rho.
    for (int mu = 0; mu < 4; ++mu)
        for (int nu = 0; nu < 4; ++nu)
            for (int rho = 0; rho < 4; ++rho)
                for (int sig = 0; sig < 4; ++sig) {
                    double val = dGamma[rho][mu][nu][sig] - dGamma[sig][mu][nu][rho];
                    for (int lam = 0; lam < 4; ++lam)
                        val += Gamma[mu][lam][rho] * Gamma[lam][nu][sig] -
                               Gamma[mu][lam][sig] * Gamma[lam][nu][rho];
                    R[mu][nu][rho][sig] = val;
                }
}

void GeodesicTracer::InitBundle(const Vec4& k, RayBundle& bundle) const {
    double nx = k(1), ny = k(2), nz = k(3);
    double nlen = std::sqrt(nx * nx + ny * ny + nz * nz);
    if (nlen < 1e-12) {
        nx = 1.0;
        ny = 0.0;
        nz = 0.0;
        nlen = 1.0;
    }
    nx /= nlen;
    ny /= nlen;
    nz /= nlen;

    double ax, ay, az, bx, by, bz;
    TransverseBasis(nx, ny, nz, ax, ay, az, bx, by, bz);

    double eps = static_cast<double>(config_.bundle_angular_size);

    bundle.xi[0] = Vec4();
    bundle.xi[1] = Vec4();
    bundle.V[0] = Vec4();
    bundle.V[1] = Vec4();

    if (config_.bundle_point_source) {
        // Pupil bundle: the rays leave a point with an angular spread, so xi = 0
        // and D xi / d lambda spans the transverse plane (P3 footprint).
        bundle.V[0](1) = eps * ax;
        bundle.V[0](2) = eps * ay;
        bundle.V[0](3) = eps * az;
        bundle.V[1](1) = eps * bx;
        bundle.V[1](2) = eps * by;
        bundle.V[1](3) = eps * bz;
    } else {
        // Parallel bundle (identity Jacobian, D xi / d lambda = 0), matching the
        // oracle's BeamStateD::Initialise so the two are comparable (P2).
        bundle.xi[0](1) = eps * ax;
        bundle.xi[0](2) = eps * ay;
        bundle.xi[0](3) = eps * az;
        bundle.xi[1](1) = eps * bx;
        bundle.xi[1](2) = eps * by;
        bundle.xi[1](3) = eps * bz;
    }
}

void GeodesicTracer::StepBundle(const Vec4& pos, const Vec4& k, double d_lambda,
                                RayBundle& bundle) {
    double Gamma[4][4][4];
    double R[4][4][4][4];
    ComputeChristoffelCart(pos, Gamma);
    ComputeRiemannCart(pos, R);

    // Deviation right-hand side, coefficients frozen over the step:
    //   d xi^mu / d lambda = V^mu - Gamma^mu_ab k^a xi^b
    //   d V^mu  / d lambda = -Gamma^mu_ab k^a V^b - R^mu_nu_rho_sig k^nu xi^rho k^sig
    auto rhs = [&](const Vec4& xi, const Vec4& V, Vec4& dxi, Vec4& dV) {
        for (int mu = 0; mu < 4; ++mu) {
            double gk_xi = 0.0, gk_V = 0.0;
            for (int a = 0; a < 4; ++a) {
                double ka = k(a);
                for (int b = 0; b < 4; ++b) {
                    gk_xi += Gamma[mu][a][b] * ka * xi(b);
                    gk_V += Gamma[mu][a][b] * ka * V(b);
                }
            }
            double r_term = 0.0;
            for (int nu = 0; nu < 4; ++nu)
                for (int rho = 0; rho < 4; ++rho)
                    for (int sig = 0; sig < 4; ++sig)
                        r_term += R[mu][nu][rho][sig] * k(nu) * xi(rho) * k(sig);
            dxi(mu) = V(mu) - gk_xi;
            dV(mu) = -gk_V - r_term;
        }
    };

    double h = d_lambda;
    for (int i = 0; i < 2; ++i) {
        Vec4 xi0 = bundle.xi[i];
        Vec4 V0 = bundle.V[i];

        Vec4 dxi1, dV1, dxi2, dV2, dxi3, dV3, dxi4, dV4;
        rhs(xi0, V0, dxi1, dV1);
        rhs(xi0 + dxi1 * (0.5 * h), V0 + dV1 * (0.5 * h), dxi2, dV2);
        rhs(xi0 + dxi2 * (0.5 * h), V0 + dV2 * (0.5 * h), dxi3, dV3);
        rhs(xi0 + dxi3 * h, V0 + dV3 * h, dxi4, dV4);

        bundle.xi[i] = xi0 + (dxi1 + dxi2 * 2.0 + dxi3 * 2.0 + dxi4) * (h / 6.0);
        bundle.V[i] = V0 + (dV1 + dV2 * 2.0 + dV3 * 2.0 + dV4) * (h / 6.0);
    }
}

Vec4 GeodesicTracer::TidalAcceleration(const Vec4& pos, const Vec4& k, const Vec4& xi) {
    double R[4][4][4][4];
    ComputeRiemannCart(pos, R);
    Vec4 accel;
    for (int mu = 0; mu < 4; ++mu) {
        double s = 0.0;
        for (int nu = 0; nu < 4; ++nu)
            for (int rho = 0; rho < 4; ++rho)
                for (int sig = 0; sig < 4; ++sig)
                    s += R[mu][nu][rho][sig] * k(nu) * xi(rho) * k(sig);
        accel(mu) = -s;  // D^2 xi^mu / d lambda^2 = -R^mu_nu_rho_sig k^nu xi^rho k^sig.
    }
    return accel;
}

double GeodesicTracer::KretschmannScalar(const Vec4& pos) {
    double R[4][4][4][4];
    ComputeRiemannCart(pos, R);

    Metric4d gm, gim;
    Tensor<Dual<double>, 4, 4, 4> dg;
    metric_->Evaluate(pos, gm, dg);
    Metric4d gim_metric = metric_->InverseMetric(pos, gim) ? gim : TensorOps::Inverse(gm);

    double g[4][4], gi[4][4];
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) {
            g[i][j] = gm(i, j).real;
            gi[i][j] = gim_metric(i, j).real;
        }

    // Fully covariant R_mu_nu_rho_sig = g_mu_alpha R^alpha_nu_rho_sig.
    double Rd[4][4][4][4];
    for (int mu = 0; mu < 4; ++mu)
        for (int nu = 0; nu < 4; ++nu)
            for (int rho = 0; rho < 4; ++rho)
                for (int sig = 0; sig < 4; ++sig) {
                    double s = 0.0;
                    for (int al = 0; al < 4; ++al) s += g[mu][al] * R[al][nu][rho][sig];
                    Rd[mu][nu][rho][sig] = s;
                }

    // Fully contravariant R^mu_nu_rho_sig = g^nu_b g^rho_c g^sig_d R^mu_b_c_d.
    double Ru[4][4][4][4];
    for (int mu = 0; mu < 4; ++mu)
        for (int nu = 0; nu < 4; ++nu)
            for (int rho = 0; rho < 4; ++rho)
                for (int sig = 0; sig < 4; ++sig) {
                    double s = 0.0;
                    for (int b = 0; b < 4; ++b)
                        for (int c = 0; c < 4; ++c)
                            for (int d = 0; d < 4; ++d)
                                s += gi[nu][b] * gi[rho][c] * gi[sig][d] * R[mu][b][c][d];
                    Ru[mu][nu][rho][sig] = s;
                }

    double K = 0.0;
    for (int mu = 0; mu < 4; ++mu)
        for (int nu = 0; nu < 4; ++nu)
            for (int rho = 0; rho < 4; ++rho)
                for (int sig = 0; sig < 4; ++sig) K += Rd[mu][nu][rho][sig] * Ru[mu][nu][rho][sig];
    return K;
}

void GeodesicTracer::FinaliseBundle(const RayBundle& bundle, const Vec4& k, double lambda,
                                    TraceResult::Beam& out) const {
    double nx = k(1), ny = k(2), nz = k(3);
    double nlen = std::sqrt(nx * nx + ny * ny + nz * nz);
    if (nlen < 1e-12) {
        nx = 1.0;
        ny = 0.0;
        nz = 0.0;
        nlen = 1.0;
    }
    nx /= nlen;
    ny /= nlen;
    nz /= nlen;

    double ex, ey, ez, fx, fy, fz;
    TransverseBasis(nx, ny, nz, ex, ey, ez, fx, fy, fz);

    // Project each deviation vector onto the transverse plane; the along-ray part
    // is a longitudinal (gauge) displacement that does not change the beam
    // cross-section. Columns of M are xi_0, xi_1 in the (e_a, e_b) basis.
    double a = bundle.xi[0](1) * ex + bundle.xi[0](2) * ey + bundle.xi[0](3) * ez;  // xi0.e_a
    double c = bundle.xi[0](1) * fx + bundle.xi[0](2) * fy + bundle.xi[0](3) * fz;  // xi0.e_b
    double b = bundle.xi[1](1) * ex + bundle.xi[1](2) * ey + bundle.xi[1](3) * ez;  // xi1.e_a
    double d = bundle.xi[1](1) * fx + bundle.xi[1](2) * fy + bundle.xi[1](3) * fz;  // xi1.e_b

    double det = a * d - b * c;
    double area = std::abs(det);
    double eps = static_cast<double>(config_.bundle_angular_size);
    double area0 = eps * eps;

    out.transverse_area = static_cast<float>(area);
    out.area_ratio = (area0 > 0.0) ? static_cast<float>(area / area0) : 0.0f;
    out.magnification = (area > 1e-30) ? static_cast<float>(area0 / area) : 1.0e12f;

    // Singular values of M give the ellipse semi-axes; the sum-of-squares
    // discriminant form is the same one the oracle uses (BeamStateD::UpdateGeometry).
    double p = a * a + b * b + c * c + d * d;
    double disc = p * p - 4.0 * det * det;
    if (disc < 0.0) disc = 0.0;
    double s = std::sqrt(disc);
    out.semi_major = static_cast<float>(std::sqrt(std::max(0.0, (p + s) / 2.0)));
    out.semi_minor = static_cast<float>(std::sqrt(std::max(0.0, (p - s) / 2.0)));
    // Output-plane orientation is the major eigenvector of M M^T. The
    // right-singular-vector expression (ab+cd) orients the input basis instead.
    double num = 2.0 * (a * c + b * d);
    double den = a * a + b * b - c * c - d * d;
    out.orientation = static_cast<float>(0.5 * std::atan2(num, den));

    // Angular footprint on the sky: transverse extent per unit affine length. For
    // the pupil bundle in flat space xi = (angular size) * lambda, so this returns
    // the pixel angular size; lensing stretches it (P3).
    double inv_lambda = (lambda > 1e-12) ? 1.0 / lambda : 0.0;
    out.footprint_major = static_cast<float>(out.semi_major * inv_lambda);
    out.footprint_minor = static_cast<float>(out.semi_minor * inv_lambda);
    out.valid = true;
}

}  // namespace sirius::backend
