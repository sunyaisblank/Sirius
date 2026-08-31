// Geodesic ray tracer implementation. The CPU path owns adaptive central-ray,
// Jacobi, transfer, and polarisation orchestration in one accepted trajectory.

#include "sirius/backend/cpu/geodesic_tracer.h"

#include "sirius/core/disk/novikov_thorne_disk.h"
#include "sirius/core/observer_frame.h"
#include "sirius/core/spectral/colour_modes.h"
#include "sirius/core/trace_boundary.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numbers>
#include <optional>

namespace sirius::backend {

using namespace sirius::core;

namespace {
struct ObserverSkySample {
    Metric4d metric;
    relativity::ObserverFrame observer;
    std::array<Vec4, 2> screen;
    Vec4 world_direction;
};

std::optional<ObserverSkySample> SampleEulerianSky(IMetric& metric_authority, const Vec4& position,
                                                   const Vec4& past_ray) {
    Metric4d metric;
    Tensor<Dual<double>, 4, 4, 4> derivatives;
    metric_authority.Evaluate(position, metric, derivatives);
    Metric4d inverse;
    if (!metric_authority.InverseMetric(position, inverse)) {
        inverse = TensorOps::Inverse(metric);
    }

    std::array<Vec4, 3> seeds;
    seeds[0](1) = 1.0;
    seeds[1](2) = 1.0;
    seeds[2](3) = 1.0;

    const auto observer = relativity::EulerianObserverFrame(metric, inverse, seeds);
    if (!observer.has_value()) return std::nullopt;
    const double frequency = TensorOps::InnerProduct(past_ray, observer->time, metric);
    if (!std::isfinite(frequency) || !(frequency > 0.0)) return std::nullopt;

    std::array<double, 3> local_direction{};
    for (std::size_t component = 0; component < local_direction.size(); ++component) {
        local_direction[component] =
            TensorOps::InnerProduct(past_ray, observer->spatial[component], metric) / frequency;
    }
    const auto screen = relativity::ObserverScreenBasis(*observer, local_direction);
    if (!screen.has_value()) return std::nullopt;

    ObserverSkySample sample{metric, *observer, *screen, Vec4{}};
    // The tetrad axes are seeded by fixed Cartesian chart axes, preserving the
    // catalogue's orientation while measuring the direction locally.
    sample.world_direction(1) = local_direction[0];
    sample.world_direction(2) = local_direction[1];
    sample.world_direction(3) = local_direction[2];
    return sample;
}

}  // namespace

bool GeodesicTracer::FindDiskIntersection(const Vec4& start_position, const Vec4& start_tangent,
                                          const Vec4& end_position, const Vec4& end_tangent,
                                          double d_lambda, float& intersection_r,
                                          float& intersection_phi, double& intersection_fraction,
                                          Vec4& intersection_position, Vec4& intersection_tangent) {
    const double start_z = start_position(3);
    const double end_z = end_position(3);
    if (!std::isfinite(start_z) || !std::isfinite(end_z) || !std::isfinite(d_lambda) ||
        !(d_lambda > 0.0)) {
        return false;
    }

    // z(s) is a cubic under the accepted-segment Hermite interpolant. Isolate
    // roots on the monotone intervals separated by z'(s)=0 rather than using
    // endpoint signs: a cubic can touch the plane or cross it twice while its
    // endpoints remain on the same side.
    const double dz_start = d_lambda * start_tangent(3);
    const double dz_end = d_lambda * end_tangent(3);
    const double coefficient_a = 2.0 * start_z - 2.0 * end_z + dz_start + dz_end;
    const double coefficient_b = -3.0 * start_z + 3.0 * end_z - 2.0 * dz_start - dz_end;
    const double coefficient_c = dz_start;
    const double coefficient_d = start_z;
    if (!std::isfinite(coefficient_a) || !std::isfinite(coefficient_b) ||
        !std::isfinite(coefficient_c) || !std::isfinite(coefficient_d)) {
        return false;
    }

    const double scale = std::max({1.0, std::abs(coefficient_a), std::abs(coefficient_b),
                                   std::abs(coefficient_c), std::abs(coefficient_d)});
    const double root_tolerance = 64.0 * std::numeric_limits<double>::epsilon() * scale;
    const double fraction_tolerance = 256.0 * std::numeric_limits<double>::epsilon();
    const auto evaluate_z = [&](double fraction) {
        return ((coefficient_a * fraction + coefficient_b) * fraction + coefficient_c) * fraction +
               coefficient_d;
    };
    const auto radius_at = [&](double fraction) {
        const AcceptedTraceSegmentSample sample = SampleAcceptedTraceSegment(
            start_position, start_tangent, end_position, end_tangent, d_lambda, fraction);
        const coordinates::Vec4Cart cart{sample.position(0), sample.position(1), sample.position(2),
                                         sample.position(3)};
        return std::pair{sample, coordinates::KerrSchildRadius(cart, cached_a_ * cached_m_)};
    };
    const auto accept_fraction = [&](double fraction) {
        const auto [sample, radius] = radius_at(fraction);
        if (!std::isfinite(radius) || radius < config_.disk_inner || radius > config_.disk_outer) {
            return false;
        }
        intersection_fraction = fraction;
        intersection_position = sample.position;
        intersection_tangent = sample.tangent;
        intersection_r = static_cast<float>(radius);
        intersection_phi = static_cast<float>(std::atan2(sample.position(2), sample.position(1)));
        return true;
    };

    if (std::abs(coefficient_a) <= root_tolerance && std::abs(coefficient_b) <= root_tolerance &&
        std::abs(coefficient_c) <= root_tolerance && std::abs(coefficient_d) <= root_tolerance) {
        // A ray contained in the equatorial plane intersects the disk where it
        // first enters the radial annulus. The governed CPU session caps an
        // accepted step at 0.1 M, far below the ISCO-to-20 M annulus width.
        if (accept_fraction(0.0)) return true;
        constexpr int kRadialBrackets = 64;
        double lower = 0.0;
        for (int interval = 1; interval <= kRadialBrackets; ++interval) {
            const double upper = static_cast<double>(interval) / kRadialBrackets;
            const double radius = radius_at(upper).second;
            const bool inside = std::isfinite(radius) && radius >= config_.disk_inner &&
                                radius <= config_.disk_outer;
            if (inside) {
                double outside_fraction = lower;
                double inside_fraction = upper;
                for (int iteration = 0; iteration < 64; ++iteration) {
                    const double midpoint = 0.5 * (outside_fraction + inside_fraction);
                    const double midpoint_radius = radius_at(midpoint).second;
                    const bool midpoint_inside = std::isfinite(midpoint_radius) &&
                                                 midpoint_radius >= config_.disk_inner &&
                                                 midpoint_radius <= config_.disk_outer;
                    if (midpoint_inside) {
                        inside_fraction = midpoint;
                    } else {
                        outside_fraction = midpoint;
                    }
                }
                return accept_fraction(inside_fraction);
            }
            lower = upper;
        }
        return false;
    }

    std::array<double, 4> breakpoints{0.0, 1.0, 1.0, 1.0};
    std::size_t breakpoint_count = 2;
    const double derivative_discriminant =
        coefficient_b * coefficient_b - 3.0 * coefficient_a * coefficient_c;
    const double discriminant_scale = std::max(
        {1.0, coefficient_b * coefficient_b, std::abs(3.0 * coefficient_a * coefficient_c)});
    const double discriminant_tolerance =
        128.0 * std::numeric_limits<double>::epsilon() * discriminant_scale;
    if (std::abs(coefficient_a) > root_tolerance &&
        derivative_discriminant >= -discriminant_tolerance) {
        const double root = std::sqrt(std::max(derivative_discriminant, 0.0));
        const double first = (-coefficient_b - root) / (3.0 * coefficient_a);
        const double second = (-coefficient_b + root) / (3.0 * coefficient_a);
        if (first > 0.0 && first < 1.0) breakpoints[breakpoint_count++] = first;
        if (second > 0.0 && second < 1.0) breakpoints[breakpoint_count++] = second;
    } else if (std::abs(coefficient_b) > root_tolerance) {
        const double critical = -coefficient_c / (2.0 * coefficient_b);
        if (critical > 0.0 && critical < 1.0) breakpoints[breakpoint_count++] = critical;
    }
    std::sort(breakpoints.begin(), breakpoints.end());

    double previous_candidate = -1.0;
    const auto consider = [&](double fraction) {
        if (previous_candidate >= 0.0 &&
            std::abs(fraction - previous_candidate) <= fraction_tolerance) {
            return false;
        }
        previous_candidate = fraction;
        return accept_fraction(fraction);
    };
    for (std::size_t interval = 0; interval + 1 < breakpoint_count; ++interval) {
        const double lower = breakpoints[interval];
        const double upper = breakpoints[interval + 1];
        const double lower_z = evaluate_z(lower);
        const double upper_z = evaluate_z(upper);
        if (std::abs(lower_z) <= root_tolerance && consider(lower)) return true;
        if ((lower_z < 0.0) != (upper_z < 0.0)) {
            double root_lower = lower;
            double root_upper = upper;
            double root_lower_z = lower_z;
            for (int iteration = 0; iteration < 64; ++iteration) {
                const double midpoint = 0.5 * (root_lower + root_upper);
                const double midpoint_z = evaluate_z(midpoint);
                if ((root_lower_z < 0.0) == (midpoint_z < 0.0)) {
                    root_lower = midpoint;
                    root_lower_z = midpoint_z;
                } else {
                    root_upper = midpoint;
                }
            }
            if (consider(0.5 * (root_lower + root_upper))) return true;
        }
    }
    const double final_z = evaluate_z(1.0);
    if (std::abs(final_z) <= root_tolerance && consider(1.0)) return true;
    return false;
}

// =============================================================================
// Initialise a Lightray from a camera ray.
// =============================================================================
Lightray GeodesicTracer::InitializeLightray(const CameraRay& camera_ray,
                                            relativity::ObserverFrame* launch_frame) {
    // Value-initialize the full device-facing record so screen coordinates and
    // alignment padding never carry indeterminate bytes into copies or hashes.
    Lightray ray{};

    // Camera ray origin is Boyer-Lindquist (t, r, theta, phi).
    double t = camera_ray.origin(0);
    double r = camera_ray.origin(1);
    double th = camera_ray.origin(2);
    double ph = camera_ray.origin(3);

    // Kerr-Schild Cartesian integration is regular on the axis. Preserve the
    // requested observer event instead of silently moving near-polar cameras.
    SIRIUS_PRE(std::isfinite(th) && th >= 0.0 && th <= std::numbers::pi);

    // Position: Boyer-Lindquist -> Kerr-Schild Cartesian using the spin-aware
    // oblate transform.
    coordinates::Vec4Bl pos_bl(t, r, th, ph);
    const double absolute_spin = cached_a_ * cached_m_;
    coordinates::Vec4Cart pos_cart = coordinates::BlToKerrSchildCart(pos_bl, absolute_spin);

    ray.position(0) = static_cast<float>(pos_cart.t);
    ray.position(1) = static_cast<float>(pos_cart.x);
    ray.position(2) = static_cast<float>(pos_cart.y);
    ray.position(3) = static_cast<float>(pos_cart.z);

    const double sin_th = std::sin(th);
    // CameraRay::direction is the screen ray in the camera rest frame, resolved
    // on the local (radial, +theta, +phi) axes.  Build those axes as coordinate
    // seeds at the actual oblate position; the metric-aware frame construction
    // below orthonormalises them before applying the observer boost.
    const double cartesian_phi = std::atan2(pos_cart.y, pos_cart.x);
    const double sin_ph = std::sin(cartesian_phi);
    const double cos_ph = std::cos(cartesian_phi);
    Metric4d g;
    Tensor<Dual<double>, 4, 4, 4> dg;

    Vec4 pos_double;
    for (int i = 0; i < 4; ++i) pos_double(i) = ray.position(i);

    metric_->Evaluate(pos_double, g, dg);
    Metric4d inverse_metric;
    if (!metric_->InverseMetric(pos_double, inverse_metric)) {
        inverse_metric = TensorOps::Inverse(g);
    }

    std::array<Vec4, 3> spatial_seeds;
    spatial_seeds[0](1) = sin_th * cos_ph;
    spatial_seeds[0](2) = sin_th * sin_ph;
    spatial_seeds[0](3) = std::cos(th);
    spatial_seeds[1](1) = std::cos(th) * cos_ph;
    spatial_seeds[1](2) = std::cos(th) * sin_ph;
    spatial_seeds[1](3) = -sin_th;
    spatial_seeds[2](1) = -sin_ph;
    spatial_seeds[2](2) = cos_ph;

    const auto reference_frame =
        relativity::EulerianObserverFrame(g, inverse_metric, spatial_seeds);
    SIRIUS_ASSERT(reference_frame.has_value());
    // Operator beta is screen-forward/up/right.  The CameraRay component basis
    // is radial/+theta/+phi, hence forward=-radial and up=-theta.
    const std::array<double, 3> local_beta{-camera_ray.beta_forward, -camera_ray.beta_up,
                                           camera_ray.beta_right};
    const auto camera_frame = relativity::BoostObserverFrame(*reference_frame, local_beta);
    SIRIUS_ASSERT(camera_frame.has_value());
    if (launch_frame != nullptr) *launch_frame = *camera_frame;
    const std::array<double, 3> rest_direction{camera_ray.direction(1), camera_ray.direction(2),
                                               camera_ray.direction(3)};
    const auto past_ray = relativity::PastDirectedCameraRay(*camera_frame, rest_direction);
    SIRIUS_ASSERT(past_ray.has_value());
    SIRIUS_ASSERT((*past_ray)(0) < 0.0);
    for (int component = 0; component < 4; ++component) {
        ray.velocity(component) = static_cast<float>((*past_ray)(component));
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
    // -ray.velocity is the physical future photon.  The frame construction
    // normalises its launch frequency -(-k).u_camera to exactly one.
    ray.ku_uobsu = 1.0f;

    return ray;
}

void GeodesicTracer::InitPolarisationFrame(const Lightray& ray,
                                           const std::array<Vec4, 2>& launch_screen,
                                           PolarisationFrame& frame) {
    Vec4 position = ray.position;
    Vec4 velocity = ray.velocity;

    frame.reference.position = position;
    frame.reference.velocity = velocity;
    frame.reference.polarisation = launch_screen[0];
    frame.reference.affine = 0.0;

    frame.perpendicular.position = position;
    frame.perpendicular.velocity = velocity;
    frame.perpendicular.polarisation = launch_screen[1];
    frame.perpendicular.affine = 0.0;

    ReconditionPolarisationFrame(frame, position, velocity);
}

void GeodesicTracer::AdvancePolarisationFrame(PolarisationFrame& frame, const Vec4& end_position,
                                              const Vec4& end_tangent, double d_lambda) {
    SIRIUS_PRE(d_lambda > 0.0);
    const Vec4 start_position = frame.reference.position;
    const Vec4 start_tangent = frame.reference.velocity;
    frame.reference.polarisation =
        ParallelTransportAlongAcceptedSegment(*metric_, start_position, start_tangent, end_position,
                                              end_tangent, d_lambda, frame.reference.polarisation);
    frame.perpendicular.polarisation = ParallelTransportAlongAcceptedSegment(
        *metric_, start_position, start_tangent, end_position, end_tangent, d_lambda,
        frame.perpendicular.polarisation);
    frame.reference.position = end_position;
    frame.reference.velocity = end_tangent;
    frame.reference.affine += d_lambda;
    frame.perpendicular.position = end_position;
    frame.perpendicular.velocity = end_tangent;
    frame.perpendicular.affine += d_lambda;
}

void GeodesicTracer::ReconditionPolarisationFrame(PolarisationFrame& frame, const Vec4& position,
                                                  const Vec4& velocity) {
    // Reset the carrier path to the accepted RK45 endpoint and project
    // numerical transport drift back into the Eulerian observer's physical
    // Sachs screen. A coordinate-time basis vector is not generally timelike
    // (notably inside the Kerr ergosphere), so it cannot define this screen.
    const auto sky = SampleEulerianSky(*metric_, position, velocity);
    SIRIUS_PRE(sky.has_value());

    frame.reference.position = position;
    frame.reference.velocity = velocity;
    const auto reference = relativity::ProjectToObserverScreen(
        sky->metric, sky->observer.time, velocity, frame.reference.polarisation);
    SIRIUS_PRE(reference.has_value());
    frame.reference.polarisation = *reference;

    frame.perpendicular.position = position;
    frame.perpendicular.velocity = velocity;
    const auto projected_perpendicular = relativity::ProjectToObserverScreen(
        sky->metric, sky->observer.time, velocity, frame.perpendicular.polarisation);
    SIRIUS_PRE(projected_perpendicular.has_value());
    Vec4 perpendicular = *projected_perpendicular;
    perpendicular -=
        frame.reference.polarisation *
        TensorOps::InnerProduct(perpendicular, frame.reference.polarisation, sky->metric);
    const double norm = TensorOps::InnerProduct(perpendicular, perpendicular, sky->metric);
    SIRIUS_PRE(norm > 0.0);
    frame.perpendicular.polarisation = perpendicular / std::sqrt(norm);
}

void GeodesicTracer::SetDiskPolarisation(const PolarisationFrame& frame,
                                         TraceResult::DiskCrossing& crossing) {
    const Vec4& position = frame.reference.position;
    const Vec4& past_wave = frame.reference.velocity;
    const Vec4 physical_wave = past_wave * -1.0;

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

    const double emitted_frequency = -TensorOps::InnerProduct(emitter, physical_wave, g);
    SIRIUS_PRE(emitted_frequency > 0.0);
    const Vec4 photon_direction = physical_wave / emitted_frequency - emitter;

    Vec4 disk_normal;
    disk_normal(3) = 1.0;
    const double normal_norm = TensorOps::InnerProduct(disk_normal, disk_normal, g);
    SIRIUS_PRE(normal_norm > 0.0);
    disk_normal = disk_normal / std::sqrt(normal_norm);

    const float inclination_cosine = static_cast<float>(
        std::clamp(std::abs(TensorOps::InnerProduct(photon_direction, disk_normal, g)), 0.0, 1.0));
    const auto atmosphere =
        polarised_emission::ChandrasekharElectronScatteringAtmosphere(inclination_cosine);
    SIRIUS_PRE(atmosphere.has_value());
    crossing.polarisation_degree = atmosphere->linear_polarisation_degree;
    crossing.polarisation_intensity_scale = atmosphere->intensity_scale;

    // The electric vector of the represented semi-infinite scattering atmosphere is
    // parallel to the disk plane, hence perpendicular to the disk normal
    // projected into the emitter's physical screen.  The emitter worldline,
    // rather than a coordinate-time basis vector, defines that local screen.
    const auto meridian_screen =
        relativity::ProjectToObserverScreen(g, emitter, past_wave, disk_normal);
    SIRIUS_PRE(meridian_screen.has_value());
    const Vec4& meridian = *meridian_screen;
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
    SIRIUS_PRE(IsRepresentedCameraRay(camera_ray) && camera_ray.active);
    TraceResult result;
    result.steps_taken = 0;
    result.numerical_failure = false;

    CacheMetricParameters();

    relativity::ObserverFrame launch_frame;
    Lightray ray = InitializeLightray(camera_ray, &launch_frame);

    if (HasInvalidState(ray)) {
        result.outcome = TraceResult::Outcome::MaxSteps;
        result.numerical_failure = true;
        return result;
    }

    float min_r = 1e10f;

    Vec4 prev_pos;
    for (int i = 0; i < 4; ++i) prev_pos(i) = ray.position(i);

    const std::array<double, 3> launch_direction{camera_ray.direction(1), camera_ray.direction(2),
                                                 camera_ray.direction(3)};
    const auto launch_screen = relativity::ObserverScreenBasis(launch_frame, launch_direction);
    SIRIUS_ASSERT(launch_screen.has_value());

    // Ray-bundle state (P2); propagated only when enabled, so the point-sampled
    // path is untouched. prev_vel and prev_pt define every accepted segment,
    // including disk and causal-boundary events when no transported feature is
    // enabled.
    RayBundle bundle;
    RayBundle previous_bundle;
    Vec4 prev_vel;
    double prev_pt = 0.0;
    if (config_.enable_ray_bundles) {
        InitBundle(*launch_screen, bundle);
    }
    PolarisationFrame polarisation_frame;
    PolarisationFrame previous_polarisation_frame;
    if (config_.enable_polarisation) {
        InitPolarisationFrame(ray, *launch_screen, polarisation_frame);
    }

    for (int step = 0; step < config_.max_steps; ++step) {
        for (int i = 0; i < 4; ++i) prev_pos(i) = ray.position(i);
        prev_vel = ray.velocity;
        prev_pt = static_cast<double>(ray.proper_time);
        if (config_.enable_polarisation) {
            previous_polarisation_frame = polarisation_frame;
        }
        if (config_.enable_ray_bundles) {
            previous_bundle = bundle;
        }

        const double pre_step_radius =
            std::sqrt(static_cast<double>(ray.position(1)) * ray.position(1) +
                      static_cast<double>(ray.position(2)) * ray.position(2) +
                      static_cast<double>(ray.position(3)) * ray.position(3));
        auto step_config = config_.integrator;
        if (config_.strong_field_radius > 0.0f && config_.strong_field_max_step > 0.0f &&
            pre_step_radius < config_.strong_field_radius) {
            step_config.max_step = std::min(step_config.max_step, config_.strong_field_max_step);
            ray.step_size = std::min(ray.step_size, step_config.max_step);
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

        double d_lambda = static_cast<double>(ray.proper_time) - prev_pt;
        const float min_radius_before_step = min_r;
        bool terminal_causal_boundary = false;

        // A finite causal boundary clips the accepted central-ray segment
        // before any coupled state or segment source advances. This makes the
        // central ray, Jacobi bundle, polarisation frame, volume transfer, and
        // disk search consume one event-synchronised affine interval.
        if (config_.finite_causal_boundary) {
            const double accepted_x = ray.position(1);
            const double accepted_y = ray.position(2);
            const double accepted_z = ray.position(3);
            const double accepted_radius = std::sqrt(
                accepted_x * accepted_x + accepted_y * accepted_y + accepted_z * accepted_z);
            if (accepted_radius > static_cast<double>(config_.escape_radius)) {
                const auto boundary =
                    FindCausalBoundaryEvent(prev_pos, prev_vel, ray.position, ray.velocity,
                                            d_lambda, static_cast<double>(config_.escape_radius));
                if (!boundary.has_value()) {
                    ray.position = prev_pos;
                    ray.velocity = prev_vel;
                    ray.proper_time = static_cast<float>(prev_pt);
                    ray.coordinate_time = static_cast<float>(prev_pos(0));
                    result.outcome = TraceResult::Outcome::MaxSteps;
                    result.numerical_failure = true;
                    break;
                }
                d_lambda *= boundary->fraction;
                ray.position = boundary->position;
                ray.velocity = boundary->tangent;
                ray.proper_time = static_cast<float>(prev_pt + d_lambda);
                ray.coordinate_time = static_cast<float>(boundary->position(0));
                terminal_causal_boundary = true;
            }
        }

        // Advance the ray bundle over the accepted affine step, sampling the
        // connection and curvature at the bundle integrator's own stages.
        if (config_.enable_ray_bundles) {
            if (d_lambda > 0.0) {
                StepBundle(prev_pos, prev_vel, ray.position, ray.velocity, d_lambda, bundle);
            }
        }
        if (config_.enable_polarisation && d_lambda > 0.0) {
            AdvancePolarisationFrame(polarisation_frame, ray.position, ray.velocity, d_lambda);
            ReconditionPolarisationFrame(polarisation_frame, ray.position, ray.velocity);
        }

        double x = ray.position(1);
        double y = ray.position(2);
        double z = ray.position(3);
        double r = std::sqrt(x * x + y * y + z * z);

        if (r < min_r) {
            min_r = static_cast<float>(r);
        }

        // 1. Volumetric disk (ray marching through a 3D disk volume).
        if (config_.enable_disk && config_.enable_volumetric) {
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

            if (was_inside || is_inside) {
                Vec4 segment_start, segment_end;
                for (int i = 0; i < 4; ++i) {
                    segment_start(i) = prev_pos(i);
                    segment_end(i) = ray.position(i);
                }

                AccumulateVolumetricEmission(prev_vel, ray.velocity, d_lambda, ray.ku_uobsu,
                                             segment_start, segment_end, result);
            }
        }

        // 2. Thin disk: terminate at the observer-nearest opaque-surface event.
        if (config_.enable_disk && !config_.enable_volumetric) {
            Vec4 curr_pos;
            curr_pos(0) = ray.position(0);
            curr_pos(1) = ray.position(1);
            curr_pos(2) = ray.position(2);
            curr_pos(3) = ray.position(3);
            float disk_r, disk_phi;

            double crossing_fraction = 0.0;
            Vec4 crossing_position;
            Vec4 ray_vel;
            if (FindDiskIntersection(prev_pos, prev_vel, curr_pos, ray.velocity, d_lambda, disk_r,
                                     disk_phi, crossing_fraction, crossing_position, ray_vel)) {
                // Each crossing is a different image order: crossing 0 primary,
                // 1 secondary, and so on; the first sets the primary outcome.
                if (result.num_disk_crossings < TraceResult::kMaxDiskCrossings) {
                    auto& crossing = result.disk_crossings[result.num_disk_crossings];
                    crossing.r = disk_r;
                    crossing.phi = disk_phi;
                    crossing.temperature = ComputeDiskTemperature(disk_r);
                    crossing.crossing_index = result.num_disk_crossings;
                    crossing.valid = true;

                    if (result.num_disk_crossings == 0) {
                        SetDiskFrequencyTransfer(disk_r, disk_phi, ray_vel, ray.ku_uobsu, result,
                                                 crossing);
                        crossing.redshift = result.redshift;

                        result.outcome = TraceResult::Outcome::DiskHit;
                        result.disk_radius = disk_r;
                        result.disk_phi = disk_phi;
                        result.disk_temperature = crossing.temperature;
                    } else {
                        const DiskFrequencySample transfer =
                            ComputeDiskFrequencySample(disk_r, disk_phi, ray_vel, ray.ku_uobsu);
                        crossing.full_redshift = static_cast<float>(transfer.transfer.full_g);
                        crossing.zamo_redshift = static_cast<float>(transfer.transfer.zamo_g);
                        crossing.redshift = config_.doppler_beaming ? crossing.full_redshift
                                                                    : crossing.zamo_redshift;
                        crossing.emission_cosine = transfer.emission_cosine;
                    }

                    if (config_.enable_polarisation) {
                        PolarisationFrame crossing_frame = previous_polarisation_frame;
                        if (crossing_fraction > 0.0 && d_lambda > 0.0) {
                            AdvancePolarisationFrame(crossing_frame, crossing_position, ray_vel,
                                                     crossing_fraction * d_lambda);
                        }
                        SetDiskPolarisation(crossing_frame, crossing);
                        polarisation_frame = crossing_frame;
                    }

                    if (config_.enable_ray_bundles) {
                        bundle = previous_bundle;
                        if (crossing_fraction > 0.0 && d_lambda > 0.0) {
                            StepBundle(prev_pos, prev_vel, crossing_position, ray_vel,
                                       crossing_fraction * d_lambda, bundle);
                        }
                    }

                    result.num_disk_crossings++;
                }

                // The disk is optically thick. The observer-nearest crossing is
                // the physical surface; higher-order images arise at other
                // screen rays, not by adding surfaces hidden behind this one.
                ray.position = crossing_position;
                ray.velocity = ray_vel;
                ray.proper_time = static_cast<float>(prev_pt + crossing_fraction * d_lambda);
                ray.coordinate_time = static_cast<float>(crossing_position(0));
                const double crossing_cartesian_radius =
                    std::sqrt(crossing_position(1) * crossing_position(1) +
                              crossing_position(2) * crossing_position(2) +
                              crossing_position(3) * crossing_position(3));
                min_r =
                    std::min(min_radius_before_step, static_cast<float>(crossing_cartesian_radius));
                result.min_radius = min_r;
                break;
            }
        }

        // 3. Horizon capture. A represented opaque-disk crossing on this
        // accepted observer-to-scene segment has already terminated above; it
        // must not be hidden merely because the step endpoint lies inside the
        // capture surface. This is the same event precedence used by Slang.
        if (metric_->InsideCaptureSurface(ray.position, config_.horizon_factor - 1.0)) {
            result.outcome = TraceResult::Outcome::Horizon;
            result.final_position(0) = ray.position(0);
            result.final_position(1) = ray.position(1);
            result.final_position(2) = ray.position(2);
            result.final_position(3) = ray.position(3);
            result.min_radius = min_r;
            break;
        }

        if (terminal_causal_boundary) {
            result.outcome = TraceResult::Outcome::Escaped;
            const auto sky = SampleEulerianSky(*metric_, ray.position, ray.velocity);
            SIRIUS_ASSERT(sky.has_value());
            if (!sky.has_value()) {
                result.outcome = TraceResult::Outcome::MaxSteps;
                result.numerical_failure = true;
                break;
            }
            result.final_direction = sky->world_direction;
            result.min_radius = min_r;
            break;
        }

        // 4. Escape to the configured outer boundary. Segment emission has
        // observer-first precedence, matching the device path. Finite causal
        // boundaries were already clipped and classified above.
        if (r > config_.escape_radius) {
            const double vx = ray.velocity(1);
            const double vy = ray.velocity(2);
            const double vz = ray.velocity(3);
            const double v_radial = (x * vx + y * vy + z * vz) / r;

            if (v_radial > 0) {
                result.outcome = TraceResult::Outcome::Escaped;
                const auto sky = SampleEulerianSky(*metric_, ray.position, ray.velocity);
                SIRIUS_ASSERT(sky.has_value());
                if (!sky.has_value()) {
                    result.outcome = TraceResult::Outcome::MaxSteps;
                    result.numerical_failure = true;
                    break;
                }
                result.final_direction = sky->world_direction;
                result.min_radius = min_r;
                break;
            }
        }
    }

    // Publish the actual terminal central-ray event for every outcome. For an
    // opaque disk this is the Hermite-rooted crossing, not the accepted RK45
    // endpoint beyond it; the bundle and polarisation frames use the same event.
    for (int component = 0; component < 4; ++component) {
        result.final_position(component) = ray.position(component);
    }
    result.affine_length = ray.proper_time;
    result.min_radius = std::min(result.min_radius, min_r);

    // Beam ellipse from the propagated bundle (P2), evaluated against the ray's
    // terminal direction. Only when ray bundles are enabled.
    if (config_.enable_ray_bundles) {
        FinaliseBundle(bundle, ray.position, ray.velocity, result.beam);
    }

    // No termination reached: result stays MaxSteps (default).
    return result;
}

// =============================================================================
// Keplerian orbital velocity.
// =============================================================================
float GeodesicTracer::ComputeOrbitalVelocity(float r) {
    const auto angular_velocity = relativity::TryKerrCircularOrbitAngularVelocity(
        cached_m_, cached_a_ * cached_m_, static_cast<double>(r));
    SIRIUS_ASSERT(angular_velocity.has_value());
    return static_cast<float>(*angular_velocity);
}

// =============================================================================
// Invariant disk frequency transfer.
// =============================================================================
GeodesicTracer::DiskFrequencySample GeodesicTracer::ComputeDiskFrequencySample(
    float r, float cartesian_phi, const Vec4& ray_vel, float observer_frequency) {
    const double mass = cached_m_;
    const double spin = cached_a_ * mass;
    const double radius = static_cast<double>(r);
    const double cylindrical_radius = std::sqrt(radius * radius + spin * spin);
    Vec4 position;
    position(1) = cylindrical_radius * std::cos(static_cast<double>(cartesian_phi));
    position(2) = cylindrical_radius * std::sin(static_cast<double>(cartesian_phi));

    Metric4d metric;
    Tensor<Dual<double>, 4, 4, 4> derivatives;
    metric_->Evaluate(position, metric, derivatives);
    // The live tracer advances the past-directed camera-to-source tangent.
    // Frequency transfer is defined for the physical future photon, -k_past.
    const Vec4 physical_photon = ray_vel * -1.0;
    const Vec4 covector = TensorOps::LowerIndex(physical_photon, metric);
    const double energy = -covector(0);
    const double angular_momentum = -position(2) * covector(1) + position(1) * covector(2);
    const auto transfer = relativity::KerrDiskTransfer(
        static_cast<double>(observer_frequency), energy, angular_momentum, mass, spin, radius);
    SIRIUS_ASSERT(transfer.has_value());

    const double emission_cosine =
        std::clamp(std::abs(covector(3)) / transfer->emitter_frequency, 0.0, 1.0);
    return DiskFrequencySample{*transfer, static_cast<float>(emission_cosine)};
}

// =============================================================================
// Store the circular-emitter physical branch and the explicit ZAMO diagnostic.
// =============================================================================
void GeodesicTracer::SetDiskFrequencyTransfer(float r, float phi, const Vec4& ray_vel,
                                              float observer_frequency, TraceResult& result,
                                              TraceResult::DiskCrossing& crossing) {
    const DiskFrequencySample sample =
        ComputeDiskFrequencySample(r, phi, ray_vel, observer_frequency);
    result.full_disk_redshift = static_cast<float>(sample.transfer.full_g);
    result.zamo_disk_redshift = static_cast<float>(sample.transfer.zamo_g);
    crossing.full_redshift = result.full_disk_redshift;
    crossing.zamo_redshift = result.zamo_disk_redshift;
    result.redshift =
        config_.doppler_beaming ? result.full_disk_redshift : result.zamo_disk_redshift;
    crossing.emission_cosine = sample.emission_cosine;
}

// =============================================================================
// Volumetric disk.
// =============================================================================
float GeodesicTracer::ComputeScaleHeight(float r) {
    if (r <= 0) return 0;

    // Phenomenological Gaussian scale height
    // H(r)/r = clamp[(H/r)_ref (r/r_ref)^h_power, 0.01, 0.5].
    // The saturation bounds are part of the represented model; the operator's
    // reference value is validated inside the same interval.
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

    // Normalise over the represented finite support, not an unrepresented
    // infinite Gaussian tail. Thus the integral from -3H to +3H is exactly
    // tau_ref (r/r_ref)^(-3/2).
    float density = static_cast<float>(core::volumetric_disk::TruncatedGaussianOpacityDensity(
        config_.volumetric_tau_midplane, r, config_.disk_inner, H, z));
    if (config_.enable_turbulence) {
        const float spherical_r = std::sqrt(r * r + z * z);
        const float theta =
            spherical_r > 0.0f ? std::acos(std::clamp(z / spherical_r, -1.0f, 1.0f)) : 0.0f;
        density *= turbulence_noise::SampleDensityPerturbation(spherical_r, theta, phi,
                                                               config_.turbulence);
    }
    return density;
}

float GeodesicTracer::ComputeVolumetricTemperature(float r, [[maybe_unused]] float z) {
    if (!IsInVolumetricDisk(r, z)) return 0;
    // Vertically isothermal closure: the represented volume borrows the
    // selected thin-disk radial temperature and makes no ungrounded atmosphere
    // gradient claim.
    return ComputeDiskTemperature(r);
}

void GeodesicTracer::AccumulateVolumetricEmission(const Vec4& entry_velocity,
                                                  const Vec4& exit_velocity, double affine_length,
                                                  float observer_frequency, const Vec4& entry_pos,
                                                  const Vec4& exit_pos, TraceResult& result) {
    if (!std::isfinite(affine_length) || affine_length <= 0.0) return;

    int N = config_.volumetric_samples;
    const double d_lambda = affine_length / N;

    // Carry the radiative-transfer solution across integration segments. The
    // previous implementation restarted from zero for every RK45 segment and
    // overwrote the result, so the reported optical depth depended on the final
    // step rather than the traversed volume.
    relativity::GreyTransferState transfer{
        {result.volumetric_emission[0], result.volumetric_emission[1],
         result.volumetric_emission[2]},
        result.optical_depth};
    double active_affine_length = 0.0;
    const double max_tau = static_cast<double>(config_.volumetric_tau_max);

    if (transfer.optical_depth >= max_tau) {
        result.optical_depth = static_cast<float>(max_tau);
        return;
    }

    for (int i = 0; i < N; i++) {
        const double fraction = (static_cast<double>(i) + 0.5) / N;
        const AcceptedTraceSegmentSample sample = SampleAcceptedTraceSegment(
            entry_pos, entry_velocity, exit_pos, exit_velocity, affine_length, fraction);
        const Vec4& position = sample.position;
        const Vec4& past_velocity = sample.tangent;
        const double x = position(1);
        const double y = position(2);
        const double z = position(3);

        const coordinates::Vec4Cart sample_cart{position(0), x, y, z};
        const double disk_r = coordinates::KerrSchildRadius(sample_cart, cached_a_ * cached_m_);
        const float phi = static_cast<float>(std::atan2(y, x));

        Metric4d metric;
        Tensor<Dual<double>, 4, 4, 4> derivatives;
        metric_->Evaluate(position, metric, derivatives);
        Metric4d inverse_metric;
        if (!metric_->InverseMetric(position, inverse_metric)) {
            inverse_metric = TensorOps::Inverse(metric);
        }
        const double inverse_g_tt = inverse_metric(0, 0).real;
        if (!std::isfinite(inverse_g_tt) || !(inverse_g_tt < 0.0)) continue;
        const double lapse = 1.0 / std::sqrt(-inverse_g_tt);
        Vec4 eulerian;
        for (int component = 0; component < 4; ++component) {
            eulerian(component) = -lapse * inverse_metric(component, 0).real;
        }
        const double eulerian_frequency = TensorOps::InnerProduct(past_velocity, eulerian, metric);

        double disk_dtau = 0.0;
        core::spectral::Rgb disk_source;
        if (config_.enable_volumetric &&
            IsInVolumetricDisk(static_cast<float>(disk_r), static_cast<float>(z))) {
            const float disk_opacity = ComputeVolumetricOpacityDensity(static_cast<float>(disk_r),
                                                                       static_cast<float>(z), phi);
            const float temperature =
                ComputeVolumetricTemperature(static_cast<float>(disk_r), static_cast<float>(z));
            if (disk_opacity > 0.0f && temperature > 0.0f) {
                const double omega = ComputeOrbitalVelocity(static_cast<float>(disk_r));
                Vec4 emitter;
                emitter(0) = 1.0;
                emitter(1) = -omega * y;
                emitter(2) = omega * x;
                const double emitter_norm = TensorOps::InnerProduct(emitter, emitter, metric);
                if (!std::isfinite(emitter_norm) || !(emitter_norm < 0.0)) continue;
                emitter = emitter / std::sqrt(-emitter_norm);
                const double emitted_frequency =
                    TensorOps::InnerProduct(past_velocity, emitter, metric);
                if (!std::isfinite(emitted_frequency) || !(emitted_frequency > 0.0) ||
                    !std::isfinite(eulerian_frequency) || !(eulerian_frequency > 0.0)) {
                    continue;
                }
                const auto disk_path =
                    relativity::ComovingPathLength(past_velocity, emitter, metric, d_lambda);
                if (!disk_path.has_value()) continue;
                disk_dtau = static_cast<double>(disk_opacity) * *disk_path;
                const float emitted_source =
                    std::pow(temperature / config_.disk_temperature_inner, 4.0f);
                const double source_frequency =
                    config_.doppler_beaming ? emitted_frequency : eulerian_frequency;
                const float g = static_cast<float>(observer_frequency / source_frequency);
                disk_source = core::color_modes::ApplyColorMode(
                    config_.color_mode, temperature, g, emitted_source, nullptr,
                    config_.disk_temperature_scale_kelvin);
            }
        }

        double dtau = disk_dtau;
        if (!std::isfinite(dtau) || !(dtau > 0.0)) continue;
        const std::array<double, 3> source{disk_source.r, disk_source.g, disk_source.b};
        const auto accepted =
            relativity::AccumulateObserverToSourceLayer(transfer, source, dtau, max_tau);
        SIRIUS_ASSERT(accepted.has_value());
        if (!accepted.has_value()) continue;
        active_affine_length += d_lambda * *accepted;

        if (transfer.optical_depth >= max_tau) break;
    }

    result.volumetric_emission[0] = static_cast<float>(transfer.observed_emission[0]);
    result.volumetric_emission[1] = static_cast<float>(transfer.observed_emission[1]);
    result.volumetric_emission[2] = static_cast<float>(transfer.observed_emission[2]);
    result.optical_depth = static_cast<float>(transfer.optical_depth);
    result.volumetric_affine_length += static_cast<float>(active_affine_length);
    result.volumetric_hit = result.volumetric_hit || (transfer.optical_depth > 0.01);
}

// =============================================================================
// Ray-bundle (geodesic deviation) machinery (P2).
// =============================================================================
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
    // Scale the central-difference stencil to the local length. The CPU path is
    // double precision, so its smaller relative step resolves the analytic
    // Schwarzschild tidal contraction below 5e-7; the fp32 device path uses the
    // same stencil with its own precision-appropriate step.
    constexpr double kRelativeRiemannStep = 2.5e-5;
    double eps = kRelativeRiemannStep * std::max(1.0, rr);

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

void GeodesicTracer::InitBundle(const std::array<Vec4, 2>& launch_screen, RayBundle& bundle) const {
    double eps = static_cast<double>(config_.bundle_angular_size);

    bundle.xi[0] = Vec4();
    bundle.xi[1] = Vec4();
    bundle.V[0] = Vec4();
    bundle.V[1] = Vec4();

    if (config_.bundle_point_source) {
        // Pupil bundle: the rays leave a point with an angular spread, so xi = 0
        // and D xi / d lambda spans the observer's Sachs screen (P3 footprint).
        bundle.V[0] = launch_screen[0] * eps;
        bundle.V[1] = launch_screen[1] * eps;
    } else {
        // Parallel bundle (identity Jacobian, D xi / d lambda = 0), matching the
        // oracle's BeamStateD::Initialise so the two are comparable (P2).
        bundle.xi[0] = launch_screen[0] * eps;
        bundle.xi[1] = launch_screen[1] * eps;
    }
}

void GeodesicTracer::StepBundle(const Vec4& start_position, const Vec4& start_tangent,
                                const Vec4& end_position, const Vec4& end_tangent, double d_lambda,
                                RayBundle& bundle) {
    // Deviation right-hand side at one central-ray stage:
    //   d xi^mu / d lambda = V^mu - Gamma^mu_ab k^a xi^b
    //   d V^mu  / d lambda = -Gamma^mu_ab k^a V^b - R^mu_nu_rho_sig k^nu xi^rho k^sig
    auto rhs = [&](const Vec4& position, const Vec4& tangent, const RayBundle& state,
                   RayBundle& derivative) {
        double Gamma[4][4][4];
        double R[4][4][4][4];
        ComputeChristoffelCart(position, Gamma);
        ComputeRiemannCart(position, R);
        for (int column = 0; column < 2; ++column) {
            const Vec4& xi = state.xi[column];
            const Vec4& V = state.V[column];
            Vec4& dxi = derivative.xi[column];
            Vec4& dV = derivative.V[column];
            for (int mu = 0; mu < 4; ++mu) {
                double gk_xi = 0.0, gk_V = 0.0;
                for (int a = 0; a < 4; ++a) {
                    const double ka = tangent(a);
                    for (int b = 0; b < 4; ++b) {
                        gk_xi += Gamma[mu][a][b] * ka * xi(b);
                        gk_V += Gamma[mu][a][b] * ka * V(b);
                    }
                }
                double r_term = 0.0;
                for (int nu = 0; nu < 4; ++nu)
                    for (int rho = 0; rho < 4; ++rho)
                        for (int sig = 0; sig < 4; ++sig)
                            r_term += R[mu][nu][rho][sig] * tangent(nu) * xi(rho) * tangent(sig);
                dxi(mu) = V(mu) - gk_xi;
                dV(mu) = -gk_V - r_term;
            }
        }
    };

    const auto advance = [](const RayBundle& state, const RayBundle& derivative, double amount) {
        RayBundle advanced;
        for (int column = 0; column < 2; ++column) {
            advanced.xi[column] = state.xi[column] + derivative.xi[column] * amount;
            advanced.V[column] = state.V[column] + derivative.V[column] * amount;
        }
        return advanced;
    };

    const double h = d_lambda;
    // Cubic Hermite central-ray midpoint. It matches both accepted endpoint
    // events and tangents and is fourth-order accurate for a smooth geodesic.
    const AcceptedTraceSegmentSample midpoint = SampleAcceptedTraceSegment(
        start_position, start_tangent, end_position, end_tangent, h, 0.5);

    RayBundle stage1, stage2, stage3, stage4;
    rhs(start_position, start_tangent, bundle, stage1);
    rhs(midpoint.position, midpoint.tangent, advance(bundle, stage1, 0.5 * h), stage2);
    rhs(midpoint.position, midpoint.tangent, advance(bundle, stage2, 0.5 * h), stage3);
    rhs(end_position, end_tangent, advance(bundle, stage3, h), stage4);

    for (int column = 0; column < 2; ++column) {
        bundle.xi[column] += (stage1.xi[column] + stage2.xi[column] * 2.0 +
                              stage3.xi[column] * 2.0 + stage4.xi[column]) *
                             (h / 6.0);
        bundle.V[column] += (stage1.V[column] + stage2.V[column] * 2.0 + stage3.V[column] * 2.0 +
                             stage4.V[column]) *
                            (h / 6.0);
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

void GeodesicTracer::FinaliseBundle(const RayBundle& bundle, const Vec4& position, const Vec4& k,
                                    TraceResult::Beam& out) const {
    const auto sky = SampleEulerianSky(*metric_, position, k);
    if (!sky.has_value()) return;

    // Project with the spacetime metric onto the terminal observer's Sachs
    // screen. Longitudinal gauge additions proportional to k vanish here.
    const double a = TensorOps::InnerProduct(bundle.xi[0], sky->screen[0], sky->metric);
    const double c = TensorOps::InnerProduct(bundle.xi[0], sky->screen[1], sky->metric);
    const double b = TensorOps::InnerProduct(bundle.xi[1], sky->screen[0], sky->metric);
    const double d = TensorOps::InnerProduct(bundle.xi[1], sky->screen[1], sky->metric);

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

    // The endpoint lies on the asymptotic celestial sphere. Its physical
    // transverse displacement divided by that sphere's radius is the angular
    // footprint; affine distance is not a coordinate-independent substitute.
    const double radius = std::sqrt(position(1) * position(1) + position(2) * position(2) +
                                    position(3) * position(3));
    const double inverse_radius = radius > 1.0e-12 ? 1.0 / radius : 0.0;
    out.footprint_major = static_cast<float>(out.semi_major * inverse_radius);
    out.footprint_minor = static_cast<float>(out.semi_minor * inverse_radius);
    out.valid = true;
}

}  // namespace sirius::backend
