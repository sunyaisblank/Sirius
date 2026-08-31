#pragma once

// Invariant frequency-transfer primitives for stationary, axisymmetric
// spacetimes.  The render path integrates in Kerr-Schild Cartesian
// coordinates, while circular disk motion is most compactly specified by the
// Boyer-Lindquist Killing fields.  Energy E=-p_t and axial angular momentum
// L_z=p_phi are invariant under the stationary/axisymmetric chart change, so
// they form the exact bridge without introducing a Euclidean orbital speed.
//
// For any stationary axisymmetric frame
// u = u^t(\partial_t + Omega \partial_phi), the measured photon frequency is
//
//   nu = -p.u = u^t(E - Omega L_z),
//
// and g=nu_observer/nu_emitter.  The locally non-rotating (ZAMO) frame supplies
// the physically timelike gravitational/frame-dragging-only branch used as an
// explicitly nonphysical diagnostic substitution when orbital Doppler transfer
// is disabled; unlike a static observer it remains valid inside the Kerr
// ergosphere.

#include "sirius/base/contracts.h"
#include "sirius/core/kerr_orbits.h"
#include "sirius/core/tensor.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <optional>

namespace sirius::core::relativity {

struct StationaryAxisymmetricFrame {
    double angular_velocity = 0.0;
    double time_component = 0.0;
};

struct KerrFrameFrequencyTransfer {
    double g = 0.0;
    double frame_frequency = 0.0;
    StationaryAxisymmetricFrame frame;
};

struct KerrDiskFrequencyTransfer {
    double full_g = 0.0;
    double zamo_g = 0.0;
    double emitter_frequency = 0.0;
    double zamo_frequency = 0.0;
    StationaryAxisymmetricFrame emitter;
    StationaryAxisymmetricFrame zamo;
};

struct GreyTransferState {
    std::array<double, 3> observed_emission{};
    double optical_depth = 0.0;
};

// Stable evaluation of the fraction absorbed/emitted by one homogeneous grey
// layer.  The direct 1-exp(-delta_tau) form loses the entire first-order term
// for sufficiently optically thin layers, even though those layers are inside
// the represented non-negative optical-depth domain.
[[nodiscard]] inline std::optional<double> GreyLayerAbsorbedFraction(double delta_tau) {
    if (!std::isfinite(delta_tau) || delta_tau < 0.0) return std::nullopt;
    return -std::expm1(-delta_tau);
}

// For a past-directed ray k and a future-directed fluid worldline u, k.u is
// the positive comoving photon frequency. With observer-normalised affine
// parameter, the proper path traversed in that fluid frame is (k.u) dlambda.
[[nodiscard]] inline std::optional<double> ComovingPathLength(const Vec4& past_ray,
                                                              const Vec4& fluid_velocity,
                                                              const Metric4d& metric,
                                                              double affine_length) {
    if (!std::isfinite(affine_length) || !(affine_length > 0.0)) return std::nullopt;
    const double norm = TensorOps::InnerProduct(fluid_velocity, fluid_velocity, metric);
    const double frequency = TensorOps::InnerProduct(past_ray, fluid_velocity, metric);
    if (!std::isfinite(norm) || !(norm < 0.0) || !std::isfinite(frequency) || !(frequency > 0.0)) {
        return std::nullopt;
    }
    return frequency * affine_length;
}

// Formal grey transfer while marching a past ray from the observer toward the
// source. A newly visited (farther) layer is attenuated by foreground tau; it
// does not attenuate emission already accumulated closer to the observer.
// Returns the accepted fraction of delta_tau when max_tau clips the layer.
[[nodiscard]] inline std::optional<double> AccumulateObserverToSourceLayer(
    GreyTransferState& state, const std::array<double, 3>& observed_source, double delta_tau,
    double max_tau) {
    if (!std::isfinite(state.optical_depth) || state.optical_depth < 0.0 ||
        !std::isfinite(delta_tau) || !(delta_tau > 0.0) || !std::isfinite(max_tau) ||
        !(max_tau > state.optical_depth)) {
        return std::nullopt;
    }
    for (const double channel : observed_source) {
        if (!std::isfinite(channel) || channel < 0.0) return std::nullopt;
    }
    const double accepted_tau = std::min(delta_tau, max_tau - state.optical_depth);
    const double accepted_fraction = accepted_tau / delta_tau;
    const auto absorbed_fraction = GreyLayerAbsorbedFraction(accepted_tau);
    if (!absorbed_fraction.has_value()) return std::nullopt;
    const double layer_weight = std::exp(-state.optical_depth) * *absorbed_fraction;
    for (std::size_t channel = 0; channel < state.observed_emission.size(); ++channel) {
        state.observed_emission[channel] += observed_source[channel] * layer_weight;
    }
    state.optical_depth += accepted_tau;
    return accepted_fraction;
}

[[nodiscard]] inline std::optional<StationaryAxisymmetricFrame> StationaryAxisymmetricFrameAt(
    double g_tt, double g_t_phi, double g_phi_phi, double angular_velocity) {
    if (!std::isfinite(g_tt) || !std::isfinite(g_t_phi) || !std::isfinite(g_phi_phi) ||
        !std::isfinite(angular_velocity) || !(g_phi_phi > 0.0)) {
        return std::nullopt;
    }
    const double norm =
        g_tt + 2.0 * angular_velocity * g_t_phi + angular_velocity * angular_velocity * g_phi_phi;
    if (!std::isfinite(norm) || !(norm < 0.0)) return std::nullopt;
    return StationaryAxisymmetricFrame{angular_velocity, 1.0 / std::sqrt(-norm)};
}

// Frequency transfer for an arbitrary-latitude stationary azimuthal frame in
// Kerr.  cos_theta is the Boyer-Lindquist polar cosine; in Kerr-Schild
// Cartesian coordinates it is exactly z/r.  The Killing quantities bridge the
// charts without importing the radial motion of a coordinate-slicing normal.
[[nodiscard]] inline std::optional<KerrFrameFrequencyTransfer>
TryKerrStationaryFrameFrequencyTransfer(double observer_frequency, double photon_energy,
                                        double photon_angular_momentum, double mass, double spin,
                                        double radius, double cos_theta, double angular_velocity) {
    if (!std::isfinite(observer_frequency) || !(observer_frequency > 0.0) ||
        !std::isfinite(photon_energy) || !(photon_energy > 0.0) ||
        !std::isfinite(photon_angular_momentum) || !std::isfinite(mass) || !(mass > 0.0) ||
        !std::isfinite(spin) || std::abs(spin) > mass || !std::isfinite(radius) ||
        !(radius > 0.0) || !std::isfinite(cos_theta) || std::abs(cos_theta) >= 1.0 ||
        !std::isfinite(angular_velocity)) {
        return std::nullopt;
    }

    const double sin_theta_squared = (1.0 - cos_theta) * (1.0 + cos_theta);
    const double sigma = radius * radius + spin * spin * cos_theta * cos_theta;
    if (!std::isfinite(sigma) || !(sigma > 0.0) || !(sin_theta_squared > 0.0)) {
        return std::nullopt;
    }
    const double g_tt = -(1.0 - 2.0 * mass * radius / sigma);
    const double g_t_phi = -2.0 * mass * spin * radius * sin_theta_squared / sigma;
    const double g_phi_phi =
        sin_theta_squared * (radius * radius + spin * spin +
                             2.0 * mass * spin * spin * radius * sin_theta_squared / sigma);
    const auto frame = StationaryAxisymmetricFrameAt(g_tt, g_t_phi, g_phi_phi, angular_velocity);
    if (!frame) return std::nullopt;

    const double frame_frequency =
        frame->time_component * (photon_energy - frame->angular_velocity * photon_angular_momentum);
    if (!std::isfinite(frame_frequency) || !(frame_frequency > 0.0)) {
        return std::nullopt;
    }
    return KerrFrameFrequencyTransfer{observer_frequency / frame_frequency, frame_frequency,
                                      *frame};
}

// Locally non-rotating Kerr frame at arbitrary latitude.  Unlike the normal
// to a horizon-penetrating Kerr-Schild time slice, this ZAMO has no radial
// Boyer-Lindquist component and therefore removes only circular-emitter motion.
[[nodiscard]] inline std::optional<KerrFrameFrequencyTransfer> KerrZamoFrequencyTransfer(
    double observer_frequency, double photon_energy, double photon_angular_momentum, double mass,
    double spin, double radius, double cos_theta) {
    if (!std::isfinite(mass) || !(mass > 0.0) || !std::isfinite(spin) || std::abs(spin) > mass ||
        !std::isfinite(radius) || !(radius > 0.0) || !std::isfinite(cos_theta) ||
        std::abs(cos_theta) >= 1.0) {
        return std::nullopt;
    }
    const double sin_theta_squared = (1.0 - cos_theta) * (1.0 + cos_theta);
    const double sigma = radius * radius + spin * spin * cos_theta * cos_theta;
    const double g_t_phi = -2.0 * mass * spin * radius * sin_theta_squared / sigma;
    const double g_phi_phi =
        sin_theta_squared * (radius * radius + spin * spin +
                             2.0 * mass * spin * spin * radius * sin_theta_squared / sigma);
    if (!std::isfinite(g_phi_phi) || !(g_phi_phi > 0.0)) return std::nullopt;
    return TryKerrStationaryFrameFrequencyTransfer(observer_frequency, photon_energy,
                                                   photon_angular_momentum, mass, spin, radius,
                                                   cos_theta, -g_t_phi / g_phi_phi);
}

// Exact equatorial Kerr transfer for a circular geodesic emitter and a ZAMO.
// M and r share geometric units; a is the signed dimensional Kerr parameter.
// observer_frequency is the positive -k.u measured at the camera.  E and L_z
// are the photon's conserved Killing quantities in the same affine scaling.
[[nodiscard]] inline std::optional<KerrDiskFrequencyTransfer> KerrDiskTransfer(
    double observer_frequency, double photon_energy, double photon_angular_momentum, double mass,
    double spin, double radius) {
    if (!std::isfinite(observer_frequency) || !(observer_frequency > 0.0) ||
        !std::isfinite(photon_energy) || !(photon_energy > 0.0) ||
        !std::isfinite(photon_angular_momentum) || !std::isfinite(mass) || !(mass > 0.0) ||
        !std::isfinite(spin) || std::abs(spin) > mass || !std::isfinite(radius) ||
        !(radius > 0.0)) {
        return std::nullopt;
    }

    const auto emitter_omega = TryKerrCircularOrbitAngularVelocity(mass, spin, radius);
    if (!emitter_omega) return std::nullopt;
    const auto emitter = TryKerrStationaryFrameFrequencyTransfer(observer_frequency, photon_energy,
                                                                 photon_angular_momentum, mass,
                                                                 spin, radius, 0.0, *emitter_omega);
    const auto zamo = KerrZamoFrequencyTransfer(observer_frequency, photon_energy,
                                                photon_angular_momentum, mass, spin, radius, 0.0);
    if (!emitter || !zamo) return std::nullopt;

    return KerrDiskFrequencyTransfer{
        emitter->g,     zamo->g,    emitter->frame_frequency, zamo->frame_frequency,
        emitter->frame, zamo->frame};
}

// Static-observer special case. Both g_tt values must lie in the timelike
// static domain; callers in an ergoregion must use an invariant observer frame.
[[nodiscard]] inline std::optional<double> StaticObserverFrequencyRatio(double g_tt_emit,
                                                                        double g_tt_observer) {
    if (!std::isfinite(g_tt_emit) || !std::isfinite(g_tt_observer) || !(g_tt_emit < 0.0) ||
        !(g_tt_observer < 0.0)) {
        return std::nullopt;
    }
    return std::sqrt(g_tt_emit / g_tt_observer);
}

}  // namespace sirius::core::relativity
