#pragma once

// Invariant frequency-transfer primitives for stationary, axisymmetric
// spacetimes.  The render path integrates in Kerr-Schild Cartesian
// coordinates, while circular disk motion is most compactly specified by the
// Boyer-Lindquist Killing fields.  Energy E=-p_t and axial angular momentum
// L_z=p_phi are invariant under the stationary/axisymmetric chart change, so
// they form the exact bridge without introducing a Euclidean orbital speed.
//
// For an equatorial circular observer u = u^t(\partial_t + Omega \partial_phi),
// the emitted photon frequency is
//
//   nu = -p.u = u^t(E - Omega L_z),
//
// and g=nu_observer/nu_emitter.  The locally non-rotating (ZAMO) frame supplies
// the physically timelike gravitational/frame-dragging-only branch used as an
// explicitly nonphysical diagnostic substitution when orbital Doppler transfer
// is disabled; unlike a static observer it remains valid inside the Kerr
// ergosphere.

#include "sirius/base/contracts.h"
#include "sirius/core/tensor.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <optional>

namespace sirius::core::relativity {

struct EquatorialCircularFrame {
    double angular_velocity = 0.0;
    double time_component = 0.0;
};

struct KerrDiskFrequencyTransfer {
    double full_g = 0.0;
    double zamo_g = 0.0;
    double emitter_frequency = 0.0;
    double zamo_frequency = 0.0;
    EquatorialCircularFrame emitter;
    EquatorialCircularFrame zamo;
};

struct GreyTransferState {
    std::array<double, 3> observed_emission{};
    double optical_depth = 0.0;
};

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
    const double layer_weight = std::exp(-state.optical_depth) * (1.0 - std::exp(-accepted_tau));
    for (std::size_t channel = 0; channel < state.observed_emission.size(); ++channel) {
        state.observed_emission[channel] += observed_source[channel] * layer_weight;
    }
    state.optical_depth += accepted_tau;
    return accepted_fraction;
}

[[nodiscard]] inline std::optional<EquatorialCircularFrame> EquatorialFrame(
    double g_tt, double g_t_phi, double g_phi_phi, double angular_velocity) {
    if (!std::isfinite(g_tt) || !std::isfinite(g_t_phi) || !std::isfinite(g_phi_phi) ||
        !std::isfinite(angular_velocity) || !(g_phi_phi > 0.0)) {
        return std::nullopt;
    }
    const double norm =
        g_tt + 2.0 * angular_velocity * g_t_phi + angular_velocity * angular_velocity * g_phi_phi;
    if (!std::isfinite(norm) || !(norm < 0.0)) return std::nullopt;
    return EquatorialCircularFrame{angular_velocity, 1.0 / std::sqrt(-norm)};
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

    const double g_tt = -(1.0 - 2.0 * mass / radius);
    const double g_t_phi = -2.0 * mass * spin / radius;
    const double g_phi_phi = radius * radius + spin * spin + 2.0 * mass * spin * spin / radius;
    const double sqrt_mass = std::sqrt(mass);
    const double emitter_omega = sqrt_mass / (std::pow(radius, 1.5) + spin * sqrt_mass);
    const double zamo_omega = -g_t_phi / g_phi_phi;
    const auto emitter = EquatorialFrame(g_tt, g_t_phi, g_phi_phi, emitter_omega);
    const auto zamo = EquatorialFrame(g_tt, g_t_phi, g_phi_phi, zamo_omega);
    if (!emitter || !zamo) return std::nullopt;

    const double emitter_frequency =
        emitter->time_component *
        (photon_energy - emitter->angular_velocity * photon_angular_momentum);
    const double zamo_frequency =
        zamo->time_component * (photon_energy - zamo->angular_velocity * photon_angular_momentum);
    if (!std::isfinite(emitter_frequency) || !(emitter_frequency > 0.0) ||
        !std::isfinite(zamo_frequency) || !(zamo_frequency > 0.0)) {
        return std::nullopt;
    }

    return KerrDiskFrequencyTransfer{observer_frequency / emitter_frequency,
                                     observer_frequency / zamo_frequency,
                                     emitter_frequency,
                                     zamo_frequency,
                                     *emitter,
                                     *zamo};
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
