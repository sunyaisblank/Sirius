#pragma once

// Single construction authority for metrics consumed by the Cartesian CPU
// geodesic tracer. Registry support and concrete construction meet here: a
// metric not advertised for the CPU path returns null, never a substituted
// spacetime. The same parameter domain is consumed by the typed session and
// direct factory callers.

#include "sirius/core/metrics/kerr_schild_family.h"
#include "sirius/core/metrics/morris_thorne_family.h"
#include "sirius/core/metrics/registry.h"
#include "sirius/core/metrics/warp_drive_family.h"

#include <cmath>
#include <memory>
#include <optional>
#include <string_view>

namespace sirius::core {

struct MetricConstructionParameters {
    double mass = 1.0;
    double dimensionless_spin = 0.0;
    double dimensionless_charge = 0.0;
    double cosmological_constant = 0.0;
    double throat_radius = kDefaultMorrisThorneThroatRadius;
    WormholeTopology wormhole_topology = WormholeTopology::OneSheetCapture;
    double warp_velocity = kDefaultAlcubierreWarpVelocity;
    double bubble_radius = kDefaultAlcubierreBubbleRadius;
    double bubble_sigma = kDefaultAlcubierreBubbleSigma;
};

// Complete parameter domain for a concrete CPU metric request. Keeping this
// beside construction lets the typed session and direct factory callers share
// one fail-closed boundary; no field absent from the selected metric can be
// silently ignored.
[[nodiscard]] inline std::optional<std::string_view> MetricConstructionParameterIssue(
    MetricId id, const MetricConstructionParameters& parameters) noexcept {
    if (const auto issue =
            MetricParameterIssue(id, parameters.dimensionless_spin, parameters.dimensionless_charge,
                                 parameters.cosmological_constant);
        issue.has_value()) {
        return issue;
    }
    if (!std::isfinite(parameters.mass) || parameters.mass < 0.0 || parameters.mass > 100.0 ||
        !std::isfinite(parameters.dimensionless_spin) || parameters.dimensionless_spin < 0.0 ||
        parameters.dimensionless_spin > 0.998 || !std::isfinite(parameters.dimensionless_charge) ||
        parameters.dimensionless_charge < 0.0 || parameters.dimensionless_charge > 0.999 ||
        !std::isfinite(parameters.cosmological_constant) ||
        std::abs(parameters.cosmological_constant) > 0.1) {
        return "metric mass, spin, charge, or lambda is outside the represented domain";
    }
    if (parameters.dimensionless_spin * parameters.dimensionless_spin +
            parameters.dimensionless_charge * parameters.dimensionless_charge >=
        0.999) {
        return "combined spin squared plus charge squared must be below 0.999";
    }
    const bool uses_mass = MetricUsesMass(id);
    if ((!uses_mass && parameters.mass != 0.0) || (uses_mass && parameters.mass < 0.1)) {
        return uses_mass ? "mass must be between 0.1 and 100 for this metric"
                         : "metrics without a mass parameter require mass to be zero";
    }
    if (const auto issue =
            MetricHorizonIssue(id, parameters.mass, parameters.cosmological_constant);
        issue.has_value()) {
        return issue;
    }
    if (!std::isfinite(parameters.throat_radius) ||
        parameters.throat_radius < kMinMorrisThorneThroatRadius ||
        parameters.throat_radius > kMaxMorrisThorneThroatRadius ||
        !std::isfinite(parameters.warp_velocity) || std::abs(parameters.warp_velocity) > 10.0 ||
        !std::isfinite(parameters.bubble_radius) || parameters.bubble_radius <= 0.0 ||
        parameters.bubble_radius > 1000.0 || !std::isfinite(parameters.bubble_sigma) ||
        parameters.bubble_sigma <= 0.0 || parameters.bubble_sigma > 1000.0) {
        return "wormhole or warp-drive parameters are outside the represented domain";
    }
    if (const auto issue = MetricSpecificParameterIssue(
            id, parameters.throat_radius, parameters.wormhole_topology, parameters.warp_velocity,
            parameters.bubble_radius, parameters.bubble_sigma);
        issue.has_value()) {
        return issue;
    }
    if (id == MetricId::Alcubierre) {
        return AlcubierreScaleIssue(parameters.bubble_radius, parameters.bubble_sigma);
    }
    return std::nullopt;
}

[[nodiscard]] inline std::unique_ptr<IMetric> CreateCpuMetric(
    MetricId id, const MetricConstructionParameters& parameters) {
    bool cpu_supported = false;
    for (const auto& info : MetricRegistry()) {
        if (info.id == id) {
            cpu_supported = info.cpu_supported;
            break;
        }
    }
    if (!cpu_supported || MetricConstructionParameterIssue(id, parameters).has_value()) {
        return nullptr;
    }

    switch (id) {
        case MetricId::Minkowski:
            return std::make_unique<KerrSchildFamily>(KerrSchildParams::Minkowski());
        case MetricId::Schwarzschild:
            return std::make_unique<KerrSchildFamily>(
                KerrSchildParams::Schwarzschild(parameters.mass));
        case MetricId::Kerr:
            return std::make_unique<KerrSchildFamily>(KerrSchildParams::Kerr(
                parameters.mass, parameters.dimensionless_spin * parameters.mass));
        case MetricId::ReissnerNordstrom:
            return std::make_unique<KerrSchildFamily>(KerrSchildParams::ReissnerNordstrom(
                parameters.mass, parameters.dimensionless_charge * parameters.mass));
        case MetricId::KerrNewman:
            return std::make_unique<KerrSchildFamily>(KerrSchildParams::KerrNewman(
                parameters.mass, parameters.dimensionless_spin * parameters.mass,
                parameters.dimensionless_charge * parameters.mass));
        case MetricId::DeSitter:
            return std::make_unique<KerrSchildFamily>(
                KerrSchildParams::DeSitter(parameters.cosmological_constant));
        case MetricId::SchwarzschildDeSitter: {
            KerrSchildParams kottler = KerrSchildParams::Schwarzschild(parameters.mass);
            kottler.Lambda = parameters.cosmological_constant;
            return std::make_unique<KerrSchildFamily>(kottler);
        }
        case MetricId::MorrisThorne:
            return std::make_unique<MorrisThorneCartesian>(
                MorrisThorneParams::Ellis(parameters.throat_radius));
        case MetricId::Alcubierre: {
            WarpDriveParams warp;
            warp.vs = parameters.warp_velocity;
            warp.R = parameters.bubble_radius;
            warp.sigma = parameters.bubble_sigma;
            return std::make_unique<WarpDriveFamily>(warp);
        }
    }
    return nullptr;
}

}  // namespace sirius::core
