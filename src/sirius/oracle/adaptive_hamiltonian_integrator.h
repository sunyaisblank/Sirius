// Adaptive double-precision Hamiltonian geodesic integrator for oracle work.
//
// Dormand-Prince 5(4) advances the canonical system
//   dq/dlambda = dH/dp,  dp/dlambda = -dH/dq
// using the analytic derivatives supplied by IMetricD. Despite its historical
// name, IMetricD::dHdq returns the already-negated canonical momentum rate
// dp/dlambda, as its interface contract states. This path is deliberately
// independent of the live RK45 state/equation implementation and of the
// fixed-step implicit-midpoint/Yoshida symplectic oracle.

#pragma once

#include "sirius/base/contracts.h"
#include "sirius/oracle/metric_interface.h"
#include "sirius/oracle/transport_types.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

namespace sirius::oracle {

class AdaptiveHamiltonianIntegratorD {
  public:
    struct Config {
        double absolute_tolerance = 1.0e-13;
        double relative_tolerance = 1.0e-13;
        double initial_step = 1.0e-2;
        double minimum_step = 1.0e-10;
        double maximum_step = 5.0e-2;
        int maximum_steps = 1000000;
    };

    [[nodiscard]] static bool IsRepresentedConfig(const Config& config) noexcept {
        return std::isfinite(config.absolute_tolerance) && config.absolute_tolerance > 0.0 &&
               std::isfinite(config.relative_tolerance) && config.relative_tolerance > 0.0 &&
               std::isfinite(config.initial_step) && std::isfinite(config.minimum_step) &&
               config.minimum_step > 0.0 && std::isfinite(config.maximum_step) &&
               config.maximum_step >= config.minimum_step &&
               config.initial_step >= config.minimum_step &&
               config.initial_step <= config.maximum_step && config.maximum_steps > 0;
    }

    struct IntegrationResult {
        HamiltonianStateD state;
        double affine_advance = 0.0;
        double maximum_scaled_error = 0.0;
        int accepted_steps = 0;
        int rejected_steps = 0;
        bool complete = false;
    };

    explicit AdaptiveHamiltonianIntegratorD(const IMetricD* metric)
        : AdaptiveHamiltonianIntegratorD(metric, Config()) {}

    AdaptiveHamiltonianIntegratorD(const IMetricD* metric, Config config)
        : metric_(metric), config_(config) {
        SIRIUS_PRE(metric != nullptr);
        SIRIUS_PRE(IsRepresentedConfig(config));
    }

    [[nodiscard]] IntegrationResult Integrate(const HamiltonianStateD& initial,
                                              double affine_span) const {
        SIRIUS_PRE(std::isfinite(affine_span) && affine_span >= 0.0);
        IntegrationResult result;
        result.state = initial;
        if (affine_span == 0.0) {
            result.complete = true;
            return result;
        }

        double step = std::clamp(config_.initial_step, config_.minimum_step, config_.maximum_step);
        while (result.affine_advance < affine_span &&
               result.accepted_steps + result.rejected_steps < config_.maximum_steps) {
            step = std::min(step, affine_span - result.affine_advance);
            if (step < config_.minimum_step) break;

            const Derivative k1 = Evaluate(result.state);
            const Derivative k2 = Evaluate(Advance(result.state, step, {{1.0 / 5.0, &k1}}));
            const Derivative k3 =
                Evaluate(Advance(result.state, step, {{3.0 / 40.0, &k1}, {9.0 / 40.0, &k2}}));
            const Derivative k4 = Evaluate(Advance(
                result.state, step, {{44.0 / 45.0, &k1}, {-56.0 / 15.0, &k2}, {32.0 / 9.0, &k3}}));
            const Derivative k5 = Evaluate(Advance(result.state, step,
                                                   {{19372.0 / 6561.0, &k1},
                                                    {-25360.0 / 2187.0, &k2},
                                                    {64448.0 / 6561.0, &k3},
                                                    {-212.0 / 729.0, &k4}}));
            const Derivative k6 = Evaluate(Advance(result.state, step,
                                                   {{9017.0 / 3168.0, &k1},
                                                    {-355.0 / 33.0, &k2},
                                                    {46732.0 / 5247.0, &k3},
                                                    {49.0 / 176.0, &k4},
                                                    {-5103.0 / 18656.0, &k5}}));
            const Derivative k7 = Evaluate(Advance(result.state, step,
                                                   {{35.0 / 384.0, &k1},
                                                    {500.0 / 1113.0, &k3},
                                                    {125.0 / 192.0, &k4},
                                                    {-2187.0 / 6784.0, &k5},
                                                    {11.0 / 84.0, &k6}}));

            const HamiltonianStateD fifth = Advance(result.state, step,
                                                    {{35.0 / 384.0, &k1},
                                                     {500.0 / 1113.0, &k3},
                                                     {125.0 / 192.0, &k4},
                                                     {-2187.0 / 6784.0, &k5},
                                                     {11.0 / 84.0, &k6}});
            const HamiltonianStateD fourth = Advance(result.state, step,
                                                     {{5179.0 / 57600.0, &k1},
                                                      {7571.0 / 16695.0, &k3},
                                                      {393.0 / 640.0, &k4},
                                                      {-92097.0 / 339200.0, &k5},
                                                      {187.0 / 2100.0, &k6},
                                                      {1.0 / 40.0, &k7}});

            const double scaled_error = ScaledError(result.state, fifth, fourth);
            if (!std::isfinite(scaled_error)) break;
            result.maximum_scaled_error = std::max(result.maximum_scaled_error, scaled_error);
            if (scaled_error <= 1.0) {
                if (!metric_->IsValid(fifth.q)) break;
                result.state = fifth;
                result.affine_advance += step;
                ++result.accepted_steps;
            } else {
                ++result.rejected_steps;
            }

            const double factor = scaled_error == 0.0
                                      ? 5.0
                                      : std::clamp(0.9 * std::pow(scaled_error, -0.2), 0.2, 5.0);
            step = std::clamp(step * factor, config_.minimum_step, config_.maximum_step);
        }

        result.complete = result.affine_advance >= affine_span;
        if (result.complete) result.state.H = metric_->Hamiltonian(result.state.q, result.state.p);
        return result;
    }

  private:
    struct Derivative {
        Vec4d q;
        Vec4d p;
    };

    struct WeightedDerivative {
        double weight;
        const Derivative* derivative;
    };

    const IMetricD* metric_;
    Config config_;

    [[nodiscard]] Derivative Evaluate(const HamiltonianStateD& state) const {
        return Derivative{metric_->dHdp(state.q, state.p), metric_->dHdq(state.q, state.p)};
    }

    [[nodiscard]] static HamiltonianStateD Advance(
        const HamiltonianStateD& state, double step,
        std::initializer_list<WeightedDerivative> terms) {
        HamiltonianStateD advanced = state;
        for (const WeightedDerivative& term : terms) {
            advanced.q += term.derivative->q * (step * term.weight);
            advanced.p += term.derivative->p * (step * term.weight);
        }
        return advanced;
    }

    [[nodiscard]] double ScaledError(const HamiltonianStateD& before,
                                     const HamiltonianStateD& fifth,
                                     const HamiltonianStateD& fourth) const {
        double maximum = 0.0;
        for (int component = 0; component < 4; ++component) {
            const double q_scale =
                config_.absolute_tolerance +
                config_.relative_tolerance *
                    std::max(std::abs(before.q[component]), std::abs(fifth.q[component]));
            const double p_scale =
                config_.absolute_tolerance +
                config_.relative_tolerance *
                    std::max(std::abs(before.p[component]), std::abs(fifth.p[component]));
            maximum =
                std::max(maximum, std::abs(fifth.q[component] - fourth.q[component]) / q_scale);
            maximum =
                std::max(maximum, std::abs(fifth.p[component] - fourth.p[component]) / p_scale);
        }
        return maximum;
    }
};

}  // namespace sirius::oracle
