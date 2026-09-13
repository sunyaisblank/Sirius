#include "sirius/backend/retained_integrator.h"

#include "sirius/core/twofold.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace sirius::backend {
namespace {
constexpr double kInvalid = std::numeric_limits<double>::infinity();
using core::CoupledStepFailure;

core::Twofold Center(const RetainedValue& value) {
    return core::Twofold(value.high) + core::Twofold(value.low) + core::Twofold(value.tail);
}

// Retained limbs can be widely separated in exponent. Accumulate their
// difference exactly before rounding, including when common leading terms
// cancel. Two doubles cannot preserve every such sparse expansion.
class CenterDifference {
  public:
    CenterDifference(const RetainedValue& first, const RetainedValue& second) {
        for (double value : {double(first.high), double(first.low), double(first.tail),
                             -double(second.high), -double(second.low), -double(second.tail)})
            Add(value);
    }

    void Add(double value) {
        std::size_t written = 0;
        for (std::size_t i = 0; i < size_; ++i) {
            const double term = terms_[i];
            const double sum = value + term;
            const double recovered = sum - value;
            const double residual = (value - (sum - recovered)) + (term - recovered);
            if (residual != 0) terms_[written++] = residual;
            value = sum;
        }
        if (value != 0) terms_[written++] = value;
        size_ = written;
    }

    double Rounded() const {
        double result = 0;
        for (std::size_t i = 0; i < size_; ++i) result += terms_[i];
        return result;
    }

    float Extract() {
        const float value = static_cast<float>(Rounded());
        Add(-double(value));
        return value;
    }

    double Radius(double first, double second) const {
        double result = first;
        const auto include = [&result](double value) {
            if (value > 0) result = std::nextafter(result + value, kInvalid);
        };
        include(second);
        for (std::size_t i = 0; i < size_; ++i) include(std::abs(terms_[i]));
        return result;
    }

  private:
    // Six input limbs and at most three extraction subtractions.
    std::array<double, 9> terms_{};
    std::size_t size_ = 0;
};

bool ValidControl(const RetainedIntervalControl& control) {
    if (!core::IsRepresentedIntegratorStepControl(control.integrator)) return false;
    for (double value : {control.length_scale, control.frequency_scale, control.tolerance})
        if (!std::isfinite(value) || value <= 0) return false;
    for (double value : control.column_scale)
        if (!std::isfinite(value) || value <= 0) return false;
    return true;
}

double EmbeddedError(const RetainedStepInput& input, const RetainedStepOutput& output,
                     const core::IntegratorConfig& config) {
    double sum = 0;
    for (std::size_t i = 0; i < 8; ++i) {
        const double scale =
            config.abs_tolerance +
            config.rel_tolerance * std::max(std::abs(Center(input.values[4 + i]).Rounded()),
                                            std::abs(Center(output.fifth[i]).Rounded()));
        const double error = std::abs(Center(output.error[i]).Rounded()) / scale;
        if (!std::isfinite(scale) || scale <= 0 || !std::isfinite(error)) return kInvalid;
        sum += error * error;
    }
    return std::sqrt(sum / 8);
}

RetainedEndpointInput EndpointInput(const RetainedIntervalInput& input,
                                    const std::array<RetainedValue, 40>& phase) {
    RetainedEndpointInput result;
    std::copy(input.metric.begin(), input.metric.end(), result.values.begin());
    std::copy(phase.begin(), phase.end(), result.values.begin() + 4);
    result.values[44] = RetainedValue::FromDouble(input.chart);
    return result;
}

std::array<RetainedValue, 20> Increments(const RetainedStepOutput& output, bool lower) {
    std::array<RetainedValue, 20> result;
    for (std::size_t group = 0; group < 5; ++group)
        for (std::size_t axis = 0; axis < 4; ++axis) {
            const auto index = group * 8 + axis;
            const auto& full = output.increment[index];
            if (!lower) {
                result[group * 4 + axis] = full;
                continue;
            }
            CenterDifference difference(full, output.error[index]);
            const float high = difference.Extract();
            const float low = difference.Extract();
            const float tail = difference.Extract();
            const double radius = difference.Radius(full.radius, output.error[index].radius);
            const float upper = radius == 0
                                    ? 0
                                    : std::nextafter(static_cast<float>(radius),
                                                     std::numeric_limits<float>::infinity());
            result[group * 4 + axis] = {high, low, tail, upper, 1};
        }
    return result;
}
}  // namespace

double RetainedPhysicalError(const std::array<RetainedValue, 40>& first,
                             const std::array<RetainedValue, 40>& second,
                             const RetainedIntervalControl& control) {
    if (!ValidControl(control)) return kInvalid;
    double sum = 0, ratio = 0;
    for (std::size_t i = 0; i < 40; ++i) {
        if (!first[i].IsRepresented() || !second[i].IsRepresented()) return kInvalid;
        const auto a = Center(first[i]), b = Center(second[i]);
        const double magnitude = std::max(std::abs(a.Rounded()), std::abs(b.Rounded()));
        const bool displacement = i % 8 < 4;
        const double unit = displacement ? control.length_scale : control.frequency_scale;
        const double scale =
            i < 8 ? control.integrator.abs_tolerance * unit +
                        control.integrator.rel_tolerance * magnitude
                  : control.tolerance * (unit * control.column_scale[(i - 8) / 8] + magnitude);
        const double error = std::abs(CenterDifference(first[i], second[i]).Rounded()) / scale;
        if (!std::isfinite(scale) || scale <= 0 || !std::isfinite(error)) return kInvalid;
        if (i < 8)
            sum += error * error;
        else
            ratio = std::max(ratio, error);
    }
    return std::max(ratio, std::sqrt(sum / 8));
}

base::Expected<std::vector<RetainedIntervalOutput>> AttemptRetainedIntervals(
    RetainedCompute& compute, std::span<const RetainedIntervalInput> inputs) {
    if (inputs.empty() || inputs.size() > compute.Capacity())
        return base::Fail(base::ErrorDomain::kDevice, "attempt retained intervals",
                          "invalid batch size");
    const auto count = inputs.size();
    std::vector<RetainedIntervalOutput> work(count);
    std::vector<bool> active(count, true);
    for (std::size_t row = 0; row < count; ++row) {
        const auto& input = inputs[row];
        bool valid = ValidControl(input.control) && input.start.valid &&
                     input.start.component < 4 && (input.chart == 1 || input.chart == -1) &&
                     std::isfinite(input.interval) &&
                     input.interval >= input.control.integrator.min_step &&
                     input.interval <= input.control.integrator.max_step;
        for (const auto& value : input.metric) valid = valid && value.IsRepresented();
        for (const auto& value : input.start.phase) valid = valid && value.IsRepresented();
        for (const auto& value : input.start.physical) valid = valid && value.IsRepresented();
        if (!valid) {
            active[row] = false;
            work[row].failure = CoupledStepFailure::InvalidState;
        }
    }

    for (std::size_t part = 0; part < 3; ++part) {
        std::vector<RetainedStepInput> packets(count);
        for (std::size_t row = 0; row < count; ++row)
            if (active[row]) {
                const auto& phase = part == 2 ? work[row].midpoint.phase : inputs[row].start.phase;
                const auto endpoint = EndpointInput(inputs[row], phase);
                std::copy(endpoint.values.begin(), endpoint.values.end(),
                          packets[row].values.begin());
                packets[row].values[45] =
                    RetainedValue::FromDouble(inputs[row].interval / (part == 0 ? 1 : 2));
            }
        const auto stepped = compute.Step(packets);
        if (!stepped) return std::unexpected(stepped.error());
        std::vector<RetainedEndpointInput> upper(count), lower(count);
        for (std::size_t row = 0; row < count; ++row)
            if (active[row]) {
                auto& state = work[row];
                const auto& step = (*stepped)[row];
                state.attempted_stages += step.stages;
                const double error =
                    step.valid ? EmbeddedError(packets[row], step, inputs[row].control.integrator)
                               : kInvalid;
                state.error_ratio = std::max(state.error_ratio, error);
                if (error > 1) {
                    active[row] = false;
                    state.failure = CoupledStepFailure::DerivativeDomain;
                    continue;
                }
                upper[row] = EndpointInput(inputs[row], step.fifth);
                lower[row] = EndpointInput(inputs[row], step.fourth);
            }
        const auto projected = compute.Endpoint(upper);
        if (!projected) return std::unexpected(projected.error());
        const auto projected_lower = compute.Endpoint(lower);
        if (!projected_lower) return std::unexpected(projected_lower.error());
        for (std::size_t row = 0; row < count; ++row)
            if (active[row]) {
                auto& state = work[row];
                const auto& high = (*projected)[row];
                const auto& low = (*projected_lower)[row];
                const double error =
                    high.valid && low.valid && high.component == low.component
                        ? RetainedPhysicalError(high.physical, low.physical, inputs[row].control)
                        : kInvalid;
                state.error_ratio = std::max(state.error_ratio, error);
                if (error > 1) {
                    active[row] = false;
                    state.failure = CoupledStepFailure::Projection;
                    continue;
                }
                if (part == 0) {
                    state.full = high;
                    state.lower = low;
                    state.full_increment = Increments((*stepped)[row], false);
                    state.lower_increment = Increments((*stepped)[row], true);
                } else if (part == 1) {
                    state.midpoint = high;
                    state.midpoint_increment = Increments((*stepped)[row], false);
                } else {
                    state.refined = high;
                    state.refined_increment = Increments((*stepped)[row], false);
                }
            }
    }

    std::vector<RetainedDenseInput> dense_inputs(count);
    for (std::size_t row = 0; row < count; ++row)
        if (active[row]) {
            auto& values = dense_inputs[row].values;
            std::copy(inputs[row].metric.begin(), inputs[row].metric.end(), values.begin());
            std::copy(inputs[row].start.physical.begin(), inputs[row].start.physical.end(),
                      values.begin() + 4);
            std::copy(work[row].full.physical.begin(), work[row].full.physical.end(),
                      values.begin() + 44);
            std::copy(work[row].full_increment.begin(), work[row].full_increment.end(),
                      values.begin() + 84);
            values[104] = RetainedValue::FromDouble(inputs[row].interval);
            values[105] = RetainedValue::FromDouble(.5);
            for (std::size_t i = 106; i < 111; ++i) values[i] = RetainedValue::FromDouble(0);
            values[111] = RetainedValue::FromDouble(inputs[row].chart);
        }
    const auto dense = compute.Dense(dense_inputs);
    if (!dense) return std::unexpected(dense.error());
    for (std::size_t row = 0; row < count; ++row) {
        auto& state = work[row];
        if (active[row]) {
            const double error =
                (*dense)[row].valid
                    ? std::max(RetainedPhysicalError((*dense)[row].physical,
                                                     state.midpoint.physical, inputs[row].control),
                               RetainedPhysicalError(state.full.physical, state.refined.physical,
                                                     inputs[row].control))
                    : kInvalid;
            state.error_ratio = std::max(state.error_ratio, error);
            state.admissible = state.error_ratio <= 1;
            for (const auto* increment : {&state.full_increment, &state.lower_increment,
                                          &state.midpoint_increment, &state.refined_increment})
                for (const auto& value : *increment)
                    state.admissible = state.admissible && value.IsRepresented();
            if (!state.admissible) state.failure = CoupledStepFailure::Interpolation;
        }
        if (!state.admissible) {
            // No earlier successful substage can escape a rejected attempt.
            const auto stages = state.attempted_stages;
            const auto error = state.error_ratio;
            const auto failure = state.failure;
            state = {};
            state.attempted_stages = stages;
            state.error_ratio = error;
            state.failure = failure;
        }
    }
    return work;
}

}  // namespace sirius::backend
