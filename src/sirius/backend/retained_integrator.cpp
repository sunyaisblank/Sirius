#include "sirius/backend/retained_integrator.h"

#include "sirius/core/twofold.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <utility>

namespace sirius::backend {
namespace {
constexpr double kInvalid = std::numeric_limits<double>::infinity();
using core::CoupledStepFailure;

void Observe(RetainedErrorObservation& observation, double ratio, std::size_t field,
             bool evaluated) {
    if (!observation.observed || ratio > observation.ratio)
        observation = {ratio, field, true, evaluated};
}

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

bool ExactZero(const RetainedValue& value) {
    return value.IsRepresented() && value.high == 0 && value.low == 0 && value.tail == 0 &&
           value.radius == 0;
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

std::array<RetainedValue, 40> PhaseIncrements(const RetainedStepOutput& output, bool lower) {
    auto result = output.increment;
    if (!lower) return result;
    for (std::size_t i = 0; i < result.size(); ++i) {
        CenterDifference difference(output.increment[i], output.error[i]);
        const float high = difference.Extract(), low = difference.Extract(),
                    tail = difference.Extract();
        const double radius = difference.Radius(output.increment[i].radius, output.error[i].radius);
        const float upper = radius == 0 ? 0
                                        : std::nextafter(static_cast<float>(radius),
                                                         std::numeric_limits<float>::infinity());
        result[i] = {high, low, tail, upper, 1};
    }
    return result;
}

RetainedDopriPhaseInput PhasePacket(const RetainedEndpointOutput& start,
                                    const RetainedStepOutput& step, double interval, bool lower) {
    RetainedDopriPhaseInput packet;
    std::copy(start.phase.begin(), start.phase.end(), packet.values.begin());
    const auto increments = PhaseIncrements(step, lower);
    std::copy(increments.begin(), increments.end(), packet.values.begin() + 40);
    for (std::size_t stage = 0; stage < 7; ++stage)
        std::copy(step.rhs[stage].begin(), step.rhs[stage].end(),
                  packet.values.begin() + 80 + stage * 40);
    packet.values[360] = RetainedValue::FromDouble(interval);
    packet.values[361] = RetainedValue::FromDouble(1);
    return packet;
}

core::DopriPositionSegment PositionCurve(const RetainedDopriPhaseInput& input,
                                         const RetainedDopriPhaseOutput& output) {
    core::DopriPositionSegment curve;
    curve.interval = Center(input.values[360]).Rounded();
    for (int i = 0; i < 4; ++i) {
        curve.origin(i) = Center(input.values[i]).Rounded();
        curve.increment(i) = Center(input.values[40 + i]).Rounded();
        curve.a(i) = Center(output.a[i]).Rounded();
        curve.b(i) = Center(output.b[i]).Rounded();
        curve.c(i) = Center(output.c[i]).Rounded();
    }
    return curve;
}
}  // namespace

double RetainedPhysicalError(const std::array<RetainedValue, 40>& first,
                             const std::array<RetainedValue, 40>& second,
                             const RetainedIntervalControl& control, std::size_t* limiting_field) {
    if (limiting_field) *limiting_field = 41;
    if (!ValidControl(control)) return kInvalid;
    double sum = 0, ratio = 0;
    std::size_t field = 40;
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
        else {
            if (error > ratio) field = i;
            ratio = std::max(ratio, error);
        }
    }
    const double central = std::sqrt(sum / 8);
    if (limiting_field && std::isfinite(central)) *limiting_field = central > ratio ? 40 : field;
    return std::max(ratio, central);
}

namespace {
base::Expected<std::vector<RetainedIntervalOutput>> AttemptIntervals(
    RetainedCompute& compute, std::span<const RetainedIntervalInput> inputs,
    std::size_t projection_row_budget, bool dopri) {
    if (inputs.empty() || inputs.size() > compute.Capacity())
        return base::Fail(base::ErrorDomain::kDevice, "attempt retained intervals",
                          "invalid batch size");
    if (projection_row_budget == 0 || projection_row_budget > compute.Capacity() ||
        inputs.size() > projection_row_budget)
        return base::Fail(base::ErrorDomain::kDevice, "attempt retained intervals",
                          "invalid projection row budget");
    const auto count = inputs.size();
    // Keep every row position, including inactive rows. A paired dispatch uses
    // the unchanged endpoint shader and capacity layout, inside the current
    // governor budget; larger batches retain two synchronous submissions.
    const bool paired_projections = count <= projection_row_budget / 2;
    std::vector<RetainedIntervalOutput> work(count);
    std::vector<std::shared_ptr<RetainedDopriInterval>> curves(count);
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

    const auto dense_packets = [&] {
        std::vector<RetainedDenseInput> packets(count);
        for (std::size_t row = 0; row < count; ++row)
            if (active[row]) {
                auto& values = packets[row].values;
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
        return packets;
    };
    std::optional<RetainedEndpointDenseOutput> shared;
    std::vector<bool> speculative_dense;

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
        std::vector<RetainedEndpointInput> upper(paired_projections ? 2 * count : count),
            lower(paired_projections ? 0 : count);
        for (std::size_t row = 0; row < count; ++row)
            if (active[row]) {
                auto& state = work[row];
                const auto& step = (*stepped)[row];
                state.attempted_stages += step.stages;
                const double error =
                    step.valid ? EmbeddedError(packets[row], step, inputs[row].control.integrator)
                               : kInvalid;
                Observe(state.error_checks[0], error, std::isfinite(error) ? 40 : 41, step.valid);
                state.error_ratio = std::max(state.error_ratio, error);
                if (error > 1) {
                    active[row] = false;
                    state.failure = CoupledStepFailure::DerivativeDomain;
                    continue;
                }
                upper[row] = EndpointInput(inputs[row], step.fifth);
                if (paired_projections)
                    upper[count + row] = EndpointInput(inputs[row], step.fourth);
                else
                    lower[row] = EndpointInput(inputs[row], step.fourth);
                if (dopri) {
                    const bool flat = ExactZero(inputs[row].metric[0]) &&
                                      ExactZero(inputs[row].metric[2]) &&
                                      ExactZero(inputs[row].metric[3]);
                    if (!flat && !step.rhs_valid) {
                        active[row] = false;
                        state.failure = CoupledStepFailure::Interpolation;
                        continue;
                    }
                    if (!flat) {
                        if (part == 0) {
                            curves[row] = std::make_shared<RetainedDopriInterval>();
                            curves[row]->metric = inputs[row].metric;
                            curves[row]->chart = inputs[row].chart;
                            curves[row]->control = inputs[row].control;
                        }
                        auto& curve = *curves[row];
                        const auto& start = part == 2 ? state.midpoint : inputs[row].start;
                        const double h = inputs[row].interval / (part == 0 ? 1 : 2);
                        const auto trial = part == 0 ? 0 : part + 1;
                        curve.starts[trial] = start;
                        curve.packets[trial] = PhasePacket(start, step, h, false);
                        if (part == 0) {
                            curve.starts[1] = start;
                            curve.packets[1] = PhasePacket(start, step, h, true);
                        }
                    }
                }
            }
        base::Expected<std::vector<RetainedEndpointOutput>> projected;
        // Dense uses only the already projected full state and original start.
        // Keep budget-one retry completely serialized; a shared wait is reducible.
        const bool legacy_dense = !dopri || std::none_of(curves.begin(), curves.end(),
                                                         [](const auto& p) { return bool(p); });
        if (legacy_dense && part == 2 && paired_projections && projection_row_budget > 1 &&
            compute.SupportsIndependentPair()) {
            const auto dense_inputs = dense_packets();
            auto paired = compute.EndpointAndDense(upper, dense_inputs);
            if (!paired) return std::unexpected(paired.error());
            shared = std::move(*paired);
            projected = std::move(shared->endpoints);
            speculative_dense = active;
        } else {
            projected = compute.Endpoint(upper);
        }
        if (!projected) return std::unexpected(projected.error());
        std::vector<RetainedEndpointOutput> projected_lower;
        if (!paired_projections) {
            auto separate_lower = compute.Endpoint(lower);
            if (!separate_lower) return std::unexpected(separate_lower.error());
            projected_lower = std::move(*separate_lower);
        }
        for (std::size_t row = 0; row < count; ++row)
            if (active[row]) {
                auto& state = work[row];
                const auto& high = (*projected)[row];
                const auto& low =
                    paired_projections ? (*projected)[count + row] : projected_lower[row];
                std::size_t field = 41;
                const bool evaluated = high.valid && low.valid && high.component == low.component;
                const double error = evaluated ? RetainedPhysicalError(high.physical, low.physical,
                                                                       inputs[row].control, &field)
                                               : kInvalid;
                Observe(state.error_checks[1], error, field, evaluated);
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
                if (curves[row]) {
                    const auto trial = part == 0 ? 0 : part + 1;
                    curves[row]->endpoints[trial] = high;
                    if (part == 0) curves[row]->endpoints[1] = low;
                }
            }
    }

    for (auto& state : work) state.embedded_projected_error_ratio = state.error_ratio;
    std::vector<RetainedDenseOutput> dense(count);
    if (shared) {
        for (std::size_t row = 0; row < count; ++row) {
            // Ignore only a speculative result whose final projection rejected.
            // Actual transfer failures and malformed original zero-input rows
            // keep their original global-error meaning.
            if (speculative_dense[row] && !active[row]) continue;
            if (!shared->dense[row]) return std::unexpected(shared->dense[row].error());
            dense[row] = *shared->dense[row];
        }
    } else {
        // Only the exact flat fast path retains Hermite sampling in the DP
        // route. Invalid placeholders preserve original row positions.
        auto packets = dense_packets();
        bool needed = !dopri;
        if (dopri)
            for (std::size_t row = 0; row < count; ++row) {
                if (curves[row]) packets[row] = {};
                if (active[row] && !curves[row]) needed = true;
            }
        if (needed) {
            const auto evaluated = compute.Dense(packets);
            if (!evaluated) return std::unexpected(evaluated.error());
            dense = std::move(*evaluated);
        }
    }
    if (dopri) {
        std::vector<RetainedEndpointInput> midpoint_packets(count);
        bool needed = false;
        // All four private trial packets already exist. Preserve trial-major
        // row positions while filling each current governor-sized dispatch;
        // budget one retains the original serialized evaluation. No packet
        // depends on another trial's interpolation result.
        const auto total = core::kCoupledTrialCount * count;
        for (std::size_t offset = 0; offset < total;) {
            bool trial_needed = false;
            for (std::size_t row = 0; row < count; ++row)
                trial_needed = trial_needed || (active[row] && curves[row]);
            if (!trial_needed) break;
            const auto rows = std::min(projection_row_budget, total - offset);
            std::vector<RetainedDopriPhaseInput> packets(rows);
            for (std::size_t item = 0; item < rows; ++item) {
                const auto trial = (offset + item) / count;
                const auto row = (offset + item) % count;
                if (active[row] && curves[row]) {
                    packets[item] = curves[row]->packets[trial];
                    if (trial == 0) packets[item].values[361] = RetainedValue::FromDouble(.5);
                }
            }
            const auto sampled = compute.DopriPhase(packets);
            // Performed transfer/malformed-result failures stay global, even
            // if an earlier private trial in this dispatch refused its row.
            if (!sampled) return std::unexpected(sampled.error());
            for (std::size_t item = 0; item < rows; ++item) {
                const auto trial = (offset + item) / count;
                const auto row = (offset + item) % count;
                if (active[row] && curves[row]) {
                    auto& curve = *curves[row];
                    curve.positions[trial] = PositionCurve(packets[item], (*sampled)[item]);
                    if (!(*sampled)[item].valid || !curve.positions[trial].IsFinite()) {
                        active[row] = false;
                        work[row].failure = CoupledStepFailure::Interpolation;
                        continue;
                    }
                    curve.bases[trial] = (*sampled)[item].basis;
                    if (trial == 0) {
                        needed = true;
                        midpoint_packets[row] = EndpointInput(inputs[row], (*sampled)[item].phase);
                    }
                }
            }
            offset += rows;
        }
        if (needed) {
            const auto sampled = compute.Endpoint(midpoint_packets);
            if (!sampled) return std::unexpected(sampled.error());
            for (std::size_t row = 0; row < count; ++row)
                if (active[row] && curves[row]) {
                    dense[row].physical = (*sampled)[row].physical;
                    dense[row].valid = (*sampled)[row].valid;
                }
        }
    }
    for (std::size_t row = 0; row < count; ++row) {
        auto& state = work[row];
        if (active[row]) {
            std::size_t midpoint_field = 41, refined_field = 41;
            const double midpoint_error =
                dense[row].valid
                    ? RetainedPhysicalError(dense[row].physical, state.midpoint.physical,
                                            inputs[row].control, &midpoint_field)
                    : kInvalid;
            const double refined_error =
                dense[row].valid
                    ? RetainedPhysicalError(state.full.physical, state.refined.physical,
                                            inputs[row].control, &refined_field)
                    : kInvalid;
            Observe(state.error_checks[2], midpoint_error, midpoint_field, dense[row].valid);
            Observe(state.error_checks[3], refined_error, refined_field, dense[row].valid);
            const double error =
                dense[row].valid ? std::max(midpoint_error, refined_error) : kInvalid;
            state.error_ratio = std::max(state.error_ratio, error);
            state.admissible = state.error_ratio <= 1;
            for (const auto* increment : {&state.full_increment, &state.lower_increment,
                                          &state.midpoint_increment, &state.refined_increment})
                for (const auto& value : *increment)
                    state.admissible = state.admissible && value.IsRepresented();
            if (!state.admissible) state.failure = CoupledStepFailure::Interpolation;
            if (state.admissible) state.dopri = std::move(curves[row]);
        }
        if (!state.admissible) {
            // No earlier successful substage can escape a rejected attempt.
            const auto stages = state.attempted_stages;
            const auto error = state.error_ratio;
            const auto embedded_projected_error = state.embedded_projected_error_ratio;
            const auto error_checks = state.error_checks;
            const auto failure = state.failure;
            state = {};
            state.attempted_stages = stages;
            state.error_ratio = error;
            state.embedded_projected_error_ratio = embedded_projected_error;
            state.error_checks = error_checks;
            state.failure = failure;
        }
    }
    return work;
}
}  // namespace

base::Expected<std::vector<RetainedIntervalOutput>> AttemptRetainedIntervals(
    RetainedCompute& compute, std::span<const RetainedIntervalInput> inputs) {
    return AttemptRetainedIntervals(compute, inputs, compute.Capacity());
}

base::Expected<std::vector<RetainedIntervalOutput>> AttemptRetainedIntervals(
    RetainedCompute& compute, std::span<const RetainedIntervalInput> inputs,
    std::size_t projection_row_budget) {
    return AttemptIntervals(compute, inputs, projection_row_budget, false);
}

base::Expected<std::vector<RetainedIntervalOutput>> AttemptRetainedDopriIntervals(
    RetainedCompute& compute, std::span<const RetainedIntervalInput> inputs) {
    return AttemptRetainedDopriIntervals(compute, inputs, compute.Capacity());
}

base::Expected<std::vector<RetainedIntervalOutput>> AttemptRetainedDopriIntervals(
    RetainedCompute& compute, std::span<const RetainedIntervalInput> inputs,
    std::size_t projection_row_budget) {
    return AttemptIntervals(compute, inputs, projection_row_budget, true);
}

namespace {
RetainedDenseInput ArrivalPacket(const RetainedDopriInterval& interval,
                                 const std::array<RetainedValue, 40>& physical,
                                 const core::Vec4& normal) {
    RetainedDenseInput packet;
    std::copy(interval.metric.begin(), interval.metric.end(), packet.values.begin());
    // At exact fraction zero Dense selects the supplied physical endpoint and
    // covariant V. Only its existing retained arrival transformation is used;
    // these identical endpoints do not refit the accepted quartic.
    std::copy(physical.begin(), physical.end(), packet.values.begin() + 4);
    std::copy(physical.begin(), physical.end(), packet.values.begin() + 44);
    for (std::size_t i = 84; i < 106; ++i) packet.values[i] = RetainedValue::FromDouble(0);
    packet.values[104] = RetainedValue::FromDouble(1);
    for (int axis = 0; axis < 4; ++axis)
        packet.values[106 + axis] = RetainedValue::FromDouble(normal(axis));
    packet.values[110] = RetainedValue::FromDouble(1);
    packet.values[111] = RetainedValue::FromDouble(interval.chart);
    return packet;
}

bool ArrivalDenominators(const RetainedDopriSampleOutput& sample, const core::Vec4& normal,
                         core::Twofold& physical, core::Twofold& polynomial) {
    double physical_sum = 0, polynomial_sum = 0;
    for (int axis = 0; axis < 4; ++axis) {
        if (!std::isfinite(normal(axis))) return false;
        const auto k = Center(sample.physical[4 + axis]) * normal(axis);
        const auto w = Center(sample.polynomial_tangent[axis]) * normal(axis);
        physical += k;
        polynomial += w;
        physical_sum += std::abs(k.Rounded());
        polynomial_sum += std::abs(w.Rounded());
    }
    const double k = physical.Rounded(), w = polynomial.Rounded();
    constexpr double cancellation = 256 * std::numeric_limits<double>::epsilon();
    return std::isfinite(k) && std::isfinite(w) && std::isfinite(physical_sum) &&
           std::isfinite(polynomial_sum) && std::abs(k) > cancellation * physical_sum &&
           std::abs(w) > cancellation * polynomial_sum && std::signbit(k) == std::signbit(w);
}

bool ArrivalAgreement(const RetainedDopriSampleOutput& fixed, const RetainedDenseOutput& physical,
                      const RetainedDenseOutput& polynomial, const RetainedDopriInterval& interval,
                      const core::Vec4& normal) {
    if (!physical.valid || !polynomial.valid) return false;
    auto comparison = physical.physical;
    for (std::size_t column = 0; column < 4; ++column)
        for (std::size_t axis = 0; axis < 4; ++axis)
            comparison[8 + 8 * column + axis] = polynomial.physical[8 + 8 * column + axis];
    if (RetainedPhysicalError(physical.physical, comparison, interval.control) > 1) return false;
    core::Twofold k_normal, w_normal;
    if (!ArrivalDenominators(fixed, normal, k_normal, w_normal)) return false;
    for (std::size_t column = 0; column < 4; ++column) {
        core::Twofold numerator;
        for (int axis = 0; axis < 4; ++axis)
            numerator += Center(fixed.physical[8 + 8 * column + axis]) * normal(axis);
        const auto difference = -numerator / k_normal + numerator / w_normal;
        for (std::size_t axis = 0; axis < 4; ++axis) {
            const auto index = 8 + 8 * column + axis;
            const double magnitude = std::max(std::abs(Center(physical.physical[index]).Rounded()),
                                              std::abs(Center(comparison[index]).Rounded()));
            const double budget =
                interval.control.tolerance *
                (interval.control.length_scale * interval.control.column_scale[column] + magnitude);
            const double shift =
                std::abs((Center(fixed.physical[4 + axis]) * difference).Rounded());
            if (!std::isfinite(budget) || budget <= 0 || !std::isfinite(shift) || shift > budget)
                return false;
        }
    }
    return true;
}
}  // namespace

base::Expected<std::vector<RetainedDopriSampleOutput>> SampleRetainedDopriIntervals(
    RetainedCompute& compute, std::span<const RetainedDopriSampleInput> inputs,
    std::size_t row_budget) {
    if (inputs.empty() || inputs.size() > compute.Capacity() || row_budget == 0 ||
        row_budget > compute.Capacity() || inputs.size() > row_budget)
        return base::Fail(base::ErrorDomain::kDevice, "sample retained DP intervals",
                          "invalid batch size or row budget");
    const auto count = inputs.size();
    std::vector<RetainedDopriSampleOutput> work(count);
    std::vector<RetainedDopriPhaseInput> packets(count);
    std::vector<std::shared_ptr<const RetainedDopriBasis>> bases(count);
    std::vector<bool> interior(count, false);
    bool sample_interior = false;
    for (std::size_t row = 0; row < count; ++row) {
        const auto& input = inputs[row];
        if (!input.interval || input.trial >= 4 || !std::isfinite(input.fraction) ||
            input.fraction < 0 || input.fraction > 1 || !ValidControl(input.interval->control))
            continue;
        const auto& curve = *input.interval;
        if (input.fraction == 0 || input.fraction == 1) {
            const auto& endpoint =
                input.fraction == 0 ? curve.starts[input.trial] : curve.endpoints[input.trial];
            if (!endpoint.valid) continue;
            work[row].physical = endpoint.physical;
            const auto stage = input.fraction == 0 ? 0 : 6;
            for (std::size_t axis = 0; axis < 4; ++axis)
                work[row].polynomial_tangent[axis] =
                    curve.packets[input.trial].values[80 + 40 * stage + axis];
            work[row].valid = true;
        } else {
            sample_interior = interior[row] = true;
            packets[row] = curve.packets[input.trial];
            packets[row].values[361] = RetainedValue::FromDouble(input.fraction);
            bases[row] = curve.bases[input.trial];
        }
    }
    if (sample_interior) {
        const auto phases = compute.DopriPhaseFromBasis(packets, bases);
        if (!phases) return std::unexpected(phases.error());
        std::vector<RetainedEndpointInput> endpoints(count);
        for (std::size_t row = 0; row < count; ++row)
            if (interior[row] && (*phases)[row].valid) {
                const auto& interval = *inputs[row].interval;
                std::copy(interval.metric.begin(), interval.metric.end(),
                          endpoints[row].values.begin());
                std::copy((*phases)[row].phase.begin(), (*phases)[row].phase.end(),
                          endpoints[row].values.begin() + 4);
                endpoints[row].values[44] = RetainedValue::FromDouble(interval.chart);
                std::copy_n((*phases)[row].derivative.begin(), 4,
                            work[row].polynomial_tangent.begin());
            }
        const auto projected = compute.Endpoint(endpoints);
        if (!projected) return std::unexpected(projected.error());
        for (std::size_t row = 0; row < count; ++row)
            if (interior[row] && (*phases)[row].valid && (*projected)[row].valid) {
                work[row].physical = (*projected)[row].physical;
                work[row].valid = true;
            }
    }
    std::vector<RetainedDenseInput> physical_packets(count), polynomial_packets(count);
    std::vector<bool> moving(count, false);
    bool sample_arrival = false;
    for (std::size_t row = 0; row < count; ++row)
        if (work[row].valid && inputs[row].normal) {
            core::Twofold k_normal, w_normal;
            if (!ArrivalDenominators(work[row], *inputs[row].normal, k_normal, w_normal)) {
                work[row] = {};
                continue;
            }
            sample_arrival = moving[row] = true;
            physical_packets[row] =
                ArrivalPacket(*inputs[row].interval, work[row].physical, *inputs[row].normal);
            auto polynomial = work[row].physical;
            std::copy(work[row].polynomial_tangent.begin(), work[row].polynomial_tangent.end(),
                      polynomial.begin() + 4);
            polynomial_packets[row] =
                ArrivalPacket(*inputs[row].interval, polynomial, *inputs[row].normal);
        }
    if (sample_arrival) {
        const auto physical = compute.Dense(physical_packets);
        if (!physical) return std::unexpected(physical.error());
        const auto polynomial = compute.Dense(polynomial_packets);
        if (!polynomial) return std::unexpected(polynomial.error());
        for (std::size_t row = 0; row < count; ++row)
            if (moving[row]) {
                if (ArrivalAgreement(work[row], (*physical)[row], (*polynomial)[row],
                                     *inputs[row].interval, *inputs[row].normal))
                    work[row].physical = (*physical)[row].physical;
                else
                    work[row] = {};
            }
    }
    for (auto& output : work) {
        if (output.valid) {
            for (const auto& value : output.physical) output.valid &= value.IsRepresented();
            for (const auto& value : output.polynomial_tangent)
                output.valid &= value.IsRepresented();
        }
        if (!output.valid) output = {};
    }
    return work;
}

}  // namespace sirius::backend
