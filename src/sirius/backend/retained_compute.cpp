#include "sirius/backend/retained_compute.h"

#include "retained_kernels.h"

#include <algorithm>
#include <bit>
#include <chrono>
#include <cmath>
#include <cstring>
#include <limits>

namespace sirius::backend {
namespace {
using base::ErrorDomain;
using base::Fail;
using namespace retained_program;
// One ray owns one workgroup; stay within Vulkan's portable X dispatch bound.
constexpr std::size_t kMaximumCapacity = 65535;

struct StageLayout {
    std::span<const std::uint32_t> program;
    std::size_t input_words;
    std::size_t row_words;
};
// Same order as KernelStage; planning and allocation share these exact spans.
constexpr std::array<StageLayout, RetainedCompute::kStageCount> kStageLayouts{{
    {kCameraProgram, 160, kCameraRowWords},
    {kTransportProgram, 230, kTransportRowWords},
    {kEndpointProgram, 225, kEndpointRowWords},
    {kDenseProgram, 560, kDenseRowWords},
    {kInitializeProgram, 225, kInitializeRowWords},
    {kRayCameraProgram, 225, kRayCameraRowWords},
    {kDopriPhaseProgram, 1810, kDopriPhaseRowWords},
}};

std::array<std::uint64_t, 2> StageBufferBytes(const StageLayout& layout, std::size_t capacity) {
    const auto rows = static_cast<std::uint64_t>(capacity);
    return {4 * (1 + rows * layout.input_words + layout.program.size()),
            4 * rows * layout.row_words};
}

RetainedValue Decode(const std::uint32_t* words) {
    return std::bit_cast<RetainedValue>(
        std::array<std::uint32_t, 5>{words[0], words[1], words[2], words[3], words[4]});
}

base::Expected<RetainedEndpointOutput> DecodeEndpoint(const std::uint32_t* words) {
    RetainedEndpointOutput output;
    if (words[0] == 0) return output;
    if (words[0] != 1 || words[1] > 3)
        return Fail(ErrorDomain::kKernel, "read retained endpoint", "invalid projection state");
    for (std::size_t i = 0; i < 40; ++i) {
        output.phase[i] = Decode(words + 104 + i * 5);
        output.physical[i] = Decode(words + 304 + i * 5);
        if (!output.phase[i].IsRepresented() || !output.physical[i].IsRepresented())
            return Fail(ErrorDomain::kKernel, "read retained endpoint", "invalid retained value");
    }
    output.component = words[1];
    output.valid = true;
    return output;
}

base::Expected<RetainedDenseOutput> DecodeDense(const std::uint32_t* words) {
    RetainedDenseOutput output;
    if (words[0] == 0) return output;
    if (words[0] != 1)
        return Fail(ErrorDomain::kKernel, "read retained dense segment",
                    "invalid completion state");
    for (std::size_t i = 0; i < 40; ++i) {
        output.physical[i] = Decode(words + 4 + i * 5);
        if (!output.physical[i].IsRepresented())
            return Fail(ErrorDomain::kKernel, "read retained dense segment",
                        "invalid retained value");
    }
    output.valid = true;
    return output;
}
}  // namespace

RetainedValue RetainedValue::FromDouble(double value) {
    if (!std::isfinite(value) || std::abs(value) > 0x1p120) return {};
    RetainedValue result;
    result.high = static_cast<float>(value);
    const double residual = value - double(result.high);
    result.low = static_cast<float>(residual);
    const double tail = residual - double(result.low);
    result.tail = static_cast<float>(tail);
    const double remainder = std::abs(tail - double(result.tail));
    result.radius = static_cast<float>(remainder);
    if (remainder != 0)
        result.radius = std::nextafter(result.radius, std::numeric_limits<float>::infinity());
    result.valid = 1;
    return result.IsRepresented() ? result : RetainedValue{};
}

bool RetainedValue::IsRepresented() const {
    const auto spacing = [](float value) {
        const float magnitude = std::abs(value);
        return double(std::nextafter(magnitude, std::numeric_limits<float>::infinity())) -
               double(magnitude);
    };
    return valid == 1 && std::isfinite(high) && std::isfinite(low) && std::isfinite(tail) &&
           std::isfinite(radius) && std::abs(high) <= 0x1p120f && radius >= 0 &&
           radius <= 0x1p120f && std::abs(double(low)) <= spacing(high) &&
           std::abs(double(tail)) <= spacing(low) && (high != 0 || (low == 0 && tail == 0));
}

base::Expected<std::unique_ptr<RetainedCompute>> RetainedCompute::Create(
    ComputeDevice& device, std::size_t capacity, bool fp64_products, double dispatch_target_ms) {
    if (capacity == 0 || capacity > kMaximumCapacity)
        return Fail(ErrorDomain::kDevice, "create retained compute stages", "invalid capacity");
    if (!std::isfinite(dispatch_target_ms) || dispatch_target_ms < 0)
        return Fail(ErrorDomain::kDevice, "create retained compute stages",
                    "invalid dispatch target");
    if (const auto issue = RetainedArithmeticIssue(device.Info(), fp64_products))
        return Fail(ErrorDomain::kDevice, "create retained compute stages", std::string(*issue));
    auto result = std::unique_ptr<RetainedCompute>(new RetainedCompute(device, capacity));
    result->dispatch_target_ms_ = dispatch_target_ms;
    const auto create = [&](Stage& stage, KernelStage kind,
                            std::span<const std::uint32_t> code) -> base::Expected<void> {
        auto kernel = device.LoadKernel(code);
        if (!kernel) return std::unexpected(kernel.error());
        const auto& layout = kStageLayouts[static_cast<std::size_t>(kind)];
        const auto bytes = StageBufferBytes(layout, capacity);
        stage.kind = kind;
        stage.kernel = *kernel;
        stage.input.resize(bytes[0] / 4);
        stage.output.resize(bytes[1] / 4);
        stage.input_words_per_row = layout.input_words;
        stage.stats.command_row_counts.resize(capacity + 1);
        stage.input[0] = static_cast<std::uint32_t>(capacity);
        std::copy(layout.program.begin(), layout.program.end(),
                  stage.input.begin() + 1 + layout.input_words * capacity);
        auto input = device.CreateBuffer(bytes[0], BufferUsage::kStorage);
        if (!input) return std::unexpected(input.error());
        auto output = device.CreateBuffer(bytes[1], BufferUsage::kStorage);
        if (!output) return std::unexpected(output.error());
        stage.buffers = {*input, *output};
        return {};
    };
    const bool portable = RetainedUsesPortableArithmetic(device.Info());
    const auto shader = [fp64_products, portable](std::span<const std::uint32_t> narrow,
                                                  std::span<const std::uint32_t> wide,
                                                  std::span<const std::uint32_t> integer_narrow,
                                                  std::span<const std::uint32_t> integer_wide) {
        if (portable) return fp64_products ? integer_wide : integer_narrow;
        return fp64_products ? wide : narrow;
    };
    auto status = create(
        result->camera_, KernelStage::kCamera,
        shader(kCameraShader, kCameraFp64Shader, kCameraPortableShader, kCameraPortableFp64Shader));
    if (!status) return std::unexpected(status.error());
    const bool fma = fp64_products && !portable && device.Info().fma_fp32_enabled &&
                     device.Info().preserves_fp32_signed_zero_inf_nan;
    // The normal-lattice transform needs RTE32 but no denormal preservation.
    // Adapters without RTE32 retain the pure-integer portable modules.
    const bool normal_sum = portable && device.Info().rounds_fp32_to_nearest;
    const std::span<const std::uint32_t> transport =
        normal_sum ? std::span(kTransportPortableNormalSumShader)
                   : (fma && kTransportFmaAvailable
                          ? std::span(kTransportFmaShader)
                          : shader(kTransportShader, kTransportFp64Shader, kTransportPortableShader,
                                   kTransportPortableFp64Shader));
    status = create(result->transport_, KernelStage::kTransport, transport);
    if (!status) return std::unexpected(status.error());
    result->transport_.paired_rows = !portable && !fp64_products;
    const std::span<const std::uint32_t> endpoint =
        normal_sum ? std::span(kEndpointPortableNormalSumShader)
                   : (fma && kEndpointFmaAvailable
                          ? std::span(kEndpointFmaShader)
                          : shader(kEndpointShader, kEndpointFp64Shader, kEndpointPortableShader,
                                   kEndpointPortableFp64Shader));
    status = create(result->endpoint_, KernelStage::kEndpoint, endpoint);
    if (!status) return std::unexpected(status.error());
    status = create(
        result->dense_, KernelStage::kDense,
        shader(kDenseShader, kDenseFp64Shader, kDensePortableShader, kDensePortableFp64Shader));
    if (!status) return std::unexpected(status.error());
    status = create(result->initialize_, KernelStage::kInitialize,
                    shader(kInitializeShader, kInitializeFp64Shader, kInitializePortableShader,
                           kInitializePortableFp64Shader));
    if (!status) return std::unexpected(status.error());
    status = create(result->ray_camera_, KernelStage::kRayCamera,
                    shader(kRayCameraShader, kRayCameraFp64Shader, kRayCameraPortableShader,
                           kRayCameraPortableFp64Shader));
    if (!status) return std::unexpected(status.error());
    // Reuse the same guarded normal-lattice sums without changing the DP
    // program or its integer fallback, including on adapters without RTE32.
    const std::span<const std::uint32_t> dopri_phase =
        normal_sum ? std::span(kDopriPhasePortableNormalSumShader)
                   : shader(kDopriPhaseShader, kDopriPhaseFp64Shader, kDopriPhasePortableShader,
                            kDopriPhasePortableFp64Shader);
    status = create(result->dopri_phase_, KernelStage::kDopriPhase, dopri_phase);
    if (!status) return std::unexpected(status.error());
    return result;
}

base::Expected<void> RetainedCompute::PrepareSoftwareRendererStages(
    PreparationStats& observation, const std::function<bool()>& should_cancel) {
    using Clock = std::chrono::steady_clock;
    observation = {};
    const auto started = Clock::now();
    const auto finish = [&](base::Expected<void> status) {
        observation.wall_ms =
            std::chrono::duration<double, std::milli>(Clock::now() - started).count();
        return status;
    };
    if (device_.Info().kind != DeviceKind::kSoftware)
        return finish(Fail(ErrorDomain::kDevice, "prepare retained renderer",
                           "zero-row initialization is restricted to software devices"));
    const auto cancelled = [&] { return should_cancel && should_cancel(); };
    const auto cancellation = [] {
        return Fail(ErrorDomain::kInternal, "prepare retained renderer", "cancelled");
    };
    const auto valid_timing = [](const DispatchTiming& timing) {
        return std::isfinite(timing.pipeline_setup_ms) && timing.pipeline_setup_ms >= 0 &&
               std::isfinite(timing.command_setup_ms) && timing.command_setup_ms >= 0 &&
               std::isfinite(timing.submit_wait_ms) && timing.submit_wait_ms >= 0 &&
               std::isfinite(timing.cleanup_ms) && timing.cleanup_ms >= 0 &&
               std::isfinite(timing.total_ms) && timing.total_ms >= 0;
    };
    for (auto* stage :
         {&ray_camera_, &initialize_, &transport_, &endpoint_, &dense_, &dopri_phase_}) {
        if (cancelled()) return finish(cancellation());
        auto& stats = observation.stages[static_cast<std::size_t>(stage->kind)];
        ++stats.attempts;
        const auto write_header = [&](std::uint32_t header) {
            ++stats.write_buffer_calls;
            const auto write_started = Clock::now();
            auto status =
                device_.WriteBuffer(stage->buffers[0], std::as_bytes(std::span(&header, 1)));
            stats.write_buffer_ms +=
                std::chrono::duration<double, std::milli>(Clock::now() - write_started).count();
            if (status) stats.write_buffer_bytes += sizeof(header);
            return status;
        };
        // Do not use the physical Dispatch owner: preparation must not upload
        // programs, alter feedback, clear outputs or count completed ray work.
        auto status = write_header(0);
        if (status) {
            if (cancelled()) {
                status = cancellation();
            } else {
                ++stats.dispatch_attempts;
                status = device_.Dispatch(stage->kernel, stage->buffers, 1, 1, 1, &stats.timing);
                if (status) {
                    if (!valid_timing(stats.timing)) {
                        status = Fail(ErrorDomain::kDevice, "prepare retained renderer",
                                      "invalid submission timing");
                    } else {
                        ++stats.completed_dispatches;
                    }
                }
            }
        }
        // Restore even after a failed write/dispatch or cancellation. Preserve
        // the original error if restoration also fails; no later stage executes.
        auto restored = write_header(static_cast<std::uint32_t>(capacity_));
        stats.header_restored = restored.has_value();
        if (stats.completed_dispatches && restored) ++stats.completed;
        if (!status) return finish(std::move(status));
        if (!restored) return finish(std::move(restored));
        if (cancelled()) return finish(cancellation());
    }
    return finish({});
}

std::uint64_t RetainedCompute::RequiredBufferBytes(std::size_t capacity) {
    if (capacity == 0 || capacity > kMaximumCapacity) return 0;
    std::uint64_t total = 0;
    for (const auto& layout : kStageLayouts)
        for (const auto bytes : StageBufferBytes(layout, capacity)) total += bytes;
    return total;
}

base::Expected<std::uint64_t> RetainedCompute::RequiredAllocationBytes(ComputeDevice& device,
                                                                       std::size_t capacity) {
    if (capacity == 0 || capacity > kMaximumCapacity)
        return Fail(ErrorDomain::kDevice, "plan retained buffers", "invalid capacity");
    std::uint64_t total = 0;
    for (const auto& layout : kStageLayouts) {
        for (const auto bytes : StageBufferBytes(layout, capacity)) {
            auto required = device.RequiredBufferAllocationBytes(bytes, BufferUsage::kStorage);
            if (!required) return std::unexpected(required.error());
            if (*required < bytes)
                return Fail(ErrorDomain::kDevice, "plan retained buffers",
                            "allocation requirement is smaller than the requested buffer");
            if (*required > std::numeric_limits<std::uint64_t>::max() - total)
                return Fail(ErrorDomain::kDevice, "plan retained buffers",
                            "allocation requirement sum overflows");
            total += *required;
        }
    }
    return total;
}

std::array<RetainedCompute::StageStats, RetainedCompute::kStageCount> RetainedCompute::Statistics()
    const {
    return {camera_.stats,     transport_.stats,  endpoint_.stats,   dense_.stats,
            initialize_.stats, ray_camera_.stats, dopri_phase_.stats};
}

void RetainedCompute::ResetStatistics() {
    for (auto* stage :
         {&camera_, &transport_, &endpoint_, &dense_, &initialize_, &ray_camera_, &dopri_phase_}) {
        auto counts = std::move(stage->stats.command_row_counts);
        std::fill(counts.begin(), counts.end(), 0);
        stage->stats = {.command_row_counts = std::move(counts)};
    }
    endpoint_dense_stats_ = {};
    submission_feedback_ = {};
}

base::Expected<void> RetainedCompute::Dispatch(Stage& stage, std::size_t active_rows,
                                               DispatchTiming* timing) {
    using Clock = std::chrono::steady_clock;
    const auto milliseconds = [](Clock::time_point started, Clock::time_point finished) {
        return std::chrono::duration<double, std::milli>(finished - started).count();
    };
    // A failed write must not leave the caller's previous dispatch observation.
    if (timing) *timing = {};
    auto input = std::span(stage.input);
    // The fixed-capacity header and program tail define shader indexing. Once
    // that immutable tail has reached the device, only requested rows change.
    if (stage.program_uploaded) input = input.first(1 + active_rows * stage.input_words_per_row);
    const auto write_started = Clock::now();
    auto status = device_.WriteBuffer(stage.buffers[0], std::as_bytes(input));
    const auto write_finished = Clock::now();
    stage.stats.write_buffer_ms += milliseconds(write_started, write_finished);
    if (!status) return status;
    stage.program_uploaded = true;
    stage.stats.write_buffer_bytes += input.size_bytes();
    DispatchTiming measured;
    auto* observed = timing ? timing : &measured;
    // Buffer strides and the immutable program retain their capacity layout.
    // Only requested rows execute; a later larger batch clears its own rows.
    status = device_.Dispatch(
        stage.kernel, stage.buffers,
        static_cast<std::uint32_t>(stage.paired_rows ? (active_rows + 1) / 2
                                                     : (active_rows + kRetainedGroupRows - 1) /
                                                           kRetainedGroupRows),
        stage.paired_rows && active_rows % 2 != 0 ? 2 : 1, 1, observed);
    if (!status) return status;
    if (!std::isfinite(observed->submit_wait_ms) || observed->submit_wait_ms < 0)
        return Fail(ErrorDomain::kDevice, "dispatch retained stage", "invalid submission timing");
    ++stage.stats.submissions;
    ++stage.stats.command_row_counts[active_rows];
    stage.stats.submit_wait_ms += observed->submit_wait_ms;
    stage.stats.maximum_submit_wait_ms =
        std::max(stage.stats.maximum_submit_wait_ms, observed->submit_wait_ms);
    stage.stats.pipeline_setup_ms += observed->pipeline_setup_ms;
    stage.stats.command_setup_ms += observed->command_setup_ms;
    stage.stats.cleanup_ms += observed->cleanup_ms;
    stage.stats.dispatch_total_ms += observed->total_ms;
    stage.stats.pipeline_creations += observed->pipeline_created ? 1 : 0;
    if (observed->submit_wait_ms >= submission_feedback_.peak_ms) {
        submission_feedback_.peak_ms = observed->submit_wait_ms;
        submission_feedback_.peak_rows = active_rows;
    }
    submission_feedback_.maximum_rows = std::max(submission_feedback_.maximum_rows, active_rows);
    if (active_rows == 1 && observed->submit_wait_ms >= submission_feedback_.maximum_one_row_ms) {
        submission_feedback_.maximum_one_row_ms = observed->submit_wait_ms;
        submission_feedback_.maximum_one_row_stage = stage.kind;
    }
    if (dispatch_target_ms_ > 0 && observed->submit_wait_ms > dispatch_target_ms_)
        ++stage.stats.target_overshoots;
    auto output = std::span(stage.output).first(active_rows * (stage.output.size() / capacity_));
    const auto read_started = Clock::now();
    status = device_.ReadBuffer(stage.buffers[1], std::as_writable_bytes(output));
    const auto read_finished = Clock::now();
    stage.stats.read_buffer_ms += milliseconds(read_started, read_finished);
    if (status) {
        stage.stats.read_buffer_bytes += output.size_bytes();
        ObserveReadback(stage, active_rows);
    }
    return status;
}

void RetainedCompute::ObserveReadback(Stage& stage, std::size_t active_rows) {
    const bool camera = stage.kind == KernelStage::kCamera || stage.kind == KernelStage::kRayCamera;
    const std::uint32_t one = camera ? std::bit_cast<std::uint32_t>(1.0f) : 1U;
    const auto stride = stage.output.size() / capacity_;
    for (std::size_t row = 0; row < active_rows; ++row) {
        const auto flag = stage.output[row * stride + (camera ? 3 : 0)];
        ++stage.stats.completion_flag_counts[flag == 0 ? 0 : (flag == one ? 1 : 2)];
    }
}

RetainedCompute::SubmissionFeedback RetainedCompute::TakeSubmissionFeedback() {
    const auto feedback = submission_feedback_;
    submission_feedback_ = {};
    return feedback;
}

double RetainedCompute::TakeSubmissionPeakMs() { return TakeSubmissionFeedback().peak_ms; }

base::Expected<std::vector<RetainedCameraOutput>> RetainedCompute::Camera(
    std::span<const RetainedCameraInput> inputs, DispatchTiming* timing) {
    if (inputs.empty() || inputs.size() > capacity_)
        return Fail(ErrorDomain::kDevice, "dispatch retained camera", "invalid batch size");
    static_assert(sizeof(RetainedCameraInput) == 160 * sizeof(std::uint32_t));
    std::fill_n(camera_.input.begin() + 1, capacity_ * 160, 0U);
    std::memcpy(camera_.input.data() + 1, inputs.data(), inputs.size_bytes());
    auto status = Dispatch(camera_, inputs.size(), timing);
    if (!status) return std::unexpected(status.error());
    std::vector<RetainedCameraOutput> result(inputs.size());
    for (std::size_t row = 0; row < inputs.size(); ++row) {
        const auto* words = camera_.output.data() + row * kCameraRowWords;
        if (words[3] == 0) continue;
        if (std::bit_cast<float>(words[3]) != 1 ||
            std::memcmp(&inputs[row], words + 8, sizeof(RetainedCameraInput)) != 0)
            return Fail(ErrorDomain::kKernel, "read retained camera", "invalid row ownership");
        for (std::size_t i = 0; i < result[row].values.size(); ++i) {
            auto& value = result[row].values[i];
            value = {std::bit_cast<float>(words[168 + 3 * i]),
                     std::bit_cast<float>(words[169 + 3 * i]), 0,
                     std::bit_cast<float>(words[170 + 3 * i]), 1};
            if (!value.IsRepresented())
                return Fail(ErrorDomain::kKernel, "read retained camera", "invalid retained value");
        }
        result[row].valid = true;
    }
    return result;
}

base::Expected<std::vector<RetainedStepOutput>> RetainedCompute::Step(
    std::span<const RetainedStepInput> inputs, DispatchTiming* timing) {
    if (inputs.empty() || inputs.size() > capacity_)
        return Fail(ErrorDomain::kDevice, "dispatch retained transport", "invalid batch size");
    static_assert(sizeof(RetainedStepInput) == 230 * sizeof(std::uint32_t));
    std::fill_n(transport_.input.begin() + 1, capacity_ * 230, 0U);
    std::memcpy(transport_.input.data() + 1, inputs.data(), inputs.size_bytes());
    auto status = Dispatch(transport_, inputs.size(), timing);
    if (!status) return std::unexpected(status.error());
    std::vector<RetainedStepOutput> result(inputs.size());
    for (std::size_t row = 0; row < inputs.size(); ++row) {
        const auto* words = transport_.output.data() + row * kTransportRowWords;
        auto& output = result[row];
        output.stages = words[1];
        if (output.stages > 7)
            return Fail(ErrorDomain::kKernel, "read retained transport", "invalid stage count");
        if (words[0] == 0) continue;
        if (words[0] != 1 || output.stages != 7 || words[2] > 1)
            return Fail(ErrorDomain::kKernel, "read retained transport",
                        "invalid completion state");
        const std::array groups{&output.fifth, &output.fourth, &output.increment, &output.error};
        for (std::size_t group = 0; group < groups.size(); ++group)
            for (std::size_t i = 0; i < 40; ++i) {
                auto& value = (*groups[group])[i];
                value = Decode(words + 4 + group * 200 + i * 5);
                if (!value.IsRepresented())
                    return Fail(ErrorDomain::kKernel, "read retained transport",
                                "invalid retained value");
            }
        output.rhs_valid = words[2] == 1;
        if (output.rhs_valid)
            for (std::size_t stage = 0; stage < 7; ++stage)
                for (std::size_t i = 0; i < 40; ++i) {
                    auto& value = output.rhs[stage][i];
                    value = Decode(words + 1004 + stage * 200 + i * 5);
                    if (!value.IsRepresented())
                        return Fail(ErrorDomain::kKernel, "read retained transport",
                                    "invalid retained RHS value");
                }
        output.valid = true;
    }
    return result;
}

base::Expected<std::vector<RetainedEndpointOutput>> RetainedCompute::Endpoint(
    std::span<const RetainedEndpointInput> inputs, DispatchTiming* timing) {
    if (inputs.empty() || inputs.size() > capacity_)
        return Fail(ErrorDomain::kDevice, "dispatch retained endpoint", "invalid batch size");
    static_assert(sizeof(RetainedEndpointInput) == 225 * sizeof(std::uint32_t));
    std::fill_n(endpoint_.input.begin() + 1, capacity_ * 225, 0U);
    std::memcpy(endpoint_.input.data() + 1, inputs.data(), inputs.size_bytes());
    auto status = Dispatch(endpoint_, inputs.size(), timing);
    if (!status) return std::unexpected(status.error());
    std::vector<RetainedEndpointOutput> result(inputs.size());
    for (std::size_t row = 0; row < inputs.size(); ++row) {
        auto output = DecodeEndpoint(endpoint_.output.data() + row * kEndpointRowWords);
        if (!output) return std::unexpected(output.error());
        result[row] = *output;
    }
    return result;
}

base::Expected<std::vector<RetainedDenseOutput>> RetainedCompute::Dense(
    std::span<const RetainedDenseInput> inputs, DispatchTiming* timing) {
    if (inputs.empty() || inputs.size() > capacity_)
        return Fail(ErrorDomain::kDevice, "dispatch retained dense segment", "invalid batch size");
    static_assert(sizeof(RetainedDenseInput) == 560 * sizeof(std::uint32_t));
    std::fill_n(dense_.input.begin() + 1, capacity_ * 560, 0U);
    std::memcpy(dense_.input.data() + 1, inputs.data(), inputs.size_bytes());
    auto status = Dispatch(dense_, inputs.size(), timing);
    if (!status) return std::unexpected(status.error());
    std::vector<RetainedDenseOutput> result(inputs.size());
    for (std::size_t row = 0; row < inputs.size(); ++row) {
        auto output = DecodeDense(dense_.output.data() + row * kDenseRowWords);
        if (!output) return std::unexpected(output.error());
        result[row] = *output;
    }
    return result;
}

base::Expected<std::vector<RetainedDopriPhaseOutput>> RetainedCompute::DopriPhase(
    std::span<const RetainedDopriPhaseInput> inputs, DispatchTiming* timing) {
    if (inputs.empty() || inputs.size() > capacity_)
        return Fail(ErrorDomain::kDevice, "dispatch retained DP phase", "invalid batch size");
    static_assert(sizeof(RetainedDopriPhaseInput) == 1810 * sizeof(std::uint32_t));
    std::fill_n(dopri_phase_.input.begin() + 1, capacity_ * 1810, 0U);
    std::memcpy(dopri_phase_.input.data() + 1, inputs.data(), inputs.size_bytes());
    auto status = Dispatch(dopri_phase_, inputs.size(), timing);
    if (!status) return std::unexpected(status.error());
    std::vector<RetainedDopriPhaseOutput> result(inputs.size());
    for (std::size_t row = 0; row < inputs.size(); ++row) {
        const auto* words = dopri_phase_.output.data() + row * kDopriPhaseRowWords;
        auto& output = result[row];
        if (words[0] == 0) continue;
        if (words[0] != 1)
            return Fail(ErrorDomain::kKernel, "read retained DP phase", "invalid completion state");
        const std::array groups{&output.phase, &output.derivative, &output.a, &output.b, &output.c};
        for (std::size_t group = 0; group < groups.size(); ++group)
            for (std::size_t i = 0; i < 40; ++i) {
                auto& value = (*groups[group])[i];
                value = Decode(words + 4 + group * 200 + i * 5);
                if (!value.IsRepresented())
                    return Fail(ErrorDomain::kKernel, "read retained DP phase",
                                "invalid retained value");
            }
        output.valid = true;
    }
    return result;
}

base::Expected<RetainedEndpointDenseOutput> RetainedCompute::EndpointAndDense(
    std::span<const RetainedEndpointInput> endpoints, std::span<const RetainedDenseInput> dense) {
    const auto rows = std::max(endpoints.size(), dense.size());
    if (!SupportsIndependentPair() || endpoints.empty() || dense.empty() || rows <= 1 ||
        rows > capacity_)
        return Fail(ErrorDomain::kDevice, "dispatch retained endpoint and dense",
                    "invalid or unsupported independent pair");
    using Clock = std::chrono::steady_clock;
    const auto milliseconds = [](auto started, auto finished) {
        return std::chrono::duration<double, std::milli>(finished - started).count();
    };
    std::fill_n(endpoint_.input.begin() + 1, capacity_ * 225, 0U);
    std::memcpy(endpoint_.input.data() + 1, endpoints.data(), endpoints.size_bytes());
    std::fill_n(dense_.input.begin() + 1, capacity_ * 560, 0U);
    std::memcpy(dense_.input.data() + 1, dense.data(), dense.size_bytes());
    const auto write = [&](Stage& stage, std::size_t active) -> base::Expected<void> {
        auto input = std::span(stage.input);
        if (stage.program_uploaded) input = input.first(1 + active * stage.input_words_per_row);
        const auto started = Clock::now();
        auto status = device_.WriteBuffer(stage.buffers[0], std::as_bytes(input));
        stage.stats.write_buffer_ms += milliseconds(started, Clock::now());
        if (status) {
            stage.program_uploaded = true;
            stage.stats.write_buffer_bytes += input.size_bytes();
        }
        return status;
    };
    auto status = write(endpoint_, endpoints.size());
    if (!status) return std::unexpected(status.error());
    status = write(dense_, dense.size());
    if (!status) return std::unexpected(status.error());
    const std::array commands{
        ComputeDispatch{endpoint_.kernel, endpoint_.buffers,
                        static_cast<std::uint32_t>((endpoints.size() + kRetainedGroupRows - 1) /
                                                   kRetainedGroupRows),
                        1, 1},
        ComputeDispatch{dense_.kernel, dense_.buffers,
                        static_cast<std::uint32_t>((dense.size() + kRetainedGroupRows - 1) /
                                                   kRetainedGroupRows),
                        1, 1}};
    IndependentPairTiming observed;
    status = device_.DispatchIndependentPair(commands, &observed);
    if (!status) return std::unexpected(status.error());
    const auto& timing = observed.combined;
    if (!std::isfinite(timing.submit_wait_ms) || timing.submit_wait_ms < 0)
        return Fail(ErrorDomain::kDevice, "dispatch retained endpoint and dense",
                    "invalid submission timing");
    ++endpoint_.stats.submissions;
    ++dense_.stats.submissions;
    ++endpoint_.stats.command_row_counts[endpoints.size()];
    ++dense_.stats.command_row_counts[dense.size()];
    auto& stats = endpoint_dense_stats_;
    ++stats.submissions;
    stats.submit_wait_ms += timing.submit_wait_ms;
    stats.maximum_submit_wait_ms = std::max(stats.maximum_submit_wait_ms, timing.submit_wait_ms);
    stats.pipeline_setup_ms += timing.pipeline_setup_ms;
    stats.command_setup_ms += timing.command_setup_ms;
    stats.cleanup_ms += timing.cleanup_ms;
    stats.dispatch_total_ms += timing.total_ms;
    stats.pipeline_creations += observed.pipeline_creations;
    if (dispatch_target_ms_ > 0 && timing.submit_wait_ms > dispatch_target_ms_)
        ++stats.target_overshoots;
    if (timing.submit_wait_ms >= submission_feedback_.peak_ms) {
        submission_feedback_.peak_ms = timing.submit_wait_ms;
        submission_feedback_.peak_rows = rows;
    }
    submission_feedback_.maximum_rows = std::max(submission_feedback_.maximum_rows, rows);
    // A pair is reducible even for one logical ray. The caller must retry it
    // through budget-one individual stages, never call it an irreducible stage.
    const auto read = [&](Stage& stage, std::size_t active) -> base::Expected<void> {
        auto output = std::span(stage.output).first(active * (stage.output.size() / capacity_));
        const auto started = Clock::now();
        auto result = device_.ReadBuffer(stage.buffers[1], std::as_writable_bytes(output));
        stage.stats.read_buffer_ms += milliseconds(started, Clock::now());
        if (result) {
            stage.stats.read_buffer_bytes += output.size_bytes();
            ObserveReadback(stage, active);
        }
        return result;
    };
    status = read(endpoint_, endpoints.size());
    if (!status) return std::unexpected(status.error());
    RetainedEndpointDenseOutput result;
    result.endpoints.reserve(endpoints.size());
    for (std::size_t row = 0; row < endpoints.size(); ++row) {
        auto output = DecodeEndpoint(endpoint_.output.data() + row * kEndpointRowWords);
        if (!output) return std::unexpected(output.error());
        result.endpoints.push_back(*output);
    }
    status = read(dense_, dense.size());
    if (!status) return std::unexpected(status.error());
    result.dense.reserve(dense.size());
    for (std::size_t row = 0; row < dense.size(); ++row)
        result.dense.push_back(DecodeDense(dense_.output.data() + row * kDenseRowWords));
    return result;
}

base::Expected<std::vector<RetainedInitializeOutput>> RetainedCompute::Initialize(
    std::span<const RetainedInitializeInput> inputs, DispatchTiming* timing) {
    if (inputs.empty() || inputs.size() > capacity_)
        return Fail(ErrorDomain::kDevice, "dispatch retained initialization", "invalid batch size");
    static_assert(sizeof(RetainedInitializeInput) == 225 * sizeof(std::uint32_t));
    std::fill_n(initialize_.input.begin() + 1, capacity_ * 225, 0U);
    std::memcpy(initialize_.input.data() + 1, inputs.data(), inputs.size_bytes());
    auto status = Dispatch(initialize_, inputs.size(), timing);
    if (!status) return std::unexpected(status.error());
    std::vector<RetainedInitializeOutput> result(inputs.size());
    for (std::size_t row = 0; row < inputs.size(); ++row) {
        const auto* words = initialize_.output.data() + row * kInitializeRowWords;
        if (words[0] == 0) continue;
        if (words[0] != 1)
            return Fail(ErrorDomain::kKernel, "read retained initialization",
                        "invalid completion state");
        for (std::size_t i = 0; i < 40; ++i) {
            result[row].phase[i] = Decode(words + 4 + i * 5);
            if (!result[row].phase[i].IsRepresented())
                return Fail(ErrorDomain::kKernel, "read retained initialization",
                            "invalid retained value");
        }
        result[row].valid = true;
    }
    return result;
}

base::Expected<std::vector<RetainedCameraOutput>> RetainedCompute::RayCamera(
    std::span<const RetainedRayCameraInput> inputs, DispatchTiming* timing) {
    if (inputs.empty() || inputs.size() > capacity_)
        return Fail(ErrorDomain::kDevice, "dispatch retained ray camera", "invalid batch size");
    static_assert(sizeof(RetainedRayCameraInput) == 225 * sizeof(std::uint32_t));
    std::fill_n(ray_camera_.input.begin() + 1, capacity_ * 225, 0U);
    std::memcpy(ray_camera_.input.data() + 1, inputs.data(), inputs.size_bytes());
    auto status = Dispatch(ray_camera_, inputs.size(), timing);
    if (!status) return std::unexpected(status.error());
    std::vector<RetainedCameraOutput> result(inputs.size());
    for (std::size_t row = 0; row < inputs.size(); ++row) {
        const auto* words = ray_camera_.output.data() + row * kRayCameraRowWords;
        if (words[3] == 0) continue;
        if (std::bit_cast<float>(words[3]) != 1 ||
            std::memcmp(&inputs[row], words + 8, sizeof(RetainedRayCameraInput)) != 0)
            return Fail(ErrorDomain::kKernel, "read retained ray camera", "invalid row ownership");
        for (std::size_t i = 0; i < 104; ++i) {
            auto& value = result[row].values[i];
            value = {std::bit_cast<float>(words[233 + 3 * i]),
                     std::bit_cast<float>(words[234 + 3 * i]), 0,
                     std::bit_cast<float>(words[235 + 3 * i]), 1};
            if (!value.IsRepresented())
                return Fail(ErrorDomain::kKernel, "read retained ray camera",
                            "invalid retained value");
        }
        result[row].valid = true;
    }
    return result;
}

}  // namespace sirius::backend
