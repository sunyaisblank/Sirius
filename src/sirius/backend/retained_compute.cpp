#include "sirius/backend/retained_compute.h"

#include "retained_kernels.h"

#include <algorithm>
#include <bit>
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

RetainedValue Decode(const std::uint32_t* words) {
    return std::bit_cast<RetainedValue>(
        std::array<std::uint32_t, 5>{words[0], words[1], words[2], words[3], words[4]});
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
    if (!device.Info().preserves_fp32_denormals || !device.Info().rounds_fp32_to_nearest)
        return Fail(ErrorDomain::kDevice, "create retained compute stages",
                    "device lacks binary32 subnormal preservation or round-to-nearest control");
    if (fp64_products && (!device.Info().supports_fp64 || !device.Info().rounds_fp64_to_nearest))
        return Fail(ErrorDomain::kDevice, "create retained compute stages",
                    "device lacks binary64 products or round-to-nearest control");
    auto result = std::unique_ptr<RetainedCompute>(new RetainedCompute(device, capacity));
    result->dispatch_target_ms_ = dispatch_target_ms;
    const auto create = [&](Stage& stage, std::span<const std::uint32_t> code,
                            std::span<const std::uint32_t> program, std::size_t input_words,
                            std::size_t row_words) -> base::Expected<void> {
        auto kernel = device.LoadKernel(code);
        if (!kernel) return std::unexpected(kernel.error());
        stage.kernel = *kernel;
        stage.input.resize(1 + input_words * capacity + program.size());
        stage.output.resize(row_words * capacity);
        stage.input[0] = static_cast<std::uint32_t>(capacity);
        std::copy(program.begin(), program.end(), stage.input.begin() + 1 + input_words * capacity);
        auto input = device.CreateBuffer(stage.input.size() * 4, BufferUsage::kStorage);
        if (!input) return std::unexpected(input.error());
        auto output = device.CreateBuffer(stage.output.size() * 4, BufferUsage::kStorage);
        if (!output) return std::unexpected(output.error());
        stage.buffers = {*input, *output};
        return {};
    };
    const auto shader = [fp64_products](std::span<const std::uint32_t> narrow,
                                        std::span<const std::uint32_t> wide) {
        return fp64_products ? wide : narrow;
    };
    auto status = create(result->camera_, shader(kCameraShader, kCameraFp64Shader), kCameraProgram,
                         160, kCameraRowWords);
    if (!status) return std::unexpected(status.error());
    status = create(result->transport_, shader(kTransportShader, kTransportFp64Shader),
                    kTransportProgram, 230, kTransportRowWords);
    if (!status) return std::unexpected(status.error());
    status = create(result->endpoint_, shader(kEndpointShader, kEndpointFp64Shader),
                    kEndpointProgram, 225, kEndpointRowWords);
    if (!status) return std::unexpected(status.error());
    status = create(result->dense_, shader(kDenseShader, kDenseFp64Shader), kDenseProgram, 560,
                    kDenseRowWords);
    if (!status) return std::unexpected(status.error());
    status = create(result->initialize_, shader(kInitializeShader, kInitializeFp64Shader),
                    kInitializeProgram, 225, kInitializeRowWords);
    if (!status) return std::unexpected(status.error());
    status = create(result->ray_camera_, shader(kRayCameraShader, kRayCameraFp64Shader),
                    kRayCameraProgram, 225, kRayCameraRowWords);
    if (!status) return std::unexpected(status.error());
    return result;
}

std::uint64_t RetainedCompute::RequiredBufferBytes(std::size_t capacity) {
    if (capacity == 0 || capacity > kMaximumCapacity) return 0;
    return 4 * (6 + kCameraProgram.size() + kTransportProgram.size() + kEndpointProgram.size() +
                kDenseProgram.size() + kInitializeProgram.size() + kRayCameraProgram.size() +
                capacity *
                    (160 + kCameraRowWords + 230 + kTransportRowWords + 225 + kEndpointRowWords +
                     560 + kDenseRowWords + 225 + kInitializeRowWords + 225 + kRayCameraRowWords));
}

std::array<RetainedCompute::StageStats, 6> RetainedCompute::Statistics() const {
    return {camera_.stats, transport_.stats,  endpoint_.stats,
            dense_.stats,  initialize_.stats, ray_camera_.stats};
}

base::Expected<void> RetainedCompute::Dispatch(Stage& stage, std::size_t active_rows,
                                               DispatchTiming* timing) {
    auto status = device_.WriteBuffer(stage.buffers[0], std::as_bytes(std::span(stage.input)));
    if (!status) return status;
    DispatchTiming measured;
    auto* observed = timing ? timing : &measured;
    // Buffer strides and the immutable program retain their capacity layout.
    // Only requested rows execute; a later larger batch clears its own rows.
    status = device_.Dispatch(
        stage.kernel, stage.buffers,
        static_cast<std::uint32_t>((active_rows + kRetainedGroupRows - 1) / kRetainedGroupRows), 1,
        1, observed);
    if (!status) return status;
    if (!std::isfinite(observed->submit_wait_ms) || observed->submit_wait_ms < 0)
        return Fail(ErrorDomain::kDevice, "dispatch retained stage", "invalid submission timing");
    ++stage.stats.submissions;
    stage.stats.submit_wait_ms += observed->submit_wait_ms;
    stage.stats.maximum_submit_wait_ms =
        std::max(stage.stats.maximum_submit_wait_ms, observed->submit_wait_ms);
    stage.stats.pipeline_setup_ms += observed->pipeline_setup_ms;
    submission_peak_ms_ = std::max(submission_peak_ms_, observed->submit_wait_ms);
    if (dispatch_target_ms_ > 0 && observed->submit_wait_ms > dispatch_target_ms_)
        ++stage.stats.target_overshoots;
    return device_.ReadBuffer(stage.buffers[1], std::as_writable_bytes(std::span(stage.output)));
}

double RetainedCompute::TakeSubmissionPeakMs() {
    const double peak = submission_peak_ms_;
    submission_peak_ms_ = 0;
    return peak;
}

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
        if (words[0] != 1 || output.stages != 7)
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
        const auto* words = endpoint_.output.data() + row * kEndpointRowWords;
        if (words[0] == 0) continue;
        if (words[0] != 1 || words[1] > 3)
            return Fail(ErrorDomain::kKernel, "read retained endpoint", "invalid projection state");
        auto& output = result[row];
        for (std::size_t i = 0; i < 40; ++i) {
            output.phase[i] = Decode(words + 104 + i * 5);
            output.physical[i] = Decode(words + 304 + i * 5);
            if (!output.phase[i].IsRepresented() || !output.physical[i].IsRepresented())
                return Fail(ErrorDomain::kKernel, "read retained endpoint",
                            "invalid retained value");
        }
        output.component = words[1];
        output.valid = true;
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
        const auto* words = dense_.output.data() + row * kDenseRowWords;
        if (words[0] == 0) continue;
        if (words[0] != 1)
            return Fail(ErrorDomain::kKernel, "read retained dense segment",
                        "invalid completion state");
        for (std::size_t i = 0; i < 40; ++i) {
            result[row].physical[i] = Decode(words + 4 + i * 5);
            if (!result[row].physical[i].IsRepresented())
                return Fail(ErrorDomain::kKernel, "read retained dense segment",
                            "invalid retained value");
        }
        result[row].valid = true;
    }
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
