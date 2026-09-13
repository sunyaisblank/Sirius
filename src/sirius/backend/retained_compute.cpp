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
constexpr std::size_t kMaximumCapacity = 65536;

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

base::Expected<std::unique_ptr<RetainedCompute>> RetainedCompute::Create(ComputeDevice& device,
                                                                         std::size_t capacity) {
    if (capacity == 0 || capacity > kMaximumCapacity)
        return Fail(ErrorDomain::kDevice, "create retained compute stages", "invalid capacity");
    if (!device.Info().preserves_fp32_denormals || !device.Info().rounds_fp32_to_nearest)
        return Fail(ErrorDomain::kDevice, "create retained compute stages",
                    "device lacks binary32 subnormal preservation or round-to-nearest control");
    auto result = std::unique_ptr<RetainedCompute>(new RetainedCompute(device, capacity));
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
    auto status = create(result->camera_, kCameraShader, kCameraProgram, 160, kCameraRowWords);
    if (!status) return std::unexpected(status.error());
    status =
        create(result->transport_, kTransportShader, kTransportProgram, 230, kTransportRowWords);
    if (!status) return std::unexpected(status.error());
    return result;
}

base::Expected<void> RetainedCompute::Dispatch(Stage& stage, DispatchTiming* timing) {
    auto status = device_.WriteBuffer(stage.buffers[0], std::as_bytes(std::span(stage.input)));
    if (!status) return status;
    status = device_.Dispatch(stage.kernel, stage.buffers,
                              static_cast<std::uint32_t>((capacity_ + 63) / 64), 1, 1, timing);
    if (!status) return status;
    return device_.ReadBuffer(stage.buffers[1], std::as_writable_bytes(std::span(stage.output)));
}

base::Expected<std::vector<RetainedCameraOutput>> RetainedCompute::Camera(
    std::span<const RetainedCameraInput> inputs, DispatchTiming* timing) {
    if (inputs.empty() || inputs.size() > capacity_)
        return Fail(ErrorDomain::kDevice, "dispatch retained camera", "invalid batch size");
    static_assert(sizeof(RetainedCameraInput) == 160 * sizeof(std::uint32_t));
    std::fill_n(camera_.input.begin() + 1, capacity_ * 160, 0U);
    std::memcpy(camera_.input.data() + 1, inputs.data(), inputs.size_bytes());
    auto status = Dispatch(camera_, timing);
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
    auto status = Dispatch(transport_, timing);
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

}  // namespace sirius::backend
