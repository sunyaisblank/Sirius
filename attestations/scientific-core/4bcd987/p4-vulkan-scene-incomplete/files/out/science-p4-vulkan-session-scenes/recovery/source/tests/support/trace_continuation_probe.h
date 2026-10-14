#pragma once

// Direct test producer for the real trace ABI. Physics and record validation
// belong to production; this helper owns only bounded phases and private rows.
#include "sirius/backend/device.h"
#include "sirius/render/dispatch_governor.h"
#include "sirius/render/trace_continuation.h"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <string>
#include <vector>

namespace sirius::test {

template <typename Real>
class TraceContinuationProbe {
  public:
    using Record = render::TraceContinuationRecord<Real>;
    static constexpr std::size_t kMaximumActive = 64;

    TraceContinuationProbe(backend::ComputeDevice& device, backend::KernelHandle kernel)
        : device_(device), kernel_(kernel) {}

    base::Expected<void> Prepare() {
        // The empty CSR covers all cells even in a point-bundle diagnostic.
        const std::array<std::size_t, 7> sizes{kMaximumActive * 4 * sizeof(float),
                                               render::kTraceParameterCount * sizeof(float),
                                               sizeof(std::uint32_t),
                                               8 * sizeof(float),
                                               (256 * 512 + 1) * sizeof(std::uint32_t),
                                               sizeof(std::uint32_t),
                                               kMaximumActive * sizeof(Record)};
        for (const auto size : sizes) {
            auto buffer = device_.CreateBuffer(size, backend::BufferUsage::kStorage);
            if (!buffer) return std::unexpected(buffer.error());
            bindings_.push_back(*buffer);
            const std::vector<std::byte> zero(size);
            auto uploaded = device_.WriteBuffer(*buffer, zero);
            if (!uploaded) return uploaded;
        }
        records_.resize(kMaximumActive);
        if (device_.Info().kind == backend::DeviceKind::kSoftware) {
            // Production's separate zero-active driver initialization. It owns
            // no rays; every real Init/Advance/Finalize below has the hard bound.
            std::array<float, render::kTraceParameterCount> empty{};
            auto uploaded = device_.WriteBuffer(bindings_[1], std::as_bytes(std::span(empty)));
            if (!uploaded) return uploaded;
            backend::DispatchTiming timing;
            auto warmed = device_.Dispatch(kernel_, bindings_, 1, 1, 1, &timing);
            if (!warmed) return warmed;
            if (!std::isfinite(timing.submit_wait_ms) || timing.submit_wait_ms <= 0.0)
                return Failure("invalid zero-active initialization timing");
        }
        return {};
    }

    base::Expected<std::vector<float>> RunRegion(std::span<const float> parameters) {
        if (parameters.size() != render::kTraceParameterCount || bindings_.size() != 7)
            return Failure("missing production parameters/bindings");
        const float sample = parameters[46];
        if (!std::isfinite(sample) || sample < 0 || sample != std::floor(sample) ||
            sample > static_cast<float>(std::numeric_limits<int>::max() / 2))
            return Failure("invalid original sample index");
        if (sample == 0.0f) {
            failed_ = false;
            next_sample_ = 0;
        } else if (failed_ || static_cast<std::uint32_t>(sample) != next_sample_) {
            return Failure("failed or incomplete sample batch requires sample-zero restart");
        }
        failed_ = true;  // Cleared only by completed, validated Finalize/readback.
        for (const auto index : {31, 32, 33, 34}) {
            if (!std::isfinite(parameters[index]) ||
                parameters[index] != std::floor(parameters[index]) ||
                parameters[index] < (index < 33 ? 0.0f : 1.0f) || parameters[index] > 100000.0f)
                return Failure("invalid active rectangle");
        }
        const auto width = static_cast<std::uint32_t>(parameters[33]);
        const auto height = static_cast<std::uint32_t>(parameters[34]);
        if (static_cast<std::uint64_t>(width) * height > kMaximumActive)
            return Failure("direct probe exceeds its active-pixel cap");
        active_ = width * height;
        if (!std::isfinite(parameters[21]) || parameters[21] < 1 || parameters[21] > 1000000.0f ||
            parameters[21] != std::floor(parameters[21]))
            return Failure("invalid attempt bound");
        maximum_attempts_ = static_cast<std::uint32_t>(parameters[21]);
        std::copy(parameters.begin(), parameters.end(), packet_.begin());
        if (++epoch_ == 0) return Failure("epoch exhausted");
        packet_[render::kTraceEpochParameter] = std::bit_cast<float>(epoch_);
        auto initial = DispatchPhase(render::TraceAction::Initialise);
        if (!initial) return std::unexpected(initial.error());
        bool terminal = *initial;
        for (std::uint32_t step = 0; !terminal && step < maximum_attempts_; ++step) {
            auto advanced = DispatchPhase(render::TraceAction::Advance);
            if (!advanced) return std::unexpected(advanced.error());
            terminal = *advanced;
        }
        if (!terminal) return Failure("unfinished trace after bounded attempts");
        auto finalized = DispatchPhase(render::TraceAction::Finalize);
        if (!finalized) return std::unexpected(finalized.error());
        auto output = ReadRadiance();
        if (!output) return output;
        for (const auto value : *output)
            if (!std::isfinite(value)) return Failure("non-finite finalized radiance");
        failed_ = false;
        next_sample_ = static_cast<std::uint32_t>(sample) + 1;
        return output;
    }

    // A negative control sends Finalize to the real failed shader state. It
    // deliberately skips success validation, then exposes actual bytes so the
    // test can require the same failed status and no radiance publication.
    base::Expected<void> FinalizeFailedStateForControl() {
        if (!failed_ || active_ == 0) return Failure("no failed state to inspect");
        auto result = DispatchPhase(render::TraceAction::Finalize, false);
        if (!result) return std::unexpected(result.error());
        return {};
    }

    base::Expected<std::vector<float>> ReadRadiance() {
        std::vector<float> values(active_ * 4);
        auto read = device_.ReadBuffer(bindings_[0], std::as_writable_bytes(std::span(values)));
        if (!read) return std::unexpected(read.error());
        return values;
    }

    [[nodiscard]] std::span<const Record> Records() const {
        return std::span(records_).first(active_);
    }
    [[nodiscard]] std::uint64_t PhysicalSubmissions() const { return physical_submissions_; }

  private:
    static std::unexpected<base::Error> Failure(const std::string& detail) {
        return base::Fail(base::ErrorDomain::kKernel, "direct trace continuation probe", detail);
    }

    base::Expected<bool> DispatchPhase(render::TraceAction action, bool validate = true) {
        packet_[render::kTraceActionParameter] = static_cast<float>(action);
        auto uploaded = device_.WriteBuffer(bindings_[1], std::as_bytes(std::span(packet_)));
        if (!uploaded) return std::unexpected(uploaded.error());
        const auto previous = records_;
        const auto width = static_cast<std::uint32_t>(packet_[33]);
        const auto height = static_cast<std::uint32_t>(packet_[34]);
        backend::DispatchTiming timing;
        ++physical_submissions_;
        auto dispatched =
            device_.Dispatch(kernel_, bindings_, (width + 7) / 8, (height + 7) / 8, 1, &timing);
        if (!dispatched) return std::unexpected(dispatched.error());
        const bool admitted = bands_.Record(active_, timing.submit_wait_ms);
        if (!admitted || !std::isfinite(timing.submit_wait_ms) || timing.submit_wait_ms <= 0.0 ||
            timing.submit_wait_ms > render::kDispatchStopMs)
            return Failure("physical phase exceeded the valid 1000 ms timing envelope");
        auto read = device_.ReadBuffer(bindings_[6], std::as_writable_bytes(std::span(records_)));
        if (!read) return std::unexpected(read.error());
        bool terminal = true;
        for (std::size_t index = 0; index < active_ && validate; ++index) {
            const std::array<std::uint32_t, 4> identity{
                static_cast<std::uint32_t>(packet_[31]) + static_cast<std::uint32_t>(index) % width,
                static_cast<std::uint32_t>(packet_[32]) + static_cast<std::uint32_t>(index) / width,
                static_cast<std::uint32_t>(packet_[46]), epoch_};
            auto checked = render::ValidateTraceContinuation(
                records_[index],
                action == render::TraceAction::Initialise ? nullptr : &previous[index], action,
                identity, maximum_attempts_);
            if (!checked) return std::unexpected(checked.error());
            terminal = terminal && *checked;
        }
        return terminal;
    }

    backend::ComputeDevice& device_;
    backend::KernelHandle kernel_;
    std::vector<backend::BufferHandle> bindings_;
    std::vector<Record> records_;
    std::array<float, render::kTraceParameterCount> packet_{};
    std::uint32_t epoch_ = 0, next_sample_ = 0, maximum_attempts_ = 0, active_ = 0;
    std::uint64_t physical_submissions_ = 0;
    render::BandController bands_{64, 0.0, 1, 64};
    bool failed_ = true;
};

template <typename Real>
base::Expected<std::vector<float>> TraceContinuationImage(backend::ComputeDevice& device,
                                                          backend::KernelHandle kernel,
                                                          std::span<const float> original,
                                                          std::uint32_t width,
                                                          std::uint32_t height) {
    TraceContinuationProbe<Real> probe(device, kernel);
    auto prepared = probe.Prepare();
    if (!prepared) return std::unexpected(prepared.error());
    if (original.size() != render::kTraceParameterCount)
        return base::Fail(base::ErrorDomain::kKernel, "direct trace image",
                          "wrong parameter count");
    std::vector<float> result(static_cast<std::size_t>(width) * height * 4);
    std::vector<float> packet(original.begin(), original.end());
    for (std::uint32_t y = 0; y < height; ++y) {
        for (std::uint32_t x = 0; x < width; x += 64) {
            const auto count = std::min(64U, width - x);
            packet[31] = original[31] + static_cast<float>(x);
            packet[32] = original[32] + static_cast<float>(y);
            packet[33] = static_cast<float>(count);
            packet[34] = 1.0f;
            auto values = probe.RunRegion(packet);
            if (!values) return std::unexpected(values.error());
            std::copy(values->begin(), values->end(),
                      result.begin() + static_cast<std::ptrdiff_t>(
                                           (static_cast<std::size_t>(y) * width + x) * 4));
        }
    }
    return result;
}

}  // namespace sirius::test
