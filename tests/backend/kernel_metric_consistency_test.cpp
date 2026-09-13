#include "sirius/backend/device.h"
#include "sirius/render/dispatch_governor.h"

#include <gtest/gtest.h>

#include "../support/metric_consistency_reference.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <span>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

namespace {

std::string MetricNumber(double value) {
    std::ostringstream stream;
    stream << std::setprecision(std::numeric_limits<double>::max_digits10) << value;
    return stream.str();
}

std::vector<std::uint32_t> ReadMetricProbe(const std::string& path) {
    std::ifstream input(path, std::ios::binary | std::ios::ate);
    if (!input) return {};
    const auto size = input.tellg();
    if (size <= 0 || static_cast<std::size_t>(size) % 4 != 0) return {};
    std::vector<std::uint32_t> words(static_cast<std::size_t>(size) / 4);
    input.seekg(0);
    input.read(reinterpret_cast<char*>(words.data()), size);
    return input ? words : std::vector<std::uint32_t>{};
}

template <class Real>
void CheckMetricConsistency(const std::string& artifact) {
    using namespace sirius::backend;
#ifdef SIRIUS_KERNEL_DIR
    const std::string directory = SIRIUS_KERNEL_DIR;
#else
    const std::string directory;
    GTEST_SKIP() << "Vulkan kernels are not configured";
#endif
    constexpr bool wide = std::is_same_v<Real, double>;
    // Independent smooth-point reference budgets, declared before device runs.
    // They do not replace trajectory, event or image acceptance tolerances.
    constexpr double metric_budget = wide ? 2e-13 : 4e-6;
    constexpr double first_budget = wide ? 2e-12 : 3e-5;
    constexpr double second_budget = wide ? 2e-11 : 2e-4;
    const auto inventory = EnumerateVulkanDevices();
    ASSERT_TRUE(inventory.has_value()) << inventory.error().Description();
    if (inventory->empty()) GTEST_SKIP() << "no Vulkan device available";
    const auto index = ResolveVulkanDeviceIndex(*inventory);
    ASSERT_TRUE(index.has_value()) << index.error().Description();
    auto opened = CreateVulkanDevice(*index);
    ASSERT_TRUE(opened.has_value()) << opened.error().Description();
    auto device = std::move(*opened);
    const auto& info = device->Info();
    ::testing::Test::RecordProperty("metric_device_index", std::to_string(*index));
    ::testing::Test::RecordProperty("metric_device_name", info.name);
    ::testing::Test::RecordProperty("metric_device_kind", ToString(info.kind));
    ::testing::Test::RecordProperty("metric_vendor_id", std::to_string(info.vendor_id));
    ::testing::Test::RecordProperty("metric_device_id", std::to_string(info.device_id));
    ::testing::Test::RecordProperty("metric_driver_id", std::to_string(info.driver_id));
    ::testing::Test::RecordProperty("metric_driver_name", info.driver_name);
    ::testing::Test::RecordProperty("metric_driver_info", info.driver_info);
    ::testing::Test::RecordProperty("metric_supports_fp64", info.supports_fp64 ? "true" : "false");
    const auto words = ReadMetricProbe(directory + "/" + artifact);
    ASSERT_FALSE(words.empty()) << artifact;
    const auto kernel = device->LoadKernel(words);
    if constexpr (wide) {
        if (!info.supports_fp64) {
            ASSERT_FALSE(kernel.has_value());
            EXPECT_EQ(kernel.error().domain(), sirius::base::ErrorDomain::kKernel);
            EXPECT_NE(kernel.error().detail().find("shaderFloat64"), std::string::npos);
            ::testing::Test::RecordProperty("metric_evidence", "unsupported_fp64_declined");
            return;
        }
    }
    ASSERT_TRUE(kernel.has_value()) << kernel.error().Description();
    constexpr std::uint64_t allocation_limit = 1024 * 1024;
    constexpr std::uint64_t requested_bytes = 428 * sizeof(Real);
    ASSERT_TRUE(device->SetBufferAllocationLimit(allocation_limit).has_value());
    const auto input_buffer = device->CreateBuffer(8 * sizeof(Real), BufferUsage::kStorage);
    const auto output_buffer = device->CreateBuffer(420 * sizeof(Real), BufferUsage::kStorage);
    ASSERT_TRUE(input_buffer.has_value());
    ASSERT_TRUE(output_buffer.has_value());
    const std::array<BufferHandle, 2> bindings{*input_buffer, *output_buffer};
    ::testing::Test::RecordProperty("metric_allocated_bytes",
                                    std::to_string(device->BufferAllocationBytes()));
    ::testing::Test::RecordProperty("metric_requested_bytes", std::to_string(requested_bytes));
    ::testing::Test::RecordProperty("metric_allocation_limit_bytes",
                                    std::to_string(allocation_limit));
    // Native adapters may require padding beyond the buffer payload.
    ASSERT_GE(device->BufferAllocationBytes(), requested_bytes);
    ASSERT_LE(device->BufferAllocationBytes(), allocation_limit);
    std::size_t ordinal = 0;
    for (const auto& reference : metric_consistency_reference::Cases()) {
        SCOPED_TRACE(reference.name);
        std::array<Real, 8> input{};
        std::transform(reference.input.begin(), reference.input.end(), input.begin(),
                       [](double x) { return Real(x); });
        std::array<Real, 420> output;
        output.fill(std::numeric_limits<Real>::quiet_NaN());
        ASSERT_TRUE(
            device->WriteBuffer(*input_buffer, std::as_bytes(std::span(input))).has_value());
        ASSERT_TRUE(
            device->WriteBuffer(*output_buffer, std::as_bytes(std::span(output))).has_value());
        DispatchTiming timing;
        const auto dispatched = device->Dispatch(*kernel, bindings, 1, 1, 1, &timing);
        ASSERT_TRUE(dispatched.has_value()) << dispatched.error().Description();
        const std::string key = "metric_" + std::to_string(ordinal++);
        ::testing::Test::RecordProperty(key + "_submit_wait_ms",
                                        MetricNumber(timing.submit_wait_ms));
        ::testing::Test::RecordProperty(key + "_pipeline_setup_ms",
                                        MetricNumber(timing.pipeline_setup_ms));
        ::testing::Test::RecordProperty(key + "_total_ms", MetricNumber(timing.total_ms));
        ASSERT_TRUE(std::isfinite(timing.submit_wait_ms));
        ASSERT_GT(timing.submit_wait_ms, 0.0);
        ASSERT_LE(timing.submit_wait_ms, sirius::render::kDispatchStopMs);
        ASSERT_TRUE(device->ReadBuffer(*output_buffer, std::as_writable_bytes(std::span(output)))
                        .has_value());
        for (std::size_t j = 0; j < output.size(); ++j)
            ASSERT_TRUE(std::isfinite(output[j])) << "unwritten/nonfinite scalar " << j;
        ASSERT_EQ(output[0], Real(1));
        std::array<double, 3> maximum_error{};
        for (std::size_t j = 0; j < reference.values.size(); ++j) {
            const double budget = j < 32 ? metric_budget : (j < 160 ? first_budget : second_budget);
            const double expected = reference.values[j];
            const double normalized =
                std::abs(double(output[j + 1]) - expected) / (1.0 + std::abs(expected));
            const std::size_t category = j < 32 ? 0 : (j < 160 ? 1 : 2);
            maximum_error[category] = std::max(maximum_error[category], normalized);
            EXPECT_LE(normalized, budget) << "reference scalar " << j;
        }
        // Algebraic identities use the actual returned operands; a host
        // long-double accumulator limits extra rounding in these checks.
        for (int mu = 0; mu < 4; ++mu)
            for (int nu = 0; nu < 4; ++nu) {
                long double product = 0, scale = 1;
                for (int rho = 0; rho < 4; ++rho) {
                    const long double term = static_cast<long double>(output[1 + 4 * mu + rho]) *
                                             output[17 + 4 * rho + nu];
                    product += term;
                    scale += std::abs(term);
                }
                EXPECT_LE(std::abs(product - (mu == nu ? 1 : 0)),
                          64 * std::numeric_limits<Real>::epsilon() * scale)
                    << "inverse " << mu << "," << nu;
                for (int column = 0; column < 4; ++column) {
                    long double compatible = 0, derivative_scale = 1;
                    for (int rho = 0; rho < 4; ++rho) {
                        const long double first =
                            static_cast<long double>(output[33 + 16 * rho + 4 * column + mu]) *
                            output[1 + 4 * rho + nu];
                        const long double second =
                            static_cast<long double>(output[33 + 16 * rho + 4 * column + nu]) *
                            output[1 + 4 * mu + rho];
                        compatible += first + second;
                        derivative_scale += std::abs(first) + std::abs(second);
                    }
                    const long double derivative = output[97 + 16 * column + 4 * mu + nu];
                    EXPECT_LE(std::abs(derivative - compatible),
                              first_budget * (derivative_scale + std::abs(derivative)))
                        << "metric compatibility " << column << "," << mu << "," << nu;
                }
            }
        EXPECT_LE(output[417], metric_budget);
        EXPECT_LE(output[418], first_budget);
        for (std::size_t j = 0; j < maximum_error.size(); ++j)
            ::testing::Test::RecordProperty(key + "_maximum_normalized_error_" + std::to_string(j),
                                            MetricNumber(maximum_error[j]));
        std::cout << "[MetricConsistency] " << reference.name
                  << " submit_wait_ms=" << timing.submit_wait_ms
                  << " pipeline_setup_ms=" << timing.pipeline_setup_ms
                  << " max_metric=" << maximum_error[0] << " max_first=" << maximum_error[1]
                  << " max_second=" << maximum_error[2] << std::endl;
    }
    ::testing::Test::RecordProperty("metric_actual_submissions", std::to_string(ordinal));
    ::testing::Test::RecordProperty("metric_evidence",
                                    "actual_metric_inverse_connection_and_four_partials");
}

TEST(KernelParity, MetricConnectionAndJetShareRepresentedParametersFp32) {
    ASSERT_NO_FATAL_FAILURE(CheckMetricConsistency<float>("metric_consistency_probe.spv"));
}

TEST(KernelParity, MetricConnectionAndJetShareRepresentedParametersCompensatedFp32) {
    ASSERT_NO_FATAL_FAILURE(CheckMetricConsistency<float>("metric_consistency_probe_fp32comp.spv"));
}

TEST(KernelParity, MetricConnectionAndJetShareRepresentedParametersFp64) {
    ASSERT_NO_FATAL_FAILURE(CheckMetricConsistency<double>("metric_consistency_probe_fp64.spv"));
}

}  // namespace
