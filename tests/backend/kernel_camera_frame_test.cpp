#include "sirius/backend/device.h"
#include "sirius/core/observer_frame.h"
#include "sirius/render/dispatch_governor.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <bit>
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
#include <utility>
#include <vector>

namespace {

std::string CameraNumber(double value) {
    std::ostringstream stream;
    stream << std::setprecision(std::numeric_limits<double>::max_digits10) << value;
    return stream.str();
}

std::vector<std::uint32_t> ReadCameraProbe(const std::string& path) {
    std::ifstream input(path, std::ios::binary | std::ios::ate);
    if (!input) return {};
    const auto size = input.tellg();
    if (size <= 0 || static_cast<std::size_t>(size) % 4 != 0) return {};
    std::vector<std::uint32_t> words(static_cast<std::size_t>(size) / 4);
    input.seekg(0);
    input.read(reinterpret_cast<char*>(words.data()), size);
    return input ? words : std::vector<std::uint32_t>{};
}

// Defining Kerr-Schild metric and inverse, evaluated independently on the host.
// Long-double accumulation reduces contraction error; correctness does not
// depend on long double being wider than double on every supported platform.
using Matrix = std::array<std::array<long double, 4>, 4>;
std::pair<Matrix, Matrix> ReferenceMetric(const std::array<double, 24>& input) {
    Matrix g{}, inverse{};
    for (int i = 0; i < 4; ++i) g[i][i] = inverse[i][i] = i == 0 ? -1 : 1;
    const long double mass = input[0], spin = input[1];
    if (mass == 0 && spin == 0) return {g, inverse};
    const long double x = input[5], y = input[6], z = input[7];
    const long double reduced = x * x + y * y + z * z - spin * spin;
    const long double r2 = (reduced + std::sqrt(reduced * reduced + 4 * spin * spin * z * z)) / 2;
    const long double r = std::sqrt(r2);
    const long double h = 2 * mass * r / (r2 + spin * spin * z * z / r2);
    const std::array<long double, 4> ell{1, (r * x + spin * y) / (r2 + spin * spin),
                                         (r * y - spin * x) / (r2 + spin * spin), z / r};
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) {
            g[i][j] += h * ell[i] * ell[j];
            inverse[i][j] -= h * ell[i] * ell[j] * (i == 0 ? -1 : 1) * (j == 0 ? -1 : 1);
        }
    return {g, inverse};
}

std::vector<std::array<double, 24>> CameraCases() {
    // Non-axis directions exercise all least-aligned choices, then the x-first
    // tie. Float ABI values are promoted exactly before the host calculation.
    constexpr std::array<std::array<float, 3>, 7> directions{
        {{0.1f, 0.7f, 0.7f},
         {0.8f, 0.1f, 0.59f},
         {0.6f, 0.7f, 0.05f},
         {1.0f, 1.0f, 1.0f},
         {std::bit_cast<float>(0x3e600009U), std::bit_cast<float>(0x3e600008U),
          std::bit_cast<float>(0x3f894e83U)},
         {std::bit_cast<float>(0x3f894e83U), std::bit_cast<float>(0x3e600009U),
          std::bit_cast<float>(0x3e600008U)},
         {std::bit_cast<float>(0x3e600008U), std::bit_cast<float>(0x3f894e83U),
          std::bit_cast<float>(0x3e600009U)}}};
    std::vector<std::array<double, 24>> cases;
    for (int family = 0; family < 3; ++family) {
        for (const auto& direction : directions) {
            std::array<double, 24> input{family == 0 ? 0.0 : 1.0,
                                         family == 2 ? double(0.9f) : 0.0,
                                         0,
                                         0,
                                         0,
                                         4,
                                         1,
                                         0.5,
                                         -1,
                                         double(0.05f),
                                         double(0.02f),
                                         double(0.01f),
                                         0,
                                         -1,
                                         0,
                                         1,
                                         double(0.03f),
                                         double(0.1f),
                                         double(0.02f),
                                         double(-0.01f),
                                         direction[0],
                                         direction[1],
                                         direction[2],
                                         0};
            cases.push_back(input);
        }
        auto stationary = cases.back();
        stationary[17] = stationary[18] = stationary[19] = 0;
        cases.push_back(stationary);
    }
    return cases;
}

template <class Real>
void CheckCameraFrame(const std::string& artifact) {
    using namespace sirius::backend;
#ifdef SIRIUS_KERNEL_DIR
    const std::string directory = SIRIUS_KERNEL_DIR;
#else
    const std::string directory;
    GTEST_SKIP() << "Vulkan kernels are not configured";
#endif
    constexpr bool wide = std::is_same_v<Real, double>;
    // Finite moderate-frame witnesses; these are not runtime admission policy.
    // The wide budget matches the existing CPU camera-frame regressions.
    constexpr double budget = wide ? 3e-13 : 4e-6;
    const auto inventory = EnumerateVulkanDevices();
    ASSERT_TRUE(inventory.has_value()) << inventory.error().Description();
    if (inventory->empty()) GTEST_SKIP() << "no Vulkan device available";
    const auto index = ResolveVulkanDeviceIndex(*inventory);
    ASSERT_TRUE(index.has_value()) << index.error().Description();
    auto opened = CreateVulkanDevice(*index);
    ASSERT_TRUE(opened.has_value()) << opened.error().Description();
    auto device = std::move(*opened);
    const auto& info = device->Info();
    ::testing::Test::RecordProperty("camera_device_index", std::to_string(*index));
    ::testing::Test::RecordProperty("camera_device_name", info.name);
    ::testing::Test::RecordProperty("camera_device_kind", ToString(info.kind));
    ::testing::Test::RecordProperty("camera_vendor_id", std::to_string(info.vendor_id));
    ::testing::Test::RecordProperty("camera_device_id", std::to_string(info.device_id));
    ::testing::Test::RecordProperty("camera_driver_id", std::to_string(info.driver_id));
    ::testing::Test::RecordProperty("camera_driver_name", info.driver_name);
    ::testing::Test::RecordProperty("camera_driver_info", info.driver_info);
    ::testing::Test::RecordProperty("camera_supports_fp64", info.supports_fp64 ? "true" : "false");
    const auto words = ReadCameraProbe(directory + "/" + artifact);
    ASSERT_FALSE(words.empty()) << artifact;
    const auto kernel = device->LoadKernel(words);
    if constexpr (wide) {
        if (!info.supports_fp64) {
            ASSERT_FALSE(kernel.has_value());
            EXPECT_EQ(kernel.error().domain(), sirius::base::ErrorDomain::kKernel);
            EXPECT_NE(kernel.error().detail().find("shaderFloat64"), std::string::npos);
            ::testing::Test::RecordProperty("camera_evidence", "unsupported_fp64_declined");
            return;
        }
    }
    ASSERT_TRUE(kernel.has_value()) << kernel.error().Description();

    constexpr std::uint64_t allocation_limit = 1024 * 1024;
    constexpr std::uint64_t requested_bytes = 63 * sizeof(Real);
    ASSERT_TRUE(device->SetBufferAllocationLimit(allocation_limit).has_value());
    const auto input_buffer = device->CreateBuffer(24 * sizeof(Real), BufferUsage::kStorage);
    const auto output_buffer = device->CreateBuffer(39 * sizeof(Real), BufferUsage::kStorage);
    ASSERT_TRUE(input_buffer.has_value());
    ASSERT_TRUE(output_buffer.has_value());
    const std::array<BufferHandle, 2> bindings{*input_buffer, *output_buffer};
    ::testing::Test::RecordProperty("camera_allocated_bytes",
                                    std::to_string(device->BufferAllocationBytes()));
    ::testing::Test::RecordProperty("camera_requested_bytes", std::to_string(requested_bytes));
    ::testing::Test::RecordProperty("camera_allocation_limit_bytes",
                                    std::to_string(allocation_limit));
    ASSERT_GE(device->BufferAllocationBytes(), requested_bytes);
    ASSERT_LE(device->BufferAllocationBytes(), allocation_limit);
    std::size_t ordinal = 0;
    for (const auto& reference_input : CameraCases()) {
        SCOPED_TRACE(ordinal);
        std::array<Real, 24> input{};
        std::transform(reference_input.begin(), reference_input.end(), input.begin(),
                       [](double x) { return Real(x); });
        std::array<Real, 39> output;
        output.fill(std::numeric_limits<Real>::quiet_NaN());
        ASSERT_TRUE(
            device->WriteBuffer(*input_buffer, std::as_bytes(std::span(input))).has_value());
        ASSERT_TRUE(
            device->WriteBuffer(*output_buffer, std::as_bytes(std::span(output))).has_value());
        DispatchTiming timing;
        const auto dispatched = device->Dispatch(*kernel, bindings, 1, 1, 1, &timing);
        ASSERT_TRUE(dispatched.has_value()) << dispatched.error().Description();
        const std::string key = "camera_" + std::to_string(ordinal++);
        ::testing::Test::RecordProperty(key + "_submit_wait_ms",
                                        CameraNumber(timing.submit_wait_ms));
        ::testing::Test::RecordProperty(key + "_pipeline_setup_ms",
                                        CameraNumber(timing.pipeline_setup_ms));
        ::testing::Test::RecordProperty(key + "_total_ms", CameraNumber(timing.total_ms));
        ASSERT_TRUE(std::isfinite(timing.submit_wait_ms));
        ASSERT_GT(timing.submit_wait_ms, 0.0);
        ASSERT_LE(timing.submit_wait_ms, sirius::render::kDispatchStopMs);
        ASSERT_TRUE(device->ReadBuffer(*output_buffer, std::as_writable_bytes(std::span(output)))
                        .has_value());
        for (std::size_t j = 0; j < output.size(); ++j)
            ASSERT_TRUE(std::isfinite(output[j])) << "unwritten/nonfinite scalar " << j;
        ASSERT_EQ(output[32], Real(1));
        for (int j = 0; j < 4; ++j) ASSERT_EQ(output[j], input[j + 4]);
        // Use the exact returned event and represented ABI parameters.
        auto actual_input = reference_input;
        for (int j = 0; j < 4; ++j) actual_input[j + 4] = double(output[j]);
        const auto [g, inverse] = ReferenceMetric(actual_input);
        const auto dot = [&](int a, int b) {
            long double value = 0;
            for (int i = 0; i < 4; ++i)
                for (int j = 0; j < 4; ++j)
                    value += g[i][j] * static_cast<long double>(output[a + i]) * output[b + j];
            return value;
        };
        long double maximum_invariant = 0;
        const auto check = [&](long double value, long double expected) {
            const long double error = std::abs(value - expected);
            maximum_invariant = std::max(maximum_invariant, error);
            EXPECT_LE(error, budget);
        };
        for (int a = 0; a < 4; ++a)
            for (int b = 0; b < 4; ++b)
                check(dot(4 + 4 * a, 4 + 4 * b), a == b ? (a == 0 ? -1 : 1) : 0);
        check(dot(20, 20), 0);
        check(dot(20, 4), 1);  // Physical future photon is minus this past ray.
        for (int screen : {24, 28}) {
            check(dot(screen, screen), 1);
            check(dot(screen, 4), 0);
            check(dot(screen, 20), 0);
        }
        check(dot(24, 28), 0);
        ASSERT_GT(output[4], Real(0));

        // Invariants alone permit a wrongly rotated or relabelled frame.
        // Compare the actual labelled vectors with the independent CPU route.
        using namespace sirius::core;
        Metric4d cpu_metric, cpu_inverse;
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j) {
                cpu_metric(i, j) = Dual<double>(double(g[i][j]));
                cpu_inverse(i, j) = Dual<double>(double(inverse[i][j]));
            }
        std::array<Vec4, 3> seeds{};
        for (int axis = 0; axis < 3; ++axis)
            for (int j = 0; j < 3; ++j)
                seeds[axis](j + 1) = reference_input[8 + 3 * axis + j] * (axis == 1 ? -1 : 1);
        const auto frame = relativity::EulerianObserverFrame(cpu_metric, cpu_inverse, seeds);
        ASSERT_TRUE(frame.has_value());
        const auto boosted = relativity::BoostObserverFrame(
            *frame, {reference_input[17], reference_input[18], reference_input[19]});
        ASSERT_TRUE(boosted.has_value());
        const std::array<double, 3> direction{reference_input[20], reference_input[21],
                                              reference_input[22]};
        const auto catalogue = relativity::MakeCelestialTangentBasis<float>(
            {float(direction[0]), float(direction[1]), float(direction[2])});
        ASSERT_TRUE(catalogue.has_value());
        for (int j = 0; j < 3; ++j) {
            EXPECT_NEAR(output[33 + j], catalogue->first[j], 2e-6);
            EXPECT_NEAR(output[36 + j], catalogue->second[j], 2e-6);
            // Project the observed physical screens back into their own frame;
            // a wrong yet orthonormal reference axis must fail this comparison.
            EXPECT_NEAR(static_cast<double>(dot(24, 8 + 4 * j)), output[33 + j], 4e-6);
            EXPECT_NEAR(static_cast<double>(dot(28, 8 + 4 * j)), output[36 + j], 4e-6);
        }
        const auto ray = relativity::PastDirectedCameraRay(*boosted, direction);
        const auto screen = relativity::ObserverScreenBasis(*boosted, direction);
        ASSERT_TRUE(ray.has_value());
        ASSERT_TRUE(screen.has_value());
        const std::array<Vec4, 7> expected{
            boosted->time, boosted->spatial[0], boosted->spatial[1], boosted->spatial[2],
            *ray,          (*screen)[0],        (*screen)[1]};
        double maximum_component = 0;
        for (int v = 0; v < 7; ++v)
            for (int j = 0; j < 4; ++j) {
                const double error = std::abs(double(output[4 + 4 * v + j]) - expected[v](j));
                maximum_component = std::max(maximum_component, error);
                EXPECT_LE(error, budget) << "labelled vector " << v << " component " << j;
            }
        ::testing::Test::RecordProperty(key + "_maximum_invariant_error",
                                        CameraNumber(double(maximum_invariant)));
        ::testing::Test::RecordProperty(key + "_maximum_component_error",
                                        CameraNumber(maximum_component));
        std::cout << "[CameraFrame] " << key << " invariant=" << double(maximum_invariant)
                  << " component=" << maximum_component
                  << " submit_wait_ms=" << timing.submit_wait_ms << std::endl;
    }
    ASSERT_EQ(ordinal, 24U);
    ::testing::Test::RecordProperty("camera_actual_submissions", std::to_string(ordinal));
    ::testing::Test::RecordProperty("camera_evidence", "actual_frame_ray_and_labelled_screens");
}

TEST(KernelParity, CameraFramePreservesActivePrecisionFp32) {
    ASSERT_NO_FATAL_FAILURE(CheckCameraFrame<float>("camera_frame_probe.spv"));
}

TEST(KernelParity, CameraFramePreservesActivePrecisionCompensatedFp32) {
    ASSERT_NO_FATAL_FAILURE(CheckCameraFrame<float>("camera_frame_probe_fp32comp.spv"));
}

TEST(KernelParity, CameraFramePreservesActivePrecisionFp64) {
    ASSERT_NO_FATAL_FAILURE(CheckCameraFrame<double>("camera_frame_probe_fp64.spv"));
}

}  // namespace
