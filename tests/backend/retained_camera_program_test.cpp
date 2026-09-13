#include "sirius/backend/device.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <span>
#include <string>
#include <vector>

#ifdef SIRIUS_RETAINED_CAMERA_TEST_DIR
#include "program_fixture.h"

namespace {
using namespace sirius::test::retained_camera;

float Value(std::uint32_t bits) { return std::bit_cast<float>(bits); }

// The independent reference's 100/180-digit gap is a finite witness, not a
// certified interval. Include conversion/accumulation rounding at this host
// comparison boundary; never discard the device's low part or error radius.
bool Agrees(const Case& fixture, const std::vector<std::uint32_t>& words) {
    if (words.size() != kOutputWords || Value(words[3]) != 1) return false;
    std::array<long double, 104> errors{};
    for (std::size_t i = 0; i < 104; ++i) {
        const long double hi = Value(words[40 + 3 * i]);
        const long double lo = Value(words[41 + 3 * i]);
        const long double radius = Value(words[42 + 3 * i]);
        if (!std::isfinite(hi) || !std::isfinite(lo) || !std::isfinite(radius) || radius < 0)
            return false;
        errors[i] = std::abs((hi + lo) - fixture.reference[i]);
        const long double rounding = 16 * std::numeric_limits<long double>::epsilon() *
                                     (std::abs(hi) + std::abs(lo) + std::abs(fixture.reference[i]));
        if (errors[i] > radius + fixture.reference_gap[i] + rounding) return false;
    }
    // A large radius cannot turn an inaccurate result into successful launch
    // admission. Apply the unchanged 1e-11 stage target separately to x, k,
    // the frame, X, coordinate K, covariant V, observer du and repeated frame.
    constexpr std::array<std::size_t, 9> boundaries{0, 4, 8, 24, 40, 56, 72, 88, 104};
    for (std::size_t group = 1; group < boundaries.size(); ++group) {
        long double scale = 0, error = 0;
        for (std::size_t i = boundaries[group - 1]; i < boundaries[group]; ++i) {
            scale = std::max(scale, std::abs(fixture.reference[i]));
            error = std::max(error, errors[i]);
        }
        if (error > 1e-11L * scale) return false;
    }
    return true;
}

void CheckProgram(bool compensated) {
    using namespace sirius::backend;
    const auto inventory = EnumerateVulkanDevices();
    ASSERT_TRUE(inventory) << inventory.error().Description();
    if (inventory->empty()) GTEST_SKIP() << "No Vulkan device available";
    const auto selected = ResolveVulkanDeviceIndex(*inventory);
    ASSERT_TRUE(selected) << selected.error().Description();
    auto opened = CreateVulkanDevice(*selected);
    ASSERT_TRUE(opened) << opened.error().Description();
    auto device = std::move(*opened);
    ::testing::Test::RecordProperty("retained_camera_device", device->Info().name);
    ::testing::Test::RecordProperty("retained_camera_driver", device->Info().driver_info);
    ::testing::Test::RecordProperty("retained_camera_registers", std::to_string(kRegisters));
    ::testing::Test::RecordProperty("retained_camera_instructions", std::to_string(kInstructions));
    const std::string shader =
        std::string(SIRIUS_RETAINED_CAMERA_TEST_DIR) +
        (compensated ? "/program_camera_probe-fp32comp.spv" : "/program_camera_probe-fp32.spv");
    std::ifstream file(shader, std::ios::binary | std::ios::ate);
    ASSERT_TRUE(file) << shader;
    const auto size = file.tellg();
    ASSERT_GT(size, 0);
    ASSERT_EQ(size % 4, 0);
    std::vector<std::uint32_t> code(static_cast<std::size_t>(size) / 4);
    file.seekg(0);
    file.read(reinterpret_cast<char*>(code.data()), size);
    ASSERT_TRUE(file);
    const auto kernel = device->LoadKernel(code);
    ASSERT_TRUE(kernel) << kernel.error().Description();
    ASSERT_TRUE(device->SetBufferAllocationLimit(256 * 1024));
    const auto in = device->CreateBuffer((32 + kProgram.size()) * 4, BufferUsage::kStorage);
    const auto out = device->CreateBuffer(kOutputWords * 4, BufferUsage::kStorage);
    ASSERT_TRUE(in) << in.error().Description();
    ASSERT_TRUE(out) << out.error().Description();
    const std::array<BufferHandle, 2> buffers{*in, *out};
    std::vector<std::uint32_t> packet(32 + kProgram.size());
    std::copy(kProgram.begin(), kProgram.end(), packet.begin() + 32);
    std::vector<std::uint32_t> result(kOutputWords);
    double largest_submit = 0, largest_preparation = 0;
    const auto dispatch = [&] {
        std::fill(result.begin(), result.end(), 0x7fc00000U);
        auto written = device->WriteBuffer(*in, std::as_bytes(std::span(packet)));
        if (!written) return false;
        written = device->WriteBuffer(*out, std::as_bytes(std::span(result)));
        if (!written) return false;
        DispatchTiming timing;
        const auto dispatched = device->Dispatch(*kernel, buffers, 1, 1, 1, &timing);
        if (!dispatched) return false;
        largest_submit = std::max(largest_submit, timing.submit_wait_ms);
        largest_preparation = std::max(largest_preparation, timing.pipeline_setup_ms);
        return device->ReadBuffer(*out, std::as_writable_bytes(std::span(result))).has_value();
    };
    std::size_t component_controls = 0, low_part_controls = 0;
    for (const auto& fixture : kCases) {
        SCOPED_TRACE(fixture.name);
        std::copy(fixture.input.begin(), fixture.input.end(), packet.begin());
        ASSERT_TRUE(dispatch());
        ASSERT_EQ(Value(result[0]), 271832);
        ASSERT_EQ(Value(result[2]), compensated ? 1 : 0);
        ASSERT_EQ(result[4], fixture.input[30]);
        ASSERT_EQ(Value(result[5]), packet.size());
        ASSERT_EQ(Value(result[6]), result.size());
        for (std::size_t i = 0; i < 32; ++i) ASSERT_EQ(result[8 + i], fixture.input[i]);
        ASSERT_TRUE(Agrees(fixture, result)) << "complete retained launch differs from reference";
        auto narrowed = result;
        for (std::size_t i = 0; i < 104; ++i) narrowed[41 + 3 * i] = 0;
        EXPECT_FALSE(Agrees(fixture, narrowed)) << "discarded low parts went undetected";
        ++low_part_controls;
        for (std::size_t i = 0; i < 104; ++i) {
            auto mutated = result;
            const float hi = Value(mutated[40 + 3 * i]);
            mutated[40 + 3 * i] =
                std::bit_cast<std::uint32_t>(hi + .01f * std::max(1.0f, std::abs(hi)));
            ASSERT_FALSE(Agrees(fixture, mutated))
                << "component " << i << " mutation was invisible";
            ++component_controls;
        }
    }
    // Invalid inputs and malformed late output ownership must leave the whole
    // scientific region empty, even after an earlier valid dispatch.
    for (int invalid = 0; invalid < 8; ++invalid) {
        SCOPED_TRACE(invalid);
        std::copy(kCases[4].input.begin(), kCases[4].input.end(), packet.begin());
        std::copy(kProgram.begin(), kProgram.end(), packet.begin() + 32);
        switch (invalid) {
            case 0:
                packet[20] = 0x7fc00000U;
                break;
            case 1:
                packet[20] = 0;
                break;
            case 2:
                packet[17] = std::bit_cast<std::uint32_t>(2.0f);
                break;
            case 3:
                packet[3] = std::bit_cast<std::uint32_t>(.001f);
                break;
            case 4:
                packet[29] = 0;
                break;  // Nonzero pupil on a pinhole.
            case 5:
                packet[31] = std::bit_cast<std::uint32_t>(1.0f);
                break;
            case 6:
                packet[34 + 103] = static_cast<std::uint32_t>(kRegisters);
                break;
            case 7:
                packet[138] = 12;
                break;
        }
        ASSERT_TRUE(dispatch());
        ASSERT_EQ(result[3], 0U);
        for (std::size_t i = 40; i < 352; ++i) ASSERT_EQ(result[i], 0U) << i;
    }
    EXPECT_EQ(component_controls, 2080U);
    EXPECT_EQ(low_part_controls, 20U);
    std::cout << "[RetainedCamera] registers=" << kRegisters << " instructions=" << kInstructions
              << " buffers=" << device->BufferAllocationBytes()
              << " max_submit_ms=" << largest_submit << " max_pipeline_ms=" << largest_preparation
              << '\n';
}
}  // namespace
#endif

TEST(RetainedCameraProgram, Fp32CompletePhysicalLaunchAndRefusals) {
#ifdef SIRIUS_RETAINED_CAMERA_TEST_DIR
    ASSERT_NO_FATAL_FAILURE(CheckProgram(false));
#else
    GTEST_SKIP() << "Retained camera compiler and SPIR-V tools are not configured";
#endif
}

TEST(RetainedCameraProgram, CompensatedCompletePhysicalLaunchAndRefusals) {
#ifdef SIRIUS_RETAINED_CAMERA_TEST_DIR
    ASSERT_NO_FATAL_FAILURE(CheckProgram(true));
#else
    GTEST_SKIP() << "Retained camera compiler and SPIR-V tools are not configured";
#endif
}
