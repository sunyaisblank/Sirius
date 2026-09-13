#include "sirius/backend/retained_compute.h"

#include "sirius/core/twofold.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstring>
#include <limits>

#if defined(SIRIUS_HAS_RETAINED_COMPUTE) && defined(SIRIUS_RETAINED_CAMERA_TEST_DIR)
#define SIRIUS_RETAINED_TESTS_AVAILABLE 1
#include "program_fixture.h"
#include "support/retained_camera/continuous_reference.h"
#include "support/retained_transport/reference_cases.h"
#endif

namespace {
using namespace sirius::backend;

class RetainedComputeTest : public ::testing::Test {
  protected:
    void SetUp() override {
#ifdef SIRIUS_RETAINED_TESTS_AVAILABLE
        const auto inventory = EnumerateVulkanDevices();
        ASSERT_TRUE(inventory) << inventory.error().Description();
        if (inventory->empty()) GTEST_SKIP() << "No Vulkan device available";
        const auto index = ResolveVulkanDeviceIndex(*inventory);
        ASSERT_TRUE(index) << index.error().Description();
        auto opened = CreateVulkanDevice(*index);
        ASSERT_TRUE(opened) << opened.error().Description();
        device = std::move(*opened);
        ASSERT_TRUE(device->SetBufferAllocationLimit(8 * 1024 * 1024));
        auto created = RetainedCompute::Create(*device, 24);
        ASSERT_TRUE(created) << created.error().Description();
        compute = std::move(*created);
#else
        GTEST_SKIP() << "Retained compute build tools unavailable";
#endif
    }
    std::unique_ptr<ComputeDevice> device;
    std::unique_ptr<RetainedCompute> compute;
};

#ifdef SIRIUS_RETAINED_TESTS_AVAILABLE
bool Encloses(const RetainedValue& value, long double expected, long double reference_gap) {
    if (!value.IsRepresented()) return false;
    const long double center = (static_cast<long double>(value.high) + value.low) + value.tail;
    const long double rounding =
        16 * std::numeric_limits<long double>::epsilon() * (std::abs(center) + std::abs(expected));
    return std::abs(center - expected) <= value.radius + reference_gap + rounding;
}

template <typename Fixture>
bool CameraAgrees(const RetainedCameraOutput& output, const Fixture& fixture) {
    if (!output.valid) return false;
    constexpr std::array<std::size_t, 9> boundaries{0, 4, 8, 24, 40, 56, 72, 88, 104};
    for (std::size_t group = 1; group < boundaries.size(); ++group) {
        long double scale = 0, error = 0;
        for (std::size_t i = boundaries[group - 1]; i < boundaries[group]; ++i) {
            if (!Encloses(output.values[i], fixture.reference[i], fixture.reference_gap[i]))
                return false;
            const long double center =
                (static_cast<long double>(output.values[i].high) + output.values[i].low) +
                output.values[i].tail;
            scale = std::max(scale, std::abs(fixture.reference[i]));
            error = std::max(error, std::abs(center - fixture.reference[i]));
        }
        if (error > 1e-11L * scale) return false;
    }
    return true;
}

bool StepAgrees(const RetainedStepOutput& output,
                const sirius::test::retained_transport::Case& fixture) {
    if (!output.valid || output.stages != 7) return false;
    const std::array records{&output.fifth, &output.fourth, &output.increment, &output.error};
    for (std::size_t record = 0; record < records.size(); ++record)
        for (std::size_t group = 0; group < 10; ++group) {
            long double scale = 0, error = 0;
            for (std::size_t axis = 0; axis < 4; ++axis) {
                const auto i = group * 4 + axis, index = record * 40 + i;
                const auto& value = (*records[record])[i];
                if (!value.IsRepresented()) return false;
                // Preserve the third device term even on platforms where
                // long double has only binary64 precision. The frozen oracle
                // supplies an independent twofold binary64 conversion.
                const auto reference = fixture.wide_reference[index];
                const sirius::core::Twofold exact(reference.high, reference.low);
                const auto observed = sirius::core::Twofold(value.high) +
                                      sirius::core::Twofold(value.low) +
                                      sirius::core::Twofold(value.tail);
                const double difference = std::abs((observed - exact).Rounded());
                constexpr double epsilon = std::numeric_limits<double>::epsilon();
                const double rounding =
                    128 * epsilon * epsilon * (std::abs(reference.high) + std::abs(value.high));
                if (difference > value.radius + fixture.precision_gap[index] + rounding)
                    return false;
                const long double center =
                    (static_cast<long double>(value.high) + value.low) + value.tail;
                error = std::max(error, std::abs(center - fixture.reference[index]));
                // The small embedded difference is checked on the increment's
                // scale. Its finite arithmetic enclosure remains checked above.
                const auto scale_index = (record == 3 ? 2 : record) * 40 + i;
                scale = std::max(scale, std::abs(fixture.reference[scale_index]));
            }
            if (error > 1e-11L * scale) return false;
        }
    return true;
}

#endif

TEST_F(RetainedComputeTest, BatchedCameraPreservesPhysicalColumnsAndRejectsInvalidRows) {
#ifdef SIRIUS_RETAINED_TESTS_AVAILABLE
    const auto allocation = device->BufferAllocationBytes();
    std::vector<RetainedCameraInput> inputs;
    for (const auto& fixture : sirius::test::retained_camera::kCases) {
        RetainedCameraInput input;
        for (std::size_t i = 0; i < 32; ++i)
            input.values[i] = RetainedValue::FromDouble(std::bit_cast<float>(fixture.input[i]));
        inputs.push_back(input);
    }
    auto outputs = compute->Camera(inputs);
    ASSERT_TRUE(outputs) << outputs.error().Description();
    ASSERT_EQ(outputs->size(), inputs.size());
    for (std::size_t row = 0; row < inputs.size(); ++row) {
        SCOPED_TRACE(row);
        const auto& fixture = sirius::test::retained_camera::kCases[row];
        ASSERT_TRUE(CameraAgrees((*outputs)[row], fixture));
        auto narrowed = (*outputs)[row];
        for (auto& value : narrowed.values) value.low = 0;
        EXPECT_FALSE(CameraAgrees(narrowed, fixture));
    }
    inputs[2].values[20].high = std::numeric_limits<float>::quiet_NaN();
    inputs[6].values[29] = RetainedValue::FromDouble(2);
    outputs = compute->Camera(inputs);
    ASSERT_TRUE(outputs) << outputs.error().Description();
    for (std::size_t row = 0; row < inputs.size(); ++row) {
        if (row == 2 || row == 6) {
            EXPECT_FALSE((*outputs)[row].valid);
            for (const auto& value : (*outputs)[row].values) EXPECT_EQ(value.valid, 0U);
        } else {
            EXPECT_TRUE(CameraAgrees((*outputs)[row], sirius::test::retained_camera::kCases[row]));
        }
    }
    inputs.resize(1);
    outputs = compute->Camera(inputs);
    ASSERT_TRUE(outputs) << outputs.error().Description();
    ASSERT_EQ(outputs->size(), 1U);
    EXPECT_TRUE(CameraAgrees(outputs->front(), sirius::test::retained_camera::kCases.front()));
    inputs.clear();
    for (const auto& fixture : sirius::test::continuous_retained_camera::cases) {
        RetainedCameraInput input;
        for (std::size_t i = 0; i < 32; ++i)
            input.values[i] = RetainedValue::FromDouble(fixture.input[i]);
        inputs.push_back(input);
    }
    ASSERT_EQ(inputs[0].values[25].high, inputs[1].values[25].high);
    ASSERT_NE(inputs[0].values[25].low, inputs[1].values[25].low);
    ASSERT_EQ(inputs[0].values[27].high, inputs[2].values[27].high);
    ASSERT_NE(inputs[0].values[27].low, inputs[2].values[27].low);
    outputs = compute->Camera(inputs);
    ASSERT_TRUE(outputs) << outputs.error().Description();
    for (std::size_t row = 0; row < inputs.size(); ++row)
        EXPECT_TRUE(
            CameraAgrees((*outputs)[row], sirius::test::continuous_retained_camera::cases[row]));
    for (auto& input : inputs)
        for (auto& value : input.values) {
            value.low = 0;
            value.tail = 0;
        }
    outputs = compute->Camera(inputs);
    ASSERT_TRUE(outputs) << outputs.error().Description();
    EXPECT_FALSE(CameraAgrees((*outputs)[1], sirius::test::continuous_retained_camera::cases[1]));
    EXPECT_FALSE(CameraAgrees((*outputs)[2], sirius::test::continuous_retained_camera::cases[2]));
    EXPECT_EQ(device->BufferAllocationBytes(), allocation);
#else
    GTEST_SKIP() << "Retained compute build tools unavailable";
#endif
}

TEST_F(RetainedComputeTest, JointRkStagesRetainCriticalIncrementsAndEmbeddedError) {
#ifdef SIRIUS_RETAINED_TESTS_AVAILABLE
    const auto allocation = device->BufferAllocationBytes();
    std::vector<RetainedStepInput> inputs;
    for (const auto& fixture : sirius::test::retained_transport::cases) {
        inputs.push_back(std::bit_cast<RetainedStepInput>(fixture.input));
    }
    const auto outputs = compute->Step(inputs);
    ASSERT_TRUE(outputs) << outputs.error().Description();
    ASSERT_EQ(outputs->size(), inputs.size());
    for (std::size_t row = 0; row < inputs.size(); ++row) {
        SCOPED_TRACE(sirius::test::retained_transport::cases[row].name);
        ASSERT_TRUE(StepAgrees((*outputs)[row], sirius::test::retained_transport::cases[row]));
        if (row + 1 < inputs.size()) {
            auto two_term = (*outputs)[row];
            for (auto* record :
                 {&two_term.fifth, &two_term.fourth, &two_term.increment, &two_term.error})
                for (auto& value : *record) value.tail = 0;
            EXPECT_FALSE(StepAgrees(two_term, sirius::test::retained_transport::cases[row]));
            auto narrowed = (*outputs)[row];
            for (auto* record :
                 {&narrowed.fifth, &narrowed.fourth, &narrowed.increment, &narrowed.error})
                for (auto& value : *record) value.low = 0;
            EXPECT_FALSE(StepAgrees(narrowed, sirius::test::retained_transport::cases[row]));
        }
    }
    EXPECT_EQ(device->BufferAllocationBytes(), allocation);
    const auto& flat = outputs->back();
    for (std::size_t i = 0; i < 40; ++i) {
        EXPECT_EQ(flat.error[i].Center(), 0);
        EXPECT_EQ(flat.fifth[i].Center(), flat.fourth[i].Center());
    }
#else
    GTEST_SKIP() << "Retained compute build tools unavailable";
#endif
}

TEST_F(RetainedComputeTest, RejectedStepRowsCannotExposeOldOrPartialCandidates) {
#ifdef SIRIUS_RETAINED_TESTS_AVAILABLE
    const auto good =
        std::bit_cast<RetainedStepInput>(sirius::test::retained_transport::cases.front().input);
    std::vector<RetainedStepInput> inputs(9, good);
    auto outputs = compute->Step(inputs);
    ASSERT_TRUE(outputs) << outputs.error().Description();
    for (const auto& output : *outputs) ASSERT_TRUE(output.valid);
    inputs[0].values[4].high = std::numeric_limits<float>::quiet_NaN();
    inputs[1].values[5].low = std::numeric_limits<float>::infinity();
    inputs[2].values[6].radius = -1;
    inputs[3].values[7].valid = 0;
    inputs[4].values[44] = RetainedValue::FromDouble(0);
    inputs[5].values[45] = RetainedValue::FromDouble(0);
    inputs[6].values[45] = RetainedValue::FromDouble(-1);
    inputs[7].values[3] = RetainedValue::FromDouble(.01);
    for (std::size_t i = 5; i < 8; ++i) inputs[8].values[i] = RetainedValue::FromDouble(0);
    outputs = compute->Step(inputs);
    ASSERT_TRUE(outputs) << outputs.error().Description();
    for (const auto& output : *outputs) {
        EXPECT_FALSE(output.valid);
        for (const auto* record : {&output.fifth, &output.fourth, &output.increment, &output.error})
            for (const auto& value : *record) EXPECT_EQ(value.valid, 0U);
    }
    EXPECT_EQ(outputs->back().stages, 1U);
    const auto recovered = compute->Step(std::span(&good, 1));
    ASSERT_TRUE(recovered) << recovered.error().Description();
    ASSERT_EQ(recovered->size(), 1U);
    EXPECT_TRUE(StepAgrees(recovered->front(), sirius::test::retained_transport::cases.front()));
#else
    GTEST_SKIP() << "Retained compute build tools unavailable";
#endif
}

TEST(RetainedValue, Binary64InputsKeepTheirRepresentationResidual) {
#ifdef SIRIUS_RETAINED_TESTS_AVAILABLE
    for (const double input : {0.0, -0.0, 1.0 / 3.0, -1.0 / 7.0, 0x1.123456789abcdp80,
                               0x1.123456789abcdp-120, 0x1p-149, 0x1p-150, 0x1p-1022}) {
        const auto value = RetainedValue::FromDouble(input);
        ASSERT_TRUE(value.IsRepresented());
        EXPECT_TRUE(Encloses(value, input, 0));
    }
    EXPECT_FALSE(
        RetainedValue::FromDouble(std::numeric_limits<double>::infinity()).IsRepresented());
    EXPECT_FALSE(RetainedValue::FromDouble(0x1p121).IsRepresented());
    EXPECT_FALSE((RetainedValue{1, 1, 0, 0, 1}.IsRepresented()));
#else
    GTEST_SKIP() << "Retained compute build tools unavailable";
#endif
}
}  // namespace
