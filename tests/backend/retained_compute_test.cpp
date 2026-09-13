#include "sirius/backend/retained_compute.h"

#include "sirius/backend/cpu/geodesic_tracer.h"
#include "sirius/backend/retained_integrator.h"
#include "sirius/backend/retained_trace_executor.h"
#include "sirius/core/twofold.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <bit>
#include <cmath>
#include <cstring>
#include <future>
#include <iomanip>
#include <limits>
#include <numbers>

#if defined(SIRIUS_HAS_RETAINED_COMPUTE) && defined(SIRIUS_RETAINED_CAMERA_TEST_DIR)
#define SIRIUS_RETAINED_TESTS_AVAILABLE 1
#include "program_fixture.h"
#include "support/cpu_critical/transport_reference.h"
#include "support/retained_camera/continuous_reference.h"
#include "support/retained_camera/ray_reference.h"
#include "support/retained_transport/dense_reference.h"
#include "support/retained_transport/endpoint_reference.h"
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
::testing::AssertionResult CameraAgrees(const RetainedCameraOutput& output,
                                        const Fixture& fixture) {
    if (!output.valid) return ::testing::AssertionFailure() << "invalid camera row";
    constexpr std::array<std::size_t, 9> boundaries{0, 4, 8, 24, 40, 56, 72, 88, 104};
    for (std::size_t group = 1; group < boundaries.size(); ++group) {
        long double scale = 0, error = 0;
        for (std::size_t i = boundaries[group - 1]; i < boundaries[group]; ++i) {
            if (!Encloses(output.values[i], fixture.reference[i], fixture.reference_gap[i]))
                return ::testing::AssertionFailure()
                       << "camera component " << i << " center=" << std::setprecision(20)
                       << output.values[i].Center() << " reference=" << fixture.reference[i]
                       << " radius=" << output.values[i].radius;
            const long double center =
                (static_cast<long double>(output.values[i].high) + output.values[i].low) +
                output.values[i].tail;
            scale = std::max(scale, std::abs(fixture.reference[i]));
            error = std::max(error, std::abs(center - fixture.reference[i]));
        }
        if (error > 1e-11L * scale) return ::testing::AssertionFailure();
    }
    return ::testing::AssertionSuccess();
}

::testing::AssertionResult StepAgrees(const RetainedStepOutput& output,
                                      const sirius::test::retained_transport::Case& fixture) {
    if (!output.valid || output.stages != 7) return ::testing::AssertionFailure();
    const std::array records{&output.fifth, &output.fourth, &output.increment, &output.error};
    for (std::size_t record = 0; record < records.size(); ++record)
        for (std::size_t group = 0; group < 10; ++group) {
            long double scale = 0, error = 0;
            for (std::size_t axis = 0; axis < 4; ++axis) {
                const auto i = group * 4 + axis, index = record * 40 + i;
                const auto& value = (*records[record])[i];
                if (!value.IsRepresented()) return ::testing::AssertionFailure();
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
                    return ::testing::AssertionFailure()
                           << "step component " << index << " difference=" << difference
                           << " radius=" << value.radius;
                const long double center =
                    (static_cast<long double>(value.high) + value.low) + value.tail;
                error = std::max(error, std::abs(center - fixture.reference[index]));
                // The small embedded difference is checked on the increment's
                // scale. Its finite arithmetic enclosure remains checked above.
                const auto scale_index = (record == 3 ? 2 : record) * 40 + i;
                scale = std::max(scale, std::abs(fixture.reference[scale_index]));
            }
            if (error > 1e-11L * scale) return ::testing::AssertionFailure();
        }
    return ::testing::AssertionSuccess();
}

#endif

TEST_F(RetainedComputeTest, BatchedCameraPreservesPhysicalColumnsAndRejectsInvalidRows) {
#ifdef SIRIUS_RETAINED_TESTS_AVAILABLE
    const auto allocation = device->BufferAllocationBytes();
    EXPECT_FALSE(RetainedCompute::Create(*device, 0));
    EXPECT_FALSE(RetainedCompute::Create(*device, 65536));
    EXPECT_EQ(device->BufferAllocationBytes(), allocation);
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

TEST_F(RetainedComputeTest, SmoothRayCameraPreservesPhysicalLensDerivatives) {
#ifdef SIRIUS_RETAINED_TESTS_AVAILABLE
    const auto allocation = device->BufferAllocationBytes();
    std::vector<RetainedRayCameraInput> inputs;
    for (const auto& fixture : sirius::test::retained_ray_camera::cases)
        inputs.push_back(std::bit_cast<RetainedRayCameraInput>(fixture.input));
    auto outputs = compute->RayCamera(inputs);
    ASSERT_TRUE(outputs) << outputs.error().Description();
    ASSERT_EQ(outputs->size(), inputs.size());
    for (std::size_t row = 0; row < inputs.size(); ++row) {
        SCOPED_TRACE(sirius::test::retained_ray_camera::cases[row].name);
        ASSERT_TRUE(CameraAgrees((*outputs)[row], sirius::test::retained_ray_camera::cases[row]));
        auto narrowed = (*outputs)[row];
        for (auto& value : narrowed.values) value.low = 0;
        EXPECT_FALSE(CameraAgrees(narrowed, sirius::test::retained_ray_camera::cases[row]));
    }
    inputs[2].values[22].high = std::numeric_limits<float>::quiet_NaN();
    outputs = compute->RayCamera(inputs);
    ASSERT_TRUE(outputs) << outputs.error().Description();
    for (std::size_t row = 0; row < inputs.size(); ++row) {
        if (row == 2) {
            EXPECT_FALSE((*outputs)[row].valid);
            for (const auto& value : (*outputs)[row].values) EXPECT_EQ(value.valid, 0U);
        } else {
            EXPECT_TRUE(
                CameraAgrees((*outputs)[row], sirius::test::retained_ray_camera::cases[row]));
        }
    }
    inputs.resize(1);
    outputs = compute->RayCamera(inputs);
    ASSERT_TRUE(outputs) << outputs.error().Description();
    EXPECT_TRUE(CameraAgrees(outputs->front(), sirius::test::retained_ray_camera::cases.front()));
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

TEST_F(RetainedComputeTest, ProjectedEndpointsKeepPhysicalColumnsAndRetainedContinuation) {
#ifdef SIRIUS_RETAINED_TESTS_AVAILABLE
    const auto allocation = device->BufferAllocationBytes();
    std::vector<RetainedEndpointInput> inputs;
    for (const auto& fixture : sirius::test::retained_endpoint::cases)
        inputs.push_back(std::bit_cast<RetainedEndpointInput>(fixture.input));
    const auto outputs = compute->Endpoint(inputs);
    ASSERT_TRUE(outputs) << outputs.error().Description();
    for (std::size_t row = 0; row < inputs.size(); ++row) {
        const auto& fixture = sirius::test::retained_endpoint::cases[row];
        SCOPED_TRACE(fixture.name);
        const auto& output = (*outputs)[row];
        ASSERT_TRUE(output.valid);
        EXPECT_EQ(output.component, fixture.component);
        for (std::size_t i = 0; i < 80; ++i) {
            SCOPED_TRACE(i);
            const auto& value = i < 40 ? output.phase[i] : output.physical[i - 40];
            const auto oracle = fixture.reference[i];
            const sirius::core::Twofold exact(oracle.high, oracle.low);
            const auto center = sirius::core::Twofold(value.high) +
                                sirius::core::Twofold(value.low) +
                                sirius::core::Twofold(value.tail);
            const double difference = std::abs((center - exact).Rounded());
            EXPECT_LE(difference, value.radius + 1e-29 * (1 + std::abs(oracle.high)));
            EXPECT_LE(difference, 1e-11 * (1 + std::abs(oracle.high)));
        }
    }
    // Reuse the accepted phase without converting through the physical output
    // or a binary64 sum. Projection must not move an already projected state.
    for (std::size_t row = 0; row < inputs.size(); ++row)
        std::copy((*outputs)[row].phase.begin(), (*outputs)[row].phase.end(),
                  inputs[row].values.begin() + 4);
    const auto repeated = compute->Endpoint(inputs);
    ASSERT_TRUE(repeated) << repeated.error().Description();
    for (std::size_t row = 0; row < inputs.size(); ++row) {
        ASSERT_TRUE((*repeated)[row].valid) << row;
        for (std::size_t i = 0; i < 40; ++i)
            EXPECT_NEAR((*repeated)[row].physical[i].Center(), (*outputs)[row].physical[i].Center(),
                        1e-11 * (1 + std::abs((*outputs)[row].physical[i].Center())));
    }
    inputs[0].values[4].valid = 0;
    inputs[1].values[8] = RetainedValue::FromDouble(100);
    const auto invalid = compute->Endpoint(inputs);
    ASSERT_TRUE(invalid) << invalid.error().Description();
    for (std::size_t row = 0; row < 2; ++row) {
        EXPECT_FALSE((*invalid)[row].valid);
        for (const auto& value : (*invalid)[row].physical) EXPECT_EQ(value.valid, 0U);
        for (const auto& value : (*invalid)[row].phase) EXPECT_EQ(value.valid, 0U);
    }
    EXPECT_EQ(device->BufferAllocationBytes(), allocation);
#else
    GTEST_SKIP() << "Retained compute build tools unavailable";
#endif
}

TEST_F(RetainedComputeTest, DenseSegmentsPreserveSmallCovariantArrivalDerivatives) {
#ifdef SIRIUS_RETAINED_TESTS_AVAILABLE
    const auto allocation = device->BufferAllocationBytes();
    const auto& cases = sirius::test::retained_dense::cases;
    for (std::size_t begin = 0; begin < cases.size(); begin += compute->Capacity()) {
        const auto count = std::min(compute->Capacity(), cases.size() - begin);
        std::vector<RetainedDenseInput> inputs;
        for (std::size_t i = 0; i < count; ++i)
            inputs.push_back(std::bit_cast<RetainedDenseInput>(cases[begin + i].input));
        const auto outputs = compute->Dense(inputs);
        ASSERT_TRUE(outputs) << outputs.error().Description();
        for (std::size_t row = 0; row < count; ++row) {
            const auto& fixture = cases[begin + row];
            SCOPED_TRACE(fixture.name);
            ASSERT_TRUE((*outputs)[row].valid);
            for (std::size_t i = 0; i < 40; ++i) {
                SCOPED_TRACE(i);
                const auto& value = (*outputs)[row].physical[i];
                const auto oracle = fixture.reference[i];
                const auto center = sirius::core::Twofold(value.high) +
                                    sirius::core::Twofold(value.low) +
                                    sirius::core::Twofold(value.tail);
                const double difference =
                    std::abs((center - sirius::core::Twofold(oracle.high, oracle.low)).Rounded());
                EXPECT_LE(difference, value.radius + 1e-28 * (1 + std::abs(oracle.high)));
                EXPECT_LE(difference, 1e-10 * (1 + std::abs(oracle.high)));
            }
        }
    }
    auto input = std::bit_cast<RetainedDenseInput>(cases[1].input);
    for (std::size_t i = 106; i < 110; ++i) input.values[i] = RetainedValue::FromDouble(0);
    auto invalid = compute->Dense(std::span(&input, 1));
    ASSERT_TRUE(invalid) << invalid.error().Description();
    EXPECT_FALSE(invalid->front().valid);
    for (const auto& value : invalid->front().physical) EXPECT_EQ(value.valid, 0U);
    input = std::bit_cast<RetainedDenseInput>(cases[0].input);
    input.values[105] = RetainedValue::FromDouble(1.01);
    invalid = compute->Dense(std::span(&input, 1));
    ASSERT_TRUE(invalid) << invalid.error().Description();
    EXPECT_FALSE(invalid->front().valid);
    EXPECT_EQ(device->BufferAllocationBytes(), allocation);
#else
    GTEST_SKIP() << "Retained compute build tools unavailable";
#endif
}

TEST_F(RetainedComputeTest, PhysicalInitializationRetainsTheHamiltonianResidual) {
#ifdef SIRIUS_RETAINED_TESTS_AVAILABLE
    const auto& initial = sirius::test::critical_fixture::transport_cases;
    std::vector<RetainedInitializeInput> inputs(initial.size());
    for (std::size_t row = 0; row < initial.size(); ++row) {
        auto& input = inputs[row];
        const auto original =
            std::bit_cast<RetainedStepInput>(sirius::test::retained_transport::cases[row].input);
        std::copy_n(original.values.begin(), 4, input.values.begin());
        for (std::size_t i = 0; i < 40; ++i)
            input.values[4 + i] = RetainedValue::FromDouble(initial[row].initial[i]);
        input.values[44] = RetainedValue::FromDouble(-1);
    }
    const auto outputs = compute->Initialize(inputs);
    ASSERT_TRUE(outputs) << outputs.error().Description();
    for (std::size_t row = 0; row < inputs.size(); ++row) {
        SCOPED_TRACE(row);
        ASSERT_TRUE((*outputs)[row].valid);
        const auto expected =
            std::bit_cast<RetainedStepInput>(sirius::test::retained_transport::cases[row].input);
        for (std::size_t i = 0; i < 40; ++i) {
            SCOPED_TRACE(i);
            const auto& actual = (*outputs)[row].phase[i];
            const auto& reference = expected.values[i + 4];
            const auto center = [](const RetainedValue& value) {
                return sirius::core::Twofold(value.high) + sirius::core::Twofold(value.low) +
                       sirius::core::Twofold(value.tail);
            };
            const double difference = std::abs((center(actual) - center(reference)).Rounded());
            EXPECT_LE(difference, double(actual.radius) + reference.radius +
                                      1e-29 * (1 + std::abs(reference.high)));
            EXPECT_LE(difference, 1e-11 * (1 + std::abs(reference.high)));
        }
    }
    inputs[0].values[16].valid = 0;
    const auto invalid = compute->Initialize(std::span(inputs.data(), 1));
    ASSERT_TRUE(invalid) << invalid.error().Description();
    EXPECT_FALSE(invalid->front().valid);
    for (const auto& value : invalid->front().phase) EXPECT_EQ(value.valid, 0U);
#else
    GTEST_SKIP() << "Retained compute build tools unavailable";
#endif
}

TEST_F(RetainedComputeTest, CoupledIntervalsRequireEmbeddedAndIndependentDenseAgreement) {
#ifdef SIRIUS_RETAINED_TESTS_AVAILABLE
    std::vector<RetainedEndpointInput> initial;
    for (const auto& fixture : sirius::test::retained_transport::cases) {
        const auto step = std::bit_cast<RetainedStepInput>(fixture.input);
        RetainedEndpointInput input;
        std::copy_n(step.values.begin(), 45, input.values.begin());
        initial.push_back(input);
    }
    const auto launches = compute->Endpoint(initial);
    ASSERT_TRUE(launches) << launches.error().Description();
    std::vector<RetainedIntervalInput> inputs(initial.size());
    for (std::size_t row = 0; row < inputs.size(); ++row) {
        ASSERT_TRUE((*launches)[row].valid) << row;
        const auto step =
            std::bit_cast<RetainedStepInput>(sirius::test::retained_transport::cases[row].input);
        auto& input = inputs[row];
        std::copy_n(step.values.begin(), 4, input.metric.begin());
        input.start = (*launches)[row];
        input.chart = step.values[44].Center();
        input.interval = step.values[45].Center();
        input.control.length_scale = 1;
        input.control.frequency_scale = 1;
        input.control.tolerance = 1e-4 / (4 * 30000);
        input.control.integrator.min_step = static_cast<float>(input.interval);
        input.control.integrator.max_step = 1;
        input.control.integrator.abs_tolerance = 1e-9f;
        input.control.integrator.rel_tolerance = 1e-9f;
    }
    const auto outputs = AttemptRetainedIntervals(*compute, inputs);
    ASSERT_TRUE(outputs) << outputs.error().Description();
    for (std::size_t row = 0; row < inputs.size(); ++row) {
        SCOPED_TRACE(sirius::test::retained_transport::cases[row].name);
        const auto& output = (*outputs)[row];
        EXPECT_TRUE(output.admissible)
            << "error=" << output.error_ratio
            << " failure=" << sirius::core::CoupledStepFailureName(output.failure);
        EXPECT_EQ(output.attempted_stages, 21U);
        EXPECT_LE(output.error_ratio, 1);
    }
    // Exact flat flow remains exact through all three independent candidates.
    const auto& flat = outputs->back();
    ASSERT_TRUE(flat.admissible);
    EXPECT_EQ(flat.error_ratio, 0);
    EXPECT_EQ(flat.full.physical[0].Center(), .75);
    EXPECT_EQ(flat.full.physical[3].Center(), 4.25);
    EXPECT_EQ(flat.midpoint.physical[0].Center(), .875);
    EXPECT_EQ(flat.refined.physical[3].Center(), 4.25);

    // Sparse retained tails must survive the lower-order increment split even
    // when the embedded error is exactly zero in flat space.
    auto sparse = inputs.back();
    sparse.start.phase[13] = sparse.start.physical[13] = {1, 0x1.000002p-35f, 0x1p-120f, 0, 1};
    const auto sparse_output = AttemptRetainedIntervals(*compute, {&sparse, 1});
    ASSERT_TRUE(sparse_output) << sparse_output.error().Description();
    ASSERT_TRUE(sparse_output->front().admissible);
    const auto& full_increment = sparse_output->front().full_increment[5];
    const auto& lower_increment = sparse_output->front().lower_increment[5];
    EXPECT_EQ(full_increment.tail, 0x1p-122f);
    EXPECT_EQ(lower_increment.high, full_increment.high);
    EXPECT_EQ(lower_increment.low, full_increment.low);
    EXPECT_EQ(lower_increment.tail, full_increment.tail);

    inputs.resize(2);
    inputs[0].control.tolerance = 1e-30;
    inputs[1].start.phase[5].valid = 0;
    const auto rejected = AttemptRetainedIntervals(*compute, inputs);
    ASSERT_TRUE(rejected) << rejected.error().Description();
    for (const auto& output : *rejected) {
        EXPECT_FALSE(output.admissible);
        EXPECT_NE(output.failure, sirius::core::CoupledStepFailure::None);
        EXPECT_FALSE(output.full.valid);
        EXPECT_FALSE(output.lower.valid);
        EXPECT_FALSE(output.midpoint.valid);
        EXPECT_FALSE(output.refined.valid);
        for (const auto& value : output.full_increment) EXPECT_EQ(value.valid, 0U);
    }
#else
    GTEST_SKIP() << "Retained compute build tools unavailable";
#endif
}

TEST_F(RetainedComputeTest, SharedTracerCompletesDeviceIntervalsAndRetainsRollbackState) {
#ifdef SIRIUS_RETAINED_TESTS_AVAILABLE
    std::atomic<bool> cancelled{false};
    RetainedTraceExecutor executor(*compute, [&] { return cancelled.load(); });
    for (const double mass : {0.0, 1.0}) {
        sirius::core::KerrSchildFamily metric(mass == 0
                                                  ? sirius::core::KerrSchildParams::Minkowski()
                                                  : sirius::core::KerrSchildParams::Kerr(1, .7));
        sirius::core::CameraConfig config;
        config.width = 32;
        config.height = 20;
        config.theta = 1.1;
        config.phi = .2;
        config.beta_x = .1;
        config.beta_y = .8;
        config.beta_z = .01;
        config.focus_distance = 50;
        sirius::core::ThinLensCamera camera(config);
        const auto film = camera.ProjectFilmForObserver(4.125, 7.5, .25f, .6f);
        ASSERT_TRUE(film);
        const auto expected = sirius::core::LaunchCameraRay(metric, mass * .7, film->ray);
        const auto actual = executor.Launch(metric, mass * .7, film->ray);
        ASSERT_TRUE(expected);
        ASSERT_TRUE(actual);
        const auto compare = [](double observed, double reference) {
            EXPECT_NEAR(observed, reference, 1e-10 * (1 + std::abs(reference)));
        };
        for (int axis = 0; axis < 4; ++axis) {
            compare(actual->position(axis), expected->position(axis));
            compare(actual->tangent(axis), expected->tangent(axis));
            compare(actual->observer.time(axis), expected->observer.time(axis));
            for (int basis = 0; basis < 3; ++basis)
                compare(actual->observer.spatial[basis](axis),
                        expected->observer.spatial[basis](axis));
            for (int column = 0; column < 4; ++column) {
                compare(actual->variations[column].displacement(axis),
                        expected->variations[column].displacement(axis));
                compare(actual->variations[column].derivative(axis),
                        expected->variations[column].derivative(axis));
            }
        }
    }
    std::vector<std::future<std::array<TraceResult, 2>>> workers;
    for (int worker = 0; worker < 4; ++worker) {
        workers.push_back(std::async(std::launch::async, [&, worker] {
            sirius::core::KerrSchildFamily metric(
                worker % 2 == 0 ? sirius::core::KerrSchildParams::Minkowski()
                                : sirius::core::KerrSchildParams::Kerr(1, .5));
            TracerConfig config;
            config.enable_disk = false;
            config.enable_ray_bundles = true;
            config.bundle_point_source = true;
            config.escape_radius = 12;
            config.max_steps = 1000;
            config.integrator.min_step = 1e-5f;
            config.integrator.max_step = 2;
            config.integrator.initial_step = 1;
            sirius::core::CameraRay ray;
            ray.origin(1) = 5;
            ray.origin(2) = std::numbers::pi / 2;
            ray.origin(3) = worker * .2;
            ray.direction(1) = 1;
            GeodesicTracer device_trace(&metric, config), reference(&metric, config);
            device_trace.SetStepExecutor(&executor);
            return std::array{device_trace.Trace(ray), reference.Trace(ray)};
        }));
    }
    for (std::size_t worker = 0; worker < workers.size(); ++worker) {
        SCOPED_TRACE(worker);
        const auto results = workers[worker].get();
        const auto& actual = results[0];
        const auto& expected = results[1];
        ASSERT_FALSE(actual.numerical_failure)
            << sirius::core::CoupledStepFailureName(actual.coupled_failure)
            << " termination=" << actual.integrator_termination
            << " attempts=" << actual.steps_taken;
        EXPECT_EQ(actual.outcome, TraceResult::Outcome::Escaped);
        EXPECT_EQ(actual.outcome, expected.outcome);
        EXPECT_GT(actual.central_stages, 0U);
        EXPECT_TRUE(actual.beam.valid);
        EXPECT_NEAR(actual.affine_length, expected.affine_length, 1e-5);
        EXPECT_NEAR(actual.redshift, expected.redshift, 1e-5);
        for (int axis = 0; axis < 4; ++axis) {
            EXPECT_NEAR(actual.final_position(axis), expected.final_position(axis), 1e-5);
            EXPECT_NEAR(actual.final_direction(axis), expected.final_direction(axis), 1e-5);
        }
    }
    ASSERT_FALSE(executor.Error());
    EXPECT_GT(executor.Statistics().interval_batches, 0U);
    EXPECT_GT(executor.Statistics().camera_batches, 0U);
    EXPECT_GT(executor.Statistics().reused_phases, 0U);
    EXPECT_EQ(executor.Statistics().initialized_phases, 4U);

    sirius::core::KerrSchildFamily flat(sirius::core::KerrSchildParams::Minkowski());
    sirius::core::Lightray ray{};
    ray.position(1) = 5;
    ray.velocity(0) = -1;
    ray.velocity(1) = 1;
    ray.step_size = 1;
    sirius::core::Rk45CoupledState coupled;
    coupled.length_scale = 1;
    coupled.frequency_scale = 1;
    coupled.tolerance = 1e-9;
    coupled.variations[0].derivative(2) = .001;
    coupled.variations[1].derivative(3) = .001;
    coupled.variations[2].displacement(2) = 1;
    coupled.variations[3].displacement(3) = 1;
    sirius::core::IntegratorConfig config;
    config.min_step = .01f;
    config.max_step = 2;
    const auto before = ray;
    const auto columns = coupled;
    sirius::core::Rk45CoupledComparison comparison;
    ASSERT_TRUE(executor.Step(ray, flat, config, coupled, comparison));
    const auto midpoint = comparison.midpoint;
    const auto reused = executor.Statistics().reused_phases;
    // Mimic a rejected localized event: restore the original public state and
    // retry half the interval. The future phase must not contaminate this retry.
    executor.RejectLastInterval();
    ray = before;
    coupled = columns;
    ray.step_size = .5f;
    ASSERT_TRUE(executor.Step(ray, flat, config, coupled, comparison));
    EXPECT_GT(executor.Statistics().reused_phases, reused);
    for (int axis = 0; axis < 4; ++axis) {
        EXPECT_EQ(ray.position(axis), midpoint.position(axis));
        EXPECT_EQ(ray.velocity(axis), midpoint.velocity(axis));
    }
    const auto phase_reuses = executor.Statistics().reused_phases;
    ray.proper_time = std::nextafter(ray.proper_time, std::numeric_limits<float>::infinity());
    ASSERT_TRUE(executor.Step(ray, flat, config, coupled, comparison));
    EXPECT_GT(executor.Statistics().reused_phases, phase_reuses);
    const auto accepted = ray;
    const auto submissions = executor.Statistics().interval_batches;
    cancelled = true;
    EXPECT_FALSE(executor.Step(ray, flat, config, coupled, comparison));
    EXPECT_EQ(executor.Statistics().interval_batches, submissions);
    for (int axis = 0; axis < 4; ++axis) {
        EXPECT_EQ(ray.position(axis), accepted.position(axis));
        EXPECT_EQ(ray.velocity(axis), accepted.velocity(axis));
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

TEST_F(RetainedComputeTest, Fp64ProductsPreserveIndependentScienceOrDeclineUnsupportedDevices) {
#ifdef SIRIUS_RETAINED_TESTS_AVAILABLE
    auto wide = RetainedCompute::Create(*device, 24, true);
    if (!device->Info().supports_fp64 || !device->Info().rounds_fp64_to_nearest) {
        ASSERT_FALSE(wide);
        EXPECT_NE(wide.error().detail().find("binary64"), std::string::npos);
        return;
    }
    ASSERT_TRUE(wide) << wide.error().Description();
    std::vector<RetainedRayCameraInput> cameras;
    for (const auto& fixture : sirius::test::retained_ray_camera::cases)
        cameras.push_back(std::bit_cast<RetainedRayCameraInput>(fixture.input));
    const auto launched = (*wide)->RayCamera(cameras);
    ASSERT_TRUE(launched) << launched.error().Description();
    for (std::size_t row = 0; row < cameras.size(); ++row)
        EXPECT_TRUE(CameraAgrees((*launched)[row], sirius::test::retained_ray_camera::cases[row]));
    std::vector<RetainedStepInput> steps;
    std::vector<RetainedEndpointInput> endpoints;
    for (const auto& fixture : sirius::test::retained_transport::cases) {
        steps.push_back(std::bit_cast<RetainedStepInput>(fixture.input));
        RetainedEndpointInput endpoint;
        std::copy_n(steps.back().values.begin(), 45, endpoint.values.begin());
        endpoints.push_back(endpoint);
    }
    const auto stepped = (*wide)->Step(steps);
    ASSERT_TRUE(stepped) << stepped.error().Description();
    for (std::size_t row = 0; row < steps.size(); ++row)
        EXPECT_TRUE(StepAgrees((*stepped)[row], sirius::test::retained_transport::cases[row]));
    const auto projected = (*wide)->Endpoint(endpoints);
    ASSERT_TRUE(projected) << projected.error().Description();
    std::vector<RetainedIntervalInput> intervals(steps.size());
    std::vector<RetainedInitializeInput> physical(steps.size());
    for (std::size_t row = 0; row < steps.size(); ++row) {
        ASSERT_TRUE((*projected)[row].valid);
        auto& input = intervals[row];
        std::copy_n(steps[row].values.begin(), 4, input.metric.begin());
        input.start = (*projected)[row];
        input.chart = steps[row].values[44].Center();
        input.interval = steps[row].values[45].Center();
        input.control.length_scale = 1;
        input.control.frequency_scale = 1;
        input.control.tolerance = 1e-4 / (4 * 30000);
        input.control.integrator.min_step = 1e-15f;
        input.control.integrator.max_step = 1;
        input.control.integrator.abs_tolerance = 1e-9f;
        input.control.integrator.rel_tolerance = 1e-9f;
        std::copy(input.metric.begin(), input.metric.end(), physical[row].values.begin());
        std::copy(input.start.physical.begin(), input.start.physical.end(),
                  physical[row].values.begin() + 4);
        physical[row].values[44] = steps[row].values[44];
    }
    const auto initialized = (*wide)->Initialize(physical);
    ASSERT_TRUE(initialized) << initialized.error().Description();
    for (const auto& row : *initialized) EXPECT_TRUE(row.valid);
    const auto narrow_intervals = AttemptRetainedIntervals(*compute, intervals);
    const auto wide_intervals = AttemptRetainedIntervals(**wide, intervals);
    ASSERT_TRUE(narrow_intervals) << narrow_intervals.error().Description();
    ASSERT_TRUE(wide_intervals) << wide_intervals.error().Description();
    for (std::size_t row = 0; row < intervals.size(); ++row) {
        SCOPED_TRACE(sirius::test::retained_transport::cases[row].name);
        ASSERT_TRUE((*wide_intervals)[row].admissible)
            << sirius::core::CoupledStepFailureName((*wide_intervals)[row].failure)
            << " error=" << (*wide_intervals)[row].error_ratio;
        ASSERT_TRUE((*narrow_intervals)[row].admissible);
        EXPECT_LT(
            RetainedPhysicalError((*wide_intervals)[row].full.physical,
                                  (*narrow_intervals)[row].full.physical, intervals[row].control),
            1e-4);
    }
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

    std::array<RetainedValue, 40> first;
    first.fill(RetainedValue::FromDouble(0));
    first[8] = {1, 0x1.000002p-35f, 0x1p-120f, 0, 1};
    ASSERT_TRUE(first[8].IsRepresented());
    auto second = first;
    second[8].tail = 0;
    RetainedIntervalControl control;
    control.length_scale = control.frequency_scale = 1;
    control.tolerance = 1e-40;
    const double expected = 0x1p-120 / (control.tolerance * (1 + first[8].Center()));
    EXPECT_DOUBLE_EQ(RetainedPhysicalError(first, second, control), expected);
    EXPECT_DOUBLE_EQ(RetainedPhysicalError(second, first, control), expected);
    EXPECT_EQ(RetainedPhysicalError(first, first, control), 0);
#else
    GTEST_SKIP() << "Retained compute build tools unavailable";
#endif
}
}  // namespace
