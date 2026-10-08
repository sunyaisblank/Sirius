#include "sirius/backend/retained_compute.h"
#include "sirius/core/twofold.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <functional>
#include <limits>
#include <memory>
#include <optional>
#include <span>
#include <string>
#include <vector>

#ifdef SIRIUS_HAS_RETAINED_COMPUTE
#include "sirius/backend/vulkan/vulkan_device.h"

#include "support/retained_transport/dopri_reference.h"
#include "support/retained_transport/reference_cases.h"
#endif

namespace {
using namespace sirius::backend;

class RetainedDopriTest : public ::testing::Test {
  protected:
    void SetUp() override {
#ifdef SIRIUS_HAS_RETAINED_COMPUTE
        const auto inventory = EnumerateVulkanDevices();
        ASSERT_TRUE(inventory) << inventory.error().Description();
        if (inventory->empty()) {
            GTEST_SKIP() << "No Vulkan device available";
        }
        const auto index = ResolveVulkanDeviceIndex(*inventory);
        ASSERT_TRUE(index) << index.error().Description();
        inventory_ = *inventory;
        device_index_ = *index;
#else
        GTEST_SKIP() << "Retained compute build tools unavailable";
#endif
    }

    void BothProductModes(const std::function<void(RetainedCompute&)>& check) {
#ifdef SIRIUS_HAS_RETAINED_COMPUTE
        std::optional<DeviceInfo> opened_identity;
        for (const bool wide : {false, true}) {
            SCOPED_TRACE(wide ? "fp64 products" : "default products");
            const auto current_inventory = EnumerateVulkanDevices();
            ASSERT_TRUE(current_inventory) << current_inventory.error().Description();
            ASSERT_EQ(*current_inventory, inventory_);
            auto opened = CreateVulkanDevice(device_index_);
            ASSERT_TRUE(opened) << opened.error().Description();
            device = std::move(*opened);
            if (opened_identity) {
                EXPECT_EQ(device->Info(), *opened_identity);
            } else {
                opened_identity = device->Info();
            }
            ASSERT_TRUE(device->SetBufferAllocationLimit(8 * 1024 * 1024));
            ASSERT_EQ(device->BufferAllocationBytes(), 0U);
            auto created = RetainedCompute::Create(*device, 24, wide);
            if (wide && (!device->Info().supports_fp64 || !device->Info().rounds_fp64_to_nearest)) {
                ASSERT_FALSE(created);
                EXPECT_NE(created.error().detail().find("binary64"), std::string::npos);
                EXPECT_EQ(device->BufferAllocationBytes(), 0U);
                device.reset();
                continue;
            }
            ASSERT_TRUE(created) << created.error().Description();
            const auto required = RetainedCompute::RequiredAllocationBytes(*device, 24);
            ASSERT_TRUE(required) << required.error().Description();
            EXPECT_EQ(device->BufferAllocationBytes(), *required);
            check(**created);
            created->reset();
            // Numeric kernel/buffer handles belong to the logical device.
            // Release its whole owner before opening the next product mode.
            device.reset();
            if (HasFatalFailure()) return;
        }
#else
        (void)check;
#endif
    }

    std::unique_ptr<ComputeDevice> device;
    std::vector<DeviceInfo> inventory_;
    std::size_t device_index_ = 0;
};

#ifdef SIRIUS_HAS_RETAINED_COMPUTE
namespace reference = sirius::test::retained_dopri;

sirius::core::Twofold Center(const RetainedValue& value) {
    return sirius::core::Twofold(value.high) + sirius::core::Twofold(value.low) +
           sirius::core::Twofold(value.tail);
}

void CheckValue(const RetainedValue& value, reference::Wide expected, double padding = 1e-28,
                double accuracy = 1e-10) {
    ASSERT_TRUE(value.IsRepresented());
    EXPECT_TRUE(std::isfinite(value.high));
    EXPECT_TRUE(std::isfinite(value.low));
    EXPECT_TRUE(std::isfinite(value.tail));
    EXPECT_TRUE(std::isfinite(value.radius));
    EXPECT_GE(value.radius, 0);
    const double difference =
        std::abs((Center(value) - sirius::core::Twofold(expected.high, expected.low)).Rounded());
    EXPECT_LE(difference,
              double(value.radius) + expected.gap + padding * (1 + std::abs(expected.high)));
    EXPECT_LE(difference, expected.gap + accuracy * (1 + std::abs(expected.high)));
}

auto Groups(const RetainedDopriPhaseOutput& output) {
    return std::array{&output.phase, &output.derivative, &output.a, &output.b, &output.c};
}

void CheckOutput(const RetainedDopriPhaseOutput& output, const reference::Case& fixture) {
    SCOPED_TRACE(fixture.name);
    ASSERT_TRUE(output.valid);
    const auto groups = Groups(output);
    for (std::size_t group = 0; group < groups.size(); ++group) {
        SCOPED_TRACE(group);
        for (std::size_t field = 0; field < 40; ++field) {
            SCOPED_TRACE(field);
            CheckValue((*groups[group])[field], fixture.reference[group * 40 + field]);
            if (group == 0) {
                const auto input_radius = std::bit_cast<float>(fixture.input[5 * field + 3]);
                EXPECT_GE(output.phase[field].radius, input_radius);
            }
        }
    }
    if (fixture.exact_phase) {
        for (std::size_t field = 0; field < 40; ++field) {
            const auto expected = fixture.reference[field];
            ASSERT_EQ(expected.gap, 0);
            // Sparse endpoint cancellation is an exact represented sum. An
            // absolute tolerance would otherwise hide a lost 2^-120 tail.
            EXPECT_EQ(
                (Center(output.phase[field]) - sirius::core::Twofold(expected.high, expected.low))
                    .Rounded(),
                0);
        }
    }
}

void SameWords(const RetainedDopriPhaseOutput& first, const RetainedDopriPhaseOutput& second) {
    ASSERT_EQ(first.valid, second.valid);
    const auto a = Groups(first), b = Groups(second);
    for (std::size_t group = 0; group < a.size(); ++group) {
        for (std::size_t field = 0; field < 40; ++field) {
            SCOPED_TRACE(group * 40 + field);
            EXPECT_EQ((std::bit_cast<std::array<std::uint32_t, 5>>((*a[group])[field])),
                      (std::bit_cast<std::array<std::uint32_t, 5>>((*b[group])[field])));
        }
    }
}

void Refused(const RetainedDopriPhaseOutput& output) {
    EXPECT_FALSE(output.valid);
    for (const auto* group : Groups(output)) {
        for (const auto& value : *group) {
            EXPECT_EQ((std::bit_cast<std::array<std::uint32_t, 5>>(value)),
                      (std::array<std::uint32_t, 5>{}));
        }
    }
}
#endif

TEST_F(RetainedDopriTest, IndependentQuarticPreservesCompletePhase) {
#ifdef SIRIUS_HAS_RETAINED_COMPUTE
    BothProductModes([&](RetainedCompute& compute) {
        const auto allocation = device->BufferAllocationBytes();
        std::vector<RetainedDopriPhaseOutput> original;
        for (std::size_t begin = 0; begin < reference::cases.size(); begin += compute.Capacity()) {
            const auto count = std::min(compute.Capacity(), reference::cases.size() - begin);
            std::vector<RetainedDopriPhaseInput> inputs;
            for (std::size_t i = 0; i < count; ++i) {
                inputs.push_back(
                    std::bit_cast<RetainedDopriPhaseInput>(reference::cases[begin + i].input));
            }
            const auto sampled = compute.DopriPhase(inputs);
            ASSERT_TRUE(sampled) << sampled.error().Description();
            ASSERT_EQ(sampled->size(), count);
            for (std::size_t i = 0; i < count; ++i) {
                CheckOutput((*sampled)[i], reference::cases[begin + i]);
                original.push_back((*sampled)[i]);
            }
            EXPECT_EQ(device->BufferAllocationBytes(), allocation);
        }

        const auto named = [](const char* name) {
            const auto found =
                std::find_if(reference::cases.begin(), reference::cases.end(),
                             [name](const auto& row) { return std::string(row.name) == name; });
            return static_cast<std::size_t>(found - reference::cases.begin());
        };
        const auto complete = named("noncanonical-1/2");
        ASSERT_LT(complete, original.size());
        for (const auto* name : {"low-limb-removed", "tail-limb-removed"}) {
            const auto removed = named(name);
            ASSERT_LT(removed, original.size());
            EXPECT_GT(
                std::abs((Center(original[complete].phase[0]) - Center(original[removed].phase[0]))
                             .Rounded()),
                1e-19);
        }

        // Reorder complete rows and shrink the active prefix. Retained words,
        // including all arithmetic radii, must remain independent of row slots.
        const std::array<std::size_t, 4> indices{reference::cases.size() - 1, 3, 0, 9};
        std::array<RetainedDopriPhaseInput, 4> reordered;
        for (std::size_t row = 0; row < indices.size(); ++row) {
            reordered[row] =
                std::bit_cast<RetainedDopriPhaseInput>(reference::cases[indices[row]].input);
        }
        const auto repeated = compute.DopriPhase(reordered);
        ASSERT_TRUE(repeated) << repeated.error().Description();
        ASSERT_EQ(repeated->size(), reordered.size());
        for (std::size_t row = 0; row < indices.size(); ++row) {
            CheckOutput((*repeated)[row], reference::cases[indices[row]]);
            SameWords((*repeated)[row], original[indices[row]]);
        }
        EXPECT_EQ(device->BufferAllocationBytes(), allocation);
        const auto stage = static_cast<std::size_t>(RetainedCompute::KernelStage::kDopriPhase);
        const auto& statistics = compute.Statistics();
        EXPECT_EQ(statistics[stage].submissions,
                  (reference::cases.size() + compute.Capacity() - 1) / compute.Capacity() + 1);
        for (std::size_t i = 0; i < statistics.size(); ++i) {
            if (i != stage) {
                EXPECT_EQ(statistics[i].submissions, 0U);
            }
        }
    });
#endif
}

TEST_F(RetainedDopriTest, MalformedRowsRefuseAndRecover) {
#ifdef SIRIUS_HAS_RETAINED_COMPUTE
    BothProductModes([&](RetainedCompute& compute) {
        // Strictly interior fraction and positive h have zero-width encoded
        // enclosures, making the boundary changes below independently known.
        const auto fixture =
            std::find_if(reference::cases.begin(), reference::cases.end(),
                         [](const auto& row) { return std::string(row.name) == "arbitrary-1/2"; });
        ASSERT_NE(fixture, reference::cases.end());
        const auto good = std::bit_cast<RetainedDopriPhaseInput>(fixture->input);
        const auto allocation = device->BufferAllocationBytes();
        const auto first = compute.DopriPhase(std::span(&good, 1));
        ASSERT_TRUE(first) << first.error().Description();
        CheckOutput(first->front(), *fixture);
        EXPECT_FALSE(compute.DopriPhase({}));
        std::vector<RetainedDopriPhaseInput> oversized(compute.Capacity() + 1, good);
        EXPECT_FALSE(compute.DopriPhase(oversized));

        std::array<RetainedDopriPhaseInput, 24> inputs;
        inputs.fill(good);
        inputs[0].values[0].high = std::numeric_limits<float>::quiet_NaN();
        inputs[1].values[40].low = std::numeric_limits<float>::infinity();
        inputs[2].values[119].tail = -std::numeric_limits<float>::infinity();
        inputs[3].values[120].valid = 0;  // k2 has zero weight, but must be validated.
        inputs[4].values[159].valid = 2;
        inputs[5].values[200].radius = -1;
        inputs[6].values[239].radius = std::numeric_limits<float>::quiet_NaN();
        inputs[7].values[319].radius = std::numeric_limits<float>::infinity();
        inputs[8].values[359].tail = std::numeric_limits<float>::quiet_NaN();
        inputs[9].values[0] = RetainedValue::FromDouble(0x1p121);
        inputs[10].values[360] = RetainedValue::FromDouble(0);
        inputs[11].values[360] = RetainedValue::FromDouble(-1);
        inputs[12].values[360] = {1, -1, -0x1p-60f, 0, 1};
        inputs[13].values[360].radius = inputs[13].values[360].high;
        inputs[14].values[361] = RetainedValue::FromDouble(-0x1p-30);
        inputs[15].values[361] = {1, 0x1p-30f, 0, 0, 1};
        inputs[16].values[361] = {1, -1, -0x1p-60f, 0, 1};
        inputs[17].values[361] = {1, 0, 0x1p-60f, 0, 1};
        inputs[18].values[361].radius = .5f;  // Enclosure touches both endpoints.
        inputs[19].values[361] = {0, 0, 0, 0x1p-80f, 1};
        inputs[20].values[361] = {1, 0, 0, 0x1p-80f, 1};
        // All operands are individually admitted; the coefficient arithmetic
        // exceeds the retained public bound and must refuse without publication.
        inputs[21].values[360] = RetainedValue::FromDouble(0x1p100);
        inputs[21].values[80] = RetainedValue::FromDouble(0x1p100);
        // Final two rows remain good, detecting cross-row refusal contamination.
        const auto refused = compute.DopriPhase(inputs);
        ASSERT_TRUE(refused) << refused.error().Description();
        ASSERT_EQ(refused->size(), inputs.size());
        for (std::size_t row = 0; row < 22; ++row) {
            SCOPED_TRACE(row);
            Refused((*refused)[row]);
        }
        for (std::size_t row = 22; row < inputs.size(); ++row) {
            CheckOutput((*refused)[row], *fixture);
            SameWords((*refused)[row], first->front());
        }
        // Exact endpoint ownership does not bypass input or coefficient
        // validation. Even the zero-weight stage must remain represented.
        std::array<RetainedDopriPhaseInput, 5> endpoint_refusals{inputs[21], inputs[21], inputs[3],
                                                                 inputs[3], good};
        for (std::size_t row = 0; row < 4; ++row) {
            endpoint_refusals[row].values[361] =
                RetainedValue::FromDouble(static_cast<double>(row & 1));
        }
        const auto endpoint_results = compute.DopriPhase(endpoint_refusals);
        ASSERT_TRUE(endpoint_results) << endpoint_results.error().Description();
        ASSERT_EQ(endpoint_results->size(), endpoint_refusals.size());
        for (std::size_t row = 0; row < 4; ++row) {
            SCOPED_TRACE(row);
            Refused((*endpoint_results)[row]);
        }
        CheckOutput(endpoint_results->back(), *fixture);
        SameWords(endpoint_results->back(), first->front());
        const auto recovered = compute.DopriPhase(std::span(&good, 1));
        ASSERT_TRUE(recovered) << recovered.error().Description();
        ASSERT_EQ(recovered->size(), 1U);
        CheckOutput(recovered->front(), *fixture);
        SameWords(recovered->front(), first->front());
        EXPECT_EQ(device->BufferAllocationBytes(), allocation);
    });
#endif
}

TEST_F(RetainedDopriTest, CurvedTransportConnectsSamplerAndProjection) {
#ifdef SIRIUS_HAS_RETAINED_COMPUTE
    BothProductModes([&](RetainedCompute& compute) {
        const auto allocation = device->BufferAllocationBytes();
        const auto& fixture = reference::connection;
        const auto input = std::bit_cast<RetainedStepInput>(fixture.input);
        const auto stepped = compute.Step(std::span(&input, 1));
        ASSERT_TRUE(stepped) << stepped.error().Description();
        ASSERT_EQ(stepped->size(), 1U);
        ASSERT_TRUE(stepped->front().valid);
        ASSERT_TRUE(stepped->front().rhs_valid);
        EXPECT_EQ(stepped->front().stages, 7U);
        RetainedDopriPhaseInput dense;
        std::copy_n(input.values.begin() + 4, 40, dense.values.begin());
        std::copy(stepped->front().increment.begin(), stepped->front().increment.end(),
                  dense.values.begin() + 40);
        for (std::size_t stage = 0; stage < 7; ++stage) {
            for (std::size_t field = 0; field < 40; ++field) {
                SCOPED_TRACE(stage * 40 + field);
                const auto value = stepped->front().rhs[stage][field];
                CheckValue(value, fixture.reference[stage * 40 + field], 1e-29, 1e-11);
                dense.values[80 + stage * 40 + field] = value;
            }
        }
        dense.values[360] = input.values[45];
        dense.values[361] = RetainedValue::FromDouble(3. / 8);
        const auto sampled = compute.DopriPhase(std::span(&dense, 1));
        ASSERT_TRUE(sampled) << sampled.error().Description();
        ASSERT_EQ(sampled->size(), 1U);
        ASSERT_TRUE(sampled->front().valid);
        const auto groups = Groups(sampled->front());
        for (std::size_t group = 0; group < groups.size(); ++group) {
            for (std::size_t field = 0; field < 40; ++field) {
                SCOPED_TRACE(group * 40 + field);
                CheckValue((*groups[group])[field], fixture.reference[280 + group * 40 + field]);
            }
        }
        RetainedEndpointInput endpoint;
        std::copy_n(input.values.begin(), 4, endpoint.values.begin());
        std::copy(sampled->front().phase.begin(), sampled->front().phase.end(),
                  endpoint.values.begin() + 4);
        endpoint.values[44] = input.values[44];
        const auto projected = compute.Endpoint(std::span(&endpoint, 1));
        ASSERT_TRUE(projected) << projected.error().Description();
        ASSERT_EQ(projected->size(), 1U);
        ASSERT_TRUE(projected->front().valid);
        EXPECT_EQ(projected->front().component, fixture.component);
        for (std::size_t field = 0; field < 40; ++field) {
            SCOPED_TRACE(field);
            CheckValue(projected->front().phase[field], fixture.reference[480 + field], 1e-29,
                       1e-11);
            CheckValue(projected->front().physical[field], fixture.reference[520 + field], 1e-29,
                       1e-11);
        }
        // The original analytic flat shortcut is valid without seven stored
        // RHS records. Its marker must not fabricate a sampler-ready capsule.
        const auto flat =
            std::bit_cast<RetainedStepInput>(sirius::test::retained_transport::cases.back().input);
        const auto flat_step = compute.Step(std::span(&flat, 1));
        ASSERT_TRUE(flat_step) << flat_step.error().Description();
        ASSERT_TRUE(flat_step->front().valid);
        EXPECT_FALSE(flat_step->front().rhs_valid);
        for (const auto& stage : flat_step->front().rhs) {
            for (const auto& value : stage) {
                EXPECT_EQ(value.valid, 0U);
            }
        }
        EXPECT_EQ(device->BufferAllocationBytes(), allocation);
    });
#endif
}
}  // namespace
