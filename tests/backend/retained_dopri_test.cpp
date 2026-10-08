#include "sirius/backend/retained_compute.h"
#include "sirius/backend/retained_integrator.h"
#include "sirius/backend/retained_trace_executor.h"
#include "sirius/core/metrics/outgoing_kerr_schild.h"
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
#include <utility>
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

void SameEndpoint(const RetainedEndpointOutput& actual, const RetainedEndpointOutput& expected) {
    ASSERT_EQ(actual.valid, expected.valid);
    EXPECT_EQ(actual.component, expected.component);
    for (std::size_t field = 0; field < 40; ++field) {
        SCOPED_TRACE(field);
        for (const auto& pair : {std::pair{&actual.phase[field], &expected.phase[field]},
                                 std::pair{&actual.physical[field], &expected.physical[field]}}) {
            EXPECT_EQ((std::bit_cast<std::array<std::uint32_t, 5>>(*pair.first)),
                      (std::bit_cast<std::array<std::uint32_t, 5>>(*pair.second)));
        }
    }
}

double IndependentPhysicalRatio(const std::array<RetainedValue, 40>& first,
                                const std::array<RetainedValue, 40>& second,
                                const RetainedIntervalControl& control) {
    double sum = 0, columns = 0;
    for (std::size_t field = 0; field < 40; ++field) {
        const auto a = Center(first[field]), b = Center(second[field]);
        const double unit = field % 8 < 4 ? control.length_scale : control.frequency_scale;
        const double magnitude = std::max(std::abs(a.Rounded()), std::abs(b.Rounded()));
        const double allowance =
            field < 8
                ? control.integrator.abs_tolerance * unit +
                      control.integrator.rel_tolerance * magnitude
                : control.tolerance * (unit * control.column_scale[(field - 8) / 8] + magnitude);
        const double ratio = std::abs((a - b).Rounded()) / allowance;
        if (field < 8) {
            sum += ratio * ratio;
        } else {
            columns = std::max(columns, ratio);
        }
    }
    return std::max(std::sqrt(sum / 8), columns);
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

TEST_F(RetainedDopriTest, ConnectedIntervalsPreserveTrialsAndIndependentBudgets) {
#ifdef SIRIUS_HAS_RETAINED_COMPUTE
    BothProductModes([&](RetainedCompute& compute) {
        const auto allocation = device->BufferAllocationBytes();
        const auto original = std::bit_cast<RetainedStepInput>(reference::connection.input);
        RetainedEndpointInput launch;
        std::copy_n(original.values.begin(), 45, launch.values.begin());
        const auto initial = compute.Endpoint(std::span(&launch, 1));
        ASSERT_TRUE(initial) << initial.error().Description();
        ASSERT_TRUE(initial->front().valid);
        RetainedIntervalInput input;
        std::copy_n(original.values.begin(), 4, input.metric.begin());
        input.start = initial->front();
        input.chart = original.values[44].Center();
        input.interval = original.values[45].Center();
        // These are the unchanged original coupled-fixture component budgets.
        input.control.length_scale = input.control.frequency_scale = 1;
        input.control.tolerance = 1e-4 / (4 * 30000);
        input.control.integrator.min_step = static_cast<float>(input.interval);
        input.control.integrator.max_step = 1;
        input.control.integrator.abs_tolerance = input.control.integrator.rel_tolerance = 1e-9f;
        const auto same_values = [](const auto& actual, const auto& expected) {
            ASSERT_EQ(actual.size(), expected.size());
            for (std::size_t field = 0; field < actual.size(); ++field) {
                SCOPED_TRACE(field);
                EXPECT_EQ((std::bit_cast<std::array<std::uint32_t, 5>>(actual[field])),
                          (std::bit_cast<std::array<std::uint32_t, 5>>(expected[field])));
            }
        };
        const auto step_control = [](const sirius::core::IntegratorConfig& control) {
            return std::array{control.abs_tolerance, control.rel_tolerance,  control.min_step,
                              control.max_step,      control.initial_step,   control.safety_factor,
                              control.step_grow_max, control.step_shrink_min};
        };
        const auto same_attempt = [&](const RetainedIntervalOutput& actual,
                                      const RetainedIntervalOutput& expected) {
            EXPECT_EQ(actual.admissible, expected.admissible);
            EXPECT_EQ(actual.failure, expected.failure);
            EXPECT_EQ(actual.attempted_stages, expected.attempted_stages);
            EXPECT_EQ(actual.error_ratio, expected.error_ratio);
            EXPECT_EQ(actual.embedded_projected_error_ratio,
                      expected.embedded_projected_error_ratio);
            for (std::size_t check = 0; check < actual.error_checks.size(); ++check) {
                SCOPED_TRACE(check);
                const auto& a = actual.error_checks[check];
                const auto& b = expected.error_checks[check];
                EXPECT_EQ(a.ratio, b.ratio);
                EXPECT_EQ(a.field, b.field);
                EXPECT_EQ(a.observed, b.observed);
                EXPECT_EQ(a.evaluated, b.evaluated);
            }
            SameEndpoint(actual.full, expected.full);
            SameEndpoint(actual.lower, expected.lower);
            SameEndpoint(actual.midpoint, expected.midpoint);
            SameEndpoint(actual.refined, expected.refined);
            same_values(actual.full_increment, expected.full_increment);
            same_values(actual.lower_increment, expected.lower_increment);
            same_values(actual.midpoint_increment, expected.midpoint_increment);
            same_values(actual.refined_increment, expected.refined_increment);
            ASSERT_EQ(bool(actual.dopri), bool(expected.dopri));
            if (!actual.admissible) {
                EXPECT_FALSE(actual.dopri);
                for (const auto* endpoint :
                     {&actual.full, &actual.lower, &actual.midpoint, &actual.refined}) {
                    SameEndpoint(*endpoint, RetainedEndpointOutput{});
                }
                const std::array<RetainedValue, 20> erased{};
                same_values(actual.full_increment, erased);
                same_values(actual.lower_increment, erased);
                same_values(actual.midpoint_increment, erased);
                same_values(actual.refined_increment, erased);
            }
            if (!actual.dopri) return;
            const auto& a = *actual.dopri;
            const auto& b = *expected.dopri;
            same_values(a.metric, b.metric);
            EXPECT_EQ(a.chart, b.chart);
            EXPECT_EQ(a.control.length_scale, b.control.length_scale);
            EXPECT_EQ(a.control.frequency_scale, b.control.frequency_scale);
            EXPECT_EQ(a.control.tolerance, b.control.tolerance);
            EXPECT_EQ(a.control.column_scale, b.control.column_scale);
            EXPECT_EQ(step_control(a.control.integrator), step_control(b.control.integrator));
            for (std::size_t trial = 0; trial < sirius::core::kCoupledTrialCount; ++trial) {
                SCOPED_TRACE(trial);
                same_values(a.packets[trial].values, b.packets[trial].values);
                SameEndpoint(a.starts[trial], b.starts[trial]);
                SameEndpoint(a.endpoints[trial], b.endpoints[trial]);
                const auto& p = a.positions[trial];
                const auto& q = b.positions[trial];
                EXPECT_EQ(p.interval, q.interval);
                EXPECT_EQ(p.parameter_limit, q.parameter_limit);
                for (int axis = 0; axis < 4; ++axis) {
                    EXPECT_EQ(p.origin(axis), q.origin(axis));
                    EXPECT_EQ(p.increment(axis), q.increment(axis));
                    EXPECT_EQ(p.a(axis), q.a(axis));
                    EXPECT_EQ(p.b(axis), q.b(axis));
                    EXPECT_EQ(p.c(axis), q.c(axis));
                }
            }
        };
        constexpr auto phase_stage =
            static_cast<std::size_t>(RetainedCompute::KernelStage::kDopriPhase);
        const auto serial_before = compute.Statistics()[phase_stage].submissions;
        const auto attempted = AttemptRetainedDopriIntervals(compute, std::span(&input, 1), 1);
        ASSERT_TRUE(attempted) << attempted.error().Description();
        ASSERT_EQ(attempted->size(), 1U);
        const auto serial_after = compute.Statistics()[phase_stage].submissions;
        EXPECT_EQ(serial_after - serial_before, 4U);
        const auto& accepted = attempted->front();
        ASSERT_TRUE(accepted.admissible) << sirius::core::CoupledStepFailureName(accepted.failure)
                                         << " " << accepted.error_ratio;
        EXPECT_EQ(accepted.failure, sirius::core::CoupledStepFailure::None);
        EXPECT_EQ(accepted.attempted_stages, 21U);
        EXPECT_LE(accepted.error_ratio, 1);
        EXPECT_GE(accepted.embedded_projected_error_ratio, 0);
        EXPECT_LE(accepted.embedded_projected_error_ratio, accepted.error_ratio);
        ASSERT_TRUE(accepted.dopri);
        const auto& capsule = *accepted.dopri;
        constexpr auto full = static_cast<std::size_t>(sirius::core::CoupledTrial::Full);
        constexpr auto lower = static_cast<std::size_t>(sirius::core::CoupledTrial::Lower);
        constexpr auto first_half = static_cast<std::size_t>(sirius::core::CoupledTrial::FirstHalf);
        constexpr auto second_half =
            static_cast<std::size_t>(sirius::core::CoupledTrial::SecondHalf);

        // Count only the connected attempt. The standalone independent trials
        // below deliberately add their own phase submissions afterwards.
        const auto packed_before = compute.Statistics()[phase_stage].submissions;
        const auto packed = AttemptRetainedDopriIntervals(compute, std::span(&input, 1));
        ASSERT_TRUE(packed) << packed.error().Description();
        ASSERT_EQ(packed->size(), 1U);
        const auto packed_after = compute.Statistics()[phase_stage].submissions;
        EXPECT_EQ(packed_after - packed_before, 1U);
        same_attempt(packed->front(), accepted);
        EXPECT_EQ(device->BufferAllocationBytes(), allocation);

        // Two different intervals cross a trial boundary inside a three-row
        // chunk: Full0, Full1, Lower0, then Lower1, FirstHalf0, FirstHalf1.
        auto shorter = input;
        shorter.interval *= .5;
        shorter.control.integrator.min_step = static_cast<float>(shorter.interval);
        const auto shorter_serial =
            AttemptRetainedDopriIntervals(compute, std::span(&shorter, 1), 1);
        ASSERT_TRUE(shorter_serial) << shorter_serial.error().Description();
        ASSERT_EQ(shorter_serial->size(), 1U);
        ASSERT_TRUE(shorter_serial->front().admissible)
            << sirius::core::CoupledStepFailureName(shorter_serial->front().failure) << " "
            << shorter_serial->front().error_ratio;
        ASSERT_TRUE(shorter_serial->front().dopri);
        EXPECT_NE(shorter_serial->front().full.physical[0].Center(),
                  accepted.full.physical[0].Center());
        const std::array distinct{input, shorter};
        const auto distinct_before = compute.Statistics()[phase_stage].submissions;
        (void)compute.TakeSubmissionFeedback();
        const auto distinct_packed = AttemptRetainedDopriIntervals(compute, distinct, 3);
        ASSERT_TRUE(distinct_packed) << distinct_packed.error().Description();
        ASSERT_EQ(distinct_packed->size(), distinct.size());
        const auto distinct_after = compute.Statistics()[phase_stage].submissions;
        EXPECT_EQ(distinct_after - distinct_before, 3U);
        EXPECT_EQ(compute.TakeSubmissionFeedback().maximum_rows, 3U);
        for (const auto& result : *distinct_packed) {
            ASSERT_TRUE(result.admissible) << sirius::core::CoupledStepFailureName(result.failure)
                                           << " " << result.error_ratio;
            ASSERT_TRUE(result.dopri);
        }
        same_attempt((*distinct_packed)[0], accepted);
        same_attempt((*distinct_packed)[1], shorter_serial->front());
        EXPECT_EQ(device->BufferAllocationBytes(), allocation);

        // Form three independent original Step/Endpoint trials, without using
        // capsule starts or increments to supply their initial values.
        std::array<RetainedStepOutput, 3> steps;
        std::array<RetainedEndpointOutput, 3> endpoints;
        RetainedEndpointOutput low_endpoint;
        for (std::size_t part = 0; part < steps.size(); ++part) {
            const auto& start = part == 2 ? endpoints[1] : input.start;
            RetainedStepInput packet;
            std::copy(input.metric.begin(), input.metric.end(), packet.values.begin());
            std::copy(start.phase.begin(), start.phase.end(), packet.values.begin() + 4);
            packet.values[44] = RetainedValue::FromDouble(input.chart);
            packet.values[45] = RetainedValue::FromDouble(input.interval / (part == 0 ? 1 : 2));
            const auto evaluated = compute.Step(std::span(&packet, 1));
            ASSERT_TRUE(evaluated) << evaluated.error().Description();
            ASSERT_TRUE(evaluated->front().valid);
            ASSERT_TRUE(evaluated->front().rhs_valid);
            EXPECT_EQ(evaluated->front().stages, 7U);
            steps[part] = evaluated->front();
            std::array<RetainedEndpointInput, 2> projection;
            for (std::size_t row = 0; row < projection.size(); ++row) {
                std::copy(input.metric.begin(), input.metric.end(), projection[row].values.begin());
                const auto& phase = row == 0 ? steps[part].fifth : steps[part].fourth;
                std::copy(phase.begin(), phase.end(), projection[row].values.begin() + 4);
                projection[row].values[44] = RetainedValue::FromDouble(input.chart);
            }
            const auto projected = compute.Endpoint(projection);
            ASSERT_TRUE(projected) << projected.error().Description();
            ASSERT_TRUE(projected->front().valid);
            ASSERT_TRUE(projected->back().valid);
            endpoints[part] = projected->front();
            if (part == 0) low_endpoint = projected->back();
        }
        SameEndpoint(accepted.full, endpoints[0]);
        SameEndpoint(accepted.lower, low_endpoint);
        SameEndpoint(accepted.midpoint, endpoints[1]);
        SameEndpoint(accepted.refined, endpoints[2]);
        const std::array trial_endpoints{&accepted.full, &accepted.lower, &accepted.midpoint,
                                         &accepted.refined};
        const std::array<std::size_t, 4> step_index{0, 0, 1, 2};
        auto endpoint_packets = capsule.packets;
        for (auto& packet : endpoint_packets) packet.values[361] = RetainedValue::FromDouble(1);
        const auto sampled_endpoints = compute.DopriPhase(endpoint_packets);
        ASSERT_TRUE(sampled_endpoints) << sampled_endpoints.error().Description();
        ASSERT_EQ(sampled_endpoints->size(), sirius::core::kCoupledTrialCount);
        const auto& coefficients = *sampled_endpoints;
        for (const auto trial : {full, lower, first_half, second_half}) {
            SCOPED_TRACE(trial);
            const auto& start = trial == second_half ? endpoints[1] : input.start;
            SameEndpoint(capsule.starts[trial], start);
            SameEndpoint(capsule.endpoints[trial], *trial_endpoints[trial]);
            const auto& packet = capsule.packets[trial];
            ASSERT_TRUE(coefficients[trial].valid);
            for (std::size_t field = 0; field < 40; ++field) {
                SCOPED_TRACE(field);
                EXPECT_EQ((std::bit_cast<std::array<std::uint32_t, 5>>(packet.values[field])),
                          (std::bit_cast<std::array<std::uint32_t, 5>>(start.phase[field])));
                const auto& step = steps[step_index[trial]];
                const auto& end_phase = trial == lower ? step.fourth : step.fifth;
                EXPECT_LE(
                    std::abs((Center(coefficients[trial].phase[field]) - Center(end_phase[field]))
                                 .Rounded()),
                    double(coefficients[trial].phase[field].radius) + end_phase[field].radius +
                        1e-29 * (1 + std::abs(Center(end_phase[field]).Rounded())));
                if (trial == lower) {
                    const auto expected = Center(step.increment[field]) - Center(step.error[field]);
                    EXPECT_LE(std::abs((Center(packet.values[40 + field]) - expected).Rounded()),
                              packet.values[40 + field].radius +
                                  1e-29 * (1 + std::abs(expected.Rounded())));
                } else {
                    EXPECT_EQ(
                        (std::bit_cast<std::array<std::uint32_t, 5>>(packet.values[40 + field])),
                        (std::bit_cast<std::array<std::uint32_t, 5>>(step.increment[field])));
                }
                for (std::size_t stage = 0; stage < 7; ++stage) {
                    EXPECT_EQ(
                        (std::bit_cast<std::array<std::uint32_t, 5>>(
                            packet.values[80 + stage * 40 + field])),
                        (std::bit_cast<std::array<std::uint32_t, 5>>(step.rhs[stage][field])));
                }
            }
            const double h = input.interval / (trial == full || trial == lower ? 1 : 2);
            EXPECT_EQ(Center(packet.values[360]).Rounded(), h);
            EXPECT_EQ(Center(packet.values[361]).Rounded(), 1);
            const auto& position = capsule.positions[trial];
            EXPECT_TRUE(position.IsFinite());
            EXPECT_EQ(position.interval, h);
            EXPECT_EQ(position.parameter_limit, 1);
            for (int axis = 0; axis < 4; ++axis) {
                EXPECT_EQ(position.origin(axis), Center(packet.values[axis]).Rounded());
                EXPECT_EQ(position.increment(axis), Center(packet.values[40 + axis]).Rounded());
                EXPECT_EQ(position.a(axis), Center(coefficients[trial].a[axis]).Rounded());
                EXPECT_EQ(position.b(axis), Center(coefficients[trial].b[axis]).Rounded());
                EXPECT_EQ(position.c(axis), Center(coefficients[trial].c[axis]).Rounded());
            }
        }
        EXPECT_EQ(capsule.chart, input.chart);
        EXPECT_EQ(capsule.control.length_scale, input.control.length_scale);
        EXPECT_EQ(capsule.control.frequency_scale, input.control.frequency_scale);
        EXPECT_EQ(capsule.control.tolerance, input.control.tolerance);
        EXPECT_EQ(capsule.control.column_scale, input.control.column_scale);
        EXPECT_EQ(step_control(capsule.control.integrator), step_control(input.control.integrator));
        for (std::size_t parameter = 0; parameter < 4; ++parameter) {
            EXPECT_EQ((std::bit_cast<std::array<std::uint32_t, 5>>(capsule.metric[parameter])),
                      (std::bit_cast<std::array<std::uint32_t, 5>>(input.metric[parameter])));
        }

        auto full_packet = capsule.packets[full];
        full_packet.values[361] = RetainedValue::FromDouble(.5);
        const auto dense_midpoint = compute.DopriPhase(std::span(&full_packet, 1));
        ASSERT_TRUE(dense_midpoint) << dense_midpoint.error().Description();
        ASSERT_TRUE(dense_midpoint->front().valid);
        RetainedEndpointInput midpoint_projection = launch;
        std::copy(dense_midpoint->front().phase.begin(), dense_midpoint->front().phase.end(),
                  midpoint_projection.values.begin() + 4);
        const auto midpoint = compute.Endpoint(std::span(&midpoint_projection, 1));
        ASSERT_TRUE(midpoint) << midpoint.error().Description();
        ASSERT_TRUE(midpoint->front().valid);
        const double midpoint_ratio = IndependentPhysicalRatio(
            midpoint->front().physical, endpoints[1].physical, input.control);
        const double refined_ratio =
            IndependentPhysicalRatio(endpoints[0].physical, endpoints[2].physical, input.control);
        EXPECT_TRUE(std::isfinite(midpoint_ratio));
        EXPECT_TRUE(std::isfinite(refined_ratio));
        EXPECT_LE(midpoint_ratio, 1);
        EXPECT_LE(refined_ratio, 1);
        EXPECT_LE(std::abs(accepted.error_checks[2].ratio - midpoint_ratio),
                  1e-12 * (1 + midpoint_ratio));
        EXPECT_LE(std::abs(accepted.error_checks[3].ratio - refined_ratio),
                  1e-12 * (1 + refined_ratio));
        for (const auto& observation : accepted.error_checks) {
            EXPECT_TRUE(observation.observed);
            EXPECT_TRUE(observation.evaluated);
            EXPECT_TRUE(std::isfinite(observation.ratio));
            EXPECT_LE(observation.ratio, accepted.error_ratio);
        }

        // The new full packet still meets the independently frozen Kerr RHS
        // and interior polynomial reference, rather than a coordinator oracle.
        for (std::size_t field = 0; field < 280; ++field) {
            CheckValue(capsule.packets[full].values[80 + field],
                       reference::connection.reference[field], 1e-29, 1e-11);
        }
        full_packet.values[361] = RetainedValue::FromDouble(3. / 8);
        const auto interior = compute.DopriPhase(std::span(&full_packet, 1));
        ASSERT_TRUE(interior) << interior.error().Description();
        ASSERT_TRUE(interior->front().valid);
        const auto groups = Groups(interior->front());
        for (std::size_t group = 0; group < groups.size(); ++group) {
            for (std::size_t field = 0; field < 40; ++field) {
                CheckValue((*groups[group])[field],
                           reference::connection.reference[280 + group * 40 + field]);
            }
        }

        auto malformed = input;
        malformed.start.phase.back().valid = 0;
        const auto flat_step =
            std::bit_cast<RetainedStepInput>(sirius::test::retained_transport::cases.back().input);
        RetainedEndpointInput flat_launch;
        std::copy_n(flat_step.values.begin(), 45, flat_launch.values.begin());
        const auto flat_start = compute.Endpoint(std::span(&flat_launch, 1));
        ASSERT_TRUE(flat_start) << flat_start.error().Description();
        ASSERT_TRUE(flat_start->front().valid);
        auto flat_input = input;
        std::copy_n(flat_step.values.begin(), 4, flat_input.metric.begin());
        flat_input.start = flat_start->front();
        flat_input.chart = flat_step.values[44].Center();
        flat_input.interval = flat_step.values[45].Center();
        flat_input.control.integrator.min_step = static_cast<float>(flat_input.interval);
        const std::array mixed{malformed, input, flat_input};
        const auto mixed_result = AttemptRetainedDopriIntervals(compute, mixed, 3);
        ASSERT_TRUE(mixed_result) << mixed_result.error().Description();
        EXPECT_FALSE(mixed_result->front().admissible);
        EXPECT_FALSE(mixed_result->front().dopri);
        EXPECT_EQ(mixed_result->front().failure, sirius::core::CoupledStepFailure::InvalidState);
        EXPECT_EQ(mixed_result->front().attempted_stages, 0U);
        ASSERT_TRUE((*mixed_result)[1].admissible);
        ASSERT_TRUE((*mixed_result)[1].dopri);
        SameEndpoint((*mixed_result)[1].full, accepted.full);
        ASSERT_TRUE(mixed_result->back().admissible);
        EXPECT_FALSE(mixed_result->back().dopri);  // Original flat Hermite route.
        EXPECT_EQ(mixed_result->back().attempted_stages, 21U);
        EXPECT_EQ(mixed_result->back().error_ratio, 0);
        EXPECT_EQ(mixed_result->back().full.physical[0].Center(), .75);
        EXPECT_EQ(mixed_result->back().full.physical[3].Center(), 4.25);
        const auto mixed_packed = AttemptRetainedDopriIntervals(compute, mixed, compute.Capacity());
        ASSERT_TRUE(mixed_packed) << mixed_packed.error().Description();
        ASSERT_EQ(mixed_packed->size(), mixed_result->size());
        for (std::size_t row = 0; row < mixed.size(); ++row) {
            SCOPED_TRACE(row);
            same_attempt((*mixed_packed)[row], (*mixed_result)[row]);
        }
        auto strict = input;
        strict.control.tolerance = 1e-30;
        const auto rejected = AttemptRetainedDopriIntervals(compute, std::span(&strict, 1), 1);
        ASSERT_TRUE(rejected) << rejected.error().Description();
        EXPECT_FALSE(rejected->front().admissible);
        EXPECT_FALSE(rejected->front().dopri);
        EXPECT_NE(rejected->front().failure, sirius::core::CoupledStepFailure::None);
        EXPECT_GT(rejected->front().error_ratio, 1);
        for (const auto* endpoint : {&rejected->front().full, &rejected->front().lower,
                                     &rejected->front().midpoint, &rejected->front().refined}) {
            EXPECT_FALSE(endpoint->valid);
        }
        const auto rejected_packed = AttemptRetainedDopriIntervals(compute, std::span(&strict, 1));
        ASSERT_TRUE(rejected_packed) << rejected_packed.error().Description();
        ASSERT_EQ(rejected_packed->size(), rejected->size());
        same_attempt(rejected_packed->front(), rejected->front());
        const auto recovered = AttemptRetainedDopriIntervals(compute, std::span(&input, 1));
        ASSERT_TRUE(recovered) << recovered.error().Description();
        ASSERT_TRUE(recovered->front().admissible);
        ASSERT_TRUE(recovered->front().dopri);
        SameEndpoint(recovered->front().full, accepted.full);
        SameEndpoint(recovered->front().lower, accepted.lower);
        SameEndpoint(recovered->front().midpoint, accepted.midpoint);
        SameEndpoint(recovered->front().refined, accepted.refined);
        SameEndpoint(capsule.starts[second_half], endpoints[1]);
        SameEndpoint(capsule.endpoints[full], endpoints[0]);
        EXPECT_EQ(device->BufferAllocationBytes(), allocation);
    });
#endif
}

TEST_F(RetainedDopriTest, SamplerPreservesPhysicalArrivalAndRejectsInconsistentRates) {
#ifdef SIRIUS_HAS_RETAINED_COMPUTE
    BothProductModes([&](RetainedCompute& compute) {
        const auto allocation = device->BufferAllocationBytes();
        const auto original = std::bit_cast<RetainedStepInput>(reference::connection.input);
        RetainedEndpointInput launch;
        std::copy_n(original.values.begin(), 45, launch.values.begin());
        const auto projected_start = compute.Endpoint(std::span(&launch, 1));
        ASSERT_TRUE(projected_start) << projected_start.error().Description();
        ASSERT_TRUE(projected_start->front().valid);
        RetainedIntervalInput input;
        std::copy_n(original.values.begin(), 4, input.metric.begin());
        input.start = projected_start->front();
        input.chart = original.values[44].Center();
        input.interval = original.values[45].Center();
        input.control.length_scale = input.control.frequency_scale = 1;
        input.control.tolerance = 1e-4 / (4 * 30000);
        input.control.integrator.min_step = static_cast<float>(input.interval);
        input.control.integrator.max_step = 1;
        input.control.integrator.abs_tolerance = input.control.integrator.rel_tolerance = 1e-9f;
        const auto interval = AttemptRetainedDopriIntervals(compute, std::span(&input, 1));
        ASSERT_TRUE(interval) << interval.error().Description();
        ASSERT_TRUE(interval->front().admissible);
        ASSERT_TRUE(interval->front().dopri);
        const auto capsule = interval->front().dopri;
        constexpr auto full = static_cast<std::size_t>(sirius::core::CoupledTrial::Full);
        std::array<RetainedDopriSampleInput, 8> endpoint_inputs;
        for (std::size_t trial = 0; trial < 4; ++trial) {
            endpoint_inputs[2 * trial] = {capsule, trial, 0, std::nullopt};
            endpoint_inputs[2 * trial + 1] = {capsule, trial, 1, std::nullopt};
        }
        const auto endpoint_samples = SampleRetainedDopriIntervals(compute, endpoint_inputs, 8);
        ASSERT_TRUE(endpoint_samples) << endpoint_samples.error().Description();
        ASSERT_EQ(endpoint_samples->size(), endpoint_inputs.size());
        for (std::size_t row = 0; row < endpoint_inputs.size(); ++row) {
            SCOPED_TRACE(row);
            ASSERT_TRUE((*endpoint_samples)[row].valid);
            const auto trial = row / 2;
            const auto& expected = row & 1 ? capsule->endpoints[trial] : capsule->starts[trial];
            for (std::size_t field = 0; field < 40; ++field) {
                EXPECT_EQ((std::bit_cast<std::array<std::uint32_t, 5>>(
                              (*endpoint_samples)[row].physical[field])),
                          (std::bit_cast<std::array<std::uint32_t, 5>>(expected.physical[field])));
            }
        }

        RetainedDopriSampleInput request{capsule, full, 3. / 8, std::nullopt};
        const auto fixed = SampleRetainedDopriIntervals(compute, std::span(&request, 1), 1);
        ASSERT_TRUE(fixed) << fixed.error().Description();
        ASSERT_TRUE(fixed->front().valid);
        auto packet = capsule->packets[full];
        packet.values[361] = RetainedValue::FromDouble(request.fraction);
        const auto phase = compute.DopriPhase(std::span(&packet, 1));
        ASSERT_TRUE(phase) << phase.error().Description();
        ASSERT_TRUE(phase->front().valid);
        RetainedEndpointInput projected_packet = launch;
        std::copy(phase->front().phase.begin(), phase->front().phase.end(),
                  projected_packet.values.begin() + 4);
        const auto direct = compute.Endpoint(std::span(&projected_packet, 1));
        ASSERT_TRUE(direct) << direct.error().Description();
        ASSERT_TRUE(direct->front().valid);
        for (std::size_t field = 0; field < 40; ++field) {
            EXPECT_EQ(
                (std::bit_cast<std::array<std::uint32_t, 5>>(fixed->front().physical[field])),
                (std::bit_cast<std::array<std::uint32_t, 5>>(direct->front().physical[field])));
            CheckValue(fixed->front().physical[field], reference::connection.reference[520 + field],
                       1e-29, 1e-11);
        }
        for (std::size_t axis = 0; axis < 4; ++axis) {
            EXPECT_EQ(
                (std::bit_cast<std::array<std::uint32_t, 5>>(
                    fixed->front().polynomial_tangent[axis])),
                (std::bit_cast<std::array<std::uint32_t, 5>>(phase->front().derivative[axis])));
            CheckValue(fixed->front().polynomial_tangent[axis],
                       reference::connection.reference[320 + axis]);
        }

        sirius::core::Vec4 normal;
        normal(0) = 1;
        normal(1) = .25;
        normal(2) = -.125;
        normal(3) = .0625;
        request.normal = normal;
        const auto moving = SampleRetainedDopriIntervals(compute, std::span(&request, 1), 1);
        ASSERT_TRUE(moving) << moving.error().Description();
        ASSERT_TRUE(moving->front().valid);
        sirius::core::Twofold denominator;
        for (int axis = 0; axis < 4; ++axis) {
            denominator += Center(fixed->front().physical[4 + axis]) * normal(axis);
            CheckValue(moving->front().physical[axis],
                       {Center(fixed->front().physical[axis]).hi,
                        Center(fixed->front().physical[axis]).lo, 0});
            CheckValue(moving->front().physical[4 + axis],
                       {Center(fixed->front().physical[4 + axis]).hi,
                        Center(fixed->front().physical[4 + axis]).lo, 0});
        }
        ASSERT_GT(std::abs(denominator.Rounded()), .1);
        for (std::size_t column = 0; column < 4; ++column) {
            sirius::core::Twofold numerator, arrived;
            double arrived_radius = 0;
            for (int axis = 0; axis < 4; ++axis) {
                numerator += Center(fixed->front().physical[8 + 8 * column + axis]) * normal(axis);
            }
            const auto shift = -numerator / denominator;
            for (int axis = 0; axis < 4; ++axis) {
                const auto displacement = Center(fixed->front().physical[8 + 8 * column + axis]) +
                                          Center(fixed->front().physical[4 + axis]) * shift;
                CheckValue(moving->front().physical[8 + 8 * column + axis],
                           {displacement.hi, displacement.lo, 0});
                const auto covariant = Center(fixed->front().physical[12 + 8 * column + axis]);
                CheckValue(moving->front().physical[12 + 8 * column + axis],
                           {covariant.hi, covariant.lo, 0});
                arrived += Center(moving->front().physical[8 + 8 * column + axis]) * normal(axis);
                arrived_radius +=
                    moving->front().physical[8 + 8 * column + axis].radius * std::abs(normal(axis));
            }
            EXPECT_LE(std::abs(arrived.Rounded()), arrived_radius + 1e-28);
            EXPECT_LE(std::abs(arrived.Rounded()), 1e-10);
        }

        // Same-sign, non-tangential denominators can still disagree beyond
        // the original X budget. At s=1 the physical endpoint is immutable;
        // changing only the supplied k7 time rate gives a controlled witness.
        const auto& terminal = (*endpoint_samples)[2 * full + 1];
        auto inconsistent = std::make_shared<RetainedDopriInterval>(*capsule);
        inconsistent->packets[full].values[320] =
            RetainedValue::FromDouble(2 * Center(terminal.polynomial_tangent[0]).Rounded());
        const auto k = Center(terminal.physical[4]);
        const auto w = Center(inconsistent->packets[full].values[320]);
        EXPECT_EQ(std::signbit(k.Rounded()), std::signbit(w.Rounded()));
        double largest_ratio = 0;
        for (std::size_t column = 0; column < 4; ++column) {
            const auto numerator = Center(terminal.physical[8 + 8 * column]);
            const auto shift_difference = -numerator / k + numerator / w;
            for (int axis = 0; axis < 4; ++axis) {
                const auto central = Center(terminal.physical[4 + axis]);
                const auto physical_x =
                    Center(terminal.physical[8 + 8 * column + axis]) - central * (numerator / k);
                const double budget =
                    capsule->control.tolerance *
                    (capsule->control.length_scale * capsule->control.column_scale[column] +
                     std::abs(physical_x.Rounded()));
                largest_ratio = std::max(largest_ratio,
                                         std::abs((central * shift_difference).Rounded()) / budget);
            }
        }
        ASSERT_GT(largest_ratio, 1);
        sirius::core::Vec4 temporal_normal;
        temporal_normal(0) = 1;
        RetainedDopriSampleInput disagreement{inconsistent, full, 1, temporal_normal};
        const auto declined = SampleRetainedDopriIntervals(compute, std::span(&disagreement, 1), 1);
        ASSERT_TRUE(declined) << declined.error().Description();
        EXPECT_FALSE(declined->front().valid);

        sirius::core::Vec4 tangent_normal;
        tangent_normal(0) = Center(terminal.physical[5]).Rounded();
        tangent_normal(1) = -Center(terminal.physical[4]).Rounded();
        auto malformed_control = std::make_shared<RetainedDopriInterval>(*capsule);
        malformed_control->control.tolerance = 0;
        auto malformed_endpoint = std::make_shared<RetainedDopriInterval>(*capsule);
        malformed_endpoint->endpoints[full].physical.back().valid = 0;
        sirius::core::Vec4 invalid_normal = normal;
        invalid_normal(2) = std::numeric_limits<double>::quiet_NaN();
        const std::array<RetainedDopriSampleInput, 10> invalid{{
            {nullptr, full, .5, std::nullopt},
            {capsule, 4, .5, std::nullopt},
            {capsule, full, -1, std::nullopt},
            {capsule, full, std::nextafter(1., 2.), std::nullopt},
            {capsule, full, std::numeric_limits<double>::quiet_NaN(), std::nullopt},
            {malformed_control, full, .5, std::nullopt},
            {malformed_endpoint, full, 1, std::nullopt},
            {capsule, full, 1, sirius::core::Vec4{}},
            {capsule, full, 1, tangent_normal},
            {capsule, full, .5, invalid_normal},
        }};
        std::vector<RetainedDopriSampleInput> mixed(invalid.begin(), invalid.end());
        mixed.push_back({capsule, full, 3. / 8, std::nullopt});
        const auto mixed_result = SampleRetainedDopriIntervals(compute, mixed, mixed.size());
        ASSERT_TRUE(mixed_result) << mixed_result.error().Description();
        ASSERT_EQ(mixed_result->size(), mixed.size());
        for (std::size_t row = 0; row < invalid.size(); ++row) {
            SCOPED_TRACE(row);
            EXPECT_FALSE((*mixed_result)[row].valid);
            for (const auto& value : (*mixed_result)[row].physical) {
                EXPECT_EQ((std::bit_cast<std::array<std::uint32_t, 5>>(value)),
                          (std::array<std::uint32_t, 5>{}));
            }
            for (const auto& value : (*mixed_result)[row].polynomial_tangent) {
                EXPECT_EQ((std::bit_cast<std::array<std::uint32_t, 5>>(value)),
                          (std::array<std::uint32_t, 5>{}));
            }
        }
        ASSERT_TRUE(mixed_result->back().valid);
        for (std::size_t field = 0; field < 40; ++field) {
            EXPECT_EQ(
                (std::bit_cast<std::array<std::uint32_t, 5>>(mixed_result->back().physical[field])),
                (std::bit_cast<std::array<std::uint32_t, 5>>(fixed->front().physical[field])));
        }
        EXPECT_FALSE(SampleRetainedDopriIntervals(compute, {}, 1));
        EXPECT_FALSE(SampleRetainedDopriIntervals(compute, mixed, 1));
        EXPECT_EQ(device->BufferAllocationBytes(), allocation);

        // Queue the same finite scene through the real executor. Rollback and
        // trace end revoke sampling; no completed render is asserted here.
        sirius::core::KerrSchildFamily family({input.metric[0].Center(), input.metric[1].Center(),
                                               input.metric[2].Center(), input.metric[3].Center()});
        sirius::core::OutgoingKerrSchild metric(family);
        sirius::core::Lightray ray{};
        ray.step_size = static_cast<float>(input.interval);
        sirius::core::Rk45CoupledState coupled;
        coupled.length_scale = input.control.length_scale;
        coupled.frequency_scale = input.control.frequency_scale;
        coupled.tolerance = input.control.tolerance;
        coupled.column_scale = input.control.column_scale;
        for (int axis = 0; axis < 4; ++axis) {
            ray.position(axis) = Center(input.start.physical[axis]).Rounded();
            ray.velocity(axis) = Center(input.start.physical[4 + axis]).Rounded();
            for (int column = 0; column < 4; ++column) {
                coupled.variations[column].displacement(axis) =
                    Center(input.start.physical[8 + 8 * column + axis]).Rounded();
                coupled.variations[column].derivative(axis) =
                    Center(input.start.physical[12 + 8 * column + axis]).Rounded();
            }
        }
        bool cancelled = false;
        RetainedTraceExecutor executor(compute, [&] { return cancelled; }, 0, 1);
        sirius::core::Rk45CoupledComparison comparison;
        executor.BeginTrace();
        const bool stepped =
            executor.Step(ray, metric, input.control.integrator, coupled, comparison);
        EXPECT_TRUE(stepped);
        if (stepped) {
            const auto sampled = executor.Sample(sirius::core::CoupledTrial::Full, 1);
            EXPECT_TRUE(sampled);
            if (sampled) {
                EXPECT_TRUE(sampled->polynomial_tangent);
                for (int axis = 0; axis < 4; ++axis) {
                    EXPECT_EQ(sampled->ray.position(axis), ray.position(axis));
                    EXPECT_EQ(sampled->ray.velocity(axis), ray.velocity(axis));
                }
            }
            EXPECT_EQ(executor.Statistics().sample_rows, 0U);
            const auto interior_sample = executor.Sample(sirius::core::CoupledTrial::Full, 3. / 8);
            EXPECT_TRUE(interior_sample);
            if (interior_sample) {
                EXPECT_TRUE(interior_sample->polynomial_tangent);
                for (std::size_t field = 0; field < 40; ++field) {
                    const auto axis = static_cast<int>(field % 4);
                    const double value =
                        field < 4   ? interior_sample->ray.position(axis)
                        : field < 8 ? interior_sample->ray.velocity(axis)
                        : (field - 8) % 8 < 4
                            ? interior_sample->variations[(field - 8) / 8].displacement(axis)
                            : interior_sample->variations[(field - 8) / 8].derivative(axis);
                    const auto expected = reference::connection.reference[520 + field];
                    EXPECT_TRUE(std::isfinite(value));
                    EXPECT_LE(std::abs(value - (expected.high + expected.low)),
                              1e-11 * (1 + std::abs(expected.high)));
                }
                if (interior_sample->polynomial_tangent) {
                    for (int axis = 0; axis < 4; ++axis) {
                        const auto expected = reference::connection.reference[320 + axis];
                        const double value = (*interior_sample->polynomial_tangent)(axis);
                        EXPECT_TRUE(std::isfinite(value));
                        EXPECT_LE(std::abs(value - (expected.high + expected.low)),
                                  1e-10 * (1 + std::abs(expected.high)));
                    }
                }
            }
        }
        // Starting a new attempt revokes the old capsule even when an early
        // cancellation performs no device work. Clear the cancellation before
        // querying, so the refusal demonstrates capsule lifetime ownership.
        cancelled = true;
        EXPECT_FALSE(executor.Step(ray, metric, input.control.integrator, coupled, comparison));
        cancelled = false;
        EXPECT_FALSE(executor.Sample(sirius::core::CoupledTrial::Full, .5));
        executor.RejectLastInterval();
        executor.RejectLastInterval();
        EXPECT_FALSE(executor.Sample(sirius::core::CoupledTrial::Full, .5));
        executor.EndTrace();
        EXPECT_FALSE(executor.Sample(sirius::core::CoupledTrial::Full, .5));
        EXPECT_FALSE(executor.Error());
        const auto statistics = executor.Statistics();
        EXPECT_EQ(statistics.accepted_intervals, 1U);
        EXPECT_EQ(statistics.sample_rows, 1U);
        EXPECT_EQ(statistics.tracer_rollbacks, 1U);
        EXPECT_EQ(coupled.central_stages, 21U);
        EXPECT_EQ(coupled.variation_stages, 21U);
        EXPECT_EQ(device->BufferAllocationBytes(), allocation);
    });
#endif
}
}  // namespace
