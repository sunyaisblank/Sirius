#include "sirius/core/first_order_number.h"
#include "sirius/core/metrics/outgoing_kerr_schild.h"
#include "sirius/core/twofold.h"

#include <gtest/gtest.h>

#include "support/cpu_critical/metric_reference.h"

#include <cmath>

namespace sirius::test {
using namespace sirius::core;

TEST(RetainedArithmetic, SmallTermsSurviveLargeSumsProductsAndDivision) {
    const double tiny = 0x1p-54;
    const auto sum = Twofold(1) + tiny;
    EXPECT_EQ(sum.hi, 1);
    EXPECT_EQ(sum.lo, tiny);
    EXPECT_EQ((sum - Twofold(1)).Rounded(), tiny);
    const auto product = Twofold::Product(1 + 0x1p-27, 1 - 0x1p-27);
    EXPECT_EQ(product.hi, 1);
    EXPECT_EQ(product.lo, -tiny);
    EXPECT_EQ((product - Twofold(1)).Rounded(), -tiny);
    const auto third = Twofold(1) / Twofold(3);
    EXPECT_EQ(third.hi, 0x1.5555555555555p-2);
    EXPECT_EQ(third.lo, 0x1.5555555555555p-56);
    const auto residual = third * 3.0 - Twofold(1);
    EXPECT_LT(std::abs(residual.hi) + std::abs(residual.lo), 0x1p-104);
}

TEST(RetainedArithmetic, CartesianGradientMatchesIndependentRationalValues) {
    using Jet = FirstOrder3<Twofold>;
    const auto x = Jet::Variable(3, 0);
    const auto y = Jet::Variable(4, 1);
    const auto z = Jet::Variable(12, 2);
    const auto value = sqrt(x * x + y * y + z * z) / (x + y);
    const std::array<Twofold, 4> expected{Twofold(13) / 7.0, Twofold(-148) / 637.0,
                                          Twofold(-141) / 637.0, Twofold(12) / 91.0};
    const std::array<Twofold, 4> actual{value.value, value.gradient[0], value.gradient[1],
                                        value.gradient[2]};
    for (int i = 0; i < 4; ++i) {
        const auto error = actual[i] - expected[i];
        EXPECT_LT(std::abs(error.hi) + std::abs(error.lo), 1e-29);
    }
}

TEST(RetainedArithmetic, CompleteCriticalMetricMatchesIndependentPrecisionWitnesses) {
    KerrSchildFamily ingoing(KerrSchildParams::Kerr(1.0, 0.998));
    OutgoingKerrSchild metric(ingoing);
    int low_part_mutants = 0;
    for (const auto& fixture : critical_fixture::metric_cases) {
        Vec4 position;
        for (int axis = 0; axis < 4; ++axis) position(axis) = fixture.position[axis];
        RetainedMetricSample sample;
        ASSERT_TRUE(metric.EvaluateRetained(position, sample));
        std::array<Twofold, 96> actual;
        int field = 0;
        for (int mu = 0; mu < 4; ++mu)
            for (int nu = 0; nu < 4; ++nu) actual[field++] = sample.metric(mu, nu);
        for (int mu = 0; mu < 4; ++mu)
            for (int nu = 0; nu < 4; ++nu) actual[field++] = sample.inverse(mu, nu);
        for (int axis = 0; axis < 4; ++axis)
            for (int mu = 0; mu < 4; ++mu)
                for (int nu = 0; nu < 4; ++nu) actual[field++] = sample.derivative(axis, mu, nu);
        for (int i = 0; i < field; ++i) {
            const auto& reference = fixture.fields[i];
            const Twofold expected(reference.hi, reference.lo);
            const double bound = 1e-27 * (1 + std::abs(reference.hi));
            const auto error = actual[i] - expected;
            EXPECT_LT(std::abs(error.hi) + std::abs(error.lo), bound) << i;
            const auto mutant = Twofold(actual[i].hi) - expected;
            if (std::abs(mutant.hi) + std::abs(mutant.lo) > bound) ++low_part_mutants;
        }
    }
    // Exact zeros do not have low parts. The remaining independent fields
    // must distinguish retained geometry from a scalar-only substitution.
    EXPECT_GT(low_part_mutants, 700);
}

}  // namespace sirius::test
