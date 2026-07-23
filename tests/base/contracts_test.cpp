// Postconditions verified here: enforce-mode contracts terminate on
// violation, satisfied contracts are silent, and axioms never evaluate.
// This target compiles with SIRIUS_CONTRACT_MODE=2 (tests/CMakeLists.txt).

#include "sirius/base/contracts.h"

#include <gtest/gtest.h>

namespace {

int Divide(int numerator, int denominator) {
    SIRIUS_PRE(denominator != 0);
    const int quotient = numerator / denominator;
    SIRIUS_POST(quotient * denominator + numerator % denominator == numerator);
    return quotient;
}

TEST(Contracts, SatisfiedContractsPassSilently) {
    EXPECT_EQ(Divide(10, 2), 5);
    EXPECT_EQ(Divide(-9, 3), -3);
}

TEST(ContractsDeathTest, EnforcedPreconditionViolationTerminates) {
    EXPECT_DEATH(Divide(1, 0), "precondition violated .* denominator != 0");
}

TEST(ContractsDeathTest, EnforcedAssertionViolationTerminates) {
    EXPECT_DEATH([] { SIRIUS_ASSERT(2 + 2 == 5); }(), "assertion violated");
}

TEST(Contracts, AxiomIsNeverEvaluated) {
    int evaluations = 0;
    SIRIUS_AXIOM(++evaluations > 0);
    EXPECT_EQ(evaluations, 0);
}

}  // namespace
