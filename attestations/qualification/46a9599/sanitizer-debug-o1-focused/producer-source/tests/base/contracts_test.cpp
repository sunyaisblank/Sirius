// Postconditions verified here: enforce-mode contracts terminate on
// violation, satisfied contracts are silent, and axioms never evaluate.
// This target compiles with SIRIUS_CONTRACT_MODE=2 (tests/CMakeLists.txt).

#include "sirius/base/contracts.h"

#include "sirius/base/language_capabilities.h"

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
    EXPECT_DEATH(
        [] {
            volatile int impossible_sum = 5;
            SIRIUS_ASSERT(2 + 2 == impossible_sum);
        }(),
        "assertion violated");
}

TEST(Contracts, AxiomIsNeverEvaluated) {
    int evaluations = 0;
    SIRIUS_AXIOM(++evaluations > 0);
    EXPECT_EQ(evaluations, 0);
}

TEST(Contracts, LanguageCapabilitiesDistinguishNativeFeaturesFromSubstitutes) {
#if defined(__cpp_contracts) && __cpp_contracts >= 202502L
    EXPECT_TRUE(sirius::base::kNativeP2900Contracts);
    EXPECT_EQ(sirius::base::kContractImplementation,
              sirius::base::ContractImplementation::NativeP2900);
#else
    EXPECT_FALSE(sirius::base::kNativeP2900Contracts);
    EXPECT_EQ(sirius::base::kContractImplementation,
              sirius::base::ContractImplementation::CheckedMacro);
#endif

    // Reflection adoption is deliberate: current bindings remain explicit even
    // on an experimental compiler so the representation cannot change merely
    // because a toolchain feature macro appeared.
    EXPECT_EQ(sirius::base::kReflectionImplementation,
              sirius::base::ReflectionImplementation::ExplicitSchema);
}

}  // namespace
