// Postconditions verified here: Expected propagates values untouched, Fail
// carries domain, operation, and detail to the caller, and Description
// formats all three for an operator.

#include "sirius/base/error.h"

#include <gtest/gtest.h>

#include <string>

namespace {

using sirius::base::ErrorDomain;
using sirius::base::Expected;
using sirius::base::Fail;

Expected<int> ParsePositive(int value) {
    if (value <= 0) {
        return Fail(ErrorDomain::kConfiguration, "parse positive integer",
                    std::to_string(value));
    }
    return value;
}

TEST(Error, ExpectedPropagatesValue) {
    const auto result = ParsePositive(7);
    ASSERT_TRUE(result.has_value());
    EXPECT_EQ(*result, 7);
}

TEST(Error, FailCarriesDomainOperationAndDetail) {
    const auto result = ParsePositive(-3);
    ASSERT_FALSE(result.has_value());
    EXPECT_EQ(result.error().domain(), ErrorDomain::kConfiguration);
    EXPECT_EQ(result.error().operation(), "parse positive integer");
    EXPECT_EQ(result.error().detail(), "-3");
}

TEST(Error, DescriptionNamesDomainOperationAndDetail) {
    const auto result = ParsePositive(0);
    ASSERT_FALSE(result.has_value());
    EXPECT_EQ(result.error().Description(), "configuration: parse positive integer: 0");
}

}  // namespace
