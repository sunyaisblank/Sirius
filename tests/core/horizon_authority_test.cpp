#include "sirius/core/metrics/kerr_schild_family.h"
#include "sirius/core/metrics/outgoing_kerr_schild.h"
#include <gtest/gtest.h>
#include <array>
#include <cmath>
#include <limits>

namespace sirius::core {
namespace {

void ExpectAbsent(const KerrSchildParams& p) {
    KerrSchildFamily metric(p);
    EXPECT_FALSE(metric.HasHorizon());
    EXPECT_EQ(metric.OuterHorizonRadius(), -1.0);
    EXPECT_EQ(metric.InnerHorizonRadius(), -1.0);
}

TEST(HorizonAuthority, TinyChargeAtEqualMassAndSpinIsExactlySuperextremal) {
    const double charge = static_cast<double>(1.0e-20f);
    // Exact represented relation: M²-a²-Q²=-Q²<0, regardless of rounded sums.
    const std::array<KerrSchildParams, 5> cases{{
        {1.0, -1.0, charge, 0.0},
        {1.0, charge, -1.0, 0.0},
        {1.0, 1.0, charge, 0.0},
        {1.0, charge, 1.0, 0.0},
        {1.0, 1.0, std::numeric_limits<double>::denorm_min(), 0.0}
    }};
    for (const auto& parameters : cases) {
        KerrSchildFamily metric(parameters);
        EXPECT_FALSE(metric.HasHorizon());
        EXPECT_EQ(metric.OuterHorizonRadius(), -1.0);
        EXPECT_EQ(metric.InnerHorizonRadius(), -1.0);
    }
}

TEST(HorizonAuthority, RoundedPythagoreanEqualityDoesNotDefineExtremality) {
    // Exact binary64 0.6²+0.8² exceeds 1 (independent rational reference).
    ExpectAbsent({1.0, 0.6, 0.8, 0.0});
    ExpectAbsent({5.0, 3.0, std::nextafter(4.0, 5.0), 0.0});
    KerrSchildFamily exact({5.0, 3.0, 4.0, 0.0});
    ASSERT_TRUE(exact.HasHorizon());
    EXPECT_EQ(exact.OuterHorizonRadius(), 5.0);
    EXPECT_EQ(exact.InnerHorizonRadius(), 5.0);
    KerrSchildFamily below({5.0, 3.0, std::nextafter(4.0, 3.0), 0.0});
    ASSERT_TRUE(below.HasHorizon());
    EXPECT_GT(below.OuterHorizonRadius(), 5.0);
    EXPECT_LT(below.InnerHorizonRadius(), 5.0);
}

TEST(HorizonAuthority, NearExtremalFloatSpinKeepsDistinctFiniteRoots) {
    const double spin = static_cast<double>(std::nextafter(1.0f, 0.0f));
    KerrSchildFamily metric({1.0, spin, 0.0, 0.0});
    ASSERT_TRUE(metric.HasHorizon());
    // Independent 100-digit decimal roots of the exact represented spin.
    EXPECT_NEAR(metric.OuterHorizonRadius(), 1.0003452669778563649, 5.0e-16);
    EXPECT_NEAR(metric.InnerHorizonRadius(), 0.9996547330221436351, 5.0e-16);
    EXPECT_GT(metric.OuterHorizonRadius(), metric.InnerHorizonRadius());
}

TEST(HorizonAuthority, ScalingAvoidsOverflowAndUnderflowOfSquaredParameters) {
    for (const int exponent : {-600, 0, 600}) {
        const double scale = std::ldexp(1.0, exponent);
        KerrSchildFamily metric({5.0 * scale, 3.0 * scale, 0.0, 0.0});
        ASSERT_TRUE(metric.HasHorizon());
        EXPECT_DOUBLE_EQ(metric.OuterHorizonRadius() / scale, 9.0);
        EXPECT_DOUBLE_EQ(metric.InnerHorizonRadius() / scale, 1.0);
        KerrSchildFamily extreme({5.0 * scale, 3.0 * scale, 4.0 * scale, 0.0});
        ASSERT_TRUE(extreme.HasHorizon());
        EXPECT_EQ(extreme.OuterHorizonRadius(), 5.0 * scale);
        EXPECT_EQ(extreme.InnerHorizonRadius(), 5.0 * scale);
    }
}

TEST(HorizonAuthority, SmallPositiveInnerRootSurvivesCancellation) {
    KerrSchildFamily metric({1.0, 0.0, 1.0e-150, 0.0});
    ASSERT_TRUE(metric.HasHorizon());
    EXPECT_EQ(metric.OuterHorizonRadius(), 2.0);
    ASSERT_GT(metric.InnerHorizonRadius(), 0.0);
    // r-=Q²/r+; dividing first avoids both square overflow and subtraction.
    EXPECT_NEAR(metric.InnerHorizonRadius() / 5.0e-301, 1.0, 8.0e-16);
}

TEST(HorizonAuthority, UnrepresentableRootPairDoesNotSelectHorizonFreeRoute) {
    for (const KerrSchildParams p : {
             KerrSchildParams{std::numeric_limits<double>::max(), 0.0, 0.0, 0.0},
             KerrSchildParams{1.0, 0.0, std::numeric_limits<double>::denorm_min(), 0.0}}) {
        KerrSchildFamily metric(p);
        ASSERT_TRUE(metric.HasHorizon());
        EXPECT_EQ(metric.OuterHorizonRadius(), -1.0);
        EXPECT_EQ(metric.InnerHorizonRadius(), -1.0);
        OutgoingKerrSchild outgoing(metric);
        Tensor<double, 4> event;
        event(0) = 0.0; event(1) = 10.0; event(2) = 0.0; event(3) = 0.0;
        EXPECT_FALSE(outgoing.FromIngoing(event).has_value());
    }
}

TEST(HorizonAuthority, SchwarzschildNativeAndKottlerContractsRemain) {
    KerrSchildFamily schwarzschild({2.0, 0.0, 0.0, 0.0});
    ASSERT_TRUE(schwarzschild.HasHorizon());
    EXPECT_EQ(schwarzschild.OuterHorizonRadius(), 4.0);
    EXPECT_EQ(schwarzschild.InnerHorizonRadius(), 0.0);
    ExpectAbsent(KerrSchildParams::Minkowski());
    ExpectAbsent(KerrSchildParams::DeSitter(0.01));
    KerrSchildFamily kottler({1.0, 0.0, 0.0, 0.01});
    ASSERT_TRUE(kottler.HasHorizon());
    EXPECT_EQ(kottler.OuterHorizonRadius(), KottlerBlackHoleHorizonRadius(1.0, 0.01));
    EXPECT_EQ(kottler.InnerHorizonRadius(), 0.0);
}

}  // namespace
}  // namespace sirius::core
