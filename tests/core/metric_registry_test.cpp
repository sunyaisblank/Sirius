// TSMT200A.cpp - Metric Identity Registry Tests
// Component ID: TSMT200A
// Tests for: PHMT200A.h
//
// The registry is the single naming authority; these tests pin the
// postconditions that ended the stringly-typed routing defects: every
// canonical name and alias parses to its identity, display names from the
// metric families round-trip, unknown names fail rather than defaulting,
// and the error text can enumerate the accepted set.

#include "sirius/core/metrics/kerr_schild_family.h"
#include "sirius/core/metrics/morris_thorne_family.h"
#include "sirius/core/metrics/registry.h"
#include "sirius/core/metrics/warp_drive_family.h"

#include <gtest/gtest.h>

using namespace sirius::core;

TEST(MetricRegistryTests, EveryCanonicalNameParsesToItsOwnId) {
    for (const auto& info : MetricRegistry()) {
        auto id = ParseMetricName(info.canonical_name);
        ASSERT_TRUE(id.has_value()) << "canonical name failed: " << info.canonical_name;
        EXPECT_EQ(*id, info.id) << "canonical name misrouted: " << info.canonical_name;
    }
}

TEST(MetricRegistryTests, EveryAliasParsesToItsOwnId) {
    for (const auto& info : MetricRegistry()) {
        for (const char* alias : info.aliases) {
            if (alias == nullptr) continue;
            auto id = ParseMetricName(alias);
            ASSERT_TRUE(id.has_value()) << "alias failed: " << alias;
            EXPECT_EQ(*id, info.id) << "alias misrouted: " << alias;
        }
    }
}

TEST(MetricRegistryTests, ParsingIsCaseInsensitive) {
    EXPECT_EQ(ParseMetricName("schwarzschild"), MetricId::Schwarzschild);
    EXPECT_EQ(ParseMetricName("KERR"), MetricId::Kerr);
    EXPECT_EQ(ParseMetricName("morristhorne"), MetricId::MorrisThorne);
    EXPECT_EQ(ParseMetricName("wormhole"), MetricId::MorrisThorne);
    EXPECT_EQ(ParseMetricName("warpdrive"), MetricId::Alcubierre);
}

TEST(MetricRegistryTests, UnknownNamesFailInsteadOfDefaulting) {
    EXPECT_FALSE(ParseMetricName("").has_value());
    EXPECT_FALSE(ParseMetricName("Schwarzchild").has_value());  // common misspelling
    EXPECT_FALSE(ParseMetricName("Godel").has_value());         // unimplemented on purpose
    EXPECT_FALSE(ParseMetricName("Taub-NUT").has_value());
    EXPECT_FALSE(ParseMetricName("Kerr ").has_value());  // trailing space is not a name
}

TEST(MetricRegistryTests, MetricInfoRoundTripsById) {
    for (const auto& info : MetricRegistry()) {
        EXPECT_EQ(MetricInfoFor(info.id).id, info.id);
        EXPECT_STREQ(MetricInfoFor(info.id).canonical_name, info.canonical_name);
    }
    EXPECT_DEATH(static_cast<void>(MetricInfoFor(static_cast<MetricId>(255))), "violated");
    EXPECT_DEATH(static_cast<void>(DiskSupportFor(static_cast<MetricId>(255))), "violated");
    EXPECT_DEATH(static_cast<void>(MetricUsesMass(static_cast<MetricId>(255))), "violated");
    EXPECT_DEATH(static_cast<void>(ToString(static_cast<DiskSupport>(255))), "violated");

    EXPECT_TRUE(MetricUsesMass(MetricId::Schwarzschild));
    EXPECT_TRUE(MetricUsesMass(MetricId::KerrNewman));
    EXPECT_TRUE(MetricUsesMass(MetricId::SchwarzschildDeSitter));
    EXPECT_FALSE(MetricUsesMass(MetricId::Minkowski));
    EXPECT_FALSE(MetricUsesMass(MetricId::DeSitter));
    EXPECT_FALSE(MetricUsesMass(MetricId::MorrisThorne));
    EXPECT_FALSE(MetricUsesMass(MetricId::Alcubierre));

    EXPECT_FALSE(MetricSpecificParameterIssue(
                     MetricId::MorrisThorne, 2.0, true, kDefaultAlcubierreWarpVelocity,
                     kDefaultAlcubierreBubbleRadius, kDefaultAlcubierreBubbleSigma)
                     .has_value());
    EXPECT_TRUE(MetricSpecificParameterIssue(
                    MetricId::Schwarzschild, 2.0, true, kDefaultAlcubierreWarpVelocity,
                    kDefaultAlcubierreBubbleRadius, kDefaultAlcubierreBubbleSigma)
                    .has_value());
    EXPECT_TRUE(MetricSpecificParameterIssue(
                    MetricId::MorrisThorne, kDefaultMorrisThorneThroatRadius, true, 1.0,
                    kDefaultAlcubierreBubbleRadius, kDefaultAlcubierreBubbleSigma)
                    .has_value());
    EXPECT_FALSE(MetricSpecificParameterIssue(
                     MetricId::Alcubierre, kDefaultMorrisThorneThroatRadius, true, 1.0, 2.0, 0.25)
                     .has_value());
    EXPECT_TRUE(
        MetricSpecificParameterIssue(static_cast<MetricId>(255), kDefaultMorrisThorneThroatRadius,
                                     true, kDefaultAlcubierreWarpVelocity,
                                     kDefaultAlcubierreBubbleRadius, kDefaultAlcubierreBubbleSigma)
            .has_value());

    EXPECT_FALSE(MetricParameterIssue(MetricId::Kerr, 0.9, 0.0, 0.0).has_value());
    EXPECT_TRUE(MetricParameterIssue(MetricId::Kerr, 0.9, 0.1, 0.0).has_value());
    EXPECT_TRUE(MetricParameterIssue(MetricId::MorrisThorne, 0.1, 0.0, 0.0).has_value());
    EXPECT_TRUE(MetricParameterIssue(MetricId::Alcubierre, 0.0, 0.1, 0.0).has_value());
    EXPECT_TRUE(MetricParameterIssue(MetricId::DeSitter, 0.0, 0.0, 0.0).has_value());
    EXPECT_FALSE(MetricParameterIssue(MetricId::DeSitter, 0.0, 0.0, 0.01).has_value());
    EXPECT_FALSE(MetricHorizonIssue(MetricId::SchwarzschildDeSitter, 1.0, 0.01).has_value());
    EXPECT_TRUE(MetricHorizonIssue(MetricId::SchwarzschildDeSitter, 2.0, 0.03).has_value());
    EXPECT_TRUE(MetricParameterIssue(static_cast<MetricId>(255), 0.0, 0.0, 0.0).has_value());
}

TEST(MetricRegistryTests, PositiveLambdaObserverAndHorizonShareOneCausalDomain) {
    constexpr double lambda = 0.001;
    const double de_sitter_horizon = KottlerCosmologicalHorizonRadius(0.0, lambda);
    ASSERT_TRUE(std::isfinite(de_sitter_horizon));
    EXPECT_NEAR(KottlerStaticLapse(0.0, lambda, de_sitter_horizon), 0.0, 3.0e-16);
    ASSERT_TRUE(MetricCosmologicalHorizonRadius(MetricId::DeSitter, 0.0, lambda).has_value());
    EXPECT_DOUBLE_EQ(*MetricCosmologicalHorizonRadius(MetricId::DeSitter, 0.0, lambda),
                     de_sitter_horizon);
    EXPECT_EQ(MetricObserverRadiusIssueFor(MetricId::DeSitter, 0.0, lambda, 50.0, 1.0, 1.0, 0.5),
              MetricObserverRadiusIssue::None);
    EXPECT_EQ(MetricObserverRadiusIssueFor(MetricId::DeSitter, 0.0, lambda, 55.0, 1.0, 1.0, 0.5),
              MetricObserverRadiusIssue::CosmologicalHorizon);

    KerrSchildFamily kottler({2.0, 0.0, 0.0, lambda});
    const double kottler_horizon = KottlerCosmologicalHorizonRadius(2.0, lambda);
    EXPECT_DOUBLE_EQ(kottler.CosmologicalHorizonRadius(), kottler_horizon);
    EXPECT_DOUBLE_EQ(kottler.OuterHorizonRadius(), KottlerBlackHoleHorizonRadius(2.0, lambda));
    EXPECT_EQ(MetricObserverRadiusIssueFor(MetricId::SchwarzschildDeSitter, 2.0, lambda, 50.0, 1.0,
                                           1.0, 0.5),
              MetricObserverRadiusIssue::None);
    EXPECT_EQ(MetricObserverRadiusIssueFor(MetricId::SchwarzschildDeSitter, 2.0, 0.02, 50.0, 1.0,
                                           1.0, 0.5),
              MetricObserverRadiusIssue::CosmologicalHorizon);
    EXPECT_FALSE(MetricCosmologicalHorizonRadius(MetricId::Kerr, 2.0, 0.0).has_value());
    EXPECT_FALSE(MetricCosmologicalHorizonRadius(MetricId::DeSitter, 0.0, 0.0).has_value());
    EXPECT_FALSE(
        MetricCosmologicalHorizonRadius(MetricId::SchwarzschildDeSitter, 2.0, 0.03).has_value());

    const double tiny_lambda_horizon =
        KottlerCosmologicalHorizonRadius(0.0, std::numeric_limits<double>::denorm_min());
    EXPECT_TRUE(std::isfinite(tiny_lambda_horizon));
    EXPECT_NEAR(
        KottlerStaticLapse(0.0, std::numeric_limits<double>::denorm_min(), tiny_lambda_horizon),
        0.0, 3.0e-16);
}

TEST(MetricRegistryTests, KnownMetricNamesListsEveryCanonicalName) {
    std::string names = KnownMetricNames();
    for (const auto& info : MetricRegistry()) {
        EXPECT_NE(names.find(info.canonical_name), std::string::npos)
            << "missing from error text: " << info.canonical_name;
    }
}

// Display names produced by the metric families must resolve through the
// registry (the GL pipeline routes IMetric::GetName() output through it).
TEST(MetricRegistryTests, FamilyDisplayNamesRoundTrip) {
    KerrSchildFamily schw(KerrSchildParams::Schwarzschild(1.0));
    EXPECT_EQ(ParseMetricName(schw.GetName()), MetricId::Schwarzschild);

    KerrSchildFamily kerr(KerrSchildParams::Kerr(1.0, 0.9));
    EXPECT_EQ(ParseMetricName(kerr.GetName()), MetricId::Kerr);

    KerrSchildFamily rn(KerrSchildParams::ReissnerNordstrom(1.0, 0.5));
    EXPECT_EQ(ParseMetricName(rn.GetName()), MetricId::ReissnerNordstrom);

    KerrSchildFamily kn(KerrSchildParams::KerrNewman(1.0, 0.5, 0.3));
    EXPECT_EQ(ParseMetricName(kn.GetName()), MetricId::KerrNewman);

    KerrSchildFamily mink(KerrSchildParams::Minkowski());
    EXPECT_EQ(ParseMetricName(mink.GetName()), MetricId::Minkowski);

    MorrisThorneFamily wormhole(MorrisThorneParams::Ellis(1.0));
    EXPECT_EQ(ParseMetricName(wormhole.GetName()), MetricId::MorrisThorne);

    WarpDriveFamily warp(WarpDriveParams::Alcubierre(1.0, 1.0));
    EXPECT_EQ(ParseMetricName(warp.GetName()), MetricId::Alcubierre);

    EXPECT_DEATH(
        {
            KerrSchildParams invalid = KerrSchildParams::Kerr(1.0, 0.5);
            invalid.Lambda = 0.01;
            KerrSchildFamily rotating_de_sitter(invalid);
        },
        "violated");
    EXPECT_DEATH(schw.SetParameter("mas", 1.0), "violated");
    EXPECT_DEATH(wormhole.SetParameter("throat", 1.0), "violated");
    EXPECT_DEATH(warp.SetParameter("speed", 1.0), "violated");
    EXPECT_DEATH(schw.SetParameter("mass", std::numeric_limits<double>::quiet_NaN()), "violated");
    EXPECT_DEATH(schw.SetParameter("mass", 1000001.0), "violated");
    EXPECT_DEATH(wormhole.SetParameter("throat_radius", 1001.0), "violated");
    EXPECT_DEATH(wormhole.SetParameter("throat_radius", 0.09), "violated");
    EXPECT_DEATH(warp.SetParameter("velocity", 11.0), "violated");
}

// Backend support flags document reality; these pin the deliberate ones so a
// flag change is a conscious decision, not a drive-by edit.
TEST(MetricRegistryTests, BackendSupportMatchesImplementations) {
    // Morris-Thorne renders on the CPU through the exact isotropic Cartesian
    // chart of the zero-tidal Ellis member.
    EXPECT_TRUE(MetricInfoFor(MetricId::MorrisThorne).cpu_supported);
    EXPECT_TRUE(MetricInfoFor(MetricId::MorrisThorne).gpu_supported);

    // Charge and lambda are not plumbed to the kernel yet.
    EXPECT_FALSE(MetricInfoFor(MetricId::ReissnerNordstrom).gpu_supported);
    EXPECT_FALSE(MetricInfoFor(MetricId::KerrNewman).gpu_supported);
    EXPECT_FALSE(MetricInfoFor(MetricId::DeSitter).gpu_supported);
    EXPECT_FALSE(MetricInfoFor(MetricId::SchwarzschildDeSitter).gpu_supported);

    EXPECT_TRUE(MetricInfoFor(MetricId::Schwarzschild).cpu_supported);
    EXPECT_TRUE(MetricInfoFor(MetricId::Schwarzschild).gpu_supported);
    EXPECT_TRUE(MetricInfoFor(MetricId::Alcubierre).cpu_supported);
}
