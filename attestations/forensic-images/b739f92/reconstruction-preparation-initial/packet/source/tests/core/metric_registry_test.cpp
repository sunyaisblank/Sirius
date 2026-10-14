// TSMT200A.cpp - Metric Identity Registry Tests
// Component ID: TSMT200A
// Tests for: PHMT200A.h
//
// The registry is the single naming authority; these tests pin the
// postconditions that ended the stringly-typed routing defects: every
// canonical name and alias parses to its identity, display names from the
// metric families round-trip, unknown names fail rather than defaulting,
// and the error text can enumerate the accepted set.

#include <gtest/gtest.h>
#include "sirius/core/metrics/registry.h"
#include "sirius/core/metrics/kerr_schild_family.h"
#include "sirius/core/metrics/morris_thorne_family.h"
#include "sirius/core/metrics/warp_drive_family.h"

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
    EXPECT_FALSE(ParseMetricName("Schwarzchild").has_value());   // common misspelling
    EXPECT_FALSE(ParseMetricName("Godel").has_value());          // unimplemented on purpose
    EXPECT_FALSE(ParseMetricName("Taub-NUT").has_value());
    EXPECT_FALSE(ParseMetricName("Kerr ").has_value());          // trailing space is not a name
}

TEST(MetricRegistryTests, MetricInfoRoundTripsById) {
    for (const auto& info : MetricRegistry()) {
        EXPECT_EQ(MetricInfoFor(info.id).id, info.id);
        EXPECT_STREQ(MetricInfoFor(info.id).canonical_name, info.canonical_name);
    }
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
}

// Backend support flags document reality; these pin the deliberate ones so a
// flag change is a conscious decision, not a drive-by edit.
TEST(MetricRegistryTests, BackendSupportMatchesImplementations) {
    // Morris-Thorne evaluates in spherical coordinates; the Cartesian CPU
    // tracer cannot drive it (see PHMT200A).
    EXPECT_FALSE(MetricInfoFor(MetricId::MorrisThorne).cpu_supported);
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
