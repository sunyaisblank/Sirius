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
#include <PHMT200A.h>
#include <PHMT100A.h>
#include <PHMT101A.h>
#include <PHMT102A.h>

using namespace Sirius;

TEST(MetricRegistryTests, EveryCanonicalNameParsesToItsOwnId) {
    for (const auto& info : metricRegistry()) {
        auto id = parseMetricName(info.canonicalName);
        ASSERT_TRUE(id.has_value()) << "canonical name failed: " << info.canonicalName;
        EXPECT_EQ(*id, info.id) << "canonical name misrouted: " << info.canonicalName;
    }
}

TEST(MetricRegistryTests, EveryAliasParsesToItsOwnId) {
    for (const auto& info : metricRegistry()) {
        for (const char* alias : info.aliases) {
            if (alias == nullptr) continue;
            auto id = parseMetricName(alias);
            ASSERT_TRUE(id.has_value()) << "alias failed: " << alias;
            EXPECT_EQ(*id, info.id) << "alias misrouted: " << alias;
        }
    }
}

TEST(MetricRegistryTests, ParsingIsCaseInsensitive) {
    EXPECT_EQ(parseMetricName("schwarzschild"), MetricId::Schwarzschild);
    EXPECT_EQ(parseMetricName("KERR"), MetricId::Kerr);
    EXPECT_EQ(parseMetricName("morristhorne"), MetricId::MorrisThorne);
    EXPECT_EQ(parseMetricName("wormhole"), MetricId::MorrisThorne);
    EXPECT_EQ(parseMetricName("warpdrive"), MetricId::Alcubierre);
}

TEST(MetricRegistryTests, UnknownNamesFailInsteadOfDefaulting) {
    EXPECT_FALSE(parseMetricName("").has_value());
    EXPECT_FALSE(parseMetricName("Schwarzchild").has_value());   // common misspelling
    EXPECT_FALSE(parseMetricName("Godel").has_value());          // unimplemented on purpose
    EXPECT_FALSE(parseMetricName("Taub-NUT").has_value());
    EXPECT_FALSE(parseMetricName("Kerr ").has_value());          // trailing space is not a name
}

TEST(MetricRegistryTests, MetricInfoRoundTripsById) {
    for (const auto& info : metricRegistry()) {
        EXPECT_EQ(metricInfo(info.id).id, info.id);
        EXPECT_STREQ(metricInfo(info.id).canonicalName, info.canonicalName);
    }
}

TEST(MetricRegistryTests, KnownMetricNamesListsEveryCanonicalName) {
    std::string names = knownMetricNames();
    for (const auto& info : metricRegistry()) {
        EXPECT_NE(names.find(info.canonicalName), std::string::npos)
            << "missing from error text: " << info.canonicalName;
    }
}

// Display names produced by the metric families must resolve through the
// registry (the GL pipeline routes IMetric::getName() output through it).
TEST(MetricRegistryTests, FamilyDisplayNamesRoundTrip) {
    KerrSchildFamily schw(KerrSchildParams::Schwarzschild(1.0));
    EXPECT_EQ(parseMetricName(schw.getName()), MetricId::Schwarzschild);

    KerrSchildFamily kerr(KerrSchildParams::Kerr(1.0, 0.9));
    EXPECT_EQ(parseMetricName(kerr.getName()), MetricId::Kerr);

    KerrSchildFamily rn(KerrSchildParams::ReissnerNordstrom(1.0, 0.5));
    EXPECT_EQ(parseMetricName(rn.getName()), MetricId::ReissnerNordstrom);

    KerrSchildFamily kn(KerrSchildParams::KerrNewman(1.0, 0.5, 0.3));
    EXPECT_EQ(parseMetricName(kn.getName()), MetricId::KerrNewman);

    KerrSchildFamily mink(KerrSchildParams::Minkowski());
    EXPECT_EQ(parseMetricName(mink.getName()), MetricId::Minkowski);

    MorrisThorneFamily wormhole(MorrisThorneParams::Ellis(1.0));
    EXPECT_EQ(parseMetricName(wormhole.getName()), MetricId::MorrisThorne);

    WarpDriveFamily warp(WarpDriveParams::Alcubierre(1.0, 1.0));
    EXPECT_EQ(parseMetricName(warp.getName()), MetricId::Alcubierre);
}

// Backend support flags document reality; these pin the deliberate ones so a
// flag change is a conscious decision, not a drive-by edit.
TEST(MetricRegistryTests, BackendSupportMatchesImplementations) {
    // Morris-Thorne evaluates in spherical coordinates; the Cartesian CPU
    // tracer cannot drive it (see PHMT200A).
    EXPECT_FALSE(metricInfo(MetricId::MorrisThorne).cpuSupported);
    EXPECT_TRUE(metricInfo(MetricId::MorrisThorne).gpuSupported);

    // Charge and lambda are not plumbed to the kernel yet.
    EXPECT_FALSE(metricInfo(MetricId::ReissnerNordstrom).gpuSupported);
    EXPECT_FALSE(metricInfo(MetricId::KerrNewman).gpuSupported);
    EXPECT_FALSE(metricInfo(MetricId::DeSitter).gpuSupported);
    EXPECT_FALSE(metricInfo(MetricId::SchwarzschildDeSitter).gpuSupported);

    EXPECT_TRUE(metricInfo(MetricId::Schwarzschild).cpuSupported);
    EXPECT_TRUE(metricInfo(MetricId::Schwarzschild).gpuSupported);
    EXPECT_TRUE(metricInfo(MetricId::Alcubierre).cpuSupported);
}
