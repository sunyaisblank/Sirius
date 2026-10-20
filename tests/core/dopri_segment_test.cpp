#include "sirius/core/dopri_segment.h"

#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <limits>

namespace sirius::test {
namespace {

using namespace sirius::core;

// Supply an independently specified spatial power polynomial to the Hairer
// representation. The constant and endpoint increment are exact dyadics in
// these witnesses, so no fitted endpoint or numerical reference is involved.
DopriPositionSegment PolynomialSegment(const std::array<double, 5>& coefficients, int axis = 1,
                                       double interval = 1.0) {
    DopriPositionSegment segment;
    segment.interval = interval;
    segment.origin(axis) = coefficients[0];
    for (int power = 1; power <= 4; ++power) segment.increment(axis) += coefficients[power];
    segment.a(axis) = coefficients[1] - segment.increment(axis);
    segment.c(axis) = coefficients[4];
    segment.b(axis) = -coefficients[3] - 2.0 * coefficients[4];
    return segment;
}

TEST(DopriPositionSegment, SamplesTheQuarticAndItsAffineDerivative) {
    // x(s) = 3 + 2s - 3s^2 + 4s^3 - s^4, affine interval two.
    auto segment = PolynomialSegment({3.0, 2.0, -3.0, 4.0, -1.0}, 1, 2.0);
    segment.origin(0) = 0x1p50;
    segment.increment(0) = 0x1p-20;
    for (double fraction : {0.0, 0.125, 0.5, 0.875, 1.0}) {
        const auto sample = segment.Sample(fraction);
        const double position =
            3.0 + fraction * (2.0 + fraction * (-3.0 + fraction * (4.0 - fraction)));
        const double tangent = (2.0 + fraction * (-6.0 + fraction * (12.0 - 4.0 * fraction))) / 2.0;
        EXPECT_DOUBLE_EQ(sample.position(1), position);
        EXPECT_DOUBLE_EQ(sample.tangent(1), tangent);
        EXPECT_DOUBLE_EQ(sample.fraction, fraction);
        EXPECT_DOUBLE_EQ(segment.Displacement(fraction)(0), fraction * 0x1p-20);
        EXPECT_DOUBLE_EQ(sample.tangent(0), 0x1p-21);
    }
    EXPECT_DOUBLE_EQ(segment.Displacement(1.0)(1), segment.increment(1));
}

TEST(DopriPositionSegment, DegreeEightIsolationCoversEveryRootWithoutGridSampling) {
    // The eight simple roots are independently prescribed at i/8, i=0..7.
    std::array<double, 9> coefficients{};
    coefficients[0] = 1.0;
    int degree = 0;
    for (int index = 0; index < 8; ++index) {
        const double root = static_cast<double>(index) / 8.0;
        for (int power = degree + 1; power >= 1; --power) {
            coefficients[power] = coefficients[power - 1] - root * coefficients[power];
        }
        coefficients[0] *= -root;
        ++degree;
    }
    const auto roots = core::detail::FindPolynomialRootsOnUnitInterval(coefficients, degree);
    ASSERT_EQ(roots.count, 8);
    for (int index = 0; index < roots.count; ++index) {
        EXPECT_NEAR(roots.values[index], static_cast<double>(index) / 8.0, 5e-10);
    }
}

TEST(DopriPositionSegment, SameSideQuarticContactsRespectDirectionAndOriginalRestriction) {
    // x(s) = 1 + 8(s-1/4)(s-3/4)(1+s^2). Both endpoints are outside.
    const auto segment = PolynomialSegment({2.5, -8.0, 9.5, -8.0, 8.0});
    EXPECT_GT(segment.Sample(0.0).position(1), 1.0);
    EXPECT_GT(segment.Sample(1.0).position(1), 1.0);
    const auto inward =
        FindSphericalBoundaryEvent(segment, 1.0, SphericalBoundarySense::DecreasingRadius);
    const auto outward =
        FindSphericalBoundaryEvent(segment, 1.0, SphericalBoundarySense::IncreasingRadius);
    ASSERT_TRUE(inward);
    ASSERT_TRUE(outward);
    EXPECT_NEAR(inward->fraction, 0.25, 2e-12);
    EXPECT_NEAR(outward->fraction, 0.75, 2e-12);
    EXPECT_NEAR(inward->tangent(1), -4.25, 2e-11);
    EXPECT_NEAR(outward->tangent(1), 6.25, 2e-11);

    const auto restricted = segment.Restricted(0.5);
    EXPECT_DOUBLE_EQ(restricted.increment(1), segment.increment(1));
    EXPECT_DOUBLE_EQ(restricted.interval, segment.interval);
    EXPECT_DOUBLE_EQ(restricted.Sample(0.25).position(1), segment.Sample(0.25).position(1));
    const auto restricted_inward =
        FindSphericalBoundaryEvent(restricted, 1.0, SphericalBoundarySense::DecreasingRadius);
    ASSERT_TRUE(restricted_inward);
    EXPECT_NEAR(restricted_inward->fraction, 0.25, 2e-12);
    EXPECT_FALSE(
        FindSphericalBoundaryEvent(restricted, 1.0, SphericalBoundarySense::IncreasingRadius));
    const auto retained_outward = FindSphericalBoundaryEvent(
        segment.Restricted(0.875), 1.0, SphericalBoundarySense::IncreasingRadius);
    ASSERT_TRUE(retained_outward);
    EXPECT_NEAR(retained_outward->fraction, 0.75, 2e-12);
}

TEST(DopriPositionSegment, TangenciesAndKerrEllipsoidUseTheSameQuartic) {
    // x(s) = 1 + (s-1/2)^2(1+s^2): an isolated tangent contact.
    auto tangent = PolynomialSegment({1.25, -1.0, 1.25, -1.0, 1.0});
    const auto contact =
        FindSphericalBoundaryEvent(tangent, 1.0, SphericalBoundarySense::AnyContact);
    ASSERT_TRUE(contact);
    EXPECT_NEAR(contact->fraction, 0.5, 2e-12);
    EXPECT_NEAR(contact->tangent(1), 0.0, 2e-12);
    EXPECT_FALSE(
        FindSphericalBoundaryEvent(tangent, 1.0, SphericalBoundarySense::IncreasingRadius));
    EXPECT_FALSE(
        FindSphericalBoundaryEvent(tangent, 1.0, SphericalBoundarySense::DecreasingRadius));

    // r=4,a=3 gives equatorial semiaxis five, rather than a sphere of radius four.
    auto oblate = PolynomialSegment({12.5, -40.0, 47.5, -40.0, 40.0});
    const auto entry =
        FindKerrEllipsoidBoundaryEvent(oblate, 4.0, 3.0, SphericalBoundarySense::DecreasingRadius);
    ASSERT_TRUE(entry);
    EXPECT_NEAR(entry->fraction, 0.25, 2e-12);
    EXPECT_NEAR(entry->position(1), 5.0, 2e-11);
    const auto spherical =
        FindSphericalBoundaryEvent(oblate, 5.0, SphericalBoundarySense::DecreasingRadius);
    const auto zero_spin =
        FindKerrEllipsoidBoundaryEvent(oblate, 5.0, 0.0, SphericalBoundarySense::DecreasingRadius);
    ASSERT_TRUE(spherical);
    ASSERT_TRUE(zero_spin);
    EXPECT_DOUBLE_EQ(zero_spin->fraction, spherical->fraction);
    const auto polar =
        FindKerrEllipsoidBoundaryEvent(PolynomialSegment({6.0, -4.0, 0.0, 0.0, 0.0}, 3), 4.0, 3.0,
                                       SphericalBoundarySense::DecreasingRadius);
    ASSERT_TRUE(polar);
    EXPECT_NEAR(polar->fraction, 0.5, 2e-12);
    EXPECT_NEAR(polar->position(3), 4.0, 2e-12);
}

TEST(DopriPositionSegment, DiskRootsIncludeTangenciesAndKeepOriginalFractions) {
    // z(s)=(s-1/8)(s-3/8)(s-5/8)(s-7/8).
    const auto disk = PolynomialSegment({105.0 / 4096.0, -11.0 / 32.0, 43.0 / 32.0, -2.0, 1.0}, 3);
    const auto roots = FindDiskPlaneRoots(disk);
    ASSERT_EQ(roots.count, 4);
    for (int index = 0; index < roots.count; ++index) {
        EXPECT_NEAR(roots.values[index], static_cast<double>(2 * index + 1) / 8.0, 2e-12);
    }
    const auto restricted = FindDiskPlaneRoots(disk.Restricted(0.5));
    ASSERT_EQ(restricted.count, 2);
    EXPECT_NEAR(restricted.values[0], 0.125, 2e-12);
    EXPECT_NEAR(restricted.values[1], 0.375, 2e-12);
    const auto closed_limit = FindDiskPlaneRoots(disk.Restricted(0.375));
    ASSERT_EQ(closed_limit.count, 2);
    EXPECT_NEAR(closed_limit.values[1], 0.375, 2e-12);
    const auto tangent = FindDiskPlaneRoots(PolynomialSegment({0.25, -1.0, 1.0, 0.0, 0.0}, 3));
    ASSERT_EQ(tangent.count, 1);
    EXPECT_NEAR(tangent.values[0], 0.5, 2e-12);
    EXPECT_EQ(FindDiskPlaneRoots(PolynomialSegment({0, 0, 0, 0, 0}, 3)).count, 0);
}

TEST(DopriPositionSegment, IncrementBoundarySurvivesLargeOriginsAndInvalidInputsFailClosed) {
    // A displacement of one reaches a radius2^40 from radius2^40-1. The
    // additional quartic vanishes at s=0,1/2,1 and does not refit the endpoint.
    constexpr double radius = 0x1p40;
    auto segment = PolynomialSegment({radius - 1.0, 2.25, -1.25, 2.0, -1.0});
    const auto event =
        FindSphericalBoundaryEvent(segment, radius, SphericalBoundarySense::IncreasingRadius);
    ASSERT_TRUE(event);
    EXPECT_NEAR(event->fraction, 0.5, 2e-12);
    EXPECT_DOUBLE_EQ(segment.Displacement(0.5)(1), 1.0);
    EXPECT_DOUBLE_EQ(event->position(1), radius);
    EXPECT_FALSE(FindSphericalBoundaryEvent(segment.Restricted(0.25), radius,
                                            SphericalBoundarySense::AnyContact));
    EXPECT_FALSE(FindSphericalBoundaryEvent(segment, -radius, SphericalBoundarySense::AnyContact));
    EXPECT_FALSE(FindKerrEllipsoidBoundaryEvent(segment, radius,
                                                std::numeric_limits<double>::infinity(),
                                                SphericalBoundarySense::AnyContact));
    segment.a(3) = std::numeric_limits<double>::quiet_NaN();
    EXPECT_FALSE(FindSphericalBoundaryEvent(segment, radius, SphericalBoundarySense::AnyContact));
    EXPECT_EQ(FindDiskPlaneRoots(segment).count, 0);
}

}  // namespace
}  // namespace sirius::test
