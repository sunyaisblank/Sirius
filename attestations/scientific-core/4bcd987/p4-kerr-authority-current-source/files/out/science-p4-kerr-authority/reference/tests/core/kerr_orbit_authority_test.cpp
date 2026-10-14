// Cross-consumer controls for the exact CPU Kerr circular-orbit authority.
// The Boyer-Lindquist oracle and Slang implementation intentionally remain
// separate mathematical comparators; this suite prevents production CPU
// consumers from growing mutually inconsistent ISCO or emitter laws.

#include "sirius/core/disk/novikov_thorne_disk.h"
#include "sirius/core/kerr_orbits.h"
#include "sirius/core/metrics/kerr_schild_family.h"
#include "sirius/core/relativistic_transfer.h"

#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <limits>

namespace {

using sirius::core::AccretionDiskD;
using sirius::core::KerrSchildFamily;
using sirius::core::KerrSchildParams;
using sirius::core::relativity::KerrDiskTransfer;
using sirius::core::relativity::TryKerrCircularOrbitAngularVelocity;
using sirius::core::relativity::TryKerrIscoRadius;

TEST(KerrOrbitAuthority, EveryCpuConsumerSharesSignedIscoAndCircularEmitterLaw) {
    struct Case {
        double mass;
        double dimensionless_spin;
        double orbit_radius_in_mass_units;
    };
    constexpr std::array cases{
        Case{0.5, -0.998, 12.0}, Case{2.0, -0.5, 12.0},   Case{3.0, -0.0009, 12.0},
        Case{0.5, 0.0, 12.0},    Case{2.0, 0.0009, 12.0}, Case{3.0, 0.5, 12.0},
        Case{2.0, 0.998, 12.0},
    };

    for (const auto& sample : cases) {
        const double spin = sample.dimensionless_spin * sample.mass;
        const double radius = sample.orbit_radius_in_mass_units * sample.mass;
        const auto isco = TryKerrIscoRadius(sample.mass, spin);
        const auto omega = TryKerrCircularOrbitAngularVelocity(sample.mass, spin, radius);
        ASSERT_TRUE(isco.has_value());
        ASSERT_TRUE(omega.has_value());

        KerrSchildFamily metric(KerrSchildParams::Kerr(sample.mass, spin));
        EXPECT_DOUBLE_EQ(metric.IscoRadius(), *isco);
        EXPECT_DOUBLE_EQ(sample.mass * AccretionDiskD::ComputeIsco(sample.dimensionless_spin),
                         *isco);
        EXPECT_DOUBLE_EQ(AccretionDiskD::AngularVelocity(radius, sample.mass, spin), *omega);

        const auto transfer = KerrDiskTransfer(1.0, 1.0, 0.0, sample.mass, spin, radius);
        ASSERT_TRUE(transfer.has_value());
        EXPECT_DOUBLE_EQ(transfer->emitter.angular_velocity, *omega);
    }

    EXPECT_DOUBLE_EQ(*TryKerrIscoRadius(2.0, 0.0), 12.0);
    EXPECT_DOUBLE_EQ(*TryKerrIscoRadius(2.0, 2.0), 2.0);
    EXPECT_DOUBLE_EQ(*TryKerrIscoRadius(2.0, -2.0), 18.0);
}

TEST(KerrOrbitAuthority, UnrepresentedInputsDeclineInsteadOfChangingTheOrbit) {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double infinity = std::numeric_limits<double>::infinity();

    EXPECT_FALSE(TryKerrIscoRadius(0.0, 0.0).has_value());
    EXPECT_FALSE(TryKerrIscoRadius(-1.0, 0.0).has_value());
    EXPECT_FALSE(TryKerrIscoRadius(1.0, 1.0001).has_value());
    EXPECT_FALSE(TryKerrIscoRadius(1.0, nan).has_value());
    EXPECT_FALSE(TryKerrIscoRadius(infinity, 0.0).has_value());

    EXPECT_FALSE(TryKerrCircularOrbitAngularVelocity(0.0, 0.0, 6.0).has_value());
    EXPECT_FALSE(TryKerrCircularOrbitAngularVelocity(1.0, 1.0001, 6.0).has_value());
    EXPECT_FALSE(TryKerrCircularOrbitAngularVelocity(1.0, 0.0, 0.0).has_value());
    EXPECT_FALSE(TryKerrCircularOrbitAngularVelocity(1.0, 0.0, nan).has_value());
    EXPECT_FALSE(TryKerrCircularOrbitAngularVelocity(1.0, -1.0, 1.0).has_value());

    EXPECT_DEATH((void)AccretionDiskD::ComputeIsco(nan), "precondition.*enforced, terminating");
    EXPECT_DEATH((void)AccretionDiskD::AngularVelocity(1.0, 1.0, -1.0),
                 "precondition.*enforced, terminating");
}

}  // namespace
