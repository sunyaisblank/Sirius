// Independent spectral-physics oracles for the production blackbody authority.

#include "sirius/core/constants.h"
#include "sirius/core/relativistic_transfer.h"
#include "sirius/core/spectral/blackbody.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numbers>

namespace sirius::test {
namespace {

using sirius::core::spectral::BlackbodyToRgb;
using sirius::core::spectral::DopplerFactor;
using sirius::core::spectral::Rgb;
using sirius::core::spectral::TotalRedshift;
using sirius::core::spectral::TryPlanckSpectralRadiancePerMetre;
using sirius::core::spectral::TryStefanBoltzmannExitance;
using sirius::core::spectral::TryWienPeakWavelength;

constexpr double kReferencePlanck = 6.62607015e-34;
constexpr double kReferenceLightSpeed = 299792458.0;
constexpr double kReferenceBoltzmann = 1.380649e-23;

double IndependentPlanckRadiance(double wavelength, double temperature) {
    const double numerator = 2.0 * kReferencePlanck * kReferenceLightSpeed * kReferenceLightSpeed;
    const double exponent =
        kReferencePlanck * kReferenceLightSpeed / (wavelength * kReferenceBoltzmann * temperature);
    return numerator / (std::pow(wavelength, 5) * std::expm1(exponent));
}

double RelativeError(double actual, double expected) {
    return std::abs(actual - expected) / std::max(std::abs(expected), 1.0e-300);
}

double IntegratedExitance(double temperature) {
    // Midpoint quadrature in log(lambda) spans all material thermal radiance.
    // B_lambda d(lambda) becomes B_exp(x) exp(x) dx.
    constexpr int kPanels = 32768;
    constexpr double kLogMinimum = -20.72326583694641;  // log(1 nm)
    constexpr double kLogMaximum = -6.907755278982137;  // log(1 mm)
    constexpr double kStep = (kLogMaximum - kLogMinimum) / kPanels;
    double integral = 0.0;
    for (int panel = 0; panel < kPanels; ++panel) {
        const double x = kLogMinimum + (panel + 0.5) * kStep;
        const double wavelength = std::exp(x);
        integral += IndependentPlanckRadiance(wavelength, temperature) * wavelength;
    }
    return std::numbers::pi * integral * kStep;
}

double NumericalPeakWavelength(double temperature) {
    // Golden-section maximisation is independent of the closed-form Wien law.
    constexpr double kPhi = 1.6180339887498948482;
    double lower = 50.0e-9;
    double upper = 20.0e-6;
    double left = upper - (upper - lower) / kPhi;
    double right = lower + (upper - lower) / kPhi;
    for (int iteration = 0; iteration < 120; ++iteration) {
        if (IndependentPlanckRadiance(left, temperature) <
            IndependentPlanckRadiance(right, temperature)) {
            lower = left;
            left = right;
            right = lower + (upper - lower) / kPhi;
        } else {
            upper = right;
            right = left;
            left = upper - (upper - lower) / kPhi;
        }
    }
    return 0.5 * (lower + upper);
}

}  // namespace

TEST(SpectralValidationTests, PlanckSpectralRadianceMatchesIndependentCodataEquation) {
    for (const double temperature : {1200.0, 3000.0, 5778.0, 10000.0, 1.0e7}) {
        for (const double wavelength : {100.0e-9, 380.0e-9, 550.0e-9, 780.0e-9, 10.0e-6}) {
            SCOPED_TRACE(::testing::Message() << "T=" << temperature << " lambda=" << wavelength);
            const double expected = IndependentPlanckRadiance(wavelength, temperature);
            const auto actual = TryPlanckSpectralRadiancePerMetre(wavelength, temperature);
            ASSERT_TRUE(actual.has_value());
            EXPECT_LT(RelativeError(*actual, expected), 2.0e-14);
        }
    }
}

TEST(SpectralValidationTests, WienPeakMatchesIndependentNumericalPlanckMaximum) {
    for (const double temperature : {1200.0, 3000.0, 5778.0, 10000.0}) {
        SCOPED_TRACE(temperature);
        const auto analytic = TryWienPeakWavelength(temperature);
        ASSERT_TRUE(analytic.has_value());
        const double numerical = NumericalPeakWavelength(temperature);
        EXPECT_LT(RelativeError(*analytic, numerical), 2.0e-8);
    }
}

TEST(SpectralValidationTests, StefanBoltzmannExitanceMatchesIndependentHemisphericPlanckIntegral) {
    for (const double temperature : {1200.0, 3000.0, 5778.0, 10000.0}) {
        SCOPED_TRACE(temperature);
        const auto analytic = TryStefanBoltzmannExitance(temperature);
        ASSERT_TRUE(analytic.has_value());
        const double quadrature = IntegratedExitance(temperature);
        EXPECT_LT(RelativeError(*analytic, quadrature), 2.0e-7);
    }
}

TEST(SpectralValidationTests, BlackbodyLawsRejectUnrepresentedDomains) {
    const double infinity = std::numeric_limits<double>::infinity();
    const double nan = std::numeric_limits<double>::quiet_NaN();
    for (const double invalid : {0.0, -1.0, infinity, nan}) {
        SCOPED_TRACE(invalid);
        EXPECT_FALSE(TryPlanckSpectralRadiancePerMetre(invalid, 5778.0).has_value());
        EXPECT_FALSE(TryPlanckSpectralRadiancePerMetre(550.0e-9, invalid).has_value());
        EXPECT_FALSE(TryWienPeakWavelength(invalid).has_value());
        EXPECT_FALSE(TryStefanBoltzmannExitance(invalid).has_value());
    }

    const auto deep_wien_tail = TryPlanckSpectralRadiancePerMetre(1.0e-20, 1.0);
    ASSERT_TRUE(deep_wien_tail.has_value());
    EXPECT_DOUBLE_EQ(*deep_wien_tail, 0.0);

    // The defining expression has an infinity-times-zero intermediate here,
    // but its Rayleigh-Jeans limit is finite and represented.
    const auto extreme_rayleigh_jeans = TryPlanckSpectralRadiancePerMetre(1.0e100, 1.0e300);
    ASSERT_TRUE(extreme_rayleigh_jeans.has_value());
    const double rayleigh_jeans = sirius::core::constants::physical::kPlanckC1 /
                                  sirius::core::constants::physical::kPlanckC2 * 1.0e-100;
    EXPECT_LT(RelativeError(*extreme_rayleigh_jeans, rayleigh_jeans), 2.0e-13);

    // This positive finite input has a mathematical radiance above DBL_MAX;
    // zero would confuse overflow with a physically underflowed Wien tail.
    const double overflow_temperature = sirius::core::constants::physical::kPlanckC2 / 1.0e-297;
    EXPECT_FALSE(TryPlanckSpectralRadiancePerMetre(1.0e-300, overflow_temperature).has_value());
    EXPECT_FALSE(TryStefanBoltzmannExitance(std::numeric_limits<double>::max()).has_value());
}

TEST(SpectralValidationTests, DopplerFactorMatchesIndependentLorentzFormula) {
    for (const double beta : {-0.95, -0.5, 0.0, 0.5, 0.95}) {
        SCOPED_TRACE(beta);
        const double expected = std::sqrt((1.0 - beta) / (1.0 + beta));
        EXPECT_NEAR(DopplerFactor(beta), expected, 2.0e-15);
    }
    EXPECT_TRUE(std::isnan(DopplerFactor(-1.0)));
    EXPECT_TRUE(std::isnan(DopplerFactor(1.0)));
    EXPECT_TRUE(std::isnan(DopplerFactor(0.5, 0.0)));
    EXPECT_TRUE(std::isnan(DopplerFactor(0.5, -1.0)));
}

TEST(SpectralValidationTests, TotalRedshiftComposesGravitationalAndDopplerFactors) {
    constexpr double kEmitGtt = -0.42;
    constexpr double kObserverGtt = -0.98;
    constexpr double kBeta = 0.37;
    const double gravitational = std::sqrt(kEmitGtt / kObserverGtt);
    const double doppler = std::sqrt((1.0 - kBeta) / (1.0 + kBeta));
    const double expected = 1.0 / (gravitational * doppler) - 1.0;
    EXPECT_NEAR(TotalRedshift(kEmitGtt, kObserverGtt, kBeta), expected, 2.0e-15);
}

TEST(SpectralValidationTests, StaticMetricRedshiftHasThePhysicalDirection) {
    constexpr double kEmitRadius = 6.0;
    constexpr double kObserverRadius = 100.0;
    constexpr double kMass = 1.0;
    const double emit_gtt = -(1.0 - 2.0 * kMass / kEmitRadius);
    const double observer_gtt = -(1.0 - 2.0 * kMass / kObserverRadius);
    const double expected_g = std::sqrt(emit_gtt / observer_gtt);
    const double redshift = TotalRedshift(emit_gtt, observer_gtt, 0.0);
    EXPECT_GT(redshift, 0.0);
    EXPECT_NEAR(1.0 / (1.0 + redshift), expected_g, 2.0e-15);
    EXPECT_TRUE(std::isnan(TotalRedshift(0.1, observer_gtt, 0.0)));
    EXPECT_TRUE(std::isnan(TotalRedshift(emit_gtt, 0.0, 0.0)));
}

TEST(SpectralValidationTests, KerrDiskTransferMatchesKillingFieldContraction) {
    constexpr double kMass = 1.0;
    constexpr double kSpin = 0.7;
    constexpr double kRadius = 6.0;
    constexpr double kEnergy = 1.0;
    constexpr double kAngularMomentum = 2.0;
    constexpr double kObserverFrequency = 0.93;
    const auto transfer = sirius::core::relativity::KerrDiskTransfer(
        kObserverFrequency, kEnergy, kAngularMomentum, kMass, kSpin, kRadius);
    ASSERT_TRUE(transfer.has_value());

    const double g_tt = -(1.0 - 2.0 * kMass / kRadius);
    const double g_t_phi = -2.0 * kMass * kSpin / kRadius;
    const double g_phi_phi =
        kRadius * kRadius + kSpin * kSpin + 2.0 * kMass * kSpin * kSpin / kRadius;
    const double omega = std::sqrt(kMass) / (std::pow(kRadius, 1.5) + kSpin * std::sqrt(kMass));
    const double u_t = 1.0 / std::sqrt(-(g_tt + 2.0 * omega * g_t_phi + omega * omega * g_phi_phi));
    const double expected_emitter_frequency = u_t * (kEnergy - omega * kAngularMomentum);
    EXPECT_NEAR(transfer->emitter_frequency, expected_emitter_frequency, 2.0e-15);
    EXPECT_NEAR(transfer->full_g, kObserverFrequency / expected_emitter_frequency, 2.0e-15);

    const double zamo_omega = -g_t_phi / g_phi_phi;
    const double zamo_u_t =
        1.0 / std::sqrt(-(g_tt + 2.0 * zamo_omega * g_t_phi + zamo_omega * zamo_omega * g_phi_phi));
    const double expected_zamo_frequency = zamo_u_t * (kEnergy - zamo_omega * kAngularMomentum);
    EXPECT_NEAR(transfer->zamo_frequency, expected_zamo_frequency, 2.0e-15);
    EXPECT_NEAR(transfer->zamo_g, kObserverFrequency / expected_zamo_frequency, 2.0e-15);
}

TEST(SpectralValidationTests, ZamoBranchRemainsTimelikeInsideTheErgosphere) {
    constexpr double kRadius = 1.5;
    const auto static_ratio =
        sirius::core::relativity::StaticObserverFrequencyRatio(1.0 / 3.0, -1.0);
    EXPECT_FALSE(static_ratio.has_value());
    const auto transfer =
        sirius::core::relativity::KerrDiskTransfer(1.0, 1.0, 1.5, 1.0, 0.998, kRadius);
    ASSERT_TRUE(transfer.has_value());
    EXPECT_GT(transfer->zamo_frequency, 0.0);
    EXPECT_TRUE(std::isfinite(transfer->zamo_g));
}

TEST(SpectralValidationTests, ComovingOpacityUsesInvariantAffinePathLength) {
    sirius::core::Metric4d metric;
    metric(0, 0) = -1.0;
    metric(1, 1) = 1.0;
    metric(2, 2) = 1.0;
    metric(3, 3) = 1.0;
    sirius::core::Vec4 past_ray;
    past_ray(0) = -1.0;
    past_ray(1) = 1.0;
    constexpr double kBeta = 0.6;
    const double gamma = 1.0 / std::sqrt(1.0 - kBeta * kBeta);
    sirius::core::Vec4 fluid;
    fluid(0) = gamma;
    fluid(1) = gamma * kBeta;

    const auto path = sirius::core::relativity::ComovingPathLength(past_ray, fluid, metric, 2.5);
    ASSERT_TRUE(path.has_value());
    EXPECT_NEAR(*path, gamma * (1.0 + kBeta) * 2.5, 2.0e-15);
}

TEST(SpectralValidationTests, ObserverToSourceTransferPreservesForegroundEmissionOrder) {
    sirius::core::relativity::GreyTransferState state;
    const double layer_tau = std::log(2.0);
    ASSERT_TRUE(sirius::core::relativity::AccumulateObserverToSourceLayer(state, {1.0, 0.0, 0.0},
                                                                          layer_tau, 10.0)
                    .has_value());
    ASSERT_TRUE(sirius::core::relativity::AccumulateObserverToSourceLayer(state, {0.0, 0.0, 1.0},
                                                                          layer_tau, 10.0)
                    .has_value());

    EXPECT_NEAR(state.observed_emission[0], 0.5, 2.0e-15);
    EXPECT_NEAR(state.observed_emission[1], 0.0, 2.0e-15);
    EXPECT_NEAR(state.observed_emission[2], 0.25, 2.0e-15);
    EXPECT_NEAR(std::exp(-state.optical_depth), 0.25, 2.0e-15);
}

TEST(SpectralValidationTests, BlackbodyColourProgressionConsumesIntegratedSpectrum) {
    const Rgb very_cold = BlackbodyToRgb(500.0);
    const Rgb former_cold_boundary = BlackbodyToRgb(1000.0);
    const Rgb cold = BlackbodyToRgb(2000.0);
    const Rgb solar = BlackbodyToRgb(5500.0);
    const Rgb hot = BlackbodyToRgb(10000.0);
    const Rgb former_hot_boundary = BlackbodyToRgb(100000.0);
    const Rgb very_hot = BlackbodyToRgb(1000000.0);
    const auto blue_to_red = [](const Rgb& colour) {
        return colour.b / std::max(colour.r, 1.0e-6f);
    };
    EXPECT_GT(cold.r, cold.b);
    EXPECT_GT(hot.b, hot.r);
    EXPECT_LT(blue_to_red(cold), blue_to_red(solar));
    EXPECT_LT(blue_to_red(solar), blue_to_red(hot));
    EXPECT_GT(std::abs(very_cold.r - former_cold_boundary.r) +
                  std::abs(very_cold.g - former_cold_boundary.g) +
                  std::abs(very_cold.b - former_cold_boundary.b),
              1.0e-3f);
    EXPECT_GT(std::abs(very_hot.r - former_hot_boundary.r) +
                  std::abs(very_hot.g - former_hot_boundary.g) +
                  std::abs(very_hot.b - former_hot_boundary.b),
              1.0e-3f);
}

}  // namespace sirius::test
