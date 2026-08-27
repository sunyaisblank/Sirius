// Independent spectral-physics oracles for the production blackbody authority.

#include "sirius/core/spectral/blackbody.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <numbers>

namespace sirius::test {
namespace {

using sirius::core::spectral::ApplyLimbDarkening;
using sirius::core::spectral::BlackbodyToRgb;
using sirius::core::spectral::DopplerFactor;
using sirius::core::spectral::LimbDarkeningCoeff;
using sirius::core::spectral::PlanckRadiance;
using sirius::core::spectral::Rgb;
using sirius::core::spectral::StefanBoltzmannRadiance;
using sirius::core::spectral::TotalRedshift;
using sirius::core::spectral::WienPeakWavelength;

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
        if (PlanckRadiance(left, temperature) < PlanckRadiance(right, temperature)) {
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

TEST(SpectralValidationTests, PlanckRadianceMatchesIndependentCodataEquation) {
    for (const double temperature : {1200.0, 3000.0, 5778.0, 10000.0, 1.0e7}) {
        for (const double wavelength : {100.0e-9, 380.0e-9, 550.0e-9, 780.0e-9, 10.0e-6}) {
            SCOPED_TRACE(::testing::Message() << "T=" << temperature << " lambda=" << wavelength);
            const double expected = IndependentPlanckRadiance(wavelength, temperature);
            EXPECT_LT(RelativeError(PlanckRadiance(wavelength, temperature), expected), 2.0e-14);
        }
    }
}

TEST(SpectralValidationTests, WienAuthorityMatchesNumericalPlanckMaximum) {
    for (const double temperature : {1200.0, 3000.0, 5778.0, 10000.0}) {
        SCOPED_TRACE(temperature);
        const double analytic = WienPeakWavelength(temperature);
        const double numerical = NumericalPeakWavelength(temperature);
        EXPECT_LT(RelativeError(analytic, numerical), 2.0e-8);
    }
}

TEST(SpectralValidationTests, StefanBoltzmannAuthorityMatchesIntegratedPlanckSpectrum) {
    for (const double temperature : {1200.0, 3000.0, 5778.0, 10000.0}) {
        SCOPED_TRACE(temperature);
        const double analytic = StefanBoltzmannRadiance(temperature);
        const double quadrature = IntegratedExitance(temperature);
        EXPECT_LT(RelativeError(analytic, quadrature), 2.0e-7);
    }
}

TEST(SpectralValidationTests, DopplerFactorMatchesIndependentLorentzFormula) {
    for (const double beta : {-0.95, -0.5, 0.0, 0.5, 0.95}) {
        SCOPED_TRACE(beta);
        const double expected = std::sqrt((1.0 - beta) / (1.0 + beta));
        EXPECT_NEAR(DopplerFactor(beta), expected, 2.0e-15);
    }
    EXPECT_DOUBLE_EQ(DopplerFactor(-1.0), 1.0);
    EXPECT_DOUBLE_EQ(DopplerFactor(1.0), 1.0);
}

TEST(SpectralValidationTests, TotalRedshiftComposesGravitationalAndDopplerFactors) {
    constexpr double kEmitGtt = -0.42;
    constexpr double kObserverGtt = -0.98;
    constexpr double kBeta = 0.37;
    const double gravitational = std::sqrt(std::abs(kObserverGtt / kEmitGtt));
    const double doppler = std::sqrt((1.0 - kBeta) / (1.0 + kBeta));
    const double expected = 1.0 / (gravitational * doppler) - 1.0;
    EXPECT_NEAR(TotalRedshift(kEmitGtt, kObserverGtt, kBeta), expected, 2.0e-15);
}

TEST(SpectralValidationTests, LimbDarkeningAuthorityMatchesEmpiricalLawAndLiveApplication) {
    const Rgb input{0.8f, 0.6f, 0.4f};
    constexpr double kCosTheta = 0.31;
    const Rgb output = ApplyLimbDarkening(input, kCosTheta);
    for (const auto [wavelength, channel_in, channel_out] : std::array{
             std::array<double, 3>{650.0, input.r, output.r},
             std::array<double, 3>{550.0, input.g, output.g},
             std::array<double, 3>{450.0, input.b, output.b},
         }) {
        const double t = (wavelength - 400.0) / 300.0;
        const double expected_u = 0.9 + t * (0.4 - 0.9);
        const double expected_channel =
            channel_in * (1.0 + expected_u * kCosTheta) / (1.0 + expected_u);
        EXPECT_NEAR(LimbDarkeningCoeff(wavelength), expected_u, 2.0e-15);
        EXPECT_NEAR(channel_out, expected_channel, 2.0e-7);
    }
    const Rgb hidden = ApplyLimbDarkening(input, 0.0);
    EXPECT_FLOAT_EQ(hidden.r, 0.0f);
    EXPECT_FLOAT_EQ(hidden.g, 0.0f);
    EXPECT_FLOAT_EQ(hidden.b, 0.0f);
}

TEST(SpectralValidationTests, BlackbodyColourProgressionConsumesIntegratedSpectrum) {
    const Rgb cold = BlackbodyToRgb(2000.0);
    const Rgb solar = BlackbodyToRgb(5500.0);
    const Rgb hot = BlackbodyToRgb(10000.0);
    const auto blue_to_red = [](const Rgb& colour) {
        return colour.b / std::max(colour.r, 1.0e-6f);
    };
    EXPECT_GT(cold.r, cold.b);
    EXPECT_GT(hot.b, hot.r);
    EXPECT_LT(blue_to_red(cold), blue_to_red(solar));
    EXPECT_LT(blue_to_red(solar), blue_to_red(hot));
}

}  // namespace sirius::test
