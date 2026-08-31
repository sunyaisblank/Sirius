#include "sirius/core/cie1931_observer.h"
#include "sirius/core/spectral/blackbody.h"
#include "sirius/core/spectral/spectral_radiance.h"

#include <gtest/gtest.h>

#include <array>
#include <limits>

namespace {

struct OfficialCieSample {
    double wavelength_nanometres;
    double x_bar;
    double y_bar;
    double z_bar;
};

// Independent oracle samples from CIE 018:2019, Table 6, DOI
// 10.25039/CIE.DS.xvudnb9b. The source CSV's CIE-published MD5 is
// 17cca777db64b17170f06f67ce9d3ab7 and its SHA-256 is
// fa663e3535a7e0763a745993a1f0a192eb0275ac46ad2d1befd7626841e713c1.
constexpr std::array<OfficialCieSample, 22> kOfficialCieSamples = {{
    {380.0, 0.001368000000, 0.0000390000000, 0.006450001000},
    {400.0, 0.014310000000, 0.0003960000000, 0.067850010000},
    {420.0, 0.134380000000, 0.0040000000000, 0.645600000000},
    {440.0, 0.348280000000, 0.0230000000000, 1.747060000000},
    {460.0, 0.290800000000, 0.0600000000000, 1.669200000000},
    {480.0, 0.095640000000, 0.1390200000000, 0.812950100000},
    {500.0, 0.004900000000, 0.3230000000000, 0.272000000000},
    {520.0, 0.063270000000, 0.7100000000000, 0.078249990000},
    {540.0, 0.290400000000, 0.9540000000000, 0.020300000000},
    {560.0, 0.594500000000, 0.9950000000000, 0.003900000000},
    {580.0, 0.916300000000, 0.8700000000000, 0.001650001000},
    {600.0, 1.062200000000, 0.6310000000000, 0.000800000000},
    {620.0, 0.854449900000, 0.3810000000000, 0.000190000000},
    {630.0, 0.642400000000, 0.2650000000000, 0.000049999990},
    {640.0, 0.447900000000, 0.1750000000000, 0.000020000000},
    {660.0, 0.164900000000, 0.0610000000000, 0.000000000000},
    {680.0, 0.046770000000, 0.0170000000000, 0.000000000000},
    {700.0, 0.011359160000, 0.0041020000000, 0.000000000000},
    {720.0, 0.002899327000, 0.0010470000000, 0.000000000000},
    {740.0, 0.000690078600, 0.0002492000000, 0.000000000000},
    {760.0, 0.000166150500, 0.0000600000000, 0.000000000000},
    {780.0, 0.000041509940, 0.0000149900000, 0.000000000000},
}};

TEST(Cie1931ObserverAuthority, FitIsBoundedAgainstOfficialCieDatasetAcrossRepresentedBand) {
    for (const OfficialCieSample& reference : kOfficialCieSamples) {
        const sirius::core::colour::CieXyzMatching actual =
            sirius::core::colour::Cie1931TwoDegreeFit(reference.wavelength_nanometres);
        EXPECT_NEAR(actual.x_bar, reference.x_bar, 0.014)
            << reference.wavelength_nanometres << " nm x-bar";
        EXPECT_NEAR(actual.y_bar, reference.y_bar, 0.014)
            << reference.wavelength_nanometres << " nm y-bar";
        EXPECT_NEAR(actual.z_bar, reference.z_bar, 0.014)
            << reference.wavelength_nanometres << " nm z-bar";
    }
}

TEST(Cie1931ObserverAuthority, HostSpectralFacadesEvaluateTheirActualWavelengths) {
    for (int bin = 0; bin < sirius::core::kNumWavelengthBins; ++bin) {
        const double wavelength = sirius::core::SpectralRadiance::Wavelength(bin);
        const auto matching = sirius::core::colour::Cie1931TwoDegreeFit(wavelength);

        sirius::core::SpectralRadiance impulse = sirius::core::SpectralRadiance::Zero();
        impulse.L[bin] = 1.0;
        const auto integrated = impulse.ToXyz();
        EXPECT_DOUBLE_EQ(integrated.X, matching.x_bar * sirius::core::kLambdaStep);
        EXPECT_DOUBLE_EQ(integrated.Y, matching.y_bar * sirius::core::kLambdaStep);
        EXPECT_DOUBLE_EQ(integrated.Z, matching.z_bar * sirius::core::kLambdaStep);

        const auto facade = sirius::core::spectral::WavelengthToXyz(wavelength);
        EXPECT_FLOAT_EQ(facade.X, static_cast<float>(matching.x_bar));
        EXPECT_FLOAT_EQ(facade.Y, static_cast<float>(matching.y_bar));
        EXPECT_FLOAT_EQ(facade.Z, static_cast<float>(matching.z_bar));
    }
}

TEST(Cie1931ObserverAuthority, UnrepresentedWavelengthsFailClosed) {
    for (const double wavelength : {
             -std::numeric_limits<double>::infinity(),
             379.999,
             780.001,
             std::numeric_limits<double>::infinity(),
             std::numeric_limits<double>::quiet_NaN(),
         }) {
        const auto matching = sirius::core::colour::Cie1931TwoDegreeFit(wavelength);
        EXPECT_DOUBLE_EQ(matching.x_bar, 0.0);
        EXPECT_DOUBLE_EQ(matching.y_bar, 0.0);
        EXPECT_DOUBLE_EQ(matching.z_bar, 0.0);
    }
    EXPECT_GT(sirius::core::colour::Cie1931TwoDegreeFit(380.0).z_bar, 0.0);
    EXPECT_GT(sirius::core::colour::Cie1931TwoDegreeFit(780.0).x_bar, 0.0);
}

}  // namespace
