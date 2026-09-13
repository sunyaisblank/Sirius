#include "sirius/core/spectral/point_source_transfer.h"

#include "sirius/core/point_source_response.h"
#include "sirius/core/starfield.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numbers>

namespace sirius::core {
namespace {

using spectral::PointSourceTransferFailure;
using spectral::TransferPointSourceBand;
using Channels = std::array<long double, 3>;

// Independent defining Planck equation, SI defining constants and long-double
// arithmetic. No call to the production Planck, colour or transfer helpers.
long double Planck(long double wavelength, long double temperature) {
    constexpr long double h = 6.62607015e-34L, c = 299792458.0L, k = 1.380649e-23L;
    return 2 * h * c * c /
           (std::pow(wavelength, 5) * std::expm1(h * c / (wavelength * k * temperature)));
}

Channels Matching(long double nm) {
    const auto gaussian = [nm](long double centre, long double left, long double right) {
        const long double t = (nm - centre) * (nm < centre ? left : right);
        return std::exp(-t * t / 2);
    };
    return {0.362L * gaussian(442, .0624L, .0374L) + 1.056L * gaussian(599.8L, .0264L, .0323L) -
                .065L * gaussian(501.1L, .0490L, .0382L),
            .821L * gaussian(568.8L, .0213L, .0247L) + .286L * gaussian(530.9L, .0613L, .0322L),
            1.217L * gaussian(437, .0845L, .0278L) + .681L * gaussian(459, .0385L, .0725L)};
}

Channels Rgb(const Channels& xyz) {
    return {
        std::max(0.0L, 12831.L / 3959 * xyz[0] - 329.L / 214 * xyz[1] - 1974.L / 3959 * xyz[2]),
        std::max(0.0L, -851781.L / 878810 * xyz[0] + 1648619.L / 878810 * xyz[1] +
                           36519.L / 878810 * xyz[2]),
        std::max(0.0L, 705.L / 12673 * xyz[0] - 2585.L / 12673 * xyz[1] + 705.L / 667 * xyz[2])};
}

// shift_wavelength uses g^5 B_lambda(g lambda,T) independently of gT.
Channels Integrate(long double temperature, long double g, int bins, bool simpson,
                   bool shift_wavelength, bool truncate_source = false) {
    Channels xyz{};
    const long double step = 400e-9L / bins;
    const int count = simpson ? bins + 1 : bins;
    for (int i = 0; i < count; ++i) {
        const long double lambda = 380e-9L + (i + (simpson ? 0.L : .5L)) * step;
        long double value = shift_wavelength ? std::pow(g, 5) * Planck(g * lambda, temperature)
                                             : Planck(lambda, g * temperature);
        if (truncate_source && (g * lambda < 380e-9L || g * lambda > 780e-9L)) value = 0;
        const long double weight = simpson ? (i == 0 || i == bins ? 1 : i % 2 ? 4 : 2) / 3.L : 1;
        const auto matching = Matching(lambda * 1e9L);
        for (int c = 0; c < 3; ++c) xyz[c] += matching[c] * value * step * weight;
    }
    return xyz;
}

Channels Reference(long double temperature, long double g) {
    const auto source = Rgb(Integrate(temperature, 1, 32, false, false));
    auto shifted = Rgb(Integrate(temperature, g, 32, false, true));
    const long double normalizer = std::max({source[0], source[1], source[2], .001L});
    for (auto& channel : shifted) channel /= normalizer;
    return shifted;
}

TEST(PointSourceTransfer, ReferenceFrequencyPreservesCatalogueColourAndFlux) {
    for (float temperature : {100.f, 300.f, 1000.f, 3000.f, 6500.f, 20000.f, 100000.f, 1e6f}) {
        StarEntry star{1, 0, 0, 10, 4.25f, .2f, temperature, 0};
        ASSERT_TRUE(IsRepresentedStarEntry(star));
        float r, g, b;
        star.ComputeColor(r, g, b);
        const std::array<float, 3> colour{r, g, b};
        const double flux = static_cast<double>(star.Intensity()) * 100.0;
        const auto unit = TransferPointSourceBand(temperature, 1, 1, 1);
        const auto weighted = TransferPointSourceBand(temperature, 1, flux, 7.25);
        ASSERT_TRUE(unit);
        ASSERT_TRUE(weighted);
        for (int c = 0; c < 3; ++c) {
            // One continuous double spectrum, original float N(T): reference
            // agreement is roundoff-limited, not a special g==1 branch.
            EXPECT_NEAR((*unit)[c], static_cast<double>(colour[c]), 2e-6);
            EXPECT_DOUBLE_EQ((*weighted)[c], (*unit)[c] * flux * 7.25);
        }
        if (temperature >= 3000 && temperature <= 20000) {
            for (double frequency : {std::nextafter(1., 0.), std::nextafter(1., 2.)}) {
                const auto adjacent = TransferPointSourceBand(temperature, frequency, 1, 1);
                ASSERT_TRUE(adjacent);
                for (int c = 0; c < 3; ++c) EXPECT_NEAR((*adjacent)[c], (*unit)[c], 3e-14);
            }
            constexpr double step = 1e-7;
            const auto plus = TransferPointSourceBand(temperature, 1 + step, 1, 1);
            const auto minus = TransferPointSourceBand(temperature, 1 - step, 1, 1);
            ASSERT_TRUE(plus);
            ASSERT_TRUE(minus);
            for (int c = 0; c < 3; ++c) {
                const double forward = ((*plus)[c] - (*unit)[c]) / step;
                const double backward = ((*unit)[c] - (*minus)[c]) / step;
                EXPECT_NEAR(forward, backward, 1e-5);
            }
        }
    }
}

TEST(PointSourceTransfer, ShiftedSpectrumMatchesIndependentFixedBandPlanckLaw) {
    for (double temperature : {3000., 6500., 20000.}) {
        for (double g : {.5, 1., 1.25, 2.}) {
            const auto actual = TransferPointSourceBand(temperature, g, 2.5, .75);
            ASSERT_TRUE(actual);
            const auto expected = Reference(temperature, g);
            const long double scale = std::max({expected[0], expected[1], expected[2]});
            for (int c = 0; c < 3; ++c)
                EXPECT_NEAR((*actual)[c], static_cast<double>(expected[c] * 2.5L * .75L),
                            static_cast<double>(2e-6L * scale * 2.5L * .75L));
            const auto shifted = Integrate(temperature, g, 8192, true, true);
            const auto direct = Integrate(temperature, g, 8192, true, false);
            const auto coarse = Integrate(temperature, g, 4096, true, false);
            for (int c = 0; c < 3; ++c) {
                EXPECT_NEAR(static_cast<double>(shifted[c] / direct[c]), 1.0, 2e-13);
                EXPECT_NEAR(static_cast<double>(coarse[c] / direct[c]), 1.0, 2e-8);
            }
        }
    }
}

TEST(PointSourceTransfer, RejectsReciprocalExtraGainAndShiftedNormalizerSubstitutions) {
    constexpr double temperature = 6500;
    for (double g : {.5, 2.}) {
        const auto actual = TransferPointSourceBand(temperature, g, 1, 1);
        ASSERT_TRUE(actual);
        const auto right = Reference(temperature, g);
        const auto reciprocal = Reference(temperature, 1 / g);
        const auto renormalized = Reference(temperature * g, 1);
        const auto truncated = Rgb(Integrate(temperature, g, 32, false, true, true));
        const auto source = Rgb(Integrate(temperature, 1, 32, false, false));
        const long double normalizer = std::max({source[0], source[1], source[2], .001L});
        for (int c = 0; c < 3; ++c) {
            const double tolerance = static_cast<double>(right[c] * 2e-6L);
            EXPECT_NEAR((*actual)[c], static_cast<double>(right[c]), tolerance);
            EXPECT_GT(std::abs((*actual)[c] - static_cast<double>(reciprocal[c])), 100 * tolerance);
            EXPECT_GT(std::abs((*actual)[c] - static_cast<double>(right[c] * std::pow(g, 4))),
                      100 * tolerance);
            EXPECT_GT(std::abs((*actual)[c] - static_cast<double>(renormalized[c])),
                      100 * tolerance);
            EXPECT_GT(std::abs((*actual)[c] - static_cast<double>(truncated[c] / normalizer)),
                      100 * tolerance);
        }
    }
}

TEST(PointSourceTransfer, MovingObserverUsesBandGainAndAngularAreaExactlyOnce) {
    constexpr double beta = .6, gamma = 1.25, sigma = .001;
    for (double source_cosine : {1., -1., 0.}) {
        // Minkowski contraction of p=(1,-N) and u=gamma(1,0,0,beta).
        const double g = gamma * (1 + beta * source_cosine);
        const double expected_g = source_cosine == 1 ? 2 : source_cosine == -1 ? .5 : 1.25;
        EXPECT_DOUBLE_EQ(g, expected_g);
        const double observed_cosine = (source_cosine + beta) / (1 + beta * source_cosine);
        const double area = (1 - beta * beta) / std::pow(1 - beta * observed_cosine, 2);
        EXPECT_NEAR(area, g * g, 2e-15);
        for (double parity : {-1., 1.}) {
            const AngularMatrix2 source_map{{{parity * g, 0}, {0, g}}};
            const AngularMatrix2 film_map{{{sigma, 0}, {0, sigma}}};
            const auto response = MakeAffinePointResponse(source_map, film_map, 1, .001);
            ASSERT_TRUE(response);
            const auto band = TransferPointSourceBand(6500, g, 1, 1);
            ASSERT_TRUE(band);
            // Integrate circular image support. Density already supplies 1/|J|.
            constexpr int intervals = 2048;
            double integral = 0;
            // Two-node Gauss quadrature uses strictly interior points, so
            // rounding of a measure-zero compact-support edge cannot bias it.
            const double interval = 4 * sigma / intervals;
            for (int i = 0; i < intervals; ++i) {
                for (double node : {-1 / std::sqrt(3.), 1 / std::sqrt(3.)}) {
                    const double r = interval * (i + .5 + .5 * node);
                    const auto density = response->Density({parity * g * r, 0});
                    ASSERT_TRUE(density);
                    integral += .5 * interval * 2 * std::numbers::pi * r * *density;
                }
            }
            EXPECT_NEAR(integral, 1 / area, 1e-10 / area);
            const auto reference = Reference(6500, g);
            for (int c = 0; c < 3; ++c)
                EXPECT_NEAR((*band)[c] * integral, static_cast<double>(reference[c]) / area,
                            static_cast<double>(reference[c]) * 2e-6 / area);
            // A separately defined bolometric source would instead yield g².
            EXPECT_NEAR(std::pow(g, 4) * integral, g * g, 1e-9);
        }
    }
}

TEST(PointSourceTransfer, ShiftedTemperatureHasNoCatalogueLookupClamp) {
    for (const auto& fixture :
         std::array<std::array<double, 2>, 3>{{{100, .5}, {1e6, 2}, {6500, 20}}}) {
        const double temperature = fixture[0], g = fixture[1];
        const auto actual = TransferPointSourceBand(temperature, g, 1, 1);
        ASSERT_TRUE(actual);
        const auto expected = Reference(temperature, g);
        const double scale = static_cast<double>(std::max({expected[0], expected[1], expected[2]}));
        for (int c = 0; c < 3; ++c)
            EXPECT_NEAR((*actual)[c], static_cast<double>(expected[c]), 2e-6 * scale);
    }
    // This band is below binary32 but not double. A represented bright source
    // recovers visible radiance; do not narrow the shifted spectrum first.
    const auto faint = TransferPointSourceBand(1000, .15, 1e38, 1);
    ASSERT_TRUE(faint);
    const auto faint_reference = Reference(1000, .15);
    EXPECT_GT(std::max({(*faint)[0], (*faint)[1], (*faint)[2]}), 0.0);
    const long double faint_scale =
        std::max({faint_reference[0], faint_reference[1], faint_reference[2]});
    for (int c = 0; c < 3; ++c)
        EXPECT_NEAR((*faint)[c], static_cast<double>(faint_reference[c] * 1e38L),
                    static_cast<double>(faint_scale * 1e38L * 2e-6L));
    const auto hot = TransferPointSourceBand(6500, 1e100, 1, 1);
    ASSERT_TRUE(hot);
    EXPECT_TRUE(std::isfinite((*hot)[0]));
    EXPECT_GT((*hot)[0], 1e90);
    const auto zero = TransferPointSourceBand(100, .001, 1, 1);
    ASSERT_TRUE(zero);
    EXPECT_EQ(*zero, (std::array<double, 3>{0, 0, 0}));
}

TEST(PointSourceTransfer, InvalidOrOverflowingTransferHasNoPublishedValue) {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();
    for (int field = 0; field < 4; ++field) {
        for (double invalid : {nan, inf, -inf, -1.}) {
            std::array<double, 4> inputs{6500, 1, 1, 1};
            inputs[field] = invalid;
            const auto result = TransferPointSourceBand(inputs[0], inputs[1], inputs[2], inputs[3]);
            ASSERT_FALSE(result);
            EXPECT_EQ(result.error(), PointSourceTransferFailure::InvalidInput);
        }
    }
    EXPECT_FALSE(TransferPointSourceBand(0, 1, 0, 0));
    EXPECT_FALSE(TransferPointSourceBand(6500, 0, 0, 0));
    for (const auto& inputs : std::array<std::array<double, 4>, 4>{
             {{1e308, 2, 1, 1},
              {std::numeric_limits<double>::denorm_min(), .5, 1, 1},
              {6500, 1e304, 1, 1},
              {6500, 2, 1e308, 1e308}}}) {
        const auto result = TransferPointSourceBand(inputs[0], inputs[1], inputs[2], inputs[3]);
        ASSERT_FALSE(result);
        EXPECT_EQ(result.error(), PointSourceTransferFailure::Arithmetic);
    }
    const auto scaled = TransferPointSourceBand(6500, 2, 1e308, 1e-300);
    const auto unit = TransferPointSourceBand(6500, 2, 1, 1);
    ASSERT_TRUE(scaled);
    ASSERT_TRUE(unit);
    for (int c = 0; c < 3; ++c) EXPECT_NEAR((*scaled)[c] / (*unit)[c], 1e8, 1e-6);
    const auto dark = TransferPointSourceBand(6500, 2, 0, 12);
    ASSERT_TRUE(dark);
    EXPECT_EQ(*dark, (std::array<double, 3>{0, 0, 0}));
    EXPECT_FALSE(TransferPointSourceBand(nan, 1, 0, 0));
    EXPECT_FALSE(TransferPointSourceBand(6500, nan, 0, 0));
}

}  // namespace
}  // namespace sirius::core
