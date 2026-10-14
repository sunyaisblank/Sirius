#pragma once

#include <algorithm>
#include <array>
#include <cmath>

namespace sirius::test::point_source_band_reference {

using Channels = std::array<long double, 3>;

// Independent defining Planck equation, SI defining constants and long-double
// arithmetic. No call to the production Planck, colour or transfer helpers.
inline long double Planck(long double wavelength, long double temperature) {
    constexpr long double h = 6.62607015e-34L, c = 299792458.0L, k = 1.380649e-23L;
    return 2 * h * c * c /
           (std::pow(wavelength, 5) * std::expm1(h * c / (wavelength * k * temperature)));
}

inline Channels Matching(long double nm) {
    const auto gaussian = [nm](long double centre, long double left, long double right) {
        const long double t = (nm - centre) * (nm < centre ? left : right);
        return std::exp(-t * t / 2);
    };
    return {0.362L * gaussian(442, .0624L, .0374L) + 1.056L * gaussian(599.8L, .0264L, .0323L) -
                .065L * gaussian(501.1L, .0490L, .0382L),
            .821L * gaussian(568.8L, .0213L, .0247L) + .286L * gaussian(530.9L, .0613L, .0322L),
            1.217L * gaussian(437, .0845L, .0278L) + .681L * gaussian(459, .0385L, .0725L)};
}

inline Channels Rgb(const Channels& xyz) {
    return {
        std::max(0.0L, 12831.L / 3959 * xyz[0] - 329.L / 214 * xyz[1] - 1974.L / 3959 * xyz[2]),
        std::max(0.0L, -851781.L / 878810 * xyz[0] + 1648619.L / 878810 * xyz[1] +
                           36519.L / 878810 * xyz[2]),
        std::max(0.0L, 705.L / 12673 * xyz[0] - 2585.L / 12673 * xyz[1] + 705.L / 667 * xyz[2])};
}

// shift_wavelength uses g^5 B_lambda(g lambda,T) independently of gT.
inline Channels Integrate(long double temperature, long double g, int bins, bool simpson,
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

inline Channels Reference(long double temperature, long double g) {
    const auto source = Rgb(Integrate(temperature, 1, 32, false, false));
    auto shifted = Rgb(Integrate(temperature, g, 32, false, true));
    const long double normalizer = std::max({source[0], source[1], source[2], .001L});
    for (auto& channel : shifted) channel /= normalizer;
    return shifted;
}

}  // namespace sirius::test::point_source_band_reference
