#pragma once

// Analytic approximation to the CIE 1931 2-degree standard-observer colour
// matching functions over 380--780 nm. This is the Wyman--Sloan--Shirley fit,
// not the official CIE tabulation; wavelength is expressed in nanometres.
// Values outside the represented band fail closed to zero.

#include <cmath>

namespace sirius::core::colour {

struct CieXyzMatching {
    double x_bar;
    double y_bar;
    double z_bar;
};

namespace detail {

inline double CieGaussian(double wavelength_nanometres, double centre, double left_width,
                          double right_width, double amplitude) {
    const double width = wavelength_nanometres < centre ? left_width : right_width;
    const double offset = (wavelength_nanometres - centre) * width;
    return amplitude * std::exp(-0.5 * offset * offset);
}

}  // namespace detail

inline CieXyzMatching Cie1931TwoDegreeFit(double wavelength_nanometres) {
    if (!std::isfinite(wavelength_nanometres) || wavelength_nanometres < 380.0 ||
        wavelength_nanometres > 780.0) {
        return {0.0, 0.0, 0.0};
    }

    const double x_bar = detail::CieGaussian(wavelength_nanometres, 442.0, 0.0624, 0.0374, 0.362) +
                         detail::CieGaussian(wavelength_nanometres, 599.8, 0.0264, 0.0323, 1.056) -
                         detail::CieGaussian(wavelength_nanometres, 501.1, 0.0490, 0.0382, 0.065);
    const double y_bar = detail::CieGaussian(wavelength_nanometres, 568.8, 0.0213, 0.0247, 0.821) +
                         detail::CieGaussian(wavelength_nanometres, 530.9, 0.0613, 0.0322, 0.286);
    const double z_bar = detail::CieGaussian(wavelength_nanometres, 437.0, 0.0845, 0.0278, 1.217) +
                         detail::CieGaussian(wavelength_nanometres, 459.0, 0.0385, 0.0725, 0.681);
    return {x_bar, y_bar, z_bar};
}

}  // namespace sirius::core::colour
