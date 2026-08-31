#include "sirius/core/constants.h"
#include "sirius/core/spectral/blackbody.h"
#include "sirius/core/spectral/spectral_radiance.h"
#include "sirius/core/srgb_transfer.h"
#include "sirius/core/xyz_srgb.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

namespace sirius::core::test {
namespace {

struct Vector3 {
    double x;
    double y;
    double z;
};

struct Matrix3 {
    double value[3][3];
};

Vector3 Multiply(const Matrix3& matrix, const Vector3& vector) {
    return {
        matrix.value[0][0] * vector.x + matrix.value[0][1] * vector.y +
            matrix.value[0][2] * vector.z,
        matrix.value[1][0] * vector.x + matrix.value[1][1] * vector.y +
            matrix.value[1][2] * vector.z,
        matrix.value[2][0] * vector.x + matrix.value[2][1] * vector.y +
            matrix.value[2][2] * vector.z,
    };
}

Matrix3 Inverse(const Matrix3& matrix) {
    const auto& m = matrix.value;
    const double determinant = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) -
                               m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
                               m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);

    return {{{(m[1][1] * m[2][2] - m[1][2] * m[2][1]) / determinant,
              (m[0][2] * m[2][1] - m[0][1] * m[2][2]) / determinant,
              (m[0][1] * m[1][2] - m[0][2] * m[1][1]) / determinant},
             {(m[1][2] * m[2][0] - m[1][0] * m[2][2]) / determinant,
              (m[0][0] * m[2][2] - m[0][2] * m[2][0]) / determinant,
              (m[0][2] * m[1][0] - m[0][0] * m[1][2]) / determinant},
             {(m[1][0] * m[2][1] - m[1][1] * m[2][0]) / determinant,
              (m[0][1] * m[2][0] - m[0][0] * m[2][1]) / determinant,
              (m[0][0] * m[1][1] - m[0][1] * m[1][0]) / determinant}}};
}

Vector3 UnitLuminancePrimary(double x, double y) { return {x / y, 1.0, (1.0 - x - y) / y}; }

Matrix3 IndependentlyDerivedXyzD65ToLinearSrgb() {
    // IEC sRGB primary chromaticities and the CIE D65 white chromaticity.
    const Vector3 red = UnitLuminancePrimary(0.640, 0.330);
    const Vector3 green = UnitLuminancePrimary(0.300, 0.600);
    const Vector3 blue = UnitLuminancePrimary(0.150, 0.060);
    const Vector3 d65 = UnitLuminancePrimary(0.3127, 0.3290);

    const Matrix3 unscaled_primaries = {
        {{red.x, green.x, blue.x}, {red.y, green.y, blue.y}, {red.z, green.z, blue.z}}};
    const Vector3 scales = Multiply(Inverse(unscaled_primaries), d65);
    const Matrix3 linear_srgb_to_xyz = {
        {{red.x * scales.x, green.x * scales.y, blue.x * scales.z},
         {red.y * scales.x, green.y * scales.y, blue.y * scales.z},
         {red.z * scales.x, green.z * scales.y, blue.z * scales.z}}};
    return Inverse(linear_srgb_to_xyz);
}

TEST(XyzSrgbAuthority, MatchesIndependentPrimaryWhiteDerivation) {
    const Matrix3 oracle = IndependentlyDerivedXyzD65ToLinearSrgb();
    constexpr std::array samples{
        Vector3{1.0, 0.0, 0.0},  Vector3{0.0, 1.0, 0.0},
        Vector3{0.0, 0.0, 1.0},  Vector3{0.9504559270516716, 1.0, 1.0890577507598784},
        Vector3{0.25, 0.4, 0.1}, Vector3{-0.1, 0.2, 1.5},
    };

    for (const Vector3& xyz : samples) {
        const Vector3 expected = Multiply(oracle, xyz);
        const colour::LinearSrgb actual = colour::XyzD65ToLinearSrgb(xyz.x, xyz.y, xyz.z);
        const double scale =
            std::max({1.0, std::abs(expected.x), std::abs(expected.y), std::abs(expected.z)});
        const double tolerance = 64.0 * constants::kEpsilonD * scale;
        EXPECT_NEAR(actual.r, expected.x, tolerance);
        EXPECT_NEAR(actual.g, expected.y, tolerance);
        EXPECT_NEAR(actual.b, expected.z, tolerance);
    }
}

TEST(XyzSrgbAuthority, PreservesExtendedGamutAndPropagatesNonfiniteInputs) {
    const colour::LinearSrgb green_axis = colour::XyzD65ToLinearSrgb(0.0, 1.0, 0.0);
    EXPECT_LT(green_axis.r, 0.0);
    EXPECT_GT(green_axis.g, 1.0);
    EXPECT_LT(green_axis.b, 0.0);

    const double nan = std::numeric_limits<double>::quiet_NaN();
    const colour::LinearSrgb malformed = colour::XyzD65ToLinearSrgb(nan, 1.0, 1.0);
    EXPECT_TRUE(std::isnan(malformed.r));
    EXPECT_TRUE(std::isnan(malformed.g));
    EXPECT_TRUE(std::isnan(malformed.b));

    const colour::LinearSrgb unbounded =
        colour::XyzD65ToLinearSrgb(std::numeric_limits<double>::infinity(), 1.0, 1.0);
    EXPECT_FALSE(std::isfinite(unbounded.r));
    EXPECT_FALSE(std::isfinite(unbounded.g));
    EXPECT_FALSE(std::isfinite(unbounded.b));
}

TEST(XyzSrgbAuthority, HostSpectralFacadesDelegateToTheExactTransform) {
    const spectral::Xyz xyz{0.25f, 0.4f, 0.1f};
    const spectral::Rgb facade = spectral::XyzToLinearRgb(xyz);
    const colour::LinearSrgb direct = colour::XyzD65ToLinearSrgb(xyz.X, xyz.Y, xyz.Z);
    EXPECT_FLOAT_EQ(facade.r, direct.r);
    EXPECT_FLOAT_EQ(facade.g, direct.g);
    EXPECT_FLOAT_EQ(facade.b, direct.b);

    SpectralRadiance spectrum = SpectralRadiance::Zero();
    spectrum.L[4] = 0.002;
    spectrum.L[12] = 0.004;
    spectrum.L[23] = 0.001;
    const SpectralRadiance::Xyz spectrum_xyz = spectrum.ToXyz();
    const colour::LinearSrgb spectrum_linear =
        colour::XyzD65ToLinearSrgb(spectrum_xyz.X, spectrum_xyz.Y, spectrum_xyz.Z);
    const SpectralRadiance::Rgb encoded = spectrum.ToSrgb();
    EXPECT_DOUBLE_EQ(encoded.r, colour::EncodeClippedSrgbChannel(spectrum_linear.r));
    EXPECT_DOUBLE_EQ(encoded.g, colour::EncodeClippedSrgbChannel(spectrum_linear.g));
    EXPECT_DOUBLE_EQ(encoded.b, colour::EncodeClippedSrgbChannel(spectrum_linear.b));
}

}  // namespace
}  // namespace sirius::core::test
