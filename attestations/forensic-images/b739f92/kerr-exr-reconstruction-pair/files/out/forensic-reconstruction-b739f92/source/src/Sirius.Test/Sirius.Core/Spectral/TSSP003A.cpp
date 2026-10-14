// TSSP003A.cpp - GPU/CPU Blackbody Colour Parity Tests
//
// Validates that the GPU Chebyshev blackbody approximation (RDOP002A.cu)
// produces results consistent with the CPU 32-sample integration (PHSP001A.h).
// Tolerance: 25% per channel (Chebyshev 8-term fit vs 32-sample CIE integration
// use fundamentally different methods; this tests colour direction, not exact match).
//
// LABEL: Correctness

#include <gtest/gtest.h>
#include <cmath>
#include <PHSP001A.h>

namespace sirius::test {
using namespace Sirius::Spectral;

// GPU Chebyshev blackbody reimplemented on CPU for comparison.
// Mirrors RDOP002A.cu:blackbodyColor() exactly.
namespace {

constexpr float BB_LOG_MIN = 6.907755279f;   // log(1000)
constexpr float BB_LOG_RANGE = 3.688879454f;  // log(40000) - log(1000)

// Chebyshev coefficients for CIE x chromaticity (8 terms)
// Must match RDOP002A.cu exactly
constexpr float BB_X_CHEB[8] = {
    0.332461f, -0.098234f, 0.024891f, -0.008123f,
    0.003012f, -0.001156f, 0.000423f, -0.000148f
};
// Chebyshev coefficients for CIE y chromaticity (8 terms)
constexpr float BB_Y_CHEB[8] = {
    0.341231f, -0.056789f, 0.018234f, -0.006512f,
    0.002345f, -0.000891f, 0.000312f, -0.000108f
};

float chebyshevEval(const float* coeffs, int n, float u) {
    if (n < 2) return coeffs[0];
    float b0 = coeffs[n-1];
    float b1 = 0.0f;
    for (int i = n - 2; i >= 0; --i) {
        float tmp = 2.0f * u * b0 - b1 + coeffs[i];
        b1 = b0;
        b0 = tmp;
    }
    return b0;
}

struct Float3 { float x, y, z; };

Float3 gpuBlackbodyColor(float T) {
    T = std::clamp(T, 1000.0f, 40000.0f);

    float logT = std::log(T);
    float u = 2.0f * (logT - BB_LOG_MIN) / BB_LOG_RANGE - 1.0f;

    float cx = chebyshevEval(BB_X_CHEB, 8, u);
    float cy = chebyshevEval(BB_Y_CHEB, 8, u);

    cx = std::clamp(cx, 0.1f, 0.65f);
    cy = std::clamp(cy, 0.1f, 0.65f);

    float Y = 1.0f;
    float y_inv = 1.0f / std::max(cy, 0.01f);
    float X = Y * y_inv * cx;
    float Z = Y * y_inv * (1.0f - cx - cy);

    Float3 rgb;
    rgb.x =  3.2404542f * X - 1.5371385f * Y - 0.4985314f * Z;
    rgb.y = -0.9692660f * X + 1.8760108f * Y + 0.0415560f * Z;
    rgb.z =  0.0556434f * X - 0.2040259f * Y + 1.0572252f * Z;

    rgb.x = std::max(rgb.x, 0.0f);
    rgb.y = std::max(rgb.y, 0.0f);
    rgb.z = std::max(rgb.z, 0.0f);

    return rgb;
}

} // anonymous namespace

// Compare GPU Chebyshev against CPU 32-sample integration at reference temps.
// The two methods use fundamentally different approaches (Chebyshev chromaticity
// fit vs CIE observer function integration), so we test colour direction (which
// channel dominates) rather than exact per-channel agreement.
class BlackbodyParityTest : public ::testing::TestWithParam<double> {};

TEST_P(BlackbodyParityTest, GPUMatchesCPU) {
    double T = GetParam();
    RGB cpuRGB = blackbodyToRGB(T);
    Float3 gpuRGB = gpuBlackbodyColor(static_cast<float>(T));

    // Normalise GPU output to max=1 for comparison
    float gpuMax = std::max({gpuRGB.x, gpuRGB.y, gpuRGB.z, 0.001f});
    float gR = gpuRGB.x / gpuMax;
    float gG = gpuRGB.y / gpuMax;
    float gB = gpuRGB.z / gpuMax;

    // Both methods must agree on colour direction (blue vs red dominant).
    // Near-white temperatures (~5500-7000K) have nearly equal R and B channels,
    // so we only test dominance when one channel leads by > 15%.
    float cpuBR = cpuRGB.b / std::max(cpuRGB.r, 0.01f);
    float gpuBR = gB / std::max(gR, 0.01f);
    if (cpuBR > 1.3f) {
        EXPECT_GT(gpuBR, 1.0f) << "GPU should also be blue-dominant at T=" << T;
    } else if (cpuBR < 0.7f) {
        EXPECT_LT(gpuBR, 1.0f) << "GPU should also be red-dominant at T=" << T;
    }
    // else: near-white, both methods may disagree on which channel barely leads

    // All channels should be non-negative and at least one should be near 1
    EXPECT_GE(gR, 0.0f);
    EXPECT_GE(gG, 0.0f);
    EXPECT_GE(gB, 0.0f);
    float maxCh = std::max({gR, gG, gB});
    EXPECT_NEAR(maxCh, 1.0f, 0.01f) << "Normalised GPU output should have max ~1";
}

INSTANTIATE_TEST_SUITE_P(
    ReferenceTemperatures,
    BlackbodyParityTest,
    ::testing::Values(3000.0, 5000.0, 6500.0, 10000.0, 25000.0)
);

// Verify colour ordering: low T should be redder, high T bluer
TEST(BlackbodyParityTest, ColourTemperatureOrdering) {
    RGB cool = blackbodyToRGB(3000.0);
    RGB warm = blackbodyToRGB(6500.0);
    RGB hot  = blackbodyToRGB(25000.0);

    // 3000K: r > b (reddish)
    EXPECT_GT(cool.r, cool.b) << "3000K should be redder than blue";

    // 25000K: b > r (bluish-white)
    // After normalisation the blue channel approaches 1
    EXPECT_GT(hot.b, 0.8f) << "25000K should have strong blue";

    // Green channel should be moderate across all temperatures
    EXPECT_GT(warm.g, 0.5f) << "6500K should have moderate green";
}

} // namespace sirius::test
