// Full-image smoke gate for the trace kernel (workstream deliverable 3): a
// 64x64 Kerr render with no disk, background shaded by escape direction, must
// dispatch on Lavapipe and produce a finite, non-constant field whose
// horizon-shadow fraction sits in a broad sane band. This is a behavioural
// smoke, not a pixel-parity gate; the physics it exercises (Kerr-Schild metric,
// Christoffels, Cartesian RK4) is pinned by kernel_parity_test.cpp.

#include "sirius/backend/device.h"

#include <gtest/gtest.h>

#include "../support/trace_continuation_probe.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <span>
#include <string>
#include <vector>

namespace {

using sirius::backend::BufferHandle;
using sirius::backend::BufferUsage;
using sirius::backend::ComputeDevice;
using sirius::backend::CreateVulkanDevice;
using sirius::backend::EnumerateVulkanDevices;
using sirius::backend::KernelHandle;
using sirius::backend::ResolveVulkanDeviceIndex;

std::vector<std::uint32_t> LoadSpirv(const std::string& path) {
    std::ifstream file(path, std::ios::binary | std::ios::ate);
    if (!file) {
        return {};
    }
    const auto size = static_cast<std::size_t>(file.tellg());
    std::vector<std::uint32_t> words(size / sizeof(std::uint32_t));
    file.seekg(0);
    file.read(reinterpret_cast<char*>(words.data()), static_cast<std::streamsize>(size));
    return words;
}

TEST(KernelTrace, KerrRenderIsFiniteNonConstantWithBoundedShadow) {
#ifndef SIRIUS_KERNEL_DIR
    GTEST_SKIP() << "kernels not compiled (slangc absent at configure time)";
#else
    const auto devices = EnumerateVulkanDevices();
    ASSERT_TRUE(devices.has_value()) << devices.error().Description();
    if (devices->empty()) {
        GTEST_SKIP() << "no Vulkan device present";
    }
    const auto selected = ResolveVulkanDeviceIndex(*devices);
    ASSERT_TRUE(selected.has_value()) << selected.error().Description();
    auto device = CreateVulkanDevice(*selected);
    ASSERT_TRUE(device.has_value()) << device.error().Description();

    const auto spirv = LoadSpirv(std::string(SIRIUS_KERNEL_DIR) + "/trace.spv");
    ASSERT_FALSE(spirv.empty()) << "trace.spv missing or empty";
    const auto kernel = (*device)->LoadKernel(spirv);
    ASSERT_TRUE(kernel.has_value()) << kernel.error().Description();

    constexpr std::uint32_t kWidth = 64;
    constexpr std::uint32_t kHeight = 64;

    // Kerr M=1, a=0.9 viewed down the spin axis from 40M. fov and distance are
    // chosen so the capture cross-section covers a modest image fraction. The
    // full image is a single tile here (tileOrigin 0, tile == image). Background
    // is the analytic gradient (starfield disabled), so no texture is needed.
    std::vector<float> params(sirius::render::kTraceParameterCount, 0.0f);
    params[44] = 0.5f;
    params[45] = 0.5f;
    params[66] = 0.5f;
    params[67] = 0.5f;
    params[51] = 1.0f;
    params[0] = kWidth;   // imageWidth
    params[1] = kHeight;  // imageHeight
    params[2] = 0.0f;     // metricId (Kerr-Schild family)
    params[3] = 1.0f;     // M
    params[4] = 0.9f;     // a
    params[5] = 0.0f;     // Q
    params[6] = 0.0f;     // Lambda
    params[7] = 0.0f;
    params[8] = 0.0f;
    params[9] = -40.0f;  // camera position
    params[10] = 0.0f;
    params[11] = 0.0f;
    params[12] = 1.0f;  // forward
    params[13] = 1.0f;
    params[14] = 0.0f;
    params[15] = 0.0f;  // right
    params[16] = 0.0f;
    params[17] = 1.0f;
    params[18] = 0.0f;            // up
    params[19] = 0.6f;            // fov (radians)
    params[20] = 1.0f;            // aspect
    params[21] = 2000.0f;         // max_steps
    params[22] = 0.06f;           // stepScale
    params[23] = 0.02f;           // min_step
    params[24] = 2.0f;            // max_step
    params[25] = 100.0f;          // escape_radius
    params[26] = 1.12f;           // captureFactor
    params[27] = 0.0f;            // disk disabled
    params[31] = 0.0f;            // tileOriginX
    params[32] = 0.0f;            // tileOriginY
    params[33] = float(kWidth);   // tileWidth
    params[34] = float(kHeight);  // tileHeight
    params[35] = 0.0f;            // starfield disabled (gradient background)
    params[36] = 1.0f;            // starfieldWidth (dummy)
    params[37] = 1.0f;            // starfieldHeight (dummy)

    const auto output =
        sirius::test::TraceContinuationImage<float>(**device, *kernel, params, kWidth, kHeight);
    ASSERT_TRUE(output.has_value()) << output.error().Description();
    const auto& radiance = *output;

    // Finiteness across every channel.
    for (float value : radiance) {
        ASSERT_TRUE(std::isfinite(value)) << "non-finite radiance sample";
    }

    // Non-constancy and shadow fraction from per-pixel luminance.
    float min_lum = std::numeric_limits<float>::max();
    float max_lum = 0.0f;
    std::size_t shadow_pixels = 0;
    for (std::uint32_t p = 0; p < kWidth * kHeight; ++p) {
        const float r = radiance[p * 4 + 0];
        const float g = radiance[p * 4 + 1];
        const float b = radiance[p * 4 + 2];
        const float lum = 0.2126f * r + 0.7152f * g + 0.0722f * b;
        min_lum = std::min(min_lum, lum);
        max_lum = std::max(max_lum, lum);
        if (r < 1e-4f && g < 1e-4f && b < 1e-4f) {
            ++shadow_pixels;
        }
    }
    EXPECT_GT(max_lum - min_lum, 1e-3f) << "radiance field is constant";

    const double fraction =
        static_cast<double>(shadow_pixels) / static_cast<double>(kWidth * kHeight);
    std::cout << "[ trace    ] shadow fraction=" << fraction << " luminance range=[" << min_lum
              << ", " << max_lum << "]\n";
    EXPECT_GE(fraction, 0.005) << "shadow fraction too small: " << fraction;
    EXPECT_LE(fraction, 0.30) << "shadow fraction too large: " << fraction;
#endif
}

#ifdef SIRIUS_KERNEL_DIR
// Dispatch one 64x64 Kerr scene (the same scene as the smoke gate above)
// through the given SPIR-V module and return the RGBA radiance field. Keep the
// direct probe under the product's 64-active-pixel work bound. This helper
// retains conservative row strips; the product packs blocks into 8x8 groups.
template <typename Real>
std::vector<float> RunKerrScene(ComputeDevice& device, const std::vector<std::uint32_t>& spirv) {
    const auto kernel = device.LoadKernel(spirv);
    if (!kernel) {
        ADD_FAILURE() << kernel.error().Description();
        return {};
    }

    constexpr std::uint32_t kWidth = 64;
    constexpr std::uint32_t kHeight = 64;

    std::vector<float> params(sirius::render::kTraceParameterCount, 0.0f);
    params[44] = 0.5f;
    params[45] = 0.5f;
    params[66] = 0.5f;
    params[67] = 0.5f;
    params[51] = 1.0f;
    params[0] = kWidth;
    params[1] = kHeight;
    params[2] = 0.0f;
    params[3] = 1.0f;
    params[4] = 0.9f;
    params[7] = 0.0f;
    params[8] = 0.0f;
    params[9] = -40.0f;
    params[10] = 0.0f;
    params[11] = 0.0f;
    params[12] = 1.0f;
    params[13] = 1.0f;
    params[14] = 0.0f;
    params[15] = 0.0f;
    params[16] = 0.0f;
    params[17] = 1.0f;
    params[18] = 0.0f;
    params[19] = 0.6f;
    params[20] = 1.0f;
    params[21] = 2000.0f;
    params[22] = 0.06f;
    params[23] = 0.02f;
    params[24] = 2.0f;
    params[25] = 100.0f;
    params[26] = 1.12f;
    params[33] = float(kWidth);
    params[34] = float(kHeight);
    params[36] = 1.0f;
    params[37] = 1.0f;

    const auto output =
        sirius::test::TraceContinuationImage<Real>(device, *kernel, params, kWidth, kHeight);
    if (!output) {
        ADD_FAILURE() << output.error().Description();
        return {};
    }
    return *output;
}
#endif

// The fp64 rung (trace_fp64.spv, same source with the Cartesian trajectory
// core widened to double) renders the same Kerr scene as the fp32 kernel and
// agrees with it closely away from discontinuous capture flips. The geometry
// is identical, while each precision policy owns its adaptive step schedule;
// invariant accuracy and physical termination are qualified by KernelParity.
// Verifies the rung is a real double-precision trajectory, not a relabelled
// fp32 module, by requiring the artefact to load, dispatch, and stay finite
// on a shaderFloat64 device (Lavapipe reports it).
TEST(KernelTrace, Fp64RungAgreesWithFp32OnKerrScene) {
#ifndef SIRIUS_KERNEL_DIR
    GTEST_SKIP() << "kernels not compiled (slangc absent at configure time)";
#else
    const auto devices = EnumerateVulkanDevices();
    ASSERT_TRUE(devices.has_value()) << devices.error().Description();
    if (devices->empty()) {
        GTEST_SKIP() << "no Vulkan device present";
    }
    const auto selected = ResolveVulkanDeviceIndex(*devices);
    ASSERT_TRUE(selected.has_value()) << selected.error().Description();
    auto device = CreateVulkanDevice(*selected);
    ASSERT_TRUE(device.has_value()) << device.error().Description();
    const auto spirv32 = LoadSpirv(std::string(SIRIUS_KERNEL_DIR) + "/trace.spv");
    const auto spirv64 = LoadSpirv(std::string(SIRIUS_KERNEL_DIR) + "/trace_fp64.spv");
    ASSERT_FALSE(spirv32.empty()) << "trace.spv missing";
    ASSERT_FALSE(spirv64.empty()) << "trace_fp64.spv missing";

    const auto r32 = RunKerrScene<float>(**device, spirv32);
    ASSERT_FALSE(r32.empty());
    if (!(*device)->Info().supports_fp64) {
        const auto refused = (*device)->LoadKernel(spirv64);
        ASSERT_FALSE(refused.has_value());
        EXPECT_EQ(refused.error().domain(), sirius::base::ErrorDomain::kKernel);
        EXPECT_NE(refused.error().detail().find("shaderFloat64"), std::string::npos);
        for (const float value : r32) EXPECT_TRUE(std::isfinite(value));
        const auto [low, high] = std::minmax_element(r32.begin(), r32.end());
        EXPECT_GT(*high - *low, 1e-3f);
        RecordProperty("fp64_evidence", "unsupported_kernel_declined");
        return;
    }
    RecordProperty("fp64_evidence", "native_comparison_executed");
    const auto r64 = RunKerrScene<double>(**device, spirv64);
    ASSERT_EQ(r32.size(), r64.size());
    ASSERT_FALSE(r64.empty());

    double sum_abs_diff = 0.0;
    double max_abs_diff = 0.0;
    float min64 = std::numeric_limits<float>::max();
    float max64 = 0.0f;
    for (std::size_t i = 0; i < r64.size(); ++i) {
        ASSERT_TRUE(std::isfinite(r64[i])) << "non-finite fp64-rung radiance";
        const double d = std::abs(double(r64[i]) - double(r32[i]));
        sum_abs_diff += d;
        max_abs_diff = std::max(max_abs_diff, d);
        min64 = std::min(min64, r64[i]);
        max64 = std::max(max64, r64[i]);
    }
    const double mean_abs_diff = sum_abs_diff / double(r64.size());

    std::size_t black32 = 0;
    std::size_t black64 = 0;
    for (std::size_t pixel = 0; pixel < r64.size() / 4; ++pixel) {
        const std::size_t base = pixel * 4;
        if (r32[base] < 1e-4f && r32[base + 1] < 1e-4f && r32[base + 2] < 1e-4f) {
            ++black32;
        }
        if (r64[base] < 1e-4f && r64[base + 1] < 1e-4f && r64[base + 2] < 1e-4f) {
            ++black64;
        }
    }

    std::cout << "[ trace64  ] mean|d|=" << mean_abs_diff << " max|d|=" << max_abs_diff
              << " range64=[" << min64 << ", " << max64 << "] black32=" << black32
              << " black64=" << black64 << "\n";

    EXPECT_GT(max64 - min64, 1e-3f) << "fp64 radiance field is constant";
    const std::size_t shadow_count_tolerance = (r64.size() / 4) / 50;  // two percent
    const std::size_t shadow_count_difference =
        black64 > black32 ? black64 - black32 : black32 - black64;
    EXPECT_LE(shadow_count_difference, shadow_count_tolerance)
        << "precision rungs disagree on the physical-termination population";
    // Same scene, independently controlled precision schedules: the fields
    // must agree closely in the mean. Individual pixels may flip across the
    // capture edge, so the max is bounded loosely by the background dynamic
    // range rather than tightly.
    EXPECT_LT(mean_abs_diff, 1e-2) << "fp64 rung diverges from fp32 in the mean";
    EXPECT_LT(max_abs_diff, 2.0) << "fp64 rung wildly diverges at a pixel";
#endif
}

// The compensated rung (trace_fp32comp.spv, Kahan state accumulation) renders
// the same scene, stays finite and non-constant, agrees with plain fp32
// closely, and — measured against the fp64 reference — tracks it at least as
// well as plain fp32 does (within slack for noise; on hardware with sloppier
// fp32 the compensation's gain is larger, which this bound also admits).
TEST(KernelTrace, CompensatedRungTracksFp64AtLeastAsWellAsFp32) {
#ifndef SIRIUS_KERNEL_DIR
    GTEST_SKIP() << "kernels not compiled (slangc absent at configure time)";
#else
    const auto devices = EnumerateVulkanDevices();
    ASSERT_TRUE(devices.has_value()) << devices.error().Description();
    if (devices->empty()) {
        GTEST_SKIP() << "no Vulkan device present";
    }
    const auto selected = ResolveVulkanDeviceIndex(*devices);
    ASSERT_TRUE(selected.has_value()) << selected.error().Description();
    auto device = CreateVulkanDevice(*selected);
    ASSERT_TRUE(device.has_value()) << device.error().Description();
    const auto spirv32 = LoadSpirv(std::string(SIRIUS_KERNEL_DIR) + "/trace.spv");
    const auto spirvC = LoadSpirv(std::string(SIRIUS_KERNEL_DIR) + "/trace_fp32comp.spv");
    const auto spirv64 = LoadSpirv(std::string(SIRIUS_KERNEL_DIR) + "/trace_fp64.spv");
    ASSERT_FALSE(spirv32.empty());
    ASSERT_FALSE(spirvC.empty()) << "trace_fp32comp.spv missing";
    ASSERT_FALSE(spirv64.empty());

    const auto r32 = RunKerrScene<float>(**device, spirv32);
    const auto rC = RunKerrScene<float>(**device, spirvC);
    ASSERT_EQ(rC.size(), r32.size());
    ASSERT_FALSE(rC.empty());

    double sum_c32 = 0.0;
    float minc = std::numeric_limits<float>::max();
    float maxc = 0.0f;
    for (std::size_t i = 0; i < rC.size(); ++i) {
        ASSERT_TRUE(std::isfinite(rC[i])) << "non-finite compensated-rung radiance";
        ASSERT_TRUE(std::isfinite(r32[i])) << "non-finite fp32-rung radiance";
        sum_c32 += std::abs(double(rC[i]) - double(r32[i]));
        minc = std::min(minc, rC[i]);
        maxc = std::max(maxc, rC[i]);
    }
    const double n = double(rC.size());
    const double mean_c32 = sum_c32 / n;
    EXPECT_GT(maxc - minc, 1e-3f) << "compensated radiance field is constant";
    EXPECT_LT(mean_c32, 1e-2) << "compensated rung diverges from fp32 in the mean";
    if (!(*device)->Info().supports_fp64) {
        const auto refused = (*device)->LoadKernel(spirv64);
        ASSERT_FALSE(refused.has_value());
        EXPECT_EQ(refused.error().domain(), sirius::base::ErrorDomain::kKernel);
        EXPECT_NE(refused.error().detail().find("shaderFloat64"), std::string::npos);
        RecordProperty("fp64_evidence", "unsupported_kernel_declined");
        return;
    }
    RecordProperty("fp64_evidence", "native_comparison_executed");
    const auto r64 = RunKerrScene<double>(**device, spirv64);
    ASSERT_EQ(rC.size(), r64.size());
    double sum_c64 = 0.0, sum_3264 = 0.0;
    for (std::size_t i = 0; i < rC.size(); ++i) {
        ASSERT_TRUE(std::isfinite(r64[i])) << "non-finite fp64-rung radiance";
        sum_c64 += std::abs(double(rC[i]) - double(r64[i]));
        sum_3264 += std::abs(double(r32[i]) - double(r64[i]));
    }
    const double mean_c64 = sum_c64 / n;
    const double mean_3264 = sum_3264 / n;
    std::cout << "[ traceC   ] mean|comp-fp64|=" << mean_c64 << " mean|fp32-fp64|=" << mean_3264
              << " mean|comp-fp32|=" << mean_c32 << "\n";
    EXPECT_LE(mean_c64, mean_3264 * 1.5 + 1e-9)
        << "compensation made the fp64 tracking worse, which defeats the rung";
#endif
}

#ifdef SIRIUS_KERNEL_DIR
// Direct production trace ABI. A one-pixel launch has seven valid bindings;
// the empty CSR covers every sky cell so the pupil-bundle diagnostic never
// relies on an out-of-bounds dummy index. No production timing/state hook.
std::array<float, sirius::render::kTraceParameterCount> FlatBoundaryTraceParams() {
    std::array<float, sirius::render::kTraceParameterCount> p{};
    p[0] = p[1] = p[33] = p[34] = 1.0f;
    p[7] = 3.0f;
    p[10] = p[14] = p[18] = 1.0f;
    p[19] = 0.5f;
    p[20] = 1.0f;
    p[21] = 128.0f;
    p[22] = 100.0f;
    p[23] = 0.001f;
    p[24] = 8.0f;
    p[25] = 5.0f;
    p[26] = 1.0f;
    p[28] = 6.0f;
    p[29] = 8.0f;
    p[30] = 1.0e7f;
    p[36] = p[37] = 1.0f;
    p[39] = p[40] = p[42] = 1.0f;
    p[44] = p[45] = p[66] = p[67] = 0.5f;
    p[50] = p[51] = 1.0f;
    p[53] = 50.0f;
    p[54] = 2.8f;
    p[55] = 10.0f;
    p[57] = 0.001f;
    p[58] = 1.0f;
    p[60] = 0.1f;
    p[61] = p[62] = 1.0f;
    p[63] = 64.0f;
    return p;
}

template <typename Real>
sirius::base::Expected<std::array<float, 4>> TracePixel(
    sirius::test::TraceContinuationProbe<Real>& probe,
    const std::array<float, sirius::render::kTraceParameterCount>& parameters) {
    auto result = probe.RunRegion(parameters);
    if (!result) return std::unexpected(result.error());
    if (result->size() != 4)
        return sirius::base::Fail(sirius::base::ErrorDomain::kKernel, "trace pixel",
                                  "wrong output extent");
    return std::array<float, 4>{(*result)[0], (*result)[1], (*result)[2], (*result)[3]};
}

#endif

TEST(KernelTrace, ActualTraceClipsJacobiAndRetainsInvalidSamplesAcrossRungs) {
#ifdef SIRIUS_KERNEL_DIR
    const auto devices = EnumerateVulkanDevices();
    ASSERT_TRUE(devices.has_value()) << devices.error().Description();
    if (devices->empty()) GTEST_SKIP() << "no Vulkan device present";
    const auto selected = ResolveVulkanDeviceIndex(*devices);
    ASSERT_TRUE(selected.has_value()) << selected.error().Description();
    auto opened = CreateVulkanDevice(*selected);
    ASSERT_TRUE(opened.has_value()) << opened.error().Description();
    auto& device = **opened;
    for (const char* name : {"trace.spv", "trace_fp32comp.spv", "trace_fp64.spv"}) {
        SCOPED_TRACE(name);
        const auto words = LoadSpirv(std::string(SIRIUS_KERNEL_DIR) + "/" + name);
        ASSERT_FALSE(words.empty());
        const auto kernel = device.LoadKernel(words);
        if (std::string(name) == "trace_fp64.spv" && !device.Info().supports_fp64) {
            ASSERT_FALSE(kernel.has_value());
            EXPECT_EQ(kernel.error().domain(), sirius::base::ErrorDomain::kKernel);
            EXPECT_NE(kernel.error().detail().find("shaderFloat64"), std::string::npos);
            RecordProperty("fp64_evidence", "unsupported_kernel_declined");
            continue;
        }
        ASSERT_TRUE(kernel.has_value());
        const auto run = [&]<typename Real>() {
            sirius::test::TraceContinuationProbe<Real> probe(device, *kernel);
            const auto prepared = probe.Prepare();
            ASSERT_TRUE(prepared.has_value()) << prepared.error().Description();
            // Independent Minkowski solution: X_direction=lambda*screen,
            // V_direction=screen, X_position=screen, V_position=0. Both labelled
            // pairs must survive regardless of which one feeds the beam image.
            const auto check_columns = [&](Real affine) {
                ASSERT_EQ(probe.Records().size(), 1U);
                const auto& record = probe.Records()[0];
                using Record = sirius::render::TraceContinuationRecord<Real>;
                const Real seed = Real(0.001f);
                const double tolerance = sizeof(Real) == 8 ? 1e-12 : 1e-8;
                const std::array<std::array<Real, 4>, 2> screens{
                    std::array<Real, 4>{0, 0, 0, -seed},
                    std::array<Real, 4>{0, 0, seed, 0}};
                EXPECT_NEAR(record.physical[Record::kIntegration][1], affine, tolerance);
                for (std::size_t column = 0; column < 4; ++column) {
                    for (std::size_t component = 0; component < 4; ++component) {
                        const Real screen = screens[column % 2][component];
                        const Real x = column < 2 ? affine * screen : screen;
                        const Real velocity = column < 2 ? screen : Real(0);
                        EXPECT_NEAR(record.physical[Record::kPositionColumns + column][component],
                                    x, tolerance) << column << component;
                        EXPECT_NEAR(record.physical[Record::kCovariantColumns + column][component],
                                    velocity, tolerance) << column << component;
                    }
                }
            };
            auto p = FlatBoundaryTraceParams();
            p[41] = p[56] = 1.0f;  // Pupil Jacobi state, empty actual catalogue.
            const auto large_step = TracePixel(probe, p);
            ASSERT_TRUE(large_step.has_value()) << large_step.error().Description();
            ASSERT_NO_FATAL_FAILURE(check_columns(Real(2)));
            p[24] = 0.75f;
            const auto small_step = TracePixel(probe, p);
            ASSERT_TRUE(small_step.has_value()) << small_step.error().Description();
            ASSERT_NO_FATAL_FAILURE(check_columns(Real(2)));
            EXPECT_GT(probe.Records()[0].control[3], 1U);
            p[56] = 0.0f;
            const auto parallel = TracePixel(probe, p);
            ASSERT_TRUE(parallel.has_value()) << parallel.error().Description();
            EXPECT_NEAR((*parallel)[3], 1.0f, 3.0e-5f);
            ASSERT_NO_FATAL_FAILURE(check_columns(Real(2)));
            p[56] = 1.0f;
            // In flat spacetime xi=lambda*epsilon at a translated observer. The
            // existing geometric alpha is lambda/R, not the angular-map Jacobian.
            // R=5, x_launch=3 => lambda=2. Overshot endpoints give a different ratio.
            for (const auto& rgba : {*large_step, *small_step}) {
                for (float value : rgba) EXPECT_TRUE(std::isfinite(value));
                EXPECT_NEAR(rgba[3], 2.0f / 5.0f, 3.0e-5f);
            }
            p[7] = 4.0f;
            p[24] = 8.0f;
            const auto translated = TracePixel(probe, p);
            ASSERT_TRUE(translated.has_value()) << translated.error().Description();
            EXPECT_NEAR((*translated)[3], 1.0f / 5.0f, 3.0e-5f);
            ASSERT_NO_FATAL_FAILURE(check_columns(Real(1)));
            p[7] = 5.0f;
            const auto at_surface = TracePixel(probe, p);
            ASSERT_TRUE(at_surface.has_value()) << at_surface.error().Description();
            for (float value : *at_surface) EXPECT_TRUE(std::isfinite(value));
            EXPECT_NEAR((*at_surface)[3], 0.0f, 3.0e-5f);
            ASSERT_NO_FATAL_FAILURE(check_columns(Real(0)));

            // An exterior outward launch has no represented source-sphere event.
            // The real typed failure must not publish a fabricated black/NaN sample.
            p[7] = 6.0f;
            const auto invalid = TracePixel(probe, p);
            ASSERT_FALSE(invalid.has_value());
            ASSERT_EQ(probe.Records().size(), 1U);
            EXPECT_EQ(probe.Records()[0].control[0],
                      static_cast<std::uint32_t>(sirius::render::TracePhase::Terminal));
            EXPECT_EQ(
                probe.Records()[0].control[1],
                static_cast<std::uint32_t>(sirius::render::TraceTermination::UnrepresentedEvent));
            const auto failed_record = probe.Records()[0];
            const auto unpublished = probe.ReadRadiance();
            ASSERT_TRUE(unpublished.has_value());
            for (const float value : *unpublished) EXPECT_FLOAT_EQ(value, 0.0f);
            const auto declined_finalize = probe.FinalizeFailedStateForControl();
            ASSERT_TRUE(declined_finalize.has_value()) << declined_finalize.error().Description();
            EXPECT_EQ(probe.Records()[0].control, failed_record.control);
            EXPECT_EQ(probe.Records()[0].identity, failed_record.identity);
            EXPECT_EQ(probe.Records()[0].physical, failed_record.physical);
            EXPECT_EQ(probe.Records()[0].source, failed_record.source);
            EXPECT_EQ(probe.Records()[0].volume, failed_record.volume);
            const auto after_finalize = probe.ReadRadiance();
            ASSERT_TRUE(after_finalize.has_value());
            EXPECT_EQ(*after_finalize, *unpublished);
            p[7] = 3.0f;
            for (int sample_index : {1, 2}) {
                p[46] = static_cast<float>(sample_index);
                const auto before_submissions = probe.PhysicalSubmissions();
                const auto retained = TracePixel(probe, p);
                EXPECT_FALSE(retained.has_value());
                EXPECT_EQ(probe.PhysicalSubmissions(), before_submissions)
                    << "failed sample batch must not dispatch a fresh later sample";
                const auto retained_output = probe.ReadRadiance();
                ASSERT_TRUE(retained_output.has_value());
                EXPECT_EQ(*retained_output, *unpublished);
            }
            // A new sample-zero accumulation must recover, proving both the
            // valid launch and prior typed failure controls actually executed.
            p[46] = 0.0f;
            const auto reset = TracePixel(probe, p);
            ASSERT_TRUE(reset.has_value()) << reset.error().Description();
            for (float value : *reset) EXPECT_TRUE(std::isfinite(value));
            EXPECT_NEAR((*reset)[3], 2.0f / 5.0f, 3.0e-5f);
            p[46] = 1.0f;
            const auto before_valid_sample = probe.PhysicalSubmissions();
            const auto continued = TracePixel(probe, p);
            ASSERT_TRUE(continued.has_value()) << continued.error().Description();
            EXPECT_GT(probe.PhysicalSubmissions(), before_valid_sample);
            EXPECT_EQ(*continued, *reset)
                << "same valid sample must average normally after an explicit restart";
        };
        if (std::string(name) == "trace_fp64.spv") {
            ASSERT_NO_FATAL_FAILURE(run.template operator()<double>());
        } else {
            ASSERT_NO_FATAL_FAILURE(run.template operator()<float>());
        }
    }
#else
    GTEST_SKIP() << "trace kernels not compiled";
#endif
}

TEST(KernelTrace, ActualTraceFiniteSphereExcludesLaterDiskAndVolumeAcrossRungs) {
#ifdef SIRIUS_KERNEL_DIR
    const auto devices = EnumerateVulkanDevices();
    ASSERT_TRUE(devices.has_value()) << devices.error().Description();
    if (devices->empty()) GTEST_SKIP() << "no Vulkan device present";
    const auto selected = ResolveVulkanDeviceIndex(*devices);
    ASSERT_TRUE(selected.has_value()) << selected.error().Description();
    auto opened = CreateVulkanDevice(*selected);
    ASSERT_TRUE(opened.has_value()) << opened.error().Description();
    auto& device = **opened;
    for (const char* name : {"trace.spv", "trace_fp32comp.spv", "trace_fp64.spv"}) {
        SCOPED_TRACE(name);
        const auto words = LoadSpirv(std::string(SIRIUS_KERNEL_DIR) + "/" + name);
        ASSERT_FALSE(words.empty());
        const auto kernel = device.LoadKernel(words);
        if (std::string(name) == "trace_fp64.spv" && !device.Info().supports_fp64) {
            ASSERT_FALSE(kernel.has_value());
            EXPECT_EQ(kernel.error().domain(), sirius::base::ErrorDomain::kKernel);
            EXPECT_NE(kernel.error().detail().find("shaderFloat64"), std::string::npos);
            RecordProperty("fp64_evidence", "unsupported_kernel_declined");
            continue;
        }
        ASSERT_TRUE(kernel.has_value());
        const auto run = [&]<typename Real>() {
            sirius::test::TraceContinuationProbe<Real> probe(device, *kernel);
            const auto prepared = probe.Prepare();
            ASSERT_TRUE(prepared.has_value()) << prepared.error().Description();
            auto p = FlatBoundaryTraceParams();
            // The material transfer authority requires positive M. This weak-field
            // Schwarzschild ray leaves R=5 before it can reach the r>=6 annulus,
            // crosses z=0 near x=7, and therefore meets material before R=12.
            p[3] = 1.0e-6f;
            p[9] = 1.0f;
            p[12] = -0.25f;
            p[16] = 0.25f;
            const auto short_clear = TracePixel(probe, p);
            ASSERT_TRUE(short_clear.has_value()) << short_clear.error().Description();
            p[25] = 12.0f;
            const auto long_clear = TracePixel(probe, p);
            ASSERT_TRUE(long_clear.has_value()) << long_clear.error().Description();
            for (int volume : {0, 1}) {
                SCOPED_TRACE(volume == 0 ? "opaque disk" : "volumetric disk");
                p[27] = 1.0f;
                p[59] = static_cast<float>(volume);
                p[25] = 5.0f;
                const auto short_material = TracePixel(probe, p);
                ASSERT_TRUE(short_material.has_value()) << short_material.error().Description();
                p[25] = 12.0f;
                const auto long_material = TracePixel(probe, p);
                ASSERT_TRUE(long_material.has_value()) << long_material.error().Description();
                float positive_control_difference = 0.0f;
                for (std::size_t channel = 0; channel < 3; ++channel) {
                    ASSERT_TRUE(std::isfinite((*short_material)[channel]));
                    ASSERT_TRUE(std::isfinite((*long_material)[channel]));
                    EXPECT_NEAR((*short_material)[channel], (*short_clear)[channel], 3.0e-5f);
                    positive_control_difference +=
                        std::abs((*long_material)[channel] - (*long_clear)[channel]);
                }
                EXPECT_GT(positive_control_difference, 0.01f)
                    << "larger sphere must expose actual material; disabling the consumer is not a "
                       "pass";
            }
        };
        if (std::string(name) == "trace_fp64.spv") {
            ASSERT_NO_FATAL_FAILURE(run.template operator()<double>());
        } else {
            ASSERT_NO_FATAL_FAILURE(run.template operator()<float>());
        }
    }
#else
    GTEST_SKIP() << "trace kernels not compiled";
#endif
}

}  // namespace
