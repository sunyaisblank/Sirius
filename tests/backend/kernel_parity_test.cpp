// Kernel-versus-core parity gate (docs/SPECIFICATION.md section 7,
// docs/ARCHITECTURE.md section 9). The parity_probe Slang kernel evaluates the
// live device physics at fixed sample points through the Vulkan adapter.
//
// Reference provenance: the values embedded below are the stdout of a
// mechanical extraction of the legacy device functions (scratchpad
// kernel-ref/ref_gen.cpp: __device__/template markers removed, make_float3 and
// sincosf shimmed, arithmetic otherwise verbatim, built with g++-14 -O2). The
// comparison is host-cmath-fp32 (reference) against Slang->SPIR-V->Lavapipe-fp32
// (kernel); tolerances below absorb the last-ULP differences between the two
// transcendental libraries and are named where they exceed the nominal 1e-6.
//
// Sample points (as the workstream brief requires): Kerr a=0.9 near the photon
// sphere, a=0 Schwarzschild at r=6M, an Ellis wormhole near its throat, and the
// Alcubierre warp bubble wall, plus a charged Kerr-Newman point for inverse
// coverage and two disk/blackbody temperatures.

#include "sirius/backend/device.h"
#include "sirius/core/disk/novikov_thorne_disk.h"
#include "sirius/core/spectral/blackbody.h"
#include "sirius/core/spectral/colour_modes.h"

#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <fstream>
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

constexpr std::uint32_t kSampleStride = 16;
constexpr std::uint32_t kResultStride = 64;

// Metric selector ids (must match gr_types.slang).
constexpr float kKerrSchild = 0.0f;
constexpr float kEllis = 1.0f;
constexpr float kWarp = 2.0f;

// Opcodes (must match parity_probe.slang).
constexpr std::uint32_t kOpMetric = 0;
constexpr std::uint32_t kOpChristoffel = 1;
constexpr std::uint32_t kOpDiskTemp = 2;
constexpr std::uint32_t kOpBlackbody = 3;
constexpr std::uint32_t kOpDeviation = 4;
constexpr std::uint32_t kOpDiskRedshift = 5;
constexpr std::uint32_t kOpLiveCartConservation = 6;
constexpr std::uint32_t kOpBeamEllipse = 7;

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

// One sample record, laid out per the parity_probe ABI.
struct Sample {
    float metric_id = 0.0f;
    float p1 = 0.0f, p2 = 0.0f, p3 = 0.0f, p4 = 0.0f;  // family params
    float c0 = 0.0f, c1 = 0.0f, c2 = 0.0f;             // position coords
    float u0 = 0.0f, u1 = 0.0f, u2 = 0.0f, u3 = 0.0f;  // 4-velocity
    float h = 0.0f, aux0 = 0.0f, aux1 = 0.0f, aux2 = 0.0f;
};

std::vector<float> Flatten(const std::vector<Sample>& samples) {
    std::vector<float> flat;
    flat.reserve(samples.size() * kSampleStride);
    for (const auto& s : samples) {
        flat.insert(flat.end(), {s.metric_id, s.p1, s.p2, s.p3, s.p4, s.c0, s.c1, s.c2, s.u0, s.u1,
                                 s.u2, s.u3, s.h, s.aux0, s.aux1, s.aux2});
    }
    return flat;
}

// Runs the probe for one opcode over a set of samples, returns the flat result
// buffer (kResultStride floats per sample).
std::vector<float> RunProbe(ComputeDevice& device, KernelHandle kernel, std::uint32_t opcode,
                            const std::vector<Sample>& samples) {
    const auto count = static_cast<std::uint32_t>(samples.size());
    const std::vector<float> flat = Flatten(samples);
    const std::vector<float> params = {static_cast<float>(opcode), static_cast<float>(count),
                                       static_cast<float>(kSampleStride),
                                       static_cast<float>(kResultStride)};
    std::vector<float> results(count * kResultStride, 0.0f);

    const auto sbuf = device.CreateBuffer(flat.size() * sizeof(float), BufferUsage::kStorage);
    const auto rbuf = device.CreateBuffer(results.size() * sizeof(float), BufferUsage::kStorage);
    const auto pbuf = device.CreateBuffer(params.size() * sizeof(float), BufferUsage::kStorage);
    EXPECT_TRUE(sbuf && rbuf && pbuf);

    EXPECT_TRUE(device.WriteBuffer(*sbuf, std::as_bytes(std::span<const float>(flat))));
    // Zero the result buffer so unused slots read back deterministically.
    EXPECT_TRUE(device.WriteBuffer(*rbuf, std::as_bytes(std::span<const float>(results))));
    EXPECT_TRUE(device.WriteBuffer(*pbuf, std::as_bytes(std::span<const float>(params))));

    const BufferHandle binding[] = {*sbuf, *rbuf, *pbuf};
    EXPECT_TRUE(device.Dispatch(kernel, binding, (count + 63) / 64, 1, 1).has_value());

    EXPECT_TRUE(device.ReadBuffer(*rbuf, std::as_writable_bytes(std::span<float>(results))));
    return results;
}

// Relative-with-floor closeness: |a-e| <= rel*|e| + abs.
::testing::AssertionResult Close(float a, float e, float rel, float abs_tol,
                                 const std::string& tag) {
    const float d = std::abs(a - e);
    if (d <= rel * std::abs(e) + abs_tol) {
        return ::testing::AssertionSuccess();
    }
    return ::testing::AssertionFailure()
           << tag << " actual=" << a << " expected=" << e << " |diff|=" << d;
}

struct Fixture {
    std::unique_ptr<ComputeDevice> device;
    KernelHandle kernel;
    bool ready = false;
};

Fixture OpenProbe() {
    Fixture f;
#ifndef SIRIUS_KERNEL_DIR
    return f;
#else
    const auto devices = EnumerateVulkanDevices();
    if (!devices.has_value() || devices->empty()) {
        return f;
    }
    const auto selected = ResolveVulkanDeviceIndex(*devices);
    if (!selected.has_value()) {
        return f;
    }
    auto device = CreateVulkanDevice(*selected);
    if (!device.has_value()) {
        return f;
    }
    const auto spirv = LoadSpirv(std::string(SIRIUS_KERNEL_DIR) + "/parity_probe.spv");
    if (spirv.empty()) {
        return f;
    }
    auto kernel = (*device)->LoadKernel(spirv);
    if (!kernel.has_value()) {
        return f;
    }
    f.device = std::move(*device);
    f.kernel = *kernel;
    f.ready = true;
    return f;
#endif
}

// --- Sample points -----------------------------------------------------------

Sample KerrSchildAt(float M, float a, float Q, float x, float y, float z) {
    Sample s;
    s.metric_id = kKerrSchild;
    s.p1 = M;
    s.p2 = a;
    s.p3 = Q;
    s.p4 = 0.0f;
    s.c0 = x;
    s.c1 = y;
    s.c2 = z;
    return s;
}

// S1 Kerr a=0.9 near photon sphere; S2 Schwarzschild r=6M; S3 Kerr-Newman.
const Sample kS1 = KerrSchildAt(1.0f, 0.9f, 0.0f, 1.6f, 0.8f, 0.3f);
const Sample kS2 = KerrSchildAt(1.0f, 0.0f, 0.0f, 4.0f, 2.0f, 4.0f);
const Sample kS3 = KerrSchildAt(1.0f, 0.9f, 0.4f, 2.0f, 1.0f, 0.6f);

// --- Embedded reference values (ref_gen.cpp stdout) --------------------------

const std::array<float, 16> kS1_g = {0.24817276f,    1.22385824f,   -0.0649836585f, 0.236396417f,
                                     1.22385824f,    2.20001745f,   -0.0637177676f, 0.231791392f,
                                     -0.0649836585f, -0.063717775f, 1.00338328f,    -0.0123075135f,
                                     0.236396417f,   0.231791407f,  -0.0123075135f, 1.04477203f};
const std::array<float, 16> kS1_ginv = {
    -2.24817276f,  1.22385824f,   -0.0649836585f, 0.236396417f, 1.22385824f,  -0.200017452f,
    0.0637177676f, -0.231791392f, -0.0649836585f, 0.063717775f, 0.996616781f, 0.0123075135f,
    0.236396417f,  -0.231791407f, 0.0123075135f,  0.955227971f};
const std::array<float, 16> kS2_g = {-0.666666627f, 0.222222239f,  0.111111119f,  0.222222239f,
                                     0.222222239f,  1.14814818f,   0.0740740821f, 0.148148164f,
                                     0.111111119f,  0.0740740821f, 1.03703701f,   0.0740740821f,
                                     0.222222239f,  0.148148164f,  0.0740740821f, 1.14814818f};
const std::array<float, 16> kS2_ginv = {
    -1.33333337f,   0.222222239f,  0.111111119f,   0.222222239f,   0.222222239f, 0.851851821f,
    -0.0740740821f, -0.148148164f, 0.111111119f,   -0.0740740821f, 0.962962985f, -0.0740740821f,
    0.222222239f,   -0.148148164f, -0.0740740821f, 0.851851821f};
const std::array<float, 16> kS3_g = {-0.1156317f,   0.84726429f,  0.0567223616f, 0.247048855f,
                                     0.84726429f,   1.81171703f,  0.054342553f,  0.236683816f,
                                     0.0567223616f, 0.054342553f, 1.00363815f,   0.0158454273f,
                                     0.247048855f,  0.236683816f, 0.0158454273f, 1.06901324f};
const std::array<float, 16> kS3_ginv = {
    -1.8843683f,   0.84726429f,   0.0567223616f,  0.247048855f,  0.84726429f,  0.188283026f,
    -0.054342553f, -0.236683816f, 0.0567223616f,  -0.054342553f, 0.996361911f, -0.0158454273f,
    0.247048855f,  -0.236683816f, -0.0158454273f, 0.930986762f};

const std::array<float, 16> kS4_g = {
    -0.186177671f, -0.902121007f, 0.0f,        0.0f, -0.902121007f, 1.0f, 0.0f, 0.0f,
    0.0f,          0.0f,          4.08212233f, 0.0f, 0.0f,          0.0f, 0.0f, 4.08212233f};
const std::array<float, 16> kS4_ginv = {
    -1.0f, -0.902121007f, 0.0f,         0.0f, -0.902121007f, 0.186177671f, 0.0f, 0.0f,
    0.0f,  0.0f,          0.244970605f, 0.0f, 0.0f,          0.0f,         0.0f, 0.244970605f};
const std::array<float, 16> kS5_g = {4.76837158e-07f,
                                     -1.00000024f,
                                     0.0f,
                                     0.0f,
                                     -1.00000024f,
                                     1.0f,
                                     0.0f,
                                     0.0f,
                                     0.0f,
                                     0.0f,
                                     1.0f,
                                     0.0f,
                                     0.0f,
                                     0.0f,
                                     0.0f,
                                     1.0f};
const std::array<float, 16> kS5_ginv = {
    -1.0f, -1.00000024f, 0.0f, 0.0f, -1.00000024f, -4.76837158e-07f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f,
    0.0f,  0.0f,         0.0f, 0.0f, 1.0f};

const std::array<float, 64> kS1_christoffel = {
    0.480513632f,   0.846556902f,    0.162684917f,   0.231745839f,    0.846556902f,
    1.16810358f,    0.100205779f,    0.442988515f,   0.162684917f,    0.100205779f,
    -0.795098364f,  0.0158032775f,   0.231745839f,   0.442988515f,    0.0158032775f,
    -0.680479527f,  -0.0957495421f,  -0.461975127f,  -0.163184449f,   -0.151768044f,
    -0.461975127f,  -0.784428596f,   -0.101850964f,  -0.360365629f,   -0.163184449f,
    -0.101850964f,  0.778974175f,    -0.0164174605f, -0.151768044f,   -0.360365629f,
    -0.0164174605f, 0.682342529f,    0.21271883f,    0.212243751f,    -0.0013024573f,
    0.129806966f,   0.212243751f,    0.210141778f,   -0.00353839644f, 0.135504469f,
    -0.0013024573f, -0.00353839644f, -0.0408864506f, -0.0053072162f,  0.129806966f,
    0.135504469f,   -0.0053072162f,  0.00243838085f, 0.0497330427f,   -0.0266991332f,
    -0.120330833f,  -0.0172360651f,  -0.0266991332f, -0.094480522f,   -0.106526606f,
    -0.0585899353f, -0.120330833f,   -0.106526606f,  0.159526646f,    -0.0199474562f,
    -0.0172360651f, -0.0585899353f,  -0.0199474562f, 0.133927062f};
const std::array<float, 64> kS2_christoffel = {
    0.00925926026f,  0.0246913582f,   0.0123456791f,   0.0246913582f,   0.0246913582f,
    -0.00205761427f, 0.0267489739f,   0.0534979478f,   0.0123456791f,   0.0267489739f,
    -0.0421810746f,  0.0267489739f,   0.0246913582f,   0.0534979478f,   0.0267489739f,
    -0.00205761427f, 0.0123456791f,   -0.00411522668f, -0.00205761334f, -0.00411522668f,
    -0.00411522668f, 0.00960219558f,  -0.0137174213f,  -0.0274348427f,  -0.00205761334f,
    -0.0137174213f,  0.0301783271f,   -0.0137174213f,  -0.00411522668f, -0.0274348427f,
    -0.0137174213f,  0.00960219558f,  0.00617283955f,  -0.00205761334f, -0.00102880667f,
    -0.00205761334f, -0.00205761334f, 0.00480109872f,  -0.00685871206f, -0.0137174223f,
    -0.00102880667f, -0.00685871206f, 0.0150891654f,   -0.00685871113f, -0.00205761334f,
    -0.0137174223f,  -0.00685871113f, 0.00480109872f,  0.0123456791f,   -0.00411522668f,
    -0.00205761334f, -0.00411522668f, -0.00411522668f, 0.00960219651f,  -0.0137174213f,
    -0.0274348427f,  -0.00205761334f, -0.0137174213f,  0.0301783271f,   -0.0137174213f,
    -0.00411522668f, -0.0274348427f,  -0.0137174213f,  0.00960219651f};
const std::array<float, 64> kS3_christoffel = {
    0.17010279f,      0.334537357f,   0.0966958255f,  0.127950162f,   0.334537357f,
    0.451506615f,     0.128602356f,   0.279217303f,   0.0966958255f,  0.128602356f,
    -0.392807484f,    0.0394485183f,  0.127950162f,   0.279217303f,   0.0394485183f,
    -0.316273749f,    0.00860523805f, -0.156128749f,  -0.0683526099f, -0.0777029321f,
    -0.156128749f,    -0.27508682f,   -0.0999394208f, -0.224506527f,  -0.0683526099f,
    -0.0999394208f,   0.378736645f,   -0.0312046371f, -0.0777029321f, -0.224506527f,
    -0.0312046371f,   0.314689457f,   0.0748754591f,  0.047447715f,   -0.000699766446f,
    0.0613086782f,    0.047447715f,   0.0243298467f,  -0.0038289465f, 0.0449797213f,
    -0.000699766446f, -0.0038289465f, 0.0255471226f,  0.00192844577f, 0.0613086782f,
    0.0449797213f,    0.00192844577f, 0.0524292588f,  0.0329135135f,  -0.0133464206f,
    -0.0674042106f,   -0.0132742766f, -0.0133464206f, -0.0464611873f, -0.0744270533f,
    -0.0556216836f,   -0.0674042106f, -0.0744270533f, 0.104218721f,   -0.0223035496f,
    -0.0132742766f,   -0.0556216836f, -0.0223035496f, 0.0946279168f};

// -----------------------------------------------------------------------------

TEST(KernelParity, KerrSchildMetricMatchesLegacyToOnePartInMillion) {
    Fixture f = OpenProbe();
    if (!f.ready) GTEST_SKIP() << "no Vulkan device or kernels absent";

    const std::vector<Sample> samples = {kS1, kS2, kS3};
    const auto r = RunProbe(*f.device, f.kernel, kOpMetric, samples);
    const std::array<const std::array<float, 16>*, 3> gref = {&kS1_g, &kS2_g, &kS3_g};
    const std::array<const std::array<float, 16>*, 3> giref = {&kS1_ginv, &kS2_ginv, &kS3_ginv};

    for (std::size_t s = 0; s < samples.size(); ++s) {
        const std::size_t base = s * kResultStride;
        for (int k = 0; k < 16; ++k) {
            EXPECT_TRUE(Close(r[base + k], (*gref[s])[k], 1e-6f, 1e-6f,
                              "g[" + std::to_string(s) + "][" + std::to_string(k) + "]"));
            EXPECT_TRUE(Close(r[base + 16 + k], (*giref[s])[k], 1e-6f, 1e-6f,
                              "ginv[" + std::to_string(s) + "][" + std::to_string(k) + "]"));
        }
        // The analytic Kerr-Schild inverse satisfies g g^-1 = I to fp32.
        EXPECT_LT(r[base + 32], 5e-6f) << "inverse-metric identity error, sample " << s;
    }
}

TEST(KernelParity, WormholeAndWarpMetricMatchLegacy) {
    Fixture f = OpenProbe();
    if (!f.ready) GTEST_SKIP() << "no Vulkan device or kernels absent";

    Sample ellis;  // Ellis drainhole m=0.5, n=1.0, near throat (r=0.6)
    ellis.metric_id = kEllis;
    ellis.p1 = 0.5f;
    ellis.p2 = 1.0f;
    ellis.c0 = 0.6f;
    ellis.c1 = 1.5707963f;
    ellis.c2 = 0.0f;

    Sample warp;  // Alcubierre bubble wall: vs=2, sigma=8, R=1 at rs=1
    warp.metric_id = kWarp;
    warp.p1 = 2.0f;
    warp.p2 = 8.0f;
    warp.p3 = 1.0f;
    warp.c0 = 1.0f;
    warp.c1 = 0.0f;
    warp.c2 = 0.0f;

    const std::vector<Sample> samples = {ellis, warp};
    const auto r = RunProbe(*f.device, f.kernel, kOpMetric, samples);
    const std::array<const std::array<float, 16>*, 2> gref = {&kS4_g, &kS5_g};
    const std::array<const std::array<float, 16>*, 2> giref = {&kS4_ginv, &kS5_ginv};

    // Loosened from the nominal 1e-6 to 1e-5 relative: the Ellis areal-radius
    // factor Rp/sqrt(1 - Fp^2) threads atan -> exp -> sqrt, and the warp shape
    // function threads tanh, so their angular components accumulate ~1-2 ULP per
    // transcendental of disagreement between Lavapipe and host glibc (observed
    // worst case 1.4e-6 relative on g_thth). The Kerr-Schild family, which is
    // algebraic apart from one sqrt for the radius, holds the strict 1e-6.
    for (std::size_t s = 0; s < samples.size(); ++s) {
        const std::size_t base = s * kResultStride;
        for (int k = 0; k < 16; ++k) {
            EXPECT_TRUE(Close(r[base + k], (*gref[s])[k], 1e-5f, 1e-6f,
                              "g[" + std::to_string(s) + "][" + std::to_string(k) + "]"));
            EXPECT_TRUE(Close(r[base + 16 + k], (*giref[s])[k], 1e-5f, 1e-6f,
                              "ginv[" + std::to_string(s) + "][" + std::to_string(k) + "]"));
        }
    }
}

TEST(KernelParity, KerrSchildChristoffelMatchesLegacyToOnePartInMillion) {
    Fixture f = OpenProbe();
    if (!f.ready) GTEST_SKIP() << "no Vulkan device or kernels absent";

    const std::vector<Sample> samples = {kS1, kS2, kS3};
    const auto r = RunProbe(*f.device, f.kernel, kOpChristoffel, samples);
    const std::array<const std::array<float, 64>*, 3> cref = {&kS1_christoffel, &kS2_christoffel,
                                                              &kS3_christoffel};

    for (std::size_t s = 0; s < samples.size(); ++s) {
        const std::size_t base = s * kResultStride;
        for (int k = 0; k < 64; ++k) {
            EXPECT_TRUE(Close(r[base + k], (*cref[s])[k], 1e-6f, 1e-6f,
                              "Gamma[" + std::to_string(s) + "][" + std::to_string(k) + "]"));
        }
    }
}

TEST(KernelParity, FullPageThorneDiskTemperatureMatchesIndependentCoreModel) {
    Fixture f = OpenProbe();
    if (!f.ready) GTEST_SKIP() << "no Vulkan device or kernels absent";

    Sample s8;  // Kerr a=0.9 at r=8; r_isco and inner temperature in aux slots
    s8.p1 = 1.0f;
    s8.p2 = 0.9f;
    s8.c0 = 8.0f;
    s8.aux0 = 2.32088304f;  // r_isco = ComputeKerrISCO(1, 0.9)
    s8.aux1 = 10000.0f;

    Sample s9;  // Schwarzschild at r=10
    s9.p1 = 1.0f;
    s9.p2 = 0.0f;
    s9.c0 = 10.0f;
    s9.aux0 = 6.0f;
    s9.aux1 = 10000.0f;

    const float retrograde_isco =
        static_cast<float>(sirius::core::AccretionDiskD::ComputeIsco(-0.7));
    Sample s10;  // Retrograde Kerr at dimensionless r/M=10 with M=2
    s10.p1 = 2.0f;
    s10.p2 = -0.7f;
    s10.c0 = 20.0f;
    s10.aux0 = 2.0f * retrograde_isco;
    s10.aux1 = 10000.0f;

    const std::vector<Sample> samples = {s8, s9, s10};
    const auto r = RunProbe(*f.device, f.kernel, kOpDiskTemp, samples);

    const auto expected_temperature = [](double spin, double radius) {
        sirius::core::AccretionDiskD::Config config;
        config.M = 1.0;
        config.a_star = spin;
        sirius::core::AccretionDiskD disk(config);
        return 10000.0 * disk.Temperature(radius) / disk.Temperature(1.5 * disk.IscoRadius());
    };

    // The kernel is fp32 quadrature while the independent Core model is double.
    EXPECT_TRUE(Close(r[0 * kResultStride + 0], static_cast<float>(expected_temperature(0.9, 8.0)),
                      2e-4f, 2e-3f, "S8 T"));
    EXPECT_GT(r[0 * kResultStride + 1], 0.0f);
    EXPECT_TRUE(Close(r[0 * kResultStride + 2], 2.32088304f, 1e-4f, 1e-6f, "S8 isco"));
    EXPECT_TRUE(Close(r[1 * kResultStride + 0], static_cast<float>(expected_temperature(0.0, 10.0)),
                      2e-4f, 2e-3f, "S9 T"));
    EXPECT_GT(r[1 * kResultStride + 1], 0.0f);
    EXPECT_TRUE(Close(r[1 * kResultStride + 2], 6.0f, 1e-4f, 1e-6f, "S9 isco"));
    EXPECT_TRUE(Close(r[2 * kResultStride + 0],
                      static_cast<float>(expected_temperature(-0.7, 10.0)), 2e-4f, 2e-3f, "S10 T"));
    EXPECT_GT(r[2 * kResultStride + 1], 0.0f);
    EXPECT_TRUE(Close(r[2 * kResultStride + 2], 2.0f * retrograde_isco, 1e-4f, 1e-6f, "S10 isco"));
}

TEST(KernelParity, BlackbodyMatchesIntegratedCoreSpectrum) {
    Fixture f = OpenProbe();
    if (!f.ready) GTEST_SKIP() << "no Vulkan device or kernels absent";

    Sample t500, t3000, t6000, t15000, t1000000;
    t500.c0 = 500.0f;
    t3000.c0 = 3000.0f;
    t6000.c0 = 6000.0f;
    t15000.c0 = 15000.0f;
    t1000000.c0 = 1000000.0f;

    const std::vector<Sample> samples = {t500, t3000, t6000, t15000, t1000000};
    const auto r = RunProbe(*f.device, f.kernel, kOpBlackbody, samples);

    for (std::size_t s = 0; s < samples.size(); ++s) {
        const auto expected = sirius::core::spectral::BlackbodyToRgb(samples[s].c0);
        const std::array<float, 3> ref = {expected.r, expected.g, expected.b};
        const std::size_t base = s * kResultStride;
        for (int k = 0; k < 3; ++k) {
            EXPECT_TRUE(Close(r[base + k], ref[k], 2e-3f, 2e-4f,
                              "bb[" + std::to_string(s) + "][" + std::to_string(k) + "]"));
        }
    }
}

TEST(KernelParity, DiskEmissionAppliesExactlyOneGFourthFactor) {
    Fixture f = OpenProbe();
    if (!f.ready) GTEST_SKIP() << "no Vulkan device or kernels absent";

    Sample receding;
    receding.c0 = 0.37f;
    receding.c1 = 0.5f;
    Sample approaching;
    approaching.c0 = 0.37f;
    approaching.c1 = 1.25f;
    const std::vector<Sample> samples = {receding, approaching};
    const auto results = RunProbe(*f.device, f.kernel, kOpDiskRedshift, samples);

    for (std::size_t i = 0; i < samples.size(); ++i) {
        const float expected =
            sirius::core::color_modes::ObservedBolometricIntensity(samples[i].c0, samples[i].c1);
        EXPECT_FLOAT_EQ(results[i * kResultStride], expected);
    }
}

TEST(KernelParity, NearExtremalKerrLiveRenderIntegratorConservesEnergyAngularMomentumAndCarter) {
    Fixture f = OpenProbe();
    if (!f.ready) GTEST_SKIP() << "no Vulkan device or kernels absent";

    // Same off-plane a=0.998 launch used by the independent CPU live-path
    // conservation gate. The shader supplies the future-directed null root,
    // then advances the ray with the production Vulkan render schedule.
    Sample near_extremal;
    near_extremal.metric_id = kKerrSchild;
    near_extremal.p1 = 1.0f;
    near_extremal.p2 = 0.998f;
    near_extremal.c0 = 12.0f;
    near_extremal.c1 = 0.0f;
    near_extremal.c2 = 3.0f;
    near_extremal.u1 = 0.08f;
    near_extremal.u2 = 1.0f;
    near_extremal.u3 = -0.03f;
    near_extremal.h = 0.08f;       // production stepScale
    near_extremal.aux0 = 3000.0f;  // production maxSteps
    near_extremal.aux1 = 0.02f;    // production minStep
    near_extremal.aux2 = 2.0f;     // production maxStep

    const auto results = RunProbe(*f.device, f.kernel, kOpLiveCartConservation, {near_extremal});
    for (int component = 0; component < 10; ++component) {
        EXPECT_TRUE(std::isfinite(results[component]))
            << "conservation result component " << component << " is non-finite";
    }
    EXPECT_GT(results[4], 100.0f) << "render integrator made too little progress";
    EXPECT_FLOAT_EQ(results[9], 2.0f)
        << "render integrator exhausted its work bound instead of escaping";
    EXPECT_GE(results[5], 200.0f) << "render integrator did not reach the escape surface";
    EXPECT_GT(std::abs(results[6]), 1.0e-6f) << "degenerate initial energy";
    EXPECT_GT(std::abs(results[7]), 1.0e-6f) << "degenerate initial angular momentum";
    EXPECT_GT(std::abs(results[8]), 1.0e-6f) << "degenerate initial Carter constant";
    EXPECT_LT(results[0], 1.0e-4f) << "energy drift";
    EXPECT_LT(results[1], 1.0e-4f) << "axial angular-momentum drift";
    EXPECT_LT(results[2], 1.0e-4f) << "Carter-constant drift";
    EXPECT_LT(results[3], 5.0e-4f) << "relative null residual";
}

TEST(KernelParity, BeamEllipseRetainsBothAxesAndOutputOrientation) {
    Fixture f = OpenProbe();
    if (!f.ready) GTEST_SKIP() << "no Vulkan device or kernels absent";

    // J = R(0.4) diag(3, 1) R(-0.2). The right rotation changes only the
    // input basis; the output ellipse must retain axes 3:1 at +0.4 radians.
    constexpr float output_angle = 0.4f;
    constexpr float input_angle = -0.2f;
    const float co = std::cos(output_angle);
    const float so = std::sin(output_angle);
    const float ci = std::cos(input_angle);
    const float si = std::sin(input_angle);
    Sample map;
    map.c0 = 3.0f * co * ci - so * si;
    map.c1 = -3.0f * co * si - so * ci;
    map.c2 = 3.0f * so * ci + co * si;
    map.u0 = -3.0f * so * si + co * ci;

    const auto results = RunProbe(*f.device, f.kernel, kOpBeamEllipse, {map});
    EXPECT_NEAR(results[0], 3.0f, 2.0e-6f);
    EXPECT_NEAR(results[1], 1.0f, 2.0e-6f);
    EXPECT_NEAR(std::cos(2.0f * results[2]), std::cos(2.0f * output_angle), 2.0e-6f);
    EXPECT_NEAR(std::sin(2.0f * results[2]), std::sin(2.0f * output_angle), 2.0e-6f);
    EXPECT_NEAR(results[3], 3.0f, 2.0e-6f);
}

TEST(KernelParity, GeodesicDeviationIsFiniteAndCurvedNearBlackHole) {
    Fixture f = OpenProbe();
    if (!f.ready) GTEST_SKIP() << "no Vulkan device or kernels absent";

    // Kerr a=0.9 at a Cartesian point ~4M out, with an infalling null tangent.
    // There is no legacy reference for this Cartesian-chart deviation (the .cu
    // evaluated it through a spherical<->Cartesian transform); the check is that
    // the full-Riemann acceleration is finite and non-zero, exercising the sole
    // represented deviation authority through the gate.
    Sample s;
    s.metric_id = kKerrSchild;
    s.p1 = 1.0f;
    s.p2 = 0.9f;
    s.p3 = 0.0f;
    s.p4 = 0.0f;
    s.c0 = 4.0f;
    s.c1 = 1.0f;
    s.c2 = 0.5f;
    s.u0 = 1.0f;
    s.u1 = -1.0f;
    s.u2 = 0.0f;
    s.u3 = 0.0f;
    s.aux0 = 0.05f;

    const std::vector<Sample> samples = {s};
    const auto r = RunProbe(*f.device, f.kernel, kOpDeviation, samples);
    for (int k = 0; k < 9; ++k) {
        EXPECT_TRUE(std::isfinite(r[k])) << "deviation component " << k << " non-finite";
    }
    EXPECT_GT(r[8], 1e-6f) << "full-Riemann tidal acceleration vanishes near the hole";
}

}  // namespace
