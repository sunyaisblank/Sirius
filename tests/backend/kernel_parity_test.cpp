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
#include "sirius/core/celestial_tangent_basis.h"
#include "sirius/core/cie1931_observer.h"
#include "sirius/core/coordinates.h"
#include "sirius/core/disk/novikov_thorne_disk.h"
#include "sirius/core/disk/volumetric_disk.h"
#include "sirius/core/metrics/kerr_schild_family.h"
#include "sirius/core/metrics/morris_thorne_family.h"
#include "sirius/core/relativistic_transfer.h"
#include "sirius/core/spectral/blackbody.h"
#include "sirius/core/spectral/colour_modes.h"
#include "sirius/core/xyz_srgb.h"
#include "sirius/oracle/kerr_boyer_lindquist.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <fstream>
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

constexpr std::uint32_t kSampleStride = 16;
constexpr std::uint32_t kResultStride = 64;

// Metric selector ids (must match gr_types.slang).
constexpr float kKerrSchild = 0.0f;
constexpr float kEllis = 1.0f;
constexpr float kWarp = 2.0f;
constexpr float kMorrisThorne = 3.0f;

// Opcodes (must match parity_probe.slang).
constexpr std::uint32_t kOpMetric = 0;
constexpr std::uint32_t kOpChristoffel = 1;
constexpr std::uint32_t kOpDiskTemp = 2;
constexpr std::uint32_t kOpBlackbody = 3;
constexpr std::uint32_t kOpDeviation = 4;
constexpr std::uint32_t kOpDiskRedshift = 5;
constexpr std::uint32_t kOpLiveCartConservation = 6;
constexpr std::uint32_t kOpBeamEllipse = 7;
constexpr std::uint32_t kOpAdaptiveEventDomain = 8;
constexpr std::uint32_t kOpVolumeOpacity = 9;
constexpr std::uint32_t kOpNullProjection = 10;
constexpr std::uint32_t kOpCelestialTangentBasis = 11;
constexpr std::uint32_t kOpJacobiRadialCongruence = 12;
constexpr std::uint32_t kOpXyzD65ToLinearSrgb = 13;
constexpr std::uint32_t kOpCie1931TwoDegreeFit = 14;
constexpr std::uint32_t kOpKerrZamoTransfer = 15;
constexpr std::uint32_t kOpKerrDiskTransfer = 16;
constexpr std::uint32_t kOpGreyLayerAbsorption = 17;
constexpr std::uint32_t kOpSphericalCaptureEvent = 18;
constexpr std::uint32_t kOpEllisTwoSheetTrace = 19;

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

bool HasSpirvCapability(const std::vector<std::uint32_t>& words, std::uint32_t capability) {
    // SPIR-V header is five words. OpCapability is opcode 17 and carries one
    // capability operand; Float64 is capability 10.
    for (std::size_t at = 5; at < words.size();) {
        const std::uint32_t instruction = words[at];
        const std::uint32_t word_count = instruction >> 16;
        const std::uint32_t opcode = instruction & 0xffffu;
        if (word_count == 0 || at + word_count > words.size()) {
            return false;
        }
        if (opcode == 17u && word_count >= 2u && words[at + 1] == capability) {
            return true;
        }
        at += word_count;
    }
    return false;
}

// One sample record, laid out per the parity_probe ABI.
struct Sample {
    float metric_id = 0.0f;
    float p1 = 0.0f, p2 = 0.0f, p3 = 0.0f, p4 = 0.0f;  // family params
    float c0 = 0.0f, c1 = 0.0f, c2 = 0.0f;             // position coords
    float u0 = 0.0f, u1 = 0.0f, u2 = 0.0f, u3 = 0.0f;  // 4-velocity
    float h = 0.0f, aux0 = 0.0f, aux1 = 0.0f, aux2 = 0.0f;
};

std::array<double, 4> SchwarzschildBlToKerrSchildVector(const sirius::oracle::Vec4d& vector,
                                                        const sirius::oracle::Vec4d& event,
                                                        double mass) {
    const double sin_theta = std::sin(event.theta);
    const double cos_theta = std::cos(event.theta);
    const double sin_phi = std::sin(event.phi);
    const double cos_phi = std::cos(event.phi);
    return {
        vector.t + (2.0 * mass / (event.r - 2.0 * mass)) * vector.r,
        sin_theta * cos_phi * vector.r + event.r * cos_theta * cos_phi * vector.theta -
            event.r * sin_theta * sin_phi * vector.phi,
        sin_theta * sin_phi * vector.r + event.r * cos_theta * sin_phi * vector.theta +
            event.r * sin_theta * cos_phi * vector.phi,
        cos_theta * vector.r - event.r * sin_theta * vector.theta,
    };
}

sirius::oracle::Vec4d KerrSchildToSchwarzschildBlVector(const std::array<float, 4>& vector,
                                                        const sirius::oracle::Vec4d& event,
                                                        double mass) {
    const double sin_theta = std::sin(event.theta);
    const double cos_theta = std::cos(event.theta);
    const double sin_phi = std::sin(event.phi);
    const double cos_phi = std::cos(event.phi);
    const double radial =
        sin_theta * cos_phi * vector[1] + sin_theta * sin_phi * vector[2] + cos_theta * vector[3];
    const double polar = (cos_theta * cos_phi * vector[1] + cos_theta * sin_phi * vector[2] -
                          sin_theta * vector[3]) /
                         event.r;
    const double azimuthal = (-sin_phi * vector[1] + cos_phi * vector[2]) / (event.r * sin_theta);
    return {vector[0] - (2.0 * mass / (event.r - 2.0 * mass)) * radial, radial, polar, azimuthal};
}

sirius::oracle::Vec4d AnalyticTidalAcceleration(sirius::oracle::KerrMetricD& metric,
                                                const sirius::oracle::Vec4d& event,
                                                const sirius::oracle::Vec4d& tangent,
                                                const sirius::oracle::Vec4d& deviation) {
    double riemann[4][4][4][4];
    metric.Riemann(event, riemann);
    sirius::oracle::Vec4d acceleration;
    for (int mu = 0; mu < 4; ++mu) {
        double contraction = 0.0;
        for (int nu = 0; nu < 4; ++nu)
            for (int rho = 0; rho < 4; ++rho)
                for (int sigma = 0; sigma < 4; ++sigma)
                    contraction +=
                        riemann[mu][nu][rho][sigma] * tangent[nu] * deviation[rho] * tangent[sigma];
        acceleration[mu] = -contraction;
    }
    return acceleration;
}

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

Fixture OpenProbe(const std::string& artefact = "parity_probe.spv") {
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
    const auto spirv = LoadSpirv(std::string(SIRIUS_KERNEL_DIR) + "/" + artefact);
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

TEST(KernelParity, RepresentedSubThresholdKerrMetricIsScaleCovariant) {
    Fixture f = OpenProbe();
    if (!f.ready) GTEST_SKIP() << "no Vulkan device or kernels absent";

    constexpr float scale = 1.0e-9f;
    const Sample unit = KerrSchildAt(1.0f, 0.05f, 0.0f, 2.0f, 1.0f, 0.5f);
    const Sample tiny =
        KerrSchildAt(scale, 0.05f * scale, 0.0f, 2.0f * scale, 1.0f * scale, 0.5f * scale);
    const auto result = RunProbe(*f.device, f.kernel, kOpMetric, {unit, tiny});

    for (int component = 0; component < 32; ++component) {
        EXPECT_TRUE(Close(result[kResultStride + component], result[component], 3.0e-5f, 3.0e-6f,
                          "scale-covariant metric component " + std::to_string(component)));
    }
    EXPECT_LT(result[32], 5.0e-6f);
    EXPECT_LT(result[kResultStride + 32], 5.0e-6f);

    const auto connection = RunProbe(*f.device, f.kernel, kOpChristoffel, {unit, tiny});
    for (int component = 0; component < 64; ++component) {
        EXPECT_TRUE(Close(scale * connection[kResultStride + component], connection[component],
                          2.0e-4f, 3.0e-5f,
                          "scale-covariant connection component " + std::to_string(component)));
    }
}

TEST(KernelParity, UnrepresentedKerrStageShrinksBeforeMetricEvaluation) {
    Fixture f = OpenProbe();
    if (!f.ready) GTEST_SKIP() << "no Vulkan device or kernels absent";

    // At a=0.5, (x,y,z)=(1,0,0) is on the differentiable positive-r
    // sheet. The RK4 midpoint for h=2 and u^x=-1 is exactly the non-unique
    // r=0 Kerr disk. The production adaptive authority must reject that
    // tableau and preserve the input state rather than sampling the metric.
    Sample crossing;
    crossing.metric_id = kKerrSchild;
    crossing.p1 = 1.0f;
    crossing.p2 = 0.5f;
    crossing.c0 = 1.0f;
    crossing.u1 = -1.0f;
    crossing.h = 2.0f;
    const auto result = RunProbe(*f.device, f.kernel, kOpAdaptiveEventDomain, {crossing});

    EXPECT_FLOAT_EQ(result[0], 1.0f) << "initial event was not represented";
    EXPECT_FLOAT_EQ(result[1], 0.0f) << "singular-sheet midpoint was admitted";
    EXPECT_FLOAT_EQ(result[2], 0.0f) << "unrepresented tableau was accepted";
    EXPECT_FLOAT_EQ(result[3], crossing.c0) << "rejected tableau mutated x";
    EXPECT_FLOAT_EQ(result[4], crossing.c1) << "rejected tableau mutated y";
    EXPECT_FLOAT_EQ(result[5], crossing.c2) << "rejected tableau mutated z";
    EXPECT_GT(result[6], 1.0f) << "rejection did not exceed the controller envelope";
    EXPECT_FLOAT_EQ(result[7], 0.5f) << "rejection did not request a smaller step";
}

TEST(KernelParity, NullProjectionPreservesConeBranchAndFailsClosed) {
    Fixture f = OpenProbe();
    if (!f.ready) GTEST_SKIP() << "no Vulkan device or kernels absent";

    Sample flat_past = KerrSchildAt(0.0f, 0.0f, 0.0f, 4.0f, 0.0f, 0.0f);
    flat_past.u0 = -0.8f;
    flat_past.u1 = 1.0f;
    Sample flat_future = flat_past;
    flat_future.u0 = 0.8f;

    // At Schwarzschild Kerr-Schild r=M, fixed outward spatial tangent has two
    // past-coordinate-time null roots, -1 and -3.  Nearest-root selection is
    // the branch witness; a sign-only selector cannot distinguish them.
    Sample inner_near = KerrSchildAt(1.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f);
    inner_near.u0 = -1.2f;
    inner_near.u1 = 1.0f;
    Sample inner_far = inner_near;
    inner_far.u0 = -2.8f;

    // A merely small nonzero g_00 remains quadratic.  Near the stationary
    // limit the second root is large but represented; collapsing this event to
    // the finite linear-limit root would switch branches discontinuously.
    Sample near_stationary = KerrSchildAt(1.0f, 0.0f, 0.0f, 1.99999f, 0.0f, 0.0f);
    near_stationary.u1 = 1.0f;
    const double h = 1.0 / static_cast<double>(near_stationary.c0);
    const double a = -1.0 + 2.0 * h;
    const double b = 4.0 * h;
    const double c = 1.0 + 2.0 * h;
    const double root_discriminant = b * b - 4.0 * a * c;
    const double first_root = (-b + std::sqrt(root_discriminant)) / (2.0 * a);
    const double second_root = (-b - std::sqrt(root_discriminant)) / (2.0 * a);
    const float large_root =
        static_cast<float>(std::abs(first_root) > std::abs(second_root) ? first_root : second_root);
    near_stationary.u0 = large_root + 1000.0f;

    // At r=2M, g_00 is exactly zero.  A radial tangent has a linear temporal
    // root; a purely transverse tangent has neither a quadratic nor a linear
    // solution and must be rejected.
    Sample linear = KerrSchildAt(1.0f, 0.0f, 0.0f, 2.0f, 0.0f, 0.0f);
    linear.u0 = 0.8f;
    linear.u1 = -1.0f;
    Sample absent = linear;
    absent.u1 = 0.0f;
    absent.u2 = 1.0f;

    const std::vector<Sample> samples = {flat_past,       flat_future, inner_near, inner_far,
                                         near_stationary, linear,      absent};
    const auto result = RunProbe(*f.device, f.kernel, kOpNullProjection, samples);
    const std::array<float, 6> expected_temporal = {-1.0f, 1.0f, -1.0f, -3.0f, large_root, 1.0f};

    for (std::size_t index = 0; index < expected_temporal.size(); ++index) {
        const std::size_t base = index * kResultStride;
        EXPECT_FLOAT_EQ(result[base], 1.0f) << "sample " << index;
        const float temporal_tolerance =
            std::max(2.0e-6f, 2.0e-5f * std::abs(expected_temporal[index]));
        EXPECT_NEAR(result[base + 1], expected_temporal[index], temporal_tolerance)
            << "temporal branch, sample " << index;
        EXPECT_FLOAT_EQ(result[base + 2], samples[index].u1);
        EXPECT_FLOAT_EQ(result[base + 3], samples[index].u2);
        EXPECT_FLOAT_EQ(result[base + 4], samples[index].u3);
        EXPECT_GT(result[base + 5], 1.0e-4f) << "sample was not detectably off-null";
        EXPECT_LT(result[base + 6], 2.0e-6f) << "projected null residual, sample " << index;
    }

    const std::size_t absent_base = expected_temporal.size() * kResultStride;
    EXPECT_FLOAT_EQ(result[absent_base], 0.0f);
    EXPECT_FLOAT_EQ(result[absent_base + 1], absent.u0);
    EXPECT_FLOAT_EQ(result[absent_base + 2], absent.u1);
    EXPECT_FLOAT_EQ(result[absent_base + 3], absent.u2);
    EXPECT_FLOAT_EQ(result[absent_base + 4], absent.u3);
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

TEST(KernelParity, IsotropicEllisMetricAndConnectionMatchCoreOnBothSheets) {
    Fixture f = OpenProbe();
    if (!f.ready) GTEST_SKIP() << "no Vulkan device or kernels absent";

    const auto sample_at = [](float b0, float x, float y, float z) {
        Sample sample;
        sample.metric_id = kMorrisThorne;
        sample.p1 = b0;
        sample.p2 = 0.0f;
        sample.p3 = 0.0f;  // exact Ellis selector
        sample.c0 = x;
        sample.c1 = y;
        sample.c2 = z;
        return sample;
    };
    const std::vector<Sample> samples = {
        sample_at(0.1f, 0.03f, 0.04f, 0.0f),  // public lower scale, exact throat
        sample_at(2.0f, 0.3f, 0.4f, 0.0f),    // second sheet, rho=b0/4
        sample_at(2.0f, 0.6f, 0.8f, 0.0f),    // exact throat, rho=b0/2
        sample_at(2.0f, 1.2f, 1.6f, 0.0f),    // exterior, rho=b0
        sample_at(1.0f, -2.0f, 1.0f, 3.0f),
        sample_at(1000.0f, 300.0f, 400.0f, 0.0f),  // public upper scale, exact throat
    };

    const auto metric_results = RunProbe(*f.device, f.kernel, kOpMetric, samples);
    const auto connection_results = RunProbe(*f.device, f.kernel, kOpChristoffel, samples);
    for (std::size_t sample_index = 0; sample_index < samples.size(); ++sample_index) {
        const auto& sample = samples[sample_index];
        sirius::core::MorrisThorneCartesian host(
            sirius::core::MorrisThorneParams::Ellis(sample.p1));
        sirius::core::Vec4 position;
        position(1) = sample.c0;
        position(2) = sample.c1;
        position(3) = sample.c2;
        sirius::core::Metric4d metric;
        sirius::core::Tensor<sirius::core::Dual<double>, 4, 4, 4> derivative;
        host.Evaluate(position, metric, derivative);
        sirius::core::Metric4d inverse;
        ASSERT_TRUE(host.InverseMetric(position, inverse));
        const auto connection = sirius::core::TensorOps::Christoffel(metric, derivative);
        const std::size_t base = sample_index * kResultStride;

        for (int row = 0; row < 4; ++row) {
            for (int column = 0; column < 4; ++column) {
                const int component = row * 4 + column;
                EXPECT_TRUE(Close(metric_results[base + component],
                                  static_cast<float>(metric(row, column).real), 8.0e-6f, 2.0e-6f,
                                  "isotropic Ellis metric"));
                EXPECT_TRUE(Close(metric_results[base + 16 + component],
                                  static_cast<float>(inverse(row, column).real), 8.0e-6f, 2.0e-6f,
                                  "isotropic Ellis inverse"));
                for (int lower = 0; lower < 4; ++lower) {
                    const int gamma_component = row * 16 + column * 4 + lower;
                    EXPECT_TRUE(Close(connection_results[base + gamma_component],
                                      static_cast<float>(connection.gamma(row, column, lower).real),
                                      2.0e-5f, 3.0e-6f, "isotropic Ellis connection"));
                }
            }
        }
        EXPECT_LT(metric_results[base + 32], 5.0e-6f);
    }
}

TEST(KernelParity, UnnormalisedOrNonEllisDeviceProfilesFailClosed) {
    Fixture f = OpenProbe();
    if (!f.ready) GTEST_SKIP() << "no Vulkan device or kernels absent";

    Sample represented;
    represented.metric_id = kMorrisThorne;
    represented.p1 = 1.0f;
    represented.p2 = 0.0f;
    represented.p3 = 0.0f;
    represented.c0 = 1.0f;

    Sample unnormalised = represented;
    unnormalised.p2 = 0.25f;
    Sample non_ellis = represented;
    non_ellis.p3 = 1.0f;

    const auto result =
        RunProbe(*f.device, f.kernel, kOpMetric, {represented, unnormalised, non_ellis});
    EXPECT_FLOAT_EQ(result[0], -1.0f);
    EXPECT_GT(result[5], 0.0f);
    for (std::size_t sample = 1; sample < 3; ++sample) {
        const std::size_t base = sample * kResultStride;
        for (std::size_t component = 0; component < 32; ++component) {
            EXPECT_FLOAT_EQ(result[base + component], 0.0f)
                << "sample=" << sample << " component=" << component;
        }
    }
}

TEST(KernelParity, SphericalCaptureEventFindsHiddenAndTangentContacts) {
    Fixture f = OpenProbe();
    if (!f.ready) GTEST_SKIP() << "no Vulkan device or kernels absent";

    const auto event_sample = [](float tangent_magnitude) {
        Sample sample;
        sample.p1 = 1.0f;                // capture radius
        sample.p2 = tangent_magnitude;   // end dx/dlambda
        sample.c0 = 1.1f;                // start x
        sample.u0 = 1.0f;                // start dt/dlambda
        sample.u1 = -tangent_magnitude;  // start dx/dlambda
        sample.h = 1.0f;                 // affine length and end t
        sample.aux0 = 1.1f;              // end x
        return sample;
    };
    const std::vector<Sample> samples = {
        event_sample(0.8f),  // enter and exit with both endpoints outside
        event_sample(0.4f),  // tangent contact at the midpoint
        event_sample(0.3f),  // entirely outside
        event_sample(std::numeric_limits<float>::quiet_NaN()),  // malformed segment
    };
    const auto result = RunProbe(*f.device, f.kernel, kOpSphericalCaptureEvent, samples);

    const float hidden_fraction = (0.8f - std::sqrt(0.32f)) / 1.6f;
    EXPECT_FLOAT_EQ(result[0], 1.0f);
    EXPECT_NEAR(result[1], hidden_fraction, 2.0e-5f);
    EXPECT_NEAR(result[2], 1.0f, 3.0e-6f);
    EXPECT_NEAR(result[5], -std::sqrt(0.32f), 3.0e-5f);
    EXPECT_NEAR(result[8], 1.0f, 3.0e-6f);

    const std::size_t tangent = kResultStride;
    EXPECT_FLOAT_EQ(result[tangent], 1.0f);
    EXPECT_NEAR(result[tangent + 1], 0.5f, 2.0e-5f);
    EXPECT_NEAR(result[tangent + 2], 1.0f, 3.0e-6f);
    EXPECT_NEAR(result[tangent + 5], 0.0f, 3.0e-5f);

    const std::size_t outside = 2 * kResultStride;
    EXPECT_FLOAT_EQ(result[outside], 0.0f);

    const std::size_t malformed = 3 * kResultStride;
    EXPECT_FLOAT_EQ(result[malformed], 0.0f);
}

TEST(KernelParity, EllisTwoSheetTraceCrossesThroatAndMapsTheOppositeSky) {
    Fixture f = OpenProbe();
    if (!f.ready) GTEST_SKIP() << "no Vulkan device or kernels absent";

    const auto central_ray = [](float b0) {
        Sample sample;
        sample.metric_id = kMorrisThorne;
        sample.p1 = b0;
        sample.p2 = 20.0f * b0;
        sample.c0 = 10.0f * b0;
        // A=1+b0^2/(4 rho^2)=1.0025 at rho=10*b0. Choosing
        // drho/dlambda=-1/A makes the past-ray frequency exactly one, so the
        // analytic proper-distance difference is also its affine length.
        sample.u1 = -1.0f / 1.0025f;
        sample.h = 0.08f;
        sample.aux0 = 20000.0f;
        sample.aux1 = 1.0e-5f * b0;
        sample.aux2 = 2.0f * b0;
        return sample;
    };
    // Cross-scale inversion identities are proved directly by the host tests;
    // one live device trajectory is sufficient to bind the shared production
    // function without tripling the software-Vulkan compute cost.
    const std::vector<Sample> samples = {central_ray(1.0f)};
    const auto result = RunProbe(*f.device, f.kernel, kOpEllisTwoSheetTrace, samples);

    for (std::size_t sample = 0; sample < samples.size(); ++sample) {
        const std::size_t base = sample * kResultStride;
        const float b0 = samples[sample].p1;
        const float opposite_radius = b0 / 80.0f;
        SCOPED_TRACE(b0);
        EXPECT_FLOAT_EQ(result[base], 2.0f);
        EXPECT_NEAR(result[base + 1], opposite_radius, 2.0e-4f * b0);
        const float expected_affine = 29.9625f * b0;
        EXPECT_NEAR(result[base + 2], expected_affine, 2.0e-2f * b0);
        EXPECT_GT(result[base + 3], 0.0f);
        EXPECT_LT(result[base + 3], samples[sample].aux0);
        EXPECT_GT(result[base + 4], 0.0f);
        EXPECT_NEAR(result[base + 5], 0.0f, 2.0e-5f * b0);
        EXPECT_NEAR(result[base + 6], 0.0f, 2.0e-5f * b0);
        EXPECT_GT(result[base + 7], 0.999f);
        EXPECT_NEAR(result[base + 8], 0.0f, 2.0e-5f);
        EXPECT_NEAR(result[base + 9], 0.0f, 2.0e-5f);
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

TEST(KernelParity, RepresentedSmallKerrSpinDoesNotAliasToSchwarzschildIsco) {
    Fixture f = OpenProbe();
    if (!f.ready) GTEST_SKIP() << "no Vulkan device or kernels absent";

    constexpr float spin = 9.0e-4f;
    Sample sample;
    sample.p1 = 1.0f;
    sample.p2 = spin;
    sample.c0 = 10.0f;
    sample.aux0 = static_cast<float>(sirius::core::AccretionDiskD::ComputeIsco(spin));
    sample.aux1 = 10000.0f;

    const auto result = RunProbe(*f.device, f.kernel, kOpDiskTemp, {sample});
    const float device_isco = result[2];
    EXPECT_LT(device_isco, 6.0f);
    EXPECT_TRUE(Close(device_isco, sample.aux0, 1.0e-4f, 1.0e-6f, "small-spin ISCO"));
}

TEST(KernelParity, ShakuraSunyaevZeroTorqueTemperatureMatchesCoreModel) {
    Fixture f = OpenProbe();
    if (!f.ready) GTEST_SKIP() << "no Vulkan device or kernels absent";

    Sample near_edge;
    near_edge.p1 = 1.0f;
    near_edge.c0 = 6.6f;
    near_edge.aux0 = 6.0f;
    near_edge.aux1 = 10000.0f;

    Sample reference = near_edge;
    reference.c0 = 9.0f;

    Sample far = near_edge;
    far.c0 = 36.0f;

    Sample edge = near_edge;
    edge.c0 = edge.aux0;

    const std::vector<Sample> samples = {near_edge, reference, far, edge};
    const auto results = RunProbe(*f.device, f.kernel, kOpDiskTemp, samples);
    for (std::size_t sample = 0; sample < samples.size(); ++sample) {
        const float expected = static_cast<float>(sirius::core::ShakuraSunyaevTemperature(
            samples[sample].aux1, samples[sample].c0, samples[sample].aux0));
        EXPECT_TRUE(Close(results[sample * kResultStride + 3], expected, 2.0e-6f, 2.0e-3f,
                          "Shakura-Sunyaev temperature"));
    }
    EXPECT_FLOAT_EQ(results[1 * kResultStride + 3], reference.aux1);
    EXPECT_GT(results[0 * kResultStride + 3], 0.0f);
    EXPECT_LT(results[2 * kResultStride + 3], reference.aux1);
    EXPECT_FLOAT_EQ(results[3 * kResultStride + 3], 0.0f);
}

TEST(KernelParity, TruncatedGaussianOpacityMatchesFiniteColumnCoreClosure) {
    Fixture f = OpenProbe();
    if (!f.ready) GTEST_SKIP() << "no Vulkan device or kernels absent";

    const auto sample = [](float radius, float height, float reference_radius, float scale_height,
                           float tau_ref) {
        Sample value;
        value.p1 = reference_radius;
        value.p2 = scale_height;
        value.p3 = tau_ref;
        value.c0 = radius;
        value.c1 = height;
        return value;
    };
    const std::vector<Sample> samples = {
        sample(6.0f, 0.0f, 6.0f, 0.6f, 2.0f),
        sample(12.0f, 1.2f, 6.0f, 1.2f, 2.0f),
        sample(12.0f, 3.01f * 1.2f, 6.0f, 1.2f, 2.0f),
    };
    const auto results = RunProbe(*f.device, f.kernel, kOpVolumeOpacity, samples);
    for (std::size_t index = 0; index < samples.size(); ++index) {
        const auto& value = samples[index];
        const float expected =
            static_cast<float>(sirius::core::volumetric_disk::TruncatedGaussianOpacityDensity(
                value.p3, value.c0, value.p1, value.p2, value.c1));
        EXPECT_TRUE(Close(results[index * kResultStride], expected, 3.0e-6f, 1.0e-7f,
                          "truncated Gaussian opacity"));
    }
    EXPECT_GT(results[0], 0.0f);
    EXPECT_GT(results[kResultStride], 0.0f);
    EXPECT_FLOAT_EQ(results[2 * kResultStride], 0.0f);
}

TEST(KernelParity, OpticallyThinGreyLayerAbsorptionMatchesHostAuthority) {
    Fixture f = OpenProbe();
    if (!f.ready) GTEST_SKIP() << "no Vulkan device or kernels absent";

    Sample zero, very_thin, thin, transition, thick;
    zero.c0 = 0.0f;
    very_thin.c0 = 1.0e-8f;
    thin.c0 = 1.0e-5f;
    transition.c0 = 1.0e-3f;
    thick.c0 = static_cast<float>(std::log(2.0));
    const std::vector<Sample> samples = {zero, very_thin, thin, transition, thick};
    const auto results = RunProbe(*f.device, f.kernel, kOpGreyLayerAbsorption, samples);

    for (std::size_t index = 0; index < samples.size(); ++index) {
        const auto expected = sirius::core::relativity::GreyLayerAbsorbedFraction(
            static_cast<double>(samples[index].c0));
        ASSERT_TRUE(expected.has_value());
        EXPECT_TRUE(Close(results[index * kResultStride], static_cast<float>(*expected), 3.0e-6f,
                          2.0e-12f, "grey layer absorbed fraction"));
    }
    EXPECT_GT(results[kResultStride], 0.0f)
        << "the fp32 transfer mirror discarded a represented optically thin layer";
    EXPECT_NEAR(results[4 * kResultStride], 0.5f, 1.0e-6f);
}

TEST(KernelParity, ArbitraryLatitudeKerrZamoMatchesHostAuthority) {
    Fixture f = OpenProbe();
    if (!f.ready) GTEST_SKIP() << "no Vulkan device or kernels absent";

    const auto sample = [](float mass, float spin, float radius, float cos_theta,
                           float observer_frequency, float energy, float angular_momentum) {
        Sample value;
        value.p1 = mass;
        value.p2 = spin;
        value.c0 = radius;
        value.c1 = cos_theta;
        value.c2 = observer_frequency;
        value.u0 = energy;
        value.u1 = angular_momentum;
        return value;
    };
    const std::vector<Sample> samples = {
        sample(1.0f, 0.7f, 6.0f, 0.35f, 0.93f, 1.0f, 1.4f),
        sample(2.0f, -1.2f, 12.0f, -0.25f, 1.1f, 1.2f, -2.0f),
        sample(1.0f, 0.998f, 1.7f, 0.2f, 1.0f, 1.0f, 0.2f),
        sample(1.0f, 0.7f, 6.0f, 1.0f, 0.93f, 1.0f, 1.4f),
    };
    const auto results = RunProbe(*f.device, f.kernel, kOpKerrZamoTransfer, samples);
    for (std::size_t index = 0; index < samples.size(); ++index) {
        const auto& value = samples[index];
        const auto expected = sirius::core::relativity::KerrZamoFrequencyTransfer(
            value.c2, value.u0, value.u1, value.p1, value.p2, value.c0, value.c1);
        const std::size_t base = index * kResultStride;
        EXPECT_EQ(results[base] != 0.0f, expected.has_value());
        if (!expected) continue;
        EXPECT_TRUE(Close(results[base + 1], static_cast<float>(expected->g), 4.0e-6f, 2.0e-7f,
                          "arbitrary-latitude ZAMO g"));
        EXPECT_TRUE(Close(results[base + 2], static_cast<float>(expected->frame_frequency), 4.0e-6f,
                          2.0e-7f, "arbitrary-latitude ZAMO frequency"));
        EXPECT_TRUE(Close(results[base + 3], static_cast<float>(expected->frame.angular_velocity),
                          4.0e-6f, 2.0e-7f, "arbitrary-latitude ZAMO omega"));
        EXPECT_TRUE(Close(results[base + 4], static_cast<float>(expected->frame.time_component),
                          4.0e-6f, 2.0e-7f, "arbitrary-latitude ZAMO u^t"));
    }
}

TEST(KernelParity, EquatorialKerrDiskTransferMatchesHostAuthority) {
    Fixture f = OpenProbe();
    if (!f.ready) GTEST_SKIP() << "no Vulkan device or kernels absent";

    const auto sample = [](float mass, float spin, float radius, float phi,
                           float observer_frequency, std::array<float, 4> past_velocity) {
        Sample value;
        value.p1 = mass;
        value.p2 = spin;
        value.c0 = radius;
        value.c1 = phi;
        value.c2 = observer_frequency;
        value.u0 = past_velocity[0];
        value.u1 = past_velocity[1];
        value.u2 = past_velocity[2];
        value.u3 = past_velocity[3];
        return value;
    };
    const std::vector<Sample> samples = {
        sample(1.0f, 0.7f, 6.0f, 0.4f, 0.93f, {-1.2f, 0.3f, -0.1f, 0.05f}),
        sample(2.0f, -1.2f, 12.0f, -0.8f, 1.1f, {-1.1f, -0.15f, 0.24f, -0.03f}),
    };
    const auto results = RunProbe(*f.device, f.kernel, kOpKerrDiskTransfer, samples);
    for (std::size_t index = 0; index < samples.size(); ++index) {
        const auto& value = samples[index];
        const double rho = std::hypot(static_cast<double>(value.c0), static_cast<double>(value.p2));
        sirius::core::Vec4 position;
        position(1) = rho * std::cos(value.c1);
        position(2) = rho * std::sin(value.c1);
        sirius::core::KerrSchildFamily metric(
            sirius::core::KerrSchildParams::Kerr(value.p1, value.p2));
        sirius::core::Metric4d local_metric;
        sirius::core::Tensor<sirius::core::Dual<double>, 4, 4, 4> derivatives;
        metric.Evaluate(position, local_metric, derivatives);
        sirius::core::Vec4 physical_photon;
        physical_photon(0) = -value.u0;
        physical_photon(1) = -value.u1;
        physical_photon(2) = -value.u2;
        physical_photon(3) = -value.u3;
        const auto covector = sirius::core::TensorOps::LowerIndex(physical_photon, local_metric);
        const double energy = -covector(0);
        const double angular_momentum = -position(2) * covector(1) + position(1) * covector(2);
        const auto expected = sirius::core::relativity::KerrDiskTransfer(
            value.c2, energy, angular_momentum, value.p1, value.p2, value.c0);
        ASSERT_TRUE(expected.has_value());
        const std::size_t base = index * kResultStride;
        EXPECT_TRUE(Close(results[base], static_cast<float>(expected->full_g), 5.0e-6f, 2.0e-7f,
                          "equatorial emitter g"));
        EXPECT_TRUE(Close(results[base + 1], static_cast<float>(expected->zamo_g), 5.0e-6f, 2.0e-7f,
                          "equatorial ZAMO g"));
    }
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

TEST(KernelParity, XyzD65ToLinearSrgbMatchesHostAuthority) {
    Fixture f = OpenProbe();
    if (!f.ready) GTEST_SKIP() << "no Vulkan device or kernels absent";

    Sample d65;
    d65.c0 = 0.9504559f;
    d65.c1 = 1.0f;
    d65.c2 = 1.0890578f;
    Sample x_basis;
    x_basis.c0 = 1.0f;
    Sample y_basis;
    y_basis.c1 = 1.0f;
    Sample z_basis;
    z_basis.c2 = 1.0f;
    Sample outside_gamut;
    outside_gamut.c0 = 0.25f;
    outside_gamut.c1 = 0.4f;
    outside_gamut.c2 = 0.1f;

    const std::vector<Sample> samples = {d65, x_basis, y_basis, z_basis, outside_gamut};
    const auto results = RunProbe(*f.device, f.kernel, kOpXyzD65ToLinearSrgb, samples);
    for (std::size_t index = 0; index < samples.size(); ++index) {
        const auto& sample_value = samples[index];
        const auto expected = sirius::core::colour::XyzD65ToLinearSrgb(
            sample_value.c0, sample_value.c1, sample_value.c2);
        const std::array reference = {expected.r, expected.g, expected.b};
        for (std::size_t channel = 0; channel < reference.size(); ++channel) {
            EXPECT_TRUE(Close(results[index * kResultStride + channel], reference[channel], 3.0e-6f,
                              2.0e-7f, "XYZ D65 to linear sRGB"));
        }
    }
}

TEST(KernelParity, Cie1931TwoDegreeFitMatchesHostAuthority) {
    Fixture f = OpenProbe();
    if (!f.ready) GTEST_SKIP() << "no Vulkan device or kernels absent";

    std::vector<Sample> samples;
    for (const float wavelength :
         {379.0f, 380.0f, 400.0f, 440.0f, 500.0f, 555.0f, 630.0f, 700.0f, 780.0f, 781.0f}) {
        Sample sample;
        sample.c0 = wavelength;
        samples.push_back(sample);
    }
    const auto results = RunProbe(*f.device, f.kernel, kOpCie1931TwoDegreeFit, samples);
    for (std::size_t index = 0; index < samples.size(); ++index) {
        const auto expected = sirius::core::colour::Cie1931TwoDegreeFit(samples[index].c0);
        const std::array reference = {expected.x_bar, expected.y_bar, expected.z_bar};
        for (std::size_t channel = 0; channel < reference.size(); ++channel) {
            EXPECT_TRUE(Close(results[index * kResultStride + channel],
                              static_cast<float>(reference[channel]), 3.0e-6f, 2.0e-7f,
                              "CIE 1931 2-degree observer fit"));
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

TEST(KernelParity, PrecisionProbeArtifactsCarryOnlyTheirDeclaredFloat64Capability) {
#ifndef SIRIUS_KERNEL_DIR
    GTEST_SKIP() << "kernels not compiled (slangc absent at configure time)";
#else
    const auto fp32 = LoadSpirv(std::string(SIRIUS_KERNEL_DIR) + "/parity_probe.spv");
    const auto compensated =
        LoadSpirv(std::string(SIRIUS_KERNEL_DIR) + "/parity_probe_fp32comp.spv");
    const auto fp64 = LoadSpirv(std::string(SIRIUS_KERNEL_DIR) + "/parity_probe_fp64.spv");
    ASSERT_FALSE(fp32.empty());
    ASSERT_FALSE(compensated.empty());
    ASSERT_FALSE(fp64.empty());
    EXPECT_FALSE(HasSpirvCapability(fp32, 10u));
    EXPECT_FALSE(HasSpirvCapability(compensated, 10u));
    EXPECT_TRUE(HasSpirvCapability(fp64, 10u));
#endif
}

TEST(KernelParity, PrecisionRungsConserveNearExtremalKerrWithoutImageComparison) {
    Fixture fp32 = OpenProbe();
    if (!fp32.ready) GTEST_SKIP() << "no Vulkan device or kernels absent";
    if (!fp32.device->Info().supports_fp64) {
        GTEST_SKIP() << "device lacks shaderFloat64";
    }
    Fixture compensated = OpenProbe("parity_probe_fp32comp.spv");
    Fixture fp64 = OpenProbe("parity_probe_fp64.spv");
    ASSERT_TRUE(compensated.ready) << "compensated precision probe unavailable";
    ASSERT_TRUE(fp64.ready) << "fp64 precision probe unavailable";

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
    near_extremal.h = 0.08f;
    near_extremal.aux0 = 3000.0f;
    near_extremal.aux1 = 0.02f;
    near_extremal.aux2 = 2.0f;

    const auto result32 =
        RunProbe(*fp32.device, fp32.kernel, kOpLiveCartConservation, {near_extremal});
    const auto result_compensated =
        RunProbe(*compensated.device, compensated.kernel, kOpLiveCartConservation, {near_extremal});
    const auto result64 =
        RunProbe(*fp64.device, fp64.kernel, kOpLiveCartConservation, {near_extremal});

    const auto validate_progress = [](const std::vector<float>& result, const char* rung) {
        for (int component = 0; component < 12; ++component) {
            EXPECT_TRUE(std::isfinite(result[component]))
                << rung << " result component " << component << " is non-finite";
        }
        EXPECT_GT(result[4], 100.0f) << rung << " made too little progress";
        EXPECT_FLOAT_EQ(result[9], 2.0f)
            << rung << " did not escape; last rejected error ratio=" << result[10]
            << ", null residual=" << result[11];
        EXPECT_GE(result[5], 200.0f) << rung << " did not reach the escape surface";
    };
    validate_progress(result32, "fp32");
    validate_progress(result_compensated, "fp32-comp");
    validate_progress(result64, "fp64");

    EXPECT_LT(result64[0], 1.0e-6f) << "fp64 energy drift";
    EXPECT_LT(result64[1], 1.0e-6f) << "fp64 angular-momentum drift";
    EXPECT_LT(result64[2], 1.0e-6f) << "fp64 Carter drift";
    EXPECT_LT(result64[3], 1.0e-8f) << "fp64 relative null residual";

    const float drift32 = result32[0] + result32[1] + result32[2] + result32[3];
    const float drift_compensated = result_compensated[0] + result_compensated[1] +
                                    result_compensated[2] + result_compensated[3];
    const float drift64 = result64[0] + result64[1] + result64[2] + result64[3];
    EXPECT_LE(drift_compensated, drift32 * 1.25f + 1.0e-7f)
        << "compensated accumulation regressed the invariant envelope";
    EXPECT_LE(drift64, drift32 + 1.0e-8f) << "fp64 failed to improve the invariant envelope";
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

TEST(KernelParity, CelestialTangentBasisIsSharedByBeamAndPointFilter) {
    Fixture f = OpenProbe();
    if (!f.ready) GTEST_SKIP() << "no Vulkan device or kernels absent";

    // Include directions for which the least-aligned axis is respectively y,
    // x, and z. The old point filter's z/near-z heuristic disagreed with the
    // beam Sachs basis on the first two fixtures.
    std::vector<Sample> samples(3);
    samples[0].c0 = 0.8f;
    samples[0].c1 = 0.1f;
    samples[0].c2 = 0.59f;
    samples[1].c0 = 0.1f;
    samples[1].c1 = 0.7f;
    samples[1].c2 = 0.7f;
    samples[2].c0 = 0.6f;
    samples[2].c1 = 0.7f;
    samples[2].c2 = 0.05f;

    const auto results = RunProbe(*f.device, f.kernel, kOpCelestialTangentBasis, samples);
    for (std::size_t sample = 0; sample < samples.size(); ++sample) {
        const std::array<float, 3> direction{samples[sample].c0, samples[sample].c1,
                                             samples[sample].c2};
        const auto expected = sirius::core::relativity::MakeCelestialTangentBasis(direction);
        ASSERT_TRUE(expected.has_value());
        const std::size_t offset = sample * kResultStride;
        for (std::size_t component = 0; component < 3; ++component) {
            EXPECT_NEAR(results[offset + component], expected->first[component], 2.0e-6f);
            EXPECT_NEAR(results[offset + 3 + component], expected->second[component], 2.0e-6f);
        }
    }
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
    s.aux0 = 1.0f;  // xi^x; h,aux0..2 carry the four deviation components.

    const std::vector<Sample> samples = {s};
    const auto r = RunProbe(*f.device, f.kernel, kOpDeviation, samples);
    for (int k = 0; k < 9; ++k) {
        EXPECT_TRUE(std::isfinite(r[k])) << "deviation component " << k << " non-finite";
    }
    EXPECT_GT(r[8], 1e-6f) << "full-Riemann tidal acceleration vanishes near the hole";
}

TEST(KernelParity, DeviceTidalContractionMatchesAnalyticSchwarzschildAtMatchedEvents) {
    Fixture f = OpenProbe();
    if (!f.ready) GTEST_SKIP() << "no Vulkan device or kernels absent";

    constexpr double mass = 1.0;
    constexpr double theta = 1.1;
    constexpr double phi = 0.37;
    constexpr double energy = 1.0;
    constexpr double polar_momentum = 1.3;
    constexpr double angular_momentum = 2.0;
    sirius::oracle::KerrMetricD oracle(mass, 0.0);

    std::vector<Sample> samples;
    std::vector<sirius::oracle::Vec4d> events;
    std::vector<sirius::oracle::Vec4d> tangents;
    std::vector<sirius::oracle::Vec4d> deviations;
    for (const double radius : {4.0, 6.0, 10.0, 20.0}) {
        const sirius::oracle::Vec4d event(0.0, radius, theta, phi);
        const double lapse = 1.0 - 2.0 * mass / radius;
        const double sin_theta = std::sin(theta);
        const double angular_norm =
            (polar_momentum * polar_momentum +
             angular_momentum * angular_momentum / (sin_theta * sin_theta)) /
            (radius * radius);
        const double radial = std::sqrt(energy * energy - lapse * angular_norm);
        const sirius::oracle::Vec4d tangent(
            energy / lapse, radial, polar_momentum / (radius * radius),
            angular_momentum / (radius * radius * sin_theta * sin_theta));
        const sirius::oracle::Vec4d deviation(0.2, 0.15, 0.07 / radius,
                                              -0.11 / (radius * sin_theta));
        const auto live_tangent = SchwarzschildBlToKerrSchildVector(tangent, event, mass);
        const auto live_deviation = SchwarzschildBlToKerrSchildVector(deviation, event, mass);

        Sample sample;
        sample.metric_id = kKerrSchild;
        sample.p1 = static_cast<float>(mass);
        sample.c0 = static_cast<float>(radius * sin_theta * std::cos(phi));
        sample.c1 = static_cast<float>(radius * sin_theta * std::sin(phi));
        sample.c2 = static_cast<float>(radius * std::cos(theta));
        sample.u0 = static_cast<float>(live_tangent[0]);
        sample.u1 = static_cast<float>(live_tangent[1]);
        sample.u2 = static_cast<float>(live_tangent[2]);
        sample.u3 = static_cast<float>(live_tangent[3]);
        sample.h = static_cast<float>(live_deviation[0]);
        sample.aux0 = static_cast<float>(live_deviation[1]);
        sample.aux1 = static_cast<float>(live_deviation[2]);
        sample.aux2 = static_cast<float>(live_deviation[3]);
        samples.push_back(sample);
        events.push_back(event);
        tangents.push_back(tangent);
        deviations.push_back(deviation);
    }

    const auto results = RunProbe(*f.device, f.kernel, kOpDeviation, samples);
    for (std::size_t sample = 0; sample < samples.size(); ++sample) {
        const std::size_t offset = sample * kResultStride;
        const std::array<float, 4> device_acceleration = {results[offset], results[offset + 1],
                                                          results[offset + 2], results[offset + 3]};
        const auto actual =
            KerrSchildToSchwarzschildBlVector(device_acceleration, events[sample], mass);
        const auto expected =
            AnalyticTidalAcceleration(oracle, events[sample], tangents[sample], deviations[sample]);

        double scale = 0.0;
        double maximum_error = 0.0;
        for (int component = 0; component < 4; ++component) {
            ASSERT_TRUE(std::isfinite(actual[component]))
                << "non-finite device contraction component " << component
                << " at r=" << events[sample].r;
            ASSERT_TRUE(std::isfinite(expected[component]))
                << "non-finite analytic contraction component " << component
                << " at r=" << events[sample].r;
            scale = std::max(scale, std::abs(expected[component]));
            maximum_error =
                std::max(maximum_error, std::abs(actual[component] - expected[component]));
        }
        // The production path is fp32 and finite-differences an analytic
        // Christoffel tensor. Its precision-matched five-point stencil stays
        // below 5e-4 here, versus a measured 1.37e-2 for the former stencil.
        EXPECT_LT(maximum_error / std::max(scale, 1.0e-30), 5.0e-4)
            << "device matched-event tidal contraction at r=" << events[sample].r;
    }
}

TEST(KernelParity, DeviceTidalContractionMatchesAnalyticKerrAtMatchedEvents) {
    Fixture f = OpenProbe();
    if (!f.ready) GTEST_SKIP() << "no Vulkan device or kernels absent";

    constexpr double mass = 1.0;
    constexpr double phi = 0.37;

    std::vector<Sample> samples;
    std::vector<sirius::core::coordinates::Vec4Cart> positions;
    std::vector<double> spins;
    std::vector<sirius::oracle::Vec4d> events;
    std::vector<sirius::oracle::Vec4d> expected_accelerations;
    // Exercise both orientations of the rotating geometry, off the equator,
    // with two independent tangent/deviation pairs. The production coordinate
    // authority maps each Cartesian vector to the matched Boyer-Lindquist
    // event; only the independent analytic oracle supplies the expected
    // Riemann contraction. This is a scalar/tiny-compute probe, not a render.
    for (const float device_spin : {-0.9f, 0.9f}) {
        const double spin = static_cast<double>(device_spin);
        sirius::oracle::KerrMetricD oracle(mass, spin);
        for (const float radius : {3.0f, 6.0f, 20.0f}) {
            for (const float theta : {0.6f, 1.1f}) {
                const sirius::core::coordinates::Vec4Bl requested_event{0.0, radius, theta, phi};
                const auto requested_position =
                    sirius::core::coordinates::BlToKerrSchildCart(requested_event, spin);
                // Round the event to the exact fp32 storage-buffer bytes, then
                // recover the matched oracle event from those bytes. The
                // oracle may not silently receive a nearby double-only event.
                const sirius::core::coordinates::Vec4Cart position{
                    0.0, static_cast<double>(static_cast<float>(requested_position.x)),
                    static_cast<double>(static_cast<float>(requested_position.y)),
                    static_cast<double>(static_cast<float>(requested_position.z))};
                const auto matched_event =
                    sirius::core::coordinates::KerrSchildCartToBl(position, spin);
                const sirius::oracle::Vec4d event(matched_event.t, matched_event.r,
                                                  matched_event.theta, matched_event.phi);
                const std::array<float, 4> tangent_components =
                    theta < 1.0f ? std::array{1.0f, -0.71f, 0.23f, -0.17f}
                                 : std::array{0.82f, 0.31f, -0.57f, 0.26f};
                const std::array<float, 4> deviation_components =
                    theta < 1.0f ? std::array{0.13f, 0.19f, -0.11f, 0.07f}
                                 : std::array{-0.09f, 0.07f, 0.17f, -0.13f};
                const sirius::core::coordinates::Vec4Cart cart_tangent = {
                    tangent_components[0], tangent_components[1], tangent_components[2],
                    tangent_components[3]};
                const sirius::core::coordinates::Vec4Cart cart_deviation = {
                    deviation_components[0], deviation_components[1], deviation_components[2],
                    deviation_components[3]};
                const auto tangent_bl =
                    sirius::core::coordinates::TransformVectorKerrSchildCartToBl(
                        cart_tangent, position, mass, spin);
                const auto deviation_bl =
                    sirius::core::coordinates::TransformVectorKerrSchildCartToBl(
                        cart_deviation, position, mass, spin);

                Sample sample;
                sample.metric_id = kKerrSchild;
                sample.p1 = static_cast<float>(mass);
                sample.p2 = device_spin;
                sample.c0 = static_cast<float>(position.x);
                sample.c1 = static_cast<float>(position.y);
                sample.c2 = static_cast<float>(position.z);
                sample.u0 = tangent_components[0];
                sample.u1 = tangent_components[1];
                sample.u2 = tangent_components[2];
                sample.u3 = tangent_components[3];
                sample.h = deviation_components[0];
                sample.aux0 = deviation_components[1];
                sample.aux1 = deviation_components[2];
                sample.aux2 = deviation_components[3];
                samples.push_back(sample);
                positions.push_back(position);
                spins.push_back(spin);
                events.push_back(event);
                const sirius::oracle::Vec4d tangent(tangent_bl.t, tangent_bl.r, tangent_bl.theta,
                                                    tangent_bl.phi);
                const sirius::oracle::Vec4d deviation(deviation_bl.t, deviation_bl.r,
                                                      deviation_bl.theta, deviation_bl.phi);
                expected_accelerations.push_back(
                    AnalyticTidalAcceleration(oracle, event, tangent, deviation));
            }
        }
    }

    const auto results = RunProbe(*f.device, f.kernel, kOpDeviation, samples);
    for (std::size_t sample = 0; sample < samples.size(); ++sample) {
        const std::size_t offset = sample * kResultStride;
        const sirius::core::coordinates::Vec4Cart device_acceleration{
            results[offset], results[offset + 1], results[offset + 2], results[offset + 3]};
        const auto actual_bl = sirius::core::coordinates::TransformVectorKerrSchildCartToBl(
            device_acceleration, positions[sample], mass, spins[sample]);
        const sirius::oracle::Vec4d actual(actual_bl.t, actual_bl.r, actual_bl.theta,
                                           actual_bl.phi);
        const auto& expected = expected_accelerations[sample];

        double scale = 0.0;
        double maximum_error = 0.0;
        for (int component = 0; component < 4; ++component) {
            ASSERT_TRUE(std::isfinite(actual[component]));
            ASSERT_TRUE(std::isfinite(expected[component]));
            scale = std::max(scale, std::abs(expected[component]));
            maximum_error =
                std::max(maximum_error, std::abs(actual[component] - expected[component]));
        }
        const double relative_error = maximum_error / std::max(scale, 1.0e-30);
        // The worst exploratory fp32 result was 7.02e-5. Keep a governed
        // 2e-4 ceiling for cross-driver arithmetic variance while remaining a
        // quantitative rotating-Kerr tensor gate.
        EXPECT_LT(relative_error, 2.0e-4)
            << "device Kerr matched-event contraction at a=" << spins[sample]
            << " r=" << events[sample].r << " theta=" << events[sample].theta
            << " relative_error=" << relative_error;
    }
}

TEST(KernelParity, DeviceRadialPointSourceCongruenceMatchesClosedForm) {
    Fixture f = OpenProbe();
    if (!f.ready) GTEST_SKIP() << "no Vulkan device or kernels absent";

    constexpr float mass = 1.0f;
    constexpr float seed = 1.0e-3f;
    std::vector<Sample> samples;
    for (const auto [start_radius, end_radius, max_step] :
         {std::array{12.0f, 10.0f, 0.2f}, std::array{20.0f, 16.0f, 0.2f},
          std::array{8.0f, 6.0f, 0.1f}, std::array{6.0f, 4.0f, 0.1f}}) {
        Sample sample;
        sample.p1 = mass;
        sample.c0 = start_radius;
        sample.c1 = end_radius;
        sample.c2 = max_step;
        sample.h = seed;
        samples.push_back(sample);
    }

    const auto results = RunProbe(*f.device, f.kernel, kOpJacobiRadialCongruence, samples);
    for (std::size_t sample = 0; sample < samples.size(); ++sample) {
        const std::size_t offset = sample * kResultStride;
        const double affine_length = samples[sample].c0 - samples[sample].c1;
        const double expected_axis = seed * affine_length;
        SCOPED_TRACE(::testing::Message() << "r=" << samples[sample].c0 << "->"
                                          << samples[sample].c1 << " h=" << samples[sample].c2);
        for (std::size_t component = 0; component < 10; ++component) {
            ASSERT_TRUE(std::isfinite(results[offset + component]))
                << "non-finite device congruence component " << component;
        }
        EXPECT_NEAR(results[offset], expected_axis, 5.0e-6 * expected_axis);
        EXPECT_NEAR(results[offset + 1], expected_axis, 5.0e-6 * expected_axis);
        for (std::size_t component = 2; component < 6; ++component) {
            EXPECT_NEAR(results[offset + component], 0.0f, 1.0e-6f * expected_axis);
        }
        EXPECT_NEAR(results[offset + 6], seed, 1.0e-6f * seed);
        EXPECT_NEAR(results[offset + 7], seed, 1.0e-6f * seed);
        EXPECT_NEAR(results[offset + 8], affine_length, 2.0e-5 * affine_length);
        EXPECT_GT(results[offset + 9], 0.0f);
        EXPECT_LT(results[offset + 9], 128.0f);
    }
}

}  // namespace
