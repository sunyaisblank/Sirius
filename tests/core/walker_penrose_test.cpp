// Live-path polarisation-transport suite (sirius::core): parallel transport of a
// polarisation vector along the Cartesian Kerr-Schild geodesics the CPU tracer
// integrates, checked against the double-precision Boyer-Lindquist oracle
// (oracle/polarisation_transport.h) through the coordinate-invariant
// Walker-Penrose constant. Validates specification E2 on the render-native path.
//
// The Walker-Penrose constant is a scalar, so it is the natural cross-chart
// comparison: the same physical ray, set up in Kerr-Schild Cartesian for the
// live path and transformed to Boyer-Lindquist for the oracle, must yield the
// same constant, and each path must conserve it. The live path meets the
// live-path conservation bar (constants::geodesic::kConservationTol = 1e-4); the
// oracle meets 1e-10. Coverage:
//  - the live path conserves the Walker-Penrose constant and preserves f.k = 0
//    and f.f = 1;
//  - live and oracle agree on the constant, which also pins the Kerr-Schild to
//    Boyer-Lindquist transform (the transported tangent stays null in the
//    oracle chart);
//  - the transport-to-EVPA map rotates Stokes Q/U through the expected angle.

#include "sirius/core/polarisation/walker_penrose.h"

#include "sirius/core/constants.h"
#include "sirius/core/coordinates.h"
#include "sirius/core/metrics/kerr_schild_family.h"
#include "sirius/oracle/polarisation_transport.h"

#include <gtest/gtest.h>

#include <cmath>
#include <complex>
#include <string>

namespace sirius::test {
using namespace sirius::core;
namespace kc = sirius::core::constants;

namespace {

// A physical ray specified once in Kerr-Schild Cartesian for the live path and
// transformed to Boyer-Lindquist for the oracle.
struct BoyerLindquistVector {
    double r, theta;
    double kt, kr, ktheta, kphi;
};

// Kerr-Schild Cartesian contravariant components -> Boyer-Lindquist contravariant
// components. Spatial part is the Kerr spheroidal transform; the t and phi
// components carry the Kerr-Schild twist dt_BL = dt_KS - (2Mr/Delta) dr,
// dphi_BL = dphi_KS - (a/Delta) dr. The twist sign is fixed by requiring the
// transported tangent to stay null in the oracle chart (asserted below).
// Reference: Visser, "The Kerr spacetime" (arXiv:0706.0622), Kerr-Schild forms.
BoyerLindquistVector KerrSchildToBoyerLindquist(double M, double a, const Vec4& X, const Vec4& V) {
    const coordinates::Vec4Cart position{X(0), X(1), X(2), X(3)};
    const coordinates::Vec4Cart vector{V(0), V(1), V(2), V(3)};
    const coordinates::Vec4Bl position_bl = coordinates::KerrSchildCartToBl(position, a);
    const coordinates::Vec4Bl vector_bl =
        coordinates::TransformVectorKerrSchildCartToBl(vector, position, M, a);
    return {position_bl.r, position_bl.theta, vector_bl.t,
            vector_bl.r,   vector_bl.theta,   vector_bl.phi};
}

std::complex<double> WalkerPenroseFromBoyerLindquist(const BoyerLindquistVector& K,
                                                     const BoyerLindquistVector& F, double a) {
    const double r = K.r, theta = K.theta;
    const double sinth = std::sin(theta), costh = std::cos(theta), sin2th = sinth * sinth;
    const double A = (K.kt * F.kr - K.kr * F.kt) + a * sin2th * (K.kr * F.kphi - K.kphi * F.kr);
    const double B = ((r * r + a * a) * (K.kphi * F.ktheta - K.ktheta * F.kphi) -
                      a * (K.kt * F.ktheta - K.ktheta * F.kt)) *
                     sinth;
    return std::complex<double>(A, -B) * std::complex<double>(r, -a * costh);
}

double CartesianRadius(const Vec4& x) { return std::sqrt(x(1) * x(1) + x(2) * x(2) + x(3) * x(3)); }

}  // namespace

//==============================================================================
// The live path conserves the Walker-Penrose constant and orthonormality
//==============================================================================

// Transporting a polarisation vector along a Cartesian Kerr-Schild geodesic
// conserves the Walker-Penrose constant (computed by transforming each sample to
// Boyer-Lindquist) and preserves f.k = 0 and f.f = 1. Bound
// constants::geodesic::kConservationTol = 1e-4, the Mandatory live-path
// conservation tolerance; the measured drift over a near-photon-sphere pass is
// ~1e-10, so the bar holds with wide margin.
TEST(WalkerPenroseLivePath, ConservesConstantAndOrthonormality) {
    const double M = 1.0, a = 0.9;
    KerrSchildFamily metric(KerrSchildParams::Kerr(M, a));

    PolarisedRay ray{};
    ray.position(1) = 12.0;
    ray.position(2) = 0.0;
    ray.position(3) = 1.5;
    ray.velocity(1) = -1.0;
    ray.velocity(2) = 0.35;
    ray.velocity(3) = -0.05;

    Metric4d g;
    Tensor<Dual<double>, 4, 4, 4> dg;
    metric.Evaluate(ray.position, g, dg);
    ray.velocity = TensorOps::NormalizeNull(ray.velocity, g);
    Vec4 trial;
    trial(3) = 1.0;
    ray.polarisation = MakeOrthonormalPolarisation(metric, ray.position, ray.velocity, trial);

    metric.Evaluate(ray.position, g, dg);
    const double fk0 = TensorOps::InnerProduct(ray.polarisation, ray.velocity, g);
    const double ff0 = TensorOps::InnerProduct(ray.polarisation, ray.polarisation, g);
    ASSERT_NEAR(fk0, 0.0, 1e-12);
    ASSERT_NEAR(ff0, 1.0, 1e-12);

    const BoyerLindquistVector K0 = KerrSchildToBoyerLindquist(M, a, ray.position, ray.velocity);
    const BoyerLindquistVector F0 =
        KerrSchildToBoyerLindquist(M, a, ray.position, ray.polarisation);
    const std::complex<double> kappa0 = WalkerPenroseFromBoyerLindquist(K0, F0, a);
    ASSERT_GT(std::abs(kappa0), 1e-6);

    double max_fk = 0, max_ff = 0, max_kappa = 0;
    bool escaped = false;
    for (int i = 0; i < 200000; ++i) {
        const double r_kerr =
            metric.ComputeKerrRadius(ray.position(1), ray.position(2), ray.position(3));
        if (r_kerr <= metric.OuterHorizonRadius() * 1.05) break;
        if (CartesianRadius(ray.position) > 55.0) {
            escaped = true;
            break;
        }
        ParallelTransportStep(metric, ray, 0.01);
        metric.Evaluate(ray.position, g, dg);
        max_fk =
            std::max(max_fk, std::abs(TensorOps::InnerProduct(ray.polarisation, ray.velocity, g)));
        max_ff = std::max(
            max_ff, std::abs(TensorOps::InnerProduct(ray.polarisation, ray.polarisation, g) - 1.0));
        const BoyerLindquistVector K = KerrSchildToBoyerLindquist(M, a, ray.position, ray.velocity);
        const BoyerLindquistVector F =
            KerrSchildToBoyerLindquist(M, a, ray.position, ray.polarisation);
        const std::complex<double> kappa = WalkerPenroseFromBoyerLindquist(K, F, a);
        max_kappa = std::max(max_kappa, std::abs(kappa - kappa0) / std::abs(kappa0));
    }
    ASSERT_TRUE(escaped);
    EXPECT_LT(max_fk, kc::geodesic::kConservationTol);
    EXPECT_LT(max_ff, kc::geodesic::kConservationTol);
    EXPECT_LT(max_kappa, kc::geodesic::kConservationTol);
}

//==============================================================================
// Live and oracle paths agree on the Walker-Penrose constant
//==============================================================================

// The same physical ray, set up in Kerr-Schild Cartesian and transformed to
// Boyer-Lindquist, produces the same Walker-Penrose constant, and both paths
// conserve it. The initial invariants (null tangent, f.k = 0, f.f = 1) carrying
// across the transform pins the Kerr-Schild to Boyer-Lindquist map, including the
// twist. Agreement bound is the live-path tolerance kConservationTol = 1e-4; the
// oracle holds its own to 1e-10.
TEST(WalkerPenroseLivePath, AgreesWithOracleAcrossCharts) {
    const double M = 1.0;
    for (const double a : {0.0, 0.9}) {
        SCOPED_TRACE("spin=" + std::to_string(a));
        KerrSchildFamily live(KerrSchildParams::Kerr(M, a));
        sirius::oracle::KerrMetricD oracle(M, a);

        // Physical ray in Kerr-Schild Cartesian.
        PolarisedRay ray{};
        ray.position(1) = 14.0;
        ray.position(2) = 0.0;
        ray.position(3) = 1.0;
        ray.velocity(1) = -1.0;
        // Schwarzschild needs the larger impact parameter to stay outside its
        // b_crit = 3 sqrt(3) M capture cone; the prograde Kerr witness remains
        // the stronger-field pass used by the original gate.
        ray.velocity(2) = a == 0.0 ? 0.45 : 0.30;
        ray.velocity(3) = -0.02;

        Metric4d g;
        Tensor<Dual<double>, 4, 4, 4> dg;
        live.Evaluate(ray.position, g, dg);
        ray.velocity = TensorOps::NormalizeNull(ray.velocity, g);
        Vec4 trial;
        trial(3) = 1.0;
        ray.polarisation = MakeOrthonormalPolarisation(live, ray.position, ray.velocity, trial);

        // Transform the shared initial data to Boyer-Lindquist for the oracle.
        const BoyerLindquistVector K0 =
            KerrSchildToBoyerLindquist(M, a, ray.position, ray.velocity);
        const BoyerLindquistVector F0 =
            KerrSchildToBoyerLindquist(M, a, ray.position, ray.polarisation);

        sirius::oracle::Vec4d x_bl(0.0, K0.r, K0.theta, 0.0);
        sirius::oracle::Vec4d k_bl(K0.kt, K0.kr, K0.ktheta, K0.kphi);
        sirius::oracle::Vec4d f_bl(F0.kt, F0.kr, F0.ktheta, F0.kphi);

        // The transform preserves the invariants: the tangent stays null and the
        // polarisation stays unit and orthogonal in the oracle chart. These pin the
        // twist (a wrong sign breaks the null condition badly, ~0.3 for Kerr).
        double gg[4][4], gg_inv[4][4];
        oracle.Evaluate(x_bl, gg, gg_inv);
        EXPECT_NEAR(sirius::oracle::InnerProductD(gg, k_bl, k_bl), 0.0, 1e-9);
        EXPECT_NEAR(sirius::oracle::InnerProductD(gg, f_bl, f_bl), 1.0, 1e-9);
        EXPECT_NEAR(sirius::oracle::InnerProductD(gg, f_bl, k_bl), 0.0, 1e-9);

        // Both charts read the same Walker-Penrose constant at the start.
        const std::complex<double> kappa_live0 = WalkerPenroseFromBoyerLindquist(K0, F0, a);
        sirius::oracle::PolarisedStateD oracle_state(x_bl, k_bl, f_bl);
        const std::complex<double> kappa_oracle0 =
            sirius::oracle::WalkerPenroseConstant(oracle_state, oracle);
        ASSERT_GT(std::abs(kappa_live0), 1e-6);
        EXPECT_LT(std::abs(kappa_live0 - kappa_oracle0) / std::abs(kappa_live0), 1e-9);

        // Advance the oracle to escape; it conserves the constant to oracle tier.
        sirius::oracle::PolarisedGeodesicIntegratorD integrator(&oracle);
        const sirius::oracle::PolarisedGeodesicIntegratorD::Result oracle_result =
            integrator.Integrate(oracle_state);
        ASSERT_TRUE(oracle_result.escaped);
        const std::complex<double> kappa_oracle_final =
            sirius::oracle::WalkerPenroseConstant(oracle_result.state, oracle);
        EXPECT_LT(std::abs(kappa_oracle_final - kappa_oracle0) / std::abs(kappa_oracle0),
                  sirius::oracle::kWalkerPenroseConservationTol);

        // Advance the live path to escape; it conserves the same constant to the
        // live-path tolerance, so the two charts stay in agreement throughout.
        std::complex<double> kappa_live_final = kappa_live0;
        bool escaped = false;
        for (int i = 0; i < 200000; ++i) {
            const double r_kerr =
                live.ComputeKerrRadius(ray.position(1), ray.position(2), ray.position(3));
            if (r_kerr <= live.OuterHorizonRadius() * 1.05) break;
            if (CartesianRadius(ray.position) > 55.0) {
                escaped = true;
                break;
            }
            ParallelTransportStep(live, ray, 0.01);
            const BoyerLindquistVector K =
                KerrSchildToBoyerLindquist(M, a, ray.position, ray.velocity);
            const BoyerLindquistVector F =
                KerrSchildToBoyerLindquist(M, a, ray.position, ray.polarisation);
            kappa_live_final = WalkerPenroseFromBoyerLindquist(K, F, a);
        }
        ASSERT_TRUE(escaped);
        EXPECT_LT(std::abs(kappa_live_final - kappa_live0) / std::abs(kappa_live0),
                  kc::geodesic::kConservationTol);

        // Live and oracle final readings of the invariant agree.
        EXPECT_LT(std::abs(kappa_live_final - kappa_oracle_final) / std::abs(kappa_live0),
                  kc::geodesic::kConservationTol);
    }
}

//==============================================================================
// Transport -> EVPA -> Stokes Q/U rotation
//==============================================================================

// A polarisation vector at angle alpha to the reference sky axis carries EVPA
// alpha, and rotating an emitted Stokes vector by it turns (Q, U) through
// 2 alpha: horizontally polarised light (Q = I) becomes (I cos 2alpha,
// I sin 2alpha). This exercises the map from the transported vector to the
// Stokes rotation, reusing the stokes.h Mueller algebra. Float tolerance 1e-6
// because the Stokes types are single precision.
TEST(WalkerPenroseLivePath, TransportedVectorRotatesStokesByEvpa) {
    // Flat Minkowski frame; e_ref = x-hat, e_perp = y-hat span the sky plane.
    Metric4d eta;
    eta(0, 0) = Dual<double>(-1.0);
    eta(1, 1) = Dual<double>(1.0);
    eta(2, 2) = Dual<double>(1.0);
    eta(3, 3) = Dual<double>(1.0);
    Vec4 e_ref;
    e_ref(1) = 1.0;
    Vec4 e_perp;
    e_perp(2) = 1.0;

    const double alpha = 0.4;  // radians
    Vec4 f;
    f(1) = std::cos(alpha);
    f(2) = std::sin(alpha);

    EXPECT_NEAR(ElectricVectorPositionAngle(f, e_ref, e_perp, eta), alpha, 1e-12);

    const StokesVector emitted = StokesVector::Horizontal(1.0f);  // (I, Q, U, V) = (1, 1, 0, 0)
    const StokesVector observed = RotateStokesToObserverFrame(emitted, f, e_ref, e_perp, eta);
    EXPECT_NEAR(observed.I, 1.0f, 1e-6f);
    EXPECT_NEAR(observed.Q, std::cos(2.0f * static_cast<float>(alpha)), 1e-6f);
    EXPECT_NEAR(observed.U, std::sin(2.0f * static_cast<float>(alpha)), 1e-6f);
    EXPECT_NEAR(observed.V, 0.0f, 1e-6f);
}

}  // namespace sirius::test
