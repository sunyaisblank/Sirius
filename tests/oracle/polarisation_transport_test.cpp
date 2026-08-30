// Oracle polarisation-transport suite (sirius::oracle): parallel transport of a
// polarisation vector along Kerr null geodesics and the conserved Walker-Penrose
// complex constant. Validates specification E2 on the double-precision reference
// path. Tolerances are drawn from the constants authority (core/constants.h)
// where one exists and from the transport header's named bar otherwise; each
// carries a comment naming why it is achievable.
//
// Coverage:
//  (a) f stays orthogonal to k and unit-normalised, and the geodesic conserves
//      E, L_z, the Carter constant Q and the null condition, under transport;
//  (b) the Walker-Penrose constant is conserved to 1e-10 relative over full
//      Kerr (a = 0.9) geodesics at several impact parameters, including a
//      near-photon-sphere pass;
//  (c) Schwarzschild equatorial limit: a polarisation vector perpendicular to
//      the orbital plane transports rigidly (f^theta * r constant, no in-plane
//      leakage), the analytic expectation;
//  (d) flat-space limit: the polarisation vector is constant in Cartesian
//      components.

#include "sirius/oracle/polarisation_transport.h"

#include "sirius/core/constants.h"

#include <gtest/gtest.h>

#include <cmath>
#include <complex>
#include <limits>
#include <numbers>
#include <optional>
#include <vector>

namespace sirius::test {
using namespace sirius::oracle;
namespace k = sirius::core::constants;

namespace {

TEST(WalkerPenrose, MalformedIntegratorStepDomainFailsClosed) {
    KerrMetricD metric(1.0, 0.0);
    PolarisedGeodesicIntegratorD::Config config;
    config.base_step = std::numeric_limits<double>::quiet_NaN();
    EXPECT_FALSE(PolarisedGeodesicIntegratorD::IsRepresentedConfig(config));
    EXPECT_DEATH((void)PolarisedGeodesicIntegratorD(&metric, config),
                 "precondition.*enforced, terminating");
}

TEST(WalkerPenrose, BoyerLindquistInitialDataAndAxisExitDeclineWithoutSubstitution) {
    KerrMetricD metric(1.0, 0.5);

    // At the equatorial stationary-limit surface g_tt=0, so the null equation
    // is linear in k^t. The constructor must solve that represented case rather
    // than divide a quadratic formula by zero.
    const Vec4d stationary_limit(0.0, 2.0, std::numbers::pi / 2.0, 0.0);
    ASSERT_TRUE(metric.IsValid(stationary_limit));
    const std::optional<Vec4d> linear_tangent =
        TryMakeNullTangent(metric, stationary_limit, 0.0, 0.0, 0.05);
    ASSERT_TRUE(linear_tangent);
    EXPECT_GT(linear_tangent->t, 0.0);
    double g[4][4], g_inv[4][4];
    metric.Evaluate(stationary_limit, g, g_inv);
    EXPECT_NEAR(InnerProductD(g, *linear_tangent, *linear_tangent), 0.0, 1.0e-13);

    // A direction with neither a quadratic nor linear temporal solution, and a
    // non-finite direction, are absent requests rather than clamped roots.
    EXPECT_FALSE(TryMakeNullTangent(metric, stationary_limit, 1.0, 0.0, 0.0));
    EXPECT_FALSE(TryMakeNullTangent(metric, stationary_limit,
                                    std::numeric_limits<double>::infinity(), 0.0, 0.05));

    // Inside the ergoregion the same fixed spatial coordinate direction can
    // produce two distinct positive k^t roots. Both are future-directed, so
    // choosing either would be an unrequested branch selection.
    const Vec4d ambiguous_event(0.0, 1.9, std::numbers::pi / 2.0, 0.0);
    ASSERT_TRUE(metric.IsValid(ambiguous_event));
    metric.Evaluate(ambiguous_event, g, g_inv);
    constexpr double ambiguous_k_phi = 0.05;
    const double quadratic_a = g[0][0];
    const double quadratic_b = 2.0 * g[0][3] * ambiguous_k_phi;
    const double quadratic_c = g[3][3] * ambiguous_k_phi * ambiguous_k_phi;
    const double discriminant = quadratic_b * quadratic_b - 4.0 * quadratic_a * quadratic_c;
    ASSERT_GT(discriminant, 0.0);
    EXPECT_GT((-quadratic_b + std::sqrt(discriminant)) / (2.0 * quadratic_a), 0.0);
    EXPECT_GT((-quadratic_b - std::sqrt(discriminant)) / (2.0 * quadratic_a), 0.0);
    EXPECT_FALSE(TryMakeNullTangent(metric, ambiguous_event, 0.0, 0.0, ambiguous_k_phi));

    const Vec4d near_axis(0.0, 10.0, 2.0 * kBoyerLindquistPoleMargin, 0.0);
    const std::optional<Vec4d> tangent = TryMakeNullTangent(metric, near_axis, 0.0, -0.01, 0.0);
    ASSERT_TRUE(tangent);
    const std::optional<Vec4d> polarisation =
        TryMakeOrthonormalPolarisation(metric, near_axis, *tangent, Vec4d(0.0, 1.0, 0.0, 0.0));
    ASSERT_TRUE(polarisation);
    EXPECT_FALSE(TryMakeOrthonormalPolarisation(metric, near_axis, Vec4d(1.0, 0.0, 0.0, 0.0),
                                                Vec4d(0.0, 1.0, 0.0, 0.0)));
    EXPECT_FALSE(TryMakeOrthonormalPolarisation(metric, near_axis, *tangent, *tangent));

    const PolarisedStateD initial(near_axis, *tangent, *polarisation);
    EXPECT_FALSE(TryParallelTransportStep(metric, initial, 1.0e-3));

    PolarisedGeodesicIntegratorD integrator(&metric);
    const PolarisedGeodesicIntegratorD::Result result = integrator.Integrate(initial);
    EXPECT_EQ(result.outcome, PolarisedGeodesicIntegratorD::Result::Outcome::kChartDomainExit);
    EXPECT_EQ(result.steps, 0);
    EXPECT_EQ(result.state.x.theta, initial.x.theta)
        << "a failed chart step must preserve the last represented state";

    PolarisedStateD invalid = initial;
    invalid.x.theta = kBoyerLindquistPoleMargin;
    const PolarisedGeodesicIntegratorD::Result invalid_result = integrator.Integrate(invalid);
    EXPECT_EQ(invalid_result.outcome,
              PolarisedGeodesicIntegratorD::Result::Outcome::kChartDomainExit);
    EXPECT_EQ(invalid_result.steps, 0);

    PolarisedStateD invalid_near_horizon = initial;
    invalid_near_horizon.x.r = metric.HorizonRadius() * 1.01;
    invalid_near_horizon.x.theta = kBoyerLindquistPoleMargin;
    const PolarisedGeodesicIntegratorD::Result invalid_near_horizon_result =
        integrator.Integrate(invalid_near_horizon);
    EXPECT_EQ(invalid_near_horizon_result.outcome,
              PolarisedGeodesicIntegratorD::Result::Outcome::kChartDomainExit)
        << "an unrepresented event cannot be promoted to physical capture by its radius";

    PolarisedStateD invalid_transport = initial;
    invalid_transport.k.r = std::numeric_limits<double>::infinity();
    const PolarisedGeodesicIntegratorD::Result invalid_transport_result =
        integrator.Integrate(invalid_transport);
    EXPECT_EQ(invalid_transport_result.outcome,
              PolarisedGeodesicIntegratorD::Result::Outcome::kConstraintFailure);
}

// E = -k_t = -(g_t nu k^nu); conserved by the t Killing vector.
double Energy(const KerrMetricD& m, const PolarisedStateD& s) {
    double g[4][4], g_inv[4][4];
    m.Evaluate(s.x, g, g_inv);
    return -(g[0][0] * s.k[0] + g[0][3] * s.k[3]);
}

// L_z = k_phi = g_phi nu k^nu; conserved by the phi Killing vector.
double AngularMomentum(const KerrMetricD& m, const PolarisedStateD& s) {
    double g[4][4], g_inv[4][4];
    m.Evaluate(s.x, g, g_inv);
    return g[3][0] * s.k[0] + g[3][3] * s.k[3];
}

// Carter constant Q = k_theta^2 + cos^2 theta (a^2(-E^2) + L_z^2/sin^2 theta)
// for a null geodesic (mu = 0). Reference: Carter (1968), Phys. Rev. 174, 1559.
double CarterConstant(const KerrMetricD& m, const PolarisedStateD& s) {
    double g[4][4], g_inv[4][4];
    m.Evaluate(s.x, g, g_inv);
    double p_theta = g[2][2] * s.k[2];  // k_theta = g_theta theta k^theta.
    double E = -(g[0][0] * s.k[0] + g[0][3] * s.k[3]);
    double Lz = g[3][0] * s.k[0] + g[3][3] * s.k[3];
    double a = m.spin();
    double costh = std::cos(s.x.theta);
    double sinth = std::sin(s.x.theta);
    return p_theta * p_theta + costh * costh * (-a * a * E * E + Lz * Lz / (sinth * sinth));
}

double NullCondition(const KerrMetricD& m, const PolarisedStateD& s) {
    double g[4][4], g_inv[4][4];
    m.Evaluate(s.x, g, g_inv);
    return InnerProductD(g, s.k, s.k);
}

double InnerProductAt(const KerrMetricD& m, const Vec4d& x, const Vec4d& u, const Vec4d& v) {
    double g[4][4], g_inv[4][4];
    m.Evaluate(x, g, g_inv);
    return InnerProductD(g, u, v);
}

// Finite-difference Christoffel from the metric's own Evaluate, the independent
// reference used to confirm the analytic connection the transport builds.
void NumericChristoffel(const KerrMetricD& m, const Vec4d& x, double G[4][4][4]) {
    double g0[4][4], g_inv[4][4];
    m.Evaluate(x, g0, g_inv);
    double dg[4][4][4];
    for (int sigma = 0; sigma < 4; ++sigma) {
        double h = 1e-6 * std::max(1.0, std::abs(x[sigma])) + 1e-8;
        Vec4d xp = x, xm = x;
        xp[sigma] += h;
        xm[sigma] -= h;
        double gp[4][4], gm[4][4], tmp[4][4];
        m.Evaluate(xp, gp, tmp);
        m.Evaluate(xm, gm, tmp);
        for (int mu = 0; mu < 4; ++mu)
            for (int nu = 0; nu < 4; ++nu) dg[sigma][mu][nu] = (gp[mu][nu] - gm[mu][nu]) / (2 * h);
    }
    for (int mu = 0; mu < 4; ++mu)
        for (int nu = 0; nu < 4; ++nu)
            for (int rho = 0; rho < 4; ++rho) {
                double s = 0.0;
                for (int sigma = 0; sigma < 4; ++sigma)
                    s += g_inv[mu][sigma] *
                         (dg[nu][sigma][rho] + dg[rho][sigma][nu] - dg[sigma][nu][rho]);
                G[mu][nu][rho] = 0.5 * s;
            }
}

}  // namespace

//==============================================================================
// The connection the transport uses is the correct one for the metric
//==============================================================================

// The analytic connection agrees with a finite-difference connection built from
// the same metric, which is why the transport uses KerrChristoffel rather than
// KerrMetricD::Christoffel (the latter disagrees with its own metric; see the
// transport header). Tolerance 1e-9: a central finite difference on a smooth
// metric resolves the connection to O(h^2) ~ 1e-12 in principle but is floored
// by cancellation near 1e-9 at these radii.
TEST(OracleConnection, AnalyticConnectionAgreesWithMetricDerivatives) {
    KerrMetricD metric(1.0, 0.998);
    const std::vector<Vec4d> points = {Vec4d(0.0, 8.0, 1.2, 0.3),
                                       Vec4d(0.0, 4.0, std::numbers::pi / 2 - 0.2, 1.0),
                                       Vec4d(0.0, 15.0, 0.7, 0.0)};
    for (const Vec4d& x : points) {
        double Ga[4][4][4], Gn[4][4][4];
        KerrChristoffel(metric, x, Ga);
        NumericChristoffel(metric, x, Gn);
        double max_diff = 0.0;
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
                for (int l = 0; l < 4; ++l)
                    max_diff = std::max(max_diff, std::abs(Ga[i][j][l] - Gn[i][j][l]));
        EXPECT_LT(max_diff, 1e-9) << "at r=" << x.r << " theta=" << x.theta;
    }
}

//==============================================================================
// The component form equals the Killing-Yano contraction (derivation check)
//==============================================================================

// Im(kappa_WP) = Y_uv k^u f^v with the Kerr Killing-Yano tensor. This pins the
// (A, B) component form to the geometric object whose conservation the header
// derives. Tolerance 1e-12: both sides are the same polynomial in double
// precision, so only round-off separates them.
TEST(WalkerPenrose, ImaginaryPartEqualsKillingYanoContraction) {
    const double a = 0.9;
    KerrMetricD metric(1.0, a);
    Vec4d x(0.0, 9.0, 1.1, 0.4);
    const std::optional<Vec4d> tangent = TryMakeNullTangent(metric, x, -1.0, 0.02, 0.05);
    ASSERT_TRUE(tangent);
    const Vec4d kk = *tangent;
    const std::optional<Vec4d> polarisation =
        TryMakeOrthonormalPolarisation(metric, x, kk, Vec4d(0.0, 0.0, 1.0, 0.0));
    ASSERT_TRUE(polarisation);
    const Vec4d ff = *polarisation;
    PolarisedStateD s(x, kk, ff);

    const double r = x.r, theta = x.theta;
    const double sinth = std::sin(theta), costh = std::cos(theta), sin2th = sinth * sinth;
    // Non-zero Kerr Killing-Yano components (lower indices), header block.
    const double Y_tr = -a * costh;
    const double Y_tth = a * r * sinth;
    const double Y_rph = -a * a * costh * sin2th;
    const double Y_thph = r * (r * r + a * a) * sinth;
    const double kY =
        Y_tr * (kk[0] * ff[1] - kk[1] * ff[0]) + Y_tth * (kk[0] * ff[2] - kk[2] * ff[0]) +
        Y_rph * (kk[1] * ff[3] - kk[3] * ff[1]) + Y_thph * (kk[2] * ff[3] - kk[3] * ff[2]);

    std::complex<double> kappa = WalkerPenroseConstant(s, a);
    EXPECT_NEAR(kappa.imag(), kY, 1e-12);
}

//==============================================================================
// (a) Orthogonality, normalisation, and geodesic invariants under transport
//==============================================================================

// Parallel transport preserves every inner product exactly, so f.k = 0 and
// f.f = 1 hold along the whole geodesic, and the tangent stays a null geodesic
// (E, L_z, Q, the null condition all conserved). Tolerance is the oracle-tier
// CPU null bar 1e-10 (constants::geodesic::kNullConditionTolCpu): these are the
// same conservation-of-a-quadratic-form conditions as the null constraint, and
// double-precision RK4 with adaptive stepping holds them to ~1e-11 here.
TEST(WalkerPenrose, PolarisationStaysOrthogonalAndNormalisedUnderTransport) {
    KerrMetricD metric(1.0, 0.998);
    Vec4d x0(0.0, 18.0, std::numbers::pi / 2 - 0.15, 0.0);
    const std::optional<Vec4d> tangent = TryMakeNullTangent(metric, x0, -1.0, 0.02, 0.04);
    ASSERT_TRUE(tangent);
    const Vec4d kk = *tangent;
    const std::optional<Vec4d> polarisation =
        TryMakeOrthonormalPolarisation(metric, x0, kk, Vec4d(0.0, 0.0, 1.0, 0.0));
    ASSERT_TRUE(polarisation);
    const Vec4d ff = *polarisation;
    PolarisedStateD s(x0, kk, ff);

    ASSERT_NEAR(InnerProductAt(metric, x0, ff, kk), 0.0, 1e-14);
    ASSERT_NEAR(InnerProductAt(metric, x0, ff, ff), 1.0, 1e-14);

    const double E0 = Energy(metric, s), Lz0 = AngularMomentum(metric, s),
                 Q0 = CarterConstant(metric, s);
    ASSERT_GT(std::abs(E0), 1e-6);

    PolarisedGeodesicIntegratorD integrator(&metric);
    double max_fk = 0, max_ff = 0, max_null = 0, max_E = 0, max_Lz = 0, max_Q = 0;
    PolarisedStateD cur = s;
    for (int i = 0; i < 400000; ++i) {
        if (cur.x.r <= metric.HorizonRadius() * 1.02 || cur.x.r > 60.0) break;
        const std::optional<PolarisedStateD> next =
            TryParallelTransportStep(metric, cur, integrator.AdaptiveStepSize(cur.x.r));
        ASSERT_TRUE(next);
        cur = *next;
        max_fk = std::max(max_fk, std::abs(InnerProductAt(metric, cur.x, cur.f, cur.k)));
        max_ff = std::max(max_ff, std::abs(InnerProductAt(metric, cur.x, cur.f, cur.f) - 1.0));
        max_null = std::max(max_null, std::abs(NullCondition(metric, cur)));
        max_E = std::max(max_E, std::abs(Energy(metric, cur) - E0) / std::abs(E0));
        max_Lz = std::max(max_Lz, std::abs(AngularMomentum(metric, cur) - Lz0) / std::abs(Lz0));
        max_Q = std::max(max_Q, std::abs(CarterConstant(metric, cur) - Q0) / std::abs(Q0));
    }
    EXPECT_LT(max_fk, k::geodesic::kNullConditionTolCpu);
    EXPECT_LT(max_ff, k::geodesic::kNullConditionTolCpu);
    EXPECT_LT(max_null, k::geodesic::kNullConditionTolCpu);
    EXPECT_LT(max_E, k::geodesic::kNullConditionTolCpu);
    EXPECT_LT(max_Lz, k::geodesic::kNullConditionTolCpu);
    EXPECT_LT(max_Q, k::geodesic::kNullConditionTolCpu);
}

//==============================================================================
// (b) Walker-Penrose constant conserved over full Kerr geodesics
//==============================================================================

// The Walker-Penrose constant is invariant along any Kerr null geodesic with a
// parallel-propagated polarisation. Tolerance kWalkerPenroseConservationTol =
// 1e-10 (the transport header's named oracle-tier bar): double-precision RK4
// with the horizon-adaptive step conserves the analytic invariant to ~1e-12,
// including the near-photon-sphere pass, so 1e-10 holds with two orders of
// margin.
TEST(WalkerPenrose, ConstantConservedAlongKerrGeodesics) {
    KerrMetricD metric(1.0, 0.9);

    struct Ray {
        const char* name;
        double r0, theta0, kr, ktheta, kphi;
    };
    // Impact parameters from weak deflection to a near-photon-sphere pass (the
    // kphi = 0.01 ray reaches r ~ 2.8, close to the a = 0.9 prograde photon
    // orbit) plus an off-equatorial ray so k^theta participates.
    const std::vector<Ray> rays = {
        {"weak prograde", 20.0, std::numbers::pi / 2, -1.0, 0.0, 0.02},
        {"near photon sphere", 20.0, std::numbers::pi / 2, -1.0, 0.0, 0.01},
        {"retrograde", 20.0, std::numbers::pi / 2, -1.0, 0.0, -0.03},
        {"off-equatorial", 18.0, std::numbers::pi / 2 - 0.2, -1.0, 0.02, 0.03},
    };

    for (const Ray& r : rays) {
        Vec4d x0(0.0, r.r0, r.theta0, 0.0);
        const std::optional<Vec4d> tangent = TryMakeNullTangent(metric, x0, r.kr, r.ktheta, r.kphi);
        ASSERT_TRUE(tangent) << r.name;
        const Vec4d kk = *tangent;
        const std::optional<Vec4d> polarisation =
            TryMakeOrthonormalPolarisation(metric, x0, kk, Vec4d(0.0, 0.0, 1.0, 0.0));
        ASSERT_TRUE(polarisation) << r.name;
        const Vec4d ff = *polarisation;
        PolarisedStateD s(x0, kk, ff);
        const std::complex<double> kappa0 = WalkerPenroseConstant(s, metric);
        ASSERT_GT(std::abs(kappa0), 1e-6) << r.name;

        PolarisedGeodesicIntegratorD integrator(&metric);
        const PolarisedGeodesicIntegratorD::Result result = integrator.Integrate(s);
        // A conserved constant is only meaningful over a completed geodesic.
        ASSERT_EQ(result.outcome, PolarisedGeodesicIntegratorD::Result::Outcome::kEscaped)
            << r.name << " did not complete at the escape surface";

        double max_drift = 0.0;
        PolarisedStateD cur = s;
        for (int i = 0; i < 400000; ++i) {
            if (cur.x.r <= metric.HorizonRadius() * 1.02 || cur.x.r > 60.0) break;
            const std::optional<PolarisedStateD> next =
                TryParallelTransportStep(metric, cur, integrator.AdaptiveStepSize(cur.x.r));
            ASSERT_TRUE(next) << r.name;
            cur = *next;
            const std::complex<double> kappa = WalkerPenroseConstant(cur, metric);
            max_drift = std::max(max_drift, std::abs(kappa - kappa0) / std::abs(kappa0));
        }
        EXPECT_LT(max_drift, kWalkerPenroseConservationTol) << r.name;
    }
}

//==============================================================================
// (c) Schwarzschild equatorial limit: rigid perpendicular transport
//==============================================================================

// In Schwarzschild an equatorial null geodesic keeps k^theta = 0, and a
// polarisation vector normal to the orbital plane (pure theta) stays normal:
// reflection symmetry theta -> pi - theta decouples it, so the in-plane
// components never grow and the covariant magnitude g_theta theta (f^theta)^2 =
// r^2 (f^theta)^2 is preserved, i.e. f^theta * r is constant. This is the
// analytic no-rotation expectation, and it also follows from Walker-Penrose
// conservation (kappa = -i r^3 k^phi f^theta with k^phi = L/r^2). In-plane
// leakage bound 1e-10 (kNullConditionTolCpu): it is zero analytically and only
// round-off drives it; f^theta * r drift bound 1e-9, the same round-off scaled
// by r.
TEST(WalkerPenrose, SchwarzschildEquatorialPerpendicularTransportsRigidly) {
    KerrMetricD metric(1.0, 0.0);
    Vec4d x0(0.0, 15.0, std::numbers::pi / 2, 0.0);
    const std::optional<Vec4d> tangent =
        TryMakeNullTangent(metric, x0, -1.0, 0.0, 0.03);  // Equatorial: k^theta = 0.
    ASSERT_TRUE(tangent);
    const Vec4d kk = *tangent;
    Vec4d ff(0.0, 0.0, 1.0 / x0.r, 0.0);  // Unit: r^2 (1/r)^2 = 1.

    ASSERT_NEAR(InnerProductAt(metric, x0, ff, kk), 0.0, 1e-14);
    ASSERT_NEAR(InnerProductAt(metric, x0, ff, ff), 1.0, 1e-14);

    PolarisedStateD cur(x0, kk, ff);
    const double f_theta_r0 = cur.f.theta * cur.x.r;
    double max_in_plane = 0.0, max_ftheta_r = 0.0;
    for (int i = 0; i < 400000; ++i) {
        if (cur.x.r <= metric.HorizonRadius() * 1.05 || cur.x.r > 50.0) break;
        const std::optional<PolarisedStateD> next = TryParallelTransportStep(metric, cur, 0.01);
        ASSERT_TRUE(next);
        cur = *next;
        max_in_plane =
            std::max({max_in_plane, std::abs(cur.f.t), std::abs(cur.f.r), std::abs(cur.f.phi)});
        max_ftheta_r = std::max(max_ftheta_r, std::abs(cur.f.theta * cur.x.r - f_theta_r0));
    }
    EXPECT_LT(max_in_plane, k::geodesic::kNullConditionTolCpu);
    EXPECT_LT(max_ftheta_r, 1e-9);
}

//==============================================================================
// (d) Flat-space limit: constant polarisation
//==============================================================================

// With M = 0 and a = 0 the Kerr chart is Minkowski in spherical coordinates.
// Parallel transport is path-independent and the polarisation is constant in
// Cartesian components, though the spherical components change with the moving
// frame. Bound 1e-9: the Cartesian components are constant analytically, and
// double-precision RK4 over the ray holds them to ~1e-11.
TEST(WalkerPenrose, FlatSpacePolarisationIsConstantInCartesian) {
    KerrMetricD metric(0.0, 0.0);
    Vec4d x0(0.0, 10.0, 1.1, 0.3);
    const std::optional<Vec4d> tangent = TryMakeNullTangent(metric, x0, -1.0, 0.02, 0.03);
    ASSERT_TRUE(tangent);
    const Vec4d kk = *tangent;
    const std::optional<Vec4d> polarisation =
        TryMakeOrthonormalPolarisation(metric, x0, kk, Vec4d(0.0, 1.0, 0.0, 0.0));
    ASSERT_TRUE(polarisation);
    const Vec4d ff = *polarisation;

    // Contravariant spherical -> Cartesian: v^i = (dx^i/dx^alpha_sph) v^alpha.
    auto to_cartesian = [](const Vec4d& X, const Vec4d& V) {
        double r = X.r, th = X.theta, ph = X.phi;
        double st = std::sin(th), ct = std::cos(th), sp = std::sin(ph), cp = std::cos(ph);
        double vx = st * cp * V.r + r * ct * cp * V.theta - r * st * sp * V.phi;
        double vy = st * sp * V.r + r * ct * sp * V.theta + r * st * cp * V.phi;
        double vz = ct * V.r - r * st * V.theta;
        return Vec4d(V.t, vx, vy, vz);
    };

    PolarisedStateD cur(x0, kk, ff);
    const Vec4d f_cart0 = to_cartesian(cur.x, cur.f);
    double max_dev = 0.0;
    for (int i = 0; i < 400000; ++i) {
        if (cur.x.r > 50.0 || cur.x.r < 0.5) break;
        const std::optional<PolarisedStateD> next = TryParallelTransportStep(metric, cur, 0.01);
        ASSERT_TRUE(next);
        cur = *next;
        const Vec4d f_cart = to_cartesian(cur.x, cur.f);
        for (int j = 0; j < 4; ++j) max_dev = std::max(max_dev, std::abs(f_cart[j] - f_cart0[j]));
    }
    EXPECT_LT(max_dev, 1e-9);
}

}  // namespace sirius::test
