// Geodesic integrators. Two live methods, exactly as the retired .cu carried
// them:
//
//   1. The time-transformed explicit symplectic integrator (TTESI) for the
//      Kerr Hamiltonian: a Yoshida 4th-order composition of leapfrog steps in
//      Preto-Saha transformed time, with the .cu's adaptive step and
//      Hamiltonian-error acceptance test. This is the conservation-preserving
//      candidate integrator (RDOP002A.cu integrateGeodesicSymplectic; the .cu
//      annotates it as not-yet-on-the-render-path, kept for parity here).
//
//   2. The Cartesian RK4 geodesic step over the Kerr-Schild family, which is
//      the render path (RDOP002A.cu integrateGeodesicRK4Cart), used by
//      trace.slang.
//
// Yoshida coefficients, the tau step, tolerance factors, and every guard
// constant are copied verbatim from RDOP002A.cu. References: Yoshida (1990)
// Phys. Lett. A 150, 262; Preto & Saha (1999) ApJ 514, 272.
module gr_integrator;

import gr_types;
import gr_metrics;

// Forest-Ruth / Yoshida composition coefficients (RDOP002A.cu YOSHIDA_*).
static const float YOSHIDA_C1 = 1.3512071919596578f;   // 1/(2 - 2^(1/3))
static const float YOSHIDA_C2 = -1.7024143839193153f;  // -2^(1/3)/(2 - 2^(1/3))

// --- Kerr Hamiltonian helpers (RDOP002A.cu computeKerrParams etc.) ------------

struct KerrHamiltonianParams {
    float Sigma;
    float Delta;
    float A_func;
    float sin_th;
    float cos_th;
    float sin2_th;
    float cos2_th;
}

KerrHamiltonianParams ComputeKerrParams(float r, float theta, float M, float a) {
    KerrHamiltonianParams kp;
    kp.sin_th = sin(theta);
    kp.cos_th = cos(theta);
    kp.sin2_th = kp.sin_th * kp.sin_th;
    kp.cos2_th = kp.cos_th * kp.cos_th;
    float r2 = r * r;
    float a2 = a * a;
    kp.Sigma = r2 + a2 * kp.cos2_th;
    kp.Delta = r2 - 2.0f * M * r + a2;
    float r2_a2 = r2 + a2;
    kp.A_func = r2_a2 * r2_a2 - a2 * kp.Delta * kp.sin2_th;
    return kp;
}

// K2: radial kinetic term, position-dependent kinetic split.
void FlowK2Radial(inout HamiltonianState state, float h, float M, float a) {
    float r = state.q[1];
    float p_r = state.p[1];
    float r2 = r * r;
    float a2 = a * a;
    float Delta = r2 - 2.0f * M * r + a2;
    float dDelta_dr = 2.0f * r - 2.0f * M;
    float Delta_mid = Delta;
    state.p[1] -= 0.5f * h * (0.5f * dDelta_dr * p_r * p_r);
    state.q[1] += h * Delta_mid * state.p[1];
    r = state.q[1];
    r2 = r * r;
    Delta = r2 - 2.0f * M * r + a2;
    dDelta_dr = 2.0f * r - 2.0f * M;
    state.p[1] -= 0.5f * h * (0.5f * dDelta_dr * state.p[1] * state.p[1]);
}

// K3: free particle in theta with trans-polar chart switch.
void FlowK3Angular(inout HamiltonianState state, float h) {
    state.q[2] += h * state.p[2];
    const float PI_VAL = 3.14159265359f;
    const float TWO_PI_VAL = 6.28318530718f;
    if (state.q[2] < 0.0f) {
        state.q[2] = -state.q[2];
        state.p[2] = -state.p[2];
        state.q[3] += PI_VAL;
    } else if (state.q[2] > PI_VAL) {
        state.q[2] = TWO_PI_VAL - state.q[2];
        state.p[2] = -state.p[2];
        state.q[3] += PI_VAL;
    }
}

// K4: effective radial-potential kick.
void FlowK4RadialPotential(inout HamiltonianState state, float h, float M, float a, float E, float Phi) {
    float r = state.q[1];
    float r2 = r * r;
    float a2 = a * a;
    float Delta = r2 - 2.0f * M * r + a2;
    float Delta_safe = max(Delta, 0.01f);
    float r2_a2 = r2 + a2;
    float term1 = E * r2_a2 - a * Phi;
    float dDelta_dr = 2.0f * r - 2.0f * M;
    float dVr_dr = -2.0f * E * 2.0f * r * term1 / Delta_safe
                  + term1 * term1 * dDelta_dr / (Delta_safe * Delta_safe);
    state.p[1] -= h * dVr_dr;
}

// K5: effective angular-potential kick (null geodesic, mu = 0).
void FlowK5AngularPotential(inout HamiltonianState state, float h, float a, float E, float Phi) {
    float theta = state.q[2];
    float sin_th = sin(theta);
    float sin2_th = sin_th * sin_th;
    float sin2_safe = max(sin2_th, 0.0001f);
    float a2 = a * a;
    float dVth_dth = 2.0f * a2 * cos(theta) * sin_th * E * E
                   - 2.0f * Phi * Phi * cos(theta) / (sin2_safe * sin_th);
    state.p[2] -= h * dVth_dth;
}

// t and phi drift from the time transformation.
void FlowTimeAzimuth(inout HamiltonianState state, float h, float M, float a) {
    float r = state.q[1];
    float theta = state.q[2];
    KerrHamiltonianParams kp = ComputeKerrParams(r, theta, M, a);
    float E = -state.p[0];
    float Phi = state.p[3];
    float a2 = a * a;
    float r2 = r * r;
    float r2_a2 = r2 + a2;
    float Delta_safe = max(kp.Delta, 0.01f);
    float sin2_safe = max(kp.sin2_th, 0.0001f);
    float dH_dpt = (r2_a2 * r2_a2 / (kp.Sigma * Delta_safe)) * E
                 - (a / sin2_safe) * Phi;
    float dH_dphi = (a * r2_a2 / (kp.Sigma * Delta_safe)) * E
                  + (1.0f / (kp.Sigma * sin2_safe)) * Phi;
    state.q[0] += h * kp.Sigma * dH_dpt;
    state.q[3] += h * kp.Sigma * dH_dphi;
}

// Second-order leapfrog S2(h): half kick, full drift, half kick, then t/phi.
void IntegrateSymplecticS2(inout HamiltonianState state, float h, float M, float a) {
    float E = -state.p[0];
    float Phi = state.p[3];
    FlowK4RadialPotential(state, 0.5f * h, M, a, E, Phi);
    FlowK5AngularPotential(state, 0.5f * h, a, E, Phi);
    FlowK2Radial(state, h, M, a);
    FlowK3Angular(state, h);
    FlowK4RadialPotential(state, 0.5f * h, M, a, E, Phi);
    FlowK5AngularPotential(state, 0.5f * h, a, E, Phi);
    FlowTimeAzimuth(state, h, M, a);
}

// Yoshida 4th-order composition S4 = S2(c1 h) S2(c2 h) S2(c1 h).
public void IntegrateSymplecticS4(inout HamiltonianState state, float h, float M, float a) {
    IntegrateSymplecticS2(state, YOSHIDA_C1 * h, M, a);
    IntegrateSymplecticS2(state, YOSHIDA_C2 * h, M, a);
    IntegrateSymplecticS2(state, YOSHIDA_C1 * h, M, a);
}

// Preto-Saha time transformation g = Sigma.
float ComputeTimeTransformG(float r, float theta, float a) {
    float a2 = a * a;
    float cos_th = cos(theta);
    return r * r + a2 * cos_th * cos_th;
}

// Adaptive coordinate-time step: dt = g(r,theta) * tau * horizon_factor.
float ComputeAdaptiveSymplecticStep(HamiltonianState state, float tau_step, float M, float a,
                                    float minStep, float maxStep) {
    float r = state.q[1];
    float theta = state.q[2];
    float g = ComputeTimeTransformG(r, theta, a);
    float r2 = r * r;
    float Delta = r2 - 2.0f * M * r + a * a;
    float horizon_factor = max(Delta / (r2 + 0.1f), 0.1f);
    float dt = g * tau_step * horizon_factor;
    return clamp(dt, minStep, maxStep);
}

// Normalised Hamiltonian-constraint residual |2 Sigma H| / (1 + E^2).
float ComputeHamiltonianError(HamiltonianState state, float M, float a) {
    float r = state.q[1];
    float theta = state.q[2];
    float p_r = state.p[1];
    float p_th = state.p[2];
    float E = -state.p[0];
    float Phi = state.p[3];
    KerrHamiltonianParams kp = ComputeKerrParams(r, theta, M, a);
    float a2 = a * a;
    float r2 = r * r;
    float r2_a2 = r2 + a2;
    float Delta_safe = max(kp.Delta, 0.01f);
    float sin2_safe = max(kp.sin2_th, 0.0001f);
    float term1 = -(r2_a2 * E - a * Phi) * (r2_a2 * E - a * Phi) / Delta_safe;
    float term2 = kp.Delta * p_r * p_r;
    float term3 = p_th * p_th;
    float term4 = (a * E - Phi / kp.sin_th) * (a * E - Phi / kp.sin_th) / sin2_safe;
    float H_2Sigma = term1 + term2 + term3 + term4;
    return abs(H_2Sigma) / (1.0f + E * E);
}

// GeodesicState (q, u) -> HamiltonianState (q, p = g u), Boyer-Lindquist Kerr.
public HamiltonianState GeodesicToHamiltonian(GeodesicState gs, float M, float a) {
    HamiltonianState hs = ZeroHamiltonian();
    hs.q[0] = gs.x.t;
    hs.q[1] = gs.x.r;
    hs.q[2] = gs.x.theta;
    hs.q[3] = gs.x.phi;

    float r = gs.x.r;
    float theta = gs.x.theta;
    float sin_th = sin(theta);
    float cos_th = cos(theta);
    float r2 = r * r;
    float a2 = a * a;

    float Sigma = r2 + a2 * cos_th * cos_th;
    float Delta = r2 - 2.0f * M * r + a2;
    float sin2_th = sin_th * sin_th;

    float g_tt = -(1.0f - 2.0f * M * r / Sigma);
    float g_rr = Sigma / max(Delta, 0.01f);
    float g_thth = Sigma;
    float g_phph = (r2 + a2 + 2.0f * M * r * a2 * sin2_th / Sigma) * sin2_th;
    float g_tph = -2.0f * M * r * a * sin2_th / Sigma;

    float u[4] = { gs.u.t, gs.u.r, gs.u.theta, gs.u.phi };
    hs.p[0] = g_tt * u[0] + g_tph * u[3];
    hs.p[1] = g_rr * u[1];
    hs.p[2] = g_thth * u[2];
    hs.p[3] = g_tph * u[0] + g_phph * u[3];
    return hs;
}

// HamiltonianState (q, p) -> GeodesicState (q, u = g^{-1} p), Boyer-Lindquist.
public GeodesicState HamiltonianToGeodesic(HamiltonianState hs, float M, float a) {
    GeodesicState gs;
    gs.x = MakeVec4(hs.q[0], hs.q[1], hs.q[2], hs.q[3]);

    float r = hs.q[1];
    float theta = hs.q[2];
    float sin_th = sin(theta);
    float cos_th = cos(theta);
    float r2 = r * r;
    float a2 = a * a;

    float Sigma = r2 + a2 * cos_th * cos_th;
    float Delta = max(r2 - 2.0f * M * r + a2, 0.01f);
    float sin2_th = max(sin_th * sin_th, 0.0001f);

    float r2_a2 = r2 + a2;
    float A = r2_a2 * r2_a2 - a2 * Delta * sin2_th;

    float g_inv_tt = -A / (Sigma * Delta);
    float g_inv_rr = Delta / Sigma;
    float g_inv_thth = 1.0f / Sigma;
    float g_inv_phph = (Delta - a2 * sin2_th) / (Sigma * Delta * sin2_th);
    float g_inv_tph = -2.0f * M * r * a / (Sigma * Delta);

    gs.u = MakeVec4(
        g_inv_tt * hs.p[0] + g_inv_tph * hs.p[3],
        g_inv_rr * hs.p[1],
        g_inv_thth * hs.p[2],
        g_inv_tph * hs.p[0] + g_inv_phph * hs.p[3]);
    return gs;
}

// Adaptive symplectic step (RDOP002A.cu integrateGeodesicSymplectic).
// h is in/out (the step-size controller), stepAccepted reports the
// Hamiltonian-conservation acceptance test.
public GeodesicState IntegrateGeodesicSymplectic(
    GeodesicState state, inout float h, float M, float a,
    float tolerance, float minStep, float maxStep, out bool stepAccepted)
{
    HamiltonianState hs = GeodesicToHamiltonian(state, M, a);
    float H_error_initial = ComputeHamiltonianError(hs, M, a);

    float tau_step = 0.1f;
    float effective_h = ComputeAdaptiveSymplecticStep(hs, tau_step, M, a, minStep, maxStep);

    HamiltonianState hs_new = hs;
    IntegrateSymplecticS4(hs_new, effective_h, M, a);

    float H_error_final = ComputeHamiltonianError(hs_new, M, a);
    float H_error_change = abs(H_error_final - H_error_initial);

    stepAccepted = (H_error_change < tolerance * 10.0f);
    if (!stepAccepted) {
        h = max(h * 0.5f, minStep);
        return state;
    }
    if (H_error_change < tolerance * 0.1f) {
        h = min(h * 1.5f, maxStep);
    }
    h = effective_h;

    GeodesicState result = HamiltonianToGeodesic(hs_new, M, a);

    const float THETA_MIN = 0.001f;
    const float THETA_MAX = PI - 0.001f;
    result.x.theta = max(THETA_MIN, min(THETA_MAX, result.x.theta));
    result.x.phi = fmod(result.x.phi, 2.0f * PI);
    if (result.x.phi < 0.0f) result.x.phi += 2.0f * PI;
    return result;
}

// --- Cartesian RK4 render integrator (RDOP002A.cu integrateGeodesicRK4Cart) --

// du^mu/dlambda = -Gamma^mu_nu_rho u^nu u^rho, Kerr-Schild Cartesian
// Christoffels. (M, a, Q, Lambda) let the render path carry charge and the
// cosmological sector; the legacy accel fixed Q = Lambda = 0.
public Vec4Cart GeodesicAccelerationCart(GeodesicStateCart state, float M, float a, float Q, float Lambda) {
    Real Gamma[4][4][4];
    GetKerrSchildFamilyChristoffel(state.x, M, a, Q, Lambda, Gamma);
    Real u[4] = { state.u.t, state.u.x, state.u.y, state.u.z };
    Real acc[4] = { 0.0f, 0.0f, 0.0f, 0.0f };
    [unroll] for (int mu = 0; mu < 4; mu++)
        [unroll] for (int nu = 0; nu < 4; nu++)
            [unroll] for (int rho = 0; rho < 4; rho++)
                acc[mu] -= Gamma[mu][nu][rho] * u[nu] * u[rho];
    return MakeVec4Cart(acc[0], acc[1], acc[2], acc[3]);
}

public GeodesicStateCart IntegrateGeodesicRK4Cart(GeodesicStateCart state, Real dt,
                                                  float M, float a, float Q, float Lambda) {
    Vec4Cart dx1 = ScaleVec4Cart(state.u, dt);
    Vec4Cart du1 = ScaleVec4Cart(GeodesicAccelerationCart(state, M, a, Q, Lambda), dt);

    GeodesicStateCart s2;
    s2.x = AddVec4Cart(state.x, ScaleVec4Cart(dx1, 0.5f));
    s2.u = AddVec4Cart(state.u, ScaleVec4Cart(du1, 0.5f));
    Vec4Cart dx2 = ScaleVec4Cart(s2.u, dt);
    Vec4Cart du2 = ScaleVec4Cart(GeodesicAccelerationCart(s2, M, a, Q, Lambda), dt);

    GeodesicStateCart s3;
    s3.x = AddVec4Cart(state.x, ScaleVec4Cart(dx2, 0.5f));
    s3.u = AddVec4Cart(state.u, ScaleVec4Cart(du2, 0.5f));
    Vec4Cart dx3 = ScaleVec4Cart(s3.u, dt);
    Vec4Cart du3 = ScaleVec4Cart(GeodesicAccelerationCart(s3, M, a, Q, Lambda), dt);

    GeodesicStateCart s4;
    s4.x = AddVec4Cart(state.x, dx3);
    s4.u = AddVec4Cart(state.u, du3);
    Vec4Cart dx4 = ScaleVec4Cart(s4.u, dt);
    Vec4Cart du4 = ScaleVec4Cart(GeodesicAccelerationCart(s4, M, a, Q, Lambda), dt);

    GeodesicStateCart result;
    float sixth = 1.0f / 6.0f;
    result.x = AddVec4Cart(state.x,
        ScaleVec4Cart(AddVec4Cart(AddVec4Cart(dx1, ScaleVec4Cart(dx2, 2.0f)),
                                  AddVec4Cart(ScaleVec4Cart(dx3, 2.0f), dx4)), sixth));
    result.u = AddVec4Cart(state.u,
        ScaleVec4Cart(AddVec4Cart(AddVec4Cart(du1, ScaleVec4Cart(du2, 2.0f)),
                                  AddVec4Cart(ScaleVec4Cart(du3, 2.0f), du4)), sixth));
    return result;
}

// --- Alcubierre warp Cartesian RK4 render integrator -------------------------
// The warp family (gr_metrics GetWarpDrive*) shares the Cartesian chart and the
// same du^mu/dlambda = -Gamma^mu_nu_rho u^nu u^rho law as the Kerr-Schild path,
// so the render path reuses the RK4 structure with the warp Christoffels. The
// bubble is centred at (xs, ys, zs); (vs, sigma, R) are the warp velocity, wall
// thickness, and bubble radius. This is lensing-only: the warp bubble has no
// horizon and no disk, so trace.slang integrates and shades the background.
public Vec4Cart GeodesicAccelerationWarpCart(GeodesicStateCart state, float vs, float sigma,
                                             float R, float xs, float ys, float zs) {
    Real Gamma[4][4][4];
    GetWarpDriveChristoffel(state.x, vs, sigma, R, xs, ys, zs, Gamma);
    Real u[4] = { state.u.t, state.u.x, state.u.y, state.u.z };
    Real acc[4] = { 0.0f, 0.0f, 0.0f, 0.0f };
    [unroll] for (int mu = 0; mu < 4; mu++)
        [unroll] for (int nu = 0; nu < 4; nu++)
            [unroll] for (int rho = 0; rho < 4; rho++)
                acc[mu] -= Gamma[mu][nu][rho] * u[nu] * u[rho];
    return MakeVec4Cart(acc[0], acc[1], acc[2], acc[3]);
}

public GeodesicStateCart IntegrateGeodesicRK4WarpCart(GeodesicStateCart state, Real dt,
                                                      float vs, float sigma, float R,
                                                      float xs, float ys, float zs) {
    Vec4Cart dx1 = ScaleVec4Cart(state.u, dt);
    Vec4Cart du1 = ScaleVec4Cart(GeodesicAccelerationWarpCart(state, vs, sigma, R, xs, ys, zs), dt);

    GeodesicStateCart s2;
    s2.x = AddVec4Cart(state.x, ScaleVec4Cart(dx1, 0.5f));
    s2.u = AddVec4Cart(state.u, ScaleVec4Cart(du1, 0.5f));
    Vec4Cart dx2 = ScaleVec4Cart(s2.u, dt);
    Vec4Cart du2 = ScaleVec4Cart(GeodesicAccelerationWarpCart(s2, vs, sigma, R, xs, ys, zs), dt);

    GeodesicStateCart s3;
    s3.x = AddVec4Cart(state.x, ScaleVec4Cart(dx2, 0.5f));
    s3.u = AddVec4Cart(state.u, ScaleVec4Cart(du2, 0.5f));
    Vec4Cart dx3 = ScaleVec4Cart(s3.u, dt);
    Vec4Cart du3 = ScaleVec4Cart(GeodesicAccelerationWarpCart(s3, vs, sigma, R, xs, ys, zs), dt);

    GeodesicStateCart s4;
    s4.x = AddVec4Cart(state.x, dx3);
    s4.u = AddVec4Cart(state.u, du3);
    Vec4Cart dx4 = ScaleVec4Cart(s4.u, dt);
    Vec4Cart du4 = ScaleVec4Cart(GeodesicAccelerationWarpCart(s4, vs, sigma, R, xs, ys, zs), dt);

    GeodesicStateCart result;
    float sixth = 1.0f / 6.0f;
    result.x = AddVec4Cart(state.x,
        ScaleVec4Cart(AddVec4Cart(AddVec4Cart(dx1, ScaleVec4Cart(dx2, 2.0f)),
                                  AddVec4Cart(ScaleVec4Cart(dx3, 2.0f), dx4)), sixth));
    result.u = AddVec4Cart(state.u,
        ScaleVec4Cart(AddVec4Cart(AddVec4Cart(du1, ScaleVec4Cart(du2, 2.0f)),
                                  AddVec4Cart(ScaleVec4Cart(du3, 2.0f), du4)), sixth));
    return result;
}
