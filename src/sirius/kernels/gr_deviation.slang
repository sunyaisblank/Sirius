// Ray-bundle geodesic deviation (DNGR's distinguishing technique). Two paths,
// exactly the two the .cu carried:
//
//   1. Full-Riemann: R^mu_nu_rho_sigma finite-differenced from the analytic
//      Christoffels, contracted D^2 xi^mu/dl^2 = R^mu_nu_rho_sigma k^nu xi^rho
//      k^sigma (RDOP002A.cu getRiemannTensor + geodesicDeviationAccel,
//      DNGR eq A.18-A.19).
//
//   2. Tidal approximation (the .cu's default fast mode): leading-order tidal
//      stretch xi *= M/r^3 * step (RDOP002A.cu render loop, "TIDAL
//      APPROXIMATION" branch).
//
// Chart choice, deliberate and documented: the .cu evaluated the Riemann path
// through a spherical<->Cartesian Christoffel transform while the render
// integrated in one chart, mixing charts within a step. This port evaluates
// deviation in the same Kerr-Schild Cartesian chart the render integrator uses
// (gr_integrator IntegrateGeodesicRK4Cart), so k, xi, and R share one chart.
// The finite-difference step and the Riemann formula are otherwise copied
// verbatim.
module gr_deviation;

import gr_types;
import gr_metrics;

// Riemann tensor from central differences of the Kerr-Schild Christoffels.
// R[mu][nu][rho][sigma]. Time derivatives vanish (stationary spacetime),
// matching the .cu.
public void GetRiemannTensorCart(
    Vec4Cart x, float M, float a, float Q, float Lambda,
    out Real R[4][4][4][4])
{
    [unroll] for (int mu = 0; mu < 4; mu++)
        [unroll] for (int nu = 0; nu < 4; nu++)
            [unroll] for (int rho = 0; rho < 4; rho++)
                [unroll] for (int sig = 0; sig < 4; sig++)
                    R[mu][nu][rho][sig] = 0.0f;

    Real Gamma[4][4][4];
    GetKerrSchildFamilyChristoffel(x, M, a, Q, Lambda, Gamma);

    Real r = ComputeKerrRadiusCart(x.x, x.y, x.z, Real(a));
    const Real eps = 0.001f * max(Real(1.0f), r);

    // dGamma[rho][mu][nu][sigma] = d_rho Gamma^mu_nu_sigma.
    Real dGamma[4][4][4][4];
    [unroll] for (int mu = 0; mu < 4; mu++)
        [unroll] for (int nu = 0; nu < 4; nu++)
            [unroll] for (int sig = 0; sig < 4; sig++)
                dGamma[0][mu][nu][sig] = 0.0f;

    for (int rho = 1; rho < 4; rho++) {
        Vec4Cart xp = x;
        Vec4Cart xm = x;
        if (rho == 1) { xp.x += eps; xm.x -= eps; }
        else if (rho == 2) { xp.y += eps; xm.y -= eps; }
        else { xp.z += eps; xm.z -= eps; }

        Real Gamma_plus[4][4][4];
        Real Gamma_minus[4][4][4];
        GetKerrSchildFamilyChristoffel(xp, M, a, Q, Lambda, Gamma_plus);
        GetKerrSchildFamilyChristoffel(xm, M, a, Q, Lambda, Gamma_minus);

        for (int mu = 0; mu < 4; mu++)
            for (int nu = 0; nu < 4; nu++)
                for (int sig = 0; sig < 4; sig++)
                    dGamma[rho][mu][nu][sig] =
                        (Gamma_plus[mu][nu][sig] - Gamma_minus[mu][nu][sig]) / (2.0f * eps);
    }

    // R^mu_nu_rho_sigma = d_rho Gamma^mu_nu_sigma - d_sigma Gamma^mu_nu_rho
    //                   + Gamma^mu_lam_rho Gamma^lam_nu_sigma
    //                   - Gamma^mu_lam_sigma Gamma^lam_nu_rho.
    for (int mu = 0; mu < 4; mu++)
        for (int nu = 0; nu < 4; nu++)
            for (int rho = 0; rho < 4; rho++)
                for (int sig = 0; sig < 4; sig++) {
                    Real val = dGamma[rho][mu][nu][sig] - dGamma[sig][mu][nu][rho];
                    for (int lam = 0; lam < 4; lam++)
                        val += Gamma[mu][lam][rho] * Gamma[lam][nu][sig]
                             - Gamma[mu][lam][sig] * Gamma[lam][nu][rho];
                    R[mu][nu][rho][sig] = val;
                }
}

// Contract a precomputed Riemann tensor with one deviation vector. A beam has
// two vectors, so separating construction from contraction prevents computing
// the same 256-component tensor twice per geodesic step.
public Real4 GeodesicDeviationAccelFromRiemann(
    Real R[4][4][4][4], Vec4Cart k, Real4 xi) {
    Real k_arr[4] = { k.t, k.x, k.y, k.z };
    Real xi_arr[4] = { xi.x, xi.y, xi.z, xi.w };
    Real accel[4] = { 0.0f, 0.0f, 0.0f, 0.0f };

    for (int mu = 0; mu < 4; mu++)
        for (int nu = 0; nu < 4; nu++)
            for (int rho = 0; rho < 4; rho++)
                for (int sig = 0; sig < 4; sig++)
                    accel[mu] += R[mu][nu][rho][sig] * k_arr[nu] * xi_arr[rho] * k_arr[sig];

    return Real4(accel[0], accel[1], accel[2], accel[3]);
}

// Full-Riemann deviation acceleration D^2 xi^mu/dl^2 (DNGR eq A.19).
public Real4 GeodesicDeviationAccelFullRiemann(
    Vec4Cart x, Vec4Cart k, Real4 xi, float M, float a, float Q, float Lambda)
{
    Real R[4][4][4][4];
    GetRiemannTensorCart(x, M, a, Q, Lambda, R);
    return GeodesicDeviationAccelFromRiemann(R, k, xi);
}

// Tidal-approximation deviation acceleration (the .cu default fast mode):
// leading-order Riemann scale M/r^3 times the step.
public float4 TidalDeviationAccel(float4 xi, float M, float r, float step) {
    float tidalFactor = M / (r * r * r);
    float tf = tidalFactor * step;
    return float4(xi.x * tf, xi.y * tf, xi.z * tf, xi.w * tf);
}

// Singular-value geometry of the propagated transverse 2x2 Jacobi map. The
// semi-axes are its singular values and the orientation is the major
// eigenvector of J J^T, so a rotation of the input basis cannot rotate the
// rendered output ellipse. Shared by trace.slang and the parity probe to keep
// this final P2 projection a single device authority.
public void BeamEllipseFromTransverseMap(float beamA, float beamB, float beamC, float beamD,
                                         out float semiMajor, out float semiMinor,
                                         out float orientation, out float determinant) {
    float mapDeterminant = beamA * beamD - beamB * beamC;
    float mapPower = beamA * beamA + beamB * beamB + beamC * beamC + beamD * beamD;
    float discriminant = sqrt(max(mapPower * mapPower -
                                  4.0f * mapDeterminant * mapDeterminant, 0.0f));
    semiMajor = sqrt(max(0.5f * (mapPower + discriminant), 0.0f));
    semiMinor = sqrt(max(0.5f * (mapPower - discriminant), 0.0f));
    orientation = 0.5f * atan2(2.0f * (beamA * beamC + beamB * beamD),
                               beamA * beamA + beamB * beamB -
                               beamC * beamC - beamD * beamD);
    determinant = mapDeterminant;
}
