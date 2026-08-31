// Ray-bundle geodesic deviation (DNGR's distinguishing technique). The sole
// represented path contracts the full Riemann tensor:
//
//   R^mu_nu_rho_sigma is finite-differenced from the analytic
//      Christoffels, contracted D^2 xi^mu/dl^2 = -R^mu_nu_rho_sigma k^nu xi^rho
//      k^sigma (RDOP002A.cu getRiemannTensor + geodesicDeviationAccel,
//      DNGR eq A.18-A.19).
//
// Chart choice, deliberate and documented: the .cu evaluated the Riemann path
// through a spherical<->Cartesian Christoffel transform while the render
// integrated in one chart, mixing charts within a step. This port evaluates
// deviation in the same Kerr-Schild Cartesian chart the render integrator uses
// (gr_integrator IntegrateGeodesicRK4Cart), so k, xi, and R share one chart.
// The Riemann formula is retained, while its connection derivative uses a
// precision-matched five-point stencil pinned against an independent analytic
// Schwarzschild tensor and closed-form radial congruence by non-image Vulkan
// parity probes.
module gr_deviation;

import gr_types;
import gr_metrics;

// Riemann tensor from fourth-order differences of the Kerr-Schild Christoffels.
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
    // A fourth-order five-point stencil uses a precision-matched step. At
    // fp32, the former 0.1% two-point subtraction lost up to 1.4% of an
    // analytic Schwarzschild tidal contraction to cancellation.
#ifdef SIRIUS_FP64
    const Real eps = Real(0.001) * max(Real(1.0), r);
#else
    const Real eps = Real(0.02) * max(Real(1.0), r);
#endif

    // dGamma[rho][mu][nu][sigma] = d_rho Gamma^mu_nu_sigma.
    Real dGamma[4][4][4][4];
    [unroll] for (int mu = 0; mu < 4; mu++)
        [unroll] for (int nu = 0; nu < 4; nu++)
            [unroll] for (int sig = 0; sig < 4; sig++)
                dGamma[0][mu][nu][sig] = 0.0f;

    for (int rho = 1; rho < 4; rho++) {
        Vec4Cart xp = x;
        Vec4Cart xm = x;
        Vec4Cart xpp = x;
        Vec4Cart xmm = x;
        if (rho == 1) {
            xp.x += eps;
            xm.x -= eps;
            xpp.x += 2.0f * eps;
            xmm.x -= 2.0f * eps;
        } else if (rho == 2) {
            xp.y += eps;
            xm.y -= eps;
            xpp.y += 2.0f * eps;
            xmm.y -= 2.0f * eps;
        } else {
            xp.z += eps;
            xm.z -= eps;
            xpp.z += 2.0f * eps;
            xmm.z -= 2.0f * eps;
        }

        Real Gamma_plus[4][4][4];
        Real Gamma_minus[4][4][4];
        GetKerrSchildFamilyChristoffel(xp, M, a, Q, Lambda, Gamma_plus);
        GetKerrSchildFamilyChristoffel(xm, M, a, Q, Lambda, Gamma_minus);

        for (int mu = 0; mu < 4; mu++)
            for (int nu = 0; nu < 4; nu++)
                for (int sig = 0; sig < 4; sig++)
                    dGamma[rho][mu][nu][sig] =
                        8.0f * (Gamma_plus[mu][nu][sig] - Gamma_minus[mu][nu][sig]);

        // Reuse the two Christoffel work arrays so the higher-order stencil
        // does not double the shader's live local-array footprint.
        GetKerrSchildFamilyChristoffel(xpp, M, a, Q, Lambda, Gamma_plus);
        GetKerrSchildFamilyChristoffel(xmm, M, a, Q, Lambda, Gamma_minus);
        for (int mu = 0; mu < 4; mu++)
            for (int nu = 0; nu < 4; nu++)
                for (int sig = 0; sig < 4; sig++)
                    dGamma[rho][mu][nu][sig] =
                        (dGamma[rho][mu][nu][sig] -
                         (Gamma_plus[mu][nu][sig] - Gamma_minus[mu][nu][sig])) /
                        (12.0f * eps);
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
                    accel[mu] -= R[mu][nu][rho][sig] * k_arr[nu] * xi_arr[rho] * k_arr[sig];

    return Real4(accel[0], accel[1], accel[2], accel[3]);
}

// Coordinate form of the covariant Jacobi system. V = D xi/dlambda is a
// vector, so both xi and V carry connection terms in coordinate components.
// Advancing xi'' = -R(k,xi)k as an ordinary Cartesian second derivative drops
// these terms and is not equivalent outside an inertial chart.
public void GeodesicDeviationRhs(
    Real Gamma[4][4][4], Real R[4][4][4][4], Vec4Cart k, Real4 xi, Real4 V,
    out Real4 dXi, out Real4 dV) {
    Real k_arr[4] = { k.t, k.x, k.y, k.z };
    Real xi_arr[4] = { xi.x, xi.y, xi.z, xi.w };
    Real V_arr[4] = { V.x, V.y, V.z, V.w };
    Real d_xi[4] = { 0.0f, 0.0f, 0.0f, 0.0f };
    Real d_v[4] = { 0.0f, 0.0f, 0.0f, 0.0f };

    for (int mu = 0; mu < 4; ++mu) {
        d_xi[mu] = V_arr[mu];
        for (int alpha = 0; alpha < 4; ++alpha) {
            for (int beta = 0; beta < 4; ++beta) {
                d_xi[mu] -= Gamma[mu][alpha][beta] * k_arr[alpha] * xi_arr[beta];
                d_v[mu] -= Gamma[mu][alpha][beta] * k_arr[alpha] * V_arr[beta];
            }
        }
        for (int nu = 0; nu < 4; ++nu)
            for (int rho = 0; rho < 4; ++rho)
                for (int sigma = 0; sigma < 4; ++sigma)
                    d_v[mu] -= R[mu][nu][rho][sigma] * k_arr[nu] * xi_arr[rho] * k_arr[sigma];
    }
    dXi = Real4(d_xi[0], d_xi[1], d_xi[2], d_xi[3]);
    dV = Real4(d_v[0], d_v[1], d_v[2], d_v[3]);
}

// Advance both columns of a ray bundle across one accepted central-ray
// segment. Connection and curvature are evaluated at the start, cubic-Hermite
// midpoint, and end so this is the single device authority used by the trace
// kernel and the non-image analytic congruence probe.
public void AdvanceJacobiPairAcceptedSegment(
    Vec4Cart startPosition, Vec4Cart startTangent,
    Vec4Cart endPosition, Vec4Cart endTangent, Real dLambda,
    float M, float a, float Q, float Lambda,
    inout Real4 xiA, inout Real4 VA, inout Real4 xiB, inout Real4 VB) {
    Real R[4][4][4][4];
    Real Gamma[4][4][4];
    Real4 k1XiA, k1VA, k2XiA, k2VA, k3XiA, k3VA, k4XiA, k4VA;
    Real4 k1XiB, k1VB, k2XiB, k2VB, k3XiB, k3VB, k4XiB, k4VB;

    GetRiemannTensorCart(startPosition, M, a, Q, Lambda, R);
    GetKerrSchildFamilyChristoffel(startPosition, M, a, Q, Lambda, Gamma);
    GeodesicDeviationRhs(Gamma, R, startTangent, xiA, VA, k1XiA, k1VA);
    GeodesicDeviationRhs(Gamma, R, startTangent, xiB, VB, k1XiB, k1VB);

    Vec4Cart midpointPosition = AddVec4Cart(
        ScaleVec4Cart(AddVec4Cart(startPosition, endPosition), Real(0.5)),
        ScaleVec4Cart(
            AddVec4Cart(startTangent, ScaleVec4Cart(endTangent, Real(-1.0))),
            dLambda / Real(8.0)));
    Vec4Cart midpointTangent = AddVec4Cart(
        ScaleVec4Cart(
            AddVec4Cart(endPosition, ScaleVec4Cart(startPosition, Real(-1.0))),
            Real(1.5) / dLambda),
        ScaleVec4Cart(AddVec4Cart(startTangent, endTangent), Real(-0.25)));
    GetRiemannTensorCart(midpointPosition, M, a, Q, Lambda, R);
    GetKerrSchildFamilyChristoffel(midpointPosition, M, a, Q, Lambda, Gamma);

    GeodesicDeviationRhs(Gamma, R, midpointTangent,
                         xiA + k1XiA * (Real(0.5) * dLambda),
                         VA + k1VA * (Real(0.5) * dLambda), k2XiA, k2VA);
    GeodesicDeviationRhs(Gamma, R, midpointTangent,
                         xiA + k2XiA * (Real(0.5) * dLambda),
                         VA + k2VA * (Real(0.5) * dLambda), k3XiA, k3VA);
    GeodesicDeviationRhs(Gamma, R, midpointTangent,
                         xiB + k1XiB * (Real(0.5) * dLambda),
                         VB + k1VB * (Real(0.5) * dLambda), k2XiB, k2VB);
    GeodesicDeviationRhs(Gamma, R, midpointTangent,
                         xiB + k2XiB * (Real(0.5) * dLambda),
                         VB + k2VB * (Real(0.5) * dLambda), k3XiB, k3VB);

    GetRiemannTensorCart(endPosition, M, a, Q, Lambda, R);
    GetKerrSchildFamilyChristoffel(endPosition, M, a, Q, Lambda, Gamma);
    GeodesicDeviationRhs(Gamma, R, endTangent,
                         xiA + k3XiA * dLambda, VA + k3VA * dLambda, k4XiA, k4VA);
    GeodesicDeviationRhs(Gamma, R, endTangent,
                         xiB + k3XiB * dLambda, VB + k3VB * dLambda, k4XiB, k4VB);

    const Real sixthStep = dLambda / Real(6.0);
    xiA += (k1XiA + k2XiA * Real(2.0) + k3XiA * Real(2.0) + k4XiA) * sixthStep;
    VA += (k1VA + k2VA * Real(2.0) + k3VA * Real(2.0) + k4VA) * sixthStep;
    xiB += (k1XiB + k2XiB * Real(2.0) + k3XiB * Real(2.0) + k4XiB) * sixthStep;
    VB += (k1VB + k2VB * Real(2.0) + k3VB * Real(2.0) + k4VB) * sixthStep;
}

// Full-Riemann deviation acceleration D^2 xi^mu/dl^2 (DNGR eq A.19).
public Real4 GeodesicDeviationAccelFullRiemann(
    Vec4Cart x, Vec4Cart k, Real4 xi, float M, float a, float Q, float Lambda)
{
    Real R[4][4][4][4];
    GetRiemannTensorCart(x, M, a, Q, Lambda, R);
    return GeodesicDeviationAccelFromRiemann(R, k, xi);
}

// Deterministic tangent basis on the unit celestial sphere. The beam ellipse
// angle is measured in this basis and point-catalogue filtering consumes that
// same angle, so the least-aligned reference-axis decision must have one device
// authority. This mirrors core/celestial_tangent_basis.h.
public void MakeCelestialTangentBasis(float3 inputDirection,
                                      out float3 first, out float3 second) {
    float3 direction = normalize(inputDirection);
    float3 reference;
    if (abs(direction.x) <= abs(direction.y) && abs(direction.x) <= abs(direction.z)) {
        reference = float3(1.0f, 0.0f, 0.0f);
    } else if (abs(direction.y) <= abs(direction.z)) {
        reference = float3(0.0f, 1.0f, 0.0f);
    } else {
        reference = float3(0.0f, 0.0f, 1.0f);
    }
    first = normalize(reference - dot(reference, direction) * direction);
    second = cross(direction, first);
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
