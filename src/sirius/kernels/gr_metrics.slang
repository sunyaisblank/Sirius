// Metric-tensor family: Kerr-Schild (M, a, Q, Lambda), Ellis drainhole,
// Morris-Thorne wormhole, and the Alcubierre warp-drive family. Each routine
// returns the covariant metric g_mu_nu, its inverse g^mu_nu, and (where the
// .cu carried a closed form) the Christoffel symbols Gamma^lambda_mu_nu.
//
// The live Kerr-Schild branch uses the exact represented parameter and chart
// domain: scale-safe radial arithmetic replaces the retired implementation's
// epsilon radii and low-spin aliases. Deliberately NOT ported: the older
// spherical getAlcubierreMetric approximation (the
// WarpDrive Cartesian family below is the render dispatch's authority and
// supersedes it) and the Godel/Taub-NUT branches (docs/ARCHITECTURE.md
// section 2 deletes them).
//
// References cited at each routine follow docs/STYLE.md section 9.
module gr_metrics;

import gr_types;

// --- Coordinate conversions (RDOP002A.cu vec4SphToCart / vec4CartToSph) -------

// Y-polar spherical->Cartesian: x = r sin th cos ph, y = r cos th (polar
// axis), z = r sin th sin ph. This is the .cu convention (note y is the pole).
public Vec4Cart Vec4SphToCart(Vec4 sph) {
    float sin_theta = sin(sph.theta);
    float cos_theta = cos(sph.theta);
    float sin_phi = sin(sph.phi);
    float cos_phi = cos(sph.phi);
    return MakeVec4Cart(
        sph.t,
        sph.r * sin_theta * cos_phi,
        sph.r * sin_theta * sin_phi,
        sph.r * cos_theta);
}

public Vec4 Vec4CartToSph(Vec4Cart cart) {
    // The spherical chart is float (see gr_types); the casts collapse to
    // no-ops on the fp32 rung and narrow deliberately on fp64, where the
    // spherical consumers are the parity oracle, not the trajectory.
    Real r = sqrt(cart.x * cart.x + cart.y * cart.y + cart.z * cart.z);
    Real r_safe = max(r, Real(1e-10f));
    float theta = acos(float(cart.z / r_safe));
    float phi = atan2(float(cart.y), float(cart.x));
    return MakeVec4(float(cart.t), float(r), theta, phi);
}

// =============================================================================
// Unified Kerr-Schild family (RDOP002A.cu getKerrSchildFamily*).
// g_mu_nu = eta_mu_nu + H l_mu l_nu ; g^mu_nu = eta^mu_nu - H l^mu l^nu.
// Covers Minkowski/Schwarzschild/Kerr/Reissner-Nordstrom/Kerr-Newman and the
// de Sitter sector via (M, a, Q, Lambda). Cartesian chart, polynomial
// Christoffels, no pole singularity.
// =============================================================================

// Kerr radius from Cartesian position: solves r^4 - (R^2 - a^2) r^2 - a^2 z^2 = 0.
public Real ComputeKerrRadiusCart(Real x, Real y, Real z, Real a) {
    Real scale = max(max(abs(x), abs(y)), max(abs(z), abs(a)));
    if (scale == Real(0.0f)) return Real(0.0f);

    Real scaledX = x / scale;
    Real scaledY = y / scale;
    Real scaledZ = z / scale;
    Real scaledA = a / scale;
    Real scaledA2 = scaledA * scaledA;
    Real reduced = scaledX * scaledX + scaledY * scaledY + scaledZ * scaledZ - scaledA2;
    Real scaledCross = Real(2.0f) * scaledA * scaledZ;
    Real discriminantRoot = sqrt(reduced * reduced + scaledCross * scaledCross);
    Real scaledRadius = Real(0.0f);
    if (reduced >= Real(0.0f)) {
        scaledRadius = sqrt(Real(0.5f) * (reduced + discriminantRoot));
    } else if (scaledA != Real(0.0f) && scaledZ != Real(0.0f)) {
        scaledRadius = sqrt(Real(2.0f)) * abs(scaledA * scaledZ) /
                       sqrt(discriminantRoot - reduced);
    }
    return scale * scaledRadius;
}

// Numerical event domain shared by the live device integrator and the
// Kerr-Schild metric. Flat space is defined at every finite event; the pure
// de-Sitter Cartesian origin has an exact analytic limit. Every other r=0
// event is either a physical spherical singularity or the non-unique rotating
// sheet boundary and must be declined before a metric/connection evaluation.
public bool IsKerrSchildFamilyEventRepresented(Vec4Cart event, float M, float a,
                                               float Q, float Lambda) {
    if (!isfinite(event.t) || !isfinite(event.x) || !isfinite(event.y) ||
        !isfinite(event.z) || !isfinite(M) || !isfinite(a) ||
        !isfinite(Q) || !isfinite(Lambda)) {
        return false;
    }
    if (M == 0.0f && a == 0.0f && Q == 0.0f && Lambda == 0.0f) return true;
    if (M == 0.0f && a == 0.0f && Q == 0.0f &&
        event.x == Real(0.0f) && event.y == Real(0.0f) && event.z == Real(0.0f)) {
        return true;
    }
    Real r = ComputeKerrRadiusCart(event.x, event.y, event.z, Real(a));
    return isfinite(r) && r > Real(0.0f);
}

public void GetKerrSchildFamilyMetric(
    Vec4Cart pos, float M, float a, float Q, float Lambda,
    out Real g[4][4], out Real g_inv[4][4])
{
    Real x = pos.x;
    Real y = pos.y;
    Real z = pos.z;

    Real a2 = a * a;
    Real Q2 = Q * Q;

    [unroll] for (int i = 0; i < 4; i++)
        [unroll] for (int j = 0; j < 4; j++) { g[i][j] = 0.0f; g_inv[i][j] = 0.0f; }

    g[0][0] = -1.0f; g_inv[0][0] = -1.0f;
    g[1][1] = 1.0f;  g_inv[1][1] = 1.0f;
    g[2][2] = 1.0f;  g_inv[2][2] = 1.0f;
    g[3][3] = 1.0f;  g_inv[3][3] = 1.0f;

    if (M == 0.0f && a == 0.0f && Q == 0.0f && Lambda == 0.0f) {
        return;
    }

    Real r = ComputeKerrRadiusCart(x, y, z, Real(a));
    Real r2 = r * r;

    if (M == 0.0f && a == 0.0f && Q == 0.0f && r == Real(0.0f)) {
        return;
    }

    Real H = Real(0.0f);
    if (M != 0.0f || Q != 0.0f) {
        Real cosine = z / r;
        Real sigma = r2 + a2 * cosine * cosine;
        H = (2.0f * M * r - Q2) / sigma;
    }
    if (a == 0.0f && Lambda != 0.0f) {
        H += Lambda * r2 / 3.0f;
    }

    Real denom = r2 + a2;
    Real l[4] = {
        1.0f,
        (r * x + a * y) / denom,
        (r * y - a * x) / denom,
        z / r
    };

    [unroll] for (int mu = 0; mu < 4; mu++)
        [unroll] for (int nu = 0; nu < 4; nu++)
            g[mu][nu] += H * l[mu] * l[nu];

    Real l_up[4] = { 1.0f, -l[1], -l[2], -l[3] };
    [unroll] for (int mu = 0; mu < 4; mu++)
        [unroll] for (int nu = 0; nu < 4; nu++)
            g_inv[mu][nu] -= H * l_up[mu] * l_up[nu];

}

// Christoffel symbols Gamma^lambda_mu_nu = (1/2) g^{lambda sigma}
// (d_mu g_sigma_nu + d_nu g_sigma_mu - d_sigma g_mu_nu), all polynomial in
// (x, y, z, r). Transcribed from RDOP002A.cu getKerrSchildFamilyChristoffel.
public void GetKerrSchildFamilyChristoffel(
    Vec4Cart pos, float M, float a, float Q, float Lambda,
    out Real Gamma[4][4][4])
{
    [unroll] for (int i = 0; i < 4; i++)
        [unroll] for (int j = 0; j < 4; j++)
            [unroll] for (int k = 0; k < 4; k++)
                Gamma[i][j][k] = 0.0f;

    if (M == 0.0f && a == 0.0f && Q == 0.0f && Lambda == 0.0f) {
        return;
    }

    Real x = pos.x;
    Real y = pos.y;
    Real z = pos.z;

    Real a2 = a * a;
    Real Q2 = Q * Q;
    Real R2 = x * x + y * y + z * z;

    Real r = ComputeKerrRadiusCart(x, y, z, Real(a));
    if (M == 0.0f && a == 0.0f && Q == 0.0f && r == Real(0.0f)) {
        return;
    }
    Real r2 = r * r;

    Real Rm2 = R2 - a2;
    Real disc = Rm2 * Rm2 + 4.0f * a2 * z * z;
    Real sqrt_disc = sqrt(disc);

    Real d_disc_dx = 4.0f * x * Rm2;
    Real d_disc_dy = 4.0f * y * Rm2;
    Real d_disc_dz = 4.0f * z * (R2 + a2);

    Real d_r2_dx = x + d_disc_dx / (4.0f * sqrt_disc);
    Real d_r2_dy = y + d_disc_dy / (4.0f * sqrt_disc);
    Real d_r2_dz = z + d_disc_dz / (4.0f * sqrt_disc);

    Real dr[4];
    dr[0] = 0.0f;
    dr[1] = d_r2_dx / (2.0f * r);
    dr[2] = d_r2_dy / (2.0f * r);
    dr[3] = d_r2_dz / (2.0f * r);

    Real denom = r2 + a2;
    Real denom2 = denom * denom;
    Real l[4] = {
        1.0f,
        (r * x + a * y) / denom,
        (r * y - a * x) / denom,
        z / r
    };
    Real l_up[4] = { 1.0f, -l[1], -l[2], -l[3] };

    Real dl[4][4];
    [unroll] for (int i = 0; i < 4; i++)
        [unroll] for (int j = 0; j < 4; j++) dl[i][j] = 0.0f;

    Real f1 = r * x + a * y;
    Real f2 = denom;

    dl[1][1] = (r + x * dr[1]) / f2 - f1 * (2.0f * r * dr[1]) / denom2;
    dl[2][1] = (x * dr[2] + a) / f2 - f1 * (2.0f * r * dr[2]) / denom2;
    dl[3][1] = (x * dr[3]) / f2 - f1 * (2.0f * r * dr[3]) / denom2;

    Real g1 = r * y - a * x;

    dl[1][2] = (y * dr[1] - a) / f2 - g1 * (2.0f * r * dr[1]) / denom2;
    dl[2][2] = (r + y * dr[2]) / f2 - g1 * (2.0f * r * dr[2]) / denom2;
    dl[3][2] = (y * dr[3]) / f2 - g1 * (2.0f * r * dr[3]) / denom2;

    dl[1][3] = -z * dr[1] / r2;
    dl[2][3] = -z * dr[2] / r2;
    dl[3][3] = 1.0f / r - z * dr[3] / r2;

    Real cosine = z / r;
    Real sigma = r2 + a2 * cosine * cosine;
    Real H = Real(0.0f);
    Real asymptoticH = Real(0.0f);
    if (M != 0.0f || Q != 0.0f) {
        asymptoticH = (2.0f * M * r - Q2) / sigma;
        H = asymptoticH;
    }
    if (a == 0.0f && Lambda != 0.0f) {
        H += Lambda * r2 / 3.0f;
    }

    Real dH[4] = { 0.0f, 0.0f, 0.0f, 0.0f };
    for (int lam = 1; lam <= 3; lam++) {
        if (M != 0.0f || Q != 0.0f) {
            Real dCosine = (lam == 3 ? Real(1.0f) / r : Real(0.0f)) -
                           z * dr[lam] / r2;
            Real dSigma = 2.0f * r * dr[lam] + 2.0f * a2 * cosine * dCosine;
            Real dNumerator = 2.0f * M * dr[lam];
            dH[lam] = (dNumerator - asymptoticH * dSigma) / sigma;
        }
        if (a == 0.0f && Lambda != 0.0f) {
            dH[lam] += (2.0f * Lambda / 3.0f) * r * dr[lam];
        }
    }

    Real dg[4][4][4];
    for (int lam = 0; lam < 4; lam++)
        for (int mu = 0; mu < 4; mu++)
            for (int nu = 0; nu < 4; nu++)
                dg[lam][mu][nu] = dH[lam] * l[mu] * l[nu]
                                + H * dl[lam][mu] * l[nu]
                                + H * l[mu] * dl[lam][nu];

    Real eta_inv[4] = { -1.0f, 1.0f, 1.0f, 1.0f };

    for (int lam = 0; lam < 4; lam++) {
        for (int mu = 0; mu < 4; mu++) {
            for (int nu = mu; nu < 4; nu++) {
                Real sum = 0.0f;
                for (int s = 0; s < 4; s++) {
                    Real g_inv_lam_s = (lam == s) ? eta_inv[lam] : Real(0.0f);
                    g_inv_lam_s -= H * l_up[lam] * l_up[s];
                    Real bracket = dg[mu][s][nu] + dg[nu][s][mu] - dg[s][mu][nu];
                    sum += g_inv_lam_s * bracket;
                }
                Gamma[lam][mu][nu] = 0.5f * sum;
                Gamma[lam][nu][mu] = Gamma[lam][mu][nu];
            }
        }
    }

}

// =============================================================================
// Ellis drainhole traversable wormhole (RDOP002A.cu getEllisDrainholeMetric).
// ds^2 = -(1 - Fp^2) dt^2 + 2 Fp dt dr + dr^2 + Rp^2 (dth^2 + sin^2 th dph^2).
// m = mass parameter, n > m ensures a valid throat at alpha = sqrt(n^2 - m^2).
// =============================================================================
public void GetEllisDrainholeMetric(
    Vec4 pos, float m, float n,
    out float g[4][4], out float g_inv[4][4])
{
    float r = pos.r;
    if (n <= m) n = m + 0.1f;

    float alpha = sqrt(n * n - m * m);

    float dx = r - m;
    float Rp2 = dx * dx + alpha * alpha;
    float Rp = sqrt(Rp2);

    float psi = (n / alpha) * (1.5707963f - atan(dx / alpha));

    float exp_term = exp(-2.0f * m * psi / n);
    float Fp_sq = 1.0f - exp_term;
    float Fp = -sqrt(max(Fp_sq, 0.0f));

    float denom = 1.0f - Fp_sq;
    if (denom < 0.001f) denom = 0.001f;
    float Rp_areal = Rp / sqrt(denom);
    float Rp_areal2 = Rp_areal * Rp_areal;

    float theta = pos.theta;
    float sin_theta = sin(theta);
    float sin2 = sin_theta * sin_theta;

    [unroll] for (int i = 0; i < 4; i++)
        [unroll] for (int j = 0; j < 4; j++) { g[i][j] = 0.0f; g_inv[i][j] = 0.0f; }

    g[0][0] = -(1.0f - Fp_sq);
    g[0][1] = Fp;
    g[1][0] = Fp;
    g[1][1] = 1.0f;
    g[2][2] = Rp_areal2;
    g[3][3] = Rp_areal2 * sin2;

    float det_tr = g[0][0] * g[1][1] - g[0][1] * g[1][0];
    if (abs(det_tr) > 0.001f) {
        g_inv[0][0] = g[1][1] / det_tr;
        g_inv[0][1] = -g[0][1] / det_tr;
        g_inv[1][0] = -g[1][0] / det_tr;
        g_inv[1][1] = g[0][0] / det_tr;
    } else {
        g_inv[0][0] = -1.0f;
        g_inv[1][1] = 1.0f;
    }
    g_inv[2][2] = 1.0f / max(Rp_areal2, 0.01f);
    g_inv[3][3] = 1.0f / max(Rp_areal2 * sin2, 0.0001f);
}

// =============================================================================
// Exact isotropic Cartesian chart for the zero-tidal Ellis member rendered by
// Sirius, mirroring MorrisThorneCartesian.  With rho=|x|,
//
//   r_areal = rho + b0^2/(4 rho),
//   g_tt = -1,
//   g_ij = A^2 delta_ij,                 A = 1 + b0^2/(4 rho^2).
//
// The throat is regular at rho=b0/2.  Positive rho below it is the exact
// second-sheet continuation needed to keep every RK stage represented; the
// trace clips the accepted segment at the throat and never exposes that sheet
// as an output topology.  rho=0 is the other asymptotic infinity and is not a
// finite chart event.
// =============================================================================
public bool IsMorrisThorneCartesianEventRepresented(
    Vec4Cart pos, float b0, float Phi0, int shapeType)
{
    // Phi0=0 is part of the live asymptotic normalisation, not merely the
    // current caller default.  A finite nonzero constant could be removed by a
    // time-coordinate rescaling, but admitting it here without rescaling the
    // observer and boundary authorities would expose a different convention.
    if (shapeType != 0 || !isfinite(b0) || !(b0 > 0.0f) || Phi0 != 0.0f ||
        !isfinite(pos.t) || !isfinite(pos.x) || !isfinite(pos.y) || !isfinite(pos.z)) {
        return false;
    }
    Real rho2 = pos.x * pos.x + pos.y * pos.y + pos.z * pos.z;
    return isfinite(rho2) && rho2 > Real(0.0f);
}

public void GetMorrisThorneCartesianMetric(
    Vec4Cart pos, float b0, float Phi0, int shapeType,
    out Real g[4][4], out Real g_inv[4][4])
{
    [unroll] for (int i = 0; i < 4; i++)
        [unroll] for (int j = 0; j < 4; j++) { g[i][j] = 0.0f; g_inv[i][j] = 0.0f; }
    if (!IsMorrisThorneCartesianEventRepresented(pos, b0, Phi0, shapeType)) return;

    Real rho2 = pos.x * pos.x + pos.y * pos.y + pos.z * pos.z;
    Real halfB0OverRho = Real(0.5f * b0) / sqrt(rho2);
    Real q = halfB0OverRho * halfB0OverRho;
    Real conformalBase = Real(1.0f) + q;
    Real conformal = conformalBase * conformalBase;
    Real inverseConformal = Real(1.0f) / conformal;
    g[0][0] = -1.0f;
    g_inv[0][0] = -1.0f;
    for (int i = 0; i < 3; i++) {
        g[i + 1][i + 1] = conformal;
        g_inv[i + 1][i + 1] = inverseConformal;
    }
}

// For g_ij=A^2 delta_ij, define w_i=d_i ln A.  The only nonzero connection is
//   Gamma^l_ij = delta^l_j w_i + delta^l_i w_j - delta_ij w_l,
//   w_i = -2 q x_i / (A rho^2).
public void GetMorrisThorneCartesianChristoffel(
    Vec4Cart pos, float b0, float Phi0, int shapeType,
    out Real Gamma[4][4][4])
{
    [unroll] for (int i = 0; i < 4; i++)
        [unroll] for (int j = 0; j < 4; j++)
            [unroll] for (int k = 0; k < 4; k++)
                Gamma[i][j][k] = 0.0f;
    if (!IsMorrisThorneCartesianEventRepresented(pos, b0, Phi0, shapeType)) return;

    Real rho2 = pos.x * pos.x + pos.y * pos.y + pos.z * pos.z;
    Real halfB0OverRho = Real(0.5f * b0) / sqrt(rho2);
    Real q = halfB0OverRho * halfB0OverRho;
    Real conformalBase = Real(1.0f) + q;
    Real coordinates[3] = { pos.x, pos.y, pos.z };
    Real w[3];
    for (int axis = 0; axis < 3; ++axis) {
        w[axis] = Real(-2.0f) * q * coordinates[axis] / (conformalBase * rho2);
    }

    for (int l = 0; l < 3; l++) {
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                Real value = ((l == j) ? w[i] : Real(0.0f))
                           + ((l == i) ? w[j] : Real(0.0f))
                           - ((i == j) ? w[l] : Real(0.0f));
                Gamma[l + 1][i + 1][j + 1] = value;
            }
        }
    }
}

// =============================================================================
// Alcubierre / warp-drive family (RDOP002A.cu getWarpDrive*).
// ds^2 = -dt^2 + (dx - vs f(rs) dt)^2 + dy^2 + dz^2, Cartesian chart.
// f(rs) = [tanh(sigma(rs+R)) - tanh(sigma(rs-R))] / [2 tanh(sigma R)].
// =============================================================================
public bool IsWarpDriveEventRepresented(
    Vec4Cart pos, float vs, float sigma, float R, float xs, float ys, float zs)
{
    if (!isfinite(pos.t) || !isfinite(pos.x) || !isfinite(pos.y) || !isfinite(pos.z) ||
        !isfinite(vs) || abs(vs) > 10.0f || !isfinite(sigma) || !(sigma > 0.0f) ||
        sigma > 1000.0f || !isfinite(R) || !(R > 0.0f) || R > 1000.0f ||
        !isfinite(xs) || !isfinite(ys) || !isfinite(zs)) {
        return false;
    }
    float sigmaR = sigma * R;
    return sigmaR >= 0.1f && sigmaR <= 100.0f;
}

public float AlcubierreShapeFunction(float rs, float sigma, float R) {
    float tanh_plus = tanh(sigma * (rs + R));
    float tanh_minus = tanh(sigma * (rs - R));
    float tanh_R = tanh(sigma * R);
    return (tanh_plus - tanh_minus) / max(2.0f * tanh_R, 1e-10f);
}

public float AlcubierreShapeDerivative(float rs, float sigma, float R) {
    float cosh_plus = cosh(sigma * (rs + R));
    float cosh_minus = cosh(sigma * (rs - R));
    float sech2_plus = 1.0f / (cosh_plus * cosh_plus);
    float sech2_minus = 1.0f / (cosh_minus * cosh_minus);
    float tanh_R = tanh(sigma * R);
    return sigma * (sech2_plus - sech2_minus) / max(2.0f * tanh_R, 1e-10f);
}

public void GetWarpDriveMetric(
    Vec4Cart pos, float vs, float sigma, float R,
    float xs, float ys, float zs,
    out Real g[4][4], out Real g_inv[4][4])
{
    [unroll] for (int i = 0; i < 4; i++)
        [unroll] for (int j = 0; j < 4; j++) { g[i][j] = 0.0f; g_inv[i][j] = 0.0f; }
    if (!IsWarpDriveEventRepresented(pos, vs, sigma, R, xs, ys, zs)) return;

    Real x = pos.x;
    Real y = pos.y;
    Real z = pos.z;

    // The bubble centre follows xs(t)=xs0+vs*t, matching
    // core::WarpDriveFamily. Freezing xs at its t=0 value solves a different,
    // stationary spacetime.
    Real dx = x - (xs + vs * pos.t);
    Real dy = y - ys;
    Real dz = z - zs;
    Real rs = sqrt(dx * dx + dy * dy + dz * dz);
    rs = max(rs, Real(1e-10f));

    // The tanh-shaped bubble wall evaluates in float: GLSL.std.450 restricts
    // tanh/cosh to 32-bit floats, and the shape function is a bounded profile
    // whose float error does not compound along the trajectory.
    float f = AlcubierreShapeFunction(float(rs), sigma, R);
    Real vsf = vs * f;
    Real vsf2 = vsf * vsf;

    g[0][0] = -(1.0f - vsf2);
    g[0][1] = -vsf;
    g[1][0] = -vsf;
    g[1][1] = 1.0f;
    g[2][2] = 1.0f;
    g[3][3] = 1.0f;

    // The represented t-x block has exact determinant -1.
    g_inv[0][0] = -1.0f;
    g_inv[0][1] = -vsf;
    g_inv[1][0] = -vsf;
    g_inv[1][1] = 1.0f - vsf2;
    g_inv[2][2] = 1.0f;
    g_inv[3][3] = 1.0f;
}

public void GetWarpDriveChristoffel(
    Vec4Cart pos, float vs, float sigma, float R,
    float xs, float ys, float zs,
    out Real Gamma[4][4][4])
{
    [unroll] for (int i = 0; i < 4; i++)
        [unroll] for (int j = 0; j < 4; j++)
            [unroll] for (int k = 0; k < 4; k++)
                Gamma[i][j][k] = 0.0f;
    if (!IsWarpDriveEventRepresented(pos, vs, sigma, R, xs, ys, zs)) return;

    Real x = pos.x;
    Real y = pos.y;
    Real z = pos.z;

    Real dx = x - (xs + vs * pos.t);
    Real dy = y - ys;
    Real dz = z - zs;
    Real rs = sqrt(dx * dx + dy * dy + dz * dz);
    rs = max(rs, Real(1e-10f));

    // Float at the transcendental boundary, as in GetWarpDriveMetric above.
    float f = AlcubierreShapeFunction(float(rs), sigma, R);
    float df_drs = AlcubierreShapeDerivative(float(rs), sigma, R);

    Real drs[4] = { -vs * dx / rs, dx / rs, dy / rs, dz / rs };
    Real df[4] = {
        df_drs * drs[0], df_drs * drs[1], df_drs * drs[2], df_drs * drs[3]
    };

    Real metric[4][4];
    Real inverseMetric[4][4];
    GetWarpDriveMetric(pos, vs, sigma, R, xs, ys, zs, metric, inverseMetric);
    Real dg[4][4][4];
    [unroll] for (int derivative = 0; derivative < 4; ++derivative) {
        [unroll] for (int mu = 0; mu < 4; ++mu) {
            [unroll] for (int nu = 0; nu < 4; ++nu) {
                dg[derivative][mu][nu] = Real(0.0f);
            }
        }
        dg[derivative][0][0] = 2.0f * vs * vs * f * df[derivative];
        dg[derivative][0][1] = -vs * df[derivative];
        dg[derivative][1][0] = dg[derivative][0][1];
    }

    [unroll] for (int alpha = 0; alpha < 4; ++alpha) {
        [unroll] for (int mu = 0; mu < 4; ++mu) {
            [unroll] for (int nu = 0; nu < 4; ++nu) {
                Real sum = Real(0.0f);
                [unroll] for (int sigmaIndex = 0; sigmaIndex < 4; ++sigmaIndex) {
                    sum += inverseMetric[alpha][sigmaIndex] *
                           (dg[mu][sigmaIndex][nu] + dg[nu][sigmaIndex][mu] -
                            dg[sigmaIndex][mu][nu]);
                }
                Gamma[alpha][mu][nu] = Real(0.5f) * sum;
            }
        }
    }
}
