// Metric-tensor family: Kerr-Schild (M, a, Q, Lambda), Ellis drainhole,
// Morris-Thorne wormhole, and the Alcubierre warp-drive family. Each routine
// returns the covariant metric g_mu_nu, its inverse g^mu_nu, and (where the
// .cu carried a closed form) the Christoffel symbols Gamma^lambda_mu_nu.
//
// Every expression, guard constant, and operation order is copied verbatim
// from the retired src/Sirius.Render/Acceleration/OptiX/RDOP002A.cu so that
// fp32 output is bit-comparable with the legacy device functions. Deliberately
// NOT ported: the older spherical getAlcubierreMetric approximation (the
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
    Real a2 = a * a;
    Real R2 = x * x + y * y + z * z;
    if (abs(a) < Real(1e-10f)) {
        return sqrt(max(R2, Real(1e-20f)));
    }
    Real Rm2 = R2 - a2;
    Real disc = Rm2 * Rm2 + 4.0f * a2 * z * z;
    Real r2 = (Rm2 + sqrt(max(disc, Real(0.0f)))) / 2.0f;
    return sqrt(max(r2, Real(1e-20f)));
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

    if (M < 1e-15f && abs(Lambda) < 1e-15f) {
        return;
    }

    Real r = ComputeKerrRadiusCart(x, y, z, Real(a));
    Real r2 = r * r;
    Real r4 = r2 * r2;

    Real f_denom = r4 + a2 * z * z;
    Real H = (2.0f * M * r - Q2) * r2 / max(f_denom, Real(1e-20f));
    if (abs(a) < 1e-15f && abs(Lambda) > 1e-15f) {
        H += Lambda * r2 / 3.0f;
    }

    Real denom = r2 + a2;
    Real l[4] = {
        1.0f,
        (r * x + a * y) / denom,
        (r * y - a * x) / denom,
        z / max(r, Real(1e-10f))
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

    if (M < 1e-15f && abs(Lambda) < 1e-15f) {
        return;
    }

    Real x = pos.x;
    Real y = pos.y;
    Real z = pos.z;

    Real a2 = a * a;
    Real Q2 = Q * Q;
    Real R2 = x * x + y * y + z * z;

    Real r = ComputeKerrRadiusCart(x, y, z, Real(a));
    Real r2 = r * r;
    Real r3 = r2 * r;
    Real r4 = r2 * r2;

    Real Rm2 = R2 - a2;
    Real disc = Rm2 * Rm2 + 4.0f * a2 * z * z;
    Real sqrt_disc = sqrt(max(disc, Real(1e-20f)));

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
    Real r_safe = max(r, Real(1e-10f));

    Real l[4] = {
        1.0f,
        (r * x + a * y) / denom,
        (r * y - a * x) / denom,
        z / r_safe
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

    Real r2_safe = r_safe * r_safe;
    dl[1][3] = -z * dr[1] / r2_safe;
    dl[2][3] = -z * dr[2] / r2_safe;
    dl[3][3] = 1.0f / r_safe - z * dr[3] / r2_safe;

    Real f_denom = r4 + a2 * z * z;
    Real numerator = (2.0f * M * r - Q2) * r2;
    Real H = numerator / max(f_denom, Real(1e-20f));
    if (abs(a) < 1e-15f && abs(Lambda) > 1e-15f) {
        H += Lambda * r2 / 3.0f;
    }

    Real dH[4] = { 0.0f, 0.0f, 0.0f, 0.0f };
    for (int lam = 1; lam <= 3; lam++) {
        Real d_num = 2.0f * r * dr[lam] * (3.0f * M * r - Q2);
        Real d_f_denom = 4.0f * r3 * dr[lam];
        if (lam == 3) d_f_denom += 2.0f * a2 * z;
        dH[lam] = (d_num * f_denom - numerator * d_f_denom) / (f_denom * f_denom);
        if (abs(a) < 1e-15f && abs(Lambda) > 1e-15f) {
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
// Morris-Thorne wormhole family (RDOP002A.cu getMorrisThorne*).
// ds^2 = -e^{2 Phi} dt^2 + dr^2/(1 - b(r)/r) + r^2 (dth^2 + sin^2 th dph^2).
// shapeType 0 = Ellis b(r) = b0^2/r ; 1 = ZeroTidal ; 2 = AbsurdlyBenign.
// =============================================================================
public Real WormholeShapeFunction(Real r, float b0, int shapeType) {
    r = max(r, Real(b0));
    if (shapeType == 1) return b0;
    if (shapeType == 2) return b0 * (2.0f - b0 / r);
    return b0 * b0 / r;  // Ellis (default)
}

public Real WormholeShapeDerivative(Real r, float b0, int shapeType) {
    r = max(r, Real(b0));
    Real r2 = r * r;
    if (shapeType == 1) return 0.0f;
    if (shapeType == 2) return b0 * b0 / r2;
    return -b0 * b0 / r2;  // Ellis (default)
}

public void GetMorrisThorneMetric(
    Vec4 pos, float b0, float Phi0, int shapeType,
    out float g[4][4], out float g_inv[4][4])
{
    float r = max(pos.r, b0 * 1.001f);
    float theta = pos.theta;

    float sin_theta = sin(theta);
    float sin2 = sin_theta * sin_theta;
    float r2 = r * r;

    float b = float(WormholeShapeFunction(Real(r), b0, shapeType));
    float one_minus_b_over_r = 1.0f - b / r;
    one_minus_b_over_r = max(one_minus_b_over_r, 0.001f);

    float exp2Phi = exp(2.0f * Phi0);

    [unroll] for (int i = 0; i < 4; i++)
        [unroll] for (int j = 0; j < 4; j++) { g[i][j] = 0.0f; g_inv[i][j] = 0.0f; }

    g[0][0] = -exp2Phi;
    g[1][1] = 1.0f / one_minus_b_over_r;
    g[2][2] = r2;
    g[3][3] = r2 * sin2;

    g_inv[0][0] = -1.0f / exp2Phi;
    g_inv[1][1] = one_minus_b_over_r;
    g_inv[2][2] = 1.0f / r2;
    g_inv[3][3] = 1.0f / max(r2 * sin2, 1e-10f);
}

// =============================================================================
// Morris-Thorne Cartesian embedding, mirroring the CPU render authority
// (core/metrics/morris_thorne_family.h MorrisThorneCartesian). With r = |x|,
// n = x/r the spatial metric is flat plus rank one:
//   g_tt = -e^{2 Phi0},  g_ij = delta_ij + f n_i n_j,  f = b/(r - b),
//   g^ij = delta_ij - (f/(1+f)) n_i n_j          (Sherman-Morrison).
// One sheet, r clamped to b0 * 1.001 exactly as the CPU class clamps; the
// throat is the capture surface and trace.slang tests it on |x| directly.
// Rational arithmetic plus sqrt only, so the fp64 rung covers it unchanged.
// =============================================================================
public void GetMorrisThorneCartesianMetric(
    Vec4Cart pos, float b0, float Phi0, int shapeType,
    out Real g[4][4], out Real g_inv[4][4])
{
    Real x = pos.x;
    Real y = pos.y;
    Real z = pos.z;
    Real r = sqrt(x * x + y * y + z * z);
    Real r_min = b0 * 1.001f;
    Real n1, n2, n3;
    if (r < Real(1e-30f)) {
        n1 = 0.0f; n2 = 0.0f; n3 = 1.0f;
        r = r_min;
    } else {
        Real inv = 1.0f / r;
        n1 = x * inv; n2 = y * inv; n3 = z * inv;
        r = max(r, r_min);
    }

    Real b = WormholeShapeFunction(r, b0, shapeType);
    Real r_minus_b = max(r - b, Real(1e-10f));
    Real f = b / r_minus_b;
    Real c = f / (1.0f + f);
    float exp2Phi = exp(2.0f * Phi0);

    Real n[3] = { n1, n2, n3 };
    [unroll] for (int i = 0; i < 4; i++)
        [unroll] for (int j = 0; j < 4; j++) { g[i][j] = 0.0f; g_inv[i][j] = 0.0f; }

    g[0][0] = -exp2Phi;
    g_inv[0][0] = -1.0f / exp2Phi;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            Real dij = (i == j) ? Real(1.0f) : Real(0.0f);
            g[i + 1][j + 1] = dij + f * n[i] * n[j];
            g_inv[i + 1][j + 1] = dij - c * n[i] * n[j];
        }
    }
}

// Christoffels of the embedding, contracted from the analytic derivative
//   d_k g_ij = f' n_k n_i n_j + (f/r)(d_ki n_j + d_kj n_i - 2 n_k n_i n_j),
//   f' = (b' r - b)/(r - b)^2,
// and the closed-form inverse; the time row and column vanish (static,
// zero-tidal Phi0). Same source formulas as MorrisThorneCartesian::Evaluate.
public void GetMorrisThorneCartesianChristoffel(
    Vec4Cart pos, float b0, float Phi0, int shapeType,
    out Real Gamma[4][4][4])
{
    [unroll] for (int i = 0; i < 4; i++)
        [unroll] for (int j = 0; j < 4; j++)
            [unroll] for (int k = 0; k < 4; k++)
                Gamma[i][j][k] = 0.0f;

    Real x = pos.x;
    Real y = pos.y;
    Real z = pos.z;
    Real r = sqrt(x * x + y * y + z * z);
    Real r_min = b0 * 1.001f;
    Real n1, n2, n3;
    if (r < Real(1e-30f)) {
        n1 = 0.0f; n2 = 0.0f; n3 = 1.0f;
        r = r_min;
    } else {
        Real inv = 1.0f / r;
        n1 = x * inv; n2 = y * inv; n3 = z * inv;
        r = max(r, r_min);
    }

    Real b = WormholeShapeFunction(r, b0, shapeType);
    Real db_dr = WormholeShapeDerivative(r, b0, shapeType);
    Real r_minus_b = max(r - b, Real(1e-10f));
    Real f = b / r_minus_b;
    Real df_dr = (db_dr * r - b) / (r_minus_b * r_minus_b);
    Real c = f / (1.0f + f);
    Real f_over_r = f / r;

    Real n[3] = { n1, n2, n3 };

    // dg[k][i][j] = d_k g_ij, spatial block only.
    Real dg[3][3][3];
    for (int k = 0; k < 3; k++) {
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                Real dki = (k == i) ? Real(1.0f) : Real(0.0f);
                Real dkj = (k == j) ? Real(1.0f) : Real(0.0f);
                dg[k][i][j] = df_dr * n[k] * n[i] * n[j]
                            + f_over_r * (dki * n[j] + dkj * n[i]
                                          - 2.0f * n[k] * n[i] * n[j]);
            }
        }
    }

    // Gamma^l_ij = 1/2 g^lm (d_i g_mj + d_j g_mi - d_m g_ij), spatial indices.
    for (int l = 0; l < 3; l++) {
        for (int i = 0; i < 3; i++) {
            for (int j = i; j < 3; j++) {
                Real sum = 0.0f;
                for (int m = 0; m < 3; m++) {
                    Real g_inv_lm = ((l == m) ? Real(1.0f) : Real(0.0f)) - c * n[l] * n[m];
                    sum += g_inv_lm * (dg[i][m][j] + dg[j][m][i] - dg[m][i][j]);
                }
                Gamma[l + 1][i + 1][j + 1] = 0.5f * sum;
                Gamma[l + 1][j + 1][i + 1] = Gamma[l + 1][i + 1][j + 1];
            }
        }
    }
}

// safe_cot: regularised cotangent, finite at the poles (RDOP002A.cu safe_cot).
float SafeCot(float theta) {
    float sin_theta = sin(theta);
    float cos_theta = cos(theta);
    const float POLE_THRESHOLD = 0.0001f;
    if (abs(sin_theta) < POLE_THRESHOLD) {
        float dist_from_pole = (theta < PI * 0.5f) ? theta : (PI - theta);
        dist_from_pole = max(dist_from_pole, 1e-6f);
        float sign = (theta < PI * 0.5f) ? 1.0f : -1.0f;
        return sign / dist_from_pole;
    }
    return cos_theta / sin_theta;
}

public void GetMorrisThorneChristoffel(
    Vec4 pos, float b0, float Phi0, int shapeType,
    out float Gamma[4][4][4])
{
    [unroll] for (int i = 0; i < 4; i++)
        [unroll] for (int j = 0; j < 4; j++)
            [unroll] for (int k = 0; k < 4; k++)
                Gamma[i][j][k] = 0.0f;

    float r = max(pos.r, b0 * 1.001f);
    float theta = pos.theta;

    float sin_theta = sin(theta);
    float cos_theta = cos(theta);
    float sin2 = sin_theta * sin_theta;
    float r2 = r * r;

    float b = float(WormholeShapeFunction(Real(r), b0, shapeType));
    float db_dr = float(WormholeShapeDerivative(Real(r), b0, shapeType));

    float one_minus_b_over_r = 1.0f - b / r;
    one_minus_b_over_r = max(one_minus_b_over_r, 0.001f);

    float d_b_over_r = db_dr / r - b / r2;

    Gamma[1][1][1] = -d_b_over_r / (2.0f * one_minus_b_over_r);
    Gamma[1][2][2] = -(r - b);
    Gamma[1][3][3] = -(r - b) * sin2;
    Gamma[2][1][2] = 1.0f / r;
    Gamma[2][2][1] = 1.0f / r;
    Gamma[2][3][3] = -sin_theta * cos_theta;
    Gamma[3][1][3] = 1.0f / r;
    Gamma[3][3][1] = 1.0f / r;

    float cot_theta = SafeCot(theta);
    Gamma[3][2][3] = cot_theta;
    Gamma[3][3][2] = cot_theta;
    // Phi0 constant -> no time-related Christoffels (matches the .cu).
}

// =============================================================================
// Alcubierre / warp-drive family (RDOP002A.cu getWarpDrive*).
// ds^2 = -dt^2 + (dx - vs f(rs) dt)^2 + dy^2 + dz^2, Cartesian chart.
// f(rs) = [tanh(sigma(rs+R)) - tanh(sigma(rs-R))] / [2 tanh(sigma R)].
// =============================================================================
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

    [unroll] for (int i = 0; i < 4; i++)
        [unroll] for (int j = 0; j < 4; j++) { g[i][j] = 0.0f; g_inv[i][j] = 0.0f; }

    g[0][0] = -(1.0f - vsf2);
    g[0][1] = -vsf;
    g[1][0] = -vsf;
    g[1][1] = 1.0f;
    g[2][2] = 1.0f;
    g[3][3] = 1.0f;

    Real det_tx = g[0][0] * g[1][1] - g[0][1] * g[1][0];
    det_tx = max(abs(det_tx), Real(1e-10f)) * (det_tx >= 0.0f ? Real(1.0f) : Real(-1.0f));

    g_inv[0][0] = g[1][1] / det_tx;
    g_inv[0][1] = -g[0][1] / det_tx;
    g_inv[1][0] = -g[1][0] / det_tx;
    g_inv[1][1] = g[0][0] / det_tx;
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
