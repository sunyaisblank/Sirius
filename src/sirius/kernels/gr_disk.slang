// Accretion-disk emission physics: the Novikov-Thorne / Page-Thorne flux
// function Q(r), the Kerr ISCO, the branch-free Chebyshev blackbody colour fit,
// and the invariant redshift g-factor. Every coefficient (the eight-term
// Chebyshev fits, the log-domain mapping, the XYZ->sRGB matrix, the clamp
// ranges) is copied verbatim from RDOP002A.cu so fp32 output matches the legacy
// device functions.
//
// References: Page & Thorne (1974) ApJ 191, 499 (flux Q); Bardeen et al.
// (1972) (ISCO); the invariant g = (k.u_obs)/(k.u_emit).
module gr_disk;

import gr_types;
import gr_metrics;

// --- Chebyshev blackbody colour fit (RDOP002A.cu BB_* / blackbodyColor) -------

// log(T) domain mapping: T in [1000 K, 40000 K].
static const float BB_LOG_MIN = 6.907755279f;    // log(1000)
static const float BB_LOG_RANGE = 3.688879454f;  // log(40000) - log(1000)

static const float BB_X_CHEB[8] = {
    0.332461f, -0.098234f, 0.024891f, -0.008123f,
    0.003012f, -0.001156f, 0.000423f, -0.000148f
};
static const float BB_Y_CHEB[8] = {
    0.341231f, -0.056789f, 0.018234f, -0.006512f,
    0.002345f, -0.000891f, 0.000312f, -0.000108f
};

// Clenshaw recurrence for sum c_i T_i(u) (RDOP002A.cu chebyshevEval).
public float ChebyshevEval(float c[8], int n, float u) {
    float b0 = 0.0f;
    float b1 = 0.0f;
    for (int i = n - 1; i >= 0; i--) {
        float b2 = b1;
        b1 = b0;
        b0 = 2.0f * u * b1 - b2 + c[i];
    }
    return b0 - u * b1;
}

// T -> linear sRGB via CIE xy chromaticity (RDOP002A.cu blackbodyColor).
public float3 BlackbodyColor(float T) {
    T = clamp(T, 1000.0f, 40000.0f);

    float logT = log(T);
    float u = 2.0f * (logT - BB_LOG_MIN) / BB_LOG_RANGE - 1.0f;

    float x = ChebyshevEval(BB_X_CHEB, 8, u);
    float y = ChebyshevEval(BB_Y_CHEB, 8, u);

    x = clamp(x, 0.1f, 0.65f);
    y = clamp(y, 0.1f, 0.65f);

    float Y = 1.0f;
    float y_inv = 1.0f / max(y, 0.01f);
    float X = Y * y_inv * x;
    float Z = Y * y_inv * (1.0f - x - y);

    float3 rgb;
    rgb.x =  3.2404542f * X - 1.5371385f * Y - 0.4985314f * Z;
    rgb.y = -0.9692660f * X + 1.8760108f * Y + 0.0415560f * Z;
    rgb.z =  0.0556434f * X - 0.2040259f * Y + 1.0572252f * Z;

    rgb.x = max(rgb.x, 0.0f);
    rgb.y = max(rgb.y, 0.0f);
    rgb.z = max(rgb.z, 0.0f);
    return rgb;
}

// --- Kerr ISCO (RDOP002A.cu computeKerrISCO), geometric units, prograde -------
public float ComputeKerrISCO(float M, float a) {
    float spin = clamp(abs(a), 0.0f, 0.998f);
    if (spin < 0.001f) {
        return 6.0f * M;
    }
    float a2 = spin * spin;
    float oneMa2_13 = pow(1.0f - a2, 1.0f / 3.0f);
    float onePa_13 = pow(1.0f + spin, 1.0f / 3.0f);
    float oneMa_13 = pow(1.0f - spin, 1.0f / 3.0f);
    float Z1 = 1.0f + oneMa2_13 * (onePa_13 + oneMa_13);
    float Z2 = sqrt(3.0f * a2 + Z1 * Z1);
    float inner = (3.0f - Z1) * (3.0f + Z1 + 2.0f * Z2);
    float r_isco = M * (3.0f + Z2 - sqrt(max(inner, 0.0f)));
    return r_isco;
}

// --- Novikov-Thorne flux Q(r) (RDOP002A.cu computeNovikovThorneQ) -------------
public float ComputeNovikovThorneQ(float r, float r_isco, float M, float a_star) {
    a_star = max(-0.998f, min(0.998f, a_star));

    float x = sqrt(r / M);
    float x_isco = sqrt(r_isco / M);

    float A = 1.0f - 2.0f / (x * x) + a_star / (x * x * x);
    float B = 1.0f - 3.0f / (x * x) + 2.0f * a_star / (x * x * x);
    float D_sq = B;
    float D = (D_sq > 0.0f) ? sqrt(D_sq) : 0.0f;

    float B_isco = 1.0f - 3.0f / (x_isco * x_isco) + 2.0f * a_star / (x_isco * x_isco * x_isco);
    float D_isco = (B_isco > 0.0f) ? sqrt(B_isco) : 0.0f;

    if (D <= 0.0f || D_isco <= 0.0f || x <= x_isco * 1.001f) {
        return 0.0f;
    }

    float x_ratio = x_isco / x;
    float base_factor = (1.0f - sqrt(x_ratio));
    float log_factor = log(max(x / x_isco, 1.001f));

    float frame_drag = 0.0f;
    if (abs(a_star) > 0.01f) {
        frame_drag = (3.0f * a_star / (2.0f * x * x * x)) * log_factor;
    }

    float redshift_factor = sqrt(A / B);
    float Q = base_factor * (1.0f - frame_drag) * min(redshift_factor, 2.0f);
    Q *= D / D_isco;
    return max(Q, 0.0f);
}

// Novikov-Thorne surface temperature: T = T_inner * Q^{1/4} (RDOP002A.cu
// sampleAccretionDiskVolumetric temperature branch).
public float NovikovThorneTemperature(float innerTemperature, float Q) {
    return innerTemperature * pow(max(Q, 0.0f), 0.25f);
}

// --- Invariant redshift g-factor (RDOP002A.cu computeRelativisticG) -----------
// g = (k.u_obs)/(k.u_emit) with a static observer at infinity u_obs=(1,0,0,0).
// The metric at the emission point is supplied by the caller.
public float ComputeRelativisticGFromMetric(float g[4][4], Vec4 k, Vec4 u_emit) {
    Vec4 u_obs = MakeVec4(1.0f, 0.0f, 0.0f, 0.0f);
    float k_dot_u_obs = Dot4D(g, k, u_obs);
    float k_dot_u_emit = Dot4D(g, k, u_emit);
    float g_factor = k_dot_u_obs / min(k_dot_u_emit, -0.001f);
    return clamp(abs(g_factor), 0.05f, 20.0f);
}

// Keplerian disk 4-velocity for an equatorial circular geodesic
// (RDOP002A.cu computeKeplerianVelocity). Metric supplied by the caller.
public Vec4 ComputeKeplerianVelocity(float r, float M, float a, float g[4][4]) {
    float r32 = pow(r, 1.5f);
    float Omega = sqrt(M) / (r32 + a * sqrt(M));
    float denom = g[0][0] + 2.0f * Omega * g[0][3] + Omega * Omega * g[3][3];
    float u_t = 1.0f / sqrt(max(-denom, 0.01f));
    float u_phi = Omega * u_t;
    return MakeVec4(u_t, 0.0f, 0.0f, u_phi);
}
