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

// Bolometric Liouville scaling for disk emission. Both thin and volumetric
// render paths call this single kernel authority so g^4 is applied once.
public float ObservedDiskBolometricSource(float emittedSource, float redshift) {
    float g2 = redshift * redshift;
    return emittedSource * g2 * g2;
}

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

// --- Full Page-Thorne profile ------------------------------------------------
// Independent fp32 port of core::AccretionDiskD::FullPageThorneFlux. The
// physical prefactor cancels when the operator temperature scale is normalised
// at 1.5 r_ISCO, leaving Q(r)/r^3.
static const float PT_GL_NODES[16] = {
    -0.989400935f, -0.944575024f, -0.865631202f, -0.755404408f,
    -0.617876244f, -0.458016778f, -0.281603551f, -0.095012510f,
     0.095012510f,  0.281603551f,  0.458016778f,  0.617876244f,
     0.755404408f,  0.865631202f,  0.944575024f,  0.989400935f
};
static const float PT_GL_WEIGHTS[16] = {
    0.027152459f, 0.062253524f, 0.095158512f, 0.124628971f,
    0.149595989f, 0.169156520f, 0.182603515f, 0.189450611f,
    0.189450611f, 0.182603515f, 0.169156520f, 0.149595989f,
    0.124628971f, 0.095158512f, 0.062253524f, 0.027152459f
};

float PtSpecificEnergy(float r, float a) {
    float sqrtR = sqrt(r);
    float root = sqrt(max(r * r - 3.0f * r + 2.0f * a * sqrtR, 1.0e-20f));
    return (r * r - 2.0f * r + a * sqrtR) / (r * root);
}

float PtSpecificAngularMomentum(float r, float a) {
    float sqrtR = sqrt(r);
    float root = sqrt(max(r * r - 3.0f * r + 2.0f * a * sqrtR, 1.0e-20f));
    return (r * r - 2.0f * a * sqrtR + a * a) / (sqrtR * root);
}

float PtAngularVelocity(float r, float a) {
    return 1.0f / (pow(r, 1.5f) + a);
}

float PtAngularMomentumDerivative(float r, float a) {
    float sqrtR = sqrt(r);
    float N = r * r - 2.0f * a * sqrtR + a * a;
    float radicand = r * r - 3.0f * r + 2.0f * a * sqrtR;
    if (radicand <= 1.0e-12f) return 0.0f;
    float aOverSqrtR = a / sqrtR;
    float dN = 2.0f * r - aOverSqrtR;
    float dRadicand = 2.0f * r - 3.0f + aOverSqrtR;
    float denominator = sqrtR * sqrt(radicand);
    float dLogDenominator = 0.5f / r + 0.5f * dRadicand / radicand;
    return (dN - N * dLogDenominator) / denominator;
}

float PtAngularVelocityDerivative(float r, float a) {
    float denominator = pow(r, 1.5f) + a;
    return -3.0f * sqrt(r) / (2.0f * denominator * denominator);
}

public float FullPageThorneFluxShape(float r, float rIsco, float aStar) {
    if (r <= rIsco * 1.001f) return 0.0f;
    float energy = PtSpecificEnergy(r, aStar);
    float angularMomentum = PtSpecificAngularMomentum(r, aStar);
    float omega = PtAngularVelocity(r, aStar);
    float halfWidth = 0.5f * (r - rIsco);
    float midpoint = 0.5f * (r + rIsco);
    float integral = 0.0f;
    for (int i = 0; i < 16; ++i) {
        float ri = midpoint + halfWidth * PT_GL_NODES[i];
        float ei = PtSpecificEnergy(ri, aStar);
        float li = PtSpecificAngularMomentum(ri, aStar);
        float oi = PtAngularVelocity(ri, aStar);
        integral += PT_GL_WEIGHTS[i] * (ei - oi * li) *
                    PtAngularMomentumDerivative(ri, aStar);
    }
    integral *= halfWidth;
    float invariant = energy - omega * angularMomentum;
    if (abs(invariant) < 1.0e-12f) return 0.0f;
    float correction =
        -PtAngularVelocityDerivative(r, aStar) * integral / (invariant * invariant);
    return max(correction / (r * r * r), 0.0f);
}

public float FullPageThorneTemperature(float temperatureScale, float r, float rIsco,
                                       float M, float aStar) {
    float radiusInM = r / M;
    float iscoInM = rIsco / M;
    float flux = FullPageThorneFluxShape(radiusInM, iscoInM, aStar);
    float reference = FullPageThorneFluxShape(1.5f * iscoInM, iscoInM, aStar);
    if (flux <= 0.0f || reference <= 0.0f) return 0.0f;
    return temperatureScale * pow(flux / reference, 0.25f);
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

// Redshift g-factor for a thin-disk crossing evaluated in the Kerr-Schild
// Cartesian render chart. This is the closed-form gravitational x Doppler factor
// the CPU tracer applies (GeodesicTracer::ComputeGFactor): the render integrator
// runs in Cartesian, where the BL-basis ComputeKeplerianVelocity above does not
// apply, so the disk redshift is authored here in the same analytic form both
// backends share. `rayVel` is the photon spatial 4-velocity at the crossing
// (Cartesian), `phi` the crossing azimuth, `a` the absolute spin a = a*/M.
//   g = grav / (gamma (1 - v_orb (A))),  A = (-n_x) sin phi - (-n_y) cos phi,
//   grav = sqrt(1 - 2 M r / (r^2 + a^2)),  v_orb = r Omega,
//   Omega = sqrt(M) / (r^{3/2} + a sqrt(M)).   (Bardeen et al. 1972; the
// clamps guard the physical singularities exactly as the CPU path does.)
public float ComputeDiskGFactorCart(float r, float phi, float3 rayVel, float M, float a) {
    float r32 = pow(r, 1.5f);
    float sqrtM = sqrt(M);
    float Omega = sqrtM / (r32 + a * sqrtM);
    float v_orb = r * Omega;

    float a2 = a * a;
    float r2 = r * r;
    float grav_factor = sqrt(max(1.0f - 2.0f * M * r / (r2 + a2), 0.0f));

    float v_mag = length(rayVel);
    if (v_mag < 1e-10f) {
        return grav_factor;
    }
    float n_x = rayVel.x / v_mag;
    float n_y = rayVel.y / v_mag;

    float gamma = 1.0f / sqrt(max(1e-10f, 1.0f - v_orb * v_orb));

    float sin_phi = sin(phi);
    float cos_phi = cos(phi);
    // n points toward the observer; negate for the outgoing-photon convention.
    float A = (-n_x) * sin_phi - (-n_y) * cos_phi;

    float doppler_denom = clamp(gamma * (1.0f - v_orb * A), 0.1f, 10.0f);
    return clamp(grav_factor / doppler_denom, 0.1f, 5.0f);
}

// Gravitational-only branch used by the artistic Doppler-suppression toggle.
// This is the same grav term as ComputeDiskGFactorCart with orbital velocity
// set to zero.
public float ComputeDiskGravitationalFactorCart(float r, float M, float a) {
    float r2 = r * r;
    float grav = sqrt(max(1.0f - 2.0f * M * r / (r2 + a * a), 0.0f));
    return clamp(grav, 0.1f, 5.0f);
}
