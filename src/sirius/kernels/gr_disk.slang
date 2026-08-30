// Accretion-disk emission physics: the Novikov-Thorne / Page-Thorne flux
// function Q(r), the Kerr ISCO, Planck/CIE blackbody colour, and the invariant
// redshift g-factor.
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

// --- Planck spectrum through the Wyman-Sloan-Shirley CIE 1931 fit -----------

float CieGaussian(float wavelengthNm, float centre, float leftWidth,
                  float rightWidth, float amplitude) {
    float width = wavelengthNm < centre ? leftWidth : rightWidth;
    float offset = (wavelengthNm - centre) * width;
    return amplitude * exp(-0.5f * offset * offset);
}

float3 CieXyz(float wavelengthNm) {
    float X = CieGaussian(wavelengthNm, 442.0f, 0.0624f, 0.0374f, 0.362f) +
              CieGaussian(wavelengthNm, 599.8f, 0.0264f, 0.0323f, 1.056f) -
              CieGaussian(wavelengthNm, 501.1f, 0.0490f, 0.0382f, 0.065f);
    float Y = CieGaussian(wavelengthNm, 568.8f, 0.0213f, 0.0247f, 0.821f) +
              CieGaussian(wavelengthNm, 530.9f, 0.0613f, 0.0322f, 0.286f);
    float Z = CieGaussian(wavelengthNm, 437.0f, 0.0845f, 0.0278f, 1.217f) +
              CieGaussian(wavelengthNm, 459.0f, 0.0385f, 0.0725f, 0.681f);
    return float3(X, Y, Z);
}

// The wavelength-independent Planck prefactor and integration interval cancel
// during chromaticity normalisation. Wavelength is in nm and c2 in nm K.
public float3 BlackbodyColor(float T) {
    if (T <= 0.0f) return float3(0.0f);
    float3 xyz = float3(0.0f);
    static const int kSamples = 32;
    static const float kMinimumNm = 380.0f;
    static const float kIntervalNm = 400.0f / float(kSamples);
    static const float kPlanckC2NmK = 1.4387768775e7f;
    for (int sample = 0; sample < kSamples; ++sample) {
        float wavelengthNm = kMinimumNm + (float(sample) + 0.5f) * kIntervalNm;
        float exponent = kPlanckC2NmK / (wavelengthNm * T);
        float relativeRadiance =
            1.0f / (pow(wavelengthNm, 5.0f) * (exp(exponent) - 1.0f));
        xyz += CieXyz(wavelengthNm) * relativeRadiance;
    }

    float3 rgb;
    rgb.x =  3.2404542f * xyz.x - 1.5371385f * xyz.y - 0.4985314f * xyz.z;
    rgb.y = -0.9692660f * xyz.x + 1.8760108f * xyz.y + 0.0415560f * xyz.z;
    rgb.z =  0.0556434f * xyz.x - 0.2040259f * xyz.y + 1.0572252f * xyz.z;
    rgb = max(rgb, float3(0.0f));
    return rgb / max(max(rgb.x, rgb.y), max(rgb.z, 1.0e-20f));
}

// --- Kerr ISCO (Bardeen et al. 1972), signed dimensionless spin ---------------
public float ComputeKerrISCO(float M, float aStar) {
    // The typed session owns the represented |a/M| <= 0.998 live domain.
    // Evaluate every admitted non-zero spin instead of substituting the
    // Schwarzschild radius throughout an arbitrary low-spin interval.
    float spin = abs(aStar);
    if (spin == 0.0f) {
        return 6.0f * M;
    }
    float a2 = spin * spin;
    float oneMa2_13 = pow(1.0f - a2, 1.0f / 3.0f);
    float onePa_13 = pow(1.0f + spin, 1.0f / 3.0f);
    float oneMa_13 = pow(1.0f - spin, 1.0f / 3.0f);
    float Z1 = 1.0f + oneMa2_13 * (onePa_13 + oneMa_13);
    float Z2 = sqrt(3.0f * a2 + Z1 * Z1);
    float inner = (3.0f - Z1) * (3.0f + Z1 + 2.0f * Z2);
    float branch = aStar >= 0.0f ? -1.0f : 1.0f;
    float r_isco = M * (3.0f + Z2 + branch * sqrt(max(inner, 0.0f)));
    return r_isco;
}

// --- Full Page-Thorne profile ------------------------------------------------
// Independent fp32 port of core::AccretionDiskD::FullPageThorneFlux. The
// physical prefactor cancels when the operator temperature scale is normalised
// at 1.5 r_ISCO, leaving q(r)/r because the Kerr equatorial determinant factor
// e^(nu+psi+mu) is r.
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
static const float DISK_TEMPERATURE_REFERENCE_RADIUS_RATIO = 1.5f;

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
    return max(correction / r, 0.0f);
}

public float FullPageThorneTemperature(float temperatureScale, float r, float rIsco,
                                       float M, float aStar) {
    float radiusInM = r / M;
    float iscoInM = rIsco / M;
    float flux = FullPageThorneFluxShape(radiusInM, iscoInM, aStar);
    float reference = FullPageThorneFluxShape(
        DISK_TEMPERATURE_REFERENCE_RADIUS_RATIO * iscoInM, iscoInM, aStar);
    if (flux <= 0.0f || reference <= 0.0f) return 0.0f;
    return temperatureScale * pow(flux / reference, 0.25f);
}

// Explicit Newtonian zero-torque standard thin-disk substitute. This is the
// Shakura-Sunyaev radial effective-temperature closure only; it does not claim
// relativistic Page-Thorne flux or alpha-disk vertical structure. As on the
// CPU, temperatureScale is the effective temperature at 1.5 rInner.
float ShakuraSunyaevFluxShape(float r, float rInner) {
    if (r <= rInner || rInner <= 0.0f) return 0.0f;
    float ratio = rInner / r;
    return ratio * ratio * ratio * max(1.0f - sqrt(ratio), 0.0f);
}

public float ShakuraSunyaevTemperature(float temperatureScale, float r,
                                       float rInner) {
    if (temperatureScale <= 0.0f) return 0.0f;
    float flux = ShakuraSunyaevFluxShape(r, rInner);
    float reference = ShakuraSunyaevFluxShape(
        DISK_TEMPERATURE_REFERENCE_RADIUS_RATIO * rInner, rInner);
    if (flux <= 0.0f || reference <= 0.0f) return 0.0f;
    return temperatureScale * pow(flux / reference, 0.25f);
}

// Exact stationary/axisymmetric transfer shared with
// core/relativistic_transfer.h.  The renderer integrates in Kerr-Schild
// Cartesian coordinates, so the photon covector is formed there.  The Killing
// quantities E=-p_t and L_z=p_phi=x p_y-y p_x are unchanged by the
// Boyer-Lindquist/Kerr-Schild chart transformation.  For every circular frame
// u=u^t(partial_t+Omega partial_phi), nu=-p.u=u^t(E-Omega L_z).
//
// x/y are reconstructed on the equatorial surface with
// rho=sqrt(r^2+a^2).  This also defines the explicit midplane projection used
// by the current axisymmetric volumetric disk model.  The return value is
// (full circular-emitter g, ZAMO g); zero means the requested frame or photon
// frequency lies outside the physical domain.  No display-oriented clamp is
// permitted in this frequency authority.
public float2 ComputeKerrDiskTransferCart(float r, float phi, Vec4Cart photonVelocity,
                                          float observerFrequency, float M, float a) {
    if (M <= 0.0f || r <= 0.0f || observerFrequency <= 0.0f || abs(a) > M) {
        return float2(0.0f);
    }

    float rho = sqrt(r * r + a * a);
    Vec4Cart position =
        MakeVec4Cart(Real(0.0f), Real(rho * cos(phi)), Real(rho * sin(phi)), Real(0.0f));
    Real metric[4][4];
    Real metricInverse[4][4];
    GetKerrSchildFamilyMetric(position, M, a, 0.0f, 0.0f, metric, metricInverse);

    // The live tracer advances the past-directed camera-to-source tangent.
    // Transfer uses the corresponding physical future photon, -k_past.
    Real k[4] = {
        -photonVelocity.t, -photonVelocity.x, -photonVelocity.y, -photonVelocity.z
    };
    Real p[4] = { Real(0.0f), Real(0.0f), Real(0.0f), Real(0.0f) };
    [unroll] for (int mu = 0; mu < 4; ++mu) {
        [unroll] for (int nu = 0; nu < 4; ++nu) {
            p[mu] += metric[mu][nu] * k[nu];
        }
    }

    float energy = -float(p[0]);
    float angularMomentum =
        -float(position.y) * float(p[1]) + float(position.x) * float(p[2]);
    if (energy <= 0.0f) return float2(0.0f);

    float gtt = -(1.0f - 2.0f * M / r);
    float gtphi = -2.0f * M * a / r;
    float gphiphi = r * r + a * a + 2.0f * M * a * a / r;
    float sqrtM = sqrt(M);
    float emitterOmega = sqrtM / (pow(r, 1.5f) + a * sqrtM);
    float zamoOmega = -gtphi / gphiphi;
    float emitterNorm =
        gtt + 2.0f * emitterOmega * gtphi + emitterOmega * emitterOmega * gphiphi;
    float zamoNorm = gtt + 2.0f * zamoOmega * gtphi + zamoOmega * zamoOmega * gphiphi;
    if (emitterNorm >= 0.0f || zamoNorm >= 0.0f) return float2(0.0f);

    float emitterFrequency =
        (energy - emitterOmega * angularMomentum) / sqrt(-emitterNorm);
    float zamoFrequency = (energy - zamoOmega * angularMomentum) / sqrt(-zamoNorm);
    if (emitterFrequency <= 0.0f || zamoFrequency <= 0.0f) return float2(0.0f);
    return float2(observerFrequency / emitterFrequency,
                  observerFrequency / zamoFrequency);
}
