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

// Stable 1-exp(-deltaTau) for the complete represented non-negative fp32
// optical-depth domain.  The series branch preserves optically thin layers
// that direct subtraction rounds to zero.
public float GreyLayerAbsorbedFraction(float deltaTau) {
    if (!isfinite(deltaTau) || deltaTau < 0.0f) return 0.0f;
    if (deltaTau <= 1.0e-2f) {
        return deltaTau * (1.0f - deltaTau * (0.5f - deltaTau / 6.0f));
    }
    return 1.0f - exp(-deltaTau);
}

// Finite-support vertical opacity for the declared phenomenological volume.
// integral[-3H,+3H] exp(-z^2/(2H^2)) dz
//   = sqrt(2 pi) erf(3/sqrt(2)) H
//   = 2.4998608894830947 H.
// The rounded fp32 constant has the same value used by the CPU closure after
// conversion to float, so tauMidplane is the represented truncated column.
static const float VOLUME_VERTICAL_TRUNCATION_SIGMA = 3.0f;
static const float TRUNCATED_GAUSSIAN_COLUMN_INTEGRAL = 2.499861f;

public float TruncatedGaussianOpacityDensity(float tauRef, float radius,
                                             float referenceRadius, float scaleHeight,
                                             float height) {
    if (!isfinite(tauRef) || tauRef < 0.0f || !isfinite(radius) || radius <= 0.0f ||
        !isfinite(referenceRadius) || referenceRadius <= 0.0f ||
        !isfinite(scaleHeight) || scaleHeight <= 0.0f || !isfinite(height)) {
        return 0.0f;
    }
    float zOverH = height / scaleHeight;
    if (!isfinite(zOverH) || abs(zOverH) > VOLUME_VERTICAL_TRUNCATION_SIGMA) {
        return 0.0f;
    }
    float radialColumn = tauRef * pow(radius / referenceRadius, -1.5f);
    float density = radialColumn * exp(-0.5f * zOverH * zOverH) /
                    (TRUNCATED_GAUSSIAN_COLUMN_INTEGRAL * scaleHeight);
    return isfinite(density) && density >= 0.0f ? density : 0.0f;
}

// --- Planck spectrum through the Wyman-Sloan-Shirley CIE 1931 fit -----------

float CieGaussian(float wavelengthNm, float centre, float leftWidth,
                  float rightWidth, float amplitude) {
    float width = wavelengthNm < centre ? leftWidth : rightWidth;
    float offset = (wavelengthNm - centre) * width;
    return amplitude * exp(-0.5f * offset * offset);
}

// Device mirror of core::colour::Cie1931TwoDegreeFit. This is an analytic
// approximation to the CIE 1931 2-degree observer, not the CIE tabulation.
public float3 Cie1931TwoDegreeFit(float wavelengthNm) {
    if (!isfinite(wavelengthNm) || wavelengthNm < 380.0f || wavelengthNm > 780.0f) {
        return float3(0.0f);
    }
    float X = CieGaussian(wavelengthNm, 442.0f, 0.0624f, 0.0374f, 0.362f) +
              CieGaussian(wavelengthNm, 599.8f, 0.0264f, 0.0323f, 1.056f) -
              CieGaussian(wavelengthNm, 501.1f, 0.0490f, 0.0382f, 0.065f);
    float Y = CieGaussian(wavelengthNm, 568.8f, 0.0213f, 0.0247f, 0.821f) +
              CieGaussian(wavelengthNm, 530.9f, 0.0613f, 0.0322f, 0.286f);
    float Z = CieGaussian(wavelengthNm, 437.0f, 0.0845f, 0.0278f, 1.217f) +
              CieGaussian(wavelengthNm, 459.0f, 0.0385f, 0.0725f, 0.681f);
    return float3(X, Y, Z);
}

// Device mirror of core::colour::XyzD65ToLinearSrgb. The ratios are the exact
// matrix implied by the sRGB primary chromaticities and D65 white point; the
// direct parity probe and repository governance prevent either mirror drifting.
public float3 XyzD65ToLinearSrgb(float3 xyz) {
    float3 rgb;
    rgb.x = (12831.0f / 3959.0f) * xyz.x - (329.0f / 214.0f) * xyz.y -
            (1974.0f / 3959.0f) * xyz.z;
    rgb.y = -(851781.0f / 878810.0f) * xyz.x + (1648619.0f / 878810.0f) * xyz.y +
            (36519.0f / 878810.0f) * xyz.z;
    rgb.z = (705.0f / 12673.0f) * xyz.x - (2585.0f / 12673.0f) * xyz.y +
            (705.0f / 667.0f) * xyz.z;
    return rgb;
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
        xyz += Cie1931TwoDegreeFit(wavelengthNm) * relativeRadiance;
    }

    float3 rgb = XyzD65ToLinearSrgb(xyz);
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

public float FullPageThorneFluxShape(float r, float rInner, float aStar) {
    if (r <= rInner * 1.001f) return 0.0f;
    float energy = PtSpecificEnergy(r, aStar);
    float angularMomentum = PtSpecificAngularMomentum(r, aStar);
    float omega = PtAngularVelocity(r, aStar);
    float halfWidth = 0.5f * (r - rInner);
    float midpoint = 0.5f * (r + rInner);
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

public float FullPageThorneTemperature(float temperatureScale, float r, float rInner,
                                       float M, float aStar) {
    float radiusInM = r / M;
    float innerInM = rInner / M;
    float flux = FullPageThorneFluxShape(radiusInM, innerInM, aStar);
    float reference = FullPageThorneFluxShape(
        DISK_TEMPERATURE_REFERENCE_RADIUS_RATIO * innerInM, innerInM, aStar);
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

// The live path advances a past-directed Kerr-Schild Cartesian tangent. Form
// the future photon covector there and extract the chart-invariant stationary
// and axial Killing quantities E=-p_t and L_z=p_phi=x p_y-y p_x.
public bool KerrPhotonKillingQuantitiesCart(Vec4Cart position, Vec4Cart pastVelocity,
                                            Real metric[4][4], out float energy,
                                            out float angularMomentum) {
    Real k[4] = {
        -pastVelocity.t, -pastVelocity.x, -pastVelocity.y, -pastVelocity.z
    };
    Real p[4] = { Real(0.0f), Real(0.0f), Real(0.0f), Real(0.0f) };
    [unroll] for (int mu = 0; mu < 4; ++mu) {
        [unroll] for (int nu = 0; nu < 4; ++nu) {
            p[mu] += metric[mu][nu] * k[nu];
        }
    }
    energy = -float(p[0]);
    angularMomentum = -float(position.y) * float(p[1]) + float(position.x) * float(p[2]);
    return isfinite(energy) && energy > 0.0f && isfinite(angularMomentum);
}

// Arbitrary-latitude Kerr stationary-frame transfer shared mathematically
// with core/relativistic_transfer.h.  cosTheta is the Boyer-Lindquist polar
// cosine, exactly z/r in Kerr-Schild Cartesian coordinates. The returned
// tuple is (g, measured frequency, Omega, u^t). Zero/false means the requested
// frame lies outside the represented timelike domain.
public bool TryKerrStationaryFrameFrequencyTransfer(
    float observerFrequency, float energy, float angularMomentum, float M, float a,
    float r, float cosTheta, float angularVelocity, out float4 transfer) {
    transfer = float4(0.0f);
    if (!isfinite(observerFrequency) || observerFrequency <= 0.0f ||
        !isfinite(energy) || energy <= 0.0f || !isfinite(angularMomentum) ||
        !isfinite(M) || M <= 0.0f || !isfinite(a) || abs(a) > M ||
        !isfinite(r) || r <= 0.0f || !isfinite(cosTheta) || abs(cosTheta) >= 1.0f ||
        !isfinite(angularVelocity)) {
        return false;
    }
    float sinThetaSquared = (1.0f - cosTheta) * (1.0f + cosTheta);
    float sigma = r * r + a * a * cosTheta * cosTheta;
    if (!isfinite(sigma) || sigma <= 0.0f || sinThetaSquared <= 0.0f) return false;
    float gtt = -(1.0f - 2.0f * M * r / sigma);
    float gtphi = -2.0f * M * a * r * sinThetaSquared / sigma;
    float gphiphi = sinThetaSquared *
        (r * r + a * a + 2.0f * M * a * a * r * sinThetaSquared / sigma);
    float norm = gtt + 2.0f * angularVelocity * gtphi +
                 angularVelocity * angularVelocity * gphiphi;
    if (!isfinite(norm) || norm >= 0.0f) return false;
    float timeComponent = 1.0f / sqrt(-norm);
    float frameFrequency = timeComponent * (energy - angularVelocity * angularMomentum);
    if (!isfinite(frameFrequency) || frameFrequency <= 0.0f) return false;
    transfer = float4(observerFrequency / frameFrequency, frameFrequency,
                      angularVelocity, timeComponent);
    return true;
}

// Exact locally non-rotating Kerr frame. This is intentionally distinct from
// the normal to a horizon-penetrating Kerr-Schild time slice, whose radial
// motion would contaminate a circular-Doppler suppression diagnostic.
public bool KerrZamoFrequencyTransfer(float observerFrequency, float energy,
                                      float angularMomentum, float M, float a, float r,
                                      float cosTheta, out float4 transfer) {
    transfer = float4(0.0f);
    if (!isfinite(M) || M <= 0.0f || !isfinite(a) || abs(a) > M ||
        !isfinite(r) || r <= 0.0f || !isfinite(cosTheta) || abs(cosTheta) >= 1.0f) {
        return false;
    }
    float sinThetaSquared = (1.0f - cosTheta) * (1.0f + cosTheta);
    float sigma = r * r + a * a * cosTheta * cosTheta;
    float gtphi = -2.0f * M * a * r * sinThetaSquared / sigma;
    float gphiphi = sinThetaSquared *
        (r * r + a * a + 2.0f * M * a * a * r * sinThetaSquared / sigma);
    if (!isfinite(gphiphi) || gphiphi <= 0.0f) return false;
    return TryKerrStationaryFrameFrequencyTransfer(
        observerFrequency, energy, angularMomentum, M, a, r, cosTheta,
        -gtphi / gphiphi, transfer);
}

// Exact equatorial circular-emitter and ZAMO transfer. x/y are reconstructed
// on the equatorial surface with rho=sqrt(r^2+a^2). No display clamp is
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

    float energy;
    float angularMomentum;
    if (!KerrPhotonKillingQuantitiesCart(position, photonVelocity, metric,
                                         energy, angularMomentum)) {
        return float2(0.0f);
    }
    float sqrtM = sqrt(M);
    float emitterOmega = sqrtM / (pow(r, 1.5f) + a * sqrtM);
    float4 emitterTransfer;
    float4 zamoTransfer;
    if (!TryKerrStationaryFrameFrequencyTransfer(
            observerFrequency, energy, angularMomentum, M, a, r, 0.0f,
            emitterOmega, emitterTransfer) ||
        !KerrZamoFrequencyTransfer(observerFrequency, energy, angularMomentum,
                                   M, a, r, 0.0f, zamoTransfer)) {
        return float2(0.0f);
    }
    return float2(emitterTransfer.x, zamoTransfer.x);
}
