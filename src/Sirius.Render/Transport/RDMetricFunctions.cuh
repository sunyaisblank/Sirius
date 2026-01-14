// RDMetricFunctions.cuh - Metric Tensor Evaluation
// Component ID: RDTR003A (Transport/MetricFunctions)

#pragma once

#include "../Integration/RDMath.cuh"
#include "../Acceleration/OptiX/RDOP003A.h" // For NumericalMetricData
#include "RDGeodesic.cuh"

namespace Sirius {

//==============================================================================
// Metric Tensor Evaluation
//==============================================================================

// Forward declarations
__device__ inline void getVolumetricSample(const Vec4& pos, const NumericalMetricData& nm, float& rho, float3& vel);


// Minkowski metric in spherical coordinates
__device__ inline void getMinkowskiMetric(const Vec4& x, float g[4][4], float g_inv[4][4]) {
    float r = x.r;
    float theta = x.theta;
    float sin_theta = sinf(theta);
    float r2 = r * r;
    
    // Initialize to zero
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            g[i][j] = 0.0f;
            g_inv[i][j] = 0.0f;
        }
    }
    
    // Metric components
    g[0][0] = -1.0f;           // g_tt
    g[1][1] = 1.0f;            // g_rr
    g[2][2] = r2;              // g_θθ
    g[3][3] = r2 * sin_theta * sin_theta;  // g_φφ
    
    // Inverse metric
    g_inv[0][0] = -1.0f;
    g_inv[1][1] = 1.0f;
    g_inv[2][2] = 1.0f / r2;
    g_inv[3][3] = 1.0f / (r2 * sin_theta * sin_theta);
}

// Schwarzschild metric
__device__ inline void getSchwarzschildMetric(const Vec4& x, float M, 
                                        float g[4][4], float g_inv[4][4]) {
    float r = x.r;
    float theta = x.theta;
    float rs = 2.0f * M;  // Schwarzschild radius
    
    // Avoid singularity
    if (r <= rs * 1.001f) {
        getMinkowskiMetric(x, g, g_inv);
        return;
    }
    
    float f = 1.0f - rs / r;  // Metric function
    float sin_theta = sinf(theta);
    float r2 = r * r;
    
    // Initialize to zero
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            g[i][j] = 0.0f;
            g_inv[i][j] = 0.0f;
        }
    }
    
    // Metric components
    g[0][0] = -f;              // g_tt
    g[1][1] = 1.0f / f;        // g_rr
    g[2][2] = r2;              // g_θθ
    g[3][3] = r2 * sin_theta * sin_theta;  // g_φφ
    
    // Inverse metric
    g_inv[0][0] = -1.0f / f;
    g_inv[1][1] = f;
    g_inv[2][2] = 1.0f / r2;
    g_inv[3][3] = 1.0f / (r2 * sin_theta * sin_theta);
}

// Kerr metric in Boyer-Lindquist coordinates
__device__ inline void getKerrMetric(const Vec4& x, float M, float a,
                               float g[4][4], float g_inv[4][4]) {
    float r = x.r;
    float theta = x.theta;
    
    float r2 = r * r;
    float a2 = a * a;
    float cos_theta = cosf(theta);
    float sin_theta = sinf(theta);
    float sin2 = sin_theta * sin_theta;
    float cos2 = cos_theta * cos_theta;
    
    float Sigma = r2 + a2 * cos2;
    float Delta = r2 - 2.0f * M * r + a2;
    
    // Check for horizon
    if (Delta <= 0.001f || Sigma <= 0.001f) {
        getMinkowskiMetric(x, g, g_inv);
        return;
    }
    
    float A = (r2 + a2) * (r2 + a2) - Delta * a2 * sin2;
    
    // Initialize to zero
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            g[i][j] = 0.0f;
            g_inv[i][j] = 0.0f;
        }
    }
    
    // Metric components
    g[0][0] = -(1.0f - 2.0f * M * r / Sigma);
    g[0][3] = -2.0f * M * r * a * sin2 / Sigma;
    g[3][0] = g[0][3];
    g[1][1] = Sigma / Delta;
    g[2][2] = Sigma;
    g[3][3] = A * sin2 / Sigma;
    
    // Inverse metric (simplified for diagonal dominance)
    float det = g[0][0] * g[3][3] - g[0][3] * g[0][3];
    g_inv[0][0] = g[3][3] / det;
    g_inv[0][3] = -g[0][3] / det;
    g_inv[3][0] = g_inv[0][3];
    g_inv[1][1] = Delta / Sigma;
    g_inv[2][2] = 1.0f / Sigma;
    g_inv[3][3] = g[0][0] / det;
}

// Reissner-Nordström metric
__device__ inline void getReissnerNordstromMetric(const Vec4& x, float M, float Q,
                                            float g[4][4], float g_inv[4][4]) {
    float r = x.r;
    float theta = x.theta;
    
    float r2 = r * r;
    float Q2 = Q * Q;
    float f = 1.0f - 2.0f * M / r + Q2 / r2;
    
    // Check for horizon
    if (f <= 0.001f || r <= 0.001f) {
        getMinkowskiMetric(x, g, g_inv);
        return;
    }
    
    float sin_theta = sinf(theta);
    
    // Initialize to zero
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            g[i][j] = 0.0f;
            g_inv[i][j] = 0.0f;
        }
    }
    
    // Metric components
    g[0][0] = -f;
    g[1][1] = 1.0f / f;
    g[2][2] = r2;
    g[3][3] = r2 * sin_theta * sin_theta;
    
    // Inverse metric
    g_inv[0][0] = -1.0f / f;
    g_inv[1][1] = f;
    g_inv[2][2] = 1.0f / r2;
    g_inv[3][3] = 1.0f / (r2 * sin_theta * sin_theta);
}

// Einstein Toolkit Numerical Metric (interpolated from HDF5 data)
__device__ inline void getNumericalMetric(const Vec4& pos, const NumericalMetricData& nm, float g[4][4], float g_inv[4][4]) {
    // 1. Convert Spherical (r, θ, φ) to Cartesian (x, y, z)
    // Sirius uses spherical coordinates for the integration
    float3 cart = make_float3(
        pos.r * sinf(pos.theta) * cosf(pos.phi),
        pos.r * sinf(pos.theta) * sinf(pos.phi),
        pos.r * cosf(pos.theta)
    );
    
    // 2. Compute texture coordinates
    // Map logical (x,y,z) to [0,1] texture space
    // x_tex = (x - x0) / (nx * dx)
    float3 uvw = make_float3(
        (cart.x - nm.origin.x) / (nm.dims.x * nm.spacing.x),
        (cart.y - nm.origin.y) / (nm.dims.y * nm.spacing.y),
        (cart.z - nm.origin.z) / (nm.dims.z * nm.spacing.z)
    );
    
    // Check bounds
    if (!nm.isLoaded || 
        uvw.x < 0.0f || uvw.x > 1.0f || 
        uvw.y < 0.0f || uvw.y > 1.0f || 
        uvw.z < 0.0f || uvw.z > 1.0f) {
        getMinkowskiMetric(pos, g, g_inv);
        return;
    }
    
    // 3. Sample grid data (Trilinear interpolation via Texture Object)
    // ADM Metric: g_ij (spatial 3x3), alpha (lapse), beta^i (shift)
    // Spacetime metric g_uv is reconstructed from 3+1 splitting:
    // g_00 = -alpha^2 + beta_k beta^k
    // g_0i = beta_i
    // g_ij = gamma_ij
    
    float gxx = tex3D<float>(nm.gxx, uvw.x, uvw.y, uvw.z);
    float gxy = tex3D<float>(nm.gxy, uvw.x, uvw.y, uvw.z);
    float gxz = tex3D<float>(nm.gxz, uvw.x, uvw.y, uvw.z);
    float gyy = tex3D<float>(nm.gyy, uvw.x, uvw.y, uvw.z);
    float gyz = tex3D<float>(nm.gyz, uvw.x, uvw.y, uvw.z);
    float gzz = tex3D<float>(nm.gzz, uvw.x, uvw.y, uvw.z);
    
    float alpha = tex3D<float>(nm.alp, uvw.x, uvw.y, uvw.z);
    float betax = tex3D<float>(nm.betax, uvw.x, uvw.y, uvw.z);
    float betay = tex3D<float>(nm.betay, uvw.x, uvw.y, uvw.z);
    float betaz = tex3D<float>(nm.betaz, uvw.x, uvw.y, uvw.z);
    
    // Handle vacuum regions (alpha -> 1, shift -> 0, g_ij -> delta_ij)
    if (alpha < 0.001f) alpha = 1.0f;
    
    // 4. Reconstruct 4-Metric (in Cartesian basis)
    float gc[4][4]; // g_cartesian
    
    // Lower shift beta_i = gamma_ij beta^j
    float beta_cov_x = gxx*betax + gxy*betay + gxz*betaz;
    float beta_cov_y = gxy*betax + gyy*betay + gyz*betaz;
    float beta_cov_z = gxz*betax + gyz*betay + gzz*betaz;
    
    // beta_sq = beta_k beta^k
    float beta_sq = betax*beta_cov_x + betay*beta_cov_y + betaz*beta_cov_z;
    
    gc[0][0] = -alpha*alpha + beta_sq;
    gc[0][1] = beta_cov_x; gc[0][2] = beta_cov_y; gc[0][3] = beta_cov_z;
    gc[1][0] = beta_cov_x; gc[2][0] = beta_cov_y; gc[3][0] = beta_cov_z;
    
    gc[1][1] = gxx; gc[1][2] = gxy; gc[1][3] = gxz;
    gc[2][1] = gxy; gc[2][2] = gyy; gc[2][3] = gyz;
    gc[3][1] = gxz; gc[3][2] = gyz; gc[3][3] = gzz;
    
    // 5. Transform to Spherical Basis: g_sph_mn = J^a_m J^b_n g_cart_ab
    // Jacobian J^a_m = dx^a / dx^m
    // x^0=t, x^1=r, x^2=theta, x^3=phi
    // y^0=t, y^1=x, y^2=y, y^3=z
    
    float st = sinf(pos.theta), ct = cosf(pos.theta);
    float sp = sinf(pos.phi), cp = cosf(pos.phi);
    float r = pos.r;
    
    // J[row][col] = J[cart][sph]
    float J[4][4] = {
        {1.0f, 0.0f,       0.0f,        0.0f},         // dt/dt, dt/dr...
        {0.0f, st*cp,      r*ct*cp,     -r*st*sp},     // dx/dt, dx/dr...
        {0.0f, st*sp,      r*ct*sp,     r*st*cp},      // dy/dt...
        {0.0f, ct,         -r*st,       0.0f}          // dz/dt...
    };
    
    // Transform g
    for (int m = 0; m < 4; m++) {
        for (int n = 0; n < 4; n++) {
            float sum = 0.0f;
            for (int a = 0; a < 4; a++) {
                for (int b = 0; b < 4; b++) {
                    sum += J[a][m] * J[b][n] * gc[a][b];
                }
            }
            g[m][n] = sum;
        }
    }
    
    // Compute inverse metric
    if (!invert4x4(g, g_inv)) {
        // Fallback if singular
        for(int i=0; i<4; i++) for(int j=0; j<4; j++) g_inv[i][j] = (i==j) ? (i==0 ? -1.0f : 1.0f) : 0.0f;
    } 
}

// Gödel Metric (Rotating Cosmology with Closed Timelike Curves)
__device__ inline void getGodelMetric(const Vec4& x, float a, 
                                float g[4][4], float g_inv[4][4]) {
    // Gödel metric: ds² = -dt² - 2ar dt dφ + dr² + (r² - a²r⁴)dφ² + dz²
    // In spherical form where theta->phi (azimuth), phi->z
    float r = fmaxf(x.r, 0.01f);  // Avoid singularity
    float r2 = r * r;
    float a2 = a * a;
    
    float g_phiphi = r2 * (1.0f - a2 * r2);  // Can become negative for CTCs
    
    // Initialize to zero
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            g[i][j] = 0.0f;
            g_inv[i][j] = 0.0f;
        }
    }
    
    // Metric components
    g[0][0] = -1.0f;           // g_tt
    g[0][2] = -a * r;          // g_tφ (using theta slot for φ)
    g[2][0] = -a * r;          // g_φt
    g[1][1] = 1.0f;            // g_rr
    g[2][2] = g_phiphi;        // g_φφ
    g[3][3] = 1.0f;            // g_zz (using phi slot for z)
    
    // Inverse metric (handle potential CTC region where g_φφ can be negative)
    float det_tc = g[0][0] * g[2][2] - g[0][2] * g[0][2];  // det of t-φ block
    if (fabsf(det_tc) > 0.001f) {
        g_inv[0][0] = g[2][2] / det_tc;
        g_inv[0][2] = -g[0][2] / det_tc;
        g_inv[2][0] = g_inv[0][2];
        g_inv[2][2] = g[0][0] / det_tc;
    } else {
        // Near singular, use Minkowski fallback
        g_inv[0][0] = -1.0f;
        g_inv[2][2] = 1.0f / fmaxf(fabsf(g_phiphi), 0.01f);
    }
    g_inv[1][1] = 1.0f;
    g_inv[3][3] = 1.0f;
}

// Taub-NUT Metric (Gravimagnetic Monopole)
__device__ inline void getTaubNUTMetric(const Vec4& x, float M, float n,
                                  float g[4][4], float g_inv[4][4]) {
    // ds² = -f(dt + 2n cosθ dφ)² + dr²/f + Σ(dθ² + sin²θ dφ²)
    // where f = (r² - 2Mr - n²)/(r² + n²), Σ = r² + n²
    
    float r = x.r;
    float theta = x.theta;
    
    float r2 = r * r;
    float n2 = n * n;
    float Sigma = r2 + n2;
    
    float numerator = r2 - 2.0f * M * r - n2;
    float f = numerator / Sigma;
    
    // Avoid horizon singularities
    float f_safe = fmaxf(fabsf(f), 0.001f);
    if (f < 0.0f) f_safe = -f_safe;
    
    float sin_theta = sinf(theta);
    float cos_theta = cosf(theta);
    float sin2 = sin_theta * sin_theta;
    float cos2 = cos_theta * cos_theta;
    
    // Handle pole singularity (Misner string)
    if (fabsf(sin_theta) < 0.01f) {
        sin_theta = (sin_theta >= 0.0f) ? 0.01f : -0.01f;
        sin2 = sin_theta * sin_theta;
    }
    
    // Initialize to zero
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            g[i][j] = 0.0f;
            g_inv[i][j] = 0.0f;
        }
    }
    
    // Metric components: g_tt, g_tφ, g_rr, g_θθ, g_φφ
    g[0][0] = -f;                                        // g_tt
    g[0][3] = -2.0f * f * n * cos_theta;                 // g_tφ
    g[3][0] = g[0][3];                                   // g_φt
    g[1][1] = 1.0f / f_safe;                             // g_rr
    g[2][2] = Sigma;                                     // g_θθ
    g[3][3] = Sigma * sin2 - 4.0f * f * n2 * cos2;       // g_φφ
    
    // Inverse metric
    // For the t-φ block
    float det_tphi = g[0][0] * g[3][3] - g[0][3] * g[0][3];
    if (fabsf(det_tphi) > 0.001f) {
        g_inv[0][0] = g[3][3] / det_tphi;
        g_inv[0][3] = -g[0][3] / det_tphi;
        g_inv[3][0] = g_inv[0][3];
        g_inv[3][3] = g[0][0] / det_tphi;
    } else {
        // Fallback
        g_inv[0][0] = -1.0f / f_safe;
        g_inv[3][3] = 1.0f / fmaxf(fabsf(g[3][3]), 0.01f);
    }
    g_inv[1][1] = f_safe;
    g_inv[2][2] = 1.0f / Sigma;
}

// Kerr-Schild Metric (Horizon-Penetrating Coordinates)
__device__ inline void getKerrSchildMetric(const Vec4& x, float M, float a,
                                     float g[4][4], float g_inv[4][4]) {
    float r = fmaxf(x.r, 0.01f);
    float theta = x.theta;
    
    float r2 = r * r;
    float a2 = a * a;
    float cos_theta = cosf(theta);
    float sin_theta = sinf(theta);
    float sin2 = sin_theta * sin_theta;
    float cos2 = cos_theta * cos_theta;
    
    // Kerr-Schild H function
    float Sigma = r2 + a2 * cos2;
    float H = 2.0f * M * r / fmaxf(Sigma, 0.01f);
    
    // Initialize to zero
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            g[i][j] = 0.0f;
            g_inv[i][j] = 0.0f;
        }
    }
    
    // Metric components (simplified Kerr-Schild in spherical-like coords)
    g[0][0] = -(1.0f - H);            // g_tt
    g[0][1] = H;                       // g_tr
    g[1][0] = H;                       // g_rt
    g[1][1] = 1.0f + H;                // g_rr
    g[2][2] = Sigma;                   // g_θθ
    g[3][3] = Sigma * sin2;            // g_φφ
    
    // Inverse metric
    float det_tr = g[0][0] * g[1][1] - g[0][1] * g[1][0];
    if (fabsf(det_tr) > 0.001f) {
        g_inv[0][0] = g[1][1] / det_tr;
        g_inv[0][1] = -g[0][1] / det_tr;
        g_inv[1][0] = -g[1][0] / det_tr;
        g_inv[1][1] = g[0][0] / det_tr;
    } else {
        g_inv[0][0] = -1.0f;
        g_inv[1][1] = 1.0f;
    }
    g_inv[2][2] = 1.0f / Sigma;
    g_inv[3][3] = 1.0f / (Sigma * sin2 + 0.0001f);
}

// Ellis Drainhole (Traversable Wormhole) Metric
__device__ inline void getEllisDrainholeMetric(const Vec4& x, float m, float n,
                                         float g[4][4], float g_inv[4][4]) {
    float r = x.r;
    
    // Ensure n > m for valid wormhole
    if (n <= m) n = m + 0.1f;
    
    // Throat radius: α = √(n² - m²)
    float alpha = sqrtf(n * n - m * m);
    
    // Proper radial coordinate
    float dx = r - m;
    float Rp2 = dx * dx + alpha * alpha;
    float Rp = sqrtf(Rp2);
    
    // Pseudo-angle ψ
    float psi = (n / alpha) * (1.5707963f - atanf(dx / alpha));  // π/2 - atan(...)
    
    // Redshift function Fp
    float exp_term = expf(-2.0f * m * psi / n);
    float Fp_sq = 1.0f - exp_term;
    float Fp = -sqrtf(fmaxf(Fp_sq, 0.0f));  // Negative for ingoing
    
    // Areal radius factor
    float denom = 1.0f - Fp_sq;
    if (denom < 0.001f) denom = 0.001f;
    float Rp_areal = Rp / sqrtf(denom);
    float Rp_areal2 = Rp_areal * Rp_areal;
    
    float theta = x.theta;
    float sin_theta = sinf(theta);
    float sin2 = sin_theta * sin_theta;
    
    // Initialize to zero
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            g[i][j] = 0.0f;
            g_inv[i][j] = 0.0f;
        }
    }
    
    // Metric components
    g[0][0] = -(1.0f - Fp_sq);
    g[0][1] = Fp;
    g[1][0] = Fp;
    g[1][1] = 1.0f;
    g[2][2] = Rp_areal2;
    g[3][3] = Rp_areal2 * sin2;
    
    // Inverse metric
    float det_tr = g[0][0] * g[1][1] - g[0][1] * g[1][0];
    if (fabsf(det_tr) > 0.001f) {
        g_inv[0][0] = g[1][1] / det_tr;
        g_inv[0][1] = -g[0][1] / det_tr;
        g_inv[1][0] = -g[1][0] / det_tr;
        g_inv[1][1] = g[0][0] / det_tr;
    } else {
        g_inv[0][0] = -1.0f;
        g_inv[1][1] = 1.0f;
    }
    g_inv[2][2] = 1.0f / Rp_areal2;
    g_inv[3][3] = 1.0f / (Rp_areal2 * sin2 + 0.0001f);
}


// Alcubierre Warp Drive Metric
// ds² = (vs²f² - 1)dt² - 2vs·f·dx·dt + dx² + dy² + dz²
__device__ inline void getAlcubierreMetric(const Vec4& x, float vs, float sigma, float R,
                                     float g[4][4], float g_inv[4][4]) {
    // For spherical → Cartesian conversion
    float r = x.r;
    float theta = x.theta;
    float t = x.t;
    
    // Bubble center moves along x-axis: xs(t) = vs * t
    float xs_t = vs * t;
    
    // Current position in Cartesian
    float cart_x = r * sinf(theta) * cosf(x.phi);
    float cart_y = r * sinf(theta) * sinf(x.phi);
    float cart_z = r * cosf(theta);
    
    // Distance from bubble center
    float dx = cart_x - xs_t;
    float rs2 = dx*dx + cart_y*cart_y + cart_z*cart_z;
    float rs = sqrtf(rs2);
    if (rs < 0.001f) rs = 0.001f;
    
    // Warp bubble shape function: f(rs) = (tanh(σ(rs+R)) - tanh(σ(rs-R))) / (2·tanh(σR))
    float tanh_plus = tanhf(sigma * (rs + R));
    float tanh_minus = tanhf(sigma * (rs - R));
    float tanh_R = tanhf(sigma * R);
    
    float f = 1.0f;
    if (fabsf(tanh_R) > 0.001f) {
        f = (tanh_plus - tanh_minus) / (2.0f * tanh_R);
    }
    
    float vs_f = vs * f;
    
    // Metric in Cartesian (simplified to spherical r direction)
    float sin_theta = sinf(theta);
    float r2 = r * r;
    
    // Initialize
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            g[i][j] = 0.0f;
            g_inv[i][j] = 0.0f;
        }
    }
    
    // Approximate effect in spherical: bubble distorts along radial direction
    g[0][0] = vs_f * vs_f - 1.0f;  // g_tt
    g[0][1] = -vs_f;               // g_tr (radial direction ~x)
    g[1][0] = -vs_f;               // g_rt
    g[1][1] = 1.0f;                // g_rr
    g[2][2] = r2;                  // g_θθ
    g[3][3] = r2 * sin_theta * sin_theta;  // g_φφ
    
    // Inverse (t-r block)
    float det_tr = g[0][0] * g[1][1] - g[0][1] * g[1][0];
    if (fabsf(det_tr) > 0.001f) {
        g_inv[0][0] = g[1][1] / det_tr;
        g_inv[0][1] = -g[0][1] / det_tr;
        g_inv[1][0] = -g[1][0] / det_tr;
        g_inv[1][1] = g[0][0] / det_tr;
    } else {
        g_inv[0][0] = -1.0f;
        g_inv[1][1] = 1.0f;
    }
    g_inv[2][2] = 1.0f / r2;
    g_inv[3][3] = 1.0f / (r2 * sin_theta * sin_theta + 0.0001f);
}

// de Sitter Metric (Cosmological Spacetime with Λ > 0)
__device__ inline void getDeSitterMetric(const Vec4& x, float Lambda,
                                   float g[4][4], float g_inv[4][4]) {
    float r = x.r;
    float theta = x.theta;
    float r2 = r * r;
    
    // f(r) = 1 - Λr²/3
    float f = 1.0f - Lambda * r2 / 3.0f;
    
    // Avoid cosmological horizon
    if (fabsf(f) < 0.001f) {
        f = (f > 0.0f) ? 0.001f : -0.001f;
    }
    
    float sin_theta = sinf(theta);
    float sin2 = sin_theta * sin_theta;
    
    // Initialize
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            g[i][j] = 0.0f;
            g_inv[i][j] = 0.0f;
        }
    }
    
    // Diagonal metric
    g[0][0] = -f;                  // g_tt
    g[1][1] = 1.0f / f;            // g_rr
    g[2][2] = r2;                  // g_θθ
    g[3][3] = r2 * sin2;           // g_φφ
    
    // Inverse
    g_inv[0][0] = -1.0f / f;
    g_inv[1][1] = f;
    g_inv[2][2] = 1.0f / r2;
    g_inv[3][3] = 1.0f / (r2 * sin2 + 0.0001f);
}

// Sample volumetric data (rho, velocity)
__device__ inline void getVolumetricSample(const Vec4& pos, const NumericalMetricData& nm, float& rho, float3& vel) {
    rho = 0.0f;
    vel = make_float3(0.0f, 0.0f, 0.0f);

    if (!nm.isLoaded) return;

    // 1. Convert Spherical (r, θ, φ) to Cartesian (x, y, z)
    // Using Y-POLAR convention: y = r cos(θ), consistent with rest of codebase
    // θ = polar angle from +Y axis, φ = azimuthal angle in XZ plane
    float3 cart = make_float3(
        pos.r * sinf(pos.theta) * cosf(pos.phi),   // x = r sinθ cosφ
        pos.r * cosf(pos.theta),                    // y = r cosθ  (polar axis)
        pos.r * sinf(pos.theta) * sinf(pos.phi)    // z = r sinθ sinφ
    );

    // 2. Compute texture coordinates
    float3 uvw = make_float3(
        (cart.x - nm.origin.x) / (nm.dims.x * nm.spacing.x),
        (cart.y - nm.origin.y) / (nm.dims.y * nm.spacing.y),
        (cart.z - nm.origin.z) / (nm.dims.z * nm.spacing.z)
    );

    // Check bounds
    if (uvw.x < 0.0f || uvw.x > 1.0f || 
        uvw.y < 0.0f || uvw.y > 1.0f || 
        uvw.z < 0.0f || uvw.z > 1.0f) {
        return;
    }

    // 3. Sample data
    rho = tex3D<float>(nm.rho, uvw.x, uvw.y, uvw.z);
    
    if (rho > 0.0001f) {
        vel.x = tex3D<float>(nm.vx, uvw.x, uvw.y, uvw.z);
        vel.y = tex3D<float>(nm.vy, uvw.x, uvw.y, uvw.z);
        vel.z = tex3D<float>(nm.vz, uvw.x, uvw.y, uvw.z);
    }
}

// Compute Kerr radius from Cartesian coordinates
__device__ inline float computeKerrRadiusCart(float x, float y, float z, float a) {
    float a2 = a * a;
    float R2 = x*x + y*y + z*z;
    
    // For a = 0 (Schwarzschild), r = R directly
    if (fabsf(a) < 1e-10f) {
        return sqrtf(fmaxf(R2, 1e-20f));
    }
    
    // Solve: r² = [(R² - a²) + √((R² - a²)² + 4a²z²)] / 2
    float Rm2 = R2 - a2;
    float disc = Rm2 * Rm2 + 4.0f * a2 * z * z;
    float r2 = (Rm2 + sqrtf(fmaxf(disc, 0.0f))) / 2.0f;
    
    return sqrtf(fmaxf(r2, 1e-20f));
}

// Unified Kerr-Schild metric in Cartesian coordinates
__device__ inline void getKerrSchildFamilyMetric(
    const Vec4Cart& pos,            // Position (t, x, y, z)
    float M, float a, float Q, float Lambda,  // Family parameters
    float g[4][4], float g_inv[4][4])
{
    float x = pos.x;
    float y = pos.y;
    float z = pos.z;
    
    float a2 = a * a;
    float Q2 = Q * Q;
    
    // Initialize to Minkowski
    #pragma unroll 4
    for (int i = 0; i < 4; i++) {
        #pragma unroll 4
        for (int j = 0; j < 4; j++) {
            g[i][j] = 0.0f;
            g_inv[i][j] = 0.0f;
        }
    }
    g[0][0] = -1.0f;  g_inv[0][0] = -1.0f;
    g[1][1] = 1.0f;   g_inv[1][1] = 1.0f;
    g[2][2] = 1.0f;   g_inv[2][2] = 1.0f;
    g[3][3] = 1.0f;   g_inv[3][3] = 1.0f;
    
    // If no mass and no cosmological constant, we're done (Minkowski)
    if (M < 1e-15f && fabsf(Lambda) < 1e-15f) {
        return;
    }
    
    // Compute Kerr radius
    float r = computeKerrRadiusCart(x, y, z, a);
    float r2 = r * r;
    float r4 = r2 * r2;
    
    // Kerr-Schild scalar function H = (2Mr - Q²)r² / (r⁴ + a²z²)
    float f_denom = r4 + a2 * z * z;
    float H = (2.0f * M * r - Q2) * r2 / fmaxf(f_denom, 1e-20f);
    
    // Null vector l_μ = (1, (rx + ay)/(r² + a²), (ry - ax)/(r² + a²), z/r)
    float denom = r2 + a2;
    float l[4] = {
        1.0f,
        (r * x + a * y) / denom,
        (r * y - a * x) / denom,
        z / fmaxf(r, 1e-10f)
    };
    
    // Add Kerr-Schild perturbation: g_μν = η_μν + H · l_μ · l_ν
    #pragma unroll 4
    for (int mu = 0; mu < 4; mu++) {
        #pragma unroll 4
        for (int nu = 0; nu < 4; nu++) {
            g[mu][nu] += H * l[mu] * l[nu];
        }
    }
    
    // Inverse metric: g^μν = η^μν - H · l^μ · l^ν
    // l^μ = (1, -l_1, -l_2, -l_3) since raising with Minkowski flips spatial
    float l_up[4] = {1.0f, -l[1], -l[2], -l[3]};
    #pragma unroll 4
    for (int mu = 0; mu < 4; mu++) {
        #pragma unroll 4
        for (int nu = 0; nu < 4; nu++) {
            g_inv[mu][nu] -= H * l_up[mu] * l_up[nu];
        }
    }
    
    // Cosmological constant contribution (de Sitter sector)
    if (fabsf(Lambda) > 1e-15f) {
        float R2 = x*x + y*y + z*z;
        float LambdaFactor = Lambda * R2 / 3.0f;
        g[0][0] -= LambdaFactor;
        g_inv[0][0] += LambdaFactor;  // Approximate for small Lambda
    }
}

// Wrapper to accept spherical Vec4 (converts internally)
__device__ inline void getKerrSchildFamilyMetricFromSpherical(
    const Vec4& pos_sph,            // Position (t, r, θ, φ)
    float M, float a, float Q, float Lambda,
    float g[4][4], float g_inv[4][4])
{
    Vec4Cart pos_cart = vec4SphToCart(pos_sph);
    getKerrSchildFamilyMetric(pos_cart, M, a, Q, Lambda, g, g_inv);
}


// Helper: Compute gradients of r (Kerr radius) and H (scalar function)
__device__ inline void computeKerrSchildGradients(
    const Vec4Cart& pos, float r, float M, float a, float Q,
    float dr_dx[4], float dH_dx[4], float dl_dx[4][4])
{
    float x = pos.x;
    float y = pos.y;
    float z = pos.z;
    float r2 = r * r;
    float a2 = a * a;
    float Q2 = Q * Q;
    
    // 1. Compute Gradient of r (dr/dx^mu)
    // Implicit r: r^4 - (R^2-a^2)r^2 - a^2 z^2 = 0
    // F(r,x,y,z) = 0 => dr/xi = - (dF/dxi) / (dF/dr)
    // dF/dr = 4r^3 - 2(R^2-a^2)r = 2r (2r^2 - R^2 + a^2)
    // dF/dx = -r^2 * 2x
    // dF/dy = -r^2 * 2y
    // dF/dz = -r^2 * 2z - 2a^2 z
    
    float R2 = x*x + y*y + z*z;
    float denom_r = 2.0f * r * (2.0f * r2 - R2 + a2);
    
    // Avoid division by zero
    float inv_denom_r = (fabsf(denom_r) > 1e-10f) ? 1.0f / denom_r : 0.0f;
    
    dr_dx[0] = 0.0f; // Static
    dr_dx[1] = (r2 * 2.0f * x) * inv_denom_r;
    dr_dx[2] = (r2 * 2.0f * y) * inv_denom_r;
    dr_dx[3] = (r2 * 2.0f * z + 2.0f * a2 * z) * inv_denom_r;
    
    // 2. Compute Gradient of H
    // H = (2Mr - Q^2) r^2 / (r^4 + a^2 z^2)
    // Let N = 2Mr^3 - Q^2 r^2
    // Let D = r^4 + a^2 z^2
    // dH = (dN * D - N * dD) / D^2
    
    float N = 2.0f * M * r * r2 - Q2 * r2;
    float D = r2 * r2 + a2 * z * z;
    float inv_D = 1.0f / fmaxf(D, 1e-20f);
    float inv_D2 = inv_D * inv_D;
    
    // Gradients of N and D
    // dN/dx^i = (6Mr^2 - 2Q^2r) * dr/dx^i
    // dD/dx^i = 4r^3 * dr/dx^i + a^2 * d(z^2)/dx^i
    
    float dN_dr = 6.0f * M * r2 - 2.0f * Q2 * r;
    float dD_dr = 4.0f * r * r2;
    
    // Calculate partials
    dH_dx[0] = 0.0f;
    
    for (int i = 1; i <= 3; i++) {
        float dN_dxi = dN_dr * dr_dx[i];
        float dD_dxi = dD_dr * dr_dx[i];
        if (i == 3) dD_dxi += 2.0f * a2 * z; // Add term for z
        
        dH_dx[i] = (dN_dxi * D - N * dD_dxi) * inv_D2;
    }
    
    // 3. Compute Gradients of l_mu
    // l = (1, (rx+ay)/(r^2+a^2), (ry-ax)/(r^2+a^2), z/r)
    // Let den = r^2 + a^2
    
    float den = r2 + a2;
    float inv_den = 1.0f / den;
    
    // Terms for derivatives
    // d(1/den)/xi = -1/den^2 * 2r * dr/xi
    float d_inv_den_dr = -inv_den * inv_den * 2.0f * r;
    
    // l_t = 1 -> grad = 0
    #pragma unroll 4
    for(int k=0; k<4; k++) dl_dx[0][k] = 0.0f;
    
    // l_x = (rx + ay) / den
    // d(lx)/xi = (d(rx+ay)/xi * den - (rx+ay) * d(den)/xi) / den^2
    //          = d(numerator)/xi * inv_den + numerator * d(inv_den)/xi
    float num_x = r*x + a*y;
    for (int i = 1; i <= 3; i++) {
        float d_num = dr_dx[i] * x + ((i==1)?r:0.0f) + ((i==2)?a:0.0f);
        dl_dx[1][i] = d_num * inv_den + num_x * d_inv_den_dr * dr_dx[i];
    }
    dl_dx[1][0] = 0.0f;
    
    // l_y = (ry - ax) / den
    float num_y = r*y - a*x;
    for (int i = 1; i <= 3; i++) {
        float d_num = dr_dx[i] * y + ((i==2)?r:0.0f) - ((i==1)?a:0.0f);
        dl_dx[2][i] = d_num * inv_den + num_y * d_inv_den_dr * dr_dx[i];
    }
    dl_dx[2][0] = 0.0f;
    
    // l_z = z / r
    for (int i = 1; i <= 3; i++) {
        float term1 = ((i==3)?1.0f:0.0f) / r;
        float term2 = -z / r2 * dr_dx[i];
        dl_dx[3][i] = term1 + term2;
    }
    dl_dx[3][0] = 0.0f;
}

// UNIFIED KERR-SCHILD BLACK HOLE FAMILY - Christoffel Symbols
__device__ inline void getKerrSchildFamilyChristoffel(
    const Vec4Cart& pos,            // Position (t, x, y, z)
    float M, float a, float Q, float Lambda,  // Family parameters
    float Gamma[4][4][4])
{
    // Initialize to zero
    #pragma unroll 4
    for (int i = 0; i < 4; i++)
        #pragma unroll 4
        for (int j = 0; j < 4; j++)
            #pragma unroll 4
            for (int k = 0; k < 4; k++)
                Gamma[i][j][k] = 0.0f;
    
    // For pure Minkowski (M=0, Λ=0), all Christoffels are zero in Cartesian
    if (M < 1e-15f && fabsf(Lambda) < 1e-15f) {
        return;
    }
    
    float x = pos.x;
    float y = pos.y;
    float z = pos.z;
    
    // Compute Kerr radius
    float r = computeKerrRadiusCart(x, y, z, a);
    float r2 = r * r;
    float a2 = a * a;
    float Q2 = Q * Q;
    float den = r2 + a2;
    
    // Metric function H and vector l
    float H = (2.0f * M * r - Q2) * r2 / (r2*r2 + a2*z*z);
    
    float l[4] = {
        1.0f,
        (r * x + a * y) / den,
        (r * y - a * x) / den,
        z / r
    };
    
    // Compute gradients
    float dr_dx[4];
    float dH_dx[4];
    float dl_dx[4][4]; // [component][derivative_direction]
    
    computeKerrSchildGradients(pos, r, M, a, Q, dr_dx, dH_dx, dl_dx);
    
    // Inverse metric (needed for raising indices)
    // g^uv = eta^uv - H l^u l^v
    // l^u = (1, -lx, -ly, -lz)
    float l_up[4] = {l[0], -l[1], -l[2], -l[3]};
    float g_inv[4][4];
    
    // Diagonal Minkowski
    for(int i=0; i<4; i++) for(int j=0; j<4; j++) g_inv[i][j] = 0.0f;
    g_inv[0][0] = -1.0f; g_inv[1][1] = 1.0f; g_inv[2][2] = 1.0f; g_inv[3][3] = 1.0f;
    
    for(int u=0; u<4; u++) {
        for(int v=0; v<4; v++) {
            g_inv[u][v] -= H * l_up[u] * l_up[v];
        }
    }
    
    // Compute lower index derivatives: dg_ab/dx^c
    // g_ab = eta_ab + H l_a l_b
    // d(g_ab)/dxc = (dH/dxc) l_a l_b + H (dl_a/dxc l_b + l_a dl_b/dxc)
    
    // Temporary storage for relevant derivatives to save registers? 
    // Just compute usage on the fly might be better but redundant.
    // Let's iterate.
    
    // Gamma^u_ab = 0.5 * g^ud * (dg_db/dx^a + dg_da/dx^b - dg_ab/dx^d)
    
    for (int u = 0; u < 4; u++) {
        for (int a_idx = 0; a_idx < 4; a_idx++) {
            for (int b_idx = a_idx; b_idx < 4; b_idx++) { // Symmetric
                
                float sum = 0.0f;
                // Contraction over d
                for (int d = 0; d < 4; d++) {
                    float G_ud = g_inv[u][d];
                    if (fabsf(G_ud) < 1e-15f) continue;
                    
                    // Term 1: dg_db/dx^a
                    float dg_db_da = dH_dx[a_idx] * l[d] * l[b_idx] + 
                                     H * (dl_dx[d][a_idx] * l[b_idx] + l[d] * dl_dx[b_idx][a_idx]);
                                     
                    // Term 2: dg_da/dx^b
                    float dg_da_db = dH_dx[b_idx] * l[d] * l[a_idx] + 
                                     H * (dl_dx[d][b_idx] * l[a_idx] + l[d] * dl_dx[a_idx][b_idx]);
                                     
                    // Term 3: -dg_ab/dx^d
                    float dg_ab_dd = dH_dx[d] * l[a_idx] * l[b_idx] + 
                                     H * (dl_dx[a_idx][d] * l[b_idx] + l[a_idx] * dl_dx[b_idx][d]);
                                     
                    sum += G_ud * (dg_db_da + dg_da_db - dg_ab_dd);
                }
                
                Gamma[u][a_idx][b_idx] = 0.5f * sum;
                if (a_idx != b_idx) Gamma[u][b_idx][a_idx] = 0.5f * sum;
            }
        }
    }
}

// Morris-Thorne wormhole Christoffel symbols
__device__ inline void getMorrisThorneChristoffel(
    const Vec4& pos,
    float b0,
    float Phi0,
    WormholeShape shapeType,
    float Gamma[4][4][4])
{
    // Initialize to zero
    #pragma unroll 4
    for (int i = 0; i < 4; i++)
        #pragma unroll 4
        for (int j = 0; j < 4; j++)
            #pragma unroll 4
            for (int k = 0; k < 4; k++)
                Gamma[i][j][k] = 0.0f;
    
    float r = fmaxf(pos.r, b0 * 1.001f);
    float theta = pos.theta;
    
    float sin_theta = sinf(theta);
    float cos_theta = cosf(theta);
    float sin2 = sin_theta * sin_theta;
    float r2 = r * r;
    
    // Shape function and derivatives
    float b = wormholeShapeFunction(r, b0, shapeType);
    float db_dr = wormholeShapeDerivative(r, b0, shapeType);
    
    float one_minus_b_over_r = 1.0f - b / r;
    one_minus_b_over_r = fmaxf(one_minus_b_over_r, 0.001f);
    
    // d(b/r)/dr = (r*db_dr - b) / r² = db_dr/r - b/r²
    float d_b_over_r = db_dr / r - b / r2;
    
    // Γ^r_rr = d(1-b/r)/dr / (2(1-b/r)) = -d_b_over_r / (2(1-b/r))
    Gamma[1][1][1] = -d_b_over_r / (2.0f * one_minus_b_over_r);
    
    // Γ^r_θθ = -r(1-b/r) = -(r - b)
    Gamma[1][2][2] = -(r - b);
    
    // Γ^r_φφ = -r sin²θ (1-b/r) = -(r - b) sin²θ
    Gamma[1][3][3] = -(r - b) * sin2;
    
    // Γ^θ_rθ = Γ^θ_θr = 1/r
    Gamma[2][1][2] = 1.0f / r;
    Gamma[2][2][1] = 1.0f / r;
    
    // Γ^θ_φφ = -sinθ cosθ
    Gamma[2][3][3] = -sin_theta * cos_theta;
    
    // Γ^φ_rφ = Γ^φ_φr = 1/r
    Gamma[3][1][3] = 1.0f / r;
    Gamma[3][3][1] = 1.0f / r;
    
    // Γ^φ_θφ = Γ^φ_φθ = cotθ
    float cot_theta = cos_theta / fmaxf(fabsf(sin_theta), 1e-10f); // safe cot
    if (sin_theta < 0.0f) cot_theta = -cot_theta; // preserve sign

    Gamma[3][2][3] = cot_theta;
    Gamma[3][3][2] = cot_theta;
    
    // If Φ ≠ 0, add time-related Christoffels
    if (fabsf(Phi0) > 1e-10f) {
        // For r-dependent Φ, we would add Γ^t_tr = dΦ/dr terms
        // With constant Φ, these are zero
    }
}

// Warp drive Christoffel symbols
__device__ inline void getWarpDriveChristoffel(
    const Vec4Cart& pos,
    float vs, float sigma, float R,
    float xs, float ys, float zs,
    float Gamma[4][4][4])
{
    // Initialize to zero
    #pragma unroll 4
    for (int i = 0; i < 4; i++)
        #pragma unroll 4
        for (int j = 0; j < 4; j++)
            #pragma unroll 4
            for (int k = 0; k < 4; k++)
                Gamma[i][j][k] = 0.0f;
    
    float x = pos.x;
    float y = pos.y;
    float z = pos.z;
    
    // Distance from bubble center
    float dx = x - xs;
    float dy = y - ys;
    float dz = z - zs;
    float rs = sqrtf(dx*dx + dy*dy + dz*dz);
    rs = fmaxf(rs, 1e-10f);
    
    // Shape function and derivative
    float f = alcubierreShapeFunction(rs, sigma, R);
    float df_drs = alcubierreShapeDerivative(rs, sigma, R);
    
    // Derivatives of rs with respect to position
    float drs[4] = {0.0f, dx / rs, dy / rs, dz / rs};
    
    // Derivatives of f with respect to position
    float df[4] = {0.0f, df_drs * drs[1], df_drs * drs[2], df_drs * drs[3]};
    
    float vsf = vs * f;
    float vsf2 = vsf * vsf;
    
    // Metric components
    // g_tt = -(1 - v²f²), g_tx = g_xt = -vf, g_xx = g_yy = g_zz = 1
    
    // Metric derivatives
    // ∂g_tt/∂x^i = 2v²f df/dx^i
    // ∂g_tx/∂x^i = -v df/dx^i
    
    float dg_tt[4], dg_tx[4];
    for (int i = 0; i < 4; i++) {
        dg_tt[i] = 2.0f * vs * vs * f * df[i];
        dg_tx[i] = -vs * df[i];
    }
    
    // Inverse metric (approximation for small perturbations)
    // g^tt ≈ -1 - v²f², g^tx ≈ -vf, g^xx = g^yy = g^zz = 1
    float g_inv_tt = -1.0f - vsf2;  // Approximate
    float g_inv_tx = -vsf;
    
    // Christoffel symbols from Γ^λ_μν = (1/2) g^λσ (∂_μ g_σν + ∂_ν g_σμ - ∂_σ g_μν)
    
    // Γ^t_tt = (1/2) g^tt ∂_t g_tt + (1/2) g^tx (2∂_t g_tx - ∂_x g_tt)
    // For static bubble: ∂_t terms = 0, so many simplify
    
    // Γ^t_xx = (1/2) g^tt (-∂_x g_tt) + (1/2) g^tx (-∂_x g_tx - ∂_x g_tx)
    //        = -(1/2) g^tt dg_tt[1] - g^tx dg_tx[1]
    Gamma[0][1][1] = -0.5f * g_inv_tt * dg_tt[1] - g_inv_tx * dg_tx[1];
    
    // Γ^t_tx = Γ^t_xt = (1/2) g^tt ∂_x g_tt
    Gamma[0][0][1] = 0.5f * g_inv_tt * dg_tt[1];
    Gamma[0][1][0] = Gamma[0][0][1];
    
    // Γ^x_tt = (1/2) g^xt ∂_t g_tt + (1/2) g^xx (2∂_t g_tx - ∂_x g_tt)
    //        = -(1/2) ∂_x g_tt (for static)
    Gamma[1][0][0] = -0.5f * dg_tt[1];
    
    // Γ^x_tx = Γ^x_xt = (1/2) g^xx ∂_x g_tt = (1/2) dg_tt[1]
    Gamma[1][0][1] = 0.5f * dg_tt[1];
    Gamma[1][1][0] = Gamma[1][0][1];
    
    // Similar for y and z components (using symmetry)
    for (int i = 2; i <= 3; i++) {
        Gamma[0][0][i] = 0.5f * g_inv_tt * dg_tt[i];
        Gamma[0][i][0] = Gamma[0][0][i];
        Gamma[0][1][i] = -g_inv_tx * dg_tx[i];
        Gamma[0][i][1] = Gamma[0][1][i];
        
        Gamma[1][0][i] = dg_tx[i];      // Γ^x_ti = ∂g_tx/∂x^i
        Gamma[1][i][0] = Gamma[1][0][i];
    }
}

// Numerical Metric Christoffel Symbols (Precomputed)
__device__ inline void getNumericalChristoffel(const Vec4& pos, const NumericalMetricData& nm, float Gamma[4][4][4]) {
    // Initialize to zero
    #pragma unroll 4
    for (int i = 0; i < 4; i++)
        #pragma unroll 4
        for (int j = 0; j < 4; j++)
            #pragma unroll 4
            for (int k = 0; k < 4; k++)
                Gamma[i][j][k] = 0.0f;

#ifdef __CUDACC__
    if (!nm.isLoaded || !nm.christoffelLoaded) return;

    // 1. Convert Spherical (r, θ, φ) to Cartesian (x, y, z) - Y-Polar
    float3 cart = make_float3(
        pos.r * sinf(pos.theta) * cosf(pos.phi),
        pos.r * cosf(pos.theta),
        pos.r * sinf(pos.theta) * sinf(pos.phi)
    );

    // 2. Compute texture coordinates
    float3 uvw = make_float3(
        (cart.x - nm.origin.x) / (nm.dims.x * nm.spacing.x),
        (cart.y - nm.origin.y) / (nm.dims.y * nm.spacing.y),
        (cart.z - nm.origin.z) / (nm.dims.z * nm.spacing.z)
    );

    // Check bounds
    if (uvw.x < 0.0f || uvw.x > 1.0f || 
        uvw.y < 0.0f || uvw.y > 1.0f || 
        uvw.z < 0.0f || uvw.z > 1.0f) {
        return;
    }

    // 3. Sample Christoffel symbols from textures
    // Indexing: Γ[mu][packed_idx] where packed_idx = ν + ρ*(ρ+1)/2 for ν ≤ ρ
    
    int idx[4][4];
    idx[0][0] = 0; idx[0][1] = 1; idx[0][2] = 2; idx[0][3] = 3;
    idx[1][0] = 1; idx[1][1] = 4; idx[1][2] = 5; idx[1][3] = 6;
    idx[2][0] = 2; idx[2][1] = 5; idx[2][2] = 7; idx[2][3] = 8;
    idx[3][0] = 3; idx[3][1] = 6; idx[3][2] = 8; idx[3][3] = 9;

    // Sample for each mu
    // mu=0 (t)
    for(int nu=0; nu<4; nu++) {
        for(int rho=nu; rho<4; rho++) {
            int p_idx = idx[nu][rho];
            float val = tex3D<float>(nm.Gamma_t[p_idx], uvw.x, uvw.y, uvw.z);
            Gamma[0][nu][rho] = val;
            if (nu != rho) Gamma[0][rho][nu] = val;
        }
    }
    // mu=1 (r/x)
    for(int nu=0; nu<4; nu++) {
        for(int rho=nu; rho<4; rho++) {
            int p_idx = idx[nu][rho];
            float val = tex3D<float>(nm.Gamma_r[p_idx], uvw.x, uvw.y, uvw.z);
            Gamma[1][nu][rho] = val;
            if (nu != rho) Gamma[1][rho][nu] = val;
        }
    }
    // mu=2 (theta/y)
    for(int nu=0; nu<4; nu++) {
        for(int rho=nu; rho<4; rho++) {
            int p_idx = idx[nu][rho];
            float val = tex3D<float>(nm.Gamma_theta[p_idx], uvw.x, uvw.y, uvw.z);
            Gamma[2][nu][rho] = val;
            if (nu != rho) Gamma[2][rho][nu] = val;
        }
    }
    // mu=3 (phi/z)
    for(int nu=0; nu<4; nu++) {
        for(int rho=nu; rho<4; rho++) {
            int p_idx = idx[nu][rho];
            float val = tex3D<float>(nm.Gamma_phi[p_idx], uvw.x, uvw.y, uvw.z);
            Gamma[3][nu][rho] = val;
            if (nu != rho) Gamma[3][rho][nu] = val;
        }
    }
#endif
}

//==============================================================================
// FIDO Basis Calculation (DNGR Paper Eq A.2-A.3)
//==============================================================================
template<Sirius::MetricType type>
__device__ Sirius::FIDOBasis computeFIDOBasis(
    const Vec4& x,
    const Sirius::MetricParams& mp)
{
    Sirius::FIDOBasis fido;
    
    float r = fmaxf(x.r, 0.01f);
    float theta = x.theta;
    // float phi = x.phi; // Unused in axisymmetric
    float M = mp.M;
    float a = mp.a;
    
    float r2 = r * r;
    float a2 = a * a;
    float cos_theta = cosf(theta);
    float sin_theta = sinf(theta);
    
    // Clamp sin_theta away from zero (poles)
    if (fabsf(sin_theta) < 0.001f) {
        sin_theta = (sin_theta >= 0.0f) ? 0.001f : -0.001f;
    }
    
    float cos2 = cos_theta * cos_theta;
    float sin2 = sin_theta * sin_theta;
    
    // Kerr metric quantities (DNGR Eq A.2)
    fido.rho = sqrtf(r2 + a2 * cos2);                              // ρ = √(r² + a²cos²θ)
    fido.Delta = r2 - 2.0f * M * r + a2;                            // Δ = r² - 2Mr + a²
    
    float r2_plus_a2 = r2 + a2;
    float Sigma_sq = r2_plus_a2 * r2_plus_a2 - a2 * fido.Delta * sin2;
    fido.Sigma = sqrtf(fmaxf(Sigma_sq, 0.01f));
    
    // Clamp Delta to avoid horizon singularity
    float Delta_safe = fmaxf(fido.Delta, 0.01f);
    
    // Lapse function: α = ρ√Δ/Σ
    fido.alpha = fido.rho * sqrtf(Delta_safe) / fido.Sigma;
    
    // Frame dragging angular velocity: ω = 2ar/Σ²
    fido.omega = 2.0f * a * r / (fido.Sigma * fido.Sigma);
    
    // Proper circumference factor: ϖ = Σsinθ/ρ
    fido.varpi = fido.Sigma * sin_theta / fido.rho;
    
    // Basis vectors (orthonormal basis vectors in coordinate basis)
    // We don't necessarily need to return vectors if we only use scalars, 
    // but the struct has float3 members. Let's fill them with placeholders
    // or Cartesian projections if needed. For now, unused by computeCameraSpeed.
    fido.e_r = make_float3(1.0f, 0.0f, 0.0f);
    fido.e_theta = make_float3(0.0f, 1.0f, 0.0f);
    fido.e_phi = make_float3(0.0f, 0.0f, 1.0f);
    
    return fido;
}

// Helper: Compute ISCO radius for Kerr metric (Bardeen et al. 1972)
__device__ inline float computeKerrISCO(float M, float a) {
    if (M < 1e-6f) return 0.0f;
    float a_star = a / M;
    a_star = (a_star > 0.999f) ? 0.999f : ((a_star < -0.999f) ? -0.999f : a_star);
    
    // Z1 = 1 + (1 - a^2)^(1/3) * ((1 + a)^(1/3) + (1 - a)^(1/3))
    float term1 = cbrtf(1.0f - a_star*a_star);
    float term2 = cbrtf(1.0f + a_star) + cbrtf(1.0f - a_star);
    float Z1 = 1.0f + term1 * term2;
    
    // Z2 = sqrt(3 a^2 + Z1^2)
    float Z2 = sqrtf(3.0f * a_star*a_star + Z1*Z1);
    
    // r_ISCO = M * (3 + Z2 - sign(a) * sqrt((3-Z1)(3+Z1+2Z2)))
    float term_sqrt = sqrtf(fmaxf((3.0f - Z1) * (3.0f + Z1 + 2.0f * Z2), 0.0f));
    float sign_a = (a_star >= 0.0f) ? 1.0f : -1.0f;
    
    return M * (3.0f + Z2 - sign_a * term_sqrt);
}

//==============================================================================
// ACCRETION DISK VOLUMETRIC SAMPLING (Accurate Radiative Transfer)
//==============================================================================
__device__ inline void sampleAccretionDiskVolumetric(
    const Vec4& pos,
    const Sirius::AccretionDiskParams& disk,
    const Sirius::MetricParams& mp,
    float stepSize,
    float beamRadius,
    float3& emission,
    float& dTau,
    float3& velocity)
{
    emission = make_float3(0.0f, 0.0f, 0.0f);
    dTau = 0.0f;
    velocity = make_float3(0.0f, 0.0f, 0.0f);
    
    if (!disk.enabled) return;

    // Convert to disk coordinates (cylindrical r, z)
    float r = pos.r;
    float theta = pos.theta;
    
    // Vertical distance from equatorial plane
    // z = r * cos(theta)
    float z = r * cosf(theta);
    float z_abs = fabsf(z);
    
    // Cylindrical radius (projected)
    // R = r * sin(theta)
    float R = r * sinf(theta);
    
    // 1. Check Radial Bounds
    float r_isco = disk.innerRadius;
    if (r_isco < 0.001f) {
         // Calculate exact ISCO for given spin
         r_isco = computeKerrISCO(mp.M, mp.a);
    }
    
    if (R < r_isco || R > disk.outerRadius) return;
    
    // 2. Check Vertical Bounds (Disk Height)
    float H = disk.heightScale * R;
    if (z_abs > 3.0f * H) return;
    
    // 3. Compute Density
    // Vertical profile: Gaussian
    float rho_vertical = expf(-0.5f * (z * z) / (H * H));
    
    // Radial tapering
    // Smoothstep implementation
    float inner_fade = 0.0f;
    float outer_fade = 0.0f;
    
    // Inner fade
    {
        float t = (R - r_isco) / (r_isco * 0.1f);
        t = (t < 0.0f) ? 0.0f : ((t > 1.0f) ? 1.0f : t);
        inner_fade = t * t * (3.0f - 2.0f * t);
    }
    
    // Outer fade
    {
         float t = (R - disk.outerRadius * 0.9f) / (disk.outerRadius * 0.1f);
         t = (t < 0.0f) ? 0.0f : ((t > 1.0f) ? 1.0f : t);
         outer_fade = 1.0f - t * t * (3.0f - 2.0f * t);
    }

    
    float rho_radial = inner_fade * outer_fade;
    float density = rho_vertical * rho_radial;
    
    if (density < 1e-4f) return;
    
    // 4. Compute Temperature and Emission
    float T = disk.innerTemperature * powf(r_isco / R, disk.temperatureExponent);
    
    // Novikov-Thorne correction Q(r)
    if (disk.temperatureModel == Sirius::TemperatureModel::NovikovThorne) {
         float x_val = sqrtf(R / mp.M);
         float x0 = sqrtf(r_isco / mp.M);
         if (x_val > x0) {
             float Q = (1.0f - x0/x_val); // Rough approx
             T *= powf(Q, 0.25f);
         } else {
             T = 0.0f;
         }
    }
    
    float3 color;
    if (disk.useSpectralColors) {
         // Simple blackbody-like ramp
         if (T > 10000.0f) color = make_float3(0.8f, 0.9f, 1.0f);
         else if (T > 5000.0f) color = make_float3(1.0f, 0.9f, 0.7f);
         else if (T > 2000.0f) color = make_float3(1.0f, 0.5f, 0.1f);
         else color = make_float3(0.8f, 0.1f, 0.0f);
    } else {
         color = make_float3(1.0f, 1.0f, 1.0f);
    }
    
    emission = color * density * disk.emissionCoefficient;
    dTau = density * disk.absorptionCoefficient * stepSize;
    
    // Velocity is handled by integrator's computeKeplerianVelocity usually
}



