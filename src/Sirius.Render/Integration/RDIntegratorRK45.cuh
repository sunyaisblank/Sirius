// RDIntegratorRK45.cuh - General Purpose RK45 Integrator
// Component ID: RDIN003A (Integration/RK45)
//
// Implements adaptive Runge-Kutta-Fehlberg (RK45) integration for
// general relativistic geodesics in arbitrary spacetimes.
// Supports both Spherical and Cartesian integration paths depending on metric.

#pragma once

#include "../Integration/RDMath.cuh"
#include "../Transport/RDMetricFunctions.cuh"

namespace Sirius {


//==============================================================================
// Coordinate Conversions (State Level)
//==============================================================================

__device__ inline GeodesicStateCart sphericalToCartesianState(const GeodesicState& s) {
    GeodesicStateCart cart;
    cart.x = vec4SphToCart(s.x);
    
    // Convert velocity using Jacobian
    float r = s.x.r;
    float theta = s.x.theta;
    float phi = s.x.phi;
    
    float st, ct; sincosf(theta, &st, &ct);
    float sp, cp; sincosf(phi, &sp, &cp);
    
    float dr = s.u.r;
    float dt = s.u.theta;
    float dp = s.u.phi;
    
    // t dot
    cart.u.t = s.u.t;
    
    // x = r st cp
    // dx = dr st cp + r ct cp dtheta - r st sp dphi
    cart.u.x = dr * st * cp + r * ct * cp * dt - r * st * sp * dp;
    
    // y = r ct (Y-Polar convention?) NO, sticking to Z-Polar for Physics Kernels
    // If using Z-Polar (z = r cos theta):
    // y = r st sp
    // dy = dr st sp + r ct sp dtheta + r st cp dphi
    cart.u.y = dr * st * sp + r * ct * sp * dt + r * st * cp * dp;
    
    // z = r ct
    // dz = dr ct - r st dtheta
    cart.u.z = dr * ct - r * st * dt;
    
    return cart;
}

__device__ inline GeodesicState cartesianToSphericalState(const GeodesicStateCart& c) {
    GeodesicState sph;
    sph.x = vec4CartToSph(c.x);
    
    // Convert velocity using inverse Jacobian
    // r = sqrt(x^2 + y^2 + z^2)
    float r = sph.x.r;
    float r2 = r * r;
    float rho = sqrtf(c.x.x*c.x.x + c.x.y*c.x.y); // xy plane distance
    
    if (r < 1e-10f) {
        sph.u = Vec4(c.u.t, 0.0f, 0.0f, 0.0f);
        return sph;
    }
    
    sph.u.t = c.u.t;
    
    // dr = (x dx + y dy + z dz) / r
    sph.u.r = (c.x.x * c.u.x + c.x.y * c.u.y + c.x.z * c.u.z) / r;
    
    // dtheta: theta = acos(z/r)
    // dtheta = -1/sqrt(1-(z/r)^2) * (dz*r - z*dr)/r^2
    //        = -r/rho * (dz*r - z*dr)/r^2
    //        = (z*dr - r*dz) / (r * rho)
    if (rho > 1e-10f) {
        sph.u.theta = (c.x.z * sph.u.r - r * c.u.z) / (r * rho);
    } else {
        sph.u.theta = 0.0f;
    }
    
    // dphi: phi = atan2(y, x)
    // dphi = (x dy - y dx) / (x^2 + y^2)
    //      = (x dy - y dx) / rho^2
    if (rho > 1e-10f) {
        sph.u.phi = (c.x.x * c.u.y - c.x.y * c.u.x) / (rho * rho);
    } else {
        sph.u.phi = 0.0f;
    }
    
    return sph;
}

//==============================================================================
// Christoffel Dispatch (General)
//==============================================================================

__device__ inline void getGeneralChristoffel(
    const Vec4& pos_sph, 
    const Vec4Cart& pos_cart,
    const LaunchParams& params,
    float Gamma[4][4][4],
    bool useCartesian) 
{
    // 1. Numerical Metric (Overrides if active)
    if (params.numericalMetric.isLoaded && params.numericalMetric.christoffelLoaded) {
        getNumericalChristoffel(pos_sph, params.numericalMetric, Gamma);
        // Note: Numerical Christoffels are typically in Coordinate Basis of the Grid
        // If Grid is Cartesian, these are Cartesian Gammas.
        // If integrating in Cartesian, this is fine.
        return;
    }

    // 2. Analytic Families
    switch(params.metricParams.family) {
        case MetricFamily::KerrSchild:
            // Always Cartesian
            getKerrSchildFamilyChristoffel(pos_cart, 
                params.metricParams.kerrSchild.M,
                params.metricParams.kerrSchild.a,
                params.metricParams.kerrSchild.Q,
                params.metricParams.kerrSchild.Lambda,
                Gamma);
            break;
            
        case MetricFamily::WarpDrive:
            // Always Cartesian
            getWarpDriveChristoffel(pos_cart,
                params.metricParams.warpDrive.vs,
                params.metricParams.warpDrive.sigma,
                params.metricParams.warpDrive.R,
                0.0f, 0.0f, 0.0f, // xs, ys, zs assumed 0 relative to current frame or handled by conversion
                Gamma);
            break;
            
        case MetricFamily::MorrisThorne:
            // Always Spherical
            getMorrisThorneChristoffel(pos_sph,
                params.metricParams.morrisThorne.b0,
                params.metricParams.morrisThorne.redshiftPhi,
                WormholeShape::Ellis, // Parameterize this?
                Gamma);
            break;
            
        default:
            // Minkowski Fallback (Zero Gamma)
            for(int i=0; i<4; i++)
                for(int j=0; j<4; j++)
                    for(int k=0; k<4; k++)
                        Gamma[i][j][k] = 0.0f;
            break;
    }
}

//==============================================================================
// Geodesic Equation (Spherical)
//==============================================================================
__device__ inline Vec4 evaluateGeodesicDerivSpherical(
    const Vec4& x, const Vec4& u, const float Gamma[4][4][4]) 
{
    float du[4] = {0.0f};
    float u_arr[4] = {u.t, u.r, u.theta, u.phi};
    
    // du^lambda/dtau = -Gamma^lambda_mu_nu * u^mu * u^nu
    for(int lam=0; lam<4; lam++) {
        float sum = 0.0f;
        for(int mu=0; mu<4; mu++) {
            float gamma_u_sum = 0.0f;
            for(int nu=0; nu<4; nu++) {
                gamma_u_sum += Gamma[lam][mu][nu] * u_arr[nu];
            }
            sum += gamma_u_sum * u_arr[mu];
        }
        du[lam] = -sum;
    }
    
    return Vec4(du[0], du[1], du[2], du[3]);
}


//==============================================================================
// Geodesic Equation (Cartesian)
//==============================================================================
__device__ inline Vec4Cart evaluateGeodesicDerivCartesian(
    const Vec4Cart& x, const Vec4Cart& u, const float Gamma[4][4][4]) 
{
    float du[4] = {0.0f};
    float u_arr[4] = {u.t, u.x, u.y, u.z};
    
    // du^lambda/dtau = -Gamma^lambda_mu_nu * u^mu * u^nu
    for(int lam=0; lam<4; lam++) {
        float sum = 0.0f;
        for(int mu=0; mu<4; mu++) {
            float gamma_u_sum = 0.0f;
            for(int nu=0; nu<4; nu++) {
                gamma_u_sum += Gamma[lam][mu][nu] * u_arr[nu];
            }
            sum += gamma_u_sum * u_arr[mu];
        }
        du[lam] = -sum;
    }
    
    return Vec4Cart(du[0], du[1], du[2], du[3]);
}

//==============================================================================
// RK45 Step (Spherical)
//==============================================================================
__device__ inline void rk45_step_spherical(
    GeodesicState& state, 
    float& h, 
    const LaunchParams& params,
    float tolerance,
    bool& stepAccepted) 
{
    float Gamma[4][4][4];
    
    // k1
    getGeneralChristoffel(state.x, Vec4Cart(), params, Gamma, false);
    Vec4 k1_x = state.u;
    Vec4 k1_u = evaluateGeodesicDerivSpherical(state.x, state.u, Gamma);
    
    // k2 (at x + h/4)
    GeodesicState s2;
    s2.x = state.x + k1_x * (h * 0.25f);
    s2.u = state.u + k1_u * (h * 0.25f);
    getGeneralChristoffel(s2.x, Vec4Cart(), params, Gamma, false); // Re-eval Gamma? Yes.
    Vec4 k2_x = s2.u;
    Vec4 k2_u = evaluateGeodesicDerivSpherical(s2.x, s2.u, Gamma);
    
    // k3 (at x + 3h/8)
    GeodesicState s3;
    s3.x = state.x + k1_x * (h * 3.0f/32.0f) + k2_x * (h * 9.0f/32.0f);
    s3.u = state.u + k1_u * (h * 3.0f/32.0f) + k2_u * (h * 9.0f/32.0f);
    getGeneralChristoffel(s3.x, Vec4Cart(), params, Gamma, false);
    Vec4 k3_x = s3.u;
    Vec4 k3_u = evaluateGeodesicDerivSpherical(s3.x, s3.u, Gamma);
    
    // k4 (at x + 12h/13)
    GeodesicState s4;
    s4.x = state.x + k1_x * (h * 1932.0f/2197.0f) - k2_x * (h * 7200.0f/2197.0f) + k3_x * (h * 7296.0f/2197.0f);
    s4.u = state.u + k1_u * (h * 1932.0f/2197.0f) - k2_u * (h * 7200.0f/2197.0f) + k3_u * (h * 7296.0f/2197.0f);
    getGeneralChristoffel(s4.x, Vec4Cart(), params, Gamma, false);
    Vec4 k4_x = s4.u;
    Vec4 k4_u = evaluateGeodesicDerivSpherical(s4.x, s4.u, Gamma);
    
    // k5 (at x + h)
    GeodesicState s5;
    s5.x = state.x + k1_x * (h * 439.0f/216.0f) - k2_x * (h * 8.0f) + k3_x * (h * 3680.0f/513.0f) - k4_x * (h * 845.0f/4104.0f);
    s5.u = state.u + k1_u * (h * 439.0f/216.0f) - k2_u * (h * 8.0f) + k3_u * (h * 3680.0f/513.0f) - k4_u * (h * 845.0f/4104.0f);
    getGeneralChristoffel(s5.x, Vec4Cart(), params, Gamma, false);
    Vec4 k5_x = s5.u;
    Vec4 k5_u = evaluateGeodesicDerivSpherical(s5.x, s5.u, Gamma);
    
    // k6 (at x + h/2)
    GeodesicState s6;
    s6.x = state.x - k1_x * (h * 8.0f/27.0f) + k2_x * (h * 2.0f) - k3_x * (h * 3544.0f/2565.0f) + k4_x * (h * 1859.0f/4104.0f) - k5_x * (h * 11.0f/40.0f);
    s6.u = state.u - k1_u * (h * 8.0f/27.0f) + k2_u * (h * 2.0f) - k3_u * (h * 3544.0f/2565.0f) + k4_u * (h * 1859.0f/4104.0f) - k5_u * (h * 11.0f/40.0f);
    getGeneralChristoffel(s6.x, Vec4Cart(), params, Gamma, false);
    Vec4 k6_x = s6.u;
    Vec4 k6_u = evaluateGeodesicDerivSpherical(s6.x, s6.u, Gamma);
    
    // 5th order solution
    Vec4 dx5 = k1_x * (16.0f/135.0f) + k3_x * (6656.0f/12825.0f) + k4_x * (28561.0f/56430.0f) - k5_x * (9.0f/50.0f) + k6_x * (2.0f/55.0f);
    Vec4 du5 = k1_u * (16.0f/135.0f) + k3_u * (6656.0f/12825.0f) + k4_u * (28561.0f/56430.0f) - k5_u * (9.0f/50.0f) + k6_u * (2.0f/55.0f);
    
    // 4th order solution (for error est)
    Vec4 dx4 = k1_x * (25.0f/216.0f) + k3_x * (1408.0f/2565.0f) + k4_x * (2197.0f/4104.0f) - k5_x * (1.0f/5.0f);
    Vec4 du4 = k1_u * (25.0f/216.0f) + k3_u * (1408.0f/2565.0f) + k4_u * (2197.0f/4104.0f) - k5_u * (1.0f/5.0f);
    
    // Error estimate
    // Simple max norms
    float error = 0.0f;
    error = fmaxf(error, fabsf(dx5.r - dx4.r));
    error = fmaxf(error, fabsf(dx5.theta - dx4.theta));
    error = fmaxf(error, fabsf(dx5.phi - dx4.phi));
    error = fmaxf(error, fabsf(du5.t - du4.t));
    
    if (error > tolerance) {
        stepAccepted = false;
    } else {
        stepAccepted = true;
        state.x = state.x + dx5 * h;
        state.u = state.u + du5 * h;
    }
}

//==============================================================================
// RK45 Step (Cartesian)
//==============================================================================
__device__ inline void rk45_step_cartesian(
    GeodesicStateCart& state, 
    float& h, 
    const LaunchParams& params,
    float tolerance,
    bool& stepAccepted) 
{
    float Gamma[4][4][4];
    
    // k1
    getGeneralChristoffel(Vec4(), state.x, params, Gamma, true);
    Vec4Cart k1_x = state.u;
    Vec4Cart k1_u = evaluateGeodesicDerivCartesian(state.x, state.u, Gamma);
    
    // k2
    GeodesicStateCart s2;
    s2.x = state.x + k1_x * (h * 0.25f);
    s2.u = state.u + k1_u * (h * 0.25f);
    getGeneralChristoffel(Vec4(), s2.x, params, Gamma, true);
    Vec4Cart k2_x = s2.u;
    Vec4Cart k2_u = evaluateGeodesicDerivCartesian(s2.x, s2.u, Gamma);
    
    // ... Assume simplified RK4 for brevity if needed, but sticking to RK45 structure ...
    // Using RK4 for Cartesian to save typing, as this is already a fallback?
    // No, be consistent.
    
    // k3
    GeodesicStateCart s3;
    s3.x = state.x + k1_x * (h * 3.0f/32.0f) + k2_x * (h * 9.0f/32.0f);
    s3.u = state.u + k1_u * (h * 3.0f/32.0f) + k2_u * (h * 9.0f/32.0f);
    getGeneralChristoffel(Vec4(), s3.x, params, Gamma, true);
    Vec4Cart k3_x = s3.u;
    Vec4Cart k3_u = evaluateGeodesicDerivCartesian(s3.x, s3.u, Gamma);
    
    // k4
    GeodesicStateCart s4;
    s4.x = state.x + k1_x * (h * 1932.0f/2197.0f) - k2_x * (h * 7200.0f/2197.0f) + k3_x * (h * 7296.0f/2197.0f);
    s4.u = state.u + k1_u * (h * 1932.0f/2197.0f) - k2_u * (h * 7200.0f/2197.0f) + k3_u * (h * 7296.0f/2197.0f);
    getGeneralChristoffel(Vec4(), s4.x, params, Gamma, true);
    Vec4Cart k4_x = s4.u;
    Vec4Cart k4_u = evaluateGeodesicDerivCartesian(s4.x, s4.u, Gamma);
    
    // k5
    GeodesicStateCart s5;
    s5.x = state.x + k1_x * (h * 439.0f/216.0f) - k2_x * (h * 8.0f) + k3_x * (h * 3680.0f/513.0f) - k4_x * (h * 845.0f/4104.0f);
    s5.u = state.u + k1_u * (h * 439.0f/216.0f) - k2_u * (h * 8.0f) + k3_u * (h * 3680.0f/513.0f) - k4_u * (h * 845.0f/4104.0f);
    getGeneralChristoffel(Vec4(), s5.x, params, Gamma, true);
    Vec4Cart k5_x = s5.u;
    Vec4Cart k5_u = evaluateGeodesicDerivCartesian(s5.x, s5.u, Gamma);
    
    // k6
    GeodesicStateCart s6;
    s6.x = state.x - k1_x * (h * 8.0f/27.0f) + k2_x * (h * 2.0f) - k3_x * (h * 3544.0f/2565.0f) + k4_x * (h * 1859.0f/4104.0f) - k5_x * (h * 11.0f/40.0f);
    s6.u = state.u - k1_u * (h * 8.0f/27.0f) + k2_u * (h * 2.0f) - k3_u * (h * 3544.0f/2565.0f) + k4_u * (h * 1859.0f/4104.0f) - k5_u * (h * 11.0f/40.0f);
    getGeneralChristoffel(Vec4(), s6.x, params, Gamma, true);
    Vec4Cart k6_x = s6.u;
    Vec4Cart k6_u = evaluateGeodesicDerivCartesian(s6.x, s6.u, Gamma);
    
    // 5th order solution
    Vec4Cart dx5 = k1_x * (16.0f/135.0f) + k3_x * (6656.0f/12825.0f) + k4_x * (28561.0f/56430.0f) - k5_x * (9.0f/50.0f) + k6_x * (2.0f/55.0f);
    Vec4Cart du5 = k1_u * (16.0f/135.0f) + k3_u * (6656.0f/12825.0f) + k4_u * (28561.0f/56430.0f) - k5_u * (9.0f/50.0f) + k6_u * (2.0f/55.0f);
    
    // 4th order solution (for error est)
    Vec4Cart dx4 = k1_x * (25.0f/216.0f) + k3_x * (1408.0f/2565.0f) + k4_x * (2197.0f/4104.0f) - k5_x * (1.0f/5.0f);
    Vec4Cart du4 = k1_u * (25.0f/216.0f) + k3_u * (1408.0f/2565.0f) + k4_u * (2197.0f/4104.0f) - k5_u * (1.0f/5.0f);
    
    // Error estimate
    float error = 0.0f;
    error = fmaxf(error, fabsf(dx5.x - dx4.x));
    error = fmaxf(error, fabsf(dx5.y - dx4.y));
    error = fmaxf(error, fabsf(dx5.z - dx4.z));
    error = fmaxf(error, fabsf(du5.t - du4.t));
    
    if (error > tolerance) {
        stepAccepted = false;
    } else {
        stepAccepted = true;
        state.x = state.x + dx5 * h;
        state.u = state.u + du5 * h;
    }
}


//==============================================================================
// Main Interface
//==============================================================================
__device__ inline GeodesicState integrateGeodesicRK45(
    const GeodesicState& state,
    float& h,
    const LaunchParams& params,
    float tolerance,
    float minStep,
    float maxStep,
    bool& stepAccepted)
{
    // Determine integration strategy
    // Numerical metrics use Cartesian grid logic -> integrate in Cartesian
    // Kerr-Schild and Warp Drive are defined in Cartesian -> integrate in Cartesian
    // Morris-Thorne is defined in Spherical -> integrate in Spherical
    
    bool useCartesian = false;
    
    if (params.numericalMetric.isLoaded) {
        useCartesian = true;
    } else if (params.metricParams.family == MetricFamily::KerrSchild ||
               params.metricParams.family == MetricFamily::WarpDrive) {
        useCartesian = true;
    }
    
    if (useCartesian) {
        GeodesicStateCart cart = sphericalToCartesianState(state);
        rk45_step_cartesian(cart, h, params, tolerance, stepAccepted);
        if (stepAccepted) {
            // Adaptive step size increase logic can be here
            return cartesianToSphericalState(cart);
        } else {
            h = fmaxf(h * 0.5f, minStep);
            return state;
        }
    } else {
        GeodesicState sph = state;
        rk45_step_spherical(sph, h, params, tolerance, stepAccepted);
        if (stepAccepted) {
            return sph;
        } else {
            h = fmaxf(h * 0.5f, minStep);
            return state;
        }
    }
}

} // namespace Sirius
