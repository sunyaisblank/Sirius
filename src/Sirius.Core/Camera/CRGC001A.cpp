// CRGC001A.cpp - Geodesic Camera Controller Implementation

#include "CRGC001A.h"

namespace Sirius {

//==============================================================================
// Constructor
//==============================================================================
GeodesicCamera::GeodesicCamera()
    : m_Mode(CameraMode::FreeWASD)
    , m_M(1.0f)
    , m_a(0.0f)
    , m_OrbitRadius(10.0f)
    , m_Prograde(true)
    , m_StepSize(0.01f)  // Initial adaptive step size
{
    // Default state: at r=50M, equatorial plane, at rest
    m_State.position = Vec4(0.0f, 50.0f, 1.5708f, 0.0f);
    m_State.velocity = Vec4(1.0f, 0.0f, 0.0f, 0.0f);
    normalizeVelocity();
}

//==============================================================================
// Initialization
//==============================================================================
void GeodesicCamera::initialize(const Vec4& position, float M, float a) {
    m_State.position = position;
    m_M = M;
    m_a = a;
    
    // Initialize at rest relative to spatial infinity
    m_State.velocity = Vec4(1.0f, 0.0f, 0.0f, 0.0f);
    normalizeVelocity();
    
    m_State.properTime = 0.0f;
    m_State.coordinateTime = 0.0f;
}

void GeodesicCamera::setMetricParams(float M, float a) {
    m_M = M;
    m_a = a;
}

//==============================================================================
// Orbit Presets
//==============================================================================
void GeodesicCamera::setCircularOrbit(float r, bool prograde) {
    m_OrbitRadius = r;
    m_Prograde = prograde;
    m_Mode = CameraMode::CircularOrbit;
    
    // Set position at equatorial plane
    m_State.position.t = 0.0f;
    m_State.position.r = r * m_M;  // r in units of M
    m_State.position.theta = 1.5707963f;  // π/2
    m_State.position.phi = 0.0f;
    
    // Compute circular orbit velocity
    m_State.velocity = computeCircularVelocity(r, prograde);
    normalizeVelocity();
    
    m_State.properTime = 0.0f;
    m_State.coordinateTime = 0.0f;
}

void GeodesicCamera::setISCO(bool prograde) {
    float r_isco = computeISCO(prograde);
    setCircularOrbit(r_isco / m_M, prograde);  // Convert back to units of M
    m_Mode = CameraMode::ISCO;
}

void GeodesicCamera::setPlungingOrbit(float r_start) {
    m_Mode = CameraMode::PlungingOrbit;
    
    // Start at given radius, equatorial plane
    m_State.position.t = 0.0f;
    m_State.position.r = r_start * m_M;
    m_State.position.theta = 1.5707963f;  // π/2
    m_State.position.phi = 0.0f;
    
    // Plunging from rest at infinity has E = 1, L = 0
    // For finite distance, we compute the appropriate radial velocity
    float r = r_start * m_M;
    float rs = 2.0f * m_M;
    
    // Energy at infinity E = 1, angular momentum L = 0
    // Radial velocity from energy conservation for Schwarzschild:
    // (dr/dτ)² = E² - (1 - rs/r)
    float energy_sq = 1.0f;  // E = 1 (dropped from rest at infinity)
    float effective_potential = 1.0f - rs / r;
    float dr_dtau_sq = energy_sq - effective_potential;
    
    if (dr_dtau_sq > 0.0f) {
        m_State.velocity.r = -sqrtf(dr_dtau_sq);  // Negative = falling inward
    } else {
        m_State.velocity.r = 0.0f;
    }
    
    // dt/dτ from energy: E = -(g_tt) * (dt/dτ) for radial motion
    // For Schwarzschild: dt/dτ = E / (1 - rs/r)
    m_State.velocity.t = 1.0f / (1.0f - rs / r);
    
    m_State.velocity.theta = 0.0f;
    m_State.velocity.phi = 0.0f;
    
    normalizeVelocity();
    
    m_State.properTime = 0.0f;
    m_State.coordinateTime = 0.0f;
}

void GeodesicCamera::setMode(CameraMode mode) {
    m_Mode = mode;
}

//==============================================================================
// Time Evolution
//==============================================================================
void GeodesicCamera::advanceByCoordinateTime(float dt) {
    if (m_Mode == CameraMode::FreeWASD) {
        return;  // Non-geodesic mode - handled elsewhere
    }
    
    // Accumulate coordinate time
    m_State.coordinateTime += dt;
    
    // Integrate geodesic equation using RK45/DOPRI adaptive step
    // The step size adapts automatically based on local truncation error
    float remaining = dt;
    while (remaining > 1e-10f) {
        float step = fminf(m_StepSize, remaining);
        rk45Step(step);
        remaining -= step;
        m_StepSize = step;  // Save adapted step for next frame
    }
    
    // Update proper time: dτ = dt * (dτ/dt)
    float time_dilation = getTimeDilationFactor();
    m_State.properTime += dt * time_dilation;
}

void GeodesicCamera::advanceByProperTime(float dtau) {
    if (m_Mode == CameraMode::FreeWASD) {
        return;
    }
    
    // Convert proper time to coordinate time: dt = dτ / (dτ/dt)
    float time_dilation = getTimeDilationFactor();
    float dt = (time_dilation > 0.001f) ? dtau / time_dilation : dtau;
    
    advanceByCoordinateTime(dt);
}

void GeodesicCamera::resetProperTime() {
    m_State.properTime = 0.0f;
    m_State.coordinateTime = 0.0f;
}

//==============================================================================
// Time Dilation
//==============================================================================
float GeodesicCamera::getTimeDilationFactor() const {
    // dτ/dt = sqrt(1 - rs/r - (r/c * dφ/dt)² + ...) simplified for Schwarzschild
    float r = m_State.position.r;
    float rs = 2.0f * m_M;
    
    // Basic gravitational time dilation
    float grav_factor = sqrtf(fmaxf(1.0f - rs / r, 0.001f));
    
    // For circular orbit, include velocity contribution
    if (m_Mode == CameraMode::CircularOrbit || m_Mode == CameraMode::ISCO) {
        float omega = getOrbitalAngularVelocity();
        float v_orbital = omega * r;  // Approximate orbital velocity
        float lorentz = sqrtf(fmaxf(1.0f - v_orbital * v_orbital, 0.001f));
        return grav_factor * lorentz;
    }
    
    return grav_factor;
}

float GeodesicCamera::getOrbitalAngularVelocity() const {
    // Keplerian angular velocity: Ω = sqrt(M) / (r^1.5 + a*sqrt(M))
    float r = m_State.position.r;
    float omega = sqrtf(m_M) / (powf(r, 1.5f) + m_a * sqrtf(m_M));
    return m_Prograde ? omega : -omega;
}

//==============================================================================
// Coordinate Conversion
//==============================================================================
Vec3 GeodesicCamera::getCartesianPosition() const {
    float r = m_State.position.r;
    float theta = m_State.position.theta;
    float phi = m_State.position.phi;
    
    // Boyer-Lindquist to Cartesian (for Kerr, approximately)
    float x = sqrtf(r * r + m_a * m_a) * sinf(theta) * cosf(phi);
    float y = sqrtf(r * r + m_a * m_a) * sinf(theta) * sinf(phi);
    float z = r * cosf(theta);
    
    return Vec3(x, y, z);
}

Vec3 GeodesicCamera::getLookDirection() const {
    // Default: look toward the black hole (origin)
    Vec3 pos = getCartesianPosition();
    float mag = sqrtf(pos.x * pos.x + pos.y * pos.y + pos.z * pos.z);
    if (mag > 0.001f) {
        return Vec3(-pos.x / mag, -pos.y / mag, -pos.z / mag);
    }
    return Vec3(0.0f, 0.0f, 1.0f);
}

//==============================================================================
// Geodesic Integration - RK45/DOPRI (Dormand-Prince) Adaptive Step
//==============================================================================

// DOPRI Butcher tableau coefficients
static constexpr float DOPRI_A21 = 1.0f/5.0f;
static constexpr float DOPRI_A31 = 3.0f/40.0f, DOPRI_A32 = 9.0f/40.0f;
static constexpr float DOPRI_A41 = 44.0f/45.0f, DOPRI_A42 = -56.0f/15.0f, DOPRI_A43 = 32.0f/9.0f;
static constexpr float DOPRI_A51 = 19372.0f/6561.0f, DOPRI_A52 = -25360.0f/2187.0f;
static constexpr float DOPRI_A53 = 64448.0f/6561.0f, DOPRI_A54 = -212.0f/729.0f;
static constexpr float DOPRI_A61 = 9017.0f/3168.0f, DOPRI_A62 = -355.0f/33.0f;
static constexpr float DOPRI_A63 = 46732.0f/5247.0f, DOPRI_A64 = 49.0f/176.0f, DOPRI_A65 = -5103.0f/18656.0f;
static constexpr float DOPRI_A71 = 35.0f/384.0f, DOPRI_A73 = 500.0f/1113.0f;
static constexpr float DOPRI_A74 = 125.0f/192.0f, DOPRI_A75 = -2187.0f/6784.0f, DOPRI_A76 = 11.0f/84.0f;

// 5th order solution coefficients (same as A7*)
static constexpr float DOPRI_B1 = 35.0f/384.0f, DOPRI_B3 = 500.0f/1113.0f;
static constexpr float DOPRI_B4 = 125.0f/192.0f, DOPRI_B5 = -2187.0f/6784.0f, DOPRI_B6 = 11.0f/84.0f;

// 4th order solution coefficients (for error estimation)
static constexpr float DOPRI_Bh1 = 5179.0f/57600.0f, DOPRI_Bh3 = 7571.0f/16695.0f;
static constexpr float DOPRI_Bh4 = 393.0f/640.0f, DOPRI_Bh5 = -92097.0f/339200.0f;
static constexpr float DOPRI_Bh6 = 187.0f/2100.0f, DOPRI_Bh7 = 1.0f/40.0f;

// Error coefficients (B - Bh)
static constexpr float DOPRI_E1 = DOPRI_B1 - DOPRI_Bh1;
static constexpr float DOPRI_E3 = DOPRI_B3 - DOPRI_Bh3;
static constexpr float DOPRI_E4 = DOPRI_B4 - DOPRI_Bh4;
static constexpr float DOPRI_E5 = DOPRI_B5 - DOPRI_Bh5;
static constexpr float DOPRI_E6 = DOPRI_B6 - DOPRI_Bh6;

void GeodesicCamera::rk45Step(float& dt) {
    const float TOL = 1.0e-6f, SAFETY = 0.9f, MIN_STEP = 1.0e-8f, MAX_STEP = 1.0f;
    
    Vec4 pos = m_State.position;
    Vec4 vel = m_State.velocity;
    
    bool accepted = false;
    for (int iter = 0; iter < 20 && !accepted; ++iter) {
        // Compute RK stages k1-k6
        Vec4 k1_pos = vel, k1_vel = computeAcceleration(pos, vel);
        
        Vec4 pos2 = pos + k1_pos * (dt * DOPRI_A21);
        Vec4 vel2 = vel + k1_vel * (dt * DOPRI_A21);
        Vec4 k2_pos = vel2, k2_vel = computeAcceleration(pos2, vel2);
        
        Vec4 pos3 = pos + k1_pos * (dt * DOPRI_A31) + k2_pos * (dt * DOPRI_A32);
        Vec4 vel3 = vel + k1_vel * (dt * DOPRI_A31) + k2_vel * (dt * DOPRI_A32);
        Vec4 k3_pos = vel3, k3_vel = computeAcceleration(pos3, vel3);
        
        Vec4 pos4 = pos + k1_pos * (dt * DOPRI_A41) + k2_pos * (dt * DOPRI_A42) + k3_pos * (dt * DOPRI_A43);
        Vec4 vel4 = vel + k1_vel * (dt * DOPRI_A41) + k2_vel * (dt * DOPRI_A42) + k3_vel * (dt * DOPRI_A43);
        Vec4 k4_pos = vel4, k4_vel = computeAcceleration(pos4, vel4);
        
        Vec4 pos5 = pos + k1_pos * (dt * DOPRI_A51) + k2_pos * (dt * DOPRI_A52) 
                       + k3_pos * (dt * DOPRI_A53) + k4_pos * (dt * DOPRI_A54);
        Vec4 vel5 = vel + k1_vel * (dt * DOPRI_A51) + k2_vel * (dt * DOPRI_A52) 
                       + k3_vel * (dt * DOPRI_A53) + k4_vel * (dt * DOPRI_A54);
        Vec4 k5_pos = vel5, k5_vel = computeAcceleration(pos5, vel5);
        
        Vec4 pos6 = pos + k1_pos * (dt * DOPRI_A61) + k2_pos * (dt * DOPRI_A62)
                       + k3_pos * (dt * DOPRI_A63) + k4_pos * (dt * DOPRI_A64) + k5_pos * (dt * DOPRI_A65);
        Vec4 vel6 = vel + k1_vel * (dt * DOPRI_A61) + k2_vel * (dt * DOPRI_A62)
                       + k3_vel * (dt * DOPRI_A63) + k4_vel * (dt * DOPRI_A64) + k5_vel * (dt * DOPRI_A65);
        Vec4 k6_pos = vel6, k6_vel = computeAcceleration(pos6, vel6);
        
        // 5th order solution
        Vec4 pos_new = pos + (k1_pos * DOPRI_B1 + k3_pos * DOPRI_B3 + k4_pos * DOPRI_B4 
                            + k5_pos * DOPRI_B5 + k6_pos * DOPRI_B6) * dt;
        Vec4 vel_new = vel + (k1_vel * DOPRI_B1 + k3_vel * DOPRI_B3 + k4_vel * DOPRI_B4 
                            + k5_vel * DOPRI_B5 + k6_vel * DOPRI_B6) * dt;
        
        // Error estimate
        Vec4 err_pos = (k1_pos * DOPRI_E1 + k3_pos * DOPRI_E3 + k4_pos * DOPRI_E4 
                      + k5_pos * DOPRI_E5 + k6_pos * DOPRI_E6) * dt;
        float rel_err = (fabsf(err_pos.r) + fabsf(vel.r) * dt) / (fabsf(pos.r) + fabsf(vel.r) * dt + 1e-10f);
        
        if (rel_err <= TOL || dt <= MIN_STEP) {
            m_State.position = pos_new;
            m_State.velocity = vel_new;
            accepted = true;
            if (rel_err > 0.0f) dt = fminf(dt * fminf(SAFETY * powf(TOL / rel_err, 0.2f), 5.0f), MAX_STEP);
        } else {
            dt = fmaxf(dt * fmaxf(SAFETY * powf(TOL / rel_err, 0.25f), 0.1f), MIN_STEP);
        }
    }
    
    // Normalize angles and velocity
    m_State.position.theta = fmaxf(0.01f, fminf(3.13159f, m_State.position.theta));
    while (m_State.position.phi < 0.0f) m_State.position.phi += 6.28318f;
    while (m_State.position.phi >= 6.28318f) m_State.position.phi -= 6.28318f;
    normalizeVelocity();
}


Vec4 GeodesicCamera::computeAcceleration(const Vec4& pos, const Vec4& vel) {
    // a^μ = -Γ^μ_αβ u^α u^β
    // For Schwarzschild metric (simplified - full Kerr would use metric interface)
    
    float r = pos.r;
    float theta = pos.theta;
    float rs = 2.0f * m_M;
    
    // Schwarzschild Christoffel symbols (non-zero components)
    float f = 1.0f - rs / r;
    float df_dr = rs / (r * r);
    
    Vec4 acc;
    
    // Γ^t_tr = (rs/2) / (r(r-rs)) * u^t * u^r
    acc.t = -df_dr / (2.0f * f) * vel.t * vel.r * 2.0f;
    
    // Γ^r_tt = (rs/2r²)(1-rs/r) * (u^t)²
    // Γ^r_rr = -(rs/2r²)/(1-rs/r) * (u^r)²
    // Γ^r_θθ = -(r-rs) * (u^θ)²
    // Γ^r_φφ = -(r-rs)sin²θ * (u^φ)²
    acc.r = -0.5f * df_dr * f * vel.t * vel.t 
            + 0.5f * df_dr / f * vel.r * vel.r
            - (r - rs) * vel.theta * vel.theta
            - (r - rs) * sinf(theta) * sinf(theta) * vel.phi * vel.phi;
    
    // Γ^θ_rθ = 1/r * u^r * u^θ
    // Γ^θ_φφ = -sinθcosθ * (u^φ)²
    acc.theta = -2.0f / r * vel.r * vel.theta
                + sinf(theta) * cosf(theta) * vel.phi * vel.phi;
    
    // Γ^φ_rφ = 1/r * u^r * u^φ
    // Γ^φ_θφ = cotθ * u^θ * u^φ
    acc.phi = -2.0f / r * vel.r * vel.phi
              - 2.0f * cosf(theta) / sinf(theta) * vel.theta * vel.phi;
    
    return acc;
}

void GeodesicCamera::normalizeVelocity() {
    // Enforce g_μν u^μ u^ν = -1 for timelike geodesic
    float r = m_State.position.r;
    float theta = m_State.position.theta;
    float rs = 2.0f * m_M;
    
    // Schwarzschild metric components
    float g_tt = -(1.0f - rs / r);
    float g_rr = 1.0f / (1.0f - rs / r);
    float g_thth = r * r;
    float g_phph = r * r * sinf(theta) * sinf(theta);
    
    // Current norm (for diagnostic purposes - should be -1 for timelike)
    [[maybe_unused]] float norm = g_tt * m_State.velocity.t * m_State.velocity.t
               + g_rr * m_State.velocity.r * m_State.velocity.r
               + g_thth * m_State.velocity.theta * m_State.velocity.theta
               + g_phph * m_State.velocity.phi * m_State.velocity.phi;

    // We need norm = -1, so solve for u^t given spatial components
    // g_tt (u^t)² + spatial_terms = -1
    // (u^t)² = (-1 - spatial_terms) / g_tt
    
    float spatial = g_rr * m_State.velocity.r * m_State.velocity.r
                  + g_thth * m_State.velocity.theta * m_State.velocity.theta
                  + g_phph * m_State.velocity.phi * m_State.velocity.phi;
    
    float ut_sq = (-1.0f - spatial) / g_tt;
    if (ut_sq > 0.0f) {
        m_State.velocity.t = sqrtf(ut_sq);  // Future-directed
    }
    
    // Compute Lorentz factor
    m_State.gamma = m_State.velocity.t * sqrtf(-g_tt);
}

//==============================================================================
// ISCO Calculation
//==============================================================================
float GeodesicCamera::computeISCO(bool prograde) const {
    // ISCO for Kerr black hole
    // r_isco = M * (3 + Z2 ∓ sqrt((3-Z1)(3+Z1+2*Z2)))
    // where Z1 = 1 + (1-a²)^(1/3) * ((1+a)^(1/3) + (1-a)^(1/3))
    //       Z2 = sqrt(3*a² + Z1²)
    // ∓ is - for prograde, + for retrograde
    
    float a_norm = m_a / m_M;  // Dimensionless spin
    
    if (fabsf(a_norm) < 0.001f) {
        // Schwarzschild: r_isco = 6M
        return 6.0f * m_M;
    }
    
    float a2 = a_norm * a_norm;
    float cbrt_1ma2 = powf(1.0f - a2, 1.0f / 3.0f);
    float cbrt_1pa = powf(1.0f + a_norm, 1.0f / 3.0f);
    float cbrt_1ma = powf(1.0f - a_norm, 1.0f / 3.0f);
    
    float Z1 = 1.0f + cbrt_1ma2 * (cbrt_1pa + cbrt_1ma);
    float Z2 = sqrtf(3.0f * a2 + Z1 * Z1);
    
    float r_isco;
    if (prograde) {
        r_isco = m_M * (3.0f + Z2 - sqrtf((3.0f - Z1) * (3.0f + Z1 + 2.0f * Z2)));
    } else {
        r_isco = m_M * (3.0f + Z2 + sqrtf((3.0f - Z1) * (3.0f + Z1 + 2.0f * Z2)));
    }
    
    return r_isco;
}

Vec4 GeodesicCamera::computeCircularVelocity(float r, bool prograde) const {
    // Circular orbit angular velocity: Ω = ±sqrt(M) / (r^1.5 ± a*sqrt(M))
    float r_phys = r * m_M;
    float omega = sqrtf(m_M) / (powf(r_phys, 1.5f) + (prograde ? 1.0f : -1.0f) * m_a * sqrtf(m_M));
    if (!prograde) omega = -omega;
    
    // For circular orbit: u^φ = Ω * u^t
    // u^t is determined by normalization constraint
    Vec4 vel;
    vel.r = 0.0f;
    vel.theta = 0.0f;
    
    // Compute u^t from circular orbit energy
    // For Schwarzschild: u^t = 1 / sqrt(1 - 3M/r) for r > 3M
    // Note: f = 1 - 2M/r is the lapse function (kept for reference)
    [[maybe_unused]] float f = 1.0f - 2.0f * m_M / r_phys;
    float orbit_factor = 1.0f - 3.0f * m_M / r_phys;
    
    if (orbit_factor > 0.0f) {
        vel.t = 1.0f / sqrtf(orbit_factor);
        vel.phi = omega * vel.t;
    } else {
        // Inside photon sphere - unstable
        vel.t = 1.0f;
        vel.phi = omega;
    }
    
    return vel;
}

} // namespace Sirius
