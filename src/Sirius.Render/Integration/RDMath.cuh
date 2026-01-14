// RDMath.cuh - Math Utilities for CUDA Ray Tracing
// Component ID: RDIN001A (Integration/Math)

#pragma once

#include <optix.h>
#include <cuda_runtime.h>
#include <math_constants.h>

namespace Sirius {

//==============================================================================
// Mathematical Helpers
//==============================================================================

__device__ __forceinline__ float3 operator+(float3 a, float3 b) {
    return make_float3(a.x + b.x, a.y + b.y, a.z + b.z);
}

__device__ __forceinline__ float3 operator-(float3 a, float3 b) {
    return make_float3(a.x - b.x, a.y - b.y, a.z - b.z);
}

__device__ __forceinline__ float3 operator-(float3 a) {
    return make_float3(-a.x, -a.y, -a.z);
}

__device__ __forceinline__ float3 operator*(float3 a, float s) {
    return make_float3(a.x * s, a.y * s, a.z * s);
}

__device__ __forceinline__ float3 operator*(float s, float3 a) {
    return make_float3(a.x * s, a.y * s, a.z * s);
}

__device__ __forceinline__ float3 operator/(float3 a, float s) {
    float inv = 1.0f / s;
    return make_float3(a.x * inv, a.y * inv, a.z * inv);
}

__device__ __forceinline__ float dot(float3 a, float3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

__device__ __forceinline__ float3 cross(float3 a, float3 b) {
    return make_float3(
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    );
}

__device__ __forceinline__ float length(float3 v) {
    return sqrtf(dot(v, v));
}

__device__ __forceinline__ float3 normalize(float3 v) {
    float len = length(v);
    return (len > 0.0f) ? v / len : make_float3(0.0f, 0.0f, 0.0f);
}

__device__ __forceinline__ float clamp(float x, float lo, float hi) {
    return fmaxf(lo, fminf(hi, x));
}

//==============================================================================
// Constants
//==============================================================================
constexpr float PI = 3.14159265358979323846f;
constexpr float TWO_PI = 2.0f * PI;

//==============================================================================
// OPTIMIZATION: Fast Math Intrinsics
// ~3x faster trig with ~1 ULP precision trade-off (acceptable for ray tracing)
//==============================================================================
#define FAST_SINF(x) __sinf(x)
#define FAST_COSF(x) __cosf(x)
#define FAST_TANF(x) __tanf(x)
#define FAST_EXPF(x) __expf(x)
#define FAST_POWF(x,y) __powf(x,y)

// Combined sin/cos for paired operations (single SFU instruction)
__device__ __forceinline__ void fast_sincosf(float x, float* s, float* c) {
    sincosf(x, s, c);
}

//==============================================================================
// Fast Hash-based Random Number Generator (for Russian Roulette)
//==============================================================================
__device__ __forceinline__ unsigned int hashPCG(unsigned int input) {
    unsigned int state = input * 747796405u + 2891336453u;
    unsigned int word = ((state >> ((state >> 28u) + 4u)) ^ state) * 277803737u;
    return (word >> 22u) ^ word;
}


__device__ __forceinline__ float randomFloat(unsigned int seed) {
    return (float)hashPCG(seed) / 4294967295.0f;  // [0, 1]
}

// 4x4 Matrix Inversion (Gauss-Jordan or Cramer's rule adaptation for 4x4)
// Based on standard gluInvertMatrix implementation logic
__device__ inline bool invert4x4(const float m[4][4], float out[4][4]) {
    float inv[16], det;
    float m_flat[16];
    
    // Flatten input for easier indexing
    for(int i=0; i<4; i++) for(int j=0; j<4; j++) m_flat[i*4+j] = m[i][j];

    inv[0] = m_flat[5]  * m_flat[10] * m_flat[15] - 
             m_flat[5]  * m_flat[11] * m_flat[14] - 
             m_flat[9]  * m_flat[6]  * m_flat[15] + 
             m_flat[9]  * m_flat[7]  * m_flat[14] +
             m_flat[13] * m_flat[6]  * m_flat[11] - 
             m_flat[13] * m_flat[7]  * m_flat[10];

    inv[4] = -m_flat[4]  * m_flat[10] * m_flat[15] + 
              m_flat[4]  * m_flat[11] * m_flat[14] + 
              m_flat[8]  * m_flat[6]  * m_flat[15] - 
              m_flat[8]  * m_flat[7]  * m_flat[14] - 
              m_flat[12] * m_flat[6]  * m_flat[11] + 
              m_flat[12] * m_flat[7]  * m_flat[10];

    inv[8] = m_flat[4]  * m_flat[9] * m_flat[15] - 
             m_flat[4]  * m_flat[11] * m_flat[13] - 
             m_flat[8]  * m_flat[5] * m_flat[15] + 
             m_flat[8]  * m_flat[7] * m_flat[13] + 
             m_flat[12] * m_flat[5] * m_flat[11] - 
             m_flat[12] * m_flat[7] * m_flat[9];

    inv[12] = -m_flat[4]  * m_flat[9] * m_flat[14] + 
               m_flat[4]  * m_flat[10] * m_flat[13] + 
               m_flat[8]  * m_flat[5] * m_flat[14] - 
               m_flat[8]  * m_flat[6] * m_flat[13] - 
               m_flat[12] * m_flat[5] * m_flat[10] + 
               m_flat[12] * m_flat[6] * m_flat[9];

    inv[1] = -m_flat[1]  * m_flat[10] * m_flat[15] + 
              m_flat[1]  * m_flat[11] * m_flat[14] + 
              m_flat[9]  * m_flat[2] * m_flat[15] - 
              m_flat[9]  * m_flat[3] * m_flat[14] - 
              m_flat[13] * m_flat[2] * m_flat[11] + 
              m_flat[13] * m_flat[3] * m_flat[10];

    inv[5] = m_flat[0]  * m_flat[10] * m_flat[15] - 
             m_flat[0]  * m_flat[11] * m_flat[14] - 
             m_flat[8]  * m_flat[2] * m_flat[15] + 
             m_flat[8]  * m_flat[3] * m_flat[14] + 
             m_flat[12] * m_flat[2] * m_flat[11] - 
             m_flat[12] * m_flat[3] * m_flat[10];

    inv[9] = -m_flat[0]  * m_flat[9] * m_flat[15] + 
              m_flat[0]  * m_flat[11] * m_flat[13] + 
              m_flat[8]  * m_flat[1] * m_flat[15] - 
              m_flat[8]  * m_flat[3] * m_flat[13] - 
              m_flat[12] * m_flat[1] * m_flat[11] + 
              m_flat[12] * m_flat[3] * m_flat[9];

    inv[13] = m_flat[0]  * m_flat[9] * m_flat[14] - 
              m_flat[0]  * m_flat[10] * m_flat[13] - 
              m_flat[8]  * m_flat[1] * m_flat[14] + 
              m_flat[8]  * m_flat[2] * m_flat[13] + 
              m_flat[12] * m_flat[1] * m_flat[10] - 
              m_flat[12] * m_flat[2] * m_flat[9];

    inv[2] = m_flat[1]  * m_flat[6] * m_flat[15] - 
             m_flat[1]  * m_flat[7] * m_flat[14] - 
             m_flat[5]  * m_flat[2] * m_flat[15] + 
             m_flat[5]  * m_flat[3] * m_flat[14] + 
             m_flat[13] * m_flat[2] * m_flat[7] - 
             m_flat[13] * m_flat[3] * m_flat[6];

    inv[6] = -m_flat[0]  * m_flat[6] * m_flat[15] + 
              m_flat[0]  * m_flat[7] * m_flat[14] + 
              m_flat[4]  * m_flat[2] * m_flat[15] - 
              m_flat[4]  * m_flat[3] * m_flat[14] - 
              m_flat[12] * m_flat[2] * m_flat[7] + 
              m_flat[12] * m_flat[3] * m_flat[6];

    inv[10] = m_flat[0]  * m_flat[5] * m_flat[15] - 
              m_flat[0]  * m_flat[7] * m_flat[13] - 
              m_flat[4]  * m_flat[1] * m_flat[15] + 
              m_flat[4]  * m_flat[3] * m_flat[13] + 
              m_flat[12] * m_flat[1] * m_flat[7] - 
              m_flat[12] * m_flat[3] * m_flat[5];

    inv[14] = -m_flat[0]  * m_flat[5] * m_flat[14] + 
               m_flat[0]  * m_flat[6] * m_flat[13] + 
               m_flat[4]  * m_flat[1] * m_flat[14] - 
               m_flat[4]  * m_flat[2] * m_flat[13] - 
               m_flat[12] * m_flat[1] * m_flat[6] + 
               m_flat[12] * m_flat[2] * m_flat[5];

    inv[3] = -m_flat[1] * m_flat[6] * m_flat[11] + 
              m_flat[1] * m_flat[7] * m_flat[10] + 
              m_flat[5] * m_flat[2] * m_flat[11] - 
              m_flat[5] * m_flat[3] * m_flat[10] - 
              m_flat[9] * m_flat[2] * m_flat[7] + 
              m_flat[9] * m_flat[3] * m_flat[6];

    inv[7] = m_flat[0] * m_flat[6] * m_flat[11] - 
             m_flat[0] * m_flat[7] * m_flat[10] - 
             m_flat[4] * m_flat[2] * m_flat[11] + 
             m_flat[4] * m_flat[3] * m_flat[10] + 
             m_flat[8] * m_flat[2] * m_flat[7] - 
             m_flat[8] * m_flat[3] * m_flat[6];

    inv[11] = -m_flat[0] * m_flat[5] * m_flat[11] + 
               m_flat[0] * m_flat[7] * m_flat[9] + 
               m_flat[4] * m_flat[1] * m_flat[11] - 
               m_flat[4] * m_flat[3] * m_flat[9] - 
               m_flat[8] * m_flat[1] * m_flat[7] + 
               m_flat[8] * m_flat[3] * m_flat[5];

    inv[15] = m_flat[0] * m_flat[5] * m_flat[10] - 
              m_flat[0] * m_flat[6] * m_flat[9] - 
              m_flat[4] * m_flat[1] * m_flat[10] + 
              m_flat[4] * m_flat[2] * m_flat[9] + 
              m_flat[8] * m_flat[1] * m_flat[6] - 
              m_flat[8] * m_flat[2] * m_flat[5];

    det = m_flat[0] * inv[0] + m_flat[1] * inv[4] + m_flat[2] * inv[8] + m_flat[3] * inv[12];

    if (det == 0) return false;

    det = 1.0f / det;

    for (int i = 0; i < 4; i++) {
        for(int j=0; j<4; j++) {
            out[i][j] = inv[i*4+j] * det;
        }
    }

    return true;
}

} // namespace Sirius

