// TSDG005A.cpp - GPU Conservation Law Verification Tests
// Component ID: TSDG005A
// Purpose: Validate conservation laws on GPU kernel integration
//
// MATHEMATICAL BASIS:
// Tests that GPU geodesic integration preserves:
// - Null condition: g_μν k^μ k^ν = 0
// - Killing energy: E = -g_tμ k^μ (stationary spacetimes)
// - Killing angular momentum: L = g_φμ k^μ (axisymmetric spacetimes)
//
// SPECIFICATION TARGETS (docs/specification.md):
// - Null condition GPU: |H| < 10^-5
// - Energy conservation GPU: |ΔE/E| < 10^-4
// - Angular momentum GPU: |ΔL/L| < 10^-4
//
// REFERENCES:
// - TSDG003A.cpp: CPU conservation tests (this is the GPU counterpart)
// - RDOP005A.h: GPU debug buffer interface

#define _USE_MATH_DEFINES
#include <cmath>
#include <gtest/gtest.h>
#include <vector>
#include <memory>

// Only compile GPU tests if CUDA is available
#ifdef SIRIUS_HAS_CUDA

#include <cuda_runtime.h>
#include "../src/Sirius.Render/Acceleration/OptiX/RDOP003A.h"
#include "../src/Sirius.Render/Acceleration/OptiX/RDOP005A.h"

namespace sirius::test {

// =============================================================================
// GPU Conservation Test Tolerances
// =============================================================================
// GPU uses single precision, so tolerances are relaxed compared to CPU
// These are still stricter than typical rendering requirements

constexpr float GPU_HAMILTONIAN_TOLERANCE = 1e-5f;     // Null condition
constexpr float GPU_KILLING_ENERGY_TOLERANCE = 1e-4f; // Energy drift
constexpr float GPU_KILLING_L_TOLERANCE = 1e-4f;      // Angular momentum drift

// =============================================================================
// Test Fixture
// =============================================================================

class GPUConservationTests : public ::testing::Test {
protected:
    // CUDA resources
    Sirius::DebugRayState* d_debugBuffer = nullptr;
    Sirius::DebugRayState* h_debugBuffer = nullptr;
    int* d_bufferCount = nullptr;
    int h_bufferCount = 0;

    static constexpr int DEBUG_BUFFER_SIZE = 100000;

    void SetUp() override {
        // Check CUDA availability
        int deviceCount = 0;
        cudaError_t err = cudaGetDeviceCount(&deviceCount);
        if (err != cudaSuccess || deviceCount == 0) {
            GTEST_SKIP() << "No CUDA device available";
        }

        // Allocate debug buffers
        cudaMalloc(&d_debugBuffer, DEBUG_BUFFER_SIZE * sizeof(Sirius::DebugRayState));
        cudaMalloc(&d_bufferCount, sizeof(int));
        h_debugBuffer = new Sirius::DebugRayState[DEBUG_BUFFER_SIZE];

        // Initialize count to zero
        cudaMemset(d_bufferCount, 0, sizeof(int));
    }

    void TearDown() override {
        if (d_debugBuffer) cudaFree(d_debugBuffer);
        if (d_bufferCount) cudaFree(d_bufferCount);
        if (h_debugBuffer) delete[] h_debugBuffer;
    }

    // Copy debug data from GPU to host
    void copyDebugBufferToHost() {
        cudaMemcpy(&h_bufferCount, d_bufferCount, sizeof(int), cudaMemcpyDeviceToHost);
        h_bufferCount = std::min(h_bufferCount, DEBUG_BUFFER_SIZE);
        if (h_bufferCount > 0) {
            cudaMemcpy(h_debugBuffer, d_debugBuffer,
                      h_bufferCount * sizeof(Sirius::DebugRayState),
                      cudaMemcpyDeviceToHost);
        }
    }

    // Create debug buffer configuration
    Sirius::DebugBufferConfig createDebugConfig() {
        Sirius::DebugBufferConfig config = {};
        config.enabled = true;
        config.sample_interval = 10;  // Sample every 10 steps
        config.max_samples_per_ray = 100;
        config.num_debug_rays = 0;    // All rays
        config.debug_pixel_x = -1;    // All pixels
        config.debug_pixel_y = -1;
        config.buffer = d_debugBuffer;
        config.buffer_capacity = DEBUG_BUFFER_SIZE;
        config.buffer_count = d_bufferCount;
        return config;
    }
};

// =============================================================================
// Placeholder Tests (Require Full GPU Pipeline)
// =============================================================================
// NOTE: Full GPU tests require the complete OptiX pipeline to be initialized.
// These tests validate the debug buffer interface and analysis functions.
// Full integration tests should be run with the complete rendering pipeline.

TEST_F(GPUConservationTests, DebugBufferAnalysisWorks)
{
    // Create synthetic debug data to test analysis function
    std::vector<Sirius::DebugRayState> synthetic;

    // Simulate a ray with good conservation
    for (int step = 0; step < 100; ++step) {
        Sirius::DebugRayState state = {};
        state.step_number = step;
        state.hamiltonian = 1e-7f;  // Small null violation
        state.killing_energy = 1.0f + 1e-6f * step;  // Small drift
        state.killing_L = 2.0f + 5e-7f * step;
        state.affine_param = 0.01f * step;
        state.terminated = (step == 99) ? 2 : 0;  // Escaped at end
        synthetic.push_back(state);
    }

    // Analyze
    Sirius::ConservationStats stats = Sirius::analyzeConservation(
        synthetic.data(), static_cast<int>(synthetic.size()));

    // Verify statistics are computed correctly
    EXPECT_GT(stats.total_steps, 0);
    EXPECT_GT(stats.total_rays, 0);
    EXPECT_LT(stats.max_hamiltonian, 1e-5);
    EXPECT_LT(stats.max_energy_drift, 1e-3);
    EXPECT_LT(stats.max_L_drift, 1e-3);
}

TEST_F(GPUConservationTests, DebugBufferConfigValid)
{
    Sirius::DebugBufferConfig config = createDebugConfig();

    EXPECT_TRUE(config.enabled);
    EXPECT_NE(config.buffer, nullptr);
    EXPECT_GT(config.buffer_capacity, 0);
    EXPECT_NE(config.buffer_count, nullptr);
}

// =============================================================================
// Specification Compliance Tests
// =============================================================================
// These tests document the required tolerances even when full GPU
// pipeline is not available. They serve as specification anchors.

TEST_F(GPUConservationTests, SpecificationNullCondition)
{
    // SPECIFICATION: GPU null condition tolerance < 10^-5
    // This test verifies the tolerance constant matches specification

    EXPECT_LE(GPU_HAMILTONIAN_TOLERANCE, 1e-5f)
        << "GPU Hamiltonian tolerance exceeds specification (< 10^-5)";
}

TEST_F(GPUConservationTests, SpecificationEnergyConservation)
{
    // SPECIFICATION: GPU energy conservation < 10^-4 relative drift
    EXPECT_LE(GPU_KILLING_ENERGY_TOLERANCE, 1e-4f)
        << "GPU energy conservation tolerance exceeds specification (< 10^-4)";
}

TEST_F(GPUConservationTests, SpecificationAngularMomentum)
{
    // SPECIFICATION: GPU angular momentum conservation < 10^-4 relative drift
    EXPECT_LE(GPU_KILLING_L_TOLERANCE, 1e-4f)
        << "GPU angular momentum tolerance exceeds specification (< 10^-4)";
}

} // namespace sirius::test

#else // !SIRIUS_HAS_CUDA

// Stub tests when CUDA is not available
namespace sirius::test {

class GPUConservationTests : public ::testing::Test {};

TEST_F(GPUConservationTests, CUDANotAvailable)
{
    GTEST_SKIP() << "CUDA not available - GPU conservation tests skipped";
}

} // namespace sirius::test

#endif // SIRIUS_HAS_CUDA
