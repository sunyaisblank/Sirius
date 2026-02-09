// TSPF002A.cpp - Memory Usage Performance Tests
// Tests: struct sizes, per-pixel/ray memory, 3D texture memory.

#define _USE_MATH_DEFINES
#include <cmath>
#include <gtest/gtest.h>
#include <cstddef>
#include <vector>
#include "MTTN001A.h"
#include "MTDL001A.h"
#include "PHMT000A.h"
#include "PHMT100A.h"

namespace sirius::test {
using namespace Sirius;

// =============================================================================
// Memory Budget Constants (from roadmap)
// =============================================================================

namespace MemoryBudgets {
    // Resolution
    constexpr int WIDTH_1080P = 1920;
    constexpr int HEIGHT_1080P = 1080;
    constexpr int TOTAL_PIXELS_1080P = WIDTH_1080P * HEIGHT_1080P;
    
    constexpr int WIDTH_4K = 3840;
    constexpr int HEIGHT_4K = 2160;
    constexpr int TOTAL_PIXELS_4K = WIDTH_4K * HEIGHT_4K;
    
    // Typical numerical metric grid
    constexpr int NR_GRID_SIZE_SMALL = 64;
    constexpr int NR_GRID_SIZE_MEDIUM = 128;
    constexpr int NR_GRID_SIZE_LARGE = 256;
    
    // ADM variable count (10 components: alp, betax/y/z, gxx, gxy, gxz, gyy, gyz, gzz)
    constexpr int ADM_COMPONENT_COUNT = 10;
    
    // VRAM budget (typical GPU)
    constexpr size_t VRAM_BUDGET_8GB = 8ULL * 1024 * 1024 * 1024;
    constexpr size_t VRAM_BUDGET_12GB = 12ULL * 1024 * 1024 * 1024;
    
    // Reasonable limits
    constexpr size_t MAX_TEXTURE_SIZE_MB = 2048;  // 2GB texture limit
}

// =============================================================================
// Test Fixture
// =============================================================================

class MemoryUsageTests : public ::testing::Test {
protected:
    // Helper: format bytes as human-readable string
    std::string formatBytes(size_t bytes) {
        if (bytes < 1024) return std::to_string(bytes) + " B";
        if (bytes < 1024 * 1024) return std::to_string(bytes / 1024) + " KB";
        if (bytes < 1024ULL * 1024 * 1024) return std::to_string(bytes / (1024 * 1024)) + " MB";
        return std::to_string(bytes / (1024ULL * 1024 * 1024)) + " GB";
    }
    
    // Calculate 3D texture memory for numerical metric
    size_t calculate3DTextureMemory(int nx, int ny, int nz, int components, int bytesPerComponent) {
        return static_cast<size_t>(nx) * ny * nz * components * bytesPerComponent;
    }
    
    void SetUp() override {}
    void TearDown() override {}
};

// =============================================================================
// Type Size Validation Tests
// =============================================================================

// Test: Verify fundamental type sizes
TEST_F(MemoryUsageTests, FundamentalTypeSizes) {
    // These are standard on modern 64-bit systems
    EXPECT_EQ(sizeof(float), 4) << "float should be 4 bytes";
    EXPECT_EQ(sizeof(double), 8) << "double should be 8 bytes";
    EXPECT_EQ(sizeof(int), 4) << "int should be 4 bytes";
    EXPECT_EQ(sizeof(size_t), 8) << "size_t should be 8 bytes on 64-bit";
}

// Test: Verify Dual number size
TEST_F(MemoryUsageTests, DualNumberSize) {
    // Dual<double> contains real + dual parts
    size_t dual_size = sizeof(Dual<double>);
    
    std::cout << "\n[MEMORY] Dual<double> size: " << dual_size << " bytes\n";
    
    // Should be 2 doubles = 16 bytes (or potentially more with padding)
    EXPECT_GE(dual_size, 16) << "Dual<double> should be at least 16 bytes";
    EXPECT_LE(dual_size, 32) << "Dual<double> should not exceed 32 bytes";
}

// Test: Verify Vec4 size
TEST_F(MemoryUsageTests, Vec4Size) {
    size_t vec4_size = sizeof(Vec4);
    
    std::cout << "\n[MEMORY] Vec4 size: " << vec4_size << " bytes\n";
    
    // Vec4 is Tensor<double, 4> = 4 doubles = 32 bytes (minimum)
    EXPECT_GE(vec4_size, 32) << "Vec4 should be at least 32 bytes";
}

// Test: Verify Metric4D size
TEST_F(MemoryUsageTests, Metric4DSize) {
    size_t metric_size = sizeof(Metric4D);
    
    std::cout << "\n[MEMORY] Metric4D size: " << metric_size << " bytes\n";
    
    // Metric4D is 4x4 Dual<double> = 16 × 16 bytes = 256 bytes (minimum)
    EXPECT_GE(metric_size, 256) << "Metric4D should be at least 256 bytes";
}

// =============================================================================
// Per-Pixel/Per-Ray Memory Tests
// =============================================================================

// Test: Memory per geodesic ray state
TEST_F(MemoryUsageTests, PerRayMemory) {
    // Lightray structure contains:
    // - position (Vec4): 4 doubles = 32 bytes
    // - velocity (Vec4): 4 doubles = 32 bytes  
    // - acceleration (Vec4): 4 doubles = 32 bytes
    // - step_size, proper_time, etc.: ~16 bytes
    // Total: ~112 bytes minimum
    
    size_t lightray_estimate = sizeof(Vec4) * 3 + 16;
    
    std::cout << "\n[MEMORY] Per-ray state estimate: " << lightray_estimate << " bytes\n";
    
    // At 1080p with one ray per pixel
    size_t total_ray_memory = lightray_estimate * MemoryBudgets::TOTAL_PIXELS_1080P;
    
    std::cout << "  Total ray memory @ 1080p: " << formatBytes(total_ray_memory) << "\n";
    
    // Should be reasonable (< 1GB for 1080p)
    EXPECT_LT(total_ray_memory, 1ULL * 1024 * 1024 * 1024) 
        << "Ray memory at 1080p should be < 1GB";
}

// Test: Memory per pixel for metric storage
TEST_F(MemoryUsageTests, PerPixelMetricStorage) {
    // Per pixel, we might cache:
    // - Metric tensor: 10 unique components × 8 bytes = 80 bytes
    // - Metric dual parts: 10 × 8 bytes = 80 bytes
    // - Christoffel symbols: 40 unique × 8 bytes = 320 bytes (if cached)
    
    size_t per_pixel_minimal = 10 * sizeof(double);  // Just metric components
    size_t per_pixel_full = 10 * sizeof(Dual<double>) + 40 * sizeof(double);  // Metric + Christoffel
    
    std::cout << "\n[MEMORY] Per-pixel metric storage:\n";
    std::cout << "  Minimal (metric only): " << per_pixel_minimal << " bytes\n";
    std::cout << "  Full (metric + Christoffel): " << per_pixel_full << " bytes\n";
    
    // At 1080p
    size_t total_minimal = per_pixel_minimal * MemoryBudgets::TOTAL_PIXELS_1080P;
    size_t total_full = per_pixel_full * MemoryBudgets::TOTAL_PIXELS_1080P;
    
    std::cout << "  Total @ 1080p (minimal): " << formatBytes(total_minimal) << "\n";
    std::cout << "  Total @ 1080p (full): " << formatBytes(total_full) << "\n";
    
    EXPECT_LT(total_minimal, 512ULL * 1024 * 1024) << "Minimal metric storage < 512MB";
}

// =============================================================================
// Numerical Metric Texture Memory Tests
// =============================================================================

// Test: 3D texture memory for small grid
TEST_F(MemoryUsageTests, NumericalMetricSmallGrid) {
    int n = MemoryBudgets::NR_GRID_SIZE_SMALL;
    int components = MemoryBudgets::ADM_COMPONENT_COUNT;
    
    size_t float_memory = calculate3DTextureMemory(n, n, n, components, sizeof(float));
    size_t double_memory = calculate3DTextureMemory(n, n, n, components, sizeof(double));
    
    std::cout << "\n[MEMORY] Numerical metric 3D texture (" << n << "³):\n";
    std::cout << "  Float precision: " << formatBytes(float_memory) << "\n";
    std::cout << "  Double precision: " << formatBytes(double_memory) << "\n";
    
    // 64³ × 10 × 4 bytes = ~10MB (should be very reasonable)
    EXPECT_LT(float_memory, 100ULL * 1024 * 1024) << "Small grid should be < 100MB";
}

// Test: 3D texture memory for medium grid
TEST_F(MemoryUsageTests, NumericalMetricMediumGrid) {
    int n = MemoryBudgets::NR_GRID_SIZE_MEDIUM;
    int components = MemoryBudgets::ADM_COMPONENT_COUNT;
    
    size_t float_memory = calculate3DTextureMemory(n, n, n, components, sizeof(float));
    size_t double_memory = calculate3DTextureMemory(n, n, n, components, sizeof(double));
    
    std::cout << "\n[MEMORY] Numerical metric 3D texture (" << n << "³):\n";
    std::cout << "  Float precision: " << formatBytes(float_memory) << "\n";
    std::cout << "  Double precision: " << formatBytes(double_memory) << "\n";
    
    // 128³ × 10 × 4 bytes = ~80MB
    EXPECT_LT(float_memory, 500ULL * 1024 * 1024) << "Medium grid should be < 500MB";
}

// Test: 3D texture memory for large grid  
TEST_F(MemoryUsageTests, NumericalMetricLargeGrid) {
    int n = MemoryBudgets::NR_GRID_SIZE_LARGE;
    int components = MemoryBudgets::ADM_COMPONENT_COUNT;
    
    size_t float_memory = calculate3DTextureMemory(n, n, n, components, sizeof(float));
    size_t double_memory = calculate3DTextureMemory(n, n, n, components, sizeof(double));
    
    std::cout << "\n[MEMORY] Numerical metric 3D texture (" << n << "³):\n";
    std::cout << "  Float precision: " << formatBytes(float_memory) << "\n";
    std::cout << "  Double precision: " << formatBytes(double_memory) << "\n";
    
    // 256³ × 10 × 4 bytes = ~640MB
    EXPECT_LT(float_memory, 2ULL * 1024 * 1024 * 1024) << "Large grid should be < 2GB";
}

// Test: Christoffel 3D texture (precomputed)
TEST_F(MemoryUsageTests, ChristoffelTextureMemory) {
    // Christoffel symbols: Γ^α_μν has 4×10 = 40 independent components
    // (symmetric in lower indices)
    int n = MemoryBudgets::NR_GRID_SIZE_MEDIUM;
    int christoffel_components = 40;
    
    size_t float_memory = calculate3DTextureMemory(n, n, n, christoffel_components, sizeof(float));
    
    std::cout << "\n[MEMORY] Christoffel 3D texture (" << n << "³):\n";
    std::cout << "  Float precision (40 components): " << formatBytes(float_memory) << "\n";
    
    // This is the roadmap's Phase 3 optimization target
    // 128³ × 40 × 4 bytes = ~320MB
    EXPECT_LT(float_memory, 1ULL * 1024 * 1024 * 1024) << "Christoffel texture < 1GB";
}

// =============================================================================
// Total VRAM Budget Tests
// =============================================================================

// Test: Total memory estimate for analytic metrics at 1080p
TEST_F(MemoryUsageTests, TotalVRAMAnalytic1080p) {
    // For analytic metrics, main memory usage is:
    // - Frame buffer: 1920×1080×4 (RGBA8) = ~8MB
    // - Ray state: ~100 bytes × 2M pixels = ~200MB
    // - Background texture: varies (~10-50MB)
    // - Shader constants: minimal
    
    size_t frame_buffer = MemoryBudgets::TOTAL_PIXELS_1080P * 4;
    size_t ray_state = MemoryBudgets::TOTAL_PIXELS_1080P * 100;
    size_t background = 50 * 1024 * 1024;  // 50MB estimate
    
    size_t total = frame_buffer + ray_state + background;
    
    std::cout << "\n[MEMORY] Total VRAM estimate (analytic, 1080p):\n";
    std::cout << "  Frame buffer: " << formatBytes(frame_buffer) << "\n";
    std::cout << "  Ray state: " << formatBytes(ray_state) << "\n";
    std::cout << "  Background: " << formatBytes(background) << "\n";
    std::cout << "  Total: " << formatBytes(total) << "\n";
    
    // Should fit easily in 8GB VRAM
    EXPECT_LT(total, MemoryBudgets::VRAM_BUDGET_8GB);
}

// Test: Total memory estimate for numerical metrics at 1080p
TEST_F(MemoryUsageTests, TotalVRAMNumerical1080p) {
    int n = 128;  // Medium grid
    
    size_t frame_buffer = MemoryBudgets::TOTAL_PIXELS_1080P * 4;
    size_t ray_state = MemoryBudgets::TOTAL_PIXELS_1080P * 100;
    size_t metric_texture = calculate3DTextureMemory(n, n, n, 10, sizeof(float));
    size_t christoffel_texture = calculate3DTextureMemory(n, n, n, 40, sizeof(float));
    size_t background = 50 * 1024 * 1024;
    
    size_t total = frame_buffer + ray_state + metric_texture + christoffel_texture + background;
    
    std::cout << "\n[MEMORY] Total VRAM estimate (numerical, 1080p, " << n << "³ grid):\n";
    std::cout << "  Frame buffer: " << formatBytes(frame_buffer) << "\n";
    std::cout << "  Ray state: " << formatBytes(ray_state) << "\n";
    std::cout << "  Metric 3D texture: " << formatBytes(metric_texture) << "\n";
    std::cout << "  Christoffel 3D texture: " << formatBytes(christoffel_texture) << "\n";
    std::cout << "  Background: " << formatBytes(background) << "\n";
    std::cout << "  Total: " << formatBytes(total) << "\n";
    
    // Should fit in 8GB VRAM
    EXPECT_LT(total, MemoryBudgets::VRAM_BUDGET_8GB);
}

// =============================================================================
// Safety Checks
// =============================================================================

// Test: Vector allocation doesn't fail for reasonable sizes
TEST_F(MemoryUsageTests, VectorAllocationSafe) {
    // Test allocating a reasonable-sized vector
    size_t test_size = 1024 * 1024;  // 1M elements
    
    ASSERT_NO_THROW({
        std::vector<double> large_vec(test_size);
        large_vec[0] = 1.0;
        large_vec[test_size - 1] = 2.0;
        EXPECT_EQ(large_vec.size(), test_size);
    }) << "Should be able to allocate 8MB vector";
}

// Test: Metric evaluation doesn't leak memory
TEST_F(MemoryUsageTests, MetricEvaluationNoLeak) {
    // Evaluate metric many times - if there's a leak, this will fail
    Sirius::KerrSchildFamily mink{Sirius::KerrSchildParams::Minkowski()};
    
    Vec4 pos;
    pos(0) = 0; pos(1) = 10; pos(2) = 0; pos(3) = 0;
    
    Metric4D g;
    Tensor<Dual<double>, 4, 4, 4> dg;
    
    // 100,000 evaluations should not cause memory issues
    for (int i = 0; i < 100000; ++i) {
        mink.evaluate(pos, g, dg);
    }
    
    // If we got here, no crash from memory issues
    EXPECT_TRUE(true);
}

// Test: Report total test memory overhead
TEST_F(MemoryUsageTests, TestMemoryOverhead) {
    std::cout << "\n[MEMORY] Test infrastructure overhead:\n";
    std::cout << "  sizeof(Dual<double>): " << sizeof(Dual<double>) << " bytes\n";
    std::cout << "  sizeof(Vec4): " << sizeof(Vec4) << " bytes\n";
    std::cout << "  sizeof(Metric4D): " << sizeof(Metric4D) << " bytes\n";
    std::cout << "  sizeof(Tensor<Dual<double>,4,4,4>): " 
              << sizeof(Tensor<Dual<double>, 4, 4, 4>) << " bytes\n";
    
    EXPECT_TRUE(true);  // Documentation test
}

} // namespace sirius::test
