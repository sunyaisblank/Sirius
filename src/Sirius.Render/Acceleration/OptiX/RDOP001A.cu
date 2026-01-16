// RDOP001A.cu - OptiX Host Code
//
// Context initialization, pipeline creation, SBT setup, CUDA/OpenGL interop.

#include "RDOP003A.h"

#include <optix.h>
#include <optix_stubs.h>
#include <optix_function_table_definition.h>

#include <cuda_runtime.h>
#include <cuda_gl_interop.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <stdexcept>
#include <cstring>

namespace Sirius {

//==============================================================================
// Error Checking Macros
//==============================================================================

#define CUDA_CHECK(call)                                                       \
    do {                                                                       \
        cudaError_t error = call;                                              \
        if (error != cudaSuccess) {                                            \
            std::cerr << "CUDA Error: " << cudaGetErrorString(error)           \
                      << " at " << __FILE__ << ":" << __LINE__ << std::endl;   \
            throw std::runtime_error("CUDA call failed");                      \
        }                                                                      \
    } while(0)

#define OPTIX_CHECK(call)                                                      \
    do {                                                                       \
        OptixResult result = call;                                             \
        if (result != OPTIX_SUCCESS) {                                         \
            std::cerr << "OptiX Error: " << optixGetErrorName(result)          \
                      << " (" << optixGetErrorString(result) << ")"            \
                      << " at " << __FILE__ << ":" << __LINE__ << std::endl;   \
            throw std::runtime_error("OptiX call failed");                     \
        }                                                                      \
    } while(0)

//==============================================================================
// OptiX Logging Callback
//==============================================================================
static void optixLogCallback(unsigned int level, const char* tag, 
                             const char* message, void* /*cbdata*/) {
    std::cerr << "[OptiX][" << level << "][" << tag << "]: " << message << std::endl;
}

//==============================================================================
// PTX Loading Utility
//==============================================================================
static std::string loadPtxFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.good()) {
        throw std::runtime_error("Failed to open PTX file: " + filename);
    }
    std::stringstream buffer;
    buffer << file.rdbuf();
    return buffer.str();
}

//==============================================================================
// Helper to create a single 3D texture
//==============================================================================
static cudaTextureObject_t create3DTexture(const float* data, int3 dims, cudaArray_t& outArray) {
    if (!data) return 0;
    
    cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
    cudaExtent extent = make_cudaExtent(dims.x, dims.y, dims.z);
    
    CUDA_CHECK(cudaMalloc3DArray(&outArray, &channelDesc, extent));
    
    cudaMemcpy3DParms copyParams = {0};
    copyParams.srcPtr = make_cudaPitchedPtr((void*)data, dims.x * sizeof(float), dims.x, dims.y);
    copyParams.dstArray = outArray;
    copyParams.extent = extent;
    copyParams.kind = cudaMemcpyHostToDevice;
    
    CUDA_CHECK(cudaMemcpy3D(&copyParams));
    
    cudaResourceDesc resDesc = {};
    resDesc.resType = cudaResourceTypeArray;
    resDesc.res.array.array = outArray;
    
    cudaTextureDesc texDesc = {};
    texDesc.addressMode[0] = cudaAddressModeClamp; // Clamp to edge
    texDesc.addressMode[1] = cudaAddressModeClamp;
    texDesc.addressMode[2] = cudaAddressModeClamp;
    texDesc.filterMode = cudaFilterModeLinear;     // Trilinear interpolation
    texDesc.readMode = cudaReadModeElementType;    // Read as float
    texDesc.normalizedCoords = 1;                  // Use [0,1] coords
    
    cudaTextureObject_t texObj = 0;
    CUDA_CHECK(cudaCreateTextureObject(&texObj, &resDesc, &texDesc, nullptr));
    return texObj;
}

static void update3DTexture(const float* data, int3 dims, cudaArray_t array) {
    if (!data || !array) return;
    
    cudaMemcpy3DParms copyParams = {0};
    copyParams.srcPtr = make_cudaPitchedPtr((void*)data, dims.x * sizeof(float), dims.x, dims.y);
    copyParams.dstArray = array;
    copyParams.extent = make_cudaExtent(dims.x, dims.y, dims.z);
    copyParams.kind = cudaMemcpyHostToDevice;
    
    CUDA_CHECK(cudaMemcpy3D(&copyParams));
}

// Container for numerical metric textures to manage lifetime
struct NumericalMetricTextureSet {
    std::vector<cudaArray_t> arrays;
    std::vector<cudaTextureObject_t> textures;
    int3 dims = {0, 0, 0};
    
    void cleanup() {
        for (auto tex : textures) {
            if (tex) cudaDestroyTextureObject(tex);
        }
        textures.clear();
        for (auto arr : arrays) {
            if (arr) cudaFreeArray(arr);
        }
        arrays.clear();
        dims = {0, 0, 0};
    }
};

// Container for Christoffel symbol textures (40 textures: 4 mu × 10 packed components)
struct ChristoffelTextureSet {
    cudaArray_t arrays[40];
    cudaTextureObject_t textures[40];
    bool initialized = false;
    
    ChristoffelTextureSet() {
        for (int i = 0; i < 40; i++) {
            arrays[i] = nullptr;
            textures[i] = 0;
        }
    }
    
    void cleanup() {
        for (int i = 0; i < 40; i++) {
            if (textures[i]) {
                cudaDestroyTextureObject(textures[i]);
                textures[i] = 0;
            }
            if (arrays[i]) {
                cudaFreeArray(arrays[i]);
                arrays[i] = nullptr;
            }
        }
        initialized = false;
    }
};

//==============================================================================
// OptixRenderer Class Implementation
//==============================================================================
class OptixRenderer {
public:
    OptixRenderer();
    ~OptixRenderer();

    // Initialization
    bool initialize(int width, int height);
    void cleanup();

    // Pipeline management
    bool createPipeline(const std::string& ptxPath);
    void buildSBT();

    // Rendering
    void launch(const LaunchParams& params);
    void resize(int width, int height);

    // Parameter updates
    void updateLaunchParams(const LaunchParams& params);
    void updateDisplay();

    // OpenGL interop
    void registerGLTexture(unsigned int glTexture);
    void unregisterGLTexture();
    float4* mapFrameBuffer();
    void unmapFrameBuffer();

    // Background texture
    bool uploadBackgroundTexture(const unsigned char* data, int width, int height);
    cudaTextureObject_t getBackgroundTexture() const { return m_BackgroundTexObj; }
    
    // Kernel Selection
    void setMetricType(int type); // 0=Minkowski, 1=Schwarzschild, etc.

    // Getters
    bool isInitialized() const { return m_Initialized; }
    int getWidth() const { return m_Width; }
    int getHeight() const { return m_Height; }

private:
    // Context
    CUcontext          m_CudaContext     = nullptr;
    OptixDeviceContext m_OptixContext    = nullptr;
    CUstream           m_Stream          = nullptr;

    // Pipeline
    OptixModule            m_Module           = nullptr;
    OptixPipeline          m_Pipeline         = nullptr;
    std::vector<OptixProgramGroup> m_RaygenPGs;
    OptixProgramGroup      m_MissPG           = nullptr;
    OptixProgramGroup      m_HitGroupPG       = nullptr;
    OptixPipelineCompileOptions m_PipelineCompileOptions = {};

    // Shader Binding Table
    std::vector<CUdeviceptr> m_RaygenRecords;
    CUdeviceptr m_MissRecord        = 0;
    CUdeviceptr m_HitGroupRecord    = 0;
    OptixShaderBindingTable m_SBT   = {};

    // Frame buffers
    CUdeviceptr m_FrameBuffer       = 0;
    CUdeviceptr m_AccumBuffer       = 0;
    CUdeviceptr m_LaunchParamsBuffer = 0;

    // Host-side frame buffer for CPU access (avoids segfault when reading GPU memory)
    std::vector<float> m_HostFrameBuffer;

    // OpenGL interop
    cudaGraphicsResource_t m_CudaGLResource = nullptr;
    bool m_GLInteropEnabled = false;

    // State
    bool m_Initialized = false;
    
    // Numerical Metric Resources
    NumericalMetricTextureSet m_NumericalMetricResources;
    
public:
    NumericalMetricData uploadNumericalMetric(const NumericalMetricHostData& hostData) {
        bool canReuse = (m_NumericalMetricResources.arrays.size() == 14) &&
                        (m_NumericalMetricResources.dims.x == hostData.dims.x) &&
                        (m_NumericalMetricResources.dims.y == hostData.dims.y) &&
                        (m_NumericalMetricResources.dims.z == hostData.dims.z);

        if (!canReuse) {
            m_NumericalMetricResources.cleanup();
            m_NumericalMetricResources.dims = hostData.dims;
            
            auto createAndStore = [&](const float* ptr) {
                cudaArray_t arr = nullptr;
                cudaTextureObject_t tex = create3DTexture(ptr, hostData.dims, arr);
                if (arr) m_NumericalMetricResources.arrays.push_back(arr);
                else m_NumericalMetricResources.arrays.push_back(nullptr); 
                
                if (tex) m_NumericalMetricResources.textures.push_back(tex); 
                else m_NumericalMetricResources.textures.push_back(0);
                
                return tex;
            };
            
            // Initial Creation (Order matters!)
            createAndStore(hostData.gxx); createAndStore(hostData.gxy); createAndStore(hostData.gxz);
            createAndStore(hostData.gyy); createAndStore(hostData.gyz); createAndStore(hostData.gzz);
            createAndStore(hostData.alp);
            createAndStore(hostData.betax); createAndStore(hostData.betay); createAndStore(hostData.betaz);
            createAndStore(hostData.rho);
            createAndStore(hostData.vx); createAndStore(hostData.vy); createAndStore(hostData.vz);
            
        } else {
            // Update Existing Arrays
            int idx = 0;
            auto update = [&](const float* ptr) {
                if (ptr && m_NumericalMetricResources.arrays[idx]) {
                    update3DTexture(ptr, hostData.dims, m_NumericalMetricResources.arrays[idx]);
                }
                idx++;
            };
            
            update(hostData.gxx); update(hostData.gxy); update(hostData.gxz);
            update(hostData.gyy); update(hostData.gyz); update(hostData.gzz);
            update(hostData.alp);
            update(hostData.betax); update(hostData.betay); update(hostData.betaz);
            update(hostData.rho);
            update(hostData.vx); update(hostData.vy); update(hostData.vz);
        }

        // Construct Data Struct
        NumericalMetricData nm = {};
        nm.dims = hostData.dims;
        nm.origin = hostData.origin;
        nm.spacing = hostData.spacing;
        nm.isLoaded = true;
        
        // Map from stored vector (Order must match creation!)
        const auto& texs = m_NumericalMetricResources.textures;
        nm.gxx = texs[0]; nm.gxy = texs[1]; nm.gxz = texs[2];
        nm.gyy = texs[3]; nm.gyz = texs[4]; nm.gzz = texs[5];
        nm.alp = texs[6];
        nm.betax = texs[7]; nm.betay = texs[8]; nm.betaz = texs[9];
        nm.rho = texs[10];
        nm.vx = texs[11]; nm.vy = texs[12]; nm.vz = texs[13];
        
        // Initialize Christoffel arrays to zero first
        for (int i = 0; i < 10; i++) {
            nm.Gamma_t[i] = 0;
            nm.Gamma_r[i] = 0;
            nm.Gamma_theta[i] = 0;
            nm.Gamma_phi[i] = 0;
        }
        nm.christoffelLoaded = false;
        
        // =====================================================================
        // PHASE 3.1: Precompute Christoffel symbols (Cartesian)
        // =====================================================================
        // The precomputed Christoffels are in Cartesian coordinates (t,x,y,z).
        // The GPU transforms them to spherical (t,r,θ,φ) at lookup time using:
        // Γ'^μ_νρ = (∂x'μ/∂xα)(∂xβ/∂x'ν)(∂xγ/∂x'ρ)Γα_βγ + (∂x'μ/∂xα)(∂²xα/∂x'ν∂x'ρ)
        // See transformChristoffelCartToSph() in RDOP002A.cu
        // =====================================================================
        precomputeChristoffelTextures(hostData, nm);
        
        return nm;
    }
    
    // =========================================================================
    // P2: Precompute Christoffel Symbols from Metric Data
    // =========================================================================
    // Computes Christoffel symbols Γ^μ_νρ via finite differences on CPU
    // and uploads them as 40 3D textures for fast GPU lookup.
    // This reduces per-integration-step texture reads from ~70 to 40.
    //
    // Call this AFTER uploadNumericalMetric() with the same hostData.
    // =========================================================================
    void precomputeChristoffelTextures(const NumericalMetricHostData& hostData,
                                        NumericalMetricData& nm) {
        if (!nm.isLoaded) {
            std::cerr << "[OptiX] Cannot precompute Christoffel: metric not loaded" << std::endl;
            return;
        }
        
        const int nx = hostData.dims.x;
        const int ny = hostData.dims.y;
        const int nz = hostData.dims.z;
        const int totalVoxels = nx * ny * nz;
        
        std::cout << "[OptiX] Precomputing Christoffel symbols (" << nx << "x" << ny << "x" << nz << ")..." << std::endl;
        
        // Allocate host memory for 40 Christoffel components
        std::vector<std::vector<float>> christoffelData(40, std::vector<float>(totalVoxels, 0.0f));
        
        // Helper to pack lower indices: when nu <= rho, idx = nu + rho*(rho+1)/2
        auto packIdx = [](int nu, int rho) -> int {
            if (nu > rho) std::swap(nu, rho);
            return nu + rho * (rho + 1) / 2;
        };
        
        // Helper to sample metric at grid position (periodic boundary)
        auto sampleMetric = [&](int ix, int iy, int iz, float g[4][4]) {
            // Clamp to grid bounds
            ix = std::max(0, std::min(nx-1, ix));
            iy = std::max(0, std::min(ny-1, iy));
            iz = std::max(0, std::min(nz-1, iz));
            
            int idx = ix + iy * nx + iz * nx * ny;
            
            // 3-metric (spatial)
            float gxx = hostData.gxx ? hostData.gxx[idx] : 1.0f;
            float gxy = hostData.gxy ? hostData.gxy[idx] : 0.0f;
            float gxz = hostData.gxz ? hostData.gxz[idx] : 0.0f;
            float gyy = hostData.gyy ? hostData.gyy[idx] : 1.0f;
            float gyz = hostData.gyz ? hostData.gyz[idx] : 0.0f;
            float gzz = hostData.gzz ? hostData.gzz[idx] : 1.0f;
            
            // Lapse and shift
            float alp = hostData.alp ? hostData.alp[idx] : 1.0f;
            float bx = hostData.betax ? hostData.betax[idx] : 0.0f;
            float by = hostData.betay ? hostData.betay[idx] : 0.0f;
            float bz = hostData.betaz ? hostData.betaz[idx] : 0.0f;
            
            // Compute 4-metric from ADM variables
            // g_00 = -α² + β_i β^i
            float beta_sq = gxx*bx*bx + gyy*by*by + gzz*bz*bz 
                          + 2.0f*(gxy*bx*by + gxz*bx*bz + gyz*by*bz);
            
            g[0][0] = -(alp*alp - beta_sq);
            g[0][1] = gxx*bx + gxy*by + gxz*bz;
            g[0][2] = gxy*bx + gyy*by + gyz*bz;
            g[0][3] = gxz*bx + gyz*by + gzz*bz;
            g[1][0] = g[0][1]; g[2][0] = g[0][2]; g[3][0] = g[0][3];
            g[1][1] = gxx; g[1][2] = gxy; g[1][3] = gxz;
            g[2][1] = gxy; g[2][2] = gyy; g[2][3] = gyz;
            g[3][1] = gxz; g[3][2] = gyz; g[3][3] = gzz;
        };
        
        // Grid spacing
        float dx = hostData.spacing.x;
        float dy = hostData.spacing.y;
        float dz = hostData.spacing.z;
        float h[4] = {1.0f, dx, dy, dz};  // Time is not differentiated
        
        // Compute Christoffel at each grid point
        for (int iz = 0; iz < nz; iz++) {
            for (int iy = 0; iy < ny; iy++) {
                for (int ix = 0; ix < nx; ix++) {
                    int voxelIdx = ix + iy * nx + iz * nx * ny;
                    
                    // Sample metric at current point and neighbors
                    float g[4][4], gp[4][4], gm[4][4];
                    float dg[4][4][4];  // dg[coord][mu][nu]
                    
                    sampleMetric(ix, iy, iz, g);
                    
                    // Compute derivatives ∂g/∂x, ∂g/∂y, ∂g/∂z via central differences
                    // For time: ∂g/∂t = 0 (static metric)
                    for (int mu = 0; mu < 4; mu++)
                        for (int nu = 0; nu < 4; nu++)
                            dg[0][mu][nu] = 0.0f;
                    
                    // ∂/∂x
                    sampleMetric(ix+1, iy, iz, gp);
                    sampleMetric(ix-1, iy, iz, gm);
                    for (int mu = 0; mu < 4; mu++)
                        for (int nu = 0; nu < 4; nu++)
                            dg[1][mu][nu] = (gp[mu][nu] - gm[mu][nu]) / (2.0f * dx);
                    
                    // ∂/∂y
                    sampleMetric(ix, iy+1, iz, gp);
                    sampleMetric(ix, iy-1, iz, gm);
                    for (int mu = 0; mu < 4; mu++)
                        for (int nu = 0; nu < 4; nu++)
                            dg[2][mu][nu] = (gp[mu][nu] - gm[mu][nu]) / (2.0f * dy);
                    
                    // ∂/∂z
                    sampleMetric(ix, iy, iz+1, gp);
                    sampleMetric(ix, iy, iz-1, gm);
                    for (int mu = 0; mu < 4; mu++)
                        for (int nu = 0; nu < 4; nu++)
                            dg[3][mu][nu] = (gp[mu][nu] - gm[mu][nu]) / (2.0f * dz);
                    
                    // =========================================================
                    // Proper 4×4 Matrix Inversion (Cramer's Rule)
                    // =========================================================
                    float g_inv[4][4];
                    
                    // Compute cofactors and determinant
                    float A2323 = g[2][2] * g[3][3] - g[2][3] * g[3][2];
                    float A1323 = g[2][1] * g[3][3] - g[2][3] * g[3][1];
                    float A1223 = g[2][1] * g[3][2] - g[2][2] * g[3][1];
                    float A0323 = g[2][0] * g[3][3] - g[2][3] * g[3][0];
                    float A0223 = g[2][0] * g[3][2] - g[2][2] * g[3][0];
                    float A0123 = g[2][0] * g[3][1] - g[2][1] * g[3][0];
                    float A2313 = g[1][2] * g[3][3] - g[1][3] * g[3][2];
                    float A1313 = g[1][1] * g[3][3] - g[1][3] * g[3][1];
                    float A1213 = g[1][1] * g[3][2] - g[1][2] * g[3][1];
                    float A2312 = g[1][2] * g[2][3] - g[1][3] * g[2][2];
                    float A1312 = g[1][1] * g[2][3] - g[1][3] * g[2][1];
                    float A1212 = g[1][1] * g[2][2] - g[1][2] * g[2][1];
                    float A0313 = g[1][0] * g[3][3] - g[1][3] * g[3][0];
                    float A0213 = g[1][0] * g[3][2] - g[1][2] * g[3][0];
                    float A0312 = g[1][0] * g[2][3] - g[1][3] * g[2][0];
                    float A0212 = g[1][0] * g[2][2] - g[1][2] * g[2][0];
                    float A0113 = g[1][0] * g[3][1] - g[1][1] * g[3][0];
                    float A0112 = g[1][0] * g[2][1] - g[1][1] * g[2][0];

                    float det = g[0][0] * (g[1][1] * A2323 - g[1][2] * A1323 + g[1][3] * A1223)
                              - g[0][1] * (g[1][0] * A2323 - g[1][2] * A0323 + g[1][3] * A0223)
                              + g[0][2] * (g[1][0] * A1323 - g[1][1] * A0323 + g[1][3] * A0123)
                              - g[0][3] * (g[1][0] * A1223 - g[1][1] * A0223 + g[1][2] * A0123);
                    
                    // Handle degenerate case
                    if (std::abs(det) < 1e-10f) {
                        // Fall back to Minkowski inverse
                        for (int i = 0; i < 4; i++)
                            for (int j = 0; j < 4; j++)
                                g_inv[i][j] = (i == j) ? (i == 0 ? -1.0f : 1.0f) : 0.0f;
                    } else {
                        float invDet = 1.0f / det;
                        
                        g_inv[0][0] =  invDet * (g[1][1] * A2323 - g[1][2] * A1323 + g[1][3] * A1223);
                        g_inv[0][1] = -invDet * (g[0][1] * A2323 - g[0][2] * A1323 + g[0][3] * A1223);
                        g_inv[0][2] =  invDet * (g[0][1] * A2313 - g[0][2] * A1313 + g[0][3] * A1213);
                        g_inv[0][3] = -invDet * (g[0][1] * A2312 - g[0][2] * A1312 + g[0][3] * A1212);
                        g_inv[1][0] = -invDet * (g[1][0] * A2323 - g[1][2] * A0323 + g[1][3] * A0223);
                        g_inv[1][1] =  invDet * (g[0][0] * A2323 - g[0][2] * A0323 + g[0][3] * A0223);
                        g_inv[1][2] = -invDet * (g[0][0] * A2313 - g[0][2] * A0313 + g[0][3] * A0213);
                        g_inv[1][3] =  invDet * (g[0][0] * A2312 - g[0][2] * A0312 + g[0][3] * A0212);
                        g_inv[2][0] =  invDet * (g[1][0] * A1323 - g[1][1] * A0323 + g[1][3] * A0123);
                        g_inv[2][1] = -invDet * (g[0][0] * A1323 - g[0][1] * A0323 + g[0][3] * A0123);
                        g_inv[2][2] =  invDet * (g[0][0] * A1313 - g[0][1] * A0313 + g[0][3] * A0113);
                        g_inv[2][3] = -invDet * (g[0][0] * A1312 - g[0][1] * A0312 + g[0][3] * A0112);
                        g_inv[3][0] = -invDet * (g[1][0] * A1223 - g[1][1] * A0223 + g[1][2] * A0123);
                        g_inv[3][1] =  invDet * (g[0][0] * A1223 - g[0][1] * A0223 + g[0][2] * A0123);
                        g_inv[3][2] = -invDet * (g[0][0] * A1213 - g[0][1] * A0213 + g[0][2] * A0113);
                        g_inv[3][3] =  invDet * (g[0][0] * A1212 - g[0][1] * A0212 + g[0][2] * A0112);
                    }
                    
                    // Compute Christoffel: Γ^μ_νρ = 0.5 g^μσ (∂_ν g_σρ + ∂_ρ g_σν - ∂_σ g_νρ)
                    for (int mu = 0; mu < 4; mu++) {
                        for (int nu = 0; nu < 4; nu++) {
                            for (int rho = nu; rho < 4; rho++) {  // Only compute nu <= rho
                                float sum = 0.0f;
                                for (int sigma = 0; sigma < 4; sigma++) {
                                    sum += g_inv[mu][sigma] * (
                                        dg[nu][sigma][rho] +
                                        dg[rho][sigma][nu] -
                                        dg[sigma][nu][rho]
                                    );
                                }
                                
                                int packedIdx = packIdx(nu, rho);
                                int dataIdx = mu * 10 + packedIdx;
                                christoffelData[dataIdx][voxelIdx] = 0.5f * sum;
                            }
                        }
                    }
                }
            }
        }
        
        std::cout << "[OptiX] Uploading Christoffel textures..." << std::endl;
        
        // Upload all 40 components as 3D textures
        // Track arrays and textures for proper cleanup
        for (int mu = 0; mu < 4; mu++) {
            cudaTextureObject_t* targetArray;
            switch (mu) {
                case 0: targetArray = nm.Gamma_t; break;
                case 1: targetArray = nm.Gamma_r; break;
                case 2: targetArray = nm.Gamma_theta; break;
                case 3: targetArray = nm.Gamma_phi; break;
            }
            
            for (int i = 0; i < 10; i++) {
                int idx = mu * 10 + i;
                cudaArray_t arr = nullptr;
                targetArray[i] = create3DTexture(christoffelData[idx].data(), hostData.dims, arr);
                // Track arrays for cleanup
                m_ChristoffelTextures.arrays[idx] = arr;
                m_ChristoffelTextures.textures[idx] = targetArray[i];
            }
        }
        
        m_ChristoffelTextures.initialized = true;
        nm.christoffelLoaded = true;
        std::cout << "[OptiX] Christoffel precomputation complete (40 textures)" << std::endl;
    }
    
    // =========================================================================
    // OptiX AI Denoiser (Phase 4 Performance Optimization)
    // =========================================================================
    bool initializeDenoiser();
    void destroyDenoiser();
    void denoise(float blendFactor = 0.0f);  // 0 = fully denoised, 1 = original
    bool isDenoiserEnabled() const { return m_DenoiserEnabled; }
    void setDenoiserEnabled(bool enabled) { m_DenoiserEnabled = enabled; }
    
    int  m_Width       = 0;
    int  m_Height      = 0;


    // Background texture
    cudaArray_t         m_BackgroundArray  = nullptr;
    cudaTextureObject_t m_BackgroundTexObj = 0;
    
    // Denoiser resources
    OptixDenoiser       m_Denoiser         = nullptr;
    CUdeviceptr         m_DenoiserState    = 0;
    CUdeviceptr         m_DenoiserScratch  = 0;
    CUdeviceptr         m_DenoiserIntensity = 0;  // Average log intensity for HDR
    CUdeviceptr         m_DenoisedBuffer   = 0;   // Output buffer
    size_t              m_DenoiserStateSize = 0;
    size_t              m_DenoiserScratchSize = 0;
    bool                m_DenoiserEnabled  = false;
    bool                m_DenoiserInitialized = false;
    
    // Christoffel symbol textures (for numerical metric cleanup)
    ChristoffelTextureSet m_ChristoffelTextures;
};

//==============================================================================
// Constructor / Destructor
//==============================================================================
OptixRenderer::OptixRenderer() = default;

OptixRenderer::~OptixRenderer() {
    cleanup();
}

//==============================================================================
// Initialize OptiX
//==============================================================================
bool OptixRenderer::initialize(int width, int height) {
    if (m_Initialized) {
        std::cerr << "[OptiX] Already initialized" << std::endl;
        return false;
    }

    m_Width = width;
    m_Height = height;

    try {
        //----------------------------------------------------------------------
        // Initialize CUDA
        //----------------------------------------------------------------------
        CUDA_CHECK(cudaFree(0)); // Force CUDA initialization
        
        int numDevices;
        CUDA_CHECK(cudaGetDeviceCount(&numDevices));
        if (numDevices == 0) {
            std::cerr << "[OptiX] No CUDA-capable devices found" << std::endl;
            return false;
        }

        // Use first CUDA device
        CUDA_CHECK(cudaSetDevice(0));
        
        // Get device properties
        cudaDeviceProp deviceProps;
        CUDA_CHECK(cudaGetDeviceProperties(&deviceProps, 0));
        std::cout << "[OptiX] Using GPU: " << deviceProps.name 
                  << " (SM " << deviceProps.major << "." << deviceProps.minor << ")"
                  << std::endl;

        // Create CUDA stream
        CUDA_CHECK(cudaStreamCreate(&m_Stream));

        //----------------------------------------------------------------------
        // Initialize OptiX
        //----------------------------------------------------------------------
        OPTIX_CHECK(optixInit());

        // Get the current CUDA context (created by the runtime API)
        CUresult cuResult = cuCtxGetCurrent(&m_CudaContext);
        if (cuResult != CUDA_SUCCESS || m_CudaContext == nullptr) {
            std::cerr << "[OptiX] Failed to get CUDA context. ";
            std::cerr << "Attempting to use primary context..." << std::endl;
            
            // Use the primary context (recommended approach for interop with runtime API)
            CUdevice device;
            cuResult = cuDeviceGet(&device, 0);
            if (cuResult != CUDA_SUCCESS) {
                std::cerr << "[OptiX] cuDeviceGet failed" << std::endl;
                throw std::runtime_error("Failed to get CUDA device");
            }
            
            // Retain the primary context (safer than creating a new one)
            cuResult = cuDevicePrimaryCtxRetain(&m_CudaContext, device);
            if (cuResult != CUDA_SUCCESS) {
                std::cerr << "[OptiX] cuDevicePrimaryCtxRetain failed" << std::endl;
                throw std::runtime_error("Failed to retain CUDA primary context");
            }
            cuResult = cuCtxSetCurrent(m_CudaContext);
            if (cuResult != CUDA_SUCCESS) {
                std::cerr << "[OptiX] cuCtxSetCurrent failed" << std::endl;
                throw std::runtime_error("Failed to set CUDA context");
            }
        }

        // Create OptiX device context
        OptixDeviceContextOptions contextOptions = {};
        contextOptions.logCallbackFunction = &optixLogCallback;
        contextOptions.logCallbackLevel = 4; // 1=fatal, 2=error, 3=warning, 4=print
        #ifdef _DEBUG
        contextOptions.validationMode = OPTIX_DEVICE_CONTEXT_VALIDATION_MODE_ALL;
        #endif

        OPTIX_CHECK(optixDeviceContextCreate(m_CudaContext, &contextOptions, &m_OptixContext));

        //----------------------------------------------------------------------
        // Allocate Frame Buffers
        //----------------------------------------------------------------------
        size_t frameBufferSize = width * height * sizeof(float4);
        CUDA_CHECK(cudaMalloc(reinterpret_cast<void**>(&m_FrameBuffer), frameBufferSize));
        CUDA_CHECK(cudaMalloc(reinterpret_cast<void**>(&m_AccumBuffer), frameBufferSize));
        CUDA_CHECK(cudaMemset(reinterpret_cast<void*>(m_FrameBuffer), 0, frameBufferSize));
        CUDA_CHECK(cudaMemset(reinterpret_cast<void*>(m_AccumBuffer), 0, frameBufferSize));

        // Allocate launch params buffer
        CUDA_CHECK(cudaMalloc(reinterpret_cast<void**>(&m_LaunchParamsBuffer), sizeof(LaunchParams)));

        std::cout << "[OptiX] Context initialized successfully" << std::endl;
        std::cout << "[OptiX] Frame buffer: " << width << "x" << height 
                  << " (" << (frameBufferSize / 1024) << " KB)" << std::endl;

        m_Initialized = true;
        return true;

    } catch (const std::exception& e) {
        std::cerr << "[OptiX] Initialization failed: " << e.what() << std::endl;
        cleanup();
        return false;
    }
}

//==============================================================================
// Cleanup
//==============================================================================
void OptixRenderer::cleanup() {
    if (m_CudaGLResource) {
        cudaGraphicsUnregisterResource(m_CudaGLResource);
        m_CudaGLResource = nullptr;
    }

    // Free Raygen Records
    for (auto rec : m_RaygenRecords) {
        if (rec) cudaFree(reinterpret_cast<void*>(rec));
    }
    m_RaygenRecords.clear();
    m_SBT.raygenRecord = 0;

    if (m_SBT.missRecordBase) {
        cudaFree(reinterpret_cast<void*>(m_SBT.missRecordBase));
        m_SBT.missRecordBase = 0;
    }
    if (m_SBT.hitgroupRecordBase) {
        cudaFree(reinterpret_cast<void*>(m_SBT.hitgroupRecordBase));
        m_SBT.hitgroupRecordBase = 0;
    }

    // Destroy Program Groups
    for (auto pg : m_RaygenPGs) {
        if (pg) optixProgramGroupDestroy(pg);
    }
    m_RaygenPGs.clear();

    if (m_MissPG) optixProgramGroupDestroy(m_MissPG);
    if (m_HitGroupPG) optixProgramGroupDestroy(m_HitGroupPG);
    if (m_Pipeline) optixPipelineDestroy(m_Pipeline);
    if (m_Module) optixModuleDestroy(m_Module);
    if (m_OptixContext) optixDeviceContextDestroy(m_OptixContext);

    if (m_FrameBuffer) cudaFree(reinterpret_cast<void*>(m_FrameBuffer));
    if (m_AccumBuffer) cudaFree(reinterpret_cast<void*>(m_AccumBuffer));
    if (m_LaunchParamsBuffer) cudaFree(reinterpret_cast<void*>(m_LaunchParamsBuffer));

    if (m_Stream) cudaStreamDestroy(m_Stream);

    // m_RaygenPG = nullptr; // Vector cleared above
    m_MissPG = nullptr;
    m_HitGroupPG = nullptr;
    m_Pipeline = nullptr;
    m_Module = nullptr;
    m_OptixContext = nullptr;
    m_FrameBuffer = 0;
    m_AccumBuffer = 0;
    m_LaunchParamsBuffer = 0;
    m_Stream = nullptr;
    m_Initialized = false;
    
    // =========================================================================
    // RESOURCE CLEANUP (Fix memory leaks)
    // =========================================================================
    
    // Destroy denoiser resources
    destroyDenoiser();
    
    // Destroy background texture
    if (m_BackgroundTexObj != 0) {
        cudaDestroyTextureObject(m_BackgroundTexObj);
        m_BackgroundTexObj = 0;
    }
    if (m_BackgroundArray != nullptr) {
        cudaFreeArray(m_BackgroundArray);
        m_BackgroundArray = nullptr;
    }
    
    // Destroy numerical metric resources
    m_NumericalMetricResources.cleanup();
    
    // Destroy Christoffel symbol textures
    m_ChristoffelTextures.cleanup();

    std::cout << "[OptiX] Cleanup complete" << std::endl;
}

//==============================================================================
// Create Pipeline
//==============================================================================
bool OptixRenderer::createPipeline(const std::string& ptxPath) {
    if (!m_OptixContext) {
        std::cerr << "[OptiX] Context not initialized" << std::endl;
        return false;
    }

    try {
        //----------------------------------------------------------------------
        // Module compilation options
        //----------------------------------------------------------------------
        OptixModuleCompileOptions moduleCompileOptions = {};
        #ifdef _DEBUG
        moduleCompileOptions.optLevel = OPTIX_COMPILE_OPTIMIZATION_LEVEL_0;
        moduleCompileOptions.debugLevel = OPTIX_COMPILE_DEBUG_LEVEL_FULL;
        #else
        // Use O2 optimization - balances performance with compilation stability
        // Note: O3 can trigger JIT bugs with very complex kernels (35k+ instructions)
        moduleCompileOptions.optLevel = OPTIX_COMPILE_OPTIMIZATION_LEVEL_2;
        moduleCompileOptions.debugLevel = OPTIX_COMPILE_DEBUG_LEVEL_NONE;
        #endif
        moduleCompileOptions.maxRegisterCount = OPTIX_COMPILE_DEFAULT_MAX_REGISTER_COUNT;

        //----------------------------------------------------------------------
        // Pipeline compile options
        //----------------------------------------------------------------------
        m_PipelineCompileOptions.usesMotionBlur = false;
        m_PipelineCompileOptions.traversableGraphFlags = OPTIX_TRAVERSABLE_GRAPH_FLAG_ALLOW_SINGLE_LEVEL_INSTANCING;
        m_PipelineCompileOptions.numPayloadValues = 3;    // RGB color payload
        m_PipelineCompileOptions.numAttributeValues = 2;  // UV attributes
        m_PipelineCompileOptions.exceptionFlags = OPTIX_EXCEPTION_FLAG_NONE;
        m_PipelineCompileOptions.pipelineLaunchParamsVariableName = "params";
        m_PipelineCompileOptions.usesPrimitiveTypeFlags = OPTIX_PRIMITIVE_TYPE_FLAGS_TRIANGLE;

        //----------------------------------------------------------------------
        // Load and compile PTX module
        //----------------------------------------------------------------------
        std::string ptxCode = loadPtxFile(ptxPath);
        
        char log[2048];
        size_t logSize = sizeof(log);

        OPTIX_CHECK(optixModuleCreate(
            m_OptixContext,
            &moduleCompileOptions,
            &m_PipelineCompileOptions,
            ptxCode.c_str(),
            ptxCode.size(),
            log, &logSize,
            &m_Module
        ));

        if (logSize > 1) {
            std::cout << "[OptiX] Module compilation log: " << log << std::endl;
        }

        //----------------------------------------------------------------------
        // Create program groups
        //----------------------------------------------------------------------
        OptixProgramGroupOptions pgOptions = {};

        // Ray generation program(s)
        const char* raygenEntries[] = {
            "__raygen__Minkowski",
            "__raygen__Schwarzschild",
            "__raygen__Kerr",
            "__raygen__ReissnerNordstrom",
            "__raygen__Godel",
            "__raygen__TaubNUT",
            "__raygen__KerrSchild",
            "__raygen__EllisDrainhole",
            "__raygen__Alcubierre",
            "__raygen__DeSitter"
        };
        
        m_RaygenPGs.resize(10);
        for(int i=0; i<10; ++i) {
             OptixProgramGroupDesc raygenPGDesc = {};
             raygenPGDesc.kind = OPTIX_PROGRAM_GROUP_KIND_RAYGEN;
             raygenPGDesc.raygen.module = m_Module;
             raygenPGDesc.raygen.entryFunctionName = raygenEntries[i];
             
             logSize = sizeof(log);
             OPTIX_CHECK(optixProgramGroupCreate(
                m_OptixContext, &raygenPGDesc, 1, &pgOptions,
                log, &logSize, &m_RaygenPGs[i]
             ));
        }

        // Miss program
        OptixProgramGroupDesc missPGDesc = {};
        missPGDesc.kind = OPTIX_PROGRAM_GROUP_KIND_MISS;
        missPGDesc.miss.module = m_Module;
        missPGDesc.miss.entryFunctionName = "__miss__background";

        logSize = sizeof(log);
        OPTIX_CHECK(optixProgramGroupCreate(
            m_OptixContext, &missPGDesc, 1, &pgOptions,
            log, &logSize, &m_MissPG
        ));

        // Hit group program (empty for geodesic raymarching - no geometry)
        OptixProgramGroupDesc hitgroupPGDesc = {};
        hitgroupPGDesc.kind = OPTIX_PROGRAM_GROUP_KIND_HITGROUP;
        hitgroupPGDesc.hitgroup.moduleCH = nullptr;
        hitgroupPGDesc.hitgroup.entryFunctionNameCH = nullptr;
        hitgroupPGDesc.hitgroup.moduleAH = nullptr;
        hitgroupPGDesc.hitgroup.entryFunctionNameAH = nullptr;

        logSize = sizeof(log);
        OPTIX_CHECK(optixProgramGroupCreate(
            m_OptixContext, &hitgroupPGDesc, 1, &pgOptions,
            log, &logSize, &m_HitGroupPG
        ));

        //----------------------------------------------------------------------
        // Create pipeline
        //----------------------------------------------------------------------
        std::vector<OptixProgramGroup> programGroups;
        programGroups.insert(programGroups.end(), m_RaygenPGs.begin(), m_RaygenPGs.end());
        programGroups.push_back(m_MissPG);
        programGroups.push_back(m_HitGroupPG);

        OptixPipelineLinkOptions pipelineLinkOptions = {};
        pipelineLinkOptions.maxTraceDepth = 1; // No recursive tracing for geodesics

        logSize = sizeof(log);
        OPTIX_CHECK(optixPipelineCreate(
            m_OptixContext,
            &m_PipelineCompileOptions,
            &pipelineLinkOptions,
            programGroups.data(), (unsigned int)programGroups.size(),
            log, &logSize,
            &m_Pipeline
        ));

        if (logSize > 1) {
            std::cout << "[OptiX] Pipeline creation log: " << log << std::endl;
        }

        //----------------------------------------------------------------------
        // Set stack sizes
        //----------------------------------------------------------------------
        OPTIX_CHECK(optixPipelineSetStackSize(
            m_Pipeline,
            2 * 1024,  // directCallableStackSizeFromTraversal
            2 * 1024,  // directCallableStackSizeFromState
            2 * 1024,  // continuationStackSize
            2          // maxTraversableGraphDepth (2 for single-level instancing: GAS + IAS)
        ));

        std::cout << "[OptiX] Pipeline created successfully" << std::endl;
        return true;

    } catch (const std::exception& e) {
        std::cerr << "[OptiX] Pipeline creation failed: " << e.what() << std::endl;
        return false;
    }
}

//==============================================================================
// Build Shader Binding Table
//==============================================================================
void OptixRenderer::buildSBT() {
    //--------------------------------------------------------------------------
    // Ray generation records - one for each metric type (must match m_RaygenPGs.size())
    //--------------------------------------------------------------------------
    const size_t numRaygenPrograms = m_RaygenPGs.size();  // Should be 7
    m_RaygenRecords.resize(numRaygenPrograms);
    for(size_t i = 0; i < numRaygenPrograms; ++i) {
        RayGenRecord raygenRecord;
        OPTIX_CHECK(optixSbtRecordPackHeader(m_RaygenPGs[i], &raygenRecord));
        
        CUdeviceptr d_record;
        size_t recSize = sizeof(RayGenRecord);
        CUDA_CHECK(cudaMalloc(reinterpret_cast<void**>(&d_record), recSize));
        CUDA_CHECK(cudaMemcpy(reinterpret_cast<void*>(d_record), &raygenRecord, recSize, cudaMemcpyHostToDevice));
        
        m_RaygenRecords[i] = d_record;
    }
    
    // Default to Minkowski
    m_SBT.raygenRecord = m_RaygenRecords[0];

    //--------------------------------------------------------------------------
    // Miss record
    //--------------------------------------------------------------------------
    MissRecord missRecord;
    OPTIX_CHECK(optixSbtRecordPackHeader(m_MissPG, &missRecord));
    missRecord.data.backgroundColor = make_float3(0.0f, 0.0f, 0.02f);

    CUDA_CHECK(cudaMalloc(reinterpret_cast<void**>(&m_MissRecord), sizeof(MissRecord)));
    CUDA_CHECK(cudaMemcpy(
        reinterpret_cast<void*>(m_MissRecord), &missRecord,
        sizeof(MissRecord), cudaMemcpyHostToDevice
    ));

    //--------------------------------------------------------------------------
    // Hit group record (empty for geodesic raymarching)
    //--------------------------------------------------------------------------
    HitGroupRecord hitgroupRecord;
    memset(&hitgroupRecord, 0, sizeof(HitGroupRecord));
    OPTIX_CHECK(optixSbtRecordPackHeader(m_HitGroupPG, &hitgroupRecord));

    CUDA_CHECK(cudaMalloc(reinterpret_cast<void**>(&m_HitGroupRecord), sizeof(HitGroupRecord)));
    CUDA_CHECK(cudaMemcpy(
        reinterpret_cast<void*>(m_HitGroupRecord), &hitgroupRecord,
        sizeof(HitGroupRecord), cudaMemcpyHostToDevice
    ));

    //--------------------------------------------------------------------------
    // Configure SBT
    //--------------------------------------------------------------------------
    // m_SBT.raygenRecord set to default (0) earlier
    m_SBT.missRecordBase = m_MissRecord;
    m_SBT.missRecordStrideInBytes = sizeof(MissRecord);
    m_SBT.missRecordCount = 1;
    m_SBT.hitgroupRecordBase = m_HitGroupRecord;
    m_SBT.hitgroupRecordStrideInBytes = sizeof(HitGroupRecord);
    m_SBT.hitgroupRecordCount = 1;

    std::cout << "[OptiX] Shader Binding Table built" << std::endl;
}

//==============================================================================
// Update Launch Parameters
//==============================================================================
void OptixRenderer::updateLaunchParams(const LaunchParams& params) {
    LaunchParams hostParams = params;
    hostParams.frameBuffer = reinterpret_cast<float4*>(m_FrameBuffer);
    hostParams.accumBuffer = reinterpret_cast<float4*>(m_AccumBuffer);
    hostParams.frameDimensions = make_int2(m_Width, m_Height);

    CUDA_CHECK(cudaMemcpyAsync(
        reinterpret_cast<void*>(m_LaunchParamsBuffer),
        &hostParams,
        sizeof(LaunchParams),
        cudaMemcpyHostToDevice,
        m_Stream
    ));
}

//==============================================================================
// Launch Ray Tracing
//==============================================================================
void OptixRenderer::launch(const LaunchParams& params) {
    updateLaunchParams(params);

    OPTIX_CHECK(optixLaunch(
        m_Pipeline,
        m_Stream,
        m_LaunchParamsBuffer,
        sizeof(LaunchParams),
        &m_SBT,
        m_Width,
        m_Height,
        1  // depth
    ));

    // Wait for completion
    CUDA_CHECK(cudaStreamSynchronize(m_Stream));
}

//==============================================================================
// Resize Frame Buffers
//==============================================================================
void OptixRenderer::resize(int width, int height) {
    if (width == m_Width && height == m_Height) return;

    m_Width = width;
    m_Height = height;

    // =========================================================================
    // P1 STABILITY: Synchronize before freeing buffers
    // =========================================================================
    // Ensure all pending GPU operations complete before freeing memory
    // This prevents race conditions during resize
    CUDA_CHECK(cudaDeviceSynchronize());

    // Reallocate frame buffers
    if (m_FrameBuffer) cudaFree(reinterpret_cast<void*>(m_FrameBuffer));
    if (m_AccumBuffer) cudaFree(reinterpret_cast<void*>(m_AccumBuffer));

    size_t frameBufferSize = width * height * sizeof(float4);
    CUDA_CHECK(cudaMalloc(reinterpret_cast<void**>(&m_FrameBuffer), frameBufferSize));
    CUDA_CHECK(cudaMalloc(reinterpret_cast<void**>(&m_AccumBuffer), frameBufferSize));
    CUDA_CHECK(cudaMemset(reinterpret_cast<void*>(m_FrameBuffer), 0, frameBufferSize));
    CUDA_CHECK(cudaMemset(reinterpret_cast<void*>(m_AccumBuffer), 0, frameBufferSize));
    
    // =========================================================================
    // DENOISER RE-INITIALIZATION (Fix issue #7)
    // =========================================================================
    // Denoiser buffers and state depend on width/height.
    // Must destroy and re-init when dimensions change.
    if (m_DenoiserInitialized) {
        bool wasEnabled = m_DenoiserEnabled;
        destroyDenoiser();
        if (wasEnabled) {
            initializeDenoiser();
            m_DenoiserEnabled = true;
        }
    }

    std::cout << "[OptiX] Resized to " << width << "x" << height << std::endl;
}

//==============================================================================
// OpenGL Interop
//==============================================================================
void OptixRenderer::registerGLTexture(unsigned int glTexture) {
    if (m_CudaGLResource) {
        unregisterGLTexture();
    }

    CUDA_CHECK(cudaGraphicsGLRegisterImage(
        &m_CudaGLResource,
        glTexture,
        GL_TEXTURE_2D,
        cudaGraphicsRegisterFlagsSurfaceLoadStore
    ));

    m_GLInteropEnabled = true;
    std::cout << "[OptiX] Registered OpenGL texture for interop" << std::endl;
}

void OptixRenderer::unregisterGLTexture() {
    if (m_CudaGLResource) {
        CUDA_CHECK(cudaGraphicsUnregisterResource(m_CudaGLResource));
        m_CudaGLResource = nullptr;
        m_GLInteropEnabled = false;
    }
}

float4* OptixRenderer::mapFrameBuffer() {
    // Copy device buffer to host for CPU access
    size_t bufferSize = static_cast<size_t>(m_Width) * m_Height * 4 * sizeof(float);
    m_HostFrameBuffer.resize(m_Width * m_Height * 4);

    cudaError_t err = cudaMemcpy(m_HostFrameBuffer.data(),
                                  reinterpret_cast<void*>(m_FrameBuffer),
                                  bufferSize,
                                  cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) {
        std::cerr << "[OptiX] Failed to copy frame buffer to host: " << cudaGetErrorString(err) << std::endl;
        return nullptr;
    }

    return reinterpret_cast<float4*>(m_HostFrameBuffer.data());
}

void OptixRenderer::unmapFrameBuffer() {
    // For direct buffer access, nothing to do
    // For texture interop, would need to unmap here
}

//==============================================================================
// Update Display (Copy buffer to GL texture)
//==============================================================================
void OptixRenderer::updateDisplay() {
    if (!m_GLInteropEnabled || !m_CudaGLResource) return;

    // Map OpenGL texture
    CUDA_CHECK(cudaGraphicsMapResources(1, &m_CudaGLResource, m_Stream));

    cudaArray_t textureArray;
    CUDA_CHECK(cudaGraphicsSubResourceGetMappedArray(&textureArray, m_CudaGLResource, 0, 0));

    // Copy linear buffer to texture array
    CUDA_CHECK(cudaMemcpy2DToArrayAsync(
        textureArray,
        0, 0,                           // dst offsets
        reinterpret_cast<void*>(m_FrameBuffer),
        m_Width * sizeof(float4),       // pitch
        m_Width * sizeof(float4),       // width in bytes
        m_Height,                       // height
        cudaMemcpyDeviceToDevice,
        m_Stream
    ));

    // Unmap
    CUDA_CHECK(cudaGraphicsUnmapResources(1, &m_CudaGLResource, m_Stream));
    
    // Ensure copy finishes before GL uses it
    // NOTE: For higher performance with fences, we could avoid full sync, 
    // but this is safe for now.
    CUDA_CHECK(cudaStreamSynchronize(m_Stream));
}

//==============================================================================
// Upload Background Texture
//==============================================================================
bool OptixRenderer::uploadBackgroundTexture(const unsigned char* data, int width, int height) {
    try {
        // Clean up existing texture
        if (m_BackgroundTexObj != 0) {
            cudaDestroyTextureObject(m_BackgroundTexObj);
            m_BackgroundTexObj = 0;
        }
        if (m_BackgroundArray != nullptr) {
            cudaFreeArray(m_BackgroundArray);
            m_BackgroundArray = nullptr;
        }

        // Create CUDA channel format (4 channels, 8-bit per channel)
        cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<uchar4>();

        // Allocate CUDA array
        CUDA_CHECK(cudaMallocArray(&m_BackgroundArray, &channelDesc, width, height));

        // Copy image data to CUDA array
        CUDA_CHECK(cudaMemcpy2DToArray(
            m_BackgroundArray,
            0, 0,                          // dst offset
            data,                          // source data
            width * 4,                     // pitch (4 bytes per pixel)
            width * 4,                     // width in bytes
            height,                        // height
            cudaMemcpyHostToDevice
        ));

        // Set up texture resource descriptor
        cudaResourceDesc resDesc = {};
        resDesc.resType = cudaResourceTypeArray;
        resDesc.res.array.array = m_BackgroundArray;

        // Set up texture descriptor
        cudaTextureDesc texDesc = {};
        texDesc.addressMode[0] = cudaAddressModeWrap;   // U wraps around
        texDesc.addressMode[1] = cudaAddressModeClamp;  // V clamps to edge
        texDesc.filterMode = cudaFilterModeLinear;      // Bilinear filtering
        texDesc.readMode = cudaReadModeNormalizedFloat; // Convert to [0,1] float
        texDesc.normalizedCoords = 1;                   // Use normalized UVs [0,1]

        // Create texture object
        CUDA_CHECK(cudaCreateTextureObject(&m_BackgroundTexObj, &resDesc, &texDesc, nullptr));

        std::cout << "[OptiX] Background texture uploaded: " << width << "x" << height << std::endl;
        return true;
    } catch (const std::exception& e) {
        std::cerr << "[OptiX] Failed to upload background texture: " << e.what() << std::endl;
        return false;
    }
    }
};

//==============================================================================
// Set Metric Type (Kernel Switch)
//==============================================================================
void Sirius::OptixRenderer::setMetricType(int type) {
    if (type >= 0 && type < m_RaygenRecords.size()) {
        m_SBT.raygenRecord = m_RaygenRecords[type];
        // Note: We don't need to rebuild SBT, just updating the host pointer struct is enough?
        // Wait, m_SBT is OptixShaderBindingTable struct passed to optixLaunch.
        // We modify it here, and verify launch uses it.
        // Yes, launch() calls optixLaunch(..., &m_SBT).
    }
}



//==============================================================================
// OptiX AI Denoiser Implementation (Phase 4 Performance)
//==============================================================================
bool Sirius::OptixRenderer::initializeDenoiser() {
    if (m_DenoiserInitialized) return true;
    if (!m_OptixContext) return false;
    
    try {
        // Create denoiser with HDR model (best for ray-traced content)
        OptixDenoiserOptions denoiserOptions = {};
        denoiserOptions.guideAlbedo = 0;   // No albedo guide
        denoiserOptions.guideNormal = 0;   // No normal guide
        
        OPTIX_CHECK(optixDenoiserCreate(
            m_OptixContext,
            OPTIX_DENOISER_MODEL_KIND_HDR,
            &denoiserOptions,
            &m_Denoiser
        ));
        
        // Compute memory requirements
        OptixDenoiserSizes denoiserSizes;
        OPTIX_CHECK(optixDenoiserComputeMemoryResources(
            m_Denoiser,
            m_Width,
            m_Height,
            &denoiserSizes
        ));
        
        m_DenoiserStateSize = denoiserSizes.stateSizeInBytes;
        m_DenoiserScratchSize = denoiserSizes.withoutOverlapScratchSizeInBytes;
        
        // Allocate buffers
        CUDA_CHECK(cudaMalloc(reinterpret_cast<void**>(&m_DenoiserState), m_DenoiserStateSize));
        CUDA_CHECK(cudaMalloc(reinterpret_cast<void**>(&m_DenoiserScratch), m_DenoiserScratchSize));
        CUDA_CHECK(cudaMalloc(reinterpret_cast<void**>(&m_DenoiserIntensity), sizeof(float)));
        CUDA_CHECK(cudaMalloc(reinterpret_cast<void**>(&m_DenoisedBuffer), m_Width * m_Height * sizeof(float4)));
        
        // Setup denoiser state
        OPTIX_CHECK(optixDenoiserSetup(
            m_Denoiser,
            m_Stream,
            m_Width,
            m_Height,
            m_DenoiserState,
            m_DenoiserStateSize,
            m_DenoiserScratch,
            m_DenoiserScratchSize
        ));
        
        m_DenoiserInitialized = true;
        std::cout << "[OptiX] AI Denoiser initialized (HDR, " 
                  << m_DenoiserStateSize / 1024 << " KB state)" << std::endl;
        return true;
        
    } catch (const std::exception& e) {
        std::cerr << "[OptiX] Denoiser init failed: " << e.what() << std::endl;
        destroyDenoiser();
        return false;
    }
}

void Sirius::OptixRenderer::destroyDenoiser() {
    if (m_Denoiser) { optixDenoiserDestroy(m_Denoiser); m_Denoiser = nullptr; }
    if (m_DenoiserState) { cudaFree(reinterpret_cast<void*>(m_DenoiserState)); m_DenoiserState = 0; }
    if (m_DenoiserScratch) { cudaFree(reinterpret_cast<void*>(m_DenoiserScratch)); m_DenoiserScratch = 0; }
    if (m_DenoiserIntensity) { cudaFree(reinterpret_cast<void*>(m_DenoiserIntensity)); m_DenoiserIntensity = 0; }
    if (m_DenoisedBuffer) { cudaFree(reinterpret_cast<void*>(m_DenoisedBuffer)); m_DenoisedBuffer = 0; }
    m_DenoiserInitialized = false;
}

void Sirius::OptixRenderer::denoise(float blendFactor) {
    if (!m_DenoiserEnabled || !m_DenoiserInitialized) return;
    
    // Setup input layer
    OptixImage2D inputLayer = {};
    inputLayer.data = m_FrameBuffer;
    inputLayer.width = m_Width;
    inputLayer.height = m_Height;
    inputLayer.rowStrideInBytes = m_Width * sizeof(float4);
    inputLayer.pixelStrideInBytes = sizeof(float4);
    inputLayer.format = OPTIX_PIXEL_FORMAT_FLOAT4;
    
    // Compute HDR intensity
    OPTIX_CHECK(optixDenoiserComputeIntensity(
        m_Denoiser, m_Stream, &inputLayer,
        m_DenoiserIntensity, m_DenoiserScratch, m_DenoiserScratchSize
    ));
    
    // Setup denoiser layer
    OptixDenoiserLayer layer = {};
    layer.input = inputLayer;
    layer.output = inputLayer;  // Same format
    layer.output.data = m_DenoisedBuffer;
    
    // Empty guide layer (required by API even when not using guides)
    OptixDenoiserGuideLayer guideLayer = {};
    
    OptixDenoiserParams denoiserParams = {};
    denoiserParams.blendFactor = blendFactor;
    denoiserParams.hdrIntensity = m_DenoiserIntensity;
    
    OPTIX_CHECK(optixDenoiserInvoke(
        m_Denoiser, m_Stream, &denoiserParams,
        m_DenoiserState, m_DenoiserStateSize,
        &guideLayer, &layer, 1, 0, 0,
        m_DenoiserScratch, m_DenoiserScratchSize
    ));
    
    // Copy result back
    CUDA_CHECK(cudaMemcpyAsync(
        reinterpret_cast<void*>(m_FrameBuffer),
        reinterpret_cast<void*>(m_DenoisedBuffer),
        m_Width * m_Height * sizeof(float4),
        cudaMemcpyDeviceToDevice, m_Stream
    ));
}

//==============================================================================
// C API for Integration with RDRT001A.cpp
//==============================================================================
extern "C" {
    
// Opaque handle type
typedef void* SiriusOptixHandle;

SiriusOptixHandle sirius_optix_create() {
    return new Sirius::OptixRenderer();
}

void sirius_optix_destroy(SiriusOptixHandle handle) {
    delete static_cast<Sirius::OptixRenderer*>(handle);
}

bool sirius_optix_initialize(SiriusOptixHandle handle, int width, int height) {
    return static_cast<Sirius::OptixRenderer*>(handle)->initialize(width, height);
}

bool sirius_optix_create_pipeline(SiriusOptixHandle handle, const char* ptxPath) {
    auto renderer = static_cast<Sirius::OptixRenderer*>(handle);
    if (!renderer->createPipeline(ptxPath)) return false;
    try {
        renderer->buildSBT();
    } catch (const std::exception& e) {
        std::cerr << "[OptiX] SBT creation failed: " << e.what() << std::endl;
        return false;
    }
    return true;
}

void sirius_optix_launch(SiriusOptixHandle handle, const Sirius::LaunchParams* params) {
    static_cast<Sirius::OptixRenderer*>(handle)->launch(*params);
}

void sirius_optix_update_display(SiriusOptixHandle handle) {
    static_cast<Sirius::OptixRenderer*>(handle)->updateDisplay();
}

void sirius_optix_resize(SiriusOptixHandle handle, int width, int height) {
    static_cast<Sirius::OptixRenderer*>(handle)->resize(width, height);
}

void sirius_optix_set_metric_type(SiriusOptixHandle handle, int type) {
    static_cast<Sirius::OptixRenderer*>(handle)->setMetricType(type);
}

void sirius_optix_cleanup(SiriusOptixHandle handle) {
    static_cast<Sirius::OptixRenderer*>(handle)->cleanup();
}

float* sirius_optix_get_frame_buffer(SiriusOptixHandle handle) {
    return reinterpret_cast<float*>(
        static_cast<Sirius::OptixRenderer*>(handle)->mapFrameBuffer()
    );
}

void sirius_optix_upload_numerical_metric(SiriusOptixHandle handle, 
                                          const Sirius::NumericalMetricHostData* hostData,
                                          Sirius::NumericalMetricData* outDeviceData) {
    if (handle && hostData && outDeviceData) {
        *outDeviceData = static_cast<Sirius::OptixRenderer*>(handle)->uploadNumericalMetric(*hostData);
    }
}



void sirius_optix_register_gl_texture(SiriusOptixHandle handle, unsigned int glTexture) {
    static_cast<Sirius::OptixRenderer*>(handle)->registerGLTexture(glTexture);
}

bool sirius_optix_is_initialized(SiriusOptixHandle handle) {
    return static_cast<Sirius::OptixRenderer*>(handle)->isInitialized();
}

bool sirius_optix_upload_background(SiriusOptixHandle handle, const unsigned char* data, int width, int height) {
    return static_cast<Sirius::OptixRenderer*>(handle)->uploadBackgroundTexture(data, width, height);
}

unsigned long long sirius_optix_get_background_texture(SiriusOptixHandle handle) {
    return static_cast<Sirius::OptixRenderer*>(handle)->getBackgroundTexture();
}

// OptiX AI Denoiser C API
bool sirius_optix_init_denoiser(SiriusOptixHandle handle) {
    return static_cast<Sirius::OptixRenderer*>(handle)->initializeDenoiser();
}

void sirius_optix_denoise(SiriusOptixHandle handle, float blendFactor) {
    static_cast<Sirius::OptixRenderer*>(handle)->denoise(blendFactor);
}

void sirius_optix_set_denoiser_enabled(SiriusOptixHandle handle, bool enabled) {
    static_cast<Sirius::OptixRenderer*>(handle)->setDenoiserEnabled(enabled);
}

} // extern "C"
