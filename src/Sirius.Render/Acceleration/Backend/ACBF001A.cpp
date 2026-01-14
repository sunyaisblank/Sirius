// ACBF001A.cpp - Accelerator Backend Factory
// Component ID: ACBF001A (Acceleration/Backend/Factory)

#include "ACIB001A.h"
#include <memory>
#include <vector>

#ifdef SIRIUS_HAS_OPTIX
#include "../OptiX/RDOX001A.h"
#endif

// Future: #include "../CUDA/RDCA001A.h"

namespace Sirius::Acceleration {

std::unique_ptr<IAccelerator> createAccelerator(BackendType type) {
    if (type == BackendType::OptiX) {
#ifdef SIRIUS_HAS_OPTIX
        return std::make_unique<OptiX::OptiXAccelerator>();
#else
        return nullptr;
#endif
    }
    
    // Future backends...
    // if (type == BackendType::CUDA) return std::make_unique<CUDA::CUDAAccelerator>();
    
    return nullptr;
}

std::vector<BackendType> getAvailableBackends() {
    std::vector<BackendType> backends;
    
#ifdef SIRIUS_HAS_OPTIX
    backends.push_back(BackendType::OptiX);
#endif

    // backends.push_back(BackendType::CUDA);
    
    return backends;
}

BackendType getBestBackend() {
    // Priority: OptiX > CUDA > Compute Shader > CPU
    
#ifdef SIRIUS_HAS_OPTIX
    return BackendType::OptiX;
#endif

    // return BackendType::CUDA;
    
    return BackendType::None;
}

} // namespace Sirius::Acceleration
