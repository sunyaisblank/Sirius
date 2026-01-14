// RDOX001A.h - OptiX Accelerator Backend
// Component ID: RDOX001A (Render/OptiX/Backend)

#pragma once

#include "../Backend/ACIB001A.h"
#include "RDOP003A.h"  // LaunchParams
#include <vector>

namespace Sirius::Acceleration::OptiX {

// Raw C API declarations (avoiding including headers that might conflict)
extern "C" {
    typedef void* SiriusOptixHandle;
    SiriusOptixHandle sirius_optix_create();
    void sirius_optix_destroy(SiriusOptixHandle handle);
    bool sirius_optix_initialize(SiriusOptixHandle handle, int width, int height);
    bool sirius_optix_create_pipeline(SiriusOptixHandle handle, const char* ptxPath);
    void sirius_optix_launch(SiriusOptixHandle handle, const Sirius::LaunchParams* params);
    void sirius_optix_cleanup(SiriusOptixHandle handle);
    float* sirius_optix_get_frame_buffer(SiriusOptixHandle handle);
    bool sirius_optix_is_initialized(SiriusOptixHandle handle);
    bool sirius_optix_upload_background(SiriusOptixHandle handle, const unsigned char* data, int width, int height);
    unsigned long long sirius_optix_get_background_texture(SiriusOptixHandle handle);
}

class OptiXAccelerator : public IAccelerator {
public:
    OptiXAccelerator();
    ~OptiXAccelerator() override;
    
    BackendType getType() const override { return BackendType::OptiX; }
    DeviceCapabilities getCapabilities() const override;
    
    bool initialise(int width, int height) override;
    bool isInitialised() const override;
    
    void launch(const LaunchConfig& config) override;
    float* getFrameBuffer() override;
    size_t getFrameBufferSize() const override { return m_Width * m_Height * 4 * sizeof(float); }
    
    bool uploadBackground(const uint8_t* data, int width, int height) override;
    void synchronise() override;
    void cleanup() override;
    std::string getLastError() const override { return m_LastError; }
    
private:
    SiriusOptixHandle m_Handle = nullptr;
    int m_Width = 0;
    int m_Height = 0;
    std::string m_LastError;
    
    // Internal conversion from generic LaunchConfig to OptiX LaunchParams
    Sirius::LaunchParams convertConfig(const LaunchConfig& config);
};

} // namespace Sirius::Acceleration::OptiX
