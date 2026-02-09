// RDOX001A.h - OptiX Accelerator Backend
// Component ID: RDOX001A (Render/OptiX/Backend)

#pragma once

#include "../Backend/ACIB001A.h"
#include "RDOP006A.h"  // Unified OptiX C API declarations
#include <vector>

namespace Sirius::Acceleration::OptiX {

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
    
    /// @brief Reset accumulation buffer (call when camera changes)
    void resetAccumulation();

private:
    SiriusOptixHandle m_Handle = nullptr;
    int m_Width = 0;
    int m_Height = 0;
    std::string m_LastError;
    unsigned int m_FrameCount = 0;  // Tracks accumulated frames

    // Internal conversion from generic LaunchConfig to OptiX LaunchParams
    Sirius::LaunchParams convertConfig(const LaunchConfig& config);
};

} // namespace Sirius::Acceleration::OptiX
