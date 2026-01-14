// ACBM001A.h - Backend Manager
// Component ID: ACBM001A (Acceleration/Backend/Manager)
//
// Manages accelerator backend selection and lifecycle.
// Singleton with lazy initialisation.

#pragma once

#include "ACIB001A.h"
#include <memory>
#include <mutex>

namespace Sirius::Acceleration {

//==============================================================================
// BackendManager - Singleton
//==============================================================================
class BackendManager {
public:
    static BackendManager& instance() {
        static BackendManager mgr;
        return mgr;
    }
    
    /// @brief Get or create accelerator with specified backend
    IAccelerator* getAccelerator(BackendType type = BackendType::None) {
        std::lock_guard<std::mutex> lock(m_Mutex);
        
        // Use best available if None specified
        if (type == BackendType::None) {
            type = getBestBackend();
        }
        
        // Return existing if type matches
        if (m_Accelerator && m_Accelerator->getType() == type) {
            return m_Accelerator.get();
        }
        
        // Create new accelerator
        m_Accelerator = createAccelerator(type);
        return m_Accelerator.get();
    }
    
    /// @brief Release current accelerator
    void release() {
        std::lock_guard<std::mutex> lock(m_Mutex);
        if (m_Accelerator) {
            m_Accelerator->cleanup();
            m_Accelerator.reset();
        }
    }
    
    /// @brief Check if accelerator is active
    bool hasAccelerator() const {
        return m_Accelerator != nullptr;
    }
    
    /// @brief Get current backend type
    BackendType getCurrentType() const {
        return m_Accelerator ? m_Accelerator->getType() : BackendType::None;
    }
    
private:
    BackendManager() = default;
    ~BackendManager() { release(); }
    
    BackendManager(const BackendManager&) = delete;
    BackendManager& operator=(const BackendManager&) = delete;
    
    std::unique_ptr<IAccelerator> m_Accelerator;
    std::mutex m_Mutex;
};

} // namespace Sirius::Acceleration
