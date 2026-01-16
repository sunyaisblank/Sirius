// ACMG001A.h - Multi-GPU Device Manager
// Component ID: ACMG001A
// GPU device enumeration, selection, and multi-device orchestration
//
// Governance: docs/specification.md

#ifndef SIRIUS_RENDER_ACMG001A_H
#define SIRIUS_RENDER_ACMG001A_H

#include <cstdint>
#include <functional>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

#ifdef SIRIUS_HAS_CUDA
#include <cuda_runtime.h>
#endif

namespace sirius::render {

//==============================================================================
// GPU Device Information
//==============================================================================

struct GPUDeviceInfo {
    int deviceId = -1;
    std::string name;
    int smMajor = 0;
    int smMinor = 0;
    size_t totalMemory = 0;      // Total global memory in bytes
    size_t freeMemory = 0;       // Available memory (approximate)
    int multiprocessorCount = 0;
    int maxThreadsPerBlock = 0;
    bool supportsOptiX = false;
    bool supportsPeerAccess = false;
    float computeScore = 0.0f;   // Relative performance metric
};

//==============================================================================
// Tile Assignment
// Represents a render tile assigned to a specific GPU
//==============================================================================

struct TileAssignment {
    int tileIndex = -1;
    int deviceId = -1;
    int x = 0, y = 0;
    int width = 0, height = 0;
};

//==============================================================================
// Multi-GPU Configuration
//==============================================================================

struct MultiGPUConfig {
    bool enabled = false;
    std::vector<int> enabledDevices;      // Empty = use all available
    bool loadBalanceByPerformance = true; // Weight tile distribution by compute score
    bool enablePeerAccess = true;         // Enable direct GPU-to-GPU transfers
    int minTilesPerDevice = 1;            // Minimum tiles before splitting across GPUs
};

//==============================================================================
// GPUDeviceManager
// Singleton for GPU enumeration and selection
//==============================================================================

class GPUDeviceManager {
public:
    static GPUDeviceManager& instance() {
        static GPUDeviceManager mgr;
        return mgr;
    }

    //--------------------------------------------------------------------------
    // Device Enumeration
    //--------------------------------------------------------------------------

    /// Refresh device list (call on startup or after GPU changes)
    bool enumerateDevices() {
#ifdef SIRIUS_HAS_CUDA
        std::lock_guard<std::mutex> lock(m_mutex);
        m_devices.clear();

        int deviceCount = 0;
        cudaError_t err = cudaGetDeviceCount(&deviceCount);
        if (err != cudaSuccess || deviceCount == 0) {
            return false;
        }

        for (int i = 0; i < deviceCount; ++i) {
            cudaDeviceProp props;
            if (cudaGetDeviceProperties(&props, i) != cudaSuccess) continue;

            GPUDeviceInfo info;
            info.deviceId = i;
            info.name = props.name;
            info.smMajor = props.major;
            info.smMinor = props.minor;
            info.totalMemory = props.totalGlobalMem;
            info.multiprocessorCount = props.multiProcessorCount;
            info.maxThreadsPerBlock = props.maxThreadsPerBlock;

            // Check OptiX compatibility (SM 5.0+)
            info.supportsOptiX = (props.major >= 5);

            // Compute performance score (rough estimate)
            // Based on SM count and clock rate
            info.computeScore = static_cast<float>(props.multiProcessorCount) *
                               (props.clockRate / 1000000.0f);

            // Check peer access capability
            info.supportsPeerAccess = (props.unifiedAddressing != 0);

            // Query free memory
            size_t freeMem = 0, totalMem = 0;
            int prevDevice;
            cudaGetDevice(&prevDevice);
            cudaSetDevice(i);
            if (cudaMemGetInfo(&freeMem, &totalMem) == cudaSuccess) {
                info.freeMemory = freeMem;
            }
            cudaSetDevice(prevDevice);

            m_devices.push_back(info);
        }

        return !m_devices.empty();
#else
        return false;
#endif
    }

    /// Get all detected devices
    const std::vector<GPUDeviceInfo>& getDevices() const {
        return m_devices;
    }

    /// Get device count
    size_t deviceCount() const {
        return m_devices.size();
    }

    /// Get device by ID
    const GPUDeviceInfo* getDevice(int deviceId) const {
        for (const auto& dev : m_devices) {
            if (dev.deviceId == deviceId) return &dev;
        }
        return nullptr;
    }

    //--------------------------------------------------------------------------
    // Device Selection
    //--------------------------------------------------------------------------

    /// Get the best device for rendering (highest compute score with OptiX)
    int getBestDevice() const {
        int best = -1;
        float bestScore = -1.0f;
        for (const auto& dev : m_devices) {
            if (dev.supportsOptiX && dev.computeScore > bestScore) {
                bestScore = dev.computeScore;
                best = dev.deviceId;
            }
        }
        return best >= 0 ? best : (m_devices.empty() ? -1 : m_devices[0].deviceId);
    }

    /// Get devices suitable for multi-GPU rendering
    std::vector<int> getMultiGPUDevices(const MultiGPUConfig& config) const {
        std::vector<int> result;

        for (const auto& dev : m_devices) {
            // Skip devices without OptiX support
            if (!dev.supportsOptiX) continue;

            // Check if device is in enabled list (empty = all)
            if (!config.enabledDevices.empty()) {
                bool found = false;
                for (int id : config.enabledDevices) {
                    if (id == dev.deviceId) { found = true; break; }
                }
                if (!found) continue;
            }

            result.push_back(dev.deviceId);
        }

        return result;
    }

    //--------------------------------------------------------------------------
    // Peer Access
    //--------------------------------------------------------------------------

    /// Enable peer access between devices (for direct memory transfers)
    bool enablePeerAccess(int srcDevice, int dstDevice) {
#ifdef SIRIUS_HAS_CUDA
        int canAccess = 0;
        cudaError_t err = cudaDeviceCanAccessPeer(&canAccess, srcDevice, dstDevice);
        if (err != cudaSuccess || !canAccess) return false;

        int prevDevice;
        cudaGetDevice(&prevDevice);
        cudaSetDevice(srcDevice);
        err = cudaDeviceEnablePeerAccess(dstDevice, 0);
        cudaSetDevice(prevDevice);

        return (err == cudaSuccess || err == cudaErrorPeerAccessAlreadyEnabled);
#else
        return false;
#endif
    }

    /// Enable peer access for all device pairs
    void enableAllPeerAccess() {
        for (size_t i = 0; i < m_devices.size(); ++i) {
            for (size_t j = 0; j < m_devices.size(); ++j) {
                if (i != j) {
                    enablePeerAccess(m_devices[i].deviceId, m_devices[j].deviceId);
                }
            }
        }
    }

private:
    GPUDeviceManager() { enumerateDevices(); }
    ~GPUDeviceManager() = default;

    GPUDeviceManager(const GPUDeviceManager&) = delete;
    GPUDeviceManager& operator=(const GPUDeviceManager&) = delete;

    mutable std::mutex m_mutex;
    std::vector<GPUDeviceInfo> m_devices;
};

//==============================================================================
// TileDistributor
// Distributes render tiles across multiple GPUs
//==============================================================================

class TileDistributor {
public:
    explicit TileDistributor(const MultiGPUConfig& config = MultiGPUConfig())
        : m_config(config) {}

    /// Distribute tiles across available GPUs
    /// @param tilesX Number of tiles horizontally
    /// @param tilesY Number of tiles vertically
    /// @param tileWidth Pixel width of each tile
    /// @param tileHeight Pixel height of each tile
    /// @return Vector of tile assignments per device
    std::vector<TileAssignment> distributeTiles(int tilesX, int tilesY,
                                                 int tileWidth, int tileHeight) {
        std::vector<TileAssignment> assignments;
        int totalTiles = tilesX * tilesY;

        // Get available devices
        auto& mgr = GPUDeviceManager::instance();
        auto devices = mgr.getMultiGPUDevices(m_config);

        if (devices.empty()) {
            // Fallback to best single device
            int best = mgr.getBestDevice();
            if (best >= 0) devices.push_back(best);
        }

        if (devices.empty() || totalTiles < m_config.minTilesPerDevice) {
            // Single device mode
            int deviceId = devices.empty() ? 0 : devices[0];
            for (int idx = 0; idx < totalTiles; ++idx) {
                TileAssignment ta;
                ta.tileIndex = idx;
                ta.deviceId = deviceId;
                ta.x = (idx % tilesX) * tileWidth;
                ta.y = (idx / tilesX) * tileHeight;
                ta.width = tileWidth;
                ta.height = tileHeight;
                assignments.push_back(ta);
            }
            return assignments;
        }

        // Calculate device weights for load balancing
        std::vector<float> weights(devices.size(), 1.0f);
        float totalWeight = 0.0f;

        if (m_config.loadBalanceByPerformance) {
            for (size_t i = 0; i < devices.size(); ++i) {
                const auto* dev = mgr.getDevice(devices[i]);
                if (dev) {
                    weights[i] = std::max(1.0f, dev->computeScore);
                }
                totalWeight += weights[i];
            }
        } else {
            totalWeight = static_cast<float>(devices.size());
        }

        // Distribute tiles based on weights
        std::vector<int> tileCounts(devices.size(), 0);
        int assigned = 0;

        for (size_t i = 0; i < devices.size(); ++i) {
            float fraction = weights[i] / totalWeight;
            int count = static_cast<int>(totalTiles * fraction);
            // Ensure at least 1 tile per device
            count = std::max(1, count);
            // Don't exceed remaining tiles
            count = std::min(count, totalTiles - assigned);
            tileCounts[i] = count;
            assigned += count;
        }

        // Distribute any remaining tiles to fastest device
        while (assigned < totalTiles) {
            tileCounts[0]++;
            assigned++;
        }

        // Create assignments
        int tileIdx = 0;
        for (size_t devIdx = 0; devIdx < devices.size(); ++devIdx) {
            for (int t = 0; t < tileCounts[devIdx]; ++t) {
                TileAssignment ta;
                ta.tileIndex = tileIdx;
                ta.deviceId = devices[devIdx];
                ta.x = (tileIdx % tilesX) * tileWidth;
                ta.y = (tileIdx / tilesX) * tileHeight;
                ta.width = tileWidth;
                ta.height = tileHeight;
                assignments.push_back(ta);
                ++tileIdx;
            }
        }

        return assignments;
    }

    /// Get tile assignments grouped by device
    std::vector<std::vector<TileAssignment>> groupByDevice(
        const std::vector<TileAssignment>& assignments) {

        std::vector<std::vector<TileAssignment>> groups;
        std::vector<int> deviceIds;

        for (const auto& ta : assignments) {
            // Find or create group for this device
            int groupIdx = -1;
            for (size_t i = 0; i < deviceIds.size(); ++i) {
                if (deviceIds[i] == ta.deviceId) {
                    groupIdx = static_cast<int>(i);
                    break;
                }
            }
            if (groupIdx < 0) {
                groupIdx = static_cast<int>(deviceIds.size());
                deviceIds.push_back(ta.deviceId);
                groups.push_back({});
            }
            groups[groupIdx].push_back(ta);
        }

        return groups;
    }

private:
    MultiGPUConfig m_config;
};

//==============================================================================
// Multi-GPU Synchronization Utilities
//==============================================================================

/// Synchronize all CUDA devices
inline void synchronizeAllDevices() {
#ifdef SIRIUS_HAS_CUDA
    int deviceCount = 0;
    cudaGetDeviceCount(&deviceCount);
    int prevDevice;
    cudaGetDevice(&prevDevice);

    for (int i = 0; i < deviceCount; ++i) {
        cudaSetDevice(i);
        cudaDeviceSynchronize();
    }

    cudaSetDevice(prevDevice);
#endif
}

/// Execute a function on a specific device
template<typename Func>
void executeOnDevice(int deviceId, Func&& func) {
#ifdef SIRIUS_HAS_CUDA
    int prevDevice;
    cudaGetDevice(&prevDevice);
    cudaSetDevice(deviceId);
    func();
    cudaSetDevice(prevDevice);
#else
    func();
#endif
}

} // namespace sirius::render

#endif // SIRIUS_RENDER_ACMG001A_H
