// SROR001A.h - Offline Renderer
// Component ID: SROR001A
// Orchestrates batch rendering with session, scheduling, and output integration
//
// Governance: docs/specification.md

#ifndef SIRIUS_RENDER_SROR001A_H
#define SIRIUS_RENDER_SROR001A_H

#include "SRRS001A.h"
#include "SRCT001A.h"
#include "SRER001A.h"
#include "SRPG001A.h"
#include "SRCK001A.h"
#include "../Scheduling/SCTM001A.h"
#include "../Output/OUMD001A.h"
#include "../Output/OUDR001A.h"
#include "../Output/OUIB001A.h"

#include <functional>
#include <memory>
#include <thread>
#include <chrono>
#include <atomic>

namespace sirius::render {

//==============================================================================
// OfflineRenderConfig
// Extended configuration for offline/batch rendering
//==============================================================================

struct OfflineRenderConfig : public RenderConfig {
    // Tile ordering
    TileOrder tileOrder = TileOrder::SPIRAL_CENTER;

    // Output settings
    std::string outputDriver = "memory";  // "null", "memory", "exr", "png"
    std::string outputPath = "output";

    // Progress settings
    bool showProgress = true;
    int progressUpdateInterval = 100;  // ms

    // Threading
    int numThreads = 1;  // For future parallel tile rendering
};

//==============================================================================
// OfflineRenderResult
// Result of a completed offline render
//==============================================================================

struct OfflineRenderResult {
    bool success = false;
    int tilesCompleted = 0;
    int tilesTotal = 0;
    double elapsedSeconds = 0.0;
    std::vector<RenderError> errors;
    std::string outputFile;

    // Access to rendered image (if memory output driver)
    std::shared_ptr<MemoryOutputDriver> imageData;
};

//==============================================================================
// TileRenderCallback
// Called to render a single tile - user provides implementation
//==============================================================================

using TileRenderCallback = std::function<bool(
    const TileInfo& tile,           // Tile to render
    std::vector<float>& pixels,     // Output RGBA buffer
    const RenderConfig& config      // Render configuration
)>;

//==============================================================================
// OfflineRenderer
// Main orchestrator for offline/batch rendering
//==============================================================================

class OfflineRenderer {
public:
    using ProgressCallback = std::function<void(const RenderProgress&)>;
    using ErrorCallback = std::function<void(const RenderError&)>;

    OfflineRenderer() = default;
    ~OfflineRenderer() { cancel(); }

    //--------------------------------------------------------------------------
    // Configuration
    //--------------------------------------------------------------------------

    void configure(const OfflineRenderConfig& config) {
        m_config = config;
        m_session = std::make_unique<RenderSession>(config);
    }

    void setTileRenderCallback(TileRenderCallback callback) {
        m_tileRenderCallback = std::move(callback);
    }

    void setProgressCallback(ProgressCallback callback) {
        m_progressCallback = std::move(callback);
    }

    void setErrorCallback(ErrorCallback callback) {
        m_errorCallback = std::move(callback);
    }

    //--------------------------------------------------------------------------
    // Execution
    //--------------------------------------------------------------------------

    /// Start rendering (blocking)
    OfflineRenderResult render() {
        if (!m_tileRenderCallback) {
            return makeErrorResult("No tile render callback set");
        }

        // Reset state
        m_cancelToken.reset();
        m_errors.clear();
        m_eta.reset();

        // Create output driver
        auto outputDriver = createOutputDriver();
        if (!outputDriver) {
            return makeErrorResult("Failed to create output driver");
        }

        // Configure and start output
        std::vector<OutputPass> passes = {OutputPass::COMBINED};
        if (!outputDriver->configure(m_config.imageWidth, m_config.imageHeight, passes)) {
            return makeErrorResult("Failed to configure output driver");
        }
        if (!outputDriver->start()) {
            return makeErrorResult("Failed to start output driver");
        }

        // Generate tile order
        int tilesX = (m_config.imageWidth + m_config.tileSize - 1) / m_config.tileSize;
        int tilesY = (m_config.imageHeight + m_config.tileSize - 1) / m_config.tileSize;
        auto tileOrder = TileScheduler::generateOrder(tilesX, tilesY, m_config.tileOrder);

        // Start session
        m_session->start();
        auto startTime = std::chrono::steady_clock::now();

        // Render tiles
        int tilesCompleted = 0;
        int tilesTotal = static_cast<int>(tileOrder.size());

        for (int tileIdx : tileOrder) {
            if (m_cancelToken.isCancelled()) break;

            const auto& tileInfo = m_session->tiles()[tileIdx];

            // Allocate pixel buffer
            std::vector<float> pixels(tileInfo.width * tileInfo.height * 4, 0.0f);

            // Render tile
            bool success = false;
            try {
                success = m_tileRenderCallback(tileInfo, pixels, m_config);
            } catch (const std::exception& e) {
                m_errors.add(RenderError::error("Tile " + std::to_string(tileIdx) +
                                                " failed: " + e.what()));
                if (m_errorCallback) {
                    m_errorCallback(m_errors.all().back());
                }
            }

            if (success) {
                // Write tile to output
                OutputTile outputTile;
                outputTile.x = tileInfo.x;
                outputTile.y = tileInfo.y;
                outputTile.width = tileInfo.width;
                outputTile.height = tileInfo.height;
                outputTile.pass = OutputPass::COMBINED;
                outputTile.pixels = std::move(pixels);

                if (!outputDriver->writeTile(outputTile)) {
                    m_errors.add(RenderError::warning("Failed to write tile " +
                                                      std::to_string(tileIdx)));
                }

                // Mark tile completed
                m_session->markTileCompleted(tileIdx);
                tilesCompleted++;

                // Update ETA
                auto now = std::chrono::steady_clock::now();
                double elapsed = std::chrono::duration<double>(now - startTime).count();
                m_eta.recordTileCompletion(elapsed);

                // Progress callback
                if (m_progressCallback) {
                    RenderProgress progress;
                    progress.tilesCompleted = tilesCompleted;
                    progress.tilesTotal = tilesTotal;
                    progress.percentComplete = 100.0 * tilesCompleted / tilesTotal;
                    progress.elapsedSeconds = elapsed;
                    progress.estimatedRemaining = m_eta.estimateRemaining(tilesCompleted, tilesTotal);
                    m_progressCallback(progress);
                }
            }
        }

        // Finish output
        auto endTime = std::chrono::steady_clock::now();
        double totalTime = std::chrono::duration<double>(endTime - startTime).count();

        OutputMetadata metadata;
        metadata.width = m_config.imageWidth;
        metadata.height = m_config.imageHeight;
        metadata.tilesRendered = tilesCompleted;
        metadata.samplesPerPixel = m_config.samplesPerPixel;
        metadata.renderTimeSeconds = totalTime;
        metadata.metricName = m_config.metricType;
        metadata.blackHoleSpin = m_config.blackHoleSpin;

        outputDriver->finish(metadata);

        // Build result
        OfflineRenderResult result;
        result.success = !m_cancelToken.isCancelled() && !m_errors.hasFatalError();
        result.tilesCompleted = tilesCompleted;
        result.tilesTotal = tilesTotal;
        result.elapsedSeconds = totalTime;
        result.errors = m_errors.all();
        result.outputFile = m_config.outputPath;

        // Store image data if memory driver
        auto memDriver = std::dynamic_pointer_cast<MemoryOutputDriver>(outputDriver);
        if (memDriver) {
            result.imageData = memDriver;
        }

        return result;
    }

    /// Start rendering asynchronously
    void renderAsync() {
        m_renderThread = std::thread([this]() {
            m_lastResult = render();
        });
    }

    /// Wait for async render to complete
    OfflineRenderResult waitForResult() {
        if (m_renderThread.joinable()) {
            m_renderThread.join();
        }
        return m_lastResult;
    }

    /// Cancel rendering
    void cancel() {
        m_cancelToken.request();
        if (m_renderThread.joinable()) {
            m_renderThread.join();
        }
    }

    //--------------------------------------------------------------------------
    // State Access
    //--------------------------------------------------------------------------

    bool isCancelled() const { return m_cancelToken.isCancelled(); }
    bool hasErrors() const { return m_errors.hasErrors(); }
    std::vector<RenderError> errors() const { return m_errors.all(); }

    const RenderSession* session() const { return m_session.get(); }

private:
    std::shared_ptr<IOutputDriver> createOutputDriver() {
        if (m_config.outputDriver == "null") {
            return std::make_shared<NullOutputDriver>();
        } else if (m_config.outputDriver == "memory") {
            return std::make_shared<MemoryOutputDriver>();
        } else if (m_config.outputDriver == "exr") {
            return std::make_shared<OpenEXROutputDriver>(m_config.outputPath);
        }
        // Default to memory
        return std::make_shared<MemoryOutputDriver>();
    }

    OfflineRenderResult makeErrorResult(const std::string& message) {
        OfflineRenderResult result;
        result.success = false;
        result.errors.push_back(RenderError::fatal(message));
        return result;
    }

    OfflineRenderConfig m_config;
    std::unique_ptr<RenderSession> m_session;
    CancellationToken m_cancelToken;
    ErrorAccumulator m_errors;
    ETACalculator m_eta;

    TileRenderCallback m_tileRenderCallback;
    ProgressCallback m_progressCallback;
    ErrorCallback m_errorCallback;

    std::thread m_renderThread;
    OfflineRenderResult m_lastResult;
};

} // namespace sirius::render

#endif // SIRIUS_RENDER_SROR001A_H
