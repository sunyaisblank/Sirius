#pragma once

#include <cstddef>
#include <string>
#include <string_view>

namespace sirius::render {

struct SessionConfig;
struct VulkanRenderStats;

inline constexpr std::string_view kSceneEvidencePrefix = "[Session] Scene evidence: ";
inline constexpr std::string_view kSourceSceneEvidencePrefix = "[Vulkan] Source scene evidence: ";
inline constexpr std::string_view kVulkanEvidencePrefix = "[Vulkan] Render evidence: ";

// Typed configuration, actual source ownership and measured completion are
// separate records. Native protocol controls bind the writers to the verifier.
[[nodiscard]] std::string SessionSceneEvidenceJson(const SessionConfig& config,
                                                   std::size_t point_star_count);
[[nodiscard]] std::string VulkanRenderEvidenceJson(const SessionConfig& config,
                                                   const VulkanRenderStats& stats);

}  // namespace sirius::render
