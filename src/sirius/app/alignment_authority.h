#pragma once

#include <cstddef>
#include <expected>
#include <filesystem>
#include <string>
#include <vector>

namespace sirius::app {

// Runtime view of the revision-bound receipt generated at configuration. The
// complete JSON is compared with the compiled authority before these summary
// fields are returned.
struct AlignmentAuthority {
    bool release_enforced = false;
    bool qualification_enforced = false;
    bool satisfied = false;
    bool source_tree_clean = false;
    std::string state;
    std::string source_revision;
    std::size_t admitted_domains = 0;
    std::size_t required_domains = 0;
    std::vector<std::string> admitted_domain_ids;
    std::vector<std::string> pending_domain_ids;
    std::vector<std::string> required_domain_ids;
};

struct BuildGateAuthority {
    bool required = false;
    std::string source_revision;
    std::size_t registered_tests = 0;
    std::size_t executed_tests = 0;
    std::size_t verified_product_artifacts = 0;
};

std::expected<AlignmentAuthority, std::string> LoadInstalledAlignmentAuthority();
std::expected<BuildGateAuthority, std::string> LoadInstalledBuildGateAuthority(
    const AlignmentAuthority& alignment);

// Validate a selected volume with the same release authority used by the
// installed loader. Keeping path selection explicit makes the verifier
// independently testable without weakening packaged resource discovery.
std::expected<BuildGateAuthority, std::string> ValidateBuildGateAuthorityForVolume(
    const AlignmentAuthority& alignment, const std::filesystem::path& executable,
    const std::filesystem::path& resource_root);

}  // namespace sirius::app
