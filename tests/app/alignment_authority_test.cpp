#include "sirius/app/alignment_authority.h"

#include "sirius/base/resource_locator.h"
#include "sirius/base/sha256.h"

#include <gtest/gtest.h>

#include "support/scoped_temporary_directory.h"
#include <nlohmann/json.hpp>

#include <algorithm>
#include <array>
#include <filesystem>
#include <fstream>
#include <string_view>
#include <vector>

namespace sirius::app {
namespace {

TEST(AlignmentAuthority, CompiledReceiptMatchesTheStagedRuntimeAuthority) {
    const auto authority = LoadInstalledAlignmentAuthority();
    if (!authority) FAIL() << authority.error();
    EXPECT_EQ(authority->required_domains, 8U);
    EXPECT_LE(authority->admitted_domains, authority->required_domains);
    const std::vector<std::string> expected_domains = {
        "macos-moltenvk",
        "macos-native-build",
        "physical-imax-5616x4096",
        "physical-radeon-780m",
        "viewer-native-window-input",
        "windows-native-build",
        "windows-native-vulkan",
        "wsl2-dozen",
    };
    EXPECT_EQ(authority->required_domain_ids, expected_domains);
    EXPECT_EQ(authority->admitted_domain_ids.size(), authority->admitted_domains);
    EXPECT_EQ(authority->pending_domain_ids.size(),
              authority->required_domains - authority->admitted_domains);
    std::vector<std::string> partition = authority->admitted_domain_ids;
    partition.insert(partition.end(), authority->pending_domain_ids.begin(),
                     authority->pending_domain_ids.end());
    std::ranges::sort(partition);
    EXPECT_EQ(partition, expected_domains);
    EXPECT_EQ(authority->satisfied, authority->source_tree_clean &&
                                        authority->admitted_domains == authority->required_domains);
    if (authority->release_enforced) {
        EXPECT_TRUE(authority->satisfied);
    }
    if (!authority->release_enforced && !authority->qualification_enforced) {
        const auto build_gate = LoadInstalledBuildGateAuthority(*authority);
        ASSERT_TRUE(build_gate.has_value()) << build_gate.error();
        EXPECT_FALSE(build_gate->required);
        EXPECT_EQ(build_gate->source_revision, authority->source_revision);
        EXPECT_EQ(build_gate->registered_tests, 0U);
        EXPECT_EQ(build_gate->verified_product_artifacts, 0U);
    } else if (const auto gate_path = base::ResolveResource("model/mandatory_gate.json");
               gate_path.has_value()) {
        auto product_executable = base::ExecutablePath().parent_path() / "sirius";
#if defined(_WIN32)
        product_executable += ".exe";
#endif
        const auto build_gate = ValidateBuildGateAuthorityForVolume(
            *authority, product_executable, gate_path->parent_path().parent_path());
        ASSERT_TRUE(build_gate.has_value()) << build_gate.error();
        EXPECT_TRUE(build_gate->required);
    } else {
        const auto build_gate = LoadInstalledBuildGateAuthority(*authority);
        ASSERT_FALSE(build_gate.has_value());
        EXPECT_NE(build_gate.error().find("missing"), std::string::npos);
    }
}

void WriteFile(const std::filesystem::path& path, std::string_view payload) {
    std::filesystem::create_directories(path.parent_path());
    std::ofstream output(path, std::ios::binary);
    ASSERT_TRUE(output.good());
    output.write(payload.data(), static_cast<std::streamsize>(payload.size()));
    ASSERT_TRUE(output.good());
}

nlohmann::json ArtifactRecord(const std::filesystem::path& path) {
    const auto digest = base::Sha256File(path);
    EXPECT_TRUE(digest.has_value()) << digest.error();
    return {
        {"root", "build"},
        {"path", path.filename().string()},
        {"bytes", std::filesystem::file_size(path)},
        {"sha256", digest.value_or(std::string(64U, '0'))},
    };
}

TEST(BuildGateAuthority, ReleaseReceiptBindsEveryInstalledProductAtInitialisation) {
    const auto staged_model = base::ResolveResource("model/operating_model.json");
    ASSERT_TRUE(staged_model.has_value());
    std::ifstream model_input(*staged_model, std::ios::binary);
    ASSERT_TRUE(model_input.good());
    const std::string model_payload{std::istreambuf_iterator<char>(model_input),
                                    std::istreambuf_iterator<char>()};

    sirius::test::ScopedTemporaryDirectory temporary("sirius-build-gate-runtime");
    const auto root = temporary.path();
    const std::array<std::pair<std::string_view, std::string_view>, 8> products = {{
        {"alignment_receipt", "model/alignment_receipt.json"},
        {"operating_model", "model/operating_model.json"},
        {"starfield", "assets/Starfield.png"},
        {"trace_fp32comp_spv", "kernels/trace_fp32comp.spv"},
        {"trace_fp64_spv", "kernels/trace_fp64.spv"},
        {"trace_spv", "kernels/trace.spv"},
        {"viewer_rdsd003a_fragment", "shaders/RDSD003A.frag"},
        {"viewer_rdsd003a_vertex", "shaders/RDSD003A.vert"},
    }};
    for (const auto& [logical_name, relative_path] : products) {
        WriteFile(root / relative_path,
                  logical_name == "operating_model" ? model_payload : std::string(logical_name));
    }

    nlohmann::json product_records = nlohmann::json::object();
    product_records["sirius"] = ArtifactRecord(base::ExecutablePath());
    for (const auto& [logical_name, relative_path] : products) {
        product_records[logical_name] = ArtifactRecord(root / relative_path);
    }
    constexpr std::array<std::string_view, 7> tested_names = {
        "sirius",
        "sirius_app_tests",
        "sirius_backend_tests",
        "sirius_base_tests",
        "sirius_core_tests",
        "sirius_oracle_tests",
        "sirius_render_tests",
    };
    nlohmann::json tested_records = nlohmann::json::object();
    for (const std::string_view name : tested_names) {
        tested_records[name] = ArtifactRecord(base::ExecutablePath());
    }
    const nlohmann::json evidence_record = ArtifactRecord(root / "assets/Starfield.png");
    nlohmann::json receipt = {
        {"schema_version", 1},
        {"status", "passed"},
        {"alignment_mode", "release"},
        {"source", {{"revision", std::string(40U, 'a')}, {"clean", true}}},
        {"ctest",
         {{"inventory_sha256", std::string(64U, 'b')},
          {"junit", evidence_record},
          {"log", evidence_record},
          {"registered", 700},
          {"executed", 700},
          {"failures", 0},
          {"errors", 0},
          {"skipped", 0}}},
        {"inputs",
         {{"operating_model_sha256", product_records["operating_model"]["sha256"]},
          {"alignment_receipt_sha256", product_records["alignment_receipt"]["sha256"]}}},
        {"tested_artifacts", tested_records},
        {"product_artifacts", product_records},
    };
    const auto gate_path = root / "model/mandatory_gate.json";
    WriteFile(gate_path, receipt.dump(2));
    const AlignmentAuthority alignment{
        .release_enforced = true,
        .satisfied = true,
        .source_tree_clean = true,
        .source_revision = std::string(40U, 'a'),
    };
    auto authority = ValidateBuildGateAuthorityForVolume(alignment, base::ExecutablePath(), root);
    ASSERT_TRUE(authority.has_value()) << authority.error();
    EXPECT_TRUE(authority->required);
    EXPECT_EQ(authority->source_revision, alignment.source_revision);
    EXPECT_EQ(authority->registered_tests, 700U);
    EXPECT_EQ(authority->executed_tests, 700U);
    EXPECT_EQ(authority->verified_product_artifacts, 9U);

    WriteFile(root / "assets/Starfield.png", "starfielD");
    authority = ValidateBuildGateAuthorityForVolume(alignment, base::ExecutablePath(), root);
    ASSERT_FALSE(authority.has_value());
    EXPECT_NE(authority.error().find("hash changed: starfield"), std::string::npos);

    WriteFile(root / "assets/Starfield.png", "starfield");
    receipt["ctest"]["skipped"] = 1;
    WriteFile(gate_path, receipt.dump(2));
    authority = ValidateBuildGateAuthorityForVolume(alignment, base::ExecutablePath(), root);
    ASSERT_FALSE(authority.has_value());
    EXPECT_NE(authority.error().find("zero-skip"), std::string::npos);

    std::filesystem::remove(gate_path);
    authority = ValidateBuildGateAuthorityForVolume(alignment, base::ExecutablePath(), root);
    ASSERT_FALSE(authority.has_value());
    EXPECT_NE(authority.error().find("missing"), std::string::npos);
}

}  // namespace
}  // namespace sirius::app
