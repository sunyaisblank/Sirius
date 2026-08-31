#include "sirius/app/alignment_authority.h"

#include "sirius/base/operating_model_embedded.h"
#include "sirius/base/resource_locator.h"
#include "sirius/base/sha256.h"

#include <nlohmann/json.hpp>

#include <algorithm>
#include <array>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <initializer_list>
#include <iterator>
#include <limits>
#include <set>
#include <span>
#include <string_view>
#include <utility>

namespace sirius::app {
namespace {

static_assert(base::kReleaseAlignmentEnforced ==
                  (static_cast<bool>(SIRIUS_RELEASE_RESOURCE_LOCKED) &&
                   !base::kQualificationAlignmentEnforced),
              "release alignment and packaged-resource authority must agree");
static_assert(!base::kQualificationAlignmentEnforced || SIRIUS_RELEASE_RESOURCE_LOCKED == 1,
              "qualification evidence must use packaged-resource authority");

constexpr std::array<std::string_view, 8> kRequiredDomains = {
    "macos-moltenvk",
    "macos-native-build",
    "physical-imax-5616x4096",
    "physical-radeon-780m",
    "viewer-native-window-input",
    "windows-native-build",
    "windows-native-vulkan",
    "wsl2-dozen",
};

constexpr std::array<std::string_view, 7> kTestedArtifactNames = {
    "sirius",
    "sirius_app_tests",
    "sirius_backend_tests",
    "sirius_base_tests",
    "sirius_core_tests",
    "sirius_oracle_tests",
    "sirius_render_tests",
};

constexpr std::array<std::pair<std::string_view, std::string_view>, 8> kResourceProducts = {{
    {"alignment_receipt", "model/alignment_receipt.json"},
    {"operating_model", "model/operating_model.json"},
    {"starfield", "assets/Starfield.png"},
    {"trace_fp32comp_spv", "kernels/trace_fp32comp.spv"},
    {"trace_fp64_spv", "kernels/trace_fp64.spv"},
    {"trace_spv", "kernels/trace.spv"},
    {"viewer_rdsd003a_fragment", "shaders/RDSD003A.frag"},
    {"viewer_rdsd003a_vertex", "shaders/RDSD003A.vert"},
}};

constexpr std::array<std::string_view, 9> kProductArtifactNames = {
    "alignment_receipt", "operating_model",          "sirius",
    "starfield",         "trace_fp32comp_spv",       "trace_fp64_spv",
    "trace_spv",         "viewer_rdsd003a_fragment", "viewer_rdsd003a_vertex",
};

bool HasExactFields(const nlohmann::json& object, std::initializer_list<std::string_view> fields) {
    if (!object.is_object() || object.size() != fields.size()) return false;
    return std::ranges::all_of(
        fields, [&object](std::string_view field) { return object.contains(field); });
}

bool NonEmptyString(const nlohmann::json& object, std::string_view field) {
    const auto iterator = object.find(field);
    return iterator != object.end() && iterator->is_string() &&
           !iterator->get_ref<const std::string&>().empty();
}

bool IsLowerHexRevision(std::string_view revision) {
    return revision.size() >= 40 && revision.size() <= 64 &&
           std::ranges::all_of(revision, [](char character) {
               return (character >= '0' && character <= '9') ||
                      (character >= 'a' && character <= 'f');
           });
}

bool IsLowerHexDigest(std::string_view digest) {
    return digest.size() == 64U && std::ranges::all_of(digest, [](char character) {
               return (character >= '0' && character <= '9') ||
                      (character >= 'a' && character <= 'f');
           });
}

bool IsArtifactRecord(const nlohmann::json& record) {
    if (!HasExactFields(record, {"root", "path", "bytes", "sha256"}) ||
        !NonEmptyString(record, "root") || !NonEmptyString(record, "path") ||
        !record["bytes"].is_number_unsigned() || record["bytes"].get<std::uint64_t>() == 0U ||
        !NonEmptyString(record, "sha256") ||
        !IsLowerHexDigest(record["sha256"].get_ref<const std::string&>())) {
        return false;
    }
    const std::string& root = record["root"].get_ref<const std::string&>();
    return root == "source" || root == "build";
}

bool HasExactArtifactNames(const nlohmann::json& records,
                           std::span<const std::string_view> expected) {
    if (!records.is_object() || records.size() != expected.size()) return false;
    return std::ranges::all_of(expected, [&records](std::string_view name) {
        const auto record = records.find(name);
        return record != records.end() && IsArtifactRecord(*record);
    });
}

std::expected<void, std::string> VerifyArtifact(const nlohmann::json& record,
                                                const std::filesystem::path& path,
                                                std::string_view logical_name) {
    std::error_code error;
    if (!std::filesystem::is_regular_file(path, error) || error) {
        return std::unexpected("build-gate product is missing: " + std::string(logical_name));
    }
    const std::uintmax_t size = std::filesystem::file_size(path, error);
    if (error || size != record["bytes"].get<std::uint64_t>()) {
        return std::unexpected("build-gate product size changed: " + std::string(logical_name));
    }
    const auto digest = base::Sha256File(path);
    if (!digest) {
        return std::unexpected("build-gate product could not be hashed: " + digest.error());
    }
    if (*digest != record["sha256"].get_ref<const std::string&>()) {
        return std::unexpected("build-gate product hash changed: " + std::string(logical_name));
    }
    return {};
}

bool ReadStringSet(const nlohmann::json& value, std::set<std::string>& result) {
    if (!value.is_array()) return false;
    for (const auto& item : value) {
        if (!item.is_string() || item.get_ref<const std::string&>().empty() ||
            !result.insert(item.get<std::string>()).second) {
            return false;
        }
    }
    return true;
}

std::expected<AlignmentAuthority, std::string> ValidateReceipt(const nlohmann::json& receipt) {
    if (!receipt.is_object() || receipt.value("schema_version", 0) != 1 ||
        receipt.value("method", "") != "revision-bound-attestation-set") {
        return std::unexpected("alignment receipt has an unsupported schema or method");
    }
    if (!NonEmptyString(receipt, "alignment_mode") || !NonEmptyString(receipt, "source_revision") ||
        !receipt.contains("source_tree_clean") || !receipt["source_tree_clean"].is_boolean()) {
        return std::unexpected("alignment receipt identity is incomplete");
    }
    const std::string& mode = receipt["alignment_mode"].get_ref<const std::string&>();
    if (mode != "development" && mode != "qualification" && mode != "release") {
        return std::unexpected("alignment receipt has an unknown mode");
    }
    if ((mode == "release") != base::kReleaseAlignmentEnforced ||
        (mode == "qualification") != base::kQualificationAlignmentEnforced) {
        return std::unexpected("installed alignment mode differs from the compiled mode");
    }
    const std::string& revision = receipt["source_revision"].get_ref<const std::string&>();
    if (!IsLowerHexRevision(revision)) {
        return std::unexpected("alignment receipt has an invalid source revision");
    }

    std::set<std::string> required;
    std::set<std::string> admitted;
    std::set<std::string> pending;
    if (!receipt.contains("required_domains") || !receipt.contains("admitted_domains") ||
        !receipt.contains("pending_domains") ||
        !ReadStringSet(receipt["required_domains"], required) ||
        !ReadStringSet(receipt["admitted_domains"], admitted) ||
        !ReadStringSet(receipt["pending_domains"], pending)) {
        return std::unexpected("alignment receipt domain sets are malformed");
    }
    const std::set<std::string> expected(kRequiredDomains.begin(), kRequiredDomains.end());
    if (required != expected) {
        return std::unexpected("alignment receipt does not name the exact required domains");
    }
    const auto operating_model = receipt.find("operating_model");
    std::set<std::string> model_domains;
    if (operating_model == receipt.end() || !operating_model->is_object() ||
        operating_model->value("schema_version", 0) != 3 ||
        !NonEmptyString(*operating_model, "sha256") ||
        (*operating_model)["sha256"].get_ref<const std::string&>() !=
            base::kEmbeddedOperatingModelSha256 ||
        !operating_model->contains("external_domains") ||
        !ReadStringSet((*operating_model)["external_domains"], model_domains) ||
        model_domains != expected) {
        return std::unexpected(
            "alignment receipt is not bound to the compiled operating-model authority");
    }
    std::set<std::string> partition = admitted;
    partition.insert(pending.begin(), pending.end());
    if (partition != expected || admitted.size() + pending.size() != expected.size()) {
        return std::unexpected(
            "alignment receipt admitted/pending domains do not partition the ideal");
    }

    const auto records = receipt.find("records");
    if (records == receipt.end() || !records->is_array()) {
        return std::unexpected("alignment receipt record ledger is malformed");
    }
    std::set<std::string> recorded_domains;
    for (const auto& record : *records) {
        if (!record.is_object() || !NonEmptyString(record, "completed_utc") ||
            !NonEmptyString(record, "source_revision") ||
            record["source_revision"].get_ref<const std::string&>() != revision ||
            !NonEmptyString(record, "sha256") ||
            record["sha256"].get_ref<const std::string&>().size() != 64) {
            return std::unexpected("an alignment receipt record is malformed");
        }
        const std::string& digest = record["sha256"].get_ref<const std::string&>();
        if (!std::ranges::all_of(digest, [](char character) {
                return (character >= '0' && character <= '9') ||
                       (character >= 'a' && character <= 'f');
            })) {
            return std::unexpected("an alignment receipt record has an invalid digest");
        }
        std::set<std::string> domains;
        if (!record.contains("domains") || !ReadStringSet(record["domains"], domains) ||
            domains.empty()) {
            return std::unexpected("an alignment receipt record has malformed domains");
        }
        for (const std::string& domain : domains) {
            if (!expected.contains(domain) || !recorded_domains.insert(domain).second) {
                return std::unexpected("alignment evidence has an unknown or duplicate authority");
            }
        }
    }
    if (recorded_domains != admitted) {
        return std::unexpected("alignment record ledger does not equal the admitted domain set");
    }

    const auto ideal = receipt.find("ultimate_ideal");
    if (ideal == receipt.end() || !ideal->is_object() || ideal->value("required", false) != true ||
        !ideal->contains("satisfied") || !(*ideal)["satisfied"].is_boolean() ||
        !NonEmptyString(*ideal, "state")) {
        return std::unexpected("alignment receipt has no explicit ultimate-ideal state");
    }
    const bool clean = receipt["source_tree_clean"].get<bool>();
    const bool satisfied = (*ideal)["satisfied"].get<bool>();
    const std::string& state = (*ideal)["state"].get_ref<const std::string&>();
    const bool derived_satisfied = clean && pending.empty();
    if (satisfied != derived_satisfied ||
        state != (derived_satisfied ? "aligned" : "attestation_required") ||
        satisfied != base::kEmbeddedAlignmentSatisfied) {
        return std::unexpected("alignment receipt state contradicts its evidence");
    }
    if constexpr (base::kReleaseAlignmentEnforced) {
        if (!satisfied) {
            return std::unexpected(
                "release initialisation requires the ultimate ideal to be aligned");
        }
    }

    return AlignmentAuthority{
        .release_enforced = base::kReleaseAlignmentEnforced,
        .qualification_enforced = base::kQualificationAlignmentEnforced,
        .satisfied = satisfied,
        .source_tree_clean = clean,
        .state = state,
        .source_revision = revision,
        .admitted_domains = admitted.size(),
        .required_domains = required.size(),
        .admitted_domain_ids = {admitted.begin(), admitted.end()},
        .pending_domain_ids = {pending.begin(), pending.end()},
        .required_domain_ids = {required.begin(), required.end()},
    };
}

std::expected<BuildGateAuthority, std::string> ValidateBuildGateReceipt(
    const nlohmann::json& receipt, const AlignmentAuthority& alignment,
    const std::filesystem::path& executable, const std::filesystem::path& resource_root) {
    if (!HasExactFields(receipt, {"schema_version", "status", "alignment_mode", "source", "ctest",
                                  "inputs", "tested_artifacts", "product_artifacts"}) ||
        receipt.value("schema_version", 0) != 1 || receipt.value("status", "") != "passed" ||
        receipt.value("alignment_mode", "") !=
            (alignment.release_enforced ? "release" : "qualification")) {
        return std::unexpected("Mandatory build-gate receipt has an unsupported schema or state");
    }
    const auto& source = receipt["source"];
    if (!HasExactFields(source, {"revision", "clean"}) || !NonEmptyString(source, "revision") ||
        source["revision"].get_ref<const std::string&>() != alignment.source_revision ||
        !source["clean"].is_boolean() || !source["clean"].get<bool>()) {
        return std::unexpected(
            "Mandatory build-gate receipt is not bound to the aligned clean revision");
    }

    const auto& ctest = receipt["ctest"];
    if (!HasExactFields(ctest, {"inventory_sha256", "junit", "log", "registered", "executed",
                                "failures", "errors", "skipped"}) ||
        !NonEmptyString(ctest, "inventory_sha256") ||
        !IsLowerHexDigest(ctest["inventory_sha256"].get_ref<const std::string&>()) ||
        !IsArtifactRecord(ctest["junit"]) || !IsArtifactRecord(ctest["log"])) {
        return std::unexpected("Mandatory build-gate CTest authority is malformed");
    }
    for (const std::string_view field :
         {"registered", "executed", "failures", "errors", "skipped"}) {
        if (!ctest[field].is_number_unsigned()) {
            return std::unexpected("Mandatory build-gate CTest counts are malformed");
        }
    }
    const std::uint64_t registered = ctest["registered"].get<std::uint64_t>();
    const std::uint64_t executed = ctest["executed"].get<std::uint64_t>();
    if (registered < 700U || executed != registered ||
        ctest["failures"].get<std::uint64_t>() != 0U ||
        ctest["errors"].get<std::uint64_t>() != 0U || ctest["skipped"].get<std::uint64_t>() != 0U ||
        registered > std::numeric_limits<std::size_t>::max()) {
        return std::unexpected(
            "Mandatory build-gate receipt is not a complete zero-skip test estate");
    }

    const auto& tested = receipt["tested_artifacts"];
    const auto& products = receipt["product_artifacts"];
    if (!HasExactArtifactNames(tested, kTestedArtifactNames) ||
        !HasExactArtifactNames(products, kProductArtifactNames)) {
        return std::unexpected(
            "Mandatory build-gate receipt does not bind the complete test/product topology");
    }
    const auto& inputs = receipt["inputs"];
    if (!HasExactFields(inputs, {"operating_model_sha256", "alignment_receipt_sha256"}) ||
        !NonEmptyString(inputs, "operating_model_sha256") ||
        !NonEmptyString(inputs, "alignment_receipt_sha256") ||
        inputs["operating_model_sha256"] != products["operating_model"]["sha256"] ||
        inputs["alignment_receipt_sha256"] != products["alignment_receipt"]["sha256"] ||
        products["operating_model"]["sha256"].get_ref<const std::string&>() !=
            base::kEmbeddedOperatingModelSha256) {
        return std::unexpected(
            "Mandatory build-gate inputs differ from the compiled alignment/model authority");
    }

    if (executable.empty()) {
        return std::unexpected(
            "running executable path is unavailable for build-gate verification");
    }
    if (auto verified = VerifyArtifact(products["sirius"], executable, "sirius"); !verified) {
        return std::unexpected(verified.error());
    }
    for (const auto& [logical_name, relative_path] : kResourceProducts) {
        const auto path = base::ResolveResourceFromRoot(resource_root, relative_path);
        if (!path) {
            return std::unexpected("build-gate product is missing: " + std::string(logical_name));
        }
        if (auto verified = VerifyArtifact(products[logical_name], *path, logical_name);
            !verified) {
            return std::unexpected(verified.error());
        }
    }

    return BuildGateAuthority{
        .required = true,
        .source_revision = alignment.source_revision,
        .registered_tests = static_cast<std::size_t>(registered),
        .executed_tests = static_cast<std::size_t>(executed),
        .verified_product_artifacts = products.size(),
    };
}

}  // namespace

std::expected<AlignmentAuthority, std::string> LoadInstalledAlignmentAuthority() {
    const auto path = base::ResolveResource("model/alignment_receipt.json");
    if (!path) return std::unexpected("installed alignment receipt is missing");
    try {
        std::ifstream input(*path, std::ios::binary);
        if (!input) return std::unexpected("installed alignment receipt is unreadable");
        std::string payload{std::istreambuf_iterator<char>(input),
                            std::istreambuf_iterator<char>()};
        const nlohmann::json installed = nlohmann::json::parse(payload);

        std::string authority_payload;
        authority_payload.reserve(base::kEmbeddedAlignmentReceiptSize);
        for (const std::string_view chunk : base::kEmbeddedAlignmentReceiptChunks) {
            authority_payload.append(chunk);
        }
        const nlohmann::json authority = nlohmann::json::parse(authority_payload);
        if (installed != authority) {
            return std::unexpected(
                "installed alignment receipt does not match the compiled authority");
        }
        return ValidateReceipt(installed);
    } catch (const std::exception& error) {
        return std::unexpected("installed alignment receipt is invalid: " +
                               std::string(error.what()));
    }
}

std::expected<BuildGateAuthority, std::string> LoadInstalledBuildGateAuthority(
    const AlignmentAuthority& alignment) {
    if (!alignment.release_enforced && !alignment.qualification_enforced) {
        return BuildGateAuthority{
            .required = false,
            .source_revision = alignment.source_revision,
        };
    }
    const auto path = base::ResolveResource("model/mandatory_gate.json");
    if (!path) return std::unexpected("installed Mandatory build-gate receipt is missing");
    return ValidateBuildGateAuthorityForVolume(alignment, base::ExecutablePath(),
                                               path->parent_path().parent_path());
}

std::expected<BuildGateAuthority, std::string> ValidateBuildGateAuthorityForVolume(
    const AlignmentAuthority& alignment, const std::filesystem::path& executable,
    const std::filesystem::path& resource_root) {
    if ((!alignment.release_enforced && !alignment.qualification_enforced) ||
        (alignment.release_enforced && !alignment.satisfied) || !alignment.source_tree_clean) {
        return std::unexpected(
            "strict build-gate verification requires clean qualification authority or "
            "satisfied clean release authority");
    }
    const auto path = base::ResolveResourceFromRoot(resource_root, "model/mandatory_gate.json");
    if (!path) return std::unexpected("installed Mandatory build-gate receipt is missing");
    try {
        std::ifstream input(*path, std::ios::binary);
        if (!input) return std::unexpected("installed Mandatory build-gate receipt is unreadable");
        const nlohmann::json receipt = nlohmann::json::parse(input);
        return ValidateBuildGateReceipt(receipt, alignment, executable, resource_root);
    } catch (const std::exception& error) {
        return std::unexpected("installed Mandatory build-gate receipt is invalid: " +
                               std::string(error.what()));
    }
}

}  // namespace sirius::app
