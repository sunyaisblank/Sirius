// The info command. Ported from CRCL003A.cpp.
//
// OptiX and the CUDA runtime are retired: system and backend reporting now name
// the CPU tracer (always present) and, when the Vulkan backend is compiled in,
// enumerate Vulkan-visible devices through the backend seam. The metric table
// prints from the single-authority registry, unchanged.

#include "sirius/app/cli/info_command.h"

#include "sirius/app/cli/cli_output.h"
#include "sirius/app/config/config_loader.h"
#include "sirius/app/platform_paths.h"
#include "sirius/base/language_capabilities.h"
#include "sirius/base/operating_model_embedded.h"
#include "sirius/base/resource_locator.h"
#include "sirius/core/metrics/registry.h"

#ifdef SIRIUS_HAS_VULKAN_BACKEND
#include "sirius/backend/device.h"
#endif

#include <nlohmann/json.hpp>

#include <algorithm>
#include <array>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <optional>
#include <string_view>

namespace sirius::app {

namespace cli = cli_output;

namespace {

constexpr std::array<std::string_view, 22> kRequiredDimensionIds = {
    "compile_time_contracts",
    "test_evidence_integrity",
    "test_registration_completeness",
    "configuration_and_input_boundary",
    "install_and_volume_initialisation",
    "operator_script_exit_and_artifact_integrity",
    "runtime_resource_integrity",
    "cpu_render_path",
    "session_lifecycle_and_cancellation",
    "vulkan_capability_and_render_path",
    "sampling_semantics",
    "device_identity_and_allocation_domain",
    "physics_oracles_and_near_extremal_boundary",
    "polarised_transport_to_film",
    "page_thorne_and_volumetric_transfer",
    "ray_bundles_and_filtered_point_catalogue",
    "camera_worldline_lens_and_film",
    "interactive_viewer_projection",
    "output_encoding",
    "memory_and_dispatch_governance",
    "kernel_portability",
    "metric_catalogue_and_declines",
};

constexpr std::array<std::string_view, 19> kRequiredCapabilityIds = {
    "polarised_thin_disk_cpu",
    "polarised_volumetric_transfer",
    "polarised_temporal_motion_blur",
    "polarised_vulkan_render",
    "cpu_scalar_motion_blur",
    "vulkan_scalar_motion_blur",
    "disk_emission_outside_schwarzschild_kerr",
    "morris_thorne_one_sheet",
    "morris_thorne_two_sheet",
    "vulkan_volumetric_sample_bound",
    "viewer_input_state_logic",
    "viewer_native_window_input",
    "native_p2900_contracts",
    "native_p2996_reflection",
    "physical_radeon_780m_runtime",
    "wsl2_dozen_runtime",
    "windows_native_build_and_runtime",
    "macos_native_build_and_moltenvk",
    "physical_imax_5616x4096",
};

constexpr std::array<std::string_view, 5> kCapabilityStates = {
    "supported", "bounded", "fail_closed", "substituted", "attestation_required",
};

bool NonEmptyString(const nlohmann::json& object, std::string_view field) {
    const auto iterator = object.find(field);
    return iterator != object.end() && iterator->is_string() &&
           !iterator->get_ref<const std::string&>().empty();
}

bool ValidateEvidence(const nlohmann::json& owner, std::string& reason) {
    const auto evidence = owner.find("evidence");
    if (evidence == owner.end() || !evidence->is_array() || evidence->empty()) {
        reason = "an operating-model entry has no evidence";
        return false;
    }
    for (const auto& item : *evidence) {
        if (!item.is_object() || !NonEmptyString(item, "kind") || !NonEmptyString(item, "name")) {
            reason = "an operating-model evidence record is malformed";
            return false;
        }
        const std::string& kind = item["kind"].get_ref<const std::string&>();
        if (kind != "gtest" && kind != "ctest") {
            reason = "an operating-model evidence record has an unknown kind";
            return false;
        }
    }
    return true;
}

template <std::size_t N>
bool ValidateExactIds(const nlohmann::json& entries,
                      const std::array<std::string_view, N>& required, std::string_view label,
                      std::string& reason) {
    if (!entries.is_array() || entries.size() != required.size()) {
        reason = "operating model does not contain exactly " + std::to_string(required.size()) +
                 " " + std::string(label) + " entries";
        return false;
    }
    for (const auto& entry : entries) {
        if (!entry.is_object() || !NonEmptyString(entry, "id")) {
            reason = "an operating-model " + std::string(label) + " id is malformed";
            return false;
        }
    }
    for (const std::string_view required_id : required) {
        std::size_t count = 0;
        for (const auto& entry : entries) {
            count += entry["id"].get_ref<const std::string&>() == required_id ? 1U : 0U;
        }
        if (count != 1) {
            reason = "operating model is missing or duplicates required " + std::string(label) +
                     " id " + std::string(required_id);
            return false;
        }
    }
    return true;
}

bool LoadAndValidateOperatingModel(const std::filesystem::path& path, nlohmann::json& model,
                                   std::string& reason) {
    std::string payload;
    try {
        std::ifstream input(path, std::ios::binary);
        if (!input) {
            reason = "operating model is unreadable";
            return false;
        }
        payload.assign(std::istreambuf_iterator<char>(input), std::istreambuf_iterator<char>());
        model = nlohmann::json::parse(payload);
        const nlohmann::json authority = nlohmann::json::parse(
            base::kEmbeddedOperatingModel.begin(), base::kEmbeddedOperatingModel.end());
        if (model != authority) {
            reason = "operating model does not match the compiled capability authority";
            return false;
        }
    } catch (const std::exception& error) {
        reason = "operating model JSON is unreadable: " + std::string(error.what());
        return false;
    }

    const auto schema = model.find("schema_version");
    const auto method = model.find("method");
    if (!model.is_object() || schema == model.end() || !schema->is_number_integer() ||
        schema->get<int>() != 2 || method == model.end() || !method->is_string() ||
        method->get_ref<const std::string&>() != "adversarial-claim-ledger" ||
        !model.contains("required_dimensions") || !model.contains("capability_contracts")) {
        reason = "operating model has an unsupported schema or method";
        return false;
    }
    if (!ValidateExactIds(model["required_dimensions"], kRequiredDimensionIds, "dimension",
                          reason) ||
        !ValidateExactIds(model["capability_contracts"], kRequiredCapabilityIds, "capability",
                          reason)) {
        return false;
    }
    for (const auto& dimension : model["required_dimensions"]) {
        if (!ValidateEvidence(dimension, reason)) return false;
    }
    for (const auto& capability : model["capability_contracts"]) {
        if (!NonEmptyString(capability, "state") || !NonEmptyString(capability, "request") ||
            !NonEmptyString(capability, "behavior")) {
            reason = "an operating-model capability contract is incomplete";
            return false;
        }
        const std::string& state = capability["state"].get_ref<const std::string&>();
        if (std::find(kCapabilityStates.begin(), kCapabilityStates.end(), state) ==
            kCapabilityStates.end()) {
            reason = "an operating-model capability contract has an unknown state";
            return false;
        }
        const auto profiles = capability.find("profiles");
        if (profiles == capability.end() || !profiles->is_array() || profiles->empty() ||
            std::any_of(profiles->begin(), profiles->end(), [](const nlohmann::json& profile) {
                return !profile.is_string() || profile.get_ref<const std::string&>().empty();
            })) {
            reason = "an operating-model capability contract has malformed profiles";
            return false;
        }
        if (!ValidateEvidence(capability, reason)) return false;
    }
    return true;
}

}  // namespace

std::string InfoCommand::Usage() const {
    return R"(Usage: sirius info [subcommand]

Display system and configuration information.

Subcommands:
  system        Display system capabilities (platform, backends)
  metrics       List available spacetime metrics
  config        Show current configuration
  backends      List available render backends
  readiness     Verify resources and live backend readiness
  capabilities  Show the installed supported/bounded/declined/attested model

Examples:
  sirius info system
  sirius info metrics
  sirius info config
  sirius info system --json
)";
}

int InfoCommand::Execute(const std::vector<std::string>& args, const GlobalOptions& globals,
                         SiriusConfig& config) {
    std::string subcommand = "system";
    if (!args.empty()) {
        if (args.size() != 1 || args[0].empty() || args[0].front() == '-') {
            cli::Error("info accepts exactly one subcommand and no additional arguments");
            std::cout << Usage() << std::endl;
            return 1;
        }
        subcommand = args[0];
    }

    if (subcommand == "system") {
        return ShowSystem(globals);
    } else if (subcommand == "metrics") {
        return ShowMetrics(globals);
    } else if (subcommand == "config") {
        return ShowConfig(globals, config);
    } else if (subcommand == "backends") {
        return ShowBackends(globals);
    } else if (subcommand == "readiness") {
        return ShowReadiness(globals);
    } else if (subcommand == "capabilities") {
        return ShowCapabilities(globals);
    } else {
        cli::Error("Unknown subcommand: " + subcommand);
        std::cout << Usage() << std::endl;
        return 1;
    }
}

int InfoCommand::ShowSystem(const GlobalOptions& globals) {
    if (globals.json_output) {
        nlohmann::json j;

        j["platform"] = PlatformPaths::PlatformName();
        j["wsl2"] = PlatformPaths::IsWsl2();
        j["language"]["p2900_contracts"]["native"] = base::kNativeP2900Contracts;
        j["language"]["p2900_contracts"]["implementation"] =
            base::kNativeP2900Contracts ? "native" : "checked-macro";
        j["language"]["p2996_reflection"]["native"] = base::kNativeP2996Reflection;
        j["language"]["p2996_reflection"]["implementation"] = "explicit-schema";

        j["backends"]["cpu"] = {{"available", true}};
#ifdef SIRIUS_HAS_INTERACTIVE_VIEWER
        j["interactive_viewer"] = {{"compiled", true}};
#else
        j["interactive_viewer"] = {
            {"compiled", false},
            {"reason", "OpenGL viewer support was unavailable or disabled at build time"},
        };
#endif

#ifdef SIRIUS_HAS_VULKAN_BACKEND
        auto devices = backend::EnumerateVulkanDevices();
        if (devices.has_value()) {
            j["backends"]["vulkan"]["compiled"] = true;
            j["backends"]["vulkan"]["available"] = !devices->empty();
            j["backends"]["vulkan"]["device_count"] = devices->size();
            nlohmann::json device_array = nlohmann::json::array();
            for (std::size_t index = 0; index < devices->size(); ++index) {
                const auto& d = (*devices)[index];
                device_array.push_back({{"name", d.name},
                                        {"driver_name", d.driver_name},
                                        {"driver_info", d.driver_info},
                                        {"index", index},
                                        {"kind", backend::ToString(d.kind)},
                                        {"vendor_id", d.vendor_id},
                                        {"device_id", d.device_id},
                                        {"api_version", d.api_version},
                                        {"driver_id", d.driver_id},
                                        {"device_local_bytes", d.device_local_bytes},
                                        {"render_memory_bytes", d.render_memory_bytes},
                                        {"supports_fp64", d.supports_fp64}});
            }
            j["backends"]["vulkan"]["devices"] = device_array;
            if (!devices->empty()) {
                if (auto selected = backend::ResolveVulkanDeviceIndex(*devices); selected) {
                    j["backends"]["vulkan"]["selected_device_index"] = *selected;
                    j["backends"]["vulkan"]["selected_device_name"] = (*devices)[*selected].name;
                } else {
                    j["backends"]["vulkan"]["selection_error"] = selected.error().Description();
                }
            }
        } else {
            j["backends"]["vulkan"]["compiled"] = true;
            j["backends"]["vulkan"]["available"] = false;
            j["backends"]["vulkan"]["error"] = devices.error().Description();
        }
#else
        j["backends"]["vulkan"]["compiled"] = false;
        j["backends"]["vulkan"]["available"] = false;
#endif

        cli::PrintJson(j.dump(2));
    } else {
        std::vector<cli::TableRow> rows;

        rows.push_back({"Platform", PlatformPaths::PlatformName(), false});
        if (PlatformPaths::IsWsl2()) {
            rows.push_back({"Environment", "WSL2", false});
        }
        rows.push_back({"P2900 contracts",
                        base::kNativeP2900Contracts ? "Native" : "Checked-macro substitute",
                        false});
        rows.push_back({"P2996 reflection",
                        base::kNativeP2996Reflection
                            ? "Compiler advertises native; Sirius remains explicit-schema"
                            : "Unavailable; explicit-schema substitute",
                        false});

        rows.push_back({"Available Backends", "", true});
        rows.push_back({"CPU", "Available (reference path)", false});

#ifdef SIRIUS_HAS_VULKAN_BACKEND
        auto devices = backend::EnumerateVulkanDevices();
        if (devices.has_value()) {
            if (devices->empty()) {
                rows.push_back({"Vulkan", "Compiled in (no devices found)", false});
            } else {
                for (const auto& d : *devices) {
                    rows.push_back(
                        {"Vulkan", d.name + " (" + backend::ToString(d.kind) + ")", false});
                }
            }
        } else {
            rows.push_back({"Vulkan", "Compiled in (enumeration declined)", false});
        }
#else
        rows.push_back({"Vulkan", "Not compiled", false});
#endif

        cli::PrintTable("System Information", rows);
    }

    return 0;
}

int InfoCommand::ShowMetrics(const GlobalOptions& globals) {
    // The identity registry is the single metric catalogue; this command has no
    // list of its own.
    const auto& registry = core::MetricRegistry();

    if (globals.json_output) {
        nlohmann::json j = nlohmann::json::array();
        for (const auto& m : registry) {
            nlohmann::json aliases = nlohmann::json::array();
            for (const char* alias : m.aliases) {
                if (alias != nullptr) aliases.push_back(alias);
            }
            j.push_back({{"name", m.canonical_name},
                         {"aliases", aliases},
                         {"parameters", m.parameters},
                         {"disk", core::ToString(core::DiskSupportFor(m.id))},
                         {"backends", {{"cpu", m.cpu_supported}, {"gpu", m.gpu_supported}}}});
        }
        cli::PrintJson(j.dump(2));
    } else {
        std::vector<cli::TableRow> rows;
        for (const auto& m : registry) {
            std::string backends =
                m.gpu_supported ? (m.cpu_supported ? "CPU, GPU" : "GPU") : "CPU only";
            rows.push_back({m.canonical_name,
                            std::string(m.parameters) + " [" + backends +
                                "; disk: " + core::ToString(core::DiskSupportFor(m.id)) + "]",
                            false});
        }
        cli::PrintTable("Available Metrics", rows);
    }

    return 0;
}

int InfoCommand::ShowConfig(const GlobalOptions& globals, const SiriusConfig& config) {
    if (globals.json_output) {
        nlohmann::json j = config;

        auto loaded_path = ConfigLoader::GetLoadedConfigPath();
        if (loaded_path.has_value()) {
            j["_source"] = loaded_path->string();
        } else {
            j["_source"] = "defaults";
        }

        cli::PrintJson(j.dump(2));
    } else {
        auto loaded_path = ConfigLoader::GetLoadedConfigPath();
        if (loaded_path.has_value()) {
            cli::Info("Configuration loaded from: " + loaded_path->string());
        } else {
            cli::Info("Using default configuration");
        }
        std::cout << std::endl;

        cli::PrintConfig(config);
    }

    return 0;
}

int InfoCommand::ShowBackends(const GlobalOptions& globals) {
    if (globals.json_output) {
        nlohmann::json j;

        j["cpu"] = {{"available", true},
                    {"description", "CPU software ray tracing"},
                    {"features", {"Multi-threaded", "Portable"}}};

#ifdef SIRIUS_HAS_VULKAN_BACKEND
        auto devices = backend::EnumerateVulkanDevices();
        j["vulkan"]["compiled"] = true;
        j["vulkan"]["available"] = devices.has_value() && !devices->empty();
        j["vulkan"]["description"] =
            "Vulkan compute backend (wired to the render session; --backend vulkan)";
        j["vulkan"]["device_count"] = devices.has_value() ? devices->size() : 0;
        if (!devices.has_value()) {
            j["vulkan"]["reason"] = devices.error().Description();
        } else if (devices->empty()) {
            j["vulkan"]["reason"] = "no Vulkan device is visible";
        }
#else
        j["vulkan"] = {{"compiled", false},
                       {"available", false},
                       {"description", "Vulkan compute backend"},
                       {"reason", "Vulkan development files not found at build time"}};
#endif

        cli::PrintJson(j.dump(2));
    } else {
        std::vector<cli::TableRow> rows;

        rows.push_back({"CPU", "Available (multi-threaded reference path)", false});

#ifdef SIRIUS_HAS_VULKAN_BACKEND
        auto devices = backend::EnumerateVulkanDevices();
        if (!devices.has_value()) {
            rows.push_back({"Vulkan",
                            "Compiled in; enumeration declined: " + devices.error().Description(),
                            false});
        } else if (devices->empty()) {
            rows.push_back({"Vulkan", "Compiled in; no device is visible", false});
        } else {
            rows.push_back({"Vulkan",
                            "Compiled in, " + std::to_string(devices->size()) +
                                " device(s) (--backend vulkan)",
                            false});
        }
#else
        rows.push_back({"Vulkan", "Not available (build without Vulkan)", false});
#endif

        cli::PrintTable("Render Backends", rows);
    }

    return 0;
}

int InfoCommand::ShowCapabilities(const GlobalOptions& globals) {
    const auto path = base::ResolveResource("model/operating_model.json");
    if (!path) {
        cli::Error("Installed operating model is missing");
        return 1;
    }
    nlohmann::json model;
    std::string reason;
    if (!LoadAndValidateOperatingModel(*path, model, reason)) {
        cli::Error("Installed " + reason);
        return 1;
    }
    if (globals.json_output) {
        nlohmann::json output{
            {"schema_version", model["schema_version"]},
            {"capability_contracts", model["capability_contracts"]},
        };
        cli::PrintJson(output.dump(2));
        return 0;
    }

    std::vector<cli::TableRow> rows;
    for (const auto& capability : model["capability_contracts"]) {
        const std::string id = capability["id"];
        const std::string state = capability["state"];
        rows.push_back({id, state + ": " + capability["behavior"].get<std::string>(), false});
    }
    cli::PrintTable("Operational Capability Contracts", rows);
    return 0;
}

int InfoCommand::ShowReadiness(const GlobalOptions& globals) {
#ifdef SIRIUS_HAS_VULKAN_BACKEND
    constexpr bool kVulkanCompiled = true;
#else
    constexpr bool kVulkanCompiled = false;
#endif
    std::vector<std::string> required = {
        "assets/Starfield.png",
        "model/operating_model.json",
    };
#ifdef SIRIUS_HAS_INTERACTIVE_VIEWER
    required.push_back("shaders/RDSD003A.vert");
    required.push_back("shaders/RDSD003A.frag");
#endif
#ifdef SIRIUS_HAS_VULKAN_BACKEND
    required.push_back("kernels/trace.spv");
    required.push_back("kernels/trace_fp32comp.spv");
    required.push_back("kernels/trace_fp64.spv");
#endif

    nlohmann::json resources = nlohmann::json::object();
    bool resources_ready = true;
    for (const std::string& relative : required) {
        const auto path = base::ResolveResource(relative);
        bool valid = path.has_value();
        std::string reason;
        if (relative == "model/operating_model.json" && path.has_value()) {
            nlohmann::json model;
            valid = LoadAndValidateOperatingModel(*path, model, reason);
        }
        resources[relative] = {
            {"available", path.has_value()},
            {"valid", valid},
            {"path", path.has_value() ? path->string() : ""},
            {"reason", reason},
        };
        resources_ready = resources_ready && valid;
    }

    bool vulkan_available = false;
    std::size_t vulkan_devices = 0;
    std::optional<std::size_t> vulkan_selected;
    std::string vulkan_reason;
#ifdef SIRIUS_HAS_VULKAN_BACKEND
    if (auto devices = backend::EnumerateVulkanDevices(); devices.has_value()) {
        vulkan_devices = devices->size();
        vulkan_available = !devices->empty();
        if (!vulkan_available) {
            vulkan_reason = "no Vulkan device is visible";
        } else if (auto selected = backend::ResolveVulkanDeviceIndex(*devices); selected) {
            vulkan_selected = *selected;
        } else {
            vulkan_available = false;
            vulkan_reason = selected.error().Description();
        }
    } else {
        vulkan_reason = devices.error().Description();
    }
#else
    vulkan_reason = "Vulkan render backend was not compiled";
#endif

    // CPU is the reference path and defines minimum system readiness. Vulkan is
    // reported separately because its device availability is environmental.
    const bool ready = resources_ready;
    if (globals.json_output) {
        nlohmann::json j;
        j["ready"] = ready;
        j["resources"] = resources;
        j["backends"]["cpu"] = {{"ready", resources_ready}};
#ifdef SIRIUS_HAS_INTERACTIVE_VIEWER
        j["interactive_viewer"] = {
            {"compiled", true},
            {"logic_ready", resources_ready},
            {"ready", false},
            {"state", "attestation_required"},
            {"reason",
             "native GLFW window creation and keyboard/pointer delivery require a "
             "viewer-native-window-input attestation"},
        };
#else
        j["interactive_viewer"] = {
            {"compiled", false},
            {"ready", false},
            {"reason", "OpenGL viewer support was unavailable or disabled at build time"},
        };
#endif
        j["backends"]["vulkan"] = {
            {"compiled", kVulkanCompiled},
            {"ready", resources_ready && vulkan_available},
            {"device_count", vulkan_devices},
            {"selected_device_index", vulkan_selected.has_value() ? nlohmann::json(*vulkan_selected)
                                                                  : nlohmann::json(nullptr)},
            {"reason", vulkan_reason},
        };
        cli::PrintJson(j.dump(2));
    } else {
        std::vector<cli::TableRow> rows;
        rows.push_back({"System", ready ? "READY" : "NOT READY", false});
        rows.push_back({"Runtime resources", resources_ready ? "Complete" : "Incomplete", false});
        rows.push_back(
            {"CPU reference", resources_ready ? "Ready" : "Blocked by resources", false});
#ifdef SIRIUS_HAS_INTERACTIVE_VIEWER
        rows.push_back({"Interactive viewer",
                        resources_ready
                            ? "Headless logic ready; native window/input requires attestation"
                            : "Blocked by resources; native window/input requires attestation",
                        false});
#else
        rows.push_back(
            {"Interactive viewer", "Not compiled (CPU/Vulkan CLI remains complete)", false});
#endif
        rows.push_back({"Vulkan",
                        vulkan_available && resources_ready
                            ? "Ready (" + std::to_string(vulkan_devices) + " device(s))"
                            : "Not ready: " + vulkan_reason,
                        false});
        cli::PrintTable("Operational Readiness", rows);
    }
    return ready ? 0 : 1;
}

}  // namespace sirius::app
