#!/usr/bin/env python3
"""Fail when a required operating-model dimension lacks build-gating evidence."""

import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
MODEL = ROOT / "tests" / "operating_model.json"
LABELS = ROOT / "tests" / "labels" / "CTestLabels.cmake"
REQUIRED_DIMENSIONS = {
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
}
REQUIRED_CAPABILITIES = {
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
}
CAPABILITY_STATES = {
    "supported",
    "bounded",
    "fail_closed",
    "substituted",
    "attestation_required",
}


def source_tests():
    pattern = re.compile(r"TEST(?:_F)?\(\s*(\w+)\s*,\s*(\w+)\s*\)")
    unsupported = re.compile(
        r"^\s*(TEST_P|TYPED_TEST|TYPED_TEST_P)\s*\(",
        re.MULTILINE,
    )
    found = set()
    for path in (ROOT / "tests").rglob("*_test.cpp"):
        content = path.read_text(encoding="utf-8")
        macro = unsupported.search(content)
        if macro:
            raise ValueError(
                f"{path.relative_to(ROOT)} uses ungoverned {macro.group(1)}"
            )
        for suite, name in pattern.findall(content):
            if not name.startswith("DISABLED_"):
                found.add(f"{suite}.{name}")
    return found


def mandatory_tests():
    content = LABELS.read_text(encoding="utf-8")
    found = set()
    block_pattern = re.compile(
        r"set_tests_properties\(\n(?P<names>.*?)"
        r'PROPERTIES LABELS "(?P<labels>[^"]+)"\n\)',
        re.DOTALL,
    )
    for match in block_pattern.finditer(content):
        if "Mandatory" not in match.group("labels").split(";"):
            continue
        for name in re.findall(r"^\s+(\w+\.\w+)\s*$", match.group("names"), re.MULTILINE):
            found.add(name)
    return found


def manual_ctest_labels():
    content = (ROOT / "tests" / "CMakeLists.txt").read_text(encoding="utf-8")
    labels = {}
    pattern = re.compile(r"set_tests_properties\((?P<body>.*?)\)", re.DOTALL)
    for match in pattern.finditer(content):
        body = match.group("body")
        before_properties, separator, properties = body.partition("PROPERTIES")
        if not separator:
            continue
        label_match = re.search(r'LABELS\s+"([^"]+)"', properties)
        if not label_match:
            continue
        for name in re.findall(r"^\s*([A-Za-z0-9_.-]+)\s*$", before_properties, re.MULTILINE):
            labels[name] = set(label_match.group(1).split(";"))
    return labels


def main():
    model = json.loads(MODEL.read_text(encoding="utf-8"))
    if model.get("schema_version") != 2:
        print("operating model: unsupported schema_version", file=sys.stderr)
        return 1

    try:
        tests = source_tests()
    except ValueError as error:
        print(f"operating model: {error}", file=sys.stderr)
        return 1
    mandatory = mandatory_tests()
    ctest_labels = manual_ctest_labels()
    failures = []

    def check_evidence(owner, evidence):
        if not evidence:
            failures.append(f"{owner}: no evidence")
        for item in evidence:
            kind = item.get("kind")
            name = item.get("name", "")
            if kind == "gtest":
                if name not in tests:
                    failures.append(f"{owner}: missing test {name}")
                elif name not in mandatory:
                    failures.append(f"{owner}: test is not Mandatory: {name}")
            elif kind == "ctest":
                if name not in ctest_labels:
                    failures.append(f"{owner}: missing CTest {name}")
                elif not {"Mandatory", "Operational"}.issubset(ctest_labels[name]):
                    failures.append(
                        f"{owner}: CTest is not Mandatory and Operational: {name}"
                    )
            else:
                failures.append(f"{owner}: unknown evidence kind {kind!r}")

    dimensions = model.get("required_dimensions", [])
    ids = [dimension.get("id") for dimension in dimensions]
    if len(ids) != len(set(ids)) or any(not item for item in ids):
        failures.append("dimension ids must be present and unique")
    missing_dimensions = REQUIRED_DIMENSIONS - set(ids)
    unexpected_dimensions = set(ids) - REQUIRED_DIMENSIONS
    if missing_dimensions:
        failures.append(f"missing required dimensions: {sorted(missing_dimensions)}")
    if unexpected_dimensions:
        failures.append(
            f"unexpected dimensions require verifier policy: {sorted(unexpected_dimensions)}"
        )

    for dimension in dimensions:
        evidence = dimension.get("evidence", [])
        check_evidence(dimension.get("id"), evidence)

    capabilities = model.get("capability_contracts", [])
    capability_ids = [capability.get("id") for capability in capabilities]
    if len(capability_ids) != len(set(capability_ids)) or any(
        not item for item in capability_ids
    ):
        failures.append("capability ids must be present and unique")
    missing_capabilities = REQUIRED_CAPABILITIES - set(capability_ids)
    unexpected_capabilities = set(capability_ids) - REQUIRED_CAPABILITIES
    if missing_capabilities:
        failures.append(
            f"missing required capability contracts: {sorted(missing_capabilities)}"
        )
    if unexpected_capabilities:
        failures.append(
            f"unexpected capability contracts require verifier policy: "
            f"{sorted(unexpected_capabilities)}"
        )
    for capability in capabilities:
        owner = f"capability {capability.get('id')}"
        if capability.get("state") not in CAPABILITY_STATES:
            failures.append(f"{owner}: invalid state {capability.get('state')!r}")
        if not isinstance(capability.get("profiles"), list) or not capability["profiles"]:
            failures.append(f"{owner}: profiles must be a non-empty list")
        for field in ("request", "behavior"):
            if not isinstance(capability.get(field), str) or not capability[field].strip():
                failures.append(f"{owner}: {field} must be a non-empty string")
        check_evidence(owner, capability.get("evidence", []))

    if failures:
        for failure in failures:
            print(f"operating model: {failure}", file=sys.stderr)
        return 1
    print(
        f"operating model covers {len(dimensions)} required dimensions and "
        f"{len(capabilities)} capability contracts with Mandatory evidence"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
