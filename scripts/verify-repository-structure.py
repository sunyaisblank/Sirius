#!/usr/bin/env python3
"""Fail when source topology or immutable build-input authority regresses."""

from __future__ import annotations

import re
import sys
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE_ROOT = ROOT / "src" / "sirius"

LAYER_DEPENDENCIES = {
    "base": {"base"},
    "core": {"base", "core"},
    "oracle": {"base", "oracle"},
    "backend": {"base", "core", "backend"},
    "render": {"base", "core", "backend", "render"},
    "app": {"base", "core", "backend", "render", "app"},
}

EXPECTED_VENDOR_DIRECTORIES = {"glad", "glfw", "stb", "tinyexr"}
RETIRED_PATHS = (
    ROOT / "lib" / "glm",
    ROOT / "lib" / "imgui",
    ROOT / "lib" / "glfw" / "examples",
    ROOT / "lib" / "glfw" / "tests",
    ROOT / "lib" / "glfw" / "deps" / "glad",
    ROOT / "lib" / "glfw" / "deps" / "getopt.c",
    ROOT / "lib" / "glfw" / "deps" / "getopt.h",
    ROOT / "lib" / "glfw" / "deps" / "linmath.h",
    ROOT / "lib" / "glfw" / "deps" / "nuklear.h",
    ROOT / "lib" / "glfw" / "deps" / "nuklear_glfw_gl2.h",
    ROOT / "lib" / "glfw" / "deps" / "stb_image_write.h",
    ROOT / "lib" / "glfw" / "deps" / "tinycthread.c",
    ROOT / "lib" / "glfw" / "deps" / "tinycthread.h",
    SOURCE_ROOT / "render" / "render_config.h",
)

CPP_TOKEN = re.compile(r'(?<![A-Za-z0-9_${])([A-Za-z0-9_./-]+\.cpp)')
PROJECT_INCLUDE = re.compile(r'^\s*#\s*include\s+"sirius/([^/]+)/[^\"]+"', re.MULTILINE)
CAMEL_IDENTIFIER = re.compile(r"\b[a-z][a-z0-9_]*[A-Z][A-Za-z0-9_]*\b")
FULL_REVISION = re.compile(r"[0-9a-f]{40}")
REMOTE_ACTION = re.compile(r"^\s*-\s+uses:\s+([^@\s]+)@([^\s#]+)", re.MULTILINE)
FETCHCONTENT_REVISION = re.compile(r"\bGIT_TAG\s+([^\s)#]+)")
CPP_NON_CODE = re.compile(
    r'"(?:\\.|[^"\\])*"|\'(?:\\.|[^\'\\])*\'|//[^\n]*|/\*.*?\*/',
    re.MULTILINE | re.DOTALL,
)
BOUNDARY_HEADERS = (
    SOURCE_ROOT / "app" / "config" / "config_schema.h",
    SOURCE_ROOT / "app" / "config" / "session_config_adapter.h",
    SOURCE_ROOT / "app" / "viewer" / "interactive_viewer.h",
    SOURCE_ROOT / "backend" / "cpu" / "geodesic_tracer.h",
    SOURCE_ROOT / "render" / "film_config.h",
    SOURCE_ROOT / "render" / "exr_writer.h",
    SOURCE_ROOT / "render" / "session" / "render_session.h",
)
STRICT_TEST_VOLUME_TARGETS = {
    ROOT / "tests" / "app" / "CMakeLists.txt": "sirius_app_tests",
    ROOT / "tests" / "render" / "CMakeLists.txt": "sirius_render_tests",
}
FULL_QUALIFICATION_JOBS = (
    "linux-gate",
    "linux-sanitizers",
    "windows-build",
    "macos-build",
)
NON_RENDER_INTEGRATION_JOBS = (
    "integration-no-render",
    "integration-windows-no-render",
    "integration-macos-no-render",
)
INTEGRATION_CONTROLS = (
    "OperationalEvidence.SourceAndIdealGovernanceRejectUncoveredClaims",
    "OperationalBuildPolicy.ReleaseCannotWeakenGates",
    "OperationalBuildGate.ReleaseInstallRequiresExactPassedEstate",
    "OperationalAttestation.FalseExternalEvidenceIsRejected",
    "OperationalAttestation.NativeRuntimeProducerRejectsWrongHostAndDevice",
    "OperationalAlignment.IncompleteOrAmbiguousAttestationSetIsRejected",
    "OperationalAlignment.InstalledReceiptTamperingBlocksReadiness",
    "AlignmentAuthority.CompiledReceiptMatchesTheStagedRuntimeAuthority",
    "BuildGateAuthority.ReleaseReceiptBindsEveryInstalledProductAtInitialisation",
)
INTEGRATION_TARGETS = (
    "SiriusAlignmentGate",
    "SiriusSourceGovernance",
    "sirius",
    "sirius_base_tests",
    "sirius_core_tests",
    "sirius_oracle_tests",
    "sirius_backend_tests",
    "sirius_app_tests",
    "sirius_render_tests",
)


def relative(path: Path) -> str:
    return path.relative_to(ROOT).as_posix()


def cmake_source_owners() -> Counter[Path]:
    owners: Counter[Path] = Counter()
    for cmake_file in SOURCE_ROOT.rglob("CMakeLists.txt"):
        cmake = re.sub(r"#[^\n]*", "", cmake_file.read_text(encoding="utf-8"))
        cmake = re.sub(r"set_source_files_properties\s*\([^)]*\)", "", cmake)
        for token in CPP_TOKEN.findall(cmake):
            candidate = (cmake_file.parent / token).resolve()
            if candidate.is_file() and candidate.is_relative_to(SOURCE_ROOT):
                owners[candidate] += 1
    return owners


def immutable_input_errors(workflow: str, dependencies: str) -> list[str]:
    errors: list[str] = []
    for action, revision in REMOTE_ACTION.findall(workflow):
        if action.startswith("./"):
            continue
        if FULL_REVISION.fullmatch(revision) is None:
            errors.append(f"remote action {action} is not pinned to a full commit")

    for revision in FETCHCONTENT_REVISION.findall(dependencies):
        if FULL_REVISION.fullmatch(revision) is None:
            errors.append(f"FetchContent input {revision} is not pinned to a full commit")

    if "releases/latest" in workflow or "/download/latest/" in workflow:
        errors.append("CI downloads an asset through a mutable latest-release route")
    if re.search(r"--output\s+['\"]?slang\.tgz", workflow):
        errors.append("CI downloads Slang into the qualification source tree")
    if re.search(r"install_(?:swiftshader|lavapipe):\s*true", workflow):
        errors.append("CI delegates a qualification rasterizer to mutable latest selection")

    shell_downloads = workflow.count("curl --fail --location")
    shell_hash_checks = workflow.count("sha256sum --check -")
    if shell_hash_checks < shell_downloads:
        errors.append("a CI shell download has no checked-in SHA-256 verification")
    powershell_downloads = workflow.count("Invoke-WebRequest")
    powershell_hash_checks = workflow.count("Get-FileHash -Algorithm SHA256")
    if powershell_hash_checks < powershell_downloads:
        errors.append("a CI PowerShell download has no checked-in SHA-256 verification")
    return errors


def workflow_job(workflow: str, name: str) -> str | None:
    match = re.search(
        rf"^  {re.escape(name)}:\s*\n(.*?)(?=^  [A-Za-z0-9_-]+:\s*\n|\Z)",
        workflow,
        re.MULTILINE | re.DOTALL,
    )
    return None if match is None else match.group(1)


def integration_boundary_errors(workflow: str) -> list[str]:
    """Keep cheap integration proof distinct from promotable qualification."""
    errors: list[str] = []
    if "workflow_dispatch:" not in workflow:
        errors.append("CI has no explicit non-render integration dispatch")
    integration_condition = (
        "if: github.event_name == 'pull_request' || "
        "github.event_name == 'workflow_dispatch'"
    )
    for name in NON_RENDER_INTEGRATION_JOBS:
        integration = workflow_job(workflow, name)
        if integration is None:
            errors.append(f"CI has no non-render integration job: {name}")
            continue
        if integration_condition not in integration:
            errors.append(f"non-render integration has an unsafe event boundary: {name}")
        if "-DSIRIUS_ALIGNMENT_MODE=qualification" not in integration:
            errors.append(f"non-render integration does not use strict topology: {name}")
        if "RunMandatoryTests" in integration or "ctest --preset" in integration:
            errors.append(f"non-render integration can execute the full estate: {name}")
        if integration.count("cmake --build") != 1 or "--target" not in integration:
            errors.append(f"non-render integration can escape compile targets: {name}")
        if integration.count("ctest ") != 1 or "--no-tests=error" not in integration:
            errors.append(f"non-render integration lacks one fail-closed CTest: {name}")
        if "-R '^(" not in integration:
            errors.append(f"non-render integration does not bound CTest controls: {name}")
        if integration.count("mandatory_gate.json") != 2:
            errors.append(f"non-render integration does not prove non-promotion: {name}")
        if (
            "verify-attestation.py --record" in integration
            or "upload-artifact@" in integration
        ):
            errors.append(f"non-render integration can publish evidence: {name}")
        for target in INTEGRATION_TARGETS:
            if (
                re.search(
                    rf"(?<![A-Za-z0-9_]){re.escape(target)}(?![A-Za-z0-9_])",
                    integration,
                )
                is None
            ):
                errors.append(f"non-render integration {name} does not compile {target}")
        for control in INTEGRATION_CONTROLS:
            if control.replace(".", r"\.") not in integration:
                errors.append(f"non-render integration {name} omits control {control}")

    windows_integration = workflow_job(workflow, "integration-windows-no-render")
    if windows_integration is not None:
        if "install_runtime: true" not in windows_integration:
            errors.append("Windows non-render integration cannot load linked Vulkan binaries")
        for driver_marker in (
            "install_swiftshader: true",
            "install_lavapipe: true",
            "VK_DRIVER_FILES=",
            "vulkaninfo",
        ):
            if driver_marker in windows_integration:
                errors.append(
                    "Windows non-render integration can select or exercise a Vulkan driver: "
                    + driver_marker
                )

    for name in FULL_QUALIFICATION_JOBS:
        job = workflow_job(workflow, name)
        if job is None:
            errors.append(f"CI full qualification job is missing: {name}")
        elif "if: github.event_name == 'push'" not in job:
            errors.append(
                f"CI full qualification job can execute outside protected pushes: {name}"
            )
    return errors


def verify_integration_boundary_policy() -> None:
    targets = " ".join(INTEGRATION_TARGETS)
    controls = "|".join(name.replace(".", r"\.") for name in INTEGRATION_CONTROLS)
    integration_condition = (
        "    if: github.event_name == 'pull_request' || "
        "github.event_name == 'workflow_dispatch'\n"
    )
    integration_body = (
        integration_condition
        + "    run: cmake -DSIRIUS_ALIGNMENT_MODE=qualification\n"
        + "    run: cmake --build --target "
        + targets
        + "\n"
        + "    run: ctest --test-dir bin/integration --no-tests=error -R '^("
        + controls
        + ")$'\n"
        + "    run: test ! -e mandatory_gate.json\n"
        + "    run: test ! -e mandatory_gate.json\n"
    )
    valid = (
        "workflow_dispatch:\n"
        "jobs:\n"
        + "".join(
            f"  {name}:\n"
            + (
                "    install_runtime: true\n"
                if name == "integration-windows-no-render"
                else ""
            )
            + integration_body
            for name in NON_RENDER_INTEGRATION_JOBS
        )
        + "".join(
            f"  {name}:\n    if: github.event_name == 'push'\n"
            for name in FULL_QUALIFICATION_JOBS
        )
    )
    if integration_boundary_errors(valid):
        raise RuntimeError("integration-boundary policy rejected the strict split")
    weakened = valid.replace(
        "github.event_name == 'pull_request'", "github.event_name == 'push'", 1
    ).replace("ctest --test-dir", "ctest --preset", 1)
    if not integration_boundary_errors(weakened):
        raise RuntimeError("integration-boundary policy accepted full-estate PR execution")
    missing_windows_loader = valid.replace("    install_runtime: true\n", "", 1)
    if not integration_boundary_errors(missing_windows_loader):
        raise RuntimeError("integration-boundary policy accepted an unloadable Windows binary")
    driver_enabled = valid.replace(
        "    install_runtime: true\n",
        "    install_runtime: true\n    install_swiftshader: true\n",
        1,
    )
    if not integration_boundary_errors(driver_enabled):
        raise RuntimeError("integration-boundary policy accepted Windows driver execution")


def verify_immutable_input_policy() -> None:
    revision = "a" * 40
    if immutable_input_errors(
        f"- uses: actions/example@{revision}\n", f"GIT_TAG {revision}\n"
    ):
        raise RuntimeError("immutable build-input policy rejected full revisions")
    weakened = immutable_input_errors(
        "- uses: actions/example@v1\n"
        "curl --fail --location https://example.invalid/download/latest/tool\n"
        "Invoke-WebRequest https://example.invalid/tool\n"
        "install_swiftshader: true\n",
        "GIT_TAG v1.2.3\n",
    )
    if len(weakened) != 6:
        raise RuntimeError("immutable build-input policy accepted a mutable input")
    dirty_checkout = immutable_input_errors(
        "curl --fail --location https://example.invalid/slang --output slang.tgz\n"
        "sha256sum --check -\n",
        f"GIT_TAG {revision}\n",
    )
    if dirty_checkout != ["CI downloads Slang into the qualification source tree"]:
        raise RuntimeError("immutable build-input policy accepted a dirty qualification tree")


def strict_test_volume_errors(documents: dict[str, tuple[str, str]]) -> list[str]:
    errors: list[str] = []
    for name, (target, document) in documents.items():
        colocated = re.search(
            rf"set_target_properties\s*\(\s*{re.escape(target)}\s+PROPERTIES\s+"
            r'RUNTIME_OUTPUT_DIRECTORY\s+"\$<TARGET_FILE_DIR:sirius>"\s*\)',
            document,
            re.DOTALL,
        )
        if colocated is None:
            errors.append(
                f"{name} must place {target} in the exact sirius executable volume"
            )
        if "SIRIUS_RESOURCE_DIR" in document:
            errors.append(
                f"{name} must not substitute an environment-selected resource volume"
            )
    return errors


def verify_strict_test_volume_policy() -> None:
    runtime_directory = (
        'set_target_properties({target} PROPERTIES\n'
        '    RUNTIME_OUTPUT_DIRECTORY "$<TARGET_FILE_DIR:sirius>")\n'
    )
    valid = {
        "app": (
            "sirius_app_tests",
            runtime_directory.format(target="sirius_app_tests"),
        ),
        "render": (
            "sirius_render_tests",
            runtime_directory.format(target="sirius_render_tests"),
        ),
    }
    if strict_test_volume_errors(valid):
        raise RuntimeError("strict test-volume policy rejected exact product colocation")

    weakened = dict(valid)
    weakened["app"] = ("sirius_app_tests", "add_executable(sirius_app_tests)\n")
    weakened["render"] = (
        "sirius_render_tests",
        valid["render"][1] + 'ENVIRONMENT "SIRIUS_RESOURCE_DIR=/tmp/forged"\n',
    )
    if len(strict_test_volume_errors(weakened)) != 2:
        raise RuntimeError("strict test-volume policy accepted a substituted volume")


def verify() -> list[str]:
    errors: list[str] = []

    verify_immutable_input_policy()
    verify_strict_test_volume_policy()
    verify_integration_boundary_policy()
    workflow = (ROOT / ".github" / "workflows" / "ci.yml").read_text(encoding="utf-8")
    dependencies = (ROOT / "cmake" / "sirius_dependencies.cmake").read_text(
        encoding="utf-8"
    )
    errors.extend(immutable_input_errors(workflow, dependencies))
    errors.extend(integration_boundary_errors(workflow))
    strict_test_documents = {
        relative(path): (target, path.read_text(encoding="utf-8"))
        for path, target in STRICT_TEST_VOLUME_TARGETS.items()
    }
    errors.extend(strict_test_volume_errors(strict_test_documents))

    for path in RETIRED_PATHS:
        if path.exists():
            errors.append(f"retired repository content returned: {relative(path)}")

    vendor_root = ROOT / "lib"
    actual_vendors = {path.name for path in vendor_root.iterdir() if path.is_dir()}
    unexpected_vendors = sorted(actual_vendors - EXPECTED_VENDOR_DIRECTORIES)
    missing_vendors = sorted(EXPECTED_VENDOR_DIRECTORIES - actual_vendors)
    if unexpected_vendors:
        errors.append("unapproved vendored directories: " + ", ".join(unexpected_vendors))
    if missing_vendors:
        errors.append("required vendored directories are missing: " + ", ".join(missing_vendors))

    owners = cmake_source_owners()
    for source in SOURCE_ROOT.rglob("*.cpp"):
        count = owners[source.resolve()]
        if count != 1:
            errors.append(
                f"{relative(source)} must have exactly one CMake owner; found {count}"
            )

    backend_cmake = (SOURCE_ROOT / "backend" / "CMakeLists.txt").read_text(encoding="utf-8")
    render_cmake = (SOURCE_ROOT / "render" / "CMakeLists.txt").read_text(encoding="utf-8")
    if "add_library(sirius_backend_cpu STATIC" not in backend_cmake:
        errors.append("the CPU tracer is not owned by the sirius_backend_cpu target")
    if "cpu/geodesic_tracer.cpp" not in backend_cmake:
        errors.append("sirius_backend_cpu does not compile cpu/geodesic_tracer.cpp")
    if "geodesic_tracer.cpp" in render_cmake:
        errors.append("the render target directly compiles a backend implementation")
    if "sirius_backend_cpu" not in render_cmake:
        errors.append("the render target has no explicit CPU-backend dependency")

    for layer, allowed in LAYER_DEPENDENCIES.items():
        for source in (SOURCE_ROOT / layer).rglob("*"):
            if source.suffix not in {".h", ".cpp"}:
                continue
            text = source.read_text(encoding="utf-8")
            for dependency in PROJECT_INCLUDE.findall(text):
                if dependency not in allowed:
                    errors.append(
                        f"{relative(source)} reaches upward or sideways into sirius::{dependency}"
                    )

    for header in BOUNDARY_HEADERS:
        text = CPP_NON_CODE.sub(" ", header.read_text(encoding="utf-8"))
        for identifier in CAMEL_IDENTIFIER.findall(text):
            if identifier.startswith("k") and len(identifier) > 1 and identifier[1].isupper():
                continue
            errors.append(
                f"{relative(header)} exposes non-snake-case identifier '{identifier}'"
            )

    return errors


def main() -> int:
    errors = verify()
    if errors:
        for error in errors:
            print(f"error: {error}", file=sys.stderr)
        return 1
    source_count = sum(1 for _ in SOURCE_ROOT.rglob("*.cpp"))
    print(
        f"repository structure owns {source_count} first-party translation units exactly once, "
        f"enforces {len(LAYER_DEPENDENCIES)} dependency layers, and permits "
        f"{len(EXPECTED_VENDOR_DIRECTORIES)} live vendor directories with immutable build inputs "
        "and strict test-volume identity"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
