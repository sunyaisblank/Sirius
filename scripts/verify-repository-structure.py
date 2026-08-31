#!/usr/bin/env python3
"""Fail when source topology or immutable build-input authority regresses."""

from __future__ import annotations

import re
import subprocess
import sys
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE_ROOT = ROOT / "src" / "sirius"
OPERATING_MODEL = ROOT / "tests" / "operating_model.json"
GIT_ATTRIBUTES = ROOT / ".gitattributes"

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
SRGB_TRANSFER_AUTHORITY = SOURCE_ROOT / "core" / "srgb_transfer.h"
SRGB_TRANSFER_CONSUMERS = {
    SOURCE_ROOT / "core" / "spectral" / "blackbody.h": "colour::EncodeSrgbChannel",
    (
        SOURCE_ROOT / "core" / "spectral" / "spectral_radiance.h"
    ): "colour::EncodeClippedSrgbChannel",
    SOURCE_ROOT / "render" / "image_buffer.h": "core::colour::TryEncodeSrgb8",
    SOURCE_ROOT / "render" / "png_writer.h": "core::colour::TryEncodeSrgb8",
    SOURCE_ROOT / "render" / "exr_writer.h": "core::colour::TryEncodeSrgb8",
}
SRGB_LINEAR_BREAKPOINT = re.compile(
    r"(?<![A-Za-z0-9_.])0\.0031308(?:[fFlL])?(?![A-Za-z0-9_.])"
)
XYZ_SRGB_AUTHORITY = SOURCE_ROOT / "core" / "xyz_srgb.h"
XYZ_SRGB_HOST_CONSUMERS = (
    SOURCE_ROOT / "core" / "spectral" / "blackbody.h",
    SOURCE_ROOT / "core" / "spectral" / "spectral_radiance.h",
)
XYZ_SRGB_DEVICE_AUTHORITY = SOURCE_ROOT / "kernels" / "gr_disk.slang"
XYZ_SRGB_PARITY_PROBE = SOURCE_ROOT / "kernels" / "parity_probe.slang"
XYZ_SRGB_COEFFICIENT_PAIRS = (
    ("12831", "3959"),
    ("329", "214"),
    ("1974", "3959"),
    ("851781", "878810"),
    ("1648619", "878810"),
    ("36519", "878810"),
    ("705", "12673"),
    ("2585", "12673"),
    ("705", "667"),
)
XYZ_SRGB_LEGACY_COEFFICIENT = re.compile(
    r"(?<![A-Za-z0-9_.])(?:3\.2404542|1\.5371385|0\.4985314|0\.9692660|"
    r"1\.8760108|0\.0415560|0\.0556434|0\.2040259|1\.0572252)"
    r"(?:[fFlL])?(?![A-Za-z0-9_.])"
)
CIE1931_OBSERVER_AUTHORITY = SOURCE_ROOT / "core" / "cie1931_observer.h"
CIE1931_OBSERVER_HOST_CONSUMERS = (
    SOURCE_ROOT / "core" / "spectral" / "blackbody.h",
    SOURCE_ROOT / "core" / "spectral" / "spectral_radiance.h",
)
CIE1931_OBSERVER_DEVICE_AUTHORITY = SOURCE_ROOT / "kernels" / "gr_disk.slang"
CIE1931_OBSERVER_PARITY_PROBE = SOURCE_ROOT / "kernels" / "parity_probe.slang"
CIE1931_GAUSSIAN_LOBES = (
    ("442.0", "0.0624", "0.0374", "0.362"),
    ("599.8", "0.0264", "0.0323", "1.056"),
    ("501.1", "0.0490", "0.0382", "0.065"),
    ("568.8", "0.0213", "0.0247", "0.821"),
    ("530.9", "0.0613", "0.0322", "0.286"),
    ("437.0", "0.0845", "0.0278", "1.217"),
    ("459.0", "0.0385", "0.0725", "0.681"),
)
CIE1931_LEGACY_FACADE = re.compile(r"\bCie[XYZ]\s*\(")
CIE1931_STALE_TABLE = re.compile(
    r"(?<![A-Za-z0-9_.])0\.014310(?:[fFlL])?.{0,160}"
    r"0\.043510(?:[fFlL])?.{0,160}0\.134380(?:[fFlL])?",
    re.DOTALL,
)
BLACKBODY_LAWS_AUTHORITY = SOURCE_ROOT / "core" / "spectral" / "blackbody_laws.h"
BLACKBODY_LAWS_CONSUMERS = (
    SOURCE_ROOT / "core" / "spectral" / "blackbody.h",
    SOURCE_ROOT / "core" / "spectral" / "spectral_radiance.h",
)
PHYSICAL_CONSTANTS_AUTHORITY = SOURCE_ROOT / "core" / "constants.h"
PLANCK_EXPM1 = re.compile(r"\bstd::expm1\s*\(")
PLANCK_CONSTANT_IDENTIFIER = re.compile(r"\bk(?:Planck(?:C[12])?|Boltzmann)\b")
WIEN_CONSTANT_LITERAL = re.compile(
    r"(?<![A-Za-z0-9_.])2\.897771955e-3(?:[fFlL])?(?![A-Za-z0-9_.])",
    re.IGNORECASE,
)
ACES_FIT_AUTHORITY = SOURCE_ROOT / "core" / "postprocess.h"
ACES_FIT_CONSUMERS = {
    SOURCE_ROOT / "app" / "config" / "config_schema.h": "kDefaultTonemapperName",
    SOURCE_ROOT / "app" / "config" / "config_loader.cpp": "ParseTonemapType",
    SOURCE_ROOT / "app" / "config" / "session_config_adapter.cpp": "ParseTonemapType",
    SOURCE_ROOT / "render" / "session" / "render_session.h": "TonemapType::AcesFit",
    SOURCE_ROOT / "render" / "session" / "render_session.cpp": "TonemapType::AcesFit",
}
ACES_FIT_CLI = SOURCE_ROOT / "app" / "cli" / "render_command.cpp"
SPECTRAL_RADIANCE_FACADE = SOURCE_ROOT / "core" / "spectral" / "spectral_radiance.h"
ACES_FIT_COEFFICIENTS = ("2.51", "0.03", "2.43", "0.59", "0.14")
LEGACY_ACES_TONEMAP_IDENTIFIER = re.compile(
    r"\bTonemapType::Aces\b|\btonemap::Aces\s*\("
)
FALSE_SPECTRAL_ACES_API = re.compile(r"\bstruct\s+Aces\b|\bToAces\s*\(")
DETACHED_XYZ_AP0_MATRIX = re.compile(
    r"1\.0498110175.{0,240}-0\.4959030231.{0,240}"
    r"1\.3733130458.{0,240}0\.9912520182",
    re.DOTALL,
)
BARE_ACES_CONFIG_LITERAL = re.compile(r'"ACES"')
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
    "OperationalAttestation.PreflightAndNativeRuntimeRejectFalseHostDevice",
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


def authority_checkout_errors(
    payload: bytes, attributes: str, attribute_source: str,
) -> list[str]:
    """Require one checkout-byte identity for the compiled model authority."""
    errors: list[str] = []
    parsed: dict[str, str] = {}
    for line in attributes.splitlines():
        parts = line.rsplit(": ", 2)
        if len(parts) == 3:
            parsed[parts[1]] = parts[2]
    if parsed.get("text") != "set" or parsed.get("eol") != "lf":
        errors.append("operating model checkout is not pinned to text eol=lf")
    if re.search(
        r"^/tests/operating_model\.json[ \t]+text[ \t]+eol=lf[ \t]*$",
        attribute_source,
        re.MULTILINE,
    ) is None:
        errors.append("repository does not own the operating-model LF attribute")
    if not payload.endswith(b"\n") or b"\r" in payload:
        errors.append("operating model is not canonical LF-terminated bytes")
    try:
        payload.decode("utf-8")
    except UnicodeDecodeError:
        errors.append("operating model is not canonical UTF-8")
    return errors


def verify_authority_checkout_policy() -> None:
    valid = (
        "tests/operating_model.json: text: set\n"
        "tests/operating_model.json: eol: lf\n"
    )
    source = "/tests/operating_model.json text eol=lf\n"
    if authority_checkout_errors(b"{}\n", valid, source):
        raise RuntimeError("authority checkout policy rejected canonical model bytes")
    if not authority_checkout_errors(b"{}\r\n", valid, source):
        raise RuntimeError("authority checkout policy accepted CRLF model bytes")
    if not authority_checkout_errors(b"{}", valid, source):
        raise RuntimeError("authority checkout policy accepted unterminated model bytes")
    if not authority_checkout_errors(b"\xff\n", valid, source):
        raise RuntimeError("authority checkout policy accepted non-UTF-8 model bytes")
    if not authority_checkout_errors(
        b"{}\n", valid.replace("lf", "unspecified"), source
    ):
        raise RuntimeError("authority checkout policy accepted unspecified line endings")
    if not authority_checkout_errors(b"{}\n", valid, ""):
        raise RuntimeError("authority checkout policy accepted an external-only attribute")


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
        expected_builds = 1 if name == "integration-no-render" else 2
        if integration.count("cmake --build") != expected_builds or "--target" not in integration:
            errors.append(f"non-render integration can escape compile targets: {name}")
        if integration.count("ctest ") != 1 or "--no-tests=error" not in integration:
            errors.append(f"non-render integration lacks one fail-closed CTest: {name}")
        if "-R '^(" not in integration:
            errors.append(f"non-render integration does not bound CTest controls: {name}")
        if integration.count("mandatory_gate.json") != 2:
            errors.append(f"non-render integration does not prove non-promotion: {name}")
        records_build = integration.count("verify-attestation.py --record-build")
        uploads = integration.count("upload-artifact@")
        if name == "integration-no-render":
            if records_build or uploads or "RunNativeBuildEvidence" in integration:
                errors.append("Linux integration can publish an unrepresented build domain")
        else:
            expected_domain = (
                "windows-native-build"
                if name == "integration-windows-no-render"
                else "macos-native-build"
            )
            if records_build != 1 or uploads != 1:
                errors.append(f"native integration lacks one build evidence producer: {name}")
            if integration.count("if: github.event_name == 'workflow_dispatch'") != 3:
                errors.append(
                    f"native build evidence is not dispatch-only inside integration: {name}"
                )
            for marker in (
                "RunNativeBuildEvidence",
                "native_build_gate_junit.xml",
                "--native-build-gate",
                "--native-build-gate-log",
                expected_domain,
            ):
                if marker not in integration:
                    errors.append(f"native integration evidence omits {marker}: {name}")
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
        for loader_marker in (
            "install_runtime: true",
            "cache: false",
            r"runtime\x64",
            "vulkan-1.dll",
            "$env:GITHUB_PATH",
        ):
            if loader_marker not in windows_integration:
                errors.append(
                    "Windows non-render integration cannot load linked Vulkan binaries: "
                    + loader_marker
                )
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
    native_evidence_body = (
        "    if: github.event_name == 'workflow_dispatch'\n"
        "    run: cmake --build --target RunNativeBuildEvidence\n"
        "    if: github.event_name == 'workflow_dispatch'\n"
        "    run: python scripts/verify-attestation.py --record-build "
        "--native-build-gate native_build_gate.json "
        "--native-build-gate-log native_build_gate_ctest.log "
        "--artifact native_build_gate_junit.xml --domain {domain}\n"
        "    if: github.event_name == 'workflow_dispatch'\n"
        "    uses: actions/upload-artifact@0123456789012345678901234567890123456789\n"
    )
    valid = (
        "workflow_dispatch:\n"
        "jobs:\n"
        + "".join(
            f"  {name}:\n"
            + (
                "    install_runtime: true\n"
                "    cache: false\n"
                "    runtime\\x64\n"
                "    vulkan-1.dll\n"
                "    $env:GITHUB_PATH\n"
                if name == "integration-windows-no-render"
                else ""
            )
            + integration_body
            + (
                native_evidence_body.format(
                    domain=(
                        "windows-native-build"
                        if name == "integration-windows-no-render"
                        else "macos-native-build"
                    )
                )
                if name != "integration-no-render"
                else ""
            )
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
    cached_without_loader = valid.replace("    cache: false\n", "    cache: true\n", 1)
    if not integration_boundary_errors(cached_without_loader):
        raise RuntimeError("integration-boundary policy accepted the runtime-blind SDK cache")
    promotable_pull_request = valid.replace(
        "    if: github.event_name == 'workflow_dispatch'\n"
        "    run: cmake --build --target RunNativeBuildEvidence\n",
        "    run: cmake --build --target RunNativeBuildEvidence\n",
        1,
    )
    if not integration_boundary_errors(promotable_pull_request):
        raise RuntimeError("integration-boundary policy accepted PR-produced build evidence")
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


def srgb_transfer_authority_errors(documents: dict[Path, str]) -> list[str]:
    errors: list[str] = []
    authority = documents.get(SRGB_TRANSFER_AUTHORITY)
    if authority is None:
        return ["the host sRGB transfer authority is missing"]

    authority_code = CPP_NON_CODE.sub(" ", authority)
    for marker in (
        "EncodeSrgbChannel",
        "EncodeClippedSrgbChannel",
        "TryEncodeSrgb8",
    ):
        if marker not in authority_code:
            errors.append(f"the host sRGB transfer authority omits {marker}")
    if len(SRGB_LINEAR_BREAKPOINT.findall(authority_code)) != 1:
        errors.append("the host sRGB transfer authority must own one IEC linear breakpoint")

    required_include = '#include "sirius/core/srgb_transfer.h"'
    for path, marker in SRGB_TRANSFER_CONSUMERS.items():
        document = documents.get(path)
        if document is None:
            errors.append(f"sRGB transfer consumer is missing: {relative(path)}")
            continue
        if required_include not in document or marker not in CPP_NON_CODE.sub(" ", document):
            errors.append(
                f"{relative(path)} does not delegate to the host sRGB transfer authority"
            )

    for path, document in documents.items():
        if path == SRGB_TRANSFER_AUTHORITY:
            continue
        if SRGB_LINEAR_BREAKPOINT.search(CPP_NON_CODE.sub(" ", document)):
            errors.append(f"{relative(path)} reimplements the IEC sRGB transfer breakpoint")
    return errors


def verify_srgb_transfer_authority_policy() -> None:
    valid = {
        SRGB_TRANSFER_AUTHORITY: (
            "EncodeSrgbChannel EncodeClippedSrgbChannel TryEncodeSrgb8 0.0031308"
        ),
        **{
            path: '#include "sirius/core/srgb_transfer.h"\n' + marker
            for path, marker in SRGB_TRANSFER_CONSUMERS.items()
        },
    }
    if srgb_transfer_authority_errors(valid):
        raise RuntimeError("sRGB transfer policy rejected the single host authority")
    duplicated = dict(valid)
    consumer = next(iter(SRGB_TRANSFER_CONSUMERS))
    duplicated[consumer] += "\n0.0031308"
    if not srgb_transfer_authority_errors(duplicated):
        raise RuntimeError("sRGB transfer policy accepted a production reimplementation")
    detached = dict(valid)
    detached[consumer] = "independent_transfer();"
    if not srgb_transfer_authority_errors(detached):
        raise RuntimeError("sRGB transfer policy accepted a detached production consumer")


def xyz_srgb_ratio(numerator: str, denominator: str) -> re.Pattern[str]:
    return re.compile(
        rf"(?<![A-Za-z0-9_.]){numerator}(?:\.0)?(?:[fFlL])?\s*/\s*"
        rf"{denominator}(?:\.0)?(?:[fFlL])?(?![A-Za-z0-9_.])"
    )


def xyz_srgb_authority_errors(documents: dict[Path, str]) -> list[str]:
    errors: list[str] = []
    authority = documents.get(XYZ_SRGB_AUTHORITY)
    if authority is None:
        return ["the host XYZ-D65 to linear-sRGB authority is missing"]

    authority_code = CPP_NON_CODE.sub(" ", authority)
    if "XyzD65ToLinearSrgb" not in authority_code:
        errors.append("the host XYZ-D65 to linear-sRGB authority omits its transform")
    for numerator, denominator in XYZ_SRGB_COEFFICIENT_PAIRS:
        if len(xyz_srgb_ratio(numerator, denominator).findall(authority_code)) != 1:
            errors.append(
                "the host XYZ-D65 to linear-sRGB authority must own exact coefficient "
                f"{numerator}/{denominator} once"
            )

    required_include = '#include "sirius/core/xyz_srgb.h"'
    for path in XYZ_SRGB_HOST_CONSUMERS:
        document = documents.get(path)
        if document is None:
            errors.append(f"XYZ-to-sRGB host consumer is missing: {relative(path)}")
            continue
        code = CPP_NON_CODE.sub(" ", document)
        if required_include not in document or "colour::XyzD65ToLinearSrgb" not in code:
            errors.append(
                f"{relative(path)} does not delegate to the host XYZ-to-sRGB authority"
            )

    device = documents.get(XYZ_SRGB_DEVICE_AUTHORITY)
    if device is None:
        errors.append("the Slang XYZ-D65 to linear-sRGB mirror is missing")
    else:
        device_code = CPP_NON_CODE.sub(" ", device)
        if device_code.count("XyzD65ToLinearSrgb") < 2:
            errors.append("the Slang blackbody path does not consume its XYZ-to-sRGB mirror")
        for numerator, denominator in XYZ_SRGB_COEFFICIENT_PAIRS:
            if len(xyz_srgb_ratio(numerator, denominator).findall(device_code)) != 1:
                errors.append(
                    "the Slang XYZ-D65 to linear-sRGB mirror must own exact coefficient "
                    f"{numerator}/{denominator} once"
                )

    parity = documents.get(XYZ_SRGB_PARITY_PROBE)
    if parity is None or "OP_XYZ_D65_TO_LINEAR_SRGB" not in parity or (
        "XyzD65ToLinearSrgb" not in CPP_NON_CODE.sub(" ", parity)
    ):
        errors.append("the direct XYZ-D65 to linear-sRGB device parity route is missing")

    allowed_ratio_paths = {XYZ_SRGB_AUTHORITY, XYZ_SRGB_DEVICE_AUTHORITY}
    for path, document in documents.items():
        code = CPP_NON_CODE.sub(" ", document)
        if XYZ_SRGB_LEGACY_COEFFICIENT.search(code):
            errors.append(f"{relative(path)} retains a rounded legacy XYZ-to-sRGB coefficient")
        if path in allowed_ratio_paths:
            continue
        for numerator, denominator in XYZ_SRGB_COEFFICIENT_PAIRS:
            if xyz_srgb_ratio(numerator, denominator).search(code):
                errors.append(
                    f"{relative(path)} reimplements exact XYZ-to-sRGB coefficient "
                    f"{numerator}/{denominator}"
                )
    return errors


def verify_xyz_srgb_authority_policy() -> None:
    host_coefficients = " ".join(
        f"{numerator}.0L / {denominator}.0L"
        for numerator, denominator in XYZ_SRGB_COEFFICIENT_PAIRS
    )
    device_coefficients = " ".join(
        f"{numerator}.0f / {denominator}.0f"
        for numerator, denominator in XYZ_SRGB_COEFFICIENT_PAIRS
    )
    valid = {
        XYZ_SRGB_AUTHORITY: "XyzD65ToLinearSrgb " + host_coefficients,
        **{
            path: '#include "sirius/core/xyz_srgb.h"\ncolour::XyzD65ToLinearSrgb'
            for path in XYZ_SRGB_HOST_CONSUMERS
        },
        XYZ_SRGB_DEVICE_AUTHORITY: (
            "public float3 XyzD65ToLinearSrgb "
            + device_coefficients
            + " BlackbodyColor XyzD65ToLinearSrgb"
        ),
        XYZ_SRGB_PARITY_PROBE: "OP_XYZ_D65_TO_LINEAR_SRGB XyzD65ToLinearSrgb",
    }
    if xyz_srgb_authority_errors(valid):
        raise RuntimeError("XYZ-to-sRGB policy rejected the governed host/device authority")

    duplicated = dict(valid)
    consumer = XYZ_SRGB_HOST_CONSUMERS[0]
    duplicated[consumer] += "\n12831.0 / 3959.0"
    if not xyz_srgb_authority_errors(duplicated):
        raise RuntimeError("XYZ-to-sRGB policy accepted a production coefficient copy")

    legacy = dict(valid)
    legacy[consumer] += "\n3.2404542f"
    if not xyz_srgb_authority_errors(legacy):
        raise RuntimeError("XYZ-to-sRGB policy accepted a rounded legacy matrix")

    detached = dict(valid)
    detached[consumer] = "independent_xyz_to_rgb();"
    if not xyz_srgb_authority_errors(detached):
        raise RuntimeError("XYZ-to-sRGB policy accepted a detached host consumer")


def cie1931_lobe(values: tuple[str, str, str, str]) -> re.Pattern[str]:
    def literal(value: str) -> str:
        whole, dot, fraction = value.partition(".")
        decimal = rf"{whole}\.{fraction}" if dot else whole
        return rf"(?<![A-Za-z0-9_.]){decimal}(?:[fFlL])?(?![A-Za-z0-9_.])"

    return re.compile(r"\s*,\s*".join(literal(value) for value in values))


def cie1931_observer_authority_errors(documents: dict[Path, str]) -> list[str]:
    errors: list[str] = []
    authority = documents.get(CIE1931_OBSERVER_AUTHORITY)
    if authority is None:
        return ["the host CIE 1931 observer-fit authority is missing"]

    authority_code = CPP_NON_CODE.sub(" ", authority)
    for marker in ("CieXyzMatching", "Cie1931TwoDegreeFit", "std::isfinite"):
        if marker not in authority_code:
            errors.append(f"the host CIE 1931 observer-fit authority omits {marker}")
    for lobe in CIE1931_GAUSSIAN_LOBES:
        if len(cie1931_lobe(lobe).findall(authority_code)) != 1:
            errors.append(
                "the host CIE 1931 observer-fit authority must own Gaussian lobe "
                + "/".join(lobe)
                + " once"
            )

    required_include = '#include "sirius/core/cie1931_observer.h"'
    for path in CIE1931_OBSERVER_HOST_CONSUMERS:
        document = documents.get(path)
        if document is None:
            errors.append(f"CIE 1931 observer-fit host consumer is missing: {relative(path)}")
            continue
        code = CPP_NON_CODE.sub(" ", document)
        if required_include not in document or "colour::Cie1931TwoDegreeFit" not in code:
            errors.append(
                f"{relative(path)} does not delegate to the host CIE 1931 observer-fit authority"
            )

    device = documents.get(CIE1931_OBSERVER_DEVICE_AUTHORITY)
    if device is None:
        errors.append("the Slang CIE 1931 observer-fit mirror is missing")
    else:
        device_code = CPP_NON_CODE.sub(" ", device)
        if device_code.count("Cie1931TwoDegreeFit") < 2:
            errors.append("the Slang blackbody path does not consume its CIE observer-fit mirror")
        for lobe in CIE1931_GAUSSIAN_LOBES:
            if len(cie1931_lobe(lobe).findall(device_code)) != 1:
                errors.append(
                    "the Slang CIE 1931 observer-fit mirror must own Gaussian lobe "
                    + "/".join(lobe)
                    + " once"
                )

    parity = documents.get(CIE1931_OBSERVER_PARITY_PROBE)
    if parity is None or "OP_CIE_1931_TWO_DEGREE_FIT" not in parity or (
        "Cie1931TwoDegreeFit" not in CPP_NON_CODE.sub(" ", parity)
    ):
        errors.append("the direct CIE 1931 observer-fit device parity route is missing")

    allowed_lobe_paths = {
        CIE1931_OBSERVER_AUTHORITY,
        CIE1931_OBSERVER_DEVICE_AUTHORITY,
    }
    for path, document in documents.items():
        code = CPP_NON_CODE.sub(" ", document)
        if CIE1931_LEGACY_FACADE.search(code):
            errors.append(f"{relative(path)} retains a detached CIE colour-matching facade")
        if CIE1931_STALE_TABLE.search(code) or re.search(r"\bk[XYZ]Bar\b", code):
            errors.append(f"{relative(path)} retains the stale wavelength-mislabeled CIE table")
        if path in allowed_lobe_paths:
            continue
        for lobe in CIE1931_GAUSSIAN_LOBES:
            if cie1931_lobe(lobe).search(code):
                errors.append(
                    f"{relative(path)} reimplements CIE observer-fit Gaussian lobe "
                    + "/".join(lobe)
                )
    return errors


def verify_cie1931_observer_authority_policy() -> None:
    coefficients = " ".join(", ".join(lobe) for lobe in CIE1931_GAUSSIAN_LOBES)
    valid = {
        CIE1931_OBSERVER_AUTHORITY: (
            "CieXyzMatching Cie1931TwoDegreeFit std::isfinite " + coefficients
        ),
        **{
            path: '#include "sirius/core/cie1931_observer.h"\ncolour::Cie1931TwoDegreeFit'
            for path in CIE1931_OBSERVER_HOST_CONSUMERS
        },
        CIE1931_OBSERVER_DEVICE_AUTHORITY: (
            "public float3 Cie1931TwoDegreeFit "
            + coefficients
            + " BlackbodyColor Cie1931TwoDegreeFit"
        ),
        CIE1931_OBSERVER_PARITY_PROBE: (
            "OP_CIE_1931_TWO_DEGREE_FIT Cie1931TwoDegreeFit"
        ),
    }
    if cie1931_observer_authority_errors(valid):
        raise RuntimeError("CIE observer policy rejected the governed host/device authority")

    consumer = CIE1931_OBSERVER_HOST_CONSUMERS[0]
    duplicated = dict(valid)
    duplicated[consumer] += "\n442.0, 0.0624, 0.0374, 0.362"
    if not cie1931_observer_authority_errors(duplicated):
        raise RuntimeError("CIE observer policy accepted a production coefficient copy")

    detached = dict(valid)
    detached[consumer] = "independent_cie_observer();"
    if not cie1931_observer_authority_errors(detached):
        raise RuntimeError("CIE observer policy accepted a detached host consumer")

    stale_table = dict(valid)
    stale_table[consumer] += "\n0.014310, 0.043510, 0.134380"
    if not cie1931_observer_authority_errors(stale_table):
        raise RuntimeError("CIE observer policy accepted the wavelength-mislabeled table")

    missing_parity = dict(valid)
    missing_parity[CIE1931_OBSERVER_PARITY_PROBE] = "unrelated_probe();"
    if not cie1931_observer_authority_errors(missing_parity):
        raise RuntimeError("CIE observer policy accepted a missing direct device parity route")


def aces_fit_sequence() -> re.Pattern[str]:
    literals = [
        rf"(?<![A-Za-z0-9_.]){re.escape(value)}(?:[fFlL])?(?![A-Za-z0-9_.])"
        for value in ACES_FIT_COEFFICIENTS
    ]
    return re.compile(r".{0,180}".join(literals), re.DOTALL)


def aces_contract_errors(documents: dict[Path, str]) -> list[str]:
    errors: list[str] = []
    authority = documents.get(ACES_FIT_AUTHORITY)
    if authority is None:
        return ["the host ACES-fit tone-map authority is missing"]

    authority_code = CPP_NON_CODE.sub(" ", authority)
    for marker in (
        "AcesFit",
        "ParseTonemapType",
        "kDefaultTonemapperName",
        "kTonemapNames",
        "SupportedTonemapperNames",
        "std::isfinite",
    ):
        if marker not in authority_code:
            errors.append(f"the host ACES-fit tone-map authority omits {marker}")
    if len(aces_fit_sequence().findall(authority_code)) != 1:
        errors.append("the host ACES-fit authority must own the five fit coefficients once")
    if "knarkowicz.wordpress.com/2016/01/06/aces-filmic-tone-mapping-curve" not in authority:
        errors.append("the ACES-fit authority does not identify the original fit source")

    for path, marker in ACES_FIT_CONSUMERS.items():
        document = documents.get(path)
        if document is None:
            errors.append(f"ACES-fit consumer is missing: {relative(path)}")
            continue
        if marker not in CPP_NON_CODE.sub(" ", document):
            errors.append(f"{relative(path)} does not consume the explicit ACES-fit authority")

    cli = documents.get(ACES_FIT_CLI)
    if cli is None:
        errors.append("the render CLI tone-map help is missing")
    elif "Tonemapper: ACESFit" not in cli or "not an ACES Output Transform" not in cli:
        errors.append("the render CLI does not distinguish ACESFit from an ACES Output Transform")

    bare_name_paths = {
        ACES_FIT_AUTHORITY,
        SOURCE_ROOT / "app" / "config" / "config_schema.h",
        SOURCE_ROOT / "app" / "config" / "config_loader.cpp",
        SOURCE_ROOT / "app" / "config" / "session_config_adapter.cpp",
    }
    coefficient_pattern = aces_fit_sequence()
    for path, document in documents.items():
        code = CPP_NON_CODE.sub(" ", document)
        if LEGACY_ACES_TONEMAP_IDENTIFIER.search(code):
            errors.append(f"{relative(path)} retains the falsely named ACES tonemapper")
        if FALSE_SPECTRAL_ACES_API.search(code):
            errors.append(
                f"{relative(path)} relabels absolute spectral tristimulus values as ACES"
            )
        if DETACHED_XYZ_AP0_MATRIX.search(code):
            errors.append(
                f"{relative(path)} retains an ungoverned XYZ-to-AP0 matrix"
            )
        if path in bare_name_paths and BARE_ACES_CONFIG_LITERAL.search(document):
            errors.append(f"{relative(path)} accepts or defaults the unrepresented bare ACES name")
        if path != ACES_FIT_AUTHORITY and coefficient_pattern.search(code):
            errors.append(f"{relative(path)} reimplements the ACES-fit coefficient set")
    return errors


def verify_aces_contract_policy() -> None:
    coefficients = " ".join(value + "f" for value in ACES_FIT_COEFFICIENTS)
    valid = {
        ACES_FIT_AUTHORITY: (
            "AcesFit ParseTonemapType kDefaultTonemapperName kTonemapNames "
            "SupportedTonemapperNames "
            "std::isfinite "
            + coefficients
            + " // knarkowicz.wordpress.com/2016/01/06/aces-filmic-tone-mapping-curve"
        ),
        **{path: marker for path, marker in ACES_FIT_CONSUMERS.items()},
        ACES_FIT_CLI: "Tonemapper: ACESFit; not an ACES Output Transform",
        SPECTRAL_RADIANCE_FACADE: "struct SpectralRadiance {};",
    }
    if aces_contract_errors(valid):
        raise RuntimeError("ACES policy rejected the explicit represented/absent contracts")

    duplicated = dict(valid)
    consumer = next(iter(ACES_FIT_CONSUMERS))
    duplicated[consumer] += "\n" + coefficients
    if not aces_contract_errors(duplicated):
        raise RuntimeError("ACES-fit policy accepted a production coefficient copy")

    detached = dict(valid)
    detached[consumer] = "independent_tonemap();"
    if not aces_contract_errors(detached):
        raise RuntimeError("ACES-fit policy accepted a detached production consumer")

    legacy = dict(valid)
    legacy[ACES_FIT_AUTHORITY] += "\nTonemapType::Aces"
    if not aces_contract_errors(legacy):
        raise RuntimeError("ACES-fit policy accepted the falsely named internal selector")

    bare_alias = dict(valid)
    bare_alias[ACES_FIT_AUTHORITY] += '\nif (name == "ACES") return AcesFit;'
    if not aces_contract_errors(bare_alias):
        raise RuntimeError("ACES-fit policy accepted bare ACES as a represented config name")

    false_help = dict(valid)
    false_help[ACES_FIT_CLI] = "Tonemapper: ACES, Reinhard"
    if not aces_contract_errors(false_help):
        raise RuntimeError("ACES-fit policy accepted CLI advertising of an ACES transform")

    false_spectral_api = dict(valid)
    false_spectral_api[SPECTRAL_RADIANCE_FACADE] += "\nstruct Aces {}; Aces ToAces() const;"
    if not aces_contract_errors(false_spectral_api):
        raise RuntimeError("ACES policy accepted absolute spectral values labeled as ACES")

    detached_ap0 = dict(valid)
    detached_ap0[SPECTRAL_RADIANCE_FACADE] += (
        "\n1.0498110175; -0.4959030231; 1.3733130458; 0.9912520182;"
    )
    if not aces_contract_errors(detached_ap0):
        raise RuntimeError("ACES policy accepted the detached XYZ-to-AP0 matrix")


def blackbody_laws_authority_errors(documents: dict[Path, str]) -> list[str]:
    errors: list[str] = []
    authority = documents.get(BLACKBODY_LAWS_AUTHORITY)
    if authority is None:
        return ["the host blackbody-laws authority is missing"]

    authority_code = CPP_NON_CODE.sub(" ", authority)
    for marker in (
        "TryPlanckSpectralRadiancePerMetre",
        "TryWienPeakWavelength",
        "TryStefanBoltzmannExitance",
        "kPlanckC1",
        "kPlanckC2",
        "kWienB",
        "kStefanBoltzmann",
    ):
        if marker not in authority_code:
            errors.append(f"the host blackbody-laws authority omits {marker}")
    if len(PLANCK_EXPM1.findall(authority_code)) != 1:
        errors.append("the host blackbody-laws authority must own one Planck expm1")

    required_include = '#include "sirius/core/spectral/blackbody_laws.h"'
    for path in BLACKBODY_LAWS_CONSUMERS:
        document = documents.get(path)
        if document is None:
            errors.append(f"blackbody-law consumer is missing: {relative(path)}")
            continue
        if (
            required_include not in document
            or "TryPlanckSpectralRadiancePerMetre"
            not in CPP_NON_CODE.sub(" ", document)
        ):
            errors.append(
                f"{relative(path)} does not delegate to the host Planck-law authority"
            )

    constants = documents.get(PHYSICAL_CONSTANTS_AUTHORITY)
    if constants is None or len(
        WIEN_CONSTANT_LITERAL.findall(CPP_NON_CODE.sub(" ", constants))
    ) != 1:
        errors.append("the physical-constants authority must own the Wien constant once")

    for path, document in documents.items():
        code = CPP_NON_CODE.sub(" ", document)
        if path != BLACKBODY_LAWS_AUTHORITY and PLANCK_EXPM1.search(code):
            errors.append(f"{relative(path)} reimplements the host Planck law")
        if (
            path not in {BLACKBODY_LAWS_AUTHORITY, PHYSICAL_CONSTANTS_AUTHORITY}
            and PLANCK_CONSTANT_IDENTIFIER.search(code)
        ):
            errors.append(f"{relative(path)} consumes Planck constants outside their authority")
        if path != PHYSICAL_CONSTANTS_AUTHORITY and WIEN_CONSTANT_LITERAL.search(code):
            errors.append(f"{relative(path)} copies the Wien displacement constant")
        if "StefanBoltzmannRadiance" in code:
            errors.append(f"{relative(path)} retains the mislabeled Stefan-Boltzmann quantity")
    return errors


def verify_blackbody_laws_authority_policy() -> None:
    valid = {
        BLACKBODY_LAWS_AUTHORITY: (
            "TryPlanckSpectralRadiancePerMetre std::expm1( kPlanckC1 kPlanckC2 "
            "TryWienPeakWavelength kWienB TryStefanBoltzmannExitance kStefanBoltzmann"
        ),
        **{
            path: (
                '#include "sirius/core/spectral/blackbody_laws.h"\n'
                "TryPlanckSpectralRadiancePerMetre"
            )
            for path in BLACKBODY_LAWS_CONSUMERS
        },
        PHYSICAL_CONSTANTS_AUTHORITY: "kWienB = 2.897771955e-3",
    }
    if blackbody_laws_authority_errors(valid):
        raise RuntimeError("blackbody-laws policy rejected the single host authority")

    consumer = BLACKBODY_LAWS_CONSUMERS[0]
    duplicated = dict(valid)
    duplicated[consumer] += "\nstd::expm1(exponent)"
    if not blackbody_laws_authority_errors(duplicated):
        raise RuntimeError("blackbody-laws policy accepted a production Planck duplicate")

    reconstructed = dict(valid)
    reconstructed[consumer] += "\nkPlanck * kSpeedOfLight / kBoltzmann"
    if not blackbody_laws_authority_errors(reconstructed):
        raise RuntimeError("blackbody-laws policy accepted reconstructed Planck constants")

    detached = dict(valid)
    detached[consumer] = "independent_planck_law();"
    if not blackbody_laws_authority_errors(detached):
        raise RuntimeError("blackbody-laws policy accepted a detached production consumer")

    copied_wien = dict(valid)
    copied_wien[consumer] += "\n2.897771955e-3"
    if not blackbody_laws_authority_errors(copied_wien):
        raise RuntimeError("blackbody-laws policy accepted a copied Wien constant")

    mislabeled = dict(valid)
    mislabeled[consumer] += "\nStefanBoltzmannRadiance(temperature)"
    if not blackbody_laws_authority_errors(mislabeled):
        raise RuntimeError("blackbody-laws policy accepted a mislabeled physical quantity")


def attestation_preflight_errors(
    preflight: str,
    reuse: str,
    native_producer: str,
    viewer_producer: str,
    tests: str,
    documentation: str,
) -> list[str]:
    errors: list[str] = []
    permitted = re.search(
        r"PERMITTED_CANDIDATE_QUERIES\s*=\s*\{(?P<body>.*?)\n\}",
        preflight,
        re.DOTALL,
    )
    if permitted is None:
        errors.append("attestation preflight has no closed candidate-query set")
    else:
        body = permitted.group("body")
        for query in ('("--json", "info", "system")',
                      '("--json", "info", "readiness")'):
            if query not in body:
                errors.append(f"attestation preflight omits permitted query {query}")
        if any(token in body for token in ('"render"', '"view"', '"test"')):
            errors.append("attestation preflight permits an execution command")

    runner = re.search(
        r"def run_candidate_json\(.*?(?=\ndef clean_source_revision\()",
        preflight,
        re.DOTALL,
    )
    if runner is None or runner.group(0).count("subprocess.run(") != 1:
        errors.append("attestation preflight candidate execution has no single authority")
    elif (
        "arguments in PERMITTED_CANDIDATE_QUERIES" not in runner.group(0)
        or "[str(candidate), *arguments]" not in runner.group(0)
    ):
        errors.append("attestation preflight candidate execution can bypass its query set")
    without_runner = preflight if runner is None else preflight.replace(runner.group(0), "")
    if without_runner.count("subprocess.run(") != 2:
        errors.append("attestation preflight contains a non-git subprocess route")

    for marker in (
        '"kind": "sirius-attestation-preflight"',
        '"status": "ready-for-external-execution"',
        '"promotable": False',
        '"external_execution_completed": False',
        '"domains" not in report',
        '"artifacts" not in report',
        "output_path_is_safe",
        'run_candidate_json(Path("unused"), ("render",), {})',
    ):
        if marker not in preflight:
            errors.append(f"attestation preflight omits fail-closed marker {marker}")
    if "from attestation_preflight import self_test as preflight_self_test" not in native_producer:
        errors.append("native producer self-test does not exercise attestation preflight")
    if "preflight_self_test()" not in native_producer:
        errors.append("native producer does not run the imported preflight controls")
    if "subprocess.run(" in reuse:
        errors.append("qualification reuse introduced an execution subprocess")
    for marker in (
        "verifier.verify_document(document, record_path)",
        '"physical-radeon-780m" in domains',
        '"viewer-native-window-input" not in domains',
        'document.get("device") == selected',
        "verify_current_volume(candidate, document, reusable, record_path)",
        '"reused-qualification-transcript.log"',
        '"promotable": False',
    ):
        if marker not in reuse:
            errors.append(f"qualification reuse omits authority marker {marker}")
    if (
        "from reuse_qualification_evidence import self_test "
        "as qualification_reuse_self_test" not in native_producer
        or "qualification_reuse_self_test()" not in native_producer
    ):
        errors.append("native producer self-test does not exercise qualification reuse")
    for marker in (
        "SIRIUS_REUSE_QUALIFICATION_ATTESTATION",
        "scripts/reuse_qualification_evidence.py",
        'cat "$OUT/reused-qualification-transcript.log"',
        'exact full estate reused from verified hardware authority',
    ):
        if marker not in viewer_producer:
            errors.append(f"viewer producer omits qualification-reuse marker {marker}")
    control = "OperationalAttestation.PreflightAndNativeRuntimeRejectFalseHostDevice"
    if tests.count(control) != 2:
        errors.append("attestation preflight control is not one exact configured CTest")
    if "scripts/attestation_preflight.py" not in documentation:
        errors.append("attestation preflight has no operator documentation")
    if "SIRIUS_REUSE_QUALIFICATION_ATTESTATION" not in documentation:
        errors.append("qualification reuse has no operator documentation")
    return errors


def verify() -> list[str]:
    errors: list[str] = []

    verify_immutable_input_policy()
    verify_strict_test_volume_policy()
    verify_integration_boundary_policy()
    verify_authority_checkout_policy()
    verify_srgb_transfer_authority_policy()
    verify_xyz_srgb_authority_policy()
    verify_cie1931_observer_authority_policy()
    verify_aces_contract_policy()
    verify_blackbody_laws_authority_policy()
    try:
        attribute_source = GIT_ATTRIBUTES.read_text(encoding="utf-8")
        attributes = subprocess.run(
            ["git", "check-attr", "text", "eol", "--", relative(OPERATING_MODEL)],
            cwd=ROOT,
            check=True,
            capture_output=True,
            text=True,
        ).stdout
    except (OSError, subprocess.CalledProcessError, UnicodeDecodeError) as error:
        errors.append(f"could not inspect operating-model checkout attributes: {error}")
    else:
        errors.extend(
            authority_checkout_errors(
                OPERATING_MODEL.read_bytes(), attributes, attribute_source
            )
        )
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
    governed_sources = {
        path: path.read_text(encoding="utf-8")
        for path in SOURCE_ROOT.rglob("*")
        if path.suffix in {".h", ".cpp", ".slang"}
    }
    errors.extend(srgb_transfer_authority_errors(governed_sources))
    errors.extend(xyz_srgb_authority_errors(governed_sources))
    errors.extend(cie1931_observer_authority_errors(governed_sources))
    errors.extend(aces_contract_errors(governed_sources))
    errors.extend(blackbody_laws_authority_errors(governed_sources))
    preflight_path = ROOT / "scripts" / "attestation_preflight.py"
    reuse_path = ROOT / "scripts" / "reuse_qualification_evidence.py"
    if not preflight_path.is_file():
        errors.append("attestation preflight producer is missing")
    elif not reuse_path.is_file():
        errors.append("qualification evidence reuse producer is missing")
    else:
        errors.extend(attestation_preflight_errors(
            preflight_path.read_text(encoding="utf-8"),
            reuse_path.read_text(encoding="utf-8"),
            (ROOT / "scripts" / "validate-native-runtime.py").read_text(
                encoding="utf-8"
            ),
            (ROOT / "scripts" / "validate-viewer.sh").read_text(encoding="utf-8"),
            (ROOT / "tests" / "CMakeLists.txt").read_text(encoding="utf-8"),
            (ROOT / "docs" / "ATTESTATION.md").read_text(encoding="utf-8"),
        ))

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
