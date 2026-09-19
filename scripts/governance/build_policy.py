"""Immutable inputs, CI execution boundaries and qualification source identity."""

from __future__ import annotations

import re

from .common import ROOT


GIT_ATTRIBUTES = ROOT / ".gitattributes"
FULL_REVISION = re.compile(r"[0-9a-f]{40}")
REMOTE_ACTION = re.compile(r"^\s*-\s+uses:\s+([^@\s]+)@([^\s#]+)", re.MULTILINE)
FETCHCONTENT_REVISION = re.compile(r"\bGIT_TAG\s+([^\s)#]+)")
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


def attestation_source_authority_errors(source: str) -> list[str]:
    """Keep standalone and producer verification bound to one selected checkout."""
    errors: list[str] = []
    required_markers = {
        "selected-root model hashing":
            "def source_operating_model_sha256(source_root=ROOT)",
        "clean Git source inspection":
            "expected_source_revision = inspect_build_source(source_root)",
        "selected-root model forwarding":
            "source_operating_model_sha256(source_root)",
        "stable source reinspection":
            "inspect_build_source(source_root, expected_source_revision)",
        "standalone CLI source forwarding":
            "verify_path(args.attestation, args.source_root or ROOT)",
        "native producer source forwarding":
            "verify_document_against_source(document, args.output, args.source_root)",
    }
    for label, marker in required_markers.items():
        if marker not in source:
            errors.append(f"attestation verifier omits {label}")
    return errors


def verify_attestation_source_authority_policy() -> None:
    source = (ROOT / "scripts" / "verify-attestation.py").read_text(encoding="utf-8")
    if attestation_source_authority_errors(source):
        raise RuntimeError("attestation source-authority policy rejected the live verifier")
    detached_cli = source.replace(
        "verify_path(args.attestation, args.source_root or ROOT)",
        "verify_path(args.attestation)",
    )
    if not attestation_source_authority_errors(detached_cli):
        raise RuntimeError("attestation policy accepted an inert --source-root argument")
    detached_model = source.replace(
        "source_operating_model_sha256(source_root)",
        "source_operating_model_sha256()",
    )
    if not attestation_source_authority_errors(detached_model):
        raise RuntimeError("attestation policy accepted verifier-checkout model substitution")


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


FULL_QUALIFICATION_CONDITION = (
    "github.event_name == 'push' || "
    "(github.event_name == 'workflow_dispatch' && inputs.full_qualification == true)"
)
GATE_DIAGNOSTIC_PATHS = (
    "bin/**/generated/sirius/*gate*.json",
    "bin/**/generated/sirius/*gate*.xml",
    "bin/**/generated/sirius/*gate*.log",
    "bin/**/Testing/Temporary/LastTest.log",
    "bin/**/Testing/Temporary/LastTestsFailed.log",
    "linux-tests.xml",
    "windows-tests.xml",
    "macos-tests.xml",
)
MACOS_RUNTIME_DIAGNOSTIC_PATHS = (
    '${{ runner.temp }}/sirius-macos-runtime/*/native-runtime-tests.xml',
    '${{ runner.temp }}/sirius-macos-runtime/*/native-runtime-transcript.log',
    '${{ runner.temp }}/sirius-macos-runtime/producer.log',
    '${{ runner.temp }}/sirius-macos-runtime/os-version.txt',
    '${{ runner.temp }}/sirius-macos-runtime/metal-displays.json',
    '${{ runner.temp }}/sirius-macos-runtime/provider-route.txt',
)
MACOS_RUNTIME_PRODUCER_STEP = r"""      - name: Exact native MoltenVK runtime evidence
        if: success()
        shell: bash
        env:
          CMAKE_BUILD_PARALLEL_LEVEL: 3
        run: |
          set -euo pipefail
          python3 scripts/validate-native-runtime.py --expected-revision "${{ github.sha }}" --output-root "$RUNNER_TEMP/sirius-macos-runtime" 2>&1 | tee "$RUNNER_TEMP/sirius-macos-runtime/producer.log"
"""
MACOS_RUNTIME_PROVENANCE_STEP = r"""      - name: Record macOS graphics route
        if: success()
        shell: bash
        run: |
          set -euo pipefail
          mkdir -p "$RUNNER_TEMP/sirius-macos-runtime"
          sw_vers | tee "$RUNNER_TEMP/sirius-macos-runtime/os-version.txt"
          system_profiler SPDisplaysDataType -json | tee "$RUNNER_TEMP/sirius-macos-runtime/metal-displays.json"
          printf 'repository=%s\nrevision=%s\nrun_id=%s\nrun_attempt=%s\njob=%s\nrunner=%s\nos=%s\narch=%s\n' "$GITHUB_REPOSITORY" "$GITHUB_SHA" "$GITHUB_RUN_ID" "$GITHUB_RUN_ATTEMPT" "$GITHUB_JOB" "$RUNNER_NAME" "$RUNNER_OS" "$RUNNER_ARCH" | tee "$RUNNER_TEMP/sirius-macos-runtime/provider-route.txt"
"""


def workflow_steps(job: str) -> list[str]:
    """Read the repository's fixed block-step layout without accepting aliases."""
    return re.findall(r"^      - .*?(?=^      - |\Z)", job, re.MULTILINE | re.DOTALL)


def workflow_step_field(step: str, field: str) -> list[str]:
    # The first mapping field can appear beside the sequence marker.
    normalized = re.sub(r"^      - ", "        ", step, count=1)
    return re.findall(rf"^        {re.escape(field)}:\s*([^\n]*)$",
                      normalized, re.MULTILINE)


def gate_diagnostic_upload_valid(step: str, *, macos_runtime: bool = False) -> bool:
    # Keep diagnostics unable to publish an attestation or broaden their file
    # scope. Unknown fields/layouts fail closed rather than becoming exceptions.
    paths = GATE_DIAGNOSTIC_PATHS + (MACOS_RUNTIME_DIAGNOSTIC_PATHS if macos_runtime else ())
    pattern = (
        r"\A      - name: Preserve gate diagnostics\n"
        r"        if: always\(\)\n"
        r"        uses: actions/upload-artifact@[0-9a-f]{40}(?: +#[^\n]*)?\n"
        r"        with:\n"
        r"          name: sirius-gate-diagnostics-\$\{\{ github\.job \}\}\n"
        r"          path: \|\n"
        + "".join(re.escape(f"            {path}\n") for path in paths)
        + r"          if-no-files-found: ignore\n\s*\Z"
    )
    return re.fullmatch(pattern, step) is not None


def macos_runtime_upload_valid(step: str) -> bool:
    pattern = (
        r"\A      - name: Preserve verified MoltenVK runtime evidence\n"
        r"        if: success\(\)\n"
        r"        uses: actions/upload-artifact@[0-9a-f]{40}(?: +#[^\n]*)?\n"
        r"        with:\n"
        r"          name: sirius-macos-moltenvk-attestation\n"
        r"          path: \$\{\{ runner\.temp \}\}/sirius-macos-runtime\n"
        r"          if-no-files-found: error\n\s*\Z"
    )
    return re.fullmatch(pattern, step) is not None


def full_qualification_input_valid(workflow: str) -> bool:
    header = workflow.split("\njobs:", 1)[0]
    dispatches = re.findall(
        r"^  workflow_dispatch:\n(.*?)(?=^  [A-Za-z0-9_-]+:|\Z)",
        header, re.MULTILINE | re.DOTALL,
    )
    if len(dispatches) != 1:
        return False
    inputs = re.findall(
        r"^      full_qualification:\n(.*?)(?=^      [A-Za-z0-9_-]+:|\Z)",
        dispatches[0], re.MULTILINE | re.DOTALL,
    )
    return (len(inputs) == 1
            and re.findall(r"^        type: ([^\n]+)$", inputs[0], re.MULTILINE)
            == ["boolean"]
            and re.findall(r"^        default: ([^\n]+)$", inputs[0], re.MULTILINE)
            == ["false"]
            and re.findall(r"^    inputs:\s*$", dispatches[0], re.MULTILINE)
            == ["    inputs:"])


def integration_boundary_errors(workflow: str) -> list[str]:
    """Keep cheap integration proof distinct from promotable qualification."""
    errors: list[str] = []
    if not full_qualification_input_valid(workflow):
        errors.append("CI full qualification dispatch must use a default-off boolean input")
    integration_condition = (
        "github.event_name == 'pull_request' || "
        "github.event_name == 'workflow_dispatch'"
    )
    for name in NON_RENDER_INTEGRATION_JOBS:
        integration = workflow_job(workflow, name)
        if integration is None:
            errors.append(f"CI has no non-render integration job: {name}")
            continue
        if re.findall(r"^    if: ([^\n]+)$", integration, re.MULTILINE) != [
            integration_condition
        ]:
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
            steps = workflow_steps(integration)
            diagnostic_steps = [step for step in steps
                                if workflow_step_field(step, "name")
                                == ["Preserve gate diagnostics"]]
            evidence_uploads = [step for step in steps
                                if "upload-artifact@" in step and step not in diagnostic_steps]
            if (records_build != 1 or uploads != 2 or len(evidence_uploads) != 1
                    or len(diagnostic_steps) != 1
                    or not gate_diagnostic_upload_valid(diagnostic_steps[0])):
                errors.append(f"native integration lacks one build evidence producer: {name}")
            promoters = [step for step in steps
                         if ("cmake --build" in step and "RunNativeBuildEvidence" in step)
                         or "verify-attestation.py --record-build" in step]
            promoters += evidence_uploads
            if (len(promoters) != 3 or any(
                workflow_step_field(step, "if")
                != ["github.event_name == 'workflow_dispatch'"] for step in promoters
            )):
                errors.append(
                    f"native build evidence is not dispatch-only inside integration: {name}"
                )
            platform = "windows" if expected_domain == "windows-native-build" else "macos"
            if (len(evidence_uploads) != 1
                    or re.findall(r"^          name: ([^\n]+)$", evidence_uploads[0],
                                  re.MULTILINE)
                    != [f"sirius-{platform}-native-attestation"]
                    or re.findall(r"^          path: ([^\n]+)$", evidence_uploads[0],
                                  re.MULTILINE) != ["attestations"]):
                errors.append(f"native build attestation upload has an unsafe scope: {name}")
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
            continue
        if re.findall(r"^    if: ([^\n]+)$", job, re.MULTILINE) != [
            FULL_QUALIFICATION_CONDITION
        ]:
            errors.append(
                f"CI full qualification requires a push or explicit opted-in dispatch: {name}"
            )
        diagnostics = [step for step in workflow_steps(job)
                       if workflow_step_field(step, "name") == ["Preserve gate diagnostics"]]
        if len(diagnostics) != 1 or not gate_diagnostic_upload_valid(
            diagnostics[0], macos_runtime=name == "macos-build"
        ):
            errors.append(f"CI full qualification does not retain bounded gate diagnostics: {name}")
        if name in {"windows-build", "macos-build"}:
            platform = name.removesuffix("-build")
            uploads = [step for step in workflow_steps(job)
                       if "upload-artifact@" in step and step not in diagnostics]
            runtime_uploads = [step for step in uploads if macos_runtime_upload_valid(step)]
            native_uploads = [step for step in uploads if step not in runtime_uploads]
            if (len(uploads) != (2 if name == "macos-build" else 1)
                    or len(runtime_uploads) != (1 if name == "macos-build" else 0)
                    or len(native_uploads) != 1
                    or re.findall(r"^          name: ([^\n]+)$", native_uploads[0], re.MULTILINE)
                    != [f"sirius-{platform}-full-native-attestation"]
                    or re.findall(r"^          path: ([^\n]+)$", native_uploads[0], re.MULTILINE)
                    != ["attestations"]):
                errors.append(f"CI full native upload must have a distinct bounded artifact: {name}")
        if name == "macos-build":
            steps = workflow_steps(job)
            producers = [step for step in steps if "validate-native-runtime.py" in step]
            provenance = [step for step in steps
                          if workflow_step_field(step, "name") == ["Record macOS graphics route"]]
            if (len(producers) != 1
                    or producers[0].rstrip() != MACOS_RUNTIME_PRODUCER_STEP.rstrip()
                    or len(provenance) != 1
                    or provenance[0].rstrip() != MACOS_RUNTIME_PROVENANCE_STEP.rstrip()
                    or steps.index(provenance[0]) >= steps.index(producers[0])):
                errors.append("macOS full qualification must run the exact runtime producer after route capture")
            if (len(producers) == 1 and len(runtime_uploads) == 1
                    and steps.index(runtime_uploads[0]) <= steps.index(producers[0])):
                errors.append("macOS runtime upload precedes its verified producer")
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
        + "    steps:\n"
        + "      - run: cmake -DSIRIUS_ALIGNMENT_MODE=qualification\n"
        + "      - run: cmake --build --target "
        + targets
        + "\n"
        + "      - run: ctest --test-dir bin/integration --no-tests=error -R '^("
        + controls
        + ")$'\n"
        + "      - run: test ! -e mandatory_gate.json\n"
        + "      - run: test ! -e mandatory_gate.json\n"
    )
    native_evidence_body = (
        "      - if: github.event_name == 'workflow_dispatch'\n"
        "        run: cmake --build --target RunNativeBuildEvidence\n"
        "      - if: github.event_name == 'workflow_dispatch'\n"
        "        run: python scripts/verify-attestation.py --record-build "
        "--native-build-gate native_build_gate.json "
        "--native-build-gate-log native_build_gate_ctest.log "
        "--artifact native_build_gate_junit.xml --domain {domain}\n"
        "      - if: github.event_name == 'workflow_dispatch'\n"
        "        uses: actions/upload-artifact@0123456789012345678901234567890123456789\n"
        "        with:\n"
        "          name: sirius-{platform}-native-attestation\n"
        "          path: attestations\n"
    )
    diagnostics = (
        "      - name: Preserve gate diagnostics\n"
        "        if: always()\n"
        "        uses: actions/upload-artifact@0123456789012345678901234567890123456789\n"
        "        with:\n"
        "          name: sirius-gate-diagnostics-${{ github.job }}\n"
        "          path: |\n"
        + "".join(f"            {path}\n" for path in GATE_DIAGNOSTIC_PATHS)
        + "          if-no-files-found: ignore\n"
    )
    macos_runtime_upload = """      - name: Preserve verified MoltenVK runtime evidence
        if: success()
        uses: actions/upload-artifact@043fb46d1a93c77aae656e7c1c64a875d1fc6a0a # v7
        with:
          name: sirius-macos-moltenvk-attestation
          path: ${{ runner.temp }}/sirius-macos-runtime
          if-no-files-found: error
"""
    macos_diagnostics = diagnostics.replace(
        "          if-no-files-found: ignore\n",
        "".join(f"            {path}\n" for path in MACOS_RUNTIME_DIAGNOSTIC_PATHS)
        + "          if-no-files-found: ignore\n",
    )
    valid = (
        "on:\n"
        "  workflow_dispatch:\n"
        "    inputs:\n"
        "      full_qualification:\n"
        "        type: boolean\n"
        "        default: false\n"
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
                    ),
                    platform="windows" if name == "integration-windows-no-render" else "macos",
                )
                + diagnostics
                if name != "integration-no-render"
                else ""
            )
            for name in NON_RENDER_INTEGRATION_JOBS
        )
        + "".join(
            f"  {name}:\n    if: {FULL_QUALIFICATION_CONDITION}\n    steps:\n"
            + (MACOS_RUNTIME_PROVENANCE_STEP + MACOS_RUNTIME_PRODUCER_STEP
               + macos_runtime_upload if name == "macos-build" else "")
            + (
                "      - uses: actions/upload-artifact@0123456789012345678901234567890123456789\n"
                "        with:\n"
                f"          name: sirius-{name.removesuffix('-build')}-full-native-attestation\n"
                "          path: attestations\n"
                if name in {"windows-build", "macos-build"} else ""
            )
            + (macos_diagnostics if name == "macos-build" else diagnostics)
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
        "      - if: github.event_name == 'workflow_dispatch'\n"
        "        run: cmake --build --target RunNativeBuildEvidence\n",
        "      - run: cmake --build --target RunNativeBuildEvidence\n",
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

    mutations = {
        "diagnostic attestation name": valid.replace(
            "name: sirius-gate-diagnostics-${{ github.job }}",
            "name: sirius-windows-native-attestation", 1),
        "diagnostic attestation scope": valid.replace(
            "            bin/**/generated/sirius/*gate*.json\n",
            "            attestations/**\n", 1),
        "diagnostic unbounded scope": valid.replace(
            "            bin/**/generated/sirius/*gate*.json\n", "            **/*\n", 1),
        "diagnostic missing failure guard": valid.replace(
            "        if: always()\n", "        if: success()\n", 1),
        "diagnostic extra execution": valid.replace(
            "      - name: Preserve gate diagnostics\n",
            "      - name: Preserve gate diagnostics\n        run: ctest --preset macos\n", 1),
        "attestation always guard": valid.replace(
            "      - if: github.event_name == 'workflow_dispatch'\n"
            "        uses: actions/upload-artifact@",
            "      - if: always()\n        uses: actions/upload-artifact@", 1),
        "default-on runtime": valid.replace("        default: false", "        default: true", 1),
        "string runtime input": valid.replace("        type: boolean", "        type: string", 1),
        "missing runtime default": valid.replace("        default: false\n", "", 1),
        "duplicated runtime default": valid.replace(
            "        default: false\n", "        default: false\n        default: true\n", 1),
    }
    for condition in (
        "true",
        "github.event_name == 'push' || true",
        "github.event_name == 'push' || inputs.full_qualification == true",
        "github.event_name == 'push' || github.event_name == 'pull_request'",
        "github.event_name == 'push' || github.event_name == 'workflow_dispatch'",
    ):
        mutations[condition] = valid.replace(FULL_QUALIFICATION_CONDITION, condition, 1)
    for platform in ("windows", "macos"):
        mutations[f"{platform} full/integration artifact collision"] = valid.replace(
            f"name: sirius-{platform}-full-native-attestation",
            f"name: sirius-{platform}-native-attestation", 1)
    runtime_mutations = {
        "missing explicit producer pipefail": (
            MACOS_RUNTIME_PRODUCER_STEP,
            MACOS_RUNTIME_PRODUCER_STEP.replace("          set -euo pipefail\n", "", 1)),
        "producer disabled guard": ("      - name: Exact native MoltenVK runtime evidence\n        if: success()", "      - name: Exact native MoltenVK runtime evidence\n        if: false"),
        "missing producer": (MACOS_RUNTIME_PRODUCER_STEP, ""),
        "producer self-test": ("validate-native-runtime.py --expected-revision", "validate-native-runtime.py --self-test --expected-revision"),
        "wrong expected revision": ('--expected-revision "${{ github.sha }}"', '--expected-revision "stale"'),
        "missing expected revision": (' --expected-revision "${{ github.sha }}"', ""),
        "wrong output root": ('--output-root "$RUNNER_TEMP/sirius-macos-runtime"', '--output-root "$RUNNER_TEMP"'),
        "unbounded parallel build": ("CMAKE_BUILD_PARALLEL_LEVEL: 3", "CMAKE_BUILD_PARALLEL_LEVEL: 8"),
        "suppressed producer failure": ('2>&1 | tee "$RUNNER_TEMP/sirius-macos-runtime/producer.log"', '2>&1 | tee "$RUNNER_TEMP/sirius-macos-runtime/producer.log" || true'),
        "producer without pipefail shell": ("        shell: bash\n        env:", "        shell: sh\n        env:"),
        "missing route capture": (MACOS_RUNTIME_PROVENANCE_STEP, ""),
        "missing runtime upload": (macos_runtime_upload, ""),
        "runtime publication on failure": ("      - name: Preserve verified MoltenVK runtime evidence\n        if: success()", "      - name: Preserve verified MoltenVK runtime evidence\n        if: always()"),
        "runtime artifact collision": ("name: sirius-macos-moltenvk-attestation", "name: sirius-macos-full-native-attestation"),
        "unbounded runtime upload": ("          path: ${{ runner.temp }}/sirius-macos-runtime", "          path: ${{ runner.temp }}"),
        "missing runtime failure JUnit": ("            " + MACOS_RUNTIME_DIAGNOSTIC_PATHS[0] + "\n", ""),
        "missing producer console log": ("            " + MACOS_RUNTIME_DIAGNOSTIC_PATHS[2] + "\n", ""),
    }
    for description, (before, after) in runtime_mutations.items():
        mutations["macOS " + description] = valid.replace(before, after, 1)
    mutations["macOS upload before producer"] = valid.replace(
        MACOS_RUNTIME_PRODUCER_STEP + macos_runtime_upload,
        macos_runtime_upload + MACOS_RUNTIME_PRODUCER_STEP, 1)
    for description, mutated in mutations.items():
        if mutated == valid or not integration_boundary_errors(mutated):
            raise RuntimeError(f"integration-boundary policy accepted {description}")


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
