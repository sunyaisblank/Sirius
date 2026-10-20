"""Immutable inputs, CI execution boundaries and qualification source identity."""

from __future__ import annotations

import ast
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


INTEGRATION_CONDITION = (
    "github.event_name == 'pull_request' || "
    "(github.event_name == 'workflow_dispatch' && inputs.macos_retained_diagnostic != true)"
)
FULL_QUALIFICATION_CONDITION = (
    "github.event_name == 'push' || "
    "(github.event_name == 'workflow_dispatch' && inputs.full_qualification == true "
    "&& inputs.macos_retained_diagnostic != true)"
)
MACOS_RETAINED_DIAGNOSTIC_CONDITION = (
    "github.event_name == 'workflow_dispatch' && inputs.macos_retained_diagnostic == true"
)
# Dispatch-only build volumes are diagnostic inputs, separate from attestations.
WINDOWS_DIAGNOSTIC_BUILD_STAGE = r"""      - name: Stage diagnostic Windows build products
        if: github.event_name == 'workflow_dispatch'
        shell: pwsh
        run: |
          $ErrorActionPreference = "Stop"
          $revision = (git rev-parse HEAD).Trim()
          if ($LASTEXITCODE -ne 0 -or $revision -ne $env:GITHUB_SHA -or (git status --porcelain)) { throw "diagnostic source identity mismatch" }
          $alignment = "bin/windows-msvc/generated/sirius/alignment_receipt.json"
          if ((Get-Content -Raw $alignment | ConvertFrom-Json).source_revision -ne $revision) { throw "diagnostic alignment revision mismatch" }
          $volumes = @("bin/windows-msvc/src/sirius/app/Release", "bin/windows-msvc/tests/backend/Release")
          $executables = @("sirius.exe", "sirius_backend_tests.exe")
          $loader = Get-Item -LiteralPath (Join-Path $env:VULKAN_SDK "runtime/x64/vulkan-1.dll")
          $cache = Get-Content -Raw "bin/windows-msvc/CMakeCache.txt"
          if ($cache -notmatch '(?m)^BUILD_SHARED_LIBS:BOOL=OFF\r?$') { throw "diagnostic volumes require the existing static-library topology" }
          if ($cache -match '(?m)^GLFW_LIBRARY_TYPE:STRING=(?!STATIC\r?$)[^\r\n]+') { throw "diagnostic GLFW DLL topology is not represented" }
          $instance = [regex]::Match($cache, '(?m)^CMAKE_GENERATOR_INSTANCE:[^=]+=(.+)$')
          if (-not $instance.Success) { throw "selected Visual Studio instance is absent" }
          $compiler = @(Get-ChildItem -Path "bin/windows-msvc/CMakeFiles/*/CMakeCXXCompiler.cmake" -File)
          if ($compiler.Count -ne 1) { throw "selected C++ compiler metadata is ambiguous" }
          $compilerVersion = [regex]::Match((Get-Content -Raw $compiler[0].FullName), 'set\(CMAKE_CXX_COMPILER_VERSION "(19\.(\d+)\.[^"]+)"\)')
          if (-not $compilerVersion.Success) { throw "selected MSVC compiler version is absent" }
          $redistVersionFile = Join-Path ($instance.Groups[1].Value.Trim()) "VC/Auxiliary/Build/Microsoft.VCRedistVersion.default.txt"
          $redistVersion = (Get-Content -Raw $redistVersionFile).Trim()
          if ($redistVersion -notmatch ('^14\.' + $compilerVersion.Groups[2].Value + '\.')) { throw "selected MSVC and release CRT versions do not match" }
          $crtDirectory = Join-Path ($instance.Groups[1].Value.Trim()) "VC/Redist/MSVC/$redistVersion/x64/Microsoft.VC143.CRT"
          $crtNames = @("concrt140.dll", "msvcp140.dll", "msvcp140_1.dll", "msvcp140_2.dll", "msvcp140_atomic_wait.dll", "msvcp140_codecvt_ids.dll", "vcruntime140.dll", "vcruntime140_1.dll", "Microsoft.VC143.CRT.manifest")
          $crt = @($crtNames | ForEach-Object { $path = Join-Path $crtDirectory $_; if (Test-Path -LiteralPath $path -PathType Leaf) { Get-Item -LiteralPath $path } })
          foreach ($name in @("msvcp140.dll", "vcruntime140.dll", "vcruntime140_1.dll")) {
            if ($name -notin $crt.Name) { throw "selected release CRT missing: $name" }
          }
          $files = @{}
          $origins = @{}
          foreach ($dll in @($loader) + @($crt)) { $origins[$dll.Name] = $dll.FullName }
          for ($index = 0; $index -lt $volumes.Count; $index++) {
            $volume = $volumes[$index]
            foreach ($dll in @($loader) + @($crt)) { Copy-Item -LiteralPath $dll.FullName -Destination $volume -Force }
            $payload = @(Get-Item -LiteralPath (Join-Path $volume $executables[$index]))
            if ($index -eq 0) {
              $payload += @(Get-Item -LiteralPath (Join-Path $volume "sirius_render_tests.exe"))
            }
            $payload += @(Get-ChildItem -Path "$volume/*.dll", "$volume/*.manifest" -File -Force)
            $payload += @(Get-ChildItem -LiteralPath "$volume/resources" -Recurse -File -Force)
            foreach ($file in $payload) {
              $relative = [IO.Path]::GetRelativePath($PWD.Path, $file.FullName).Replace('\', '/')
              $files[$relative] = @{ bytes = $file.Length; sha256 = (Get-FileHash -Algorithm SHA256 -LiteralPath $file.FullName).Hash.ToLowerInvariant() }
            }
          }
          $compilerCopy = "bin/windows-msvc/generated/sirius/windows-diagnostic-CMakeCXXCompiler.cmake"
          Copy-Item -LiteralPath $compiler[0].FullName -Destination $compilerCopy
          foreach ($path in @($alignment, "bin/windows-msvc/CMakeCache.txt", $compilerCopy)) {
            $file = Get-Item -LiteralPath $path
            $files[$path] = @{ bytes = $file.Length; sha256 = (Get-FileHash -Algorithm SHA256 -LiteralPath $path).Hash.ToLowerInvariant() }
          }
          $registrationRoot = "bin/windows-msvc/generated/sirius/windows-mandatory-registration"
          if (Test-Path -LiteralPath $registrationRoot) { throw "registration export already exists" }
          $registrationFiles = @(Get-ChildItem -LiteralPath "bin/windows-msvc" -Recurse -File | Where-Object {
            $_.Name -in @("CTestTestfile.cmake", "cmake_install.cmake") -or $_.Name -match '_(tests|include)(-[^.]+)?\.cmake$'
          })
          foreach ($target in @("base", "core", "oracle", "backend", "app", "render")) {
            foreach ($kind in @("tests", "include")) {
              if (@($registrationFiles | Where-Object { $_.Name -match "^sirius_${target}_tests\[\d+\]_${kind}\.cmake$" }).Count -ne 1) {
                throw "configured GoogleTest discovery is incomplete: $target/$kind"
              }
            }
          }
          $registrationFiles += @(Get-Item -LiteralPath @(
            "bin/windows-msvc/CTestTestfile.cmake", "bin/windows-msvc/cmake_install.cmake",
            "bin/windows-msvc/tests/evidence_governance.cmake",
            "bin/windows-msvc/tests/alignment_receipt_integrity.cmake",
            "bin/windows-msvc/tests/volume_initialization.cmake",
            "bin/windows-msvc/verify_sirius_install.cmake",
            "bin/windows-msvc/tests/source-governance-ctest.json",
            "tests/operational/required_skip_rejection.cmake",
            "tests/operational/vulkan_runtime_profile.cmake", "tests/operational/p1_burn_in.cmake"))
          $configuredRoots = @{}
          foreach ($key in @("CMAKE_HOME_DIRECTORY", "CMAKE_CACHEFILE_DIR", "CMAKE_COMMAND", "CMAKE_CTEST_COMMAND", "CMAKE_ROOT", "_Python3_EXECUTABLE")) {
            $entry = [regex]::Match($cache, "(?m)^${key}:[^=]+=(.+?)\r?$")
            if (-not $entry.Success) { throw "configured root/tool is absent: $key" }
            $configuredRoots[$key] = $entry.Groups[1].Value
          }
          if ([IO.Path]::GetFullPath($configuredRoots.CMAKE_HOME_DIRECTORY) -ne $PWD.Path -or
              [IO.Path]::GetFullPath($configuredRoots.CMAKE_CACHEFILE_DIR) -ne (Resolve-Path "bin/windows-msvc").Path) {
            throw "configured registration roots differ from the exported checkout/build"
          }
          $configuredTools = @{}
          foreach ($key in @("CMAKE_COMMAND", "CMAKE_CTEST_COMMAND", "_Python3_EXECUTABLE")) {
            $tool = Get-Item -LiteralPath $configuredRoots[$key]
            $configuredTools[$key] = @{ path = $tool.FullName; bytes = $tool.Length; sha256 = (Get-FileHash -Algorithm SHA256 -LiteralPath $tool.FullName).Hash.ToLowerInvariant() }
          }
          $registration = @{ schema = "sirius-windows-configured-mandatory-registration-v1"; source_revision = $revision; configuration = "Release"; configured_roots = $configuredRoots; configured_tools = $configuredTools; files = @{}; execution_claimed = $false; scope = "Unmodified configured registration/install recipes and the preceding governance inventory; no Mandatory execution, qualification, or tool substitution admission." }
          foreach ($file in @($registrationFiles | Sort-Object -Property FullName -Unique)) {
            $relative = [IO.Path]::GetRelativePath($PWD.Path, $file.FullName).Replace('\', '/')
            $destination = "$registrationRoot/$relative"
            New-Item -ItemType Directory -Path (Split-Path -Parent $destination) -Force | Out-Null
            Copy-Item -LiteralPath $file.FullName -Destination $destination
            $identity = @{ bytes = $file.Length; sha256 = (Get-FileHash -Algorithm SHA256 -LiteralPath $file.FullName).Hash.ToLowerInvariant() }
            if ((Get-FileHash -Algorithm SHA256 -LiteralPath $destination).Hash.ToLowerInvariant() -ne $identity.sha256) { throw "registration copy changed: $relative" }
            $files[$destination] = $identity
            $registration.files[$relative] = @{ source_path = $file.FullName; artifact_path = $destination; bytes = $identity.bytes; sha256 = $identity.sha256 }
          }
          $registrationMetadata = "$registrationRoot/registration-metadata.json"
          [IO.File]::WriteAllText((Join-Path $PWD.Path $registrationMetadata), ($registration | ConvertTo-Json -Depth 8) + "`n", [Text.UTF8Encoding]::new($false))
          $metadataFile = Get-Item -LiteralPath $registrationMetadata
          $files[$registrationMetadata] = @{ bytes = $metadataFile.Length; sha256 = (Get-FileHash -Algorithm SHA256 -LiteralPath $registrationMetadata).Hash.ToLowerInvariant() }
          $record = @{ schema = "sirius-windows-diagnostic-build-products-v1"; source_revision = $revision; repository = $env:GITHUB_REPOSITORY; run_id = $env:GITHUB_RUN_ID; run_attempt = $env:GITHUB_RUN_ATTEMPT; msvc_version = $compilerVersion.Groups[1].Value; release_crt_version = $redistVersion; runtime_dll_origins = $origins; files = $files; diagnostic_build_only = $true; qualification_claimed = $false; scope = "Strict compiled CLI/backend-test volumes and configured CTest/install recipes only; no native runtime, rendering, UI, scientific, release or native-build attestation claim. Windows OS/UCRT and the physical graphics driver remain host dependencies." }
          $output = "bin/windows-msvc/generated/sirius/windows-diagnostic-build-products.json"
          [IO.File]::WriteAllText((Join-Path $PWD.Path $output), ($record | ConvertTo-Json -Depth 8) + "`n", [Text.UTF8Encoding]::new($false))
"""
WINDOWS_DIAGNOSTIC_BUILD_UPLOAD = r"""      - name: Preserve diagnostic Windows build products
        if: github.event_name == 'workflow_dispatch'
        uses: actions/upload-artifact@043fb46d1a93c77aae656e7c1c64a875d1fc6a0a # v7
        with:
          name: sirius-windows-diagnostic-build-products
          path: |
            bin/windows-msvc/src/sirius/app/Release/sirius.exe
            bin/windows-msvc/src/sirius/app/Release/*.dll
            bin/windows-msvc/src/sirius/app/Release/*.manifest
            bin/windows-msvc/src/sirius/app/Release/resources
            bin/windows-msvc/tests/backend/Release/sirius_backend_tests.exe
            bin/windows-msvc/tests/backend/Release/*.dll
            bin/windows-msvc/tests/backend/Release/*.manifest
            bin/windows-msvc/tests/backend/Release/resources
            bin/windows-msvc/generated/sirius/alignment_receipt.json
            bin/windows-msvc/CMakeCache.txt
            bin/windows-msvc/generated/sirius/windows-diagnostic-CMakeCXXCompiler.cmake
            bin/windows-msvc/generated/sirius/windows-diagnostic-build-products.json
            bin/windows-msvc/generated/sirius/windows-mandatory-registration
          include-hidden-files: true
          if-no-files-found: error
"""

# This one-off diagnostic reuses an immutable native export. Keep its entire
# job bounded: one provider build, four original controls and diagnostic upload.
MACOS_RETAINED_DIAGNOSTIC_JOB = r"""    if: github.event_name == 'workflow_dispatch' && inputs.macos_retained_diagnostic == true
    runs-on: macos-15
    timeout-minutes: 50
    permissions:
      contents: read
      actions: read
    steps:
      - uses: actions/checkout@3d3c42e5aac5ba805825da76410c181273ba90b1 # v7
      - name: Fetch the exported executable's source identity
        run: git fetch --no-tags --depth=1 origin 4bcd9875cd9308deabed06a5fa2f17d2321b0cbc
      - name: Vulkan SDK + MoltenVK diagnostic runtime
        uses: jakoch/install-vulkan-sdk-action@37effcfa045411f8bfbbda26df2fd1b3bf3436fa
        with:
          vulkan_version: 1.4.357.0
          cache: true
          github_token: ${{ secrets.GITHUB_TOKEN }}
      - name: Fetch pinned native export and diagnostic provider sources
        env:
          GH_TOKEN: ${{ secrets.GITHUB_TOKEN }}
        run: |
          mkdir -p out/macos-moltenvk-integration
          gh api repos/sunyaisblank/Sirius/actions/artifacts/11778215779/zip \
            > out/macos-moltenvk-integration/export.zip
          git fetch --no-tags origin 21f87e4806c36cd1ae9e45f912a5817d6878b774
      - name: Build diagnostic MoltenVK and exercise original retained controls
        run: |
          export VK_DRIVER_FILES="$VULKAN_SDK/share/vulkan/icd.d/MoltenVK_icd.json"
          python3 -B scripts/macos_moltenvk_integration.py \
            --archive out/macos-moltenvk-integration/export.zip \
            --output out/macos-moltenvk-integration/observation
      - name: Preserve diagnostic inputs and observations
        if: always()
        uses: actions/upload-artifact@043fb46d1a93c77aae656e7c1c64a875d1fc6a0a # v7
        with:
          name: sirius-macos-retained-diagnostic
          path: |
            out/macos-moltenvk-integration/observation
            bin/macos/retained-diagnostic
            bin/macos/macos-moltenvk-integration
          if-no-files-found: error
"""
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
WINDOWS_VULKAN_PREFLIGHT_DIAGNOSTIC_PATHS = (
    "bin/windows-msvc/generated/sirius/windows_vulkan_preflight_inventory.json",
    "bin/windows-msvc/generated/sirius/windows_vulkan_preflight_identity.json",
    "bin/windows-msvc/generated/sirius/windows_vulkan_preflight.xml",
    "bin/windows-msvc/generated/sirius/windows_vulkan_preflight.log",
)
MACOS_RUNTIME_DIAGNOSTIC_PATHS = (
    '${{ runner.temp }}/sirius-macos-runtime/*/native-runtime-tests.xml',
    '${{ runner.temp }}/sirius-macos-runtime/*/native-runtime-transcript.log',
    '${{ runner.temp }}/sirius-macos-runtime/*/native-provider-modules.json',
    '${{ runner.temp }}/sirius-macos-runtime/*/native-provider-bindings.json',
    'bin/macos/tests/backend/sirius_backend_tests',
    '${{ runner.temp }}/sirius-macos-runtime/producer.log',
    '${{ runner.temp }}/sirius-macos-runtime/os-version.txt',
    '${{ runner.temp }}/sirius-macos-runtime/metal-displays.json',
    '${{ runner.temp }}/sirius-macos-runtime/provider-route.txt',
    'out/macos-moltenvk-integration/observation',
    'bin/macos/macos-moltenvk-integration',
)
MACOS_PROVIDER_PROVISION_STEP = r"""      - name: Provision pinned custom MoltenVK for current-source qualification
        if: success()
        shell: bash
        run: |
          set -euo pipefail
          git fetch --no-tags origin 21f87e4806c36cd1ae9e45f912a5817d6878b774
          python3 -B scripts/macos_moltenvk_provider.py --expected-revision "${{ github.sha }}"
"""
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


def gate_diagnostic_upload_valid(step: str, *, macos_runtime: bool = False,
                                 windows_preflight: bool = False) -> bool:
    # Keep diagnostics unable to publish an attestation or broaden their file
    # scope. Unknown fields/layouts fail closed rather than becoming exceptions.
    paths = GATE_DIAGNOSTIC_PATHS + (MACOS_RUNTIME_DIAGNOSTIC_PATHS if macos_runtime else ())
    if windows_preflight:
        paths += WINDOWS_VULKAN_PREFLIGHT_DIAGNOSTIC_PATHS
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


def default_off_boolean_input_valid(workflow: str, name: str) -> bool:
    header = workflow.split("\njobs:", 1)[0]
    dispatches = re.findall(
        r"^  workflow_dispatch:\n(.*?)(?=^  [A-Za-z0-9_-]+:|\Z)",
        header, re.MULTILINE | re.DOTALL,
    )
    if len(dispatches) != 1:
        return False
    inputs = re.findall(
        rf"^      {re.escape(name)}:\n(.*?)(?=^      [A-Za-z0-9_-]+:|\Z)",
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
    for name in ("full_qualification", "macos_retained_diagnostic"):
        if not default_off_boolean_input_valid(workflow, name):
            errors.append(f"CI {name} dispatch must use a default-off boolean input")
    diagnostic = workflow_job(workflow, "macos-retained-diagnostic")
    if (len(re.findall(r"^  macos-retained-diagnostic:\s*$", workflow, re.MULTILINE)) != 1
            or diagnostic is None
            or diagnostic.rstrip() != MACOS_RETAINED_DIAGNOSTIC_JOB.rstrip()):
        errors.append("macOS retained diagnostic must keep its exact manual-only export, owner and non-promoting upload boundary")
    for name in NON_RENDER_INTEGRATION_JOBS:
        integration = workflow_job(workflow, name)
        if integration is None:
            errors.append(f"CI has no non-render integration job: {name}")
            continue
        if re.findall(r"^    if: ([^\n]+)$", integration, re.MULTILINE) != [
            INTEGRATION_CONDITION
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
            product_steps = [step for step in steps if workflow_step_field(step, "name") in
                             (["Stage diagnostic Windows build products"],
                              ["Preserve diagnostic Windows build products"])]
            windows_products = name == "integration-windows-no-render"
            if windows_products:
                expected_products = [WINDOWS_DIAGNOSTIC_BUILD_STAGE, WINDOWS_DIAGNOSTIC_BUILD_UPLOAD]
                if ([step.rstrip() for step in product_steps]
                        != [step.rstrip() for step in expected_products]
                        or not (integration.find("ctest ") < integration.find(WINDOWS_DIAGNOSTIC_BUILD_STAGE)
                                < integration.find(WINDOWS_DIAGNOSTIC_BUILD_UPLOAD)
                                < integration.find("RunNativeBuildEvidence"))):
                    errors.append("Windows diagnostic products exceed their exact dispatch-only build scope")
            elif product_steps:
                errors.append("Windows diagnostic products escaped their Windows integration job")
            evidence_uploads = [step for step in steps
                                if "upload-artifact@" in step and step not in diagnostic_steps
                                and step not in product_steps]
            if (records_build != 1 or uploads != (3 if windows_products else 2) or len(evidence_uploads) != 1
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
            diagnostics[0], macos_runtime=name == "macos-build",
            windows_preflight=name == "windows-build"
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
            provision = [step for step in steps if "macos_moltenvk_provider.py" in step]
            provenance = [step for step in steps
                          if workflow_step_field(step, "name") == ["Record macOS graphics route"]]
            if (len(provision) != 1
                    or provision[0].rstrip() != MACOS_PROVIDER_PROVISION_STEP.rstrip()
                    or len(producers) != 1
                    or producers[0].rstrip() != MACOS_RUNTIME_PRODUCER_STEP.rstrip()
                    or len(provenance) != 1
                    or provenance[0].rstrip() != MACOS_RUNTIME_PROVENANCE_STEP.rstrip()
                    or steps.index(provenance[0]) >= steps.index(producers[0])
                    or steps.index(provision[0]) >= steps.index(producers[0])):
                errors.append("macOS full qualification must provision the pinned custom provider and capture its route before the exact current-source runtime producer")
            if (len(producers) == 1 and len(runtime_uploads) == 1
                    and steps.index(runtime_uploads[0]) <= steps.index(producers[0])):
                errors.append("macOS runtime upload precedes its verified producer")
    return errors


def verify_integration_boundary_policy() -> None:
    targets = " ".join(INTEGRATION_TARGETS)
    controls = "|".join(name.replace(".", r"\.") for name in INTEGRATION_CONTROLS)
    integration_body = (
        f"    if: {INTEGRATION_CONDITION}\n"
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
    windows_diagnostics = diagnostics.replace(
        "          if-no-files-found: ignore\n",
        "".join(f"            {path}\n" for path in WINDOWS_VULKAN_PREFLIGHT_DIAGNOSTIC_PATHS)
        + "          if-no-files-found: ignore\n",
    )
    valid = (
        "on:\n"
        "  workflow_dispatch:\n"
        "    inputs:\n"
        "      full_qualification:\n"
        "        type: boolean\n"
        "        default: false\n"
        "      macos_retained_diagnostic:\n"
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
            + (WINDOWS_DIAGNOSTIC_BUILD_STAGE + WINDOWS_DIAGNOSTIC_BUILD_UPLOAD
               if name == "integration-windows-no-render" else "")
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
            + (MACOS_PROVIDER_PROVISION_STEP + MACOS_RUNTIME_PROVENANCE_STEP + MACOS_RUNTIME_PRODUCER_STEP
               + macos_runtime_upload if name == "macos-build" else "")
            + (
                "      - uses: actions/upload-artifact@0123456789012345678901234567890123456789\n"
                "        with:\n"
                f"          name: sirius-{name.removesuffix('-build')}-full-native-attestation\n"
                "          path: attestations\n"
                if name in {"windows-build", "macos-build"} else ""
            )
            + (macos_diagnostics if name == "macos-build" else
               windows_diagnostics if name == "windows-build" else diagnostics)
            for name in FULL_QUALIFICATION_JOBS
        )
        + "  macos-retained-diagnostic:\n" + MACOS_RETAINED_DIAGNOSTIC_JOB
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
        "Windows missing diagnostic staging": valid.replace(WINDOWS_DIAGNOSTIC_BUILD_STAGE, "", 1),
        "Windows diagnostic upload before staging": valid.replace(
            WINDOWS_DIAGNOSTIC_BUILD_STAGE + WINDOWS_DIAGNOSTIC_BUILD_UPLOAD,
            WINDOWS_DIAGNOSTIC_BUILD_UPLOAD + WINDOWS_DIAGNOSTIC_BUILD_STAGE, 1),
        "Windows diagnostic PR upload": valid.replace(WINDOWS_DIAGNOSTIC_BUILD_UPLOAD,
            WINDOWS_DIAGNOSTIC_BUILD_UPLOAD.replace("        if: github.event_name == 'workflow_dispatch'\n", ""), 1),
        "Windows diagnostic always upload": valid.replace(WINDOWS_DIAGNOSTIC_BUILD_UPLOAD,
            WINDOWS_DIAGNOSTIC_BUILD_UPLOAD.replace("if: github.event_name == 'workflow_dispatch'", "if: always()"), 1),
        "Windows diagnostic unbounded volume": valid.replace(WINDOWS_DIAGNOSTIC_BUILD_UPLOAD,
            WINDOWS_DIAGNOSTIC_BUILD_UPLOAD.replace("bin/windows-msvc/tests/backend/Release/resources", "bin/windows-msvc/**"), 1),
        "Windows diagnostic missing identity": valid.replace(WINDOWS_DIAGNOSTIC_BUILD_UPLOAD,
            WINDOWS_DIAGNOSTIC_BUILD_UPLOAD.replace("            bin/windows-msvc/generated/sirius/windows-diagnostic-build-products.json\n", ""), 1),
        "Windows diagnostic hidden files omitted": valid.replace(WINDOWS_DIAGNOSTIC_BUILD_UPLOAD,
            WINDOWS_DIAGNOSTIC_BUILD_UPLOAD.replace("include-hidden-files: true", "include-hidden-files: false"), 1),
        "Windows diagnostic product promotion": valid.replace(WINDOWS_DIAGNOSTIC_BUILD_UPLOAD,
            WINDOWS_DIAGNOSTIC_BUILD_UPLOAD.replace("sirius-windows-diagnostic-build-products", "sirius-windows-native-attestation"), 1),
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
        "missing retained diagnostic": valid.replace(
            "  macos-retained-diagnostic:\n" + MACOS_RETAINED_DIAGNOSTIC_JOB, "", 1),
        "duplicated retained diagnostic": valid + "  macos-retained-diagnostic:\n" + MACOS_RETAINED_DIAGNOSTIC_JOB,
    }
    for condition in (
        "true",
        "github.event_name == 'push' || true",
        "github.event_name == 'push' || inputs.full_qualification == true",
        "github.event_name == 'push' || github.event_name == 'pull_request'",
        "github.event_name == 'push' || github.event_name == 'workflow_dispatch'",
        "github.event_name == 'push' || (github.event_name == 'workflow_dispatch' && inputs.macos_retained_diagnostic != true)",
        "github.event_name == 'push' || (github.event_name == 'workflow_dispatch' && inputs.full_qualification == true)",
    ):
        mutations[condition] = valid.replace(FULL_QUALIFICATION_CONDITION, condition, 1)
    for condition in (
        "true",
        "github.event_name == 'pull_request' || github.event_name == 'workflow_dispatch'",
        "github.event_name == 'pull_request' && inputs.macos_retained_diagnostic != true",
    ):
        mutations["integration boundary " + condition] = valid.replace(
            INTEGRATION_CONDITION, condition, 1)
    diagnostic_input = (
        "      macos_retained_diagnostic:\n"
        "        type: boolean\n"
        "        default: false\n"
    )
    for description, replacement in {
        "default on": diagnostic_input.replace("default: false", "default: true"),
        "string type": diagnostic_input.replace("type: boolean", "type: string"),
        "missing default": diagnostic_input.replace("        default: false\n", ""),
        "duplicate default": diagnostic_input + "        default: false\n",
        "missing input": "",
        "duplicate input": diagnostic_input * 2,
    }.items():
        mutations["retained diagnostic input " + description] = valid.replace(
            diagnostic_input, replacement, 1)
    diagnostic_mutations = {
        "unconditional dispatch": (MACOS_RETAINED_DIAGNOSTIC_CONDITION, "github.event_name == 'workflow_dispatch'"),
        "non-manual escape": (MACOS_RETAINED_DIAGNOSTIC_CONDITION, "inputs.macos_retained_diagnostic == true"),
        "PR execution": (MACOS_RETAINED_DIAGNOSTIC_CONDITION, "github.event_name == 'pull_request' || inputs.macos_retained_diagnostic == true"),
        "wrong runner": ("    runs-on: macos-15", "    runs-on: macos-15-large"),
        "unbounded timeout": ("    timeout-minutes: 50", "    timeout-minutes: 180"),
        "write permission": ("      actions: read", "      actions: write"),
        "wrong source": ("origin 4bcd9875cd9308deabed06a5fa2f17d2321b0cbc", "origin main"),
        "wrong export": ("artifacts/11778215779/zip", "artifacts/1/zip"),
        "wrong provider source": ("origin 21f87e4806c36cd1ae9e45f912a5817d6878b774", "origin main"),
        "native compiler proxy": ("scripts/macos_moltenvk_integration.py", "scripts/diagnose-macos-retained.py"),
        "different provider entry point": ("scripts/macos_moltenvk_integration.py", "scripts/macos_metal_probe.py"),
        "wrong SDK": ("vulkan_version: 1.4.357.0", "vulkan_version: latest"),
        "preparation only": ("python3 -B scripts/macos_moltenvk_integration.py", "python3 -B scripts/macos_moltenvk_integration.py --prepare-only"),
        "rebuild": ("    steps:\n", "    steps:\n      - run: cmake --build --preset macos\n"),
        "full estate": ("    steps:\n", "    steps:\n      - run: ctest --preset macos\n"),
        "promotion": ("    steps:\n", "    steps:\n      - run: python3 scripts/verify-attestation.py --record-build\n"),
        "attestation artifact": ("name: sirius-macos-retained-diagnostic", "name: sirius-macos-moltenvk-attestation"),
        "unbounded upload": ("            bin/macos/retained-diagnostic", "            attestations/**"),
        "failure data discarded": ("        if: always()", "        if: success()"),
        "missing data accepted": ("if-no-files-found: error", "if-no-files-found: ignore"),
        "extra execution": ("          if-no-files-found: error\n", "          if-no-files-found: error\n      - run: echo bypass\n"),
    }
    for description, (before, after) in diagnostic_mutations.items():
        mutations["retained diagnostic " + description] = valid.replace(
            MACOS_RETAINED_DIAGNOSTIC_JOB,
            MACOS_RETAINED_DIAGNOSTIC_JOB.replace(before, after, 1), 1)
    for platform in ("windows", "macos"):
        mutations[f"{platform} full/integration artifact collision"] = valid.replace(
            f"name: sirius-{platform}-full-native-attestation",
            f"name: sirius-{platform}-native-attestation", 1)
    runtime_mutations = {
        "missing custom provider provision": (MACOS_PROVIDER_PROVISION_STEP, ""),
        "wrong custom provider source": (MACOS_PROVIDER_PROVISION_STEP,
            MACOS_PROVIDER_PROVISION_STEP.replace("21f87e4806c36cd1ae9e45f912a5817d6878b774", "0" * 40)),
        "export reuse as current-source provider": (MACOS_PROVIDER_PROVISION_STEP,
            MACOS_PROVIDER_PROVISION_STEP.replace("macos_moltenvk_provider.py", "macos_moltenvk_integration.py")),
        "provider SDK mutation": (MACOS_PROVIDER_PROVISION_STEP,
            MACOS_PROVIDER_PROVISION_STEP + '          cp custom.dylib "$VULKAN_SDK/lib/libMoltenVK.dylib"\n'),
        "missing custom provider evidence": ("            out/macos-moltenvk-integration/observation\n", ""),
        "missing explicit producer pipefail": (
            MACOS_RUNTIME_PRODUCER_STEP,
            MACOS_RUNTIME_PRODUCER_STEP.replace("          set -euo pipefail\n", "", 1)),
        "producer disabled guard": ("      - name: Exact native MoltenVK runtime evidence\n        if: success()", "      - name: Exact native MoltenVK runtime evidence\n        if: false"),
        "missing producer": (MACOS_RUNTIME_PRODUCER_STEP, ""),
        "producer self-test": ("validate-native-runtime.py --expected-revision", "validate-native-runtime.py --self-test --expected-revision"),
        "wrong expected revision": (MACOS_RUNTIME_PRODUCER_STEP,
            MACOS_RUNTIME_PRODUCER_STEP.replace('--expected-revision "${{ github.sha }}"', '--expected-revision "stale"', 1)),
        "missing expected revision": (MACOS_RUNTIME_PRODUCER_STEP,
            MACOS_RUNTIME_PRODUCER_STEP.replace(' --expected-revision "${{ github.sha }}"', "", 1)),
        "provision wrong expected revision": (MACOS_PROVIDER_PROVISION_STEP,
            MACOS_PROVIDER_PROVISION_STEP.replace('--expected-revision "${{ github.sha }}"', '--expected-revision "stale"', 1)),
        "provision missing expected revision": (MACOS_PROVIDER_PROVISION_STEP,
            MACOS_PROVIDER_PROVISION_STEP.replace(' --expected-revision "${{ github.sha }}"', "", 1)),
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
        "missing producer console log": ("            " + MACOS_RUNTIME_DIAGNOSTIC_PATHS[5] + "\n", ""),
    }
    for description, (before, after) in runtime_mutations.items():
        mutations["macOS " + description] = valid.replace(before, after, 1)
    for path in WINDOWS_VULKAN_PREFLIGHT_DIAGNOSTIC_PATHS:
        mutations["Windows missing preflight " + path] = valid.replace(
            "            " + path + "\n", "", 1)
    mutations["Windows unbounded preflight diagnostics"] = valid.replace(
        "            " + WINDOWS_VULKAN_PREFLIGHT_DIAGNOSTIC_PATHS[0] + "\n",
        "            bin/windows-msvc/generated/sirius/**\n", 1)
    mutations["preflight diagnostics in integration job"] = valid.replace(
        diagnostics, windows_diagnostics, 1)
    mutations["macOS provision after current-source runtime"] = valid.replace(
        MACOS_PROVIDER_PROVISION_STEP, "", 1).replace(
        MACOS_RUNTIME_PRODUCER_STEP, MACOS_RUNTIME_PRODUCER_STEP + MACOS_PROVIDER_PROVISION_STEP, 1)
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


def native_runtime_preflight_errors(source: str) -> list[str]:
    """Keep complete compilation, exact-candidate preflight and Mandatory in order."""
    errors: list[str] = []
    try:
        tree = ast.parse(source)
        function = next(node for node in tree.body if isinstance(node, ast.FunctionDef)
                        and node.name == "write_runtime_attestation")
        target_values = [ast.literal_eval(node.value) for node in tree.body
                         if isinstance(node, ast.Assign) and any(
                             isinstance(target, ast.Name) and target.id == "NATIVE_COMPILE_TARGETS"
                             for target in node.targets)]
    except (SyntaxError, StopIteration, ValueError, TypeError):
        return ["native runtime compile/preflight authority is missing or malformed"]
    if target_values != [INTEGRATION_TARGETS]:
        errors.append("native runtime preflight compile must retain the complete integration topology")

    def same(node: ast.AST, expression: str) -> bool:
        return ast.dump(node) == ast.dump(ast.parse(expression, mode="eval").body)

    run_calls = sorted((node for node in ast.walk(function)
                        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
                        and node.func.id == "run"), key=lambda node: (node.lineno, node.col_offset))
    expected = (
        "run(configure_command)",
        'run(["cmake", "--build", "--preset", profile.preset, "--target", *NATIVE_COMPILE_TARGETS])',
        'run(["cmake", "--build", "--preset", profile.preset], environment=environment)',
    )
    direct_calls = [node.value for node in function.body if isinstance(node, ast.Expr)]
    if (len(run_calls) < 3 or any(not same(node, expression) or node not in direct_calls
                                 for node, expression in zip(run_calls, expected))):
        errors.append("native runtime must compile explicit targets before preflight and retain its full Mandatory build")

    assignments = [node for node in function.body if isinstance(node, ast.Assign)]

    def assigned(name: str) -> list[ast.Assign]:
        return [node for node in assignments if len(node.targets) == 1
                and isinstance(node.targets[0], ast.Name) and node.targets[0].id == name]

    binary, preflight, environment = (assigned(name) for name in ("binary", "preflight", "environment"))
    selection = [node for node in assignments if len(node.targets) == 1
                 and isinstance(node.targets[0], ast.Subscript)
                 and same(node.targets[0].value, "environment")
                 and same(node.targets[0].slice, '"SIRIUS_VULKAN_DEVICE"')
                 and same(node.value, 'str(preflight["device"]["index"])')]
    if (len(binary) != 1 or not same(binary[0].value, "executable_for(profile)")
            or len(preflight) != 1 or not same(preflight[0].value,
                "run_preflight(argparse.Namespace(profile=profile.domain, candidate=binary, "
                "expected_revision=revision, device_pattern=args.device_pattern))")
            or len(environment) != 1 or not same(environment[0].value, "os.environ.copy()")
            or len(selection) != 1 or len(run_calls) < 3
            or not (run_calls[1].lineno < binary[0].lineno < preflight[0].lineno
                    < environment[0].lineno < selection[0].lineno < run_calls[2].lineno)):
        errors.append("native runtime must directly preflight its exact candidate/profile/revision and bind the selection before Mandatory")
    for name in ("inspect_candidate_authority", "run_preflight"):
        if not any(isinstance(node, ast.ImportFrom) and node.module == "attestation_preflight"
                   and any(alias.name == name and alias.asname is None for alias in node.names)
                   for node in tree.body):
            errors.append("native runtime does not use the existing preflight authority: " + name)

    def checks(expression: str) -> list[ast.Call]:
        return [node for node in direct_calls if isinstance(node, ast.Call)
                and isinstance(node.func, ast.Name) and node.func.id == "require"
                and node.args and same(node.args[0], expression)]

    authority = checks('{"path": str(binary.resolve()), **inspect_candidate_authority(binary, revision)} == preflight["candidate"]')
    copies = [node for node in direct_calls if same(node, "shutil.copy2(binary, qualification_executable_path)")]
    final_selection = [node for node in assigned("selected") if same(node.value,
        "select_physical_device(inventory, profile, args.device_pattern)")]
    recordings = [node for node in direct_calls if isinstance(node, ast.Call)
                  and same(node.func, "inventory_path.write_text")]
    devices = checks('selected == preflight["device"]')
    if (len(authority) != 1 or len(copies) != 1 or len(run_calls) < 3
            or not run_calls[2].lineno < authority[0].lineno < copies[0].lineno):
        errors.append("native runtime must recheck the full existing candidate authority after Mandatory before copying it")
    if (len(final_selection) != 1 or len(recordings) != 1 or len(devices) != 2
            or not any(final_selection[0].lineno < node.lineno < recordings[0].lineno
                       for node in devices)):
        errors.append("native runtime must bind its final inventory to the preflight device before recording it")
    return errors


def verify_native_runtime_preflight_policy(source: str) -> None:
    """Source-only ordering controls; no compiler, candidate or provider execution."""
    if native_runtime_preflight_errors(source):
        raise RuntimeError("native runtime preflight policy rejected its valid producer")
    function = next(node for node in ast.parse(source).body if isinstance(node, ast.FunctionDef)
                    and node.name == "write_runtime_attestation")
    lines = source.splitlines(keepends=True)

    def block(node: ast.AST) -> str:
        return "".join(lines[node.lineno - 1:node.end_lineno])

    runs = [node for node in function.body if isinstance(node, ast.Expr)
            and isinstance(node.value, ast.Call) and isinstance(node.value.func, ast.Name)
            and node.value.func.id == "run"]
    compile_block, mandatory_block = block(runs[1]), block(runs[2])
    preflight_block = block(next(node for node in function.body if isinstance(node, ast.Assign)
                                and any(isinstance(target, ast.Name) and target.id == "preflight"
                                        for target in node.targets)))
    authority_block = block(next(node for node in function.body if isinstance(node, ast.Expr)
                                 and isinstance(node.value, ast.Call)
                                 and any(isinstance(child, ast.Name) and child.id == "inspect_candidate_authority"
                                         for child in ast.walk(node))))
    final_device_block = block(next(node for node in reversed(function.body)
                                    if isinstance(node, ast.Expr) and isinstance(node.value, ast.Call)
                                    and node.value.args and ast.dump(node.value.args[0])
                                    == ast.dump(ast.parse('selected == preflight["device"]', mode="eval").body)))
    mutations = {
        "early default build": source.replace(compile_block,
            '    run(["cmake", "--build", "--preset", profile.preset])\n', 1),
        "early Mandatory build": source.replace(compile_block, mandatory_block + compile_block, 1),
        "missing Mandatory build": source.replace(mandatory_block, "", 1),
        "missing preflight": source.replace(preflight_block, "", 1),
        "late preflight": source.replace(preflight_block, "", 1).replace(
            mandatory_block, mandatory_block + preflight_block, 1),
        "conditional preflight": source.replace(preflight_block,
            "    if False:\n" + "".join("    " + line for line in preflight_block.splitlines(keepends=True)), 1),
        "missing test executable": source.replace('    "sirius_backend_tests",\n', "", 1),
        "Mandatory target in compile stage": source.replace('    "SiriusAlignmentGate",\n',
            '    "RunMandatoryTests",\n', 1),
        "wrong profile": source.replace("profile=profile.domain", 'profile="macos-moltenvk"', 1),
        "wrong candidate": source.replace("candidate=binary", 'candidate=Path("other")', 1),
        "unbound revision": source.replace("expected_revision=revision", "expected_revision=None", 1),
        "ignored requested device": source.replace("device_pattern=args.device_pattern", "device_pattern=None", 1),
        "unbound Mandatory device": source.replace(
            'str(preflight["device"]["index"])', '"0"', 1),
        "unbound Mandatory environment": source.replace(mandatory_block,
            '    run(["cmake", "--build", "--preset", profile.preset])\n', 1),
        "missing post-Mandatory authority": source.replace(authority_block, "", 1),
        "early post-Mandatory authority": source.replace(authority_block, "", 1).replace(
            mandatory_block, authority_block + mandatory_block, 1),
        "unbound post-Mandatory authority": source.replace(
            "inspect_candidate_authority(binary, revision)", "inspect_candidate_authority(binary, None)", 1),
        "missing final device continuity": source.replace(final_device_block, "", 1),
    }
    for label, invalid in mutations.items():
        if invalid == source or not native_runtime_preflight_errors(invalid):
            raise RuntimeError("native runtime preflight policy accepted " + label)


def native_provider_binding_errors(source: str) -> list[str]:
    """Require actual current-source/provider binding around the original gates."""
    try:
        tree = ast.parse(source)
        function = next(node for node in tree.body if isinstance(node, ast.FunctionDef)
                        and node.name == "write_runtime_attestation")
    except (SyntaxError, StopIteration):
        return ["native provider binding producer is malformed"]
    same = lambda node, expression: ast.dump(node) == ast.dump(ast.parse(expression, mode="eval").body)
    statements = function.body
    guards = [node for node in statements if isinstance(node, ast.If)]
    def guarded(condition: str, marker: str) -> list[ast.If]:
        return [node for node in guards if same(node.test, condition)
                and marker in ast.unparse(node)]
    setup = guarded('profile.platform == "macos"', "native_provider.load(revision)")
    binding = guarded("provider is not None", "native_provider.bind_candidate(")
    authority = guarded("provider is not None", "native_provider.verify_binding(provider, provider_binding, revision)")
    frame = guarded("provider is not None", "native_provider.bind_frame(provider, provider_binding, revision, render_output, output)")
    runs = [node for node in statements if isinstance(node, ast.Expr) and isinstance(node.value, ast.Call)
            and isinstance(node.value.func, ast.Name) and node.value.func.id == "run"]
    errors = []
    if (len(setup) != 1 or len(binding) != 1 or len(authority) != 2 or len(frame) != 1
            or len(runs) < 3 or not setup[0].lineno < runs[0].lineno < runs[1].lineno
            < binding[0].lineno < runs[2].lineno < authority[0].lineno < frame[0].lineno < authority[1].lineno):
        errors.append("native current-source provider setup/binding/continuity must surround all original runtime gates")
    if len(setup) == 1:
        expected = ast.parse('''if profile.platform == "macos":
    import macos_moltenvk_provider as native_provider
    provider = native_provider.load(revision)
    native_provider.activate(provider)
''').body[0]
        if ast.dump(setup[0]) != ast.dump(expected):
            errors.append("native provider must activate its verified current-revision receipt on macOS only")
    if len(binding) == 1:
        expected = ast.parse('''if provider is not None:
    provider_binding = native_provider.bind_candidate(provider, revision, binary,
        ROOT / "bin/macos/tests/backend/sirius_backend_tests", preflight, environment, output)
''').body[0]
        if ast.dump(binding[0]) != ast.dump(expected):
            errors.append("native provider must bind the current app/backend and original exact preflight")
    calls = [node for node in ast.walk(function) if isinstance(node, ast.Call)]
    frame_runs = [node for node in calls if isinstance(node.func, ast.Name) and node.func.id == "run"
                  and node.args and isinstance(node.args[0], ast.List)
                  and any(isinstance(item, ast.Constant) and item.value == "render" for item in node.args[0].elts)]
    frame_env = guarded("provider is not None", "render_environment['DYLD_PRINT_LIBRARIES'] = '1'")
    if (len(frame_runs) != 1 or not any(keyword.arg == "environment" and same(keyword.value, "render_environment")
                                       for keyword in frame_runs[0].keywords)
            or len(frame_env) != 1):
        errors.append("original native frame must retain actual provider loaded-image observation")
    if ('native_provider.bundle_proofs(provider, provider_binding, revision, output)' not in source
            or '*provider_artifacts,' not in source):
        errors.append("native provider auxiliary proof bundle must be bound by the existing artifacts hash map")
    if "native_provider_self_test()" not in source or "verify_native_provider_policy(" not in source:
        errors.append("native producer omits portable provider/source controls")
    return errors


def native_provider_source_errors(provider: str, owner: str) -> list[str]:
    """Small source contracts for shared owner reuse and exported-product separation."""
    errors = []
    try:
        tree, owner_tree = ast.parse(provider), ast.parse(owner)
        functions = {node.name: node for node in tree.body if isinstance(node, ast.FunctionDef)}
        owner_functions = {node.name: node for node in owner_tree.body if isinstance(node, ast.FunctionDef)}
        owner_call = owner_functions["run_owned"]
        owner_test = owner_functions["owner_self_test"]
        source_state = owner_functions["source_state"]
        strict = functions["strict_source_snapshot"]
        imports = [node for node in tree.body if isinstance(node, (ast.Import, ast.ImportFrom))]
        require_lazy = functions["diagnostic"]
    except (SyntaxError, KeyError):
        return ["native provider or mature owner source authority is missing/malformed"]
    for function in (owner_call, owner_test):
        keywords = dict(zip((node.arg for node in function.args.kwonlyargs), function.args.kw_defaults))
        default = keywords.get("source_snapshot")
        if not isinstance(default, ast.Constant) or default.value is not None:
            errors.append("shared owner source callback must default to the original diagnostic predicate")
    if ('observe_source = source_state if source_snapshot is None else source_snapshot' not in owner
            or 'source = observe_source()' not in owner
            or 'unchanged = before == after and observe_source() == source' not in owner
            or 'source_snapshot=source_snapshot' not in ast.unparse(owner_test)):
        errors.append("shared owner does not retain before/after default and explicit source checks")
    whitelist = {".github/workflows/ci.yml", "scripts/diagnose-macos-retained.py",
                 "scripts/governance/build_policy.py", "scripts/macos_metal_probe.py",
                 "scripts/macos_moltenvk_integration.py", "scripts/macos_vulkan_identity.py",
                 "scripts/validate-native-runtime.py", "scripts/macos_moltenvk_provider.py", "README.md"}
    sets = [node for node in ast.walk(source_state) if isinstance(node, ast.Set)]
    if (len(sets) != 1 or ast.literal_eval(sets[0]) != whitelist
            or 'changes <= ' not in ast.unparse(source_state)
            or 'candidate_revision' not in ast.unparse(source_state)
            or 'SOURCE' not in ast.unparse(source_state)):
        errors.append("exported-product diagnostic source scope exceeds the explicit nine runner/documentation paths")
    if any("macos_moltenvk_integration" in ast.unparse(node) or "resource" in ast.unparse(node) for node in imports):
        errors.append("native provider imports the Darwin owner outside its lazy macOS path")
    if ('importlib.import_module(\'macos_moltenvk_integration\')' not in ast.unparse(require_lazy)
            or 'sys.platform' not in ast.unparse(require_lazy)):
        errors.append("native provider lacks its bounded lazy native import")
    if ('revision == expected_revision and status == \'\'' not in ast.unparse(strict)
            or '--untracked-files=normal' not in ast.unparse(strict)):
        errors.append("production owner source predicate does not bind exact clean current HEAD")
    for marker in ('source_snapshot=lambda: strict_source_snapshot(revision)',
                   'payload, payload_inputs = d.stage_payload(OUTPUT)',
                   'd.verify_payload(payload)', 'd.tool_inputs(tools)',
                   'verify_pins(receipt["sealed_inputs"])', 'd.loaded_join(inventory_folder, driver, loader)',
                   'd.common.macho_modules(backend)', 'custom["device"]["pipeline_cache_uuid"].startswith(d.CONTENT_ID[:8])',
                   'stock_sdk_modified=False', 'content_id_is_git_provenance=False'):
        if marker not in provider:
            errors.append("native provider omits existing source/tool/load/cache boundary: " + marker)
    if provider.count("source_snapshot=lambda: strict_source_snapshot(revision)") != 2:
        errors.append("both production owner routes must use the exact clean current-source callback")
    if any(marker in provider for marker in ('common.prepare(', 'common.source_state(', 'd.main(', 'common.CASES')):
        errors.append("current-source provisioner reuses the diagnostic executable/source predicate or control loop")
    return errors


def verify_native_provider_policy(producer: str, provider: str, owner: str) -> None:
    if native_provider_binding_errors(producer) or native_provider_source_errors(provider, owner):
        raise RuntimeError("native provider source policy rejected its valid current-source route")
    mutations = (
        (producer.replace("*provider_artifacts,", "", 1), provider, owner),
        (producer.replace('native_provider.load(revision)', 'native_provider.load(None)', 1), provider, owner),
        (producer.replace('native_provider.activate(provider)', 'pass', 1), provider, owner),
        (producer.replace('preflight, environment, output,', 'None, environment, output,', 1), provider, owner),
        (producer.replace('native_provider.verify_binding(provider, provider_binding, revision)', 'pass', 1), provider, owner),
        (producer.replace('native_provider.bind_frame(provider, provider_binding, revision, render_output, output)', 'pass', 1), provider, owner),
        (producer, provider.replace('revision == expected_revision and status == ""', 'True', 1), owner),
        (producer, provider.replace('payload, payload_inputs = d.stage_payload(OUTPUT)', 'payload, payload_inputs = d.common.prepare(OUTPUT)', 1), owner),
        (producer, provider.replace('source_snapshot=lambda: strict_source_snapshot(revision)', 'source_snapshot=None', 1), owner),
        (producer, provider, owner.replace('seconds=300, source_snapshot=None', 'seconds=300, source_snapshot=source_state', 1)),
        (producer, provider, owner.replace('observe_source = source_state if source_snapshot is None else source_snapshot',
                                            'observe_source = source_snapshot', 1)),
        (producer, provider, owner.replace('"README.md"},', '"README.md", "src/sirius/core/camera.h"},', 1)),
    )
    for changed_producer, changed_provider, changed_owner in mutations:
        if ((changed_producer, changed_provider, changed_owner) == (producer, provider, owner)
                or not (native_provider_binding_errors(changed_producer)
                        or native_provider_source_errors(changed_provider, changed_owner))):
            raise RuntimeError("native provider policy accepted an invalid source/provider/owner boundary")


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
    errors.extend(native_runtime_preflight_errors(native_producer))
    errors.extend(native_provider_binding_errors(native_producer))
    if "verify_native_runtime_preflight_policy(Path(__file__).read_text(" not in native_producer:
        errors.append("native producer self-test omits source-only preflight ordering controls")
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
