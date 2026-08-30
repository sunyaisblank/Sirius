#!/usr/bin/env python3
"""Prove that an external-attestation host is ready without rendering.

This producer is intentionally non-promoting.  It executes only Sirius's
machine-readable ``info system`` and ``info readiness`` queries, validates the
exact qualification candidate and host/device boundary, and emits no
``domains`` or ``artifacts`` object accepted by the attestation verifier.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
GIB = 1024 * 1024 * 1024
MINIMUM_OUTPUT_FREE_BYTES = 2 * GIB
MINIMUM_IMAX_DEVICE_BYTES = 2 * GIB
PERMITTED_CANDIDATE_QUERIES = {
    ("--json", "info", "system"),
    ("--json", "info", "readiness"),
}
QUALIFICATION_RESOURCES = (
    "assets/Starfield.png",
    "kernels/trace.spv",
    "kernels/trace_fp32comp.spv",
    "kernels/trace_fp64.spv",
    "model/alignment_receipt.json",
    "model/operating_model.json",
    "shaders/RDSD003A.frag",
    "shaders/RDSD003A.vert",
    "shaders/RDSD004A.frag",
    "shaders/RDSD004A.vert",
    "shaders/RDSD005A.frag",
    "shaders/RDSD005A.vert",
)
READINESS_RESOURCES = (
    "assets/Starfield.png",
    "kernels/trace.spv",
    "kernels/trace_fp32comp.spv",
    "kernels/trace_fp64.spv",
    "model/alignment_receipt.json",
    "model/operating_model.json",
    "shaders/RDSD003A.frag",
    "shaders/RDSD003A.vert",
)


@dataclass(frozen=True)
class Profile:
    key: str
    host_platform: str
    reported_platform: str
    wsl2: bool
    device_pattern: str
    required_driver_tokens: tuple[str, ...]
    forbidden_driver_tokens: tuple[str, ...]
    planned_domains: tuple[str, ...]
    require_viewer: bool = False
    require_dxg: bool = False
    require_imax_memory: bool = False


PROFILES = {
    "wsl2-radeon-viewer": Profile(
        key="wsl2-radeon-viewer",
        host_platform="linux",
        reported_platform="linux",
        wsl2=True,
        device_pattern="Radeon 780M",
        required_driver_tokens=("dozen", "dzn"),
        forbidden_driver_tokens=(),
        planned_domains=(
            "physical-imax-5616x4096",
            "physical-radeon-780m",
            "viewer-native-window-input",
            "wsl2-dozen",
        ),
        require_viewer=True,
        require_dxg=True,
        require_imax_memory=True,
    ),
    "windows-native-vulkan": Profile(
        key="windows-native-vulkan",
        host_platform="win32",
        reported_platform="windows",
        wsl2=False,
        device_pattern="",
        required_driver_tokens=(),
        forbidden_driver_tokens=("dozen", "dzn"),
        planned_domains=("windows-native-vulkan",),
    ),
    "macos-moltenvk": Profile(
        key="macos-moltenvk",
        host_platform="darwin",
        reported_platform="macos",
        wsl2=False,
        device_pattern="",
        required_driver_tokens=("moltenvk",),
        forbidden_driver_tokens=(),
        planned_domains=("macos-moltenvk",),
    ),
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def normalized_host_platform(platform: str) -> str:
    if platform.startswith("linux"):
        return "linux"
    return platform.casefold()


def driver_identity(device: dict) -> str:
    return " ".join(
        str(device.get(field, "")) for field in ("driver_name", "driver_info")
    ).casefold()


def has_any_token(text: str, tokens: tuple[str, ...]) -> bool:
    return any(token in text for token in tokens)


def validate_device_shape(device: dict) -> None:
    require(device.get("kind") in {"integrated", "discrete"},
            "preflight requires an integrated or discrete Vulkan device")
    require(isinstance(device.get("name"), str) and device["name"].strip(),
            "selected Vulkan device has no name")
    require(isinstance(device.get("driver_name"), str)
            and device["driver_name"].strip(),
            "selected Vulkan device has no driver name")
    for field in (
        "index", "vendor_id", "device_id", "api_version", "driver_id",
        "device_local_bytes", "render_memory_bytes",
    ):
        require(type(device.get(field)) is int and device[field] >= 0,
                f"selected Vulkan device has no valid {field}")
    require(device["vendor_id"] > 0 and device["api_version"] > 0
            and device["device_local_bytes"] > 0
            and device["render_memory_bytes"] > 0,
            "selected Vulkan device reports an empty physical capability")
    require(type(device.get("supports_fp64")) is bool,
            "selected Vulkan device has no fp64 capability declaration")


def select_profile_device(
    inventory: dict,
    profile: Profile,
    *,
    host_platform: str,
    device_pattern: str | None,
    tools: dict[str, bool],
    display: str,
    dxg_present: bool,
) -> dict:
    require(normalized_host_platform(host_platform) == profile.host_platform,
            f"profile {profile.key} cannot run on host {host_platform}")
    require(str(inventory.get("platform", "")).casefold()
            == profile.reported_platform,
            "candidate inventory platform does not match the preflight profile")
    require(inventory.get("wsl2") is profile.wsl2,
            "candidate WSL2 identity does not match the preflight profile")

    vulkan = inventory.get("backends", {}).get("vulkan", {})
    require(vulkan.get("compiled") is True and vulkan.get("available") is True,
            "candidate Vulkan backend is not compiled and available")
    devices = vulkan.get("devices")
    require(isinstance(devices, list), "candidate inventory has no Vulkan device list")
    pattern = profile.device_pattern if device_pattern is None else device_pattern
    pattern = pattern.casefold().strip()
    matching = [
        device for device in devices
        if isinstance(device, dict)
        and device.get("kind") in {"integrated", "discrete"}
        and (not pattern or pattern in str(device.get("name", "")).casefold())
    ]
    require(len(matching) == 1,
            "physical Vulkan selection must resolve exactly one device; "
            f"matched {len(matching)}")
    selected = matching[0]
    validate_device_shape(selected)
    identity = driver_identity(selected)
    if profile.required_driver_tokens:
        require(has_any_token(identity, profile.required_driver_tokens),
                "selected Vulkan driver lacks the profile-required identity")
    require(not has_any_token(identity, profile.forbidden_driver_tokens),
            "selected Vulkan driver has an identity forbidden by the profile")
    if profile.require_imax_memory:
        require(selected["render_memory_bytes"] >= MINIMUM_IMAX_DEVICE_BYTES,
                "selected device cannot expose the governed 2048 MiB IMAX budget")
    if profile.require_dxg:
        require(dxg_present, "WSL2 Radeon preflight requires /dev/dxg")
    if profile.require_viewer:
        require(inventory.get("interactive_viewer", {}).get("compiled") is True,
                "consolidated WSL2 preflight requires a compiled native viewer")
        require(bool(display), "consolidated WSL2 viewer preflight requires DISPLAY")
        for command in ("timeout", "xdotool"):
            require(tools.get(command) is True,
                    f"consolidated WSL2 viewer preflight requires {command}")
    return selected


def run_candidate_json(
    candidate: Path,
    arguments: tuple[str, ...],
    environment: dict[str, str],
) -> dict:
    require(arguments in PERMITTED_CANDIDATE_QUERIES,
            "attestation preflight may execute only non-render info queries")
    process = subprocess.run(
        [str(candidate), *arguments],
        cwd=ROOT,
        env=environment,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    require(process.returncode in {0, 1},
            f"candidate info query failed with exit code {process.returncode}: "
            f"{process.stderr.strip()}")
    document = json.loads(process.stdout)
    require(isinstance(document, dict), "candidate info query did not return an object")
    return document


def clean_source_revision(expected_revision: str | None) -> str:
    revision = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, check=True,
        capture_output=True, text=True,
    ).stdout.strip()
    status = subprocess.run(
        ["git", "status", "--porcelain", "--untracked-files=normal"],
        cwd=ROOT, check=True, capture_output=True, text=True,
    ).stdout.strip()
    require(len(revision) in range(40, 65)
            and all(character in "0123456789abcdef" for character in revision),
            "source revision is not a full lowercase hexadecimal revision")
    require(status == "", "attestation preflight requires a clean source tree")
    if expected_revision is not None:
        require(revision == expected_revision,
                "expected revision differs from the preflight source worktree")
    return revision


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def inspect_candidate_authority(candidate: Path, revision: str) -> dict:
    require(candidate.is_file(), f"qualification candidate is missing: {candidate}")
    resource_root = candidate.parent / "resources"
    for relative in QUALIFICATION_RESOURCES:
        path = resource_root / relative
        require(path.is_file() and path.stat().st_size > 0,
                f"qualification resource is missing or empty: {relative}")
    receipt_path = resource_root / "model" / "alignment_receipt.json"
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    require(receipt.get("method") == "revision-bound-attestation-set",
            "candidate alignment receipt has the wrong authority method")
    require(receipt.get("alignment_mode") == "qualification",
            "external preflight requires a qualification candidate")
    require(receipt.get("source_revision") == revision
            and receipt.get("source_tree_clean") is True,
            "candidate alignment receipt is not bound to this clean revision")
    ideal = receipt.get("ultimate_ideal", {})
    require(ideal.get("required") is True,
            "candidate alignment receipt does not require the ultimate ideal")
    model_path = resource_root / "model" / "operating_model.json"
    require(receipt.get("operating_model", {}).get("sha256") == sha256(model_path),
            "candidate operating model differs from its compiled alignment receipt")
    return {
        "bytes": candidate.stat().st_size,
        "sha256": sha256(candidate),
        "alignment_receipt_sha256": sha256(receipt_path),
        "operating_model_sha256": sha256(model_path),
    }


def inspect_readiness(readiness: dict, selected_index: int, revision: str) -> dict:
    ideal = readiness.get("ultimate_ideal", {})
    require(ideal.get("source_revision") == revision
            and ideal.get("source_tree_clean") is True,
            "readiness is not bound to the preflight source revision")
    require(ideal.get("alignment_mode") == "qualification",
            "readiness does not belong to a qualification candidate")
    require(readiness.get("ready") is False and readiness.get("release_ready") is False,
            "a non-promoting preflight must not report release readiness")
    vulkan = readiness.get("backends", {}).get("vulkan", {})
    require(vulkan.get("compiled") is True,
            "readiness does not report a compiled Vulkan backend")
    require(vulkan.get("selected_device_index") == selected_index,
            "readiness selected a different Vulkan device")
    resources = readiness.get("resources", {})
    for relative in READINESS_RESOURCES:
        resource = resources.get(relative, {})
        require(resource.get("available") is True and resource.get("valid") is True,
                f"readiness rejects qualification resource: {relative}")
    gate = readiness.get("mandatory_build_gate", {})
    return {
        "evidence_generation_ready": readiness.get("evidence_generation_ready") is True,
        "mandatory_gate_valid": gate.get("valid") is True,
        "mandatory_gate_executed_tests": gate.get("executed_tests", 0),
        "mandatory_gate_registered_tests": gate.get("registered_tests", 0),
        "vulkan_ready": vulkan.get("ready") is True,
    }


def output_path_is_safe(path: Path) -> bool:
    resolved = path.resolve()
    return not resolved.is_relative_to(ROOT.resolve())


def make_report(
    *,
    profile: Profile,
    revision: str,
    candidate: Path,
    candidate_authority: dict,
    inventory: dict,
    selected: dict,
    readiness: dict,
    tools: dict[str, bool],
    output_free_bytes: int,
) -> dict:
    return {
        "schema_version": 1,
        "kind": "sirius-attestation-preflight",
        "status": "ready-for-external-execution",
        "completed_utc": datetime.now(timezone.utc).isoformat(),
        "source_revision": revision,
        "profile": profile.key,
        "planned_domains": list(profile.planned_domains),
        "promotable": False,
        "external_execution_completed": False,
        "candidate": {
            "path": str(candidate.resolve()),
            **candidate_authority,
        },
        "platform": {
            "os": inventory["platform"],
            "wsl2": inventory["wsl2"],
        },
        "device": selected,
        "prerequisites": {
            "tools": tools,
            "output_free_bytes": output_free_bytes,
            "explicit_vulkan_driver_manifest": bool(
                os.environ.get("VK_DRIVER_FILES") or os.environ.get("VK_ICD_FILENAMES")
            ),
        },
        "readiness": readiness,
        "queries_executed": ["info system", "info readiness"],
        "pending_actions": [
            "complete the domain-specific full test authority",
            "execute the domain-specific physical runtime and image witnesses",
            "record and independently verify the revision-bound attestation",
        ],
    }


def run_preflight(args: argparse.Namespace) -> dict:
    profile = PROFILES[args.profile]
    revision = clean_source_revision(args.expected_revision)
    candidate = args.candidate.resolve()
    candidate_authority = inspect_candidate_authority(candidate, revision)
    required_tools = {"cmake", "ctest", "git"}
    if profile.require_viewer:
        required_tools.update({"timeout", "xdotool"})
    tools = {command: shutil.which(command) is not None for command in sorted(required_tools)}
    missing = [command for command, available in tools.items() if not available]
    require(not missing, f"preflight commands are unavailable: {', '.join(missing)}")
    output_free_bytes = shutil.disk_usage(candidate.parent).free
    require(output_free_bytes >= MINIMUM_OUTPUT_FREE_BYTES,
            "candidate volume has less than 2 GiB free for external evidence")

    environment = os.environ.copy()
    initial_inventory = run_candidate_json(
        candidate, ("--json", "info", "system"), environment
    )
    selected = select_profile_device(
        initial_inventory,
        profile,
        host_platform=sys.platform,
        device_pattern=args.device_pattern,
        tools=tools,
        display=environment.get("DISPLAY", ""),
        dxg_present=Path("/dev/dxg").exists(),
    )
    environment["SIRIUS_VULKAN_DEVICE"] = str(selected["index"])
    inventory = run_candidate_json(
        candidate, ("--json", "info", "system"), environment
    )
    selected = select_profile_device(
        inventory,
        profile,
        host_platform=sys.platform,
        device_pattern=args.device_pattern,
        tools=tools,
        display=environment.get("DISPLAY", ""),
        dxg_present=Path("/dev/dxg").exists(),
    )
    require(inventory["backends"]["vulkan"].get("selected_device_index")
            == selected["index"],
            "candidate did not preserve the selected Vulkan device")
    readiness_document = run_candidate_json(
        candidate, ("--json", "info", "readiness"), environment
    )
    readiness = inspect_readiness(readiness_document, selected["index"], revision)
    require(clean_source_revision(args.expected_revision) == revision,
            "source identity changed during attestation preflight")
    return make_report(
        profile=profile,
        revision=revision,
        candidate=candidate,
        candidate_authority=candidate_authority,
        inventory=inventory,
        selected=selected,
        readiness=readiness,
        tools=tools,
        output_free_bytes=output_free_bytes,
    )


def self_test() -> None:
    base_device = {
        "index": 0,
        "name": "Microsoft Direct3D12 (AMD Radeon 780M Graphics)",
        "kind": "integrated",
        "driver_name": "Dozen",
        "driver_info": "Mesa dzn",
        "vendor_id": 0x1002,
        "device_id": 0x1900,
        "api_version": 1,
        "driver_id": 23,
        "device_local_bytes": 4 * GIB,
        "render_memory_bytes": 4 * GIB,
        "supports_fp64": True,
    }

    def inventory(platform: str, wsl2: bool, device: dict, viewer: bool = True) -> dict:
        return {
            "platform": platform,
            "wsl2": wsl2,
            "interactive_viewer": {"compiled": viewer},
            "backends": {
                "vulkan": {
                    "compiled": True,
                    "available": True,
                    "selected_device_index": device["index"],
                    "devices": [device],
                }
            },
        }

    viewer_tools = {"cmake": True, "ctest": True, "git": True,
                    "timeout": True, "xdotool": True}
    wsl_profile = PROFILES["wsl2-radeon-viewer"]
    require(
        select_profile_device(
            inventory("Linux", True, base_device), wsl_profile,
            host_platform="linux", device_pattern=None, tools=viewer_tools,
            display=":0", dxg_present=True,
        ) == base_device,
        "WSL2 Radeon/viewer positive control failed",
    )
    windows_device = {**base_device, "name": "Physical Windows GPU",
                      "driver_name": "AMD Vulkan", "driver_info": "native"}
    macos_device = {**base_device, "name": "Apple GPU",
                    "driver_name": "MoltenVK", "driver_info": "native"}
    require(
        select_profile_device(
            inventory("Windows", False, windows_device),
            PROFILES["windows-native-vulkan"], host_platform="win32",
            device_pattern=None, tools={}, display="", dxg_present=False,
        ) == windows_device,
        "native Windows positive control failed",
    )
    require(
        select_profile_device(
            inventory("macOS", False, macos_device), PROFILES["macos-moltenvk"],
            host_platform="darwin", device_pattern=None, tools={}, display="",
            dxg_present=False,
        ) == macos_device,
        "MoltenVK positive control failed",
    )

    controls = (
        ("software device", inventory("Linux", True, {**base_device, "kind": "software"}),
         wsl_profile, "linux", viewer_tools, ":0", True),
        ("wrong WSL driver", inventory(
            "Linux", True,
            {**base_device, "driver_name": "RADV", "driver_info": "Mesa native"},
        ),
         wsl_profile, "linux", viewer_tools, ":0", True),
        ("missing dxg", inventory("Linux", True, base_device),
         wsl_profile, "linux", viewer_tools, ":0", False),
        ("missing display", inventory("Linux", True, base_device),
         wsl_profile, "linux", viewer_tools, "", True),
        ("missing input tool", inventory("Linux", True, base_device),
         wsl_profile, "linux", {**viewer_tools, "xdotool": False}, ":0", True),
        ("missing viewer", inventory("Linux", True, base_device, viewer=False),
         wsl_profile, "linux", viewer_tools, ":0", True),
        ("insufficient IMAX memory", inventory(
            "Linux", True, {**base_device, "render_memory_bytes": GIB}
        ), wsl_profile, "linux", viewer_tools, ":0", True),
        ("native Windows Dozen", inventory("Windows", False, base_device),
         PROFILES["windows-native-vulkan"], "win32", {}, "", False),
        ("macOS without MoltenVK", inventory("macOS", False, windows_device),
         PROFILES["macos-moltenvk"], "darwin", {}, "", False),
        ("wrong executing host", inventory("Windows", False, windows_device),
         PROFILES["windows-native-vulkan"], "linux", {}, "", False),
    )
    for label, document, profile, host, tools, display, dxg in controls:
        try:
            select_profile_device(
                document, profile, host_platform=host, device_pattern=None,
                tools=tools, display=display, dxg_present=dxg,
            )
        except ValueError:
            continue
        raise ValueError(f"attestation preflight accepted false control: {label}")

    ambiguous = inventory("Windows", False, windows_device)
    ambiguous["backends"]["vulkan"]["devices"].append(
        {**windows_device, "index": 1}
    )
    try:
        select_profile_device(
            ambiguous, PROFILES["windows-native-vulkan"], host_platform="win32",
            device_pattern=None, tools={}, display="", dxg_present=False,
        )
    except ValueError:
        pass
    else:
        raise ValueError("attestation preflight accepted an ambiguous device selection")

    try:
        run_candidate_json(Path("unused"), ("render",), {})
    except ValueError:
        pass
    else:
        raise ValueError("attestation preflight accepted a render command")

    report = make_report(
        profile=wsl_profile,
        revision="a" * 40,
        candidate=Path("candidate"),
        candidate_authority={"bytes": 1, "sha256": "b" * 64,
                             "alignment_receipt_sha256": "c" * 64,
                             "operating_model_sha256": "d" * 64},
        inventory=inventory("Linux", True, base_device),
        selected=base_device,
        readiness={"evidence_generation_ready": False,
                   "mandatory_gate_valid": False, "vulkan_ready": False},
        tools=viewer_tools,
        output_free_bytes=4 * GIB,
    )
    require(report.get("promotable") is False
            and report.get("external_execution_completed") is False
            and report.get("status") != "pass"
            and "domains" not in report and "artifacts" not in report,
            "preflight report could be confused with a promotable attestation")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--profile", choices=tuple(PROFILES))
    parser.add_argument("--candidate", type=Path)
    parser.add_argument("--device-pattern")
    parser.add_argument("--expected-revision")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    try:
        if args.self_test:
            self_test()
            print("attestation preflight rejected false host, device, and promotion controls")
            return
        require(args.profile is not None and args.candidate is not None,
                "--profile and --candidate are required")
        if args.output is not None:
            require(output_path_is_safe(args.output),
                    "preflight output must remain outside the source worktree")
        report = run_preflight(args)
        payload = json.dumps(report, indent=2) + "\n"
        if args.output is None:
            print(payload, end="")
        else:
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(payload, encoding="utf-8")
            print(f"attestation preflight ready: {args.output}")
    except (OSError, ValueError, json.JSONDecodeError, subprocess.SubprocessError) as error:
        parser.exit(1, f"attestation preflight: {error}\n")


if __name__ == "__main__":
    main()
