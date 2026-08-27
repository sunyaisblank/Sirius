#!/usr/bin/env python3
"""Produce native Windows-Vulkan or macOS-MoltenVK runtime evidence.

This is an external-domain runbook, not a local substitute.  It runs only on
the native host being attested, binds every artifact to a clean Git revision,
and delegates final admission to verify-attestation.py.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
import shutil
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
VERIFIER_PATH = ROOT / "scripts" / "verify-attestation.py"
QUALIFICATION_RUNTIME_PRODUCTS = {
    "operating_model": Path("model/operating_model.json"),
    "starfield": Path("assets/Starfield.png"),
    "trace_fp32comp_spv": Path("kernels/trace_fp32comp.spv"),
    "trace_fp64_spv": Path("kernels/trace_fp64.spv"),
    "trace_spv": Path("kernels/trace.spv"),
    "viewer_rdsd003a_fragment": Path("shaders/RDSD003A.frag"),
    "viewer_rdsd003a_vertex": Path("shaders/RDSD003A.vert"),
    "viewer_rdsd004a_fragment": Path("shaders/RDSD004A.frag"),
    "viewer_rdsd004a_vertex": Path("shaders/RDSD004A.vert"),
    "viewer_rdsd005a_fragment": Path("shaders/RDSD005A.frag"),
    "viewer_rdsd005a_vertex": Path("shaders/RDSD005A.vert"),
}


@dataclass(frozen=True)
class RuntimeProfile:
    domain: str
    platform: str
    preset: str


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def profile_for_host(host: str) -> RuntimeProfile:
    if host == "win32":
        return RuntimeProfile("windows-native-vulkan", "windows", "windows-msvc")
    if host == "darwin":
        return RuntimeProfile("macos-moltenvk", "macos", "macos")
    raise ValueError(
        "native runtime attestation must run on Windows or macOS, not " + host
    )


def load_verifier():
    spec = importlib.util.spec_from_file_location("sirius_attestation", VERIFIER_PATH)
    require(spec is not None and spec.loader is not None,
            "could not load the attestation verifier")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def run(
    command: list[str],
    *,
    environment: dict[str, str] | None = None,
    accepted_returncodes: frozenset[int] = frozenset({0}),
    print_output: bool = True,
) -> str:
    print("+ " + " ".join(command), flush=True)
    process = subprocess.run(
        command,
        cwd=ROOT,
        env=environment,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    if print_output:
        print(process.stdout, end="", flush=True)
    if process.returncode not in accepted_returncodes:
        raise ValueError(
            f"command failed with exit code {process.returncode}: {' '.join(command)}"
        )
    return process.stdout


def clean_source_revision(expected_revision: str | None = None) -> str:
    try:
        revision = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=ROOT,
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()
        status = subprocess.run(
            ["git", "status", "--porcelain", "--untracked-files=normal"],
            cwd=ROOT,
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()
    except (OSError, subprocess.CalledProcessError) as error:
        raise ValueError(f"could not inspect the native source worktree: {error}") from error
    require(len(revision) in range(40, 65) and all(c in "0123456789abcdef" for c in revision),
            "native source revision is not a full lowercase hexadecimal revision")
    require(status == "", "native runtime attestation requires a clean source tree")
    if expected_revision is not None:
        require(revision == expected_revision,
                "expected revision differs from the native source worktree")
    return revision


def select_physical_device(
    inventory: dict, profile: RuntimeProfile, device_pattern: str | None
) -> dict:
    require(str(inventory.get("platform", "")).casefold() == profile.platform,
            "binary platform does not match the native attestation profile")
    require(inventory.get("wsl2") is False,
            "native Windows/macOS evidence may not be produced through WSL2")
    vulkan = inventory.get("backends", {}).get("vulkan", {})
    require(vulkan.get("compiled") is True and vulkan.get("available") is True,
            "the native Vulkan backend is not compiled and available")
    devices = vulkan.get("devices")
    require(isinstance(devices, list), "native Vulkan inventory has no device list")
    physical = [
        device for device in devices
        if isinstance(device, dict) and device.get("kind") in {"integrated", "discrete"}
    ]
    if device_pattern:
        pattern = device_pattern.casefold()
        physical = [
            device for device in physical
            if pattern in str(device.get("name", "")).casefold()
        ]
    require(len(physical) == 1,
            "native runtime selection must resolve exactly one physical Vulkan device")
    selected = physical[0]
    driver_identity = " ".join(
        str(selected.get(field, ""))
        for field in ("driver_name", "driver_info")
    ).casefold()
    if profile.domain == "windows-native-vulkan":
        require("dozen" not in driver_identity and "dzn" not in driver_identity,
                "native Windows Vulkan evidence may not use Dozen/dzn")
    else:
        require("moltenvk" in driver_identity,
                "macOS native runtime evidence requires MoltenVK identity")
    require(type(selected.get("index")) is int and selected["index"] >= 0,
            "selected Vulkan device has no valid index")
    return selected


def executable_for(profile: RuntimeProfile) -> Path:
    base = ROOT / "bin" / profile.preset / "src" / "sirius" / "app"
    if profile.platform == "windows":
        return base / "Release" / "sirius.exe"
    return base / "sirius"


def artifact_record(path: Path) -> dict:
    payload = path.read_bytes()
    return {
        "path": path.name,
        "bytes": len(payload),
        "sha256": hashlib.sha256(payload).hexdigest(),
    }


def copy_qualification_products(resource_root: Path, output: Path) -> list[Path]:
    bundled = []
    for logical_name, relative_path in QUALIFICATION_RUNTIME_PRODUCTS.items():
        source = resource_root / relative_path
        require(source.is_file(), f"qualification product is missing: {source}")
        destination = output / f"qualification-product-{logical_name}"
        shutil.copy2(source, destination)
        bundled.append(destination)
    return bundled


def write_runtime_attestation(args: argparse.Namespace) -> Path:
    profile = profile_for_host(sys.platform)
    for command in ("cmake", "ctest", "git"):
        require(shutil.which(command) is not None,
                f"required native-attestation command is unavailable: {command}")
    revision = clean_source_revision(args.expected_revision)

    run_id = datetime.now(timezone.utc).strftime("%Y%m%d-%H%M%S-%f")
    output = args.output_root.resolve() / run_id
    output.mkdir(parents=True, exist_ok=False)
    inventory_path = output / "device-inventory.json"
    report_path = output / "native-runtime-tests.xml"
    ctest_inventory_path = output / "ctest-inventory.json"
    transcript_path = output / "native-runtime-transcript.log"
    frame_path = output / "native-runtime-frame.png"
    attestation_path = output / "attestation.json"
    alignment_receipt_path = output / "alignment_receipt.json"
    mandatory_gate_path = output / "mandatory_gate.json"
    qualification_executable_path = output / "qualification-sirius.bin"
    qualification_gate_junit_path = output / "qualification-gate-junit"
    qualification_gate_log_path = output / "qualification-gate-log"

    configure_command = [
        "cmake",
        "--preset",
        profile.preset,
        "-DSIRIUS_ALIGNMENT_MODE=qualification",
        "-DSIRIUS_REQUIRE_VULKAN_RUNTIME=ON",
    ]
    if profile.platform == "windows":
        configure_command.append("-DCMAKE_CONFIGURATION_TYPES=Release")
    run(configure_command)
    run(["cmake", "--build", "--preset", profile.preset])
    binary = executable_for(profile)
    require(binary.is_file(), f"native Sirius executable is missing: {binary}")
    shutil.copy2(binary, qualification_executable_path)
    qualification_products = copy_qualification_products(binary.parent / "resources", output)
    configured_receipt = (
        ROOT / "bin" / profile.preset / "generated" / "sirius" / "alignment_receipt.json"
    )
    require(configured_receipt.is_file(),
            f"compiled-authority alignment receipt is missing: {configured_receipt}")
    shutil.copy2(configured_receipt, alignment_receipt_path)
    configured_gate = (
        ROOT / "bin" / profile.preset / "generated" / "sirius" / "mandatory_gate.json"
    )
    require(configured_gate.is_file(),
            f"qualification build-gate receipt is missing: {configured_gate}")
    shutil.copy2(configured_gate, mandatory_gate_path)
    configured_gate_junit = configured_gate.with_name("mandatory_gate_junit.xml")
    configured_gate_log = configured_gate.with_name("mandatory_gate_ctest.log")
    require(configured_gate_junit.is_file() and configured_gate_log.is_file(),
            "qualification build-gate JUnit or log artifact is missing")
    shutil.copy2(configured_gate_junit, qualification_gate_junit_path)
    shutil.copy2(configured_gate_log, qualification_gate_log_path)

    environment = os.environ.copy()
    initial_inventory = json.loads(
        run([str(binary), "--json", "info", "system"], environment=environment)
    )
    selected = select_physical_device(initial_inventory, profile, args.device_pattern)
    environment["SIRIUS_VULKAN_DEVICE"] = str(selected["index"])

    inventory = json.loads(
        run([str(binary), "--json", "info", "system"], environment=environment)
    )
    selected = select_physical_device(inventory, profile, args.device_pattern)
    require(
        inventory["backends"]["vulkan"].get("selected_device_index")
        == selected["index"],
        "selected device did not remain authoritative in the bound inventory",
    )
    inventory_path.write_text(json.dumps(inventory, indent=2) + "\n", encoding="utf-8")

    readiness = json.loads(
        run(
            [str(binary), "--json", "info", "readiness"],
            environment=environment,
            accepted_returncodes=frozenset({0, 1}),
        )
    )
    vulkan_readiness = readiness.get("backends", {}).get("vulkan", {})
    require(
        readiness.get("evidence_generation_ready") is True
        and vulkan_readiness.get("ready") is True
        and vulkan_readiness.get("selected_device_index") == selected["index"],
        "native readiness did not bind the selected physical Vulkan device",
    )

    ctest_command = [
        "ctest",
        "--test-dir",
        str(ROOT / "bin" / profile.preset),
        "--output-on-failure",
        "--output-junit",
        str(report_path),
    ]
    if profile.platform == "windows":
        ctest_command[3:3] = ["-C", "Release"]
    inventory_command = [
        "ctest",
        "--test-dir",
        str(ROOT / "bin" / profile.preset),
        "--show-only=json-v1",
    ]
    if profile.platform == "windows":
        inventory_command[3:3] = ["-C", "Release"]
    ctest_inventory_path.write_text(
        run(inventory_command, environment=environment, print_output=False),
        encoding="utf-8",
    )
    run(ctest_command, environment=environment)
    verifier = load_verifier()
    report = verifier.inspect_junit(report_path)
    require(report["skipped"] == 0,
            "native runtime evidence requires a non-skipping test estate")

    render_output = run(
        [
            str(binary),
            "render",
            "--backend",
            "vulkan",
            "-m",
            "Kerr",
            "-a",
            "0.9",
            "-w",
            "512",
            "-h",
            "512",
            "-s",
            "1",
            "-d",
            "20",
            "-i",
            "64",
            "--no-disk",
            "-o",
            str(frame_path),
        ],
        environment=environment,
    )
    require(frame_path.is_file() and frame_path.stat().st_size > 0,
            "native Vulkan probe emitted no PNG")

    transcript = (
        f"== source revision: {revision}\n"
        f"== [1/6] configure + build ({profile.preset})\n"
        "== readiness: evidence_generation=true vulkan=true "
        f"selected_device_index={selected['index']}\n"
        f"100% tests passed, 0 tests failed out of {report['cases']}\n"
        "== native Vulkan probe output\n"
        + render_output
    )
    transcript_path.write_text(transcript, encoding="utf-8")

    # Output is ignored by the default runbook location; a custom in-tree path
    # is accepted only when it likewise leaves the tested source identity clean.
    require(clean_source_revision(args.expected_revision) == revision,
            "source revision changed during native runtime validation")

    artifacts = {
        path.name: artifact_record(path)
        for path in (
            frame_path,
            report_path,
            inventory_path,
            transcript_path,
            alignment_receipt_path,
            mandatory_gate_path,
            ctest_inventory_path,
            qualification_executable_path,
            qualification_gate_junit_path,
            qualification_gate_log_path,
            *qualification_products,
        )
    }
    document = {
        "schema_version": 1,
        "status": "pass",
        "domains": [profile.domain],
        "source_revision": revision,
        "completed_utc": datetime.now(timezone.utc).isoformat(),
        "platform": {"os": profile.platform, "wsl2": False},
        "device": selected,
        "claims": {
            "preset": profile.preset,
            "test_estate_passed": True,
            "test_report": report,
            "qualification_executable": {
                "artifact": qualification_executable_path.name,
                "bytes": qualification_executable_path.stat().st_size,
                "sha256": hashlib.sha256(
                    qualification_executable_path.read_bytes()
                ).hexdigest(),
            },
            "runtime_ready": True,
        },
        "artifacts": artifacts,
    }
    attestation_path.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")
    verifier.verify_document(document, attestation_path)
    require(clean_source_revision(args.expected_revision) == revision,
            "source identity changed while issuing native runtime evidence")
    return attestation_path


def self_test() -> None:
    windows = profile_for_host("win32")
    macos = profile_for_host("darwin")
    base_device = {
        "index": 0,
        "name": "Physical GPU",
        "kind": "integrated",
        "driver_name": "Native Vulkan",
        "driver_info": "runtime",
    }

    def inventory(platform: str, device: dict) -> dict:
        return {
            "platform": platform,
            "wsl2": False,
            "backends": {
                "vulkan": {
                    "compiled": True,
                    "available": True,
                    "devices": [device],
                }
            },
        }

    require(select_physical_device(inventory("Windows", base_device), windows, None)
            == base_device, "native Windows positive control failed")
    moltenvk_device = {**base_device, "driver_name": "MoltenVK"}
    require(select_physical_device(inventory("macOS", moltenvk_device), macos, None)
            == moltenvk_device, "MoltenVK positive control failed")
    controls = [
        ("unsupported host", lambda: profile_for_host("linux")),
        (
            "Windows through Dozen",
            lambda: select_physical_device(
                inventory("Windows", {**base_device, "driver_name": "Mesa Dozen"}),
                windows,
                None,
            ),
        ),
        (
            "macOS without MoltenVK",
            lambda: select_physical_device(inventory("macOS", base_device), macos, None),
        ),
        (
            "software Vulkan",
            lambda: select_physical_device(
                inventory("Windows", {**base_device, "kind": "software"}),
                windows,
                None,
            ),
        ),
        (
            "ambiguous physical devices",
            lambda: select_physical_device(
                {
                    **inventory("Windows", base_device),
                    "backends": {
                        "vulkan": {
                            "compiled": True,
                            "available": True,
                            "devices": [base_device, {**base_device, "index": 1}],
                        }
                    },
                },
                windows,
                None,
            ),
        ),
    ]
    for label, control in controls:
        try:
            control()
        except ValueError:
            continue
        raise ValueError(f"native runtime negative control was accepted: {label}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output-root",
        type=Path,
        default=ROOT / "renders" / "native-runtime-validation",
    )
    parser.add_argument("--device-pattern")
    parser.add_argument("--expected-revision")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    try:
        if args.self_test:
            self_test()
            print("native runtime producer rejected wrong-host and false-device controls")
            return
        attestation = write_runtime_attestation(args)
        print(f"native runtime attestation recorded and verified: {attestation}")
    except (OSError, ValueError, json.JSONDecodeError) as error:
        parser.exit(1, f"native runtime attestation: {error}\n")


if __name__ == "__main__":
    main()
