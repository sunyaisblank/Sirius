#!/usr/bin/env python3
"""Reuse one verified WSL2 hardware qualification estate for viewer evidence.

The helper never creates an attestation and never runs tests, rendering, or a
viewer.  It verifies an already complete hardware record, proves that its exact
candidate/device/platform equal the live WSL2 viewer environment, and copies
only the qualification authority artifacts consumed by validate-viewer.sh.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
import shutil
import sys
from pathlib import Path

from attestation_preflight import (
    PROFILES,
    run_candidate_json,
    select_profile_device,
)


ROOT = Path(__file__).resolve().parents[1]
VERIFIER_PATH = ROOT / "scripts" / "verify-attestation.py"
PRODUCT_ARTIFACTS = {
    "qualification-product-operating_model": "model/operating_model.json",
    "qualification-product-starfield": "assets/Starfield.png",
    "qualification-product-trace_fp32comp_spv": "kernels/trace_fp32comp.spv",
    "qualification-product-trace_fp64_spv": "kernels/trace_fp64.spv",
    "qualification-product-trace_spv": "kernels/trace.spv",
    "qualification-product-viewer_rdsd003a_fragment": "shaders/RDSD003A.frag",
    "qualification-product-viewer_rdsd003a_vertex": "shaders/RDSD003A.vert",
}
FIXED_ARTIFACTS = {
    "hardware-tests.xml": "viewer-tests.xml",
    "alignment_receipt.json": "alignment_receipt.json",
    "ctest-inventory.json": "ctest-inventory.json",
    "mandatory_gate.json": "mandatory_gate.json",
    "qualification-sirius.bin": "qualification-sirius.bin",
    "qualification-gate-junit": "qualification-gate-junit",
    "qualification-gate-log": "qualification-gate-log",
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_verifier():
    spec = importlib.util.spec_from_file_location("sirius_attestation", VERIFIER_PATH)
    require(spec is not None and spec.loader is not None,
            "could not load the attestation verifier")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def classify_reusable_record(document: dict) -> dict[str, dict]:
    domains = document.get("domains")
    require(isinstance(domains, list)
            and "physical-radeon-780m" in domains
            and "viewer-native-window-input" not in domains,
            "viewer reuse requires an upstream physical-Radeon record")
    platform = document.get("platform", {})
    require(str(platform.get("os", "")).casefold() == "linux"
            and platform.get("wsl2") is True,
            "viewer reuse requires a WSL2 hardware record")
    require(document.get("claims", {}).get("preset") == "linux-ci"
            and document.get("claims", {}).get("runtime_ready") is True,
            "viewer reuse requires a ready linux-ci hardware record")
    device = document.get("device", {})
    require(device.get("kind") in {"integrated", "discrete"}
            and "radeon 780m" in str(device.get("name", "")).casefold(),
            "viewer reuse requires the exact physical Radeon 780M")
    identity = " ".join(
        str(device.get(field, "")) for field in ("driver_name", "driver_info")
    ).casefold()
    require("dozen" in identity or "dzn" in identity,
            "consolidated WSL2 reuse requires Dozen/dzn identity")

    artifacts = document.get("artifacts")
    require(isinstance(artifacts, dict), "upstream record has no artifact map")
    required = set(FIXED_ARTIFACTS) | set(PRODUCT_ARTIFACTS)
    missing = sorted(required - set(artifacts))
    require(not missing,
            "upstream hardware record lacks reusable authority artifacts: "
            + ", ".join(missing))
    return {name: artifacts[name] for name in required}


def resolved_artifact(record_path: Path, record: dict) -> Path:
    relative = record.get("path")
    require(isinstance(relative, str) and relative,
            "upstream artifact has no relative path")
    path = (record_path.parent / relative).resolve()
    require(path.is_relative_to(record_path.parent.resolve()) and path.is_file(),
            "upstream artifact escapes or is missing from its record bundle")
    require(type(record.get("bytes")) is int and path.stat().st_size == record["bytes"],
            "upstream reusable artifact byte count changed")
    require(sha256(path) == record.get("sha256"),
            "upstream reusable artifact digest changed")
    return path


def find_device_inventory(document: dict, record_path: Path) -> Path:
    matches = []
    for record in document["artifacts"].values():
        path = resolved_artifact(record_path, record)
        if path.suffix.casefold() != ".json":
            continue
        try:
            value = json.loads(path.read_text(encoding="utf-8"))
        except (UnicodeDecodeError, json.JSONDecodeError):
            continue
        if isinstance(value, dict) and "backends" in value:
            matches.append((path, value))
    require(len(matches) == 1,
            "upstream hardware record must contain one device inventory")
    path, inventory = matches[0]
    vulkan = inventory.get("backends", {}).get("vulkan", {})
    require(vulkan.get("selected_device_index") == document["device"].get("index")
            and sum(device == document["device"]
                    for device in vulkan.get("devices", [])) == 1,
            "upstream inventory does not bind its attested device")
    return path


def verify_live_identity(
    document: dict,
    expected_revision: str,
    inventory: dict,
    selected: dict,
) -> None:
    require(document.get("source_revision") == expected_revision,
            "upstream qualification record names a different source revision")
    require(document.get("device") == selected,
            "live viewer device differs from the upstream hardware device")
    upstream_platform = document.get("platform", {})
    require(str(upstream_platform.get("os", "")).casefold()
            == str(inventory.get("platform", "")).casefold()
            and upstream_platform.get("wsl2") is inventory.get("wsl2"),
            "live viewer platform differs from the upstream hardware platform")


def find_execution_transcript(document: dict, record_path: Path) -> Path:
    revision = document["source_revision"]
    expected_lines = {
        f"== source revision: {revision}",
        "== [1/6] configure + build (linux-ci)",
        "== readiness: evidence_generation=true vulkan=true "
        f"selected_device_index={document['device']['index']}",
    }
    matches = []
    for record in document["artifacts"].values():
        path = resolved_artifact(record_path, record)
        if path.suffix.casefold() not in {".log", ".txt"}:
            continue
        lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
        has_test_summary = any(
            line.startswith("100% tests passed, 0 tests failed out of ")
            for line in lines
        )
        if expected_lines <= set(lines) and has_test_summary:
            matches.append(path)
    require(len(matches) == 1,
            "upstream hardware record must contain one complete execution transcript")
    return matches[0]


def verify_current_volume(
    candidate: Path,
    document: dict,
    reusable: dict[str, dict],
    record_path: Path,
) -> None:
    claim = document.get("claims", {}).get("qualification_executable", {})
    require(candidate.is_file() and candidate.stat().st_size == claim.get("bytes")
            and sha256(candidate) == claim.get("sha256"),
            "live viewer candidate differs from the upstream qualification candidate")
    candidate_record = reusable["qualification-sirius.bin"]
    require(candidate_record.get("bytes") == claim.get("bytes")
            and candidate_record.get("sha256") == claim.get("sha256"),
            "upstream record does not bind its qualification candidate")
    resource_root = candidate.parent / "resources"
    comparisons = {
        "alignment_receipt.json": "model/alignment_receipt.json",
        "mandatory_gate.json": "model/mandatory_gate.json",
        **PRODUCT_ARTIFACTS,
    }
    for artifact_name, relative in comparisons.items():
        current = resource_root / relative
        source = resolved_artifact(record_path, reusable[artifact_name])
        require(current.is_file() and current.stat().st_size == source.stat().st_size
                and sha256(current) == sha256(source),
                f"live viewer resource differs from upstream authority: {relative}")


def prepare_reuse(
    record_path: Path,
    candidate: Path,
    output_dir: Path,
    expected_revision: str,
    device_pattern: str | None,
) -> dict:
    record_path = record_path.resolve()
    candidate = candidate.resolve()
    output_dir = output_dir.resolve()
    document = json.loads(record_path.read_text(encoding="utf-8"))
    verifier = load_verifier()
    verifier.verify_document(document, record_path)
    reusable = classify_reusable_record(document)
    verify_current_volume(candidate, document, reusable, record_path)

    environment = os.environ.copy()
    inventory = run_candidate_json(
        candidate, ("--json", "info", "system"), environment
    )
    tools = {
        command: shutil.which(command) is not None
        for command in ("cmake", "ctest", "git", "timeout", "xdotool")
    }
    selected = select_profile_device(
        inventory,
        PROFILES["wsl2-radeon-viewer"],
        host_platform=sys.platform,
        device_pattern=device_pattern,
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
        PROFILES["wsl2-radeon-viewer"],
        host_platform=sys.platform,
        device_pattern=device_pattern,
        tools=tools,
        display=environment.get("DISPLAY", ""),
        dxg_present=Path("/dev/dxg").exists(),
    )
    require(inventory.get("backends", {}).get("vulkan", {}).get(
        "selected_device_index") == selected["index"],
        "live viewer candidate did not preserve the selected device")
    verify_live_identity(document, expected_revision, inventory, selected)

    output_dir.mkdir(parents=True, exist_ok=True)
    for source_name, destination_name in FIXED_ARTIFACTS.items():
        destination = output_dir / destination_name
        require(not destination.exists(),
                f"viewer reuse destination already exists: {destination.name}")
        shutil.copy2(resolved_artifact(record_path, reusable[source_name]), destination)
    for artifact_name in PRODUCT_ARTIFACTS:
        destination = output_dir / artifact_name
        require(not destination.exists(),
                f"viewer reuse destination already exists: {destination.name}")
        shutil.copy2(resolved_artifact(record_path, reusable[artifact_name]), destination)
    inventory_source = find_device_inventory(document, record_path)
    shutil.copy2(inventory_source, output_dir / "reused-device-inventory.json")
    transcript_source = find_execution_transcript(document, record_path)
    shutil.copy2(
        transcript_source, output_dir / "reused-qualification-transcript.log"
    )
    return {
        "schema_version": 1,
        "kind": "sirius-qualification-evidence-reuse",
        "source_revision": expected_revision,
        "upstream_domains": document["domains"],
        "selected_device_index": selected["index"],
        "copied_authority_artifacts": len(FIXED_ARTIFACTS) + len(PRODUCT_ARTIFACTS) + 2,
        "promotable": False,
    }


def self_test() -> None:
    device = {
        "index": 0,
        "name": "Microsoft Direct3D12 (AMD Radeon 780M Graphics)",
        "kind": "integrated",
        "driver_name": "Dozen",
        "driver_info": "Mesa dzn",
    }
    artifacts = {
        name: {"path": name, "bytes": 1, "sha256": "0" * 64}
        for name in set(FIXED_ARTIFACTS) | set(PRODUCT_ARTIFACTS)
    }
    valid = {
        "source_revision": "a" * 40,
        "domains": ["physical-radeon-780m", "physical-imax-5616x4096", "wsl2-dozen"],
        "platform": {"os": "linux", "wsl2": True},
        "device": device,
        "claims": {"preset": "linux-ci", "runtime_ready": True},
        "artifacts": artifacts,
    }
    require(len(classify_reusable_record(valid)) == len(artifacts),
            "qualification reuse positive control failed")
    live_inventory = {"platform": "Linux", "wsl2": True}
    verify_live_identity(valid, "a" * 40, live_inventory, device)
    controls = (
        ("viewer circularity", {**valid, "domains": ["viewer-native-window-input"]}),
        ("native host", {**valid, "platform": {"os": "linux", "wsl2": False}}),
        ("wrong driver", {**valid, "device": {**device, "driver_name": "RADV",
                                               "driver_info": "native"}}),
        ("missing gate", {**valid, "artifacts": {
            name: value for name, value in artifacts.items()
            if name != "mandatory_gate.json"
        }}),
    )
    for label, document in controls:
        try:
            classify_reusable_record(document)
        except ValueError:
            continue
        raise ValueError(f"qualification reuse accepted false control: {label}")

    identity_controls = (
        ("stale revision", valid, "b" * 40, live_inventory, device),
        ("cross-device", valid, "a" * 40, live_inventory,
         {**device, "device_id": 0x1901}),
        ("cross-platform", valid, "a" * 40,
         {"platform": "Linux", "wsl2": False}, device),
    )
    for label, document, revision, inventory, selected in identity_controls:
        try:
            verify_live_identity(document, revision, inventory, selected)
        except ValueError:
            continue
        raise ValueError(f"qualification reuse accepted false identity: {label}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--record", type=Path)
    parser.add_argument("--candidate", type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--expected-revision")
    parser.add_argument("--device-pattern")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    try:
        if args.self_test:
            self_test()
            print("qualification reuse rejected stale, circular, and cross-device controls")
            return
        require(args.record is not None and args.candidate is not None
                and args.output_dir is not None and args.expected_revision is not None,
                "--record, --candidate, --output-dir, and --expected-revision are required")
        result = prepare_reuse(
            args.record, args.candidate, args.output_dir,
            args.expected_revision, args.device_pattern,
        )
        print(json.dumps(result, sort_keys=True))
    except (OSError, ValueError, json.JSONDecodeError) as error:
        parser.exit(1, f"qualification evidence reuse: {error}\n")


if __name__ == "__main__":
    main()
