#!/usr/bin/env python3
"""Validate external-domain attestations without converting them into local evidence."""

import argparse
import hashlib
import json
import math
import re
import shutil
import struct
import tempfile
import xml.etree.ElementTree as ET
from datetime import datetime, timezone
from pathlib import Path


KNOWN_DOMAINS = {
    "physical-radeon-780m",
    "wsl2-dozen",
    "windows-native-build",
    "windows-native-vulkan",
    "macos-native-build",
    "macos-moltenvk",
    "physical-imax-5616x4096",
    "viewer-native-window-input",
}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def text_contains(device, *needles):
    text = " ".join(
        str(device.get(field, ""))
        for field in ("name", "driver_name", "driver_info")
    ).casefold()
    return any(needle.casefold() in text for needle in needles)


def require_positive_integer(value, field):
    require(type(value) is int and value > 0,
            f"runtime device {field} must be a positive integer")


def png_dimensions(path):
    with path.open("rb") as stream:
        header = stream.read(24)
    require(len(header) == 24 and header[:8] == b"\x89PNG\r\n\x1a\n"
            and header[12:16] == b"IHDR",
            f"artifact is not a PNG with an IHDR header: {path.name}")
    return struct.unpack(">II", header[16:24])


def verify_document(data, location):
    require(data.get("schema_version") == 1, "schema_version must be 1")
    require(data.get("status") == "pass", "status must be pass")
    domains = data.get("domains")
    require(isinstance(domains, list) and domains, "domains must be a non-empty list")
    require(all(isinstance(domain, str) and domain for domain in domains),
            "every domain must be a non-empty string")
    require(len(domains) == len(set(domains)), "domains must be unique")
    unknown = set(domains) - KNOWN_DOMAINS
    require(not unknown, f"unknown attestation domains: {sorted(unknown)}")

    source_revision = data.get("source_revision")
    require(isinstance(source_revision, str)
            and re.fullmatch(r"[0-9a-f]{40,64}", source_revision) is not None,
            "source_revision must be a full lowercase hexadecimal revision")
    completed = data.get("completed_utc")
    require(isinstance(completed, str), "completed_utc is required")
    completed_at = datetime.fromisoformat(completed.replace("Z", "+00:00"))
    require(completed_at.utcoffset() == timezone.utc.utcoffset(completed_at),
            "completed_utc must carry an explicit UTC offset")

    platform = data.get("platform", {})
    device = data.get("device", {})
    claims = data.get("claims", {})
    require(claims.get("test_estate_passed") is True, "test estate must have passed")

    physical_domains = {
        "physical-radeon-780m",
        "wsl2-dozen",
        "windows-native-vulkan",
        "macos-moltenvk",
        "physical-imax-5616x4096",
    }
    if physical_domains.intersection(domains):
        require(device.get("kind") in {"integrated", "discrete"},
                "runtime domains require a physical Vulkan device")
        require(isinstance(device.get("name"), str) and device["name"].strip(),
                "runtime domains require a non-empty device name")
        require(isinstance(device.get("driver_name"), str) and device["driver_name"].strip(),
                "runtime domains require a non-empty driver name")
        require_positive_integer(device.get("vendor_id"), "vendor_id")
        require_positive_integer(device.get("device_id"), "device_id")
        require_positive_integer(device.get("api_version"), "api_version")
        require(type(device.get("driver_id")) is int and device["driver_id"] >= 0,
                "runtime device driver_id must be a non-negative integer")
        require_positive_integer(device.get("device_local_bytes"), "device_local_bytes")
        require_positive_integer(device.get("render_memory_bytes"), "render_memory_bytes")
        require(type(device.get("supports_fp64")) is bool,
                "runtime device supports_fp64 must be a boolean")
        require(claims.get("runtime_ready") is True,
                "runtime domains require a successful readiness probe")

    if "physical-radeon-780m" in domains:
        require("radeon 780m" in str(device.get("name", "")).casefold(),
                "Radeon 780M domain requires an exact device-name witness")
    if "wsl2-dozen" in domains:
        require(platform.get("wsl2") is True, "Dozen domain requires WSL2")
        require(text_contains(device, "dozen", "dzn"),
                "Dozen domain requires driver identity containing Dozen/dzn")
    if "windows-native-build" in domains or "windows-native-vulkan" in domains:
        require(str(platform.get("os", "")).casefold() == "windows",
                "Windows domains require a native Windows host")
    if "windows-native-vulkan" in domains:
        require(not text_contains(device, "dozen", "dzn"),
                "native Windows Vulkan may not be attested through Dozen")
    if "macos-native-build" in domains or "macos-moltenvk" in domains:
        require(str(platform.get("os", "")).casefold() == "macos",
                "macOS domains require a native macOS host")
    if "macos-moltenvk" in domains:
        require(text_contains(device, "moltenvk"),
                "MoltenVK domain requires MoltenVK driver identity")
    if "viewer-native-window-input" in domains:
        viewer = claims.get("viewer", {})
        require(viewer.get("window_created") is True, "viewer window was not created")
        require(viewer.get("keyboard_event_observed") is True,
                "no native keyboard event was observed")
        require(viewer.get("pointer_event_observed") is True,
                "no native pointer event was observed")
    if "physical-imax-5616x4096" in domains:
        render = claims.get("imax_render", {})
        require(render.get("width") == 5616 and render.get("height") == 4096,
                "IMAX domain requires an exact 5616x4096 render")
        require(render.get("memory_budget_mib") == 2048,
                "IMAX domain requires the 2048 MiB governed budget")
        require(type(render.get("wall_seconds")) in (int, float)
                and math.isfinite(render["wall_seconds"]) and render["wall_seconds"] > 0,
                "IMAX domain requires a positive measured wall time")

    artifacts = data.get("artifacts")
    require(isinstance(artifacts, dict) and artifacts, "artifacts must be a non-empty object")
    base = location.parent
    artifact_paths = []
    resolved_artifact_paths = set()
    for name, record in artifacts.items():
        require(isinstance(name, str) and name, "artifact names must be non-empty strings")
        require(isinstance(record, dict), f"artifact {name} must be an object")
        relative = record.get("path", name)
        require(isinstance(relative, str) and relative,
                f"artifact {name} path must be a non-empty string")
        path = Path(relative)
        require(not path.is_absolute() and ".." not in path.parts,
                f"artifact {name} path must stay inside the attestation bundle")
        path = base / path
        resolved_path = path.resolve()
        try:
            resolved_path.relative_to(base.resolve())
        except ValueError as error:
            raise ValueError(f"artifact {name} escapes the attestation bundle") from error
        require(resolved_path not in resolved_artifact_paths,
                f"artifact {name} aliases another artifact path")
        resolved_artifact_paths.add(resolved_path)
        require(path.is_file(), f"artifact {name} is missing: {path}")
        payload = path.read_bytes()
        recorded_bytes = record.get("bytes")
        require(type(recorded_bytes) is int and recorded_bytes > 0,
                f"artifact {name} byte count must be a positive integer")
        require(len(payload) == recorded_bytes, f"artifact {name} size mismatch")
        recorded_digest = record.get("sha256")
        require(isinstance(recorded_digest, str)
                and len(recorded_digest) == 64
                and all(character in "0123456789abcdef" for character in recorded_digest),
                f"artifact {name} SHA-256 must be 64 lowercase hexadecimal characters")
        digest = hashlib.sha256(payload).hexdigest()
        require(digest == recorded_digest, f"artifact {name} hash mismatch")
        artifact_paths.append(path)

    build_domains = {"windows-native-build", "macos-native-build"}
    if build_domains.intersection(domains):
        reports = [path for path in artifact_paths if path.suffix.casefold() == ".xml"]
        require(len(reports) == 1, "native build attestation requires exactly one JUnit XML report")
        report = inspect_junit(reports[0])
        require(claims.get("test_report") == report,
                "native build JUnit summary does not match the attested claim")
    if physical_domains.intersection(domains):
        pngs = [path for path in artifact_paths if path.suffix.casefold() == ".png"]
        transcripts = [
            path for path in artifact_paths if path.suffix.casefold() in {".log", ".txt"}
        ]
        require(pngs, "runtime attestation requires a hashed PNG output")
        require(transcripts, "runtime attestation requires a hashed execution transcript")
    if "viewer-native-window-input" in domains:
        transcripts = [
            path for path in artifact_paths if path.suffix.casefold() in {".log", ".txt"}
        ]
        require(transcripts, "native viewer attestation requires a hashed input transcript")
    if "physical-imax-5616x4096" in domains:
        pngs = [path for path in artifact_paths if path.suffix.casefold() == ".png"]
        require(pngs, "IMAX domain requires a hashed PNG artifact")
        require(any(png_dimensions(path) == (5616, 4096) for path in pngs),
                "no hashed PNG artifact has the attested 5616x4096 dimensions")


def verify_path(path):
    data = json.loads(path.read_text(encoding="utf-8"))
    verify_document(data, path)

def inspect_junit(path):
    try:
        root = ET.parse(path).getroot()
    except ET.ParseError as error:
        raise ValueError(f"build evidence is not valid JUnit XML: {error}") from error
    cases = list(root.iter("testcase"))
    require(cases, "build evidence contains no JUnit test cases")
    failures = sum(1 for case in cases if case.find("failure") is not None)
    errors = sum(1 for case in cases if case.find("error") is not None)
    skipped = sum(1 for case in cases if case.find("skipped") is not None)
    require(failures == 0 and errors == 0,
            f"build evidence contains {failures} failure(s) and {errors} error(s)")
    return {"cases": len(cases), "failures": failures, "errors": errors, "skipped": skipped}


def self_test():
    with tempfile.TemporaryDirectory(prefix="sirius-attestation-") as directory:
        root = Path(directory)
        artifact = root / "frame.png"
        artifact.write_bytes(
            b"\x89PNG\r\n\x1a\n"
            + struct.pack(">I", 13)
            + b"IHDR"
            + struct.pack(">II", 5616, 4096)
        )
        transcript = root / "transcript.log"
        transcript.write_text("self-test physical runtime transcript\n", encoding="utf-8")
        digest = hashlib.sha256(artifact.read_bytes()).hexdigest()
        transcript_digest = hashlib.sha256(transcript.read_bytes()).hexdigest()
        valid = {
            "schema_version": 1,
            "status": "pass",
            "domains": ["physical-radeon-780m", "physical-imax-5616x4096"],
            "source_revision": "0" * 40,
            "completed_utc": "2026-07-31T00:00:00+00:00",
            "platform": {"os": "linux", "wsl2": False},
            "device": {
                "name": "AMD Radeon 780M",
                "kind": "integrated",
                "driver_name": "RADV",
                "driver_info": "self-test",
                "vendor_id": 4098,
                "device_id": 5567,
                "api_version": 4206592,
                "driver_id": 3,
                "device_local_bytes": 8589934592,
                "render_memory_bytes": 2147483648,
                "supports_fp64": True,
            },
            "claims": {
                "test_estate_passed": True,
                "runtime_ready": True,
                "imax_render": {
                    "width": 5616,
                    "height": 4096,
                    "memory_budget_mib": 2048,
                    "wall_seconds": 1.0,
                },
            },
            "artifacts": {
                "frame.png": {
                    "path": "frame.png",
                    "bytes": artifact.stat().st_size,
                    "sha256": digest,
                },
                "transcript.log": {
                    "path": "transcript.log",
                    "bytes": transcript.stat().st_size,
                    "sha256": transcript_digest,
                },
            },
        }
        verify_document(valid, root / "attestation.json")

        corruptions = [
            ("software device", lambda doc: doc["device"].update(kind="software")),
            ("wrong dimensions", lambda doc: doc["claims"]["imax_render"].update(width=4096)),
            ("false test state", lambda doc: doc["claims"].update(test_estate_passed=False)),
            ("bad hash", lambda doc: doc["artifacts"]["frame.png"].update(sha256="0" * 64)),
            ("zero artifact size", lambda doc: doc["artifacts"]["frame.png"].update(bytes=0)),
            ("missing driver identity", lambda doc: doc["device"].update(driver_name="")),
            ("boolean vendor identity", lambda doc: doc["device"].update(vendor_id=True)),
            ("zero device memory", lambda doc: doc["device"].update(render_memory_bytes=0)),
            (
                "aliased artifact",
                lambda doc: doc["artifacts"].update(
                    alias={
                        "path": "frame.png",
                        "bytes": artifact.stat().st_size,
                        "sha256": digest,
                    }
                ),
            ),
            ("ambiguous source revision", lambda doc: doc.update(source_revision="main")),
            ("naive completion time", lambda doc: doc.update(completed_utc="2026-07-31T00:00:00")),
            ("non-UTC completion time",
             lambda doc: doc.update(completed_utc="2026-07-31T10:00:00+10:00")),
        ]
        for label, corrupt in corruptions:
            candidate = json.loads(json.dumps(valid))
            corrupt(candidate)
            try:
                verify_document(candidate, root / "attestation.json")
            except ValueError:
                continue
            raise ValueError(f"negative control accepted: {label}")

        artifact.write_bytes(
            b"\x89PNG\r\n\x1a\n"
            + struct.pack(">I", 13)
            + b"IHDR"
            + struct.pack(">II", 4096, 4096)
        )
        mismatched_image = json.loads(json.dumps(valid))
        payload = artifact.read_bytes()
        mismatched_image["artifacts"]["frame.png"].update(
            bytes=len(payload), sha256=hashlib.sha256(payload).hexdigest()
        )
        try:
            verify_document(mismatched_image, root / "attestation.json")
        except ValueError:
            pass
        else:
            raise ValueError("negative control accepted: PNG dimensions contradict claim")
        artifact.write_bytes(
            b"\x89PNG\r\n\x1a\n"
            + struct.pack(">I", 13)
            + b"IHDR"
            + struct.pack(">II", 5616, 4096)
        )

        domain_controls = [
            (
                "Dozen outside WSL2",
                ["wsl2-dozen"],
                {"os": "linux", "wsl2": False},
                {"name": "AMD Radeon 780M", "kind": "integrated",
                 "driver_name": "Mesa Dozen", "driver_info": "dzn",
                 "vendor_id": 4098, "device_id": 5567, "api_version": 4206592,
                 "driver_id": 3, "device_local_bytes": 8589934592,
                 "render_memory_bytes": 2147483648, "supports_fp64": True},
            ),
            (
                "Dozen called native Windows",
                ["windows-native-vulkan"],
                {"os": "windows", "wsl2": False},
                {"name": "AMD Radeon 780M", "kind": "integrated",
                 "driver_name": "Mesa Dozen", "driver_info": "dzn",
                 "vendor_id": 4098, "device_id": 5567, "api_version": 4206592,
                 "driver_id": 3, "device_local_bytes": 8589934592,
                 "render_memory_bytes": 2147483648, "supports_fp64": True},
            ),
            (
                "non-MoltenVK macOS runtime",
                ["macos-moltenvk"],
                {"os": "macos", "wsl2": False},
                {"name": "Apple GPU", "kind": "integrated",
                 "driver_name": "unknown", "driver_info": "unknown",
                 "vendor_id": 4203, "device_id": 1, "api_version": 4206592,
                 "driver_id": 14, "device_local_bytes": 8589934592,
                 "render_memory_bytes": 2147483648, "supports_fp64": False},
            ),
        ]
        for label, domains, platform, device in domain_controls:
            candidate = json.loads(json.dumps(valid))
            candidate["domains"] = domains
            candidate["platform"] = platform
            candidate["device"] = device
            try:
                verify_document(candidate, root / "attestation.json")
            except ValueError:
                continue
            raise ValueError(f"negative control accepted: {label}")

        junit = root / "tests.xml"
        junit.write_text(
            '<testsuite tests="1" failures="0"><testcase name="pass"/></testsuite>\n',
            encoding="utf-8",
        )
        require(inspect_junit(junit)["cases"] == 1, "valid JUnit control was not counted")
        junit.write_text(
            '<testsuite tests="1" failures="1"><testcase name="fail">'
            '<failure message="intentional"/></testcase></testsuite>\n',
            encoding="utf-8",
        )
        try:
            inspect_junit(junit)
        except ValueError:
            pass
        else:
            raise ValueError("negative control accepted: failing JUnit build evidence")


def record_build(args):
    require(args.domain in {"windows-native-build", "macos-native-build"},
            "record-build accepts only native build domains")
    require(args.source_revision, "source revision is required")
    require(args.artifact.is_file(), f"build evidence is missing: {args.artifact}")
    report = inspect_junit(args.artifact)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    bundled = args.output.parent / args.artifact.name
    if args.artifact.resolve() != bundled.resolve():
        shutil.copy2(args.artifact, bundled)
    payload = bundled.read_bytes()
    document = {
        "schema_version": 1,
        "status": "pass",
        "domains": [args.domain],
        "source_revision": args.source_revision,
        "completed_utc": datetime.now(timezone.utc).isoformat(),
        "platform": {"os": args.platform, "wsl2": False},
        "device": {},
        "claims": {
            "test_estate_passed": True,
            "runtime_ready": False,
            "test_report": report,
        },
        "artifacts": {
            bundled.name: {
                "path": bundled.name,
                "bytes": len(payload),
                "sha256": hashlib.sha256(payload).hexdigest(),
            }
        },
    }
    args.output.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")
    verify_document(document, args.output)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("attestation", nargs="?", type=Path)
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--record-build", action="store_true")
    parser.add_argument("--domain")
    parser.add_argument("--platform")
    parser.add_argument("--source-revision")
    parser.add_argument("--artifact", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    try:
        if args.self_test:
            self_test()
            print("attestation verifier rejected every false-evidence control")
        elif args.record_build:
            require(args.platform in {"windows", "macos"},
                    "record-build platform must be windows or macos")
            require(args.artifact is not None and args.output is not None,
                    "record-build requires --artifact and --output")
            record_build(args)
            print(f"build attestation recorded and verified: {args.output}")
        elif args.attestation is not None:
            verify_path(args.attestation)
            print(f"attestation verified: {args.attestation}")
        else:
            parser.error("provide an attestation path or --self-test")
    except (OSError, ValueError, json.JSONDecodeError) as error:
        parser.exit(1, f"attestation: {error}\n")


if __name__ == "__main__":
    main()
