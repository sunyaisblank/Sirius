#!/usr/bin/env python3
"""Run, record, and verify the artifact-bound Sirius Mandatory build gate."""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import os
import re
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET
from pathlib import Path


SCHEMA_VERSION = 1
TEST_FLOOR = 700
NATIVE_BUILD_GATE_KIND = "sirius-native-build-gate"
NATIVE_BUILD_EVIDENCE_TESTS = (
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
TESTED_ARTIFACTS = {
    "sirius",
    "sirius_app_tests",
    "sirius_backend_tests",
    "sirius_base_tests",
    "sirius_core_tests",
    "sirius_oracle_tests",
    "sirius_render_tests",
}
BASE_PRODUCTS = {
    "alignment_receipt",
    "operating_model",
    "sirius",
    "starfield",
}
RELEASE_PRODUCTS = BASE_PRODUCTS | {
    "trace_fp32comp_spv",
    "trace_fp64_spv",
    "trace_spv",
    "viewer_rdsd003a_fragment",
    "viewer_rdsd003a_vertex",
}
INSTALLED_PRODUCTS = {
    "alignment_receipt": "share/sirius/model/alignment_receipt.json",
    "operating_model": "share/sirius/model/operating_model.json",
    "starfield": "share/sirius/assets/Starfield.png",
    "trace_fp32comp_spv": "share/sirius/kernels/trace_fp32comp.spv",
    "trace_fp64_spv": "share/sirius/kernels/trace_fp64.spv",
    "trace_spv": "share/sirius/kernels/trace.spv",
    "viewer_rdsd003a_fragment": "share/sirius/shaders/RDSD003A.frag",
    "viewer_rdsd003a_vertex": "share/sirius/shaders/RDSD003A.vert",
}
LOGICAL_NAME = re.compile(r"^[a-z][a-z0-9_]*$")
REVISION = re.compile(r"^[0-9a-f]{40}$")


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def sha256_file(path: Path) -> str:
    return sha256_bytes(path.read_bytes())


def canonical_json(document: object) -> bytes:
    return json.dumps(document, sort_keys=True, separators=(",", ":")).encode("utf-8")


def parse_bool(value: str) -> bool:
    if value == "true":
        return True
    if value == "false":
        return False
    raise ValueError(f"expected true or false, got {value!r}")


def parse_artifacts(values: list[str], option: str) -> dict[str, Path]:
    artifacts: dict[str, Path] = {}
    for value in values:
        name, separator, raw_path = value.partition("=")
        require(separator == "=" and LOGICAL_NAME.fullmatch(name) is not None,
                f"{option} requires logical_name=path, got {value!r}")
        require(name not in artifacts, f"{option} repeats {name}")
        path = Path(raw_path).resolve()
        require(path.is_file(), f"{option} artifact is missing: {path}")
        require(path.stat().st_size > 0, f"{option} artifact is empty: {path}")
        artifacts[name] = path
    return artifacts


def locate_under_roots(path: Path, source_root: Path, build_dir: Path) -> tuple[str, str]:
    path = path.resolve()
    for name, unresolved_root in (("build", build_dir), ("source", source_root)):
        root = unresolved_root.resolve()
        try:
            relative = path.relative_to(root)
        except ValueError:
            continue
        require(relative.parts and ".." not in relative.parts,
                f"artifact path escapes the {name} root: {path}")
        return name, relative.as_posix()
    raise ValueError(f"artifact is outside the source and build roots: {path}")


def artifact_records(artifacts: dict[str, Path], source_root: Path,
                     build_dir: Path) -> dict[str, dict[str, object]]:
    records = {}
    for name, path in sorted(artifacts.items()):
        root, relative = locate_under_roots(path, source_root, build_dir)
        records[name] = {
            "root": root,
            "path": relative,
            "bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        }
    return records


def resolve_record(record: dict[str, object], source_root: Path, build_dir: Path) -> Path:
    require(set(record) == {"root", "path", "bytes", "sha256"},
            "artifact receipt has an unsupported shape")
    root_name = record["root"]
    relative = record["path"]
    byte_count = record["bytes"]
    digest = record["sha256"]
    require(root_name in {"source", "build"} and isinstance(relative, str),
            "artifact receipt has an invalid root or path")
    require(type(byte_count) is int and byte_count > 0,
            "artifact receipt has an invalid byte count")
    require(isinstance(digest, str) and re.fullmatch(r"[0-9a-f]{64}", digest) is not None,
            "artifact receipt has an invalid SHA-256")
    relative_path = Path(relative)
    require(not relative_path.is_absolute() and relative_path.parts and
            ".." not in relative_path.parts,
            "artifact receipt path escapes its declared root")
    root = (source_root if root_name == "source" else build_dir).resolve()
    path = (root / relative_path).resolve()
    try:
        path.relative_to(root)
    except ValueError as error:
        raise ValueError("artifact receipt resolves outside its declared root") from error
    require(path.is_file(), f"recorded build-gate artifact is missing: {path}")
    require(path.stat().st_size == byte_count,
            f"recorded build-gate artifact size changed: {path}")
    require(sha256_file(path) == digest,
            f"recorded build-gate artifact hash changed: {path}")
    return path


def inspect_inventory(document: object) -> tuple[set[str], str]:
    require(isinstance(document, dict) and document.get("kind") == "ctestInfo" and
            document.get("version", {}).get("major") == 1,
            "CTest inventory has an unsupported schema")
    tests = document.get("tests")
    require(isinstance(tests, list) and len(tests) >= TEST_FLOOR,
            f"CTest inventory must contain at least {TEST_FLOOR} tests")
    names: list[str] = []
    for test in tests:
        require(isinstance(test, dict) and isinstance(test.get("name"), str) and
                test["name"], "CTest inventory contains a malformed test")
        names.append(test["name"])
        properties = test.get("properties")
        require(isinstance(properties, list),
                f"CTest test {test['name']} has no property inventory")
        labels: list[str] = []
        for prop in properties:
            if isinstance(prop, dict) and prop.get("name") == "LABELS":
                value = prop.get("value")
                if isinstance(value, list):
                    labels = [item for item in value if isinstance(item, str)]
        require("Mandatory" in labels,
                f"CTest registration escaped the Mandatory gate: {test['name']}")
    require(len(names) == len(set(names)), "CTest inventory contains duplicate test identities")
    return set(names), sha256_bytes(canonical_json(document))


def read_inventory(ctest: Path, build_dir: Path, config: str) -> tuple[dict, set[str], str]:
    command = [str(ctest), "--test-dir", str(build_dir), "--show-only=json-v1"]
    if config:
        command.extend(("-C", config))
    try:
        result = subprocess.run(command, check=True, capture_output=True, text=True)
        document = json.loads(result.stdout)
    except (OSError, subprocess.CalledProcessError, json.JSONDecodeError) as error:
        raise ValueError(f"could not inventory live CTest registration: {error}") from error
    names, digest = inspect_inventory(document)
    return document, names, digest


def inspect_junit(path: Path) -> tuple[set[str], dict[str, int]]:
    try:
        root = ET.parse(path).getroot()
    except (OSError, ET.ParseError) as error:
        raise ValueError(f"Mandatory gate JUnit is invalid: {error}") from error
    cases = list(root.iter("testcase"))
    require(cases, "Mandatory gate JUnit contains no test cases")
    names = [case.get("name") for case in cases]
    require(all(isinstance(name, str) and name for name in names),
            "Mandatory gate JUnit contains an unnamed test case")
    require(len(names) == len(set(names)),
            "Mandatory gate JUnit contains duplicate test identities")
    failures = sum(case.find("failure") is not None for case in cases)
    errors = sum(case.find("error") is not None for case in cases)
    skipped = sum(case.find("skipped") is not None for case in cases)
    require(failures == 0 and errors == 0 and skipped == 0,
            "Mandatory gate JUnit is not a zero-failure, zero-error, zero-skip estate")
    return set(names), {
        "executed": len(cases),
        "failures": failures,
        "errors": errors,
        "skipped": skipped,
    }


def git_identity(source_root: Path) -> tuple[str, bool]:
    try:
        revision = subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=source_root, check=True,
            capture_output=True, text=True).stdout.strip()
        status = subprocess.run(
            ["git", "status", "--porcelain", "--untracked-files=normal"],
            cwd=source_root, check=True, capture_output=True, text=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError) as error:
        raise ValueError(f"could not establish build-gate source identity: {error}") from error
    require(REVISION.fullmatch(revision) is not None, "Git returned an invalid source revision")
    return revision, not status


def validate_document(document: object) -> dict:
    require(isinstance(document, dict) and set(document) == {
        "schema_version", "status", "alignment_mode", "source", "ctest",
        "inputs", "tested_artifacts", "product_artifacts",
    }, "Mandatory gate receipt has an unsupported shape")
    require(document["schema_version"] == SCHEMA_VERSION and document["status"] == "passed",
            "Mandatory gate receipt is not a passing schema-v1 record")
    mode = document["alignment_mode"]
    require(mode in {"development", "qualification", "release"},
            "Mandatory gate receipt has an invalid mode")
    source = document["source"]
    require(isinstance(source, dict) and set(source) == {"revision", "clean"} and
            isinstance(source["clean"], bool) and
            isinstance(source["revision"], str) and
            REVISION.fullmatch(source["revision"]) is not None,
            "Mandatory gate receipt has an invalid source identity")
    if mode in {"qualification", "release"}:
        require(source["clean"], f"{mode} Mandatory gate receipt names a dirty source tree")

    ctest = document["ctest"]
    require(isinstance(ctest, dict) and set(ctest) == {
        "inventory_sha256", "junit", "log", "registered", "executed",
        "failures", "errors", "skipped",
    }, "Mandatory gate receipt has an invalid CTest summary")
    require(type(ctest["registered"]) is int and ctest["registered"] >= TEST_FLOOR and
            ctest["executed"] == ctest["registered"] and ctest["failures"] == 0 and
            ctest["errors"] == 0 and ctest["skipped"] == 0,
            "Mandatory gate receipt is not a complete zero-skip estate")
    require(isinstance(ctest["inventory_sha256"], str) and
            re.fullmatch(r"[0-9a-f]{64}", ctest["inventory_sha256"]) is not None,
            "Mandatory gate receipt has an invalid inventory digest")

    inputs = document["inputs"]
    require(isinstance(inputs, dict) and set(inputs) == {
        "operating_model_sha256", "alignment_receipt_sha256",
    }, "Mandatory gate receipt has an invalid input summary")
    for digest in inputs.values():
        require(isinstance(digest, str) and re.fullmatch(r"[0-9a-f]{64}", digest) is not None,
                "Mandatory gate receipt has an invalid input digest")

    tested = document["tested_artifacts"]
    products = document["product_artifacts"]
    require(isinstance(tested, dict) and set(tested) == TESTED_ARTIFACTS,
            "Mandatory gate receipt does not bind every test executable")
    require(isinstance(products, dict) and BASE_PRODUCTS <= set(products) <= RELEASE_PRODUCTS,
            "Mandatory gate receipt has an invalid product artifact set")
    if mode in {"qualification", "release"}:
        require(set(products) == RELEASE_PRODUCTS,
                f"{mode} Mandatory gate receipt does not bind the complete product")
    return document


def validate_native_build_document(document: object) -> dict:
    """Validate the non-promoting gate for a native compilation claim.

    This gate deliberately proves less than the Mandatory release gate: it
    proves the complete strict product/test topology and the fixed non-render
    authority estate, but makes no runtime, image, or full-estate claim.
    """
    require(isinstance(document, dict) and set(document) == {
        "schema_version", "kind", "status", "alignment_mode", "source", "ctest",
        "inputs", "tested_artifacts", "product_artifacts",
    }, "native build gate receipt has an unsupported shape")
    require(document["schema_version"] == SCHEMA_VERSION and
            document["kind"] == NATIVE_BUILD_GATE_KIND and
            document["status"] == "passed",
            "native build gate receipt is not a passing schema-v1 record")
    require(document["alignment_mode"] == "qualification",
            "native build evidence must come from qualification mode")
    source = document["source"]
    require(isinstance(source, dict) and set(source) == {"revision", "clean"} and
            source.get("clean") is True and
            isinstance(source.get("revision"), str) and
            REVISION.fullmatch(source["revision"]) is not None,
            "native build gate does not name a clean exact source revision")

    ctest = document["ctest"]
    require(isinstance(ctest, dict) and set(ctest) == {
        "inventory_sha256", "junit", "log", "registered", "selected", "executed",
        "failures", "errors", "skipped", "selection",
    }, "native build gate has an invalid CTest summary")
    selection = ctest.get("selection")
    require(selection == sorted(NATIVE_BUILD_EVIDENCE_TESTS),
            "native build gate does not name the exact non-render authority estate")
    require(type(ctest.get("registered")) is int and ctest["registered"] >= TEST_FLOOR and
            ctest.get("selected") == ctest.get("executed") == len(selection) and
            ctest.get("failures") == ctest.get("errors") == ctest.get("skipped") == 0,
            "native build gate is not a complete zero-skip authority selection")
    require(isinstance(ctest.get("inventory_sha256"), str) and
            re.fullmatch(r"[0-9a-f]{64}", ctest["inventory_sha256"]) is not None,
            "native build gate has an invalid inventory digest")
    for key in ("junit", "log"):
        record = ctest.get(key)
        require(isinstance(record, dict) and
                set(record) == {"root", "path", "bytes", "sha256"},
                f"native build gate has an invalid {key} artifact record")

    inputs = document["inputs"]
    require(isinstance(inputs, dict) and set(inputs) == {
        "operating_model_sha256", "alignment_receipt_sha256",
    }, "native build gate has an invalid input summary")
    for digest in inputs.values():
        require(isinstance(digest, str) and re.fullmatch(r"[0-9a-f]{64}", digest),
                "native build gate has an invalid input digest")

    tested = document["tested_artifacts"]
    products = document["product_artifacts"]
    require(isinstance(tested, dict) and set(tested) == TESTED_ARTIFACTS,
            "native build gate does not bind every compiled test executable")
    require(isinstance(products, dict) and set(products) == RELEASE_PRODUCTS,
            "native build gate does not bind the complete strict product")
    return document


def verify_recorded_files(document: dict, source_root: Path, build_dir: Path) -> None:
    for collection_name in ("tested_artifacts", "product_artifacts"):
        collection = document[collection_name]
        for record in collection.values():
            require(isinstance(record, dict), "artifact receipt entry is not an object")
            resolve_record(record, source_root, build_dir)
    products = document["product_artifacts"]
    require(document["inputs"]["operating_model_sha256"] ==
            products["operating_model"]["sha256"],
            "Mandatory gate operating-model input is not the product model")
    require(document["inputs"]["alignment_receipt_sha256"] ==
            products["alignment_receipt"]["sha256"],
            "Mandatory gate alignment input is not the product receipt")
    for key in ("junit", "log"):
        record = document["ctest"][key]
        require(isinstance(record, dict), f"Mandatory gate {key} record is invalid")
        resolve_record(record, source_root, build_dir)


def installed_product_path(root: Path, name: str, executable_suffix: str) -> Path:
    if name == "sirius":
        return root / "bin" / f"sirius{executable_suffix}"
    return root / INSTALLED_PRODUCTS[name]


def verify_installed(document: dict, installed_root: Path, stamp: Path,
                     executable_suffix: str) -> None:
    for name, record in document["product_artifacts"].items():
        path = installed_product_path(installed_root, name, executable_suffix)
        require(path.is_file(), f"installed product is missing: {path}")
        require(path.stat().st_size == record["bytes"],
                f"installed product size differs from the tested artifact: {path}")
        require(sha256_file(path) == record["sha256"],
                f"installed product differs from the tested artifact: {path}")
    installed_stamp = installed_root / "share/sirius/model/mandatory_gate.json"
    require(installed_stamp.is_file() and installed_stamp.read_bytes() == stamp.read_bytes(),
            "installed Mandatory gate receipt differs from the verified build receipt")


def atomic_write_json(path: Path, document: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(document, indent=2, sort_keys=True) + "\n"
    descriptor, temporary = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as stream:
            stream.write(payload)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def run_gate(args: argparse.Namespace) -> None:
    source_root = args.source_root.resolve()
    build_dir = args.build_dir.resolve()
    stamp = args.stamp.resolve()
    junit = stamp.with_name("mandatory_gate_junit.xml")
    log = stamp.with_name("mandatory_gate_ctest.log")
    stamp.parent.mkdir(parents=True, exist_ok=True)
    for stale in (stamp, junit, log):
        stale.unlink(missing_ok=True)

    configured_clean = parse_bool(args.source_tree_clean)
    revision, clean = git_identity(source_root)
    require(revision == args.source_revision and clean == configured_clean,
            "source identity changed after Sirius configuration")
    if args.alignment_mode in {"qualification", "release"}:
        require(clean, f"{args.alignment_mode} Mandatory gate requires a clean source tree")

    _, registered_names, inventory_digest = read_inventory(args.ctest, build_dir, args.config)
    command = [
        str(args.ctest), "--test-dir", str(build_dir), "--output-on-failure",
        "--no-tests=error", "-L", "Mandatory", "--output-junit", str(junit),
    ]
    if args.config:
        command.extend(("-C", args.config))
    try:
        result = subprocess.run(command, capture_output=True, text=True)
    except OSError as error:
        raise ValueError(f"could not execute the Mandatory CTest gate: {error}") from error
    log.write_text(
        f"command: {' '.join(command)}\nreturn_code: {result.returncode}\n"
        f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}",
        encoding="utf-8",
    )
    require(result.returncode == 0, "Mandatory CTest gate failed; no receipt was written")
    executed_names, summary = inspect_junit(junit)
    require(executed_names == registered_names,
            "Mandatory JUnit identities do not equal live CTest registration")

    tested = parse_artifacts(args.tested_artifact, "--tested-artifact")
    products = parse_artifacts(args.product_artifact, "--product-artifact")
    require(set(tested) == TESTED_ARTIFACTS,
            "Mandatory gate invocation omitted a test executable")
    require(BASE_PRODUCTS <= set(products) <= RELEASE_PRODUCTS,
            "Mandatory gate invocation has an invalid product set")
    if args.alignment_mode in {"qualification", "release"}:
        require(set(products) == RELEASE_PRODUCTS,
                f"{args.alignment_mode} Mandatory gate invocation omitted a product artifact")

    tested_records = artifact_records(tested, source_root, build_dir)
    product_records = artifact_records(products, source_root, build_dir)
    junit_record = artifact_records({"junit": junit}, source_root, build_dir)["junit"]
    log_record = artifact_records({"log": log}, source_root, build_dir)["log"]
    document = {
        "schema_version": SCHEMA_VERSION,
        "status": "passed",
        "alignment_mode": args.alignment_mode,
        "source": {"revision": revision, "clean": clean},
        "ctest": {
            "inventory_sha256": inventory_digest,
            "junit": junit_record,
            "log": log_record,
            "registered": len(registered_names),
            **summary,
        },
        "inputs": {
            "operating_model_sha256": product_records["operating_model"]["sha256"],
            "alignment_receipt_sha256": product_records["alignment_receipt"]["sha256"],
        },
        "tested_artifacts": tested_records,
        "product_artifacts": product_records,
    }
    validate_document(document)
    atomic_write_json(stamp, document)
    print(f"Mandatory gate receipt bound {len(registered_names)} zero-skip tests to "
          f"{len(products)} product artifacts")


def run_native_build_gate(args: argparse.Namespace) -> None:
    source_root = args.source_root.resolve()
    build_dir = args.build_dir.resolve()
    stamp = args.stamp.resolve()
    junit = stamp.with_name("native_build_gate_junit.xml")
    log = stamp.with_name("native_build_gate_ctest.log")
    stamp.parent.mkdir(parents=True, exist_ok=True)
    for stale in (stamp, junit, log):
        stale.unlink(missing_ok=True)

    require(args.alignment_mode == "qualification",
            "native build gate requires qualification mode")
    configured_clean = parse_bool(args.source_tree_clean)
    revision, clean = git_identity(source_root)
    require(revision == args.source_revision and clean == configured_clean and clean,
            "native build source identity is dirty, changed, or not exact")

    _, registered_names, inventory_digest = read_inventory(args.ctest, build_dir, args.config)
    selected_names = set(NATIVE_BUILD_EVIDENCE_TESTS)
    require(selected_names <= registered_names,
            "native build authority selection is absent from live CTest registration")
    expression = "^(" + "|".join(re.escape(name) for name in
                                     NATIVE_BUILD_EVIDENCE_TESTS) + ")$"
    command = [
        str(args.ctest), "--test-dir", str(build_dir), "--output-on-failure",
        "--no-tests=error", "-R", expression, "--output-junit", str(junit),
    ]
    if args.config:
        command.extend(("-C", args.config))
    try:
        result = subprocess.run(command, capture_output=True, text=True)
    except OSError as error:
        raise ValueError(f"could not execute the native build CTest gate: {error}") from error
    log.write_text(
        f"command: {' '.join(command)}\nreturn_code: {result.returncode}\n"
        f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}",
        encoding="utf-8",
    )
    require(result.returncode == 0,
            "native build authority gate failed; no receipt was written")
    executed_names, summary = inspect_junit(junit)
    require(executed_names == selected_names,
            "native build JUnit differs from the exact non-render authority selection")

    tested = parse_artifacts(args.tested_artifact, "--tested-artifact")
    products = parse_artifacts(args.product_artifact, "--product-artifact")
    require(set(tested) == TESTED_ARTIFACTS,
            "native build gate invocation omitted a compiled test executable")
    require(set(products) == RELEASE_PRODUCTS,
            "native build gate invocation omitted a strict product artifact")
    tested_records = artifact_records(tested, source_root, build_dir)
    product_records = artifact_records(products, source_root, build_dir)
    document = {
        "schema_version": SCHEMA_VERSION,
        "kind": NATIVE_BUILD_GATE_KIND,
        "status": "passed",
        "alignment_mode": "qualification",
        "source": {"revision": revision, "clean": clean},
        "ctest": {
            "inventory_sha256": inventory_digest,
            "junit": artifact_records({"junit": junit}, source_root, build_dir)["junit"],
            "log": artifact_records({"log": log}, source_root, build_dir)["log"],
            "registered": len(registered_names),
            "selected": len(selected_names),
            "selection": sorted(selected_names),
            **summary,
        },
        "inputs": {
            "operating_model_sha256": product_records["operating_model"]["sha256"],
            "alignment_receipt_sha256": product_records["alignment_receipt"]["sha256"],
        },
        "tested_artifacts": tested_records,
        "product_artifacts": product_records,
    }
    validate_native_build_document(document)
    atomic_write_json(stamp, document)
    print(f"native build gate bound {len(tested_records)} compiled executables and "
          f"{len(product_records)} products to {len(selected_names)} non-render controls")


def check_gate(args: argparse.Namespace) -> None:
    source_root = args.source_root.resolve()
    build_dir = args.build_dir.resolve()
    stamp = args.stamp.resolve()
    require(stamp.is_file(), "Mandatory gate receipt is missing; run the complete build gate")
    try:
        document = validate_document(json.loads(stamp.read_text(encoding="utf-8")))
    except json.JSONDecodeError as error:
        raise ValueError(f"Mandatory gate receipt is invalid JSON: {error}") from error
    require(document["alignment_mode"] == args.alignment_mode,
            "Mandatory gate receipt belongs to a different alignment mode")
    revision, clean = git_identity(source_root)
    require(document["source"] == {"revision": revision, "clean": clean},
            "Mandatory gate receipt is stale for the current source identity")
    if args.source_revision is not None:
        require(revision == args.source_revision,
                "Mandatory gate receipt differs from the configured source revision")
    _, inventory_names, inventory_digest = read_inventory(args.ctest, build_dir, args.config)
    require(document["ctest"]["inventory_sha256"] == inventory_digest and
            document["ctest"]["registered"] == len(inventory_names),
            "Mandatory gate receipt is stale for live CTest registration")
    verify_recorded_files(document, source_root, build_dir)
    junit_path = resolve_record(document["ctest"]["junit"], source_root, build_dir)
    junit_names, summary = inspect_junit(junit_path)
    require(junit_names == inventory_names and
            all(document["ctest"][key] == value for key, value in summary.items()),
            "Mandatory gate JUnit no longer proves the registered zero-skip estate")
    if args.installed_root is not None:
        verify_installed(document, args.installed_root.resolve(), stamp,
                         args.executable_suffix)
    print(f"Mandatory gate receipt verifies {len(inventory_names)} tests and the exact product")


def check_native_build_gate(args: argparse.Namespace) -> None:
    source_root = args.source_root.resolve()
    build_dir = args.build_dir.resolve()
    stamp = args.stamp.resolve()
    require(stamp.is_file(), "native build gate receipt is missing")
    try:
        document = validate_native_build_document(
            json.loads(stamp.read_text(encoding="utf-8"))
        )
    except json.JSONDecodeError as error:
        raise ValueError(f"native build gate receipt is invalid JSON: {error}") from error
    revision, clean = git_identity(source_root)
    require(document["source"] == {"revision": revision, "clean": clean},
            "native build gate receipt is stale for the current source identity")
    if args.source_revision is not None:
        require(revision == args.source_revision,
                "native build gate differs from the configured source revision")
    _, inventory_names, inventory_digest = read_inventory(args.ctest, build_dir, args.config)
    require(document["ctest"]["inventory_sha256"] == inventory_digest and
            document["ctest"]["registered"] == len(inventory_names),
            "native build gate is stale for live CTest registration")
    verify_recorded_files(document, source_root, build_dir)
    junit_path = resolve_record(document["ctest"]["junit"], source_root, build_dir)
    junit_names, summary = inspect_junit(junit_path)
    require(junit_names == set(NATIVE_BUILD_EVIDENCE_TESTS) and
            all(document["ctest"][key] == value for key, value in summary.items()),
            "native build JUnit no longer proves the exact authority selection")
    print("native build gate verifies the complete compiled topology and exact "
          "non-render authority estate")


def expect_rejection(callback, description: str) -> None:
    try:
        callback()
    except ValueError:
        return
    raise ValueError(f"negative control accepted {description}")


def self_test() -> None:
    with tempfile.TemporaryDirectory(prefix="sirius-build-gate-") as temporary:
        root = Path(temporary)
        source = root / "source"
        build = root / "build"
        installed = root / "installed"
        source.mkdir()
        build.mkdir()
        installed.mkdir()

        tested_paths = {}
        for name in TESTED_ARTIFACTS:
            path = build / "tests" / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(f"tested:{name}\n".encode())
            tested_paths[name] = path
        product_paths = {}
        for name in RELEASE_PRODUCTS:
            path = build / "products" / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(f"product:{name}\n".encode())
            product_paths[name] = path

        names = set(NATIVE_BUILD_EVIDENCE_TESTS)
        names.update(
            f"GateFixture.Case{index}"
            for index in range(TEST_FLOOR - len(NATIVE_BUILD_EVIDENCE_TESTS))
        )
        inventory = {
            "kind": "ctestInfo",
            "version": {"major": 1, "minor": 0},
            "tests": [
                {"name": name, "properties": [{"name": "LABELS", "value": ["Mandatory"]}]}
                for name in sorted(names)
            ],
        }
        inventory_names, inventory_digest = inspect_inventory(inventory)
        junit = build / "mandatory_gate_junit.xml"
        junit.write_text(
            "<testsuite>" + "".join(f'<testcase name="{name}"/>' for name in sorted(names))
            + "</testsuite>\n", encoding="utf-8")
        log = build / "mandatory_gate_ctest.log"
        log.write_text("self-test pass\n", encoding="utf-8")
        junit_names, summary = inspect_junit(junit)
        require(junit_names == inventory_names, "self-test JUnit/inventory control diverged")
        tested = artifact_records(tested_paths, source, build)
        products = artifact_records(product_paths, source, build)
        document = {
            "schema_version": SCHEMA_VERSION,
            "status": "passed",
            "alignment_mode": "release",
            "source": {"revision": "a" * 40, "clean": True},
            "ctest": {
                "inventory_sha256": inventory_digest,
                "junit": artifact_records({"junit": junit}, source, build)["junit"],
                "log": artifact_records({"log": log}, source, build)["log"],
                "registered": len(names),
                **summary,
            },
            "inputs": {
                "operating_model_sha256": products["operating_model"]["sha256"],
                "alignment_receipt_sha256": products["alignment_receipt"]["sha256"],
            },
            "tested_artifacts": tested,
            "product_artifacts": products,
        }
        validate_document(document)
        verify_recorded_files(document, source, build)
        qualification = copy.deepcopy(document)
        qualification["alignment_mode"] = "qualification"
        validate_document(qualification)

        native_junit = build / "native_build_gate_junit.xml"
        native_junit.write_text(
            "<testsuite>" + "".join(
                f'<testcase name="{name}"/>'
                for name in sorted(NATIVE_BUILD_EVIDENCE_TESTS)
            ) + "</testsuite>\n",
            encoding="utf-8",
        )
        native_log = build / "native_build_gate_ctest.log"
        native_log.write_text("self-test native build pass\n", encoding="utf-8")
        native_names, native_summary = inspect_junit(native_junit)
        require(native_names == set(NATIVE_BUILD_EVIDENCE_TESTS),
                "self-test native authority selection diverged")
        native = {
            "schema_version": SCHEMA_VERSION,
            "kind": NATIVE_BUILD_GATE_KIND,
            "status": "passed",
            "alignment_mode": "qualification",
            "source": {"revision": "a" * 40, "clean": True},
            "ctest": {
                "inventory_sha256": inventory_digest,
                "junit": artifact_records(
                    {"junit": native_junit}, source, build
                )["junit"],
                "log": artifact_records({"log": native_log}, source, build)["log"],
                "registered": len(names),
                "selected": len(native_names),
                "selection": sorted(native_names),
                **native_summary,
            },
            "inputs": {
                "operating_model_sha256": products["operating_model"]["sha256"],
                "alignment_receipt_sha256": products["alignment_receipt"]["sha256"],
            },
            "tested_artifacts": tested,
            "product_artifacts": products,
        }
        validate_native_build_document(native)
        verify_recorded_files(native, source, build)
        weakened_native = copy.deepcopy(native)
        weakened_native["ctest"]["selection"].pop()
        weakened_native["ctest"]["selected"] -= 1
        weakened_native["ctest"]["executed"] -= 1
        expect_rejection(lambda: validate_native_build_document(weakened_native),
                         "weakened native build authority selection")
        runtime_claim = copy.deepcopy(native)
        runtime_claim["kind"] = "sirius-native-runtime-gate"
        expect_rejection(lambda: validate_native_build_document(runtime_claim),
                         "native build gate relabelled as runtime evidence")
        dirty_qualification = copy.deepcopy(qualification)
        dirty_qualification["source"]["clean"] = False
        expect_rejection(lambda: validate_document(dirty_qualification),
                         "dirty qualification source")
        incomplete_qualification = copy.deepcopy(qualification)
        incomplete_qualification["product_artifacts"].pop("trace_spv")
        expect_rejection(lambda: validate_document(incomplete_qualification),
                         "incomplete qualification product")

        for key, value in (
            ("schema_version", 2),
            ("status", "claimed"),
            ("alignment_mode", "development"),
        ):
            mutated = copy.deepcopy(document)
            mutated[key] = value
            if key == "alignment_mode":
                mutated["source"]["clean"] = False
                mutated["product_artifacts"].pop("trace_spv")
                validate_document(mutated)
            else:
                expect_rejection(lambda item=mutated: validate_document(item), f"mutated {key}")
        dirty = copy.deepcopy(document)
        dirty["source"]["clean"] = False
        expect_rejection(lambda: validate_document(dirty), "dirty release source")
        missing_test = copy.deepcopy(document)
        missing_test["tested_artifacts"].pop("sirius_core_tests")
        expect_rejection(lambda: validate_document(missing_test), "missing test executable")
        traversal = copy.deepcopy(document)
        traversal["tested_artifacts"]["sirius"]["path"] = "../escape"
        expect_rejection(lambda: verify_recorded_files(traversal, source, build),
                         "artifact path traversal")

        original = product_paths["sirius"].read_bytes()
        product_paths["sirius"].write_bytes(b"tampered\n")
        expect_rejection(lambda: verify_recorded_files(document, source, build),
                         "tampered product")
        product_paths["sirius"].write_bytes(original)

        skipped = build / "skipped.xml"
        skipped.write_text('<testsuite><testcase name="skip"><skipped/></testcase></testsuite>\n',
                           encoding="utf-8")
        expect_rejection(lambda: inspect_junit(skipped), "skipped Mandatory test")

        stamp = build / "mandatory_gate.json"
        atomic_write_json(stamp, document)
        for name, record in products.items():
            destination = installed_product_path(installed, name, "")
            destination.parent.mkdir(parents=True, exist_ok=True)
            source_path = resolve_record(record, source, build)
            destination.write_bytes(source_path.read_bytes())
        installed_stamp = installed / "share/sirius/model/mandatory_gate.json"
        installed_stamp.parent.mkdir(parents=True, exist_ok=True)
        installed_stamp.write_bytes(stamp.read_bytes())
        verify_installed(document, installed, stamp, "")
        (installed / "share/sirius/kernels/trace.spv").write_bytes(b"stale\n")
        expect_rejection(lambda: verify_installed(document, installed, stamp, ""),
                         "stale installed product")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument(
        "--action",
        choices=("run", "check", "run-native-build", "check-native-build"),
    )
    parser.add_argument("--stamp", type=Path)
    parser.add_argument("--alignment-mode", choices=("development", "qualification", "release"))
    parser.add_argument("--source-root", type=Path)
    parser.add_argument("--build-dir", type=Path)
    parser.add_argument("--source-revision")
    parser.add_argument("--source-tree-clean", choices=("true", "false"))
    parser.add_argument("--ctest", type=Path)
    parser.add_argument("--config", default="")
    parser.add_argument("--tested-artifact", action="append", default=[])
    parser.add_argument("--product-artifact", action="append", default=[])
    parser.add_argument("--installed-root", type=Path)
    parser.add_argument("--executable-suffix", default="")
    args = parser.parse_args()
    try:
        if args.self_test:
            self_test()
            print("build-gate receipt rejected missing, stale, skipped, and tampered evidence")
            return 0
        required = (
            args.action, args.stamp, args.alignment_mode, args.source_root,
            args.build_dir, args.ctest,
        )
        if any(value is None for value in required):
            raise ValueError("--action, --stamp, --alignment-mode, --source-root, "
                             "--build-dir, and --ctest are required")
        if args.action in {"run", "run-native-build"}:
            if args.source_revision is None or args.source_tree_clean is None:
                raise ValueError("run requires --source-revision and --source-tree-clean")
            if args.action == "run":
                run_gate(args)
            else:
                run_native_build_gate(args)
        elif args.action == "check-native-build":
            check_native_build_gate(args)
        else:
            check_gate(args)
        return 0
    except ValueError as error:
        parser.exit(1, f"build gate: {error}\n")


if __name__ == "__main__":
    sys.exit(main())
