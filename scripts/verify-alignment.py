#!/usr/bin/env python3
"""Build or check the revision-bound Sirius alignment receipt.

The individual attestation verifier proves one evidence bundle.  This tool
closes the next authority boundary: it admits exactly one witness for every
external operating domain, requires all witnesses to name the source revision
being built, and emits the deterministic receipt embedded in the executable.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import re
import subprocess
import sys
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
ATTESTATION_VERIFIER = ROOT / "scripts" / "verify-attestation.py"
REQUIRED_DOMAINS = (
    "macos-moltenvk",
    "macos-native-build",
    "physical-imax-5616x4096",
    "physical-radeon-780m",
    "viewer-native-window-input",
    "windows-native-build",
    "windows-native-vulkan",
    "wsl2-dozen",
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def load_attestation_verifier():
    spec = importlib.util.spec_from_file_location(
        "sirius_attestation_verifier", ATTESTATION_VERIFIER
    )
    require(spec is not None and spec.loader is not None,
            "could not load the attestation verifier")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def inspect_operating_model(path: Path) -> tuple[str, tuple[str, ...]]:
    """Bind admission to the external profiles declared by the model itself."""
    require(path.is_file(), f"operating model is missing: {path}")
    payload = path.read_bytes()
    require(
        payload.endswith(b"\n") and b"\r" not in payload,
        "operating model must use canonical LF-terminated UTF-8 bytes",
    )
    try:
        model = json.loads(payload.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ValueError(f"operating model JSON is unreadable: {error}") from error
    require(
        isinstance(model, dict)
        and model.get("schema_version") == 3
        and model.get("method") == "adversarial-claim-ledger",
        "operating model has an unsupported schema or method",
    )
    require(model.get("source_available_test_floor") == 700,
            "operating model has an ungoverned source-available test floor")
    capabilities = model.get("capability_contracts")
    require(isinstance(capabilities, list),
            "operating model capability contracts are malformed")
    domains: set[str] = set()
    for capability in capabilities:
        require(isinstance(capability, dict),
                "operating model contains a malformed capability")
        if capability.get("state") != "attestation_required":
            continue
        profiles = capability.get("profiles")
        require(isinstance(profiles, list) and profiles,
                "attestation-required capability has no external profiles")
        for profile in profiles:
            require(isinstance(profile, str) and profile,
                    "attestation-required capability has a malformed profile")
            require(profile not in domains,
                    f"external profile has multiple capability authorities: {profile}")
            domains.add(profile)
    derived = tuple(sorted(domains))
    require(
        derived == REQUIRED_DOMAINS,
        "operating-model external profiles differ from release admission policy: "
        f"model={list(derived)}, admission={list(REQUIRED_DOMAINS)}",
    )
    return hashlib.sha256(payload).hexdigest(), derived


def require_source_operating_model(
    source_root: Path, expected_sha256: str
) -> None:
    source_sha256, _ = inspect_operating_model(
        source_root / "tests" / "operating_model.json"
    )
    require(source_sha256 == expected_sha256,
            "selected operating model differs from the exact source model")


def discover_attestations(root: Path | None) -> list[Path]:
    if root is None:
        return []
    require(root.is_dir(), f"attestation root is not a directory: {root}")
    attestations = []
    for path in sorted(root.rglob("*.json")):
        try:
            candidate = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as error:
            raise ValueError(f"JSON evidence is unreadable: {path}: {error}") from error
        if (
            isinstance(candidate, dict)
            and candidate.get("schema_version") == 1
            and "domains" in candidate
            and "source_revision" in candidate
        ):
            attestations.append(path)
    return attestations


def verified_records(
    paths: list[Path], source_revision: str, operating_model_sha256: str
) -> list[dict]:
    verifier = load_attestation_verifier()
    records = []
    for path in paths:
        document = json.loads(path.read_text(encoding="utf-8"))
        verifier.verify_document(document, path, operating_model_sha256)
        require(
            document["source_revision"] == source_revision,
            f"attestation {path} names revision {document['source_revision']}, "
            f"not {source_revision}",
        )
        records.append(
            {
                "domains": sorted(document["domains"]),
                "completed_utc": document["completed_utc"],
                "source_revision": document["source_revision"],
                "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            }
        )
    return records


def build_receipt(
    records: list[dict], source_revision: str, source_tree_clean: bool, mode: str,
    operating_model_sha256: str,
) -> dict:
    require(mode in {"development", "qualification", "release"}, "unknown alignment mode")
    require(
        re.fullmatch(r"[0-9a-f]{40,64}", source_revision) is not None,
        "source revision must be a full lowercase hexadecimal revision",
    )
    require(
        re.fullmatch(r"[0-9a-f]{64}", operating_model_sha256) is not None,
        "operating model must have a valid SHA-256",
    )
    admitted = []
    domain_counts = {domain: 0 for domain in REQUIRED_DOMAINS}
    for record in records:
        domains = record.get("domains")
        require(isinstance(domains, list) and domains,
                "an admitted record has no domains")
        require(len(domains) == len(set(domains)),
                "an admitted record duplicates a domain")
        require(
            isinstance(record.get("completed_utc"), str)
            and record["completed_utc"],
            "an admitted record has no completion time",
        )
        require(
            record.get("source_revision") == source_revision,
            "an admitted record names a different source revision",
        )
        require(
            isinstance(record.get("sha256"), str)
            and re.fullmatch(r"[0-9a-f]{64}", record["sha256"]) is not None,
            "an admitted record has no valid SHA-256",
        )
        for domain in domains:
            require(domain in domain_counts, f"unknown alignment domain: {domain}")
            domain_counts[domain] += 1
            require(
                domain_counts[domain] == 1,
                f"alignment domain has more than one authority: {domain}",
            )
            admitted.append(domain)

    admitted = sorted(admitted)
    pending = sorted(set(REQUIRED_DOMAINS) - set(admitted))
    satisfied = source_tree_clean and not pending
    if mode in {"qualification", "release"}:
        require(source_tree_clean, f"{mode} alignment requires a clean source tree")
    if mode == "release":
        require(
            not pending,
            "release alignment is missing required domains: " + ", ".join(pending),
        )

    canonical_records = sorted(
        records, key=lambda record: (record["domains"], record["sha256"])
    )
    return {
        "schema_version": 1,
        "method": "revision-bound-attestation-set",
        "alignment_mode": mode,
        "source_revision": source_revision,
        "source_tree_clean": source_tree_clean,
        "operating_model": {
            "schema_version": 3,
            "sha256": operating_model_sha256,
            "external_domains": list(REQUIRED_DOMAINS),
        },
        "required_domains": list(REQUIRED_DOMAINS),
        "admitted_domains": admitted,
        "pending_domains": pending,
        "records": canonical_records,
        "ultimate_ideal": {
            "required": True,
            "satisfied": satisfied,
            "state": "aligned" if satisfied else "attestation_required",
        },
    }


def render_receipt(receipt: dict) -> str:
    return json.dumps(receipt, indent=2, sort_keys=True) + "\n"


def render_receipt_bytes(receipt: dict) -> bytes:
    """Emit one cross-platform byte identity for a deterministic receipt."""
    return render_receipt(receipt).encode("utf-8")


def inspect_git_source(source_root: Path) -> tuple[str, bool]:
    require((source_root / ".git").exists(),
            "alignment source root is not a Git worktree")
    try:
        revision = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=source_root,
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()
        status = subprocess.run(
            ["git", "status", "--porcelain", "--untracked-files=normal"],
            cwd=source_root,
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()
    except (OSError, subprocess.CalledProcessError) as error:
        raise ValueError(f"could not inspect alignment source worktree: {error}") from error
    return revision, status == ""


def self_test() -> None:
    revision = "0" * 40
    model_digest, model_domains = inspect_operating_model(
        ROOT / "tests" / "operating_model.json"
    )
    require(model_domains == REQUIRED_DOMAINS,
            "operating-model domain authority was not derived exactly")
    require_source_operating_model(ROOT, model_digest)
    try:
        require_source_operating_model(ROOT, "f" * 64)
    except ValueError:
        pass
    else:
        raise ValueError("alignment accepted a substituted source model")
    complete = [
        {
            "domains": list(REQUIRED_DOMAINS),
            "completed_utc": "2026-08-27T00:00:00+00:00",
            "source_revision": revision,
            "sha256": "1" * 64,
        }
    ]
    aligned = build_receipt(complete, revision, True, "release", model_digest)
    require(aligned["ultimate_ideal"]["satisfied"] is True,
            "complete clean evidence did not align")

    controls = [
        ([], revision, True, "release"),
        (complete, revision, False, "release"),
        (
            complete
            + [{
                "domains": [REQUIRED_DOMAINS[0]],
                "completed_utc": "2026-08-27T00:00:00+00:00",
                "source_revision": revision,
                "sha256": "2" * 64,
            }],
            revision,
            True,
            "release",
        ),
        (
            [{**complete[0], "source_revision": "2" * 40}],
            revision,
            True,
            "release",
        ),
    ]
    for records, candidate_revision, clean, mode in controls:
        try:
            build_receipt(records, candidate_revision, clean, mode, model_digest)
        except ValueError:
            continue
        raise ValueError("alignment negative control was accepted")

    development = build_receipt([], revision, False, "development", model_digest)
    require(
        development["pending_domains"] == list(REQUIRED_DOMAINS)
        and development["ultimate_ideal"]["satisfied"] is False,
        "development receipt concealed pending external domains",
    )
    qualification = build_receipt([], revision, True, "qualification", model_digest)
    require(
        qualification["pending_domains"] == list(REQUIRED_DOMAINS)
        and qualification["ultimate_ideal"]["satisfied"] is False,
        "qualification receipt concealed pending external domains",
    )
    try:
        build_receipt([], revision, False, "qualification", model_digest)
    except ValueError:
        pass
    else:
        raise ValueError("qualification alignment accepted a dirty source tree")
    with tempfile.TemporaryDirectory(prefix="sirius-alignment-") as directory:
        receipt_path = Path(directory) / "receipt.json"
        receipt_payload = render_receipt_bytes(development)
        require(
            receipt_payload.endswith(b"\n") and b"\r" not in receipt_payload,
            "alignment receipt is not canonical LF-terminated UTF-8",
        )
        receipt_path.write_bytes(receipt_payload)
        require(
            receipt_path.read_bytes() == receipt_payload,
            "alignment receipt is not deterministic",
        )
        source_model_path = ROOT / "tests" / "operating_model.json"
        source_model_payload = source_model_path.read_bytes()
        noncanonical_models = {
            "CRLF": source_model_payload.replace(b"\n", b"\r\n"),
            "unterminated": source_model_payload.rstrip(b"\n"),
            "non-UTF-8": b"\xff\n",
        }
        for label, payload in noncanonical_models.items():
            invalid_model_path = Path(directory) / f"operating_model-{label}.json"
            invalid_model_path.write_bytes(payload)
            try:
                inspect_operating_model(invalid_model_path)
            except ValueError:
                continue
            raise ValueError(f"alignment accepted noncanonical {label} model identity")
        drifted_model_path = Path(directory) / "operating_model.json"
        drifted_model = json.loads(
            source_model_path.read_text(encoding="utf-8")
        )
        for capability in drifted_model["capability_contracts"]:
            if capability.get("id") == "viewer_native_window_input":
                capability["profiles"] = ["viewer-native-window"]
        drifted_model_path.write_bytes(
            (json.dumps(drifted_model, indent=2, sort_keys=True) + "\n").encode("utf-8")
        )
        try:
            inspect_operating_model(drifted_model_path)
        except ValueError:
            pass
        else:
            raise ValueError("alignment negative control accepted model/admission drift")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=("development", "qualification", "release"))
    parser.add_argument("--source-revision")
    parser.add_argument("--source-tree-clean", choices=("true", "false"))
    parser.add_argument("--operating-model", type=Path)
    parser.add_argument("--attestation-root", type=Path)
    parser.add_argument("--source-root", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    try:
        if args.self_test:
            self_test()
            print(
                "alignment admission rejected incomplete, dirty, ambiguous, and "
                "model-drifted evidence"
            )
            return
        require(args.mode is not None, "--mode is required")
        require(args.source_revision is not None, "--source-revision is required")
        require(args.source_tree_clean is not None, "--source-tree-clean is required")
        require(args.operating_model is not None, "--operating-model is required")
        require(args.output is not None, "--output is required")
        model_digest, _ = inspect_operating_model(args.operating_model)
        if args.mode in {"qualification", "release"}:
            require(args.source_root is not None,
                    f"{args.mode} alignment requires --source-root")
        if args.source_root is not None:
            require_source_operating_model(args.source_root, model_digest)
            actual_revision, actual_clean = inspect_git_source(args.source_root)
            require(
                actual_revision == args.source_revision,
                "source revision changed after configuration; reconfigure the build",
            )
            require(
                actual_clean == (args.source_tree_clean == "true"),
                "source tree cleanliness changed after configuration; reconfigure the build",
            )
        paths = discover_attestations(args.attestation_root)
        records = verified_records(paths, args.source_revision, model_digest)
        receipt = build_receipt(
            records,
            args.source_revision,
            args.source_tree_clean == "true",
            args.mode,
            model_digest,
        )
        payload = render_receipt_bytes(receipt)
        if args.check:
            require(args.output.is_file(), f"alignment receipt is missing: {args.output}")
            require(
                args.output.read_bytes() == payload,
                "alignment receipt is stale; reconfigure the build",
            )
            print(
                f"alignment receipt is current: {len(receipt['admitted_domains'])}/"
                f"{len(REQUIRED_DOMAINS)} domains admitted"
            )
            return
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_bytes(payload)
        print(
            f"alignment receipt generated: {len(receipt['admitted_domains'])}/"
            f"{len(REQUIRED_DOMAINS)} domains admitted"
        )
    except (OSError, ValueError, json.JSONDecodeError) as error:
        parser.exit(1, f"alignment: {error}\n")


if __name__ == "__main__":
    main()
