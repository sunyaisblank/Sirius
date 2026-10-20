"""Preserve one CLOSED rounded-error outcome. ROOT alone invokes after final readers.

Copies owned evidence once; references private sources in immutable Git and
restores the producer stage from its single original ZIP. Never launches a
workload, performs Git operations or deletes anything.
"""
from pathlib import Path
import argparse
import hashlib
import json
import os
import shutil
import stat
import zipfile

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
TARGET = ROOT / "attestations/qualification/2c8ef086/native-coefficient-sharing-controls-benefit-outcome"
PRODUCER = "e1c618f83861c478fe18c5a22f5d9aa27f6280fd"
CANONICAL = "2c8ef08693d46f57d5dea14efdbddd4874dac684"
ZIP_SHA = "d2f372ee86508c1e6cefb8986ac4667cf9b867d62efb460e285f58b6526978de"
NAMES = ["ci-publication", "pe-audit", "publication", "runtime-source", "runtime-timestamp",
         "timestamp-evidence", "runtime-abba", "ABBA-evidence", "runtime-remaining", "all8-evidence"]
OWNED = ["out/native-coefficient-sharing-" + name for name in NAMES] + [
    "renders/native-coefficient-sharing-runtime-" + phase for phase in ("timestamp", "abba", "remaining")]
OWNED += ["out/native-coefficient-sharing-preservation"]
OWNED += ["out/native-coefficient-sharing-readers", "out/native-coefficient-sharing-precontext-stop", "out/native-coefficient-cost"]


def require(value, message):
    if not value:
        raise ValueError(message)


def local(path):
    path = Path(path)
    require(not path.is_absolute() and ".." not in path.parts, "Expected checkout-relative path")
    return ROOT / path


def identity(path):
    info = path.lstat()
    require(stat.S_ISREG(info.st_mode) and not path.is_symlink(), "Not regular: " + str(path))
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return {"bytes": info.st_size, "sha256": digest.hexdigest(), "mode": oct(stat.S_IMODE(info.st_mode))}


def document(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def json_write(path, value):
    with path.open("x", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2)
        stream.write("\n")
    path.chmod(0o644)


def inventory(excluded):
    files, roots = {}, {}
    for relative in OWNED:
        base = local(relative)
        require(not base.is_symlink(), "Aliased owned root")
        roots[relative] = {"exists": base.exists(), "files": 0}
        if not base.exists():
            continue
        require(base.is_dir() and not base.is_symlink(), "Aliased/non-directory owned root")
        for directory, names, leafs in os.walk(base, followlinks=False):
            require(all(not (Path(directory) / name).is_symlink() for name in names), "Linked owned directory")
            for name in leafs:
                path = Path(directory) / name
                if path == excluded:
                    continue
                key = str(path.relative_to(ROOT))
                require(key not in files, "Duplicate owned file")
                files[key] = identity(path)
                roots[relative]["files"] += 1
    return files, roots


def accepted(path, schema):
    value = document(path)
    require(value["schema"] == schema and value["source_revision"] == PRODUCER, "Reader source/schema mismatch")
    require(value["passed"] is True and value["control8_passed"] is True, "Control reader not positive")
    return value


def main(args):
    output = local(args.output)
    require(output.parent == HERE and not output.exists(), "Preparation output must be new in preservation directory")
    require(not TARGET.exists(), "Refuse existing preservation target")
    control_path, benefit_path = local(args.control8), local(args.benefit)
    control = accepted(control_path, "sirius-native-packed-original-controls-acceptance-v1")
    require(control["original_controls_passed"] == 1, "Expected first control8 reader")
    benefit = document(benefit_path)
    require(benefit["schema"] == "sirius-native-packed-original22-benefit-v1" and
            benefit["source_revision"] == PRODUCER and benefit["audit_passed"] is True,
            "No completed original22 audit")
    require(all(type(benefit[key]) is bool for key in ("passed", "benefit_passed", "all22_numeric_pass")) and
            benefit["passed"] == benefit["benefit_passed"] ==
            benefit["all22_numeric_pass"] and type(benefit["checks_total"]) is int and
            benefit["checks_total"] == 22 and type(benefit["checks_passed"]) is int and
            0 <= benefit["checks_passed"] <= 22 and
            benefit["passed"] == (benefit["checks_passed"] == 22), "Inconsistent original22 outcome")
    require(benefit["fixed_order"] == ["a1", "b1", "b2", "a2"] and
            benefit["processes_per_slot"] == 1 and benefit["outer_repeats_or_warmup"] == 0,
            "Wrong original trial accounting")
    remaining = ROOT / "out/native-coefficient-sharing-runtime-remaining/native"
    remaining_owners = sorted(remaining.glob("*/owner.json")) if remaining.exists() else []
    all8_path, all8 = (local(args.all8), None) if args.all8 else (None, None)
    if not benefit["passed"]:
        require(not remaining.exists() and all8_path is None, "Remaining controls after benefit refusal")
    if all8_path:
        require(benefit["passed"], "All8 requires original22 PASS")
        all8 = accepted(all8_path, "sirius-native-packed-original-controls-acceptance-v1")
        require(all8["original_controls_passed"] == 8 and
                all8["all8_original_controls_and_both_complete_readers_passed"] is True,
                "All8 reader not complete")
        benefit_row = all8["benefit_acceptance"]
        require(benefit_row["path"] == str(benefit_path.relative_to(ROOT)) and
                {k: benefit_row[k] for k in ("bytes", "sha256")} ==
                {k: identity(benefit_path)[k] for k in ("bytes", "sha256")} and len(remaining_owners) == 7,
                "All8 benefit proof/remaining owner count mismatch")
    download_path = ROOT / "out/native-coefficient-sharing-ci-publication/actual-download.json"
    download = document(download_path)
    require(download["producer_source_revision"] == PRODUCER and download["run_id"] == 39393948651 and
            download["artifact"]["id"] == 12250389795 and
            download["API_digest_CRC_safe_members_and_manifest_exact"] is True, "Wrong producer download")
    archive = local(download["archive"])
    require(identity(archive)["bytes"] == 25862941 and identity(archive)["sha256"] == ZIP_SHA and
            download["archive_sha256"] == ZIP_SHA and download["artifact"]["digest"] == "sha256:" + ZIP_SHA,
            "Producer ZIP/API digest mismatch")
    stage = local(download["stage"])
    require(stage == ROOT / "bin/windows-msvc/native-coefficient-sharing-matched-current", "Unexpected stage")
    stage_files = {str(p.relative_to(stage)) for p in stage.rglob("*") if p.is_file()}
    require(stage_files == set(download["members"]), "Stage/download member inventory differs")
    backed = {}
    with zipfile.ZipFile(archive) as handle:
        members = handle.infolist()
        require(len(members) == len(download["members"]) and
                {item.filename for item in members} == stage_files, "ZIP member inventory differs")
        for item in members:
            name = item.filename
            relative = Path(name)
            require(not relative.is_absolute() and ".." not in relative.parts and "\\" not in name and ":" not in name,
                    "Unsafe ZIP member")
            row = download["members"][name]
            current = identity(stage / relative)
            require(current == {k: row[k] for k in ("bytes", "sha256", "mode")}, "Stage byte/mode mismatch: " + name)
            require(item.external_attr == row["zip_external_attr"] and
                    hashlib.sha256(handle.read(item)).hexdigest() == row["sha256"] and
                    item.file_size == row["bytes"], "ZIP byte/mode metadata mismatch: " + name)
            backed[str((stage / relative).relative_to(ROOT))] = {
                **current, "archive_payload": "files/" + download["archive"], "member": name,
                "zip_external_attr": item.external_attr, "restore_mode": current["mode"]}
    before, roots = inventory(output)
    for proof in [control_path, benefit_path, download_path] + ([all8_path] if all8_path else []):
        require(str(proof.relative_to(ROOT)) in before, "Acceptance proof is outside owned preservation roots")
    archive_key = str(archive.relative_to(ROOT))
    require(archive_key in before and sum(row["sha256"] == ZIP_SHA for row in before.values()) == 1,
            "Original ZIP must be copied exactly once through owned roots")
    TARGET.mkdir(parents=True)
    copies = {}
    for relative, row in sorted(before.items()):
        destination = TARGET / "files" / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(local(relative), destination)
        require(identity(destination) == row, "Copy byte/mode mismatch: " + relative)
        copies["files/" + relative] = {"source": relative, **row}
    after, current_roots = inventory(output)
    require(before == after and roots == current_roots, "Owned source files/membership changed during preservation")
    require({str(p.relative_to(stage)) for p in stage.rglob("*") if p.is_file()} == stage_files,
            "Stage member inventory changed during preservation")
    for relative, row in backed.items():
        require(identity(local(relative)) == {k: row[k] for k in ("bytes", "sha256", "mode")},
                "Stage bytes/modes changed during preservation")
    actual_payloads = {str(p.relative_to(TARGET)) for p in TARGET.rglob("*") if p.is_file()}
    require(actual_payloads == set(copies), "Unexpected/missing copied payload")
    outcome = {"control8_passed": True, "original22_audit_passed": True,
               "original22_benefit_passed": benefit["passed"], "checks_passed": benefit["checks_passed"],
               "checks_total": 22, "original_controls_passed": 8 if all8 else 1,
               "all8_reader_passed": bool(all8), "remaining_native_owner_records": len(remaining_owners),
               "completed_control8_instances_timestamp_plus_sole_ABBA": 5,
               "original_distinct_control_cases_in_timestamp_and_ABBA": 1,
               "remaining_seven_state": "all8-positive" if all8 else "incomplete-or-unread" if remaining.exists() else "not-run",
               "source_adopted": False, "full_scientific_or_release_qualification_claimed": False}
    manifest = {"schema": "sirius-native-coefficient-sharing-controls-benefit-outcome-manifest-v1",
                "canonical_source_revision": CANONICAL, "producer_source_revision": PRODUCER,
                "outcome": outcome, "payloads": copies, "owned_roots": roots,
                "archive_backed_stage": backed,
                "external_Git_sources_not_copied": [
                    {"commit": "1b02b28a3a29e30f0463fa9ac1578b92f9d6e30e", "prefix": "attestations/qualification/2c8ef086/native-coefficient-sharing-source-compiler-driver", "scope": "Complete new source/single-module compiler/exact clone and ABI proof/actual zero-dispatch driver query and prior rejection cleanup. Corrected current runtime/reader/precontext-source records copied here; no compiler/private93 payload duplication."},
                    {"commit": PRODUCER, "prefixes": ["out/native-coefficient-sharing-source", "out/native-coefficient-sharing-producer-inputs", "out/native-coefficient-sharing-windows-source"],
                     "scope": "Exact already-published private candidate/arm headers35SPVs/observer/hook/CI source; recovery uses named Git, not another copy."},
                    {"commit": "1cb6e6aaff04ddd54ff7ee133b773257338f2c52", "scope": "Saved original native driver baseline only; current changed zero-dispatch query is preserved in referenced1b02 milestone, with no timing/adoption claim."},
                    {"commit": "fd530d957d2e0737f5453d88d78d979f8ed17607", "scope": "Prior constructor8+22 separate positive controls; not current frame/science proof."},
                    {"commit": "5b65b542fcc017c99503418cda889f45b201649f", "scope": "Prior adopted constructor canonical GCC source/build proof; not repeated here."},
                    {"commit": "b8ecca0ff4d7fe399db5f37922416f589e822e98", "prefix": "attestations/qualification/2c8ef086/native-constant-platform-builds-and-frame", "scope": "Prior incomplete original frame; no frame success inferred."}],
                "limits": "Snapshot of closed evidence only; host hash-only runtime dependencies are not copied or claimed backed. ZIP supplies complete staged producer recovery; manifests do not claim source adoption/full qualification."}
    json_write(TARGET / "MANIFEST.json", manifest)
    preservation = {"schema": "sirius-native-coefficient-sharing-controls-benefit-outcome-preservation-v1",
                    "preservation_audit_completed": True, "manifest": identity(TARGET / "MANIFEST.json"),
                    "outcome": outcome, "payload_count": len(copies),
                    "payload_bytes": sum(row["bytes"] for row in copies.values()),
                    "archive_backed_stage_files": len(backed), "producer_ZIP_stored_once": True,
                    "source_inventory_before_equals_after": True, "stage_inventory_bytes_modes_before_equals_after": True,
                    "copy_bytes_modes_verified": True,
                    "intended_evidence_parent": "1b02b28a3a29e30f0463fa9ac1578b92f9d6e30e",
                    "Git_commit_push_remote_verification_pending": True, "deletion_or_hold_release_admitted": False}
    json_write(TARGET / "PRESERVATION.json", preservation)
    json_write(output, {"schema": "sirius-native-coefficient-sharing-outcome-bundle-preparation-v1", "passed": True,
                        "bundle": str(TARGET.relative_to(ROOT)), "manifest": identity(TARGET / "MANIFEST.json"),
                        "preservation": identity(TARGET / "PRESERVATION.json"), **outcome,
                        "Git_or_cleanup_performed": False})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--control8", required=True)
    parser.add_argument("--benefit", required=True)
    parser.add_argument("--all8")
    parser.add_argument("--output", default="out/native-coefficient-sharing-preservation/root-preparation.json")
    main(parser.parse_args())
