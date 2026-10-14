"""Derive a fresh outer-budget packet from the held immutable numerical packet.

Only byte/source preparation runs here; no CTest/compiler/provider/validator.
"""
from pathlib import Path
import difflib
import hashlib
import json
import shutil
import subprocess

ROOT = Path.cwd().resolve()
WORK = Path(__file__).resolve().parent
COMMIT = "c901f1bcdf71533efbe61143b9437588ab37efdd"
TAG = "evidence/issue-53-p4-opacity-setup-incomplete"
BASE = ROOT / "attestations/scientific-core/4bcd987/p4-opacity-setup-incomplete"
OLD_NAME = "out/science-p4-ordinary-fp32-parity-cases5-9"
OLD = BASE / "files" / OLD_NAME
RUNNER = "f5d649bf5067d1d68dc9215462f37830122f0540"
LIMITS = {"seconds": 600, "rss_kib": 4194304, "reserve_kib": 2097152}


def identity(path):
    path = Path(path)
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return {"bytes": path.stat().st_size, "sha256": digest.hexdigest(),
            "resolved_path": str(path.resolve(strict=True))}


def key(path):
    return str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path)


def git(*args):
    return subprocess.check_output(["git", *args], cwd=ROOT, timeout=30)


def write(path, value):
    path.write_text(json.dumps(value, sort_keys=True, indent=2) + "\n")


def main():
    assert WORK == ROOT / "out/science-p4-ordinary-fp32-cases5-9-outer600"
    assert not (WORK / "preparation.json").exists()
    assert git("rev-parse", "HEAD").decode().strip() == RUNNER and not git("status", "--porcelain").strip()
    assert git("rev-parse", TAG + "^{commit}").decode().strip() == COMMIT
    entries = {}
    for row in git("ls-tree", "-rz", COMMIT + ":" + key(BASE)).split(b"\0"):
        if row:
            metadata, name = row.split(b"\t", 1)
            mode, kind, blob = metadata.split()
            assert kind == b"blob" and mode in (b"100644", b"100755")
            entries[name.decode()] = blob.decode()
    archived = {str(path.relative_to(BASE)): path for path in BASE.rglob("*") if path.is_file()}
    assert set(archived) == set(entries)
    for name, path in archived.items():
        assert not path.is_symlink()
        raw = path.read_bytes()
        assert hashlib.sha1(b"blob " + str(len(raw)).encode() + b"\0" + raw).hexdigest() == entries[name]
    manifest = json.loads((BASE / "manifest.json").read_text())
    for name, expected in manifest["files"].items():
        actual = identity(BASE / name)
        assert all(actual[field] == expected[field] for field in ("bytes", "sha256"))
    original = json.loads((OLD / "preparation.json").read_text())
    original_actions = json.loads((OLD / "actions.json").read_text())
    assert len(original["pins"]) == 470
    assert original["original_case_ordinals"] == list(range(5, 10))
    assert len(original["cases"]) == 5 and len(original_actions) == 5
    assert original["limits"] == {**LIMITS, "seconds": 300}
    result = json.loads((BASE / "result.json").read_text())
    assert result["held_cases_5_to_9"]["prepared_only"] and result["held_cases_5_to_9"]["runtime_actions"] == 0
    assert not any(name.startswith("files/" + OLD_NAME + "/parity-") for name in archived)
    recovery = WORK / "recovery"
    recovery.mkdir()
    assets = set(archived.values())
    relocations = []
    for name, expected in original["pins"].items():
        path = BASE / "files" / name if name.startswith(OLD_NAME + "/") else ROOT / name
        actual = identity(path)
        assert all(actual[field] == expected[field] for field in ("bytes", "sha256")), name
        if not name.startswith(OLD_NAME + "/"):
            assert actual == expected, name
        assets.add(path)
        relocations.append({"original": name, "current": key(path), "identity": actual})
    write(recovery / "held470-input-relocations.json", relocations)
    for name in ("run_owned.py", "owned_guard.py"):
        shutil.copy2(OLD / name, WORK / name)
    reader = (OLD / "run_parity.py").read_text()
    (recovery / "held300-run_parity.py").write_text(reader)
    anchor = 'LIMITS = {"seconds": 300, "rss_kib": 4194304, "reserve_kib": 2097152}'
    assert reader.count(anchor) == 1
    changed = reader.replace(anchor, anchor.replace('"seconds": 300', '"seconds": 600'), 1)
    (WORK / "run_parity.py").write_text(changed)
    (recovery / "outer-budget-reader.patch").write_text("".join(difflib.unified_diff(reader.splitlines(True), changed.splitlines(True),
         fromfile="held300/run_parity.py", tofile="outer600/run_parity.py")))
    context = {name: value for name, value in original.items() if name != "pins"}
    context["limits"] = LIMITS
    context["finite_scope"] = "Previously unlaunched original ordinary-FP32 P4 cases5..9, separate600s outer diagnostic owners only; no complete P4/science/Mandatory/platform/release or performance claim."
    context["scheduling"]["case4_pass_before5"] = "Superseded root scheduling decision: cases3/4 remain incomplete; root now starts previously unlaunched5 under a fresh600s outer diagnostic allowance. No fabricated case4 acceptance."
    context["scheduling"]["change"] = "Only imported outer diagnostic allowance300s to600s in a fresh packet; original selected5..9 sequencing, all numerical/provider/source/owner checks and one-second actual submissions remain unchanged."
    context["held300_baseline"] = {"commit": COMMIT, "tag": TAG, "bundle": key(BASE),
        "manifest": identity(BASE / "manifest.json"), "raw_Git_files_joined": len(archived),
        "manifest_files_joined": len(manifest["files"]), "original_pins_relocated": 470,
        "cases5_through9_runtime_actions": 0, "old300_packet_modified": False}
    context["outer_diagnostic_contract"] = {"seconds": 600, "rss_kib": 4194304, "reserve_kib": 2097152,
        "whole_case_300s_is_original_test_timeout": False, "actual_sample_submit_ms_maximum": 1000,
        "zero_record_initialization_separate_from_actual_submission": True,
        "whole_child_wall_report": "Every parity-N/owner.json reports wall_seconds including first initialization, actual sample work, natural teardown and retirement; input sealing and XML postprocessing are separate.",
        "initial_conditions": "Stock full module; unchanged accepted backend; unset FP32 precision; explicit LVP index0; Mesa cacheoff; original call order, failure boundaries and allocation inequality; no prewarming or omitted first-call work.",
        "hard_cold_Kerr_frame_300s_unchanged": True,
        "case3_or4_retry": False, "loop_policy_change": False,
        "performance_or_completion_deadline_claimed": False}
    write(WORK / "context.json", context)
    assets.update(path for path in WORK.rglob("*") if path.is_file())
    frozen = sorted({key(path) for path in assets} | {key(WORK / "preparation.json")})
    actions = {}
    for label, original_action in original_actions.items():
        action = dict(original_action)
        assert action["limits"] == {**LIMITS, "seconds": 300}
        action["limits"] = LIMITS
        action["scope"] = context["finite_scope"] + " Case:" + next(case["name"] for case in context["cases"] if case["label"] == label)
        action["frozen_inputs"] = frozen
        actions[label] = action
    write(WORK / "actions.json", actions)
    assets.add(WORK / "actions.json")
    context["pins"] = {key(path): identity(path) for path in sorted(assets)}
    write(WORK / "preparation.json", context)
    print(json.dumps({"prepared_only": True, "pins": len(context["pins"]),
                      "preparation": identity(WORK / "preparation.json"), "owner_limits": LIMITS}))


if __name__ == "__main__":
    main()
