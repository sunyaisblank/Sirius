#!/usr/bin/env python3
"""Prepare separated original CPU controls; execute no tests."""
from pathlib import Path
import ast
import datetime
import hashlib
import json
import os
import re
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[2]
WORK = Path(__file__).resolve().parent
PRIOR = ROOT / "attestations/scientific-core/4bcd987/p6-cpu-writer-diagnostic-limit"
PRIOR_MANIFEST_SHA = "8b5b860301a6b697ce2505ab1ac141353840abe2589a75e9403ca3db041cab97"
ORIGINAL = PRIOR / "files/out/p6-cpu-writer-current-source"
KERR = "RenderSessionProbe.CpuKerrRenderProducesValidPpmThroughTheOwnedWriter"


def identity(path):
    raw = path.read_bytes()
    return {"bytes": len(raw), "sha256": hashlib.sha256(raw).hexdigest()}


def git(*args):
    return subprocess.check_output(["git", *args], cwd=ROOT).decode().strip()


def key(path):
    return str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path)


def write(name, value):
    (WORK / name).write_text(json.dumps(value, indent=2) + "\n")


assert not (WORK / "preparation.json").exists()
assert git("status", "--porcelain", "--untracked-files=normal") == ""
assert identity(PRIOR / "manifest.json")["sha256"] == PRIOR_MANIFEST_SHA
manifest = json.loads((PRIOR / "manifest.json").read_text())
records = {row["path"]: {key: row[key] for key in ("bytes", "sha256")} for row in manifest["files"]}
base = json.loads((ORIGINAL / "preparation.json").read_text())
original_owner_path = ORIGINAL / "original-five-p6-cpu-writer-controls/owner.json"
original_owner_pin = identity(original_owner_path)
original_owner = json.loads(original_owner_path.read_text())
assert original_owner["stop_reason"] == "diagnostic_timeout" and original_owner["child_exit"] == -15
assert original_owner["cleanup_errors"] == [] and original_owner["remaining_owned_processes"] == []
assert original_owner["source_unchanged"] and original_owner["whole_inputs_unchanged"]
assert not (ORIGINAL / "original-five-p6-cpu-writer-controls/ctest.xml").exists()
runner = git("rev-parse", "HEAD")
assert runner == base["runner_source_revision"]
assert set(git("diff", "--name-only", base["product_source_revision"], "HEAD").splitlines()) == set(base["diagnostic_only_difference"])
pins = {name: pin for name, pin in base["pins"].items() if not name.startswith("out/p6-cpu-writer-current-source/")}
for name, pin in pins.items():
    assert identity(ROOT / name) == pin, name
for name in ("run_owned.py", "owned_guard.py"):
    path = ORIGINAL / name
    assert identity(path) == records[str(path.relative_to(PRIOR))]
    shutil.copyfile(path, WORK / name)
reference_files = sorted((ORIGINAL / "reference").iterdir())
(WORK / "reference").mkdir()
for path in reference_files:
    assert path.is_file() and identity(path) == records[str(path.relative_to(PRIOR))]
    shutil.copyfile(path, WORK / "reference" / path.name)
all_names = base["source_and_registration_join"]["expected_names"]
assert len(all_names) == 5 and KERR in all_names
short_names = [name for name in all_names if name != KERR]
variants = {
    "short": {"case": "original-four-short-p6-cpu-writer-controls", "names": short_names, "seconds": 300},
    "kerr": {"case": "original-kerr-ppm-p6-cpu-writer-control", "names": [KERR], "seconds": 900},
}
registration = ROOT / base["source_and_registration_join"]["original_registration_file"]
registration_lines = registration.read_text().splitlines()
for line in base["source_and_registration_join"]["original_five_registration_lines"]:
    assert line in registration_lines
for name in all_names:
    properties = [line for line in registration_lines if line.startswith("set_tests_properties(") and name in line]
    assert len(properties) == 1 and "TIMEOUT" not in properties[0]
local = [WORK / name for name in ("prepare.py", "run_control.py", "run_owned.py", "owned_guard.py")]
local += [WORK / "reference" / path.name for path in reference_files]
prior_inputs = [PRIOR / "manifest.json", PRIOR / "disposition.json", ORIGINAL / "preparation.json", original_owner_path]
for path in prior_inputs + local:
    pins[key(path)] = identity(path)
ctest = base["source_and_registration_join"]["cache_required"]
cache_path = ROOT / "bin/linux-gcc/CMakeCache.txt"
match = re.search(r"^CMAKE_CTEST_COMMAND:INTERNAL=(.+)$", cache_path.read_text(), re.M)
assert match
ctest_path = match.group(1)
actions = {}
for selection, variant in variants.items():
    actions[variant["case"]] = {
        "revision": runner, "kind": "check", "expected_exit": 0,
        "limits": {"seconds": variant["seconds"], "rss_kib": 4194304, "reserve_kib": 2097152},
        "scope": "Original CPU " + selection + " writer/film selection with unchanged assertions, fixtures, producer, resources and cleanup; the outer budget is diagnostic, not a scientific threshold or the 300s Vulkan frame criterion.",
        "argv": [ctest_path, "--test-dir", "bin/linux-gcc", "-V", "--output-on-failure", "--no-tests=error",
                 "--output-junit", "{folder}/ctest.xml", "--output-log", "{folder}/ctest.log",
                 "-R", "^(" + "|".join(re.escape(name) for name in variant["names"]) + ")$"],
        "frozen_inputs": list(pins) + [key(WORK / "actions.json"), key(WORK / "preparation.json")],
    }
write("actions.json", actions)
pins[key(WORK / "actions.json")] = identity(WORK / "actions.json")
write("preparation.json", {
    "prepared_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
    "product_source_revision": base["product_source_revision"], "runner_source_revision": runner,
    "diagnostic_only_difference": base["diagnostic_only_difference"], "executed": False,
    "scope": "Two separate one-shot original CPU selections; preserve prior incomplete attempt, no full P6, forensic identity, retained release-image or Mandatory qualification.",
    "variants": variants,
    "diagnostic_changes_only": {"prior_shared_budget_seconds": 300, "short_outer_seconds": 300, "kerr_outer_seconds": 900,
        "original_cpu_test_TIMEOUT_added": False, "original_assertions_or_scene_changed": False,
        "vulkan_frame_acceptance_seconds": 300, "vulkan_frame_acceptance_changed": False,
        "logging": "standard CTest -V plus --output-log", "image_capture_hooks": False},
    "prior_incomplete_attempt": {"archive": key(PRIOR), "manifest_sha256": PRIOR_MANIFEST_SHA,
        "raw_owner": key(original_owner_path), "raw_owner_identity": original_owner_pin,
        "stop_reason": original_owner["stop_reason"], "wall_seconds": original_owner["wall_seconds"],
        "child_exit": original_owner["child_exit"], "scientific_assertion_verdict": "unavailable"},
    "finite_acceptance_scope": {
        "short": {name: value for name, value in base["finite_acceptance_scope"].items() if name != "kerr_ppm_scene"},
        "kerr": {"kerr_ppm_scene": base["finite_acceptance_scope"]["kerr_ppm_scene"]},
    },
    "reference": {**base["reference"], "source_snapshots": {key(WORK / "reference" / path.name): identity(WORK / "reference" / path.name) for path in reference_files}},
    "private_tmpdirs": {selection: {"path": str(WORK / ("temporary-" + selection)), "created_only_at_execution": True,
        "mode": "0700", "owner_uid": os.getuid(), "cleanup_interception": False} for selection in variants},
    "source_and_registration_join": base["source_and_registration_join"],
    "owner_reuse": {name: {"source": key(ORIGINAL / name), **identity(WORK / name)} for name in ("run_owned.py", "owned_guard.py")},
    "pins": pins,
})
for path in local:
    if path.suffix == ".py":
        ast.parse(path.read_text(), filename=str(path))
assert identity(original_owner_path) == original_owner_pin
print(json.dumps({"prepared_without_execution": True, "pins": len(pins), "preparation": identity(WORK / "preparation.json"),
                  "original_owner_unchanged": original_owner_pin,
                  "commands": [["python3", "-B", key(WORK / "run_control.py"), selection] for selection in variants]}, indent=2))
