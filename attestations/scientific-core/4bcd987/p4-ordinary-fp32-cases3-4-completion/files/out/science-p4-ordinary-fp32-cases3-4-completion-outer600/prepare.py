"""Prepare two original incomplete P4 witnesses; no build/test/device execution."""
from pathlib import Path
import ast
import copy
import difflib
import hashlib
import json
import re
import shutil
import subprocess

ROOT = Path.cwd().resolve()
WORK = Path(__file__).resolve().parent
RUNNER = "f5d649bf5067d1d68dc9215462f37830122f0540"
SPECS = [
    ("p4-ordinary-fp32-isolated-partial", "0d8beec5fdb1406d35add385cfa6f2d9f62848d2", "evidence/issue-53-p4-ordinary-fp32-isolated-partial"),
    ("p4-opacity-setup-incomplete", "c901f1bcdf71533efbe61143b9437588ab37efdd", "evidence/issue-53-p4-opacity-setup-incomplete"),
    ("p4-ordinary-fp32-cases5-9", "76a7d0281c1a47844d86170914dc18f72e27769a", "evidence/issue-53-p4-ordinary-fp32-cases5-9"),
]
BASE = ROOT / "attestations/scientific-core/4bcd987"
OLD_PREFIX = "out/science-p4-ordinary-fp32-cases5-9-outer600"
OLD = BASE / SPECS[2][0] / "files" / OLD_PREFIX
LIMITS = {"seconds": 600, "rss_kib": 4194304, "reserve_kib": 2097152}


def require(ok, message):
    if not ok:
        raise ValueError(message)


def git(*args):
    return subprocess.check_output(["git", *args], cwd=ROOT, timeout=30)


def key(path):
    path = Path(path)
    return str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path)


def identity(path):
    path = Path(path)
    with path.open("rb") as stream:
        sha = hashlib.file_digest(stream, "sha256").hexdigest()
    return {"bytes": path.stat().st_size, "sha256": sha,
            "resolved_path": str(path.resolve(strict=True))}


def write(path, value):
    path.write_text(json.dumps(value, sort_keys=True, indent=2) + "\n")


def relocate(name):
    return key(OLD) + name[len(OLD_PREFIX):] if name.startswith(OLD_PREFIX + "/") else name


def main():
    require(WORK == ROOT / "out/science-p4-ordinary-fp32-cases3-4-completion-outer600", "fixed fresh packet root")
    require(not (WORK / "preparation.json").exists(), "already sealed")
    require(git("rev-parse", "HEAD").decode().strip() == RUNNER and not git("status", "--porcelain").strip(), "clean pinned runner")
    archive_files = {}
    archive_joins = {}
    for folder, commit, tag in SPECS:
        bundle = BASE / folder
        require(git("rev-parse", tag + "^{commit}").decode().strip() == commit, "immutable evidence tag changed")
        entries = {}
        for row in git("ls-tree", "-rz", commit + ":" + key(bundle)).split(b"\0"):
            if row:
                metadata, name = row.split(b"\t", 1)
                mode, kind, blob = metadata.split()
                require(kind == b"blob" and mode in (b"100644", b"100755"), "unexpected archive entry")
                entries[name.decode()] = blob.decode()
        files = {str(p.relative_to(bundle)): p for p in bundle.rglob("*") if p.is_file()}
        require(set(files) == set(entries), "whole raw archive/Git coverage differs")
        manifest = json.loads((bundle / "manifest.json").read_text())
        for name, path in files.items():
            require(not path.is_symlink(), "unexpected archive symlink")
            raw = path.read_bytes()
            require(hashlib.sha1(b"blob " + str(len(raw)).encode() + b"\0" + raw).hexdigest() == entries[name], "immutable Git bytes differ:" + name)
            pin = identity(path)
            if name in manifest["files"]:
                require(all(pin[k] == manifest["files"][name][k] for k in ("bytes", "sha256")), "whole-byte manifest differs")
            archive_files[key(path)] = pin
        archive_joins[folder] = {"revision": commit, "tag": tag, "bundle": key(bundle),
             "whole_raw_files": len(files), "manifest_files": len(manifest["files"]),
             "git_blobs": {key(bundle / name): blob for name, blob in entries.items()},
             "manifest": identity(bundle / "manifest.json")}

    original = json.loads((OLD / "preparation.json").read_text())
    original_actions = json.loads((OLD / "actions.json").read_text())
    require(len(original["pins"]) == 510 and original["limits"] == LIMITS, "published five-case baseline")
    require(json.loads((BASE / SPECS[2][0] / "result.json").read_text())["original_cases_passed"] == 5, "material five-case evidence")
    pins = dict(archive_files)
    relocations = {}
    for name, expected in original["pins"].items():
        target = relocate(name)
        require(not target.startswith("out/"), "no live or deleted scratch dependency")
        actual = identity(ROOT / target)
        require(all(actual[k] == expected[k] for k in ("bytes", "sha256")), "published input changed:" + name)
        if target == name:
            require(actual == expected, "resolved input identity changed:" + name)
        else:
            require(target in archive_files, "relocated input lacks immutable Git bytes")
            relocations[name] = target
        pins[target] = actual
    first = BASE / SPECS[0][0] / "files/out/science-p4-ordinary-fp32-parity-isolated"
    original_nine = json.loads((first / "preparation.json").read_text())
    cases = copy.deepcopy(original_nine["cases"][2:4])
    require([c["label"] for c in cases] == ["parity-3", "parity-4"], "only original cases3/4")
    require([(c["opcode"], c["records"], c["minimum_allocated_bytes"]) for c in cases] == [(2, 4, 1296), (9, 3, 976)], "original data/allocation selection")
    source = (ROOT / "tests/backend/kernel_parity_test.cpp").read_text()
    registration = (ROOT / "bin/linux-gcc/tests/backend/sirius_backend_tests[1]_tests.cmake").read_text()
    properties = {}
    for case in cases:
        require(case["original_registration"] in registration, "original CTest registration changed")
        require(case["name"] in (ROOT / "tests/operating_model.json").read_text(), "original P4 model witness")
        body = source.split("TEST(KernelParity, " + case["name"].split(".")[1] + ") {", 1)[1].split("\nTEST(", 1)[0]
        require("RunProbe(*f.device, f.kernel" in body and "OpenProbe()" in body, "original ordinary module/helper")
        line = next(line for line in registration.splitlines() if line.startswith("set_tests_properties(") and case["name"] in line)
        require("TIMEOUT" not in line and "sirius_vulkan_device" in line, "original lock/no whole-case timeout")
        properties[case["name"]] = line
    require("ASSERT_LE(timing.submit_wait_ms, sirius::render::kDispatchStopMs);" in source, "original actual-submission gate")
    for filename, digest in (("run_owned.py", "492c75609a0dc8ff5eaa0299f7206820bb8a42f85c3ecc5b9b0af98f2496ba9f"),
                             ("owned_guard.py", "747fe1e0bcaa0cad8fb867f22a30935b1c367dc0e2c0ad8204a87f6593899056")):
        require(identity(OLD / filename)["sha256"] == digest, "mature owner/guard identity")
        shutil.copy2(OLD / filename, WORK / filename)
    reader = (OLD / "run_parity.py").read_text()
    changed = reader
    for before, after in (
        ('range(5, 10)], "select one original ordinal 5..9"', 'range(3, 5)], "select one original ordinal 3..4"'),
        ('index = int(sys.argv[1]) - 5', 'index = int(sys.argv[1]) - 3'),
    ):
        require(changed.count(before) == 1, "private ordinal selection locus")
        changed = changed.replace(before, after)
    ast.parse(changed)
    (WORK / "run_parity.py").write_text(changed)
    recovery = WORK / "recovery"
    recovery.mkdir()
    (recovery / "selection-reader.patch").write_text("".join(difflib.unified_diff(reader.splitlines(True), changed.splitlines(True),
        fromfile="published-cases5-9/run_parity.py", tofile="completion-cases3-4/run_parity.py")))
    write(recovery / "archive-bindings.json", {"archives": archive_joins, "published510_input_relocations": relocations})

    context = copy.deepcopy({k: v for k, v in original.items() if k != "pins"})
    context.update(cases=cases, original_case_ordinals=[3, 4], original_registration_properties=properties,
        prepared_only=True, finite_scope="Two original ordinary-FP32 P4 witnesses, cases3/4, one new completion attempt each after prior setup-incomplete300s observations; separate600s outer owners. No fullP4/science/Mandatory/domain/release/performance/adoption claim.")
    context["finite_references"] = {case["name"]: {"records": case["records"], "reference": case["reference"],
         "source_line": case["source_line"], "requested_minimum_bytes": case["minimum_allocated_bytes"],
         "allocation_maximum_bytes": case["maximum_allocated_bytes"], "actual_submissions": 1,
         "actual_dispatch_ms_maximum": 1000} for case in cases}
    context["scheduling"] = {"selected_original_ordinals": [3, 4], "separate_serial_owners": True,
         "stop_on_first_incomplete_or_failure": True, "preceding_selected_case_required_for4": 3,
         "no_preceding_old_case_acceptance_fabricated": True, "already_passed_cases_not_repeated": [1, 2, 5, 6, 7, 8, 9],
         "after_live_scene2_retirement_only": True, "one_new_attempt_per_case": True,
         "old_incomplete_records_preserved_unchanged": True,
         "root_plan_revision": "Root explicitly will prospectively revise issue53 prior no-retry plan on the original source-deadline audit plus five cold76a7d028 completions before execution. Issue prose itself is not approval; no optimization trial or altered-loop-policy retry."}
    context["outer_diagnostic_contract"]["case3_or4_retry"] = True
    context["outer_diagnostic_contract"]["reason"] = "Second attempts of the same original scientific inputs:300s was imported private diagnostic allowance, not an ordinary parity source deadline;600s is finite operational headroom backed by exact stock five-case cold completion observations, not a completion forecast. Original sample1000ms, allocation/refusal/settings/precision/provider/cache/owner gates unchanged."
    context["completion_provenance"] = {"archives": {k: {kk: vv for kk, vv in v.items() if kk != "git_blobs"} for k, v in archive_joins.items()},
         "bindings": key(recovery / "archive-bindings.json"), "published510_inputs_byte_joined": 510,
         "prior_cases3_and4_numerical_verdicts": None, "old_300s_attempts_modified": False,
         "baseline_reader": key(OLD / "run_parity.py"), "private_reader_change": "Only selected ordinal range5..9 to3..4 and index offset5 to3; all original gates/reader logic and600s allowance identical.",
         "owner_and_guard_byte_identical": True,
         "cold_fivecase_whole_wall_seconds": [268.71312268000474, 261.8229029830036, 265.4531771979964, 264.50430519600195, 256.1827670040002],
         "no_causal_speedup_or_case3_4_completion_forecast": True}
    write(WORK / "context.json", context)
    frozen = sorted(set(pins) | {key(p) for p in WORK.rglob("*") if p.is_file()}
                    | {key(WORK / filename) for filename in ("preparation.json", "actions.json")})
    actions = {}
    for case in cases:
        template = copy.deepcopy(original_actions["parity-5"])
        template.update(scope=context["finite_scope"] + " Original case:" + case["name"], frozen_inputs=frozen)
        template["argv"][template["argv"].index("-R") + 1] = "^" + re.escape(case["name"]) + "$"
        require(template["limits"] == LIMITS and template["software_provider"] is True, "original diagnostic/provider controls")
        actions[case["label"]] = template
    write(WORK / "actions.json", actions)
    pins.update({key(p): identity(p) for p in WORK.rglob("*") if p.is_file()})
    context["pins"] = dict(sorted(pins.items()))
    write(WORK / "preparation.json", context)
    require(all(set(a["frozen_inputs"]) == set(pins) | {key(WORK / "preparation.json")} for a in actions.values()), "complete owner frozen closure")
    require(not any(name.startswith("out/") and not name.startswith(key(WORK) + "/") for name in pins), "no other scratch dependency")
    print(json.dumps({"prepared_only": True, "runtime_actions": 0, "pins": len(pins),
         "preparation": identity(WORK / "preparation.json"), "reader": identity(WORK / "run_parity.py"),
         "actions": identity(WORK / "actions.json"), "whole_archive_files": sum(v["whole_raw_files"] for v in archive_joins.values())}, indent=2))


if __name__ == "__main__":
    main()
