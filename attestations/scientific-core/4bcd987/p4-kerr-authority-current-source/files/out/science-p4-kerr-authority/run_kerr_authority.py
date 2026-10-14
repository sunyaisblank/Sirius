#!/usr/bin/env python3
"""One frozen four-case original P4 Kerr authority action; no build/device/frame work."""
from pathlib import Path
import hashlib
import json
import os
import runpy
import subprocess
import sys
import xml.etree.ElementTree as ET

sys.dont_write_bytecode = True
WORK = Path(__file__).resolve().parent
ROOT = WORK.parents[1]
CASE = "original-four-kerr-authority-controls"
NAMES = {'KerrOrbitAuthority.UnrepresentedInputsDeclineInsteadOfChangingTheOrbit', 'KerrOrbitAuthority.EveryCpuConsumerSharesSignedIscoAndCircularEmitterLaw', 'KerrZamoTransfer.KerrSchildSlicingNormalIsNotTheOffEquatorialZamo', 'KerrZamoTransfer.ArbitraryLatitudeMatchesIndependentBoyerLindquistContraction'}
LIMITS = {"seconds": 300, "rss_kib": 4194304, "reserve_kib": 2097152}


def require(ok, message):
    if not ok:
        raise ValueError(message)


def identity(path):
    raw = Path(path).read_bytes()
    return {"bytes": len(raw), "sha256": hashlib.sha256(raw).hexdigest()}


def git(*args):
    return subprocess.check_output(["git", *args], cwd=ROOT, timeout=30).decode().strip()


def write_result(value):
    (WORK / "result.json").write_text(json.dumps(value, indent=2) + "\n")


def temporary_state(folder):
    if not folder.exists():
        return {"exists": False}
    require(folder.is_dir() and not folder.is_symlink(), "TMPDIR is not an ordinary directory")
    residual = {}
    for path in sorted(folder.rglob("*")):
        require(not path.is_symlink(), "unexpected TMPDIR symlink")
        if path.is_file():
            residual[str(path.relative_to(WORK))] = identity(path)
    return {"path": str(folder), "uid": folder.stat().st_uid,
            "mode": oct(folder.stat().st_mode & 0o777), "residual_files": residual}


def main():
    require(Path.cwd().resolve() == ROOT, "run from the Sirius checkout")
    require(sys.flags.optimize == 0 and not os.environ.get("PYTHONOPTIMIZE"), "owner assertions must be enabled")
    temporary = WORK / "temporary"
    require(not (WORK / CASE).exists() and not (WORK / "result.json").exists() and not temporary.exists(),
            "one-shot action already attempted; retain its original observation")
    prepared = json.loads((WORK / "preparation.json").read_text())
    result = {"passed": False, "product_source_revision": prepared["product_source_revision"],
              "runner_source_revision": prepared["runner_source_revision"], "expected_names": sorted(NAMES),
              "scope": prepared["scope"], "finite_acceptance_scope": prepared["finite_acceptance_scope"],
              "qualification_claimed": False, "p4_complete": False, "mandatory_receipt_created": False,
              "original_assertions_changed": False, "outer_diagnostic_limits": LIMITS,
              "declared_output_environment": {"TMPDIR": str(temporary), "GTEST_OUTPUT": ""}}
    write_result(result)
    try:
        require(git("rev-parse", "HEAD") == prepared["runner_source_revision"] and not git("status", "--porcelain"),
                "runner source is not the pinned clean revision")
        require(set(git("diff", "--name-only", prepared["product_source_revision"], "HEAD").splitlines())
                == set(prepared["diagnostic_only_difference"]), "product/science source difference changed")
        require(set(prepared["expected_names"]) == NAMES and len(NAMES) == 4, "original model selection differs")
        for name, pin in prepared["pins"].items():
            require(identity(ROOT / name) == pin, "prepared input changed: " + name)
        actions = json.loads((WORK / "actions.json").read_text())
        require(set(actions) == {CASE} and actions[CASE]["limits"] == LIMITS, "owner action/limits differ")
        available = dict(line.split(":", 1) for line in Path("/proc/meminfo").read_text().splitlines())
        require(int(available["MemAvailable"].split()[0]) >= LIMITS["reserve_kib"], "host reserve is below guard")
        temporary.mkdir(mode=0o700)
        require(temporary.stat().st_uid == os.getuid() and temporary.stat().st_mode & 0o777 == 0o700,
                "private TMPDIR ownership/mode differs")
        os.environ["TMPDIR"] = str(temporary)
        original_argv = sys.argv
        try:
            sys.argv = [str(WORK / "run_owned.py"), CASE]
            runpy.run_path(str(WORK / "run_owned.py"), run_name="__main__")
        finally:
            sys.argv = original_argv
        owner = json.loads((WORK / CASE / "owner.json").read_text())
        require(owner["passed"] is True and owner["child_exit"] == 0 and owner["stop_reason"] is None
                and owner["cleanup_errors"] == [] and owner["remaining_owned_processes"] == []
                and owner["source_unchanged"] and owner["whole_inputs_unchanged"]
                and owner["observed_births"] and all(row["absent"] for row in owner["observed_births"]),
                "owned completion/input/retirement observation is incomplete")
        xml = WORK / CASE / "ctest.xml"
        document = ET.parse(xml).getroot()
        cases = document.findall(".//testcase")
        require(len(cases) == 4 and {row.attrib["name"] for row in cases} == NAMES, "JUnit case set differs")
        for row in cases:
            require(row.get("status", "run") == "run" and not any(row.findall(tag) for tag in ("failure", "error", "skipped")),
                    "an original case failed, errored, skipped, or was disabled")
        for suite in [document, *document.findall(".//testsuite")]:
            require(all(int(suite.get(key, "0")) == 0 for key in ("failures", "errors", "skipped", "disabled")),
                    "JUnit suite records a non-passing case")
        log = WORK / CASE / "ctest.log"
        require(log.is_file(), "standard CTest output log is missing")
        result.update(passed=True, executed_names=sorted(NAMES), failures=0, errors=0, skips=0, disabled=0,
                      junit=identity(xml), ctest_output_log=identity(log), owner=identity(WORK / CASE / "owner.json"),
                      source_and_registration_join=prepared["source_and_registration_join"],
                      reference=prepared["reference"], temporary_output=temporary_state(temporary))
        write_result(result)
        print(json.dumps(result, indent=2), flush=True)
        return 0
    except BaseException as error:
        result["error"] = repr(error)
        try:
            result["temporary_output"] = temporary_state(temporary)
        except BaseException as output_error:
            result["temporary_output_error"] = repr(output_error)
        write_result(result)
        raise


if __name__ == "__main__":
    raise SystemExit(main())
