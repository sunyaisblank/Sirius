#!/usr/bin/env python3
"""One frozen five-case CPU writer/film action under the preserved owner.

Original assertions, fixtures and cleanup remain intact. Successful images are
deleted by the tests; preserve logs/XML and natural residuals, not forensic tapes.
The 300 s / 4 GiB RSS / 2 GiB reserve limits are outer diagnostic guards.
"""
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
TEMP = WORK / "temporary"
CASE = "original-five-p6-cpu-writer-controls"
NAMES = {
    "PNGWriterTest.DecodeRoundTripMatchesSRGBEncoding",
    "EXRRoundTripTests.HDRGradientSurvivesWriteAndRead",
    "EXRRoundTripTests.PpmRgbaBoundaryAppliesExactlyOneSrgbEncode",
    "RenderSessionProbe.CpuKerrRenderProducesValidPpmThroughTheOwnedWriter",
    "RenderSessionProbe.FilmAffectsDisplayOutputButNeverLinearExr",
}


def require(ok, message):
    if not ok:
        raise ValueError(message)


def identity(path):
    path = Path(path)
    return {"bytes": path.stat().st_size,
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def git(*args):
    return subprocess.check_output(["git", *args], cwd=ROOT, timeout=30).decode().strip()


def write_result(value):
    (WORK / "result.json").write_text(json.dumps(value, indent=2) + "\n")


def temporary_state():
    if not TEMP.exists():
        return {"path": str(TEMP), "exists": False}
    require(TEMP.is_dir() and not TEMP.is_symlink(), "private TMPDIR is not an ordinary directory")
    status = TEMP.stat()
    residual = {}
    for path in sorted(TEMP.rglob("*")):
        require(not path.is_symlink(), "unexpected symlink in private TMPDIR")
        if path.is_file():
            residual[str(path.relative_to(WORK))] = identity(path)
    return {"path": str(TEMP), "exists": True, "uid": status.st_uid,
            "mode": oct(status.st_mode & 0o777), "residual_files": residual,
            "original_tests_delete_successful_images": True,
            "successful_image_payload_preservation_claimed": False}


def main():
    require(Path.cwd().resolve() == ROOT, "run this action from the Sirius checkout")
    require(sys.flags.optimize == 0 and not os.environ.get("PYTHONOPTIMIZE"),
            "preserved owner's assertions must be enabled")
    require(not (WORK / CASE).exists() and not (WORK / "result.json").exists() and not TEMP.exists(),
            "one-shot action already attempted; preserve the original observation")
    prepared = json.loads((WORK / "preparation.json").read_text())
    result = {"passed": False, "product_source_revision": prepared["product_source_revision"],
              "runner_source_revision": prepared["runner_source_revision"],
              "expected_names": sorted(NAMES), "qualification_claimed": False,
              "p6_complete": False, "mandatory_receipt_created": False,
              "forensic_pixel_identity_claimed": False,
              "retained_release_image_qualification_claimed": False,
              "scope": "Five original CPU writer/film controls only; no GPU dispatch, full P6 or final integration claim.",
              "finite_acceptance_scope": prepared["finite_acceptance_scope"],
              "declared_output_environment": {"TMPDIR": str(TEMP), "GTEST_OUTPUT": ""}}
    write_result(result)
    try:
        require(git("rev-parse", "HEAD") == prepared["runner_source_revision"], "runner revision changed")
        require(not git("status", "--porcelain", "--untracked-files=normal"), "tracked source is not clean")
        changed = set(git("diff", "--name-only", prepared["product_source_revision"], "HEAD").splitlines())
        require(changed == set(prepared["diagnostic_only_difference"]), "product/science source difference changed")
        for name, expected in prepared["pins"].items():
            require(identity(ROOT / name) == expected, "prepared input changed: " + name)
        actions = json.loads((WORK / "actions.json").read_text())
        require(set(actions) == {CASE}, "unexpected action selection")
        available = dict(line.split(":", 1) for line in Path("/proc/meminfo").read_text().splitlines())
        require(int(available["MemAvailable"].split()[0]) >= actions[CASE]["limits"]["reserve_kib"],
                "host reserve is below the original diagnostic guard")
        TEMP.mkdir(mode=0o700)
        require(TEMP.stat().st_uid == os.getuid() and TEMP.stat().st_mode & 0o777 == 0o700,
                "private TMPDIR ownership/mode differs")
        os.environ["TMPDIR"] = str(TEMP)
        original_argv = sys.argv
        try:
            sys.argv = [str(WORK / "run_owned.py"), CASE]
            runpy.run_path(str(WORK / "run_owned.py"), run_name="__main__")
        finally:
            sys.argv = original_argv
        owner = json.loads((WORK / CASE / "owner.json").read_text())
        require(owner["passed"] is True and owner["child_exit"] == 0
                and owner["stop_reason"] is None and owner["cleanup_errors"] == []
                and owner["remaining_owned_processes"] == []
                and owner["source_unchanged"] is True and owner["whole_inputs_unchanged"] is True
                and owner["observed_births"] and all(x["absent"] for x in owner["observed_births"]),
                "owned completion or input/retirement observation is incomplete")
        xml = WORK / CASE / "ctest.xml"
        document = ET.parse(xml).getroot()
        cases = document.findall(".//testcase")
        require(len(cases) == 5 and {x.attrib["name"] for x in cases} == NAMES,
                "JUnit does not contain exactly the five original cases")
        for case in cases:
            require(not any(case.findall(tag) for tag in ("failure", "error", "skipped"))
                    and case.get("status", "run") not in {"notrun", "disabled", "skipped"},
                    "an original case failed, errored, skipped, or was disabled")
        for suite in [document, *document.findall(".//testsuite")]:
            require(all(int(suite.get(key, "0")) == 0 for key in ("failures", "errors", "skipped", "disabled")),
                    "JUnit suite records a non-passing case")
        result.update(passed=True, executed_names=sorted(NAMES), failures=0, errors=0, skips=0, disabled=0,
                      junit=identity(xml), owner=identity(WORK / CASE / "owner.json"),
                      reference=prepared["reference"], source_and_registration_join=prepared["source_and_registration_join"],
                      temporary_output=temporary_state())
        write_result(result)
        print(json.dumps(result, indent=2), flush=True)
        return 0
    except BaseException as error:
        result["error"] = repr(error)
        try:
            result["temporary_output"] = temporary_state()
        except BaseException as output_error:
            result["temporary_output_error"] = repr(output_error)
        write_result(result)
        raise


if __name__ == "__main__":
    raise SystemExit(main())
