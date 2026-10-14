#!/usr/bin/env python3
"""Separate one-shot original CPU writer controls with bounded ownership.

Short controls retain 300 s; Kerr alone has a 900 s outer diagnostic budget.
Original scenes/assertions and the goal's 300 s Vulkan frame criterion remain
unchanged. No image capture hooks: preserve logs/XML and natural TMPDIR residuals.
"""
from pathlib import Path
import argparse
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
VARIANTS = {
    "short": ("original-four-short-p6-cpu-writer-controls", {
        "PNGWriterTest.DecodeRoundTripMatchesSRGBEncoding",
        "EXRRoundTripTests.HDRGradientSurvivesWriteAndRead",
        "EXRRoundTripTests.PpmRgbaBoundaryAppliesExactlyOneSrgbEncode",
        "RenderSessionProbe.FilmAffectsDisplayOutputButNeverLinearExr",
    }, 300),
    "kerr": ("original-kerr-ppm-p6-cpu-writer-control", {
        "RenderSessionProbe.CpuKerrRenderProducesValidPpmThroughTheOwnedWriter",
    }, 900),
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


def write_result(path, value):
    path.write_text(json.dumps(value, indent=2) + "\n")


def temporary_state(folder):
    if not folder.exists():
        return {"path": str(folder), "exists": False}
    require(folder.is_dir() and not folder.is_symlink(), "private TMPDIR is not an ordinary directory")
    status = folder.stat()
    residual = {}
    for path in sorted(folder.rglob("*")):
        require(not path.is_symlink(), "unexpected symlink in private TMPDIR")
        if path.is_file():
            residual[str(path.relative_to(WORK))] = identity(path)
    return {"path": str(folder), "exists": True, "uid": status.st_uid,
            "mode": oct(status.st_mode & 0o777), "residual_files": residual,
            "original_tests_delete_successful_images": True,
            "successful_image_payload_preservation_claimed": False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("selection", choices=["short", "kerr"])
    args = parser.parse_args()
    case_name, names, seconds = VARIANTS[args.selection]
    temporary = WORK / ("temporary-" + args.selection)
    result_path = WORK / (args.selection + "-result.json")
    require(Path.cwd().resolve() == ROOT, "run this action from the Sirius checkout")
    require(sys.flags.optimize == 0 and not os.environ.get("PYTHONOPTIMIZE"),
            "preserved owner's assertions must be enabled")
    require(not (WORK / case_name).exists() and not result_path.exists() and not temporary.exists(),
            "this one-shot action was already attempted; preserve its original observation")
    prepared = json.loads((WORK / "preparation.json").read_text())
    result = {"passed": False, "selection": args.selection,
              "product_source_revision": prepared["product_source_revision"],
              "runner_source_revision": prepared["runner_source_revision"], "expected_names": sorted(names),
              "qualification_claimed": False, "p6_complete": False, "mandatory_receipt_created": False,
              "forensic_pixel_identity_claimed": False, "retained_release_image_qualification_claimed": False,
              "scope": "Original CPU writer/film assertions only; no GPU, full P6 or final integration claim.",
              "finite_acceptance_scope": prepared["finite_acceptance_scope"][args.selection],
              "outer_diagnostic_seconds": seconds, "original_scientific_acceptance_changed": False,
              "vulkan_frame_acceptance_seconds": 300,
              "prior_incomplete_attempt": prepared["prior_incomplete_attempt"],
              "declared_output_environment": {"TMPDIR": str(temporary), "GTEST_OUTPUT": ""}}
    write_result(result_path, result)
    try:
        require(git("rev-parse", "HEAD") == prepared["runner_source_revision"], "runner revision changed")
        require(not git("status", "--porcelain", "--untracked-files=normal"), "tracked source is not clean")
        changed = set(git("diff", "--name-only", prepared["product_source_revision"], "HEAD").splitlines())
        require(changed == set(prepared["diagnostic_only_difference"]), "product/science source difference changed")
        for name, expected in prepared["pins"].items():
            require(identity(ROOT / name) == expected, "prepared input changed: " + name)
        actions = json.loads((WORK / "actions.json").read_text())
        require(set(actions) == {row[0] for row in VARIANTS.values()}, "unexpected action selection")
        require(actions[case_name]["limits"] == {"seconds": seconds, "rss_kib": 4194304, "reserve_kib": 2097152},
                "declared outer diagnostic limits differ")
        available = dict(line.split(":", 1) for line in Path("/proc/meminfo").read_text().splitlines())
        require(int(available["MemAvailable"].split()[0]) >= 2097152,
                "host reserve is below the original diagnostic guard")
        temporary.mkdir(mode=0o700)
        require(temporary.stat().st_uid == os.getuid() and temporary.stat().st_mode & 0o777 == 0o700,
                "private TMPDIR ownership/mode differs")
        os.environ["TMPDIR"] = str(temporary)
        original_argv = sys.argv
        try:
            sys.argv = [str(WORK / "run_owned.py"), case_name]
            runpy.run_path(str(WORK / "run_owned.py"), run_name="__main__")
        finally:
            sys.argv = original_argv
        owner = json.loads((WORK / case_name / "owner.json").read_text())
        require(owner["passed"] is True and owner["child_exit"] == 0
                and owner["stop_reason"] is None and owner["cleanup_errors"] == []
                and owner["remaining_owned_processes"] == []
                and owner["source_unchanged"] is True and owner["whole_inputs_unchanged"] is True
                and owner["observed_births"] and all(row["absent"] for row in owner["observed_births"]),
                "owned completion or input/retirement observation is incomplete")
        xml = WORK / case_name / "ctest.xml"
        document = ET.parse(xml).getroot()
        cases = document.findall(".//testcase")
        require(len(cases) == len(names) and {row.attrib["name"] for row in cases} == names,
                "JUnit does not contain exactly the selected original cases")
        for case in cases:
            require(not any(case.findall(tag) for tag in ("failure", "error", "skipped"))
                    and case.get("status", "run") not in {"notrun", "disabled", "skipped"},
                    "an original case failed, errored, skipped, or was disabled")
        for suite in [document, *document.findall(".//testsuite")]:
            require(all(int(suite.get(key, "0")) == 0 for key in ("failures", "errors", "skipped", "disabled")),
                    "JUnit suite records a non-passing case")
        log = WORK / case_name / "ctest.log"
        require(log.is_file(), "standard CTest output log is missing")
        result.update(passed=True, executed_names=sorted(names), failures=0, errors=0, skips=0, disabled=0,
                      junit=identity(xml), ctest_output_log=identity(log), owner=identity(WORK / case_name / "owner.json"),
                      reference=prepared["reference"], source_and_registration_join=prepared["source_and_registration_join"],
                      temporary_output=temporary_state(temporary))
        write_result(result_path, result)
        print(json.dumps(result, indent=2), flush=True)
        return 0
    except BaseException as error:
        result["error"] = repr(error)
        try:
            result["temporary_output"] = temporary_state(temporary)
        except BaseException as output_error:
            result["temporary_output_error"] = repr(output_error)
        write_result(result_path, result)
        raise


if __name__ == "__main__":
    raise SystemExit(main())
