#!/usr/bin/env python3
"""One unchanged original P4 Vulkan scene case under the preserved owner."""
from pathlib import Path
import hashlib
import json
import math
import os
import runpy
import struct
import subprocess
import sys
import xml.etree.ElementTree as ET

sys.dont_write_bytecode = True
WORK = Path(__file__).resolve().parent
ROOT = WORK.parents[1]
LIMITS = {"seconds": 300, "rss_kib": 4194304, "reserve_kib": 2097152}


def require(ok, message):
    if not ok:
        raise ValueError(message)


def identity(path):
    path = Path(path)
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return {"bytes": path.stat().st_size, "sha256": digest.hexdigest(), "resolved_path": str(path.resolve(strict=True))}


def git(*args):
    return subprocess.check_output(["git", *args], cwd=ROOT, timeout=30).decode().strip()


def exr_header(path):
    # Header-only transport integrity. Original in-memory GTest assertions remain
    # the numerical authority; this does not decode or replace those assertions.
    raw = path.read_bytes()
    require(len(raw) > 8 and struct.unpack_from("<I", raw)[0] == 20000630, "EXR magic/length")
    at, attributes = 8, {}
    while at < min(len(raw), 65536) and raw[at]:
        end = raw.index(0, at, min(len(raw), 65536))
        name = raw[at:end].decode("ascii")
        at = end + 1
        end = raw.index(0, at, min(len(raw), 65536))
        kind = raw[at:end].decode("ascii")
        at = end + 1
        size = struct.unpack_from("<I", raw, at)[0]
        at += 4
        require(at + size <= min(len(raw), 65536) and name not in attributes, "EXR header bounds/duplicate")
        attributes[name] = (kind, raw[at:at + size])
        at += size
    require(at < min(len(raw), 65536) and raw[at] == 0, "EXR header terminator")
    require(attributes.get("dataWindow") == ("box2i", struct.pack("<4i", 0, 0, 39, 23)), "original EXR dimensions")
    require(attributes.get("siriusSamplesPerPixel") == ("string", b"1"), "original EXR sample metadata")
    return {"width": 40, "height": 24, "samples_per_pixel": 1, "identity": identity(path)}


def main():
    require(Path.cwd().resolve() == ROOT and len(sys.argv) == 2 and sys.argv[1] in ("1", "2"), "run from checkout with ordinal1 or2")
    require(sys.flags.optimize == 0 and not os.environ.get("PYTHONOPTIMIZE"), "owner assertions must be enabled")
    prepared = json.loads((WORK / "preparation.json").read_text())
    index = int(sys.argv[1]) - 1
    case = prepared["cases"][index]
    label = case["label"]
    folder = WORK / label
    temporary = WORK / (label + "-temporary")
    result_path = WORK / (label + "-result.json")
    require(not any(p.exists() for p in (folder, temporary, result_path)), "one-shot action already attempted")
    for prior in prepared["cases"][:index]:
        result = json.loads((WORK / (prior["label"] + "-result.json")).read_text())
        require(result.get("passed") is True and result["name"] == prior["name"], "preceding original scene not accepted")
        for filename, field in (("owner.json", "owner"), ("gtest.xml", "gtest"), ("ctest.xml", "ctest")):
            require(identity(WORK / prior["label"] / filename) == result[field], "preceding original acceptance changed")
        for name, pin in result["output_files"].items():
            require(identity(WORK / name) == pin, "preceding original image changed")
    result = {"passed": False, "name": case["name"], "original_assertions_changed": False,
              "original_precision_request": "unset, original defaultFP32", "outer_diagnostic_limits": LIMITS,
              "qualification_claimed": False, "p4_complete": False}
    def save():
        result_path.write_text(json.dumps(result, indent=2) + "\n")
    save()
    try:
        require(git("rev-parse", "HEAD") == prepared["runner_source_revision"] and not git("status", "--porcelain"), "runner source changed")
        require(set(git("diff", "--name-only", prepared["product_source_revision"], "HEAD").splitlines()) == set(prepared["diagnostic_only_difference"]), "product difference changed")
        require("SIRIUS_PRECISION" not in os.environ, "retain original unset precision request")
        for name, pin in prepared["pins"].items():
            require(identity(ROOT / name) == pin, "prepared input changed:" + name)
        action = json.loads((WORK / "actions.json").read_text())[label]
        require(action["limits"] == LIMITS and action["software_provider"] is True, "diagnostic limits/provider changed")
        temporary.mkdir(mode=0o700)
        require(temporary.stat().st_uid == os.getuid() and temporary.stat().st_mode & 0o777 == 0o700, "TMPDIR identity")
        os.environ["TMPDIR"] = str(temporary)
        old_argv = sys.argv
        try:
            sys.argv = [str(WORK / "run_owned.py"), label]
            runpy.run_path(str(WORK / "run_owned.py"), run_name="__main__")
        finally:
            sys.argv = old_argv
        owner = json.loads((folder / "owner.json").read_text())
        require(owner["passed"] is True and owner["status"] == "terminal" and owner["child_exit"] == 0
                and owner["stop_reason"] is None and owner["cleanup_errors"] == []
                and owner["source_unchanged"] and owner["whole_inputs_unchanged"]
                and owner["remaining_owned_processes"] == [] and owner["observed_births"]
                and all(row["absent"] for row in owner["observed_births"]), "owner not accepted")
        mappings = {}
        for observation in owner["provider_mappings"]:
            require(observation["owned_session"] == owner["child_pid"], "mapping ownership changed")
            for line in observation["maps"]:
                fields = line.split(None, 5)
                require(len(fields) == 6 and fields[5] in prepared["required_provider_mappings"], "unexpected provider mapping")
                path = fields[5]
                info = Path(path).stat()
                require(tuple(int(x, 16) for x in fields[3].split(":")) == (os.major(info.st_dev), os.minor(info.st_dev)) and int(fields[4]) == info.st_ino, "mapped module device/inode changed")
                mappings[path] = identity(path)
        require(set(mappings) == set(prepared["required_provider_mappings"]), "actual loader/LVP/LLVM mapping incomplete")
        require(all(pin == prepared["pins"][path] for path, pin in mappings.items()), "mapped provider bytes changed")
        gtest = ET.parse(folder / "gtest.xml").getroot()
        rows = gtest.findall(".//testcase")
        require(len(rows) == 1 and rows[0].get("classname") + "." + rows[0].get("name") == case["name"]
                and rows[0].get("status") == "run" and rows[0].get("result") == "completed", "original GTest selection/completion differs")
        ctest = ET.parse(folder / "ctest.xml").getroot()
        rows = ctest.findall(".//testcase")
        require(len(rows) == 1 and rows[0].get("name") == case["name"], "registered CTest selection differs")
        for tree in (gtest, ctest):
            for suite in [tree, *tree.findall(".//testsuite")]:
                require(all(int(suite.get(k, "0")) == 0 for k in ("failures", "errors", "skipped", "disabled")), "XML has nonpass")
            for row in tree.findall(".//testcase"):
                require(not any(row.findall(tag) for tag in ("failure", "error", "skipped")), "original assertion/skip failure")
        prefix = "[Vulkan] Render evidence: "
        frames = [json.loads(line.split(prefix, 1)[1]) for line in (folder / "stdout.log").read_text().splitlines() if prefix in line]
        require(len(frames) == len(case["outputs"]), "original complete frame count differs")
        for frame in frames:
            require(frame["schema"] == "sirius-vulkan-render-v1" and frame["route"] == "retained"
                    and frame["source_owner"] == "host" and frame["width"] == 40 and frame["height"] == 24
                    and frame["precision"] == "fp32" and frame["device_index"] == 0
                    and "llvmpipe" in frame["device_name"].lower(), "actual render route/device/scene differs")
            require(frame["dispatches"] > 0 and frame["accepted_intervals"] > 0, "no physical retained work")
            require(all(math.isfinite(frame[k]) and frame[k] >= 0 for k in ("wall_seconds", "dispatch_seconds", "maximum_dispatch_ms", "initialization_seconds", "initialization_submit_wait_ms")), "nonfinite phase record")
            preparation = frame["retained_preparation"]["stages"]
            active = [row for row in preparation if row["attempts"]]
            require({row["stage"] for row in active} == set(prepared["software_preparation_stages"]), "software preparation stages differ")
            require(all(row["attempts"] == row["dispatch_attempts"] == row["completed_dispatches"] == row["completed"] == 1 and row["header_restored"] is True for row in active), "zero-active preparation incomplete")
            require(frame["initialization_dispatches"] == 6, "software preparation dispatch count differs")
        entries = list(temporary.iterdir())
        require({p.name for p in entries} == set(case["outputs"]) and all(p.is_file() and not p.is_symlink() for p in entries), "owned image output set differs")
        outputs = {str(p.relative_to(WORK)): identity(p) for p in sorted(entries)}
        headers = {p.name: exr_header(p) for p in sorted(entries)}
        result.update(passed=True, owner=identity(folder / "owner.json"), gtest=identity(folder / "gtest.xml"),
                      ctest=identity(folder / "ctest.xml"), output_files=outputs, exr_headers=headers,
                      provider_mappings=mappings, complete_original_frames=frames,
                      finite_scope=prepared["finite_scope"], no_performance_or_cache_warmness_claim=True)
        save()
        print(json.dumps(result, indent=2), flush=True)
        return 0
    except BaseException as error:
        result["error"] = repr(error)
        save()
        raise


if __name__ == "__main__":
    raise SystemExit(main())
