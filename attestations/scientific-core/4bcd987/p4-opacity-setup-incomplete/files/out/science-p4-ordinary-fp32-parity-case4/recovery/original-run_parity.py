#!/usr/bin/env python3
"""One original ordinary-FP32 P4 parity case; root operates sequentially."""
from pathlib import Path
import hashlib
import json
import math
import os
import re
import runpy
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
    return {"bytes": path.stat().st_size, "sha256": digest.hexdigest(),
            "resolved_path": str(path.resolve(strict=True))}


def git(*args):
    return subprocess.check_output(["git", *args], cwd=ROOT, timeout=30).decode().strip()


def main():
    require(Path.cwd().resolve() == ROOT, "run from checkout")
    require(len(sys.argv) == 2 and sys.argv[1] in [str(i) for i in range(1, 10)], "select one ordinal 1..9")
    require(sys.flags.optimize == 0 and not os.environ.get("PYTHONOPTIMIZE"), "assertions must be enabled")
    prepared = json.loads((WORK / "preparation.json").read_text())
    index = int(sys.argv[1]) - 1
    case = prepared["cases"][index]
    label = case["label"]
    result_path = WORK / (label + "-result.json")
    folder = WORK / label
    temporary = WORK / (label + "-temporary")
    require(not any(p.exists() for p in (result_path, folder, temporary)), "one-shot case already attempted")
    for prior in prepared["cases"][:index]:
        row = json.loads((WORK / (prior["label"] + "-result.json")).read_text())
        require(row.get("passed") is True and row["name"] == prior["name"], "preceding original case not accepted")
        for suffix, key in (("owner.json", "owner"), ("gtest.xml", "gtest"), ("ctest.xml", "ctest")):
            require(identity(WORK / prior["label"] / suffix) == row[key], "preceding raw acceptance changed")
    result = {"passed": False, "name": case["name"], "ordinary_fp32": True,
              "original_assertions_changed": False, "qualification_claimed": False,
              "p4_complete": False, "outer_diagnostic_limits": LIMITS}
    def save():
        result_path.write_text(json.dumps(result, indent=2) + "\n")
    save()
    try:
        require(git("rev-parse", "HEAD") == prepared["runner_source_revision"] and not git("status", "--porcelain"),
                "source differs from pinned clean runner")
        require(set(git("diff", "--name-only", prepared["product_source_revision"], "HEAD").splitlines())
                == set(prepared["diagnostic_only_difference"]), "product difference changed")
        for name, pin in prepared["pins"].items():
            require(identity(ROOT / name) == pin, "frozen input changed: " + name)
        actions = json.loads((WORK / "actions.json").read_text())
        action = actions[label]
        require(action["limits"] == LIMITS and action["software_provider"] is True, "limits/provider policy changed")
        temporary.mkdir(mode=0o700)
        require(temporary.stat().st_uid == os.getuid() and temporary.stat().st_mode & 0o777 == 0o700, "TMPDIR identity")
        os.environ["TMPDIR"] = str(temporary)
        original_argv = sys.argv
        try:
            sys.argv = [str(WORK / "run_owned.py"), label]
            runpy.run_path(str(WORK / "run_owned.py"), run_name="__main__")
        finally:
            sys.argv = original_argv
        owner = json.loads((folder / "owner.json").read_text())
        require(owner["passed"] is True and owner["status"] == "terminal" and owner["child_exit"] == 0
                and owner["stop_reason"] is None and owner["cleanup_errors"] == []
                and owner["remaining_owned_processes"] == [] and owner["source_unchanged"]
                and owner["whole_inputs_unchanged"] and owner["observed_births"]
                and all(row["absent"] for row in owner["observed_births"]), "owner acceptance incomplete")
        mapped = {}
        for observation in owner["provider_mappings"]:
            require(observation["owned_session"] == owner["child_pid"], "mapping ownership differs")
            for line in observation["maps"]:
                fields = line.split(None, 5)
                require(len(fields) == 6 and not fields[5].endswith(" (deleted)"), "malformed/deleted provider mapping")
                path = fields[5]
                require(path in prepared["required_provider_mappings"], "unexpected collected provider")
                info = Path(path).stat()
                dev = tuple(int(piece, 16) for piece in fields[3].split(":"))
                require(dev == (os.major(info.st_dev), os.minor(info.st_dev)) and int(fields[4]) == info.st_ino,
                        "provider mapping device/inode differs: " + path)
                mapped[path] = identity(path)
        require(set(mapped) == set(prepared["required_provider_mappings"]), "provider/loader/LLVM maps incomplete")
        for path, pin in mapped.items():
            require(pin == prepared["pins"][path], "mapped module bytes differ")
        document = ET.parse(folder / "gtest.xml").getroot()
        cases = document.findall(".//testcase")
        require(len(cases) == 1 and cases[0].get("classname") + "." + cases[0].get("name") == case["name"], "GTest selection differs")
        require(cases[0].get("status") == "run" and cases[0].get("result") == "completed", "GTest did not complete")
        for tree in (document, ET.parse(folder / "ctest.xml").getroot()):
            for suite in [tree, *tree.findall(".//testsuite")]:
                require(all(int(suite.get(k, "0")) == 0 for k in ("failures", "errors", "skipped", "disabled")), "XML nonpass")
            for row in tree.findall(".//testcase"):
                require(not any(row.findall(tag) for tag in ("failure", "error", "skipped")), "original assertion/skip failure")
        ctest_cases = ET.parse(folder / "ctest.xml").getroot().findall(".//testcase")
        require(len(ctest_cases) == 1 and ctest_cases[0].get("name") == case["name"], "CTest registration selection differs")
        properties = {row.get("name"): row.get("value") for row in cases[0].findall("./properties/property")}
        prefix = "parity_opcode_" + str(case["opcode"])
        require(properties.get(prefix + "_completed_records") == str(case["records"]), "original record count differs")
        require(properties.get(prefix + "_actual_submissions") == "1", "original submission count differs")
        allocation_text = properties.get(prefix + "_allocated_bytes", "")
        require(re.fullmatch(r"[0-9]+", allocation_text) is not None, "malformed original allocation")
        allocation = int(allocation_text)
        require(case["minimum_allocated_bytes"] <= allocation <= case["maximum_allocated_bytes"],
                "original allocation bounds differ")
        maximum = float(properties[prefix + "_maximum_submit_ms"])
        require(math.isfinite(maximum) and 0 < maximum <= 1000, "original actual-dispatch bound fails")
        stdout = (folder / "stdout.log").read_text()
        phases = re.findall(r"\[ParityProbe\] phase=(\w+) opcode=(\d+) offset=(\d+) count=(\d+) submit_wait_ms=([^ ]+) pipeline_setup_ms=([^ ]+) total_ms=([^ ]+) success=(\d+)", stdout)
        require(len(phases) == 2, "software initialization/actual dispatch log incomplete")
        for row, phase, count in zip(phases, ("initialization", "actual"), (0, case["records"])):
            require(row[:4] == (phase, str(case["opcode"]), "0", str(count)) and row[-1] == "1", "dispatch provenance differs")
            timing = [float(x) for x in row[4:7]]
            require(all(math.isfinite(x) and x >= 0 for x in timing) and timing[0] > 0, "malformed timing")
        require(float(phases[1][4]) <= 1000, "actual sample work exceeded original hard limit")
        result.update(passed=True, owner=identity(folder / "owner.json"), gtest=identity(folder / "gtest.xml"),
                      ctest=identity(folder / "ctest.xml"), provider_mappings=mapped, original_properties=properties,
                      original_dispatch_phases=phases, reference=case["reference"],
                      provider_scope=prepared["provider_policy"], no_performance_or_cache_warmness_claim=True)
        save()
        print(json.dumps(result, indent=2), flush=True)
        return 0
    except BaseException as error:
        result["error"] = repr(error)
        save()
        raise


if __name__ == "__main__":
    raise SystemExit(main())
