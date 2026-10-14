"""Narrow P6 preparation/result checks. No provider is launched here."""
from pathlib import Path
import hashlib
import json
import math
import os
import re
import subprocess
import sys
import xml.etree.ElementTree as ET

sys.dont_write_bytecode = True
WORK = Path(__file__).resolve().parent
ROOT = WORK.parents[1]
CONTEXT = json.loads((WORK / "context.json").read_text())
ACTIONS = json.loads((WORK / "actions.json").read_text())


def identity(path):
    path = Path(path)
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return {"resolved_path": str(path.resolve(strict=True)),
            "bytes": path.stat().st_size, "sha256": h.hexdigest()}


def preflight(case):
    assert case in ACTIONS
    assert Path.cwd() == ROOT
    assert subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip() == CONTEXT["runner_revision"]
    assert subprocess.check_output(["git", "status", "--porcelain"], cwd=ROOT) == b""
    pins = json.loads((WORK / "input-pins.json").read_text())
    for path, expected in pins.items():
        assert identity(ROOT / path) == expected, ("changed preparation input", path)
    for path, expected in CONTEXT["original_inputs"].items():
        assert identity(ROOT / path) == expected, ("changed original input", path)
    assert ACTIONS[case]["revision"] == CONTEXT["runner_revision"]


def validate(case):
    preflight(case)
    folder = WORK / case
    owner = json.loads((folder / "owner.json").read_text())
    action = ACTIONS[case]
    assert owner["status"] == "terminal" and owner["passed"] is True
    assert owner["child_exit"] == 0 and owner["stop_reason"] is None
    assert owner["cleanup_errors"] == [] and owner["remaining_owned_processes"] == []
    assert owner["source"] == {"revision": CONTEXT["runner_revision"], "status": ""}
    assert owner["source_unchanged"] and owner["whole_inputs_unchanged"]
    assert owner["argv"] == action["argv"]
    assert owner["limits"] == action["limits"]
    assert owner["declared_output_environment"] == {
        "GTEST_OUTPUT": "xml:" + str(folder / "gtest.xml") if action["kind"] == "test" else "",
        **action["provider_environment"],
    }
    before = json.loads((folder / "inputs-before.json").read_text())
    after = json.loads((folder / "inputs-after.json").read_text())
    assert before == after
    pins = json.loads((WORK / "input-pins.json").read_text())
    for path, expected in {**CONTEXT["original_inputs"], **pins}.items():
        assert before[path] == expected, ("owner did not consume prepared input", path)
    assert before[str((WORK / "input-pins.json").relative_to(ROOT))] == identity(WORK / "input-pins.json")
    assert owner["observed_births"] and all(b["absent"] is True for b in owner["observed_births"])
    assert owner["boot_id"] == Path("/proc/sys/kernel/random/boot_id").read_text().strip()
    # Current observed births must still be absent; PID reuse is distinguished.
    for birth in owner["observed_births"]:
        p = Path(f'/proc/{birth["pid"]}/stat')
        try:
            ticks = p.read_text().split(")", 1)[1].split()[19]
        except FileNotFoundError:
            continue
        assert ticks != str(birth["saved_start_ticks"]), birth
    mapped = set()
    for observation in owner["provider_mappings"]:
        assert observation["owned_session"] == owner["child_pid"]
        assert any(b["pid"] == observation["pid"] and str(b["saved_start_ticks"]) == str(observation["start_ticks"])
                   for b in owner["observed_births"])
        for line in observation["maps"]:
            parts = line.split(None, 5)
            assert len(parts) == 6 and not parts[5].endswith(" (deleted)")
            path = parts[5]
            assert path in action["provider_modules"]
            stat = Path(path).stat()
            assert int(parts[4]) == stat.st_ino
            assert tuple(int(x, 16) for x in parts[3].split(":")) == (os.major(stat.st_dev), os.minor(stat.st_dev))
            mapped.add(path)
    assert set(action["required_provider_modules"]) <= mapped, ("missing actual provider maps", mapped)
    result = {"case": case, "passed": True, "scope": action["scope"],
              "owner": identity(folder / "owner.json"), "observed_provider_modules": sorted(mapped)}
    if case == "inventory":
        data = json.loads((folder / "stdout.log").read_text())
        vulkan = data["backends"]["vulkan"]
        assert data["wsl2"] is True and vulkan["compiled"] and vulkan["available"]
        assert vulkan["device_count"] == 1 and vulkan["selected_device_index"] == 0
        assert "selection_error" not in vulkan
        device = vulkan["devices"][0]
        for key, expected in CONTEXT["required_device"].items():
            assert device[key] == expected, (key, device)
        assert device["api_version"] >= (1 << 22) + (2 << 12)
        assert vulkan["selected_device_name"] == device["name"]
        result["selected_device"] = device
    else:
        validate("inventory")
        xml = ET.parse(folder / "gtest.xml").getroot()
        for node in [xml, *xml.findall("testsuite")]:
            assert node.attrib["tests"] == "1"
            for key in ["failures", "errors", "disabled", "skipped"]:
                assert int(node.attrib.get(key, "0")) == 0
        cases = xml.findall(".//testcase")
        assert len(cases) == 1
        test = cases[0]
        assert test.attrib["classname"] + "." + test.attrib["name"] == action["filter"]
        assert test.attrib["status"] == "run" and test.attrib["result"] == "completed"
        assert not test.findall("failure") and not test.findall("error") and not test.findall("skipped")
        props = {p.attrib["name"]: p.attrib["value"] for p in test.findall("properties/property")}
        spec = CONTEXT["cases"][case]
        prefix = f'parity_opcode_{spec["opcode"]}_'
        assert props[prefix + "completed_records"] == str(spec["records"])
        assert props[prefix + "actual_submissions"] == "1"
        allocated = int(props[prefix + "allocated_bytes"])
        assert spec["requested_buffer_bytes"] <= allocated <= 1024 * 1024
        timing = float(props[prefix + "maximum_submit_ms"])
        assert math.isfinite(timing) and 0 < timing <= 1000
        stdout = (folder / "stdout.log").read_text()
        dispatches = re.findall(r"\[ParityProbe\] phase=(\w+) opcode=(\d+) offset=(\d+) count=(\d+) submit_wait_ms=([^ ]+) pipeline_setup_ms=([^ ]+) total_ms=([^ ]+) success=(\d+)", stdout)
        # Dozen is an integrated device: original helper performs one actual dispatch.
        assert len(dispatches) == 1, dispatches
        phase, opcode, offset, count, submit, setup, total, success = dispatches[0]
        assert (phase, int(opcode), int(offset), int(count), success) == ("actual", spec["opcode"], 0, spec["records"], "1")
        assert 0 < float(submit) <= 1000 and all(math.isfinite(float(x)) for x in [submit, setup, total])
        assert f'[ParityProbe] completed opcode={spec["opcode"]} records={spec["records"]} submissions=1 allocated_bytes={allocated}' in stdout
        result.update(test=action["filter"], xml=identity(folder / "gtest.xml"), properties=props,
                      meaning="Finite original host-authority/kernel parity on this exact WSL Dozen provider; no full P6, native-platform, scientific-estate or release verdict.")
    return result


if __name__ == "__main__":
    assert len(sys.argv) == 3 and sys.argv[1] in {"--preflight", "--accept"}
    if sys.argv[1] == "--preflight":
        preflight(sys.argv[2])
        print("Prepared input identities match; no executable or provider launched.")
    else:
        result = validate(sys.argv[2])
        target = WORK / ("accepted-" + sys.argv[2] + ".json")
        with target.open("x") as stream:
            json.dump(result, stream, indent=2, sort_keys=True)
            stream.write("\n")
        print(json.dumps(result, indent=2, sort_keys=True))
