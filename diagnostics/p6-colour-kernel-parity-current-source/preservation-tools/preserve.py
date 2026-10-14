"""Preserve only the completed, original two-case P6 parity diagnostic."""
from pathlib import Path
import hashlib
import json
import os
import shutil
import stat
import subprocess
import time
import xml.etree.ElementTree as ET
import zlib

ROOT = Path(__file__).resolve().parents[2]
WORK = ROOT / "out/p6-colour-kernel-parity-current-source"
PREP = Path(__file__).resolve().parent
ARCHIVE = ROOT / "attestations/scientific-core/4bcd987/p6-colour-kernel-parity-current-source"
REVISION = "a1d6e75299a2bf250e5ce7d44ad5965f8acb4299"


def identity(path):
    path = Path(path)
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return {"resolved_path": str(path.resolve(strict=True)), "bytes": path.stat().st_size,
            "sha256": h.hexdigest()}


def write(path, value):
    with path.open("x") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")


def ticks(pid):
    try:
        return Path(f"/proc/{pid}/stat").read_text().split(")", 1)[1].split()[19]
    except FileNotFoundError:
        return None


def tree():
    files, directories = {}, {}
    for path in [WORK, *sorted(WORK.rglob("*"))]:
        mode = path.lstat().st_mode
        assert not stat.S_ISLNK(mode), path
        relative = str(path.relative_to(ROOT))
        if stat.S_ISDIR(mode):
            directories[relative] = stat.S_IMODE(mode)
        else:
            assert stat.S_ISREG(mode), path
            files[relative] = {**identity(path), "mode": stat.S_IMODE(mode)}
    return files, directories


def main():
    assert Path.cwd() == ROOT
    assert subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip() == REVISION
    assert subprocess.check_output(["git", "status", "--porcelain"]) == b""
    assert not ARCHIVE.exists() and not (PREP / "preservation-result.json").exists()
    boot = Path("/proc/sys/kernel/random/boot_id").read_text().strip()
    expected, births, summaries = {}, {}, {}
    paths = [WORK / "inventory", WORK / "shader-phase/cie", WORK / "shader-phase/xyz"]
    for folder in paths:
        owner = json.loads((folder / "owner.json").read_text())
        assert owner["status"] == "terminal" and owner["passed"] is True
        assert owner["child_exit"] == 0 and owner["stop_reason"] is None
        assert owner["cleanup_errors"] == [] and owner["remaining_owned_processes"] == []
        assert owner["source"] == {"revision": REVISION, "status": ""}
        assert owner["source_unchanged"] and owner["whole_inputs_unchanged"]
        assert owner["boot_id"] == boot
        before = json.loads((folder / "inputs-before.json").read_text())
        after = json.loads((folder / "inputs-after.json").read_text())
        assert before == after
        for path, value in before.items():
            assert path not in expected or expected[path] == value, path
            expected[path] = value
        for birth in [*owner["observed_births"], {
            "pid": owner["controller"]["pid"],
            "saved_start_ticks": owner["controller"]["start_ticks"],
        }]:
            saved = str(birth["saved_start_ticks"])
            current = ticks(birth["pid"])
            assert current != saved, birth
            births[(birth["pid"], saved)] = {"pid": birth["pid"], "saved_start_ticks": saved,
                                             "current_start_ticks": current, "absent": True}
        summaries[folder.name] = {key: owner[key] for key in [
            "passed", "child_exit", "wall_seconds", "peak_sampled_rss_kib", "limits"]}
        summaries[folder.name]["owner"] = identity(folder / "owner.json")
    for path, value in expected.items():
        assert identity(ROOT / path) == value, ("changed consumed input", path)
    analysis_path = WORK / "post-owner-overlay-review/analysis-inventory.json"
    analysis = json.loads(analysis_path.read_text())
    assert identity(analysis_path)["sha256"] == "0718dce91d3f115fd4775afc2625afb08baa864fee5e1ca6236ff73410f141fb"
    assert analysis["passed"] and analysis["phase"] == "separate_post_owner_overlay_analysis"
    assert analysis["original_checker_reported_failure_preserved"] and not analysis["qualification_claimed"]
    assert not (WORK / "accepted-inventory.json").exists()
    cases = {}
    for case, records, opcode in [("cie", 10, 14), ("xyz", 5, 13)]:
        folder = WORK / "shader-phase" / case
        accepted = json.loads((folder.parent / ("accepted-" + case + ".json")).read_text())
        assert accepted["passed"] and not accepted["qualification_claimed"]
        assert accepted["owner"] == identity(folder / "owner.json")
        xml = ET.parse(folder / "gtest.xml").getroot()
        tests = xml.findall(".//testcase")
        assert len(tests) == 1 and not xml.findall(".//failure") and not xml.findall(".//error") and not xml.findall(".//skipped")
        test = tests[0]
        assert test.attrib["status"] == "run" and test.attrib["result"] == "completed"
        assert test.attrib["classname"] + "." + test.attrib["name"] == accepted["test"]
        properties = {p.attrib["name"]: p.attrib["value"] for p in test.findall("properties/property")}
        assert properties == accepted["properties"]
        assert properties[f"parity_opcode_{opcode}_completed_records"] == str(records)
        assert properties[f"parity_opcode_{opcode}_actual_submissions"] == "1"
        assert 0 < float(properties[f"parity_opcode_{opcode}_maximum_submit_ms"]) <= 1000
        cases[case] = {"test": accepted["test"], "properties": properties,
                       "accepted_result": identity(folder.parent / ("accepted-" + case + ".json")),
                       "xml": identity(folder / "gtest.xml")}
    files, directories = tree()
    ARCHIVE.mkdir(parents=True)
    for relative, mode in directories.items():
        destination = ARCHIVE / "files" / relative
        destination.mkdir(parents=True, exist_ok=True)
        destination.chmod(mode)
    relocations = {}
    for relative, value in files.items():
        source = ROOT / relative
        destination = ARCHIVE / "files" / relative
        shutil.copy2(source, destination)
        archived = identity(destination)
        assert (archived["bytes"], archived["sha256"]) == (value["bytes"], value["sha256"])
        assert stat.S_IMODE(destination.stat().st_mode) == value["mode"]
        relocations[relative] = {"origin": value, "archive_relative": "files/" + relative,
                                 "archived": archived}
    assert tree() == (files, directories), "live parity root changed during preservation"
    for path, value in expected.items():
        assert identity(ROOT / path) == value, ("changed consumed input after copy", path)
    for item in births.values():
        current = ticks(item["pid"])
        assert current != item["saved_start_ticks"]
        item["current_start_ticks"] = current
    manifest = {"schema_version": 1, "scope": "Complete byte-exact finite P6 colour parity diagnostic; no qualification receipt.",
                "origin_root": str(WORK), "product_revision": "4bcd9875cd9308deabed06a5fa2f17d2321b0cbc",
                "runner_revision": REVISION,
                "files": {"files/" + relative: {"bytes": value["bytes"], "sha256": value["sha256"], "mode": value["mode"]}
                          for relative, value in files.items()},
                "directories": {"files/" + relative: mode for relative, mode in directories.items()}}
    write(ARCHIVE / "manifest.json", manifest)
    shutil.copy2(ARCHIVE / "manifest.json", ARCHIVE / "whole-byte-manifest.json")
    write(ARCHIVE / "relocation-map.json", relocations)
    process = {"boot_id": boot, "known_owner_children_and_controllers": list(births.values()),
               "all_known_births_absent": True, "scope": "Current inventory and two shader owners; no universal process absence claim."}
    write(ARCHIVE / "terminal-process-checks.json", process)
    disposition = {"decision": "TWO_ORIGINAL_P6_COLOUR_KERNEL_PARITY_CASES_PASSED",
                   "qualification_claimed": False, "full_P6_claimed": False,
                   "integrated_main_science_platform_Mandatory_release_claimed": False,
                   "scope": "Only the two unchanged original ordinary-fp32 kernel/host parity cases on exact WSL Dozen Radeon 780M.",
                   "provider": analysis["selected_device"], "owners": summaries, "cases": cases,
                   "original_postcheck_failure_preserved": True,
                   "separate_overlay_analysis": identity(analysis_path),
                   "consumed_input_bindings_unchanged": len(expected),
                   "cache_policy": "Ordinary caches and complete outputs retained; no cache-hit, cold-cost or speed claim.",
                   "timing_scope": "Owner wall time includes the complete process; actual submission milliseconds are the original test's distinct dispatch measurement.",
                   "2048_MiB_scene_qualification_claimed": False,
                   "raw_files": len(files), "raw_bytes": sum(v["bytes"] for v in files.values())}
    write(ARCHIVE / "disposition.json", disposition)
    result = {"archive": str(ARCHIVE), "manifest": identity(ARCHIVE / "manifest.json"),
              "whole_byte_manifest": identity(ARCHIVE / "whole-byte-manifest.json"),
              "relocation_map": identity(ARCHIVE / "relocation-map.json"),
              "disposition": identity(ARCHIVE / "disposition.json"),
              "raw_files": len(files), "raw_bytes": disposition["raw_bytes"],
              "raw_independent_zlib6_bytes": sum(len(zlib.compress((ROOT / relative).read_bytes(), 6)) for relative in files),
              "known_births_checked": len(births), "protected_bindings_checked": len(expected),
              "preparation_controller": {"pid": os.getpid(), "start_ticks": ticks(os.getpid()), "boot_id": boot},
              "git_changes_made": False, "deleted": False, "qualification_claimed": False}
    write(PREP / "preservation-result.json", result)
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
