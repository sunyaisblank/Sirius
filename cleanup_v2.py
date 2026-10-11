"""Preserve the finished compiler-isolation batch and delete only its recoverable scratch."""
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys

sys.dont_write_bytecode = True
ROOT = Path("/home/astra/.project/Sirius")
WORK = ROOT / "out/literal-transport-compiler-isolation"
HEAD = "b81e061fa081cdaed37eca3add950bd254bf296d"
TAG = "evidence/issue-44-transport-compiler-isolation-b81e061"
ARCHIVE = ROOT / "attestations/review/b81e061/literal-transport-compiler-isolation"
CLEANUP = ROOT / "attestations/cleanup/b81e061/literal-transport-isolation-finished"


def run(command, data=None):
    return subprocess.check_output(command, cwd=ROOT, input=data, text=True).strip()


def seal(path):
    data = Path(path).read_bytes()
    return {"bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}


def write(path, value):
    path.write_text(json.dumps(value, indent=2) + "\n")


def tracked():
    assert WORK == ROOT / "out/literal-transport-compiler-isolation"
    assert run(["git", "rev-parse", "HEAD"]) == HEAD
    assert not run(["git", "status", "--porcelain"])
    assert run(["git", "rev-list", "--left-right", "--count", "HEAD...@{upstream}"]).split() == ["0", "0"]
    assert run(["git", "worktree", "list", "--porcelain"]).count("worktree ") == 1


def phase_directories():
    # Every executed phase has its own immutable before/commands/result receipt.
    names = json.loads((WORK / "terminal-disposition.json").read_text())["phase_directories"]
    assert names in ([".", "frontend", "passprobe"], [".", "frontend", "passprobe", "minprobe"])
    discovered = {str(path.parent.relative_to(WORK)) for path in WORK.rglob("*.json")
                  if path.name in ("before.json", "commands.json", "result.json")}
    assert discovered == set(names), (discovered, names)
    for name in names:
        directory = WORK / name
        expected = ("compile-before.json", "commands.json", "compile-result.json") if name == "." else ("before.json", "commands.json", "result.json")
        assert all((directory / item).is_file() for item in expected), name
    return [WORK / name for name in names]


def process_scan():
    owners, groups = [], set()
    for directory in phase_directories():
        before_name = "compile-before.json" if directory == WORK else "before.json"
        owners.append(json.loads((directory / before_name).read_text())["owner"])
        commands = json.loads((directory / "commands.json").read_text())
        assert len(commands) == (2 if directory == WORK else 1)
        for record in commands:
            assert record["disposition"] in ("completed", "timeout")
            if record["disposition"] == "completed":
                assert isinstance(record["returncode"], int)
                if directory in (WORK, WORK / "frontend"):
                    assert record["returncode"] == 0
            else:
                assert record["returncode"] in (-15, -9)
            cleanup = record["cleanup"]
            assert not cleanup["errors"] and cleanup["known_leader_reaped"]
            assert cleanup["terminal_group_scan"] is not None
            assert not cleanup["terminal_group_scan"]["members"]
            assert not cleanup["forced_cleanup"] or record["disposition"] == "timeout"
            owners.append({"pid": record["pid"], "birth": record["birth"]})
            groups.add(record["pid"])
    matches, unreadable, inspected = [], 0, 0
    for directory in Path("/proc").iterdir():
        if not directory.name.isdecimal():
            continue
        try:
            fields = (directory / "stat").read_text().rsplit(")", 1)[1].split()
            pid, birth, group = int(directory.name), fields[19], int(fields[2])
            inspected += 1
            if group in groups or any(pid == row["pid"] and birth == row["birth"] for row in owners):
                matches.append({"pid": pid, "birth": birth, "group": group})
        except (OSError, ValueError, IndexError):
            unreadable += 1
    assert not matches, matches
    return {"scope": "point-in-time accessible Linux /proc: recorded leader PID/birth identities and recorded new-session process groups; excludes escaped groups, hidden or scheduled work",
            "leaders": owners, "groups": sorted(groups), "accessible_processes": inspected,
            "unreadable_or_raced": unreadable, "matches": matches}


def check_inputs():
    before = json.loads((WORK / "compile-before.json").read_text())
    result = json.loads((WORK / "compile-result.json").read_text())
    assert result["pass"] is False and result["error_type"] == "TimeoutExpired"
    assert result["source_unchanged"] and result["stage_inputs_unchanged"]
    refs = [before[key] for key in ("source", "emitter", "runner", "program_metadata", "emission_metadata", "schedule", "header", "builder", "baseline", "baseline_copy")]
    refs += before["stage_inputs"] + list(before["tools"].values())
    emission = json.loads((WORK / "emission.json").read_text())
    refs += [{**record, "path": str(ROOT / record["path"])} for record in emission["inputs"] + [emission["header"]]]
    for directory in phase_directories()[1:]:
        phase_before = json.loads((directory / "before.json").read_text())
        phase_result = json.loads((directory / "result.json").read_text())
        assert phase_result["input_seals_unchanged"]
        refs += phase_before["input_seals"]
    for directory in phase_directories():
        for record in json.loads((directory / "commands.json").read_text()):
            refs += [record["tool"], record["stdout_seal"], record["stderr_seal"]]
    frontend = json.loads((WORK / "frontend/result.json").read_text())
    assert frontend["frontend_completed"] and not frontend["output_exists"]
    for row in refs:
        assert seal(row["path"]) == {key: row[key] for key in ("bytes", "sha256")}, row["path"]
    assert not list(WORK.glob("literal-transport*.spv"))
    assert not (WORK / "literal-transport.spvasm").exists()
    assert not (WORK / "compiled-inspection.json").exists()
    return {"checked_reference_count": len(refs), "all_recorded_input_seals_unchanged": True,
            "full_builder_candidate_module_exists": False, "inspector_main_executed": False,
            "scope": "recorded source/imports/helpers/builder/tools/header/baseline and all diagnostic phase inputs; no universal or runtime-dependency claim"}


def files():
    paths = []
    for path in sorted(WORK.rglob("*")):
        assert not path.is_symlink(), path
        if path.is_file():
            paths.append(path)
        else:
            assert path.is_dir(), path
    return paths


def preserve():
    tracked()
    assert not ARCHIVE.exists()
    post = {"scope": "unadopted compiler-only complete-source isolation; full production compile timed out, front end completed; no numerical/device/benefit/adoption or full qualification result",
            "source_revision": HEAD, "input_checks": check_inputs(), "process_scan": process_scan()}
    write(WORK / "terminal-contract.json", post)
    chosen = ["emit.py", "compile.py", "inspect_module.py", "frontscreen.py", "passscreen.py", "finish.py",
              "emission.json", "programs.json", "schedule.json", "layer-boundary-check.json",
              "source-independent-review.json", "apparatus-independent-review.json",
              "passprobe-independent-review.json", "terminal-contract.json", "compile-result.json",
              "terminal-disposition.json", "cause-and-next-lead-review.json"]
    if WORK / "minprobe" in phase_directories():
        chosen += ["minscreen.py", "minprobe-independent-review.json"]
    assert subprocess.run(["git", "show-ref", "--verify", "--quiet", "refs/tags/" + TAG], cwd=ROOT).returncode == 1
    entries, refs = [], []
    for name in sorted(chosen):
        path = WORK / name
        value = seal(path)
        blob = run(["git", "hash-object", "-w", "--stdin"], path.read_text())
        assert subprocess.check_output(["git", "cat-file", "blob", blob], cwd=ROOT) == path.read_bytes()
        entries.append(f"100644 blob {blob}\t{name}\n")
        refs.append({"name": name, "blob": blob, **value})
    tree = run(["git", "mktree"], "".join(entries))
    commit = run(["git", "commit-tree", tree], "Preserve literal Transport compiler isolation tools and findings (#44)\n\nSource-only diagnostic snapshot; not a product patch or qualification pass.\n")
    run(["git", "tag", TAG, commit])
    run(["git", "push", "origin", "refs/tags/" + TAG])
    assert run(["git", "ls-remote", "--tags", "origin", "refs/tags/" + TAG]).split() == [commit, "refs/tags/" + TAG]
    write(WORK / "source-preservation.json", {"tag": TAG, "commit": commit, "source_only": True,
                                              "remote_exact": True, "files": refs})
    ARCHIVE.mkdir(parents=True)
    entries = []
    for path in files():
        relative = path.relative_to(WORK)
        destination = ARCHIVE / "diagnostic" / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        expected = seal(path)
        shutil.copyfile(path, destination)
        assert seal(destination) == expected
        entries.append({"path": str(destination.relative_to(ARCHIVE)), "scratch_path": str(relative), **expected})
    manifest = {"source_revision": HEAD, "scope": post["scope"], "source_tag": TAG,
                "source_snapshot": commit, "files": entries,
                "payload_count": len(entries), "payload_bytes": sum(entry["bytes"] for entry in entries)}
    write(ARCHIVE / "manifest.json", manifest)
    print(json.dumps({"archive": str(ARCHIVE.relative_to(ROOT)), "manifest": seal(ARCHIVE / "manifest.json"),
                      "source_snapshot": commit, "source_tag": TAG,
                      "payload_count": manifest["payload_count"], "payload_bytes": manifest["payload_bytes"]}))


def cleanup():
    tracked()
    assert not CLEANUP.exists()
    approval = ROOT / "attestations/review/b81e061/literal-transport-isolation-publication/cleanup-independent-review.json"
    review = json.loads(approval.read_text())
    assert review["accepted"] is True
    assert review["manifest"] == seal(ARCHIVE / "manifest.json")
    assert review["reviewed_finisher"] == seal(Path(__file__))
    manifest = json.loads((ARCHIVE / "manifest.json").read_text())
    source = json.loads((WORK / "source-preservation.json").read_text())
    assert source["tag"] == TAG and source["commit"] == manifest["source_snapshot"]
    assert run(["git", "ls-remote", "--tags", "origin", "refs/tags/" + TAG]).split() == [source["commit"], "refs/tags/" + TAG]
    assert sorted(str(path.relative_to(WORK)) for path in files()) == sorted(row["scratch_path"] for row in manifest["files"])
    for row in manifest["files"]:
        expected = {key: row[key] for key in ("bytes", "sha256")}
        assert seal(WORK / row["scratch_path"]) == expected
        assert seal(ARCHIVE / row["path"]) == expected
    for row in source["files"]:
        data = subprocess.check_output(["git", "show", source["commit"] + ":" + row["name"]], cwd=ROOT)
        assert {"bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()} == {key: row[key] for key in ("bytes", "sha256")}
    inputs, processes = check_inputs(), process_scan()
    protected = [ROOT / "bin/linux-gcc/src/sirius/backend/retained/retained_kernels.h",
                 *sorted((ROOT / "bin/linux-gcc/src/sirius/backend/retained").glob("*.spv"))]
    before = {str(path.relative_to(ROOT)): seal(path) for path in protected}
    shutil.rmtree(WORK)
    assert not WORK.exists()
    assert before == {str(path.relative_to(ROOT)): seal(path) for path in protected}
    tracked()
    CLEANUP.mkdir(parents=True)
    receipt = {"source_revision": HEAD, "removed_only": str(WORK.relative_to(ROOT)),
               "removed_regular_files": len(manifest["files"]), "removed_bytes": manifest["payload_bytes"],
               "exact_recoveries": manifest["files"], "manifest": seal(ARCHIVE / "manifest.json"),
               "source_tag": TAG, "source_snapshot": source["commit"], "inputs": inputs,
               "process_scan": processes, "protected_build_payloads_unchanged": before,
               "tracked_tree_clean_upstream_zero_zero_single_worktree": True}
    write(CLEANUP / "cleanup.json", receipt)
    print(json.dumps({"cleanup": str((CLEANUP / "cleanup.json").relative_to(ROOT)),
                      "removed_regular_files": receipt["removed_regular_files"],
                      "removed_bytes": receipt["removed_bytes"], "receipt": seal(CLEANUP / "cleanup.json")}))


if __name__ == "__main__":
    assert len(sys.argv) == 2 and sys.argv[1] in ("--preserve", "--cleanup")
    (preserve if sys.argv[1] == "--preserve" else cleanup)()
