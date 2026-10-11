"""Preserve this failed compiler screen, then remove only its verified scratch."""
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys

sys.dont_write_bytecode = True
WORK = Path(__file__).resolve().parent
ROOT = WORK.parents[1]
HEAD = "b81e061fa081cdaed37eca3add950bd254bf296d"
TAG = "evidence/issue-44-literal-transport-compiler-b81e061"
ARCHIVE = ROOT / "attestations/review/b81e061/literal-transport-compiler-screen"
CLEANUP = ROOT / "attestations/cleanup/b81e061/literal-transport-finished"


def run(command, data=None):
    return subprocess.check_output(command, cwd=ROOT, input=data, text=True).strip()


def seal(path):
    data = path.read_bytes()
    return {"bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}


def write(path, value):
    path.write_text(json.dumps(value, indent=2) + "\n")


def tracked():
    assert ROOT == Path("/home/astra/.project/Sirius")
    assert WORK == ROOT / "out/literal-transport-feasibility"
    assert run(["git", "rev-parse", "HEAD"]) == HEAD
    assert not run(["git", "status", "--porcelain"])
    assert run(["git", "rev-list", "--left-right", "--count", "HEAD...@{upstream}"]).split() == ["0", "0"]
    assert run(["git", "worktree", "list", "--porcelain"]).count("worktree ") == 1


def process_scan():
    before = json.loads((WORK / "compile-before.json").read_text())
    commands = json.loads((WORK / "commands.json").read_text())
    assert len(commands) == 2
    assert commands[0]["disposition"] == "completed" and commands[0]["returncode"] == 0
    assert commands[1]["disposition"] == "timeout" and commands[1]["returncode"] == -15
    owners = [before["owner"], *[{"pid": record["pid"], "birth": record["birth"]} for record in commands]]
    groups = {record["pid"] for record in commands}
    matches, unreadable, inspected = [], 0, 0
    for directory in Path("/proc").iterdir():
        if not directory.name.isdecimal():
            continue
        try:
            fields = (directory / "stat").read_text().rsplit(")", 1)[1].split()
            birth, group = fields[19], int(fields[2])
            inspected += 1
            pid = int(directory.name)
            if group in groups or any(pid == owner["pid"] and birth == owner["birth"] for owner in owners):
                matches.append({"pid": pid, "birth": birth, "group": group})
        except (OSError, ValueError, IndexError):
            unreadable += 1
    assert not matches, matches
    return {"scope": "point-in-time accessible Linux /proc: recorded leader PID/birth identities and two recorded new-session process groups; not escaped groups, hidden or scheduled work",
            "leaders": owners, "groups": sorted(groups), "accessible_processes": inspected,
            "unreadable_or_raced": unreadable, "matches": matches}


def check_inputs():
    before = json.loads((WORK / "compile-before.json").read_text())
    emission = json.loads((WORK / "emission.json").read_text())
    compile_result = json.loads((WORK / "compile-result.json").read_text())
    assert compile_result["pass"] is False and compile_result["error_type"] == "TimeoutExpired"
    assert compile_result["source_unchanged"] and compile_result["stage_inputs_unchanged"]
    refs = [before["source"], before["emitter"], before["builder"], before["baseline"], before["baseline_copy"],
            *before["stage_inputs"], *before["tools"].values()]
    refs += [{**record, "path": str(ROOT / record["path"])} for record in emission["inputs"] + [emission["header"]]]
    runtime = json.loads((WORK / "compiler-runtime-readback.json").read_text())
    refs += runtime["files"]
    for record in refs:
        assert seal(Path(record["path"])) == {key: record[key] for key in ("bytes", "sha256")}, record["path"]
    assert not list(WORK.glob("literal-transport*.spv"))
    assert not (WORK / "literal-transport.spvasm").exists()
    assert not (WORK / "compiled-inspection.json").exists()
    return {"checked_reference_count": len(refs), "all_recorded_input_seals_unchanged": True,
            "candidate_module_exists": False, "inspector_main_executed": False,
            "postcheck_scope": "recorded source/imports/emitter/builder/tools/header/baseline and current compiler mapped-file readback"}


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
    post = {"scope": "failed compiler-only complete emission; no module, device workload, numerical comparison, benefit result, adoption or full qualification",
            "source_revision": HEAD, "input_checks": check_inputs(), "process_scan": process_scan()}
    write(WORK / "terminal-contract.json", post)
    chosen = ["emit.py", "compile.py", "inspect_module.py", "finish.py", "emission.json",
              "source-independent-review.json", "apparatus-independent-review.json",
              "terminal-run-readback.json", "terminal-contract.json", "compile-result.json",
              "compiler-independent-disposition.json"]
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
    commit = run(["git", "commit-tree", tree], "Preserve failed literal Transport compiler screen tools (#44)\n\nSource-only diagnostic snapshot; not a product patch or qualification pass.\n")
    run(["git", "tag", TAG, commit])
    run(["git", "push", "origin", "refs/tags/" + TAG])
    remote = run(["git", "ls-remote", "--tags", "origin", "refs/tags/" + TAG]).split()
    assert remote == [commit, "refs/tags/" + TAG]
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
        entries.append({"path": str(destination.relative_to(ARCHIVE)),
                        "scratch_path": str(relative), **expected})
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
    manifest = json.loads((ARCHIVE / "manifest.json").read_text())
    source = json.loads((WORK / "source-preservation.json").read_text())
    assert source["tag"] == TAG and source["commit"] == manifest["source_snapshot"]
    assert run(["git", "ls-remote", "--tags", "origin", "refs/tags/" + TAG]).split() == [source["commit"], "refs/tags/" + TAG]
    expected_names = sorted(entry["scratch_path"] for entry in manifest["files"])
    assert [str(path.relative_to(WORK)) for path in files()] == expected_names
    for record in manifest["files"]:
        expected = {key: record[key] for key in ("bytes", "sha256")}
        assert seal(WORK / record["scratch_path"]) == expected
        assert seal(ARCHIVE / record["path"]) == expected
    inputs = check_inputs()
    processes = process_scan()
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
