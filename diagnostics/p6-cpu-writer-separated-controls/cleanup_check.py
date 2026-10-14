#!/usr/bin/env python3
"""Read-only pre-removal check for the two fixed CPU writer scratch roots.

Requires complete raw/Git preservation and published tags. Never deletes files;
process absence/consumer checks are sampled observations, not a future lock.
"""
from pathlib import Path
import argparse
import hashlib
import json
import os
import stat
import subprocess
import sys

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[4]
BASE = ROOT / "attestations/scientific-core/4bcd987"
PARENT = "a1d6e75299a2bf250e5ce7d44ad5965f8acb4299"
OLD_MANIFEST = "8b5b860301a6b697ce2505ab1ac141353840abe2589a75e9403ca3db041cab97"


def require(ok, message):
    if not ok:
        raise ValueError(message)


def identity(path):
    raw = path.read_bytes()
    return {"bytes": len(raw), "sha256": hashlib.sha256(raw).hexdigest()}


def git(*args):
    return subprocess.check_output(["git", *args], cwd=ROOT, timeout=30)


def snapshot(folder):
    result = {}
    require(folder.is_dir() and not folder.is_symlink(), "scratch/archive root is not an ordinary directory")
    for path in [folder, *sorted(folder.rglob("*"))]:
        status = path.lstat()
        require(stat.S_ISDIR(status.st_mode) or stat.S_ISREG(status.st_mode), "symlink or special file: " + str(path))
        result[str(path.relative_to(folder))] = {
            "mode": stat.S_IMODE(status.st_mode), "uid": status.st_uid, "gid": status.st_gid,
            "kind": "directory" if path.is_dir() else "file",
            **({} if path.is_dir() else identity(path)),
        }
    return result


def birth(pid):
    try:
        return (Path("/proc") / str(pid) / "stat").read_text().rsplit(")", 1)[1].split()[19]
    except FileNotFoundError:
        return None


def verify_archive(name, scratch_name, tag, expected_manifest):
    archive = BASE / name
    scratch = ROOT / "out" / scratch_name
    raw = archive / "files/out" / scratch_name
    require(identity(archive / "manifest.json")["sha256"] == expected_manifest, "manifest pin differs")
    manifest = json.loads((archive / "manifest.json").read_text())
    recorded = {row["path"]: row for row in manifest["files"]}
    actual = {str(path.relative_to(archive)) for path in archive.rglob("*") if path.is_file()}
    require(actual == set(recorded) | {"manifest.json"}, "archive file set differs")
    for name, row in recorded.items():
        require(identity(archive / name) == {key: row[key] for key in ("bytes", "sha256")}, "archive payload differs: " + name)
    require(snapshot(scratch) == snapshot(raw), "scratch bytes/modes/ownership differ from raw archive")
    reference = "refs/tags/" + tag
    local_object = git("rev-parse", reference).decode().strip()
    commit = git("rev-parse", reference + "^{commit}").decode().strip()
    require(git("show", "-s", "--format=%P", commit).decode().strip() == PARENT, "artifact parent differs")
    remote = git("ls-remote", "--exit-code", "origin", reference).decode().splitlines()
    require(remote == [local_object + "\t" + reference], "published tag differs or is missing")
    prefix = "diagnostics/" + archive.name + "/"
    names = git("ls-tree", "-r", "--name-only", commit).decode().splitlines()
    require(set(names) == {prefix + name for name in actual}, "artifact is not the complete isolated evidence tree")
    for name in actual:
        blob = git("show", commit + ":" + prefix + name)
        require({"bytes": len(blob), "sha256": hashlib.sha256(blob).hexdigest()} == identity(archive / name), "Git payload differs: " + name)
    for path in scratch.rglob("owner.json"):
        owner = json.loads(path.read_text())
        require(owner["status"] == "terminal" and owner["cleanup_errors"] == []
                and owner["remaining_owned_processes"] == [], "owner retirement is incomplete")
        members = [(row["pid"], row["saved_start_ticks"]) for row in owner["observed_births"]]
        members.append((owner["controller"]["pid"], owner["controller"]["start_ticks"]))
        require(all(birth(pid) != saved for pid, saved in members), "an original process birth remains")
    return scratch, {"scratch": str(scratch.relative_to(ROOT)), "tag": tag,
                     "commit": commit, "manifest_sha256": expected_manifest, "archive_files": len(actual)}


def consumers(roots):
    prefixes = tuple(str(path) for path in roots)
    found = []
    inaccessible = []
    def touches(value):
        return any(value == prefix or value.startswith(prefix + "/")
                   or value == prefix + " (deleted)" for prefix in prefixes)
    for process in Path("/proc").glob("[0-9]*"):
        try:
            saved = birth(process.name)
            if saved is None:
                continue
            links = [process / "cwd", process / "root", process / "exe", *list((process / "fd").iterdir())]
            for link in links:
                try:
                    target = os.readlink(link)
                except FileNotFoundError:
                    continue
                if touches(target):
                    found.append({"pid": int(process.name), "link": str(link), "target": target})
            for line in (process / "maps").read_text().splitlines():
                fields = line.split(None, 5)
                if len(fields) == 6 and touches(fields[5]):
                    found.append({"pid": int(process.name), "mapping": fields[5]})
            require(birth(process.name) in (None, saved), "process birth changed during consumer census")
        except PermissionError as error:
            inaccessible.append({"pid": int(process.name), "path": error.filename,
                                 "reason": "permission_denied"})
        except (FileNotFoundError, ProcessLookupError):
            continue
    return found, inaccessible


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--completed-tag", required=True)
    parser.add_argument("--manifest-sha256", required=True)
    args = parser.parse_args()
    require(Path.cwd().resolve() == ROOT, "run from the Sirius checkout")
    require(git("rev-parse", "HEAD").decode().strip() == PARENT and not git("status", "--porcelain").strip(), "source state differs")
    old, old_result = verify_archive("p6-cpu-writer-diagnostic-limit", "p6-cpu-writer-current-source",
                                    "evidence/issue-53-p6-writer-diagnostic-limit", OLD_MANIFEST)
    new, new_result = verify_archive("p6-cpu-writer-separated-controls", "p6-cpu-writer-separated-controls",
                                    args.completed_tag, args.manifest_sha256)
    live, inaccessible = consumers([old, new])
    require(not live, "scratch has a live consumer: " + json.dumps(live))
    print(json.dumps({"sampled_cleanup_preconditions_passed": True, "deleted": False,
                      "roots": [old_result, new_result], "live_consumers": live,
                      "inaccessible_processes": inaccessible,
                      "boundary": "Read-only original-owner birth and accessible-process checks only; permission-limited processes are listed, so global consumer absence is not claimed. Root performs any later removal."}, indent=2))


if __name__ == "__main__":
    main()
