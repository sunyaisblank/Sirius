"""Read-only pre-cleanup audit for these exact two completed scratch roots.

Never deletes, changes permissions, signals a process, or creates a receipt.
Permission/race gaps are listed; no global consumer-absence claim is made.
"""
from pathlib import Path
import argparse
import hashlib
import json
import os
import stat
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
PREP = Path(__file__).resolve().parent
WORK = ROOT / "out/p6-colour-kernel-parity-current-source"
ARCHIVE = ROOT / "attestations/scientific-core/4bcd987/p6-colour-kernel-parity-current-source"
MANIFEST_SHA256 = "19e706c79e47170b0d87a13aae6261a62bdc156217dd68464c245d40b4db284d"
PREFIX = "diagnostics/p6-colour-kernel-parity-current-source"


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def identity(path):
    path = Path(path)
    return {"resolved_path": str(path.resolve(strict=True)), "bytes": path.stat().st_size,
            "sha256": digest(path)}


def query(argv):
    return subprocess.check_output(["git", *argv], cwd=ROOT)


def birth(pid):
    try:
        return Path(f"/proc/{pid}/stat").read_text().split(")", 1)[1].split()[19]
    except FileNotFoundError:
        return None


def within(target):
    target = target.removesuffix(" (deleted)")
    return any(target == str(root) or target.startswith(str(root) + "/") for root in [WORK, PREP])


def consumers():
    hits, gaps, checked = [], [], []
    for directory in sorted(Path("/proc").glob("[0-9]*"), key=lambda p: int(p.name)):
        pid = int(directory.name)
        if pid == os.getpid():
            continue
        try:
            first = birth(pid)
            if first is None:
                continue
            # proc-directory ownership may change for a non-dumpable process;
            # use the process's reported real/effective/saved/filesystem UIDs.
            uid_line = next(line for line in (directory / "status").read_text().splitlines()
                            if line.startswith("Uid:"))
            uids = [int(value) for value in uid_line.split()[1:]]
            uid = os.getuid() if os.getuid() in uids else uids[0]
        except FileNotFoundError:
            continue
        except (PermissionError, ProcessLookupError) as error:
            gaps.append({"pid": pid, "uid": None, "operation": "process identity", "error": repr(error)})
            continue
        def gap(operation, error):
            gaps.append({"pid": pid, "uid": uid, "operation": operation, "error": repr(error)})
        for name in ["cwd", "exe"]:
            try:
                target = os.readlink(directory / name)
                if within(target):
                    hits.append({"pid": pid, "start_ticks": first, "kind": name, "target": target})
            except (FileNotFoundError, ProcessLookupError):
                pass
            except PermissionError as error:
                gap(name, error)
        try:
            lines = (directory / "maps").read_text().splitlines()
            for line in lines:
                fields = line.split(None, 5)
                if len(fields) == 6 and within(fields[5]):
                    hits.append({"pid": pid, "start_ticks": first, "kind": "map", "target": fields[5]})
        except (FileNotFoundError, ProcessLookupError):
            pass
        except PermissionError as error:
            gap("maps", error)
        try:
            descriptors = list((directory / "fd").iterdir())
        except (FileNotFoundError, ProcessLookupError):
            descriptors = []
        except PermissionError as error:
            gap("fd directory", error)
            descriptors = []
        for descriptor in descriptors:
            try:
                target = os.readlink(descriptor)
                if within(target):
                    hits.append({"pid": pid, "start_ticks": first, "kind": "fd", "target": target})
            except (FileNotFoundError, ProcessLookupError):
                pass
            except PermissionError as error:
                gap("fd link " + descriptor.name, error)
        try:
            last = birth(pid)
            if last not in [first, None]:
                gap("birth changed during census", RuntimeError((first, last)))
        except PermissionError as error:
            gap("final birth", error)
        checked.append({"pid": pid, "uid": uid, "reported_uids": uids, "saved_start_ticks": first})
    return {"scope": "Accessible Linux cwd/exe/maps/fd census, all known owned births, and complete current-UID access where available. Foreign/unknown permission gaps remain explicit.",
            "excluded_audit_pid": os.getpid(), "checked_processes": checked,
            "consumer_hits": hits, "access_gaps": gaps,
            "current_uid_or_unknown_gaps": [g for g in gaps if g["uid"] in [os.getuid(), None]],
            "global_consumer_absence_claimed": False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evidence-ref", required=True, help="Exact pushed evidence tag (refs/tags/...).")
    parser.add_argument("--git-artifact-manifest-sha256", required=True)
    args = parser.parse_args()
    report = {"deleted": False, "read_only": True, "qualification_claimed": False,
              "candidate_roots": [str(WORK), str(PREP)], "checks_passed": False,
              "raw_git_and_known_birth_checks_passed": False,
              "errors": []}
    try:
        assert Path.cwd() == ROOT
        assert query(["rev-parse", "HEAD"]).decode().strip() == "a1d6e75299a2bf250e5ce7d44ad5965f8acb4299"
        assert query(["status", "--porcelain"]) == b""
        assert WORK.resolve(strict=True) == WORK and PREP.resolve(strict=True) == PREP
        assert ARCHIVE.resolve(strict=True) == ARCHIVE
        assert digest(ARCHIVE / "manifest.json") == MANIFEST_SHA256
        assert digest(ARCHIVE / "whole-byte-manifest.json") == MANIFEST_SHA256
        manifest = json.loads((ARCHIVE / "manifest.json").read_text())
        expected_live = {}
        for relative, item in manifest["files"].items():
            archived = ARCHIVE / relative
            live = ROOT / relative.removeprefix("files/")
            for path in [archived, live]:
                assert stat.S_ISREG(path.lstat().st_mode) and not path.is_symlink(), path
                assert path.stat().st_size == item["bytes"] and digest(path) == item["sha256"], path
                assert stat.S_IMODE(path.stat().st_mode) == item["mode"], path
            expected_live[str(live.relative_to(WORK))] = item
        actual_live = {}
        actual_dirs = {"files/" + str(WORK.relative_to(ROOT))}
        for path in WORK.rglob("*"):
            assert not path.is_symlink(), path
            if path.is_file():
                actual_live[str(path.relative_to(WORK))] = path
            else:
                assert path.is_dir(), path
                actual_dirs.add("files/" + str(path.relative_to(ROOT)))
        assert set(actual_live) == set(expected_live)
        expected_dirs = manifest["directories"]
        assert actual_dirs == set(expected_dirs)
        for relative, mode in expected_dirs.items():
            live = ROOT / relative.removeprefix("files/")
            archived = ARCHIVE / relative
            for path in [live, archived]:
                assert stat.S_ISDIR(path.lstat().st_mode) and stat.S_IMODE(path.stat().st_mode) == mode
        payload = ARCHIVE / "files" / WORK.relative_to(ROOT)
        assert {str(path.relative_to(payload)) for path in payload.rglob("*") if path.is_file()} == set(expected_live)
        archived_dirs = {str(payload.relative_to(ARCHIVE)), *[str(path.relative_to(ARCHIVE)) for path in payload.rglob("*") if path.is_dir()]}
        assert archived_dirs == set(expected_dirs)
        relocations = json.loads((ARCHIVE / "relocation-map.json").read_text())
        assert set(relocations) == {relative.removeprefix("files/") for relative in manifest["files"]}
        for relative, item in relocations.items():
            expected = manifest["files"]["files/" + relative]
            assert item["archive_relative"] == "files/" + relative
            assert item["origin"] == {**identity(ROOT / relative), "mode": expected["mode"]}
            assert item["archived"] == identity(ARCHIVE / "files" / relative)
        for folder in [WORK / "inventory", WORK / "shader-phase/cie", WORK / "shader-phase/xyz"]:
            owner = json.loads((folder / "owner.json").read_text())
            assert owner["status"] == "terminal" and owner["passed"] and owner["child_exit"] == 0
            assert owner["stop_reason"] is None and owner["cleanup_errors"] == [] and owner["remaining_owned_processes"] == []
            assert owner["source_unchanged"] and owner["whole_inputs_unchanged"]
            before = json.loads((folder / "inputs-before.json").read_text())
            assert before == json.loads((folder / "inputs-after.json").read_text())
            for path, value in before.items():
                assert identity(ROOT / path) == value, ("changed original input", path)
        index = PREP / "git-artifacts-manifest.json"
        assert digest(index) == args.git_artifact_manifest_sha256
        assert digest(ARCHIVE / "git-artifacts-manifest.json") == args.git_artifact_manifest_sha256
        selection = json.loads(index.read_text())
        artifacts = PREP / "git-artifacts"
        assert {str(path.relative_to(artifacts)) for path in artifacts.rglob("*") if path.is_file()} == set(selection["files"])
        assert {str(path.relative_to(artifacts)) for path in artifacts.rglob("*") if path.is_dir()} == set(selection["directories"])
        for path in artifacts.rglob("*"):
            assert not path.is_symlink(), path
            if path.is_dir():
                assert stat.S_IMODE(path.stat().st_mode) == selection["directories"][str(path.relative_to(artifacts))]
        reference = args.evidence_ref
        assert reference.startswith("refs/tags/evidence/")
        commit = query(["rev-parse", "--verify", reference + "^{commit}"]).decode().strip()
        remote = dict((line.split()[1], line.split()[0]) for line in query(["ls-remote", "origin", reference, reference + "^{}"]).decode().splitlines())
        assert remote.get(reference + "^{}", remote.get(reference)) == commit, "Evidence tag is not exact on origin"
        git_entries = {}
        for record in query(["ls-tree", "-rz", commit, "--", PREFIX]).split(b"\0"):
            if record:
                header, path = record.split(b"\t", 1)
                mode, kind, blob = header.decode().split()
                assert kind == "blob"
                git_entries[os.fsdecode(path)] = (mode, blob)
        assert set(git_entries) == set(selection["files"]), "Git artifact tree differs from the prepared full selection"
        for path, item in selection["files"].items():
            live = PREP / "git-artifacts" / path
            recovered = ARCHIVE / item["archive_relative"]
            for source in [live, recovered]:
                assert not source.is_symlink() and source.stat().st_size == item["bytes"] and digest(source) == item["sha256"]
            mode, blob = git_entries[path]
            assert mode == item["git_mode"]
            data = query(["cat-file", "blob", blob])
            assert len(data) == item["bytes"] and hashlib.sha256(data).hexdigest() == item["sha256"]
        # Every top-level preservation helper/result has raw recovery too.
        top = {p.name for p in PREP.iterdir()}
        assert top == {"preserve.py", "audit_cleanup.py", "preservation-result.json", "git-artifacts-manifest.json", "git-artifacts"}
        for name, relative in [("preserve.py", "preservation-tools/preserve.py"),
                               ("audit_cleanup.py", "preservation-tools/audit_cleanup.py"),
                               ("preservation-result.json", "preservation-result.json")]:
            assert digest(PREP / name) == digest(ARCHIVE / relative)
        boot = Path("/proc/sys/kernel/random/boot_id").read_text().strip()
        processes = json.loads((ARCHIVE / "terminal-process-checks.json").read_text())
        assert processes["boot_id"] == boot
        known = processes["known_owner_children_and_controllers"]
        preparer = json.loads((PREP / "preservation-result.json").read_text())["preparation_controller"]
        assert preparer["boot_id"] == boot
        known = [*known, {"pid": preparer["pid"], "saved_start_ticks": preparer["start_ticks"]}]
        checks = []
        for item in known:
            current = birth(item["pid"])
            assert current != str(item["saved_start_ticks"]), item
            checks.append({"pid": item["pid"], "saved_start_ticks": str(item["saved_start_ticks"]),
                           "current_start_ticks": current, "absent": True})
        report.update(raw_files_verified=len(expected_live), git_files_verified=len(git_entries),
                      raw_manifest_sha256=MANIFEST_SHA256,
                      git_artifact_manifest_sha256=args.git_artifact_manifest_sha256,
                      evidence_ref=reference, evidence_commit=commit, exact_pushed_tag=True,
                      known_birth_checks=checks, raw_git_and_known_birth_checks_passed=True)
    except Exception as error:
        report["errors"].append(repr(error))
    try:
        census = consumers()
        report["accessible_consumer_census"] = census
        report["consumer_census_complete"] = not census["access_gaps"]
        report["accessible_consumers_clear"] = not census["consumer_hits"]
        if census["consumer_hits"]:
            report["errors"].append("Accessible consumers still reference a candidate root")
        # Unrelated inaccessible processes are observations outside this scope,
        # not proof of absence. Every known owned birth was checked strictly
        # above; any access failure for a known birth fails those checks.
    except Exception as error:
        report["errors"].append("Consumer census failed: " + repr(error))
    report["checks_passed"] = not report["errors"]
    report["meaning"] = "Read-only recovery/known-birth/accessible-consumer audit; no deletion, permission change, signal, global absence or qualification verdict."
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0 if report["checks_passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
