"""Fixed five-root retirement; root operates only after its focused runtime ends.

Default writes an audit. --remove requires the exact prior audit hash and repeats
the checks. No test, build, provider, Windows process query, Git write or tag.
"""
from pathlib import Path, PurePosixPath
import argparse
import hashlib
import json
import os
import shutil
import stat
import tempfile
import zipfile

import mature_cleanup as mature

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
CONFIG = json.loads((HERE / "configuration.json").read_text())
TARGETS = tuple(ROOT / name for name in CONFIG["targets"])
RECEIPT = ROOT / "attestations/cleanup/2026-10-15/native-terminal-scratch"
require = mature.require
identity = mature.identity


def load(path):
    return json.loads(Path(path).read_text(encoding="utf-8-sig"))


def inside(path):
    return any(path == target or path.is_relative_to(target) for target in TARGETS)


def inventory(root, maximum=10000):
    require(root in TARGETS and root.is_dir() and not root.is_symlink()
            and root.resolve() == root, "fixed owned root is missing or aliased")
    files, directories, visited = {}, {".": stat.S_IMODE(root.stat().st_mode)}, 0
    for base, dirs, names in os.walk(root, followlinks=False):
        for name in [*dirs, *names]:
            path = Path(base) / name
            visited += 1
            require(visited <= maximum, "fixed root exceeds its explicit entry bound")
            info = path.lstat()
            require(not stat.S_ISLNK(info.st_mode), "symlink in owned root: " + str(path))
            key = path.relative_to(root).as_posix()
            if stat.S_ISREG(info.st_mode):
                files[key] = path
            else:
                require(stat.S_ISDIR(info.st_mode), "special entry: " + str(path))
                directories[key] = stat.S_IMODE(info.st_mode)
    return files, directories


def compare(root, expected, maximum=10000):
    files, directories = inventory(root, maximum)
    require(set(files) == set(expected), "unaccounted/missing files under " + str(root))
    records = {}
    for name, path in files.items():
        actual = identity(path)
        require(actual == {k: expected[name][k] for k in ("bytes", "sha256")},
                "preserved/member byte mismatch: " + str(path))
        records[name] = dict(actual, mode=stat.S_IMODE(path.stat().st_mode))
    return {"root": str(root.relative_to(ROOT)), "files": records,
            "directories": directories, "regular_files": len(files),
            "regular_bytes": sum(row["bytes"] for row in records.values())}


def direct(bundle, relative):
    prefix = "files/" + relative + "/"
    return {name[len(prefix):]: pin for name, pin in load(bundle / "manifest.json")["files"].items()
            if name.startswith(prefix)}


def verify_zip(bundle, prefix, expected_zip, members, count):
    """Read immutable Git-backed parts, never extract or execute an archive."""
    manifest = load(bundle / "manifest.json")["files"]
    parts = sorted(name for name in manifest if name.startswith(prefix + "/"))
    require(parts == [prefix + "/part-" + str(i).zfill(3) for i in range(len(parts))]
            and len(parts) in (2, 3), "exact contiguous ZIP parts required")
    digest, size = hashlib.sha256(), 0
    with tempfile.TemporaryFile(dir=HERE) as joined:
        for name in parts:
            path = bundle / name
            require(identity(path) == manifest[name], "archive part changed")
            with path.open("rb") as stream:
                for chunk in iter(lambda: stream.read(1024 * 1024), b""):
                    joined.write(chunk)
                    digest.update(chunk)
                    size += len(chunk)
                    require(size <= expected_zip["bytes"], "ZIP exceeds exact byte limit")
        require({"bytes": size, "sha256": digest.hexdigest()} == expected_zip,
                "whole ZIP differs from original preserved transport")
        joined.seek(0)
        with zipfile.ZipFile(joined) as archive:
            infos = archive.infolist()
            require(len(infos) == count == len(members), "complete ZIP member count differs")
            require(len({item.filename for item in infos}) == count, "duplicate ZIP member")
            require({item.filename for item in infos} == set(members), "complete member inventory differs")
            total = 0
            for item in infos:
                name = item.filename
                path = PurePosixPath(name)
                require(not item.is_dir() and not item.flag_bits & 1 and name == str(path)
                        and not path.is_absolute() and ".." not in path.parts
                        and "\\" not in name and "\x00" not in name,
                        "unsafe/directory/encrypted member: " + name)
                require(stat.S_IFMT(item.external_attr >> 16) in (0, stat.S_IFREG),
                        "nonregular ZIP member: " + name)
                wanted = members[name]
                require(item.file_size == wanted["bytes"], "member size changed")
                h, length = hashlib.sha256(), 0
                with archive.open(item) as stream:
                    for chunk in iter(lambda: stream.read(1024 * 1024), b""):
                        h.update(chunk)
                        length += len(chunk)
                        require(length <= wanted["bytes"], "member exceeded exact bound")
                require({"bytes": length, "sha256": h.hexdigest()} == wanted,
                        "complete ZIP member/CRC mismatch: " + name)
                total += length
            require(total == sum(row["bytes"] for row in members.values()), "expanded total differs")
    return {"parts": {name: manifest[name] for name in parts}, "zip": expected_zip,
            "members_verified": count, "expanded_bytes": total, "extracted": False}


def preservation():
    # Ordinary raw/Git/remote checks retain the mature preserved() unchanged.
    evidence = {name: mature.preserved(row) for name, row in CONFIG["evidence"].items()}
    bundles = {name: ROOT / row[2] for name, row in CONFIG["evidence"].items()}
    roots = {}
    for name in ("recovery", "linux"):
        for relative in CONFIG["evidence"][name][4]:
            roots[relative] = compare(ROOT / relative, direct(bundles[name], relative))
    win = TARGETS[0]
    recovery = load(bundles["windows"] / "files" / win.relative_to(ROOT) / "export-recovery.json")
    require(recovery["whole_zip"] == CONFIG["windows_zip"]
            and recovery["member_count"] == CONFIG["windows_members"], "Windows original export changed")
    win_members = recovery["member_identities"]
    win_zip = verify_zip(bundles["export"], "windows.zip.parts", CONFIG["windows_zip"],
                         win_members, CONFIG["windows_members"])
    win_files = direct(bundles["windows"], str(win.relative_to(ROOT)))
    require("export.zip" not in win_files, "Windows archive exclusion changed")
    win_files["export.zip"] = CONFIG["windows_zip"]
    roots[str(win.relative_to(ROOT))] = compare(win, win_files)
    stage = TARGETS[1]
    stage_files = {"payload/" + name: value for name, value in win_members.items()}
    mapped = {}
    for row in recovery["mapping"]:
        path = ROOT / row["staged_path"]
        require(path.is_relative_to(stage), "Windows staged mapping escaped fixed bin root")
        value = {k: row[k] for k in ("bytes", "sha256")}
        require(value == win_members[row["export_member"]], "staged mapping does not join original member")
        key = path.relative_to(stage).as_posix()
        require(key not in mapped or mapped[key] == value, "conflicting duplicated mapping")
        mapped[key] = value
    require(len(mapped) == CONFIG["windows_unique_staged_paths"] == 46
            and len(recovery["mapping"]) == 47, "Windows staged mapping count changed")
    require(not set(mapped) & set(stage_files), "staged/payload collision")
    stage_files.update(mapped)
    require(len(stage_files) == CONFIG["windows_bin_files"] == 86, "Windows bin inventory changed")
    roots[str(stage.relative_to(ROOT))] = compare(stage, stage_files)
    mac = TARGETS[4]
    member_path = bundles["macos"] / "files" / mac.relative_to(ROOT) / "artifact-members.json"
    members = load(member_path)["members"]
    mac_zip = verify_zip(bundles["macos"], "artifact.zip.parts", CONFIG["macos_zip"],
                         members, CONFIG["macos_members"])
    expected = direct(bundles["macos"], str(mac.relative_to(ROOT)))
    extras = {"artifact/" + name: value for name, value in members.items()}
    extras.update({"archive-parts/" + Path(name).name: value for name, value in mac_zip["parts"].items()})
    extras["gate-diagnostics.zip"] = CONFIG["macos_zip"]
    require(not set(expected) & set(extras), "direct/ZIP-expanded file collision")
    expected.update(extras)
    roots[str(mac.relative_to(ROOT))] = compare(mac, expected, CONFIG["macos_special_entry_bound"])
    require(set(roots) == set(CONFIG["targets"]), "fixed retirement set incomplete")
    return {"evidence": evidence, "roots": roots, "Windows_ZIP": win_zip, "macOS_ZIP": mac_zip}


def recheck_linux_birth(pid, saved):
    require(type(pid) is int and pid > 0 and str(saved).isdigit(), "invalid Linux birth")
    try:
        current = mature.birth(pid)
    except (FileNotFoundError, ProcessLookupError):
        current = None
    require(current != int(saved), "an exact observed Linux birth remains")
    return {"pid": pid, "saved_start_ticks": str(saved), "current_start_ticks": current,
            "original_birth_absent": True}


def exact_negative_timeout(path, row):
    member = "Sirius/Sirius/out/macos-moltenvk-integration/observation/owner-self-test/timeout-tree/owner.json"
    expected_path = TARGETS[4] / "artifact" / member
    pin = {"bytes": 4311, "sha256": "a8ee7a0ce145a0e46697be9d57d6e784d316d5b01240073ce8570003cdd498ef"}
    bundle = ROOT / CONFIG["evidence"]["macos"][2]
    archived = bundle / "selected-artifact" / member
    require(path == expected_path and identity(path) == pin and identity(archived) == pin
            and load(bundle / "manifest.json")["files"]["selected-artifact/" + member] == pin,
            "unrecognized Darwin observation-incomplete owner")
    require(row["observation_complete"] is False and row["qualification_claimed"] is False
            and row["child_exit"] == -15 and row["stop_reason"] == "deadline"
            and row["limits"] == {"seconds": 1, "rss_kib": 4194304, "reserve_kib": 2097152}
            and row["source"] == {"runner_revision": "5b5fad91c08c0533e29c4f5c7ae249725162e405",
                                  "candidate_revision": "5b5fad91c08c0533e29c4f5c7ae249725162e405", "status": ""}
            and len(row["observed_births"]) == 2,
            "exact negative timeout self-test fields changed")
    return {"member": member, "identity": pin, "intentional_negative_control": True,
            "observation_complete": False, "qualification_claimed": False,
            "observation_gap": "Preserved one-second timeout self-test reports incomplete observation; "
                               "its two recorded births are absent, without a current remote-host or global census."}


def terminal_records():
    result = []
    boot = Path("/proc/sys/kernel/random/boot_id").read_text().strip()
    for root in TARGETS:
        for path in sorted(root.rglob("owner.json")):
            row = load(path)
            require(row.get("cleanup_errors") == [], "owner has unresolved cleanup errors")
            item = {"path": str(path.relative_to(ROOT)), "identity": identity(path),
                    "status": row.get("status"), "child_exit": row.get("child_exit", row.get("returncode")),
                    "stop_reason": row.get("stop_reason"), "passed": row.get("passed")}
            if row.get("kind") == "bounded-native-Windows-control":
                require(row["status"] == "completed" and row["owned_process_absent"] is True
                        and row["test_executable_unchanged"] and row["source_build_gate_unchanged"],
                        "Windows recorded owner is not retired")
                item["scope"] = "Preserved native Windows PID+UTC owner retirement; no fresh/global NT census from WSL."
            elif "retirement_census_complete" in row:
                require(row["status"] == "terminal" and row["whole_inputs_unchanged"]
                        and row["containment_complete"]
                        and row["retirement_census_complete"] and row["remaining_observed"] == []
                        and row["observed_births"] and all(b["absent"] is True for b in row["observed_births"]),
                        "Darwin recorded owner is not sealed/retired")
                if row["observation_complete"] is not True:
                    item["negative_control"] = exact_negative_timeout(path, row)
                item["scope"] = "Preserved Darwin terminal lstart/PGID census; no current remote-host process census."
            else:
                require(row["status"] == "terminal" and row["source_unchanged"]
                        and row["whole_inputs_unchanged"] and row["remaining_owned_processes"] == []
                        and row["observed_births"] and all(b["absent"] is True for b in row["observed_births"]),
                        "Linux recorded owner is not sealed/retired")
                item["same_boot"] = row.get("boot_id") == boot
                if item["same_boot"]:
                    item["fresh_birth_checks"] = [recheck_linux_birth(b["pid"], b["saved_start_ticks"])
                                                  for b in row["observed_births"]]
            before, after = path.parent / "inputs-before.json", path.parent / "inputs-after.json"
            if before.exists() and after.exists():
                require(load(before) == load(after), "owner's before/after input rows differ")
            result.append(item)
    win = TARGETS[0]
    accepted = load(win / "terminal-inventory-accepted.json")
    require(accepted["pass"] is True, "Windows original terminal reader did not accept")
    bridges = []
    for name in ("inventory-bridge.json", "terminal-inventory-bridge.json"):
        row = load(win / name)
        require(row["returncode"] == 0 and row["linux_child_birth_absent"] is True
                and not row.get("outer_stop") and not row.get("errors"), "Windows bridge incomplete")
        bridges.append({"path": str((win / name).relative_to(ROOT)), "identity": identity(win / name),
                        "fresh_linux_birth": recheck_linux_birth(row["linux_child_pid"], row["linux_child_birth"])})
    finite = load(ROOT / CONFIG["evidence"]["windows"][2] / "result.json")
    require(finite["passed"] and finite["whole_inputs_unchanged"]
            and all(x["owned_birth_absent"] for x in finite["native_owner_and_child_retired"]),
            "Windows finite terminal preservation changed")
    jobs = []
    for name in ("linux", "macos"):
        root = TARGETS[3] if name == "linux" else TARGETS[4]
        row = load(root / "collection.json")
        require(row["status"] == "completed" and row["conclusion"] == "cancelled", "original CI stop changed")
        jobs.append({"path": str((root / "collection.json").relative_to(ROOT)), "identity": identity(root / "collection.json"),
                     "conclusion": row["conclusion"], "job_id": row["job_id"]})
    return {"owners": result, "Windows_bridges": bridges, "CI_stops": jobs,
            "Windows_original_terminal_visibility": finite["accessible_terminal_census"],
            "global_NT_or_Darwin_process_absence_claimed": False}


def future_source(source):
    require(mature.git("cat-file", "-e", CONFIG["canonical_CI_commit"] + "^{commit}") == b"",
            "recorded CI source commit is unavailable")
    records = {}
    for name in CONFIG["canonical_future_paths"]:
        path = ROOT / name
        text = path.read_text()
        require(all(relative not in text for relative in CONFIG["targets"]),
                "canonical future consumer directly references retired scratch")
        records[name] = identity(path)
    return {"expected_source": source, "recorded_CI_source_commit": CONFIG["canonical_CI_commit"],
            "canonical_consumers": records,
            "floor_input_scope_from_root": "immutable43ed/11d17 archives and canonical only",
            "scope": "Fixed canonical source check plus fresh accessible /proc census; not a generic dependency proof."}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--expected-source", required=True)
    parser.add_argument("--remove", action="store_true")
    parser.add_argument("--accepted-audit-sha256")
    args = parser.parse_args()
    require(len(args.expected_source) == 40 and all(c in "0123456789abcdef" for c in args.expected_source),
            "root must pass the actual full clean final source commit")
    for name, expected in load(HERE / "source-freeze.json")["pins"].items():
        require(identity(Path(name)) == expected, "retirement source changed")
    mature.ROOT, mature.SOURCE, mature.TARGETS = ROOT, args.expected_source, TARGETS
    mature.prefix_target = lambda path, targets=TARGETS: any(path == t or path.is_relative_to(t) for t in targets)
    mature.clean_source()
    previous = None
    if args.remove:
        previous = load(HERE / "audit.json")
        require(identity(HERE / "audit.json")["sha256"] == args.accepted_audit_sha256
                and previous["status"] == "verified-not-removed"
                and previous["source"] == args.expected_source, "exact reviewed audit required")
        require(not RECEIPT.exists(), "fresh cleanup receipt directory required")
    else:
        require(args.accepted_audit_sha256 is None and not (HERE / "audit.json").exists(), "fresh audit required")
    proof = preservation()
    dependencies = future_source(args.expected_source)
    terminal = terminal_records()
    census = mature.consumer_census()
    document = {"schema": "sirius-fixed-native-terminal-retirement-v1", "status": "verified-not-removed",
                "source": args.expected_source, "preservation": proof, "dependencies": dependencies,
                "terminal_records": terminal, "initial_census": census, "removed": [],
                "qualification_claimed": False, "global_consumer_absence_claimed": False}
    if not args.remove:
        mature.write(HERE / "audit.json", document)
        print(json.dumps({"audit": identity(HERE / "audit.json"), "roots": CONFIG["targets"],
                          "permission_gaps": len(census["permission_gaps"])}))
        return
    require(previous["preservation"]["roots"] == proof["roots"], "root bytes/modes/layout changed since reviewed audit")
    mature.clean_source()
    document["final_census"] = mature.consumer_census()
    RECEIPT.mkdir(parents=True)
    for name in ("retire.py", "mature_cleanup.py", "configuration.json", "source-freeze.json", "commands.json", "audit.json"):
        shutil.copy2(HERE / name, RECEIPT / name)
    mature.write(RECEIPT / "receipt.json", document)
    for root in TARGETS:
        require(root in TARGETS and root.is_dir() and not root.is_symlink() and root.resolve() == root,
                "fixed root changed immediately before removal")
        relative = str(root.relative_to(ROOT))
        bound = CONFIG["macos_special_entry_bound"] if root == TARGETS[4] else 10000
        require(compare(root, proof["roots"][relative]["files"], bound) == proof["roots"][relative],
                "fixed root bytes/modes/layout changed immediately before removal")
        shutil.rmtree(root)
        require(not root.exists(), "retired root remains")
        document["removed"].append(str(root.relative_to(ROOT)))
        mature.write(RECEIPT / "receipt.json", document)
    mature.clean_source()
    document["status"] = "preserved-terminal-roots-removed"
    document["source_unchanged"] = True
    mature.write(RECEIPT / "receipt.json", document)
    mature.write(RECEIPT / "manifest.json", {"files": {p.name: identity(p) for p in sorted(RECEIPT.iterdir()) if p.is_file()}})
    print(json.dumps({"receipt": identity(RECEIPT / "receipt.json"), "removed": document["removed"],
                      "Git_commit_tag_push_performed": False}))


if __name__ == "__main__":
    main()
