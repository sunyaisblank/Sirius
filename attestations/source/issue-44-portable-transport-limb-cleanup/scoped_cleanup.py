"""Remove only this finished trial after exact recovery and a fresh local census."""
from pathlib import Path
import hashlib
import json
import os
import shutil
import stat
import subprocess
import time

root = Path.cwd()
archive = Path(__file__).resolve().parent
manifest = json.loads((archive / "whole-byte-manifest.json").read_text())
targets = [root / name for name in manifest["roots"]]
assert [str(p.relative_to(root)) for p in targets] == [
    "out/portable-transport-limb-accumulation", "bin/linux-gcc/portable-transport-limb-accumulation"
]
assert all(p.is_dir() and not p.is_symlink() for p in targets)
root_births = {str(p): (p.stat().st_dev, p.stat().st_ino) for p in targets}


def topology():
    entries = {}
    for target in targets:
        assert not target.is_symlink() and target.is_dir()
        assert (target.stat().st_dev, target.stat().st_ino) == root_births[str(target)]
        for path in target.rglob("*"):
            mode = path.lstat().st_mode
            assert not stat.S_ISLNK(mode), path
            assert stat.S_ISREG(mode) or stat.S_ISDIR(mode), path
            entries[str(path.relative_to(root))] = "file" if stat.S_ISREG(mode) else "directory"
    return entries


initial_topology = topology()
assert not subprocess.check_output(["git", "status", "--porcelain"])
assert subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip() == (
    "a05763a68fa371e927eda1b7f83cc97c98465323"
)
assert not subprocess.check_output(["git", "diff", "c8199bd", "HEAD"])
assert json.loads((archive / "recorded-input-recovery.json").read_text())["passed"]
bindings = json.loads((archive / "source-git-bindings.json").read_text())
remote = subprocess.check_output([
    "git", "ls-remote", "--tags", "origin", "refs/tags/" + bindings["tag"] + "^{}"
], text=True)
assert remote.split()[0] == bindings["commit"]

assert subprocess.check_output(["git", "rev-list", "--left-right", "--count", "HEAD...@{upstream}"]).strip() == b"0\t0"
branch = subprocess.check_output(["git", "branch", "--show-current"], text=True).strip()
remote_head = subprocess.check_output(["git", "ls-remote", "--heads", "origin", "refs/heads/" + branch], text=True)
assert remote_head.split()[0] == "a05763a68fa371e927eda1b7f83cc97c98465323"
guard_bindings = json.loads((archive / "cleanup-source-git-bindings.json").read_text())
guard_record = next(v for v in guard_bindings["aliases"] if v["alias"].endswith("/scoped_cleanup.py"))
guard_data = Path(__file__).read_bytes()
assert len(guard_data) == guard_record["bytes"] and hashlib.sha256(guard_data).hexdigest() == guard_record["sha256"]
remote_guard = subprocess.check_output(["git", "ls-remote", "--tags", "origin", "refs/tags/" + guard_bindings["tag"] + "^{}"], text=True)
assert remote_guard.split()[0] == guard_bindings["commit"]


def identity(path):
    data = path.read_bytes()
    return {"bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}


files = {str(p.relative_to(root)): p for t in targets for p in t.rglob("*") if p.is_file()}
assert set(files) == set(manifest["files"])
for name, value in manifest["files"].items():
    assert not files[name].is_symlink()
    assert identity(files[name]) == value, name
    assert identity(archive / "files" / name) == value, name
protected = json.loads((targets[0] / "protected-before.json").read_text())
for name, value in protected.items():
    assert identity(root / name) == value, name
windows = [
    "bin/windows-msvc/native-b81e061-executed-cost",
    "bin/windows-msvc/native-24a8e14-rhs-layer-clock",
]
assert all((root / p).is_dir() for p in windows)
births = set()


def collect(value):
    if isinstance(value, dict):
        ticks = value.get("start_ticks", value.get("saved_start_ticks"))
        if value.get("pid") is not None and ticks is not None:
            births.add((int(value["pid"]), str(ticks)))
        if value.get("child_pid") is not None and value.get("child_start_ticks") is not None:
            births.add((int(value["child_pid"]), str(value["child_start_ticks"])))
        for nested in value.values():
            collect(nested)
    elif isinstance(value, list):
        for nested in value:
            collect(nested)


boot = Path("/proc/sys/kernel/random/boot_id").read_text().strip()
for receipt in targets[0].rglob("owner.json"):
    value = json.loads(receipt.read_text())
    assert value["boot_id"] == boot
    assert value["status"] == "terminal" and not value["remaining_owned_processes"]
    collect(value)
for sample in targets[0].rglob("samples.jsonl"):
    for line in sample.read_text().splitlines():
        if line:
            collect(json.loads(line))
fresh = []
for pid, saved in sorted(births):
    try:
        proc_stat_text = Path(f"/proc/{pid}/stat").read_text()
        current = proc_stat_text[proc_stat_text.rfind(")") + 2:].split()[19]
    except (FileNotFoundError, ProcessLookupError):
        current = None
    assert current != saved, (pid, saved, "saved birth still exists")
    fresh.append({"pid": pid, "saved_start_ticks": saved,
                  "current_start_ticks": current, "absent": True})
ignored = {os.getpid(), os.getppid()}
consumers = []
inaccessible = []
examined = 0
names = [str(p).encode() for p in targets]
for proc in Path("/proc").glob("[0-9]*"):
    pid = int(proc.name)
    if pid in ignored:
        continue
    fields = {}
    for field in ("cmdline", "maps", "cwd", "exe", "fd"):
        try:
            if field in ("cmdline", "maps"):
                fields[field] = (proc / field).read_bytes()
            elif field in ("cwd", "exe"):
                fields[field] = os.readlink(proc / field).encode()
            else:
                descriptors = []
                for fd in (proc / "fd").iterdir():
                    try:
                        descriptors.append(os.readlink(fd).encode())
                    except (FileNotFoundError, ProcessLookupError):
                        pass
                    except PermissionError:
                        inaccessible.append({"pid": pid, "field": "fd-link"})
                fields[field] = b"\n".join(descriptors)
        except (FileNotFoundError, ProcessLookupError):
            pass
        except PermissionError:
            inaccessible.append({"pid": pid, "field": field})
    if fields:
        examined += 1
        hits = [key for key, data in fields.items() if any(name in data for name in names)]
        if hits:
            consumers.append({"pid": pid, "fields": hits})
preflight = {"checked_at_epoch": time.time(), "exact_raw_file_set_and_hashes": True,
             "raw_files": len(files), "protected_files_exact": len(protected),
             "fresh_saved_birth_fields": fresh, "accessible_process_fields_examined": examined,
             "inaccessible_process_fields": inaccessible, "ignored_current_cleanup_pids": sorted(ignored),
             "exact_target_consumers": consumers, "windows_stages_keep": windows,
             "scope": "Fresh saved births and accessible exact-target fields; no global, hidden or scheduled consumer absence claim."}
(archive / "cleanup-preflight.json").write_text(json.dumps(preflight, indent=2) + "\n")
assert not consumers, consumers
# Close the gap between recovery verification and the first unlink. The owned
# executions are terminal; reject any replacement or intervening local change.
assert topology() == initial_topology
final_files = {name for name, kind in initial_topology.items() if kind == "file"}
assert final_files == set(manifest["files"])
for name, value in manifest["files"].items():
    assert identity(root / name) == value, name
    assert identity(archive / "files" / name) == value, name
    assert oct((root / name).stat().st_mode & 0o777) == manifest["preserved_file_modes"][name]
for name, value in protected.items():
    assert identity(root / name) == value, name
assert not subprocess.check_output(["git", "status", "--porcelain"])
assert subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip() == "a05763a68fa371e927eda1b7f83cc97c98465323"
assert not subprocess.check_output(["git", "diff", "c8199bd", "HEAD"])
assert all((root / p).is_dir() for p in windows)
counts = {str(t.relative_to(root)): sum(p.is_file() for p in t.rglob("*")) for t in targets}
for target in targets:
    shutil.rmtree(target)
assert all(not target.exists() for target in targets)
for name, value in protected.items():
    assert identity(root / name) == value, name
assert all((root / p).is_dir() for p in windows)
receipt = {"passed": True, "removed_roots": counts,
           "immediate_pre_unlink_topology_bytes_modes_protected_source_recheck": True,
           "fresh_saved_birth_fields_absent": len(fresh), "protected_files_exact": len(protected),
           "out_empty": not any((root / "out").iterdir()), "windows_stages_preserved": windows,
           "archive_manifest_sha256": identity(archive / "whole-byte-manifest.json")["sha256"]}
(archive / "cleanup-receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
print(json.dumps(receipt))
