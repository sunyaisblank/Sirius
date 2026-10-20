"""Retire exactly four published/quiescent diagnostic input roots.

ROOT explicitly authorized this finite retirement. Only local Git blob reads,
one fresh ls-remote, file verification, and listed-file/empty-dir removal occur.
No branch/index/ref writes, workloads, provider or process queries.
"""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import os
import stat
import subprocess
import time

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
SOURCE = "5d993d12b3fde7a65eaa5f7988b73b86b3a9a291"
DRIVER = "1cb6e6aaff04ddd54ff7ee133b773257338f2c52"
SOURCE_ROOTS = ["out/native-wave32-source", "out/native-wave32-producer-inputs", "out/native-wave32-windows-source"]
DRIVER_ROOT = "out/native-shader-driver-introspection"
ROOTS = SOURCE_ROOTS + [DRIVER_ROOT]
BUNDLE = "attestations/qualification/2c8ef086/native-wave32-driver-insight"
SOURCE_REF = "refs/heads/development/native-wave32-diagnostic-2026-10-18"
DRIVER_REF = "refs/tags/evidence/native-wave32-driver-insight-2c8ef086"
RECEIPT = OUT / "published-inputs-retirement.json"


def require(value, message):
    if not value:
        raise ValueError(message)


def utc():
    return datetime.now(timezone.utc).isoformat()


def git(*arguments):
    return subprocess.check_output(["git", *arguments], cwd=ROOT)


def dump(path, value):
    with path.open("x", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2)
        stream.write("\n")
    path.chmod(0o644)


def identity(path):
    info = path.lstat()
    require(stat.S_ISREG(info.st_mode) and not path.is_symlink(), "Nonregular/symlink file: " + str(path))
    value = hashlib.sha256(path.read_bytes()).hexdigest()
    return {"bytes": info.st_size, "sha256": value, "mode": format(stat.S_IMODE(info.st_mode), "04o"),
            "device": info.st_dev, "inode": info.st_ino, "mtime_ns": info.st_mtime_ns,
            "hardlink_count": info.st_nlink}


def physical(relative):
    base = ROOT / relative
    require(base.exists() and base.is_dir() and not base.is_symlink(), "Missing/aliased root: " + relative)
    require(base.resolve() == ROOT.resolve() / relative, "Resolved root differs: " + relative)
    files, directories = {}, {}
    for directory, names, leaves in os.walk(base, followlinks=False):
        parent = Path(directory)
        info = parent.lstat()
        require(stat.S_ISDIR(info.st_mode) and not parent.is_symlink(), "Non-directory/symlink parent")
        directories[str(parent.relative_to(ROOT))] = {
            "mode": format(stat.S_IMODE(info.st_mode), "04o"), "device": info.st_dev, "inode": info.st_ino}
        for name in names:
            child = parent / name
            require(not child.is_symlink() and child.is_dir(), "Symlink/non-directory child")
        for name in leaves:
            path = parent / name
            files[str(path.relative_to(ROOT))] = identity(path)
    return files, directories


def tree(revision, prefix):
    result = {}
    for record in git("ls-tree", "-r", "-z", revision, "--", prefix).split(b"\0"):
        if not record:
            continue
        header, raw_path = record.split(b"\t", 1)
        mode, kind, blob = header.decode("ascii").split()
        require(kind == "blob" and mode in ("100644", "100755"), "Nonregular named Git entry")
        result[raw_path.decode("utf-8")] = {"Git_mode": mode, "Git_blob": blob}
    require(result, "Empty named Git prefix: " + prefix)
    return result


def main():
    require(not RECEIPT.exists(), "Refuse existing retirement receipt")
    started = utc()
    inventories, states = {}, {}
    manifest_bytes = git("show", DRIVER + ":" + BUNDLE + "/MANIFEST.json")
    manifest = json.loads(manifest_bytes)
    manifest_file = OUT / "driver-published-manifest.json"
    require(not manifest_file.exists(), "Refuse prior manifest copy")
    manifest_file.write_bytes(manifest_bytes)
    manifest_file.chmod(0o644)
    driver_tree = tree(DRIVER, BUNDLE)
    for relative in ROOTS:
        record = {"root": relative, "verified": False, "removed": False, "mismatch": None}
        files, directories = {}, {}
        try:
            files, directories = physical(relative)
            if relative in SOURCE_ROOTS:
                revision, expected = SOURCE, tree(SOURCE, relative)
            else:
                revision = DRIVER
                expected = {}
                for path, row in manifest["payload_files"].items():
                    if not path.startswith(DRIVER_ROOT + "/"):
                        continue
                    payload = row["payload"]
                    require(not Path(payload).is_absolute() and ".." not in Path(payload).parts, "Unsafe published payload")
                    stored = driver_tree[BUNDLE + "/" + payload]
                    expected[path] = {**stored, "manifest_row": row, "Git_path": BUNDLE + "/" + payload}
            record.update({"file_count": len(files), "bytes": sum(row["bytes"] for row in files.values()),
                           "directory_count": len(directories), "named_Git_revision": revision})
            require(set(files) == set(expected), "File inventory mismatch; missing=" +
                    repr(sorted(set(expected) - set(files))) + "; extra=" + repr(sorted(set(files) - set(expected))))
            for path, row in files.items():
                authority = expected[path]
                data = git("cat-file", "blob", authority["Git_blob"])
                require(len(data) == row["bytes"] and hashlib.sha256(data).hexdigest() == row["sha256"],
                        "Physical bytes differ from actual named Git blob: " + path)
                require(int(row["mode"], 8) & 0o111 == int(authority["Git_mode"], 8) & 0o111,
                        "Physical executable bits differ from named Git mode: " + path)
                if relative == DRIVER_ROOT:
                    m = authority["manifest_row"]
                    require(row["bytes"] == m["bytes"] and row["sha256"] == m["sha256"] and
                            int(row["mode"], 8) == int(m["mode"], 8), "Driver full mode/identity differs from published manifest")
                row.update({"named_Git_revision": revision, "Git_blob": authority["Git_blob"],
                            "Git_mode": authority["Git_mode"], "Git_path": authority.get("Git_path", path),
                            "actual_Git_blob_bytes_matched": True})
            record["verified"] = True
            inventories[relative] = {"files": files, "directories": directories}
        except (ValueError, KeyError, OSError, subprocess.CalledProcessError) as error:
            record["mismatch"] = str(error)
            inventories[relative] = {"files": locals().get("files", {}), "directories": locals().get("directories", {})}
        states[relative] = record
        dump(OUT / (relative.split("/")[-1] + "-physical-inventory.json"), inventories[relative])
    remote_started = utc()
    begin = time.monotonic()
    process = subprocess.run(["git", "ls-remote", "--refs", "origin", SOURCE_REF, DRIVER_REF],
                             cwd=ROOT, capture_output=True, text=True, timeout=60)
    rows = {}
    for line in process.stdout.splitlines():
        value, name = line.split()
        require(name not in rows, "Duplicate remote ref")
        rows[name] = value
    remote = {"started_utc": remote_started, "finished_utc": utc(), "elapsed_seconds": time.monotonic() - begin,
              "argv": ["git", "ls-remote", "--refs", "origin", SOURCE_REF, DRIVER_REF],
              "exit_code": process.returncode, "stdout": process.stdout, "stderr": process.stderr,
              "refs": rows, "passed": process.returncode == 0 and rows == {SOURCE_REF: SOURCE, DRIVER_REF: DRIVER}}
    dump(OUT / "published-inputs-remote-verification.json", remote)
    if remote["passed"]:
        for relative in ROOTS:
            record = states[relative]
            if not record["verified"]:
                continue
            try:
                before = inventories[relative]
                files, directories = physical(relative)
                require(set(files) == set(before["files"]), "Fresh file membership changed before deletion")
                expected_files = {key: {field: row[field] for field in files[key]}
                                  for key, row in before["files"].items()}
                require(files == expected_files and directories == before["directories"], "Fresh physical identity/membership changed before deletion")
                record["removed_file_count"] = 0
                record["removed_byte_count"] = 0
                for path in sorted(files):
                    require(identity(ROOT / path) == files[path], "File changed at removal boundary: " + path)
                    (ROOT / path).unlink()
                    record["removed_file_count"] += 1
                    record["removed_byte_count"] += files[path]["bytes"]
                for path in sorted(directories, key=lambda name: (name.count("/"), name), reverse=True):
                    (ROOT / path).rmdir()
                require(not (ROOT / relative).exists(), "Retired root remains")
                record["removed"] = True
            except (ValueError, KeyError, OSError) as error:
                record["mismatch"] = str(error)
    result = {"schema": "sirius-native-wave32-published-inputs-retirement-v1",
              "passed": remote["passed"] and all(row["verified"] and row["removed"] for row in states.values()),
              "started_utc": started, "finished_utc": utc(), "roots": states,
              "removed_files": sum(row.get("removed_file_count", 0) for row in states.values()),
              "removed_bytes": sum(row.get("removed_byte_count", 0) for row in states.values()),
              "remote_proof": "out/native-wave32-closeout/published-inputs-remote-verification.json",
              "cleanup_source": {"path": str(Path(__file__).relative_to(ROOT)), **identity(Path(__file__))},
              "canonical_new_candidate_and_other_roots_mutated": False,
              "already_retired_320_files_132789368_bytes_b567_repeated": False,
              "Git_branch_index_ref_or_source_mutation": False,
              "ownership_basis": "ROOT confirmed all completed runtimes closed, no live consumers, and new candidate consumers use immutable named Git; no new process/provider query performed.",
              "scope": "Four explicitly named published input roots only. Physical bytes/modes/inventory and actual named Git blob conservation plus fresh remote refs verified before removal; no source adoption or qualification claim."}
    dump(RECEIPT, result)
    print(json.dumps({"passed": result["passed"], "removed_files": result["removed_files"], "removed_bytes": result["removed_bytes"],
                      "receipt": str(RECEIPT.relative_to(ROOT)), "sha256": identity(RECEIPT)["sha256"]}))
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
