"""Check fixed preserved scratch; --remove also retire it and record an evidence tag.

Default is read-only. Never runs a test, compiler, validator or provider. The
/proc census covers accessible local observations and explicitly retains gaps.
"""
from pathlib import Path, PurePosixPath
import argparse
import hashlib
import json
import os
import re
import shutil
import stat
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parent.parent
SOURCE = "f5d649bf5067d1d68dc9215462f37830122f0540"
BUNDLE = ROOT / "attestations/cleanup/2026-10-15/terminal-isolated-parity-preflight-hints"
TAG = "evidence/issue-53-terminal-parity-preflight-hints-cleanup"
EVIDENCE = (
    ("evidence/issue-53-p4-ordinary-fp32-isolated-partial", "0d8beec5fdb1406d35add385cfa6f2d9f62848d2",
     "attestations/scientific-core/4bcd987/p4-ordinary-fp32-isolated-partial",
     "e5dac01afef4b81465f9d5e552de70cacb8e168d34d486e9537dad14d664ac78",
     ("out/science-p4-ordinary-fp32-parity-isolated",)),
    ("evidence/issue-48-native-preflight-order", "15366ab2b1039094a41ce517618fe120fb3fb22c",
     "attestations/native-vulkan/preflight-order-54c0f3b",
     "356748f9b6cde132452b4bba2f79b55b10f5231e626b5b1b01d44e7a17081fd4",
     ("out/native-preflight-order",)),
    ("evidence/issue-53-p4-unroll-hint-clear-setup", "d623208d82a0b9f739beaf69efa916382e366cbe",
     "attestations/scientific-core/4bcd987/p4-unroll-hint-clear-setup",
     "ecc099dbf33dff5e741740780e9fce64d0b1481cf270182c004e08b4ab56013f",
     ("out/p4-parity-unroll-hint-setup", "bin/linux-gcc/p4-parity-unroll-hint-setup")),
    ("evidence/issue-53-p4-dont-unroll-setup", "5bca4cc9bce771bda5702e69d381fdef94f5023b",
     "attestations/scientific-core/4bcd987/p4-dont-unroll-setup",
     "82b78eefbf38411b243405d78d2bbddefde953d8ec1a1557ab76cc0ed5de5a55",
     ("out/p4-parity-dont-unroll-setup", "bin/linux-gcc/p4-parity-dont-unroll-setup")),
)
TARGETS = tuple(ROOT / name for row in EVIDENCE for name in row[4])
FUTURE_PINS = {'out/native-provider-integration/.github/workflows/ci.yml': {'bytes': 33642, 'sha256': '8c629e05c0f166bc8b3795478abd808c587f28cfb8dfd09f6d431a7e7eb19146'}, 'out/native-provider-integration/README.md': {'bytes': 31137, 'sha256': '753020dd0db592897cf36eabce8cb42c15a656c1e6b8ae149c8ca230805993c4'}, 'out/native-provider-integration/candidate.json': {'bytes': 3703, 'sha256': '19d15ebd4b8767b37198cf21360b0ab121d2001ef5f288ac0fa598a0d2c3feea'}, 'out/native-provider-integration/commit_source.py': {'bytes': 2744, 'sha256': 'c0cb4747cbfb566d1a1c7cb989c74245e509e6b63b21117dc677439c7ea43cd6'}, 'out/native-provider-integration/scripts/diagnose-macos-retained.py': {'bytes': 46569, 'sha256': '8998d73985e3c330f4d7b55d5a421c33c06ceff5cfa160cda6ea54bb97cd9034'}, 'out/native-provider-integration/scripts/governance/build_policy.py': {'bytes': 71323, 'sha256': '841560051d11252e75aa1c0e3e2973e401edb1456eb53437af7df1bdf500bc67'}, 'out/native-provider-integration/scripts/macos_moltenvk_provider.py': {'bytes': 22818, 'sha256': '6736af1ebf985e2c72d0a616acfa9bd9bc37a528f86ecd7d370efbd8d2ed330c'}, 'out/native-provider-integration/scripts/validate-native-runtime.py': {'bytes': 22648, 'sha256': '1e57273bb962e94783c77d39a29d535d32b5252728de4b8fb6e34ce68fb25871'}, 'out/native-provider-integration/source-commit.json': {'bytes': 5698, 'sha256': '1fd20581a5e7af7435ff7ca1217db6fd80575309ac509b2e2848233369686d40'}}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def git(*arguments, environment=None):
    return subprocess.check_output(["git", *arguments], cwd=ROOT, env=environment, timeout=30)


def identity(path):
    before = path.stat()
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    after = path.stat()
    require((before.st_size, before.st_mtime_ns) == (after.st_size, after.st_mtime_ns), "input changed while hashing")
    return {"bytes": after.st_size, "sha256": digest.hexdigest()}


def clean_source():
    require(Path.cwd().resolve() == ROOT and not os.environ.get("GIT_INDEX_FILE"), "canonical checkout and default index required")
    require(git("rev-parse", "HEAD").decode().strip() == SOURCE
            and not git("status", "--porcelain").strip(), "pinned clean source changed")


def files_under(directory):
    require(directory.is_dir() and not directory.is_symlink() and directory.resolve() == directory,
            "owned/preserved root missing, aliased or replaced: " + str(directory))
    files = {}
    visited = 0
    for base, directories, names in os.walk(directory, followlinks=False):
        for name in [*directories, *names]:
            path = Path(base) / name
            visited += 1
            require(visited <= 10000, "fixed cleanup inventory exceeds its entry bound")
            info = path.lstat()
            require(not stat.S_ISLNK(info.st_mode), "unexpected symlink: " + str(path))
            require(stat.S_ISREG(info.st_mode) or stat.S_ISDIR(info.st_mode), "special entry: " + str(path))
            if stat.S_ISREG(info.st_mode):
                files[str(path.relative_to(directory))] = path
    return files


def remote_tag(tag, expected):
    require(git("rev-parse", "refs/tags/" + tag + "^{commit}").decode().strip() == expected, "local evidence tag changed")
    rows = git("ls-remote", "--tags", "origin", "refs/tags/" + tag, "refs/tags/" + tag + "^{}").decode().splitlines()
    refs = dict(row.split()[::-1] for row in rows)
    actual = refs.get("refs/tags/" + tag + "^{}", refs.get("refs/tags/" + tag))
    require(actual == expected, "pushed evidence tag missing/different: " + tag)
    return {"tag": tag, "commit": expected, "remote_origin_tag_verified": True}


def preserved(row):
    tag, commit, relative, manifest_sha, targets = row
    remote = remote_tag(tag, commit)
    bundle = ROOT / relative
    manifest_path = bundle / "manifest.json"
    require(identity(manifest_path)["sha256"] == manifest_sha, "preservation manifest changed")
    raw = files_under(bundle)
    entries = {}
    for item in git("ls-tree", "-rz", commit + ":" + relative).split(b"\0"):
        if item:
            metadata, name = item.split(b"\t", 1)
            mode, kind, blob = metadata.split()
            require(kind == b"blob" and mode in (b"100644", b"100755"), "unexpected preserved Git entry")
            entries[os.fsdecode(name)] = blob.decode()
    require(set(raw) == set(entries), "complete raw/Git preservation inventory differs")
    for name, path in raw.items():
        data = path.read_bytes()
        require(hashlib.sha1(b"blob " + str(len(data)).encode() + b"\0" + data).hexdigest() == entries[name],
                "preserved raw bytes differ from immutable Git: " + name)
    manifest = json.loads(manifest_path.read_text())
    for name, pin in manifest["files"].items():
        parts = PurePosixPath(name).parts
        require(parts and not name.startswith("/") and all(part not in (".", "..", ".git") for part in parts)
                and name in raw, "unsafe/missing manifest file")
        actual = identity(raw[name])
        require(all(actual[key] == pin[key] for key in ("bytes", "sha256")), "whole-byte manifest differs: " + name)
    scratch = []
    for name in targets:
        directory = ROOT / name
        local = files_under(directory)
        prefix = "files/" + name + "/"
        expected = {path[len(prefix):]: pin for path, pin in manifest["files"].items() if path.startswith(prefix)}
        require(local and set(local) == set(expected), "complete scratch/preserved file inventory differs: " + name)
        for filename, path in local.items():
            actual = identity(path)
            require(all(actual[key] == expected[filename][key] for key in ("bytes", "sha256")),
                    "scratch differs from preserved bytes: " + str(path))
        scratch.append({"path": name, "regular_files": len(local),
                        "regular_bytes": sum(pin["bytes"] for pin in expected.values())})
    return {**remote, "bundle": relative, "manifest": identity(manifest_path),
            "whole_manifest_files_verified": len(manifest["files"]), "whole_raw_Git_files_verified": len(raw),
            "scratch": scratch}


def prefix_target(path, targets=TARGETS):
    return any(path == target or path.is_relative_to(target) for target in targets)


def declared_dependencies():
    # The 759-input DontUnroll packet is now terminal and is itself being
    # retired. Only its own two roots are exempt; its inputs must avoid every
    # other cleanup root. Historical metadata strings are not live inputs.
    work = ROOT / "out/p4-parity-dont-unroll-setup"
    prepared = json.loads((work / "preparation.json").read_text())
    actions = json.loads((work / "actions.json").read_text())
    require(len(prepared["pins"]) == 759, "reviewed DontUnroll input set changed")
    other_targets = tuple(path for path in TARGETS if path not in (work, ROOT / "bin/linux-gcc/p4-parity-dont-unroll-setup"))
    inputs = list(prepared["pins"])
    inputs += [pin["resolved_path"] for pin in prepared["pins"].values()]
    inputs += [name for action in actions.values() for name in (*action["frozen_inputs"], *action["argv"])]
    for name in inputs:
        path = Path(name)
        path = path if path.is_absolute() else ROOT / path
        require(not prefix_target(path, other_targets) and not prefix_target(path.resolve(), other_targets),
                "DontUnroll declared input depends on another cleanup root: " + name)
    # Bind the reviewed future native source, not its changing CI/API output.
    # The new provider uses committed preflight code and Git parent identity;
    # no live out/native-preflight-order path is consumed.
    for name, expected in FUTURE_PINS.items():
        path = ROOT / name
        require(identity(path) == expected, "reviewed future native source changed: " + name)
        text = path.read_text()
        for target in TARGETS:
            relative = str(target.relative_to(ROOT))
            require(re.search(re.escape(relative) + r"(?:[/\s\"']|$)", text) is None,
                    "future native source refers to a cleanup root: " + name)
    return {"terminal_DontUnroll_pins": 759, "terminal_DontUnroll_actions": sorted(actions),
            "terminal_packet_self_roots_retired": True, "other_cleanup_inputs_absent": True,
            "reviewed_future_native_source": FUTURE_PINS, "native_preflight_scratch_is_live_input": False,
            "scope": "Declared pin/action paths plus exact reviewed native source; not a generic program-dependency proof"}


def birth(pid):
    raw = (Path("/proc") / str(pid) / "stat").read_text()
    return int(raw[raw.rindex(")") + 2:].split()[19])


def retired_owners():
    boot = Path("/proc/sys/kernel/random/boot_id").read_text().strip()
    result = []
    for target in TARGETS:
        for path in sorted(target.rglob("owner.json")):
            owner = json.loads(path.read_text())
            require(owner["status"] == "terminal" and owner["source_unchanged"] and owner["whole_inputs_unchanged"]
                    and not owner["cleanup_errors"] and owner["remaining_owned_processes"] == []
                    and owner["observed_births"] and all(row["absent"] is True for row in owner["observed_births"]),
                    "owner is not terminal/sealed/retired: " + str(path))
            for row in owner["observed_births"]:
                saved = row.get("saved_start_ticks")
                require(type(row.get("pid")) is int and row["pid"] > 0
                        and ((type(saved) is int and saved >= 0)
                             or (type(saved) is str and re.fullmatch(r"[0-9]+", saved) is not None)),
                        "owned birth PID/start-ticks token is malformed")
                saved = int(saved)
                if owner.get("boot_id") == boot:
                    try:
                        current = birth(row["pid"])
                    except (FileNotFoundError, ProcessLookupError):
                        continue
                    require(current != saved, "an original owned birth is still present")
            result.append({"path": str(path.relative_to(ROOT)), "identity": identity(path),
                           "stop_reason": owner.get("stop_reason"), "child_exit": owner.get("child_exit"),
                           "passed": owner.get("passed"), "same_boot_births_rechecked": owner.get("boot_id") == boot})
    return result


def consumer_census():
    consumers, gaps = [], []
    absolute = [re.compile(re.escape(os.fsencode(path)) + rb"(?:/|\x00|\s|$|[\"'])") for path in TARGETS]
    relative = [re.compile(rb"(?<![A-Za-z0-9_./-])(?:\./)*" + re.escape(os.fsencode(path.relative_to(ROOT)))
                           + rb"(?:/|\x00|\s|$|[\"'])") for path in TARGETS]
    def hits(data, include_relative=False):
        return any(pattern.search(data) for pattern in absolute + (relative if include_relative else []))
    for entry in Path("/proc").glob("[0-9]*"):
        pid = int(entry.name)
        if pid == os.getpid():
            continue
        observed, values = [], {}
        def read(label, operation):
            try:
                value = operation()
                values[label] = value
                if hits(value, label == "cmdline"):
                    observed.append(label)
            except (FileNotFoundError, ProcessLookupError):
                pass
            except PermissionError:
                gaps.append({"pid": pid, "field": label, "error": "PermissionError"})
        read("cmdline", lambda: (entry / "cmdline").read_bytes())
        read("maps", lambda: (entry / "maps").read_bytes())
        for label in ("cwd", "exe"):
            read(label, lambda label=label: os.fsencode(os.readlink(entry / label)))
        try:
            for descriptor in (entry / "fd").iterdir():
                read("fd/" + descriptor.name, lambda descriptor=descriptor: os.fsencode(os.readlink(descriptor)))
        except (FileNotFoundError, ProcessLookupError):
            pass
        except PermissionError:
            gaps.append({"pid": pid, "field": "fd", "error": "PermissionError"})
        # Resolve standalone path arguments against the observed cwd as well:
        # ./out, ../Sirius/out and symlink spellings must not evade the census.
        if "cwd" in values and "cmdline" in values:
            cwd = Path(os.fsdecode(values["cwd"]))
            for argument in values["cmdline"].split(b"\0"):
                token = os.fsdecode(argument).split("=", 1)[-1]
                if token and not token.startswith("-") and not any(char.isspace() for char in token):
                    candidate = Path(token)
                    candidate = candidate if candidate.is_absolute() else cwd / candidate
                    if prefix_target(candidate.resolve()):
                        observed.append("resolved-cmdline-path")
        if observed:
            consumers.append({"pid": pid, "fields": sorted(set(observed))})
    require(not consumers, "accessible local consumers still reference cleanup roots: " + repr(consumers))
    return {"accessible_local_consumers": consumers, "permission_gaps": gaps,
            "global_consumer_absence_claimed": False,
            "scope": "Sampled accessible local cmdline (absolute/relative prefixes), maps, cwd, exe and fd links; permission/race/remote-host visibility is not global absence"}


def write(path, document):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(document, sort_keys=True, indent=2) + "\n")
    temporary.replace(path)


def record_tag():
    fd, index = tempfile.mkstemp(dir=ROOT / "out", prefix="terminal-cleanup-index-")
    os.close(fd)
    os.unlink(index)
    environment = dict(os.environ, GIT_INDEX_FILE=index)
    try:
        subprocess.run(["git", "read-tree", "--empty"], cwd=ROOT, env=environment, check=True, timeout=30)
        rows = []
        for path in sorted(BUNDLE.iterdir()):
            blob = git("hash-object", "-w", "--", str(path)).decode().strip()
            rows.append(f"100644 {blob}\t{path.relative_to(ROOT)}\0")
        subprocess.run(["git", "update-index", "-z", "--index-info"], cwd=ROOT, env=environment,
                       input="".join(rows).encode(), check=True, timeout=30)
        tree = git("write-tree", environment=environment).decode().strip()
        commit = subprocess.check_output(["git", "commit-tree", tree, "-p", SOURCE], cwd=ROOT, env=environment,
                    input=b"Record cleanup of preserved terminal parity and preflight diagnostics\n", timeout=30).decode().strip()
        subprocess.run(["git", "tag", TAG, commit], cwd=ROOT, check=True, timeout=30)
        for path in BUNDLE.iterdir():
            require(git("show", commit + ":" + str(path.relative_to(ROOT))) == path.read_bytes(), "cleanup Git evidence differs")
        return {"tag": TAG, "commit": commit, "tree": tree, "pushed": False}
    finally:
        if Path(index).exists():
            Path(index).unlink()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--remove", action="store_true", help="retire only the verified fixed roots and create the local isolated evidence tag")
    args = parser.parse_args()
    clean_source()
    require(not BUNDLE.exists(), "cleanup receipt bundle already exists")
    require(subprocess.run(["git", "show-ref", "--verify", "--quiet", "refs/tags/" + TAG], cwd=ROOT).returncode == 1,
            "cleanup evidence tag already exists or could not be inspected")
    preservation = [preserved(row) for row in EVIDENCE]
    dependencies = declared_dependencies()
    owners = retired_owners()
    census = consumer_census()
    receipt = {"schema": "sirius-preserved-terminal-scratch-cleanup-v2", "runner_revision": SOURCE,
               "status": "verified-not-removed", "preservation": preservation, "dependencies": dependencies,
               "terminal_owners": owners, "initial_consumer_census": census, "removed": [],
               "qualification_claimed": False, "global_consumer_absence_claimed": False}
    if not args.remove:
        print(json.dumps(receipt, sort_keys=True))
        return
    clean_source()
    receipt["final_consumer_census"] = consumer_census()
    BUNDLE.mkdir(parents=True)
    shutil.copy2(Path(__file__).resolve(), BUNDLE / Path(__file__).name)
    write(BUNDLE / "receipt.json", receipt)
    for row in preservation:
        for item in row["scratch"]:
            target = ROOT / item["path"]
            require(target in TARGETS and target.is_dir() and not target.is_symlink() and target.resolve() == target,
                    "fixed cleanup root changed before removal")
            shutil.rmtree(target)
            require(not target.exists(), "cleanup root remains")
            receipt["removed"].append(item)
            write(BUNDLE / "receipt.json", receipt)
    clean_source()
    receipt["status"] = "preserved-terminal-roots-removed"
    receipt["source_unchanged"] = True
    write(BUNDLE / "receipt.json", receipt)
    write(BUNDLE / "manifest.json", {"qualification_claimed": False,
          "files": {path.name: identity(path) for path in sorted(BUNDLE.iterdir()) if path.is_file()}})
    recorded = record_tag()
    print(json.dumps({**recorded, "receipt": identity(BUNDLE / "receipt.json"),
                      "removed": receipt["removed"], "permission_gap_count": len(census["permission_gaps"])}))


if __name__ == "__main__":
    main()
