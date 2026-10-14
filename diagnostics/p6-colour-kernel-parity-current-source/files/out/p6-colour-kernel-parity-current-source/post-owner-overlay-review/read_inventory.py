"""Separate post-owner reader for this one preserved OverlayFS join failure.

No query, provider, test or original helper is run or rewritten. Original
validation logic is evaluated with one explicitly scoped device-join change.
"""
from pathlib import Path
import hashlib
import json
import os
import sys

sys.dont_write_bytecode = True
REVIEW = Path(__file__).resolve().parent
WORK = REVIEW.parent
ROOT = WORK.parents[1]
ORIGINAL = REVIEW / "original"
MANIFEST_SHA256 = "1723dc5d1faa20a40b7e10e467993ca1737be34a41eae0a8367046c236a4dd52"
JOIN_SHA256 = "5b66baf1f250c2548a67d4deda40961e11f172bf47461439052e492fc0846fd9"
CHECKER_SHA256 = "302c65a9955acc1484711b847790b89fd043c9908258fcb380d59bea7103fbb1"
THREE = {
    "/usr/lib/wsl/lib/libd3d12.so": 562949954106596,
    "/usr/lib/wsl/lib/libd3d12core.so": 562949954106597,
    "/usr/lib/wsl/lib/libdxcore.so": 562949954106598,
}


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def identity(path):
    path = Path(path)
    return {"resolved_path": str(path.resolve(strict=True)),
            "bytes": path.stat().st_size, "sha256": digest(path)}


def preserved_inputs():
    assert Path.cwd() == ROOT
    assert digest(ORIGINAL / "manifest.json") == MANIFEST_SHA256
    assert digest(REVIEW / "overlay-join.json") == JOIN_SHA256
    manifest = json.loads((ORIGINAL / "manifest.json").read_text())
    for relative, expected in manifest["files"].items():
        path = ORIGINAL / "files" / relative
        assert not path.is_symlink() and path.is_file()
        assert path.stat().st_size == expected["bytes"] and digest(path) == expected["sha256"]
        assert path.stat().st_mode & 0o7777 == expected["mode"]
        # Pin every original query/raw file and immutable helper at its origin.
        if not relative.startswith("cache/"):
            origin = WORK / relative
            assert identity(origin) == {"resolved_path": expected["original_resolved_path"],
                                        "bytes": expected["bytes"], "sha256": expected["sha256"]}
    pins = json.loads((REVIEW / "reader-pins.json").read_text())
    for relative, expected in pins.items():
        assert identity(REVIEW / relative) == expected
    primary = json.loads((REVIEW / "primary-sources.json").read_text())
    assert primary["commit"] == "427645e3db3a8896714f22a3d3fe0c3f7b317ad4"
    for item in primary["sources"]:
        path = REVIEW / "primary" / item["source_path"]
        assert path.stat().st_size == item["bytes"] and digest(path) == item["sha256"]
    join = json.loads((REVIEW / "overlay-join.json").read_text())
    assert join["boot_id"] == Path("/proc/sys/kernel/random/boot_id").read_text().strip()
    assert join["kernel_release"] == Path("/proc/sys/kernel/osrelease").read_text().strip() == "6.6.87.2-microsoft-standard-WSL2"
    assert join["overlay_xino_auto"] == Path("/sys/module/overlay/parameters/xino_auto").read_text().strip() == "N"
    assert digest(REVIEW / "mountinfo.txt") == join["mountinfo_sha256"]
    mounted = [line for line in Path("/proc/self/mountinfo").read_text().splitlines()
               if line.split()[4] == "/usr/lib/wsl/lib"]
    assert mounted == [join["overlay_mount_line"]]
    assert " - overlay " in mounted[0]
    directory = Path("/usr/lib/wsl/lib").stat()
    assert [os.major(directory.st_dev), os.minor(directory.st_dev)] == join["overlay_directory_stat_device"] == [0, 100]
    assert set(join["files"]) == set(THREE)
    return join


def main():
    assert len(sys.argv) == 1
    join = preserved_inputs()
    substitutions = []

    def device_identity_matches(path, parts, stat):
        mapped = tuple(int(x, 16) for x in parts[3].split(":"))
        visible = (os.major(stat.st_dev), os.minor(stat.st_dev))
        if mapped == visible:
            return True
        # Only this observed pseudo-device/backing-device pair, exact paths,
        # exact inodes, unchanged whole bytes and frozen read-only mmap join.
        assert path in THREE, ("unexplained device mismatch", path, mapped, visible)
        item = join["files"][path]
        assert mapped == tuple(item["backing_mapping_device"]) == (0, 35)
        assert visible == tuple(item["visible_stat_device"]) == (0, 102)
        assert stat.st_ino == int(parts[4]) == item["visible_inode"] == THREE[path]
        assert identity(path) == item["visible_identity"]
        assert item["mapping_sha256"] == item["visible_identity"]["sha256"]
        assert item["mapping_closed_and_absent"] is True
        observed = item["readonly_cpu_mapping"].split(None, 5)
        assert observed[5] == path and observed[1] == "r--p" and int(observed[2], 16) == 0
        assert tuple(int(x, 16) for x in observed[3].split(":")) == mapped
        assert int(observed[4]) == stat.st_ino
        substitutions.append({"path": path, "maps_device": list(mapped),
                              "visible_stat_device": list(visible), "inode": stat.st_ino})
        return True

    original = ORIGINAL / "files/check.py"
    assert digest(original) == CHECKER_SHA256
    source = original.read_text()
    needle = '            assert tuple(int(x, 16) for x in parts[3].split(":")) == (os.major(stat.st_dev), os.minor(stat.st_dev))'
    replacement = '            assert _device_identity_matches(path, parts, stat)'
    assert source.count(needle) == 1
    adapted = source.replace(needle, replacement)
    namespace = {"__name__": "preserved_original_p6_reader", "__file__": str(WORK / "check.py"),
                 "_device_identity_matches": device_identity_matches}
    exec(compile(adapted, str(original), "exec"), namespace)
    # The original owner is already terminal. All other original reader checks,
    # source/pin joins, actual maps, raw XML assertions and birth checks remain.
    result = namespace["validate"]("inventory")
    assert set(x["path"] for x in substitutions) == set(THREE)
    preserved_inputs()
    result.update(phase="separate_post_owner_overlay_analysis", qualification_claimed=False,
                  original_checker_reported_failure_preserved=True,
                  original_checker_sha256=CHECKER_SHA256,
                  original_preservation_manifest_sha256=MANIFEST_SHA256,
                  readonly_backing_join_sha256=JOIN_SHA256,
                  evaluated_reader_sha256=hashlib.sha256(adapted.encode()).hexdigest(),
                  exact_device_join_substitutions=substitutions,
                  scope="Inventory identity/containment and narrowly explained OverlayFS mapping join only; no P6 numerical execution or provider qualification.")
    target = REVIEW / "analysis-inventory.json"
    with target.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
