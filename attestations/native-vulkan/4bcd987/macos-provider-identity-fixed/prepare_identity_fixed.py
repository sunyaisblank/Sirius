"""Prepare one metadata-only diagnostic provider payload; no compiler execution."""
from pathlib import Path
import difflib
import hashlib
import importlib.util
import json
import os
import shutil
import stat
import sys
import tarfile

sys.dont_write_bytecode = True
WORK = Path(__file__).resolve().parent
REPO = WORK.parents[1]
assert WORK == REPO / "out/macos-moltenvk-identity-fixed"
assert {p.name for p in WORK.iterdir()} == {"prepare_identity_fixed.py"}
spec = importlib.util.spec_from_file_location("original_payload_reader", REPO / "scripts/macos_moltenvk_integration.py")
reader = importlib.util.module_from_spec(spec)
spec.loader.exec_module(reader)

def encoded(value):
    return (json.dumps(value, sort_keys=True, indent=2) + "\n").encode()

def identity(path):
    data = Path(path).read_bytes()
    return {"bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}

def write(path, value):
    path.write_bytes(encoded(value))

original, unused_inputs = reader.stage_payload(WORK)
payload = WORK / "payload"
original.rename(payload)
old_freeze = json.loads((payload / "input-freeze.json").read_text())
old_content = json.loads((payload / "content-id.json").read_text())
historical = payload / "historical/original"
historical.mkdir(parents=True)
device = "sources/moltenvk-patched/MoltenVK/MoltenVK/GPUObjects/MVKDevice.mm"
patch = "patches/MoltenVK-diagnostic-identity.patch"
changed = [device, patch, "build-profile.json", "content-id.json"]
for name in ["input-freeze.json", "preparation.json", *changed]:
    target = historical / name
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(payload / name, target)

# Legacy preparation/runner files remain preserved evidence, not a runnable
# recipe for this new identity. Translation receipts are retained in place.
moved = ["preparation.json"]
(payload / "preparation.json").unlink()
for name in ["prepare.py", "check_translation.py", "run_root_translation.py",
             "root_translation_phase.py", "root_candidate_phase.py", "run_owned.py",
             "owned_guard.py", "actions.json"]:
    if (payload / name).exists():
        shutil.move(payload / name, historical / name)
        moved.append(name)

old_line = '\tstrcpy(supportedProps12.driverInfo, MVK_VERSION_STRING "+sirius-issue48-unforced." MVK_DIAGNOSTIC_CONTENT_ID);'
new_line = ('\tstd::snprintf(supportedProps12.driverInfo, sizeof(supportedProps12.driverInfo),\n'
            '\t\t\t\t  "%s+sirius-issue48-unforced.%s", MVK_VERSION_STRING, MVK_DIAGNOSTIC_CONTENT_ID);')
text = (payload / device).read_text()
assert text.count(old_line) == 1 and text.count("#include <cmath>\n") == 1
assert "#include <cstdio>" not in text
text = text.replace("#include <cmath>\n", "#include <cmath>\n#include <cstdio>\n")
text = text.replace(old_line, new_line)
(payload / device).write_text(text)

# Recreate the complete identity patch against the exact upstream archive,
# rather than treating the earlier broken patch as an upstream source base.
archive_pin = old_content["preimage"]["archives"]["MoltenVK"]
archive = next(path for path in (payload / "archives").iterdir()
               if path.is_file() and identity(path) == {k: archive_pin[k] for k in ("bytes", "sha256")})
diff = []
with tarfile.open(archive, "r:gz") as bundle:
    for relative in ["MoltenVK/CMakeLists.txt", "MoltenVK/MoltenVK/GPUObjects/MVKDevice.mm"]:
        members = [m for m in bundle.getmembers() if m.isfile() and m.name.split("/", 1)[-1] == relative]
        assert len(members) == 1
        before = bundle.extractfile(members[0]).read().decode()
        after = (payload / "sources/moltenvk-patched" / relative).read_text()
        diff.extend(difflib.unified_diff(before.splitlines(keepends=True), after.splitlines(keepends=True),
                                         fromfile="a/" + relative, tofile="b/" + relative))
(payload / patch).write_text("".join(diff))

profile = json.loads((historical / "build-profile.json").read_text())
profile.update(schema=2, cmake_generator="Unix Makefiles",
    native_toolchain="Root must separately seal selected Xcode, macOS SDK, clang, CMake and xcrun-selected make before build",
    inherited_python_environment={"PYTHONDONTWRITEBYTECODE": "1"},
    build_only_environment_removed=["DYLD_PRINT_LIBRARIES"],
    runtime_environment={"DYLD_PRINT_LIBRARIES": "1"},
    cache_scope={"TMPDIR": "fresh owned observation/tmp", "XDG_CACHE_HOME": "fresh owned observation/cache",
                 "CLANG_MODULE_CACHE_PATH": "fresh owned observation/module-cache",
                 "CPM_SOURCE_CACHE": "fresh owned observation/cpm-cache with pinned CPM0.40.8",
                 "claim": "Declared diagnostic output placement; no cache-hit, warmness or performance claim"},
    metadata_format="Bounded std::snprintf; parenthesized MVK_VERSION_STRING is a %s argument")
write(payload / "build-profile.json", profile)
preimage = dict(old_content["preimage"])
preimage["build_profile"] = profile
preimage["patches"] = {p.name: identity(p) for p in sorted((payload / "patches").iterdir())}
content_id = hashlib.sha256(encoded(preimage)).hexdigest()
assert content_id[:8] not in {"00000000", "db445ff2", old_content["content_id"][:8]}
label = "1.4.1+sirius-issue48-unforced." + content_id
assert len(label.encode()) + 1 <= 256
new_content = dict(old_content, content_id=content_id, preimage=preimage,
                   driver_info_expected=label, cache_uuid_expected_first_four_bytes_hex=content_id[:8])
write(payload / "content-id.json", new_content)

# Exact whitelist against the recovered original map. Unmodified translation
# sources, dependencies, driver/options, captured inputs and MSL remain exact.
translation_prefixes = ("sources/cross-stock/", "sources/cross-patched/", "sources/Vulkan-Headers/",
                        "sources/SPIRV-Headers/", "sources/SPIRV-Tools/", "sources/cereal/", "captured/")
translation_records = {}
differences = {}
for name, pin in old_freeze["files"].items():
    current = historical / name if name in moved else payload / name
    expected = {k: pin[k] for k in ("bytes", "sha256")}
    actual = identity(current)
    if name in changed:
        assert actual != expected, name
        assert identity(historical / name) == expected
        differences[name] = {"original": expected, "current": actual}
    else:
        assert actual == expected and stat.S_IMODE(current.stat().st_mode) == pin["mode"], name
    if name.startswith(translation_prefixes) or name == "translate_initialize.cpp":
        translation_records[name] = expected
assert set(differences) == set(changed)
for name, target in old_freeze["symlinks"].items():
    path = payload / name
    assert path.is_symlink() and os.readlink(path) == target and path.resolve().is_relative_to(payload)
owner_reuse = {}
for phase in ("stock-translation", "candidate-translation"):
    path = payload / phase / "owner.json"
    row = json.loads(path.read_text())
    assert row["status"] == "terminal" and row["passed"] and row["child_exit"] == 0
    assert row["whole_inputs_unchanged"] and row["source_unchanged"]
    assert row["stop_reason"] is None and not row["cleanup_errors"] and not row["remaining_owned_processes"]
    assert all(b["absent"] is True for b in row["observed_births"])
    owner_reuse[phase] = {"path": str(path.relative_to(payload)), **identity(path)}
acceptance = payload / "observations/translation-byte-acceptance.json"
receipt = json.loads(acceptance.read_text())
assert receipt["passed"] and receipt["input_freeze"] == identity(historical / "input-freeze.json")
assert identity(payload / "observations/stock.metal") == {k: receipt["stock"][k] for k in ("bytes", "sha256")}
assert identity(payload / "observations/candidate.metal") == {k: receipt["candidate"][k] for k in ("bytes", "sha256")}
reuse = {"schema": 1, "decision": "REUSE_EXACT_ORIGINAL_TRANSLATION_ONLY_NATIVE_BUILD_PENDING",
    "original_evidence_commit": reader.EVIDENCE_COMMIT, "original_payload_tree": reader.PAYLOAD_TREE,
    "original_input_freeze": {"path": "historical/original/input-freeze.json", **identity(historical / "input-freeze.json")},
    "original_translation_receipt": {"path": str(acceptance.relative_to(payload)), **identity(acceptance)},
    "original_translation_owners": owner_reuse, "compiler_launched_for_new_payload": False,
    "provider_launched_for_new_payload": False, "fresh_translation_claimed": False,
    "qualification_claimed": False, "metadata_only_input_differences": differences,
    "unchanged_translation_input_files": len(translation_records),
    "unchanged_translation_input_map_sha256": hashlib.sha256(encoded(translation_records)).hexdigest(),
    "stock_msl": receipt["stock"], "candidate_msl": receipt["candidate"],
    "relocated_legacy_files": {name: "historical/original/" + name for name in moved},
    "scope": "The four whitelisted source/patch/profile identity files changed; original preparation metadata is archived and new preparation metadata is derived separately. The shared non-entry emitter, complete Cross sources/dependencies, translation API driver/options, captured SPV and both original MSL outputs remain exact. Original translation owners/receipt belong to their original freeze, not the new payload."}
write(payload / "translation-reuse.json", reuse)
write(payload / "preparation.json", {"schema": 2, "status": "prepared_without_build_or_execution",
    "product_source": reader.common.SOURCE, "payload_origin_commit": reader.EVIDENCE_COMMIT,
    "payload_origin_tree": reader.PAYLOAD_TREE, "root": str(payload), "content_id": content_id,
    "driver_info_expected": label, "driver_info_chars": len(label), "driver_info_bytes_including_nul": len(label.encode()) + 1,
    "identity_patch": identity(payload / patch), "translation_reuse": identity(payload / "translation-reuse.json"),
    "original_source_map": identity(historical / "input-freeze.json"),
    "new_compile_or_provider_execution": False, "qualification_claimed": False,
    "source_changes": {device: "Explicit cstdio include and bounded snprintf driverInfo only; no shader or numerical code changed"},
    "future_native_owner": {"build_seconds": 900, "runtime_case_seconds": 300,
                             "rss_kib": 4194304, "reserve_kib": 2097152, "file_bytes": 67108864},
    "root_operator_obligations": ["Pin new evidence commit/tree/freeze/content ID and verify this historical translation reuse", "Seal selected native tool footprint and honor build profile", "Query actual stock/diagnostic UUIDs and loaded provider identities", "Run unchanged original numerical controls only after completed provider identity acceptance"]})

# The new freeze includes every file and link in the new payload, including all
# historical receipts and maps. It excludes only itself to avoid self-reference.
files, links = {}, {}
for base, directories, names in os.walk(payload, followlinks=False):
    for name in [*directories, *names]:
        path = Path(base) / name
        relative = str(path.relative_to(payload))
        if relative == "input-freeze.json":
            continue
        if path.is_symlink():
            assert path.resolve().is_relative_to(payload)
            links[relative] = os.readlink(path)
        elif path.is_file():
            files[relative] = {**identity(path), "mode": stat.S_IMODE(path.stat().st_mode)}
write(payload / "input-freeze.json", {"schema": 2, "files": files, "symlinks": links,
    "scope": "Complete new prepared payload including historical evidence; no compiler/provider execution. Only this freeze file is excluded from its own map."})
write(WORK / "prepared-payload.json", {"payload": str(payload), "content_id": content_id,
    "input_freeze": identity(payload / "input-freeze.json"), "files": len(files), "symlinks": len(links),
    "preparation_recipe": identity(Path(__file__)), "original_reader": {"path": str(reader.HERE), **identity(reader.HERE)},
    "no_execution": True, "qualification_claimed": False})
print((WORK / "prepared-payload.json").read_text())
