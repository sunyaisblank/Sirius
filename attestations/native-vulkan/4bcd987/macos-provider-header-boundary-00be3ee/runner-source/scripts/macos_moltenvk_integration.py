#!/usr/bin/env python3
"""One pinned diagnostic MoltenVK build and four unchanged Vulkan controls.

Uses the existing Darwin owner and original executable preparation. This never
writes a qualification receipt or presents the custom provider as stock SDK.
"""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path, PurePosixPath
import platform
import posixpath
import re
import shutil
import stat
import struct
import subprocess
import sys
import tarfile

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve()
REPO = HERE.parent.parent if HERE.parent.name == "scripts" else HERE.parents[2]
sys.path.insert(0, str(REPO / "scripts"))
spec = importlib.util.spec_from_file_location("sirius_existing_darwin_owner", REPO / "scripts/diagnose-macos-retained.py")
common = importlib.util.module_from_spec(spec)
spec.loader.exec_module(common)

EVIDENCE_COMMIT = "21f87e4806c36cd1ae9e45f912a5817d6878b774"
PAYLOAD_TREE = "c2c8004cd36b91b34d4259f335b3bc6cc2adf941"
PAYLOAD_PATH = "attestations/native-vulkan/4bcd987/macos-provider-identity-fixed/payload"
FREEZE_SHA = "a0b5986daccff1453e11361feb8b8a5a020a6e83c1943deec4da943f3e4ca8d4"
CONTENT_ID = "76cb8b212bc049261146b526d75ebddafe017b18fe1f4c2c7065c3024c4ffe9e"
ORIGINAL_EVIDENCE_COMMIT = "e2422dd92017879a452377a1122f1838e6610b3c"
ORIGINAL_PAYLOAD_TREE = "8b7d54847738c427954c7c3b070708d60fe4c8a8"
ORIGINAL_FREEZE_SHA = "71ca1547b3fb616ac99b9e8871a7e2410f25443131a6652b8b1f1c3464f526c9"
REUSE_SHA = "18b318c2793c2c688fd6405cddb863a6f3d67c7fbc4d0f6148180637bd07c62c"
IDENTITY_PATCH_SHA = "4e5150ecebeda5584e42ea98f3b168d45752781726dc9114c6bac2e63009097c"
DIAGNOSTIC_INFO = "1.4.1+sirius-issue48-unforced." + CONTENT_ID
CMAKE_ROOT_PROBE = 'message(STATUS "SIRIUS_CMAKE_ROOT=${CMAKE_ROOT}")\n'
STOCK_DRIVER = {"bytes": 10553504, "sha256": "3e2ba682b6006d20f6fd7233cffb826951fef17b98c22e7dedc0fff93f3fba9e"}
STOCK_LOADER = {"bytes": 1488704, "sha256": "0ed7750a2ee514df7a0daae66e7f2a98bfedca48b98c8f475b2ca8cc5f142094"}
INITIALIZE_MSL = {"bytes": 82487, "sha256": "80b33059934cc77a685c5738b9413c94379e45f1b18d0e6635a6432ff0407921"}
EXPECTED_MODULES = (
    "retained_program::kInitializePortableShader",
    "retained_program::kEndpointPortableNormalSumShader",
    "retained_program::kTransportPortableNormalSumShader",
    "retained_program::kDensePortableShader",
)
require, identity, write = common.require, common.identity, common.write


def git(*arguments):
    return subprocess.check_output(["git", *arguments], cwd=REPO, text=True).strip()


def stage_payload(output):
    require(re.fullmatch(r"[0-9a-f]{40}", EVIDENCE_COMMIT)
            and re.fullmatch(r"[0-9a-f]{40}", PAYLOAD_TREE)
            and not PAYLOAD_PATH.startswith("PENDING_"), "new provider evidence pins remain unresolved")
    require(git("cat-file", "-t", EVIDENCE_COMMIT) == "commit", "pinned evidence commit is unavailable")
    require(git("rev-parse", EVIDENCE_COMMIT + ":" + PAYLOAD_PATH) == PAYLOAD_TREE,
            "provider tree is not joined to the pinned evidence commit")
    archive = output / "provider-source.tar"
    with archive.open("xb") as stream:
        subprocess.run(["git", "archive", "--format=tar", PAYLOAD_TREE], cwd=REPO, stdout=stream, check=True)
    require(archive.stat().st_size < 512 * 1024**2, "provider source export exceeds bound")
    target = output / "provider-source"
    target.mkdir()
    with tarfile.open(archive, "r:") as bundle:
        members = bundle.getmembers()
        require(len(members) < 30000, "provider tree member count exceeds bound")
        paths, links = set(), {}
        for member in members:
            name = member.name.rstrip("/")
            parts = PurePosixPath(name).parts
            require(parts and not name.startswith("/") and "\\" not in name
                    and all(part not in ("..", ".", ".git") for part in parts), "unsafe Git export path")
            require(name not in paths and (member.isdir() or member.isfile() or member.issym()), "duplicate/special Git export member")
            paths.add(name)
            if member.issym():
                value = member.linkname
                joined = posixpath.normpath(posixpath.join(posixpath.dirname(name), value))
                require(value and not value.startswith("/") and "\\" not in value
                        and joined != ".." and not joined.startswith("../"), "escaping Git export link")
                links[name] = value
        require(all(not any(str(parent) in links for parent in PurePosixPath(name).parents) for name in paths),
                "Git export member descends through a symlink")
        for member in members:
            path = target / member.name
            if member.isdir():
                path.mkdir(parents=True, exist_ok=True)
            elif member.isfile():
                path.parent.mkdir(parents=True, exist_ok=True)
                with bundle.extractfile(member) as source, path.open("xb") as destination:
                    shutil.copyfileobj(source, destination)
                path.chmod(member.mode & 0o777)
        for name, value in links.items():
            path = target / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.symlink_to(value)
            require(path.resolve().is_relative_to(target), "resolved Git export link escapes")
    freeze = json.loads((target / "input-freeze.json").read_text())
    require(identity(target / "input-freeze.json")["sha256"] == FREEZE_SHA, "provider input freeze differs")
    for name, pin in freeze["files"].items():
        path = target / name
        require(path.resolve().is_relative_to(target) and path.is_file() and not path.is_symlink(), "frozen file missing/replaced")
        require(identity(path) == {key: pin[key] for key in ("bytes", "sha256")}, "Git payload byte join failed: " + name)
        # Git preserves executable status, not all POSIX bits. Restore only the
        # exact modes in the immutable original input map, then verify them.
        path.chmod(pin["mode"])
    source_files = verify_payload(target)
    reuse = verify_translation_reuse(target)
    write(output / "provider-source-join.json", {"evidence_commit": EVIDENCE_COMMIT, "payload_tree": PAYLOAD_TREE,
          "payload_path": PAYLOAD_PATH, "export": identity(archive), "input_freeze": identity(target / "input-freeze.json"),
          "mode_restoration": "Original frozen file modes restored from the pinned map after Git blob export",
          "translation_reuse": reuse,
          "scope": "New full payload freeze; old translation receipts retain their original freeze and no fresh translation is claimed",
          "content_id_meaning": "Diagnostic source/configuration cache discriminator, not Git or native-build provenance"})
    return target, source_files


def verify_payload(target):
    freeze = json.loads((target / "input-freeze.json").read_text())
    require(identity(target / "input-freeze.json")["sha256"] == FREEZE_SHA, "provider freeze changed")
    paths = [target / "input-freeze.json"]
    for name, pin in freeze["files"].items():
        path = target / name
        require(path.is_file() and not path.is_symlink() and path.resolve().is_relative_to(target), "frozen provider file missing/replaced")
        require(identity(path) == {key: pin[key] for key in ("bytes", "sha256")}
                and stat.S_IMODE(path.stat().st_mode) == pin["mode"], "frozen provider input changed: " + name)
        paths.append(path)
    for name, value in freeze["symlinks"].items():
        path = target / name
        require(path.is_symlink() and os.readlink(path) == value and path.resolve().is_relative_to(target),
                "provider source link changed")
    files, links, directories = set(), set(), set()
    for base, child_dirs, names in os.walk(target, followlinks=False):
        for name in [*child_dirs, *names]:
            path = Path(base) / name
            relative = str(path.relative_to(target))
            if path.is_symlink():
                links.add(relative)
            elif path.is_file():
                files.add(relative)
            else:
                require(path.is_dir(), "special provider input entry: " + relative)
                directories.add(relative)
    expected_files = set(freeze["files"]) | {"input-freeze.json"}
    expected_links = set(freeze["symlinks"])
    expected_dirs = {str(parent) for name in expected_files | expected_links
                     for parent in PurePosixPath(name).parents if str(parent) != "."}
    require(files == expected_files and links == expected_links and directories == expected_dirs,
            "complete provider payload gained/lost files, links or directories")
    return paths


def verify_translation_reuse(target):
    canonical = lambda value: (json.dumps(value, sort_keys=True, indent=2) + "\n").encode()
    historical = target / "historical/original"
    old_map_path = historical / "input-freeze.json"
    require(identity(old_map_path)["sha256"] == ORIGINAL_FREEZE_SHA, "original translation freeze differs")
    old_map = json.loads(old_map_path.read_text())
    reuse_path = target / "translation-reuse.json"
    require(identity(reuse_path)["sha256"] == REUSE_SHA, "translation reuse statement differs")
    reuse = json.loads(reuse_path.read_text())
    require(reuse["original_evidence_commit"] == ORIGINAL_EVIDENCE_COMMIT
            and reuse["original_payload_tree"] == ORIGINAL_PAYLOAD_TREE
            and reuse["original_input_freeze"] == {"path": "historical/original/input-freeze.json", **identity(old_map_path)}
            and reuse["compiler_launched_for_new_payload"] is False
            and reuse["provider_launched_for_new_payload"] is False
            and reuse["fresh_translation_claimed"] is False and reuse["qualification_claimed"] is False,
            "historical translation scope/origin differs")
    moved = ("preparation.json", "prepare.py", "check_translation.py", "run_root_translation.py",
             "root_translation_phase.py", "root_candidate_phase.py", "run_owned.py", "owned_guard.py", "actions.json")
    require(reuse["relocated_legacy_files"] == {name: "historical/original/" + name for name in moved},
            "legacy preparation/runner relocation differs")
    device = "sources/moltenvk-patched/MoltenVK/MoltenVK/GPUObjects/MVKDevice.mm"
    patch = "patches/MoltenVK-diagnostic-identity.patch"
    allowed = {device, patch, "build-profile.json", "content-id.json"}
    require(set(reuse["metadata_only_input_differences"]) == allowed, "metadata-only difference whitelist differs")
    unchanged = {}
    prefixes = ("sources/cross-stock/", "sources/cross-patched/", "sources/Vulkan-Headers/",
                "sources/SPIRV-Headers/", "sources/SPIRV-Tools/", "sources/cereal/", "captured/")
    for name, pin in old_map["files"].items():
        expected = {key: pin[key] for key in ("bytes", "sha256")}
        path = (historical if name in moved else target) / name
        if name in allowed:
            require(identity(historical / name) == expected
                    and reuse["metadata_only_input_differences"][name]
                    == {"original": expected, "current": identity(path)}, "metadata preimage/difference differs: " + name)
        else:
            require(identity(path) == expected and stat.S_IMODE(path.stat().st_mode) == pin["mode"],
                    "original non-metadata input changed: " + name)
        if name.startswith(prefixes) or name == "translate_initialize.cpp":
            unchanged[name] = expected
    for name, link in old_map["symlinks"].items():
        path = target / name
        require(path.is_symlink() and os.readlink(path) == link and path.resolve().is_relative_to(target),
                "original source link differs")
    require(len(unchanged) == 11475 and reuse["unchanged_translation_input_files"] == len(unchanged)
            and hashlib.sha256(canonical(unchanged)).hexdigest() == reuse["unchanged_translation_input_map_sha256"],
            "complete unchanged translation input join differs")
    before = (historical / device).read_text()
    old_line = '\tstrcpy(supportedProps12.driverInfo, MVK_VERSION_STRING "+sirius-issue48-unforced." MVK_DIAGNOSTIC_CONTENT_ID);'
    new_line = ('\tstd::snprintf(supportedProps12.driverInfo, sizeof(supportedProps12.driverInfo),\n'
                '\t\t\t\t  "%s+sirius-issue48-unforced.%s", MVK_VERSION_STRING, MVK_DIAGNOSTIC_CONTENT_ID);')
    require(before.count(old_line) == 1 and before.count("#include <cmath>\n") == 1,
            "original device source guards differ")
    expected_source = before.replace("#include <cmath>\n", "#include <cmath>\n#include <cstdio>\n").replace(old_line, new_line)
    require((target / device).read_text() == expected_source
            and identity(target / patch)["sha256"] == IDENTITY_PATCH_SHA,
            "metadata source change exceeds explicit include/snprintf replacement")
    content = json.loads((target / "content-id.json").read_text())
    profile = json.loads((target / "build-profile.json").read_text())
    old_content = json.loads((historical / "content-id.json").read_text())
    require(content["content_id"] == CONTENT_ID and content["driver_info_expected"] == DIAGNOSTIC_INFO
            and content["cache_uuid_expected_first_four_bytes_hex"] == CONTENT_ID[:8]
            and hashlib.sha256(canonical(content["preimage"])).hexdigest() == CONTENT_ID
            and content["preimage"]["archives"] == old_content["preimage"]["archives"]
            and content["preimage"]["build_profile"] == profile
            and content["preimage"]["patches"] == {path.name: identity(path) for path in sorted((target / "patches").iterdir())}
            and len(DIAGNOSTIC_INFO.encode()) + 1 <= 256,
            "canonical diagnostic content/patch/archive/profile identity differs")
    require(profile["cmake_generator"] == "Unix Makefiles" and profile["build_parallelism"] == 1
            and profile["inherited_python_environment"] == {"PYTHONDONTWRITEBYTECODE": "1"}
            and profile["build_only_environment_removed"] == ["DYLD_PRINT_LIBRARIES"]
            and profile["runtime_environment"] == {"DYLD_PRINT_LIBRARIES": "1"}
            and profile["public_metal_api_only"] is True and profile["upstream_git_attestation"] is False,
            "declared diagnostic build/runtime boundary differs")
    receipt_path = target / "observations/translation-byte-acceptance.json"
    receipt = json.loads(receipt_path.read_text())
    require(reuse["original_translation_receipt"] == {"path": str(receipt_path.relative_to(target)), **identity(receipt_path)}
            and receipt["passed"] is True and receipt["input_freeze"] == identity(old_map_path)
            and receipt["candidate"]["sha256"] == INITIALIZE_MSL["sha256"], "original translation receipt differs")
    for phase, key in (("stock-translation", "stock"), ("candidate-translation", "candidate")):
        folder = target / phase
        owner = json.loads((folder / "owner.json").read_text())
        require(owner["status"] == "terminal" and owner["passed"] is True and owner["child_exit"] == 0
                and owner["whole_inputs_unchanged"] and owner["source_unchanged"]
                and owner["stop_reason"] is None and not owner["remaining_owned_processes"]
                and not owner["cleanup_errors"] and owner["observed_births"]
                and all(row["absent"] is True for row in owner["observed_births"]), "historical translation owner incomplete")
        require(reuse["original_translation_owners"][phase]
                == {"path": phase + "/owner.json", **identity(folder / "owner.json")}, "historical owner pin differs")
        inputs = json.loads((folder / "inputs-before.json").read_text())
        require(inputs == json.loads((folder / "inputs-after.json").read_text()), "historical translation inputs changed")
        sealed_map = inputs["out/macos-moltenvk-unforced/input-freeze.json"]
        require({field: sealed_map[field] for field in ("bytes", "sha256")} == identity(old_map_path),
                "historical owner is not joined to original freeze")
        expected_output = {field: receipt[key][field] for field in ("bytes", "sha256")}
        require(identity(target / ("observations/" + key + ".metal")) == expected_output
                and reuse[key + "_msl"] == receipt[key], "original translation output differs")
    return {"statement": identity(reuse_path), "original_freeze": identity(old_map_path),
            "unchanged_translation_inputs": len(unchanged), "fresh_translation_claimed": False}


def complete_owner(owner, require_zero=True):
    require(owner["status"] == "terminal" and (not require_zero or owner["child_exit"] == 0)
            and owner["containment_complete"] and owner["observation_complete"]
            and owner["whole_inputs_unchanged"] and owner["stop_reason"] is None,
            "owned action failed/incomplete; sequence stopped")


def clean_environment(output, sdk):
    environment = {key: value for key, value in os.environ.items()
        if not key.startswith(("GTEST_", "SIRIUS_", "MVK_CONFIG_", "CPM_", "CMAKE_", "GIT_", "VK_LOADER_", "DYLD_"))
        and key not in {"VK_ICD_FILENAMES", "VK_ADD_DRIVER_FILES", "VK_INSTANCE_LAYERS", "VK_LAYER_PATH",
                        "VK_ADD_LAYER_PATH", "DYLD_INSERT_LIBRARIES", "CC", "CXX", "CFLAGS", "CXXFLAGS",
                        "CPPFLAGS", "LDFLAGS", "OBJCFLAGS", "OBJCXXFLAGS", "CPATH", "C_INCLUDE_PATH",
                        "CPLUS_INCLUDE_PATH", "OBJC_INCLUDE_PATH", "LIBRARY_PATH", "SDKROOT", "MACOSX_DEPLOYMENT_TARGET",
                        "MAKEFLAGS", "MFLAGS", "MAKELEVEL", "MAKEFILES", "GNUMAKEFLAGS"}}
    for name in ("tmp", "cache", "module-cache", "cpm-cache"):
        (output / name).mkdir()
    environment.update(LC_ALL="C", PYTHONDONTWRITEBYTECODE="1", DYLD_PRINT_LIBRARIES="1", MVK_CONFIG_DEBUG="0",
        MVK_CONFIG_TRACE_VULKAN_CALLS="3", VK_DRIVER_FILES=str(sdk / "share/vulkan/icd.d/MoltenVK_icd.json"),
        DYLD_LIBRARY_PATH=os.pathsep.join(str(path) for path in (sdk / "lib/VulkanLoader/lib", sdk / "lib") if path.is_dir()),
        TMPDIR=str(output / "tmp"), XDG_CACHE_HOME=str(output / "cache"),
        CLANG_MODULE_CACHE_PATH=str(output / "module-cache"), CPM_SOURCE_CACHE=str(output / "cpm-cache"),
        CMAKE_BUILD_PARALLEL_LEVEL="1", GIT_CEILING_DIRECTORIES=str(REPO))
    return environment


def toolchain_child():
    require(sys.platform == "darwin" and platform.machine() == "arm64", "native arm64 Darwin required")
    def capture(command):
        return subprocess.check_output(command, text=True, timeout=20).strip()
    paths = {"python": str(Path(sys.executable).resolve()), "xcrun": "/usr/bin/xcrun", "xcodebuild": "/usr/bin/xcodebuild"}
    for name in ("cmake",):
        value = shutil.which(name)
        require(value, "existing build tool missing: " + name)
        paths[name] = str(Path(value).resolve())
    for name in ("clang", "clang++", "ld", "make"):
        paths[name] = capture(["/usr/bin/xcrun", "--no-cache", "--sdk", "macosx", "--find", name])
    sdk = capture(["/usr/bin/xcrun", "--no-cache", "--sdk", "macosx", "--show-sdk-path"])
    resource = capture([paths["clang"], "-print-resource-dir"])
    # Ask the selected executable for its module root; package layouts differ.
    probe = REPO / "out/macos-moltenvk-integration/observation/discover-cmake-root.cmake"
    require(Path.cwd().resolve() == REPO / "out/macos-moltenvk-integration/observation/toolchain",
            "CMake discovery must stay in the owned observation directory")
    require(probe.read_text() == CMAKE_ROOT_PROBE, "CMake discovery script differs")
    answer = capture([paths["cmake"], "-P", str(probe)])
    marker = "-- SIRIUS_CMAKE_ROOT="
    require(answer.startswith(marker) and "\n" not in answer, "CMake module root query is ambiguous")
    cmake_root = Path(answer[len(marker):])
    require(cmake_root.is_absolute() and (cmake_root / "Modules").is_dir(), "CMake reported module root is absent")
    toolchain_usr = Path(paths["clang"]).parent.parent
    roots = [Path(sdk).resolve(), Path(resource).resolve(), toolchain_usr / "include", toolchain_usr / "lib", cmake_root]
    require(all(path.is_dir() for path in roots), "selected SDK/toolchain input root missing")
    print(json.dumps({"paths": paths, "sdk": str(Path(sdk).resolve()), "input_roots": list(dict.fromkeys(str(p.resolve()) for p in roots)),
          "cmake_root_query": {"script": {"path": str(probe), **identity(probe)}, "stdout": answer},
          "versions": {"xcode": capture(["/usr/bin/xcodebuild", "-version"]),
                       "macos_sdk": capture(["/usr/bin/xcrun", "--no-cache", "--sdk", "macosx", "--show-sdk-version"]),
                       "clang": capture([paths["clang"], "--version"]),
                       "cmake": capture([paths["cmake"], "--version"]),
                       "make": capture([paths["make"], "--version"])}}))


def tool_inputs(toolchain):
    files = {Path(value).resolve() for value in toolchain["paths"].values()}
    links = {}
    for root_name in toolchain["input_roots"]:
        root = Path(root_name)
        for base, directories, names in os.walk(root, followlinks=False, onerror=lambda error: (_ for _ in ()).throw(error)):
            for name in [*directories, *names]:
                path = Path(base) / name
                if path.is_symlink():
                    links[str(path)] = {"target": os.readlink(path), "resolved": str(path.resolve()), "exists": path.exists()}
                    if path.is_file():
                        files.add(path.resolve())
                elif path.is_file():
                    files.add(path.resolve())
                require(len(files) < 150000, "selected SDK/toolchain input count exceeds diagnostic bound")
    require(all(path.is_file() for path in files), "selected tool input missing")
    return sorted(files), links


def build_child(context_path, expected_sha):
    require(identity(context_path)["sha256"] == expected_sha, "build context differs")
    context = json.loads(context_path.read_text())
    payload, output = Path(context["payload"]), Path(context["output"])
    source = payload / "sources/moltenvk-patched"
    require(source.resolve().is_relative_to(REPO / "out") and context["content_id"] == CONTENT_ID,
            "build source/content identity differs")
    tools, build = context["toolchain"], Path(context["build"])
    require(build == REPO / "bin/macos/macos-moltenvk-integration" and not build.exists(), "build directory must be fresh")
    # Both cmake configure and build inherit a source cwd below out/. The outer
    # mature owner intentionally keeps its usual observation-folder cwd.
    # Single-job noninteractive Makefiles avoid Ninja's per-command groups;
    # the owner's strict escaped-descendant refusal remains unchanged.
    os.chdir(source)
    command = [tools["paths"]["cmake"], "-S", str(source), "-B", str(build), "-G", "Unix Makefiles",
        "-DCMAKE_BUILD_TYPE=Release", "-DCMAKE_OSX_ARCHITECTURES=arm64",
        "-DCMAKE_OSX_SYSROOT=" + tools["sdk"], "-DCMAKE_OSX_DEPLOYMENT_TARGET=11.0",
        "-DCMAKE_C_COMPILER=" + tools["paths"]["clang"], "-DCMAKE_CXX_COMPILER=" + tools["paths"]["clang++"],
        "-DCMAKE_OBJC_COMPILER=" + tools["paths"]["clang"], "-DCMAKE_OBJCXX_COMPILER=" + tools["paths"]["clang++"],
        "-DCMAKE_LINKER=" + tools["paths"]["ld"], "-DCMAKE_MAKE_PROGRAM=" + tools["paths"]["make"],
        "-DMOLTEN_VK_WITH_CCACHE=OFF", "-DMVK_BUILD_SHADER_CONVERTER_TOOL=OFF",
        "-DMVK_EXCLUDE_SPIRV_TOOLS=OFF", "-DMVK_EXCLUDE_CEREAL=OFF", "-DMVK_USE_METAL_PRIVATE_API=OFF",
        "-DMVK_CONFIG_LOG_LEVEL=error", "-DSPIRV_SKIP_EXECUTABLES=ON", "-DSPIRV_SKIP_TESTS=ON",
        "-DBUILD_TESTS=OFF", "-DSKIP_PORTABILITY_TEST=ON", "-DFETCHCONTENT_FULLY_DISCONNECTED=ON",
        "-DCPM_SOURCE_CACHE=" + str(output / "cpm-cache"), "-DCPM_USE_LOCAL_PACKAGES=OFF",
        "-DMVK_DIAGNOSTIC_CONTENT_ID=" + CONTENT_ID]
    for name, relative in (("SPIRV-Cross", "cross-patched"), ("Vulkan-Headers", "Vulkan-Headers"),
                           ("SPIRV-Headers", "SPIRV-Headers"), ("SPIRV-Tools", "SPIRV-Tools"), ("cereal", "cereal")):
        command.append("-DCPM_" + name + "_SOURCE=" + str(payload / "sources" / relative))
    print(json.dumps({"cwd": str(source), "configure_argv": command}), flush=True)
    subprocess.run(command, cwd=source, check=True)
    subprocess.run([tools["paths"]["cmake"], "--build", str(build), "--target", "MoltenVK", "--parallel", "1"],
                   cwd=source, check=True)


def loaded_join(folder, driver, loader):
    raw = (folder / "stderr.log").read_text(errors="strict")
    found = {}
    for name in ("libMoltenVK", "libvulkan"):
        paths = {Path(value).resolve() for value in re.findall(r"(/[^\n]*?/" + name + r"[^/\s]*\.dylib)", raw)}
        expected = driver if name == "libMoltenVK" else loader
        require(paths == {expected.resolve()}, "actual case/query loaded image differs: " + name)
        found[name] = {"path": str(expected.resolve()), **identity(expected)}
    return found


def shader_join(folder, modules, expected_module):
    generic = common.dump_join(folder / "shaders", modules)
    require(generic["all_dumped_spirv_bound"], "dumped SPIR-V missing/unbound")
    paired, found = [], False
    for path in sorted((folder / "shaders").glob("*.spv")):
        digest = identity(path)["sha256"]
        aliases = [name for name, pin in modules.items() if pin["sha256"] == digest]
        found |= expected_module in aliases
        msl = path.with_suffix(".metal")
        require(msl.is_file(), "dumped SPIR-V lacks its MSL pair")
        text = msl.read_text(errors="strict")
        exact_initialize = "retained_program::kInitializePortableShader" in aliases
        if exact_initialize:
            require(identity(msl) == INITIALIZE_MSL, "Initialize MSL does not equal the sealed annotation control")
        paired.append({"spirv": {"path": str(path), **identity(path)}, "embedded_aliases": aliases,
                       "msl": {"path": str(msl), **identity(msl)}, "exact_initialize_control": exact_initialize,
                       "forced_inline_declaration_count": text.count("static inline __attribute__((always_inline))"),
                       "static_inline_declaration_count": text.count("static inline"),
                       "scope": "Only Initialize has an exact prior-byte comparison; other stages join original SPIR-V, frozen emitter source and actual custom-provider capture. Generated custom helpers retain their original policy."})
    require(found, "expected original stage was not captured: " + expected_module)
    return {"generic": generic, "pairs": paired, "expected_module": expected_module}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--prepare-only", action="store_true")
    args = parser.parse_args()
    output = args.output.resolve()
    require(output == REPO / "out/macos-moltenvk-integration/observation" and not output.exists(), "fresh exact output root required")
    common.source_state()
    require(args.prepare_only or (sys.platform == "darwin" and platform.machine() == "arm64"), "native arm64 Darwin required")
    output.mkdir(parents=True)
    result = {"status": "incomplete", "candidate_source": common.SOURCE, "runner_source": common.source_state(),
              "provider_source_evidence_commit": EVIDENCE_COMMIT, "provider_payload_tree": PAYLOAD_TREE,
              "provider_content_discriminator": CONTENT_ID, "provider_scope": "Custom diagnostic MoltenVK; unchanged SDK is not qualified",
              "qualification_claimed": False, "production_adopted": False, "cases": []}
    payload = None
    try:
        stage, app, backend, modules = common.prepare(args.archive.resolve(), output)
        payload, source_inputs = stage_payload(output)
        header_path = payload / "sources/Vulkan-Headers/include/vulkan/vulkan_core.h"
        source_header_patches = re.findall(r"^#define\s+VK_HEADER_VERSION\s+(\d+)\s*$", header_path.read_text(), re.MULTILINE)
        require(source_header_patches == ["334"], "pinned source Vulkan header patch differs")
        stock_api = (1 << 22) | (3 << 12) | 334
        diagnostic_api = (1 << 22) | (3 << 12) | int(source_header_patches[0])
        result["expected_api_versions"] = {"stock": "1.3.334", "diagnostic": "1.3.334",
            "source_header": {"path": str(header_path), **identity(header_path)},
            "scope": "Original app and public query both request Vulkan1.3; pinned source and stock header patch are both334"}
        if args.prepare_only:
            result["status"] = "prepared_without_execution"
            return 0
        sdk = Path(os.environ["VULKAN_SDK"]).resolve()
        require("1.4.357.0" in sdk.parts, "pinned Vulkan SDK version required")
        stock_manifest = sdk / "share/vulkan/icd.d/MoltenVK_icd.json"
        stock_driver = (stock_manifest.parent / json.loads(stock_manifest.read_text())["ICD"]["library_path"]).resolve()
        loader = (sdk / "lib/libvulkan.1.dylib").resolve()
        require(identity(stock_driver) == STOCK_DRIVER and identity(loader) == STOCK_LOADER, "stock SDK binaries differ")
        environment = clean_environment(output, sdk)
        scripts = [HERE, HERE.with_name("macos_vulkan_identity.py"), REPO / "scripts/diagnose-macos-retained.py",
                   REPO / "scripts/runtime_identity.py", REPO / "scripts/macos_metal_probe.py",
                   REPO / "scripts/governance/build_policy.py", REPO / ".github/workflows/ci.yml"]
        python_inputs = [Path(module.__file__).resolve() for module in list(sys.modules.values())
                         if getattr(module, "__file__", None) and Path(module.__file__).is_file()]
        inputs = list(dict.fromkeys([*scripts, *source_inputs, *python_inputs, Path(sys.executable).resolve(),
                    args.archive.resolve(), output / "preparation.json", output / "embedded-modules.json",
                    output / "provider-source-join.json", stock_manifest, stock_driver, loader,
                    *[path for path in stage.rglob("*") if path.is_file()]]))
        (output / "owner-self-test").mkdir()
        common.memory_self_test()
        common.owner_self_test(output / "owner-self-test", environment, inputs)

        probe = output / "discover-cmake-root.cmake"
        with probe.open("x") as stream:
            stream.write(CMAKE_ROOT_PROBE)
        inputs += [probe]
        tool_folder = output / "toolchain"
        owner = common.run_owned([sys.executable, "-B", str(HERE), "--toolchain-child"], environment, tool_folder, inputs, seconds=60)
        complete_owner(owner)
        tools = json.loads((tool_folder / "stdout.log").read_text())
        selected_tools, tool_links = tool_inputs(tools)
        write(output / "toolchain-input-layout.json", {"files": [str(path) for path in selected_tools], "symlinks": tool_links,
              "scope": "Complete regular files under the selected SDK/compiler library/include/module roots plus selected executables. System dyld-cache images are observed, not a hermetic OS attestation."})
        tool_receipts = [output / "toolchain-input-layout.json", tool_folder / "stdout.log", tool_folder / "owner.json",
                         probe]

        def inventory(label, env, driver, expected_info):
            folder = output / label
            owner = common.run_owned([str(app), "--json", "info", "system"], env, folder, inputs + [driver], seconds=60)
            complete_owner(owner)
            loaded = loaded_join(folder, driver, loader)
            document = json.loads((folder / "stdout.log").read_text())
            require(document["platform"] == "macOS" and document["wsl2"] is False, "inventory is not native macOS")
            _, selected = common.selected_inventory(document, True)
            require(selected["kind"] in ("integrated", "discrete") and selected["driver_name"] == "MoltenVK"
                    and selected["driver_info"] == expected_info, "actual selected provider label differs")
            return selected, loaded

        stock_selected, stock_loaded = inventory("stock-inventory", environment, stock_driver, "1.4.1")
        require(stock_selected["api_version"] == stock_api, "stock reported Vulkan API differs from original1.3.334")
        environment["SIRIUS_VULKAN_DEVICE"] = str(stock_selected["index"])

        def uuid_query(label, env, driver, selected, expected_info):
            folder = output / label
            command = [sys.executable, "-B", str(HERE.with_name("macos_vulkan_identity.py")), "--loader", str(loader),
                "--driver", str(driver), "--index", str(selected["index"]), "--expected-name", selected["name"],
                "--expected-driver-info", expected_info]
            owner = common.run_owned(command, env, folder, inputs + [driver], seconds=60)
            complete_owner(owner)
            data = json.loads((folder / "stdout.log").read_text())
            require(data["passed"] is True and data["qualification_claimed"] is False, "public Vulkan UUID query failed")
            require(all(data["device"][name] == selected[name] for name in ("name", "vendor_id", "device_id", "driver_id", "api_version", "driver_name", "driver_info")),
                    "public Vulkan properties do not join selected app inventory")
            return data

        stock_uuid = uuid_query("stock-uuid", environment, stock_driver, stock_selected, "1.4.1")
        cpm = output / "cpm-cache/cpm"
        cpm.mkdir()
        shutil.copy2(payload / "cpm-cache/cpm/CPM_0.40.8.cmake", cpm / "CPM_0.40.8.cmake")
        context_path = output / "build-context.json"
        build = REPO / "bin/macos/macos-moltenvk-integration"
        build_environment = dict(environment)
        build_environment.pop("DYLD_PRINT_LIBRARIES", None)
        write(context_path, {"payload": str(payload), "output": str(output), "build": str(build),
                            "toolchain": tools, "content_id": CONTENT_ID, "content_id_is_git_provenance": False,
                            "build_environment_difference": {"removed": ["DYLD_PRINT_LIBRARIES"],
                                "scope": "Build tools retain full input seals and the64MiB output cap; runtime provider inventories/queries/controls retain loaded-image logging"}})
        build_inputs = list(dict.fromkeys([*inputs, *selected_tools, *tool_receipts,
                                          context_path, cpm / "CPM_0.40.8.cmake"]))
        build_folder = output / "provider-build"
        owner = common.run_owned([sys.executable, "-B", str(HERE), "--build-child", str(context_path), identity(context_path)["sha256"]],
                                 build_environment, build_folder, build_inputs, seconds=900)
        complete_owner(owner)
        verify_payload(payload)
        after_tools, after_links = tool_inputs(tools)
        require(after_tools == selected_tools and after_links == tool_links, "tool/SDK input layout changed")
        require(identity(stock_driver) == STOCK_DRIVER and identity(loader) == STOCK_LOADER, "stock SDK changed")
        diagnostic_driver = (build / "MoltenVK/libMoltenVK.1.4.1.dylib").resolve(strict=True)
        require(diagnostic_driver.is_relative_to(build), "built provider escapes owned build root")
        header = diagnostic_driver.read_bytes()[:32]
        require(len(header) == 32 and struct.unpack("<8I", header)[:2] == (0xFEEDFACF, 0x0100000C)
                and struct.unpack("<8I", header)[3] == 6, "expected thin arm64 Mach-O dylib")
        generated = (build / "MoltenVK/mvkGitRevDerived.h").read_bytes()
        expected_header = (payload / "sources/moltenvk-patched/Templates/cmake/mvkGitRevDerived.h.in").read_bytes().replace(b"@MVK_GIT_REV@", CONTENT_ID.encode())
        require(generated == expected_header,
                "generated legacy revision field does not contain the declared diagnostic discriminator")
        icd = output / "MoltenVK-diagnostic_icd.json"
        write(icd, {"file_format_version": "1.0.0", "ICD": {"library_path": str(diagnostic_driver),
                    "api_version": "1.4.0", "is_portability_driver": True}})
        environment["VK_DRIVER_FILES"] = str(icd)
        environment["DYLD_LIBRARY_PATH"] = str(diagnostic_driver.parent) + os.pathsep + environment["DYLD_LIBRARY_PATH"]
        # Runtime consumes the completed build chain, not the SDK/compiler
        # input trees. Their complete byte seals remain in this one build owner.
        inputs += [context_path, *tool_receipts, build_folder / "owner.json",
                   build_folder / "inputs-before.json", build_folder / "inputs-after.json",
                   cpm / "CPM_0.40.8.cmake", diagnostic_driver, icd, build / "MoltenVK/mvkGitRevDerived.h"]
        selected, diagnostic_loaded = inventory("diagnostic-inventory", environment, diagnostic_driver, DIAGNOSTIC_INFO)
        require(selected["api_version"] == diagnostic_api, "rebuilt reported API does not match source-derived1.3.334")
        require({k: v for k, v in selected.items() if k not in ("driver_info", "api_version")}
                == {k: v for k, v in stock_selected.items() if k not in ("driver_info", "api_version")},
                "rebuilt provider changed selected device/capability inventory")
        diagnostic_uuid = uuid_query("diagnostic-uuid", environment, diagnostic_driver, selected, DIAGNOSTIC_INFO)
        require(diagnostic_uuid["device"]["pipeline_cache_uuid"] != stock_uuid["device"]["pipeline_cache_uuid"]
                and diagnostic_uuid["device"]["pipeline_cache_uuid"].startswith(CONTENT_ID[:8]), "diagnostic cache UUID aliases stock or lacks declared prefix")
        write(output / "provider-identities.json", {"stock": stock_uuid, "diagnostic": diagnostic_uuid,
              "app_loaded_images": {"stock": stock_loaded, "diagnostic": diagnostic_loaded},
              "expected_api_versions": result["expected_api_versions"],
              "provider_inputs": common.provider_inputs(app, environment, selected, True),
              "qualification_claimed": False, "stock_provider_qualified": False})
        inputs += [output / "provider-identities.json"]
        result["provider_identity_accepted"] = True
        for index, case in enumerate(common.CASES):
            folder = output / (str(index + 1) + "-" + case)
            local = dict(environment, GTEST_OUTPUT="xml:" + str(folder / "gtest.xml"), MVK_CONFIG_SHADER_DUMP_DIR=str(folder / "shaders"))
            full_name = "RetainedComputeTest." + case
            owner = common.run_owned([str(backend), "--gtest_filter=" + full_name], local, folder, inputs, seconds=300)
            complete_owner(owner, require_zero=False)
            require(common.STOP is None, "diagnostic interrupted")
            observed = common.case_result(folder / "gtest.xml", full_name)
            row = {"case": full_name, "owner_complete": True, "test": observed,
                   "exit": owner["child_exit"], "shader_capture_bound": False}
            result["cases"].append(row)
            if not observed["passed"] or owner["child_exit"] != 0:
                result["status"] = "original_control_failed"
                raise RuntimeError("original numerical case failed or skipped; no following case launched")
            loaded = loaded_join(folder, diagnostic_driver, loader)
            shaders = shader_join(folder, modules, EXPECTED_MODULES[index])
            write(folder / "provider-shader-joins.json", {"loaded": loaded, "shaders": shaders})
            row["shader_capture_bound"] = True
            write(output / "integration.json", result)
        result["status"] = "four_original_controls_passed_custom_provider"
        result["all_four_controls_passed"] = True
        result["numerics_executed"] = True
        return 0
    except Exception as error:
        result["error"] = f"{type(error).__name__}: {error}"
        return 1
    finally:
        if common.STOP is not None:
            result["status"] = "incomplete"
            result["error"] = "integration interrupted by signal " + str(common.STOP)
        if payload is not None:
            try:
                verify_payload(payload)
                result["provider_source_inputs_unchanged"] = True
            except Exception as error:
                result["status"] = "incomplete"
                result["provider_source_inputs_unchanged"] = False
                result["provider_source_verification_error"] = str(error)
                result.setdefault("error", str(error))
        write(output / "integration.json", result)
        print(json.dumps(result), flush=True)
        if common.STOP is not None or result.get("provider_source_inputs_unchanged") is False:
            return 1


if __name__ == "__main__":
    if sys.argv[1:] == ["--toolchain-child"]:
        toolchain_child()
    elif len(sys.argv) == 4 and sys.argv[1] == "--build-child":
        build_child(Path(sys.argv[2]), sys.argv[3])
    else:
        raise SystemExit(main())
