#!/usr/bin/env python3
"""Provision the pinned custom MoltenVK for a clean current-source native run.

Reuses the diagnostic's source/tool/build helpers and mature Darwin owner, never
its exported executable, diagnostic source predicate, or numerical-control loop.
Provision and property receipts do not admit an external qualification domain.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib
import json
import os
from pathlib import Path
import platform
import re
import shutil
import struct
import subprocess
import sys
import zipfile

ROOT = Path(__file__).resolve().parent.parent
OUTPUT = ROOT / "out/macos-moltenvk-integration/observation"
BUILD = ROOT / "bin/macos/macos-moltenvk-integration"
RECEIPT = OUTPUT / "native-provider.json"
_DARWIN = None


def require(condition, message):
    if not condition:
        raise ValueError(message)


def pin(path):
    path = Path(path).resolve(strict=True)
    data = path.read_bytes()
    return {"path": str(path), "bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}


def write(path, value):
    path.write_text(json.dumps(value, sort_keys=True, indent=2) + "\n", encoding="utf-8")


def strict_source_snapshot(expected_revision, query=None):
    """Record only Git identity/cleanliness, never the inherited environment."""
    require(re.fullmatch(r"[0-9a-f]{40}", expected_revision or "") is not None,
            "full current-source revision required")
    if query is None:
        def query(arguments):
            return subprocess.check_output(["git", *arguments], cwd=ROOT, text=True).strip()
    revision = query(["rev-parse", "HEAD"])
    status = query(["status", "--porcelain", "--untracked-files=normal"])
    require(revision == expected_revision and status == "", "current native source changed or is dirty")
    return {"runner_revision": revision, "candidate_revision": revision, "status": status}


def diagnostic():
    # resource and the Darwin owner are unavailable on Windows. Import only in
    # the actual macOS provisioning/binding path; metadata controls stay portable.
    global _DARWIN
    require(sys.platform == "darwin" and platform.machine() == "arm64", "native arm64 macOS required")
    if _DARWIN is None:
        _DARWIN = importlib.import_module("macos_moltenvk_integration")
    return _DARWIN


def source_inputs():
    names = subprocess.check_output(["git", "ls-files", "-z"], cwd=ROOT).decode().split("\0")
    paths = [ROOT / name for name in names if name]
    require(paths and all(path.is_file() and path.resolve().is_relative_to(ROOT) for path in paths),
            "tracked native source contains missing/escaping inputs")
    python = [Path(module.__file__).resolve() for module in list(sys.modules.values())
              if getattr(module, "__file__", None) and Path(module.__file__).is_file()]
    return list(dict.fromkeys([*paths, *python, Path(sys.executable).resolve()]))


def owned(d, revision, argv, environment, folder, inputs, seconds):
    owner = d.common.run_owned(argv, environment, folder, inputs, seconds=seconds,
                              source_snapshot=lambda: strict_source_snapshot(revision))
    d.complete_owner(owner)
    return owner


def verify_pins(pins):
    require(isinstance(pins, list) and pins, "native provider has no sealed inputs")
    for record in pins:
        require(pin(record["path"]) == record, "native provider input changed: " + record["path"])


def validate_metadata(receipt, revision, expected):
    require(receipt.get("schema_version") == 1 and receipt.get("status") == "built-pinned-custom-provider"
            and receipt.get("source") == {"runner_revision": revision, "candidate_revision": revision, "status": ""}
            and receipt.get("qualification_claimed") is False
            and receipt.get("stock_sdk_modified") is False, "native provider receipt is incomplete/stale")
    require(all(receipt.get(key) == value for key, value in expected.items()),
            "native provider source/content metadata differs")
    driver = (BUILD / "MoltenVK/libMoltenVK.1.4.1.dylib").resolve()
    environment = receipt["runtime_environment"]
    require(receipt["driver"]["path"] == str(driver)
            and set(environment) == {"VK_DRIVER_FILES", "DYLD_LIBRARY_PATH", "PYTHONDONTWRITEBYTECODE"}
            and environment["VK_DRIVER_FILES"] == str(OUTPUT / "MoltenVK-diagnostic_icd.json")
            and environment["DYLD_LIBRARY_PATH"] == os.pathsep.join([str(driver.parent), *receipt["sdk_library_paths"]])
            and receipt["sdk_library_paths"]
            and environment["PYTHONDONTWRITEBYTECODE"] == "1",
            "native provider selection differs from its private driver and unchanged SDK loader route")


def verify_receipt(receipt, revision):
    d = diagnostic()
    validate_metadata(receipt, revision, {"evidence_commit": d.EVIDENCE_COMMIT, "payload_tree": d.PAYLOAD_TREE,
                                        "content_id": d.CONTENT_ID, "driver_info": d.DIAGNOSTIC_INFO})
    verify_pins(receipt["sealed_inputs"])
    d.verify_payload(OUTPUT / "provider-source")
    require(pin(receipt["stock_driver"]["path"]) == receipt["stock_driver"]
            and pin(receipt["loader"]["path"]) == receipt["loader"], "stock SDK changed")
    strict_source_snapshot(revision)
    return receipt


def load(revision):
    return verify_receipt(json.loads(RECEIPT.read_text()), revision)


def activate(receipt):
    """Select the private ICD and the unchanged SDK loader for all native gates."""
    for key in list(os.environ):
        if key.startswith(("DYLD_", "MVK_CONFIG_", "VK_LOADER_")) or key in {
            "VK_ICD_FILENAMES", "VK_DRIVER_FILES", "VK_ADD_DRIVER_FILES", "VK_INSTANCE_LAYERS",
            "VK_LAYER_PATH", "VK_ADD_LAYER_PATH"}:
            del os.environ[key]
    os.environ.update(receipt["runtime_environment"])


def provision(revision):
    d = diagnostic()
    strict_source_snapshot(revision)
    require(not OUTPUT.exists() and not BUILD.exists(), "native provider roots must be fresh")
    OUTPUT.mkdir(parents=True)
    result = {"schema_version": 1, "status": "incomplete", "source": strict_source_snapshot(revision),
              "qualification_claimed": False, "stock_sdk_modified": False,
              "scope": "Custom annotation-only SPIRV-Cross provider with its declared diagnostic identity; current-source qualification remains separate"}
    try:
        payload, payload_inputs = d.stage_payload(OUTPUT)
        sdk = Path(os.environ["VULKAN_SDK"]).resolve(strict=True)
        require("1.4.357.0" in sdk.parts, "pinned Vulkan SDK version required")
        stock_manifest = sdk / "share/vulkan/icd.d/MoltenVK_icd.json"
        stock_driver = (stock_manifest.parent / json.loads(stock_manifest.read_text())["ICD"]["library_path"]).resolve(strict=True)
        loader = (sdk / "lib/libvulkan.1.dylib").resolve(strict=True)
        require(d.identity(stock_driver) == d.STOCK_DRIVER and d.identity(loader) == d.STOCK_LOADER,
                "stock SDK binary pins differ")
        environment = d.clean_environment(OUTPUT, sdk)
        inputs = list(dict.fromkeys([*source_inputs(), *payload_inputs, OUTPUT / "provider-source-join.json",
                                     stock_manifest, stock_driver, loader]))
        (OUTPUT / "owner-self-test").mkdir()
        d.common.memory_self_test()
        d.common.owner_self_test(OUTPUT / "owner-self-test", environment, inputs,
                                source_snapshot=lambda: strict_source_snapshot(revision))
        probe = OUTPUT / "discover-cmake-root.cmake"
        probe.write_text(d.CMAKE_ROOT_PROBE)
        inputs.append(probe)
        tool_folder = OUTPUT / "toolchain"
        owned(d, revision, [sys.executable, "-B", str(d.HERE), "--toolchain-child"], environment,
              tool_folder, inputs, 60)
        tools = json.loads((tool_folder / "stdout.log").read_text())
        selected_tools, links = d.tool_inputs(tools)
        layout = OUTPUT / "toolchain-input-layout.json"
        write(layout, {"files": [str(path) for path in selected_tools], "symlinks": links,
                      "scope": "Selected complete SDK/compiler/CMake trees; system dyld cache is not a hermetic OS attestation"})
        tool_receipts = [layout, tool_folder / "stdout.log", tool_folder / "owner.json", probe]
        cpm = OUTPUT / "cpm-cache/cpm"
        cpm.mkdir()
        shutil.copy2(payload / "cpm-cache/cpm/CPM_0.40.8.cmake", cpm / "CPM_0.40.8.cmake")
        context = OUTPUT / "build-context.json"
        write(context, {"payload": str(payload), "output": str(OUTPUT), "build": str(BUILD),
                        "toolchain": tools, "content_id": d.CONTENT_ID, "content_id_is_git_provenance": False,
                        "build_environment_difference": {"removed": ["DYLD_PRINT_LIBRARIES"],
                            "scope": "Build tools retain complete input seals; runtime property observations retain loaded-image logging"}})
        build_inputs = list(dict.fromkeys([*inputs, *selected_tools, *tool_receipts, context, cpm / "CPM_0.40.8.cmake"]))
        build_environment = dict(environment)
        build_environment.pop("DYLD_PRINT_LIBRARIES", None)
        build_folder = OUTPUT / "provider-build"
        owned(d, revision, [sys.executable, "-B", str(d.HERE), "--build-child", str(context),
                           d.identity(context)["sha256"]], build_environment, build_folder, build_inputs, 900)
        d.verify_payload(payload)
        after_files, after_links = d.tool_inputs(tools)
        require(after_files == selected_tools and after_links == links, "selected tool layout changed")
        require(d.identity(stock_driver) == d.STOCK_DRIVER and d.identity(loader) == d.STOCK_LOADER, "stock SDK changed")
        driver = (BUILD / "MoltenVK/libMoltenVK.1.4.1.dylib").resolve(strict=True)
        require(driver.is_relative_to(BUILD), "native provider escapes the private build")
        header = driver.read_bytes()[:32]
        require(len(header) == 32 and struct.unpack("<8I", header)[:2] == (0xFEEDFACF, 0x0100000C)
                and struct.unpack("<8I", header)[3] == 6, "native provider is not a thin arm64 dylib")
        generated = BUILD / "MoltenVK/mvkGitRevDerived.h"
        expected = (payload / "sources/moltenvk-patched/Templates/cmake/mvkGitRevDerived.h.in").read_bytes().replace(
            b"@MVK_GIT_REV@", d.CONTENT_ID.encode()) + b"\n"
        require(generated.read_bytes() == expected, "generated provider identity differs")
        icd = OUTPUT / "MoltenVK-diagnostic_icd.json"
        write(icd, {"file_format_version": "1.0.0", "ICD": {"library_path": str(driver),
                    "api_version": "1.4.0", "is_portability_driver": True}})
        sdk_library_paths = environment["DYLD_LIBRARY_PATH"].split(os.pathsep)
        environment["VK_DRIVER_FILES"] = str(icd)
        environment["DYLD_LIBRARY_PATH"] = str(driver.parent) + os.pathsep + environment["DYLD_LIBRARY_PATH"]
        snapshots = OUTPUT / "provider-binaries"
        snapshots.mkdir()
        binary_snapshots = {}
        for name, original in (("custom-MoltenVK.dylib", driver), ("SDK-vulkan-loader.dylib", loader)):
            copied = snapshots / name
            shutil.copy2(original, copied)
            require(d.identity(copied) == d.identity(original), "provider/loader snapshot byte join failed")
            binary_snapshots[name] = pin(copied)
        sealed = list(dict.fromkeys([*inputs, *tool_receipts, context, cpm / "CPM_0.40.8.cmake", driver, icd,
            generated, *[Path(row["path"]) for row in binary_snapshots.values()], build_folder / "owner.json", build_folder / "inputs-before.json", build_folder / "inputs-after.json"]))
        result.update(status="built-pinned-custom-provider", evidence_commit=d.EVIDENCE_COMMIT,
            payload_tree=d.PAYLOAD_TREE, content_id=d.CONTENT_ID, content_id_is_git_provenance=False,
            driver_info=d.DIAGNOSTIC_INFO, binary_snapshots=binary_snapshots, driver=pin(driver), stock_driver=pin(stock_driver), loader=pin(loader),
            sdk_library_paths=sdk_library_paths, runtime_environment={key: environment[key] for key in ("VK_DRIVER_FILES", "DYLD_LIBRARY_PATH", "PYTHONDONTWRITEBYTECODE")},
            sealed_inputs=[pin(path) for path in sealed],
            tool_input_scope="Complete build inputs are preserved by the completed owner, not rehashed as runtime dependencies",
            source=strict_source_snapshot(revision))
        verify_receipt(result, revision)
        return result
    except Exception as error:
        result.update(status="incomplete", error=f"{type(error).__name__}: {error}")
        raise
    finally:
        write(RECEIPT, result)


def loaded_paths(raw, driver, loader):
    found = {}
    for name, expected in (("libMoltenVK", Path(driver)), ("libvulkan", Path(loader))):
        paths = {Path(value).resolve() for value in re.findall(r"(/[^\n]*?/" + name + r"[^/\s]*\.dylib)", raw)}
        require(paths == {expected.resolve()}, "actual native loaded image differs: " + name)
        found[name] = pin(expected)
    return found


def bind_candidate(receipt, revision, binary, backend, preflight, environment, output):
    d = diagnostic()
    verify_receipt(receipt, revision)
    require(binary.resolve() == ROOT / "bin/macos/src/sirius/app/sirius"
            and backend.resolve() == ROOT / "bin/macos/tests/backend/sirius_backend_tests",
            "native provider must bind current preset executables")
    require(pin(binary)["sha256"] == preflight["candidate"]["sha256"], "current candidate differs from preflight")
    modules, _ = d.common.macho_modules(backend)
    write(output / "native-provider-modules.json", {"backend": pin(backend), "modules": modules,
          "msl_capture_claimed": False,
          "scope": "35 named current-backend embedded modules only; no exported executable reuse or new MSL capture claim"})
    folder = OUTPUT / "current-candidate"
    folder.mkdir()
    inputs = list(dict.fromkeys([Path(record["path"]) for record in receipt["sealed_inputs"]] +
                               [RECEIPT, binary, backend, output / "native-provider-modules.json"]))
    local = dict(environment, DYLD_PRINT_LIBRARIES="1")
    selected = preflight["device"]
    driver, loader = Path(receipt["driver"]["path"]), Path(receipt["loader"]["path"])
    inventory_folder = folder / "inventory"
    owned(d, revision, [str(binary), "--json", "info", "system"], local, inventory_folder, inputs, 60)
    inventory = json.loads((inventory_folder / "stdout.log").read_text())
    _, actual = d.common.selected_inventory(inventory, True)
    require(actual == selected and actual["driver_info"] == d.DIAGNOSTIC_INFO
            and actual["driver_name"] == "MoltenVK", "current candidate provider/device differs from preflight")
    loaded = d.loaded_join(inventory_folder, driver, loader)
    def uuid(label, env, selected_driver, info):
        subfolder = folder / label
        command = [sys.executable, "-B", str(ROOT / "scripts/macos_vulkan_identity.py"),
                   "--loader", str(loader), "--driver", str(selected_driver), "--index", str(selected["index"]),
                   "--expected-name", selected["name"], "--expected-driver-info", info]
        owned(d, revision, command, env, subfolder, inputs, 60)
        data = json.loads((subfolder / "stdout.log").read_text())
        require(data["passed"] is True and data["qualification_claimed"] is False, "native UUID property observation failed")
        require(all(data["device"][key] == selected[key] for key in (
            "name", "vendor_id", "device_id", "driver_id", "api_version", "driver_name")),
            "native UUID device differs from current candidate")
        d.loaded_join(subfolder, selected_driver, loader)
        return data
    custom = uuid("custom-uuid", local, driver, d.DIAGNOSTIC_INFO)
    stock = dict(local)
    stock["VK_DRIVER_FILES"] = str(Path(os.environ["VULKAN_SDK"]).resolve() / "share/vulkan/icd.d/MoltenVK_icd.json")
    stock["DYLD_LIBRARY_PATH"] = os.pathsep.join(part for part in local["DYLD_LIBRARY_PATH"].split(os.pathsep)
                                               if Path(part).resolve() != driver.parent)
    original = uuid("stock-uuid", stock, Path(receipt["stock_driver"]["path"]), "1.4.1")
    require(custom["device"]["pipeline_cache_uuid"].startswith(d.CONTENT_ID[:8])
            and custom["device"]["pipeline_cache_uuid"] != original["device"]["pipeline_cache_uuid"],
            "custom cache UUID aliases stock or lacks its declared content prefix")
    result = {"source_revision": revision, "candidate": pin(binary), "backend": pin(backend),
              "provision_receipt": pin(RECEIPT), "modules": pin(output / "native-provider-modules.json"),
              "custom": custom, "stock": original, "loaded_images": loaded,
              "qualification_claimed": False, "scope": "Current source build and actual custom provider/property joins; all original qualification gates remain required"}
    write(output / "native-provider-bindings.json", result)
    return result


def verify_binding(receipt, binding, revision):
    verify_receipt(receipt, revision)
    verify_pins([binding[name] for name in ("candidate", "backend", "provision_receipt", "modules")])


def bind_frame(receipt, binding, revision, raw, output):
    verify_binding(receipt, binding, revision)
    write(output / "native-provider-frame-images.json", {
        "loaded_images": loaded_paths(raw, receipt["driver"]["path"], receipt["loader"]["path"]),
        "candidate": binding["candidate"], "qualification_claimed": False,
        "scope": "Actual loader/driver images in the unchanged original native Vulkan frame transcript"})


def bundle_proofs(receipt, binding, revision, output):
    """Hash auxiliary JSON as an explicit bundle, preserving core JSON roles."""
    verify_binding(receipt, binding, revision)
    copied = output / "native-provider-provision.json"
    shutil.copy2(RECEIPT, copied)
    require({key: pin(copied)[key] for key in ("bytes", "sha256")}
            == {key: binding["provision_receipt"][key] for key in ("bytes", "sha256")},
            "copied provision receipt differs from the consumed receipt")
    files = [output / name for name in ("native-provider-bindings.json", "native-provider-modules.json",
                                        "native-provider-frame-images.json", "native-provider-provision.json")]
    archive = output / "native-provider-proof.zip"
    with zipfile.ZipFile(archive, "x", compression=zipfile.ZIP_STORED) as bundle:
        for path in files:
            member = zipfile.ZipInfo(path.name, date_time=(1980, 1, 1, 0, 0, 0))
            member.create_system, member.external_attr = 3, 0o100644 << 16
            bundle.writestr(member, path.read_bytes())
    with zipfile.ZipFile(archive, "r") as bundle:
        require(bundle.namelist() == [path.name for path in files]
                and len(bundle.infolist()) == 4
                and all(bundle.read(path.name) == path.read_bytes() for path in files),
                "auxiliary provider proof archive differs from its four exact sidecars")
    return archive


def self_test():
    imported = _DARWIN
    revision = "a" * 40
    def query(arguments):
        return revision if arguments[0] == "rev-parse" else ""
    require(strict_source_snapshot(revision, query)["candidate_revision"] == revision, "valid clean source rejected")
    for expected, reader in (("stale", query), ("b" * 40, query),
                             (revision, lambda args: revision if args[0] == "rev-parse" else " M scripts/source.py")):
        try:
            strict_source_snapshot(expected, reader)
        except ValueError:
            pass
        else:
            raise ValueError("stale/dirty current source was accepted")
    # Source-only controls run through the existing native producer self-test;
    # none of these checks import the Darwin owner or execute a tool/provider.
    driver = str((BUILD / "MoltenVK/libMoltenVK.1.4.1.dylib").resolve())
    expected = {"evidence_commit": "b" * 40, "payload_tree": "c" * 40,
                "content_id": "d" * 64, "driver_info": "explicit custom identity"}
    valid = dict(expected, schema_version=1, status="built-pinned-custom-provider",
                 source=strict_source_snapshot(revision, query), qualification_claimed=False,
                 stock_sdk_modified=False, driver={"path": driver}, sdk_library_paths=["/sdk/lib"],
                 runtime_environment={"VK_DRIVER_FILES": str(OUTPUT / "MoltenVK-diagnostic_icd.json"),
                     "DYLD_LIBRARY_PATH": str(Path(driver).parent) + os.pathsep + "/sdk/lib",
                     "PYTHONDONTWRITEBYTECODE": "1"})
    validate_metadata(valid, revision, expected)
    mutations = []
    for key, value in (("status", "incomplete"), ("source", {"runner_revision": "e" * 40}),
                       ("qualification_claimed", True), ("stock_sdk_modified", True),
                       ("content_id", "wrong"), ("driver", {"path": "/sdk/lib/libMoltenVK.dylib"})):
        mutations.append(dict(valid, **{key: value}))
    for key, value in (("VK_DRIVER_FILES", "/sdk/stock.json"), ("DYLD_LIBRARY_PATH", "/sdk/lib"),
                       ("PYTHONDONTWRITEBYTECODE", "0"), ("DYLD_INSERT_LIBRARIES", "/injected.dylib")):
        mutations.append(dict(valid, runtime_environment=dict(valid["runtime_environment"], **{key: value})))
    for invalid in mutations:
        try:
            validate_metadata(invalid, revision, expected)
        except ValueError:
            pass
        else:
            raise ValueError("stale/promoting/stock/injected provider metadata accepted")
    require(_DARWIN is imported, "portable provider metadata controls imported native helpers")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--expected-revision")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return 0
    receipt = provision(args.expected_revision)
    # Only these nonsecret, explicit provider-selection values enter later job
    # steps; the stock SDK remains untouched and its loader is sealed above.
    if os.environ.get("GITHUB_ENV"):
        with Path(os.environ["GITHUB_ENV"]).open("a") as stream:
            for key, value in receipt["runtime_environment"].items():
                require("\n" not in value and "\r" not in value, "unsafe provider environment value")
                stream.write(key + "=" + value + "\n")
    print(json.dumps({"provider_receipt": pin(RECEIPT), "qualification_claimed": False}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
