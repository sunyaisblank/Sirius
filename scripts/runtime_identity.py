#!/usr/bin/env python3
"""Capture actual non-render inventory before numerics; never issue a gate receipt.

Provider files are hashed loader inputs/candidates, not proof that every listed
ICD was loaded or executed. Device and float controls come from the real app.
"""
import argparse
import base64
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import tempfile

CONTROLS = ("supports_fp64", "preserves_fp32_denormals",
            "rounds_fp32_to_nearest", "rounds_fp64_to_nearest")
ENVIRONMENT = ("VK_DRIVER_FILES", "VK_ICD_FILENAMES", "VK_ADD_DRIVER_FILES",
               "VK_LOADER_DRIVERS_SELECT", "VK_LOADER_DRIVERS_DISABLE",
               "VULKAN_SDK", "LD_LIBRARY_PATH", "DYLD_LIBRARY_PATH", "PATH", "SystemRoot",
               "SIRIUS_VULKAN_DEVICE", "SIRIUS_PRECISION", "LP_PERF",
               "LP_NATIVE_VECTOR_WIDTH", "GALLIVM_PERF", "MESA_SHADER_CACHE_DISABLE")


def require(condition, message):
    if not condition:
        raise ValueError(message)


def identity(path):
    path = Path(path).resolve(strict=True)
    require(path.is_file(), f"not a regular identity input: {path}")
    before = path.stat()
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    after = path.stat()
    require((before.st_size, before.st_mtime_ns) == (after.st_size, after.st_mtime_ns),
            f"identity input changed while hashing: {path}")
    return {"path": str(path), "bytes": after.st_size, "sha256": digest.hexdigest()}


def selected_inventory(inventory, required):
    vk = inventory.get("backends", {}).get("vulkan", {})
    require(type(vk.get("compiled")) is bool, "inventory lacks actual compiled-backend state")
    if not vk["compiled"]:
        require(not required, "required Vulkan identity cannot use a CPU-only product")
        return "cpu_only", None
    devices = vk.get("devices", [])
    require(isinstance(devices, list), "Vulkan inventory devices must be a list")
    require(all(isinstance(device, dict) and type(device.get("index")) is int
                and device["index"] >= 0 for device in devices), "invalid actual device indices")
    if not devices:
        require(not required, "required Vulkan identity has no actual visible device")
        return "vulkan_without_visible_device", None
    index = vk.get("selected_device_index")
    require(type(index) is int and index >= 0, "Vulkan inventory lacks an actual selected index")
    matches = [device for device in devices if device.get("index") == index]
    require(len(matches) == 1, "selected Vulkan device is missing or ambiguous")
    device = matches[0]
    for field in CONTROLS:
        require(type(device.get(field)) is bool, f"missing actual device control: {field}")
    for field in ("name", "kind"):
        require(isinstance(device.get(field), str) and device[field], f"missing device identity: {field}")
    for field in ("driver_name", "driver_info"):
        require(isinstance(device.get(field), str), f"missing reported device identity: {field}")
    for field in ("vendor_id", "device_id", "driver_id", "api_version"):
        require(type(device.get(field)) is int and device[field] >= 0, f"invalid device identity: {field}")
    return "vulkan_selected", device


def search_directories(executable, environment):
    directories = [Path(executable).parent]
    for key in ("PATH", "LD_LIBRARY_PATH", "DYLD_LIBRARY_PATH"):
        directories.extend(Path(value) for value in environment.get(key, "").split(os.pathsep) if value)
    if environment.get("VULKAN_SDK"):
        sdk = Path(environment["VULKAN_SDK"])
        directories.extend(sdk / child for child in ("bin", "lib", "Bin", "Lib", "runtime/x64"))
    if sys.platform == "win32":
        directories.append(Path(environment.get("SystemRoot") or "C:/Windows") / "System32")
    else:
        directories.extend(Path(value) for value in ("/lib", "/usr/lib", "/usr/local/lib"))
        for base in (Path("/lib"), Path("/usr/lib")):
            directories.extend(sorted(base.glob("*-linux-gnu")))
    return list(dict.fromkeys(path.resolve() for path in directories))


def library_candidates(name, manifest, directories):
    path = Path(name)
    if path.is_absolute():
        candidates = [path]
    elif "/" in name or "\\" in name:
        candidates = [manifest.parent / path]
    else:
        # A bare soname uses the system loader search. Preserve resolvable
        # candidates; never claim this list establishes which image was loaded.
        candidates = [directory / name for directory in directories]
    return [identity(path) for path in dict.fromkeys(path.resolve() for path in candidates)
            if path.is_file()]


def provider_inputs(executable, environment, selected, required):
    effective = "VK_DRIVER_FILES" if environment.get("VK_DRIVER_FILES") else "VK_ICD_FILENAMES"
    explicit = environment.get(effective, "")
    if explicit:
        paths = [Path(value) for value in explicit.split(os.pathsep) if value]
        selector = {"name": effective, "value": explicit, "kind": "explicit_manifest_list"}
    else:
        # Default catalogues are input candidates, not a reconstructed loader
        # selection. Full hosted jobs use an explicit provider manifest.
        roots = [Path("/etc/vulkan/icd.d"), Path("/usr/share/vulkan/icd.d"),
                 Path("/usr/local/share/vulkan/icd.d")]
        sdk = environment.get("VULKAN_SDK")
        if sdk:
            roots.append(Path(sdk) / "share/vulkan/icd.d")
        paths = sorted({path for root in roots for path in root.glob("*.json")}) if selected is not None else []
        selector = {"name": None, "value": None, "kind": "default_catalogue_candidates"}
    if selected is None:
        return {"scope": "hashed input candidates; not actual loaded-module or per-ICD execution proof",
                "selector": selector, "manifests": [], "loader_candidates": [],
                "unresolved_optional_manifests": [], "no_selected_device": True}
    directories = search_directories(executable, environment)
    manifests = []
    unresolved = []
    for path in paths:
        if not path.is_file() and not required:
            unresolved.append(str(path))
            continue
        record = identity(path)
        raw = path.read_bytes()
        require(len(raw) == record["bytes"] and hashlib.sha256(raw).hexdigest() == record["sha256"],
                f"provider manifest changed while reading: {path}")
        document = json.loads(raw.decode("utf-8-sig"))
        library = document.get("ICD", {}).get("library_path")
        require(isinstance(library, str) and library, f"ICD lacks library_path: {path}")
        candidates = library_candidates(library, path, directories)
        if required:
            require(candidates, f"required provider library input is not resolvable: {library}")
        manifests.append({"selector_path": str(path), "manifest": record,
                          "contents_base64": base64.b64encode(raw).decode("ascii"), "library_path": library,
                          "driver_candidates": candidates})
    loader_names = ("vulkan-1.dll",) if sys.platform == "win32" else (
        ("libvulkan.1.dylib", "libvulkan.dylib") if sys.platform == "darwin" else ("libvulkan.so.1",))
    loaders = [record for name in loader_names
               for record in library_candidates(name, Path(executable), directories)]
    loaders = list({record["path"]: record for record in loaders}.values())
    if required and selected is not None:
        require(manifests and loaders, "required Vulkan provider/loader input identity is unresolved")
    return {"scope": "hashed input candidates; not actual loaded-module or per-ICD execution proof",
            "selector": selector, "manifests": manifests, "loader_candidates": loaders,
            "unresolved_optional_manifests": unresolved}


def validate_record(record, revision, executable_record, inventory=None, required=True, clean=None):
    require(type(record.get("schema_version")) is int and record["schema_version"] == 1
            and record.get("kind") == "sirius-runtime-identity"
            and record.get("status") == "captured", "runtime identity is not a captured diagnostic")
    require(record.get("qualification_claimed") is False and record.get("numerics_executed") is False,
            "identity query must not claim numerical or release qualification")
    require(record.get("phase") in {"before_mandatory_numerics", "before_native_runtime_estate"},
            "identity was not captured before its numerical estate")
    require(record.get("source", {}).get("revision") == revision, "identity source differs")
    require(type(record.get("source", {}).get("clean")) is bool, "identity source clean state is missing")
    if clean is not None:
        require(record["source"]["clean"] is clean, "identity source clean state differs")
    ci = record.get("ci", {})
    require(type(ci.get("hosted")) is bool, "identity execution context is missing")
    ci_fields = ("repository", "run_id", "run_attempt", "job", "runner", "os", "arch")
    require(all(value in ci and (ci[value] is None or isinstance(ci[value], str)) for value in ci_fields),
            "identity run/job metadata is malformed")
    if ci["hosted"]:
        require(all(ci[value] for value in ci_fields), "hosted identity lacks actual run/job metadata")
        require(all(ci[value].isdigit() and int(ci[value]) > 0 for value in ("run_id", "run_attempt")),
                "hosted identity run/attempt is invalid")
    executable = record.get("executable", {})
    require(all(executable.get(key) == executable_record[key] for key in ("bytes", "sha256")),
            "identity executable differs")
    observed = record.get("inventory")
    require(isinstance(observed, dict), "identity inventory is missing")
    if inventory is not None:
        require(observed == inventory, "identity differs from the bound runtime inventory")
    classification, selected = selected_inventory(observed, required)
    require(record.get("classification") == classification and record.get("selected_device") == selected,
            "identity selected device/classification differs from actual inventory")
    query = record.get("query", {})
    require(query.get("argv") == ["--json", "info", "system"] and type(query.get("exit_code")) is int
            and query["exit_code"] == 0,
            "identity is not the successful non-render query")
    for stream in ("stdout", "stderr"):
        payload = query.get(stream + "_base64")
        require(isinstance(payload, str), f"identity raw query {stream} payload is missing")
        raw = base64.b64decode(payload, validate=True)
        require(hashlib.sha256(raw).hexdigest() == query.get(stream + "_sha256"),
                f"identity query {stream} digest differs")
        if stream == "stdout":
            require(json.loads(raw.decode("utf-8-sig")) == observed,
                    "raw actual query differs from the recorded inventory")
    inputs = record.get("provider_inputs", {})
    require(inputs.get("scope") == "hashed input candidates; not actual loaded-module or per-ICD execution proof",
            "provider identity scope was overstated")
    environment = record.get("environment", {})
    require(isinstance(environment, dict) and all(isinstance(value, str) for value in environment.values()),
            "identity environment is missing or malformed")
    effective = "VK_DRIVER_FILES" if environment.get("VK_DRIVER_FILES") else "VK_ICD_FILENAMES"
    selector = inputs.get("selector", {})
    explicit = environment.get(effective, "")
    if explicit:
        require(selector == {"name": effective, "value": explicit, "kind": "explicit_manifest_list"},
                "effective provider selector differs from the captured environment")
    else:
        require(selector == {"name": None, "value": None, "kind": "default_catalogue_candidates"},
                "default provider catalogue was misclassified as explicit selection")
    manifests = inputs.get("manifests", [])
    loaders = inputs.get("loader_candidates", [])
    require(isinstance(manifests, list) and isinstance(loaders, list), "provider inputs must be lists")
    require(all(isinstance(value, dict) for value in manifests + loaders), "provider input records must be objects")
    if required:
        require(manifests and loaders, "required provider input identity was omitted")
        require(all(value.get("driver_candidates") for value in manifests),
                "required provider driver inputs were omitted")
    if explicit and selected is not None:
        host_platform = record.get("host", {}).get("platform")
        require(host_platform in ("win32", "darwin", "linux"), "provider input platform is missing")
        separator = ";" if host_platform == "win32" else ":"
        paths = [value for value in explicit.split(separator) if value]
        require([value.get("selector_path") for value in manifests] == paths,
                "provider manifests differ from the effective explicit selector")
    artifacts = list(loaders)
    for value in manifests:
        require(isinstance(value, dict) and isinstance(value.get("library_path"), str) and value["library_path"],
                "provider library path is missing")
        raw = base64.b64decode(value.get("contents_base64", ""), validate=True)
        manifest = value.get("manifest", {})
        require(len(raw) == manifest.get("bytes") and hashlib.sha256(raw).hexdigest() == manifest.get("sha256"),
                "provider manifest bytes/digest differ")
        require(json.loads(raw.decode("utf-8-sig")).get("ICD", {}).get("library_path") == value["library_path"],
                "provider library path differs from the bound manifest")
        require(isinstance(value.get("driver_candidates"), list), "provider driver candidates must be a list")
        artifacts.append(value.get("manifest", {}))
        artifacts.extend(value.get("driver_candidates", []))
    for value in artifacts:
        require(isinstance(value, dict), "provider artifact record must be an object")
        require(isinstance(value.get("path"), str) and value["path"], "provider input path is missing")
        require(type(value.get("bytes")) is int and value["bytes"] > 0, "provider input byte count is invalid")
        digest = value.get("sha256")
        require(isinstance(digest, str) and len(digest) == 64 and
                all(character in "0123456789abcdef" for character in digest), "provider input digest is invalid")
    return record


def capture(source_root, revision, clean, executable, required, phase, inventory=None,
            environment=None, output=None):
    source_root = Path(source_root).resolve()
    actual_revision = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=source_root, text=True).strip()
    actual_clean = not subprocess.check_output(["git", "status", "--porcelain"], cwd=source_root)
    require(actual_revision == revision and actual_clean == clean, "identity source changed after configuration")
    executable_before = identity(executable)
    effective_environment = dict(os.environ if environment is None else environment)
    environment = {key: effective_environment.get(key, "") for key in ENVIRONMENT}
    context = {"schema_version": 1, "kind": "sirius-runtime-identity", "status": "capture_failed",
               "qualification_claimed": False, "numerics_executed": False, "phase": phase,
               "source": {"revision": revision, "clean": clean}, "executable": executable_before,
               "environment": environment, "ci": {"hosted": effective_environment.get("GITHUB_ACTIONS") == "true",
                   **{key: effective_environment.get(name) for key, name in (
                   ("repository", "GITHUB_REPOSITORY"), ("run_id", "GITHUB_RUN_ID"),
                   ("run_attempt", "GITHUB_RUN_ATTEMPT"), ("job", "GITHUB_JOB"),
                   ("runner", "RUNNER_NAME"), ("os", "RUNNER_OS"), ("arch", "RUNNER_ARCH"))}}}
    if effective_environment.get("GITHUB_ACTIONS") == "true":
        require(all(isinstance(value, str) and value for key, value in context["ci"].items() if key != "hosted"),
                "hosted identity lacks actual run/job/runner metadata")
    if output is not None:
        Path(output).write_text(json.dumps(context, indent=2) + "\n")
    process = subprocess.run([str(executable), "--json", "info", "system"],
                             capture_output=True, timeout=60, check=False, env=effective_environment)
    query = {"argv": ["--json", "info", "system"], "exit_code": process.returncode}
    for stream in ("stdout", "stderr"):
        raw = getattr(process, stream)
        query[stream + "_base64"] = base64.b64encode(raw).decode("ascii")
        query[stream + "_sha256"] = hashlib.sha256(raw).hexdigest()
    if output is not None:
        Path(output).write_text(json.dumps({**context, "query": query}, indent=2) + "\n")
    require(process.returncode == 0, f"actual inventory query failed: {process.returncode}: {process.stderr.decode(errors='replace')}")
    actual_inventory = json.loads(process.stdout.decode("utf-8-sig"))
    if inventory is not None:
        require(actual_inventory == inventory, "actual identity query differs from producer's bound inventory")
    inventory = actual_inventory
    classification, selected = selected_inventory(inventory, required)
    inputs = provider_inputs(executable, environment, selected, required)
    require(identity(executable) == executable_before, "identity executable changed during query")
    require(subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=source_root, text=True).strip() == revision
            and (not subprocess.check_output(["git", "status", "--porcelain"], cwd=source_root)) == clean,
            "identity source changed during query")
    record = {"schema_version": 1, "kind": "sirius-runtime-identity", "status": "captured",
              "qualification_claimed": False, "numerics_executed": False, "phase": phase,
              "captured_utc": datetime.now(timezone.utc).isoformat(),
              "source": {"revision": revision, "clean": clean}, "executable": executable_before,
              "ci": context["ci"],
              "host": {"platform": sys.platform, "machine": platform.machine()},
              "environment": environment, "query": query, "inventory": inventory,
              "classification": classification, "selected_device": selected, "provider_inputs": inputs}
    return validate_record(record, revision, executable_before, inventory, required, clean)


def self_test():
    device = {"index": 0, "name": "software witness", "driver_name": "provider",
              "driver_info": "version", "kind": "software", "vendor_id": 1,
              "device_id": 0, "driver_id": 13, "api_version": 1,
              **dict.fromkeys(CONTROLS, False)}
    inventory = {"backends": {"vulkan": {"compiled": True, "devices": [device],
                                           "selected_device_index": 0}}}
    require(selected_inventory(inventory, True)[1] == device, "portable false controls were rejected")
    empty_driver = json.loads(json.dumps(inventory))
    empty_driver["backends"]["vulkan"]["devices"][0]["driver_info"] = ""
    require(selected_inventory(empty_driver, True)[1]["driver_info"] == "", "reported empty driver_info was rejected")
    require(selected_inventory({"backends": {"vulkan": {"compiled": False}}}, False)[0] == "cpu_only",
            "legitimate CPU-only classification failed")
    rejected = 0
    for change in ("missing_control", "ambiguous", "invalid_index", "boolean_index", "cpu_required", "device_required"):
        value = json.loads(json.dumps(inventory))
        vk = value["backends"]["vulkan"]
        if change == "missing_control": del vk["devices"][0][CONTROLS[1]]
        if change == "ambiguous": vk["devices"].append(dict(device))
        if change == "invalid_index": vk["selected_device_index"] = -1
        if change == "boolean_index": vk["devices"][0]["index"] = False
        if change == "cpu_required": vk["compiled"] = False
        if change == "device_required": vk["devices"] = []
        try: selected_inventory(value, True)
        except ValueError: rejected += 1
        else: raise ValueError(f"identity accepted negative control: {change}")
    temporary_root = Path(__file__).resolve().parent.parent / "out"
    temporary_root.mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="sirius-identity-", dir=temporary_root) as temporary:
        directory = Path(temporary)
        library = directory / "driver.bin"
        library.write_bytes(b"synthetic provider input; never loaded")
        manifest = directory / "icd.json"
        manifest.write_text(json.dumps({"ICD": {"library_path": "./driver.bin"}}))
        candidates = library_candidates("./driver.bin", manifest, [])
        require(candidates == [identity(library)], "relative provider input did not bind bytes")
        require(not library_candidates("./absent.bin", manifest, []), "absent provider was fabricated")
        executable = identity(library)
        raw = json.dumps(inventory).encode()
        raw_stderr = b"\xff\xfe retained raw stderr"
        record = {"schema_version": 1, "kind": "sirius-runtime-identity", "status": "captured",
                  "qualification_claimed": False, "numerics_executed": False,
                  "phase": "before_native_runtime_estate", "source": {"revision": "synthetic", "clean": True},
                  "executable": executable, "inventory": inventory, "classification": "vulkan_selected",
                  "selected_device": device, "environment": {"VK_DRIVER_FILES": str(manifest)},
                  "host": {"platform": sys.platform},
                  "ci": {"hosted": False, **dict.fromkeys(("repository", "run_id", "run_attempt", "job", "runner", "os", "arch"))},
                  "query": {"argv": ["--json", "info", "system"], "exit_code": 0,
                            "stdout_base64": base64.b64encode(raw).decode(),
                            "stdout_sha256": hashlib.sha256(raw).hexdigest(),
                            "stderr_base64": base64.b64encode(raw_stderr).decode(),
                            "stderr_sha256": hashlib.sha256(raw_stderr).hexdigest()},
                  "provider_inputs": {"scope": "hashed input candidates; not actual loaded-module or per-ICD execution proof",
                       "selector": {"name": "VK_DRIVER_FILES", "value": str(manifest), "kind": "explicit_manifest_list"},
                       "manifests": [{"selector_path": str(manifest), "manifest": identity(manifest),
                                      "contents_base64": base64.b64encode(manifest.read_bytes()).decode(), "library_path": "./driver.bin",
                                      "driver_candidates": candidates}], "loader_candidates": candidates}}
        validate_record(record, "synthetic", executable, inventory, True, True)
        for change in ("raw_digest", "raw_inventory", "selector", "source_clean", "executable",
                       "driver_digest", "loader_record", "missing_driver", "missing_inputs", "phase",
                       "manifest_selector", "manifest_digest", "library_path", "hosted_metadata"):
            value = json.loads(json.dumps(record))
            if change == "raw_digest": value["query"]["stdout_sha256"] = "0" * 64
            if change == "raw_inventory":
                changed = b"{}"
                value["query"]["stdout_base64"] = base64.b64encode(changed).decode()
                value["query"]["stdout_sha256"] = hashlib.sha256(changed).hexdigest()
            if change == "selector": value["environment"]["VK_DRIVER_FILES"] = "different.json"
            if change == "source_clean": value["source"]["clean"] = False
            if change == "executable": value["executable"]["sha256"] = "0" * 64
            if change == "driver_digest": value["provider_inputs"]["manifests"][0]["driver_candidates"][0]["sha256"] = "invalid"
            if change == "loader_record": value["provider_inputs"]["loader_candidates"] = [{"anything": True}]
            if change == "missing_driver": value["provider_inputs"]["manifests"][0]["driver_candidates"] = []
            if change == "missing_inputs": value["provider_inputs"]["manifests"] = []
            if change == "phase": value["phase"] = "after_numerics"
            if change == "manifest_selector": value["provider_inputs"]["manifests"][0]["selector_path"] = "different.json"
            if change == "manifest_digest": value["provider_inputs"]["manifests"][0]["manifest"]["sha256"] = "0" * 64
            if change == "library_path": value["provider_inputs"]["manifests"][0]["library_path"] = "wrong.dll"
            if change == "hosted_metadata": value["ci"]["hosted"] = True
            try: validate_record(value, "synthetic", executable, inventory, True, True)
            except ValueError: rejected += 1
            else: raise ValueError(f"identity accepted negative sidecar control: {change}")
        require(provider_inputs(library, {"VK_DRIVER_FILES": "unrelated-absent.json"}, None, False)["manifests"] == [],
                "CPU-only capture inspected an unrelated provider file")
    return {"passed": True, "rejected": rejected, "actual_query_executed": False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path)
    parser.add_argument("--source-revision")
    parser.add_argument("--source-tree-clean", choices=("true", "false"))
    parser.add_argument("--executable", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--require-vulkan", choices=("true", "false"))
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    try:
        if args.self_test:
            print(json.dumps(self_test()))
            return
        require(all((args.source_root, args.source_revision, args.source_tree_clean,
                     args.executable, args.output, args.require_vulkan)), "identity capture inputs are incomplete")
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.unlink(missing_ok=True)
        try:
            record = capture(args.source_root, args.source_revision, args.source_tree_clean == "true",
                             args.executable, args.require_vulkan == "true", "before_mandatory_numerics",
                             output=args.output)
        except Exception as error:
            failed = json.loads(args.output.read_text()) if args.output.is_file() else {"schema_version": 1, "kind": "sirius-runtime-identity",
                "status": "capture_failed", "qualification_claimed": False,
                "numerics_executed": False}
            failed["error"] = str(error)
            args.output.write_text(json.dumps(failed, indent=2) + "\n")
            raise
        args.output.write_text(json.dumps(record, indent=2) + "\n")
        print(f"Runtime identity captured ({record['classification']}): {args.output}")
    except (OSError, ValueError, KeyError, TypeError, subprocess.SubprocessError) as error:
        parser.exit(1, f"runtime identity capture: {error}\n")


if __name__ == "__main__":
    main()
