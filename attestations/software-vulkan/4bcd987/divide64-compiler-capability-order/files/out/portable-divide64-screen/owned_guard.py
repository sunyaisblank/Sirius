#!/usr/bin/env python3
"""One owned execution of the genuine linux-gcc RunMandatoryTests target.

Only this diagnostic's output directory is written directly. The unchanged
CMake target owns its real gate outputs. No partial or replacement receipt is
manufactured. Requires root's completed explicit artifact build beforehand.
"""
import argparse
import ctypes
import datetime
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import signal
import subprocess
import sys
import time

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[3]
BUILD = ROOT / "bin/linux-gcc"
ICD = Path("/usr/share/vulkan/icd.d/lvp_icd.json")
LVP = Path("/usr/lib/x86_64-linux-gnu/libvulkan_lvp.so")
LIMIT_SECONDS = 21600
RSS_KIB = 12 * 1024 * 1024
RESERVE_KIB = 2 * 1024 * 1024
SAMPLE_SECONDS = 2
STOP = None


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    digest_value = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest_value.update(chunk)
    return digest_value.hexdigest()


def identity(path):
    path = Path(path)
    return {"resolved_path": str(path.resolve(strict=True)), "bytes": path.stat().st_size,
            "sha256": digest(path)}


def write(path, value):
    temporary = path.with_suffix(".tmp")
    with temporary.open("w") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    temporary.replace(path)


def query(command, env=None):
    return subprocess.check_output(command, cwd=ROOT, env=env, timeout=60)


def source_identity():
    return {"revision": query(["git", "rev-parse", "HEAD"]).decode().strip(),
            "status": query(["git", "status", "--porcelain", "--untracked-files=normal"]).decode()}


def proc(pid):
    stat = Path(f"/proc/{pid}/stat").read_text()
    fields = stat[stat.rfind(")") + 2:].split()
    status = dict(line.split(":", 1) for line in Path(f"/proc/{pid}/status").read_text().splitlines()
                  if ":" in line)
    return {"pid": pid, "state": fields[0], "ppid": int(fields[1]),
            "pgrp": int(fields[2]), "session": int(fields[3]), "start_ticks": fields[19],
            "rss_kib": int(status.get("VmRSS", "0 kB").split()[0])}


def group_members(pgid):
    members = []
    for entry in Path("/proc").glob("[0-9]*"):
        try:
            item = proc(int(entry.name))
            if item["pgrp"] == pgid and item["session"] == pgid:
                members.append(item)
        except (FileNotFoundError, ProcessLookupError, PermissionError):
            pass
    return members


def available_kib():
    values = dict(line.split(":", 1) for line in Path("/proc/meminfo").read_text().splitlines())
    return int(values["MemAvailable"].split()[0])


def signal_stop(number, _frame):
    global STOP
    STOP = f"supervisor_signal_{number}"


def terminate_group(child, start_ticks, reaped):
    """Never signal a reused leader identity or a process outside our session."""
    for number, grace in ((signal.SIGTERM, 10), (signal.SIGKILL, 10)):
        members = group_members(child.pid)
        require(not any(p["pid"] == child.pid and p["start_ticks"] != start_ticks for p in members),
                "owned process-group leader PID was reused; refusing to signal")
        if members:
            os.killpg(child.pid, number)
        deadline = time.monotonic() + grace
        while time.monotonic() < deadline:
            # waitpid(-pgid) can also select the leader. Do not use it until
            # Popen has already reaped that exact child and saved its status.
            if child.poll() is not None:
                while True:
                    try:
                        pid, status = os.waitpid(-child.pid, os.WNOHANG)
                    except ChildProcessError:
                        break
                    if pid == 0:
                        break
                    reaped.append({"pid": pid, "wait_status": status})
            if child.poll() is not None and not group_members(child.pid):
                return
            time.sleep(0.1)
    raise RuntimeError("owned process group did not terminate within bounded cleanup grace")


def load_gate():
    spec = importlib.util.spec_from_file_location("local_mandatory_gate", ROOT / "scripts/verify-build-gate.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def cache():
    values = {}
    for line in (BUILD / "CMakeCache.txt").read_text().splitlines():
        match = re.match(r"([^#/:][^:]*):[^=]*=(.*)$", line)
        if match:
            values[match[1]] = match[2]
    expected = {"BUILD_TESTS": "ON", "SIRIUS_MANDATORY_TESTS": "ON",
                "SIRIUS_REQUIRE_VULKAN_RUNTIME": "ON", "SIRIUS_ALIGNMENT_MODE": "qualification",
                "SIRIUS_CONTRACT_MODE": "2", "SIRIUS_WERROR": "ON",
                "SIRIUS_BUILD_VIEWER": "ON", "SIRIUS_SANITIZERS": "none", "CMAKE_BUILD_TYPE": "Release"}
    require(all(values.get(key) == value for key, value in expected.items()), "CMake qualification/runtime/Mandatory flags differ")
    require(values.get("CMAKE_HOME_DIRECTORY") == str(ROOT), "CMake cache names another source checkout")
    return values


def authority(expected_revision, gate):
    lines = [line for line in (BUILD / "build.ninja").read_text().splitlines()
             if "  COMMAND = " in line and str(ROOT / "scripts/verify-build-gate.py") in line and "--action run --stamp" in line]
    require(len(lines) == 1 and "$" not in lines[0], "generated Mandatory invocation is not unique/literal")
    words = shlex.split(lines[0].split(" = ", 1)[1])
    first = words.index(str(ROOT / "scripts/verify-build-gate.py"))
    arguments = words[first + 1:words.index("&&", first)]

    def values(flag):
        return [arguments[index + 1] for index, word in enumerate(arguments) if word == flag]

    expected = {"--action": "run", "--alignment-mode": "qualification", "--source-root": str(ROOT),
                "--build-dir": str(BUILD), "--source-revision": expected_revision, "--source-tree-clean": "true",
                "--stamp": str(BUILD / "generated/sirius/mandatory_gate.json"), "--config": "Release"}
    require(all(values(flag) == [value] for flag, value in expected.items()), "generated gate invocation is stale or altered")
    collections = {key: gate.parse_artifacts(values(flag), flag) for key, flag in (
        ("tested_artifacts", "--tested-artifact"), ("product_artifacts", "--product-artifact"),
        ("test_input_artifacts", "--test-input-artifact"))}
    require(set(collections["tested_artifacts"]) == gate.TESTED_ARTIFACTS and
            set(collections["product_artifacts"]) == gate.RELEASE_PRODUCTS,
            "generated gate omitted complete executable/product topology")
    gate.require_test_input_set(collections["test_input_artifacts"], "qualification", collections["product_artifacts"])
    records = {key: gate.artifact_records(paths, ROOT, BUILD) for key, paths in collections.items()}
    gate.validate_test_input_records(records["test_input_artifacts"], "qualification", records["product_artifacts"])
    return arguments, collections, records, Path(values("--ctest")[0]), Path(words[first - 1])


def environment():
    env = os.environ.copy()
    allowed = {"VK_DRIVER_FILES": str(ICD), "SIRIUS_VULKAN_DEVICE": "0"}
    forbidden = {key: value for key, value in env.items() if value and
                 (key.startswith(("SIRIUS_", "VK_", "LP_", "GALLIVM_", "MESA_", "GTEST_", "CTEST_"))
                  or key in {"LD_PRELOAD", "LD_AUDIT", "LD_LIBRARY_PATH"}) and allowed.get(key) != value}
    require(not forbidden, f"inherited test/precision/driver/tuning overrides refused: {sorted(forbidden)}")
    env.update(allowed)
    return env


def inventory(ctest, env, destination, gate):
    raw = query([str(ctest), "--test-dir", str(BUILD), "-C", "Release", "--show-only=json-v1"], env)
    destination.write_bytes(raw)
    document = json.loads(raw)
    names, digest_value = gate.inspect_inventory(document)
    return names, digest_value


def run(expected_revision, name):
    require(re.fullmatch(r"[0-9a-f]{40,64}", expected_revision), "expected revision must be full lowercase hexadecimal")
    require(re.fullmatch(r"[A-Za-z0-9_-]+", name), "unsafe diagnostic destination name")
    folder = Path(__file__).parent / name
    folder.mkdir(exist_ok=False)  # Refuse overwrite/relaunch of this observation.
    child = None
    start_ticks = None
    reaped = []
    report = {"scope": "complete existing Mandatory estate; no overall release/ultimate qualification claim",
              "status": "preflight", "expected_revision": expected_revision,
              "supervisor": proc(os.getpid()), "boot_id": Path("/proc/sys/kernel/random/boot_id").read_text().strip(),
              "runner": identity(__file__), "started_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
              "timeout_seconds": LIMIT_SECONDS, "maximum_rss_kib": RSS_KIB,
              "minimum_available_kib": RESERVE_KIB, "sample_seconds": SAMPLE_SECONDS,
              "rss_scope": "supervisor plus owned child process group", "passed": False}
    write(folder / "progress.json", report)
    try:
        env = environment()
        source = source_identity()
        require(source == {"revision": expected_revision, "status": ""}, "source must be clean at the expected revision")
        cache_values = cache()
        gate = load_gate()
        invocation, collections, records, ctest, producer_python = authority(expected_revision, gate)
        require(str(ctest) == cache_values["CMAKE_CTEST_COMMAND"], "gate and cache CTest differ")
        cmake = Path(shutil.which("cmake")).resolve(strict=True)
        require(str(cmake) == cache_values["CMAKE_COMMAND"], "selected CMake differs from configured CMake")
        names, inventory_digest = inventory(ctest, env, folder / "inventory-before.json", gate)
        ctypes.CDLL(str(LVP))  # Resolve/bind real stock LVP/LLVM without creating a Vulkan instance.
        mapped = {Path(line.split()[-1]).resolve(strict=True) for line in Path("/proc/self/maps").read_text().splitlines()
                  if len(line.split()) >= 6 and ("libvulkan_lvp.so" in line or "libLLVM.so" in line)}
        require(LVP.resolve() in mapped, "stock LVP mapping is absent")
        icd_document = json.loads(ICD.read_text())
        require(icd_document["ICD"]["library_path"] in {LVP.name, str(LVP)}, "stock ICD does not select stock LVP")
        tracked = [ROOT / os.fsdecode(value) for value in query(["git", "ls-files", "-z"]).split(b"\0") if value]
        modules = sorted((BUILD / "src/sirius/backend/retained").glob("retained_*.spv"))
        require(len(modules) == 24, "expected all 24 retained modules")
        files = [*tracked, *mapped, ICD, cmake, ctest, producer_python, Path(sys.executable), BUILD / "build.ninja", BUILD / "CMakeCache.txt",
                 BUILD / "compile_commands.json", BUILD / "generated/sirius/base/operating_model_embedded.h",
                 BUILD / "src/sirius/backend/retained/retained_kernels.h", *modules,
                 *BUILD.rglob("CTestTestfile.cmake"), *BUILD.rglob("*_tests.cmake"), *BUILD.rglob("*_include.cmake")]
        for paths in collections.values():
            files.extend(paths.values())
        exe = collections["tested_artifacts"]["sirius"]
        resources = {key: exe.parent / "resources" / Path(path).relative_to("share/sirius")
                     for key, path in gate.INSTALLED_PRODUCTS.items()}
        files.extend(resources.values())
        before = {str(path): identity(path) for path in dict.fromkeys(files)}
        preflight = {"source": source, "cache_required": {key: cache_values[key] for key in (
            "BUILD_TESTS", "SIRIUS_MANDATORY_TESTS", "SIRIUS_REQUIRE_VULKAN_RUNTIME", "SIRIUS_ALIGNMENT_MODE",
            "SIRIUS_CONTRACT_MODE", "SIRIUS_WERROR", "SIRIUS_BUILD_VIEWER", "SIRIUS_SANITIZERS", "CMAKE_BUILD_TYPE")},
            "gate_invocation": invocation, "gate_artifact_records": records, "bound_files": before,
            "registered": len(names), "inventory_sha256": inventory_digest,
            "environment": {key: value for key, value in env.items() if key.startswith(("SIRIUS_", "VK_", "LP_", "GALLIVM_", "MESA_", "GTEST_", "CTEST_"))
                            or key in {"PATH", "LANG", "LC_ALL", "DISPLAY", "WAYLAND_DISPLAY", "XDG_RUNTIME_DIR", "LD_LIBRARY_PATH"}}}
        write(folder / "preflight.json", preflight)
        require(available_kib() >= RESERVE_KIB, "host memory reserve is already below guard")
        libc = ctypes.CDLL(None, use_errno=True)
        require(libc.prctl(36, 1, 0, 0, 0) == 0, "cannot become child subreaper")
        command = [str(cmake), "--build", "--preset", "linux-gcc", "--target", "RunMandatoryTests", "-j2"]
        report.update({"status": "running", "command": command, "preflight_sha256": digest(folder / "preflight.json"),
                       "peak_sampled_rss_kib": 0, "minimum_sampled_available_kib": None, "sampled_process_identities": {}})
        started = time.monotonic()
        for number in (signal.SIGTERM, signal.SIGINT, signal.SIGHUP):
            signal.signal(number, signal_stop)
        with (folder / "target.log").open("wb") as output, (folder / "samples.jsonl").open("w") as samples:
            child = subprocess.Popen(command, cwd=ROOT, env=env, stdout=output, stderr=subprocess.STDOUT, start_new_session=True)
            start_ticks = proc(child.pid)["start_ticks"]
            report.update({"child_pid": child.pid, "child_start_ticks": start_ticks, "owned_pgid": child.pid})
            reason = None
            while child.poll() is None:
                members = group_members(child.pid)
                available = available_kib()
                rss = sum(item["rss_kib"] for item in members) + proc(os.getpid())["rss_kib"]
                elapsed = time.monotonic() - started
                report["peak_sampled_rss_kib"] = max(report["peak_sampled_rss_kib"], rss)
                report["minimum_sampled_available_kib"] = min(report["minimum_sampled_available_kib"] or available, available)
                for item in members:
                    report["sampled_process_identities"][f'{item["pid"]}:{item["start_ticks"]}'] = item
                sample = {"elapsed_seconds": elapsed, "available_kib": available, "rss_kib": rss, "members": members}
                samples.write(json.dumps(sample) + "\n")
                samples.flush()
                report["wall_seconds"] = elapsed
                write(folder / "progress.json", report)
                reason = STOP or ("diagnostic_timeout" if elapsed >= LIMIT_SECONDS else
                                  "diagnostic_process_rss_guard" if rss > RSS_KIB else
                                  "diagnostic_host_memory_guard" if available < RESERVE_KIB else None)
                if reason:
                    break
                time.sleep(SAMPLE_SECONDS)
            remaining = group_members(child.pid)
            if reason is None and any(item["state"] != "Z" for item in remaining):
                reason = "unexpected_owned_descendants_after_target_exit"
            report["stop_reason"] = reason
            write(folder / "progress.json", report)
            if reason or remaining:
                terminate_group(child, start_ticks, reaped)
            report["child_exit"] = child.wait(timeout=10)
            report["wall_seconds"] = time.monotonic() - started
        report["remaining_owned_processes"] = group_members(child.pid)
        report["source_after"] = source_identity()
        after = {path: identity(path) for path in before}
        write(folder / "bindings-after.json", after)
        report["source_unchanged"] = report["source_after"] == source
        report["artifacts_unchanged"] = after == before
        post_names, post_digest = inventory(ctest, env, folder / "inventory-after.json", gate)
        report["inventory_unchanged"] = post_names == names and post_digest == inventory_digest
        require(report["child_exit"] == 0 and reason is None and not report["remaining_owned_processes"], "actual target did not finish successfully without guard/interruption")
        require(report["source_unchanged"] and report["artifacts_unchanged"] and report["inventory_unchanged"], "bound source/artifacts/registration changed")
        stamp = BUILD / "generated/sirius/mandatory_gate.json"
        receipt = gate.validate_document(json.loads(stamp.read_text()))
        require(receipt["source"] == {"revision": expected_revision, "clean": True} and receipt["alignment_mode"] == "qualification", "producer receipt source/mode differ")
        require(all(receipt[key] == value for key, value in records.items()), "producer receipt differs from preflight original artifacts")
        gate.check_gate(argparse.Namespace(source_root=ROOT, build_dir=BUILD, stamp=stamp, alignment_mode="qualification",
                                          source_revision=expected_revision, ctest=ctest, config="Release", installed_root=None))
        staged = exe.parent / "resources/model/mandatory_gate.json"
        require(staged.read_bytes() == stamp.read_bytes(), "actual target did not stage its exact successful receipt")
        for key, path in resources.items():
            actual = identity(path)
            require(all(actual[field] == receipt["product_artifacts"][key][field] for field in ("bytes", "sha256")), f"deployed runtime product differs: {key}")
        runtime_path = BUILD / "generated/sirius/mandatory_runtime_identity_gate.json"
        runtime_spec = importlib.util.spec_from_file_location("local_runtime_identity", ROOT / "scripts/runtime_identity.py")
        runtime_validator = importlib.util.module_from_spec(runtime_spec)
        runtime_spec.loader.exec_module(runtime_validator)
        runtime_record = runtime_validator.validate_record(
            json.loads(runtime_path.read_text()), expected_revision,
            records["tested_artifacts"]["sirius"], required=True, clean=True)
        require(all(runtime_record["environment"].get(key) == env[key]
                    for key in ("VK_DRIVER_FILES", "SIRIUS_VULKAN_DEVICE")),
                "runtime identity differs from the supervised provider selection")
        report["runtime_identity"] = identity(runtime_path)
        report.update({"passed": True, "ctest": receipt["ctest"], "producer_receipt": identity(stamp), "staged_receipt": identity(staged)})
    except BaseException as error:
        report["error"] = f"{type(error).__name__}: {error}"
        if child is not None:
            try:
                terminate_group(child, start_ticks, reaped)
                report["child_exit"] = child.wait(timeout=10)
                report["remaining_owned_processes"] = group_members(child.pid)
            except BaseException as cleanup_error:
                report["cleanup_error"] = f"{type(cleanup_error).__name__}: {cleanup_error}"
    finally:
        report["adopted_children_reaped"] = reaped
        report["status"] = "terminal"
        report["finished_utc"] = datetime.datetime.now(datetime.timezone.utc).isoformat()
        report["retained_outputs"] = {}
        report["preserved_producer_outputs"] = {}
        for path in (
            folder / "target.log", BUILD / "generated/sirius/mandatory_gate.json",
            BUILD / "generated/sirius/mandatory_runtime_identity_gate.json",
            BUILD / "generated/sirius/mandatory_gate_junit.xml", BUILD / "generated/sirius/mandatory_gate_ctest.log"):
            try:
                if path.is_file():
                    report["retained_outputs"][str(path)] = identity(path)
                    if path.parent != folder:
                        preserved = folder / ("producer-" + path.name)
                        shutil.copyfile(path, preserved)
                        report["preserved_producer_outputs"][str(path)] = identity(preserved)
            except OSError as error:
                report["passed"] = False
                report.setdefault("output_binding_errors", []).append(f"{path}: {error}")
        write(folder / "report.json", report)
        write(folder / "progress.json", report)
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--expected-revision", required=True)
    parser.add_argument("--name", required=True)
    options = parser.parse_args()
    raise SystemExit(run(options.expected_revision, options.name))
