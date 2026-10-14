#!/usr/bin/env python3
"""Root-operated one-shot phases; postprocessing occurs after owner terminal."""
from pathlib import Path
import collections
import hashlib
import json
import math
import os
import re
import runpy
import shutil
import subprocess
import sys

sys.dont_write_bytecode = True
WORK = Path(__file__).resolve().parent
ROOT = WORK.parents[1]
BUILD = ROOT / "bin/linux-gcc/p4-parity-setup-diagnostic"
PHASES = ["compile-device", "compile-driver", "link", "run"]


def require(ok, message):
    if not ok:
        raise ValueError(message)


def identity(path):
    path = Path(path)
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return {"bytes": path.stat().st_size, "sha256": digest.hexdigest(), "resolved_path": str(path.resolve(strict=True))}


def key(path):
    return str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path)


def query(argv):
    return subprocess.check_output(argv, cwd=ROOT, timeout=30, text=True)


def accepted_owner(path):
    owner = json.loads(path.read_text())
    require(owner["passed"] is True and owner["status"] == "terminal" and owner["child_exit"] == 0
            and owner["stop_reason"] is None and owner["cleanup_errors"] == []
            and owner["source_unchanged"] and owner["whole_inputs_unchanged"]
            and owner["remaining_owned_processes"] == [] and owner["observed_births"]
            and all(row["absent"] for row in owner["observed_births"]), "original owner completion/seals/retirement incomplete")
    return owner


def check_dependencies(path, prepared):
    text = path.read_text().replace("\\\n", "")
    require(":" in text, "missing compiler dependency rule")
    import shlex
    dependencies = [(Path(p) if Path(p).is_absolute() else ROOT / p).resolve() for p in shlex.split(text.split(":", 1)[1])]
    require(dependencies, "empty actual compiler footprint")
    for dependency in dependencies:
        name = key(dependency)
        require(name in prepared["pins"], "actual compiler dependency was not frozen:" + name)
        require(identity(dependency) == prepared["pins"][name], "actual compiler dependency changed:" + name)
    return [key(p) for p in dependencies]


def strong_symbols(path):
    text = query(["/usr/bin/nm", "-g", "--defined-only", str(path)])
    return {row.split()[2] for row in text.splitlines() if len(row.split()) == 3 and row.split()[1] in ("T", "D", "B", "R")}, text


def review_link():
    required = set(json.loads((WORK / "original-symbols.json").read_text())["required_strong_sirius_definitions"])
    object_symbols, object_nm = strong_symbols(BUILD / "instrumented-vulkan_device.o")
    executable_symbols, executable_nm = strong_symbols(BUILD / "setup_probe")
    require(required <= object_symbols and required <= executable_symbols, "one complete replacement TU does not provide every original strong Sirius definition")
    link_map = (BUILD / "setup_probe.map").read_text()
    original = ROOT / "bin/linux-gcc/src/sirius/backend/libsirius_backend.a"
    require(str(original) + "(vulkan_device.cpp.o)" not in link_map, "old archive Vulkan-device member was extracted")
    direct = "LOAD " + str(BUILD / "instrumented-vulkan_device.o")
    archive = "LOAD " + str(original)
    require(direct in link_map and archive in link_map and link_map.index(direct) < link_map.index(archive), "complete private TU must precede unchanged backend archive")
    dynamic = query(["/usr/bin/readelf", "-d", str(BUILD / "setup_probe")])
    require("libvulkan.so.1" in dynamic and not any(flag in dynamic for flag in ("(RPATH)", "(RUNPATH)")), "unexpected loader search wiring")
    report = {"passed": True, "required_strong_sirius_definitions": sorted(required), "original_vulkan_archive_member_extracted": False,
              "instrumented_object_before_original_archive": True, "instrumented_object_nm": object_nm,
              "executable_nm": executable_nm, "dynamic_section": dynamic,
              "link_map": identity(BUILD / "setup_probe.map"), "executable": identity(BUILD / "setup_probe")}
    (WORK / "link-review.json").write_text(json.dumps(report, indent=2) + "\n")
    return identity(WORK / "link-review.json")


def read_markers(folder):
    markers, stack = [], []
    for line in (folder / "stderr.log").read_text().splitlines():
        if not line.startswith("[SetupPhase] "):
            continue
        marker = json.loads(line[len("[SetupPhase] "):])
        require(marker["event"] in ("begin", "end", "failure") and isinstance(marker["monotonic_ns"], int), "malformed marker")
        require(not markers or marker["monotonic_ns"] >= markers[-1]["monotonic_ns"], "marker time order changed")
        if marker["event"] == "begin":
            stack.append(marker["phase"])
        elif marker["event"] == "end":
            require(stack and stack[-1] == marker["phase"], "marker scope nesting changed")
            stack.pop()
        markers.append(marker)
    return {"markers": markers, "unclosed_phases": stack, "last_marker": markers[-1] if markers else None}


def review_run(folder, owner, prepared):
    mappings = {}
    for observation in owner["provider_mappings"]:
        require(observation["owned_session"] == owner["child_pid"], "map ownership mismatch")
        for line in observation["maps"]:
            fields = line.split(None, 5)
            require(len(fields) == 6 and fields[5] in prepared["required_provider_mappings"], "unexpected provider map")
            path = fields[5]
            info = Path(path).stat()
            require(tuple(int(n, 16) for n in fields[3].split(":")) == (os.major(info.st_dev), os.minor(info.st_dev)) and int(fields[4]) == info.st_ino, "mapped module dev/inode changed")
            mappings[path] = identity(path)
    require(set(mappings) == set(prepared["required_provider_mappings"]) and all(pin == prepared["pins"][path] for path, pin in mappings.items()), "actual loader/LVP/LLVM join incomplete")
    observation = read_markers(folder)
    require(not observation["unclosed_phases"], "setup phase has no terminal end")
    markers = observation["markers"]
    began = collections.Counter(row["phase"] for row in markers if row["event"] == "begin")
    for phase in ("EnumerateVulkanDevices", "ResolveVulkanDeviceIndex", "CreateVulkanDevice", "InitialisePipelineCache", "LoadKernel", "vkCreateShaderModule", "SetBufferAllocationLimit", "DispatchCommands", "GetOrCreatePipeline", "vkCreateComputePipelines", "vkQueueSubmit", "vkQueueWaitIdle", "ReadBuffer", "VulkanDevice.destructor", "vkDeviceWaitIdle"):
        require(began[phase] == 1, "public/native call count differs:" + phase)
    require(began["CreateBuffer"] == 3 and began["WriteBuffer"] == 2, "original three-buffer/two-write sequence differs")
    native = [row for row in markers if row["event"] == "end" and row["phase"].startswith("vk")]
    require(len(native) == 5 and all(row.get("has_result") is True and row["result"] == 0 for row in native), "an original native call failed")
    lines = [line for line in (folder / "stdout.log").read_text().splitlines() if line.startswith("{")]
    require(len(lines) == 1, "setup summary absent/ambiguous")
    summary = json.loads(lines[0])
    expected = {"schema": "sirius-parity-setup-only-v1", "setup_completed": True, "active_records": 0,
                "fixture_records": 4, "flat_words": 64, "completed_words": 256, "result_words": 256,
                "module_words": 1390559, "device_index": 0,
                "numerical_acceptance": False, "qualification_claimed": False}
    require(all(summary.get(k) == value for k, value in expected.items()), "setup shape/scope differs")
    require(1296 <= summary["allocated_bytes"] <= 1024 * 1024, "original actual-allocation bounds failed")
    require("llvmpipe" in summary["device_name"].lower() and summary["driver_name"] == "llvmpipe", "software provider identity differs")
    require(all(math.isfinite(summary[k]) and summary[k] >= 0 for k in ("whole_setup_seconds_including_teardown", "pipeline_setup_ms", "submit_wait_ms", "dispatch_total_ms")) and summary["submit_wait_ms"] > 0, "original initialization timing observations invalid")
    return {"summary": summary, "phase_observation": observation, "actual_provider_mappings": mappings}


def main():
    require(Path.cwd().resolve() == ROOT and len(sys.argv) == 2 and sys.argv[1] in PHASES, "select one declared phase from checkout")
    require(sys.flags.optimize == 0 and not os.environ.get("PYTHONOPTIMIZE"), "owner assertions must be enabled")
    phase = sys.argv[1]
    prepared = json.loads((WORK / "preparation.json").read_text())
    result_path = WORK / (phase + "-result.json")
    folder = WORK / phase
    require(not folder.exists() and not result_path.exists(), "one-shot phase already attempted")
    result = {"passed": False, "phase": phase, "scope": prepared["scope"], "numerical_acceptance": False, "qualification_claimed": False}
    def save():
        result_path.write_text(json.dumps(result, indent=2) + "\n")
    save()
    try:
        require(query(["git", "rev-parse", "HEAD"]).strip() == prepared["runner_source_revision"] and not query(["git", "status", "--porcelain"]).strip(), "runner source changed")
        require(set(query(["git", "diff", "--name-only", prepared["product_source_revision"], "HEAD"]).splitlines()) == set(prepared["diagnostic_only_difference"]), "product difference changed")
        for path, pin in prepared["pins"].items():
            require(identity(ROOT / path) == pin, "prepared input changed:" + path)
        transform = json.loads((WORK / "instrumentation-transform.json").read_text())
        text = (WORK / "recovery/original-vulkan_device.cpp").read_text()
        for row in transform["changes"]:
            require(text.count(row["before"]) == 1, "source instrumentation anchor changed")
            text = text.replace(row["before"], row["after"], 1)
        require(text == (WORK / "instrumented-vulkan_device.cpp").read_text(), "instrumented TU contains another change")
        for prior in PHASES[:PHASES.index(phase)]:
            row = json.loads((WORK / (prior + "-result.json")).read_text())
            require(row["passed"] is True, "preceding phase incomplete:" + prior)
            accepted_owner(WORK / prior / "owner.json")
            require(identity(WORK / prior / "owner.json") == row["owner"], "preceding owner changed")
            for path, pin in row.get("products", {}).items():
                require(identity(ROOT / path) == pin, "preceding product changed:" + path)
        if phase == "compile-device":
            require(not BUILD.exists(), "private product root is not fresh")
            (BUILD / "resources/kernels").mkdir(parents=True)
            shutil.copyfile(ROOT / prepared["staged_shader"], ROOT / prepared["runtime_shader"])
        require(BUILD.is_dir() and not BUILD.is_symlink(), "private product root invalid")
        source_pin = identity(ROOT / prepared["staged_shader"])
        runtime_pin = identity(ROOT / prepared["runtime_shader"])
        require(all(source_pin[k] == runtime_pin[k] for k in ("bytes", "sha256")), "runtime shader differs from full stock")
        if phase == "run":
            link = json.loads((WORK / "link-result.json").read_text())
            require(identity(WORK / "link-review.json") == link["link_review"], "terminal link review changed")
            require("SIRIUS_PRECISION" not in os.environ, "preserve unset original precision environment")
        temporary = WORK / (phase + "-temporary")
        require(not temporary.exists(), "TMPDIR already exists")
        temporary.mkdir(mode=0o700)
        require(temporary.stat().st_uid == os.getuid() and temporary.stat().st_mode & 0o777 == 0o700, "TMPDIR ownership/mode")
        os.environ["TMPDIR"] = str(temporary)
        old = sys.argv
        try:
            sys.argv = [str(WORK / "run_owned.py"), phase]
            runpy.run_path(str(WORK / "run_owned.py"), run_name="__main__")
        finally:
            sys.argv = old
        owner = accepted_owner(folder / "owner.json")
        products = {path: identity(ROOT / path) for path in prepared["products"].get(phase, [])}
        if phase.startswith("compile-"):
            result["actual_compile_dependencies"] = check_dependencies(ROOT / prepared["products"][phase][1], prepared)
        if phase == "link":
            result["link_review"] = review_link()
        if phase == "run":
            result["setup_observation"] = review_run(folder, owner, prepared)
        result.update(passed=True, owner=identity(folder / "owner.json"), products=products, runtime_shader=runtime_pin)
        save()
        print(json.dumps({k: result[k] for k in ("passed", "phase", "scope", "numerical_acceptance", "qualification_claimed")}, indent=2), flush=True)
        return 0
    except BaseException as error:
        result["error"] = repr(error)
        if phase == "run" and (folder / "stderr.log").exists():
            try:
                result["partial_setup_observation"] = read_markers(folder)
            except BaseException as marker_error:
                result["marker_read_error"] = repr(marker_error)
        save()
        raise


if __name__ == "__main__":
    raise SystemExit(main())
