#!/usr/bin/env python3
"""Source-only preparation; never invokes a compiler or provider."""
from pathlib import Path
import difflib
import hashlib
import json
import re
import shlex
import shutil
import subprocess
import sys

sys.dont_write_bytecode = True
ROOT = Path.cwd().resolve()
WORK = Path(__file__).resolve().parent
BUILD = ROOT / "bin/linux-gcc/p4-parity-setup-diagnostic"
REVISION = "f5d649bf5067d1d68dc9215462f37830122f0540"
PRODUCT = "4bcd9875cd9308deabed06a5fa2f17d2321b0cbc"
SOURCE = ROOT / "src/sirius/backend/vulkan/vulkan_device.cpp"
ORIGINAL_OBJECT = ROOT / "bin/linux-gcc/src/sirius/backend/CMakeFiles/sirius_backend.dir/vulkan/vulkan_device.cpp.o"
ARCHIVES = [ROOT / "bin/linux-gcc/src/sirius/backend/libsirius_backend.a", ROOT / "bin/linux-gcc/src/sirius/core/libsirius_core.a", ROOT / "bin/linux-gcc/src/sirius/libsirius_base.a"]
BASE = ROOT / "attestations/alignment/4bcd987/gcc-revision-synchronisation"
PRODUCER = BASE / "files/out/gcc-revision-synchronisation"
OLD = ROOT / "attestations/scientific-core/4bcd987/p4-cpu-radiation-current-source/files/out/p4-cpu-radiation-current-source"
MAPPINGS = ["/usr/lib/x86_64-linux-gnu/libvulkan_lvp.so", "/usr/lib/x86_64-linux-gnu/libLLVM.so.20.1", "/usr/lib/x86_64-linux-gnu/libvulkan.so.1.3.275"]


def identity(path):
    path = Path(path)
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return {"bytes": path.stat().st_size, "sha256": digest.hexdigest(), "resolved_path": str(path.resolve(strict=True))}


def key(path):
    return str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path)


def write(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def main():
    assert WORK == ROOT / "out/p4-parity-setup-diagnostic" and not BUILD.exists()
    assert not (WORK / "preparation.json").exists()
    assert subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip() == REVISION
    assert not subprocess.check_output(["git", "status", "--porcelain"], text=True).strip()
    recovery = WORK / "recovery"
    recovery.mkdir()
    raw = SOURCE.read_text()
    assert subprocess.check_output(["git", "show", PRODUCT + ":" + key(SOURCE)]) == raw.encode()
    (recovery / "original-vulkan_device.cpp").write_text(raw)
    changes = []
    def replace(before, after):
        nonlocal raw
        assert raw.count(before) == 1, before
        changes.append({"before": before, "after": after})
        raw = raw.replace(before, after, 1)
    replace("#include <cstdlib>", "#include <cstdlib>\n#include <cstdio>")
    helper = '''// Diagnostic markers only; no Vulkan call or error handling is reordered.
struct SetupPhaseScope {
    const char* phase;
    bool has_result = false;
    int result = 0;
    static void Print(const char* event, const char* label, bool has_result, int result) {
        const auto ticks = std::chrono::duration_cast<std::chrono::nanoseconds>(
            std::chrono::steady_clock::now().time_since_epoch()).count();
        std::fprintf(stderr, "[SetupPhase] {\\"event\\":\\"%s\\",\\"phase\\":\\"%s\\",\\"monotonic_ns\\":%lld,\\"has_result\\":%s,\\"result\\":%d}\\n",
                     event, label, static_cast<long long>(ticks), has_result ? "true" : "false", result);
        std::fflush(stderr);
    }
    explicit SetupPhaseScope(const char* value) : phase(value) { Print("begin", phase, false, 0); }
    ~SetupPhaseScope() { Print("end", phase, has_result, result); }
};
template <typename Function>
auto SetupMeasuredCall(const char* phase, Function&& function) {
    SetupPhaseScope observation(phase);
    const auto result = function();
    observation.has_result = true;
    observation.result = static_cast<int>(result);
    return result;
}
'''
    replace("namespace {\n\nusing base::ErrorDomain;", "namespace {\n\n" + helper + "\nusing base::ErrorDomain;")
    scopes = [
      ("Expected<std::vector<DeviceInfo>> EnumerateVulkanDevices() {", "EnumerateVulkanDevices"),
      ("Expected<std::size_t> ResolveVulkanDeviceIndex(std::span<const DeviceInfo> devices) {", "ResolveVulkanDeviceIndex"),
      ("Expected<std::unique_ptr<ComputeDevice>> CreateVulkanDevice(std::size_t index) {", "CreateVulkanDevice"),
      ("VulkanDevice::~VulkanDevice() {", "VulkanDevice.destructor"),
      ("void VulkanDevice::InitialisePipelineCache() {", "InitialisePipelineCache"),
      ("Expected<KernelHandle> VulkanDevice::LoadKernel(std::span<const std::uint32_t> spirv) {", "LoadKernel"),
      ("Expected<void> VulkanDevice::SetBufferAllocationLimit(std::uint64_t bytes) {", "SetBufferAllocationLimit"),
      ("Expected<BufferHandle> VulkanDevice::CreateBuffer(std::uint64_t size_bytes, BufferUsage usage) {", "CreateBuffer"),
      ("Expected<void> VulkanDevice::WriteBuffer(BufferHandle handle, std::span<const std::byte> data) {", "WriteBuffer"),
      ("Expected<void> VulkanDevice::ReadBuffer(BufferHandle handle, std::span<std::byte> out) {", "ReadBuffer"),
      ("    KernelHandle kernel, std::span<const BufferHandle> buffers, bool* created) {", "GetOrCreatePipeline"),
      ("                                              std::uint32_t* pipeline_creations) {", "DispatchCommands"),
    ]
    for anchor, phase in scopes:
        replace(anchor, anchor + '\n    SetupPhaseScope setup_phase_scope("' + phase + '");')
    calls = [
      ("vkCreateShaderModule(device_, &create_info, nullptr, &kernel.module)", "vkCreateShaderModule"),
      ("vkCreateComputePipelines(device_, pipeline_cache_, 1, &pipeline_info,\n                                                    nullptr, &pipeline.pipeline)", "vkCreateComputePipelines"),
      ("vkQueueSubmit(queue_, 1, &submit_info, VK_NULL_HANDLE)", "vkQueueSubmit"),
      ("vkQueueWaitIdle(queue_)", "vkQueueWaitIdle"),
      ("vkDeviceWaitIdle(device_)", "vkDeviceWaitIdle"),
    ]
    for call, phase in calls:
        replace(call, 'SetupMeasuredCall("' + phase + '", [&] { return ' + call + '; })')
    instrumented = WORK / "instrumented-vulkan_device.cpp"
    instrumented.write_text(raw)
    original = (recovery / "original-vulkan_device.cpp").read_text()
    inverse = raw
    for change in reversed(changes):
        assert inverse.count(change["after"]) == 1
        inverse = inverse.replace(change["after"], change["before"], 1)
    assert inverse == original
    write(WORK / "instrumentation-transform.json", {"original": identity(SOURCE), "instrumented": identity(instrumented), "changes": changes, "reverse_transform_recovers_original_exactly": True})
    (WORK / "instrumentation.patch").write_text("".join(difflib.unified_diff(original.splitlines(True), raw.splitlines(True), fromfile=key(SOURCE), tofile=key(instrumented))))
    commands = json.loads((ROOT / "bin/linux-gcc/compile_commands.json").read_text())
    command = next(row for row in commands if row["file"] == str(SOURCE))
    tokens = shlex.split(command["command"])
    assert tokens[-2:] == ["-c", str(SOURCE)] and tokens[-4] == "-o"
    flags = tokens[:-4]
    assert all(word in flags for word in ("-O3", "-DNDEBUG", "-std=c++2c", "-Werror", "-DSIRIUS_CONTRACT_MODE=2"))
    assert flags[0] == "/usr/bin/g++-14"
    member = subprocess.check_output(["/usr/bin/ar", "p", str(ARCHIVES[0]), "vulkan_device.cpp.o"])
    assert member == ORIGINAL_OBJECT.read_bytes()
    symbols_text = subprocess.check_output(["/usr/bin/nm", "-g", "--defined-only", str(ORIGINAL_OBJECT)], text=True)
    strong = sorted(row.split()[2] for row in symbols_text.splitlines() if len(row.split()) == 3 and row.split()[1] in ("T", "D", "B", "R") and "sirius" in row.split()[2])
    assert len(strong) >= 20
    write(WORK / "original-symbols.json", {"object": identity(ORIGINAL_OBJECT), "archive": identity(ARCHIVES[0]), "member": "vulkan_device.cpp.o", "member_sha256": hashlib.sha256(member).hexdigest(), "required_strong_sirius_definitions": strong, "nm_stdout": symbols_text})
    # Ninja's read-only deps tool reads existing metadata; it does not build or
    # run any target. Preserve the actual prior-TU header footprint for resealing.
    dependencies = set()
    dependency_text = {}
    for target in (command["output"], "tests/backend/CMakeFiles/sirius_backend_tests.dir/kernel_parity_test.cpp.o"):
        output = subprocess.check_output(["/usr/bin/ninja", "-C", str(ROOT / "bin/linux-gcc"), "-t", "deps", target], text=True)
        assert "deps not found" not in output and "#deps" in output, target
        dependency_text[target] = output
        dependencies.update(Path(line.strip()).resolve() for line in output.splitlines() if line.startswith("    /"))
    assert SOURCE in dependencies and ROOT / "tests/support/test_resource.h" in dependencies
    write(WORK / "original-dependencies.json", dependency_text)
    sources = {p for p in dependencies if p.is_relative_to(ROOT) and not p.is_relative_to(ROOT / "bin")}
    sources.update(ROOT / name for name in ("tests/backend/kernel_parity_test.cpp", "tests/support/test_resource.h", "src/sirius/base/resource_locator.cpp", "src/sirius/backend/CMakeLists.txt", "tests/backend/CMakeLists.txt"))
    before = json.loads((PRODUCER / "build/inputs-before.json").read_text())
    after = json.loads((PRODUCER / "build/inputs-after.json").read_text())
    source_records = []
    for path in sorted(sources):
        name = key(path)
        content = subprocess.check_output(["git", "show", PRODUCT + ":" + name])
        assert path.read_bytes() == content and before[name] == after[name] == identity(path), name
        destination = recovery / "source" / name
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(content)
        source_records.append({"path": name, "revision": PRODUCT, "sha256": hashlib.sha256(content).hexdigest(), "recovery": key(destination)})
    shader = ROOT / "bin/linux-gcc/tests/backend/resources/kernels/parity_probe.spv"
    assert identity(shader)["sha256"] == "19991852c2bf8025aaa30b4cb3029e7bb980425d74de515efc4879df95eb0891"
    staged = WORK / "resources/kernels/parity_probe.spv"
    staged.parent.mkdir(parents=True)
    shutil.copyfile(shader, staged)
    owner = (OLD / "run_owned.py").read_text()
    assert hashlib.sha256(owner.encode()).hexdigest() == "b47e4bb22e63c16244802a5d121160a63b942f5c37614d26ded319fe4aaea57d"
    owner = owner.replace("if allowed and not report['provider_mappings']:", "if allowed and {line.split(None, 5)[5] for e in report['provider_mappings'] for line in e['maps']} != set(" + repr(MAPPINGS) + "):")
    owner = owner.replace("lines = [line for line in Path(f\"/proc/{item['pid']}/maps\").read_text().splitlines() if provider in line]", "required_mappings = " + repr(MAPPINGS) + "\n                        lines = [line for line in Path(f\"/proc/{item['pid']}/maps\").read_text().splitlines()\n                                 if len(line.split(None, 5)) == 6 and line.split(None, 5)[5] in required_mappings]")
    assert hashlib.sha256(owner.encode()).hexdigest() == "492c75609a0dc8ff5eaa0299f7206820bb8a42f85c3ecc5b9b0af98f2496ba9f"
    (WORK / "run_owned.py").write_text(owner)
    shutil.copyfile(OLD / "owned_guard.py", WORK / "owned_guard.py")
    prior = ROOT / "out/science-p4-ordinary-fp32-parity"
    historical = recovery / "allocation-reader-correction"
    historical.mkdir()
    preserved = []
    for filename in ("run_parity.py", "preparation.json", "actions.json", "parity-1/owner.json"):
        origin = prior / filename
        target = historical / filename
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(origin, target)
        require_same = identity(origin)
        assert identity(target)["sha256"] == require_same["sha256"]
        preserved.append({"historical_origin": key(origin), "origin_identity": require_same, "recovery": key(target), "scope": "Historical exact copy, not a live dependency or altered prior input."})
    old_owner = json.loads((historical / "parity-1/owner.json").read_text())
    assert old_owner["stop_reason"] == "diagnostic_timeout" and not old_owner["passed"]
    original_test = (ROOT / "tests/backend/kernel_parity_test.cpp").read_text()
    assert "ASSERT_GE(device.BufferAllocationBytes(), resident_before + requested);" in original_test and "ASSERT_LE(device.BufferAllocationBytes(), allocation_limit);" in original_test
    write(historical / "correction.json", {"preserved_originals": preserved,
      "source_authority": "tests/backend/kernel_parity_test.cpp:221-231; src/sirius/backend/vulkan/vulkan_device.cpp:796-845",
      "original_relation": "actual resident bytes >= resident_before+requested and <=1MiB; driver VkMemoryRequirements.size is charged.",
      "earlier_extra_equality": "Prior private reader expects completed_records*320+16 exactly. This is stronger than the original allocation assertions and is not retained in the new setup reader.",
      "first_attempt_scope": "First ordinary attempt ended at outer diagnostic_timeout before any returned ParityProbe phase/numerical verdict. No prior acceptance is promoted or rewritten.",
      "new_diagnostic_relation": "Fresh zero-active four-record setup retains original1296byte requested minimum and1MiB cap, recording actual charged bytes without exact-equality assumption.",
      "old_inputs_mutated": False})
    assets = set(dependencies) | sources | set(ARCHIVES) | {ORIGINAL_OBJECT, shader, ROOT / "bin/linux-gcc/kernels/parity_probe.spv", BASE / "current-versions/19991852c2bf8025aaa30b4cb3029e7bb980425d74de515efc4879df95eb0891"}
    assets.update(ROOT / "bin/linux-gcc" / name for name in ("compile_commands.json", "CMakeCache.txt", "build.ninja", ".ninja_deps"))
    assets.update(PRODUCER / name for name in ("build/inputs-before.json", "build/inputs-after.json", "build/owner.json", "accepted-result.json"))
    assert json.loads((PRODUCER / "build/owner.json").read_text())["passed"] and json.loads((PRODUCER / "accepted-result.json").read_text())["passed"]
    assets.update(Path(p) for p in [*MAPPINGS, "/usr/share/vulkan/icd.d/lvp_icd.json", "/usr/lib/x86_64-linux-gnu/libvulkan.so", "/usr/bin/python3.12", "/usr/bin/g++-14", "/usr/bin/ld", "/usr/bin/ar", "/usr/bin/nm", "/usr/bin/readelf", "/usr/bin/ninja", "/usr/libexec/gcc/x86_64-linux-gnu/14/cc1plus", "/usr/libexec/gcc/x86_64-linux-gnu/14/collect2", "/usr/bin/x86_64-linux-gnu-as"])
    pending = list(assets)
    while pending:
        path = pending.pop()
        with path.open("rb") as stream:
            magic = stream.read(4)
        if magic != b"\x7fELF":
            continue
        dynamic = subprocess.check_output(["/usr/bin/readelf", "-d", str(path)], stderr=subprocess.DEVNULL, text=True)
        for library in re.findall(r'\(NEEDED\).*\[([^\]]+)\]', dynamic):
            found = next((folder / library for folder in (path.resolve().parent, Path("/usr/lib/x86_64-linux-gnu"), Path("/lib/x86_64-linux-gnu")) if (folder / library).is_file()), None)
            assert found is not None, (path, library)
            if found not in assets:
                assets.add(found)
                pending.append(found)
    products = {"compile-device": [BUILD / "instrumented-vulkan_device.o", BUILD / "instrumented-vulkan_device.d"], "compile-driver": [BUILD / "setup_probe.o", BUILD / "setup_probe.d"], "link": [BUILD / "setup_probe", BUILD / "setup_probe.map"]}
    phases = ["compile-device", "compile-driver", "link", "run"]
    context = {"prepared_only": True, "product_source_revision": PRODUCT, "runner_source_revision": REVISION,
      "diagnostic_only_difference": subprocess.check_output(["git", "diff", "--name-only", PRODUCT, "HEAD"], text=True).splitlines(),
      "owned_source_root": key(WORK), "owned_product_root": key(BUILD), "phases": phases,
      "original_compile_command": command, "complete_source_bindings": source_records,
      "original_archive_membership": key(WORK / "original-symbols.json"), "required_provider_mappings": MAPPINGS,
      "products": {name: [key(p) for p in paths] for name, paths in products.items()},
      "staged_shader": key(staged), "runtime_shader": key(BUILD / "resources/kernels/parity_probe.spv"),
      "scope": "One full ordinary module and exact four-record original zero-active setup shape/public call order; diagnosis only. No active record dispatch, changed shader, optimization, numerical/scientific acceptance or qualification.",
      "public_failure_boundaries": "Original Vulkan TU bodies and error branches remain unchanged; RAII markers and wrappers only surround original calls at the original sites. No eager pipeline, deferred construction, prewarming, preload/layer, specialization or global driver-control flag.",
      "cost_scope": "Owner wall and driver whole_setup include construction/module read/load/buffers/writes/first pipeline and submit/wait/readback/destruction. Marker overhead remains included; no warm/cold or performance conclusion.",
      "call_contract": {"fixture_records": 4, "opcode": 2, "active_records": 0, "sample_stride": 16, "result_stride": 64, "buffer_bytes_in_order": [256, 1024, 16], "buffer_usage": "storage", "initial_writes_in_order": ["params16bytes", "result1024zero_bytes"], "sample_buffer_write": False, "dispatch": [1, 1, 1], "readback_bytes": 1024, "actual_record_dispatches": 0},
      "limits": {"seconds": 300, "rss_kib": 4194304, "reserve_kib": 2097152},
      "owner_reuse": {"source": key(OLD / "run_owned.py"), "guard": key(OLD / "owned_guard.py"), "change": "Same actual loader/LVP/LLVM collection extension as released parity owner; process/source/input/cleanup limits unchanged."}}
    write(WORK / "context.json", context)
    argv = {
      "compile-device": [*flags, "-MD", "-MF", str(products["compile-device"][1]), "-o", str(products["compile-device"][0]), "-c", str(instrumented)],
      "compile-driver": [*flags, "-I" + str(ROOT / "tests"), "-MD", "-MF", str(products["compile-driver"][1]), "-o", str(products["compile-driver"][0]), "-c", str(WORK / "setup_probe.cpp")],
      "link": [flags[0], "-O3", "-DNDEBUG", str(products["compile-driver"][0]), str(products["compile-device"][0]), str(ARCHIVES[0]), "/usr/lib/x86_64-linux-gnu/libvulkan.so", "-ldl", str(ARCHIVES[1]), str(ARCHIVES[2]), "-Wl,-Map," + str(products["link"][1]) + ",--cref", "-o", str(products["link"][0])],
      "run": [str(products["link"][0])],
    }
    frozen = sorted({key(p) for p in assets} | {key(p) for p in WORK.rglob("*") if p.is_file()} | {key(WORK / "preparation.json"), context["runtime_shader"]})
    actions = {}
    for index, phase in enumerate(phases):
        prior_products = [path for name in phases[:index] for path in context["products"].get(name, [])]
        prior_receipts = [key(WORK / (name + "-result.json")) for name in phases[:index]]
        actions[phase] = {"kind": "diagnostic" if phase == "run" else "build", "revision": REVISION, "expected_exit": 0, "software_provider": phase == "run", "scope": context["scope"] + " Phase:" + phase,
                          "limits": context["limits"], "argv": argv[phase], "frozen_inputs": sorted(set(frozen + prior_products + prior_receipts + ([key(WORK / "link-review.json")] if phase == "run" else [])))}
    write(WORK / "actions.json", actions)
    assets.update(p for p in WORK.rglob("*") if p.is_file())
    context["pins"] = {key(p): identity(p) for p in sorted(assets)}
    write(WORK / "preparation.json", context)
    print(json.dumps({"prepared_only": True, "pins": len(context["pins"]), "source_bindings": len(source_records), "original_strong_symbols": len(strong), "preparation": identity(WORK / "preparation.json")}, indent=2))


if __name__ == "__main__":
    main()
