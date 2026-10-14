#!/usr/bin/env python3
"""Prepare two existing P4 Vulkan scenes, without compilation or execution."""
from pathlib import Path
import difflib
import hashlib
import json
import re
import shutil
import struct
import subprocess
import sys

sys.dont_write_bytecode = True
ROOT = Path.cwd().resolve()
WORK = Path(__file__).resolve().parent
PRODUCT = "4bcd9875cd9308deabed06a5fa2f17d2321b0cbc"
RUNNER = "f5d649bf5067d1d68dc9215462f37830122f0540"
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
    path.write_text(json.dumps(value, sort_keys=True, indent=2) + "\n")


def main():
    assert WORK == ROOT / "out/science-p4-vulkan-session-scenes"
    assert not (WORK / "preparation.json").exists(), "already sealed"
    assert subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip() == RUNNER
    assert not subprocess.check_output(["git", "status", "--porcelain"], text=True).strip()
    original = (OLD / "run_owned.py").read_text()
    assert hashlib.sha256(original.encode()).hexdigest() == "b47e4bb22e63c16244802a5d121160a63b942f5c37614d26ded319fe4aaea57d"
    patched = original.replace("lines = [line for line in Path(f\"/proc/{item['pid']}/maps\").read_text().splitlines() if provider in line]",
        "required_mappings = " + repr(MAPPINGS) + "\n                        lines = [line for line in Path(f\"/proc/{item['pid']}/maps\").read_text().splitlines()\n                                 if len(line.split(None, 5)) == 6 and line.split(None, 5)[5] in required_mappings]")
    patched = patched.replace("if allowed and not report['provider_mappings']:", "if allowed and {line.split(None, 5)[5] for e in report['provider_mappings'] for line in e['maps']} != set(" + repr(MAPPINGS) + "):")
    assert hashlib.sha256(patched.encode()).hexdigest() == "492c75609a0dc8ff5eaa0299f7206820bb8a42f85c3ecc5b9b0af98f2496ba9f"
    (WORK / "run_owned.py").write_text(patched)
    shutil.copyfile(OLD / "owned_guard.py", WORK / "owned_guard.py")
    recovery = WORK / "recovery"
    recovery.mkdir()
    (recovery / "original-run_owned.py").write_text(original)
    (WORK / "owner-provider-maps.patch").write_text("".join(difflib.unified_diff(original.splitlines(True), patched.splitlines(True), fromfile="original-run_owned.py", tofile="run_owned.py")))
    cases = [
      {"label": "scene-1", "name": "VulkanRenderSession.ProceduralVolumetricTurbulenceReachesLiveKernel",
       "outputs": ["sirius_vk_volume_baseline.exr", "sirius_vk_volume_turbulence.exr"],
       "original_frames": [{"volumetric": True, "turbulence": False, "doppler": True}, {"volumetric": True, "turbulence": True, "doppler": True}],
       "assertions": "Original baseline nonempty, represented same length, every represented channel finite, sum absolute linear-radiance difference>1e-5."},
      {"label": "scene-2", "name": "VulkanRenderSession.ThinAndVolumetricDopplerSuppressionAffectLiveEmission",
       "outputs": ["sirius_vk_doppler_thin_on.exr", "sirius_vk_doppler_thin_off.exr", "sirius_vk_doppler_volume_on.exr", "sirius_vk_doppler_volume_off.exr"],
       "original_frames": [{"volumetric": False, "doppler": True}, {"volumetric": False, "doppler": False}, {"volumetric": True, "doppler": True}, {"volumetric": True, "doppler": False}],
       "assertions": "Original thin/volume pairs nonempty and same size, both members finite, separate sum absolute linear-radiance difference>1e-5 for each pair."},
    ]
    source = (ROOT / "tests/render/vulkan_render_test.cpp").read_text()
    registered = (ROOT / "bin/linux-gcc/tests/render/sirius_render_tests[1]_tests.cmake").read_text()
    for case in cases:
        name = case["name"].split(".", 1)[1]
        body = source.split("TEST(VulkanRenderSession, " + name + ") {", 1)[1].split("\nTEST(", 1)[0]
        assert "SIRIUS_PRECISION" not in body and "40, 24" in body and "1.0e-5" in body
        assert case["name"] in (ROOT / "tests/operating_model.json").read_text()
        case["source_line"] = source[:source.index("TEST(VulkanRenderSession, " + name)].count("\n") + 1
        case["original_registration"] = next(line for line in registered.splitlines() if line.startswith("add_test(") and "[=[" + case["name"] + "]=]" in line)
    binary = ROOT / "bin/linux-gcc/src/sirius/app/sirius_render_tests"
    assert identity(binary)["sha256"] == "2489bc56bd0c1c53ab01da021119d6292d8fdb6f51037b82819d3dcdf6e87ab4"
    executable = binary.read_bytes()
    arrays = []
    retained = ROOT / "bin/linux-gcc/src/sirius/backend/retained"
    for filename in ("retained_kernels.h", "retained_owner_kernels.h"):
        for match in re.finditer(r"inline constexpr std::array<std::uint32_t, (\d+)> (\w+)\{\{(.*?)\}\};", (retained / filename).read_text(), re.S):
            words = [int(n) for n in re.findall(r"(\d+)u", match[3])]
            assert len(words) == int(match[1])
            raw = struct.pack("<" + "I" * len(words), *words)
            assert raw in executable, match[2]
            arrays.append({"header": key(retained / filename), "name": match[2], "words": len(words), "sha256": hashlib.sha256(raw).hexdigest(), "whole_bytes_present_in_current_render_elf": True})
    assert len(arrays) == 42 and len([r for r in arrays if r["name"].endswith("Program")]) == 7
    seeds = ["tests/render/vulkan_render_test.cpp", "tests/operating_model.json", "tests/render/CMakeLists.txt", "tests/CMakeLists.txt", "cmake/sirius_testing.cmake", "src/sirius/backend/CMakeLists.txt", "scripts/build-retained-kernels.py", "src/sirius/kernels/retained_program.py", "src/sirius/render/session/render_session.cpp", "src/sirius/render/session/pixel_shading.cpp", "src/sirius/render/vulkan_renderer.cpp", "src/sirius/render/render_evidence.cpp", "src/sirius/render/exr_writer.cpp", "src/sirius/render/stb_impl.cpp", "src/sirius/backend/retained_compute.cpp", "src/sirius/backend/retained_trace_executor.cpp", "src/sirius/backend/retained_integrator.cpp", "src/sirius/backend/cpu/geodesic_tracer.cpp", "src/sirius/backend/vulkan/vulkan_device.cpp", "src/sirius/base/resource_locator.cpp", "lib/tinyexr/miniz.c", "lib/tinyexr/miniz.h", "lib/tinyexr/tinyexr.h", "lib/stb/stb_image.h"]
    seeds += [key(p) for p in (ROOT / "src/sirius/kernels").glob("retained_*.slang")]
    seeds += ["src/sirius/kernels/portable_binary32.h"]
    pending = [ROOT / name for name in seeds]
    paths = set()
    while pending:
        path = pending.pop().resolve()
        if path in paths:
            continue
        paths.add(path)
        if path.suffix in (".cpp", ".h", ".c", ".slang"):
            for name in re.findall(r'^\s*#\s*include\s*"([^\"]+)"', path.read_text(), re.M):
                search = [path.parent / name, ROOT / "src" / name, ROOT / "tests" / name, ROOT / "lib/tinyexr" / name, ROOT / "lib/stb" / name]
                target = next((p.resolve() for p in search if p.is_file()), None)
                if target is not None and target.is_relative_to(ROOT) and not target.is_relative_to(ROOT / "bin"):
                    pending.append(target)
    before = json.loads((PRODUCER / "build/inputs-before.json").read_text())
    after = json.loads((PRODUCER / "build/inputs-after.json").read_text())
    sources = []
    for path in sorted(paths):
        name = key(path)
        raw = subprocess.check_output(["git", "show", PRODUCT + ":" + name])
        assert path.read_bytes() == raw and before[name] == after[name] == identity(path), name
        dest = recovery / "source" / name
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_bytes(raw)
        sources.append({"path": name, "revision": PRODUCT, "git_blob": subprocess.check_output(["git", "rev-parse", PRODUCT + ":" + name], text=True).strip(), "sha256": hashlib.sha256(raw).hexdigest(), "recovery": key(dest)})
    owner = json.loads((PRODUCER / "build/owner.json").read_text())
    assert owner["passed"] and owner["child_exit"] == 0 and owner["source_unchanged"] and owner["whole_inputs_unchanged"]
    assert json.loads((PRODUCER / "accepted-result.json").read_text())["passed"]
    command = next(r for r in json.loads((ROOT / "bin/linux-gcc/compile_commands.json").read_text()) if r["file"].endswith("/tests/render/vulkan_render_test.cpp"))
    assert all(t in command["command"] for t in ("SIRIUS_CONTRACT_MODE=2", "SIRIUS_TEST_HAS_KERNELS=1", "SIRIUS_TEST_REQUIRE_VULKAN_RUNTIME=1", "SIRIUS_HAS_RETAINED_COMPUTE=1", "-O3 -DNDEBUG", "-std=c++2c", "-Werror", "tests/support/fail_on_skip.h"))
    generated = sorted(p for p in retained.rglob("*") if p.is_file())
    protected = []
    for path in generated:
        name = key(path)
        if name in before:
            assert before[name] == after[name] == identity(path), name
            protected.append(name)
        dest = recovery / "generated" / path.relative_to(retained)
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(path, dest)
    assert len(protected) == 34
    resource_root = binary.parent / "resources"
    resources = sorted(p for p in resource_root.rglob("*") if p.is_file())
    prior_p6 = ROOT / "attestations/scientific-core/4bcd987/p6-cpu-writer-separated-controls/files/out/p6-cpu-writer-separated-controls/preparation.json"
    historical_volume = json.loads(prior_p6.read_text())["source_and_registration_join"]["runtime_resource_files"]
    for path in resources:
        pin = identity(path)
        assert {k: pin[k] for k in ("bytes", "sha256")} == historical_volume[key(path)], path
    # Reuse review is source/read-only. Older cold timeouts do not become a new
    # scene observation or a current cache-off setup-duration guarantee.
    prior_setup = ROOT / "attestations/software-vulkan/8d5b89d/retained-preparation-count"
    assert (prior_setup / "executed-inputs/bin/linux-gcc/src/sirius/backend/retained/retained_kernels.h").read_bytes() == (retained / "retained_kernels.h").read_bytes()
    prior_setup_owner = json.loads((prior_setup / "controls/renderer/owner.json").read_text())
    assert prior_setup_owner["returncode"] == 0 and prior_setup_owner["source_and_inputs_unchanged"] and prior_setup_owner["owned_birth_absent"] and prior_setup_owner["owned_group_absent"]
    review_inputs = [prior_setup / "controls/renderer/owner.json", prior_setup / "controls/renderer/stdout.log", prior_setup / "controls/renderer/tests.xml", prior_setup / "executed-inputs/bin/linux-gcc/src/sirius/backend/retained/retained_kernels.h"]
    cold = ROOT / "attestations/software-vulkan/997fe85/cold-operational-timeout"
    review_inputs += [cold / name for name in ("classification.json", "preflight.json", "runtime-identity.json", "ctest.log", "report.json")]
    assets = set(paths) | set(generated) | set(resources) | set(review_inputs) | {binary, BASE / "current-versions" / identity(binary)["sha256"], prior_p6}
    assert identity(BASE / "current-versions" / identity(binary)["sha256"])["sha256"] == identity(binary)["sha256"]
    assets.update(PRODUCER / n for n in ("accepted-result.json", "current-product-authority.json", "source-bindings.json", "build/owner.json", "build/inputs-before.json", "build/inputs-after.json"))
    assets.add(BASE / "manifest.json")
    assets.update(ROOT / "bin/linux-gcc" / n for n in ("CMakeCache.txt", "build.ninja", "compile_commands.json", "CTestTestfile.cmake", "tests/render/CTestTestfile.cmake", "tests/render/sirius_render_tests[1]_include.cmake", "tests/render/sirius_render_tests[1]_tests.cmake"))
    assets.update(Path(p) for p in MAPPINGS)
    assets.update(Path(p) for p in ("/usr/share/vulkan/icd.d/lvp_icd.json", "/usr/bin/ctest", "/usr/bin/python3.12", "/usr/bin/g++-14", "/usr/bin/gcc-14", "/usr/bin/ld", "/usr/bin/readelf", "/usr/bin/spirv-opt", "/opt/slang/bin/slangc", "/opt/slang/lib/libslang-compiler.so", "/opt/slang/lib/libslang-llvm.so"))
    assets.update(ROOT / "bin/linux-gcc/qualification-tools/spirv-tools-v2026.3/build/tools" / name for name in ("spirv-as", "spirv-dis", "spirv-val"))
    pending = list(assets)
    while pending:
        path = pending.pop()
        with path.open("rb") as stream:
            magic = stream.read(4)
        if magic != b"\x7fELF":
            continue
        dynamic = subprocess.check_output(["/usr/bin/readelf", "-d", str(path)], stderr=subprocess.DEVNULL, text=True)
        for lib in re.findall(r'\(NEEDED\).*\[([^\]]+)\]', dynamic):
            target = next((directory / lib for directory in (path.resolve().parent, Path("/opt/slang/lib"), Path("/usr/lib/x86_64-linux-gnu"), Path("/lib/x86_64-linux-gnu")) if (directory / lib).is_file()), None)
            assert target is not None, (path, lib)
            if target not in assets:
                assets.add(target)
                pending.append(target)
    provider = {"VK_DRIVER_FILES": "/usr/share/vulkan/icd.d/lvp_icd.json", "SIRIUS_VULKAN_DEVICE": "0", "MESA_SHADER_CACHE_DISABLE": "true", "precision_request": "SIRIUS_PRECISION absent, unchanged defaultFP32", "scope": "Existing explicit llvmpipe/Lavapipe software route only; cache disabled by mature owner. No Dozen/native/performance/cache-warmness/qualification claim."}
    context = {"schema_version": 1, "prepared_only": True, "product_source_revision": PRODUCT, "runner_source_revision": RUNNER,
      "diagnostic_only_difference": subprocess.check_output(["git", "diff", "--name-only", PRODUCT, "HEAD"], text=True).splitlines(),
      "cases": cases, "source_recovery": sources, "compile_command": command, "required_provider_mappings": MAPPINGS,
      "provider_policy": provider, "software_preparation_stages": ["ray_camera", "initialize", "transport", "endpoint", "dense", "dopri_phase"],
      "retained_payload": {"all_embedded_arrays": arrays, "whole_generated_tree": [key(p) for p in generated], "producer_before_after_protected": protected, "scope": "All42 full shader/program arrays occur byte-exactly in current accepted renderELF;7 complete programs, unchanged layouts/ABI/mode-selection source. No regenerated or reduced program."},
      "runtime_resource_volume": {key(p): identity(p) for p in resources},
      "source_and_producer_join": {"accepted_result": key(PRODUCER / "accepted-result.json"), "owner": key(PRODUCER / "build/owner.json"), "render": identity(binary), "scope": "Current source/ELF/resource producer reuse; current compiler/tool bytes pinned now, not retroactive historical compiler-byte attestation."},
      "original_internal_controls": {"session": "Execute calls Start and WaitForCompletion; original thread join/condition waits have no wall deadline.", "driver": "vkQueueSubmit followed by vkQueueWaitIdle; no host timeout parameter. Outer owner supplies diagnostic containment.", "renderer": "Six zero-active software preparation stages restore capacity headers; these phase times are excluded from governed physical-work feedback and included in frame wall. Film Camera is not prepared.", "physical_dispatch": "Unchanged kDispatchStopMs=1000 supplied to RetainedTraceExecutor; feedback may subdivide reducible work and refuses irreducible single-row submissions above1000ms. This is observed post-submit safety, not preemption.", "scheduling": "Original1ms coalescing wait,10ms completion polling and30s progress reporting; original adaptive dispatch target and arithmetic admission unchanged.", "outer": "300s/4GiB sampled process-tree RSS/2GiB MemAvailable reserve per separate case; diagnostic admission only, no completion forecast."},
      "arithmetic_admission": "Original RetainedArithmeticIssue and RetainedUsesPortableArithmetic on actual device; denorm/RTE32 select native versus integer portable, normal-sum requiresRTE32, no forced capability or mode override. FP32 request remains unset; no weaker precision substitution.",
      "prior_setup_review": {"same_current_embedded_header": key(prior_setup / "executed-inputs/bin/linux-gcc/src/sirius/backend/retained/retained_kernels.h"), "finite_prior_control": key(prior_setup / "controls/renderer/owner.json"), "prior_setup_wall_seconds_as_logged": 20.2377, "prior_control_wall_seconds": prior_setup_owner["elapsed_seconds"], "prior_control_peak_sampled_rss_kib": prior_setup_owner["peak_sampled_rss_kib"], "cache_policy": "Historical Mesa cache enabled; no cache-hit/warmness attribution. Not current cache-off cost evidence.", "cold_incomplete_comparison": key(cold / "classification.json"), "comparison": "997fe85 completed five preparations218.3s then reached2666 accepted physical intervals before original300s timeout; current main shader header/modes and six-phase preparation differ. No unchanged current six-phase setup-failure record established by bounded cold-owner inventory. Both scenes still share unresolved retained whole-frame costs; no completion forecast or launch recommendation follows from this source-only comparison.", "parity_distinction": "These scenes load retained embedded stage modules/programs, not the ordinary monolithic5.56MB parity_probe.spv whose first P4 case hit outer300s before any phase record."},
      "finite_scope": "Only two named original P4 full Vulkan-session witnesses with original2/4 complete40x24 frames and four volumetric samples where enabled; no fullP4/science/integrated/platform/Mandatory/release qualification.",
      "owner_reuse": {"original": key(OLD / "run_owned.py"), "guard": key(OLD / "owned_guard.py"), "private_change": "Only exact actual loader/LVP/LLVM path collection; same process/source/input/limit/cleanup rules."}}
    write(WORK / "context.json", context)
    frozen = sorted({key(p) for p in assets} | {key(p) for p in WORK.rglob("*") if p.is_file()} | {key(WORK / "preparation.json")})
    actions = {}
    for case in cases:
        actions[case["label"]] = {"kind": "test", "revision": RUNNER, "expected_exit": 0, "software_provider": True, "scope": context["finite_scope"] + " Case:" + case["name"],
          "limits": {"seconds": 300, "rss_kib": 4194304, "reserve_kib": 2097152},
          "argv": ["/usr/bin/ctest", "--test-dir", str(ROOT / "bin/linux-gcc"), "-R", "^" + re.escape(case["name"]) + "$", "--output-on-failure", "--output-junit", "{folder}/ctest.xml", "-V", "--output-log", "{folder}/ctest.log"], "frozen_inputs": frozen}
    write(WORK / "actions.json", actions)
    assets.update(p for p in WORK.rglob("*") if p.is_file())
    context["pins"] = {key(p): identity(p) for p in sorted(assets)}
    write(WORK / "preparation.json", context)
    print(json.dumps({"prepared_only": True, "cases": 2, "source_files": len(sources), "pins": len(context["pins"]), "whole_embedded_arrays": len(arrays), "preparation": identity(WORK / "preparation.json")}, indent=2))


if __name__ == "__main__":
    main()
