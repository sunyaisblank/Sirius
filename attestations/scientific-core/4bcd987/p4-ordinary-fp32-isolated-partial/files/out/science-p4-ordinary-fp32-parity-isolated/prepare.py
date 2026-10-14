#!/usr/bin/env python3
"""Deterministic source/producer/resource preparation, never executes a test."""
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
PREFIX = str(WORK.relative_to(ROOT))
BASE = ROOT / "attestations/alignment/4bcd987/gcc-revision-synchronisation"
PRODUCER = BASE / "files/out/gcc-revision-synchronisation"
OLD = ROOT / "attestations/scientific-core/4bcd987/p4-cpu-radiation-current-source/files/out/p4-cpu-radiation-current-source"
P6 = ROOT / "attestations/scientific-core/4bcd987/p6-colour-kernel-parity-current-source/files/out/p6-colour-kernel-parity-current-source"
MAPPINGS = ["/usr/lib/x86_64-linux-gnu/libvulkan_lvp.so", "/usr/lib/x86_64-linux-gnu/libLLVM.so.20.1", "/usr/lib/x86_64-linux-gnu/libvulkan.so.1.3.275"]


def identity(path):
    path = Path(path)
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return {"bytes": path.stat().st_size, "sha256": digest.hexdigest(), "resolved_path": str(path.resolve(strict=True))}


def write(path, value):
    path.write_text(json.dumps(value, sort_keys=True, indent=2) + "\n")


def key(path):
    return str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path)


def main():
    assert WORK == ROOT / "out/science-p4-ordinary-fp32-parity-isolated"
    assert not (WORK / "preparation.json").exists(), "preparation already sealed"
    assert subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip() == RUNNER
    assert not subprocess.check_output(["git", "status", "--porcelain"], text=True).strip()
    original = (OLD / "run_owned.py").read_text()
    assert hashlib.sha256(original.encode()).hexdigest() == "b47e4bb22e63c16244802a5d121160a63b942f5c37614d26ded319fe4aaea57d"
    patched = original.replace("lines = [line for line in Path(f\"/proc/{item['pid']}/maps\").read_text().splitlines() if provider in line]",
        "required_mappings = " + repr(MAPPINGS) + "\n                        lines = [line for line in Path(f\"/proc/{item['pid']}/maps\").read_text().splitlines()\n                                 if len(line.split(None, 5)) == 6 and line.split(None, 5)[5] in required_mappings]")
    # Keep sampling until all three actual module paths have been observed.
    patched = patched.replace("if allowed and not report['provider_mappings']:",
        "if allowed and {line.split(None, 5)[5] for e in report['provider_mappings'] for line in e['maps']} != set(" + repr(MAPPINGS) + "):")
    assert patched != original
    (WORK / "run_owned.py").write_text(patched)
    shutil.copyfile(OLD / "owned_guard.py", WORK / "owned_guard.py")
    recovery = WORK / "recovery"
    recovery.mkdir(exist_ok=True)
    (recovery / "original-run_owned.py").write_text(original)
    (WORK / "owner-provider-maps.patch").write_text("".join(difflib.unified_diff(original.splitlines(True), patched.splitlines(True), fromfile="original-run_owned.py", tofile="run_owned.py")))
    details = [
      ("FullPageThorneDiskTemperatureMatchesIndependentCoreModel", 2, 4, "Double AccretionDiskD Page–Thorne temperature/ISCO; prograde, Schwarzschild, retrograde M=2, and truncated inner edge; original tolerances and positivity/diagnostics."),
      ("RepresentedSmallKerrSpinDoesNotAliasToSchwarzschildIsco", 2, 1, "Represented spin9e-4: original device ISCO<6 and Core ComputeIsco agreement."),
      ("ShakuraSunyaevZeroTorqueTemperatureMatchesCoreModel", 2, 4, "Core ShakuraSunyaevTemperature near/reference/far/exact edge; unchanged tolerances, exact reference/edge and ordering."),
      ("TruncatedGaussianOpacityMatchesFiniteColumnCoreClosure", 9, 3, "Core TruncatedGaussianOpacityDensity at midplane, finite height and outside3H; positive inner samples/exact zero exterior."),
      ("OpticallyThinGreyLayerAbsorptionMatchesHostAuthority", 17, 5, "Core GreyLayerAbsorbedFraction from zero to ln2 incl represented1e-8/1e-5; original tight absolute tolerance, positivity and half absorption."),
      ("ArbitraryLatitudeKerrZamoMatchesHostAuthority", 15, 4, "Core KerrZamoFrequencyTransfer including latitude/mass/spin and polar refusal; validity plus g/frequency/omega/u^t unchanged."),
      ("EquatorialKerrDiskTransferMatchesHostAuthority", 16, 2, "Core metric/covector Killing E,L and KerrDiskTransfer for two mass/spin/orientation samples; full and ZAMO g unchanged."),
      ("BlackbodyMatchesIntegratedCoreSpectrum", 3, 5, "Core spectral::BlackbodyToRgb integrated spectrum at500,3000,6000,15000,1000000K; original RGB tolerances."),
      ("DiskEmissionAppliesExactlyOneGFourthFactor", 5, 2, "Core color_modes::ObservedBolometricIntensity for receding/approaching g; original EXPECT_FLOAT_EQ."),
    ]
    source = (ROOT / "tests/backend/kernel_parity_test.cpp").read_text()
    registration = (ROOT / "bin/linux-gcc/tests/backend/sirius_backend_tests[1]_tests.cmake").read_text()
    model = (ROOT / "tests/operating_model.json").read_text()
    cases = []
    for ordinal, (name, opcode, records, reference) in enumerate(details, 1):
        full = "KernelParity." + name
        body = source.split("TEST(KernelParity, " + name + ") {", 1)[1].split("\nTEST(", 1)[0]
        assert "Fixture f = OpenProbe();" in body and 'OpenProbe("' not in body and full in model
        line = next(line for line in registration.splitlines() if line.startswith("add_test(") and "[=[" + full + "]=]" in line)
        assert "--gtest_filter=" + full in line and "--gtest_also_run_disabled_tests" in line
        cases.append({"label": f"parity-{ordinal}", "name": full, "opcode": opcode, "records": records,
                      "minimum_allocated_bytes": records * 320 + 16, "maximum_allocated_bytes": 1024 * 1024, "reference": reference,
                      "original_registration": line, "source_line": source[:source.index("TEST(KernelParity, " + name)].count("\n") + 1})
    p6 = json.loads((P6 / "context.json").read_text())
    source_paths = {row["path"] for row in p6["source_recovery"]}
    # Include the entire current CPU authority include closure rooted at the exact original TU.
    pending = [ROOT / "tests/backend/kernel_parity_test.cpp"]
    visited = set()
    while pending:
        path = pending.pop()
        if path in visited:
            continue
        visited.add(path)
        source_paths.add(key(path))
        for name in re.findall(r'^\s*#\s*include\s*"([^\"]+)"', path.read_text(), re.M):
            candidates = [path.parent / name, ROOT / "src" / name, ROOT / "tests" / name]
            target = next((p.resolve() for p in candidates if p.is_file()), None)
            if target is not None and target.is_relative_to(ROOT / "src") or target is not None and target.is_relative_to(ROOT / "tests"):
                pending.append(target)
    sources = []
    before = json.loads((PRODUCER / "build/inputs-before.json").read_text())
    after = json.loads((PRODUCER / "build/inputs-after.json").read_text())
    for name in sorted(source_paths):
        raw = subprocess.check_output(["git", "show", PRODUCT + ":" + name])
        assert (ROOT / name).read_bytes() == raw
        assert before[name] == after[name] == identity(ROOT / name), name
        dest = recovery / "source" / name
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_bytes(raw)
        sources.append({"path": name, "revision": PRODUCT, "git_blob": subprocess.check_output(["git", "rev-parse", PRODUCT + ":" + name], text=True).strip(), "sha256": hashlib.sha256(raw).hexdigest(), "recovery": key(dest)})
    owner = json.loads((PRODUCER / "build/owner.json").read_text())
    accepted = json.loads((PRODUCER / "accepted-result.json").read_text())
    assert owner["passed"] and owner["child_exit"] == 0 and owner["whole_inputs_unchanged"] and owner["source_unchanged"] and accepted["passed"]
    compiler = next(row for row in json.loads((ROOT / "bin/linux-gcc/compile_commands.json").read_text()) if row["file"].endswith("/tests/backend/kernel_parity_test.cpp"))
    assert compiler["command"] == p6["compile_command"]
    assert all(t in compiler["command"] for t in ("/usr/bin/g++-14", "SIRIUS_CONTRACT_MODE=2", "SIRIUS_TEST_HAS_KERNELS=1", "-O3 -DNDEBUG", "-std=c++2c", "-Werror", "tests/support/fail_on_skip.h"))
    shader_command = next(block for block in (ROOT / "bin/linux-gcc/build.ninja").read_text().split("\n\n") if block.startswith("build kernels/parity_probe.spv "))
    assert "-target spirv -profile spirv_1_5 -entry ComputeMain -stage compute" in shader_command
    assert "SIRIUS_REAL_" not in shader_command and "SIRIUS_FP" not in shader_command
    modes = {}
    products = []
    for filename in ("parity_probe.spv", "parity_probe_fp32comp.spv", "parity_probe_fp64.spv"):
        path = ROOT / "bin/linux-gcc/kernels" / filename
        pin = identity(path)
        staged = ROOT / "bin/linux-gcc/tests/backend/resources/kernels" / filename
        assert identity(staged)["sha256"] == pin["sha256"]
        words = struct.unpack("<" + "I" * (pin["bytes"] // 4), path.read_bytes())
        caps = []
        at = 5
        while at < len(words):
            count, op = words[at] >> 16, words[at] & 65535
            assert count > 0 and at + count <= len(words)
            if op == 17:
                caps.append(words[at + 1])
            at += count
        modes[filename] = {"sha256": pin["sha256"], "bytes": pin["bytes"], "capabilities": caps, "selected": filename == "parity_probe.spv"}
        for origin in (path, staged):
            products.append(origin)
        recovered = BASE / "current-versions" / pin["sha256"]
        assert recovered.is_file() and identity(recovered)["sha256"] == pin["sha256"]
        products.append(recovered)
    assert modes["parity_probe.spv"]["sha256"] == "19991852c2bf8025aaa30b4cb3029e7bb980425d74de515efc4879df95eb0891"
    assert modes["parity_probe.spv"]["capabilities"] == [1] and len({v["sha256"] for v in modes.values()}) == 3
    binary = ROOT / "bin/linux-gcc/tests/backend/sirius_backend_tests"
    assert identity(binary)["sha256"] == "e5f3b0f120ec6a12490e3e9b6e491d19eb12bfad62c1cdce11722c091a152149"
    products += [binary, BASE / "current-versions" / identity(binary)["sha256"]]
    assets = set(products)
    assets.update(ROOT / name for name in source_paths)
    assets.update(PRODUCER / name for name in ("accepted-result.json", "current-product-authority.json", "source-bindings.json", "build/inputs-before.json", "build/inputs-after.json", "build/owner.json"))
    assets.add(BASE / "manifest.json")
    assets.update(ROOT / "bin/linux-gcc" / name for name in ("CMakeCache.txt", "build.ninja", "compile_commands.json", "tests/backend/CTestTestfile.cmake", "tests/backend/sirius_backend_tests[1]_tests.cmake", "tests/backend/sirius_backend_tests[1]_include.cmake", "CTestTestfile.cmake"))
    assets.update(Path(p) for p in MAPPINGS)
    assets.update(Path(p) for p in ("/usr/share/vulkan/icd.d/lvp_icd.json", "/usr/bin/ctest", "/usr/bin/python3.12", "/usr/bin/g++-14", "/usr/bin/gcc-14", "/usr/bin/ld", "/usr/bin/readelf", "/opt/slang/bin/slangc", "/opt/slang/lib/libslang-compiler.so", "/opt/slang/lib/libslang-llvm.so"))
    # Read ELF NEEDED entries only; never use ldd or invoke a compiler/provider.
    pending = list(assets)
    while pending:
        path = pending.pop()
        if path.read_bytes()[:4] != b"\x7fELF":
            continue
        result = subprocess.check_output(["/usr/bin/readelf", "-d", str(path)], stderr=subprocess.DEVNULL, text=True)
        for lib in re.findall(r'\(NEEDED\).*\[([^\]]+)\]', result):
            search = [path.resolve().parent, Path("/opt/slang/lib"), Path("/usr/lib/x86_64-linux-gnu"), Path("/lib/x86_64-linux-gnu")]
            target = next((directory / lib for directory in search if (directory / lib).is_file()), None)
            assert target is not None, (path, lib)
            if target not in assets:
                assets.add(target)
                pending.append(target)
    prior_lvp = json.loads((ROOT / "attestations/software-vulkan/4bcd987/guarded-product-high-controls/files/out/guarded-product-high-controls/control-1/inputs-before.json").read_text())
    for name in (MAPPINGS[0], MAPPINGS[1], "/usr/share/vulkan/icd.d/lvp_icd.json"):
        assert prior_lvp[name] == identity(name)
    provider = {"VK_DRIVER_FILES": "/usr/share/vulkan/icd.d/lvp_icd.json", "SIRIUS_VULKAN_DEVICE": "0", "MESA_SHADER_CACHE_DISABLE": "true",
                "scope": "Explicit existing Mesa llvmpipe/Lavapipe route only; actual loader/LVP/LLVM mappings required. No Dozen, native platform, performance or cache warmness claim.",
                "cache_policy": "Original mature software-provider policy disables Mesa shader cache; no new cache output or cold-performance gate."}
    context = {"schema_version": 1, "prepared_only": True, "product_source_revision": PRODUCT, "runner_source_revision": RUNNER,
      "diagnostic_only_difference": subprocess.check_output(["git", "diff", "--name-only", PRODUCT, "HEAD"], text=True).splitlines(),
      "cases": cases, "source_recovery": sources, "compile_command": compiler, "shader_compile_rule": shader_command,
      "spirv_modes": modes, "required_provider_mappings": MAPPINGS, "provider_policy": provider,
      "producer": {"accepted": key(PRODUCER / "accepted-result.json"), "build_owner": key(PRODUCER / "build/owner.json"), "scope": "Previously accepted current backend/source/resource producer; current compiler/tool bytes pinned now, not retroactive historical compiler-byte attestation."},
      "batch_decision": "Nine separate one-case300s owners; no measured same-module llvmpipe total-cost bound supports a single300s batch. No forecast of completion.",
      "finite_scope": "Only nine named original P4 ordinary-FP32 numerical witnesses; no full P4/scientific/integrated/platform/Mandatory/release qualification.",
      "excluded_full_sessions": ["VulkanRenderSession.ProceduralVolumetricTurbulenceReachesLiveKernel", "VulkanRenderSession.ThinAndVolumetricDopplerSuppressionAffectLiveEmission"],
      "original_owner": {"source": key(OLD / "run_owned.py"), "guard": key(OLD / "owned_guard.py"), "private_change": "Collect exact loader/LVP/LLVM paths until observed; ownership, source/input/limits/cleanup acceptance unchanged."}}
    premise_path = WORK / "recovery/scheduling-premise.json"
    premise = json.loads(premise_path.read_text())
    assert premise["baseline"]["all_four_compile_link_run_owners_terminal_pass"]
    assert premise["baseline"]["numerical_acceptance"] is False
    assert premise["scheduling"]["numerical_evidence_reused"] is False
    assert premise["scheduling"]["causal_speedup_claimed"] is False
    assert all(case["records"] <= 64 for case in cases)
    for row in premise["historical_relocations"]:
        current = identity(ROOT / row["recovery"])
        assert all(current[field] == row[field] for field in ("bytes", "sha256"))
    context["scheduling_premise"] = {"path": key(premise_path), "identity": identity(premise_path),
        "numerical_evidence_reused": False, "causal_speedup_claimed": False,
        "root_runs_first_case_only_then_stops_on_incomplete": True,
        "original_300s_rss_reserve_provider_cache_dispatch_gates_unchanged": True}
    context["allocation_contract"] = {"source": "tests/backend/kernel_parity_test.cpp:201-224",
        "reader": "requested minimum <= actual allocation <= 1MiB; original allocation padding is permitted",
        "historical_exact_equality_reader_unchanged_in_recovery": True}
    context["batch_decision"] = "Nine separate one-case300s owners. Exact complete zero-active setup passed under the same guard with isolated scheduling; no original numerical case bound or completion forecast follows. Root runs first only, serially, and stops on incomplete."
    write(WORK / "context.json", context)
    actions = {}
    frozen = sorted({key(p) for p in assets} | {key(p) for p in WORK.rglob("*") if p.is_file()} | {PREFIX + "/preparation.json"})
    for case in cases:
        actions[case["label"]] = {"kind": "test", "revision": RUNNER, "expected_exit": 0, "software_provider": True,
          "scope": context["finite_scope"] + " Case:" + case["name"], "limits": {"seconds": 300, "rss_kib": 4194304, "reserve_kib": 2097152},
          "argv": ["/usr/bin/ctest", "--test-dir", str(ROOT / "bin/linux-gcc"), "-R", "^" + re.escape(case["name"]) + "$", "--output-on-failure", "--output-junit", "{folder}/ctest.xml", "-V", "--output-log", "{folder}/ctest.log"], "frozen_inputs": frozen}
    write(WORK / "actions.json", actions)
    assets.update(p for p in WORK.rglob("*") if p.is_file())
    context["pins"] = {key(p): identity(p) for p in sorted(assets)}
    write(WORK / "preparation.json", context)
    print(json.dumps({"prepared_only": True, "cases": len(cases), "frozen_pins": len(context["pins"]), "preparation": identity(WORK / "preparation.json")}, indent=2))


if __name__ == "__main__":
    main()
