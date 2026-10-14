#!/usr/bin/env python3
"""CPU-only byte/source preparation. Never run a compiler, validator or provider."""
from pathlib import Path
import hashlib
import json
import re
import shutil
import subprocess

from hint_transform import derive

ROOT = Path.cwd().resolve()
WORK = Path(__file__).resolve().parent
BUILD = ROOT / "bin/linux-gcc/p4-parity-unroll-hint-setup"
BASE = ROOT / "attestations/scientific-core/4bcd987/p4-original-setup-diagnostic"
BASE_OUT = BASE / "files/out/p4-parity-setup-diagnostic"
BASE_BIN = BASE / "files/bin/linux-gcc/p4-parity-setup-diagnostic"
RUNNER = "f5d649bf5067d1d68dc9215462f37830122f0540"
PRODUCT = "4bcd9875cd9308deabed06a5fa2f17d2321b0cbc"
VALIDATOR = ROOT / "bin/linux-gcc/qualification-tools/spirv-tools-v2026.3/build/tools/spirv-val"
LIMITS = {"seconds": 300, "rss_kib": 4194304, "reserve_kib": 2097152}


def identity(path):
    path = Path(path)
    return {"bytes": path.stat().st_size, "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "resolved_path": str(path.resolve(strict=True))}


def key(path):
    path = Path(path)
    return str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path)


def write(path, value):
    path.write_text(json.dumps(value, sort_keys=True, indent=2) + "\n")


def git(*arguments):
    return subprocess.check_output(["git", *arguments], cwd=ROOT)


def main():
    assert WORK == ROOT / "out/p4-parity-unroll-hint-setup"
    assert not (WORK / "preparation.json").exists() and not BUILD.exists()
    assert git("rev-parse", "HEAD").decode().strip() == RUNNER
    assert not git("status", "--porcelain").strip()
    evidence_commit = git("rev-parse", "4632e77^{commit}").decode().strip()
    manifest_payload = (BASE / "manifest.json").read_bytes()
    assert git("show", evidence_commit + ":" + key(BASE / "manifest.json")) == manifest_payload
    manifest = json.loads(manifest_payload)
    assets = {BASE / "manifest.json"}
    for relative, expected in manifest["files"].items():
        path = BASE / relative
        actual = identity(path)
        assert all(actual[field] == expected[field] for field in ("bytes", "sha256"))
        assert git("show", evidence_commit + ":" + key(path)) == path.read_bytes()
        assets.add(path)
    baseline = json.loads((BASE_OUT / "preparation.json").read_text())
    history = []
    for name, expected in baseline["pins"].items():
        path = BASE / "files" / name if name.startswith("out/p4-parity-setup-diagnostic/") else ROOT / name
        actual = identity(path)
        assert all(actual[field] == expected[field] for field in ("bytes", "sha256")), name
        assets.add(path)
        history.append({"original": name, "current": key(path), "bytes": actual["bytes"], "sha256": actual["sha256"]})
    for phase in ("compile-device", "compile-driver", "link", "run"):
        owner = json.loads((BASE_OUT / phase / "owner.json").read_text())
        result = json.loads((BASE_OUT / (phase + "-result.json")).read_text())
        assert result["passed"] and owner["passed"] and owner["status"] == "terminal" and owner["child_exit"] == 0
        assert owner["source_unchanged"] and owner["whole_inputs_unchanged"] and owner["stop_reason"] is None
        assert not owner["cleanup_errors"] and not owner["remaining_owned_processes"]
        assert owner["observed_births"] and all(row["absent"] for row in owner["observed_births"])
        assert identity(BASE_OUT / phase / "inputs-before.json")["sha256"] == identity(BASE_OUT / phase / "inputs-after.json")["sha256"]
        assert all(identity(BASE_OUT / phase / "owner.json")[field] == result["owner"][field] for field in ("bytes", "sha256"))
    for row in baseline["complete_source_bindings"]:
        assert identity(ROOT / row["path"])["sha256"] == row["sha256"]
    recovery = WORK / "recovery"
    recovery.mkdir()
    write(recovery / "historical-input-relocations.json", history)
    (WORK / "resources").mkdir()
    stock = (BASE_BIN / "resources/kernels/parity_probe.spv").read_bytes()
    (WORK / "resources/stock-parity_probe.spv").write_bytes(stock)
    candidate, proof = derive(stock)
    (WORK / "resources/unroll-cleared-parity_probe.spv").write_bytes(candidate)
    write(WORK / "unroll-bit-diff.json", proof)
    (BUILD / "resources/kernels").mkdir(parents=True)
    shutil.copy2(BASE_BIN / "setup_probe", BUILD / "setup_probe")
    assert identity(BUILD / "setup_probe")["sha256"] == "1349129193c71b9577cb12d2ecd52e498425adb54244340bfd9589e12b57af01"
    (BUILD / "resources/kernels/parity_probe.spv").write_bytes(candidate)
    for name in ("run_owned.py", "owned_guard.py"):
        shutil.copy2(BASE_OUT / name, WORK / name)
    local_mesa = []
    for name, lines in (("nir/nir_opt_loop_unroll.c", [896, 924]), ("spirv/vtn_structured_cfg.c", [1458, 1474])):
        origin = ROOT / "bin/toolchains/mesa-25.2.8/src/compiler" / name
        copied = recovery / "mesa-25.2.8" / name
        copied.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(origin, copied)
        local_mesa.append({"path": key(copied), "identity": identity(copied), "original": key(origin), "lines": lines})
    isolated = ROOT / "out/science-p4-ordinary-fp32-parity-isolated"
    parity_observations = []
    for ordinal, passed in ((1, True), (2, True), (3, False)):
        owner_path = isolated / f"parity-{ordinal}/owner.json"
        owner = json.loads(owner_path.read_text())
        assert owner["status"] == "terminal" and owner["passed"] is passed
        assert owner["whole_inputs_unchanged"] and owner["source_unchanged"] and not owner["cleanup_errors"]
        assert all(row["absent"] for row in owner["observed_births"])
        if not passed:
            assert owner["stop_reason"] == "diagnostic_timeout"
        for relative in (f"parity-{ordinal}/owner.json", f"parity-{ordinal}-result.json"):
            copied = recovery / "ordinary-parity-observations" / relative
            copied.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(isolated / relative, copied)
        parity_observations.append({"ordinal": ordinal, "passed": passed, "wall_seconds": owner["wall_seconds"],
                                    "stop_reason": owner["stop_reason"], "numerical_verdict_reused_for_hint_candidate": False})
    assert not any((isolated / f"parity-{i}").exists() for i in range(4, 10))
    assets.add(VALIDATOR)
    # Read ELF dynamic metadata only; never execute validator, compiler or provider.
    pending = [VALIDATOR]
    while pending:
        path = pending.pop()
        dynamic = subprocess.check_output(["/usr/bin/readelf", "-d", str(path)], text=True)
        for library in re.findall(r"\(NEEDED\).*\[([^\]]+)\]", dynamic):
            found = next((p / library for p in (path.resolve().parent, Path("/usr/lib/x86_64-linux-gnu"), Path("/lib/x86_64-linux-gnu")) if (p / library).is_file()), None)
            assert found is not None, library
            if found not in assets:
                assets.add(found)
                pending.append(found)
    baseline_result = json.loads((BASE_OUT / "run-result.json").read_text())
    context = {"schema_version": 1, "prepared_only": True, "runner_source_revision": RUNNER,
        "product_source_revision": PRODUCT,
        "diagnostic_only_difference": git("diff", "--name-only", PRODUCT, "HEAD").decode().splitlines(),
        "scope": "One hint-only full-module setup diagnosis using the unchanged archived instrumented ELF and original zero-record public call order. No numerical, performance, adoption or qualification verdict.",
        "owned_source_root": key(WORK), "owned_product_root": key(BUILD), "limits": LIMITS,
        "stock_shader": key(WORK / "resources/stock-parity_probe.spv"),
        "runtime_shader": key(BUILD / "resources/kernels/parity_probe.spv"),
        "baseline_reader": key(BASE_OUT / "operate.py"),
        "baseline": {"evidence_commit": evidence_commit, "manifest": key(BASE / "manifest.json"),
            "manifest_identity": identity(BASE / "manifest.json"), "elf": identity(BUILD / "setup_probe"),
            "all_111_archived_files_byte_joined_to_git": len(manifest["files"]) == 111,
            "all_639_historical_preparation_inputs_current_or_archived_byte_joined": len(history) == 639,
            "compile_link_owners_and_source_tool_archive_joins_reused": True,
            "setup_observation": baseline_result["setup_observation"]["summary"],
            "fresh_compile_claimed": False, "numerical_acceptance_reused": False},
        "original_public_call_contract": baseline["call_contract"],
        "public_failure_boundaries": baseline["public_failure_boundaries"],
        "cost_scope": baseline["cost_scope"], "whole_setup_and_teardown_costs_retained": True,
        "hint_change": {"proof": key(WORK / "unroll-bit-diff.json"), "changed_words": 2494,
            "only_bit_zero_cleared": True, "dont_unroll_added": False, "all_other_words_identical": True,
            "no_specialization_precision_mode_arithmetic_or_source_change": True},
        "primary_semantics": {"loop_control": "https://registry.khronos.org/SPIR-V/specs/unified1/SPIRV.html#_loop_control",
            "terms": "https://registry.khronos.org/SPIR-V/specs/unified1/SPIRV.html#_terms",
            "unroll_is_a_performance_hint": True, "hints_do_not_affect_program_semantics": True,
            "numerical_or_driver_correctness_demonstrated_by_this_preparation": False},
        "mesa_source_cost_premise": {"source_version": "25.2.8", "sources": local_mesa,
            "forced_unroll_return_precedes_ordinary_cost_and_iteration_restrictions": True,
            "actual_loaded_driver_internal_stage_observed": False,
            "ordinary_heuristics_may_still_unroll": True, "runtime_benefit_claimed": False},
        "ordinary_case_observations": parity_observations, "ordinary_cases_4_through_9_unlaunched": True,
        "no_unchanged_retry_or_cap_increase": True,
        "required_provider_mappings": baseline["required_provider_mappings"],
        "provider_policy": {"VK_DRIVER_FILES": "/usr/share/vulkan/icd.d/lvp_icd.json", "SIRIUS_VULKAN_DEVICE": "0",
            "MESA_SHADER_CACHE_DISABLE": "true", "scope": "Original explicit LVP/llvmpipe cache-off route; actual loader/LVP/LLVM dev/inode/hash joins required."},
        "phases": ["validate", "run"], "numerical_acceptance": False, "qualification_claimed": False}
    write(WORK / "context.json", context)
    assets.update(p for p in WORK.rglob("*") if p.is_file())
    assets.update(p for p in BUILD.rglob("*") if p.is_file())
    frozen = sorted({key(p) for p in assets} | {key(WORK / "preparation.json")})
    actions = {}
    for phase, argv in (("validate", [str(VALIDATOR), str(BUILD / "resources/kernels/parity_probe.spv")]), ("run", [str(BUILD / "setup_probe")])):
        actions[phase] = {"kind": "probe", "revision": RUNNER, "expected_exit": 0, "software_provider": phase == "run",
            "scope": context["scope"], "limits": LIMITS, "argv": argv, "frozen_inputs": frozen}
    write(WORK / "actions.json", actions)
    assets.add(WORK / "actions.json")
    context["pins"] = {key(p): identity(p) for p in sorted(assets)}
    write(WORK / "preparation.json", context)
    print(json.dumps({"prepared_only": True, "pins": len(context["pins"]), "candidate": proof["candidate"],
        "changed_words": len(proof["changes"]), "preparation": identity(WORK / "preparation.json")}))


if __name__ == "__main__":
    main()
