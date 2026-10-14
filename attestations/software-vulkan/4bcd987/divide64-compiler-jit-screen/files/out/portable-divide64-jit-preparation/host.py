"""Narrow root-operated host build, runtime input pin and post-owner readers."""
from pathlib import Path
import hashlib
import importlib.util
import json
import re
import shlex
import subprocess
import sys

sys.dont_write_bytecode = True
WORK = Path(__file__).resolve().parent
ROOT = WORK.parents[1]
BUILD = ROOT / "bin/linux-gcc/portable-divide64-jit"
SCREEN = ROOT / "out/portable-divide64-screen-capability-order"


def record(path):
    path = Path(path)
    with path.open("rb") as stream:
        digest = hashlib.file_digest(stream, "sha256").hexdigest()
    return {"bytes": path.stat().st_size, "sha256": digest,
            "resolved_path": str(path.resolve(strict=True))}


def write(name, value):
    (WORK / name).write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value


def verify_static():
    frozen = json.loads((WORK / "source-freeze.json").read_text())
    for path, expected in frozen["pins"].items():
        assert record(path) == expected, path


def accepted(label):
    folder = WORK / label
    result = json.loads((folder / "owner.json").read_text())
    assert result["status"] == "terminal" and result["passed"] and result["child_exit"] == 0
    assert result["stop_reason"] is None and result["cleanup_errors"] == []
    assert result["remaining_owned_processes"] == [] and all(x["absent"] for x in result["observed_births"])
    assert result["source_unchanged"] and result["whole_inputs_unchanged"]
    assert (folder / "inputs-before.json").read_bytes() == (folder / "inputs-after.json").read_bytes()
    return result


def compile_output_paths(phase):
    assert phase in ("compile-device", "compile-probe")
    stem = "int64-vulkan_device" if phase == "compile-device" else "short_input_jit_probe"
    return [BUILD / (stem + suffix) for suffix in (".o", ".d")]


def verify_compile_outputs(phase):
    accepted(phase)
    seal = json.loads((WORK / phase / "output-seal.json").read_text())
    assert seal["phase"] == phase
    assert seal["argv"] == json.loads((WORK / "host-recipes.json").read_text())[phase]
    assert seal["source_freeze"] == record(WORK / "source-freeze.json")
    assert seal["operator_context"] == record(WORK / "operator-context.json")
    assert seal["inputs_before"] == record(WORK / phase / "inputs-before.json")
    expected = {str(path): record(path) for path in compile_output_paths(phase)}
    assert seal["outputs"] == expected
    return record(WORK / phase / "output-seal.json")


def pin_context():
    verify_static()
    assert not BUILD.exists() and not (WORK / "operator-context.json").exists()
    assert not subprocess.check_output(["git", "status", "--porcelain=v1"], cwd=ROOT)
    revision = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    assert (WORK / "build-actions.json").read_bytes() == (WORK / "build-actions.unpinned.json").read_bytes()
    actions = json.loads((WORK / "build-actions.json").read_text())
    for row in actions.values():
        assert row["revision"] == "PENDING_OPERATOR_CLEAN_CONTEXT"
        row["revision"] = revision
    write("build-actions.json", actions)
    write("operator-context.json", {"revision": revision, "source_freeze": record(WORK / "source-freeze.json"),
          "preparation": record(WORK / "preparation.json"), "build_actions": record(WORK / "build-actions.json")})
    print("Pinned source-only host context", revision, flush=True)


def build(phase):
    verify_static()
    context = json.loads((WORK / "operator-context.json").read_text())
    assert context["build_actions"] == record(WORK / "build-actions.json")
    recipe = json.loads((WORK / "host-recipes.json").read_text())
    assert phase in ("compile-device", "compile-probe", "link")
    if phase == "compile-device":
        assert not BUILD.exists()
        BUILD.mkdir(parents=True)
    else:
        accepted("compile-device")
        if phase == "link":
            for producer in ("compile-device", "compile-probe"):
                verify_compile_outputs(producer)
    if phase != "link":
        assert not (WORK / phase / "output-seal.json").exists()
        assert all(not path.exists() for path in compile_output_paths(phase))
    subprocess.run(recipe[phase], cwd=ROOT, check=True)
    if phase != "link":
        # The producing child seals its complete outputs before returning to
        # the still-running owner. No terminal-owner verdict is needed here.
        write(phase + "/output-seal.json", {"phase": phase, "argv": recipe[phase],
              "source_freeze": record(WORK / "source-freeze.json"),
              "operator_context": record(WORK / "operator-context.json"),
              "inputs_before": record(WORK / phase / "inputs-before.json"),
              "outputs": {str(path): record(path) for path in compile_output_paths(phase)}})


def check_build():
    verify_static()
    for phase in ("compile-device", "compile-probe", "link"):
        accepted(phase)
    output_seals = {phase: verify_compile_outputs(phase)
                    for phase in ("compile-device", "compile-probe")}
    known = json.loads((WORK / "source-freeze.json").read_text())["pins"]
    for name in ("int64-vulkan_device.d", "short_input_jit_probe.d"):
        text = (BUILD / name).read_text().replace("\\\n", "")
        for item in shlex.split(text.split(":", 1)[1]):
            path = Path(item).resolve(strict=True)
            assert str(path) in known and record(path) == known[str(path)], path
    needed = json.loads((WORK / "recovery/original-symbols.json").read_text())["required_strong_sirius_definitions"]
    strong = lambda text: {row.split()[2] for row in text.splitlines()
                           if len(row.split()) == 3 and row.split()[1] in ("T", "D", "B", "R")}
    symbols = {}
    for name in ("int64-vulkan_device.o", "short_input_jit_probe"):
        text = subprocess.check_output(["/usr/bin/nm", "-g", "--defined-only", str(BUILD / name)], text=True)
        assert set(needed) <= strong(text)
        symbols[name] = text
    link_map = (BUILD / "short_input_jit_probe.map").read_text()
    assert "libsirius_backend.a(vulkan_device.cpp.o)" not in link_map
    recipe = json.loads((WORK / "host-recipes.json").read_text())["link"]
    assert recipe.index(str(BUILD / "int64-vulkan_device.o")) < recipe.index(
        str(ROOT / "bin/linux-gcc/src/sirius/backend/libsirius_backend.a"))
    write("link-review.json", {"passed": True, "required_strong_sirius_definitions": needed,
          "all_required_definitions_in_private_object_and_executable": True,
          "original_vulkan_archive_member_extracted": False,
          "private_object_before_original_archive": True, "nm_outputs": symbols,
          "compile_output_seals": output_seals,
          "link_map": record(BUILD / "short_input_jit_probe.map"),
          "executable": record(BUILD / "short_input_jit_probe")})
    print("Host compile/link/closure review passed; no provider was executed.", flush=True)


def pin_runtime():
    verify_static()
    assert not (WORK / "runtime-context.json").exists() and not (WORK / "module-inputs").exists()
    link_review = json.loads((WORK / "link-review.json").read_text())
    assert link_review["passed"]
    assert link_review["executable"] == record(BUILD / "short_input_jit_probe")
    assert link_review["link_map"] == record(BUILD / "short_input_jit_probe.map")
    for phase in ("compile-device", "compile-probe", "link"):
        accepted(phase)
    output_seals = {phase: verify_compile_outputs(phase)
                    for phase in ("compile-device", "compile-probe")}
    assert link_review["compile_output_seals"] == output_seals
    compiler_owner = json.loads((SCREEN / "compile-endpoint-normal/owner.json").read_text())
    assert compiler_owner["status"] == "terminal" and compiler_owner["passed"]
    assert compiler_owner["child_exit"] == 0 and compiler_owner["stop_reason"] is None
    assert compiler_owner["source_unchanged"] and compiler_owner["whole_inputs_unchanged"]
    assert compiler_owner["cleanup_errors"] == [] and compiler_owner["remaining_owned_processes"] == []
    assert all(row["absent"] for row in compiler_owner["observed_births"])
    acceptance = json.loads((SCREEN / "compiler-acceptance.json").read_text())
    assert acceptance["passed"] and acceptance["compiler_owner"]["sha256"] == record(
        SCREEN / "compile-endpoint-normal/owner.json")["sha256"]
    assert acceptance["compiler_result"]["sha256"] == record(SCREEN / "compiler-result.json")["sha256"]
    compiled = json.loads((SCREEN / "compiler-result.json").read_text())
    inputs = WORK / "module-inputs"
    inputs.mkdir()
    modules = {}
    for owner in ("stock", "candidate"):
        original = Path(compiled["modules"][owner]["module"]["resolved_path"])
        assert record(original)["sha256"] == compiled["modules"][owner]["module"]["sha256"]
        target = inputs / (owner + "-endpoint-normal.spv")
        target.write_bytes(original.read_bytes())
        modules[owner] = record(target)
    context = json.loads((WORK / "operator-context.json").read_text())
    assert not subprocess.check_output(["git", "status", "--porcelain=v1"], cwd=ROOT)
    assert subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip() == context["revision"]
    mutable = json.loads((WORK / "capture-actions.unpinned.json").read_text())
    complete = [p for p in BUILD.iterdir() if p.is_file()]
    complete += [WORK / "link-review.json", *[WORK / phase / item
                  for phase in ("compile-device", "compile-probe", "link")
                  for item in ("owner.json", "inputs-before.json", "inputs-after.json")]]
    complete += [SCREEN / name for name in ("compiler-result.json", "compiler-acceptance.json",
                 "compile-endpoint-normal/owner.json", "compile-endpoint-normal/inputs-before.json",
                 "compile-endpoint-normal/inputs-after.json")]
    complete += [Path(row["resolved_path"]) for row in modules.values()]
    complete += [WORK / "runtime-context.json"]
    complete += [Path(row["resolved_path"]) for row in output_seals.values()]
    for action in mutable.values():
        assert action["revision"] == "PENDING_OPERATOR_CLEAN_CONTEXT"
        action["revision"] = context["revision"]
        action["frozen_inputs"] += [str(path) for path in complete]
    write("capture-actions.json", mutable)
    write("runtime-context.json", {"revision": context["revision"], "modules": modules,
          "host": record(BUILD / "short_input_jit_probe"), "link_review": record(WORK / "link-review.json"),
          "compiler_acceptance": record(SCREEN / "compiler-acceptance.json"),
          "compile_output_seals": output_seals,
          "capture_actions": record(WORK / "capture-actions.json"),
          "scope": "Frozen runtime inputs only; no Vulkan or GDB process was launched by pin-runtime."})
    print("Pinned separately owned stock/candidate capture inputs; no device launched.", flush=True)


def review_capture(label):
    verify_static()
    assert label in ("capture-stock-endpoint-normal", "capture-candidate-endpoint-normal")
    accepted(label)
    # The original analyzer requires a terminal passing owner: never invoke it
    # inside run_capture.py or while the enclosing capture owner is running.
    original = module("original_jit_analysis", WORK / "analyze_capture.py")
    report = original.inspect(WORK / label)
    features = module("observed_device_features", WORK / "read_features.py").inspect(WORK / label)
    target = WORK / label / "actual-code-review.json"
    assert not target.exists()
    target.write_text(json.dumps(report, indent=2) + "\n")
    write(label + "/feature-review.json", features)
    print("Terminal original IR/JIT/relocation and actual feature/device/module joins accepted.", flush=True)


def review_pair():
    verify_static()
    assert not (WORK / "pair-review.json").exists()
    layout = json.loads((WORK / "feature-layout.json").read_text())
    assert len(layout["features_fields"]) == 55 and len(set(layout["features_fields"])) == 55
    original = module("original_pair_jit_analysis", WORK / "analyze_capture.py")
    reader = module("paired_device_features", WORK / "read_features.py")
    normalized, captures = {}, {}
    for label in ("capture-stock-endpoint-normal", "capture-candidate-endpoint-normal"):
        accepted(label)
        folder = WORK / label
        code = json.loads((folder / "actual-code-review.json").read_text())
        features = json.loads((folder / "feature-review.json").read_text())
        # Both separately completed readers must still reproduce exactly from
        # the original raw captures and terminal owners.
        assert original.inspect(folder) == code
        assert reader.inspect(folder) == features
        assert set(features["features"]) == set(layout["features_fields"])
        features["features"] = {name: features["features"][name] for name in layout["features_fields"]}
        assert features["enabled_extension_count"] == len(features["enabled_extensions"])
        normalized[label] = features
        files = sorted(path for path in folder.rglob("*") if path.is_file())
        assert files and all(not path.is_symlink() for path in files)
        captures[label] = {str(path.relative_to(folder)): record(path) for path in files}
    stock = normalized["capture-stock-endpoint-normal"]
    assert stock == normalized["capture-candidate-endpoint-normal"]
    write("pair-review.json", {"passed": True, "normalized_feature_report": stock,
          "complete_55_enabled_fields_queue_extensions_order_and_feature_chain_equal": True,
          "captures": captures,
          "scope": "Two accepted original short-input owners/readers and identical observed feature policy; raw handles/addresses remain separately bound. No retained arithmetic, timing benefit or qualification verdict."})
    print("Both terminal owners/readers and complete normalized feature reports agree.", flush=True)


if __name__ == "__main__":
    command = sys.argv[1]
    if command == "build":
        build(sys.argv[2])
    elif command == "review-capture":
        review_capture(sys.argv[2])
    else:
        {"pin-context": pin_context, "check-build": check_build,
         "pin-runtime": pin_runtime, "review-pair": review_pair}[command]()
