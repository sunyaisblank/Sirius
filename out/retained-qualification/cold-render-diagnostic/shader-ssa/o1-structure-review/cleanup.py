"""Apply only the existing product SSA/dead pipeline to exact reproduced O1."""
import importlib.util
import json
from pathlib import Path
import re
import subprocess
import sys

sys.dont_write_bytecode = True
OWN = Path(__file__).resolve().parent
BASE = OWN.parent
FLAGS = ["--target-env=vulkan1.2", "--ssa-rewrite", "--eliminate-dead-code-aggressive",
         "--preserve-bindings", "--preserve-interface"]


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    loaded = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(loaded)
    return loaded


def main():
    reproduction = module("o1_reproduction", OWN / "reproduce.py")
    inspection = module("ssa_inspection", BASE / "scalar-replacement-review/review.py")
    recorded = json.loads((OWN / "report.json").read_text())
    source = OWN / "ray-camera-o1-reproduced.spv"
    assert reproduction.sha(source) == recorded["expected_original_module_sha256"]
    optimizer = Path("/usr/bin/spirv-opt")
    expected_optimizer = json.loads((BASE / "ssa-dead-result.json").read_text())["optimizer_sha256"]
    assert reproduction.sha(optimizer) == expected_optimizer
    optimizer_before = reproduction.binding(optimizer)
    candidate = OWN / "ray-camera-o1-ssa-dead.spv"
    assembly = candidate.with_suffix(".spvasm")
    commands = [[str(optimizer), *FLAGS, str(source), "-o", str(candidate)],
        ["/usr/bin/spirv-val", "--target-env", "vulkan1.2", str(candidate)],
        ["/usr/bin/spirv-dis", str(candidate), "-o", str(assembly)]]
    for command in commands:
        subprocess.run(command, check=True, timeout=60)
    original = source.with_suffix(".spvasm").read_text()
    cleaned = assembly.read_text()
    inspection.constraints(cleaned)
    assert inspection.surface(original) == inspection.surface(cleaned)
    loops = lambda text: sorted(re.findall(r"OpLoopMerge %\S+ %\S+ (.*)", text))
    assert loops(original) == loops(cleaned)
    assert reproduction.binding(optimizer) == optimizer_before
    assert reproduction.sha(source) == recorded["expected_original_module_sha256"]
    report = {"scope": "one offline O1 plus existing SSA/dead cleanup candidate; no GPU/numerical/timing verdict",
        "reproduction_report": reproduction.binding(OWN / "report.json"),
        "review_script": reproduction.binding(Path(__file__)),
        "optimizer": optimizer_before,
        "input": reproduction.binding(source), "output": reproduction.binding(candidate),
        "flags": FLAGS, "commands": commands,
        "input_structure": reproduction.structure(source.with_suffix(".spvasm")),
        "output_structure": reproduction.structure(assembly),
        "validation_capability_layout_interface_loop_controls": "PASS"}
    (OWN / "cleanup-report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"input": report["input"], "output": report["output"],
        "input_structure": {key: report["input_structure"][key] for key in
            ["functions", "function_calls", "loops", "loop_controls", "function_locals",
             "function_local_access_chains", "loads", "stores"]},
        "output_structure": {key: report["output_structure"][key] for key in
            ["functions", "function_calls", "loops", "loop_controls", "function_locals",
             "function_local_access_chains", "loads", "stores"]},
        "validation_capability_layout_interface_loop_controls": "PASS"}, indent=2))


if __name__ == "__main__":
    main()
