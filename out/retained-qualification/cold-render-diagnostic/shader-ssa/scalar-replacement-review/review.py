"""One offline aggregate-scalarization candidate; no compiler or device execution."""
import hashlib
import json
from pathlib import Path
import re
import subprocess

OWN = Path(__file__).resolve().parent
BASE = OWN.parent
FLAGS = ["--target-env=vulkan1.2", "--scalar-replacement=100", "--ssa-rewrite",
         "--eliminate-dead-code-aggressive", "--preserve-bindings", "--preserve-interface"]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def inventory(text):
    definitions = {m.group(1): m.group(2).split() for m in re.finditer(
        r"^\s*(%\S+)\s*=\s*(Op\S+.*)$", text, re.M)}
    locals_ = {key: value for key, value in definitions.items()
               if value[0] == "OpVariable" and value[2] == "Function"}
    kinds = {}
    aggregate_shapes = {}
    for key, value in locals_.items():
        pointee = definitions[value[1]][2]
        kind = definitions[pointee][0]
        kinds[kind] = kinds.get(kind, 0) + 1
        if kind in {"OpTypeStruct", "OpTypeArray", "OpTypeVector"}:
            shape = " ".join(definitions[pointee])
            aggregate_shapes[shape] = aggregate_shapes.get(shape, 0) + 1
    roots = {key: key for key in locals_}
    changed = True
    while changed:
        changed = False
        for key, value in definitions.items():
            if value[0] in {"OpAccessChain", "OpInBoundsAccessChain", "OpCopyObject"}:
                base = value[2]
                if base in roots and key not in roots:
                    roots[key] = roots[base]
                    changed = True
    calls = [value for value in definitions.values() if value[0] == "OpFunctionCall"]
    escaped = {roots[arg] for value in calls for arg in value[3:] if arg in roots}
    direct = {arg for value in calls for arg in value[3:] if arg in locals_}
    aggregate_members = sum(value[0] in {"OpAccessChain", "OpInBoundsAccessChain"}
                            and key in roots for key, value in definitions.items())
    return {
        "function_locals": len(locals_), "pointee_kinds": kinds,
        "aggregate_shapes": aggregate_shapes,
        "locals_passed_to_calls_directly": len(direct),
        "locals_passed_to_calls_via_any_tracked_pointer": len(escaped),
        "function_local_access_chains": aggregate_members,
        "function_calls": len(calls), "loops": text.count("OpLoopMerge "),
        "loads": text.count("OpLoad "), "stores": text.count("OpStore "),
    }


def surface(text):
    decorations = sorted(line.strip() for line in text.splitlines() if re.search(
        r"Op(?:Member)?Decorate .* (?:Binding|DescriptorSet|Offset|ArrayStride) ", line))
    variables = {}
    for storage in ["Input", "StorageBuffer", "Workgroup"]:
        variables[storage] = len(re.findall(r"OpVariable %\S+ " + storage + r"\b", text))
    return {"decorations": decorations, "global_variable_counts": variables,
            "execution_modes": sorted(re.findall(r"OpExecutionMode .*", text)),
            "entry_points": sorted(re.findall(r"OpEntryPoint .*", text))}


def constraints(text):
    assert re.findall(r"OpCapability (\S+)", text) == ["Shader"]
    assert set(re.findall(r"OpTypeInt (\d+) [01]", text)) == {"32"}
    assert "OpTypeFloat" not in text and "OpControlBarrier" not in text
    assert not re.search(r"OpExecutionMode\S* .* (?:Denorm|RoundingMode|SignedZeroInfNan)", text)
    entry = re.search(r"OpEntryPoint GLCompute (%\S+)", text)[1]
    assert f"OpExecutionMode {entry} LocalSize 1 1 1" in text


def main():
    old = json.loads((BASE / "six-stage-ssa-result.json").read_text())
    optimizer = Path("/usr/bin/spirv-opt")
    assert digest(optimizer) == old["optimizer_sha256"]
    report = {"scope": "one offline scalar-replacement candidate, no numerical or timing verdict",
              "flags": FLAGS, "optimizer_sha256": digest(optimizer),
              "review_script_sha256": digest(Path(__file__)), "stages": []}
    for previous in old["stages"]:
        stage = previous["stage"]
        source = BASE / (stage + "-o0.spv")
        baseline = BASE / (stage + "-ssa-dead.spv")
        assert digest(source) == previous["input_sha256"]
        assert digest(baseline) == previous["output_sha256"]
        destination = OWN / (stage + "-scalar-ssa.spv")
        assembly = destination.with_suffix(".spvasm")
        subprocess.run([str(optimizer), *FLAGS, str(source), "-o", str(destination)],
                       check=True, timeout=30)
        subprocess.run(["spirv-val", "--target-env", "vulkan1.2", str(destination)],
                       check=True, timeout=30)
        subprocess.run(["spirv-dis", str(destination), "-o", str(assembly)],
                       check=True, timeout=30)
        current = assembly.read_text()
        prior = (BASE / (stage + "-ssa-dead.spvasm")).read_text()
        constraints(current)
        assert surface(current) == surface(prior)
        assert current.count("OpLoopMerge ") == prior.count("OpLoopMerge ")
        report["stages"].append({"stage": stage, "input_sha256": digest(source),
            "baseline_sha256": digest(baseline), "candidate_sha256": digest(destination),
            "baseline_bytes": baseline.stat().st_size, "candidate_bytes": destination.stat().st_size,
            "baseline": inventory(prior), "candidate": inventory(current),
            "validation_surface_and_loop_checks": "PASS"})
    (OWN / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    for stage in report["stages"]:
        print(json.dumps({key: stage[key] for key in ["stage", "baseline_bytes", "candidate_bytes"]} |
                         {"baseline": stage["baseline"], "candidate": stage["candidate"]}))


if __name__ == "__main__":
    main()
