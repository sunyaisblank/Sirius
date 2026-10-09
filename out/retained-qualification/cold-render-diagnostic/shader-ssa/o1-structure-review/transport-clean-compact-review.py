"""One authorized Transport-only cheap cleanup/compact/SSA candidate, offline only."""
from collections import Counter
import hashlib
import importlib.util
import json
from pathlib import Path
import re
import resource
import shlex
import subprocess
import sys
import time

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
OWN = HERE / "six-stage/transport-clean-compact"
ROOT = HERE.parents[4]
FLAGS = ["--target-env=vulkan1.2", "--eliminate-local-single-block",
         "--eliminate-local-single-store", "--eliminate-dead-code-aggressive",
         "--compact-ids", "--ssa-rewrite", "--eliminate-dead-code-aggressive",
         "--preserve-bindings", "--preserve-interface"]


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def normalized_surface(text):
    # Describe declaration meaning, retaining declared names and storage classes,
    # while recursively replacing result IDs by their type/constant definitions.
    lines = [shlex.split(line) for line in text.splitlines()
             if line.strip() and not line.lstrip().startswith(";")]
    definitions = {line[0]: line[2:] for line in lines if len(line) > 2 and line[1] == "="}
    names = {line[1]: line[2] for line in lines if line[0] == "OpName"}
    entries = {line[2]: line[3] for line in lines if line[0] == "OpEntryPoint"}
    memo = {}
    active = set()

    def token(value):
        if not value.startswith("%"):
            return value
        if value in entries:
            return ["entry", entries[value]]
        if value in memo:
            return memo[value]
        assert value not in active, "recursive interface type unsupported by this finite check"
        active.add(value)
        definition = definitions[value]
        assert definition[0].startswith("OpType") or definition[0].startswith("OpConstant") or definition[0] in {"OpVariable", "OpSpecConstant", "OpSpecConstantComposite"}, definition[0]
        result = [names.get(value), *[token(part) for part in definition]]
        active.remove(value)
        memo[value] = result
        return result

    selected = {"OpDecorate", "OpMemberDecorate", "OpExecutionMode", "OpExecutionModeId", "OpEntryPoint", "OpMemoryModel", "OpCapability", "OpExtension"}
    declarations = [[token(part) for part in line] for line in lines if line[0] in selected]
    # Exact all non-Function global variable declarations, including names,
    # type shape, initializers and storage; do not rely only on descriptor counts.
    globals_ = [[names.get(key), *[token(part) for part in value]]
                for key, value in definitions.items() if value[0] == "OpVariable" and value[2] != "Function"]
    stable = lambda items: sorted(items, key=lambda item: json.dumps(item, sort_keys=True))
    return {"declarations": stable(declarations), "globals": stable(globals_),
            "global_storage_counts": dict(Counter(value[2] for value in definitions.values()
                if value[0] == "OpVariable" and value[2] != "Function"))}


def main():
    OWN.mkdir(parents=True, exist_ok=True)
    h = load("binding_helpers", HERE / "reproduce.py")
    checks = load("module_constraints", HERE.parent / "scalar-replacement-review/review.py")
    stats = load("module_structure", HERE / "six-stage-review.py")
    prior_paths = [HERE / "six-stage/report.json", HERE / "six-stage/transport-2g/report.json"]
    prior = json.loads(prior_paths[0].read_text())
    source = HERE / "six-stage/transport-o1.spv"
    expected = json.loads((ROOT / "attestations/software-vulkan/850d9e7/ray-camera-first-and-repeat/identity.json").read_text())["artifacts"]["bin/linux-gcc/src/sirius/backend/retained/retained_transport_portable.spv"]
    assert h.sha(source) == expected["sha256"] and source.stat().st_size == expected["bytes"]
    tools = [Path("/usr/bin/spirv-opt"), Path("/usr/bin/spirv-val"), Path("/usr/bin/spirv-dis")]
    before = {str(path): h.binding(path) for path in tools}
    assert all(identity == prior["tools_before"][name] for name, identity in before.items())
    source_before = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, check=True, capture_output=True, text=True).stdout.strip()
    assert source_before == prior["source_before"]
    assert not subprocess.run(["git", "status", "--porcelain"], cwd=ROOT, check=True, capture_output=True, text=True).stdout.strip()
    report = {"scope": "one offline exact Transport-only cheap cleanup/compact/SSA candidate; no numerical/GPU/performance verdict", "status": "running",
        "preserved_prior_receipts": [h.binding(path) for path in prior_paths],
        "input": h.binding(source), "input_matches_original_850": True,
        "flags": FLAGS, "limits": {"per_tool_wall_seconds": 60, "per_tool_address_space_bytes": 2147483648},
        "tools_before": before, "source_before": source_before, "commands": [], "script": h.binding(Path(__file__))}
    output = OWN / "transport-o1-clean-compact-ssa.spv"
    assembly = output.with_suffix(".spvasm")

    def limits():
        resource.setrlimit(resource.RLIMIT_AS, (2147483648, 2147483648))

    def run(name, command):
        stdout = OWN / (name + "-stdout.log")
        stderr = OWN / (name + "-stderr.log")
        entry = {"argv": [str(part) for part in command]}
        report["commands"].append(entry)
        start = time.monotonic()
        with stdout.open("wb") as out, stderr.open("wb") as err:
            try:
                process = subprocess.Popen(command, stdout=out, stderr=err, preexec_fn=limits)
                entry["pid"] = process.pid
                try:
                    entry["exit_code"] = process.wait(timeout=60)
                except subprocess.TimeoutExpired:
                    process.kill()
                    entry["exit_code"] = process.wait()
                    entry["status"] = "timeout-killed-reaped"
                    raise
                if entry["exit_code"]:
                    raise subprocess.CalledProcessError(entry["exit_code"], command)
            finally:
                entry["wall_seconds"] = time.monotonic() - start
                entry["stdout"] = h.binding(stdout)
                entry["stderr"] = h.binding(stderr)

    try:
        run("optimize", [str(tools[0]), *FLAGS, str(source), "-o", str(output)])
        run("validate", [str(tools[1]), "--target-env", "vulkan1.2", str(output)])
        run("disassemble", [str(tools[2]), str(output), "-o", str(assembly)])
        original = source.with_suffix(".spvasm").read_text()
        candidate = assembly.read_text()
        checks.constraints(candidate)
        old_surface = normalized_surface(original)
        new_surface = normalized_surface(candidate)
        assert old_surface == new_surface, "normalized interface/layout/global declarations differ"
        loops = lambda text: sorted(re.findall(r"OpLoopMerge %\S+ %\S+ (.*)", text))
        assert loops(original) == loops(candidate), "loop/control masks differ"
        report["surface_normalization"] = "Result IDs replaced recursively by named declarations, type/constant structure and entry names; all global variables and all decorations retained."
        report["surface_sha256"] = hashlib.sha256(json.dumps(old_surface, sort_keys=True).encode()).hexdigest()
        report["input_structure"] = stats.structure(source.with_suffix(".spvasm"))
        report["output_structure"] = stats.structure(assembly)
        report["output"] = h.binding(output)
        report["constraints_normalized_layout_interface_globals_loops_controls"] = "PASS"
        report["status"] = "passed"
    except Exception as error:
        report["status"] = "incomplete"
        report["error_type"] = type(error).__name__
        report["error"] = str(error)
    finally:
        assert [h.binding(path) for path in prior_paths] == report["preserved_prior_receipts"]
        assert h.binding(source) == report["input"]
        report["tools_after"] = {str(path): h.binding(path) for path in tools}
        assert report["tools_after"] == before
        report["source_after"] = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, check=True, capture_output=True, text=True).stdout.strip()
        assert report["source_after"] == source_before
        assert not subprocess.run(["git", "status", "--porcelain"], cwd=ROOT, check=True, capture_output=True, text=True).stdout.strip()
        (OWN / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return 0 if report["status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
