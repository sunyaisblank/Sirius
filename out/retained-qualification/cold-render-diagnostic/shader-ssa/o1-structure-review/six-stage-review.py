"""One bounded offline six-stage O1/SSA candidate, never runtime evidence."""
from collections import Counter
import difflib
import importlib.util
import json
from pathlib import Path
import re
import resource
import signal
import subprocess
import sys
import time

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
OWN = HERE / "six-stage"
ROOT = HERE.parents[4]
REVISION = "850d9e7858b9b3b589126f0441808582c62f2ba6"
FLAGS = ["--target-env=vulkan1.2", "--ssa-rewrite", "--eliminate-dead-code-aggressive",
         "--preserve-bindings", "--preserve-interface"]
TOOLS = {"compiler": "/opt/slang/bin/slangc", "assembler": "/usr/bin/spirv-as",
         "disassembler": "/usr/bin/spirv-dis", "validator": "/usr/bin/spirv-val",
         "optimizer": "/usr/bin/spirv-opt"}


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def structure(path):
    text = path.read_text()
    opcodes = Counter(re.findall(r"^\s*(?:%\S+\s*=\s*)?(Op\w+)\b", text, re.M))
    return {"bytes": path.with_suffix(".spv").stat().st_size,
        "functions": opcodes["OpFunction"], "function_calls": opcodes["OpFunctionCall"],
        "function_locals": len(re.findall(r"OpVariable %\S+ Function\b", text)),
        "access_chains": opcodes["OpAccessChain"] + opcodes["OpInBoundsAccessChain"],
        "loads": opcodes["OpLoad"], "stores": opcodes["OpStore"],
        "loops": opcodes["OpLoopMerge"],
        "loop_controls": dict(Counter(re.findall(r"OpLoopMerge %\S+ %\S+ (.*)", text)))}


def main():
    OWN.mkdir(exist_ok=True)
    reproduction = load("reproduction_helpers", HERE / "reproduce.py")
    inspection = load("surface_helpers", HERE.parent / "scalar-replacement-review/review.py")
    original_run = subprocess.run
    start = time.monotonic()
    commands = []

    def limits():
        resource.setrlimit(resource.RLIMIT_AS, (1024 * 1024 * 1024, 1024 * 1024 * 1024))

    def run(command, **kwargs):
        remaining = 600 - (time.monotonic() - start)
        if remaining <= 0:
            raise TimeoutError("600-second total offline bound reached")
        entry = {"argv": [str(part) for part in command]}
        commands.append(entry)
        try:
            result = original_run(command, timeout=min(60, remaining), preexec_fn=limits, **kwargs)
            entry["exit_code"] = result.returncode
            return result
        except subprocess.CalledProcessError as error:
            entry["exit_code"] = error.returncode
            raise
        except subprocess.TimeoutExpired:
            entry["disposition"] = "timeout"
            raise

    expected = json.loads((ROOT / "attestations/software-vulkan/850d9e7/ray-camera-first-and-repeat/identity.json").read_text())
    recorded = expected["artifacts"]
    report = {"scope": "one offline six-stage O1 reproduction and existing SSA candidate; no modes/native/GPU/numerical/timing qualification",
        "status": "preparing", "historical_revision": REVISION, "flags": FLAGS,
        "limits": {"total_wall_seconds": 600, "per_tool_wall_seconds": 60,
                   "per_tool_address_space_bytes": 1073741824},
        "commands": commands, "stages": [], "script": reproduction.binding(Path(__file__))}
    try:
        assert expected["source"]["revision"] == REVISION
        report["source_before"] = original_run(["git", "rev-parse", "HEAD"], cwd=ROOT,
            check=True, capture_output=True, text=True).stdout.strip()
        assert report["source_before"] == "9c9ac14f63fbfb8157924b89e05cea72fc5aa6e4"
        assert not original_run(["git", "status", "--porcelain"], cwd=ROOT,
            check=True, capture_output=True, text=True).stdout.strip()
        sources = {}
        for name, prior in recorded.items():
            if name.startswith("src/sirius/kernels/"):
                current = reproduction.binding(ROOT / name)
                assert current["bytes"] == prior["bytes"] and current["sha256"] == prior["sha256"], name
                sources[name] = current
        report["bound_kernel_inputs"] = sources
        tool_paths = [Path(path) for path in TOOLS.values()] + sorted(Path("/opt/slang/lib").glob("libslang*.so*"))
        report["tools_before"] = {str(path): reproduction.binding(path) for path in tool_paths}
        prior_tools = json.loads((HERE / "report.json").read_text())["current_toolchain_before"]
        for path, prior in prior_tools.items():
            assert report["tools_before"][path] == prior
        assert report["tools_before"][TOOLS["optimizer"]]["sha256"] == json.loads(
            (HERE.parent / "ssa-dead-result.json").read_text())["optimizer_sha256"]
        builder_path = OWN / "build-retained-kernels-850d9e7.py"
        builder_path.write_bytes(original_run(["git", "show", REVISION + ":scripts/build-retained-kernels.py"],
            cwd=ROOT, check=True, capture_output=True).stdout)
        assert reproduction.sha(builder_path) == recorded["scripts/build-retained-kernels.py"]["sha256"]
        report["historical_builder"] = reproduction.binding(builder_path)
        candidate_source = (ROOT / "scripts/build-retained-kernels.py").read_text()
        needle = '    subprocess.run([compiler, str(source), *definitions, "-I", str(source.parent), "-O0",\n'
        assert candidate_source.count(needle) == 1
        replacement = ('    # Inline portable helpers before the existing SSA cleanup; preserve native -O0.\n'
                       '    optimization = "-O1" if portable else "-O0"\n'
                       '    subprocess.run([compiler, str(source), *definitions, "-I", str(source.parent), optimization,\n')
        modified = candidate_source.replace(needle, replacement)
        patch = ''.join(difflib.unified_diff(candidate_source.splitlines(True), modified.splitlines(True),
            fromfile="a/scripts/build-retained-kernels.py", tofile="b/scripts/build-retained-kernels.py"))
        (OWN / "product-o1-ssa-candidate.patch").write_text(patch)
        (OWN / "build-retained-kernels-candidate.py").write_text(modified)
        original_run(["git", "apply", "--check", str(OWN / "product-o1-ssa-candidate.patch")],
            cwd=ROOT, check=True)
        report["candidate_patch"] = reproduction.binding(OWN / "product-o1-ssa-candidate.patch")
        builder = load("historical_six_stage_builder", builder_path)
        programs = load("recorded_six_stage_programs", ROOT / "src/sirius/kernels/retained_program.py")
        builder.subprocess.run = run
        for stage, build in [("camera", programs.build_camera_program),
                             ("transport", programs.build_transport_program),
                             ("endpoint", programs.build_endpoint_program),
                             ("dense", programs.build_dense_program),
                             ("initialize", programs.build_initialize_program),
                             ("ray_camera", programs.build_ray_camera_program)]:
            program = build(parallel=True)
            input_path = OWN / (stage + "-o1.spv")
            output_path = OWN / (stage + "-o1-ssa.spv")
            builder.compile_shader(ROOT / ("src/sirius/kernels/retained_" + stage + ".slang"), input_path,
                TOOLS["compiler"], TOOLS["assembler"], TOOLS["disassembler"], TOOLS["validator"],
                program["registers"], 4 if stage in ["camera", "ray_camera"] else 5,
                len(program["layer_offsets"]) - 1, program.get("prefix_instructions", 0),
                fp64=False, portable=True)
            prior = recorded["bin/linux-gcc/src/sirius/backend/retained/retained_" + stage + "_portable.spv"]
            assert reproduction.sha(input_path) == prior["sha256"] and input_path.stat().st_size == prior["bytes"], stage
            run([TOOLS["disassembler"], str(input_path), "-o", str(input_path.with_suffix(".spvasm"))], check=True)
            run([TOOLS["optimizer"], *FLAGS, str(input_path), "-o", str(output_path)], check=True)
            run([TOOLS["validator"], "--target-env", "vulkan1.2", str(output_path)], check=True)
            run([TOOLS["disassembler"], str(output_path), "-o", str(output_path.with_suffix(".spvasm"))], check=True)
            source_text = input_path.with_suffix(".spvasm").read_text()
            output_text = output_path.with_suffix(".spvasm").read_text()
            inspection.constraints(output_text)
            assert inspection.surface(source_text) == inspection.surface(output_text), stage
            loops = lambda text: sorted(re.findall(r"OpLoopMerge %\S+ %\S+ (.*)", text))
            assert loops(source_text) == loops(output_text), stage
            entry = {"stage": stage, "input": reproduction.binding(input_path),
                "output": reproduction.binding(output_path), "input_matches_original_850": True,
                "input_structure": structure(input_path.with_suffix(".spvasm")),
                "output_structure": structure(output_path.with_suffix(".spvasm")),
                "validation_capability_layout_interface_loops_controls": "PASS"}
            report["stages"].append(entry)
            (OWN / "report.json").write_text(json.dumps(report, indent=2) + "\n")
            print(json.dumps({"stage": stage, "input": entry["input_structure"],
                              "output": entry["output_structure"], "output_sha256": entry["output"]["sha256"]}), flush=True)
            del source_text, output_text, program
        report["tools_after"] = {str(path): reproduction.binding(path) for path in tool_paths}
        assert report["tools_before"] == report["tools_after"]
        for name, prior in sources.items():
            assert reproduction.binding(ROOT / name) == prior
        report["source_after"] = original_run(["git", "rev-parse", "HEAD"], cwd=ROOT,
            check=True, capture_output=True, text=True).stdout.strip()
        assert report["source_before"] == report["source_after"]
        assert not original_run(["git", "status", "--porcelain"], cwd=ROOT,
            check=True, capture_output=True, text=True).stdout.strip()
        report["status"] = "passed"
        report["historical_compiler_hash_limit"] = "Original 850 receipt lacks compiler binary hashes; current tools are bound/stable and every reproduced output matches its original recorded hash."
    except Exception as error:
        report["status"] = "incomplete"
        report["error"] = str(error)
        raise
    finally:
        (OWN / "report.json").write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
