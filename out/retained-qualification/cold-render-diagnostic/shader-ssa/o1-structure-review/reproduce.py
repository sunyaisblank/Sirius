"""Offline reproduction, distinct from original 850 numerical/runtime evidence."""
from collections import Counter
import hashlib
import importlib.util
import json
from pathlib import Path
import re
import subprocess
import sys

sys.dont_write_bytecode = True
OWN = Path(__file__).resolve().parent
BASE = OWN.parent
ROOT = OWN.parents[4]
REVISION = "850d9e7858b9b3b589126f0441808582c62f2ba6"
EXPECTED = "cdcbbe033fb74318b4486277e8827197b2bb67fefab8b2342c92190978e111bb"


def sha(path):
    result = hashlib.sha256()
    with path.open("rb") as stream:
        while chunk := stream.read(1024 * 1024):
            result.update(chunk)
    return result.hexdigest()


def binding(path):
    path = path.resolve()
    return {"path": str(path), "bytes": path.stat().st_size, "sha256": sha(path)}


def structure(path):
    text = path.read_text()
    opcodes = Counter(re.findall(r"^\s*(?:%\S+\s*=\s*)?(Op\w+)\b", text, re.M))
    pointer_types = {m.group(1): (m.group(2), m.group(3)) for m in re.finditer(
        r"(%\S+) = OpTypePointer (\S+) (%\S+)", text)}
    locals_ = {m.group(1) for m in re.finditer(
        r"(%\S+) = OpVariable (%\S+) Function\b", text)}
    defs = {m.group(1): m.group(2).split() for m in re.finditer(
        r"^\s*(%\S+)\s*=\s*(Op\S+.*)$", text, re.M)}
    roots = {key: key for key in locals_}
    changed = True
    while changed:
        changed = False
        for key, value in defs.items():
            if value[0] in {"OpAccessChain", "OpInBoundsAccessChain", "OpCopyObject"}:
                if value[2] in roots and key not in roots:
                    roots[key] = roots[value[2]]
                    changed = True
    local_chains = sum(value[0] in {"OpAccessChain", "OpInBoundsAccessChain"}
                       and key in roots for key, value in defs.items())
    functions = [line for line in text.splitlines() if re.search(r"= OpFunction ", line)]
    return {"functions": opcodes["OpFunction"], "function_calls": opcodes["OpFunctionCall"],
            "function_definitions": functions, "loops": opcodes["OpLoopMerge"],
            "loop_controls": dict(Counter(re.findall(r"OpLoopMerge %\S+ %\S+ (.*)", text))),
            "function_locals": len(locals_), "function_local_access_chains": local_chains,
            "loads": opcodes["OpLoad"], "stores": opcodes["OpStore"],
            "opcode_histogram": dict(sorted(opcodes.items()))}


def main():
    identity_path = ROOT / "attestations/software-vulkan/850d9e7/ray-camera-first-and-repeat/identity.json"
    identity = json.loads(identity_path.read_text())
    assert identity["source"]["revision"] == REVISION
    recorded = identity["artifacts"]
    sources = {}
    for name, previous in recorded.items():
        if name.startswith("src/sirius/kernels/"):
            current = binding(ROOT / name)
            assert current["bytes"] == previous["bytes"] and current["sha256"] == previous["sha256"], name
            sources[name] = current
    builder_path = OWN / "build-retained-kernels-850d9e7.py"
    builder_bytes = subprocess.run(["git", "show", REVISION + ":scripts/build-retained-kernels.py"],
        cwd=ROOT, check=True, capture_output=True, timeout=10).stdout
    builder_path.write_bytes(builder_bytes)
    assert sha(builder_path) == recorded["scripts/build-retained-kernels.py"]["sha256"]
    compiler = Path("/opt/slang/bin/slangc")
    tools = {"compiler": compiler, "assembler": Path("/usr/bin/spirv-as"),
             "disassembler": Path("/usr/bin/spirv-dis"), "validator": Path("/usr/bin/spirv-val")}
    libraries = sorted(Path("/opt/slang/lib").glob("libslang*.so*"))
    before = {str(path): binding(path) for path in [*tools.values(), *libraries]}
    version_result = subprocess.run([str(compiler), "-version"], check=True,
                                   capture_output=True, text=True, timeout=10)
    version = (version_result.stdout + version_result.stderr).strip()
    assert version == "2026.12.0.1"
    spec = importlib.util.spec_from_file_location("historical_retained_builder", builder_path)
    builder = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(builder)
    spec = importlib.util.spec_from_file_location("recorded_retained_program", ROOT / "src/sirius/kernels/retained_program.py")
    program_module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(program_module)
    program = program_module.build_ray_camera_program(parallel=True)
    assert (program["instructions"], program["registers"]) == (6115, 660)
    commands = []
    original_run = subprocess.run

    def limited_run(command, **kwargs):
        commands.append([str(part) for part in command])
        return original_run(command, timeout=60, **kwargs)

    builder.subprocess.run = limited_run
    output = OWN / "ray-camera-o1-reproduced.spv"
    builder.compile_shader(ROOT / "src/sirius/kernels/retained_ray_camera.slang", output,
        *[str(tools[key]) for key in ["compiler", "assembler", "disassembler", "validator"]],
        program["registers"], 4, len(program["layer_offsets"]) - 1,
        program.get("prefix_instructions", 0), fp64=False, portable=True)
    assert sha(output) == EXPECTED, "Reproduction differs: structural attribution is forbidden"
    assert output.stat().st_size == 1091920
    after = {str(path): binding(path) for path in [*tools.values(), *libraries]}
    assert before == after
    for name, previous in sources.items():
        assert binding(ROOT / name) == previous
    assembly = output.with_suffix(".spvasm")
    original_run([str(tools["disassembler"]), str(output), "-o", str(assembly)], check=True, timeout=30)
    report = {"scope": "offline exact-byte O1 reproduction and structural comparison; not original runtime evidence",
        "historical_revision": REVISION, "original_identity": binding(identity_path),
        "builder": binding(builder_path), "bound_kernel_inputs": sources,
        "current_toolchain_before": before, "current_toolchain_after": after,
        "compiler_version": version,
        "compiler_version_raw": {"stdout": version_result.stdout, "stderr": version_result.stderr},
        "historical_compiler_binary_hash_available": False,
        "compiler_identity_limit": "Original 850 receipt did not record Slang binary/library hashes; current full toolchain is bound and stable, and reproduced module is byte-identical to the original recorded output.",
        "commands": commands, "reproduced_module": binding(output),
        "expected_original_module_sha256": EXPECTED, "exact_output_match": True,
        "comparison": {"o0": structure(BASE / "ray_camera-o0.spvasm"),
                       "ssa": structure(BASE / "ray_camera-ssa-dead.spvasm"),
                       "reproduced_o1": structure(assembly)},
        "review_script": binding(Path(__file__))}
    (OWN / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"exact_output_match": True, "module_sha256": sha(output),
        "compiler_version": version, "bound_kernel_inputs": len(sources),
        "comparison": {key: {name: value[name] for name in ["functions", "function_calls", "loops",
          "loop_controls", "function_locals", "function_local_access_chains", "loads", "stores"]}
          for key, value in report["comparison"].items()}}, indent=2))


if __name__ == "__main__":
    main()
