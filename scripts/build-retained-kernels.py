"""Build and embed the bounded retained arithmetic stages and their programs."""

import argparse
import sys
sys.dont_write_bytecode = True

import importlib.util
from pathlib import Path
import re
import struct
import subprocess


WORKGROUP_ROWS = 1
WORKGROUP_LANES = 64


def compile_shader(source, destination, compiler, assembler, disassembler, validator, registers, terms, layers, prefix=0, fp64=False, portable=False, optimizer=None):
    raw = destination.with_suffix(".compiler.spv")
    assembly = destination.with_suffix(".spvasm")
    definitions = ["-DSIRIUS_RETAINED_FP64=1"] if fp64 else []
    if portable:
        definitions.append("-DSIRIUS_RETAINED_PORTABLE=1")
    if registers * terms * 4 + 8 > 16384:
        raise ValueError("retained program exceeds the portable shared-memory bound")
    definitions += [f"-DSIRIUS_RETAINED_REGISTERS={registers}",
                    f"-DSIRIUS_RETAINED_TERMS={terms}",
                    f"-DSIRIUS_RETAINED_LANES={WORKGROUP_LANES}",
                    f"-DSIRIUS_RETAINED_EXECUTION_LANES={1 if portable else WORKGROUP_LANES}",
                    f"-DSIRIUS_RETAINED_LAYERS={layers}",
                    f"-DSIRIUS_RETAINED_PREFIX={prefix}"]
    float_controls = [] if portable else ["-denorm-mode-fp32", "preserve"]
    subprocess.run([compiler, str(source), *definitions, "-I", str(source.parent), "-O0",
                    "-target", "spirv", "-profile", "spirv_1_5", "-entry", "ComputeMain",
                    "-stage", "compute", *float_controls, "-o", str(raw)],
                   check=True)
    if portable:
        if not optimizer:
            raise ValueError("portable retained shaders require spirv-opt")
        # Promote local state without the default inlining/unrolling expansion.
        optimized = destination.with_suffix(".optimized.spv")
        subprocess.run([optimizer, "--target-env=vulkan1.2", "--ssa-rewrite",
                        "--eliminate-dead-code-aggressive", "--preserve-bindings",
                        "--preserve-interface", str(raw), "-o", str(optimized)], check=True)
        optimized.replace(raw)
    subprocess.run([disassembler, str(raw), "-o", str(assembly)], check=True)
    text = assembly.read_text()
    if portable:
        # Inspect the exact module being embedded, not an adjacent cached dump.
        # Software RTE/gradual underflow must never depend on native float modes.
        capabilities = re.findall(r"OpCapability (\S+)", text)
        integer_widths = re.findall(r"OpTypeInt (\d+) [01]", text)
        if capabilities != ["Shader"] or not integer_widths or set(integer_widths) != {"32"}:
            raise ValueError("portable retained stage introduced an optional capability or non-32-bit integer")
        if "OpTypeFloat" in text or re.search(r"OpExecutionMode\S* .* (?:Denorm|RoundingMode|SignedZeroInfNan)", text):
            raise ValueError("portable retained stage depends on native floating arithmetic")
        entry = re.search(r"OpEntryPoint GLCompute (%\S+)", text)[1]
        if f"OpExecutionMode {entry} LocalSize 1 1 1" not in text:
            raise ValueError("portable retained stage lost its one-invocation-per-ray layout")
        if "OpControlBarrier" in text:
            raise ValueError("serial portable stage retained an inter-invocation barrier")
        subprocess.run([validator, "--target-env", "vulkan1.2", str(raw)], check=True)
        raw.replace(destination)
        data = destination.read_bytes()
        assembly.unlink()
        return struct.unpack("<" + str(len(data) // 4) + "I", data)
    if ("OpCapability Float64" in text) != fp64 or "OpTypeInt 64" in text or " Fma " in text:
        raise ValueError("retained stage introduced wide arithmetic or contraction")
    operations = re.findall(r"(%\S+) = OpF(?:Add|Sub|Mul) ", text)
    decorated = set(re.findall(r"OpDecorate (%\S+) NoContraction", text))
    if not operations or any(operation not in decorated for operation in operations):
        raise ValueError("retained arithmetic lost NoContraction")
    entry = re.search(r"OpEntryPoint GLCompute (%\S+)", text)[1]
    if f"OpExecutionMode {entry} DenormPreserve 32" not in text:
        raise ValueError("retained arithmetic lost subnormal preservation")
    text = text.replace("OpCapability Shader", "OpCapability Shader\n"
                        "               OpCapability RoundingModeRTE", 1)
    local_size = re.search(r"^.*OpExecutionMode " + re.escape(entry) + r" LocalSize.*$",
                           text, re.M)[0]
    text = text.replace(local_size, local_size + "\n               OpExecutionMode " +
                        entry + " RoundingModeRTE 32" +
                        ("\n               OpExecutionMode " + entry + " RoundingModeRTE 64" if fp64 else ""), 1)
    assembly.write_text(text)
    subprocess.run([assembler, "--target-env", "spv1.5", str(assembly), "-o", str(destination)],
                   check=True)
    subprocess.run([validator, "--target-env", "vulkan1.2", str(destination)], check=True)
    data = destination.read_bytes()
    raw.unlink()
    assembly.unlink()
    return struct.unpack("<" + str(len(data) // 4) + "I", data)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    for name in ("compiler", "assembler", "disassembler", "validator", "optimizer"):
        parser.add_argument("--" + name, required=True)
    args = parser.parse_args()
    source = Path(__file__).resolve().parents[1] / "src/sirius/kernels"
    spec = importlib.util.spec_from_file_location("retained_program", source / "retained_program.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    lines = ["// Generated retained programs and validated shaders; do not edit.",
             "#pragma once", "#include <array>", "#include <cstdint>",
             "namespace sirius::backend::retained_program {"]
    lines.append(f"inline constexpr std::size_t kRetainedGroupRows = {WORKGROUP_ROWS};")

    def array(name, values):
        lines.append(f"inline constexpr std::array<std::uint32_t, {len(values)}> {name}{{{{")
        for start in range(0, len(values), 16):
            lines.append(",".join(str(v) + "u" for v in values[start:start + 16]) + ",")
        lines.append("}};")

    for kind, build in (("Camera", module.build_camera_program),
                         ("Transport", module.build_transport_program),
                         ("Endpoint", module.build_endpoint_program),
                         ("Dense", module.build_dense_program),
                         ("Initialize", module.build_initialize_program),
                         ("RayCamera", module.build_ray_camera_program)):
        program = build(parallel=True)
        prefix = [program["instructions"], program["registers"]]
        if kind not in ("Camera", "RayCamera"):
            prefix.append(len(program["outputs"]))
        array("k" + kind + "Program", prefix + program["outputs"] + program["operations"] + program["layer_offsets"])
        stem = "retained_" + ("ray_camera" if kind == "RayCamera" else kind.lower())
        terms = 4 if kind in ("Camera", "RayCamera") else 5
        sizes = []
        for suffix, name, wide, portable in (("", "", False, False),
                                             ("_fp64", "Fp64", True, False),
                                             ("_portable", "Portable", False, True),
                                             ("_portable_fp64", "PortableFp64", True, True)):
            code = compile_shader(source / (stem + ".slang"),
                                  args.output.parent / (stem + suffix + ".spv"),
                                  args.compiler, args.assembler, args.disassembler, args.validator,
                                  program["registers"], terms, len(program["layer_offsets"])-1,
                                  program.get("prefix_instructions", 0), fp64=wide, portable=portable, optimizer=args.optimizer)
            array("k" + kind + name + "Shader", code)
            sizes.append(len(code) * 4)
        words = ((512 if kind == "Camera" else 576) + 4 * program["registers"] if kind in ("Camera", "RayCamera")
                 else {"Transport":2404, "Endpoint":769, "Dense":204, "Initialize":204}[kind] + 5 * program["registers"])
        lines.append(f"inline constexpr std::size_t k{kind}RowWords = {words};")
        print(kind, program["instructions"], "instructions;", program["registers"],
              "registers;", sizes, "native/native-wide/portable/portable-wide shader bytes", flush=True)
    lines.append("} // namespace sirius::backend::retained_program")
    args.output.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
