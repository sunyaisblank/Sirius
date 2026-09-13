"""Build and embed the bounded retained arithmetic stages and their programs."""

import argparse
import sys
sys.dont_write_bytecode = True

import importlib.util
from pathlib import Path
import re
import struct
import subprocess


def compile_shader(source, destination, compiler, assembler, disassembler, validator):
    raw = destination.with_suffix(".compiler.spv")
    assembly = destination.with_suffix(".spvasm")
    subprocess.run([compiler, str(source), "-I", str(source.parent), "-O0",
                    "-target", "spirv", "-profile", "spirv_1_5", "-entry", "ComputeMain",
                    "-stage", "compute", "-denorm-mode-fp32", "preserve", "-o", str(raw)],
                   check=True)
    subprocess.run([disassembler, str(raw), "-o", str(assembly)], check=True)
    text = assembly.read_text()
    if "OpCapability Float64" in text or "OpTypeInt 64" in text or " Fma " in text:
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
                        entry + " RoundingModeRTE 32", 1)
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
    for name in ("compiler", "assembler", "disassembler", "validator"):
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

    def array(name, values):
        lines.append(f"inline constexpr std::array<std::uint32_t, {len(values)}> {name}{{{{")
        for start in range(0, len(values), 16):
            lines.append(",".join(str(v) + "u" for v in values[start:start + 16]) + ",")
        lines.append("}};")

    for kind, build in (("Camera", module.build_camera_program),
                         ("Transport", module.build_transport_program),
                         ("Endpoint", module.build_endpoint_program),
                         ("Dense", module.build_dense_program)):
        program = build()
        prefix = [program["instructions"], program["registers"]]
        if kind != "Camera":
            prefix.append(len(program["outputs"]))
        array("k" + kind + "Program", prefix + program["outputs"] + program["operations"])
        stem = "retained_" + kind.lower()
        code = compile_shader(source / (stem + ".slang"),
                              args.output.parent / (stem + ".spv"),
                              args.compiler, args.assembler, args.disassembler, args.validator)
        array("k" + kind + "Shader", code)
        words = (512 + 4 * program["registers"] if kind == "Camera"
                 else {"Transport":2404, "Endpoint":769, "Dense":204}[kind] + 5 * program["registers"])
        lines.append(f"inline constexpr std::size_t k{kind}RowWords = {words};")
        print(kind, program["instructions"], "instructions;", program["registers"],
              "registers;", len(code) * 4, "shader bytes")
    lines.append("} // namespace sirius::backend::retained_program")
    args.output.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
