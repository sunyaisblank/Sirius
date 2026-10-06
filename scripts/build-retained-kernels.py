"""Build and embed the bounded retained arithmetic stages and their programs."""

import argparse
import ast
from fractions import Fraction
import sys
sys.dont_write_bytecode = True

import importlib.util
from pathlib import Path
import re
import struct
import subprocess
import tempfile
import shutil


WORKGROUP_ROWS = 1
WORKGROUP_LANES = 64


def native_controls(text, fma=False):
    if "OpCapability Float64" in text or "OpTypeFloat 64" in text or "OpTypeInt 64" in text or " Fma " in text:
        raise ValueError("retained stage introduced wide arithmetic or contraction")
    if re.search(r"\b(?:FPFastMathMode|FPFastMathDefault|RelaxedPrecision)\b", text):
        raise ValueError("retained stage introduced relaxed arithmetic")
    operations = re.findall(r"(%\S+) = OpF(?:Add|Sub|Mul) ", text)
    fmas = re.findall(r"(%\S+) = OpFmaKHR ", text)
    decorated = set(re.findall(r"OpDecorate (%\S+) NoContraction", text))
    if (not operations and not fmas) or any(operation not in decorated for operation in operations + fmas):
        raise ValueError("retained arithmetic lost NoContraction")
    if bool(fmas) != fma:
        raise ValueError("retained FMA instructions disagree with the selected variant")
    entry = re.search(r"OpEntryPoint GLCompute (%\S+)", text)[1]
    if f"OpExecutionMode {entry} DenormPreserve 32" not in text:
        raise ValueError("retained arithmetic lost subnormal preservation")
    capabilities = "OpCapability Shader\n               OpCapability RoundingModeRTE"
    if fma:
        capabilities += "\n               OpCapability FMAKHR\n               OpCapability SignedZeroInfNanPreserve"
    text = text.replace("OpCapability Shader", capabilities, 1)
    if fma:
        end = list(re.finditer(r"^.*OpCapability .*?$", text, re.M))[-1].end()
        text = text[:end] + '\n               OpExtension "SPV_KHR_fma"' + text[end:]
    local_size = re.search(r"^.*OpExecutionMode " + re.escape(entry) + r" LocalSize.*$", text, re.M)[0]
    modes = "\n               OpExecutionMode " + entry + " RoundingModeRTE 32"
    if fma:
        modes += "\n               OpExecutionMode " + entry + " SignedZeroInfNanPreserve 32"
    return text.replace(local_size, local_size + modes, 1)


def supports_fma32(directory, compiler, assembler, disassembler, validator):
    # Probe only toolset capability. Errors in the real candidate remain fatal.
    # Ordinary GLSL Fma does not guarantee the fused, correctly rounded residual.
    with tempfile.TemporaryDirectory(prefix="fma-tool-probe-", dir=directory) as temporary:
        root = Path(temporary)
        source = root / "probe.slang"
        raw, assembly, final = (root / name for name in ("raw.spv", "probe.spvasm", "probe.spv"))
        source.write_text('''[[vk::binding(0,0)]] RWStructuredBuffer<float> values;
[shader("compute")][numthreads(1,1,1)] void ComputeMain() {
 float a=values[0],b=values[1],c=values[2];
 values[3]=spirv_asm {
  OpFmaKHR $$float %v $a $b $c;
  OpDecorate %v NoContraction;
  OpCopyObject $$float result %v;
 };
}
''')
        try:
            subprocess.run([compiler, str(source), "-O0", "-target", "spirv", "-profile", "spirv_1_5",
                "-entry", "ComputeMain", "-stage", "compute", "-denorm-mode-fp32", "preserve", "-o", str(raw)],
                capture_output=True, check=True)
            subprocess.run([disassembler, str(raw), "-o", str(assembly)], capture_output=True, check=True)
            assembly.write_text(native_controls(assembly.read_text(), fma=True))
            subprocess.run([assembler, "--target-env", "spv1.5", str(assembly), "-o", str(final)], capture_output=True, check=True)
            subprocess.run([validator, "--target-env", "vulkan1.2", str(final)], capture_output=True, check=True)
        except subprocess.CalledProcessError as error:
            print("Optional retained FMA32 unavailable in configured toolset:",
                  Path(error.cmd[0]).name, "exit", error.returncode,
                  error.stderr.decode(errors="replace").strip(), flush=True)
            return False
    return True


def validate_portable_coefficients(source):
    """Check stored words against the exact tableau, without shader arithmetic."""
    text = re.sub(r"/\*.*?\*/|//[^\n]*", "", source.read_text(), flags=re.S)
    def table(name):
        match = re.search(r"static const int " + name + r"\[7\](?:\[7\])? = (\{.*?\});", text, re.S)
        if not match:
            raise ValueError(f"missing retained tableau {name}")
        return ast.literal_eval(match[1].replace("{", "[").replace("}", "]"))

    cn, cd, an, ad, en, ed = map(table, ("cn", "cd", "an", "ad", "en", "ed"))
    if (any(type(values) is not list or len(values) != 7 or
            any(type(value) is not int for value in values) for values in (cn, cd, en, ed)) or
            any(type(matrix) is not list or len(matrix) != 7 or
                any(type(row) is not list or len(row) != 7 or
                    any(type(value) is not int for value in row) for row in matrix)
                for matrix in (an, ad))):
        raise ValueError("retained tableau requires integer 7 / 7-by-7 arrays")
    pairs = [(cn[stage], cd[stage]) for stage in range(1, 7)]
    for stage in range(2, 7):
        for previous in range(1, stage):
            if an[stage][previous] != 0:
                consumer = 6 + (stage - 2) * (stage - 1) // 2 + previous - 1
                if stage == 6:
                    consumer -= 1
                if consumer != len(pairs):
                    raise ValueError("retained correction producer/consumer slots disagree")
                pairs.append((an[stage][previous], ad[stage][previous]))
    if len(pairs) != 20 or en[1] != 0:
        raise ValueError("retained embedded-error coefficient layout changed")
    pairs += [(en[stage], ed[stage]) for stage in range(2, 7)]
    function = re.search(r"\bvoid\s+CacheCoefficient\s*\(\s*uint\s+lane\s*\)\s*\{\s*"
        r"if\s*\(\s*lane\s*>=\s*SIRIUS_RETAINED_COEFFICIENTS\s*\)\s*return\s*;\s*"
        r"#ifdef\s+SIRIUS_RETAINED_PORTABLE\s*(.*?)\s*#else\b", text, re.S)
    if not function:
        raise ValueError("retained literal-cache branch is missing")
    branch = function[1]
    base = re.match(r"\s*uint\s+base\s*=\s*lane\s*\*\s*5\s*;\s*", branch)
    if not base:
        raise ValueError("retained literal-cache base changed")
    cursor, blocks = base.end(), []
    pattern = re.compile(r"\s*if\s*\(\s*lane\s*==\s*(\d+)u\s*\)\s*\{(.*?)\}\s*", re.S)
    while cursor < len(branch):
        match = pattern.match(branch, cursor)
        if not match:
            raise ValueError("unexpected statement in retained literal-cache branch")
        blocks.append((match[1], match[2]))
        cursor = match.end()
    if len(blocks) != 25 or sorted(int(slot) for slot, _ in blocks) != list(range(25)):
        raise ValueError("retained tableau literal slots do not match the coefficient cache")

    def number(word):
        exponent = (word >> 23) & 255
        if exponent == 255:
            raise ValueError("retained tableau literal is not finite")
        significand = word & 0x7fffff
        power = -149 if exponent == 0 else exponent - 150
        if exponent:
            significand += 1 << 23
        scale = Fraction(2**power) if power >= 0 else Fraction(1, 2**-power)
        return (-1 if word >> 31 else 1) * significand * scale

    for slot, block in blocks:
        stores, cursor = [], 0
        store = re.compile(r"\s*retainedCoefficients\s*\[\s*base\s*\+\s*(\d+)u\s*\]\s*"
                           r"=\s*0x([0-9a-fA-F]{8})u\s*;\s*")
        for _ in range(5):
            match = store.match(block, cursor)
            if not match:
                raise ValueError(f"retained tableau slot {slot} store sequence changed")
            stores.append((match[1], match[2]))
            cursor = match.end()
        if not re.fullmatch(r"\s*return\s*;\s*", block[cursor:]):
            raise ValueError(f"unexpected statement in retained tableau slot {slot}")
        if len(stores) != 5 or sorted(int(index) for index, _ in stores) != list(range(5)):
            raise ValueError(f"retained tableau slot {slot} does not store five distinct words")
        words = [word for _, word in sorted((int(index), int(word, 16)) for index, word in stores)]
        hi, lo, tail, radius = map(number, words[:4])
        spacing = lambda word: number((word & 0x7fffffff) + 1) - abs(number(word))
        if (words[4] != 1 or words[3] & 0x80000000 or
                max(abs(hi), abs(lo), abs(tail), radius) > 2**120 or
                abs(lo) > spacing(words[0]) or abs(tail) > spacing(words[1]) or
                (hi == 0 and (lo != 0 or tail != 0))):
            raise ValueError(f"retained tableau slot {slot} is not represented")
        numerator, denominator = pairs[int(slot)]
        if denominator <= 0 or abs(hi + lo + tail - Fraction(numerator, denominator)) > radius:
            raise ValueError(f"retained tableau slot {slot} does not enclose its exact rational")


def compile_shader(source, destination, compiler, assembler, disassembler, validator, registers, terms, layers, prefix=0, fp64=False, portable=False, optimizer=None, fma=False):
    raw = destination.with_suffix(".compiler.spv")
    assembly = destination.with_suffix(".spvasm")
    definitions = ["-DSIRIUS_RETAINED_FP64=1"] if fp64 else []
    if fma:
        if not fp64 or portable or source.stem not in ("retained_transport", "retained_endpoint"):
            raise ValueError("FMA32 is qualified only for native-wide retained Transport and Endpoint")
        definitions.append("-DSIRIUS_RETAINED_FMA32=1")
    if portable:
        definitions.append("-DSIRIUS_RETAINED_PORTABLE=1")
        if source.stem == "retained_transport":
            validate_portable_coefficients(source)
    original_inputs = ({"retained_camera": 32, "retained_ray_camera": 45}.get(source.stem, 0)
                       if portable else 0)
    coefficients = 25 if source.stem == "retained_transport" else 0
    if coefficients:
        definitions.append(f"-DSIRIUS_RETAINED_COEFFICIENTS={coefficients}")
    if (registers * terms + original_inputs * 4 + coefficients * 5) * 4 + 8 > 16384:
        raise ValueError("retained program exceeds the portable shared-memory bound")
    definitions += [f"-DSIRIUS_RETAINED_REGISTERS={registers}",
                    f"-DSIRIUS_RETAINED_TERMS={terms}",
                    f"-DSIRIUS_RETAINED_LANES={WORKGROUP_LANES}",
                    f"-DSIRIUS_RETAINED_EXECUTION_LANES={1 if portable else WORKGROUP_LANES}",
                    f"-DSIRIUS_RETAINED_LAYERS={layers}",
                    f"-DSIRIUS_RETAINED_PREFIX={prefix}"]
    float_controls = [] if portable else ["-denorm-mode-fp32", "preserve"]
    # Bound inlining to the two portable camera stages.
    optimization = ("-O1" if portable and source.stem in ("retained_camera", "retained_ray_camera")
                    else "-O0")
    subprocess.run([compiler, str(source), *definitions, "-I", str(source.parent), optimization,
                    "-target", "spirv", "-profile", "spirv_1_5", "-entry", "ComputeMain",
                    "-stage", "compute", *float_controls, "-o", str(raw)],
                   check=True)
    if portable:
        if not optimizer:
            raise ValueError("portable retained shaders require spirv-opt")
        # Promote local temporaries before driver compilation.
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
    text = native_controls(text, fma=fma)
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
    fma_available = supports_fma32(args.output.parent, args.compiler, args.assembler,
                                  args.disassembler, args.validator)
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
        if kind in ("Transport", "Endpoint"):
            destination = args.output.parent / (stem + "_fma.spv")
            lines.append(f"inline constexpr bool k{kind}FmaAvailable = {'true' if fma_available else 'false'};")
            if fma_available:
                code = compile_shader(source / (stem + ".slang"), destination,
                    args.compiler, args.assembler, args.disassembler, args.validator,
                    program["registers"], terms, len(program["layer_offsets"])-1,
                    program.get("prefix_instructions", 0), fp64=True, fma=True)
                array("k" + kind + "FmaShader", code)
            else:
                shutil.copyfile(args.output.parent / (stem + "_fp64.spv"), destination)
                lines.append(f"inline constexpr auto& k{kind}FmaShader = k{kind}Fp64Shader;")
            print(kind + " optional FMA32:", "validated" if fma_available else "integer fallback", flush=True)
        words = ((512 if kind == "Camera" else 576) + 4 * program["registers"] if kind in ("Camera", "RayCamera")
                 else {"Transport":2404, "Endpoint":769, "Dense":204, "Initialize":204}[kind] + 5 * program["registers"])
        lines.append(f"inline constexpr std::size_t k{kind}RowWords = {words};")
        print(kind, program["instructions"], "instructions;", program["registers"],
              "registers;", sizes, "native/native-wide/portable/portable-wide shader bytes", flush=True)
    lines.append("} // namespace sirius::backend::retained_program")
    args.output.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
