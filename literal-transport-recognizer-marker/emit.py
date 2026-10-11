"""Compiler-only screen; never imports this diagnostic into product builds."""
import hashlib
import importlib.util
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[2]
WORK = Path(__file__).resolve().parent
SOURCE = ROOT / "src/sirius/kernels"
STAGE = WORK / "kernels"


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def seal(path):
    data = path.read_bytes()
    return {"path": str(path.relative_to(ROOT)), "bytes": len(data),
            "sha256": hashlib.sha256(data).hexdigest()}


def words(program):
    return [program["instructions"], program["registers"], len(program["outputs"]),
            *program["outputs"], *program["operations"], *program["layer_offsets"]]


def read_literal(slot, variable):
    base = slot * 5
    return [f"RetainedTriple {variable};",
            *[f"{variable}.{field} = retainedScratch[{base + offset}u];"
              for offset, field in enumerate(("hi", "lo", "tail", "error"))],
            f"{variable}.valid = retainedScratch[{base + 4}u];"]


def node_code(index, operation):
    op, destination, a, b, c = operation
    lines = [f"// node {index}: {op} {destination} {a} {b} {c}", "{",
             "RetainedTriple r = RTInvalid();"]
    if op == 0:
        lines += [f"r = RTConstant(RPFromBits({a}u));"]
    elif op == 1:
        lines += [f"r = ReadStored(rowBase + {804 + 5 * (a - 4)}u);" if 4 <= a < 44
                  else f"r = ReadOriginal({a}u, row);"]
    else:
        lines += read_literal(a, "av")
        if 2 <= op <= 6:
            lines += ["RetainedTriple bv = RTInvalid();"] if op == 6 else read_literal(b, "bv")
            lines += [f"r = RTEvaluateArithmetic(av, bv, {op}u);"]
        elif op == 7:
            lines += ["r = RTNegate(av);"]
        elif op == 8:
            lines += ["RPScalar bound = RTAbsoluteUpper(av);",
                      "if (av.valid != 0 && RPScalarGreater(bound, RPZero()))",
                      "    r = RTConstant(RPPowerOfTwo(RPFloorExponent(bound) + 1));"]
        else:
            lines += read_literal(b, "bv")
            if op == 9:
                lines += ["if (av.valid != 0 && bv.valid != 0)",
                          "    r = RTConstant(RPScalarMax(RTAbsoluteUpper(av), RTAbsoluteUpper(bv)));"]
            else:
                assert op in (10, 11)
                lines += read_literal(c, "cv")
                predicate = ("RPScalarGreaterEqual(cv.hi, RPZero())" if op == 10 else
                             "(RPScalarEqual(cv.hi, RPZero()) && RPScalarEqual(cv.lo, RPZero()) && "
                             "RPScalarEqual(cv.tail, RPZero()) && RPScalarEqual(cv.error, RPZero()))")
                lines += ["if (cv.valid != 0) {", f"    bool choose = {predicate};",
                          "    r = choose ? av : bv;", "}"]
    lines += [f"retainedScratch[{destination * 5 + offset}u] = r.{field};"
              for offset, field in enumerate(("hi", "lo", "tail", "error"))]
    lines += [f"retainedScratch[{destination * 5 + 4}u] = r.valid;", "}"]
    return lines


def literal_program(name, program):
    operations = [program["operations"][i:i + 5] for i in range(0, len(program["operations"]), 5)]
    offsets = program["layer_offsets"]
    lines = []
    schedule = []
    for layer, (start, end) in enumerate(zip(offsets, offsets[1:])):
        lines += [f"void Literal{name}Layer{layer:03d}(uint lane, uint row, uint rowBase) {{",
                  "switch (lane) {"]
        for lane in range(min(16, end - start)):
            indices = list(range(start + lane, end, 16))
            schedule.append({"family": name, "layer": layer, "lane": lane, "nodes": indices})
            lines += [f"case {lane}u:"]
            for index in indices:
                lines += node_code(index, operations[index])
            lines += ["break;"]
        lines += ["}", "}"]
    lines += [f"void Literal{name}Layer(uint layer, uint lane, uint row, uint rowBase) {{",
              "switch (layer) {"]
    for layer in range(len(offsets) - 1):
        lines += [f"case {layer}u:", f"Literal{name}Layer{layer:03d}(lane, row, rowBase);", "break;"]
    lines += ["}", "}",
              f"bool Literal{name}Program(uint row, uint rowBase, uint stage, uint lane) {{",
              "if (lane == 0) retainedStatus = 0;", "GroupMemoryBarrierWithGroupSync();",
              f"for (uint layer = 0; layer < {len(offsets) - 1}u; ++layer) {{",
              f"Literal{name}Layer(layer, lane, row, rowBase);",
              "GroupMemoryBarrierWithGroupSync();",
              "if (retainedStatus != 0) return false;", "}", "switch (lane) {"]
    for lane in range(16):
        lines += [f"case {lane}u:"]
        for output in range(lane, 40, 16):
            lines += [f"if (retainedScratch[{program['outputs'][output] * 5 + 4}u] == 0)",
                      "    RejectRetainedProgram();"]
        lines += ["break;"]
    lines += ["}", "GroupMemoryBarrierWithGroupSync();",
              "if (retainedStatus != 0) return false;", "switch (lane) {"]
    for lane in range(16):
        lines += [f"case {lane}u:"]
        for output in range(lane, 40, 16):
            lines += ["{", *read_literal(program["outputs"][output], "r"),
                      f"WriteStored(rowBase + {1004 + output * 5}u + stage * 200u, r);", "}"]
        lines += ["break;"]
    lines += ["}", "AllMemoryBarrierWithGroupSync();", "return true;", "}"]
    return "\n".join(lines), schedule


def main():
    STAGE.mkdir(parents=True, exist_ok=True)
    for path in sorted(SOURCE.glob("*.slang")):
        shutil.copyfile(path, STAGE / path.name)
    shutil.copyfile(SOURCE / "portable_binary32.h", STAGE / "portable_binary32.h")
    program_module = load("literal_screen_program", SOURCE / "retained_program.py")
    programs = {"General": program_module.build_transport_program(parallel=True),
                "Schwarzschild": program_module.build_schwarzschild_transport_program(parallel=True)}
    encoded = {name: words(program) for name, program in programs.items()}
    header = ROOT / "bin/linux-gcc/src/sirius/backend/retained/retained_kernels.h"
    table = re.search(r"kTransportProgram\{\{(.*?)\}\};", header.read_text(), re.S)
    cached = [int(value, 16) if value.startswith("0x") else int(value)
              for value in re.findall(r"(0x[0-9a-fA-F]+|\d+)u?", table[1])]
    assert encoded["General"] + encoded["Schwarzschild"] == cached
    assert len(cached) == 15967
    coverage = []
    generated = []
    for name, program in programs.items():
        assert program["registers"] == 459 and len(program["outputs"]) == 40
        assert program["layer_offsets"][0] == 0
        assert program["layer_offsets"][-1] == program["instructions"]
        for start, end in zip(program["layer_offsets"], program["layer_offsets"][1:]):
            assert 0 < end - start <= 64
            destinations, inputs = set(), set()
            for index in range(start, end):
                op, destination, a, b, c = program["operations"][5 * index:5 * index + 5]
                assert 0 <= op <= 11 and 0 <= destination < 459
                assert destination not in destinations
                destinations.add(destination)
                if op == 1:
                    assert a < 45
                elif op >= 2:
                    assert a < 459
                    inputs.add(a)
                    if op not in (6, 7, 8):
                        assert b < 459
                        inputs.add(b)
                    if op in (10, 11):
                        assert c < 459
                        inputs.add(c)
            assert not destinations.intersection(inputs)
        literal, schedule = literal_program(name, program)
        coverage += schedule
        emitted_nodes = sorted(index for item in schedule for index in item["nodes"])
        assert emitted_nodes == list(range(program["instructions"]))
        values = encoded[name]
        generated += [f"static uint Literal{name}Table[{len(values)}] = {{",
                      *[", ".join(f"{value}u" for value in values[i:i + 16]) + ","
                        for i in range(0, len(values), 16)], "};", literal]
    generated += ["groupshared uint literalTransportMatch;",
                  "bool RecognizeLiteralTransport(uint programBase, bool schwarzschild, uint lane) {",
                  "if (lane == 0) literalTransportMatch = 1;",
                  "GroupMemoryBarrierWithGroupSync();",
                  "if (schwarzschild) {",
                  f"for (uint i = lane; i < {len(encoded['Schwarzschild'])}u; i += 16u)",
                  "    if (inputs[programBase + i] != LiteralSchwarzschildTable[i])",
                  "        InterlockedAnd(literalTransportMatch, 0u);", "} else {",
                  f"for (uint i = lane; i < {len(encoded['General'])}u; i += 16u)",
                  "    if (inputs[programBase + i] != LiteralGeneralTable[i])",
                  "        InterlockedAnd(literalTransportMatch, 0u);", "}",
                  "GroupMemoryBarrierWithGroupSync();", "return literalTransportMatch != 0;", "}"]
    original = (SOURCE / "retained_transport.slang").read_text()
    assert original.count("bool EvaluateProgram(") == 1
    modified = original.replace("bool EvaluateProgram(", "bool EvaluateProgramGeneric(", 1)
    marker = "RetainedTriple Rational(int numerator, int denominator) {"
    assert modified.count(marker) == 1
    modified = modified.replace(marker, "\n".join(generated) + "\n" + marker, 1)
    marker = "#if !defined(SIRIUS_RETAINED_PORTABLE) || defined(SIRIUS_RETAINED_PARALLEL_TRANSPORT)\n    // The immutable plan leaves never-assigned slots invalid across all stages."
    assert modified.count(marker) == 1
    modified = modified.replace(marker,
        "bool literalSchwarzschild = programBase != 1u + 230u * rows;\n"
        "    bool literalProgram = RecognizeLiteralTransport(programBase, literalSchwarzschild, lane);\n"
        "    if (lane == 0) outputs[rowBase + 3] = literalProgram ? 0u : 1u;\n" + marker, 1)
    marker = "if (!EvaluateProgram(row, rowBase, programBase, count, registers, layers, stage, lane)) return;"
    assert modified.count(marker) == 1
    modified = modified.replace(marker,
        "if (!EvaluateProgramGeneric(row, rowBase, programBase, count, registers, layers, stage, lane)) return;", 1)
    modified = ("#if !defined(SIRIUS_RETAINED_PORTABLE) || !defined(SIRIUS_RETAINED_NORMAL_SUM32) || "
                "!defined(SIRIUS_RETAINED_PARALLEL_TRANSPORT) || SIRIUS_RETAINED_EXECUTION_LANES != 16\n"
                "#error This compiler-only diagnostic requires the existing portable normal Transport layout.\n"
                "#endif\n" + modified)
    (STAGE / "retained_transport.slang").write_text(modified)
    (WORK / "programs.json").write_text(json.dumps(programs, indent=2) + "\n")
    (WORK / "schedule.json").write_text(json.dumps(coverage, indent=2) + "\n")
    facts = {"scope": "compiler-only; no runtime, device, performance, or product adoption",
             "revision": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
             "header": seal(header), "inputs": [seal(path) for path in sorted(SOURCE.glob("*.slang"))] +
             [seal(SOURCE / "retained_program.py"), seal(ROOT / "scripts/build-retained-kernels.py")],
             "programs": {name: {"words": len(encoded[name]), "nodes": program["instructions"],
                         "layers": len(program["layer_offsets"]) - 1, "roots": 40, "registers": 459}
                         for name, program in programs.items()},
             "constructor_words_equal_current_embedded_table": True,
             "recognition": "diagnostic isolation only: compare all selected words after original guards and flat bypass; write reserved diagnostic word3 as zero on canonical match and one on mismatch; always original generic interpreter; noncanonical raw outputs intentionally differ",
             "emission": seal(STAGE / "retained_transport.slang"),
             "schedule": seal(WORK / "schedule.json")}
    facts["inputs"].append(seal(SOURCE / "portable_binary32.h"))
    (WORK / "emission.json").write_text(json.dumps(facts, indent=2) + "\n")
    print(json.dumps(facts["programs"]))


if __name__ == "__main__":
    main()
