"""Emit the fixed arithmetic program and independently computed camera witnesses."""

import sys
sys.dont_write_bytecode = True

import argparse
import json
from pathlib import Path
import struct

from program import build_program


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    program = build_program()
    source = Path(__file__).resolve().parent
    references = json.loads((source / "reference_cases.json").read_text())
    assert references["schema"] == "sirius-retained-camera-reference-cases-v1"
    lines = ["// Generated from program.py and reference_cases.json; do not edit.",
             "#pragma once", "#include <array>", "#include <cstdint>",
             "namespace sirius::test::retained_camera {",
             "struct Case { const char* name; std::array<std::uint32_t, 32> input; "
             "std::array<long double, 104> reference, reference_gap; };",
             f"inline constexpr std::size_t kInstructions = {program['instructions']};",
             f"inline constexpr std::size_t kRegisters = {program['registers']};",
             f"inline constexpr std::size_t kOutputWords = {program['output_words']};"]
    code = [program["instructions"], program["registers"],
            *program["outputs"], *program["operations"]]
    lines.append(f"inline constexpr std::array<std::uint32_t, {len(code)}> kProgram{{{{")
    for offset in range(0, len(code), 20):
        lines.append(",".join(str(value) + "u" for value in code[offset:offset + 20]) + ",")
    lines.append("}};")
    lines.append(f"inline constexpr std::array<Case, {len(references['cases'])}> kCases{{{{")
    for case in references["cases"]:
        assert len(case["input"]) == 32
        assert len(case["scientific"]) == len(case["reference_gap"]) == 104
        inputs = struct.unpack("<32I", struct.pack("<32f", *case["input"]))
        lines.append("{" + json.dumps(case["name"]) + ", {{" +
                     ",".join(str(value) + "u" for value in inputs) + "}}, {{")
        lines.append(",\n".join(value + "L" for value in case["scientific"]))
        lines.append("}}, {{")
        lines.append(",\n".join(value + "L" for value in case["reference_gap"]))
        lines.append("}}},")
    lines.extend(["}};", "} // namespace sirius::test::retained_camera"])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
