"""Freeze independent 75/105-digit Hamiltonian RK witnesses. Requires mpmath.

Geometry uses the defining metric, generic inversion and numerical high
precision differentiation from cpu_critical/reference.py. It never imports the
device's arithmetic graph, symbolic differentiation or retained pair operations.
Ordinary builds consume the checked-in header without regenerating this oracle.
"""
import sys
sys.dont_write_bytecode = True

import importlib.util
import json
from pathlib import Path
import struct

import mpmath as mp

folder = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("independent_geometry",
                                             folder.parent / "cpu_critical/reference.py")
ref = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ref)


def float32(value):
    return struct.unpack("<f", struct.pack("<f", float(value)))[0]


def pair(value):
    high = float32(value)
    low = float32(value - mp.mpf(high))
    tail = float32(value - mp.mpf(high) - mp.mpf(low))
    remaining = abs(value - mp.mpf(high) - mp.mpf(low) - mp.mpf(tail))
    bits = struct.unpack("<I", struct.pack("<f", float(remaining)))[0]
    radius = struct.unpack("<f", struct.pack("<I", bits + int(remaining != 0)))[0]
    return struct.unpack("<5I", struct.pack("<ffffI", high, low, tail, radius, 1))


def integrate(phase, h):
    rational = lambda a, b: mp.mpf(a) / b
    table = [[], [rational(1, 5)], [rational(3, 40), rational(9, 40)],
             [rational(44, 45), rational(-56, 15), rational(32, 9)],
             [rational(19372, 6561), rational(-25360, 2187), rational(64448, 6561),
              rational(-212, 729)],
             [rational(9017, 3168), rational(-355, 33), rational(46732, 5247),
              rational(49, 176), rational(-5103, 18656)],
             [rational(35, 384), 0, rational(500, 1113), rational(125, 192),
              rational(-2187, 6784), rational(11, 84)]]
    low = [rational(5179, 57600), 0, rational(7571, 16695), rational(393, 640),
           rational(-92097, 339200), rational(187, 2100), rational(1, 40)]
    slopes = []
    for stage in range(7):
        current = phase + h * sum((coefficient * slope for coefficient, slope
                                  in zip(table[stage], slopes)), mp.zeros(40, 1))
        slopes.append(ref.rhs(current))
    lower = phase + h * sum((coefficient * slope for coefficient, slope in zip(low, slopes)),
                           mp.zeros(40, 1))
    return [value for group in [current, lower, current - phase, current - lower]
            for value in group]


def calculate(row, incoming=False):
    ref.reflection = [-1, 1, -1, 1]
    phase = ref.initialize(row["previous"])
    if incoming:
        for i in range(40):
            phase[i] *= ref.reflection[i % 4]
        ref.reflection = [1, 1, 1, 1]
    h = mp.mpf(row["h"])
    inputs = [*ref.params, *phase, 1 if incoming else -1, h]
    return [word for value in inputs for word in pair(value)], integrate(phase, h)


def main():
    rows = json.loads((folder.parent / "cpu_critical/projected_inputs.json").read_text())["states"]
    cases = []
    for index, (row, incoming) in enumerate([(row, False) for row in rows] + [(rows[0], True)]):
        samples = []
        packets = []
        for precision in [75, 105]:
            mp.mp.dps = precision
            packet, values = calculate(row, incoming)
            packets.append(packet)
            samples.append(values)
        assert packets[0] == packets[1]
        gaps = [abs(low - high) for low, high in zip(*samples)]
        assert max(gaps) < mp.mpf("1e-55")
        cases.append({"name": "incoming-reflection" if incoming else "critical-" + str(index),
                      "input": packets[1], "reference": [str(v) for v in samples[1]],
                      "precision_gap": [str(v) for v in gaps]})
        print(cases[-1]["name"], "precision gap", max(gaps), flush=True)
    # Analytic flat-space control: every derivative and increment is dyadic.
    mp.mp.dps = 105
    phase = mp.matrix([1, 2, 3, 4, 1, 0, 0, 1,
                       0, 0, 0, 0, 0, mp.mpf(1) / 1024, 0, 0,
                       0, 0, 0, 0, 0, 0, mp.mpf(1) / 1024, 0,
                       0, 1, 0, 0, 0, 0, 0, 0,
                       0, 0, 1, 0, 0, 0, 0, 0])
    slope = mp.zeros(40, 1)
    for start in [0, 8, 16, 24, 32]:
        for axis in range(4):
            slope[start + axis] = phase[start + 4 + axis] * (-1 if axis == 0 else 1)
    h = mp.mpf(1) / 4
    finish = phase + h * slope
    cases.append({"name": "analytic-flat", "input": [word for value in [0, 0, 0, 0, *phase, 1, h]
                                                         for word in pair(mp.mpf(value))],
                  "reference": [str(v) for group in [finish, finish, h * slope, mp.zeros(40, 1)]
                                for v in group], "precision_gap": ["0"] * 160})
    document = {"schema": "sirius-retained-transport-reference-v2", "precision": [75, 105],
                "description": "Finite precision witnesses, not an interval certificate or trajectory admission.",
                "cases": cases}
    (folder / "reference_cases.json").write_text(json.dumps(document, indent=2) + "\n")
    lines = ["// Independent defining-metric and analytic-flat witnesses. See README.md.",
             "// clang-format off", "#pragma once", "#include <array>", "#include <cstdint>",
             "namespace sirius::test::retained_transport {",
             "struct WideReference {double high,low;};",
             "struct Case {const char* name; std::array<std::uint32_t,230> input; "
             "std::array<long double,160> reference,precision_gap; "
             "std::array<WideReference,160> wide_reference;};",
             f"inline constexpr std::array<Case,{len(cases)}> cases{{{{"]
    for case in cases:
        lines.append("{" + json.dumps(case["name"]) + ", {{")
        for offset in range(0, 230, 16):
            lines.append(",".join(str(v) + "u" for v in case["input"][offset:offset + 16]) + ",")
        lines.append("}}, {{")
        lines.append(",\n".join(v + "L" for v in case["reference"]))
        lines.append("}}, {{")
        lines.append(",\n".join(v + "L" for v in case["precision_gap"]))
        lines.append("}}, {{")
        for text in case["reference"]:
            value = mp.mpf(text)
            high = float(value)
            low = float(value - mp.mpf(high))
            lines.append("{" + high.hex() + "," + low.hex() + "},")
        lines.append("}}},")
    lines.extend(["}};", "} // namespace sirius::test::retained_transport", "// clang-format on", ""])
    (folder / "reference_cases.h").write_text("\n".join(lines))


if __name__ == "__main__":
    main()
