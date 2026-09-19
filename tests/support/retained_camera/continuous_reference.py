"""Freeze smooth camera inputs below binary32 film/pupil spacing (mpmath)."""
import sys
sys.dont_write_bytecode = True

import importlib.util
import json
from pathlib import Path
import struct
import mpmath as mp

folder = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("independent_camera", folder / "reference.py")
reference = importlib.util.module_from_spec(spec)
spec.loader.exec_module(reference)


def main():
    _, rows = reference.fixtures()
    base = rows[4][:]
    base[20:23] = [5616, 4096, reference.f32(5616 / 4096)]
    base[25:27] = [4095.75, 2048.5]
    film = base[:]
    film[25] += 2**-20
    film[26] -= 2**-21
    pupil = base[:]
    pupil[27] += 2**-30
    pupil[28] -= 2**-31
    cases = []
    for name, row in [("base", base), ("film-below-float-spacing", film),
                      ("pupil-below-float-spacing", pupil)]:
        low = reference.calculate(row, 100)["scientific"]
        high = reference.calculate(row, 180)["scientific"]
        with mp.workdps(180):
            gaps = [str(abs(mp.mpf(a) - mp.mpf(b))) for a, b in zip(low, high)]
            assert max(map(mp.mpf, gaps)) < mp.mpf("1e-90")
        cases.append({"name": name, "input": row, "reference": high, "reference_gap": gaps})
        print(name, flush=True)
    for field in [25, 26]:
        assert reference.f32(base[field]) == reference.f32(film[field])
    for field in [27, 28]:
        assert reference.f32(base[field]) == reference.f32(pupil[field])
    document = {"schema": "sirius-continuous-retained-camera-v1", "precision": [100, 180],
                "cases": cases}
    (folder / "continuous_reference.json").write_text(json.dumps(document, indent=2) + "\n")
    lines = ["// Independent complete-family differentiation. See README.md.",
             "// clang-format off", "#pragma once", "#include <array>",
             "namespace sirius::test::continuous_retained_camera {",
             "struct Case {const char* name; std::array<double,32> input; "
             "std::array<long double,104> reference,reference_gap;};",
             "inline constexpr std::array<Case,3> cases{{"]
    for case in cases:
        lines.append("{" + json.dumps(case["name"]) + ", {{" +
                     ",".join(float(v).hex() for v in case["input"]) + "}}, {{")
        lines.append(",\n".join(v + "L" for v in case["reference"]))
        lines.append("}}, {{")
        lines.append(",\n".join(v + "L" for v in case["reference_gap"]))
        lines.append("}}},")
    lines.extend(["}};", "} // namespace sirius::test::continuous_retained_camera",
                  "// clang-format on", ""])
    (folder / "continuous_reference.h").write_text("\n".join(lines))


if __name__ == "__main__":
    main()
