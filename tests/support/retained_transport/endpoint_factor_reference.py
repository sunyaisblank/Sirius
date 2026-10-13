"""Independent 100-root Endpoint witnesses for the factored continuation.

Requires mpmath. Ordinary builds consume frozen fixtures; regeneration is an
explicit operation. Geometry comes from the defining metric, generic inverse
and mp.diff in the unchanged endpoint_reference.py dependency chain. No
production graph, retained arithmetic or observed device output is imported.
The 75/105-digit agreement is a finite refinement witness, not a certified
input-ball, projection or trajectory enclosure.
"""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import struct
import sys

sys.dont_write_bytecode = True
FOLDER = Path(__file__).resolve().parent
ROOT = FOLDER.parents[2]
PRECISIONS = (75, 105)


def load(path):
    spec = importlib.util.spec_from_file_location("independent_endpoint", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def identity(path):
    data = path.read_bytes()
    return {"bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}


def decode(mp, packet):
    assert len(packet) == 225
    return [sum((mp.mpf(struct.unpack("<f", struct.pack("<I", word))[0])
                 for word in packet[offset:offset + 3]), mp.mpf(0))
            for offset in range(0, 225, 5)]


def profiles(mp):
    zero, one = mp.mpf(0), mp.mpf(1)
    return [
        ("weak-kerr", [one, mp.mpf(".7"), zero, zero], [50, 2, 10]),
        ("near-extremal-strong", [one, mp.mpf(".998"), zero, zero],
         [mp.mpf("2.2"), mp.mpf(".3"), -mp.mpf(".4")]),
        ("charged-kerr", [one, mp.mpf(".7"), mp.mpf(".25"), zero],
         [mp.mpf(7) / 4, mp.mpf(1) / 2, mp.mpf(3) / 8]),
        ("positive-lambda", [one, zero, zero, mp.mpf(".0001")], [6, .5, -.25]),
        ("negative-lambda", [one, zero, zero, -mp.mpf(".0001")], [6, .5, -.25]),
        ("kerr-axis", [one, mp.mpf(".7"), zero, zero], [0, 0, mp.mpf("3.3")]),
        ("scale-below-4", [one, mp.mpf(".7"), zero, zero],
         [4 - mp.power(2, -40), 0, 0]),
        ("scale-above-4", [one, mp.mpf(".7"), zero, zero],
         [4 + mp.power(2, -40), 0, 0]),
        # Small mass weakens derivatives while keeping the unchanged metric
        # prefix's radius-squared division inside its multiplication domain.
        ("weak-large-gradient", [mp.power(2, -40), mp.mpf(".7"), zero, zero],
         [8, 4, 2]),
        ("small-radius-large-X", [mp.power(2, -80), zero, zero, zero],
         [mp.power(2, -20), mp.power(2, -21), mp.power(2, -22)]),
    ]


def make_packet(endpoint, profile_index, chart):
    mp, ref = endpoint.mp, endpoint.ref
    name, parameters, position = profiles(mp)[profile_index]
    ref.params = parameters
    ref.reflection = [chart, 1, chart, 1]
    momentum_scale = mp.mpf(1)
    if name in ("weak-large-gradient", "small-radius-large-X"):
        momentum_scale = mp.power(2, 25 if name == "weak-large-gradient" else 10)
        position_scale = mp.power(2, 40 if name == "weak-large-gradient" else 45)
        columns = [[mp.mpf(node ** field) * position_scale for field in range(4)] +
                   [mp.mpf(node ** field) * momentum_scale for field in range(4)]
                   for node in (1, 2, 3, 4)]
        determinant = 12 * position_scale ** 4
    else:
        columns = [[mp.mpf(node ** field) / mp.power(2, 8 + field)
                    for field in range(8)] for node in (1, 2, 3, 4)]
        determinant = mp.mpf(12) / mp.power(2, 38)
    # Projection leaves X unchanged. Its Vandermonde determinant establishes
    # four independent coupled columns even after correcting their momenta.
    assert mp.det(mp.matrix([[columns[c][f] for c in range(4)] for f in range(4)])) == determinant
    momentum = [value * momentum_scale for value in
                (mp.mpf(7) / 8, mp.mpf(1) / 4, -mp.mpf(3) / 8, mp.mpf(1) / 2)]
    phase = [mp.mpf(1) / 8, *map(mp.mpf, position), *momentum,
             *(value for column in columns for value in column)]
    _, projected = endpoint.project(mp.matrix(phase))
    assert len(projected) == 80
    # Launch from the independent null projection, then evaluate the actual
    # decoded packet centres rather than the pre-packing high-precision state.
    return [word for value in [*parameters, *projected[:40], mp.mpf(chart)]
            for word in endpoint.transport.pair(value)]


def forced_projection(endpoint, phase, component):
    """The same coordinate-differentiated equations with a specified component.

    Only the existing guarded-temporal-root witness uses this function. Its two
    accepted spatial alternatives are independent of a device's tie decision.
    """
    mp, ref = endpoint.mp, endpoint.ref
    x, momentum, columns = ref.unpack(phase)
    g, inverse, dg, _, connection = ref.geometry(list(x), False)
    k = inverse * momentum
    a = g[component, component]
    b = 2 * sum(g[component, i] * k[i] for i in range(4) if i != component)
    d = sum(g[i, j] * k[i] * k[j] for i in range(4) for j in range(4)
            if i != component and j != component)
    assert a != 0 and b * b - 4 * a * d >= 0
    roots = [(-b + sign * mp.sqrt(b * b - 4 * a * d)) / (2 * a) for sign in (-1, 1)]
    projected = k.copy()
    projected[component] = min(roots, key=lambda value: abs(value - k[component]))
    physical, projected_columns = [*x, *projected], []
    for X, P in columns:
        coordinate = inverse * (P - sum((dg[a] * k * X[a] for a in range(4)),
                                        mp.zeros(4, 1)))
        numerator = sum((projected.T * dg[a] * projected)[0] * X[a] for a in range(4)) + (
            2 * sum(g[a, b] * projected[a] * coordinate[b]
                    for a in range(4) for b in range(4) if b != component))
        denominator = 2 * sum(g[component, a] * projected[a] for a in range(4))
        assert denominator != 0
        coordinate[component] = -numerator / denominator
        V = coordinate + mp.matrix([
            sum(connection[m][a][b] * projected[a] * X[b] for a in range(4) for b in range(4))
            for m in range(4)])
        next_P = g * coordinate + sum((dg[a] * projected * X[a] for a in range(4)),
                                      mp.zeros(4, 1))
        physical.extend([*X, *V])
        projected_columns.append((X, next_P))
    return [*ref.pack(x, g * projected, projected_columns), *physical]


def calculate(endpoint, packet, guarded=False, rank4=True):
    mp, ref = endpoint.mp, endpoint.ref
    values = decode(mp, packet)
    ref.params = values[:4]
    ref.reflection = [values[44], 1, values[44], 1]
    phase = mp.matrix(values[4:44])
    x, momentum, columns = ref.unpack(phase)
    if rank4:
        seeds = mp.matrix([[list(X)[i] if i < 4 else list(P)[i - 4]
                            for X, P in columns] for i in range(8)])
        assert mp.det(seeds.T * seeds) > 0
    g, inverse, _, _, _ = ref.geometry(list(x), False)
    tangent = inverse * momentum
    prefix = [g[i, j] for i in range(4) for j in range(4)] + list(tangent)
    if not guarded:
        component, output = endpoint.project(phase)
        return {component: prefix + output}
    # Existing dyadic witness: M=r=1, chart=+1, k=(-1,.5,.5-2^-90,0).
    # Both real temporal roots fail the unchanged relative denominator guard.
    # No floating-device word or observed selection participates in this test.
    a = g[0, 0]
    b = 2 * sum(g[0, i] * tangent[i] for i in range(1, 4))
    d = sum(g[i, j] * tangent[i] * tangent[j] for i in range(1, 4) for j in range(1, 4))
    discriminant = b * b - 4 * a * d
    assert discriminant > 0
    for sign in (-1, 1):
        root = (-b + sign * mp.sqrt(discriminant)) / (2 * a)
        terms = [g[0, i] * (root if i == 0 else tangent[i]) for i in range(4)]
        assert abs(sum(terms)) < mp.power(2, -44) * sum(abs(value) for value in terms)
    return {component: prefix + forced_projection(endpoint, phase, component)
            for component in (1, 2)}


def guarded_packet():
    # Byte-for-byte input from the independent dyadic control in
    # ProjectedEndpointsKeepPhysicalColumnsAndRetainedContinuation.
    packet = [word for _ in range(45) for word in (0, 0, 0, 0, 1)]
    for index, high in ((0, 0x3F800000), (5, 0x3F800000), (9, 0xBF000000), (44, 0x3F800000)):
        packet[5 * index] = high
    packet[50:55] = [0x3F000000, 0x92800000, 0, 0, 1]
    return packet


def freeze_case(endpoint, name, make_input, binding, guarded=False, rank4=True):
    mp = endpoint.mp
    packets, samples = [], []
    for precision in PRECISIONS:
        with mp.workdps(precision):
            packet = make_input()
            packets.append(packet)
            samples.append(calculate(endpoint, packet, guarded, rank4))
    assert packets[0] == packets[1]
    assert samples[0].keys() == samples[1].keys()
    alternatives = []
    with mp.workdps(105):
        for component in samples[1]:
            low, high = samples[0][component], samples[1][component]
            assert len(low) == len(high) == 100 and all(mp.isfinite(v) for v in high)
            gaps = [abs(a - b) for a, b in zip(low, high)]
            assert max(gaps) < mp.mpf("1e-55")
            alternatives.append({"component": component,
                                 "reference": [mp.nstr(v, 105) for v in high],
                                 "precision_gap": [mp.nstr(v, 105) for v in gaps]})
    payload = struct.pack("<225I", *packets[1])
    return {"name": name, "input": packets[1], "input_sha256": hashlib.sha256(payload).hexdigest(),
            "input_binding": binding, "rank4_coupled_columns_required": rank4,
            "alternatives": alternatives}


def header(cases, mp):
    lines = ["// Independent exact-decoded-centre witnesses. See endpoint_factor_reference.py.",
             "// clang-format off", "#pragma once", "#include <array>", "#include <cstddef>",
             "#include <cstdint>", "namespace sirius::test::retained_endpoint_factor {",
             "struct Wide {double high,low;};",
             "struct Expected {std::uint32_t component;std::array<Wide,100> reference;std::array<double,100> gap;};",
             "struct Case {const char* name;std::array<std::uint32_t,225> input;std::size_t alternatives;std::array<Expected,2> expected;};",
             f"inline constexpr std::array<Case,{len(cases)}> cases{{{{"]
    with mp.workdps(105):
        for case in cases:
            lines.append("{" + json.dumps(case["name"]) + ",{{" +
                         ",".join(str(v) + "u" for v in case["input"]) + "}}," +
                         str(len(case["alternatives"])) + ",{{")
            for alternative in case["alternatives"]:
                lines.append("{" + str(alternative["component"]) + ",{{")
                for text in alternative["reference"]:
                    value = mp.mpf(text)
                    high = float(value)
                    low = float(value - mp.mpf(high))
                    lines.append("{" + high.hex() + "," + low.hex() + "},")
                lines.append("}},{{" + ",".join(float(v).hex() for v in
                                                alternative["precision_gap"]) + "}}},")
            if len(case["alternatives"]) == 1:
                lines.append("{},")
            lines.append("}}},")
    lines.extend(["}};", "} // namespace sirius::test::retained_endpoint_factor",
                  "// clang-format on", ""])
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--generate", action="store_true", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    output = args.output_dir.resolve()
    output.relative_to(ROOT)
    paths = [output / ("endpoint_factor_reference." + suffix) for suffix in ("json", "h")]
    assert not any(path.exists() for path in paths), "preserve existing outputs before regeneration"
    authorities = [Path(__file__).resolve(), FOLDER / "endpoint_reference.py",
                   FOLDER / "reference.py", FOLDER.parent / "cpu_critical/reference.py",
                   FOLDER.parent / "retained_camera/reference.py",
                   FOLDER / "schwarzschild_reference_cases.json", FOLDER / "endpoint_reference.json",
                   ROOT / "tests/backend/retained_compute_test.cpp"]
    bindings = {str(path.relative_to(ROOT)): identity(path) for path in authorities}
    endpoint = load(FOLDER / "endpoint_reference.py")
    mp = endpoint.mp
    cases = []
    for index, (name, _, _) in enumerate(profiles(mp)):
        for chart in (-1, 1):
            label = name + ("-outgoing" if chart == -1 else "-ingoing")
            cases.append(freeze_case(endpoint, label,
                                    lambda i=index, c=chart: make_packet(endpoint, i, c),
                                    {"authority": "tests/support/retained_transport/endpoint_factor_reference.py",
                                     "profile": name, "chart": chart,
                                     "construction": "independent null projection before packing"}))
    frozen = json.loads((FOLDER / "schwarzschild_reference_cases.json").read_text())
    physical = [item for item in frozen["cases"] if item["name"].startswith(
        ("null_weak-", "public_mass_min-", "public_mass_max-"))]
    assert len(physical) == 6
    for item in physical:
        packet = item["input"][:225]
        cases.append(freeze_case(endpoint, "frozen-" + item["name"], lambda p=packet: list(p),
                                {"authority": "tests/support/retained_transport/schwarzschild_reference_cases.json",
                                 "case": item["name"], "source_words": 230,
                                 "retained_prefix_words": 225, "omitted_affine_words": 5}))
    original = json.loads((FOLDER / "endpoint_reference.json").read_text())
    for name in ("exact-linear-horizon", "spatial-ergoregion-root"):
        item = next(case for case in original["cases"] if case["name"] == name)
        packet = item["input"]
        cases.append(freeze_case(endpoint, "frozen-" + name, lambda p=packet: list(p),
                                {"authority": "tests/support/retained_transport/endpoint_reference.json", "case": name}))
    cases.append(freeze_case(endpoint, "guarded-temporal-roots-use-spatial-fallback", guarded_packet,
                             {"authority": "tests/backend/retained_compute_test.cpp",
                              "test": "ProjectedEndpointsKeepPhysicalColumnsAndRetainedContinuation",
                              "control": "guarded temporal roots use spatial fallback"},
                             guarded=True, rank4=False))
    assert len(cases) == 29
    assert bindings == {str(path.relative_to(ROOT)): identity(path) for path in authorities}
    document = {"classification": "finite independent exact-decoded-centre Endpoint witnesses",
                "precisions": list(PRECISIONS), "reference_order":
                "metric[16],unprojected tangent[4],projected phase[40],physical[40]",
                "case_count": len(cases), "input_words": 225, "roots_per_alternative": 100,
                "input_bindings": bindings, "source_and_inputs_unchanged": True,
                "mpmath_version": mp.__version__, "cases": cases,
                "limits": ["Raw arithmetic witnesses, not physical trajectories or device acceptance.",
                           "Observed precision agreement is not a certified reference enclosure.",
                           "Original seventeen fixture files and expectations are unchanged.",
                           "Large-gradient and small-radius/large-X families are finite raw-arithmetic witnesses, not typed physical mass or trajectory claims.",
                           "Guarded spatial alternatives express the existing component1/2 control, not observed selection.",
                           "Input radii/zero-word/refusal behavior requires separate connected device controls."]}
    output.mkdir(parents=True, exist_ok=True)
    paths[0].write_text(json.dumps(document, indent=2, sort_keys=True) + "\n")
    paths[1].write_text(header(cases, mp))
    print(json.dumps({"cases": len(cases),
                      "root_alternatives": sum(len(case["alternatives"]) for case in cases),
                      "device_calls": 0,
                      "outputs": {str(path): identity(path) for path in paths}}), flush=True)


if __name__ == "__main__":
    main()
