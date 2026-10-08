"""Independent rational DP continuous-extension witnesses.

Ordinary builds consume the checked-in header. Regeneration uses exact Fraction
arithmetic for the represented input problem, an expanded quartic differentiated
symbolically, and a conventional seven-weight dense coefficient. The separate
Kerr connection uses the existing defining-metric, generic-inverse, numerical
differentiation oracle at 75/105 digits; it imports no production arithmetic.
Coefficients: https://www.unige.ch/~hairer/prog/nonstiff/dopri5.f (CDOPRI).
Run: python3 tests/support/retained_transport/dopri_reference.py --generate
"""

import argparse
from collections import Counter
from decimal import Decimal, localcontext
from fractions import Fraction as F
import importlib.util
import json
import math
from pathlib import Path
import struct
import sys

sys.dont_write_bytecode = True
FOLDER = Path(__file__).resolve().parent
D = (F(-12715105075, 11282082432), F(0), F(87487479700, 32700410799),
     F(-10690763975, 1880347072), F(701980252875, 199316789632),
     F(-1453857185, 822651844), F(69997945, 29380423))
TABLE = ((), (F(1, 5),), (F(3, 40), F(9, 40)),
         (F(44, 45), F(-56, 15), F(32, 9)),
         (F(19372, 6561), F(-25360, 2187), F(64448, 6561), F(-212, 729)),
         (F(9017, 3168), F(-355, 33), F(46732, 5247), F(49, 176), F(-5103, 18656)),
         (F(35, 384), F(0), F(500, 1113), F(125, 192), F(-2187, 6784), F(11, 84)))
NODES = (F(0), F(1, 5), F(3, 10), F(4, 5), F(8, 9), F(1), F(1))
assert sum(D) == 0
assert tuple(sum(row) for row in TABLE) == NODES


def power2(exponent):
    return F(2 ** exponent) if exponent >= 0 else F(1, 2 ** -exponent)


def decoded(word):
    sign = -1 if word >> 31 else 1
    exponent, fraction = (word >> 23) & 255, word & 0x7fffff
    assert exponent != 255
    return sign * (F(fraction) * power2(-149) if exponent == 0 else
                   F(0x800000 + fraction) * power2(exponent - 150))


def nearest_even(numerator, denominator):
    q, r = divmod(numerator, denominator)
    return q + int(2 * r > denominator or (2 * r == denominator and q & 1))


def round32(value):
    """Exact nearest-even binary32 projection, including subnormal results."""
    value = F(value)
    if not value:
        return 0
    sign = 0x80000000 if value < 0 else 0
    value = abs(value)
    exponent = value.numerator.bit_length() - value.denominator.bit_length()
    if value < power2(exponent):
        exponent -= 1
    if exponent < -126:
        scaled = value / power2(-149)
        q = nearest_even(scaled.numerator, scaled.denominator)
        return sign | q
    scaled = value / power2(exponent - 23)
    q = nearest_even(scaled.numerator, scaled.denominator)
    if q == 1 << 24:
        exponent, q = exponent + 1, 1 << 23
    assert exponent <= 127
    return sign | ((exponent + 127) << 23) | (q - (1 << 23))


def retained(value, radius=F(0)):
    remaining, words = F(value), []
    for _ in range(3):
        word = round32(remaining)
        words.append(word)
        remaining -= decoded(word)
    bound = abs(remaining) + radius
    word = round32(bound)
    if decoded(word) < bound:
        word += 1
    return [*words, word, 1]


def raw(high, low=F(0), tail=F(0), radius=F(0)):
    values = [F(high), F(low), F(tail), F(radius)]
    words = [round32(value) for value in values]
    assert all(decoded(word) == value for word, value in zip(words, values))
    return [*words, 1]


def centers(words):
    return [sum(decoded(word) for word in words[i:i + 3]) for i in range(0, len(words), 5)]


def expansion(values):
    """Conventional weighted coefficient and expanded power-basis derivative."""
    y0, delta = values[:40], values[40:80]
    slopes = [values[80 + 40 * j:120 + 40 * j] for j in range(7)]
    h, s = values[360:]
    groups = [[], [], [], [], []]
    for i in range(40):
        a = h * slopes[0][i] - delta[i]
        b = 2 * delta[i] - h * (slopes[0][i] + slopes[6][i])
        c = h * sum(d * slope[i] for d, slope in zip(D, slopes))
        coefficients = (delta[i] + a, -a + b + c, -b - 2 * c, c)
        phase = y0[i] + sum(coefficient * s ** degree
                            for degree, coefficient in enumerate(coefficients, 1))
        derivative = sum(degree * coefficient * s ** (degree - 1)
                         for degree, coefficient in enumerate(coefficients, 1)) / h
        # A second exact evaluation checks the power expansion, not a rounded
        # implementation order. k2's coefficient vanishes, but its input exists.
        assert phase == y0[i] + s * (delta[i] + (1 - s) * (a + s * (b + (1 - s) * c)))
        if s == 0:
            assert phase == y0[i] and derivative == slopes[0][i]
        if s == 1:
            assert phase == y0[i] + delta[i] and derivative == slopes[6][i]
        for group, value in zip(groups, (phase, derivative, a, b, c)):
            group.append(value)
    return [value for group in groups for value in group]


def text(value):
    with localcontext() as context:
        context.prec = 110
        return str(Decimal(value.numerator) / Decimal(value.denominator))


def wide(value):
    high = float(value)
    low = float(value - F(high))
    gap = abs(value - F(high) - F(low))
    upper = float(gap)
    if F(upper) < gap:
        upper = math.nextafter(upper, math.inf)
    return high, low, upper


def case(name, kind, words, exact_phase=False, analytic=None):
    assert len(words) == 1810
    values = centers(words)
    result = expansion(values)
    entry = {"name": name, "kind": kind, "input": words,
             "exact_phase": exact_phase, "reference": [text(v) for v in result],
             "rational_reference": [[str(v.numerator), str(v.denominator)] for v in result]}
    if analytic is not None:
        gaps = [abs(a - b) for a, b in zip(result[:80], analytic)]
        assert max(gaps) < power2(-70)
        entry["analytic_encoding_gap"] = [text(v) for v in gaps]
    return entry


def fixtures():
    cases = []
    fractions = (F(0), F(1, 7), F(1, 3), F(1, 2), F(13, 16), F(1))
    for kind in ("arbitrary", "polynomial-ode", "large-origin", "noncanonical"):
        for s in fractions:
            h = F(3, 16) if kind == "arbitrary" else F(1, 8)
            y0, delta, slopes = [], [], [[] for _ in range(7)]
            analytic_phase, analytic_rate = [], []
            for i in range(40):
                sign = -1 if i & 1 else 1
                if kind == "arbitrary":
                    y0.append(retained(sign * F(i + 3, 8)))
                    delta.append(retained(sign * F(2 * i + 1, 128)))
                    for j in range(7):
                        slopes[j].append(retained(sign * F((j + 2) ** 2 * (i + 5), 256)))
                elif kind == "polynomial-ode":
                    # y'=alpha+beta*lambda+gamma*lambda^2+eta*lambda^3;
                    # DP tableau nodes and quadrature weights independently
                    # reproduce its defining integral and derivative exactly.
                    coefficients = [sign * F(i + n + 1, 32 * (n + 1)) for n in range(4)]
                    rhs = [sum(v * (node * h) ** degree
                               for degree, v in enumerate(coefficients)) for node in NODES]
                    increment = h * sum(w * k for w, k in zip(TABLE[-1], rhs))
                    integral = lambda t: sum(v * t ** (n + 1) / (n + 1)
                                              for n, v in enumerate(coefficients))
                    assert increment == integral(h)
                    origin = sign * F(i + 1, 4)
                    y0.append(retained(origin))
                    delta.append(retained(increment))
                    for j in range(7):
                        slopes[j].append(retained(rhs[j]))
                    # First establish exact extension versus analytic ODE before
                    # input encoding; afterward explicitly bound encoding error.
                    exact = [origin] * 40 + [increment] * 40
                    exact += [k for k in rhs for _ in range(40)] + [h, s]
                    reference = expansion(exact)
                    assert reference[0] == origin + integral(h * s)
                    assert reference[40] == sum(v * (h * s) ** n for n, v in enumerate(coefficients))
                    sampled_s = centers(retained(s))[0]
                    analytic_phase.append(origin + integral(h * sampled_s))
                    analytic_rate.append(sum(v * (h * sampled_s) ** n
                                             for n, v in enumerate(coefficients)))
                elif kind == "large-origin":
                    increment = sign * (i + 1) * power2(-20)
                    y0.append(retained(sign * power2(60)))
                    delta.append(retained(increment))
                    for j in range(7):
                        slopes[j].append(retained(increment / h))
                else:
                    # Raw opposing limbs promote a tiny retained tail. Other
                    # rows overlap limbs; none require a canonical input sum.
                    y0.append(raw(sign, -sign, sign * (i + 1) * power2(-60)))
                    delta.append(raw(sign * F(1, 2), -sign * F(1, 2),
                                     sign * (i + 1) * power2(-65)))
                    for j in range(7):
                        slopes[j].append(raw(sign * F(j + 1, 16), -sign * F(j + 1, 16),
                                             sign * (i + j + 1) * power2(-60)))
            words = [word for v in [*y0, *delta, *(v for stage in slopes for v in stage),
                                    retained(h), retained(s)] for word in v]
            analytic = analytic_phase + analytic_rate if analytic_phase else None
            cases.append(case(f"{kind}-{s}", kind, words, analytic=analytic))

    original = next(row for row in cases if row["name"] == "noncanonical-1/2")
    original_values = expansion(centers(original["input"]))
    for limb, name in ((1, "low-limb-removed"), (2, "tail-limb-removed")):
        words = list(original["input"])
        for offset in range(limb, len(words), 5):
            words[offset] = 0
        changed = case(name, "limb-negative-control", words)
        separation = abs(expansion(centers(words))[0] - original_values[0])
        assert separation > F(1, 10 ** 19)
        changed["phase_zero_separation_from_complete_limbs"] = text(separation)
        cases.append(changed)

    # Entire complete phase is reduced to sparse terms at the exact right
    # endpoint; losing low or tail limbs cannot hide behind relative tolerance.
    values = [raw((-1) ** i, (-1) ** i * power2(-35), (-1) ** i * power2(-120))
              for i in range(40)]
    values += [raw(-(-1) ** i, -(-1) ** i * power2(-35)) for i in range(40)]
    values += [retained(F((j + 1) * (i + 1), 256)) for j in range(7) for i in range(40)]
    values += [retained(1), retained(1)]
    cases.append(case("sparse-endpoint-cancellation", "sparse-endpoint",
                      [word for v in values for word in v], exact_phase=True))

    # h and fraction themselves contain cancelling/overlapping limbs. Origin
    # uncertainty belongs to phase, not its affine derivative.
    for name, h, s, radius in (
        ("noncanonical-h-and-fraction", raw(1, F(-3, 4), power2(-35)),
         raw(1, F(-3, 4), power2(-30)), F(0)),
        ("origin-radius", retained(F(1, 4)), retained(F(1, 2)), power2(-40)),
        ("tiny-positive-h", retained(power2(-60)), retained(F(5, 8)), F(0)),
        ("near-left-endpoint", retained(F(1, 8)), retained(power2(-30)), F(0)),
        ("near-right-endpoint", retained(F(1, 8)), raw(1, -power2(-30)), F(0)),
    ):
        hs = centers(h)[0]
        values = [retained(F((-1) ** i * (i + 1), 8), radius) for i in range(40)]
        values += [retained(hs * F(i + 1, 256)) for i in range(40)]
        values += [retained(F(i + 1, 256)) for _ in range(7) for i in range(40)]
        values += [h, s]
        cases.append(case(name, name, [word for v in values for word in v]))
    return cases


def kerr_connection():
    import mpmath as mp
    spec = importlib.util.spec_from_file_location("independent_dopri_projection",
                                                FOLDER / "endpoint_reference.py")
    endpoint = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(endpoint)
    original = next(row for row in json.loads((FOLDER / "reference_cases.json").read_text())["cases"]
                    if row["name"] == "ordinary-outgoing-.5")
    exact_input = centers(original["input"])
    samples = []
    for precision in (75, 105):
        mp.mp.dps = precision
        convert = lambda v: mp.mpf(v.numerator) / v.denominator
        values = list(map(convert, exact_input))
        endpoint.ref.params = values[:4]
        endpoint.ref.reflection = [values[44], 1, values[44], 1]
        y0, h, s = mp.matrix(values[4:44]), values[45], mp.mpf(3) / 8
        slopes = []
        for row in TABLE:
            state = y0 + h * sum((convert(weight) * slope for weight, slope in zip(row, slopes)),
                                  mp.zeros(40, 1))
            slopes.append(endpoint.ref.rhs(state))
        increment = state - y0
        a = h * slopes[0] - increment
        b = 2 * increment - h * (slopes[0] + slopes[6])
        c = h * sum((convert(weight) * slope for weight, slope in zip(D, slopes)), mp.zeros(40, 1))
        p1, p2, p3, p4 = increment + a, -a + b + c, -b - 2 * c, c
        phase = y0 + p1 * s + p2 * s ** 2 + p3 * s ** 3 + p4 * s ** 4
        rate = (p1 + 2 * p2 * s + 3 * p3 * s ** 2 + 4 * p4 * s ** 3) / h
        component, projected = endpoint.project(phase)
        samples.append([v for group in [*slopes, phase, rate, a, b, c] for v in group] + projected)
    gaps = [abs(a - b) for a, b in zip(*samples)]
    assert max(gaps) < mp.mpf("1e-55")
    return {"name": original["name"], "input": original["input"], "fraction": "3/8",
            "component": component, "reference": [mp.nstr(v, 110) for v in samples[-1]],
            "precision_gap": [mp.nstr(v, 110) for v in gaps]}


def render_header(cases, connected):
    lines = ["// Independent rational quartic and defining-metric DP witnesses. See dopri_reference.py.",
             "// clang-format off", "#pragma once", "#include <array>", "#include <cstdint>",
             "namespace sirius::test::retained_dopri {",
             "struct Wide {double high,low,gap;};",
             "struct Case {const char* name; bool exact_phase; std::array<std::uint32_t,1810> input; std::array<Wide,200> reference;};",
             f"inline constexpr std::array<Case,{len(cases)}> cases{{{{"]
    def words(values):
        for begin in range(0, len(values), 16):
            lines.append(",".join(str(v) + "u" for v in values[begin:begin + 16]) + ",")
    def references(values):
        for value in values:
            lines.append("{" + ",".join(v.hex() for v in wide(value)) + "},")
    for entry in cases:
        lines.append("{" + json.dumps(entry["name"]) + "," + str(entry["exact_phase"]).lower() + ",{{")
        words(entry["input"])
        lines.append("}},{{")
        references([F(int(n), int(d)) for n, d in entry["rational_reference"]])
        lines.append("}}},")
    lines += ["}};", "struct Connection {std::uint32_t component; std::array<std::uint32_t,230> input; std::array<Wide,560> reference;};",
              "inline constexpr Connection connection{" + str(connected["component"]) + ",{{"]
    words(connected["input"])
    lines.append("}},{{")
    # These high-precision values include 280 RHS + 200 sampler + 80 projection
    # fields. Incorporate the independently measured precision gap explicitly.
    for value, gap in zip(connected["reference"], connected["precision_gap"]):
        high, low, rounding = wide(F(value))
        upper = math.nextafter(rounding + float(gap), math.inf)
        lines.append("{" + ",".join(v.hex() for v in (high, low, upper)) + "},")
    lines += ["}}};", "} // namespace sirius::test::retained_dopri", "// clang-format on", ""]
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--generate", action="store_true", required=True)
    parser.parse_args()
    cases = fixtures()
    connected = kerr_connection()
    record = {"scope": "Finite represented-input quartic arithmetic and one independent Kerr connection; no event locator, trajectory or truncation certification.",
              "coefficient_source": "https://www.unige.ch/~hairer/prog/nonstiff/dopri5.f",
              "precision": [75, 105], "input_values": 362, "outputs_per_case": 200,
              "case_count": len(cases), "kinds": dict(Counter(row["kind"] for row in cases)),
              "coefficient_sum_exact": "0", "enclosure_padding": "1e-28*(1+abs(reference.high))",
              "accuracy_bound": "1e-10*(1+abs(reference.high))", "cases": cases, "connection": connected}
    (FOLDER / "dopri_reference.json").write_text(json.dumps(record, indent=2) + "\n")
    (FOLDER / "dopri_reference.h").write_text(render_header(cases, connected))
    print(json.dumps({"cases": len(cases), "rational_values": 200 * len(cases),
                      "connection_values": len(connected["reference"]), "precision": [75, 105]}))


if __name__ == "__main__":
    main()
