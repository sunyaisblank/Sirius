"""Independent projected endpoints, using coordinate differentiation at 75/105 digits."""
import sys
sys.dont_write_bytecode = True
import importlib.util
import json
from pathlib import Path
import struct
import mpmath as mp

folder = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("transport_reference", folder / "reference.py")
transport = importlib.util.module_from_spec(spec)
spec.loader.exec_module(transport)
ref = transport.ref


def project(phase):
    x, momentum, columns = ref.unpack(phase)
    g, inverse, dg, _, connection = ref.geometry(list(x), False)
    k = inverse * momentum
    options = []
    for c in range(4):
        a = g[c, c]
        b = 2 * sum(g[c, i] * k[i] for i in range(4) if i != c)
        d = sum(g[i, j] * k[i] * k[j] for i in range(4) for j in range(4) if i != c and j != c)
        if a == 0:
            roots = [-d/b] if b != 0 else []
        elif b*b-4*a*d >= 0:
            roots = [(-b+s*mp.sqrt(b*b-4*a*d))/(2*a) for s in [-1, 1]]
        else:
            roots = []
        options.extend((abs(root-k[c])/(1+abs(k[c])), c, root) for root in roots)
        if c == 0 and roots:
            break
    _, component, root = min(options)
    projected = k.copy()
    projected[component] = root
    physical = [*x, *projected]
    projected_columns = []
    for X, P in columns:
        coordinate = inverse * (P-sum((dg[a]*k*X[a] for a in range(4)), mp.zeros(4, 1)))
        numerator = sum((projected.T*dg[a]*projected)[0]*X[a] for a in range(4)) + 2*sum(
            g[a, b]*projected[a]*coordinate[b] for a in range(4) for b in range(4) if b != component)
        coordinate[component] = -numerator/(2*sum(g[component, a]*projected[a] for a in range(4)))
        V = coordinate + mp.matrix([sum(connection[m][a][b]*projected[a]*X[b]
                                       for a in range(4) for b in range(4)) for m in range(4)])
        next_P = g*coordinate + sum((dg[a]*projected*X[a] for a in range(4)), mp.zeros(4, 1))
        physical.extend([*X, *V])
        projected_columns.append((X, next_P))
    return component, [*ref.pack(x, g*projected, projected_columns), *physical]


def main():
    mp.mp.dps = 105
    transport_cases = json.loads((folder / "reference_cases.json").read_text())["cases"]
    columns = [0,0,0,0,0,0,mp.mpf(1)/1024,0,
               0,0,0,0,0,0,0,mp.mpf(1)/1024,
               0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0]
    for name, central in [
        ("exact-linear-horizon", [0,2,0,0,1,1,0,0]),
        ("spatial-ergoregion-root", [0,1,0,0,0,mp.mpf('-0.5'),mp.mpf('0.5000000001'),0])]:
        phase = [*central, *columns]
        transport_cases.append({"name":name, "reference":[str(v) for v in phase],
                                "input":[word for v in [1,0,0,0,*phase,1,1]
                                         for word in transport.pair(mp.mpf(v))]})
    cases = []
    for original in transport_cases:
        values = []
        for i in range(46):
            h, l, t, _, _ = struct.unpack("<ffffI", struct.pack("<5I", *original["input"][5*i:5*i+5]))
            values.append(mp.mpf(h)+mp.mpf(l)+mp.mpf(t))
        parameters = values[:4]
        chart = values[44]
        initial = list(map(mp.mpf, original["reference"][:40]))
        samples = []
        for precision in [75, 105]:
            mp.mp.dps = precision
            ref.params = parameters
            ref.reflection = [chart, 1, chart, 1]
            component, output = project(mp.matrix(initial))
            samples.append(output)
        gaps = [abs(a-b) for a, b in zip(*samples)]
        assert max(gaps) < mp.mpf("1e-55")
        cases.append({"name":original["name"], "component":component,
                      "input":[word for v in [*parameters, *initial, chart] for word in transport.pair(v)],
                      "reference":[str(v) for v in samples[1]], "precision_gap":[str(v) for v in gaps]})
        print(original["name"], "component", component, "gap", max(gaps), flush=True)
    (folder / "endpoint_reference.json").write_text(json.dumps({"precision":[75,105], "cases":cases}, indent=2)+"\n")
    lines = ["// Independent coordinate-differentiated projection. See endpoint_reference.py.",
             "// clang-format off", "#pragma once", "#include <array>", "#include <cstdint>",
             "namespace sirius::test::retained_endpoint {",
             "struct Wide {double high,low;};",
             "struct Case {const char* name; std::uint32_t component; std::array<std::uint32_t,225> input; std::array<Wide,80> reference;};",
             f"inline constexpr std::array<Case,{len(cases)}> cases{{{{"]
    for case in cases:
        lines.append("{"+json.dumps(case["name"])+","+str(case["component"])+",{{")
        for i in range(0, 225, 16):
            lines.append(",".join(str(v)+"u" for v in case["input"][i:i+16])+",")
        lines.append("}},{{")
        for text in case["reference"]:
            v = mp.mpf(text)
            high = float(v)
            low = float(v-mp.mpf(high))
            lines.append("{"+high.hex()+","+low.hex()+"},")
        lines.append("}}},")
    lines.extend(["}};", "} // namespace sirius::test::retained_endpoint", "// clang-format on", ""])
    (folder / "endpoint_reference.h").write_text("\n".join(lines))


if __name__ == "__main__":
    main()
