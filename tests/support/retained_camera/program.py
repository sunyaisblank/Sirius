"""Compile a fixed retained camera arithmetic graph into a bounded instruction stream.

This is a development prototype, not the production trace launch. Python emits
operations and register references only; every scientific value, derivative,
low part, error radius and validity decision is evaluated on the device.
"""
import json, struct
from pathlib import Path
ops = []
cache = {}

class P:

    def __init__(self, i):
        self.i = i

    def __add__(a, b):
        return node(2, a, b)

    def __radd__(a, b):
        return p(b) + a

    def __sub__(a, b):
        return node(3, a, b)

    def __rsub__(a, b):
        return p(b) - a

    def __mul__(a, b):
        return node(4, a, b)

    def __rmul__(a, b):
        return p(b) * a

    def __truediv__(a, b):
        return node(5, a, b)

    def __rtruediv__(a, b):
        return p(b) / a

    def __neg__(a):
        return node(7, a)

    def sqrt(a):
        return node(6, a)

def p(x):
    if isinstance(x, P):
        return x
    key = (0, struct.unpack('<I', struct.pack('<f', x))[0], 0, 0)
    if key not in cache:
        cache[key] = len(ops)
        ops.append(key)
    return P(cache[key])

def node(op, a, b=0, c=0):
    key = (op, p(a).i, p(b).i, p(c).i)
    if key not in cache:
        cache[key] = len(ops)
        ops.append(key)
    return P(cache[key])

def inp(i):
    key = (1, i, 0, 0)
    cache[key] = len(ops)
    ops.append(key)
    return P(cache[key])

class J:

    def __init__(self, v, d=None):
        self.v = p(v)
        self.d = [p(0)] * 4 if d is None else d

    def __add__(a, b):
        b = j(b)
        return J(a.v + b.v, [x + y for x, y in zip(a.d, b.d)])

    def __radd__(a, b):
        return a + j(b)

    def __neg__(a):
        return J(-a.v, [-x for x in a.d])

    def __sub__(a, b):
        return a + -j(b)

    def __rsub__(a, b):
        return j(b) + -a

    def __mul__(a, b):
        b = j(b)
        return J(a.v * b.v, [x * b.v + a.v * y for x, y in zip(a.d, b.d)])

    def __rmul__(a, b):
        return a * j(b)

    def __truediv__(a, b):
        b = j(b)
        v = a.v / b.v
        return J(v, [(x - v * y) / b.v for x, y in zip(a.d, b.d)])

    def __rtruediv__(a, b):
        return j(b) / a

    def sqrt(a):
        v = a.v.sqrt()
        return J(v, [x / (2 * v) for x in a.d])

def j(x):
    return x if isinstance(x, J) else J(x)

def coord(x, c):
    return J(x, [p(int(i == c)) for i in range(4)])

def choose(a, b, c, op=10):
    a, b = (j(a), j(b))
    return J(node(op, a.v, b.v, c), [node(op, x, y, c) for x, y in zip(a.d, b.d)])

def metric(position, row):
    M, a, Q, L = map(j, row[:4])
    xx, yy, zz = [coord(position[i], i) for i in range(1, 4)]
    bound = node(9, node(9, position[1], position[2]), node(9, position[3], a.v))
    scale = J(node(8, bound))
    x, y, z = (xx / scale, yy / scale, zz / scale)
    s = a / scale
    a2 = s * s
    z2 = z * z
    R2 = x * x + y * y + z2
    reduced = R2 - a2
    root = (reduced * reduced + 4 * a2 * z2).sqrt()
    r2 = choose(0.5 * (reduced + root), 2 * a2 * z2 / (root - reduced), reduced.v)
    r = scale * r2.sqrt()
    den = r * r + a * a
    ell = [j(1), (r * xx + a * yy) / den, (r * yy - a * xx) / den, zz / r]
    H = (2 * M * r - Q * Q) / (r * r + a * a * ell[3] * ell[3]) + L * r * r / 3
    position_bound = node(9, node(9, position[1], position[2]), position[3])
    flat = node(9, node(9, M.v, a.v), node(9, Q.v, L.v * position_bound))
    g = []
    inv = []
    for i in range(4):
        gr = []
        ir = []
        for k in range(4):
            eta = -1 if i == k == 0 else int(i == k)
            correction = H * ell[i] * ell[k]
            gr.append(choose(eta, eta + correction, flat, 11))
            ir.append(choose(eta, eta + correction if (i == 0) != (k == 0) else eta - correction, flat, 11))
        g.append(gr)
        inv.append(ir)
    return (g, inv)

def frame(g, inv, row):
    lapse = 1 / (-inv[0][0]).sqrt()
    ref = [[-lapse * inv[i][0] for i in range(4)]]

    def dot(a, b):
        return sum((g[i][k] * a[i] * b[k] for i in range(1, 4) for k in range(1, 4)))
    for axis in range(3):
        e = [j(0)] + [j(row[8 + 3 * axis + i]) * (-1 if axis == 1 else 1) for i in range(3)]
        for old in ref[1:]:
            d = dot(e, old)
            e = [x - d * y for x, y in zip(e, old)]
        norm = dot(e, e).sqrt()
        ref.append([v / norm for v in e])
    beta = row[17:20]
    b2 = sum((v * v for v in beta))
    s = (1 - b2).sqrt()
    gamma = 1 / s
    coef = 1 / (s * (1 + s))
    bv = [sum((ref[a + 1][mu] * beta[a] for a in range(3))) for mu in range(4)]
    return [[(ref[0][mu] + bv[mu]) * gamma for mu in range(4)]] + [[ref[a + 1][mu] + bv[mu] * (coef * beta[a]) + ref[0][mu] * (gamma * beta[a]) for mu in range(4)] for a in range(3)]

def build():
    ops.clear()
    cache.clear()
    row = [inp(i) for i in range(32)]
    position = row[4:8]
    X = [[p(0) for _ in range(4)] for _ in range(4)]
    for phase in range(2):
        cg, ci = metric(position, row)

        def contract(g):
            return [[J(v.v, [sum((v.d[mu] * X[c][mu] for mu in range(4))) for c in range(4)]) for v in gr] for gr in g]
        f = frame(contract(cg), contract(ci), row)
        if phase == 0:
            X[2] = [v.v * row[29] for v in f[3]]
            X[3] = [v.v * row[29] for v in f[2]]
            position = [x + X[3][mu] * row[28] + X[2][mu] * row[27] for mu, x in enumerate(position)]
    z = [coord(row[25 + i], i) for i in range(4)]
    fx, fy, pr, pu = z
    D = 1 + row[29] * (row[24] - 1)
    T = row[23]
    q = [j(D), j(D * T) * (1 - 2 * fy / row[21]) - pu * row[29], j(D * T) * (2 * fx / row[20] - 1) * row[22] - pr * row[29]]
    norm = sum((v * v for v in q)).sqrt()
    n = [v / norm for v in q]
    k = [-f[0][mu] + sum((f[a + 1][mu] * n[a] for a in range(3))) for mu in range(4)]
    G = [[[sum((ci[mu][ss].v * (cg[ss][b].d[a] + cg[ss][a].d[b] - cg[a][b].d[ss]) for ss in range(4))) * 0.5 for b in range(4)] for a in range(4)] for mu in range(4)]
    K = [[k[mu].d[c] for mu in range(4)] for c in range(4)]
    V = [[K[c][mu] + sum((G[mu][a][b] * k[a].v * X[c][b] for a in range(4) for b in range(4))) for mu in range(4)] for c in range(4)]
    du = [[f[0][mu].d[c] for mu in range(4)] for c in range(4)]
    framevalues = [v.v for axis in f for v in axis]
    scientific = position + [v.v for v in k] + framevalues + [v for group in [X, K, V, du] for col in group for v in col] + framevalues
    assert len(scientific) == 104
    return [v.i for v in scientific]

def build_program():
    outputs = build()
    live = set()

    def visit(i):
        if i in live:
            return
        live.add(i)
        op, a, b, c = ops[i]
        if op not in (0, 1):
            visit(a)
            if op not in (6, 7, 8):
                visit(b)
            if op in (10, 11):
                visit(c)
    for i in outputs:
        visit(i)
    sequence = sorted(live)
    last_use = {i: i for i in sequence}
    for old in sequence:
        op, a, b, c = ops[old]
        dependencies = [] if op in (0, 1) else [a] + ([] if op in (6, 7, 8) else [b]) + ([c] if op in (10, 11) else [])
        for dependency in dependencies:
            last_use[dependency] = old
    for output in outputs:
        last_use[output] = len(ops)
    slots = {}
    free = []
    registers = 0
    program = []
    for old in sequence:
        op, a, b, c = ops[old]
        dependencies = [] if op in (0, 1) else [a] + ([] if op in (6, 7, 8) else [b]) + ([c] if op in (10, 11) else [])
        if op not in (0, 1):
            a = slots[a]
            b = slots[b] if op not in (6, 7, 8) else 0
            c = slots[c] if op in (10, 11) else 0
        for dependency in sorted(set(dependencies)):
            if last_use[dependency] == old:
                free.append(slots[dependency])
        if free:
            destination = free.pop()
        else:
            destination = registers
            registers += 1
        slots[old] = destination
        program.extend((op, destination, a, b, c))
    outputs = [slots[i] for i in outputs]
    return {'instructions': len(sequence), 'registers': registers, 'input_words': 138 + len(program), 'output_words': 384 + 4 * registers, 'outputs': outputs, 'operations': program}


def write_packets(program, source, destination):
    data = Path(source).read_bytes()
    if not data or len(data) % 128:
        raise ValueError("expected complete 32-word input records")
    rows = struct.iter_unpack("<32I", data)
    payload = b"".join(struct.pack(
        "<" + str(program["input_words"]) + "I", *row,
        program["instructions"], program["registers"],
        *program["outputs"], *program["operations"]) for row in rows)
    Path(destination).write_bytes(payload)


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path)
    parser.add_argument("--packet-output", type=Path)
    args = parser.parse_args()
    if bool(args.inputs) != bool(args.packet_output):
        parser.error("--inputs and --packet-output must be supplied together")
    program = build_program()
    output = Path(__file__).resolve().parents[3] / "out" / "retained-camera-build"
    output.mkdir(parents=True, exist_ok=True)
    (output / "camera-program.json").write_text(json.dumps(program, indent=2) + "\n")
    if args.inputs:
        write_packets(program, args.inputs, args.packet_output)
    print(f"{program['instructions']} instructions; {program['registers']} live registers")


if __name__ == "__main__":
    main()
