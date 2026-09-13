"""Emit bounded retained camera and coupled Hamiltonian arithmetic programs.

Only opcodes and register ownership are generated on the host. The device
executes every operation, retaining each low part, arithmetic radius and status.
"""
import struct
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

def metric_profile(position, row):
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
    return H, ell, flat


def metric(position, row):
    H, ell, flat = metric_profile(position, row)
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

def compile_program(outputs):
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
    return {'instructions': len(sequence), 'registers': registers, 'outputs': outputs, 'operations': program}



def build_camera_program():
    program = compile_program(build())
    program['input_words'] = 138 + len(program['operations'])
    program['output_words'] = 384 + 4 * program['registers']
    return program


def differentiate(value, seeds):
    memo = {}
    zero = p(0)
    def add(a, b):
        return b if a.i == zero.i else a if b.i == zero.i else a + b
    def multiply(a, b):
        return zero if a.i == zero.i or b.i == zero.i else a * b
    def negate(a):
        return zero if a.i == zero.i else -a
    def divide(a, b):
        return zero if a.i == zero.i else a / b
    def d(index):
        if index in memo:
            return memo[index]
        op, a, b, c = ops[index]
        x, y = P(a), P(b)
        if op == 0:
            result = p(0)
        elif op == 1:
            result = seeds.get(a, p(0))
        elif op == 2:
            result = add(d(a), d(b))
        elif op == 3:
            result = add(d(a), negate(d(b)))
        elif op == 4:
            result = add(multiply(d(a), y), multiply(x, d(b)))
        elif op == 5:
            result = divide(add(d(a), negate(multiply(P(index), d(b)))), y)
        elif op == 6:
            result = divide(d(a), 2 * P(index))
        elif op == 7:
            result = negate(d(a))
        elif op in (8, 9):
            result = p(0)  # computational power-of-two scale, not a physical coordinate
        elif op in (10, 11):
            result = zero if d(a).i == zero.i and d(b).i == zero.i else node(op, d(a), d(b), P(c))
        else:
            raise ValueError(op)
        memo[index] = result
        return result
    return [d(v.i) for v in value]


def build_hamiltonian_rhs():
    ops.clear()
    cache.clear()
    row = [inp(i) for i in range(45)]
    M, spin, charge, cosmological = row[:4]
    # Record: x[4], p[4], followed by (X[4], P[4]) for each column.
    position, momentum = row[4:8], row[8:12]
    reflection = [row[44], p(1), row[44], p(1)]
    reflected = [v * s for v, s in zip(position, reflection)]
    H, ell, flat = metric_profile(reflected, row)
    raised = [ell[i].v * reflection[i] * (-1 if i == 0 else 1) for i in range(4)]
    contraction = sum(raised[i] * momentum[i] for i in range(4))
    tangent = [node(11, momentum[i] * (-1 if i == 0 else 1),
                    momentum[i] * (-1 if i == 0 else 1) - H.v * raised[i] * contraction,
                    flat) for i in range(4)]
    correction = node(11, p(0), H.v * contraction * contraction, flat)
    hamiltonian = (sum(momentum[i] * momentum[i] * (-1 if i == 0 else 1)
                       for i in range(4)) - correction) / 2
    force = [-differentiate([hamiltonian], {4+i: p(1)})[0] for i in range(4)]
    central = tangent + force
    outputs = central[:]
    for column in range(4):
        seeds = {4+i: row[12 + 8*column+i] for i in range(8)}
        outputs.extend(differentiate(central, seeds))
    return [v.i for v in outputs]



def build_transport_program():
    return compile_program(build_hamiltonian_rhs())


def build_endpoint_program():
    """Metric and projected physical columns from the complete phase expansion.

    The endpoint kernel first consumes the metric/tangent outputs to select a
    null root, then evaluates the same program with that retained root and its
    component selector. No rounded coordinate variation is subtracted from a
    rounded connection to recover a small covariant variation.
    """
    ops.clear()
    cache.clear()
    row = [inp(i) for i in range(53)]
    g, inverse = chart_geometry(row[4:8], row, row[44])
    tangent = [sum(inverse[i][k] * row[8+k] for k in range(4)) for i in range(4)]
    projected = row[45:49]
    selected = row[49:53]
    delta = [projected[i] - tangent[i] for i in range(4)]
    covector = [sum(g[i][k].v * projected[k] for k in range(4)) for i in range(4)]
    denominator = sum(covector[i] * selected[i] for i in range(4))

    def first(a, b, c):
        return (g[a][c].d[b] + g[a][b].d[c] - g[b][c].d[a]) / 2

    phase = row[4:8] + covector
    physical = row[4:8] + projected
    for column in range(4):
        X = row[12+8*column:16+8*column]
        P = row[16+8*column:20+8*column]
        lowered = [P[mu] + sum((first(mu, a, b) * delta[a] -
                                first(a, mu, b) * tangent[a]) * X[b]
                               for a in range(4) for b in range(4)) for mu in range(4)]
        raw = [sum(inverse[mu][nu] * lowered[nu] for nu in range(4)) for mu in range(4)]
        numerator = sum(covector[i] * raw[i] * (1-selected[i]) for i in range(4))
        solved = -numerator / denominator
        V = [(1-selected[i]) * raw[i] + selected[i] * solved for i in range(4)]
        # Preserve P itself. Reconstructing it from g*V and a large connection
        # contraction would discard the very residual retained by transport.
        corrected = [P[mu] + sum(g[mu][a].d[b] * delta[a] * X[b]
                                 for a in range(4) for b in range(4)) +
                     sum(g[mu][a].v * selected[a] * (V[a]-raw[a]) for a in range(4))
                     for mu in range(4)]
        phase.extend(X + corrected)
        physical.extend(X + V)
    outputs = [v.v for line in g for v in line] + tangent + phase + physical
    assert len(outputs) == 100
    return compile_program([v.i for v in outputs])


def chart_geometry(position, row, chart):
    reflection = [chart, p(1), chart, p(1)]
    reflected = [position[i] * reflection[i] for i in range(4)]
    g, inverse = metric(reflected, row)
    g = [[J(g[i][k].v * reflection[i] * reflection[k],
            [g[i][k].d[a] * reflection[i] * reflection[k] * reflection[a]
             for a in range(4)]) for k in range(4)] for i in range(4)]
    inverse = [[inverse[i][k].v * reflection[i] * reflection[k]
                for k in range(4)] for i in range(4)]
    return g, inverse


def build_dense_program():
    """Retained Hermite flow and physical arrival derivatives on one segment."""
    ops.clear()
    cache.clear()
    row = [inp(i) for i in range(112)]
    start, end, increment = row[4:44], row[44:84], row[84:104]
    h, fraction = row[104:106]
    normal, moving, chart = row[106:110], row[110], row[111]

    def interpolate(first, last, slope0, slope1, delta):
        secant = delta/h
        quadratic = 3*(secant-slope0)-(slope1-slope0)
        cubic = (slope1-slope0)-2*(secant-slope0)
        value = first + (slope0+(quadratic+cubic*fraction)*fraction)*(h*fraction)
        slope = slope0 + (2*quadratic+3*cubic*fraction)*fraction
        # Fixed endpoints already own their covariant columns; preserve the
        # actual retained increments while avoiding connection round trips.
        value = node(11, first, node(11, first+delta, value, fraction-1), fraction)
        slope = node(11, slope0, node(11, slope1, slope, fraction-1), fraction)
        return value, slope

    central = [interpolate(start[i], end[i], start[i+4], end[i+4], increment[i]) for i in range(4)]
    position, tangent = [v[0] for v in central], [v[1] for v in central]
    geometries = [chart_geometry(x, row, chart) for x in [start[:4], end[:4], position]]

    def contraction(geometry, k, X):
        g, inverse = geometry
        lower = [sum((g[mu][b].d[a]+g[mu][a].d[b]-g[a][b].d[mu])*k[a]*X[b]
                     for a in range(4) for b in range(4))/2 for mu in range(4)]
        return [sum(inverse[mu][nu]*lower[nu] for nu in range(4)) for mu in range(4)]

    denominator = sum(normal[i]*tangent[i] for i in range(4))
    outputs = position + tangent
    for column in range(4):
        X0, V0 = start[8+8*column:12+8*column], start[12+8*column:16+8*column]
        X1, V1 = end[8+8*column:12+8*column], end[12+8*column:16+8*column]
        C0 = contraction(geometries[0], start[4:8], X0)
        C1 = contraction(geometries[1], end[4:8], X1)
        values = [interpolate(X0[i], X1[i], V0[i]-C0[i], V1[i]-C1[i],
                              increment[4+4*column+i]) for i in range(4)]
        X, K = [v[0] for v in values], [v[1] for v in values]
        C = contraction(geometries[2], tangent, X)
        V = [node(11, V0[i], node(11, V1[i], K[i]+C[i], fraction-1), fraction) for i in range(4)]
        shift = node(11, p(0), -sum(normal[i]*X[i] for i in range(4))/denominator, moving)
        # Dk/dlambda=0 cancels the geodesic-flow acceleration and connection
        # terms exactly. Arrival changes X by k*shift and leaves covariant V.
        outputs.extend([X[i]+tangent[i]*shift for i in range(4)] + V)
    return compile_program([v.i for v in outputs])


def build_initialize_program():
    """Convert physical camera/trace columns to retained Hamiltonian phase."""
    ops.clear()
    cache.clear()
    row = [inp(i) for i in range(45)]
    position, tangent = row[4:8], row[8:12]
    g, _ = chart_geometry(position, row, row[44])
    momentum = [sum(g[mu][nu].v*tangent[nu] for nu in range(4)) for mu in range(4)]
    outputs = position + momentum
    for column in range(4):
        X, V = row[12+8*column:16+8*column], row[16+8*column:20+8*column]
        P = [sum(g[mu][nu].v*V[nu] for nu in range(4)) +
             sum((g[a][b].d[mu]+g[a][mu].d[b]-g[mu][b].d[a])*tangent[a]*X[b]
                 for a in range(4) for b in range(4))/2 for mu in range(4)]
        outputs.extend(X+P)
    return compile_program([v.i for v in outputs])


def build_ray_camera_program():
    """Metric launch from the smooth lens's physical direction and four seeds."""
    ops.clear()
    cache.clear()
    row = [inp(i) for i in range(45)]
    position = row[4:8]
    X = [[p(0) for _ in range(4)] for _ in range(4)]
    for phase in range(2):
        cg, ci = metric(position,row)
        def contract(g):
            return [[J(v.v,[sum(v.d[mu]*X[c][mu] for mu in range(4)) for c in range(4)])
                     for v in line] for line in g]
        f = frame(contract(cg),contract(ci),row)
        if phase == 0:
            X = [[f[3][mu].v*row[37+c]+f[2][mu].v*row[41+c] for mu in range(4)] for c in range(4)]
            position = [x+f[3][mu].v*row[20]+f[2][mu].v*row[21] for mu,x in enumerate(position)]
    n = [J(row[22+a],row[25+4*a:29+4*a]) for a in range(3)]
    norm = sum(v*v for v in n).sqrt()
    n = [v/norm for v in n]
    k = [-f[0][mu]+sum(f[a+1][mu]*n[a] for a in range(3)) for mu in range(4)]
    G = [[[sum(ci[mu][s].v*(cg[s][b].d[a]+cg[s][a].d[b]-cg[a][b].d[s])
                for s in range(4))/2 for b in range(4)] for a in range(4)] for mu in range(4)]
    K = [[k[mu].d[c] for mu in range(4)] for c in range(4)]
    V = [[K[c][mu]+sum(G[mu][a][b]*k[a].v*X[c][b] for a in range(4) for b in range(4))
          for mu in range(4)] for c in range(4)]
    du = [[f[0][mu].d[c] for mu in range(4)] for c in range(4)]
    framevalues = [v.v for axis in f for v in axis]
    outputs = position+[v.v for v in k]+framevalues+[v for group in [X,K,V,du] for column in group for v in column]+framevalues
    return compile_program([v.i for v in outputs])
