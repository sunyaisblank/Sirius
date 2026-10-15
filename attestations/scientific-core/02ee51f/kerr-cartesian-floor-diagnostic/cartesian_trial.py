"""Independent ideal Cartesian RK4 comparator; no GPU trial or projection replay."""
from decimal import Decimal, localcontext
from pathlib import Path
import hashlib
import json
import struct
import sys
import time


SOURCE = '02ee51f85e45e8342923bae1d9f26d50f1d39e6d'
EXPECTED_INPUT = {'bytes': 65734,
                  'sha256': 'e93e6f2e1925d086238f550dd849da13f1e7b7dc1ae465a3fadcf3118c2275b1'}
EPSILONS = ('1e-18', '1e-22')
PRECISION = 80
MAX_SECONDS = 60
MAX_RHS = 64
MAX_METRICS = 1000
MAX_STEPS = 16


def require(condition, message):
    if not condition:
        raise ValueError(message)


def pin(raw):
    return {'bytes': len(raw), 'sha256': hashlib.sha256(raw).hexdigest()}


def finite(values):
    require(all(value.is_finite() for value in values), 'Nonfinite Decimal state')


def exact(value):
    result = Decimal.from_float(value)
    finite((result,))
    return result


def binary32(bits):
    return exact(struct.unpack('<f', struct.pack('<I', bits))[0])


def raw_hex(text, size):
    require(type(text) is str and len(text) == 2*size, 'Wrong raw byte count')
    raw = bytes.fromhex(text)
    require(raw.hex() == text, 'Noncanonical raw bytes')
    return raw


def inverse(metric):
    # Generic partial-pivot Gauss-Jordan; no null-vector or mass-shell assumption.
    augmented = [list(row) + [Decimal(int(i == j)) for j in range(4)]
                 for i, row in enumerate(metric)]
    for column in range(4):
        pivot = max(range(column, 4), key=lambda row: abs(augmented[row][column]))
        require(augmented[pivot][column] != 0, 'Singular defining metric')
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        scale = augmented[column][column]
        augmented[column] = [value/scale for value in augmented[column]]
        for row in range(4):
            if row != column:
                scale = augmented[row][column]
                augmented[row] = [value-scale*other
                                  for value, other in zip(augmented[row], augmented[column])]
    result = [row[4:] for row in augmented]
    finite(value for row in result for value in row)
    residual = max(abs(sum(metric[i][k]*result[k][j] for k in range(4)) - int(i == j))
                   for i in range(4) for j in range(4))
    return result, residual


class Cartesian:
    def __init__(self, mass, spin, epsilon, work):
        self.mass, self.spin, self.epsilon = mass, spin, epsilon
        self.work = work

    def tick(self, kind):
        self.work[kind] += 1
        limit = {'rhs': MAX_RHS, 'metric': MAX_METRICS, 'step': MAX_STEPS}[kind]
        require(self.work[kind] <= limit, 'Cartesian reference work bound: ' + kind)
        require(time.monotonic() - self.work['started'] <= MAX_SECONDS,
                'Cartesian reference diagnostic wall bound')

    def metric(self, q):
        self.tick('metric')
        finite(q)
        x, y, z = q[1:]
        a2 = self.spin*self.spin
        reduced = x*x + y*y + z*z - a2
        root = (reduced*reduced + 4*a2*z*z).sqrt()
        r2 = ((reduced+root)/2 if reduced >= 0 else 2*a2*z*z/(root-reduced))
        require(r2 > 0, 'Nonunique/singular Kerr radial chart')
        r = r2.sqrt()
        h = 2*self.mass*r/(r2 + a2*z*z/r2)
        ell = (Decimal(1), (r*x+self.spin*y)/(r2+a2),
               (r*y-self.spin*x)/(r2+a2), z/r)
        eta = (-1, 1, 1, 1)
        metric = [[Decimal(eta[i] if i == j else 0) + h*ell[i]*ell[j]
                   for j in range(4)] for i in range(4)]
        finite(value for row in metric for value in row)
        return metric, r

    def rhs(self, state):
        self.tick('rhs')
        finite(state)
        q, u = state[:4], state[4:]
        metric, _ = self.metric(q)
        reciprocal, residual = inverse(metric)
        self.work['maximum_inverse_identity_residual'] = max(
            self.work['maximum_inverse_identity_residual'], residual)
        dg = [[[Decimal(0) for _ in range(4)] for _ in range(4)] for _ in range(4)]
        # The defining metric is stationary, so dg[0] is identically zero.
        for axis in range(1, 4):
            above, below = list(q), list(q)
            above[axis] += self.epsilon
            below[axis] -= self.epsilon
            require(above[axis] != below[axis], 'Metric displacement rounded away')
            plus, _ = self.metric(above)
            minus, _ = self.metric(below)
            dg[axis] = [[(plus[i][j]-minus[i][j])/(2*self.epsilon)
                         for j in range(4)] for i in range(4)]
        # Euler-Lagrange for L=.5*g_ab*u^a*u^b, valid at nonzero norm:
        # g_mu,nu*a^nu = .5*dg[mu,a,b]*u^a*u^b - dg[k,mu,j]*u^k*u^j.
        c = [sum(dg[mu][a][b]*u[a]*u[b] for a in range(4) for b in range(4))/2 -
             sum(dg[k][mu][j]*u[k]*u[j] for k in range(4) for j in range(4))
             for mu in range(4)]
        acceleration = tuple(sum(reciprocal[i][mu]*c[mu] for mu in range(4))
                             for i in range(4))
        finite(acceleration)
        return (*u, *acceleration)

    def rk4(self, state, step):
        self.tick('step')
        require(step.is_finite() and step > 0, 'Invalid affine interval')
        def shifted(rate, interval):
            return tuple(value+interval*change for value, change in zip(state, rate))
        a = self.rhs(state)
        b = self.rhs(shifted(a, step/2))
        c = self.rhs(shifted(b, step/2))
        d = self.rhs(shifted(c, step))
        result = tuple(value+step*(aa+2*bb+2*cc+dd)/6
                       for value, aa, bb, cc, dd in zip(state, a, b, c, d))
        finite(result)
        return result

    def summary(self, state, initial_norm, capture_radius, null_scale_floor, null_tolerance):
        q, u = state[:4], state[4:]
        metric, radius = self.metric(q)
        terms = [metric[i][j]*u[i]*u[j] for i in range(4) for j in range(4)]
        norm = sum(terms)
        term_sum = sum(abs(term) for term in terms)
        relative_null = abs(norm)/max(term_sum, null_scale_floor)
        p = [sum(metric[i][j]*u[j] for j in range(4)) for i in range(4)]
        return {
            'q': [str(x) for x in q], 'u': [str(x) for x in u],
            'p_covariant': [str(x) for x in p], 'gkk': str(norm),
            'norm_drift_from_captured_initial': str(norm-initial_norm),
            'absolute_metric_tangent_term_sum': str(term_sum),
            'controller_style_unprojected_null_residual': str(relative_null),
            'diagnostic_null_ratio_to_binary32_5e_4': str(relative_null/null_tolerance),
            'diagnostic_null_comparator_le_one': relative_null <= null_tolerance,
            'energy_minus_p_t': str(-p[0]),
            'axial_angular_momentum': str(q[1]*p[2]-q[2]*p[1]),
            'radius': str(radius), 'radius_minus_independent_buffered_capture': str(radius-capture_radius),
            'mass_shell_projection_performed': False,
        }


def trial_comparator(start, full, two_half, tolerance):
    # This is the source controller's scale only, not independent acceptance.
    scales = [max(Decimal(1), abs(x), abs(z)) for x, z in zip(start, two_half)]
    components = [abs(z-y)/scale for y, z, scale in zip(full, two_half, scales)]
    error = max(components)/15
    return {
        'component_order': ['t', 'x', 'y', 'z', 'u_t_contravariant', 'u_x', 'u_y', 'u_z'],
        'component_scales': [str(x) for x in scales],
        'component_scaled_differences_before_divide15': [str(x) for x in components],
        'relative_step_difference_divide15': str(error),
        'diagnostic_ratio_to_binary32_2e_5': str(error/tolerance),
        'diagnostic_ratio_le_one': error <= tolerance,
        'source_controller_scaling_comparator_only': True,
        'actual_GPU_trial_or_controller_acceptance': False,
    }


def endpoint_difference(summary, fine):
    q = [Decimal(value) for value in fine['q']]
    u = [Decimal(value) for value in fine['u_contravariant']]
    return {
        'q_maximum_absolute_difference': str(max(abs(Decimal(x)-y)
                                               for x, y in zip(summary['q'], q))),
        'u_maximum_absolute_difference': str(max(abs(Decimal(x)-y)
                                               for x, y in zip(summary['u'], u))),
        'norm_absolute_difference': str(abs(Decimal(summary['gkk'])-Decimal(fine['norm_2H']))),
        'reference_is_computed_canonical128_not_an_exact_or_certified_solution': True,
    }


def main():
    require(not sys.flags.optimize and len(sys.argv) == 3,
            'usage: cartesian_trial.py exact-local-reference-result.json new-output.json')
    input_path, output_path = map(Path, sys.argv[1:])
    original, script = input_path.read_bytes(), Path(__file__).read_bytes()
    require(pin(original) == EXPECTED_INPUT, 'Not the exact accepted local reference result')
    reference = json.loads(original)
    require(reference['schema'] == 'sirius-kerr-local-off-null-hamiltonian-reference-v1' and
            reference['analysis_completed'] is True and reference['captured_source_revision'] == SOURCE and
            reference['numerical_test_passed'] is False and
            reference['mass_shell_projection_performed'] is False and
            reference['null_only_Carter_formula_used'] is False, 'Wrong canonical reference authority')
    failure = json.loads(reference['raw_context_json'])
    require(failure['schema'] == 'sirius-trace-continuation-failure-v1' and
            failure['scalar_bytes'] == failure['packet_scalar_bytes'] == 4 and
            failure['record_bytes'] == 272 and failure['byte_order'] == 'little' and
            failure['action'] == 1 and failure['active'] == 64, 'Wrong captured layout/action')
    previous = raw_hex(failure['previous_record_hex'], 272)
    current = raw_hex(failure['current_record_hex'], 272)
    packet_raw = raw_hex(failure['packet_hex'], 280)
    packet, words = struct.unpack('<70f', packet_raw), struct.unpack('<70I', packet_raw)
    expected = failure['expected_identity']
    before, after = struct.unpack_from('<4I', previous, 240), struct.unpack_from('<4I', current, 240)
    require(list(struct.unpack_from('<4I', previous, 256)) ==
            list(struct.unpack_from('<4I', current, 256)) == expected and
            expected[2] == 0 and expected[3] == words[69] > 0 and
            before[:2] == (1, 0) and after[:2] == (2, 5) and
            after[2] == before[2]+1 and after[3] == before[3], 'Captured transition/identity mismatch')
    require(packet[2] == packet[5] == packet[6] == packet[46] == 0 and packet[3] == 1 and
            0 < packet[4] < 1 and packet[21] == 2000 and packet[68] == 1,
            'Only the exact captured asymptotically flat Kerr sector is represented')
    work = {'started': time.monotonic(), 'rhs': 0, 'metric': 0, 'step': 0,
            'maximum_inverse_identity_residual': Decimal(0)}
    with localcontext() as context:
        context.prec, context.rounding = PRECISION, 'ROUND_HALF_EVEN'
        start = tuple(exact(value) for value in struct.unpack_from('<8f', previous))
        require([str(x) for x in start[:4]] == reference['initial']['q'] and
                [str(x) for x in start[4:]] == reference['initial_exact_decoded_u'],
                'Initial q/u are not the same exact decoded binary32 state')
        mass, spin, floor = exact(packet[3]), exact(packet[4]), exact(packet[23])
        require(floor == Decimal(reference['exact_affine_interval']) ==
                exact(struct.unpack_from('<f', previous, 192)[0]) ==
                Decimal('0.0199999995529651641845703125'), 'Wrong represented affine floor')
        # Exact IEEE binary32 encodings of the unchanged shader literals.
        integration_tolerance = binary32(0x37a7c5ac)  # 2.0e-5f
        null_tolerance = binary32(0x3a03126f)         # 5.0e-4f
        null_scale_floor = binary32(0x1e3ce508)      # 1.0e-20f
        capture_radius = Decimal(reference['independent_buffered_capture_radius'])
        first = Cartesian(mass, spin, Decimal(EPSILONS[0]), work)
        initial_metric, _ = first.metric(start[:4])
        initial_norm = sum(initial_metric[i][j]*start[i+4]*start[j+4]
                           for i in range(4) for j in range(4))
        canonical_norm = Decimal(reference['initial']['norm_2H'])
        require(initial_norm != 0 and canonical_norm != 0, 'Off-null initial state was lost')
        algebra_join_tolerance = Decimal('1e-60')*max(Decimal(1), abs(initial_norm))
        require(abs(initial_norm-canonical_norm) <= algebra_join_tolerance,
                'Independent initial metric norm does not join the canonical initial norm')
        runs = {}
        for epsilon in EPSILONS:
            diagnostic = Cartesian(mass, spin, Decimal(epsilon), work)
            full = diagnostic.rk4(start, floor)
            half = diagnostic.rk4(start, floor/2)
            two_half = diagnostic.rk4(half, floor/2)
            summaries = {label: diagnostic.summary(state, initial_norm, capture_radius,
                                                  null_scale_floor, null_tolerance)
                         for label, state in [('full', full), ('two_half', two_half)]}
            fine = reference['runs'][epsilon]['128']
            require(fine['substeps'] == 128 and Decimal(fine['finite_difference_epsilon']) ==
                    Decimal(epsilon) and Decimal(fine['affine_interval']) == floor,
                    'Wrong canonical refinement endpoint')
            runs[epsilon] = {
                **summaries, 'source_scaling_comparator': trial_comparator(
                    start, full, two_half, integration_tolerance),
                'difference_from_matching_epsilon_canonical128': {
                    label: endpoint_difference(summary, fine['endpoint'])
                    for label, summary in summaries.items()},
            }
        sensitivity = {}
        for label in ('full', 'two_half'):
            a, b = (runs[epsilon][label] for epsilon in EPSILONS)
            sensitivity[label] = {
                key+'_maximum_absolute_difference': str(max(abs(Decimal(x)-Decimal(y))
                                                          for x, y in zip(a[key], b[key])))
                for key in ('q', 'u')}
            sensitivity[label]['gkk_absolute_difference'] = str(abs(Decimal(a['gkk'])-Decimal(b['gkk'])))
        result = {
            'schema': 'sirius-kerr-ideal-cartesian-floor-comparator-v1',
            'analysis_completed': True, 'captured_source_revision': SOURCE,
            'numerical_test_passed': False, 'actual_gpu_trial_reproduced': False,
            'product_or_numerical_settings_modified': False,
            'input': {'path': str(input_path), **pin(original)},
            'script': {'path': str(Path(__file__)), **pin(script)},
            'raw_context_json': reference['raw_context_json'],
            'decimal_precision': PRECISION, 'decimal_rounding': context.rounding,
            'mass': str(mass), 'spin': str(spin), 'affine_floor': str(floor),
            'initial_q': [str(x) for x in start[:4]], 'initial_u': [str(x) for x in start[4:]],
            'initial_independent_metric_norm': str(initial_norm),
            'initial_canonical_norm_2H': str(canonical_norm),
            'initial_norm_join_absolute_difference': str(abs(initial_norm-canonical_norm)),
            'algebra_join_absolute_tolerance_not_a_product_tolerance': str(algebra_join_tolerance),
            'binary32_comparator_literals': {
                'integration': {'bits': '37a7c5ac', 'exact_decimal': str(integration_tolerance)},
                'raw_null': {'bits': '3a03126f', 'exact_decimal': str(null_tolerance)},
                'null_scale_floor': {'bits': '1e3ce508', 'exact_decimal': str(null_scale_floor)}},
            'runs': runs, 'finite_difference_epsilon_sensitivity': sensitivity,
            'maximum_generic_inverse_identity_residual': str(work['maximum_inverse_identity_residual']),
            'work_counts': {key: work[key] for key in ('rhs', 'metric', 'step')},
            'bounds': {'maximum_seconds': MAX_SECONDS, 'maximum_rhs': MAX_RHS,
                       'maximum_metrics': MAX_METRICS, 'maximum_RK4_steps': MAX_STEPS},
            'elapsed_reference_seconds': time.monotonic()-work['started'],
            'scope': 'Independent defining-metric Euler-Lagrange Cartesian acceleration and ideal Decimal RK4 only, starting from the exact saved nonzero-norm q/u. No production metric, connection, projector or controller imports and no projection or null-only formula. Source-style step/null ratios are comparators, not an independent acceptance oracle, actual GPU binary32 trials or evidence identifying any rejected GPU conjunct. Fixed finite-difference sensitivity and canonical128 differences are not certified error bounds. No physics fix, product floor/tolerance changes or full-ray/scientific qualification.',
        }
    require(input_path.read_bytes() == original and Path(__file__).read_bytes() == script,
            'Input/source changed during comparison')
    with output_path.open('x', encoding='utf-8') as stream:
        stream.write(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'analysis_completed': True, 'numerical_test_passed': False,
                      'output': str(output_path)}))


if __name__ == '__main__':
    main()
