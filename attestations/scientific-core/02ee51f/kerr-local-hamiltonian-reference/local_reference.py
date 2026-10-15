"""Local canonical Kerr reference for one exact failed capture; no null projection."""
from decimal import Decimal, localcontext
from pathlib import Path
import hashlib
import json
import struct
import sys
import time


SOURCE = '02ee51f85e45e8342923bae1d9f26d50f1d39e6d'
INPUT_IDENTITY = {
    'bytes': 28046,
    'sha256': 'aa2eef9d357e9c0afad76dce0021eb0b65dd804b41b3f41d4b87b602fd644327',
}
PRECISION = 80
SUBSTEPS = (32, 64, 128)
EPSILONS = ('1e-18', '1e-22')
MAX_SECONDS = 240
MAX_RHS = 20000
MAX_GEOMETRIES = 200000
MAX_RK4 = 5000
MAX_BISECTIONS = 96
EVENT_RELATIVE_AFFINE_WIDTH = Decimal('1e-28')


def require(condition, message):
    if not condition:
        raise ValueError(message)


def identity(raw):
    return {'bytes': len(raw), 'sha256': hashlib.sha256(raw).hexdigest()}


def raw_hex(value, size):
    require(type(value) is str and len(value) == size * 2, 'Wrong raw byte count')
    raw = bytes.fromhex(value)
    require(raw.hex() == value, 'Noncanonical raw hex')
    return raw


def exact(value):
    result = Decimal.from_float(value)
    require(result.is_finite(), 'Nonfinite represented input')
    return result


def finite(values):
    require(all(value.is_finite() for value in values), 'Nonfinite reference state')


class Work:
    def __init__(self):
        self.started = time.monotonic()
        self.counts = {'rhs': 0, 'geometry': 0, 'rk4': 0}

    def tick(self, kind):
        self.counts[kind] += 1
        limit = {'rhs': MAX_RHS, 'geometry': MAX_GEOMETRIES, 'rk4': MAX_RK4}[kind]
        require(self.counts[kind] <= limit, 'Reference work bound exceeded: ' + kind)
        require(time.monotonic() - self.started <= MAX_SECONDS,
                'Reference diagnostic wall bound exceeded')


class Reference:
    def __init__(self, mass, spin, epsilon, work):
        self.mass = mass
        self.spin = spin
        self.epsilon = epsilon
        self.work = work

    def geometry(self, q):
        self.work.tick('geometry')
        finite(q)
        x, y, z = q[1:]
        a2 = self.spin * self.spin
        reduced = x*x + y*y + z*z - a2
        root = (reduced*reduced + 4*a2*z*z).sqrt()
        radius2 = ((reduced + root) / 2 if reduced >= 0
                   else 2*a2*z*z / (root - reduced))
        require(radius2 > 0, 'Reference encountered the nonunique/singular radial chart')
        radius = radius2.sqrt()
        sigma = radius2 + a2*z*z/radius2
        h = 2*self.mass*radius/sigma
        ell = ((radius*x + self.spin*y)/(radius2 + a2),
               (radius*y - self.spin*x)/(radius2 + a2), z/radius)
        finite((radius, h, *ell))
        return radius, h, ell

    def hamiltonian(self, q, p):
        finite(p)
        _, h, ell = self.geometry(q)
        w = -p[0] + sum(ell[i] * p[i+1] for i in range(3))
        result = (-p[0]*p[0] + sum(value*value for value in p[1:]) - h*w*w) / 2
        finite((result,))
        return result

    def velocity(self, q, p):
        # Analytic dH/dp for g^mu_nu=eta^mu_nu-h*l^mu*l^nu,
        # l^mu=(-1, ell). No mass-shell constraint is imposed.
        _, h, ell = self.geometry(q)
        w = -p[0] + sum(ell[i] * p[i+1] for i in range(3))
        result = (-p[0] + h*w,
                  *(p[i+1] - h*w*ell[i] for i in range(3)))
        finite(result)
        return result

    def initial_momentum(self, q, u):
        # Defining covariant g=eta+h*(1,ell)*(1,ell), not a null-cone solve.
        _, h, ell = self.geometry(q)
        contraction = u[0] + sum(ell[i] * u[i+1] for i in range(3))
        p = (-u[0] + h*contraction,
             *(u[i+1] + h*ell[i]*contraction for i in range(3)))
        finite(p)
        return p

    def rhs(self, state):
        self.work.tick('rhs')
        finite(state)
        q, p = state[:4], state[4:]
        dq = self.velocity(q, p)
        dp = [Decimal(0)]  # Stationary metric: p_t is constant, including off null.
        for axis in range(1, 4):
            above, below = list(q), list(q)
            above[axis] += self.epsilon
            below[axis] -= self.epsilon
            require(above[axis] != below[axis], 'Finite-difference displacement rounded away')
            dp.append(-(self.hamiltonian(above, p) - self.hamiltonian(below, p)) /
                      (2*self.epsilon))
        result = (*dq, *dp)
        finite(result)
        return result

    def rk4(self, state, delta):
        self.work.tick('rk4')
        require(delta > 0 and delta.is_finite(), 'Invalid affine reference interval')
        def shifted(rate, scale):
            return tuple(value + scale*change for value, change in zip(state, rate))
        k1 = self.rhs(state)
        k2 = self.rhs(shifted(k1, delta/2))
        k3 = self.rhs(shifted(k2, delta/2))
        k4 = self.rhs(shifted(k3, delta))
        result = tuple(value + delta*(a + 2*b + 2*c + d)/6
                       for value, a, b, c, d in zip(state, k1, k2, k3, k4))
        finite(result)
        require(result[4] == state[4], 'Stationary momentum changed')
        return result

    def summary(self, state, capture_radius, initial):
        q, p = state[:4], state[4:]
        radius, _, _ = self.geometry(q)
        h = self.hamiltonian(q, p)
        energy = -p[0]
        angular = q[1]*p[2] - q[2]*p[1]
        return {
            'q': [str(value) for value in q], 'p_covariant': [str(value) for value in p],
            'u_contravariant': [str(value) for value in self.velocity(q, p)],
            'hamiltonian': str(h), 'norm_2H': str(2*h),
            'energy_minus_p_t': str(energy), 'axial_angular_momentum': str(angular),
            'hamiltonian_drift': str(h - initial['H']),
            'norm_drift': str(2*(h - initial['H'])),
            'energy_drift': str(energy - initial['E']),
            'angular_momentum_drift': str(angular - initial['Lz']),
            'radius': str(radius), 'radius_minus_buffered_capture': str(radius - capture_radius),
        }

    def refine_event(self, start, finish, affine_start, step, floor,
                     capture_radius, initial):
        # The bracket endpoint at step is exactly the original RK4 step. Each
        # midpoint is a new partial integration from the fixed bracket start;
        # neither q nor p is linearly interpolated between saved endpoints.
        lo, hi = Decimal(0), step
        outside, inside = start, finish
        require(self.geometry(outside[:4])[0] > capture_radius and
                self.geometry(inside[:4])[0] <= capture_radius,
                'Invalid sampled capture bracket')
        width_goal = floor * EVENT_RELATIVE_AFFINE_WIDTH
        iterations = 0
        for _ in range(MAX_BISECTIONS):
            if hi - lo <= width_goal:
                break
            mid = (lo + hi)/2
            candidate = self.rk4(start, mid)
            if self.geometry(candidate[:4])[0] > capture_radius:
                lo, outside = mid, candidate
            else:
                hi, inside = mid, candidate
            iterations += 1
        require(hi - lo <= width_goal, 'Capture bracket refinement did not meet its bound')
        return {
            'status': 'first_sampled_outside_to_inside_bracket_refined',
            'original_affine_bracket': [str(affine_start), str(affine_start + step)],
            'refined_affine_bracket': [str(affine_start + lo), str(affine_start + hi)],
            'refined_affine_width': str(hi - lo), 'width_goal': str(width_goal),
            'bisection_iterations': iterations,
            'method': 'Fresh single partial RK4 from the fixed original bracket-start state',
            'outside_state': self.summary(outside, capture_radius, initial),
            'inside_state': self.summary(inside, capture_radius, initial),
            'continuous_first_crossing_certified': False,
            'width_is_numerical_bracket_not_a_geodesic_error_bound': True,
        }

    def integrate(self, state, floor, count, capture_radius, initial):
        step = floor / count
        require(step*count == floor, 'Affine substeps do not sum to the exact represented floor')
        distance = self.geometry(state[:4])[0] - capture_radius
        crossing = None
        if distance <= 0:
            crossing = {'status': 'initial_state_at_or_inside_independent_buffered_capture',
                        'refined_affine_bracket': ['0', '0'],
                        'state': self.summary(state, capture_radius, initial)}
        max_h = max_lz = Decimal(0)
        radial_increases = 0
        for index in range(count):
            previous = state
            state = self.rk4(previous, step)
            next_distance = self.geometry(state[:4])[0] - capture_radius
            radial_increases += int(next_distance > distance)
            if crossing is None and distance > 0 and next_distance <= 0:
                crossing = self.refine_event(previous, state, step*index, step,
                                            floor, capture_radius, initial)
            summary = self.summary(state, capture_radius, initial)
            max_h = max(max_h, abs(Decimal(summary['hamiltonian_drift'])))
            max_lz = max(max_lz, abs(Decimal(summary['angular_momentum_drift'])))
            distance = next_distance
        return {
            'substeps': count, 'finite_difference_epsilon': str(self.epsilon),
            'affine_interval': str(floor), 'affine_step': str(step),
            'endpoint': self.summary(state, capture_radius, initial),
            'maximum_sampled_absolute_H_drift': str(max_h),
            'maximum_sampled_absolute_Lz_drift': str(max_lz),
            'sampled_radial_increases': radial_increases,
            'capture_event': crossing if crossing is not None else {
                'status': 'no_sampled_capture_crossing_in_this_affine_interval'},
        }


def compare(first, second):
    a, b = first['endpoint'], second['endpoint']
    result = {key: str(abs(Decimal(a[key]) - Decimal(b[key]))) for key in
              ('hamiltonian', 'norm_2H', 'energy_minus_p_t',
               'axial_angular_momentum', 'radius', 'radius_minus_buffered_capture')}
    for key in ('q', 'p_covariant', 'u_contravariant'):
        result[key + '_maximum_absolute_difference'] = str(max(
            abs(Decimal(x) - Decimal(y)) for x, y in zip(a[key], b[key])))
    ea, eb = first['capture_event'], second['capture_event']
    if (ea['status'] == eb['status'] ==
            'first_sampled_outside_to_inside_bracket_refined'):
        result['inside_event_affine_difference'] = str(abs(
            Decimal(ea['refined_affine_bracket'][1]) - Decimal(eb['refined_affine_bracket'][1])))
        result['inside_event_q_maximum_absolute_difference'] = str(max(
            abs(Decimal(x) - Decimal(y))
            for x, y in zip(ea['inside_state']['q'], eb['inside_state']['q'])))
    else:
        result['event_comparison'] = [ea['status'], eb['status']]
    return result


def main():
    require(not sys.flags.optimize, 'Assertions/normal interpreter mode required')
    require(len(sys.argv) == 3, 'usage: local_reference.py exact-analysis-result.json new-output.json')
    input_path, output_path = map(Path, sys.argv[1:])
    original = input_path.read_bytes()
    require(identity(original) == INPUT_IDENTITY, 'Input is not the exact reviewed Decimal analysis')
    script = Path(__file__).read_bytes()
    analysis = json.loads(original)
    require(analysis['schema'] == 'sirius-kerr-failure-independent-decimal-state-analysis-v1' and
            analysis['analysis_completed'] is True and analysis['source_revision'] == SOURCE and
            analysis['numerical_test_passed'] is False, 'Wrong input analysis authority')
    require(analysis['script_sha256'] ==
            '5c97a89bf27325bfcdcd85af3c1a83ccdff7f211810f3b1fe67b5a7aeb0d5e31',
            'Wrong independent state-analysis source')
    failure = json.loads(analysis['raw_context_json'])
    require(failure['schema'] == 'sirius-trace-continuation-failure-v1' and
            failure['scalar_bytes'] == failure['packet_scalar_bytes'] == 4 and
            failure['record_bytes'] == 272 and failure['byte_order'] == 'little' and
            failure['action'] == 1 and failure['active'] == 64, 'Wrong captured record layout')
    raw = raw_hex(failure['previous_record_hex'], 272)
    current_raw = raw_hex(failure['current_record_hex'], 272)
    packet_raw = raw_hex(failure['packet_hex'], 280)
    packet = struct.unpack('<70f', packet_raw)
    packet_words = struct.unpack('<70I', packet_raw)
    before_control = list(struct.unpack_from('<4I', raw, 240))
    after_control = list(struct.unpack_from('<4I', current_raw, 240))
    expected_identity = failure['expected_identity']
    require(list(struct.unpack_from('<4I', raw, 256)) ==
            list(struct.unpack_from('<4I', current_raw, 256)) == expected_identity and
            expected_identity[2] == 0 and expected_identity[3] == packet_words[69] > 0,
            'Capture identity/epoch mismatch')
    require(before_control == analysis['previous']['control'] and
            after_control == analysis['current']['control'] and
            before_control[:2] == [1, 0] and after_control[:2] == [2, 5] and
            after_control[2] == before_control[2] + 1 and
            after_control[3] == before_control[3], 'Wrong observed failed transition')
    require(packet[2] == packet[5] == packet[6] == packet[46] == 0.0 and packet[3] == 1.0 and
            0 < packet[4] < packet[3] and packet[21] == 2000.0 and packet[68] == 1.0,
            'Reference only represents this captured asymptotically flat Kerr sector')
    work = Work()
    with localcontext() as context:
        context.prec = PRECISION
        context.rounding = 'ROUND_HALF_EVEN'
        values = tuple(exact(value) for value in struct.unpack_from('<8f', raw))
        q, u = values[:4], values[4:]
        mass, spin = exact(packet[3]), exact(packet[4])
        floor = exact(packet[23])
        require(floor == Decimal('0.0199999995529651641845703125'), 'Unexpected exact affine floor')
        require(exact(struct.unpack_from('<f', raw, 192)[0]) == floor,
                'Captured previous adaptive step does not equal the represented floor')
        horizon = mass + (mass*mass - spin*spin).sqrt()
        capture_radius = horizon * exact(packet[26])
        reference = Reference(mass, spin, Decimal(EPSILONS[0]), work)
        p = reference.initial_momentum(q, u)
        state = (*q, *p)
        h0 = reference.hamiltonian(q, p)
        require(h0 != 0, 'Measured off-null norm was lost at the reference precision')
        initial = {'H': h0, 'E': -p[0], 'Lz': q[1]*p[2] - q[2]*p[1]}
        initial_summary = reference.summary(state, capture_radius, initial)
        recovered_u = reference.velocity(q, p)
        runs = {}
        for epsilon in EPSILONS:
            reference = Reference(mass, spin, Decimal(epsilon), work)
            runs[epsilon] = {str(count): reference.integrate(
                state, floor, count, capture_radius, initial) for count in SUBSTEPS}
        convergence = {
            epsilon: {str(low) + '_to_' + str(high): compare(runs[epsilon][str(low)],
                                                            runs[epsilon][str(high)])
                      for low, high in zip(SUBSTEPS, SUBSTEPS[1:])}
            for epsilon in EPSILONS}
        sensitivity = {str(count): compare(runs[EPSILONS[0]][str(count)],
                                          runs[EPSILONS[1]][str(count)]) for count in SUBSTEPS}
        result = {
            'schema': 'sirius-kerr-local-off-null-hamiltonian-reference-v1',
            'analysis_completed': True, 'captured_source_revision': SOURCE,
            'numerical_test_passed': False, 'production_or_shader_modified': False,
            'input': {'path': str(input_path), **identity(original)},
            'script': {'path': str(Path(__file__)), **identity(script)},
            'raw_context_json': analysis['raw_context_json'],
            'decimal_precision': PRECISION, 'decimal_rounding': context.rounding,
            'mass': str(mass), 'spin': str(spin), 'charge': '0', 'Lambda': '0',
            'exact_affine_interval': str(floor), 'independent_outer_horizon_radius': str(horizon),
            'independent_buffered_capture_radius': str(capture_radius),
            'initial': initial_summary, 'initial_exact_decoded_u': [str(value) for value in u],
            'initial_p_to_u_max_absolute_residual': str(max(
                abs(a-b) for a, b in zip(u, recovered_u))),
            'initial_direct_metric_norm_from_state_analysis':
                analysis['calculations']['80']['previous']['gkk'],
            'initial_2H_minus_direct_metric_norm': str(2*h0 -
                Decimal(analysis['calculations']['80']['previous']['gkk'])),
            'mass_shell_projection_performed': False, 'null_only_Carter_formula_used': False,
            'canonical_rhs': 'qdot=dH/dp analytic; pdot_i=central finite difference of -H at fixed p; p_tdot=0',
            'runs': runs, 'substep_refinement_absolute_differences': convergence,
            'finite_difference_epsilon_absolute_sensitivity': sensitivity,
            'work_counts': dict(work.counts),
            'bounds': {'maximum_seconds': MAX_SECONDS, 'maximum_rhs': MAX_RHS,
                       'maximum_geometries': MAX_GEOMETRIES, 'maximum_rk4_steps': MAX_RK4,
                       'maximum_bisections_per_event': MAX_BISECTIONS},
            'elapsed_reference_seconds': time.monotonic() - work.started,
            'scope': 'Independent local Decimal canonical RK4 only, preserving the captured nonzero norm. Fixed CPU substeps refine this one represented affine-floor interval; no product floor or controller changes. First sampled capture brackets use fresh partial integration and are not proofs of the first continuous event. Refinement differences, finite-difference sensitivity, invariant drifts and numerical event brackets are not certified error bounds. No production metric/connection/projector imports, hidden GPU trial candidates, full-ray qualification, or inference identifying the actual rejected GPU conjunct.',
        }
    require(input_path.read_bytes() == original, 'Input changed during reference')
    require(Path(__file__).read_bytes() == script, 'Reference source changed during execution')
    with output_path.open('x', encoding='utf-8') as stream:
        stream.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'analysis_completed': True, 'numerical_test_passed': False,
                      'output': str(output_path), 'work_counts': work.counts}))


if __name__ == '__main__':
    main()
