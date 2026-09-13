#!/usr/bin/env python3
"""Independent Decimal metric/differentiation oracle. No product imports/tools."""
from decimal import Decimal as D, localcontext
from pathlib import Path
import hashlib
import json
import struct
import time
import argparse

# These scientific inputs and diagnostic budgets are fixed before device observation.
CASES = [
    ('flat_off_axis', (0, 0, 0, 0, 0.25, 6, 2, 3)),
    ('schwarzschild_near', (1, 0, 0, 0, 0.25, 2, 0.375, 0.5)),
    ('kerr_positive_near', (1, 0.9, 0, 0, 0.25, 2, 0.375, 0.5)),
    ('kerr_negative_near', (1, -0.9, 0, 0, 0.25, 2, 0.375, 0.5)),
    ('kerr_positive_off_axis', (1, 0.9, 0, 0, 0.25, 6, 2, 3)),
    ('charged_kerr_off_axis', (1, 0.7, 0.2, 0, 0.25, 3, -0.75, 1.25)),
    ('schwarzschild_de_sitter', (1, 0, 0, 0.001, 0.25, 6, 2, 3)),
    ('pure_de_sitter', (0, 0, 0, 0.001, 0.25, 2, -0.375, 0.5)),
]
FAMILIES = {'metric': (0, 16), 'inverse': (16, 32), 'connection': (32, 96),
            'metric_partial': (96, 160), 'connection_partial': (160, 416)}
# Absolute + relative*abs(reference), all finite components, including zeros.
BUDGETS = {
    'fp64': {'metric': [2e-13, 2e-13], 'inverse': [2e-13, 2e-13],
             'connection': [2e-12, 2e-12], 'metric_partial': [2e-12, 2e-12],
             'connection_partial': [2e-11, 2e-11]},
    'fp32_and_fp32comp': {'metric': [4e-6, 4e-6], 'inverse': [4e-6, 4e-6],
                         'connection': [3e-5, 3e-5], 'metric_partial': [3e-5, 3e-5],
                         'connection_partial': [2e-4, 2e-4]},
}


def exact_inputs(raw):
    # Float multiplication is not used to construct physical squared parameters.
    parameters = [D.from_float(struct.unpack('<f', struct.pack('<f', v))[0]) for v in raw[:4]]
    return parameters + [D.from_float(float(v)) for v in raw[4:]]


def metric(parameters, position):
    mass, spin, charge, cosmological = parameters
    _, x, y, z = position
    signs = [-1, 1, 1, 1]
    eta = [[D(signs[i] if i == j else 0) for j in range(4)] for i in range(4)]
    if all(p == 0 for p in parameters):
        return [v for row in eta for v in row]
    square = spin * spin
    R2 = x*x + y*y + z*z
    # Independent unscaled positive algebraic root, at the declared finite fixtures.
    r2 = ((R2-square) + ((R2-square)**2 + 4*square*z*z).sqrt()) / 2
    if r2 <= 0:
        raise ValueError('reference fixture outside unique positive radius domain')
    radius = r2.sqrt()
    covector = [D(1), (radius*x+spin*y)/(r2+square),
                (radius*y-spin*x)/(r2+square), z/radius]
    H = (2*mass*radius-charge*charge)/(r2+square*z*z/r2)
    if cosmological:
        if spin != 0 or charge != 0:
            raise ValueError('unsupported rotating/charged cosmological configuration')
        H += cosmological*r2/3
    return [eta[i][j] + H*covector[i]*covector[j] for i in range(4) for j in range(4)]


def inverse(matrix):
    # General pivoted matrix inversion, deliberately not the rank-one identity.
    rows = [[matrix[4*i+j] for j in range(4)] + [D(i == j) for j in range(4)] for i in range(4)]
    for col in range(4):
        pivot = max(range(col, 4), key=lambda row: abs(rows[row][col]))
        rows[col], rows[pivot] = rows[pivot], rows[col]
        scale = rows[col][col]
        if scale == 0:
            raise ValueError('singular reference metric')
        rows[col] = [v/scale for v in rows[col]]
        for row in range(4):
            if row != col:
                scale = rows[row][col]
                rows[row] = [a-scale*b for a, b in zip(rows[row], rows[col])]
    return [rows[i][j+4] for i in range(4) for j in range(4)]


def partial(evaluate, point, column, step):
    # Fourth-order central finite difference, including independent t evaluations.
    samples = []
    for offset in [-2, -1, 1, 2]:
        p = list(point)
        p[column] += offset*step
        samples.append(evaluate(tuple(p)))
    return [(a-8*b+8*c-d)/(12*step) for a, b, c, d in zip(*samples)]


def calculate(raw, precision, step_string):
    with localcontext() as context:
        context.prec = precision
        inputs = exact_inputs(raw)
        parameters, point = tuple(inputs[:4]), tuple(inputs[4:])
        step = D(step_string)
        cache = {}
        def met(p):
            if p not in cache:
                cache[p] = metric(parameters, p)
            return cache[p]
        def geometry(p):
            g = met(p)
            inv = inverse(g)
            dg = [partial(met, p, c, step) for c in range(4)]
            gamma = [sum(inv[mu*4+s]*(dg[nu][s*4+rho]+dg[rho][s*4+nu]-dg[s][nu*4+rho])
                         for s in range(4))/2
                     for mu in range(4) for nu in range(4) for rho in range(4)]
            return g, inv, dg, gamma
        g, inv, dg, gamma = geometry(point)
        dgamma = [partial(lambda p: geometry(p)[3], point, c, step) for c in range(4)]
        values = g + inv + gamma + [v for row in dg for v in row] + [v for row in dgamma for v in row]
        assert len(values) == 416 and all(v.is_finite() for v in values)
        identity = max(abs(sum(g[4*i+k]*inv[4*k+j] for k in range(4))-D(i == j))
                       for i in range(4) for j in range(4))
        # Levi-Civita symmetry and metric compatibility use independent finite dg.
        compatibility = max(abs(dg[c][4*i+j]-sum(g[4*k+j]*gamma[16*k+4*c+i]+
                            g[4*i+k]*gamma[16*k+4*c+j] for k in range(4)))
                            for c in range(4) for i in range(4) for j in range(4))
        return values, {'matrix_inverse_residual': str(identity),
                        'metric_compatibility_residual': str(compatibility),
                        'unique_metric_samples': len(cache)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, required=True)
    ROOT = parser.parse_args().output_dir
    ROOT.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    schedule = [(100, '1e-8'), (100, '5e-9'), (100, '2.5e-9'),
                (80, '2.5e-9'), (120, '1.25e-9')]
    rows = []
    fixtures = []
    for name, raw in CASES:
        evaluations = [calculate(raw, precision, step) for precision, step in schedule]
        final = evaluations[-1][0]
        families = {}
        for family, (lo, hi) in FAMILIES.items():
            def difference(a, b):
                with localcontext() as context:
                    context.prec = 120
                    return max(abs(x-y) for x, y in zip(evaluations[a][0][lo:hi], evaluations[b][0][lo:hi]))
            d01, d12, d24 = difference(0, 1), difference(1, 2), difference(2, 4)
            precision_effect = difference(2, 3)
            assert d24 < D('1e-24'), (name, family, d24)
            assert precision_effect < D('1e-55'), (name, family, precision_effect)
            families[family] = {'step_difference_h_h2': str(d01), 'step_difference_h2_h4': str(d12),
                                'step_difference_h4_h8': str(d24),
                                'first_difference_ratio': str(d01/d12) if d12 else None,
                                'second_difference_ratio': str(d12/d24) if d24 else None,
                                'precision80_vs100_at_h4': str(precision_effect),
                                'final_to_binary64_max_error': str(max(abs(D.from_float(float(v))-v) for v in final[lo:hi]))}
        inputs = exact_inputs(raw)
        rows.append({'name': name, 'inputs_decimal': [str(v) for v in inputs],
                     'inputs_binary64_hex': [float(v).hex() for v in inputs],
                     'families': families, 'algebraic_checks': [v[1] for v in evaluations],
                     'reference_decimal': [str(v) for v in final],
                     'reference_binary64_hex': [float(v).hex() for v in final]})
        fixtures.append((name, inputs, final))
        print(name, 'complete', flush=True)
    header = ['#pragma once', '// Generated by generate_metric_consistency_reference.py; see its refinement checks.',
              '#include <array>', '#include <string_view>', '', 'namespace metric_consistency_reference {',
              'struct Case {', '    std::string_view name;', '    std::array<double, 8> input;',
              '    std::array<double, 416> values;', '};',
              'inline const std::array<Case, 8>& Cases() {',
              '    static const std::array<Case, 8> cases{{']
    for name, inputs, values in fixtures:
        header += ['        Case{"'+name+'",', '            {'+', '.join(float(v).hex() for v in inputs)+'},', '            {']
        for i in range(0, 416, 4):
            header.append('                '+', '.join(float(v).hex() for v in values[i:i+4])+',')
        header += ['            }},']
    header += ['    }};', '    return cases;', '}', '}  // namespace metric_consistency_reference', '']
    (ROOT/'metric_consistency_reference.h').write_text('\n'.join(header))
    report = {'scope': 'Independent Decimal defining-metric and numerical-derivative reference; not device observation or formal interval certificate.',
              'schedule': schedule, 'family_offsets': FAMILIES, 'diagnostic_comparison_budgets': BUDGETS,
              'reference_acceptance': {'finest_step_change_below': '1e-24', 'precision_change_below': '1e-55'},
              'rows': rows, 'wall_seconds': time.monotonic()-started,
              'generator_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'header_sha256': hashlib.sha256((ROOT/'metric_consistency_reference.h').read_bytes()).hexdigest()}
    (ROOT/'reference.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({'success': True, 'cases': len(rows), 'wall_seconds': report['wall_seconds'],
                      'header_sha256': report['header_sha256']}, indent=2))

if __name__ == '__main__':
    main()
