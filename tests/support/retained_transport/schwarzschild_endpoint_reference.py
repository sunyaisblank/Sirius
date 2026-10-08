"""Additional defining-metric Schwarzschild Endpoint witnesses, at 75/105 digits.

This adds only missing signed-axis, horizon, signed-mass and scale/domain
boundaries to the unchanged factor29/original17 authorities. No specialized
H=2M/r expression, production graph or observed device result is imported.
"""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import struct
import sys

sys.dont_write_bytecode = True
W = Path(__file__).resolve().parent
R = W.parents[2]
F = R / 'tests/support/retained_transport'


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def identity(path):
    data = path.read_bytes()
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}


def specifications(mp):
    cases = []
    for axis in range(3):
        for sign in (-1, 1):
            position = [mp.mpf(0)] * 3
            position[axis] = mp.mpf(4 * sign)
            cases.append((f'axis-{axis + 1}-' + ('negative' if sign < 0 else 'positive'),
                          mp.mpf(1), position, True, (-1, 1)))
    for sign in (-1, 1):
        cases.append(('horizon-' + ('below' if sign < 0 else 'above'), mp.mpf(1),
                      [2 + sign * mp.power(2, -35), mp.mpf(0), mp.mpf(0)], True, (-1, 1)))
    # The ingoing exact linear horizon is already in factor29/original17.
    cases.append(('exact-linear-horizon', mp.mpf(1),
                  [mp.mpf(2), mp.mpf(0), mp.mpf(0)], True, (-1,)))
    cases.append(('negative-mass', mp.mpf(-1), list(map(mp.mpf, (8, 4, 2))), True, (-1, 1)))
    for sign in (-1, 1):
        cases.append(('scale-4-' + ('below' if sign < 0 else 'above'), mp.mpf(1),
                      [4 + sign * mp.power(2, -40), mp.mpf(0), mp.mpf(0)], True, (-1, 1)))
    # These exact dyadic points bracket the largest r whose square's high limb
    # is at the original 2^60 product operand bound. A finite real reference
    # does not require the retained baseline or candidate to admit the row.
    for name, radius in [('below', mp.power(2, 30) - mp.power(2, 6)),
                         ('at', mp.power(2, 30))]:
        cases.append(('prefix-domain-' + name, mp.mpf(1),
                      [radius, mp.mpf(0), mp.mpf(0)], False, (-1, 1)))
    return cases


def make_packet(endpoint, specification, chart):
    mp, ref = endpoint.mp, endpoint.ref
    _, mass, position, _, _ = specification
    ref.params = [mass, mp.mpf(0), mp.mpf(0), mp.mpf(0)]
    ref.reflection = [chart, 1, chart, 1]
    columns = [[mp.mpf(node ** field) / mp.power(2, 8 + field)
                for field in range(8)] for node in (1, 2, 3, 4)]
    assert mp.det(mp.matrix([[columns[c][f] for c in range(4)] for f in range(4)])) == mp.mpf(12) / mp.power(2, 38)
    phase = [mp.mpf(1) / 8, *position, mp.mpf(7) / 8, mp.mpf(1) / 4,
             -mp.mpf(3) / 8, mp.mpf(1) / 2,
             *(value for column in columns for value in column)]
    _, projected = endpoint.project(mp.matrix(phase))
    assert len(projected) == 80
    return [word for value in [*ref.params, *projected[:40], mp.mpf(chart)]
            for word in endpoint.transport.pair(value)]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--generate', action='store_true', required=True)
    parser.add_argument('--output-dir', type=Path, default=W)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    paths = [args.output_dir / ('schwarzschild_endpoint_reference.' + suffix) for suffix in ('json', 'h')]
    assert not any(path.exists() for path in paths), 'Refusing to overwrite frozen references'
    authorities = [Path(__file__).resolve(), F / 'endpoint_factor_reference.py',
                   F / 'endpoint_factor_reference.json', F / 'endpoint_factor_reference.h',
                   F / 'endpoint_reference.py', F / 'endpoint_reference.json', F / 'endpoint_reference.h',
                   F / 'reference.py', F.parent / 'cpu_critical/reference.py',
                   F.parent / 'retained_camera/reference.py', F / 'schwarzschild_reference_cases.json']
    bindings = {str(path.relative_to(R)): identity(path) for path in authorities}
    factor = load('independent_factor', F / 'endpoint_factor_reference.py')
    endpoint = load('independent_endpoint', F / 'endpoint_reference.py')
    mp = endpoint.mp
    cases = []
    # Specifications must be reconstructed under each working precision,
    # including exact dyadic offsets; do not construct them at the default dps.
    for index, specification in enumerate(specifications(mp)):
        name, _, _, required, charts = specification
        for chart in charts:
            label = name + ('-outgoing' if chart == -1 else '-ingoing')
            case = factor.freeze_case(
                endpoint, label,
                lambda i=index, c=chart: make_packet(endpoint, specifications(mp)[i], c),
                {'authority': str(Path(__file__).resolve().relative_to(R)),
                 'profile': name, 'chart': chart,
                 'construction': 'generic defining-metric null projection before three-limb packing'})
            case['admission_required'] = required
            case['domain'] = ('finite raw arithmetic witness' if required else
                              'representation/admission-boundary probe; check all100 roots only if admitted')
            cases.append(case)
            print(json.dumps({'case': label, 'admission_required': required, 'passed': True}), flush=True)
    assert len(cases) == 27 and sum(c['admission_required'] for c in cases) == 23
    # None of these packets replaces an existing factor29 witness.
    old = json.loads((F / 'endpoint_factor_reference.json').read_text())['cases']
    assert not {c['input_sha256'] for c in cases} & {c['input_sha256'] for c in old}
    for case in cases:
        if case['name'].startswith('prefix-domain-'):
            values = factor.decode(mp, case['input'])
            target = mp.power(2, 30) - mp.power(2, 6) if '-below-' in case['name'] else mp.power(2, 30)
            assert values[5] == target and values[6] == values[7] == 0
            assert values[1] == values[2] == values[3] == 0
            assert case['input'][5 * 5 + 1:5 * 5 + 4] == [0, 0, 0]
    assert bindings == {str(path.relative_to(R)): identity(path) for path in authorities}
    document = {
        'classification': 'additional finite independent exact-decoded-centre Schwarzschild Endpoint witnesses',
        'precisions': [75, 105], 'reference_order': 'metric[16],unprojected tangent[4],projected phase[40],physical[40]',
        'case_count': len(cases), 'required_admission_cases': 23, 'admission_boundary_probes': 4,
        'input_words': 225, 'roots_per_alternative': 100, 'input_bindings': bindings,
        'source_and_inputs_unchanged': True, 'mpmath_version': mp.__version__, 'cases': cases,
        'reused_unchanged_authorities': ['factor29 retains weak/public-mass/tiny-radius families and ingoing horizon/spatial fallback',
                                       'original17 retains flat and Kerr continuation criteria'],
        'limits': ['Observed75/105 agreement below unchanged1e-55 is not a certified reference enclosure.',
                   'Negative mass and prefix-domain cases are raw representation witnesses, not typed physical-mass acceptance.',
                   'Finite real boundary references do not require baseline or candidate numerical admission.',
                   'Every admitted row must still meet original1e-29 enclosure slack and1e-11 accuracy criteria.',
                   'No device, trajectory, timing, native hardware, full-frame or release acceptance.']}
    text = factor.header(cases, mp).replace('retained_endpoint_factor', 'retained_schwarzschild_endpoint')
    text = text.replace('See endpoint_factor_reference.py.', 'See schwarzschild_endpoint_reference.py and its unchanged generic authorities.')
    text = text.replace('std::array<Expected,2> expected;};', 'std::array<Expected,2> expected;bool admission_required;};')
    lines = text.splitlines()
    ends = [i for i, line in enumerate(lines) if line == '}}},']
    assert len(ends) == len(cases)
    for index, line in enumerate(ends):
        lines[line] = '}},' + ('true' if cases[index]['admission_required'] else 'false') + '},'
    paths[0].write_text(json.dumps(document, indent=2, sort_keys=True) + '\n')
    paths[1].write_text('\n'.join(lines) + '\n')
    print(json.dumps({'cases': len(cases), 'required_admission_cases': 23, 'admission_boundary_probes': 4,
                      'root_alternatives': sum(len(c['alternatives']) for c in cases), 'device_calls': 0,
                      'outputs': {str(path): identity(path) for path in paths}}), flush=True)


if __name__ == '__main__':
    main()
