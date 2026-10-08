"""Reproduce the frozen finite Schwarzschild RK witnesses. Requires mpmath.

Run explicitly with Python; ordinary builds consume the checked-in header.
The unchanged reference.py uses the defining metric, generic inverse and mp.diff.
These witnesses are evaluated at exact decoded three-limb centres at 75/105
decimal digits; their agreement is a reference refinement check, not a certified
input-ball or trajectory enclosure. No production graph or private helper is imported.
"""
import importlib.util
import json
import struct
import sys
from pathlib import Path

sys.dont_write_bytecode = True
W = Path(__file__).resolve().parent


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


oracle = load('independent_rk', W / 'reference.py')
mp = oracle.mp


def decode(packet):
    return [sum((mp.mpf(struct.unpack('<f', struct.pack('<I', word))[0])
                 for word in packet[offset:offset + 3]), mp.mpf(0))
            for offset in range(0, len(packet), 5)]


def make_input(name, chart):
    if name in ('null_weak', 'public_mass_min', 'public_mass_max'):
        mass = {'null_weak': mp.mpf(1), 'public_mass_min': mp.mpf('.1'),
                'public_mass_max': mp.mpf(100)}[name]
        oracle.ref.params = [mass, mp.mpf(0), mp.mpf(0), mp.mpf(0)]
        oracle.ref.reflection = [chart, 1, chart, 1]
        row = {'x': [0, 50 * mass, 2 * mass, 10 * mass],
               'k': [-1, mp.mpf('.2'), mp.mpf('.3'), -mp.mpf('.7')],
               'columns': [
                   {'X': [0, 0, 0, 0], 'V': [0, mp.mpf('.003'), 0, 0]},
                   {'X': [0, 0, 0, 0], 'V': [0, 0, mp.mpf('.004'), 0]},
                   {'X': [0, mass, 0, 0], 'V': [0, 0, 0, 0]},
                   {'X': [0, 0, mass, 0], 'V': [0, 0, 0, 0]}]}
        phase = oracle.ref.initialize(oracle.ref.project(oracle.ref.initialize(row)))
        return [*oracle.ref.params, *phase, mp.mpf(chart), mass / 2]
    mass = mp.mpf(1)
    positions = {'weak': ['50', '2', '10'], 'strong': ['2.2', '.3', '-.4'],
                 'axis': ['0', '0', '3.3'], 'interior': ['.125', '.1875', '-.0625']}
    if name in positions:
        position = list(map(mp.mpf, positions[name]))
    elif name in ('scale_below', 'scale_above'):
        position = [mp.mpf(4) + (1 if name == 'scale_above' else -1) * mp.power(2, -40),
                    mp.mpf(0), mp.mpf(0)]
    elif name in ('small_mass', 'large_mass'):
        mass = mp.power(2, -20 if name == 'small_mass' else 20)
        position = [mass * v for v in (mp.mpf(7) / 4, mp.mpf(1) / 2, mp.mpf(3) / 8)]
    else:
        raise ValueError(name)
    phase = [mp.mpf(1) / 8, *position, mp.mpf(7) / 8,
             mp.mpf(1) / 4, -mp.mpf(3) / 8, mp.mpf(1) / 2]
    columns = [[mp.mpf(node ** field) / mp.power(2, 8 + field)
                for field in range(8)] for node in (1, 2, 3, 4)]
    assert mp.det(mp.matrix([[columns[c][f] for c in range(4)] for f in range(4)])) == mp.mpf(12) / mp.power(2, 38)
    phase.extend(value for column in columns for value in column)
    assert len(phase) == 40
    return [mass, mp.mpf(0), mp.mpf(0), mp.mpf(0), *phase,
            mp.mpf(chart), mass / (8192 if name == 'interior' else 128)]


def calculate(packet):
    inputs = decode(packet)
    oracle.ref.params = inputs[:4]
    oracle.ref.reflection = [inputs[44], 1, inputs[44], 1]
    rhs_records = []
    original_rhs = oracle.ref.rhs

    def record_rhs(phase):
        result = original_rhs(phase)
        rhs_records.extend(result)
        return result

    oracle.ref.rhs = record_rhs
    try:
        final = oracle.integrate(mp.matrix(inputs[4:44]), inputs[45])
    finally:
        oracle.ref.rhs = original_rhs
    assert len(final) == 160 and len(rhs_records) == 280
    return final + rhs_records


def main():
    cases = []
    names = ('weak', 'strong', 'axis', 'interior', 'scale_below', 'scale_above',
             'small_mass', 'large_mass', 'null_weak', 'public_mass_min', 'public_mass_max')
    for name in names:
        for chart in (-1, 1):
            packets, samples = [], []
            for precision in (75, 105):
                with mp.workdps(precision):
                    packet = [word for value in make_input(name, chart) for word in oracle.pair(value)]
                    packets.append(packet)
                    samples.append(calculate(packet))
            assert packets[0] == packets[1]
            with mp.workdps(105):
                gaps = [abs(a - b) for a, b in zip(*samples)]
                assert max(gaps) < mp.mpf('1e-55')
                cases.append({'name': name + ('-outgoing' if chart == -1 else '-ingoing'),
                              'input': packets[1],
                              'reference': [mp.nstr(v, 105) for v in samples[1]],
                              'precision_gap': [mp.nstr(v, 105) for v in gaps]})
            print(json.dumps({'case': cases[-1]['name'], 'passed': True}), flush=True)
    frozen = json.loads((W / 'reference_cases.json').read_text())
    flat = next(case for case in frozen['cases'] if case['name'] == 'analytic-flat')
    cases.append({'name': 'frozen-analytic-flat', 'input': flat['input'],
                  'reference': flat['reference'] + ['0'] * 280,
                  'precision_gap': flat['precision_gap'] + ['0'] * 280})
    document = {'classification': 'Finite independent arithmetic RK witnesses; no trajectory, interval certificate, device or performance acceptance.',
                'reference_order': 'fifth[40],fourth[40],increment[40],error[40],seven RHS[40]',
                'precisions': [75, 105], 'curved_cases': 22, 'cases': cases,
                'input_reference_binding': 'RK is evaluated at exact decoded three-term packet centres, not at pre-packing decimal values.',
                'limits': ['Observed precision agreement is not a certified reference enclosure.',
                           'Input-ball and trajectory/truncation certificates are not established.',
                           'Eight arbitrary-covector families use four independent coupled columns; three physical families use independent metric null projection with angular and pupil columns.',
                           'Tiny/large raw masses are raw arithmetic controls, beyond typed physical mass range.']}
    (W / 'schwarzschild_reference_cases.json').write_text(json.dumps(document, indent=2, sort_keys=True) + '\n')
    lines = ['// Frozen defining-metric RK witnesses; regenerate with schwarzschild_reference.py.',
             '// Curved cases retain all 160 finals and 280 RHS fields; the flat case retains its original oracle.',
             '#pragma once', '#include <array>', '#include <cstdint>',
             '// clang-format off', 'namespace sirius::test::retained_schwarzschild {', 'struct Pair {double high,low;};',
             'struct Case {const char* name;std::array<std::uint32_t,230> input;std::array<Pair,440> reference;std::array<double,440> gap;};',
             f'inline constexpr std::array<Case,{len(cases)}> cases{{{{']
    with mp.workdps(105):
        for case in cases:
            lines.append('{' + json.dumps(case['name']) + ',{{' + ','.join(str(v) + 'u' for v in case['input']) + '}},{{')
            for text in case['reference']:
                value = mp.mpf(text)
                high = float(value)
                low = float(value - mp.mpf(high))
                lines.append('{' + high.hex() + ',' + low.hex() + '},')
            lines.append('}},{{' + ','.join(float(v).hex() for v in case['precision_gap']) + '}}},')
    lines.extend(['}};', '}  // namespace sirius::test::retained_schwarzschild', '// clang-format on'])
    (W / 'schwarzschild_reference_cases.h').write_text('\n'.join(lines) + '\n')
    print(json.dumps({'curved_cases': 22, 'flat_cases': 1, 'reference_values': 22 * 440 + 160}), flush=True)


if __name__ == '__main__':
    main()
