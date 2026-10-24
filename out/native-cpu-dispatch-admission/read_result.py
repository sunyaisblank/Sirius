"""Retain full original C8 checks and bind all three pure-admission returns."""
from pathlib import Path
import hashlib

basis = Path.cwd() / 'out/native-cpu-active-prefix/read_result.py'
assert hashlib.sha256(basis.read_bytes()).hexdigest() == '05926a8e76d0786ac344b811c6b170378c79ec7104beb1981549d8c69639e87d'
source = basis.read_text()
proof = '''
assert all(row['groups'] == 1 for row in observation['calls'])
assert len(observation['preflight_calls']) == 3
assert {x['function'] for x in observation['preflight_calls']} == {
    'sirius_native_' + s + '_admit_packet' for s in ('transport', 'endpoint', 'dopri_phase')}
for row in observation['preflight_calls']:
    assert row['is_preflight'] and row['returned'] and row['handled']
    assert row['library'] == action['runtime_provider']
    assert row['thread'] == row['return_thread']
    assert row['mxcsr'] == row['return_mxcsr']
    assert row['first_x'] == row['y'] == row['z'] == 0
    assert row['grid_y'] == row['grid_z'] == row['groups'] == 1
    assert 1 <= row['grid_x'] <= row['capacity']
    bases = {int(line.split()[0].split('-')[0], 16) - int(line.split()[2], 16)
             for line in row['provider_maps']}
    assert row['pc'] - symbols[row['function']] in bases
    worker, = [x for x in observation['calls'] if
               x['function'] == row['function'].replace('_admit_packet', '_try_packet')]
    assert not worker['is_preflight']
    for key in ('pid', 'input', 'input_bytes', 'output', 'output_bytes',
                'capacity', 'grid_x', 'grid_y', 'grid_z', 'groups'):
        assert row[key] == worker[key], (row['function'], key)
'''
changes = [
    ('providerB', 'providerF'),
    ("words[2].endswith('_try_packet')", "words[2].endswith(('_try_packet', '_admit_packet'))"),
    ('endpoint, = [row for row in observation[\'calls\']', proof + '\nendpoint, = [row for row in observation[\'calls\']'),
    ("'native_first_packets_handled': 3,", "'native_first_packets_handled': 3, 'native_pure_preflights_handled': 3,"),
    ('patch changes objectpaths/cacheUUID only',
     'patch changes admission/scheduling/cache UUID, with unchanged property/feature implementation'),
]
for before, after in changes:
    assert before in source, before
    source = source.replace(before, after)
exec(compile(source, str(basis), 'exec'))
