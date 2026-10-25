"""Keep the complete original C8 reader; add Q's passive route proof only."""
from pathlib import Path
import hashlib

basis = Path.cwd() / 'out/native-cpu-active-prefix/read_result.py'
assert hashlib.sha256(basis.read_bytes()).hexdigest() == '05926a8e76d0786ac344b811c6b170378c79ec7104beb1981549d8c69639e87d'
source = basis.read_text().replace('_try_packet', '_rows_try_packet').replace('providerB', 'providerAA')
before = "words[2].startswith('sirius_native_') and words[2].endswith('_rows_try_packet')"
after = "(words[2].startswith('sirius_native_') and words[2].endswith(('_try_packet', '_admit_packet')) or words[2].startswith('sirius_simd_') and words[2].endswith('_packet') or words[2] in ('sirius_native_simd_supported', 'sirius_native_transport_simd_policy'))"
assert source.count(before) == 1
source = source.replace(before, after)
proof = '''
slugs = ('transport', 'endpoint', 'dopri_phase')
capability, = observation['capability_calls']
symbol_sizes = {}
for line in subprocess.check_output(['/usr/bin/nm', '-S', '-a', action['runtime_provider']], text=True).splitlines():
    words = line.split()
    if len(words) == 4 and words[3] in symbols:
        symbol_sizes[words[3]] = int(words[1], 16)

def checked_entry(row, pure=False):
    assert row['returned'] and row['library'] == action['runtime_provider']
    assert row['thread'] == row['return_thread']
    assert row['entry_sequence'] < row['return_sequence']
    if pure:
        assert row['mxcsr'] == row['return_mxcsr']
    else:
        assert row['mxcsr'] & 0xe040 == row['return_mxcsr'] & 0xe040 == 0
    bases = {int(line.split()[0].split('-')[0], 16) - int(line.split()[2], 16)
             for line in row['provider_maps']}
    base = row['pc'] - symbols[row['function']]
    assert base in bases
    return base

checked_entry(capability, pure=True)
assert isinstance(capability['supported'], bool)
policies = {}
for row in observation['policy_calls']:
    assert row['category'] == 'policy_calls'
    checked_entry(row, pure=True)
    assert capability['supported'] and capability['return_sequence'] < row['entry_sequence']
    assert isinstance(row['selected'], bool)
    assert row['selected'] == (row['workers'] == 0 or row['rows'] > row['workers'])
    assert 1 <= row['rows'] <= 65536 and row['rows'] not in policies
    policies[row['rows']] = row
if capability['supported']:
    assert policies and len({row['workers'] for row in policies.values()}) == 1
else:
    assert not policies

def packet_policy(row):
    if not capability['supported']:
        return False
    if 'transport' not in row['function']:
        return True
    policy = policies[row['grid_x']]
    assert policy['return_sequence'] < row['entry_sequence']
    return policy['selected']

def checked_range(row):
    assert 1 <= row['grid_x'] <= row['capacity'] <= 65536
    assert row['y'] == row['z'] == 0 and row['grid_y'] == row['grid_z'] == 1
    full_rows = row['grid_x'] // 8 * 8
    if packet_policy(row) and row['first_x'] < full_rows:
        assert row['first_x'] % 8 == 0 and row['groups'] == 8
    else:
        assert row['groups'] == 1
    assert row['first_x'] + row['groups'] <= row['grid_x']

assert len(observation['preflight_calls']) == 3
assert {row['function'] for row in observation['preflight_calls']} == {
    'sirius_native_' + slug + '_admit_packet' for slug in slugs}
for row in observation['calls']:
    assert row['category'] == 'calls' and row['handled']
    checked_entry(row)
    checked_range(row)
    assert capability['return_sequence'] < row['entry_sequence']
for row in observation['preflight_calls']:
    assert row['category'] == 'preflight_calls' and row['handled']
    checked_entry(row, pure=True)
    assert row['first_x'] == row['y'] == row['z'] == 0
    assert row['groups'] == row['grid_y'] == row['grid_z'] == 1
    assert 1 <= row['grid_x'] <= row['capacity']
    worker, = [item for item in observation['calls'] if
               item['function'] == row['function'].replace('_admit_packet', '_rows_try_packet')]
    for key in ('pid', 'input', 'input_bytes', 'output', 'output_bytes',
                'capacity', 'grid_x', 'grid_y', 'grid_z'):
        assert row[key] == worker[key], (row['function'], key)
    assert row['return_sequence'] < worker['entry_sequence']

expected_simd_slugs = ({'endpoint', 'dopri_phase'} if capability['supported'] else set())
if capability['supported'] and any(row['selected'] and row['rows'] >= 8 for row in policies.values()):
    expected_simd_slugs.add('transport')
expected_simd = len(expected_simd_slugs)
assert len(observation['simd_calls']) == expected_simd
assert {row['function'] for row in observation['simd_calls']} == {
    'sirius_simd_' + slug + '_packet' for slug in expected_simd_slugs}
assert {row['function'] for row in observation['simd_parent_calls']} == {
    'sirius_native_' + slug + '_rows_try_packet' for slug in expected_simd_slugs}
parents = {row['entry_sequence']: row for row in observation['simd_parent_calls']}
assert len(parents) == len(observation['simd_parent_calls'])
for row in parents.values():
    assert row['category'] == 'simd_parent_calls' and row['handled'] and row['groups'] == 8
    checked_entry(row)
    checked_range(row)
for row in observation['simd_calls']:
    assert row['category'] == 'simd_calls' and row['completed'] and row['groups'] == 8
    base = checked_entry(row)
    parent = parents[row['parent_entry_sequence']]
    assert parent['contains_first_simd']
    assert parent['function'] == row['caller_function'] == row['function'].replace('sirius_simd_', 'sirius_native_').replace('_packet', '_rows_try_packet')
    for key in ('pid', 'thread', 'input', 'output', 'first_x', 'capacity', 'groups'):
        assert row[key] == parent[key], (row['function'], key)
    assert row['input_words'] * 4 == parent['input_bytes']
    assert row['output_words'] * 4 == parent['output_bytes']
    assert parent['entry_sequence'] < row['entry_sequence'] < row['return_sequence'] < parent['return_sequence']
    assert symbols[parent['function']] <= row['caller_pc'] - base < symbols[parent['function']] + symbol_sizes[parent['function']]
    stage = row['function'][len('sirius_simd_'):-len('_packet')]
    input_stride, output_stride, program_words = {
        'transport': (230,4699,15967), 'endpoint': (225,3829,21143),
        'dopri_phase': (1810,4149,10689)}[stage]
    assert row['input_words'] == 1 + row['capacity'] * input_stride + program_words
    assert row['output_words'] == row['capacity'] * output_stride

scalar_parents = {row['entry_sequence']: row for row in observation['scalar_parent_calls']}
assert len(scalar_parents) == len(observation['scalar_parent_calls'])
for row in scalar_parents.values():
    assert row['function'] == 'sirius_native_transport_rows_try_packet' and row['handled']
    checked_entry(row)
    checked_range(row)
seen_scalar_grids = set()
for row in observation['scalar_calls']:
    assert row['category'] == 'scalar_calls' and row['function'] == 'sirius_native_transport_try_packet'
    assert row['handled'] and row['groups'] == 1 and row['grid_x'] not in seen_scalar_grids
    seen_scalar_grids.add(row['grid_x'])
    base = checked_entry(row)
    checked_range(row)
    parent = scalar_parents[row['parent_entry_sequence']]
    assert parent['contains_first_scalar']
    for key in ('pid', 'thread', 'input', 'input_bytes', 'output', 'output_bytes',
                'first_x', 'y', 'z', 'grid_x', 'grid_y', 'grid_z', 'capacity', 'groups'):
        assert row[key] == parent[key], (row['function'], key)
    assert parent['entry_sequence'] < row['entry_sequence']
    if row['entry_rsp'] == parent['entry_rsp']:
        assert row['entry_return_pc'] == parent['entry_return_pc']
        # Tail jump shares the caller return site. Finish callbacks can run in
        # either order; both must observe true after the captured F entry.
        assert row['entry_sequence'] < parent['return_sequence']
    else:
        assert row['entry_rsp'] < parent['entry_rsp']
        assert row['caller_function'] == parent['function']
        for pc in (row['caller_pc'], row['entry_return_pc']):
            assert symbols[parent['function']] <= pc - base < symbols[parent['function']] + symbol_sizes[parent['function']]
        assert row['return_sequence'] < parent['return_sequence']
scalar12, = [row for row in observation['scalar_calls'] if row['grid_x'] == 12]
assert scalar12['input_bytes'] == 4 * (1 + scalar12['capacity'] * 230 + 15967)
assert scalar12['output_bytes'] == scalar12['capacity'] * 4699 * 4
selected_transport12 = 'SIMD8 plus scalar tails' if packet_policy(scalar12) else 'F scalar one-row tasks'
for grid, policy in policies.items():
    if not policy['selected']:
        assert grid in seen_scalar_grids, 'every observed scalar-selected grid needs its first F return'

# Changed R Transport SIMD has fresh component execution and exact new-code joins.
compatibility_namespace = {}
compatibility_source = work / 'compatibility.py'
exec(compile(compatibility_source.read_text(), str(compatibility_source), 'exec'), compatibility_namespace)
transport_component_scope = compatibility_namespace['check'](root)['scope']
'''
before = "endpoint, = [row for row in observation['calls']"
assert source.count(before) == 1
source = source.replace(before, proof + '\n' + before)
before = "assert endpoint['first_x'] + endpoint['groups'] <= 3"
assert source.count(before) == 1
source = source.replace(before, before + "\nassert 0 <= endpoint['first_x'] < 3 and endpoint['groups'] == 1, 'first Endpoint retains three independent scalar tasks'")
source = source.replace("'native_first_packets_handled': 3,", "'native_first_packets_handled': 3, 'pure_preflights_handled': 3, 'first_raw_SIMD_entries_completed': expected_simd, 'runtime_SIMD_supported': capability['supported'], 'observed_transport_policy': observation['policy_calls'], 'observed_Transport12_implementation': selected_transport12, 'transport_component_scope': transport_component_scope,")
source = source.replace('patch changes objectpaths/cacheUUID only',
    'AA propagates family-guarded canonical register counts as compile-time uniform metadata; original Z math/scratch ABI and Y scheduling retained; full9380tuples and generic fallback unchanged')
source = source.replace('three actual first native packets returnedtrue,51markers',
    'three actual first rows packets and pure preflights returnedtrue; policy arguments/returns and actual selected implementation PCs/first returns joined passively;51markers')
exec(compile(source, str(basis), 'exec'))
