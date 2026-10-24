"""Retain the entire original C8 reader; join changed invocation metadata."""
from pathlib import Path
import hashlib

basis = Path.cwd() / 'out/native-cpu-active-prefix/read_result.py'
assert hashlib.sha256(basis.read_bytes()).hexdigest() == '05926a8e76d0786ac344b811c6b170378c79ec7104beb1981549d8c69639e87d'
source = basis.read_text().replace('_try_packet', '_rows_try_packet').replace('providerB', 'providerG')
before = "words[2].startswith('sirius_native_') and words[2].endswith('_rows_try_packet')"
after = "(words[2].startswith('sirius_native_') and words[2].endswith(('_rows_try_packet', '_admit_packet')) or words[2].startswith('sirius_simd_') and words[2].endswith('_packet'))"
assert source.count(before) == 1
source = source.replace(before, after)
proof = '''
slugs = ('transport', 'endpoint', 'dopri_phase')
assert len(observation['preflight_calls']) == len(observation['simd_calls']) == 3
assert {row['function'] for row in observation['preflight_calls']} == {
    'sirius_native_' + slug + '_admit_packet' for slug in slugs}
assert {row['function'] for row in observation['simd_calls']} == {
    'sirius_simd_' + slug + '_packet' for slug in slugs}
assert len(observation['simd_parent_calls']) >= 3
assert {row['function'] for row in observation['simd_parent_calls']} == {
    'sirius_native_' + slug + '_rows_try_packet' for slug in slugs}
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
        assert row['handled'] and row['mxcsr'] == row['return_mxcsr']
    else:
        assert row['mxcsr'] & 0xe040 == row['return_mxcsr'] & 0xe040 == 0
    bases = {int(line.split()[0].split('-')[0], 16) - int(line.split()[2], 16)
             for line in row['provider_maps']}
    base = row['pc'] - symbols[row['function']]
    assert base in bases
    return base

def checked_range(row):
    assert 1 <= row['grid_x'] <= row['capacity'] <= 65536
    assert row['y'] == row['z'] == 0 and row['grid_y'] == row['grid_z'] == 1
    assert row['first_x'] % 8 == 0
    assert row['groups'] == min(8, row['grid_x'] - row['first_x'])
    assert row['first_x'] + row['groups'] <= row['grid_x']

for row in observation['calls']:
    assert row['category'] == 'calls'
    checked_entry(row)
    checked_range(row)
for row in observation['preflight_calls']:
    assert row['category'] == 'preflight_calls'
    checked_entry(row, pure=True)
    assert row['first_x'] == row['y'] == row['z'] == 0
    assert row['groups'] == row['grid_y'] == row['grid_z'] == 1
    assert 1 <= row['grid_x'] <= row['capacity']
    worker, = [item for item in observation['calls'] if
               item['function'] == row['function'].replace('_admit_packet', '_rows_try_packet')]
    # The pure predicate checks one initial group. Worker packing preserves
    # the original grid/capacity but can start at another concurrently queued packet.
    for key in ('pid', 'input', 'input_bytes', 'output', 'output_bytes',
                'capacity', 'grid_x', 'grid_y', 'grid_z'):
        assert row[key] == worker[key], (row['function'], key)
    assert row['return_sequence'] < worker['entry_sequence']
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
'''
before = "endpoint, = [row for row in observation['calls']"
assert source.count(before) == 1
source = source.replace(before, proof + '\n' + before)
before = "assert endpoint['first_x'] + endpoint['groups'] <= 3"
assert source.count(before) == 1
source = source.replace(before, before + "\nassert endpoint['first_x'] == 0 and endpoint['groups'] == 3, 'first Endpoint must retain scalar tail3'")
source = source.replace("'native_first_packets_handled': 3,", "'native_first_packets_handled': 3, 'pure_preflights_handled': 3, 'first_raw_SIMD_entries_completed': 3,")
source = source.replace('patch changes objectpaths/cacheUUID only',
    'patch changes packed native scheduling/entry routing/cache UUID with unchanged property/feature implementation')
source = source.replace('three actual first native packets returnedtrue,51markers',
    'three actual first rows packets and pure preflights returnedtrue; first three actual rawSIMD entries completed with exact concurrent parent/word/range joins;51markers')
exec(compile(source, str(basis), 'exec'))
