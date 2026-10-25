"""Original sole C8, with passive selected-route metadata checks only."""
from pathlib import Path
import hashlib

basis = Path.cwd() / 'out/native-cpu-active-prefix/run_case.py'
assert hashlib.sha256(basis.read_bytes()).hexdigest() == 'a5b98807d45afbc56435e553367e8baed9d7ce6628b0d1c4735713f0494748a9'
source = basis.read_text().replace('_try_packet', '_rows_try_packet')
before = "assert (folder / 'gtest.xml').is_file()"
after = """capability, = record['capability_calls']
assert capability['returned'] and isinstance(capability['supported'], bool)
assert all(row['returned'] and row['selected'] == (row['workers'] == 0 or row['rows'] > row['workers'])
           for row in record['policy_calls'])
expected_slugs = {'endpoint', 'dopri_phase'} if capability['supported'] else set()
if capability['supported'] and any(row['selected'] and row['rows'] >= 8 for row in record['policy_calls']):
    expected_slugs.add('transport')
assert len(record['preflight_calls']) == 3 and len(record['simd_calls']) == len(expected_slugs)
assert {row['function'] for row in record['preflight_calls']} == {
    'sirius_native_' + slug + '_admit_packet' for slug in ('transport', 'endpoint', 'dopri_phase')}
assert {row['function'] for row in record['simd_calls']} == {
    'sirius_simd_' + slug + '_packet' for slug in expected_slugs}
assert all(row['returned'] and row['handled'] for row in record['preflight_calls'])
assert all(row['returned'] and row['completed'] for row in record['simd_calls'])
assert all(row['returned'] and row['handled'] for row in record['scalar_calls'])
assert len([row for row in record['scalar_calls'] if row['grid_x'] == 12]) == 1
assert (folder / 'gtest.xml').is_file()"""
assert source.count(before) == 1
source = source.replace(before, after)
exec(compile(source, str(basis), 'exec'))
