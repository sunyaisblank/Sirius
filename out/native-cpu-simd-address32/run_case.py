"""Execute the unchanged sole C8 with passive first-call metadata capture."""
from pathlib import Path
import hashlib

basis = Path.cwd() / 'out/native-cpu-active-prefix/run_case.py'
assert hashlib.sha256(basis.read_bytes()).hexdigest() == 'a5b98807d45afbc56435e553367e8baed9d7ce6628b0d1c4735713f0494748a9'
source = basis.read_text().replace('_try_packet', '_rows_try_packet')
before = "assert (folder / 'gtest.xml').is_file()"
after = """assert len(record['preflight_calls']) == len(record['simd_calls']) == 3
assert {row['function'] for row in record['preflight_calls']} == {
    'sirius_native_' + slug + '_admit_packet' for slug in ('transport', 'endpoint', 'dopri_phase')}
assert {row['function'] for row in record['simd_calls']} == {
    'sirius_simd_' + slug + '_packet' for slug in ('transport', 'endpoint', 'dopri_phase')}
assert all(row['returned'] and row['handled'] for row in record['preflight_calls'])
assert all(row['returned'] and row['completed'] for row in record['simd_calls'])
assert (folder / 'gtest.xml').is_file()"""
assert source.count(before) == 1
source = source.replace(before, after)
exec(compile(source, str(basis), 'exec'))
