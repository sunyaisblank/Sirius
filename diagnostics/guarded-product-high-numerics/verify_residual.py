"""Verify every actual residual word against the immutable independent corpus."""
from pathlib import Path
import hashlib
import json
import struct
import sys

sys.dont_write_bytecode = True
root = Path.cwd().resolve()
work = Path(__file__).resolve().parent
stage = root / 'bin/linux-gcc/guarded-product-high-numerics'
mode = sys.argv[1]
assert mode in ['core', 'rte32']
folder = work / f'gpu-residual-candidate-{mode}'
owner = json.loads((folder / 'owner.json').read_text())
assert owner['passed'] and owner['status'] == 'terminal' and owner['child_exit'] == 0
assert not owner['cleanup_errors'] and not owner['remaining_owned_processes']
assert owner['whole_inputs_unchanged'] and owner['source_unchanged']
assert all(record['absent'] for record in owner['observed_births'])
fixture = json.loads((work / 'residual-fixture.json').read_text())
inputs = (work / 'residual-input.bin').read_bytes()
expected = (work / 'residual-expected.bin').read_bytes()
actual = (folder / 'actual.bin').read_bytes()
for name, raw in [('residual-input.bin', inputs), ('residual-expected.bin', expected)]:
    record = fixture['files'][name]
    assert len(raw) == record['bytes'] and hashlib.sha256(raw).hexdigest() == record['sha256']
assert struct.unpack_from('<I', inputs)[0] == fixture['records'] == 53528
assert len(actual) == len(expected) == 4 * fixture['complete_output_words']
assert actual == expected, 'Actual residual changed a complete independent result word'
rows = list(struct.iter_unpack('<4I', actual))
assert all(row[2] in (0, 1) and row[3] == 0 for row in rows)
assert len(rows) == fixture['records']
logs = [json.loads(line) for line in (folder / 'stdout.log').read_text().splitlines() if line.startswith('{')]
assert any(row.get('driver_id') == 13 and row.get('driver') == 'llvmpipe' and row.get('RTE32') for row in logs)
done = next(row for row in logs if row.get('completed_dispatches') == 1)
assert done['cases'] == fixture['records'] and done['observed_words'] == fixture['complete_output_words']
assert done['mismatches'] == done['untouched_complement_words'] == 0
assert done['requested_bytes'] == done['actual_resident_bytes'] == fixture['requested_device_bytes'] <= 8388608
compile_record = json.loads((work / f'compiled-residual-candidate-{mode}.json').read_text())
module = root / compile_record['actual_module']['path']
assert hashlib.sha256(module.read_bytes()).hexdigest() == compile_record['actual_module']['sha256']
assert compile_record['complete_actual_module_validated'] and compile_record['private_guarded_fmul_allowance']
report = {'passed': True, 'mode': mode, 'records': len(rows), 'original_ordered_records': 53524,
    'signed_binade_witnesses': 4, 'complete_words': fixture['complete_output_words'],
    'complete_independent_result_words_equal': True, 'untouched_complement_words': 0,
    'actual_module': compile_record['actual_module'], 'output_bytes': len(actual),
    'output_sha256': hashlib.sha256(actual).hexdigest(), 'actual_device_observation': done,
    'hardware_estimate_input_present': False,
    'scope': 'Actual shader residual on the pinned RTE32-capable software provider. Real FMul estimate only; does not exercise both possible non-RTE neighbor estimates, prove all platforms, or satisfy connected/integrated/performance/science/release qualification.'}
target = work / f'residual-check-candidate-{mode}.json'
assert not target.exists()
target.write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps(report))
