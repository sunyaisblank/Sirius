"""Check independent quotient intervals, then full original/candidate words."""
from pathlib import Path
import hashlib
import importlib.util
import json
import struct
import sys

sys.dont_write_bytecode = True
root = Path.cwd()
work = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('independent_reference', root / 'tests/support/portable_binary32/reference.py')
ref = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ref)
fixture = json.loads((work / 'triple-fixture.json').read_text())
assert len(sys.argv) == 3
producer, mode = sys.argv[1:]
assert producer in ('baseline', 'candidate') and mode in ('integer', 'normal', 'normal_transport')
folder = work / ('gpu-triple-' + producer + '-' + mode)
data = (folder / 'actual.bin').read_bytes()
assert len(data) == fixture['complete_words_each'] * 4
rows = list(struct.iter_unpack('<8I', data))
assert len(rows) == fixture['records']
BOUND = 0x7b800000


def transform(a, b):
    # Pre-existing ordered FastTwoSum contract, each step rounded on the
    # independent exact Fraction lattice; no production word algorithm.
    if a & ref.MASK < b & ref.MASK:
        a, b = b, a
    am, bm = a & ref.MASK, b & ref.MASK
    if mode == 'normal_transport' and bm == 0 and am < ref.INF:
        return (a if am else a & b), b
    high = ref.add(a, b)
    displaced = ref.add(high, a ^ ref.SIGN)
    low = ref.add(b, displaced ^ ref.SIGN)
    if mode != 'integer' and 0x0c000000 <= am <= 0x7e800000 and 0x0c000000 <= bm <= 0x7e800000:
        high = 0 if high & ref.MASK == 0 else high
        low = 0 if low & ref.MASK == 0 else low
    return high, low


def unit_pack(a):
    hi, lo, tail, error = a
    x, y = transform(lo, tail)
    z, w = transform(hi, x)
    x, carry = transform(w, y)
    high, w = transform(z, x)
    low, tail = transform(w, carry)
    assert ref.exact(high) + ref.exact(low) + ref.exact(tail) == sum(ref.exact(v) for v in a[:3])
    return (high, low, tail, error, 1, 0, 0, 0)


valid, refusals, known_valid, known_refusals, unit_identities = 0, 0, 0, 0, 0
for index, (case, row) in enumerate(zip(fixture['cases'], rows, strict=True)):
    assert row[4] in (0, 1) and row[5:] == (0, 0, 0), (index, 'status/padding')
    required = case['required_status']
    if required >= 0:
        assert row[4] == required, (index, case['family'], 'required admission/refusal', row)
        known_valid += required
        known_refusals += 1 - required
    if row[4] == 0:
        refusals += 1
        assert row == (0,) * 8, (index, 'incomplete refusal')
        continue
    valid += 1
    assert all(v & ref.MASK <= BOUND for v in row[:4]) and not row[3] & ref.SIGN
    a, b = case['a'], case['b']
    assert all(v & ref.MASK <= BOUND for v in (*a, *b))
    assert not (a[3] | b[3]) & ref.SIGN
    ac = sum(ref.exact(v) for v in a[:3])
    bc = sum(ref.exact(v) for v in b[:3])
    ar, br = ref.exact(a[3]), ref.exact(b[3])
    al, ah, bl, bh = ac - ar, ac + ar, bc - br, bc + br
    assert not bl <= 0 <= bh, (index, 'admitted denominator interval contains zero')
    corners = [x / y for x in (al, ah) for y in (bl, bh)]
    center = sum(ref.exact(v) for v in row[:3])
    radius = ref.exact(row[3])
    assert center - radius <= min(corners) and max(corners) <= center + radius, (
        index, case['family'], 'independent quotient interval not enclosed', row)
    if case['contract'] == 'unit':
        assert row == unit_pack(a), (index, 'unit denominator complete pack identity')
        unit_identities += 1
    elif case['contract'] == 'zero':
        assert center == radius == 0, (index, 'exact zero numerator')
    elif case['contract'] == 'self':
        assert center - radius <= 1 <= center + radius

assert known_valid > 0 and known_refusals > 0 and valid > 0 and refusals > 0
markers = [rows[i] for i in fixture['recovery_indices']]
assert len(markers) > 1 and all(row == markers[0] for row in markers)
if producer == 'candidate':
    baseline = work / ('gpu-triple-baseline-' + mode) / 'actual.bin'
    assert data == baseline.read_bytes(), 'Candidate changed a complete original result word'
owner = json.loads((folder / 'owner.json').read_text())
assert owner['passed'] and owner['status'] == 'terminal' and not owner['cleanup_errors'] and not owner['remaining_owned_processes']
logs = [json.loads(line) for line in (folder / 'stdout.log').read_text().splitlines() if line.startswith('{')]
assert any(r.get('driver_id') == 13 and r.get('driver') == 'llvmpipe' for r in logs)
done = next(r for r in logs if r.get('completed_dispatches') == 1)
assert done['actual_resident_bytes'] == fixture['requested_device_bytes'] <= 8388608
assert done['untouched_complement_words'] == done['malformed_wire_words'] == 0
assert done['observed_words'] == fixture['complete_words_each']
assert done['capture_only'] == (producer == 'baseline')
assert done['mismatches'] == (None if producer == 'baseline' else 0)
report = {'passed': True, 'producer': producer, 'mode': mode, 'records': len(rows),
          'complete_words': fixture['complete_words_each'], 'valid_results': valid, 'complete_zero_refusals': refusals,
          'required_valid_cases': known_valid, 'required_refusals': known_refusals,
          'unit_complete_pack_identities': unit_identities, 'recovery_markers_complete_words_exact': len(markers),
          'every_valid_result_encloses_independent_exact_input_quotient_interval': True,
          'complete_original_candidate_words_equal': producer == 'candidate',
          'output_bytes': len(data), 'output_sha256': hashlib.sha256(data).hexdigest(),
          'scope': 'Finite full-division/refusal/recovery preservation plus independent mathematical enclosures. '
                   'Original capture is a word comparison target, not an independent oracle; integrated/performance/full qualification remain separate.'}
(work / ('triple-check-' + producer + '-' + mode + '.json')).write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps(report))
