"""Check complete native sidecars against the reviewed layout and original test receipts."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import collections, hashlib, json, struct, xml.etree.ElementTree as ET

WORK = Path(__file__).resolve().parent
assert len(sys.argv) == 2 and sys.argv[1] in ('device-timestamps', 'mixed-transport', 'factored-endpoints')
case = sys.argv[1]
directory = WORK / case
def read(path): return json.loads(path.read_text(encoding='utf-8-sig'))
plan = read(WORK / 'execution-plan.json')
bridge = read(directory / 'native-control-bridge.json')
retirement = read(directory / 'retirement-run.json')
assert bridge['accepted'] and retirement['accepted']
assert read(directory / 'native-retirement.json')['pass']
owner = read(directory / 'run/owner.json')
assert owner['numerical_control_passed'] and owner['source_revision'] == plan['source_revision']
xml = ET.parse(directory / 'run/gtest.xml').getroot()
tests = list(xml.iter('testcase'))
assert len(tests) == 1 and tests[0].get('classname') + '.' + tests[0].get('name') == plan['cases'][case]
assert not list(xml.iter('failure')) and not list(xml.iter('error')) and not list(xml.iter('skipped'))
properties = {p.get('name'): p.get('value') for p in tests[0].iter('property')}
capture_properties = {k: json.loads(v) for k, v in properties.items() if k.startswith('shader_clock_capture_')}
paths = sorted((directory / 'run/shader-clock-observations').glob('*.json'))
assert paths and len(paths) == len(capture_properties)
assert len(list((directory / 'run/shader-clock-observations').glob('*.bin'))) == len(paths)
rows = []
ordinals = set()
xml_ordinals = {int(name.removeprefix('shader_clock_capture_')) for name in capture_properties}
for path in paths:
    metadata = read(path)
    ordinal = metadata['ordinal']
    assert type(ordinal) is int and ordinal >= 0 and ordinal not in ordinals
    ordinals.add(ordinal)
    assert type(metadata['paired']) is bool
    assert metadata == capture_properties['shader_clock_capture_' + str(ordinal)]
    payload = path.with_suffix('.bin').read_bytes()
    assert len(payload) == metadata['bytes'] == metadata['capacity'] * 36 * 4
    assert hashlib.sha256(payload).hexdigest() == metadata['sha256']
    assert metadata['capacity'] in (24, 32) and 0 < metadata['active_rows'] <= metadata['capacity']
    assert metadata['groups_y'] == metadata['groups_z'] == 1
    assert metadata['slot'] in (0, 1) and (metadata['paired'] or metadata['slot'] == 0)
    assert metadata['sidecar'] not in (metadata['input'], metadata['output'])
    assert 0 < metadata['sidecar_allocation_bytes'] <= 8 * 1024 * 1024
    stride = 36 if metadata['stage'] == 'transport' else 12
    assert metadata['stage'] in ('transport', 'endpoint') and metadata['words_per_row'] == stride
    slots = (stride - 4) // 2
    words = struct.unpack('<%dI' % (len(payload) // 4), payload)
    assert not any(words[metadata['capacity'] * stride:]), 'Endpoint maximum-buffer padding was not cleared'
    assert not any(words[metadata['active_rows'] * stride:metadata['capacity'] * stride]), 'Inactive rows retained earlier observations'
    masks = collections.Counter()
    completed = 0
    observations = []
    for row in range(metadata['active_rows']):
        record = words[row * stride:(row + 1) * stride]
        magic, observed_row, mask, completion = record[:4]
        masks[mask] += 1
        assert completion in (0, 1)
        if mask == 0:
            assert not any(record), 'Unobserved bypass/refusal row contains stale words'
            continue
        assert magic == 0x53434c4b and observed_row == row and mask < (1 << slots)
        assert mask & (mask + 1) == 0, 'Observed cuts are not a reached execution prefix'
        if metadata['stage'] == 'endpoint' and mask == (1 << slots) - 1:
            assert completion == 1, 'Successful Endpoint root publication lacks completion'
        if completion:
            assert mask == (1 << slots) - 1
            completed += 1
        ticks = []
        for slot in range(slots):
            tick = record[4 + 2 * slot] | record[5 + 2 * slot] << 32
            if mask & (1 << slot): ticks.append(tick)
            else: assert tick == 0, 'Unreached cut retained a clock value'
        observations.append({'row': row, 'mask': mask, 'completion': completion, 'raw_ticks': ticks,
            'adjacent_modulo64_differences': [(b - a) % (1 << 64) for a, b in zip(ticks, ticks[1:])]})
    rows.append({'capture': metadata, 'mask_counts': dict(masks), 'completed_rows': completed, 'observed_rows': observations})
assert ordinals == xml_ordinals, 'JSON and XML capture ordinal sets differ'
stages = {row['capture']['stage'] for row in rows}
assert stages == {'device-timestamps': {'transport', 'endpoint'}, 'mixed-transport': {'transport'}, 'factored-endpoints': {'endpoint'}}[case]
assert any(row['completed_rows'] for row in rows)
result = {'accepted': True, 'source_revision': plan['source_revision'], 'case': plan['cases'][case],
    'original_test_passed': True, 'captures': len(rows), 'stages': sorted(stages),
    'records': sorted(rows, key=lambda row: row['capture']['ordinal']),
    'scope': 'Complete diagnostic sidecars and all assertions of the named original test passed. All-zero rows are unobserved; sidecar bytes alone do not identify flat bypass versus prior refusal. Raw ticks and modulo64 differences only; clock rate and exclusive costs unestablished; perturbed observers; no ordinary rendering benefit or full science/release claim.'}
(directory / 'observation-analysis.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps({k: result[k] for k in ('accepted', 'case', 'captures', 'stages')}))
