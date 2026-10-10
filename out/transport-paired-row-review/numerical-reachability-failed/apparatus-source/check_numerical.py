"""Check actual full buffers and unchanged independent finite witnesses."""
from fractions import Fraction
from pathlib import Path
import hashlib
import json
import math
import struct

TASK = Path(__file__).resolve().parent
GATE = TASK / 'numerical'
ROW_WORDS = 4699
CAPACITY = 24
EPSILON = Fraction(1, 1 << 52)


def read(name):
    return json.loads((GATE / name).read_text(encoding='utf-8-sig'))


def sha(data):
    return hashlib.sha256(data).hexdigest()


def binary32(word):
    value = struct.unpack('<f', struct.pack('<I', word))[0]
    assert math.isfinite(value), f'nonfinite retained word {word:08x}'
    return Fraction.from_float(value)


def represented(words):
    high, low, tail, radius = map(binary32, words[:4])
    assert words[4] == 1
    assert abs(high) <= 1 << 120 and 0 <= radius <= 1 << 120
    high_spacing = binary32((words[0] & 0x7fffffff) + 1) - abs(high)
    low_spacing = binary32((words[1] & 0x7fffffff) + 1) - abs(low)
    assert abs(low) <= high_spacing and abs(tail) <= low_spacing
    assert high != 0 or (low == 0 and tail == 0)
    return high + low + tail, radius, high


fixtures = {'general': read('reference_cases.json')['cases'],
            'schwarzschild': read('schwarzschild_reference_cases.json')['cases']}
sequences = read('sequences.json')
result, bridge, owner = read('query.json'), read('bridge-query.json'), read('run-query/owner.json')
assert bridge['completed'] and bridge['returncode'] == 0 and bridge['all_selected_inputs_unchanged']
assert not bridge['outer_stop'] and not bridge['bridge_error'] and not bridge['cleanup_error']
assert bridge['linux_bridge_child_birth_absent']
assert owner['status'] == 'completed' and owner['screen_execution_completed']
assert owner['returncode'] == 0 and not owner['outer_stop'] and owner['owned_process_absent']
assert not owner['cleanup_errors']
assert result['completed'] and result['capacity'] == CAPACITY
assert result['completed_dispatches'] == len(sequences) * 2 == 54
assert len(result['samples']) == 54
assert result['observed_words_per_dispatch'] == CAPACITY * ROW_WORDS
assert result['module_sha256'] == 'b52e6373f6b30aac6c7fe877f6e07f3c47c7ac9eafa1671840eacd19960d7cb1'
assert result['alternative_module_sha256'] == '721be7625f8eebbf722adf3e79a89050bffe9bf9b7260b8f5120b0a85ccd7f27'
assert result['shader_float64_enabled'] and result['shader_fma32_enabled']
assert result['actual_resident_bytes'] == sum(a['requirement_bytes'] for a in result['allocations'])
assert result['actual_resident_bytes'] <= 8 * 1024 * 1024
assert len(result['allocations']) == 2
samples = {(s['sequence'], s['variant']): s for s in result['samples']}
assert len(samples) == len(result['samples'])
shadow = bytearray((GATE / 'initial.bin').read_bytes())
records = []
fields = curved_rhs_fields = refused_rows = late_rows = valid_rows = 0
represented_rhs_fields = 0
coverage = {family: set() for family in fixtures}
for sequence in sequences:
    name, active = sequence['name'], sequence['active_rows']
    update = (GATE / 'inputs' / (name + '.bin')).read_bytes()
    assert sha(update) == sequence['input_sha256'] and len(update) == sequence['upload_bytes']
    shadow[:len(update)] = update
    variants = []
    for variant in ('baseline', 'candidate'):
        sample = samples[(name, variant)]
        path = GATE / 'observations' / ('output.' + name + '.' + variant + '.bin')
        data = path.read_bytes()
        assert sha(data) == sample['output_sha256']
        assert len(data) == CAPACITY * ROW_WORDS * 4
        assert sample['input_shadow_sha256'] == sha(shadow)
        assert sample['whole_capacity_word_mismatches'] == sample['inactive_word_changes'] == 0
        assert sample['active_rows'] == active and sample['uploaded_bytes'] == len(update)
        assert sample['groups_x'] == (active if variant == 'baseline' else (active + 1) // 2)
        assert sample['groups_y'] == (2 if variant == 'candidate' and active % 2 else 1)
        assert data[active * ROW_WORDS * 4:] == b'\xff' * ((CAPACITY - active) * ROW_WORDS * 4)
        variants.append(data)
    assert variants[0] == variants[1], name
    words = struct.unpack('<' + str(CAPACITY * ROW_WORDS) + 'I', variants[1])
    row_records = []
    for index, row in enumerate(sequence['rows']):
        at = index * ROW_WORDS
        output = words[at:at + ROW_WORDS]
        valid, stages, rhs = output[:3]
        assert valid in (0, 1) and 0 <= stages <= 7 and rhs in (0, 1)
        assert output[3] == 0 and not any(output[2404:]), (name, index, 'reserved scratch')
        expected = row['expected']
        if expected != 'valid':
            assert valid == rhs == 0, (name, index, 'refusal completion')
            if expected == 'late-refused':
                assert stages >= 2, (name, index, 'late-refusal path not reached')
                late_rows += 1
            elif expected == 'refused':
                assert stages == 0, (name, index, 'early input refusal')
                assert not any(output), (name, index, 'early refusal not complete')
            else:
                assert expected == 'table-refused' and stages <= 1
            refused_rows += 1
            row_records.append({'name': row['name'], 'expected': expected, 'stages': stages})
            continue
        fixture = fixtures[row['family']][row['index']]
        flat = fixture['name'] in ('analytic-flat', 'frozen-analytic-flat')
        assert valid == 1 and stages == 7 and rhs == int(not flat)
        count = 160 if row['family'] == 'general' or flat else 440
        reference = list(map(Fraction, fixture['reference']))
        gap = list(map(Fraction, fixture['precision_gap']))
        errors = []
        for field in range(count):
            offset = 4 + field * 5 if field < 160 else 1004 + (field - 160) * 5
            center, radius, high = represented(output[offset:offset + 5])
            difference = abs(center - reference[field])
            reference_high = Fraction.from_float(float(fixture['reference'][field]))
            rounding = 128 * EPSILON * EPSILON * (abs(reference_high) + abs(high))
            assert difference <= radius + gap[field] + rounding, (name, index, field, 'independent enclosure')
            errors.append(difference)
        if not flat:
            # The production decoder validates every returned RHS value, even
            # when the frozen general oracle covers only the final 160 fields.
            for field in range(280):
                offset = 1004 + field * 5
                represented(output[offset:offset + 5])
            represented_rhs_fields += 280
        for begin in range(0, count, 4):
            indices = [i - 40 if 120 <= i < 160 else i for i in range(begin, begin + 4)]
            scale = max(abs(reference[i]) for i in indices)
            assert max(errors[begin:begin + 4]) <= Fraction(1, 10**11) * scale, (name, index, begin, 'original relative accuracy')
        if flat:
            assert not any(output[1004:2404])
        valid_rows += 1
        fields += count
        curved_rhs_fields += count - 160
        coverage[row['family']].add(row['index'])
        row_records.append({'name': row['name'], 'expected': expected, 'stages': stages,
                            'independent_fields': count, 'rhs_valid': bool(rhs)})
    records.append({'sequence': name, 'active_rows': active, 'whole_capacity_equal': True,
                    'inactive_canary_unchanged': True, 'actual_input_shadow_equal': True,
                    'rows': row_records})
assert coverage['general'] == set(range(15)) and coverage['schwarzschild'] == set(range(23))
assert late_rows > 0
receipt = {'pass': True, 'capacity': CAPACITY, 'sequences': len(sequences), 'dispatches': 54,
           'whole_capacity_words_per_pair': CAPACITY * ROW_WORDS,
           'valid_rows': valid_rows, 'refused_rows': refused_rows, 'late_refused_rows': late_rows,
           'independent_fields': fields, 'independent_curved_rhs_fields': curved_rhs_fields,
           'represented_curved_rhs_fields': represented_rhs_fields,
           'all_original_general_fixtures': 15, 'all_original_schwarzschild_fixtures': 23,
           'actual_resident_bytes': result['actual_resident_bytes'], 'records': records,
           'scope': 'Exact-rational decoding against unchanged finite independent witnesses, original1e-11 group accuracy and enclosure/gap/binary64-conversion allowance. Same-provider full-capacity raw/actual-input/inactive-canary identity. No timing, full integration/stream/rollback, original frame, full scientific/release or adoption verdict.'}
(GATE / 'reference-check.json').write_text(json.dumps(receipt, indent=2) + '\n')
print(json.dumps({k: v for k, v in receipt.items() if k not in ('records', 'scope')}, indent=2))
