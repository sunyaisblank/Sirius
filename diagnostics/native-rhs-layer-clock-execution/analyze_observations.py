"""Validate exact captures and bind reached layer ranges to actual input tables."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import collections, hashlib, json, statistics, struct, xml.etree.ElementTree as ET

WORK = Path(__file__).resolve().parent
assert len(sys.argv) == 2 and sys.argv[1] in ('device-timestamps', 'mixed-transport', 'factored-endpoints')
case = sys.argv[1]
directory = WORK / case
MOD = 1 << 64
def read(path): return json.loads(path.read_text(encoding='utf-8-sig'))
def seal(path):
    data = path.read_bytes()
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}
def unpack(data): return struct.unpack('<%dI' % (len(data) // 4), data)
def ticks(words, count): return [words[4 + 2*i] | words[5 + 2*i] << 32 for i in range(count)]
def summary(values):
    return {'count': len(values), 'minimum': min(values), 'median': statistics.median(values), 'maximum': max(values)}

plan = read(WORK / 'execution-plan.json')
assert read(directory / 'native-control-bridge.json')['accepted']
assert read(directory / 'retirement-run.json')['accepted'] and read(directory / 'native-retirement.json')['pass']
owner = read(directory / 'run/owner.json')
assert owner['numerical_control_passed'] and owner['source_revision'] == plan['source_revision']
assert read(WORK / 'placement/placement-query-bridge.json')['accepted']
assert read(WORK / 'placement/retirement-run.json')['accepted']
placement = read(WORK / 'placement-analysis.json')['stages']['Transport']
assert placement['all_printed_instruction_bytes_match_ELF_prefix']
assert placement['module'] == read(WORK / 'compiled-review.json')['module']
xml = ET.parse(directory / 'run/gtest.xml').getroot()
tests = list(xml.iter('testcase'))
assert len(tests) == 1 and tests[0].get('classname') + '.' + tests[0].get('name') == plan['cases'][case]
assert not list(xml.iter('failure')) and not list(xml.iter('error')) and not list(xml.iter('skipped'))
property_nodes = list(tests[0].iter('property'))
assert len({p.get('name') for p in property_nodes}) == len(property_nodes)
properties = {p.get('name'): p.get('value') for p in property_nodes}
capture_properties = {k: json.loads(v) for k,v in properties.items() if k.startswith('shader_clock_capture_')}
capture_directory = directory / 'run/shader-clock-observations'
paths = sorted(capture_directory.glob('*.json'))
assert paths and len(paths) == len(capture_properties)
assert {p.stem for p in paths} == {p.stem for p in capture_directory.glob('*.bin') if not p.name.endswith('.input.bin')}
ordinals = set()
records = []
input_files = set()
table_records = {}
layer_samples = collections.defaultdict(lambda: {'body': [], 'join': [], 'span': []})
phase_bands = collections.defaultdict(list)

def tables_from_input(metadata, path):
    source = path.with_suffix('.input.bin')
    input_files.add(source.name)
    data = source.read_bytes()
    assert len(data) == metadata['input_bytes'] == (1 + metadata['capacity']*230 + 15967)*4
    assert hashlib.sha256(data).hexdigest() == metadata['input_sha256']
    words = unpack(data)
    capacity = words[0]
    assert capacity == metadata['stored_capacity']
    # Preserve rejected layouts instead of strengthening the shader's packet admission.
    base = 1 + capacity*230
    if capacity == 0 or capacity > 65536 or base + 15967 != len(words):
        return {}, {'stored_capacity': capacity, 'layout_available': False}
    tables = {}
    for name, offset, length, count, layers in (('general', 0, 13096, 2587, 117), ('schwarzschild', 13096, 2871, 557, 42)):
        program = words[base + offset:base + offset + length]
        key = name + ':' + hashlib.sha256(struct.pack('<%dI' % len(program), *program)).hexdigest()
        shape = program[:3] == (count, 459, 40)
        table = {'key': key, 'family': name, 'shape_valid': shape, 'header': list(program[:3])}
        if shape:
            layer_base = 43 + 5*count
            offsets = program[layer_base:layer_base + layers + 1]
            assert len(offsets) == layers + 1
            table.update(count=count, layers=layers, offsets=list(offsets),
                         roots=list(program[3:43]), nodes=[list(program[43+5*i:48+5*i]) for i in range(count)])
            table['range_plan_valid'] = offsets[0] == 0 and offsets[-1] == count and all(0 < b-a <= 64 for a,b in zip(offsets, offsets[1:]))
        if key in table_records: assert table_records[key] == table
        else: table_records[key] = table
        if shape: tables[name] = table
    return tables, {'stored_capacity': capacity, 'layout_available': True, 'tables': {n:t['key'] for n,t in tables.items()}}

for path in paths:
    metadata = read(path)
    ordinal = metadata['ordinal']
    assert type(ordinal) is int and ordinal >= 0 and ordinal not in ordinals
    ordinals.add(ordinal)
    assert type(metadata['paired']) is bool
    assert metadata == capture_properties['shader_clock_capture_' + str(ordinal)]
    payload = path.with_suffix('.bin').read_bytes()
    assert len(payload) == metadata['bytes'] == metadata['capacity']*8226*4
    assert hashlib.sha256(payload).hexdigest() == metadata['sha256']
    assert metadata['capacity'] in (24, 32) and 0 < metadata['active_rows'] <= metadata['capacity']
    assert metadata['groups_y'] == metadata['groups_z'] == 1
    assert metadata['slot'] in (0,1) and (metadata['paired'] or metadata['slot'] == 0)
    assert metadata['sidecar'] not in (metadata['input'], metadata['output'])
    assert 0 < metadata['sidecar_allocation_bytes'] <= 8*1024*1024
    assert metadata['stage'] in ('transport', 'endpoint')
    transport = metadata['stage'] == 'transport'
    stride = 8226 if transport else 12
    assert metadata['words_per_row'] == stride
    if transport:
        tables, input_record = tables_from_input(metadata, path)
    else:
        assert metadata['input_bytes'] == metadata['stored_capacity'] == 0 and metadata['input_sha256'] == ''
        assert not path.with_suffix('.input.bin').exists()
        tables, input_record = {}, None
    words = unpack(payload)
    assert not any(words[metadata['capacity']*stride:]), 'Endpoint maximum-buffer padding was not cleared'
    assert not any(words[metadata['active_rows']*stride:metadata['capacity']*stride]), 'Inactive sidecar rows retained earlier observations'
    masks = collections.Counter()
    completed = 0
    observations = []
    for row in range(metadata['active_rows']):
        record = words[row*stride:(row+1)*stride]
        magic, observed_row, mask, completion = record[:4]
        masks[mask] += 1
        assert completion in (0,1)
        if mask == 0:
            assert not any(record), 'Unobserved row contains stale words'
            continue
        slots = 16 if transport else 4
        assert magic == 0x53434c4b and observed_row == row and 0 < mask < (1 << slots)
        assert mask & (mask+1) == 0, 'Phase cuts do not form a reached prefix'
        if not transport and mask == (1 << slots)-1: assert completion == 1
        if completion:
            assert mask == (1 << slots)-1
            completed += 1
        phase_ticks = ticks(record, slots)
        for slot,tick in enumerate(phase_ticks):
            if not mask & (1 << slot): assert tick == 0
        observation = {'row': row, 'mask': mask, 'completion': completion,
                       'raw_ticks': phase_ticks[:mask.bit_count()],
                       'adjacent_modulo64_differences': [(b-a)%MOD for a,b in zip(phase_ticks, phase_ticks[1:mask.bit_count()])]}
        if transport:
            candidates = set(tables)
            reached = []
            for stage in range(7):
                entries = []
                gap = False
                for layer in range(117):
                    entry = record[36+(stage*117+layer)*10:46+(stage*117+layer)*10]
                    start,end,layer_mask,reserved = entry[:4]
                    assert reserved == 0
                    if layer_mask == 0:
                        assert not any(entry)
                        gap = True
                        continue
                    # Both body and join cuts precede every ordinary layer refusal.
                    assert layer_mask == 7 and not gap, 'Incomplete layer observation in a completed synchronous command'
                    assert mask & (1 << (1+2*stage)), 'Layer precedes its reached stage-input cut'
                    values = ticks(entry, 3)
                    candidates = {n for n in candidates if layer < tables[n]['layers'] and tables[n]['offsets'][layer:layer+2] == [start,end]}
                    entries.append([start,end,*values])
                published = bool(mask & (1 << (2+2*stage)))
                if published:
                    assert entries, 'Published RHS stage has no reached layers'
                    candidates = {n for n in candidates if len(entries) == tables[n]['layers'] and tables[n]['range_plan_valid']}
                reached.append({'stage': stage, 'published': published, 'entries': entries})
            if any(s['entries'] for s in reached): assert candidates, 'Observed layer ranges match no captured table family'
            else: candidates = set()
            observation['family_candidates'] = sorted(candidates)
            observation['layer_stages'] = reached
            # Rank complete published RHS stages only. Refused entries remain raw records.
            for stage_record in reached:
                if not stage_record['published']: continue
                stage = stage_record['stage']
                a,b = phase_ticks[1+2*stage],phase_ticks[2+2*stage]
                band = (b-a)%MOD
                phase_bands[','.join(sorted(candidates))].append(band)
                if len(candidates) != 1: continue
                previous_join = 0
                for layer,entry in enumerate(stage_record['entries']):
                    start,end,t0,t1,t2 = entry
                    # This checks observed modular nesting, without converting ticks to time.
                    relative = [(t-a)%MOD for t in (t0,t1,t2)]
                    assert previous_join <= relative[0] <= relative[1] <= relative[2] <= band
                    previous_join = relative[2]
                    for name in candidates:
                        table = tables[name]
                        assert 0 <= start < end <= table['count']
                        composition = tuple(sorted(collections.Counter(n[0] for n in table['nodes'][start:end]).items()))
                        key = (table['key'],layer,start,end,composition)
                        sample = layer_samples[key]
                        sample['body'].append((t1-t0)%MOD)
                        sample['join'].append((t2-t1)%MOD)
                        sample['span'].append((t2-t0)%MOD)
        observations.append(observation)
    records.append({'capture': metadata, 'input_layout': input_record, 'mask_counts': dict(masks),
                    'completed_rows': completed, 'observed_rows': observations})

assert ordinals == {int(n.removeprefix('shader_clock_capture_')) for n in capture_properties}
assert input_files == {p.name for p in capture_directory.glob('*.input.bin')}
stages = {r['capture']['stage'] for r in records}
assert stages == {'device-timestamps': {'transport','endpoint'}, 'mixed-transport': {'transport'}, 'factored-endpoints': {'endpoint'}}[case]
assert any(r['completed_rows'] for r in records)
layer_summary = []
for (key,layer,start,end,composition),samples in layer_samples.items():
    layer_summary.append({'table': key, 'layer': layer, 'start': start, 'end': end,
        'opcode_composition': dict(composition), 'lane_zero_before_to_after_body': summary(samples['body']),
        'lane_zero_after_body_to_after_original_join': summary(samples['join']),
        'whole_observed_layer': summary(samples['span'])})
result = {'accepted': True, 'source_revision': plan['source_revision'], 'case': plan['cases'][case],
          'original_test_passed': True, 'captures': len(records), 'stages': sorted(stages),
          'layer_entry_columns': ['start','end','before_body_tick','after_body_tick','after_join_tick'],
          'records': sorted(records, key=lambda r:r['capture']['ordinal']), 'actual_tables': table_records,
          'published_rhs_bands': {k:summary(v) for k,v in phase_bands.items()}, 'published_layer_summary': layer_summary,
          'scope': 'Original full test assertions pass; complete sidecars/input snapshots hash/XML bound. Layer masks7 denote reached joins, including refusal; publication is separate. Unobserved/ambiguous families stay unclassified. Actual recorded opcode mixtures are not assigned proportional costs. Raw modulo64 tick intervals only, inclusive of scheduling/waits/observer traffic, perturbed; rate/exclusive arithmetic or barrier cost/ordinary rendering benefit/full science or release unestablished.'}
target = directory / 'observation-analysis.json'
assert not target.exists()
target.write_text(json.dumps(result, separators=(',',':')) + '\n')
print(json.dumps({k:result[k] for k in ('accepted','case','captures','stages')}))
