"""Validate division captures and bind selected records to actual complete input tables."""
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
assert placement['module'] == read(WORK / 'compiler.json')['product']['spv']
assert read(WORK / 'actual-placement-review.json')['accepted']
preparation = read(WORK / 'preparation.json')
selections = preparation['selections']
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
division_samples = collections.defaultdict(lambda: collections.defaultdict(list))
terminal_counts = collections.Counter()
unclassified_records = 0
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
                for slot in range(12):
                    entry = record[36+(stage*12+slot)*56:92+(stage*12+slot)*56]
                    node,count,reached_mask,reason,valid,entry_magic = entry[:6]
                    if reached_mask == 0:
                        assert not any(entry), 'Unobserved division contains stale metadata/operands/ticks'
                        continue
                    assert entry_magic == 0x4456434b and 0 < reached_mask < (1 << 20)
                    assert reached_mask & (reached_mask+1) == 0, 'Division cuts are not a reached prefix'
                    assert mask & (1 << (1+2*stage)), 'Division precedes reached stage input'
                    expected = [name for name,selection in selections.items()
                                if selection['program_words'] == count and slot < len(selection['nodes'])
                                and selection['nodes'][slot]['node'] == node]
                    assert len(expected) == 1, 'Metadata does not join the actual compiled slot selector'
                    # A returned validity word and a full tick mask are separate observations.
                    if reason in (2,3,4): assert reached_mask == 1
                    elif reason == 5: assert reached_mask in tuple((1 << (4+5*i))-1 for i in range(3))
                    elif reason == 6: assert reached_mask in tuple((1 << (7+5*i))-1 for i in range(3))
                    elif reason == 7: assert reached_mask == (1 << 19)-1
                    elif reason in (1,8): assert reached_mask == (1 << 20)-1
                    else: raise AssertionError(('Unknown terminal reason',reason))
                    if reason in (2,3,5,6,7,8): assert valid == 0
                    else: assert valid != 0
                    operands = [list(entry[6:11]),list(entry[11:16])]
                    if reason == 2: assert operands[0][4] == 0 or operands[1][4] == 0
                    else: assert operands[0][4] != 0 and operands[1][4] != 0
                    if reason == 4: assert valid == operands[0][4]
                    values = [entry[16+2*i] | entry[17+2*i] << 32 for i in range(20)]
                    for cut,tick in enumerate(values):
                        if not reached_mask & (1 << cut): assert tick == 0
                    values = values[:reached_mask.bit_count()]
                    terminal_counts[(reason,valid)] += 1
                    matches = []
                    for name,table in tables.items():
                        if table['count'] != count or not (node < len(table['nodes']) and table['nodes'][node][0] == 5): continue
                        layers = [i for i,(start,end) in enumerate(zip(table['offsets'],table['offsets'][1:])) if start <= node < end]
                        matches.append({'family':name,'table':table['key'],'node_words':table['nodes'][node],
                                        'actual_layers':layers,'range_plan_valid':table['range_plan_valid']})
                    # Modified/refused input tables are retained, never relabelled as canonical.
                    if matches: candidates &= {m['family'] for m in matches}
                    else: candidates.clear()
                    entry_record = {'slot':slot,'node':node,'program_count':count,'mask':reached_mask,
                                    'terminal_reason':reason,'returned_valid_word':valid,'a_words':operands[0],
                                    'b_words':operands[1],'raw_ticks':values,'actual_table_matches':matches,
                                    'adjacent_modulo64_differences':[(y-x)%MOD for x,y in zip(values,values[1:])]}
                    entries.append(entry_record)
                published = bool(mask & (1 << (2+2*stage)))
                reached.append({'stage':stage,'published':published,'entries':entries})
            assert not any(record[4740:]), 'Unused suffix of Transport row was not cleared'
            observation['family_candidates'] = sorted(candidates) if any(s['entries'] for s in reached) else []
            observation['division_stages'] = reached
            for stage_record in reached:
                stage = stage_record['stage']
                a = phase_ticks[1+2*stage]
                published = stage_record['published']
                end_index = 2+2*stage if published else mask.bit_count()-1
                if stage_record['entries']:
                    # Refused records must remain within the reached command phase prefix.
                    assert end_index >= 1+2*stage
                b = phase_ticks[end_index]
                band = (b-a)%MOD
                if published and stage_record['entries']:
                    phase_bands[','.join(sorted({m['family'] for e in stage_record['entries'] for m in e['actual_table_matches']}))].append(band)
                previous_end = 0
                for entry_record in stage_record['entries']:
                    relative = [(t-a)%MOD for t in entry_record['raw_ticks']]
                    assert all(x <= y for x,y in zip(relative,relative[1:])), 'Division clocks are not modularly ordered'
                    assert previous_end <= relative[0] and relative[-1] <= band
                    previous_end = relative[-1]
                    entry_record['relative_to_stage_input'] = relative
                    matches = entry_record['actual_table_matches']
                    if len(matches) != 1:
                        unclassified_records += 1
                        continue
                    match = matches[0]
                    # Published, valid full records may be summarized per actual table/node/stage.
                    # Unit/refused/partial records remain explicit, rather than receiving missing costs.
                    if not published or entry_record['terminal_reason'] != 1: continue
                    key = (match['table'],entry_record['node'],stage)
                    sample = division_samples[key]
                    sample['whole_observed_division'].append(relative[-1]-relative[0])
                    for cut,value in enumerate(entry_record['adjacent_modulo64_differences']):
                        sample[preparation['cut_layout'][cut]+'__to__'+preparation['cut_layout'][cut+1]].append(value)
        observations.append(observation)
    records.append({'capture': metadata, 'input_layout': input_record, 'mask_counts': dict(masks),
                    'completed_rows': completed, 'observed_rows': observations})

assert ordinals == {int(n.removeprefix('shader_clock_capture_')) for n in capture_properties}
assert input_files == {p.name for p in capture_directory.glob('*.input.bin')}
stages = {r['capture']['stage'] for r in records}
assert stages == {'device-timestamps': {'transport','endpoint'}, 'mixed-transport': {'transport'}, 'factored-endpoints': {'endpoint'}}[case]
assert any(r['completed_rows'] for r in records)
division_summary = []
for (key,node,stage),samples in sorted(division_samples.items()):
    division_summary.append({'table':key,'node':node,'stage':stage,
                             'inclusive_raw_clock_bands':{name:summary(values) for name,values in samples.items()}})
result = {'accepted':True,'source_revision':plan['source_revision'],'case':plan['cases'][case],
          'original_test_passed':True,'captures':len(records),'stages':sorted(stages),
          'records':sorted(records,key=lambda r:r['capture']['ordinal']),'actual_tables':table_records,
          'published_rhs_bands':{k:summary(v) for k,v in phase_bands.items()},
          'division_terminal_counts':[{'reason':k[0],'valid_word':k[1],'count':v} for k,v in sorted(terminal_counts.items())],
          'unclassified_division_records':unclassified_records,'published_valid_division_summary':division_summary,
          'scope':'Original full test assertions pass; complete sidecars/input snapshots hash/XML bound. All termination paths/full-vs-valid distinctions preserved. Literal count selects slots only; every classification joins captured full tables and actual opcode5 node. Modified/refused/unclassified tables remain raw. Modulo64 raw-clock bands only, inclusive of scheduling/control/observer traffic, perturbed; no calibrated time/exclusive primitive or total division cost/ordinary rendering benefit/full science or release verdict.'}
target = directory / 'observation-analysis.json'
assert not target.exists()
target.write_text(json.dumps(result, separators=(',',':')) + '\n')
print(json.dumps({k:result[k] for k in ('accepted','case','captures','stages')}))
