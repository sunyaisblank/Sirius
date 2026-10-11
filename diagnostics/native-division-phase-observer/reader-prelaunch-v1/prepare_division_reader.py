"""Prepare a strict captured-record reader; no shader or test execution."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import hashlib,json
WORK=Path(__file__).resolve().parent
OLD=Path.cwd()/'attestations/native-vulkan/24a8e14/native-rhs-layer-clock-execution/files/analyze_observations.py'
s=OLD.read_text()
def one(a,b):
    global s
    assert s.count(a)==1,(a,s.count(a))
    s=s.replace(a,b)
one('Validate exact captures and bind reached layer ranges to actual input tables.', 'Validate division captures and bind selected records to actual complete input tables.')
one("assert placement['module'] == read(WORK / 'compiled-review.json')['module']", "assert placement['module'] == read(WORK / 'compiler.json')['product']['spv']\nassert read(WORK / 'actual-placement-review.json')['accepted']\npreparation = read(WORK / 'preparation.json')\nselections = preparation['selections']")
one("layer_samples = collections.defaultdict(lambda: {'body': [], 'join': [], 'span': []})", "division_samples = collections.defaultdict(lambda: collections.defaultdict(list))\nterminal_counts = collections.Counter()\nunclassified_records = 0")
a=s.index('            candidates = set(tables)\n')
b=s.index('        observations.append(observation)\n',a)
s=s[:a]+'''            candidates = set(tables)
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
                # A refused RHS has no after-publication cut. Preserve that absence;
                # its input cut orders records but supplies no enclosing upper bound.
                band = (phase_ticks[2+2*stage]-a)%MOD if published else None
                if published and stage_record['entries']:
                    phase_bands[','.join(sorted({m['family'] for e in stage_record['entries'] for m in e['actual_table_matches']}))].append(band)
                previous_end = 0
                for entry_record in stage_record['entries']:
                    relative = [(t-a)%MOD for t in entry_record['raw_ticks']]
                    assert all(x <= y for x,y in zip(relative,relative[1:])), 'Division clocks are not modularly ordered'
                    assert previous_end <= relative[0]
                    if band is not None: assert relative[-1] <= band
                    previous_end = relative[-1]
                    entry_record['enclosing_after_publication_cut_available'] = published
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
''' +s[b:]
a=s.index('layer_summary = []\n')
b=s.index('target = directory / ',a)
s=s[:a]+'''division_summary = []
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
''' +s[b:]
p=WORK/'analyze_division_observations.py';assert not p.exists();p.write_text(s)
print(json.dumps({'output':str(p),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'scope':'Reader source preparation only; not run.'}))
