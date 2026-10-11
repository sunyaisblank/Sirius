"""Summarize original controls and inclusive layer bands without opcode cost allocation."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import collections, hashlib, json, statistics, xml.etree.ElementTree as ET
WORK = Path(__file__).resolve().parent
ROOT = Path.cwd()
MOD = 1 << 64
def read(p): return json.loads(p.read_text(encoding='utf-8-sig'))
def seal(p):
    data = p.read_bytes()
    return {'bytes':len(data), 'sha256':hashlib.sha256(data).hexdigest()}
def properties(p):
    xml = ET.parse(p).getroot()
    assert len(list(xml.iter('testcase'))) == 1 and not any(list(xml.iter(x)) for x in ('failure','error','skipped'))
    nodes = list(xml.iter('property'))
    assert len({n.get('name') for n in nodes}) == len(nodes)
    return {n.get('name'):n.get('value') for n in nodes}
def marker(value): return dict(p.split('=',1) for p in value.split(';'))
def median_record(values): return {'rows':len(values), 'median':statistics.median(values)}
plan = read(WORK / 'execution-plan.json')
result = {'accepted':True, 'source_revision':plan['source_revision'], 'cases':{},
          'scope':'All original test assertions and exact capture validators must pass. Raw inclusive diagnostic ticks only; body/join/other bands include observer and scheduling effects. Median row fractions cannot be added as an exact median partition. No exclusive opcode/barrier cost, calibrated time, ordinary rendering benefit or full scientific/release qualification.'}
for case in ('device-timestamps','mixed-transport','factored-endpoints'):
    w = WORK / case
    bridge = read(w / 'native-control-bridge.json')
    owner = read(w / 'run/owner.json')
    analysis = read(w / 'observation-analysis.json')
    assert bridge['accepted'] and analysis['accepted'] and owner['numerical_control_passed']
    assert read(w / 'retirement-run.json')['accepted'] and read(w / 'native-retirement.json')['pass']
    assert bridge['source_revision'] == analysis['source_revision'] == owner['source_revision'] == plan['source_revision']
    total_masks = collections.defaultdict(collections.Counter)
    totals = collections.Counter()
    distributions = collections.defaultdict(lambda: collections.defaultdict(list))
    for capture in analysis['records']:
        stage = capture['capture']['stage']
        totals[stage+'_captures'] += 1
        totals[stage+'_active_rows'] += capture['capture']['active_rows']
        totals[stage+'_completed_rows'] += capture['completed_rows']
        total_masks[stage].update({int(k):v for k,v in capture['mask_counts'].items()})
        for row in capture['observed_rows']:
            if stage != 'transport': continue
            totals['transport_reached_layers'] += sum(len(s['entries']) for s in row['layer_stages'])
            totals['transport_published_stages'] += sum(s['published'] for s in row['layer_stages'])
            if not row['completion']: continue
            assert len(row['family_candidates']) == 1
            family = row['family_candidates'][0]
            ticks = row['raw_ticks']
            rhs = sum((ticks[2+2*s]-ticks[1+2*s])%MOD for s in range(7))
            body = join = 0
            for s in row['layer_stages']:
                assert s['published']
                for start,end,a,b,c in s['entries']:
                    body += (b-a)%MOD
                    join += (c-b)%MOD
            assert rhs > 0 and body+join <= rhs
            values = {'rhs_ticks':rhs, 'layer_body_ticks':body, 'layer_join_ticks':join,
                      'other_inside_rhs_ticks':rhs-body-join, 'body_fraction_of_rhs':body/rhs,
                      'join_fraction_of_rhs':join/rhs, 'other_fraction_of_rhs':(rhs-body-join)/rhs}
            for name,value in values.items(): distributions[family][name].append(value)
    item = {'original_case':analysis['case'], 'native_elapsed_seconds':owner['elapsed_seconds'],
            'bridge_elapsed_seconds':bridge['elapsed_seconds'], 'peak_sampled_rss_bytes':owner['peak_sampled_rss_bytes'],
            'observations':dict(totals), 'phase_masks':{s:dict(c) for s,c in total_masks.items()},
            'complete_row_rhs_distribution':{f:{k:median_record(v) for k,v in d.items()} for f,d in distributions.items()},
            'analysis':seal(w / 'observation-analysis.json')}
    if case == 'device-timestamps':
        current = properties(w / 'run/gtest.xml')
        old_path = ROOT / 'attestations/native-vulkan/8f1118d/native-shader-clock-execution/files/device-timestamps/run/gtest.xml'
        prior = properties(old_path)
        readbacks = {k:v for k,v in current.items() if '_readback_' in k}
        assert len(readbacks) == 32 and readbacks == {k:v for k,v in prior.items() if '_readback_' in k}
        modules = {k:v for k,v in current.items() if k.startswith('selected_module_')}
        assert len(modules) == 7 and set(modules) == {k for k in prior if k.startswith('selected_module_')}
        expected = 'stage=transport;bytes=87928;sha256=e3a1589d852a4254312bb0e9e4b0343ab6aea32855c02c4f0f93b9fe1e18e5a3'
        assert modules['selected_module_1'] == expected
        assert all(v == prior[k] for k,v in modules.items() if k != 'selected_module_1')
        markers = {k:marker(v) for k,v in current.items() if k.startswith('device_observation_')}
        assert len(markers) == 51 and set(markers) == {k for k in prior if k.startswith('device_observation_')}
        for k,v in markers.items():
            previous = marker(prior[k])
            assert set(v) == set(previous)
            assert all(v[n] == previous[n] for n in ('phase','stage','available0','available1'))
            assert v['available0'] == v['available1'] == '1'
        item['original_timestamp_stream'] = {'unchanged_full_buffer_digests':32,
            'selected_modules':7, 'changed_module':'Exact new diagnostic Transport; six other modules unchanged',
            'ordered_available_markers':51, 'prior_xml':seal(old_path)}
    result['cases'][case] = item
target = WORK / 'results.json'
assert not target.exists()
target.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result))
