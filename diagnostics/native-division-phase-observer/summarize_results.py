"""Join original numerical/stream controls and inclusive division observations."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import collections, hashlib, json, xml.etree.ElementTree as ET
ROOT = Path.cwd(); WORK = Path(__file__).resolve().parent

def read(path): return json.loads(path.read_text(encoding='utf-8-sig'))
def seal(path):
    data = path.read_bytes()
    return {'bytes':len(data), 'sha256':hashlib.sha256(data).hexdigest()}
def properties(path):
    xml = ET.parse(path).getroot()
    cases = list(xml.iter('testcase'))
    assert len(cases) == 1 and not any(list(xml.iter(x)) for x in ('failure','error','skipped'))
    nodes = list(cases[0].iter('property'))
    assert len({n.get('name') for n in nodes}) == len(nodes)
    return {n.get('name'):n.get('value') for n in nodes}
def marker(value): return dict(p.split('=',1) for p in value.split(';'))

prior_xml = ROOT/'attestations/native-vulkan/24a8e14/native-rhs-layer-clock-execution/files/device-timestamps/run/gtest.xml'
plan = read(WORK/'execution-plan.json')
review = read(WORK/'summary-source-review.json')
assert review['accepted'] and review['source_revision'] == plan['source_revision']
required = (Path(__file__),WORK/'execution-plan.json',WORK/'reader-source-review.json',
            WORK/'analyze_division_observations.py',WORK/'candidate-source.json',
            WORK/'host-source-review.json',prior_xml)
assert all(str(path) in review['reviewed_seals'] for path in required)
assert review['reviewed_seals'][str(Path(__file__))] == seal(Path(__file__))
for path, expected in review['reviewed_seals'].items(): assert seal(Path(path)) == expected
inputs = {path:seal(Path(path)) for path in review['reviewed_seals']}
inputs[str(WORK/'summary-source-review.json')] = seal(WORK/'summary-source-review.json')
inputs[str(prior_xml)] = seal(prior_xml)
result = {'accepted':True,'source_revision':plan['source_revision'],'cases':{},
    'scope':'Original assertions, complete output/stream joins and capture validators. Raw inclusive diagnostic ticks are perturbed by observation and SIMT scheduling; no calibrated time, exclusive primitive/total division cost, ordinary rendering benefit or full scientific/release verdict.'}
for case in ('device-timestamps','mixed-transport','factored-endpoints'):
    directory = WORK/case
    bridge = read(directory/'native-control-bridge.json'); owner = read(directory/'run/owner.json')
    analysis = read(directory/'observation-analysis.json')
    retirement = read(directory/'retirement-run.json'); native = read(directory/'native-retirement.json')
    assert bridge['accepted'] and analysis['accepted'] and owner['numerical_control_passed']
    assert retirement['accepted'] and native['pass']
    assert bridge['source_revision'] == analysis['source_revision'] == owner['source_revision'] == plan['source_revision']
    assert analysis['case'] == plan['cases'][case]
    for path in directory.rglob('*'):
        if path.is_file(): inputs[str(path)] = seal(path)
    masks = collections.defaultdict(collections.Counter); totals = collections.Counter()
    for capture in analysis['records']:
        stage = capture['capture']['stage']
        totals[stage+'_captures'] += 1
        totals[stage+'_active_rows'] += capture['capture']['active_rows']
        totals[stage+'_completed_rows'] += capture['completed_rows']
        masks[stage].update({int(k):v for k,v in capture['mask_counts'].items()})
        for row in capture['observed_rows']:
            if stage != 'transport': continue
            totals['transport_reached_divisions'] += sum(len(s['entries']) for s in row['division_stages'])
            totals['transport_published_rhs_stages'] += sum(s['published'] for s in row['division_stages'])
    summaries = analysis['published_valid_division_summary']
    measurement_available = bool(summaries) and all(
        sample['inclusive_raw_clock_bands']['whole_observed_division']['count'] > 0 for sample in summaries)
    # Endpoints have no selected division observer; availability belongs to Transport.
    if case != 'factored-endpoints': assert measurement_available
    item = {'original_case':analysis['case'],'native_elapsed_seconds':owner['elapsed_seconds'],
        'bridge_elapsed_seconds':bridge['elapsed_seconds'],'peak_sampled_rss_bytes':owner['peak_sampled_rss_bytes'],
        'observations':dict(totals),'phase_masks':{s:dict(c) for s,c in masks.items()},
        'actual_full_tables':len(analysis['actual_tables']),
        'division_terminal_counts':analysis['division_terminal_counts'],
        'unclassified_division_records':analysis['unclassified_division_records'],
        'measurement_available':measurement_available,'published_valid_table_node_stage_groups':len(summaries),
        'published_valid_division_summary':summaries,'analysis':seal(directory/'observation-analysis.json')}
    if case == 'device-timestamps':
        current = properties(directory/'run/gtest.xml'); prior = properties(prior_xml)
        readbacks = {k:v for k,v in current.items() if '_readback_' in k}
        assert len(readbacks) == 32 and readbacks == {k:v for k,v in prior.items() if '_readback_' in k}
        modules = {k:v for k,v in current.items() if k.startswith('selected_module_')}
        assert len(modules) == 7 and set(modules) == {k for k in prior if k.startswith('selected_module_')}
        observer = plan['exact_observer_modules']['retained_transport.spv']
        expected = 'stage=transport;bytes='+str(observer['bytes'])+';sha256='+observer['sha256']
        assert modules['selected_module_1'] == expected
        assert all(v == prior[k] for k,v in modules.items() if k != 'selected_module_1')
        markers = {k:marker(v) for k,v in current.items() if k.startswith('device_observation_')}
        assert len(markers) == 51 and set(markers) == {k for k in prior if k.startswith('device_observation_')}
        for k,v in markers.items():
            previous = marker(prior[k]); assert set(v) == set(previous)
            assert all(v[n] == previous[n] for n in ('phase','stage','available0','available1'))
            assert v['available0'] == v['available1'] == '1'
        item['original_timestamp_stream'] = {'unchanged_full_buffer_digests':32,'selected_modules':7,
            'changed_module':'Exact division diagnostic Transport; six other modules unchanged',
            'ordered_available_markers':51,'prior_xml':seal(prior_xml)}
    result['cases'][case] = item
assert {path:seal(Path(path)) for path in inputs} == inputs
result['input_seals'] = inputs
path = WORK/'results.json'; assert not path.exists()
path.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'accepted':result['accepted'],'source_revision':result['source_revision'],
    'cases':{name:{k:v for k,v in data.items() if k not in ('published_valid_division_summary',)} for name,data in result['cases'].items()}}))
