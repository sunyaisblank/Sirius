import json
from pathlib import Path
import re
import sys
import xml.etree.ElementTree as ET

baseline, candidate, output, route_count = sys.argv[1:]
baseline, candidate, output = Path(baseline), Path(candidate), Path(output)
route_count = int(route_count)

def properties(folder, case):
    root = ET.parse(folder/case/'gtest.xml').getroot()
    assert int(root.get('tests')) == 1
    assert all(int(root.get(k, 0)) == 0 for k in ['failures', 'errors', 'disabled'])
    assert all(tc.find('skipped') is None for tc in root.iter('testcase'))
    return {p.get('name'): p.get('value') for p in root.iter('property')}

def endpoint_timing(folder):
    logs = ''.join(p.read_text() for p in (folder/'factored').glob('*.log'))
    records = re.findall(r'returned: call=(\d+);success=(\d+);pipeline_ms=([\d.]+);submit_wait_ms=([\d.]+);total_ms=([\d.]+)', logs)
    assert len(records) == 14 * route_count
    rows = [{'call': int(a), 'success': int(b), 'pipeline_ms': float(c), 'wait_ms': float(d), 'total_ms': float(e)} for a,b,c,d,e in records]
    for i in range(route_count):
        assert [r['call'] for r in rows[14*i:14*i+14]] == list(range(14))
    assert all(r['success'] == 1 for r in rows)
    return [{'route': i, 'initial_pipeline_ms': rows[14*i]['pipeline_ms'],
             'warm13_wait_ms': sum(r['wait_ms'] for r in rows[14*i+1:14*i+14]),
             'all14_calls': rows[14*i:14*i+14]} for i in range(route_count)]

def markers(properties):
    result = []
    for i in range(51):
        row = dict(field.split('=', 1) for field in properties['device_observation_'+str(i)].split(';'))
        assert row['available0'] == row['available1'] == '1'
        result.append(row)
    return result

def phase_totals(rows):
    totals = {}
    for row in rows:
        phase = totals.setdefault(row['phase'], {'device_ms': 0, 'host_wait_ms': 0, 'endpoint_device_ms': 0, 'endpoint_wait_ms': 0})
        for field in ['device_ms', 'host_wait_ms']:
            phase[field] += float(row[field])
        if row['stage'].startswith('endpoint'):
            phase['endpoint_device_ms'] += float(row['device_ms'])
            phase['endpoint_wait_ms'] += float(row['host_wait_ms'])
    return totals

a,b = (properties(folder,'factored') for folder in [baseline,candidate])
ka = {k for k in a if k.startswith('factored_endpoint_readback_')}
assert ka == {k for k in b if k.startswith('factored_endpoint_readback_')} and len(ka) == 14*route_count
c,d = (properties(folder,'timestamps') for folder in [baseline,candidate])
kc = {k for k in c if k.startswith('coupled_readback_')}
assert kc == {k for k in d if k.startswith('coupled_readback_')} and len(kc) == 11
ma,mb = markers(c),markers(d)
assert [(r['phase'],r['stage']) for r in ma] == [(r['phase'],r['stage']) for r in mb]
report = {'baseline_controls': str(baseline), 'candidate_controls': str(candidate),
          'original_endpoint_and_complete_coupled_controls_pass': True,
          'factored_complete_buffer_hashes_equal': sum(a[k] == b[k] for k in ka),
          'factored_complete_buffer_hashes_changed': sum(a[k] != b[k] for k in ka),
          'coupled_complete_buffer_hashes_equal': sum(c[k] == d[k] for k in kc),
          'coupled_complete_buffer_hashes_changed': sum(c[k] != d[k] for k in kc),
          'all51_markers_available_in_original_phase_stage_order': True,
          'baseline_endpoint_routes': endpoint_timing(baseline),
          'candidate_endpoint_routes': endpoint_timing(candidate),
          'baseline_coupled_phases': phase_totals(ma),
          'candidate_coupled_phases': phase_totals(mb),
          'cold_pipeline_cache_state_unverified': True,
          'qualification_claimed': False}
output.write_text(json.dumps(report, indent=2)+'\n')
print(json.dumps({k:v for k,v in report.items() if k not in ['baseline_endpoint_routes','candidate_endpoint_routes']}))
print(json.dumps({'baseline_warm13_ms': [r['warm13_wait_ms'] for r in report['baseline_endpoint_routes']],
                  'candidate_warm13_ms': [r['warm13_wait_ms'] for r in report['candidate_endpoint_routes']]}))
