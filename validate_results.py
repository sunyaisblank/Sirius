"""Check the original finite restoration controls and whole input seals."""
from pathlib import Path
import hashlib, json, xml.etree.ElementTree as ET

WORK = Path(__file__).resolve().parent
actions = json.loads((WORK / 'actions.json').read_text())
results = []
def readbacks(path):
    props = {p.get('name'): p.get('value') for p in ET.parse(path).findall('.//property')}
    return {k: v for k, v in props.items() if '_readback_' in k}
reference = readbacks(WORK.parent / 'submission-fence-review/candidate-timestamps/gtest.xml')
for label, action in actions.items():
    folder = WORK / label
    owner = json.loads((folder / 'owner.json').read_text())
    assert owner['passed'] and owner['child_exit'] == 0 and not owner['stop_reason']
    assert not owner['cleanup_errors'] and owner['remaining_owned_processes'] == []
    assert all(v['absent'] for v in owner['observed_births'])
    assert (folder / 'inputs-before.json').read_bytes() == (folder / 'inputs-after.json').read_bytes()
    result = {'case': label, 'seconds': owner['wall_seconds'], 'input_count': len(json.loads((folder / 'inputs-before.json').read_text())),
              'owner_sha256': hashlib.sha256((folder / 'owner.json').read_bytes()).hexdigest()}
    if action['kind'] == 'test':
        root = ET.parse(folder / 'gtest.xml').getroot()
        assert root.get('tests') == '1' and all(root.get(k, '0') == '0' for k in ('failures', 'errors', 'disabled'))
        cases = root.findall('.//testcase')
        assert len(cases) == 1
        actual = cases[0].get('classname') + '.' + cases[0].get('name')
        assert actual == action['expected_test_name']
        assert not root.findall('.//failure') and not root.findall('.//error') and not root.findall('.//skipped')
        assert cases[0].get('status') == 'run' and cases[0].get('result') == 'completed'
        result.update(test=actual, xml_sha256=hashlib.sha256((folder / 'gtest.xml').read_bytes()).hexdigest())
        if label == 'candidate-timestamps':
            assert readbacks(folder / 'gtest.xml') == reference
            markers = [p.get('value') for p in root.findall('.//property') if p.get('name').startswith('device_observation_')]
            assert len(markers) == 51 and all('submission_fence=' not in v and 'available0=1;available1=1' in v for v in markers)
            result.update(ordered_markers=51, same_provider_complete_readback_records=len(reference), readbacks_equal_trial=True)
    results.append(result)
(WORK / 'software-controls.json').write_text(json.dumps({'pass': True, 'results': results,
    'scope': 'Exact accepted source tree restoration. Strict configure/build and four original finite controls; historical CPU33 evidence remains at312bda3. No full-frame/science/runtime/release admission.'}, indent=2) + '\n')
print(json.dumps({'pass': True, 'executions': len(results), 'original_test_cases': 4, 'readback_records': len(reference)}))
