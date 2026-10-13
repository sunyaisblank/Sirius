from pathlib import Path
import hashlib, json, xml.etree.ElementTree as ET

R = Path.cwd()
W = Path(__file__).resolve().parent
actions = json.loads((W / 'actions.json').read_text())
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
rows = []
for i in range(1, 12):
    label = 'control-' + str(i)
    folder = W / label
    owner = json.loads((folder / 'owner.json').read_text())
    assert owner['status'] == 'terminal' and owner['passed'] and owner['child_exit'] == 0
    assert not owner['stop_reason'] and not owner['cleanup_errors']
    assert owner['remaining_owned_processes'] == []
    assert all(birth['absent'] for birth in owner['observed_births'])
    assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
    assert owner['source'] == {'revision': actions[label]['revision'], 'status': ''}
    assert owner['argv'] == actions[label]['argv']
    assert owner['provider_mappings']
    before = json.loads((folder / 'inputs-before.json').read_text())
    assert before == json.loads((folder / 'inputs-after.json').read_text())
    xml = ET.parse(folder / 'gtest.xml').getroot()
    assert int(xml.get('tests', '-1')) == 1
    assert not any(xml.findall('.//' + name) for name in ['failure', 'error', 'skipped'])
    for element in [xml, *xml.findall('.//testsuite')]:
        for name in ['failures', 'errors', 'skipped', 'disabled']:
            assert int(element.get(name, '0')) == 0
    cases = xml.findall('.//testcase')
    assert len(cases) == 1
    case = cases[0]
    assert case.get('status') == 'run' and case.get('result') == 'completed'
    name = case.get('classname') + '.' + case.get('name')
    assert actions[label]['argv'][1] == '--gtest_filter=' + name
    rows.append({'label': label, 'case': name, 'wall_seconds': owner['wall_seconds'],
                 'owner_sha256': sha(folder / 'owner.json'), 'xml_sha256': sha(folder / 'gtest.xml'),
                 'recorded_births': len(owner['observed_births'])})
result = {'passed': True, 'candidate_revision': actions['control-1']['revision'],
          'binding_controls': rows[:3], 'original_retained_controls': rows[3:],
          'scope': 'Eleven finite executed cases with zero failures/errors/skips/disabled. Original eight correctness prerequisites pass; no performance, complete scientific or release qualification verdict.'}
destination = W / 'control-acceptance.json'
assert not destination.exists()
destination.write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps({'passed': True, 'binding_controls': 3, 'original_controls': 8,
                  'wall_seconds': sum(x['wall_seconds'] for x in rows)}))
