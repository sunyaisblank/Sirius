"""Read complete original execution receipts without launching tests."""
from pathlib import Path
import hashlib, importlib.util, json, xml.etree.ElementTree as ET

R = Path.cwd().resolve()
W = Path(__file__).resolve().parent
actions = json.loads((W / 'actions.json').read_text())
def identity(path):
    data = path.read_bytes()
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}
rows = []
for index in range(1, 9):
    label = 'control-' + str(index)
    folder = W / label
    owner = json.loads((folder / 'owner.json').read_text())
    assert owner['status'] == 'terminal' and owner['passed'] and owner['child_exit'] == 0
    assert not owner['stop_reason'] and not owner['cleanup_errors'] and owner['remaining_owned_processes'] == []
    assert owner['provider_mappings'] and owner['observed_births'] and all(b['absent'] for b in owner['observed_births'])
    assert owner['whole_inputs_unchanged'] and owner['source_unchanged']
    assert owner['argv'] == actions[label]['argv']
    before = json.loads((folder / 'inputs-before.json').read_text())
    assert before == json.loads((folder / 'inputs-after.json').read_text())
    for name, pin in before.items():
        assert identity(R / name) == {k: pin[k] for k in ['bytes','sha256']}, name
    xml = ET.parse(folder / 'gtest.xml').getroot()
    assert int(xml.get('tests', '-1')) == 1
    assert not any(xml.findall('.//' + tag) for tag in ['failure', 'error', 'skipped'])
    for node in [xml, *xml.findall('.//testsuite')]:
        for key in ['failures', 'errors', 'skipped', 'disabled']:
            assert int(node.get(key, '0')) == 0
    cases = xml.findall('.//testcase')
    assert len(cases) == 1 and cases[0].get('status') == 'run' and cases[0].get('result') == 'completed'
    full_name = cases[0].get('classname') + '.' + cases[0].get('name')
    assert actions[label]['argv'][1] == '--gtest_filter=' + full_name
    rows.append({'label': label, 'case': full_name, 'wall_seconds': owner['wall_seconds'],
                 'peak_sampled_rss_kib': owner['peak_sampled_rss_kib'], 'owner': identity(folder / 'owner.json'),
                 'xml': identity(folder / 'gtest.xml'), 'recorded_births': len(owner['observed_births'])})
spec = importlib.util.spec_from_file_location('immutable_original_reader', W / 'contract/software_contracts.py')
reader = importlib.util.module_from_spec(spec)
spec.loader.exec_module(reader)
props = reader.xml_properties(W / 'control-8/gtest.xml', 'timestamps')
timestamp = reader.timestamp_result(props, 'candidate')
baseline_props = reader.xml_properties(W / 'original-baseline-timestamps.xml', 'timestamps')
baseline_readbacks = {key: reader.fields(baseline_props[key]) for key in timestamp['readbacks']}
assert timestamp['readbacks'] == baseline_readbacks
assert len(timestamp['readbacks']) == 32 and len(timestamp['executed_modules']) == 54
protected = json.loads((W / 'protected-before.json').read_text())
for name, pin in protected.items():
    assert identity(R / name) == pin, name
result = {'passed': True, 'source_revision': actions['control-1']['revision'], 'original_controls': rows,
          'timestamp_schema': {'complete_readbacks': 32, 'ordered_markers': 51, 'actual_module_command_joins': 54,
                               'all_readbacks_exact_to_repaired_original_baseline': True,
                               'selected_modules': timestamp['selected_modules'],
                               'inclusive_API_schema_accepted': True},
          'protected_payloads_unchanged': len(protected), 'original_reader': identity(W / 'contract/software_contracts.py'),
          'scope': 'Eight original integrated software-provider controls, complete assertions and unchanged mathematical fixtures. No performance acceptance, Endpoint-core execution, production adoption or full scientific/rendering/release qualification verdict.'}
(W / 'control-acceptance.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps({'passed': True, 'original_controls': 8, 'complete_readbacks': 32, 'ordered_markers': 51,
                  'actual_module_command_joins': 54, 'wall_seconds': sum(row['wall_seconds'] for row in rows)}))
