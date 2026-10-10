from pathlib import Path
import hashlib, json, math
import xml.etree.ElementTree as ET

work = Path(__file__).resolve().parent
actions = json.loads((work / 'actions-final-tests.json').read_text())
results = {}
for key in ('final-regression', 'final-shared-tracer', 'final-cpu-science'):
    folder = work / key
    owner = json.loads((folder / 'owner.json').read_text())
    before = json.loads((folder / 'inputs-before.json').read_text())
    after = json.loads((folder / 'inputs-after.json').read_text())
    assert before == after
    assert owner['source'] == {'revision': actions[key]['revision'], 'status': ''}
    assert owner['passed'] and owner['child_exit'] == 0
    assert owner['stop_reason'] is None and not owner['cleanup_errors']
    assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
    assert owner['remaining_owned_processes'] == []
    assert owner['observed_births'] and all(x['absent'] for x in owner['observed_births'])
    expected = actions[key]['expected_test_name']
    documents = {}
    for filename in ('gtest.xml', 'ctest.xml'):
        path = folder / filename
        root = ET.parse(path).getroot()
        assert root.get('tests') == '1'
        for suite in [root, *root.findall('.//testsuite')]:
            assert all(int(suite.get(field, '0')) == 0
                       for field in ('failures', 'errors', 'skipped', 'disabled'))
        cases = root.findall('.//testcase')
        assert len(cases) == 1
        case = cases[0]
        name = case.get('classname') + '.' + case.get('name') if filename == 'gtest.xml' else case.get('name')
        assert name == expected and case.get('status') == 'run'
        assert not any(root.findall('.//' + tag) for tag in ('failure', 'error', 'skipped'))
        if filename == 'gtest.xml':
            assert case.get('result') == 'completed'
        documents[filename] = {'bytes': path.stat().st_size,
                               'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}
    results[key] = {'case': expected, 'passed': True, 'sealed_inputs': len(before),
                    'xml': documents, 'owner_wall_seconds': owner['wall_seconds'],
                    'peak_sampled_rss_kib': owner['peak_sampled_rss_kib']}

cpu = ET.parse(work / 'final-cpu-science/gtest.xml').getroot()
names = ('flat_moving_pupil', 'schwarzschild_capture', 'schwarzschild_turn_mass_01',
         'schwarzschild_critical_plus_005', 'kerr_moving_off_plane', 'kerr_0998_inner_turn',
         'kerr_minus_0998_mass_100', 'kerr_extremal_capture', 'kerr_minus_extremal_outward',
         'schwarzschild_first_disk', 'kerr_0998_moving_first_disk')
properties = cpu.findall('.//property')
expected = {name + '_r' + str(level) for name in names for level in range(3)}
assert len(properties) == 33 and {p.get('name') for p in properties} == expected
for prop in properties:
    values = dict(field.split('=', 1) for field in prop.get('value').split(';'))
    assert all(math.isfinite(float(value)) for value in values.values())
    assert float(values['attempts']) > 0
results['final-cpu-science']['independent_witness_refinement_records'] = 33
bindings = json.loads((work / 'linux-readonly-payloads.json').read_text())
assert bindings['source_revision'] == actions['final-regression']['revision'] and bindings['pass_']
assert len(bindings['actual_consumers']) == 4
assert all(len(consumer['arrays']) == 40 for consumer in bindings['actual_consumers'].values())
for relative, consumer in bindings['actual_consumers'].items():
    assert hashlib.sha256(Path(relative).read_bytes()).hexdigest() == consumer['sha256']
record = {'source_revision': actions['final-regression']['revision'], 'passed': True,
          'results': results, 'unchanged_complete_kernel_payload_joins': 160,
          'scope': 'Focused integrated issue-91 acceptance only; no full frame, scientific estate, platform/runtime or release admission.'}
(work / 'payload-validation.json').write_text(json.dumps(record, indent=2) + '\n')
print(json.dumps(record, indent=2))
