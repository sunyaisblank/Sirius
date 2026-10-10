"""Preserved original native-owner validation, scoped to FP64 science."""
from pathlib import Path
import hashlib, json, subprocess, xml.etree.ElementTree as ET
ROOT = Path(__file__).resolve().parent.parent.parent
WORK = Path(__file__).resolve().parent
CASES = {'science': 'RetainedComputeTest.Fp64ProductsPreserveIndependentScienceOrDeclineUnsupportedDevices'}
def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def identity(path):
    return {'path': str(path.relative_to(ROOT)), 'bytes': path.stat().st_size, 'sha256': sha(path)}

def read(path):
    assert path.is_file() and path.stat().st_size <= 16 * 1024 * 1024, str(path)
    return path.read_text(encoding='utf-8-sig')

def document(path):
    def invalid(value):
        raise ValueError('Nonfinite JSON constant: ' + value)
    return json.loads(read(path), parse_constant=invalid)

def xml_properties(path, case):
    xml = ET.fromstring(read(path))
    assert int(xml.get('tests', '-1')) == 1
    for node in [xml, *xml.findall('.//testsuite')]:
        for key in ('failures', 'errors', 'skipped', 'disabled'):
            assert int(node.get(key, '0')) == 0
    assert not any(xml.findall('.//' + tag) for tag in ('failure', 'error', 'skipped'))
    cases = xml.findall('.//testcase')
    assert len(cases) == 1 and cases[0].get('classname') + '.' + cases[0].get('name') == CASES[case]
    assert cases[0].get('status') == 'run' and cases[0].get('result') == 'completed'
    props = {}
    for item in xml.findall('.//property'):
        key = item.get('name')
        assert key and key not in props
        props[key] = item.get('value')
    return props

def read_case(mode, label, case, bindings):
    assert (mode, label, case) == ('controls', 'original', 'science')
    assert subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip() == bindings['candidate_revision']
    assert not subprocess.check_output(['git', 'status', '--porcelain'], cwd=ROOT)
    for relative, seal in bindings['files'].items():
        path = Path(relative)
        if not path.is_absolute():
            path = ROOT / path
        assert path.stat().st_size == seal['bytes'] and sha(path) == seal['sha256'], relative
    producer = 'baseline' if mode == 'baseline' or (mode == 'matched' and label.startswith('a')) else 'candidate'
    folder = WORK / 'native' / mode / label / case
    owner = document(folder / 'owner.json')
    bridge_path = WORK / ('native-' + mode + '-' + label + '-' + case + '-bridge.json')
    bridge = document(bridge_path)
    assert bridge['completed'] and bridge['accepted_owner_execution'] and bridge['returncode'] == 0
    binding_path = WORK / ('baseline-native-bindings.json' if mode == 'baseline' else 'execution-bindings.json')
    assert bridge['execution_bindings_sha256'] == sha(binding_path)
    assert bridge['input_bindings'] == bindings['files']
    assert bridge['producer'] == producer and bridge['case'] == case and bridge['mode'] == mode and bridge['label'] == label
    assert bridge['input_seals_current_and_exact_at_end'] and bridge['execution_bindings_unchanged_at_end']
    assert bridge['source_tree_clean_at_start'] and bridge['source_tree_clean_at_end']
    assert bridge['live_source_revision'] == bridge['live_source_revision_at_end'] == bindings['candidate_revision']
    assert owner['status'] == 'completed' and owner['returncode'] == 0 and not owner['outer_stop']
    assert not owner['cleanup_errors']
    for key in ('owned_process_absent', 'numerical_control_passed', 'test_executable_unchanged', 'loaded_driver_modules_recorded', 'source_build_gate_unchanged'):
        assert owner[key] is True
    assert owner['expected_cases'] == [CASES[case]]
    assert owner['outer_guard_seconds'] == 300 and owner['rss_guard_bytes'] == 4294967296
    assert owner['source_revision'] == bridge['source_revision'] == bindings['producers'][producer]['source_revision']
    assert owner['live_source_revision'] == bindings['candidate_revision']
    assert owner['source_build_gate_sha256'] == bindings['producers'][producer]['gate_sha256']
    exe = bindings['producers'][producer]['stage'] + '/tests/backend/Release/sirius_backend_tests.exe'
    assert owner['test_executable_sha256'] == bindings['files'][exe]['sha256']
    for key, word in owner['environment'].items():
        if key not in ('TEMP', 'TMP'):
            assert word in (None, '')
    for key in ('owner_pid', 'test_pid'):
        assert type(owner[key]) is int and owner[key] > 0
    for key in ('owner_start_time_utc', 'test_start_time_utc'):
        assert owner[key]
    modules = document(folder / 'driver-modules.json')
    assert isinstance(modules, list) and len(modules) == 2
    module_identities = []
    bound_paths = {path.casefold(): path for path in bindings['files'] if path.startswith('/mnt/c/')}
    assert {module['name'].lower() for module in modules} == {'vulkan-1.dll', 'amdvlk64.dll'}
    for module in modules:
        path = '/mnt/' + module['path'][0].lower() + module['path'][2:].replace('\\', '/')
        assert {'bytes': module['bytes'], 'sha256': module['sha256']} == bindings['files'][bound_paths[path.casefold()]]
        module_identities.append({k: module[k].casefold() if k in ('name', 'path') else module[k] for k in ('name', 'path', 'bytes', 'sha256', 'version')})
    props = xml_properties(folder / 'gtest.xml', case)
    plan = document(WORK / 'native-plan.json')
    assert plan['source_revision'] == bindings['candidate_revision']
    for key in ('selected_device', 'selected_driver'):
        assert props[key] == plan[key], (key, props[key])
    for key, value in plan['required_capabilities'].items():
        assert props[key] in ('0', '1', 'false', 'true')
        assert (props[key] in ('1', 'true')) is value, (key, props[key])
    record = {'case': case, 'label': label, 'source_revision': owner['source_revision'],
              'owner': owner, 'modules': sorted(module_identities, key=lambda x: x['name'].lower()),
              'artifacts': [identity(folder / name) for name in ('owner.json', 'gtest.xml', 'stdout.log', 'stderr.log', 'driver-modules.json', 'memory.json')] + [identity(bridge_path)]}
    return record
