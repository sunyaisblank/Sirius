"""Read completed current AP controls before preparing or scheduling any ABBA slot."""
from pathlib import Path
import hashlib
import importlib.util
import json


def check():
    root = Path.cwd()
    current = root / 'out/native-cpu-endpoint-header-cohort'
    controls = current / 'controls'
    runtime = controls / 'runtime'
    def doc(path):
        return json.loads(path.read_text())
    def pin(path):
        with path.open('rb') as stream:
            return dict(bytes=path.stat().st_size, sha256=hashlib.file_digest(stream, 'sha256').hexdigest())
    c8 = doc(current / 'control-8-acceptance.json')
    assert c8['passed'] and c8['optional_math']['all_selected_optional_math_observed']
    candidate = doc(current / 'actions.json')['control-8']
    actions = doc(runtime / 'actions.json')
    assert set(actions) == {'control-' + str(i) for i in range(1, 8)}
    expected_identities = actions['control-1']['expected_identities']
    spec = importlib.util.spec_from_file_location('AP_original_control_XML', controls / 'readers/native_evidence.py')
    original = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(original)
    inputs = [runtime / 'actions.json', current / 'actions.json', current / 'control-8-acceptance.json']
    for i in range(1, 8):
        label = 'control-' + str(i)
        folder, action = runtime / label, actions[label]
        owner = doc(folder / 'owner.json')
        assert owner['passed'] and owner['status'] == 'terminal' and owner['child_exit'] == 0
        assert owner['stop_reason'] is None and owner['source_unchanged'] and owner['whole_inputs_unchanged']
        assert owner['source'] == dict(revision=action['revision'], status='')
        assert not owner['remaining_owned_processes'] and not owner['cleanup_errors']
        assert owner['observed_births'] and all(row['absent'] for row in owner['observed_births'])
        assert action['limits'] == owner['limits'] == dict(seconds=300, rss_kib=4194304, reserve_kib=2097152)
        assert owner['wall_seconds'] <= 300 and owner['peak_sampled_rss_kib'] <= 4194304
        assert owner['minimum_sampled_available_kib'] >= 2097152
        assert action['prerequisites'] == ['control-8'] + ['control-' + str(n) for n in range(1, i)]
        assert action['case'] == original.CASES[label]
        for key in ('executable', 'runtime_provider', 'runtime_icd'):
            assert action[key] == candidate[key]
        assert owner['declared_output_environment'] == {
            'GTEST_OUTPUT': 'xml:' + str(folder / 'gtest.xml'),
            'VK_DRIVER_FILES': candidate['runtime_icd'], 'SIRIUS_VULKAN_DEVICE': '0',
            'MESA_SHADER_CACHE_DISABLE': 'true'}
        assert owner['provider_mappings'] and all(row['owned_session'] == owner['child_pid']
            and row['maps'] and all(candidate['runtime_provider'] in line for line in row['maps'])
            for row in owner['provider_mappings'])
        assert doc(folder / 'inputs-before.json') == doc(folder / 'inputs-after.json')
        original.xml_properties(folder / 'gtest.xml', label)
        assert action['expected_identities'] == expected_identities
        inputs += [folder / name for name in ('owner.json', 'gtest.xml', 'inputs-before.json', 'inputs-after.json')]
    for name, expected in expected_identities.items():
        assert pin(root / name) == expected, name
    inputs += [path for path in controls.rglob('*.py')]
    return inputs
