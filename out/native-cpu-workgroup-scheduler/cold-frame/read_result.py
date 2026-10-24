"""Read the sole completed original cold-frame owner and original test XML."""
from pathlib import Path
import hashlib
import importlib.util
import json

root = Path.cwd()
work = Path(__file__).resolve().parent
current = root / 'out/native-cpu-workgroup-scheduler'


def doc(path):
    return json.loads(path.read_text())


def pin(path):
    return {'bytes': path.stat().st_size,
            'sha256': hashlib.file_digest(path.open('rb'), 'sha256').hexdigest()}


controls = doc(current / 'controls/readers/controls-1-8-acceptance.json')
assert controls['all8_original_controls_pass'] and controls['all22_original_benefit_checks_pass']
benefit = doc(current / 'benefit/readers/benefit.json')
decimal = doc(current / 'benefit/readers/decimal-agreement.json')
assert benefit['audit_passed'] and benefit['all22_numeric_pass'] and benefit['checks_passed'] == 22
assert decimal['agrees_on_all22'] and decimal['all22_pass']
assert doc(current / 'control-8-acceptance.json')['passed']
folder = work / 'frame'
owner = doc(folder / 'owner.json')
action = doc(work / 'actions.json')['frame']
assert owner['status'] == 'terminal' and owner['source'] == {'revision': action['revision'], 'status': ''}
assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
assert not owner['cleanup_errors'] and owner['remaining_owned_processes'] == []
assert owner['observed_births'] and all(x['absent'] for x in owner['observed_births'])
before = doc(folder / 'inputs-before.json')
assert before == doc(folder / 'inputs-after.json')
for name, value in before.items():
    assert pin(root / name) == {k: value[k] for k in ('bytes', 'sha256')}, name
limits = {'seconds': 300, 'rss_kib': 4194304, 'reserve_kib': 2097152}
assert owner['limits'] == action['limits'] == limits
provider = doc(current / 'benefit/providers.json')['providers']['candidate']
for key in ('library', 'ICD'):
    assert pin(Path(provider[key])) == provider[key + '_identity']
assert owner['declared_output_environment'] == {
    'GTEST_OUTPUT': 'xml:' + str(folder / 'gtest.xml'), 'VK_DRIVER_FILES': provider['ICD'],
    'SIRIUS_VULKAN_DEVICE': '0', 'MESA_SHADER_CACHE_DISABLE': 'true'}
assert owner['provider_mappings'] and all(x['owned_session'] == owner['child_pid']
    and x['maps'] and all(provider['library'] in line for line in x['maps'])
    for x in owner['provider_mappings'])
build = doc(work / 'build-result.json')
assert build['passed'] and build['new_CPP_compilations'] == 0
executable = root / build['executable']
assert pin(executable) == build['executable_identity']
assert action['argv'] == ['/usr/bin/cmake', '-DSIRIUS_RENDER_TESTS=' + str(executable),
                          '-P', str(root / 'tests/operational/vulkan_runtime_profile.cmake')]
assert action['runtime_provider'] == provider['library'] and action['runtime_icd'] == provider['ICD']
assert set(build['module_joins']) == {'Transport', 'Endpoint', 'DopriPhase'}
modules = doc(root / 'out/native-cpu-row-ownership/compiler-result.json')['modules']
data = executable.read_bytes()
for kind, row in build['module_joins'].items():
    original = modules[kind]
    assert row == {'module': original['module'], 'bytes': original['bytes'],
                   'sha256': original['sha256'], 'complete_occurrences': 1}
    module = root / row['module']
    assert pin(module) == {k: row[k] for k in ('bytes', 'sha256')}
    assert data.count(module.read_bytes()) == 1
accepted = owner['passed']
if accepted:
    assert owner['child_exit'] == 0 and owner['stop_reason'] is None
    assert owner['wall_seconds'] <= 300 and owner['peak_sampled_rss_kib'] <= 4194304
    assert owner['minimum_sampled_available_kib'] >= 2097152
    basis = root / 'out/native-cpu-arithmetic-sharing/readers/native_evidence.py'
    assert pin(basis)['sha256'] == '1b44f0678dbfc7ee455b3d49e170be884b294882896fd1ac3706ff5ebded20ab'
    spec = importlib.util.spec_from_file_location('original_frame_xml', basis)
    original = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(original)
    original.CASES['frame'] = 'RenderCommandParse.ExplicitGpuRequestRunsVulkanWhenDevicePresent'
    original.xml_properties(folder / 'gtest.xml', 'frame')
else:
    assert owner['stop_reason'] is not None or owner['child_exit'] != 0
result = {'attempted': True, 'accepted': accepted, 'owner': pin(folder / 'owner.json'),
    'XML': pin(folder / 'gtest.xml') if (folder / 'gtest.xml').exists() else None,
    'stop_reason': owner['stop_reason'], 'child_exit': owner['child_exit'],
    'wall_seconds': owner['wall_seconds'], 'peak_sampled_rss_kib': owner['peak_sampled_rss_kib'],
    'build': pin(work / 'build-result.json'), 'module_joins': build['module_joins'],
    'controls': pin(current / 'controls/readers/controls-1-8-acceptance.json'),
    'benefit': pin(current / 'benefit/readers/benefit.json'),
    'decimal': pin(current / 'benefit/readers/decimal-agreement.json'),
    'scope': 'One original cold128x128/1SPP Kerr runtime profile, unchanged defaults/governors, cache disabled. Original test validates and removes its PPM. No full scientific, Mandatory or release qualification claim.'}
destination = work / 'result.json'
assert not destination.exists()
destination.write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps({'accepted': accepted, 'owner_terminal': True, 'complete_module_joins': 3}))
