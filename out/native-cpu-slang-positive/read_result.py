"""Original complete C8 parser plus current native-entry/owner association."""
from pathlib import Path
from types import SimpleNamespace
import hashlib
import importlib.util
import json
import subprocess

root = Path.cwd()
work = Path(__file__).resolve().parent
folder = work / 'control-8'


def doc(path):
    return json.loads(path.read_text())


def pin(path):
    data = path.read_bytes()
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}


assert not (work / 'control-8-acceptance.json').exists()
action = doc(work / 'actions.json')['control-8']
owner = doc(folder / 'owner.json')
assert owner['passed'] and owner['status'] == 'terminal' and owner['child_exit'] == 0
assert owner['source'] == {'revision': action['revision'], 'status': ''}
assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
assert owner['limits'] == {'seconds': 300, 'rss_kib': 4194304, 'reserve_kib': 2097152}
assert owner['stop_reason'] is None and not owner['remaining_owned_processes'] and not owner['cleanup_errors']
assert owner['wall_seconds'] <= 300 and owner['peak_sampled_rss_kib'] <= 4194304
assert owner['minimum_sampled_available_kib'] >= 2097152
assert owner['observed_births'] and all(row['absent'] for row in owner['observed_births'])
before = doc(folder / 'inputs-before.json')
assert before == doc(folder / 'inputs-after.json')
for name, expected in before.items():
    assert pin(root / name) == {key: expected[key] for key in ('bytes', 'sha256')}, name
environment = owner['declared_output_environment']
assert environment == {'GTEST_OUTPUT': 'xml:' + str(folder / 'gtest.xml'),
    'VK_DRIVER_FILES': action['runtime_icd'], 'SIRIUS_VULKAN_DEVICE': '0',
    'MESA_SHADER_CACHE_DISABLE': 'true'}
assert owner['provider_mappings'] and all(action['runtime_provider'] in line
    for row in owner['provider_mappings'] for line in row['maps'])
observation = doc(folder / 'native-observation.json')
assert not observation['failures'] and observation['inferior_exit'] == 0
assert len(observation['calls']) == 3
assert {row['function'] for row in observation['calls']} == {
    'sirius_native_' + slug + '_try_packet' for slug in ('transport', 'endpoint', 'dopri_phase')}
symbols = {}
for line in subprocess.check_output(['/usr/bin/nm', '-a', action['runtime_provider']], text=True).splitlines():
    words = line.split()
    if len(words) == 3 and words[2].startswith('sirius_native_') and words[2].endswith('_try_packet'):
        symbols[words[2]] = int(words[0], 16)
for row in observation['calls']:
    assert row['returned'] and row['handled'] and row['library'] == action['runtime_provider']
    assert row['thread'] == row['return_thread']
    assert row['mxcsr'] & 0xe040 == row['return_mxcsr'] & 0xe040 == 0
    bases = {int(line.split()[0].split('-')[0], 16) - int(line.split()[2], 16)
             for line in row['provider_maps']}
    assert row['pc'] - symbols[row['function']] in bases

source = root / 'out/native-cpu-arithmetic-sharing/readers/native_evidence.py'
assert pin(source)['sha256'] == '1b44f0678dbfc7ee455b3d49e170be884b294882896fd1ac3706ff5ebded20ab'
spec = importlib.util.spec_from_file_location('original_timestamp_reader', source)
old = importlib.util.module_from_spec(spec)
spec.loader.exec_module(old)
old.WORK = source.parent
names = dict(zip(old.STAGES, ['k' + prefix + 'Shader' for prefix in old.PREFIXES]))
for stage in ('transport', 'endpoint', 'dopri_phase'):
    names[stage] = 'k' + dict(zip(old.STAGES, old.PREFIXES))[stage] + 'NativeOwnedShader'
old.selected_names = lambda context, arm, wide=False: names
modules = doc(root / 'out/native-cpu-row-ownership/compiler-result.json')['modules']
arrays = {}
binary = Path(action['executable']).read_bytes()
for stage, name in names.items():
    path = (root / modules[{'transport': 'Transport', 'endpoint': 'Endpoint', 'dopri_phase': 'DopriPhase'}[stage]]['module']
            if stage in ('transport', 'endpoint', 'dopri_phase') else
            root / 'bin/linux-gcc/src/sirius/backend/retained' / ('retained_' + stage + '.spv'))
    assert binary.count(path.read_bytes()) == 1
    arrays[name] = pin(path)
history = root / 'out/native-original-input-ownership/provider-history/candidate-default/stdout.log'
info, = [json.loads(line) for line in history.read_text().splitlines()
         if line.startswith('{') and 'resident_driver_path' in line]
creation, = doc(root / 'out/native-cpu-slang-bridge/ordinary/transport/device-features.json')['device_creations']
properties = doc(root / 'out/native-cpu-slang-bridge/ordinary/transport/float-controls-properties.json')[0]['values']
fma = any(item['fma16_32_64'][1] == 1 for item in creation['feature_chain'])
hardware = {'selected_device': info['device'], 'selected_driver': info['driver'] + ': ' + info['driver_info'],
    'retained_scalar_route': 'native binary32', 'timestamp_valid_bits': '64', 'timestamp_period_ns': '1',
    'fma_fp32_enabled': str(int(fma)),
    'preserves_fp32_signed_zero_inf_nan': str(properties['shaderSignedZeroInfNanPreserveFloat32'])}
context = SimpleNamespace(value={'native_identity': hardware}, arms={'candidate': {'arrays': arrays}})
record = old.timestamp_result(old.xml_properties(folder / 'gtest.xml', 'control-8'), 'candidate', context)
reference = doc(root / 'out/native-cpu-row-ownership/previous-control-8-acceptance.json')['record']['readbacks']
assert record['readbacks'] == reference
result = {'passed': True, 'all32_readbacks_match_preserved_originals': True,
    'native_first_packets_handled': 3, 'native_observation': pin(folder / 'native-observation.json'),
    'owner': pin(folder / 'owner.json'), 'record': record,
    'scope': 'Sole original C8, three actual first native packets returnedtrue,51markers54modulejoins32complete-bufferdigests exact. Passive debugger timing is not benefit evidence. No ABBA/all8/frame/scientific/release/adoption claim.'}
(work / 'control-8-acceptance.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps({'passed': True, 'complete_readbacks': 32, 'native_stages_handled': 3}))
