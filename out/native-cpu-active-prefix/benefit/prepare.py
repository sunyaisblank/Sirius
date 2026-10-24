"""Prepare exactly the original ABBA after current native B C8 admission.

No build, provider load, warmup, test or benefit execution occurs here.
"""
from pathlib import Path
import hashlib
import json
import re
import subprocess

root = Path.cwd()
work = Path(__file__).resolve().parent
current = work.parent
runtime = work / 'runtime'


def doc(path):
    return json.loads(path.read_text())


def pin(path):
    data = path.read_bytes()
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}


assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
assert not (runtime / 'actions.json').exists()
assert not any((runtime / label).exists() for label in ('a1', 'b1', 'b2', 'a2'))
assert doc(current / 'control-8-acceptance.json')['passed']
control_action = doc(current / 'actions.json')['control-8']
control_owner = doc(current / 'control-8/owner.json')
assert control_owner['passed'] and control_owner['source_unchanged'] and control_owner['whole_inputs_unchanged']
assert doc(current / 'provider/linked-native-proof.json')['passed']
build = doc(current / 'provider/provider-build-result.json')
assert build['passed']
baseline = doc(root / 'out/native-cpu-arithmetic-sharing/providers.json')['providers']['baseline']
for key in ('library', 'ICD'):
    assert pin(Path(baseline[key])) == baseline[key + '_identity']
assert baseline['library_identity']['sha256'] == 'cde2d0f8f243448a3f3b28408121412fb03a8bd9c01f5731f9e87a107150e952'
candidate = {'library': str(root / build['library']), 'ICD': str(root / build['ICD']),
    'library_identity': build['library_identity'], 'ICD_identity': build['ICD_identity'],
    'historical_library': '/home/astra/.project/Sirius/bin/linux-gcc/lavapipe-ieee32-staged/candidate/libvulkan_lvp.so',
    'scope': 'Historical unchanged provider properties supply expectations; current B identity and successful native admission are separately bound by original C8. No new property query is claimed.'}
providers = {'baseline': baseline, 'candidate': candidate}
for value in providers.values():
    assert pin(Path(value['library'])) == value['library_identity']
    assert pin(Path(value['ICD'])) == value['ICD_identity']
(work / 'providers.json').write_text(json.dumps({'providers': providers}, indent=2) + '\n')
executable = Path(control_action['executable'])
assert pin(executable)['sha256'] == '7f2d1fa88953b954ab116233458a4cb962470f863fda9bbb3a1cad7d4a7d6a26'
inputs = set(control_action['frozen_inputs'])
inputs |= {str(p.relative_to(root)) for p in work.rglob('*') if p.is_file()}
inputs |= {str(current / 'control-8-acceptance.json'), str(current / 'control-8/owner.json'),
           str(current / 'control-8/gtest.xml'), str(current / 'control-8/native-observation.json'),
           'out/native-original-input-ownership/readers/read_results.py',
           'out/native-original-input-ownership/readers/decimal_check.py',
           'out/native-original-input-ownership/readers/run_abba.py',
           'out/native-cpu-arithmetic-sharing/readers/native_evidence.py',
           'out/native-cpu-arithmetic-sharing/readers/numeric_gates.py',
           'out/native-cpu-row-ownership/compiler-result.json'}
history = root / 'out/native-original-input-ownership/provider-history'
for label in ('baseline-poison-v21', 'candidate-default'):
    for name in ('stdout.log', 'float-controls-properties.json', 'device-features.json'):
        inputs.add(str((history / label / name).relative_to(root)))
for value in providers.values():
    inputs |= {value['library'], value['ICD']}
    for name in re.findall(r'(/[^\s()]+)', subprocess.check_output(['/usr/bin/ldd', value['library']], text=True)):
        if Path(name).is_file():
            inputs.add(str(Path(name).resolve(strict=True)))
assert all((root / name).is_file() for name in inputs)
expected = {str(executable): pin(executable)}
for value in providers.values():
    expected[value['library']] = value['library_identity']
    expected[value['ICD']] = value['ICD_identity']
original = doc(root / 'out/native-original-input-ownership/runtime/actions.json')
revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
actions = {'control-8': control_action}
for label in ('a1', 'b1', 'b2', 'a2'):
    old = original[label]
    action = {key: old[key] for key in ('arm', 'case', 'prerequisites', 'kind', 'expected_exit', 'limits')}
    value = providers[action['arm']]
    assert action['limits'] == {'seconds': 300, 'rss_kib': 4194304, 'reserve_kib': 2097152}
    action.update(revision=revision,
        argv=['/usr/bin/python3', '-B', str(runtime / 'run_case.py'), label, '{folder}'],
        frozen_inputs=sorted(inputs), expected_identities=expected, executable=str(executable),
        software_provider=True, runtime_icd=value['ICD'], runtime_provider=value['library'],
        scope='Original sole ABBA22 slot ' + label + ': same exact7f2d1fa8 host, actual registered stock baseline and repaired native B candidate. Original3internalrepeats, all32readbacks/51markers/54joins, allocations/counts/cache-disable/governors/300s4GiB2GiB and assertions unchanged. No debugger, warmup, resampling, pooling or cost subtraction. Fourwhole>=5%,16inclusive-stages>=-5%,2provisioning-inclusiveAPI>=0 with independentDecimal; remaining7/frame only after22PASS.')
    actions[label] = action
(runtime / 'actions.json').write_text(json.dumps(actions, indent=2) + '\n')
print(json.dumps({'actions': 4, 'order': ['a1', 'b1', 'b2', 'a2'], 'frozen_inputs': len(inputs),
                  'same_executable_sha256': expected[str(executable)]['sha256'], 'executed': False}))
