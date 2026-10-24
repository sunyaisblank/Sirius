"""Seal the original remaining7 only after the sole original22+Decimal pass."""
from pathlib import Path
import hashlib
import json
import subprocess

root = Path.cwd()
work = Path(__file__).resolve().parent
current = work.parent
runtime = work / 'runtime'
benefit = current / 'benefit'


def doc(path):
    return json.loads(path.read_text())


def pin(path):
    data = path.read_bytes()
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}


assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
assert not (runtime / 'actions.json').exists()
assert not any((runtime / ('control-' + str(i))).exists() for i in range(1, 8))
assert doc(current / 'control-8-acceptance.json')['passed']
accepted = doc(benefit / 'readers/benefit.json')
decimal = doc(benefit / 'readers/decimal-agreement.json')
assert accepted['audit_passed'] and accepted['all22_numeric_pass'] and accepted['checks_passed'] == 22
assert decimal['agrees_on_all22'] and decimal['all22_pass']
providers = doc(benefit / 'providers.json')['providers']
provider = providers['candidate']
for key in ('library', 'ICD'):
    assert pin(Path(provider[key])) == provider[key + '_identity']
cohort = doc(benefit / 'runtime/actions.json')
candidate = cohort['b2']
executable = Path(candidate['executable'])
assert pin(executable)['sha256'] == '7f2d1fa88953b954ab116233458a4cb962470f863fda9bbb3a1cad7d4a7d6a26'
inputs = set(candidate['frozen_inputs'])
inputs |= {str(p.relative_to(root)) for p in work.rglob('*') if p.is_file()}
inputs |= {str(benefit / 'readers/benefit.json'), str(benefit / 'readers/decimal-agreement.json'),
           str(benefit / 'runtime/actions.json'), str(benefit / 'providers.json'),
           str(current / 'control-8-acceptance.json'),
           'out/native-original-input-ownership/readers/read_controls.py',
           'out/native-cpu-arithmetic-sharing/readers/native_evidence.py'}
for label in ('a1', 'b1', 'b2', 'a2'):
    owner = doc(benefit / 'runtime' / label / 'owner.json')
    assert owner['passed'] and owner['status'] == 'terminal'
    assert not owner['remaining_owned_processes'] and not owner['cleanup_errors']
    for name in ('owner.json', 'gtest.xml'):
        inputs.add(str((benefit / 'runtime' / label / name).relative_to(root)))
assert all((root / name).is_file() for name in inputs)
expected = {str(executable): pin(executable), provider['library']: provider['library_identity'],
            provider['ICD']: provider['ICD_identity']}
original = doc(root / 'out/native-original-input-ownership/runtime/actions.json')
revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
actions = {}
for i in range(1, 8):
    label = 'control-' + str(i)
    old = original[label]
    action = {key: old[key] for key in ('arm', 'case', 'prerequisites', 'kind', 'expected_exit', 'limits')}
    assert action['arm'] == 'candidate'
    assert action['limits'] == {'seconds': 300, 'rss_kib': 4194304, 'reserve_kib': 2097152}
    action.update(revision=revision,
        argv=['/usr/bin/python3', '-B', str(runtime / 'run_case.py'), label, '{folder}'],
        frozen_inputs=sorted(inputs), expected_identities=expected, executable=str(executable),
        software_provider=True, runtime_icd=provider['ICD'], runtime_provider=provider['library'],
        scope='Original ' + label + ' only after current C8 and soleABBA22+Decimal PASS. Same exact7f2host/Nprovider, originalcases/inputs/allocations/cache/governors/limits/assertions. OriginalC1 retains two actualprecisionroutes,28complete readbacks and wholeinactiveEndpointallocation comparisons. No retry or scientific/frame/release claim.')
    actions[label] = action
    if i == 4:
        action['scope'] += (' Original C4 runs once under a passive debugger with bounded first-entry/return records and one fixed45s stack sample if still running. All observation overhead remains inside the original300s/RSS envelope; no performance-benefit evidence.')
        action['frozen_inputs'] = sorted(set(action['frozen_inputs']) | {'/usr/bin/gdb'})
(runtime / 'actions.json').write_text(json.dumps(actions, indent=2) + '\n')
print(json.dumps({'original_controls_prepared': 7, 'executed': False, 'frozen_inputs': len(inputs)}))
