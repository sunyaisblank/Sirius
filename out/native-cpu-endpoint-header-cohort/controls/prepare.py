"""Seal the original seven controls after C8; leave final all8+22 admission pending."""
from pathlib import Path
import hashlib
import json
import subprocess
from control_authority import ao, authority_files

root = Path.cwd()
work = Path(__file__).resolve().parent
current = work.parent
runtime = work / 'runtime'
provider_work = root / 'out/native-cpu-endpoint-header-provider/provider'


def doc(path):
    return json.loads(path.read_text())


def pin(path):
    with path.open('rb') as stream:
        return dict(bytes=path.stat().st_size, sha256=hashlib.file_digest(stream, 'sha256').hexdigest())


ao('prepare.py')
assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
assert not (runtime / 'actions.json').exists()
assert not any((runtime / ('control-' + str(i))).exists() for i in range(1, 8))
accepted = doc(current / 'control-8-acceptance.json')
assert accepted['passed'] and accepted['optional_math']['all_selected_optional_math_observed']
owner = doc(current / 'control-8/owner.json')
assert owner['passed'] and owner['status'] == 'terminal' and owner['child_exit'] == 0
assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
assert not owner['remaining_owned_processes'] and not owner['cleanup_errors']
candidate = doc(current / 'actions.json')['control-8']
build = doc(provider_work / 'provider-build-result.json')
contract = doc(provider_work / 'source-contract.json')
proof = doc(provider_work / 'linked-native-proof.json')
assert build['passed'] and proof['passed']
assert build['native_objects'] == contract['native_objects'] == [row['object'] for row in proof['objects']]
assert len(build['native_objects']) == 8
provider = dict(library=candidate['runtime_provider'], ICD=candidate['runtime_icd'])
for key in ('library', 'ICD'):
    assert provider[key] == str(root / build[key])
    provider[key + '_identity'] = build[key + '_identity']
    assert pin(Path(provider[key])) == provider[key + '_identity']
executable = Path(candidate['executable'])
host = doc(root / 'out/native-cpu-owned-lower-increments/host-build-result.json')
assert host['passed'] and str(executable) == str(root / host['executable'])
assert pin(executable) == host['executable_identity']
inputs = set(candidate['frozen_inputs']) | set(contract['source_inputs'])
inputs |= {str(path.relative_to(root)) for path in work.rglob('*') if path.is_file()}
inputs |= {str(path.relative_to(root)) for path in authority_files()}
inputs |= {'out/native-original-input-ownership/readers/read_controls.py',
           'out/native-cpu-arithmetic-sharing/readers/native_evidence.py'}
for path in (current / 'actions.json', current / 'control-8-acceptance.json',
             current / 'control-8/owner.json', current / 'control-8/gtest.xml',
             current / 'control-8/native-observation.json',
             provider_work / 'source-contract.json', provider_work / 'provider-build-result.json',
             provider_work / 'linked-native-proof.json'):
    inputs.add(str(path.relative_to(root)))
assert all((root / name).is_file() for name in inputs)
expected = dict(candidate['expected_identities'])
expected.update({str(executable): host['executable_identity'],
                 provider['library']: provider['library_identity'], provider['ICD']: provider['ICD_identity']})
for name, identity in expected.items():
    assert pin(root / name) == identity, name
original = doc(root / 'out/native-original-input-ownership/runtime/actions.json')
revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
actions = {}
for i in range(1, 8):
    label = 'control-' + str(i)
    old = original[label]
    action = {key: old[key] for key in ('arm', 'case', 'prerequisites', 'kind', 'expected_exit', 'limits')}
    assert action['arm'] == 'candidate'
    assert action['limits'] == dict(seconds=300, rss_kib=4194304, reserve_kib=2097152)
    action['prerequisites'] = [name for name in action['prerequisites'] if name not in ('a1', 'b1', 'b2', 'a2')]
    assert action['prerequisites'] == ['control-8'] + ['control-' + str(n) for n in range(1, i)]
    action.update(revision=revision,
        argv=['/usr/bin/python3', '-B', str(runtime / 'run_case.py'), label, '{folder}'],
        frozen_inputs=sorted(inputs), expected_identities=expected, executable=str(executable),
        software_provider=True, runtime_icd=provider['ICD'], runtime_provider=provider['library'],
        scope='Original ' + label + ' after current AP C8, once in original control order before timing study. Same genuine H host, AP Endpoint header admission provider, original cases/inputs/allocations/cache/governors/300s4GiB2GiB/assertions. C1 retains two actual precision routes,28 full readbacks and whole inactive Endpoint allocation comparisons. Only ABBA scheduling dependencies are deferred; the final all8+22 reader still requires fresh ABBA and independent Decimal acceptance. No prior timing acceptance, retry, frame/scientific/release claim.')
    actions[label] = action
(runtime / 'actions.json').write_text(json.dumps(actions, indent=2) + '\n')
print(json.dumps(dict(original_controls_prepared=7, executed=False, frozen_inputs=len(inputs),
                     final_all8_plus22_reader_pending=True)))
