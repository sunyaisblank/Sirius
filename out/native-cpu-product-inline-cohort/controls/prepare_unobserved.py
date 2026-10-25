"""Complete only C4-C7 after the copied C4 observer harness aborted."""
from pathlib import Path
import hashlib
import json
import subprocess

R = Path.cwd(); W = Path(__file__).resolve().parent; runtime = W / 'runtime-unobserved'
doc = lambda p: json.loads(p.read_text())
pin = lambda p: {'bytes': p.stat().st_size, 'sha256': hashlib.file_digest(p.open('rb'), 'sha256').hexdigest()}
assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
old = doc(W / 'runtime/actions.json')
for i in range(1, 4):
    owner = doc(W / 'runtime' / ('control-' + str(i)) / 'owner.json')
    assert owner['passed'] and owner['source_unchanged'] and owner['whole_inputs_unchanged']
failed = W / 'runtime/control-4'
owner = doc(failed / 'owner.json')
assert not owner['passed'] and owner['child_exit'] == 1 and owner['stop_reason'] is None
assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
assert not owner['remaining_owned_processes'] and not owner['cleanup_errors']
assert all(row['absent'] for row in owner['observed_births'])
assert not (failed / 'gtest.xml').exists()
assert 'control4_observer.py: No such file' in (failed / 'stderr.log').read_text()
assert not any((W / 'runtime' / ('control-' + str(i))).exists() for i in range(5, 8))
assert not (runtime / 'actions.json').exists()
inputs = {name for a in old.values() for name in a['frozen_inputs']}
inputs |= {str(p.relative_to(R)) for p in failed.rglob('*') if p.is_file()}
for i in range(1, 4):
    inputs |= {str(p.relative_to(R)) for p in (W / 'runtime' / ('control-' + str(i))).rglob('*') if p.is_file()}
for root in (runtime, W / 'readers-unobserved'):
    inputs |= {str(p.relative_to(R)) for p in root.rglob('*') if p.is_file()}
inputs.add(str(Path(__file__).resolve().relative_to(R)))
assert all((R / name).is_file() for name in inputs)
revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
actions = {}
for i in range(4, 8):
    label = 'control-' + str(i)
    action = dict(old[label])
    action.update(revision=revision, frozen_inputs=sorted(inputs),
        argv=['/usr/bin/python3', '-B', str(runtime / 'run_case.py'), label, '{folder}'],
        scope='Original C' + str(i) + ' uninstrumented after the earlier C4 harness aborted for its missing observer before any XML/verdict. Original C1-C3 accepted executions remain; no repeat. Same W provider/host/inputs/assertions/cache/governors. Earlier aborted prefix preserved separately; no inferred numerical rejection or qualification from it.')
    actions[label] = action
(runtime / 'actions.json').write_text(json.dumps(actions, indent=2) + '\n')
print(json.dumps({'remaining_uninstrumented_controls': [4, 5, 6, 7], 'earlier_C4_abort_preserved': True}))
