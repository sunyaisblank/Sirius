"""Seal the narrow read-only completion or first workload after completed proof."""
from pathlib import Path
import hashlib
import json
import subprocess
import sys

R = Path.cwd()
W = Path(__file__).resolve().parent


def pin(path):
    with path.open('rb') as stream:
        return {'bytes': path.stat().st_size, 'sha256': hashlib.file_digest(stream, 'sha256').hexdigest()}


mode, = sys.argv[1:]
assert mode in ('proof', 'run')
assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
actions = json.loads((W / 'actions.json').read_text())
prepared = json.loads((W / 'prepared.json').read_text())
inputs = dict(prepared['source_inputs'])
paths = [W / 'prepared.json', W / 'build-result.json', W / 'build/owner.json',
         W / 'complete_build.py', W / 'elf_proof.py', W / 'prepare_completion.py', W / 'run_completed.py',
         R / prepared['executable'], (R / prepared['executable']).parent / 'shim.o']
if mode == 'run':
    owner = json.loads((W / 'proof/owner.json').read_text())
    build = json.loads((W / 'build-completion.json').read_text())
    assert owner['passed'] and owner['whole_inputs_unchanged'] and owner['source_unchanged']
    assert not owner['remaining_owned_processes'] and not owner['cleanup_errors'] and build['passed']
    inputs.update(build['runtime_libraries'])
    paths += [W / 'proof/owner.json', W / 'build-completion.json', W / 'link-proof.json']
for path in paths:
    inputs[str(path.relative_to(R))] = pin(path)
for name, expected in inputs.items():
    assert pin(R / name) == expected, name
assert mode not in actions and not (W / mode).exists()
actions[mode] = {**actions['build'],
    'revision': subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
    'argv': ['/usr/bin/python3', '-B', str(W / ('complete_build.py' if mode == 'proof' else 'run_completed.py')), mode, '{folder}'],
    'frozen_inputs': sorted(inputs),
    'scope': ('Read-only exact code association of already successfully compiled/linked new field Transport. No compile, link or workload.' if mode == 'proof' else prepared['scope'])}
(W / 'actions.json').write_text(json.dumps(actions, indent=2) + '\n')
print(json.dumps({'prepared': mode, 'compilers_or_workload_executed': False}))
