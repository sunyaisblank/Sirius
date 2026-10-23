"""Identify the changed outputs once; seal all reused private build inputs."""
from pathlib import Path
import hashlib
import json
import re
import shlex
import subprocess

R = Path.cwd()
W = Path(__file__).resolve().parent
P = W / 'provider'
S = R / 'out/lavapipe-acqrel-ordering/upstream/mesa-25.2.8'
B = R / 'bin/linux-gcc/lavapipe-acqrel-ordering'


def pin(path):
    return {'bytes': path.stat().st_size,
            'sha256': hashlib.file_digest(path.open('rb'), 'sha256').hexdigest()}


assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
assert json.loads((P / 'configure/owner.json').read_text())['passed']
assert json.loads((P / 'configure-result.json').read_text())['exit'] == 0
actions = json.loads((P / 'actions.json').read_text())
assert 'build' not in actions
dry = subprocess.check_output(['ninja', '-C', str(B), '-n', '-v',
                              'src/gallium/targets/lavapipe/libvulkan_lvp.so'], text=True)
commands = [match[1] for line in dry.splitlines()
            if (match := re.match(r'^\[\d+/\d+\] (.*)$', line))]
assert commands
mutable = set()
for command in commands:
    words = shlex.split(command)
    outputs = {words[i+1] for i, word in enumerate(words[:-1])
               if word in {'-o', '--output', 'csr', 'csrD', 'csrDT', 'crs', 'crsD', 'rcs', 'rcsD'}}
    assert outputs, ('Review unknown build command before executing', command)
    for output in outputs:
        path = (B / output).resolve()
        assert path.is_relative_to(B)
        mutable.add(str(path.relative_to(R)))
target = B / 'src/gallium/targets/lavapipe/libvulkan_lvp.so'
assert str(target.relative_to(R)) in mutable
reused = {str(path.relative_to(R)): pin(path) for path in B.rglob('*') if path.is_file()
          and str(path.relative_to(R)) not in mutable
          and path.name not in {'.ninja_log', '.ninja_deps', '.ninja_lock'}}
plan = {'dry_run_commands': commands, 'mutable_outputs': sorted(mutable),
        'reused_build_inputs': reused, 'previous_target': pin(target),
        'configuration_owner': pin(P / 'configure/owner.json')}
assert not (P / 'incremental-plan.json').exists()
(P / 'incremental-plan.json').write_text(json.dumps(plan, indent=2) + '\n')
(P / 'incremental-dry-run.txt').write_text(dry)
config = actions['configure']
inputs = {name for name in config['frozen_inputs'] if not (R / name).is_relative_to(B)}
inputs |= set(reused)
inputs |= {str(path.relative_to(R)) for path in P.rglob('*') if path.is_file()
           and path.name != 'actions.json'}
inputs |= {str(W / 'build_provider.py'), str(W / 'prepare_provider_build.py')}
action = {**config,
          'revision': subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
          'kind': 'build', 'limits': {'seconds': 300, 'rss_kib': 4194304, 'reserve_kib': 2097152},
          'argv': ['/usr/bin/python3', '-B', str(W / 'build_provider.py')],
          'frozen_inputs': sorted(inputs),
          'scope': 'One incremental private Mesa provider build linking3 sealed same-source native CPU objects; same configuration/toolchain and predeclared dirty outputs only. No provider loading/numerical/performance/science/release acceptance.'}
actions['build'] = action
(P / 'actions.json').write_text(json.dumps(actions, indent=2) + '\n')
print(json.dumps({'pending_commands': len(commands), 'mutable_outputs': len(mutable),
                  'sealed_reused_files': len(reused), 'sealed_inputs': len(inputs)}))
