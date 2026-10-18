"""Run the six source-frozen matched C++ producer commands, once."""
from pathlib import Path
import hashlib
import json
import shutil
import subprocess
import sys

sys.dont_write_bytecode = True
W = Path(__file__).resolve().parent
R = W.parents[1]
T = R / 'bin/linux-gcc/packed-retained-scratch-controls'
assert Path.cwd().resolve() == R
recipes = json.loads((W / 'recipes.pending.json').read_text())
labels = ['compile-baseline-compute', 'compile-baseline-observer', 'link-baseline',
          'compile-candidate-compute', 'compile-candidate-observer', 'link-candidate']
assert set(recipes['producer_actions']) == set(labels)
assert not T.exists()
T.mkdir()
record_path = T / 'producer-command-results.json'
record = {'schema': 'sirius-packed-scratch-matched-command-results-v1',
          'source_revision': recipes['producer_actions'][labels[0]]['revision'],
          'passed': False, 'commands': [], 'scope': 'C++ producer only; no acceptance or runtime.'}
def write():
    temporary = record_path.with_suffix('.tmp')
    temporary.write_text(json.dumps(record, indent=2) + '\n')
    temporary.replace(record_path)
def identity(path):
    data = path.read_bytes()
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}
write()
shutil.copytree(R / 'bin/linux-gcc/tests/backend/resources', T / 'resources')
for label in labels:
    action = recipes['producer_actions'][label]
    assert action['kind'] == 'diagnostic' and action['expected_exit'] == 0
    assert action['limits'] == {'seconds': 300, 'rss_kib': 4194304, 'reserve_kib': 2097152}
    row = {'label': label, 'argv': action['argv'], 'started': True, 'exit': None}
    record['commands'].append(row)
    write()
    print(json.dumps({'starting': label, 'argv': action['argv']}), flush=True)
    row['exit'] = subprocess.run(action['argv'], cwd=R).returncode
    write()
    if row['exit'] != 0:
        raise SystemExit(row['exit'])
outputs = [T / (arm + '-' + kind + suffix)
           for arm in ('baseline', 'candidate') for kind in ('compute', 'observer')
           for suffix in ('.o', '.d')]
outputs += [T / (arm + '.map') for arm in ('baseline', 'candidate')]
outputs += [T / ('original_' + arm + '_controls') for arm in ('baseline', 'candidate')]
assert len(outputs) == 12 and all(path.is_file() for path in outputs)
record['products'] = {str(path.relative_to(R)): identity(path) for path in outputs}
record['resources'] = {str(path.relative_to(R)): identity(path)
                       for path in sorted((T / 'resources').rglob('*')) if path.is_file()}
record['passed'] = True
write()
print(json.dumps({'producer_commands_passed': 6, 'products': 12,
                  'runtime_executed': False}), flush=True)
