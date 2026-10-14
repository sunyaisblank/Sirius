"""Pin prepared private actions to the operator's current clean checkout."""
from pathlib import Path
import hashlib
import json
import subprocess
W = Path(__file__).resolve().parent
R = Path.cwd().resolve()
assert not subprocess.check_output(['git', 'status', '--porcelain=v1'])
assert not (R / 'bin/linux-gcc/computed-carry-exit-original-prerequisites').exists()
assert not (W / 'operator-context.json').exists()
def identity(path):
    with Path(path).open('rb') as stream:
        digest = hashlib.file_digest(stream, 'sha256').hexdigest()
    return {'bytes': Path(path).stat().st_size, 'sha256': digest}
for path, pin in json.loads((W / 'source-freeze.json').read_text()).items():
    assert identity(path) == pin, path
revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
for folder in [W, W / 'arithmetic']:
    actions = json.loads((folder / 'actions.json').read_text())
    for name, action in actions.items():
        assert not (folder / name).exists() and action['revision'] == 'PENDING_OPERATOR_CLEAN_CONTEXT'
        action['revision'] = revision
    (folder / 'actions.json').write_text(json.dumps(actions, indent=2) + '\n')
(W / 'operator-context.json').write_text(json.dumps({'revision': revision, 'status': '',
    'source_freeze': identity(W / 'source-freeze.json'),
    'initial_actions': identity(W / 'actions.json'),
    'arithmetic_actions': identity(W / 'arithmetic/actions.json')}, indent=2) + '\n')
print('Pinned original-prerequisite actions to', revision)
