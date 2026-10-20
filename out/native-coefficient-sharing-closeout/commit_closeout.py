import hashlib
import json
import os
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
WORK = Path(__file__).resolve().parent
PARENT = 'b9a5703a18f06070fe809a9bfd14a414c4f2f906'
BRANCH = 'refs/heads/development/native-fp32-packed-scratch-runtime-2026-10-18'


def run(*args, data=None, env=None):
    return subprocess.run(args, cwd=ROOT, input=data, stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE, check=True, env=env).stdout


assert not run('git', 'status', '--porcelain')
assert run('git', 'rev-parse', BRANCH).decode().strip() == PARENT
assert json.loads((WORK / 'retirement.json').read_text())['passed']
assert json.loads((WORK / 'published-input-retirement.json').read_text())['passed']
paths = sorted(path for source in [WORK, ROOT / 'out/native-coefficient-sharing-preservation']
               for path in source.rglob('*') if path.is_file())
assert all(not path.is_symlink() for path in paths)
index = WORK / 'closeout.index'
assert not index.exists()
env = dict(os.environ, GIT_INDEX_FILE=str(index))
records = {}
try:
    run('git', 'read-tree', PARENT, env=env)
    rows = []
    for path in paths:
        data = path.read_bytes()
        blob = run('git', 'hash-object', '-w', '--stdin', data=data).decode().strip()
        mode = '100755' if path.stat().st_mode & 0o111 else '100644'
        name = str(path.relative_to(ROOT))
        records[name] = {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest(),
                         'Git_mode': mode, 'mode': oct(path.stat().st_mode & 0o7777)}
        rows.append(f'{mode} {blob}\t{name}\0'.encode())
    run('git', 'update-index', '-z', '--index-info', data=b''.join(rows), env=env)
    tree = run('git', 'write-tree', env=env).decode().strip()
    commit = run('git', 'commit-tree', tree, '-p', PARENT,
                 data=b'Conserve coefficient-sharing rejection cleanup and exact issue/PR status\n').decode().strip()
    for name, row in records.items():
        data = run('git', 'show', commit + ':' + name)
        assert len(data) == row['bytes'] and hashlib.sha256(data).hexdigest() == row['sha256']
    run('git', 'update-ref', BRANCH, commit, PARENT)
    receipt = {'schema': 'sirius-native-coefficient-sharing-closeout-Git-v1', 'commit': commit,
               'parent': PARENT, 'branch': BRANCH, 'files': records,
               'full_Git_bytes_verified': True, 'pushed': False}
    (WORK / 'root-Git.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({key: value for key, value in receipt.items() if key != 'files'}))
finally:
    index.unlink(missing_ok=True)
