import hashlib
import json
import os
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
WORK = Path(__file__).resolve().parent
PARENT = '1cb6e6aaff04ddd54ff7ee133b773257338f2c52'
BRANCH = 'refs/heads/development/native-fp32-packed-scratch-runtime-2026-10-18'
TAG = 'refs/tags/evidence/native-wave32-rejected-2c8ef086'
BUNDLE = ROOT / 'attestations/qualification/2c8ef086/native-wave32-controls-benefit-outcome'


def run(*argv, data=None, env=None):
    return subprocess.run(argv, cwd=ROOT, input=data, stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE, check=True, env=env).stdout


assert run('git', 'rev-parse', 'HEAD').decode().strip() == '2c8ef08693d46f57d5dea14efdbddd4874dac684'
assert not run('git', 'status', '--porcelain')
assert run('git', 'rev-parse', BRANCH).decode().strip() == PARENT
manifest = json.loads((BUNDLE / 'MANIFEST.json').read_text())
preservation = json.loads((BUNDLE / 'PRESERVATION.json').read_text())
assert preservation['outcome']['original22_benefit_passed'] is False
assert preservation['outcome']['checks_passed'] == 20 and preservation['payload_count'] == 199
paths = sorted(path for path in BUNDLE.rglob('*') if path.is_file())
assert len(paths) == 201
paths += [Path(__file__).resolve(), WORK / 'root-preparation.json']
index = WORK / 'preservation.index'
assert not index.exists()
env = dict(os.environ, GIT_INDEX_FILE=str(index))
records = {}
try:
    run('git', 'read-tree', PARENT, env=env)
    rows = []
    for path in paths:
        assert not path.is_symlink()
        data = path.read_bytes()
        blob = run('git', 'hash-object', '-w', '--stdin', data=data).decode().strip()
        mode = '100755' if path.stat().st_mode & 0o111 else '100644'
        name = str(path.relative_to(ROOT))
        records[name] = {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest(),
                         'Git_mode': mode, 'physical_mode': oct(path.stat().st_mode & 0o7777), 'blob': blob}
        rows.append(f'{mode} {blob}\t{name}\0'.encode())
    run('git', 'update-index', '-z', '--index-info', data=b''.join(rows), env=env)
    tree = run('git', 'write-tree', env=env).decode().strip()
    commit = run('git', 'commit-tree', tree, '-p', PARENT,
                 data=b'Preserve native Transport wave32 original benefit rejection\n').decode().strip()
    changes = run('git', 'diff-tree', '--no-commit-id', '--name-only', '-r', commit).decode().splitlines()
    assert set(changes) == set(records)
    for name, row in records.items():
        data = run('git', 'show', commit + ':' + name)
        assert len(data) == row['bytes'] and hashlib.sha256(data).hexdigest() == row['sha256']
    run('git', 'update-ref', BRANCH, commit, PARENT)
    run('git', 'update-ref', TAG, commit, '0' * 40)
    receipt = {'schema': 'sirius-native-wave32-rejection-Git-v1', 'commit': commit,
               'parent': PARENT, 'branch': BRANCH, 'tag': TAG,
               'bundle': str(BUNDLE.relative_to(ROOT)), 'files': records,
               'full_stored_Git_bytes_verified': True, 'Git_exec_modes_verified': True,
               'full_POSIX_modes_restored_by_manifest': True, 'pushed': False}
    (WORK / 'root-Git.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({key: value for key, value in receipt.items() if key != 'files'}))
finally:
    index.unlink(missing_ok=True)
