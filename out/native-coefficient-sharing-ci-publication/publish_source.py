import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
WORK = Path(__file__).resolve().parent
BASE = '2c8ef08693d46f57d5dea14efdbddd4874dac684'
BRANCH = 'refs/heads/development/native-coefficient-sharing-diagnostic-2026-10-20'


def run(*argv, data=None, env=None):
    return subprocess.run(argv, cwd=ROOT, input=data, stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE, env=env, check=True).stdout


assert run('git', 'rev-parse', 'HEAD').decode().strip() == BASE
assert not run('git', 'status', '--porcelain')
assert subprocess.run(['git', 'show-ref', '--verify', '--quiet', BRANCH], cwd=ROOT).returncode == 1
assert len(sys.argv) == 2
admission = ROOT / sys.argv[1]
assert admission.resolve(strict=True) == admission and admission.is_relative_to(ROOT / 'out')
composition_path = ROOT / 'out/native-coefficient-sharing-producer-inputs/root-composed-inputs.json'
composition = json.loads(composition_path.read_text())
assert composition['schema'] == 'sirius-native-coefficient-sharing-diagnostic-input-composition-v1'
assert composition['passed'] is True and composition['source_revision'] == BASE and len(composition['files']) == 77
assert composition['runtime_released'] is False and composition['adoption_admitted'] is False
admission_bytes = admission.read_bytes()
assert composition['emitted_admission'] == {'path': str(admission.relative_to(ROOT)), 'bytes': len(admission_bytes), 'sha256': hashlib.sha256(admission_bytes).hexdigest()}
assert json.loads(admission_bytes)['passed'] is True
assert not (WORK / 'source-Git.json').exists()
paths = []
for root in ['out/native-coefficient-sharing-source', 'out/native-coefficient-sharing-windows-source',
             'out/native-coefficient-sharing-producer-inputs']:
    paths.extend(p for p in (ROOT / root).rglob('*') if p.is_file())
paths += [admission, Path(__file__).resolve()]
sources = {str(path.relative_to(ROOT)): path for path in paths}
sources['.github/workflows/native-coefficient-sharing-diagnostic.yml'] = ROOT / 'out/native-coefficient-sharing-windows-source/native-coefficient-sharing-diagnostic.yml'
index = WORK / 'source.index'
assert not index.exists()
env = dict(os.environ, GIT_INDEX_FILE=str(index))
records = {}
try:
    run('git', 'read-tree', BASE, env=env)
    rows = []
    for name, path in sorted(sources.items()):
        assert not path.is_symlink()
        data = path.read_bytes()
        blob = run('git', 'hash-object', '-w', '--stdin', data=data).decode().strip()
        mode = '100755' if path.stat().st_mode & 0o111 else '100644'
        rows.append(f'{mode} {blob}\t{name}\0'.encode())
        records[name] = {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest(), 'blob': blob, 'mode': mode}
    run('git', 'update-index', '-z', '--index-info', data=b''.join(rows), env=env)
    tree = run('git', 'write-tree', env=env).decode().strip()
    commit = run('git', 'commit-tree', tree, '-p', BASE,
                 data=b'Prepare matched native coefficient-call sharing producer\n').decode().strip()
    delta = run('git', 'diff-tree', '--no-commit-id', '--name-only', '-r', commit).decode().splitlines()
    assert set(delta) == set(sources)
    for name, identity in records.items():
        data = run('git', 'show', commit + ':' + name)
        assert len(data) == identity['bytes'] and hashlib.sha256(data).hexdigest() == identity['sha256']
    assert run('git', 'rev-parse', 'HEAD').decode().strip() == BASE
    assert not run('git', 'status', '--porcelain')
    run('git', 'update-ref', BRANCH, commit, '0' * 40)
    receipt = {'schema': 'sirius-native-coefficient-sharing-private-source-publication-v1', 'commit': commit,
               'base': BASE, 'branch': BRANCH, 'files': records,
               'all_canonical_files_unchanged': True, 'stored_full_bytes_verified': True,
               'pushed': False, 'compile_only': True, 'runtime_admitted': False}
    (WORK / 'source-Git.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({key: value for key, value in receipt.items() if key != 'files'}))
finally:
    index.unlink(missing_ok=True)
