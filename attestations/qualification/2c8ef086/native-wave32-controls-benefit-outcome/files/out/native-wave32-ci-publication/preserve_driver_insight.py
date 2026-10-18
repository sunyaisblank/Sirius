import hashlib
import json
import os
import shutil
import stat
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
WORK = Path(__file__).resolve().parent
PARENT = 'b8ecca0ff4d7fe399db5f37922416f589e822e98'
SOURCE = '2c8ef08693d46f57d5dea14efdbddd4874dac684'
BRANCH = 'refs/heads/development/native-fp32-packed-scratch-runtime-2026-10-18'
TAG = 'refs/tags/evidence/native-wave32-driver-insight-2c8ef086'
BUNDLE = ROOT / 'attestations/qualification/2c8ef086/native-wave32-driver-insight'


def run(*argv, data=None, env=None):
    return subprocess.run(argv, cwd=ROOT, input=data, stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE, check=True, env=env).stdout


def identity(path):
    data = path.read_bytes()
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest(),
            'mode': format(stat.S_IMODE(path.stat().st_mode), '04o')}


assert run('git', 'rev-parse', 'HEAD').decode().strip() == SOURCE
assert not run('git', 'status', '--porcelain')
assert run('git', 'rev-parse', BRANCH).decode().strip() == PARENT
assert not BUNDLE.exists()
source_root = ROOT / 'out/native-shader-driver-introspection'
sources = sorted(p for p in source_root.rglob('*') if p.is_file())
sources += [ROOT / 'out/native-constant-preservation/root-platform-frame-Git.json',
            ROOT / 'out/native-constant-preservation/independent-platform-frame-preservation.json',
            WORK / 'source-Git.json', WORK / 'runs-first.json', Path(__file__).resolve()]
records = {}
for source in sources:
    assert not source.is_symlink()
    relative = source.relative_to(ROOT)
    destination = BUNDLE / 'files' / relative
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)
    expected = identity(source)
    assert identity(destination) == expected
    records[str(relative)] = dict(expected, payload=str(destination.relative_to(BUNDLE)))
manifest = {
    'schema': 'sirius-native-wave32-driver-insight-preservation-v1',
    'source_revision': SOURCE, 'payload_files': records,
    'owned_roots': {'out/native-shader-driver-introspection': sorted(str(p.relative_to(ROOT)) for p in source_root.rglob('*') if p.is_file())},
    'scope': 'Two bounded no-dispatch diagnostic pipeline queries, raw statistics/code objects/disassembly, sources, reviews and root findings; no runtime benefit or scientific/release qualification.',
    'canonical_source': SOURCE,
    'private_candidate_source': '5d993d12b3fde7a65eaa5f7988b73b86b3a9a291',
    'external_payload_references': {
        'same_current_Transport_SPIRV_and_main_header': 'fd530d957d2e0737f5453d88d78d979f8ed17607:attestations/qualification/5b2df58b/native-constant-controls-benefit-pass/files/out/native-constant-producer-inputs/candidate',
        'actual_staged_SDK_loader': 'b8ecca0ff4d7fe399db5f37922416f589e822e98:attestations/qualification/2c8ef086/native-constant-platform-builds-and-frame/MANIFEST.json archive_backed_stage',
        'actual_host_AMD_driver': 'Host dependency; actual loaded path/size/SHA/version preserved in each query result, not a complete host backup.'},
    'numerical_performance_or_adoption_claimed': False}
(BUNDLE / 'MANIFEST.json').write_text(json.dumps(manifest, indent=2) + '\n')
for source in sources:
    assert identity(source) == {k: records[str(source.relative_to(ROOT))][k] for k in ['bytes', 'sha256', 'mode']}
index = WORK / 'driver-insight.index'
assert not index.exists()
env = dict(os.environ, GIT_INDEX_FILE=str(index))
stored = {}
try:
    run('git', 'read-tree', PARENT, env=env)
    rows = []
    for path in sorted(p for p in BUNDLE.rglob('*') if p.is_file()):
        data = path.read_bytes()
        blob = run('git', 'hash-object', '-w', '--stdin', data=data).decode().strip()
        name = str(path.relative_to(ROOT))
        mode = '100755' if path.stat().st_mode & 0o111 else '100644'
        rows.append(f'{mode} {blob}\t{name}\0'.encode())
        stored[name] = dict(identity(path), Git_mode=mode, blob=blob)
    run('git', 'update-index', '-z', '--index-info', data=b''.join(rows), env=env)
    tree = run('git', 'write-tree', env=env).decode().strip()
    commit = run('git', 'commit-tree', tree, '-p', PARENT,
                 data=b'Preserve actual driver wave64 and wave32 shader insight\n').decode().strip()
    changed = run('git', 'diff-tree', '--no-commit-id', '--name-only', '-r', commit).decode().splitlines()
    assert set(changed) == set(stored)
    for name, row in stored.items():
        data = run('git', 'show', commit + ':' + name)
        assert len(data) == row['bytes'] and hashlib.sha256(data).hexdigest() == row['sha256']
    run('git', 'update-ref', BRANCH, commit, PARENT)
    run('git', 'update-ref', TAG, commit, '0' * 40)
    receipt = {'schema': 'sirius-native-wave32-driver-insight-Git-v1', 'commit': commit,
               'parent': PARENT, 'branch': BRANCH, 'tag': TAG, 'bundle': str(BUNDLE.relative_to(ROOT)),
               'payload_files': len(records), 'payload_bytes': sum(row['bytes'] for row in records.values()),
               'stored_files': stored, 'stored_full_bytes_verified': True,
               'all_source_copies_and_modes_equal': True, 'pushed': False}
    (WORK / 'driver-insight-Git.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({key: value for key, value in receipt.items() if key != 'stored_files'}))
finally:
    index.unlink(missing_ok=True)
